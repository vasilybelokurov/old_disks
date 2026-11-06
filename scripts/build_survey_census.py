"""Count velocity-selected candidate members of Local Group dwarfs in public WSDB spectroscopic surveys.

Candidate member: elliptical radius < 3 R_half and |v - v_sys| < 3 sqrt(sigma^2 + e_v^2)
(sigma = 10 km/s where LVDB has none). No proper-motion or CMD cut: counts are upper limits
for galaxies whose v_sys is near the Milky Way foreground.

Output: data/survey_census.csv, one row per (galaxy, survey) with at least one candidate,
and all stars inside 3 R_half (column `candidate`) in ~/data/old_disks/survey_members.parquet.
"""
import numpy as np
import pandas as pd
import sqlutilpy as sqlutil

from olddisks.paths import LVDB, SURVEY_MEMBERS
from olddisks.geometry import elliptical_radius

EXCLUDE = ["LMC", "SMC", "Sagittarius", "M 33"]  # LMC/SMC: sibling projects; Sgr: tidal stream; M33: not a dwarf
NRHALF = 3

# sdss5_dr20.allstar (APOGEE) is omitted: the WSDB copy has 390k rows, 258 at dec < -20, none in any dwarf.
# Each query returns id, ra, dec, v, ev, feh (NaN where absent); {cone} is the q3c predicate.
SURVEYS = {
    "apogee_dr17": """select apogee_id as id, ra, dec, vhelio_avg as v, verr as ev, fe_h as feh
        from apogee_dr17.allstar where {cone} and (extratarg & 16) = 0 and vhelio_avg is not null""",
    "sdss5_dr20_boss": """select sdss_id as id, ra, dec, v_rad as v, e_v_rad as ev, fe_h as feh
        from sdss5_dr20.astraallstarbossnet_boss where {cone} and v_rad is not null""",
    "desi_dr1": """select targetid as id, target_ra as ra, target_dec as dec, vrad as v, vrad_err as ev, feh
        from desi_dr1.stellar_rvs where {cone} and is_primary = 1 and rvs_warn = 0 and success = 1""",
}


def galaxies():
    d = pd.read_csv(LVDB)
    keep = (d.confirmed_galaxy == 1) & (d.distance < 1500) & (d.M_V < -7) & d.vlos_systemic.notna()
    return d[keep & ~d.name.isin(EXCLUDE)].reset_index(drop=True)


def members(g, survey):
    cone = f"q3c_radial_query(ra, dec, {g.ra}, {g.dec}, {NRHALF * g.rhalf / 60})"
    if survey == "desi_dr1":
        cone = cone.replace("(ra, dec", "(target_ra, target_dec")
    s = pd.DataFrame(sqlutil.get(SURVEYS[survey].format(cone=cone), asDict=True))
    if s.empty:
        return s
    s["id"] = s["id"].astype(str)
    s = s.drop_duplicates("id")
    ell = 0.0 if np.isnan(g.ellipticity) else g.ellipticity
    pa = 0.0 if np.isnan(g.position_angle) else g.position_angle
    s["r_ell"] = elliptical_radius(s.ra.values, s.dec.values, g.ra, g.dec, pa, ell) * 60 / g.rhalf
    sig = 10.0 if np.isnan(g.vlos_sigma) else g.vlos_sigma
    s = s[s.r_ell < NRHALF].copy()
    ev = s.ev.astype(float).fillna(0)
    s["candidate"] = np.abs(s.v - g.vlos_systemic) < 3 * np.sqrt(sig**2 + ev**2)
    s["galaxy"], s["survey"] = g["name"], survey
    return s


def main():
    gal = galaxies()
    out = pd.concat([members(g, sv) for _, g in gal.iterrows() for sv in SURVEYS], ignore_index=True)
    out["feh"] = out.feh.astype(float)
    out.to_parquet(SURVEY_MEMBERS)
    cand = out[out.candidate]
    census = (
        cand.assign(has_feh=np.isfinite(cand.feh))
        .groupby(["galaxy", "survey"])
        .agg(n_v=("v", "size"), n_v_feh=("has_feh", "sum"), med_ev=("ev", "median"), max_r_ell=("r_ell", "max"))
        .reset_index()
    )
    census = census.merge(gal[["name", "host", "M_V", "distance", "vlos_systemic", "vlos_sigma"]],
                          left_on="galaxy", right_on="name").drop(columns="name")
    census.round(2).sort_values(["M_V", "survey"]).to_csv("data/survey_census.csv", index=False)
    print(f"{len(gal)} galaxies queried; {len(census)} (galaxy, survey) pairs with candidates")
    print(census.sort_values("n_v", ascending=False).to_string(index=False))


if __name__ == "__main__":
    main()
