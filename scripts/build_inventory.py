"""Per-galaxy inventory of stars with line-of-sight velocity and [Fe/H] from published catalogues and WSDB surveys.

Member: the paper's flag where given (member_paper); otherwise |v - v_sys| < 3 sqrt(sigma^2 + e_v^2) (sigma = 10 km/s if
LVDB has none). Elliptical radius R_ell in units of the LVDB R_half.
Union over sources: a star within 1 arcsec of one already counted (in another source or the same one) is a duplicate; sources are added in decreasing order
of members with [Fe/H].

Outputs: data/inventory.csv (one row per galaxy and source, plus source = "UNION"),
         ~/data/old_disks/inventory_stars.parquet (all members, with `in_union`).
"""
import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord

from olddisks.catalogues import load
from olddisks.geometry import elliptical_radius
from olddisks.paths import DATA, LVDB, SURVEY_MEMBERS


def survey_stars():
    s = pd.read_parquet(SURVEY_MEMBERS)
    s = s.rename(columns={"id": "star_id"}).assign(source=s.survey, catalogue="WSDB " + s.survey, instrument=s.survey,
                                                   efeh=np.nan, member_paper=pd.NA)
    return s[["source", "catalogue", "galaxy", "instrument", "star_id", "ra", "dec", "v", "ev", "feh", "efeh",
              "member_paper"]]


def members(stars, lvdb):
    g = lvdb.loc[stars.galaxy].reset_index(drop=True)
    s = stars.reset_index(drop=True)
    ell, pa = g.ellipticity.fillna(0).to_numpy(), g.position_angle.fillna(0).to_numpy()
    s["r_ell"] = elliptical_radius(s.ra.to_numpy(), s.dec.to_numpy(), g.ra.to_numpy(), g.dec.to_numpy(), pa, ell) \
        * 60 / g.rhalf.to_numpy()
    sig = g.vlos_sigma.fillna(10).to_numpy()
    window = np.abs(s.v - g.vlos_systemic.to_numpy()) < 3 * np.sqrt(sig**2 + s.ev.fillna(0) ** 2)
    s["member"] = s.member_paper.fillna(pd.Series(window, dtype="boolean")).astype(bool)
    return s[s.member & np.isfinite(s.v)]


def union(m):
    """Flag in_union: first occurrence of each star (1 arcsec) when sources are taken in order of n with [Fe/H]."""
    m = m.copy()
    m["in_union"] = False
    for gal, s in m.groupby("galaxy"):
        order = s.assign(f=np.isfinite(s.feh)).groupby("source").f.sum().sort_values(ascending=False).index
        kept_ra, kept_dec = np.empty(0), np.empty(0)
        for src in order:
            idx = s.index[s.source == src]
            ra, dec = m.loc[idx, "ra"].to_numpy(), m.loc[idx, "dec"].to_numpy()
            new = np.ones(len(idx), bool)
            c = SkyCoord(ra * u.deg, dec * u.deg)
            i, j, _, _ = c.search_around_sky(c, 1 * u.arcsec)
            new[np.unique(j[i < j])] = False  # repeats within the source
            if len(kept_ra):
                sep = SkyCoord(ra * u.deg, dec * u.deg).match_to_catalog_sky(SkyCoord(kept_ra * u.deg, kept_dec * u.deg))[1]
                new = sep.arcsec > 1.0
            m.loc[idx[new], "in_union"] = True
            kept_ra, kept_dec = np.r_[kept_ra, ra[new]], np.r_[kept_dec, dec[new]]
    return m


def summary(m):
    agg = dict(n_v=("v", "size"), n_v_feh=("has_feh", "sum"), n_beyond_2rh=("far", "sum"),
               max_r_ell=("r_ell", "max"), med_ev=("ev", "median"), med_efeh=("efeh", "median"))
    m = m.assign(has_feh=np.isfinite(m.feh), far=m.r_ell > 2)
    per_src = m.groupby(["galaxy", "source", "instrument"]).agg(**agg).reset_index()
    tot = m[m.in_union].groupby("galaxy").agg(**agg).reset_index().assign(source="UNION", instrument="")
    return pd.concat([tot, per_src]).round(2)


def main():
    lvdb = pd.read_csv(LVDB).drop_duplicates("name").set_index("name")
    stars = pd.concat([load(), survey_stars()], ignore_index=True)
    m = union(members(stars, lvdb))
    m.to_parquet(DATA / "inventory_stars.parquet")
    inv = summary(m).merge(lvdb[["host", "M_V", "distance"]], left_on="galaxy", right_index=True)
    inv = inv.sort_values(["M_V", "galaxy", "n_v_feh"], ascending=[True, True, False])
    inv.to_csv("data/inventory.csv", index=False)
    print(inv[inv.source == "UNION"].sort_values("n_v_feh", ascending=False).to_string(index=False))


if __name__ == "__main__":
    main()
