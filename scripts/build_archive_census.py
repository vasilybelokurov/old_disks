"""Census of archival spectroscopy of Local Group dwarfs in the ESO archive and the Keck Observatory Archive (KOA).

Footprint: a box of half-size min(3 R_half, 2 deg) around each LVDB centre (box, not circle: ESO raw table rejects
circles because some rows carry invalid coordinates). Counts science frames, programmes and exposure hours per
(galaxy, archive, instrument), and lists programme IDs and PIs so the data can be traced to publications.
For "ESO phase3", n_frames counts reduced products: one per star for GIRAFFE/UVES fibres, so `hours` there is summed
over stars, not telescope time.

Output: data/archive_census.csv
"""
import numpy as np
import pandas as pd
import pyvo

from olddisks.paths import LVDB

ESO = pyvo.dal.TAPService("https://archive.eso.org/tap_obs")
KOA = pyvo.dal.TAPService("https://koa.ipac.caltech.edu/TAP")
ESO_SPECTROGRAPHS = "'GIRAFFE','UVES','FORS2','FORS1','MUSE','XSHOOTER','VIMOS','FLAMES'"
EXCLUDE = ["LMC", "SMC", "Sagittarius", "M 33"]


def galaxies():
    d = pd.read_csv(LVDB)
    keep = (d.confirmed_galaxy == 1) & (d.distance < 1500) & (d.M_V < -7)
    return d[keep & ~d.name.isin(EXCLUDE)].reset_index(drop=True)


def box(g, ra="ra", dec="dec"):
    h = min(3 * g.rhalf / 60, 2.0)
    dra = h / np.cos(np.radians(g.dec))
    return f"{dec} between {g.dec - h} and {g.dec + h} and {ra} between {g.ra - dra} and {g.ra + dra}"


def summarise(t, g, archive):
    if len(t) == 0:
        return []
    t = t.to_pandas()
    out = []
    for inst, s in t.groupby("instrument"):
        progs = s.groupby("prog").agg(pi=("piname", "first"), hr=("texp", lambda x: x.sum() / 3600)).sort_values("hr")[::-1]
        out.append(dict(galaxy=g["name"], host=g.host, M_V=g.M_V, distance=g.distance, archive=archive, instrument=inst,
                        n_frames=len(s), n_prog=len(progs), hours=round(s.texp.sum() / 3600, 1),
                        last_release=str(s.rel.max())[:10],
                        programmes="; ".join(f"{p} ({r.pi}, {r.hr:.1f} h)" for p, r in progs.iterrows())))
    return out


def eso(g):
    q = f"""select instrument, prog_id as prog, pi_coi as piname, exposure as texp, release_date as rel
        from dbo.raw where dp_cat = 'SCIENCE' and instrument in ({ESO_SPECTROGRAPHS}) and {box(g)}"""
    return summarise(ESO.search(q, maxrec=200000).to_table(), g, "ESO raw")


def eso_phase3(g):
    q = f"""select instrument_name as instrument, proposal_id as prog, obs_creator_name as piname, t_exptime as texp,
        obs_release_date  as rel from ivoa.ObsCore
        where dataproduct_type in ('spectrum', 'cube') and instrument_name in ({ESO_SPECTROGRAPHS})
        and {box(g, 's_ra', 's_dec')}"""
    return summarise(ESO.search(q, maxrec=200000).to_table(), g, "ESO phase3")


def koa(g):
    q = f"""select progid as prog, progpi as piname, elaptime as texp, date_obs, propint
        from koa_deimos where koaimtyp = 'object' and obsmode in ('mos', 'longslit') and {box(g)}"""
    t = KOA.search(q, maxrec=200000).to_table()
    if len(t):
        t["instrument"] = "DEIMOS"
        # propint is the proprietary period in months
        t["rel"] = [str(pd.Timestamp(str(d)) + pd.DateOffset(months=int(m))) for d, m in zip(t["date_obs"], t["propint"])]
    return summarise(t, g, "KOA")


def main():
    rows = []
    for _, g in galaxies().iterrows():
        for f in (eso, eso_phase3, koa):
            try:
                rows += f(g)
            except Exception as e:  # report and continue: one failed service must not lose the rest
                print(f"FAILED {f.__name__} {g['name']}: {str(e)[:200]}")
    out = pd.DataFrame(rows).sort_values(["M_V", "archive", "hours"], ascending=[True, True, False])
    out.to_csv("data/archive_census.csv", index=False)
    print(out.drop(columns="programmes").to_string(index=False))


if __name__ == "__main__":
    main()
