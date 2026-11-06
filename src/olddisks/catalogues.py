"""Harmonise published per-star catalogues (VizieR copies in DATA/vizier) into one table.

Columns of `load()`: source, catalogue, galaxy (LVDB name), star_id, instrument, ra, dec, v, ev [km/s], feh, efeh [dex],
member_paper (True/False from the paper's flag or probability > 0.5; None where the paper gives none).
"""
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table

from olddisks.paths import DATA, LVDB

ROMAN = {1: "I", 2: "II", 3: "III", 4: "IV", 5: "V", 6: "VI", 7: "VII", 8: "VIII", 9: "IX", 10: "X", 11: "XI",
         12: "XII", 13: "XIII", 14: "XIV", 15: "XV", 16: "XVI", 17: "XVII", 18: "XVIII", 19: "XIX", 20: "XX",
         21: "XXI", 22: "XXII", 23: "XXIII", 24: "XXIV", 25: "XXV"}
KIRBY_NAMES = {"NGC6822": "NGC 6822", "IC1613": "IC 1613", "Cetus": "Cetus", "Pegasus": "Pegasus dIrr",
               "Peg": "Pegasus dIrr", "VV124": "UGC 4879", "LeoA": "Leo A", "Aquarius": "Aquarius",
               "SagDIG": "Sagittarius dIrr"}
MMFS_NAMES = {"Car": "Carina", "For": "Fornax", "Scl": "Sculptor", "Sex": "Sextans"}


def table(name):
    cat = name.rsplit("/", 1)[0]
    return Table.read(DATA / "vizier" / cat.replace("/", "_") / (name.replace("/", "_") + ".ecsv")).to_pandas()


def walker_name(target, lvdb_names):
    """'Leo_1' -> 'Leo I', 'Carina_1' -> 'Carina', 'Antlia_2' -> 'Antlia II'; None if not in LVDB."""
    base, n = target.rsplit("_", 1)
    base = base.replace("_", " ")
    for cand in (f"{base} {ROMAN[int(n)]}", base if n == "1" else None):
        if cand in lvdb_names:
            return cand
    return None


def frame(source, catalogue, t, galaxy, v, ev, feh=None, efeh=None, member=None, star_id=None, instrument=""):
    """Build the standard columns from table t; galaxy is a name or an array, member a boolean array or None."""
    n = len(t)
    col = lambda c: t[c].to_numpy(float) if c else np.full(n, np.nan)
    ra = t["RAJ2000"] if "RAJ2000" in t else t["RA_ICRS"]
    dec = t["DEJ2000"] if "DEJ2000" in t else t["DE_ICRS"]
    if not pd.api.types.is_numeric_dtype(ra):  # sexagesimal strings
        c = SkyCoord(ra.astype(str), dec.astype(str), unit=(u.hourangle, u.deg))
        ra, dec = pd.Series(c.ra.deg), pd.Series(c.dec.deg)
    return pd.DataFrame(dict(
        source=source, catalogue=catalogue, galaxy=galaxy, instrument=instrument,
        star_id=t[star_id].astype(str).to_numpy() if star_id else t["recno"].astype(str).to_numpy(),
        ra=ra.to_numpy(float), dec=dec.to_numpy(float), v=col(v), ev=col(ev), feh=col(feh), efeh=col(efeh),
        member_paper=pd.array(member, dtype="boolean") if member is not None else pd.array([pd.NA] * n, dtype="boolean")))


def kirby2014():
    """Velocities from Kirby+2014 joined to [Fe/H] from Kirby+2013 table3 by galaxy and star name."""
    v = table("J/MNRAS/439/1015/table2")
    f = table("J/ApJ/779/102/table3")
    v["gal"] = v.Galaxy.map(KIRBY_NAMES)
    f["gal"] = f.Gal.map(KIRBY_NAMES)
    t = v.merge(f[["gal", "Name", "[Fe/H]", "e_[Fe/H]"]].dropna(subset=["gal"]), on=["gal", "Name"], how="left")
    return frame("Kirby+2014", "J/MNRAS/439/1015 + J/ApJ/779/102", t, t.gal, "HRV", "e_HRV", "[Fe/H]", "e_[Fe/H]",
                 t.Mmb == "Y", "Name", "DEIMOS")


def kirby2017():
    v = table("J/ApJ/834/9/table2")
    f = table("J/ApJ/834/9/table3")
    t = v.merge(f[["Gal", "Name", "[Fe/H]", "e_[Fe/H]"]], on=["Gal", "Name"], how="left")
    return frame("Kirby+2017", "J/ApJ/834/9", t, t.Gal.map(KIRBY_NAMES), "HRV", "e_HRV", "[Fe/H]", "e_[Fe/H]",
                 t.Mm == "Y", "Name", "DEIMOS")


def walker2023(lvdb_names):
    """One row per star: the per-star averages (VlosAvg, FeHAvg) repeated on every epoch row."""
    out = []
    for tab, inst in (("hectocat", "Hectochelle"), ("m2fshi", "M2FS HiRes"), ("m2fsmed", "M2FS MedRes")):
        t = table(f"J/ApJS/268/19/{tab}")
        t = t[t.goodNobs > 0].drop_duplicates("ObjID") if "goodNobs" in t else t.drop_duplicates("ObjID")
        t = t.assign(gal=[walker_name(x, lvdb_names) for x in t.Target]).dropna(subset=["gal"])
        t = t.reset_index(drop=True)
        out.append(frame("Walker+2023", "J/ApJS/268/19", t, t.gal.to_numpy(), "VlosAvg", "e_VlosAvg", "FeHAvg",
                         "e_FeHAvg", None, "ObjID", inst)
                   .assign(pmra=t.pmRA.to_numpy(float), pmdec=t.pmDE.to_numpy(float),
                           epmra=t.e_pmRA.to_numpy(float), epmdec=t.e_pmDE.to_numpy(float)))
    return pd.concat(out)


def ho2015():
    """Ho+2015 CaT [Fe/H] (members) joined within 1 arcsec to Geha+2010 (NGC 147, 185) and Geha+2006 (NGC 205)
    velocities; the nearest Ho star wins when two match one velocity star."""
    h = table("J/ApJ/798/77/table2")
    g10, g06 = table("J/ApJ/711/361/veloc"), table("J/AJ/131/332/table3")
    out = []
    for code, gal, v, vid in (("N147", "NGC 147", g10[g10.NGC == 147], "Name"),
                              ("N185", "NGC 185", g10[g10.NGC == 185], "Name"),
                              ("N205", "NGC 205", g06, "[GGR2006]")):
        hh, v = h[h.Gal == code].reset_index(drop=True), v.reset_index(drop=True)
        ch = SkyCoord(hh.RAJ2000.astype(str), hh.DEJ2000.astype(str), unit=(u.hourangle, u.deg))
        cv = SkyCoord(v.RAJ2000.astype(str), v.DEJ2000.astype(str), unit=(u.hourangle, u.deg))
        idx, sep, _ = ch.match_to_catalog_sky(cv)
        m = pd.DataFrame(dict(i=np.arange(len(hh)), j=idx, sep=sep.arcsec))
        m = m[m.sep < 1].sort_values("sep").drop_duplicates("j")
        t = pd.concat([v.loc[m.j].reset_index(drop=True)[["RAJ2000", "DEJ2000", "HRV", "e_HRV", vid]],
                       hh.loc[m.i].reset_index(drop=True)[["[Fe/H]", "e_[Fe/H]"]]], axis=1)
        out.append(frame("Ho+2015", "J/ApJ/798/77 + J/ApJ/711/361, J/AJ/131/332", t, gal, "HRV", "e_HRV", "[Fe/H]",
                         "e_[Fe/H]", np.ones(len(t), bool), vid, "DEIMOS"))
    return pd.concat(out)


def load():
    lvdb_names = set(pd.read_csv(LVDB).name)
    t = table
    parts = [
        frame("Swan+2016", "J/MNRAS/456/4315", t("J/MNRAS/456/4315/tablea2"), "NGC 6822", "RV", "e_RV", "[Fe/H]",
              "e_[Fe/H]", None, "Star", "FORS2"),
        kirby2014(),
        kirby2017(),
        ho2015(),
        frame("Hermosa Munoz+2020", "J/A+A/634/A10", x := t("J/A+A/634/A10/tablea2"), "Aquarius", "HV", "e_HV",
              "[Fe/H]", "e_[Fe/H]", x.Memb == "Y", None, "FORS2"),
        frame("Taibi+2018", "J/A+A/618/A122", x := t("J/A+A/618/A122/table5"), "Cetus", "HV", "e_HV", "[Fe/H]",
              "e_[Fe/H]", x.Memb == "Y", None, "FORS2"),
        frame("Taibi+2020", "J/A+A/635/A152", x := t("J/A+A/635/A152/tablec1"), "Tucana", "HV", "e_HV", "[Fe/H]",
              "e_[Fe/H]", x.P1 > 0.5, None, "FORS2"),
        frame("Taibi+2020", "J/A+A/635/A152", x := t("J/A+A/635/A152/tablec2"), "Tucana", "HV", "e_HV", "[Fe/H]",
              "e_[Fe/H]", x.P1 > 0.5, None, "FORS2"),
        frame("Taibi+2020", "J/A+A/635/A152", x := t("J/A+A/635/A152/tablec3"), "Tucana", "HV", "e_HV", None, None,
              x.P1 > 0.5, None, "GIRAFFE"),
        frame("Kacharov+2017", "J/MNRAS/466/2006", x := t("J/MNRAS/466/2006/phxcat"), "Phoenix", "HRV", "e_HRV",
              "[Fe/H]", "e_[Fe/H]", x.Memb == "yes", "Name", "FORS2"),
        frame("Tolstoy+2023", "J/A+A/675/A49", x := t("J/A+A/675/A49/tablee1").pipe(  # q=1: low-S/N [Fe/H], unused
                  lambda d: d.assign(**{"[Fe/H]": d["[Fe/H]"].where(d["q_[Fe/H]"] == 0)})), "Sculptor", "Vlos", "e_Vlos",
              "[Fe/H]", "e_[Fe/H]", x.Mem == "m", "GaiaDR3", "GIRAFFE"),
        frame("Walker+2009", "J/AJ/137/3100", x := t("J/AJ/137/3100/stars"), x.Target.str[:3].map(MMFS_NAMES),
              "<HV>", "e_<HV>", None, None, x.Mmb > 0.5, "Target", "MMFS"),
        frame("Pace+2021", "J/ApJ/923/77", x := t("J/ApJ/923/77/table3"), "Fornax", "VLOS", "e_VLOS", "Z", "e_Z",
              x.Pdsph > 0.5, "Gaia", "M2FS+MMFS"),
        frame("Battaglia+2011", "J/MNRAS/411/1013", t("J/MNRAS/411/1013/table3"), "Sextans", "HV", "e_HV", "[Fe/H]",
              "e_[Fe/H]", None, "Seq", "GIRAFFE"),
        frame("Tolstoy+2025", "J/A+A/698/A53", x := t("J/A+A/698/A53/tabled1"), "Sextans", "Vel", "e_Vel",
              "[Fe/H]", "e_[Fe/H]", x.Pmem > 0.5, "GaiaDR3", "GIRAFFE"),
        frame("Walker+2015", "J/MNRAS/448/2717", t("J/MNRAS/448/2717/table5"), "Draco", "vlos", "e_vlos", "[Fe/H]",
              "e_[Fe/H]", None, "All", "Hectochelle"),
        frame("Pace+2020", "J/MNRAS/495/3022", x := t("J/MNRAS/495/3022/table2"), "Ursa Minor", "vLOS", "e_vLOS",
              "[Fe/H]", "e_[Fe/H]", x.PdSph > 0.5, "ID", "DEIMOS"),
        frame("Spencer+2017", "J/ApJ/836/202", x := t("J/ApJ/836/202/table3"), "Leo II", "RVel", "e_RVel", "[Fe/H]",
              "e_[Fe/H]", x.Mm == "Y", None, "Hectochelle"),
        frame("Ji+2021", "J/ApJ/921/32", x := t("J/ApJ/921/32/table4"), "Antlia II", "HRV", "e_HRV", "[Fe/H]1",
              "e_[Fe/H]1", x.Mm > 0.5, "Gaia", "AAOmega"),  # table4: Antlia 2
        frame("Ji+2021", "J/ApJ/921/32", x := t("J/ApJ/921/32/table5"), "Crater II", "HRV", "e_HRV", "[Fe/H]1",
              "e_[Fe/H]1", x.Mm > 0.5, "Gaia", "AAOmega"),  # table5: Crater 2
        frame("Jenkins+2021", "J/ApJ/920/92", x := t("J/ApJ/920/92/table6"), "Bootes I", "RVel", "e_RVel", None, None,
              x.Member.astype(str) == "1", "ID", "AAOmega"),
        frame("Longeard+2022", "J/MNRAS/516/2348", x := t("J/MNRAS/516/2348/table1"), "Bootes I", "RV", "e_RV",
              "[Fe/H]G", "e_[Fe/H]G", x.Member.str.startswith("Y"), None, "AAOmega"),
        walker2023(lvdb_names),
        frame("Geha+2010", "J/ApJ/711/361", x := t("J/ApJ/711/361/veloc"), "NGC " + x.NGC.astype(str), "HRV",
              "e_HRV", None, None, None, "Name", "DEIMOS"),
        frame("Tollerud+2012", "J/ApJ/752/45", x := t("J/ApJ/752/45/table2"),
              "Andromeda " + x.Sat.astype(int).map(ROMAN), "HRV", "e_HRV", None, None, x.Mmb > 0.5, None, "DEIMOS"),
        frame("Kirby+2020", "J/AJ/159/46", x := t("J/AJ/159/46/table2"), x.dSph, "HRV", "e_HRV", "[Fe/H]",
              "e_[Fe/H]", x.Mm == "Y", "Name", "DEIMOS"),
    ]
    out = pd.concat(parts, ignore_index=True)
    return out[out.galaxy.isin(lvdb_names)].reset_index(drop=True)
