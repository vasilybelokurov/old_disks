"""Can the existing samples measure v/sigma versus [Fe/H]? Mocks at the real positions and velocity errors.

For each galaxy (union members with v and [Fe/H], data/inventory.csv), split at the median [Fe/H] into a metal-poor (MP)
and a metal-rich (MR) half. Mock velocities from olddisks.kinematics with sigma = LVDB sigma_v (10 km/s if missing),
kinematic axis = LVDB position angle, and
  'equal': V/sigma = 0.75 in both halves;
  'trend': V/sigma = 0.5 (MP) and 1.0 (MR);
  'null' : V = 0 in both halves (sets the bias floor of the amplitude estimator).
Each mock is refitted per half (maximum likelihood). Power = fraction of 'trend' mocks whose recovered
Delta = (V/sigma)_MR - (V/sigma)_MP exceeds the 95th percentile of Delta in 'equal' mocks.

Outputs: results/feasibility.csv, figures/feasibility.png
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from olddisks.geometry import tangent_plane
from olddisks.kinematics import fit_rotation, mock_velocities
from olddisks.paths import DATA, LVDB

GALAXIES = ["NGC 6822", "IC 1613", "Phoenix", "Leo A", "Pegasus dIrr", "Andromeda VII", "Andromeda V", "Aquarius",
            "Cetus", "Tucana", "UGC 4879", "Sagittarius dIrr", "Andromeda I", "Sculptor"]
SCENARIOS = {"null": (0.0, 0.0), "equal": (0.75, 0.75), "trend": (0.5, 1.0)}
NMOCK = 300


def halves(s, g):
    xi, eta = tangent_plane(s.ra.to_numpy(), s.dec.to_numpy(), g.ra, g.dec)
    s = s.assign(phi=np.arctan2(xi, eta))
    med = s.feh.median()
    return s[s.feh < med], s[s.feh >= med]


def run(galaxy, s, g, rng):
    sigma = 10.0 if np.isnan(g.vlos_sigma) else g.vlos_sigma
    theta = np.radians(np.nan_to_num(g.position_angle))
    rows = []
    for scen, ratios in SCENARIOS.items():
        for k in range(NMOCK):
            rec = []
            for h, ratio in zip(halves(s, g), ratios):
                r, phi, ev = h.r_ell.to_numpy(), h.phi.to_numpy(), h.ev.fillna(h.ev.median()).to_numpy()
                v = mock_velocities(r, phi, ev, g.vlos_systemic, ratio * sigma, theta, sigma, rng)
                f = fit_rotation(r, phi, v, ev)
                rec.append(f["vrot"] / f["sigma"])
            rows.append(dict(galaxy=galaxy, scenario=scen, mock=k, vs_mp=rec[0], vs_mr=rec[1]))
    return rows


def main():
    rng = np.random.default_rng(42)
    lvdb = pd.read_csv(LVDB).drop_duplicates("name").set_index("name")
    m = pd.read_parquet(DATA / "inventory_stars.parquet")
    m = m[m.in_union & np.isfinite(m.feh)]
    mocks = pd.DataFrame([row for gal in GALAXIES for row in run(gal, m[m.galaxy == gal], lvdb.loc[gal], rng)])
    mocks["delta"] = mocks.vs_mr - mocks.vs_mp
    out = []
    for gal in GALAXIES:
        s, g = m[m.galaxy == gal], lvdb.loc[gal]
        q = {sc: mocks[(mocks.galaxy == gal) & (mocks.scenario == sc)] for sc in SCENARIOS}
        thr = q["equal"].delta.quantile(0.95)
        out.append(dict(galaxy=gal, n=len(s), n_half=len(s) // 2, med_ev=s.ev.median(), sigma_v=g.vlos_sigma,
                        r_ell_p90=s.r_ell.quantile(0.9), feh_mp=halves(s, g)[0].feh.median(),
                        feh_mr=halves(s, g)[1].feh.median(),
                        null_vs_med=q["null"].vs_mp.median(),
                        mp_vs_med=q["trend"].vs_mp.median(), mp_vs_sd=q["trend"].vs_mp.std(),
                        mr_vs_med=q["trend"].vs_mr.median(), mr_vs_sd=q["trend"].vs_mr.std(),
                        power=(q["trend"].delta > thr).mean()))
    res = pd.DataFrame(out).round(3)
    res.to_csv("results/feasibility.csv", index=False)
    print(res.to_string(index=False))

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
    y = np.arange(len(GALAXIES))
    for sc, col, off in (("trend", "C0", -0.15), ("null", "0.5", 0.15)):
        for half, mk, truth in (("vs_mp", "o", 0.5), ("vs_mr", "s", 1.0)):
            p = mocks[mocks.scenario == sc].groupby("galaxy")[half].quantile([0.16, 0.5, 0.84]).unstack().reindex(GALAXIES)
            dy = off + (0.07 if half == "vs_mr" else -0.07)
            ax[0].errorbar(p[0.5], y + dy, xerr=[p[0.5] - p[0.16], p[0.84] - p[0.5]], fmt=mk, color=col, ms=4,
                           mfc="white" if half == "vs_mp" else col, lw=1,
                           label=f"{sc}: {'MP' if half == 'vs_mp' else 'MR'} (truth {truth if sc == 'trend' else 0})")
    ax[0].axvline(0.5, color="C0", ls=":", lw=0.8), ax[0].axvline(1.0, color="C0", ls="--", lw=0.8)
    ax[0].set_yticks(y, [f"{g} (N={n})" for g, n in zip(res.galaxy, res.n)], fontsize=8), ax[0].invert_yaxis()
    ax[0].set_xlabel("recovered V/σ (median, 16–84%) per [Fe/H] half"), ax[0].legend(fontsize=6, loc="lower right")
    ax[1].barh(y, res.power, color="C0")
    ax[1].axvline(0.8, color="k", ls=":", lw=0.8)
    ax[1].set_yticks(y, []), ax[1].invert_yaxis(), ax[1].set_xlim(0, 1)
    ax[1].set_xlabel("power to detect ΔV/σ = 0.5 (MP 0.5 → MR 1.0) at 95%")
    fig.tight_layout()
    fig.savefig("figures/feasibility.png", dpi=120)


if __name__ == "__main__":
    main()
