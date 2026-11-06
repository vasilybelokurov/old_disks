"""Rotation (velocity gradient) versus [Fe/H] in MW dSphs, pre-registered in docs/plans/de_model.md (bottom section).

Per galaxy and catalogue (Walker+2023 baseline; flagged-membership catalogue as robustness): sample_dsph (perspective-corrected,
PM-cleaned members within 6 R_half), [Fe/H]-blind fit with the kinematic axis free and R0 from the grid, R* = median R of the
stars with [Fe/H]; T_rot and T_sigma with 10^4 within-R permutations (alpha = 0.05/11); tercile contrast of V/sigma at R*
calibrated by NMOCK mocks per slope (0 and 1 dex^-1, k = 1, 5 % contaminants) at the real positions and errors.
Tucana II (15 Walker+2023 members) is skipped.
Outputs: results/dsph_fits.csv, figures/dsph_fit.png.
"""
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from olddisks.defit import blind_fit, make_mock, sample_dsph, statistics

RUNS = [("Sculptor", "Walker+2023"), ("Sculptor", "Tolstoy+2023"), ("Draco", "Walker+2023"), ("Ursa Minor", "Walker+2023"),
        ("Ursa Minor", "Pace+2020"), ("Carina", "Walker+2023"), ("Sextans", "Walker+2023"), ("Sextans", "Tolstoy+2025"),
        ("Fornax", "Walker+2023"), ("Fornax", "Pace+2021"), ("Antlia II", "Walker+2023"), ("Antlia II", "Ji+2021"),
        ("Crater II", "Walker+2023"), ("Crater II", "Ji+2021"), ("Leo II", "Walker+2023"), ("Leo II", "Spencer+2017"),
        ("Bootes I", "Walker+2023")]
R0_GRID = [0.5, 1.0, 2.0, 4.0]
ALPHA = 0.05 / 11
NMOCK = 100


@lru_cache(maxsize=None)  # one sample and blind fit per (galaxy, catalogue) per process
def setup(gal, src):
    s = sample_dsph(gal, src)
    Rstar = float(np.median(s.R[np.isfinite(s.feh)]))
    b = min((blind_fit(s, R0, Rstar, fit_theta=True) for R0 in R0_GRID), key=lambda b: b["nll"])
    return s, b


def mock(args):
    gal, src, slope, seed = args
    s, b = setup(gal, src)
    rng = np.random.default_rng(seed)
    m = make_mock(s, dict(b, f=max(b["f"], 0.05)), slope, 1.0, rng)
    st = statistics(m, b["R0"], b["Rstar"], nperm=500, rng=rng, fit_theta=True)
    return dict(galaxy=gal, source=src, slope=slope, p_rot=st["p_rot"], delta=st["delta"])


def real(args):
    gal, src = args
    s, b = setup(gal, src)
    st = statistics(s, b["R0"], b["Rstar"], nperm=10000, rng=np.random.default_rng(7), fit_theta=True)
    bb = st.pop("blind")
    t = s[np.isfinite(s.feh)]
    res = t.feh - np.polyval(np.polyfit(t.R, t.feh, 1), t.R)
    lam = (np.var(res) - np.mean(t.efeh**2)) / np.var(res)  # fraction of the [Fe/H] variance (about the radial trend) not error
    return dict(galaxy=gal, source=src, n=len(s), n_feh=len(t), feh_med=float(t.feh.median()), efeh_med=float(t.efeh.median()),
                lam=lam,
                Rstar=bb["Rstar"], R0=bb["R0"], dtheta_deg=np.degrees(bb["dtheta"]), V=bb["V"], sigma0=bb["sigma0"],
                t_rot=st["t_rot"], p_rot=st["p_rot"], t_sig=st["t_sig"], p_sig=st["p_sig"], delta=st["delta"],
                ratio_mp=st["ratio"][0], ratio_mid=st["ratio"][1], ratio_mr=st["ratio"][2],
                z_mp=st["zmean"][0], z_mid=st["zmean"][1], z_mr=st["zmean"][2])


def main():
    jobs = [(g, s, sl, 10000 * i + 100 * j + n) for i, (g, s) in enumerate(RUNS) for j, sl in enumerate((0.0, 1.0))
            for n in range(NMOCK)]
    with ProcessPoolExecutor() as ex:
        mocks = pd.DataFrame(list(ex.map(mock, jobs, chunksize=10)))
        r = pd.DataFrame(list(ex.map(real, RUNS)))
    mocks.to_csv("results/dsph_mocks.csv", index=False)
    cal = mocks.groupby(["galaxy", "source", "slope"]).agg(d=("delta", "mean"), sd=("delta", "std"),
                                                            rate=("p_rot", lambda p: (p < ALPHA).mean())).unstack("slope")
    r = r.set_index(["galaxy", "source"])
    r["C"] = cal["d"][1.0] - cal["d"][0.0]
    r["C_err"] = np.sqrt(cal["sd"][0.0] ** 2 + cal["sd"][1.0] ** 2) / np.sqrt(NMOCK)
    r["sd0"] = cal["sd"][0.0]
    ok = r.C > 3 * r.C_err  # the tercile contrast responds to a slope: amplitude calibratable
    r["slope"] = (r.delta / r.C).where(ok)
    r["slope_sd"] = (r.sd0 / r.C).where(ok)
    r["fpr_mock"], r["power_1dex"] = cal["rate"][0.0], cal["rate"][1.0]
    r["detected_rot"] = r.p_rot < ALPHA
    r = r.reset_index()
    r.round(4).to_csv("results/dsph_fits.csv", index=False)
    pd.set_option("display.width", 250)
    print(r[["galaxy", "source", "n_feh", "feh_med", "efeh_med", "lam", "C", "C_err", "Rstar", "dtheta_deg", "V", "sigma0", "p_rot", "p_sig", "ratio_mp", "ratio_mr",
             "slope", "slope_sd", "fpr_mock", "power_1dex", "detected_rot"]].round(3).to_string(index=False))
    figure(r)


def figure(r):
    fig, ax = plt.subplots(1, 2, figsize=(14, 7), gridspec_kw=dict(width_ratios=[1.1, 1]))
    y = np.arange(len(r))[::-1]
    walker = (r.source == "Walker+2023").to_numpy()
    for sel, mk, lab in ((walker, "o", "Walker+2023 (baseline)"), (~walker, "s", "flagged-membership catalogue")):
        ax[0].errorbar(r.slope[sel], y[sel], xerr=1.96 * r.slope_sd[sel], fmt=mk, color="k", mfc="k" if mk == "o" else "white",
                       label=lab)
    for yi, (_, row) in zip(y, r.iterrows()):
        ax[0].text(3.05, yi, f"p={row.p_rot:.3f} N={row.n_feh} λ={row.lam:.2f}", va="center", fontsize=7)
        if not np.isfinite(row.slope):
            ax[0].text(0, yi, "no [Fe/H] leverage (C ≈ 0)", va="center", ha="center", fontsize=7, color="C3")
    ax[0].axvline(0, color="k", lw=0.5)
    ax[0].set_yticks(y, [f"{g} ({s.split('+')[0]})" for g, s in zip(r.galaxy, r.source)], fontsize=8)
    ax[0].set_xlim(-3, 4.2)
    ax[0].set_xlabel(r"d(V/σ)/d[Fe/H] at R* [dex$^{-1}$], 95 % (k = 1, mock-calibrated)")
    ax[0].legend(fontsize=7, loc="lower right")
    ax[0].set_title(f"velocity gradient (rotation or tides) vs [Fe/H]; α = {ALPHA:.4f}", fontsize=9)
    for i, (_, row) in enumerate(r[walker].iterrows()):
        col = plt.cm.tab20(i % 20)
        ax[1].plot([row.z_mp, row.z_mid, row.z_mr], [row.ratio_mp, row.ratio_mid, row.ratio_mr], "-o", color=col, ms=4, label=row.galaxy)
    ax[1].axhline(0, color="k", lw=0.5)
    ax[1].set_xlabel("mean observed [Fe/H] of tercile (Walker+2023)"), ax[1].set_ylabel("V/σ at R* (projected, LOS)")
    ax[1].legend(fontsize=7, ncol=2), ax[1].set_title("tercile V/σ (observed [Fe/H]: diluted)", fontsize=9)
    fig.tight_layout()
    fig.savefig("figures/dsph_fit.png", dpi=110)


if __name__ == "__main__":
    main()
