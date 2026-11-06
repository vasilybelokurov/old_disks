"""Pre-registered real-data tests for NGC 185 (primary: T_rot, alpha 0.01) and NGC 147 (alpha 0.0025), plan in
docs/plans/de_model.md. Baseline: R0 from the [Fe/H]-blind fit (grid for NGC 185, linear for NGC 147), LVDB PA, all members.
Variants (robustness, not extra tests): TRGB cut (drop I < I_TRGB + 0.3, I_TRGB = (m-M) - 4.05, Rizzi+2007 2007ApJ...661..815R;
assumes Geha I is extinction-corrected, UNVERIFIED), PA +-10 deg, each R0 of the grid (NGC 185).
Amplitude: slope = Delta(V/sigma) / C(k), sd from the null mocks (results/de_mocks_summary.csv), k = 1 (headline) and 0.5.
Outputs: results/de_fits.csv, figures/de_fit.png.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from olddisks.defit import blind_fit, h, permutation_test, sample, statistics, tercile_fit
from olddisks.paths import LVDB

SETUP = {"NGC 185": dict(R0s=[0.5, 1.0, 2.0, 4.0], Rstar=2.5, alpha=0.01),
         "NGC 147": dict(R0s=[None], Rstar=1.5, alpha=0.0025)}
NPERM = 10000


def run(gal, s, R0, label, rng):
    st = statistics(s, R0, SETUP[gal]["Rstar"], NPERM, rng)
    b = st.pop("blind")
    return dict(galaxy=gal, variant=label, n=len(s), n_feh=int(np.isfinite(s.feh).sum()), R0=R0, vsys=b["vsys"], V=b["V"],
                sigma0=b["sigma0"], f=b["f"], **{k: v for k, v in st.items() if k not in ("ratio", "zmean")},
                ratio_mp=st["ratio"][0], ratio_mid=st["ratio"][1], ratio_mr=st["ratio"][2])


def main():
    rng = np.random.default_rng(2026)
    lvdb = pd.read_csv(LVDB).drop_duplicates("name").set_index("name")
    mocks = pd.read_csv("results/de_mocks_summary.csv")
    rows = []
    for gal, cfg in SETUP.items():
        s = sample(gal)
        best = min(cfg["R0s"], key=lambda R0: blind_fit(s, R0, cfg["Rstar"])["nll"])
        rows.append(run(gal, s, best, "baseline", rng))
        itrgb = lvdb.loc[gal, "distance_modulus"] - 4.05
        rows.append(run(gal, s[s.imag >= itrgb + 0.3], best, f"TRGB cut (I >= {itrgb + 0.3:.2f})", rng))
        for dpa in (-10, 10):
            rows.append(run(gal, sample(gal, dpa), best, f"PA {dpa:+d}", rng))
        for R0 in cfg["R0s"]:
            if R0 != best:
                rows.append(run(gal, s, R0, f"R0 = {R0}", rng))
    r = pd.DataFrame(rows)
    for k in (1.0, 0.5):
        cal = mocks[(mocks.slope == 0) & (mocks.k == k)].set_index("galaxy")
        r[f"slope_k{k}"] = r.delta / r.galaxy.map(cal.C)
        r[f"slope_sd_k{k}"] = r.galaxy.map(cal.slope_sd)
    r["alpha"] = r.galaxy.map({g: c["alpha"] for g, c in SETUP.items()})
    r["detected_rot"] = r.p_rot < r.alpha
    r.round(4).to_csv("results/de_fits.csv", index=False)
    pd.set_option("display.width", 250)
    print(r[["galaxy", "variant", "n_feh", "R0", "V", "sigma0", "t_rot", "p_rot", "t_sig", "p_sig", "ratio_mp", "ratio_mid",
             "ratio_mr", "delta", "slope_k1.0", "slope_sd_k1.0", "slope_k0.5", "slope_sd_k0.5", "detected_rot"]].round(3).to_string(index=False))
    figure(r, rng)


def figure(r, rng):
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    for ax, gal in zip(axes, SETUP):
        base = r[(r.galaxy == gal) & (r.variant == "baseline")].iloc[0]
        s = sample(gal)
        R0, Rstar = base.R0 if np.isfinite(base.R0) else None, SETUP[gal]["Rstar"]
        b = blind_fit(s, R0, Rstar)
        t = s[np.isfinite(s.feh)]
        q = pd.qcut(t.feh, 3, labels=False)
        y = (t.v - b["vsys"]) * np.sign(t.x) * np.sign(b["V"])
        Rg = np.linspace(0.05, t.R.max(), 100)
        for i, col, lab in ((0, "C0", "metal-poor third"), (2, "C3", "metal-rich third")):
            sel = q == i
            ax[0].scatter(t.R[sel], y[sel], s=8, c=col, alpha=0.5, label=f"{lab} ([Fe/H] med {t.feh[sel].median():.2f})")
            vs = [base.ratio_mp, base.ratio_mr][i // 2] * base.sigma0
            ax[0].plot(Rg, vs * h(Rg, R0, Rstar) * 0.93, color=col, lw=2)
        ax[0].plot(Rg, abs(b["V"]) * h(Rg, R0, Rstar) * 0.93, "k--", lw=1, label="all stars ([Fe/H]-blind)")
        ax[0].axhline(0, color="k", lw=0.5), ax[0].axvline(Rstar, color="0.6", ls=":", lw=1)
        ax[0].set_ylim(-80, 80), ax[0].set_xlabel(r"$R/R_\mathrm{half}$"), ax[0].set_ylabel(r"$(v - v_\mathrm{sys})\,\mathrm{sign}(x)$ [km/s]")
        ax[0].set_title(f"{gal}: rotation by [Fe/H] tercile (lines: tercile V/σ × σ0, ⟨|x|/R⟩ = 0.93)", fontsize=9)
        ax[0].legend(fontsize=7, loc="lower right")
        yrot = (s.v - b["vsys"]) * np.sign(s.x) * np.sign(b["V"])
        pt = permutation_test(yrot, s.R, s.feh, NPERM, rng)
        ax[1].hist(pt["null"], bins=60, color="0.7", density=True, label="within-R permutations")
        ax[1].axvline(base.t_rot, color="C3", lw=2, label=f"observed T_rot = {base.t_rot:.3f}, p = {base.p_rot:.3f}")
        ax[1].set_xlabel("T_rot (mean within-bin Spearman ρ)"), ax[1].legend(fontsize=7)
        ax[1].set_title(f"pre-registered α = {SETUP[gal]['alpha']}", fontsize=9)
        v = r[r.galaxy == gal]
        yy = np.arange(len(v))
        ax[2].errorbar(v["slope_k1.0"], yy, xerr=1.96 * v["slope_sd_k1.0"], fmt="o", color="k", label="k = 1 (headline), 95 %")
        ax[2].errorbar(v["slope_k0.5"], yy + 0.2, xerr=1.96 * v["slope_sd_k0.5"], fmt="s", color="0.5", mfc="white", label="k = 0.5")
        ax[2].axvline(0, color="k", lw=0.5)
        ax[2].set_yticks(yy, v.variant, fontsize=7), ax[2].invert_yaxis()
        ax[2].set_xlabel(r"d(V/σ)/d[Fe/H] [dex$^{-1}$] (tercile contrast / mock dilution)"), ax[2].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig("figures/de_fit.png", dpi=110)


if __name__ == "__main__":
    main()
