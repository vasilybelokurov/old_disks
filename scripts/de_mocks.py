"""Mocks for the pre-registered dE tests (docs/plans/de_model.md, agreed plan).

For NGC 147 (linear h, Rstar 1.5 R_half) and NGC 185 (R0 chosen from the grid by the [Fe/H]-blind likelihood, Rstar 2.5 R_half):
truth = the [Fe/H]-blind fit to the real velocities. Mocks at the real positions, velocity and [Fe/H] errors, with the real
missing-[Fe/H] pattern, for slopes d(V/sigma)/d[Fe/H] in SLOPES and error scales k in KS, NMOCK each, 1000 permutations.
Outputs: results/de_mocks.csv (one row per mock), results/de_mocks_summary.csv, figures/de_mocks.png.
"""
from concurrent.futures import ProcessPoolExecutor

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from olddisks.defit import blind_fit, make_mock, sample, statistics

SLOPES = [0.0, 0.3, 0.5, 1.0]
KS = [1.0, 0.5]
NMOCK = 200
R0_GRID = [0.5, 1.0, 2.0, 4.0]
SETUP = {"NGC 147": dict(R0s=[None], Rstar=1.5), "NGC 185": dict(R0s=R0_GRID, Rstar=2.5)}
ALPHA = {("NGC 185", "rot"): 0.01}  # pre-registered primary; every other test at 0.0025
ALPHA_OTHER = 0.0025


def setup(gal):
    s = sample(gal)
    fits = [blind_fit(s, R0, SETUP[gal]["Rstar"]) for R0 in SETUP[gal]["R0s"]]
    return s, min(fits, key=lambda b: b["nll"])


def one(args):
    gal, slope, k, seed = args
    s, truth = setup(gal)
    truth = dict(truth, f=max(truth["f"], 0.05))  # pre-registered 5 % contaminants
    rng = np.random.default_rng(seed)
    st = statistics(make_mock(s, truth, slope, k, rng), truth["R0"], truth["Rstar"], nperm=1000, rng=rng)
    return dict(galaxy=gal, slope=slope, k=k, seed=seed, p_rot=st["p_rot"], p_sig=st["p_sig"], t_rot=st["t_rot"],
                delta=st["delta"])


def main():
    jobs = [(g, sl, k, (SLOPES.index(sl) * 10 + KS.index(k)) * 100000 + n + (0 if g == "NGC 147" else 50000)) for g in SETUP for sl in SLOPES for k in KS for n in range(NMOCK)]
    with ProcessPoolExecutor() as ex:
        rows = list(ex.map(one, jobs, chunksize=20))
    m = pd.DataFrame(rows)
    m.to_csv("results/de_mocks.csv", index=False)

    out = []
    for (g, k, sl), d in m.groupby(["galaxy", "k", "slope"]):
        a_rot = ALPHA.get((g, "rot"), ALPHA_OTHER)
        out.append(dict(galaxy=g, k=k, slope=sl, rate_rot=(d.p_rot < a_rot).mean(), alpha_rot=a_rot,
                        rate_sig=(d.p_sig < ALPHA_OTHER).mean(), delta_mean=d.delta.mean(), delta_sd=d.delta.std()))
    summ = pd.DataFrame(out)
    cal = summ[summ.slope == 1.0].set_index(["galaxy", "k"]).delta_mean - summ[summ.slope == 0].set_index(["galaxy", "k"]).delta_mean
    summ["C"] = [cal.loc[(g, k)] for g, k in zip(summ.galaxy, summ.k)]
    summ["slope_sd"] = summ.delta_sd / summ.C
    summ.round(4).to_csv("results/de_mocks_summary.csv", index=False)
    print(summ.round(3).to_string(index=False))

    fig, ax = plt.subplots(1, 3, figsize=(14, 4))
    for (g, k), d in summ.groupby(["galaxy", "k"]):
        ls = "-" if k == 1 else "--"
        col = "C0" if g == "NGC 185" else "C1"
        ax[0].plot(d.slope, d.rate_rot, ls, marker="o", color=col, label=f"{g}, k={k}")
        ax[1].plot(d.slope, d.delta_mean, ls, marker="o", color=col)
        ax[1].fill_between(d.slope, d.delta_mean - d.delta_sd, d.delta_mean + d.delta_sd, color=col, alpha=0.1)
    ax[0].axhline(0.8, color="k", ls=":", lw=0.8)
    ax[0].set_xlabel("true d(V/σ)/d[Fe/H] [dex⁻¹]"), ax[0].set_ylabel("rate p_rot < α (α = 0.01 NGC 185, 0.0025 NGC 147)")
    ax[0].legend(fontsize=7), ax[0].set_title("detection rate of T_rot (slope 0: false-positive rate)", fontsize=9)
    ax[1].plot([0, 1], [0, 1], "k:", lw=0.8)
    ax[1].set_xlabel("true slope"), ax[1].set_ylabel("tercile contrast Δ(V/σ) (mean ± sd)")
    ax[1].set_title("dilution: Δ against true slope", fontsize=9)
    for g, col in (("NGC 185", "C0"), ("NGC 147", "C1")):
        d = m[(m.galaxy == g) & (m.slope == 0) & (m.k == 1)]
        ax[2].hist(d.p_rot, bins=20, range=(0, 1), histtype="step", color=col, label=f"{g} null (k=1)")
    ax[2].set_xlabel("p_rot under the null"), ax[2].legend(fontsize=7), ax[2].set_title("null p-values (should be flat)", fontsize=9)
    fig.tight_layout()
    fig.savefig("figures/de_mocks.png", dpi=120)


if __name__ == "__main__":
    main()
