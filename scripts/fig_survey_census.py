"""Diagnostic figure for the WSDB survey census (scripts/build_survey_census.py).

One row per (galaxy, survey) with >= 25 candidates:
  left   - line-of-sight velocity minus v_sys of all stars within 3 R_half; candidate window shaded;
  middle - sky positions (tangent plane, units of R_half along the major axis), candidates coloured by [Fe/H];
  right  - [Fe/H] distribution of candidates.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from olddisks.paths import LVDB, SURVEY_MEMBERS
from olddisks.geometry import tangent_plane

stars = pd.read_parquet(SURVEY_MEMBERS)
census = pd.read_csv("data/survey_census.csv")
lvdb = pd.read_csv(LVDB).set_index("name")
rows = census[census.n_v >= 25].sort_values("n_v", ascending=False)

fig, axes = plt.subplots(len(rows), 3, figsize=(10, 2.6 * len(rows)), squeeze=False)
for ax, (_, r) in zip(axes, rows.iterrows()):
    g = lvdb.loc[r.galaxy]
    s = stars[(stars.galaxy == r.galaxy) & (stars.survey == r.survey)]
    c = s[s.candidate]
    dv = s.v - g.vlos_systemic
    ax[0].hist(dv, bins=np.arange(-150, 151, 4), color="0.5")
    ax[0].hist(c.v - g.vlos_systemic, bins=np.arange(-150, 151, 4), color="C0")
    ax[0].set_xlim(-150, 150)
    ax[0].set_xlabel(r"$v_\mathrm{los} - v_\mathrm{sys}$ [km/s]")
    ax[0].set_title(f"{r.galaxy}, {r.survey}: {len(c)} cand. / {len(s)} in cone ({(np.abs(dv) > 150).sum()} off-axis)",
                    fontsize=8, loc="left")
    xi, eta = tangent_plane(c.ra.values, c.dec.values, g.ra, g.dec)
    sc = ax[1].scatter(xi * 60 / g.rhalf, eta * 60 / g.rhalf, c=c.feh, s=6, cmap="viridis", vmin=-3, vmax=-0.5)
    ax[1].set_aspect("equal")
    ax[1].set_xlim(3, -3), ax[1].set_ylim(-3, 3)
    ax[1].set_xlabel(r"$\xi / R_\mathrm{half}$ (East left)"), ax[1].set_ylabel(r"$\eta / R_\mathrm{half}$")
    fig.colorbar(sc, ax=ax[1], label="[Fe/H]", fraction=0.046)
    ax[2].hist(c.feh.dropna(), bins=np.arange(-4, 0.01, 0.1), color="C0")
    ax[2].set_xlabel("[Fe/H] (survey pipeline)")
fig.tight_layout()
fig.savefig("figures/survey_census.png", dpi=110)
print("figures/survey_census.png")
