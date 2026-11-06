"""Heatmap of archival spectroscopic telescope time per galaxy and instrument (data/archive_census.csv).

Uses ESO raw and KOA rows only (ESO phase3 hours are summed over stars). Cell text: hours / programmes.
Galaxies ordered by M_V (brightest at top); only galaxies with any spectroscopy are shown.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

c = pd.read_csv("data/archive_census.csv")
c = c[c.archive != "ESO phase3"]
hours = c.pivot_table(index="galaxy", columns="instrument", values="hours", aggfunc="sum")
progs = c.pivot_table(index="galaxy", columns="instrument", values="n_prog", aggfunc="sum")
order = c.drop_duplicates("galaxy").set_index("galaxy").M_V.sort_values().index
cols = ["DEIMOS", "GIRAFFE", "UVES", "FORS2", "FORS1", "VIMOS", "MUSE", "XSHOOTER"]
cols = [k for k in cols if k in hours.columns]
hours, progs = hours.reindex(index=order, columns=cols), progs.reindex(index=order, columns=cols)

fig, ax = plt.subplots(figsize=(1.0 * len(cols) + 3, 0.28 * len(order) + 1.2))
im = ax.imshow(np.log10(hours.values), cmap="Blues", vmin=-0.5, vmax=2.7, aspect="auto")
for i, j in zip(*np.where(np.isfinite(hours.values))):
    h = hours.values[i, j]
    ax.text(j, i, f"{h:.0f}/{progs.values[i, j]:.0f}", ha="center", va="center", fontsize=6,
            color="white" if np.log10(h) > 1.6 else "black")
ax.set_xticks(range(len(cols)), cols, rotation=45, ha="right")
mv = c.groupby("galaxy").M_V.first()
ax.set_yticks(range(len(order)), [f"{g} ({mv[g]:.1f})" for g in order], fontsize=7)
ax.set_title("Archival spectroscopy within 3 R_half: hours / programmes (ESO raw + KOA)", fontsize=9, loc="left")
fig.colorbar(im, ax=ax, label="log10 hours", fraction=0.04)
fig.tight_layout()
fig.savefig("figures/archive_census.png", dpi=120)
print("figures/archive_census.png")
