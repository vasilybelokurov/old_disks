"""First look at the M31 dwarf ellipticals (Ho+2015 [Fe/H] joined to Geha velocities). No model.

One row per galaxy (NGC 147, NGC 185, NGC 205):
  (a) sky, R_half units (LVDB ellipse at 1 and 3 R_half): all velocity members (grey) and the Ho+2015 subset coloured
      by v - v_sys (shows rotation);
  (b) v - v_sys against the coordinate along the photometric major axis, Ho+2015 stars split at their median [Fe/H]:
      metal-poor (blue) and metal-rich (red), with running medians (bins of 0.5 R_half, >= 8 stars);
  (c) [Fe/H] histogram of the Ho+2015 subset, and the median [Fe/H] error as a bar.
Output: figures/de_data.png
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from olddisks.geometry import tangent_plane
from olddisks.paths import DATA, LVDB

GALAXIES = ["NGC 147", "NGC 185", "NGC 205"]
lvdb = pd.read_csv(LVDB).drop_duplicates("name").set_index("name")
m = pd.read_parquet(DATA / "inventory_stars.parquet")

fig, axes = plt.subplots(len(GALAXIES), 3, figsize=(13, 4 * len(GALAXIES)))
t = np.linspace(0, 2 * np.pi, 200)
for ax, gal in zip(axes, GALAXIES):
    g = lvdb.loc[gal]
    allv = m[(m.galaxy == gal) & m.in_union]
    ho = m[(m.galaxy == gal) & (m.source == "Ho+2015")]
    p, e = np.radians(g.position_angle), g.ellipticity

    def xy(s):
        xi, eta = tangent_plane(s.ra.to_numpy(), s.dec.to_numpy(), g.ra, g.dec)
        return xi * 60 / g.rhalf, eta * 60 / g.rhalf

    x, y = xy(allv)
    ax[0].scatter(x, y, s=3, c="0.75", label=f"all v members ({len(allv)})")
    x, y = xy(ho)
    dv = ho.v - g.vlos_systemic
    sc = ax[0].scatter(x, y, s=10, c=dv, cmap="RdBu_r", vmin=-50, vmax=50, label=f"Ho+2015 ({len(ho)})")
    for r in (1, 3):
        a, b = r * np.cos(t), r * (1 - e) * np.sin(t)
        ax[0].plot(a * np.sin(p) + b * np.cos(p), a * np.cos(p) - b * np.sin(p), "k-", lw=0.6)
    ax[0].set_aspect("equal"), ax[0].invert_xaxis()
    ax[0].set_xlabel(r"$\xi/R_\mathrm{half}$ (East left)"), ax[0].set_ylabel(r"$\eta/R_\mathrm{half}$")
    ax[0].set_title(f"{gal}  PA={g.position_angle:.0f}°  e={e:.2f}  σ_LVDB={g.vlos_sigma:.0f} km/s", fontsize=9)
    ax[0].legend(fontsize=7, loc="lower left")
    fig.colorbar(sc, ax=ax[0], label=r"$v - v_\mathrm{sys}$ [km/s]", fraction=0.046)

    major = x * np.sin(p) + y * np.cos(p)
    med = ho.feh.median()
    for sel, col, lab in ((ho.feh < med, "C0", f"[Fe/H] < {med:.2f}"), (ho.feh >= med, "C3", f"[Fe/H] ≥ {med:.2f}")):
        sel = sel.to_numpy()
        ax[1].scatter(major[sel], dv[sel], s=6, c=col, alpha=0.5, label=f"{lab} ({sel.sum()})")
        bins = np.arange(-6, 6.01, 0.5)
        k = np.digitize(major[sel], bins)
        cen = [(bins[i - 1] + bins[i]) / 2 for i in np.unique(k) if (k == i).sum() >= 8 and 0 < i < len(bins)]
        val = [np.median(dv[sel][k == i]) for i in np.unique(k) if (k == i).sum() >= 8 and 0 < i < len(bins)]
        ax[1].plot(cen, val, "-o", color=col, ms=4)
    ax[1].axhline(0, color="k", lw=0.5)
    ax[1].set_ylim(-120, 120)
    ax[1].set_xlabel(r"major-axis coordinate [$R_\mathrm{half}$]"), ax[1].set_ylabel(r"$v - v_\mathrm{sys}$ [km/s]")
    ax[1].legend(fontsize=7)

    ax[2].hist(ho.feh, bins=np.arange(-3, 0.81, 0.1), color="0.5")
    ax[2].errorbar(med, ax[2].get_ylim()[1] * 0.9, xerr=ho.efeh.median(), fmt="o", color="k",
                   label=f"median error {ho.efeh.median():.2f} dex")
    ax[2].set_xlabel("[Fe/H] (Ho+2015, CaT)"), ax[2].legend(fontsize=7)
fig.tight_layout()
fig.savefig("figures/de_data.png", dpi=110)
print("figures/de_data.png")
