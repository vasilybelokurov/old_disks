"""Inventory figures (scripts/build_inventory.py).

figures/inventory_counts.png: members with v and [Fe/H] per galaxy (union over sources, deduplicated at 1 arcsec), stacked
    by source, galaxies with >= 40 such members ordered by M_V.
figures/inventory_maps.png: members with [Fe/H] on the sky (tangent plane, R_half units; LVDB ellipse at 1 and 3 R_half),
    coloured by [Fe/H]; marker shape = whether membership is the paper's flag (o) or our velocity window (x).
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from olddisks.geometry import tangent_plane
from olddisks.paths import DATA, LVDB

m = pd.read_parquet(DATA / "inventory_stars.parquet")
lvdb = pd.read_csv(LVDB).drop_duplicates("name").set_index("name")
u = m[m.in_union & np.isfinite(m.feh)]
n = u.groupby("galaxy").size()
gals = lvdb.loc[n[n >= 40].index].M_V.sort_values().index

# counts
tab = u[u.galaxy.isin(gals)].groupby(["galaxy", "source"]).size().unstack(fill_value=0).reindex(gals)
tab = tab[tab.sum().sort_values(ascending=False).index]
fig, ax = plt.subplots(figsize=(9, 6))
left = np.zeros(len(tab))
colors = plt.cm.tab20(np.linspace(0, 1, 20))
for i, src in enumerate(tab.columns):
    ax.barh(range(len(tab)), tab[src], left=left, color=colors[i % 20] if i < 19 else "0.6", label=src, height=0.8)
    left += tab[src].to_numpy()
for i, tot in enumerate(left):
    ax.text(tot + 20, i, f"{tot:.0f}", va="center", fontsize=7)
ax.set_yticks(range(len(tab)), [f"{g} ({lvdb.M_V[g]:.1f})" for g in tab.index], fontsize=8)
ax.invert_yaxis()
ax.set_xlabel("members with $v_\\mathrm{los}$ and [Fe/H] (union, 1\" dedup)")
ax.legend(fontsize=6, ncol=2, loc="lower right")
fig.tight_layout()
fig.savefig("figures/inventory_counts.png", dpi=120)

# maps
ncol = 5
nrow = int(np.ceil(len(gals) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(3.0 * ncol, 3.0 * nrow), squeeze=False)
t = np.linspace(0, 2 * np.pi, 200)
for ax, gal in zip(axes.flat, gals):
    g, s = lvdb.loc[gal], u[u.galaxy == gal]
    xi, eta = tangent_plane(s.ra.to_numpy(), s.dec.to_numpy(), g.ra, g.dec)
    xi, eta = xi * 60 / g.rhalf, eta * 60 / g.rhalf
    flag = s.member_paper.notna().to_numpy()
    for mk, sel in (("o", flag), ("x", ~flag)):
        ax.scatter(xi[sel], eta[sel], c=s.feh[sel], s=3, marker=mk, cmap="viridis", vmin=-3, vmax=-0.5, lw=0.5)
    p, e = np.radians(np.nan_to_num(g.position_angle)), np.nan_to_num(g.ellipticity)
    for r in (1, 3):
        a, b = r * np.cos(t), r * (1 - e) * np.sin(t)
        ax.plot(a * np.sin(p) + b * np.cos(p), a * np.cos(p) - b * np.sin(p), "k-", lw=0.6)
    lim = 6
    ax.set_xlim(lim, -lim), ax.set_ylim(-lim, lim), ax.set_aspect("equal")
    out = (np.hypot(xi, eta) > lim).sum()
    ax.set_title(f"{gal}: {len(s)} ({out} outside)", fontsize=8)
    ax.tick_params(labelsize=6)
for ax in axes.flat[len(gals):]:
    ax.axis("off")
fig.supxlabel(r"$\xi / R_\mathrm{half}$ (East left)"), fig.supylabel(r"$\eta / R_\mathrm{half}$")
fig.tight_layout()
fig.savefig("figures/inventory_maps.png", dpi=110)
print("figures/inventory_counts.png figures/inventory_maps.png")
