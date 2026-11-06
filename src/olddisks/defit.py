"""Rotation versus [Fe/H] in the M31 dwarf ellipticals NGC 147 and NGC 185 (plan: docs/plans/de_model.md).

sample(): all Geha+2010 velocity members with Ho+2015 [Fe/H] attached where matched (catalogues.ho2015), in R_half units
along the LVDB photometric axes: x (major-axis coordinate), R (elliptical radius).
permutation_test(): the model-free primary statistic. T_rot = Spearman rho of (v - v_sys) sign(x) with [Fe/H], computed within
radial bins (equal numbers) and averaged; its null distribution from permutations of [Fe/H] within the same bins, which keeps
the radial [Fe/H] gradient and the rising rotation curve from producing a signal.
"""
import numpy as np
import pandas as pd
from scipy.stats import rankdata

from olddisks.catalogues import ho2015, table
from olddisks.geometry import tangent_plane
from olddisks.paths import LVDB

EV_FLOOR = 1.5  # km/s; 41 Geha+2010 rows have e_HRV = 0


def sample(galaxy):
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    g = pd.read_csv(LVDB).drop_duplicates("name").set_index("name").loc[galaxy]
    t = table("J/ApJ/711/361/veloc")
    t = t[t.NGC == int(galaxy.split()[1])].reset_index(drop=True)
    ho = ho2015()
    ho = ho[ho.galaxy == galaxy].set_index("star_id")
    c = SkyCoord(t.RAJ2000.astype(str), t.DEJ2000.astype(str), unit=(u.hourangle, u.deg))
    xi, eta = tangent_plane(c.ra.deg, c.dec.deg, g.ra, g.dec)
    p = np.radians(g.position_angle)
    xmaj = (xi * np.sin(p) + eta * np.cos(p)) * 60 / g.rhalf
    xmin = (xi * np.cos(p) - eta * np.sin(p)) * 60 / g.rhalf
    name = t.Name.astype(str)
    return pd.DataFrame(dict(star_id=name, x=xmaj, R=np.hypot(xmaj, xmin / (1 - g.ellipticity)), v=t.HRV.to_numpy(float),
                             ev=np.maximum(t.e_HRV.to_numpy(float), EV_FLOOR), imag=t.Imag.to_numpy(float),
                             feh=name.map(ho.feh).to_numpy(float), efeh=name.map(ho.efeh).to_numpy(float)))


def _binned_rho(y, z, bins):
    """Mean over bins of the Spearman rho of y and z (z shape (..., n) allowed for permutations)."""
    out = []
    for b in bins:
        ry = rankdata(y[b])
        rz = rankdata(z[..., b], axis=-1)
        ry, rz = ry - ry.mean(), rz - rz.mean(axis=-1, keepdims=True)
        out.append((rz * ry).sum(-1) / np.sqrt((rz**2).sum(-1) * (ry**2).sum()))
    return np.mean(out, axis=0)


def permutation_test(x, R, v, feh, vsys, nbins=5, nperm=10000, rng=None):
    """T_rot and its two-sided permutation p-value (stars with finite [Fe/H] only)."""
    rng = np.random.default_rng() if rng is None else rng
    ok = np.isfinite(feh)
    x, R, v, feh = x[ok], R[ok], v[ok], feh[ok]
    edges = np.quantile(R, np.linspace(0, 1, nbins + 1))
    k = np.clip(np.searchsorted(edges, R, side="right") - 1, 0, nbins - 1)
    bins = [np.flatnonzero(k == i) for i in range(nbins)]
    y = (v - vsys) * np.sign(x)
    t = _binned_rho(y, feh, bins)
    perm = np.tile(feh, (nperm, 1))
    for b in bins:
        perm[:, b] = rng.permuted(perm[:, b], axis=1)
    null = _binned_rho(y, perm, bins)
    return dict(t_rot=t, p_rot=(np.sum(np.abs(null) >= abs(t)) + 1) / (nperm + 1), null=null)
