"""Sky geometry of dwarf galaxies."""
import numpy as np


def tangent_plane(ra, dec, ra0, dec0):
    """Gnomonic projection (deg) about (ra0, dec0); xi points East, eta North."""
    a, d, a0, d0 = map(np.radians, (ra, dec, ra0, dec0))
    cosc = np.sin(d0) * np.sin(d) + np.cos(d0) * np.cos(d) * np.cos(a - a0)
    xi = np.cos(d) * np.sin(a - a0) / cosc
    eta = (np.cos(d0) * np.sin(d) - np.sin(d0) * np.cos(d) * np.cos(a - a0)) / cosc
    return np.degrees(xi), np.degrees(eta)


def elliptical_radius(ra, dec, ra0, dec0, pa, ell):
    """Elliptical radius (deg, along the major axis); pa East of North (deg), ell = 1 - b/a."""
    xi, eta = tangent_plane(ra, dec, ra0, dec0)
    p = np.radians(pa)
    x_major = xi * np.sin(p) + eta * np.cos(p)
    x_minor = xi * np.cos(p) - eta * np.sin(p)
    return np.hypot(x_major, x_minor / (1 - ell))


def perspective_vlos(ra, dec, ra0, dec0, vsys, pmra, pmdec, D):
    """Line-of-sight velocity [km/s] at each (ra, dec) of a body moving as a whole with the 3D heliocentric velocity of its
    centre (ra0, dec0): line-of-sight vsys, proper motion (pmra = mu_alpha*, pmdec) [mas/yr] at distance D [kpc]."""
    def basis(a, d):
        a, d = np.radians(a), np.radians(d)
        n = np.stack([np.cos(d) * np.cos(a), np.cos(d) * np.sin(a), np.sin(d)], -1)
        ea = np.stack([-np.sin(a), np.cos(a), np.zeros_like(a)], -1)
        ed = np.stack([-np.sin(d) * np.cos(a), -np.sin(d) * np.sin(a), np.cos(d)], -1)
        return n, ea, ed
    n0, ea0, ed0 = basis(np.asarray(ra0, float), np.asarray(dec0, float))
    V = vsys * n0 + 4.740470 * D * (pmra * ea0 + pmdec * ed0)
    return basis(np.asarray(ra, float), np.asarray(dec, float))[0] @ V
