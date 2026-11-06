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
