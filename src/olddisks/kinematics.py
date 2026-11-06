"""Line-of-sight rotation + dispersion model for a dwarf galaxy.

Mean velocity v_sys + V f(r) cos(phi - theta), f(r) = min(r, 1) (solid body inside R_half, flat outside), where r is
the elliptical radius in R_half and phi the position angle on the sky (rad, East of North). Gaussian dispersion sigma,
convolved with each star's error. V is the projected amplitude (V sin i for a disc).
"""
import numpy as np
from scipy.optimize import minimize


def mean_velocity(r, phi, vsys, vrot, theta):
    return vsys + vrot * np.minimum(r, 1.0) * np.cos(phi - theta)


def mock_velocities(r, phi, ev, vsys, vrot, theta, sigma, rng):
    return mean_velocity(r, phi, vsys, vrot, theta) + rng.normal(0, np.sqrt(sigma**2 + ev**2))


def fit_rotation(r, phi, v, ev):
    """Maximum-likelihood (vsys, vrot >= 0, theta, sigma)."""
    def nll(p):
        vsys, a, b, lns = p  # rotation as (a, b) = V (cos theta, sin theta): no degeneracy at V = 0
        mu = vsys + np.minimum(r, 1.0) * (a * np.cos(phi) + b * np.sin(phi))
        var = np.exp(2 * lns) + ev**2
        return 0.5 * np.sum((v - mu) ** 2 / var + np.log(var))
    p0 = [np.median(v), 0.0, 0.0, np.log(np.std(v))]
    vsys, a, b, lns = minimize(nll, p0, method="Nelder-Mead", options=dict(maxiter=4000, xatol=1e-4, fatol=1e-6)).x
    return dict(vsys=vsys, vrot=np.hypot(a, b), theta=np.arctan2(b, a), sigma=np.exp(lns))
