import numpy as np

from olddisks.geometry import elliptical_radius


def test_major_and_minor_axis():
    # Galaxy at (150, 30) with PA = 90 deg (major axis along East), ellipticity 0.5.
    r = 0.1  # deg
    ra_major = 150 + r / np.cos(np.radians(30))
    R = elliptical_radius(np.array([ra_major, 150.0]), np.array([30.0, 30 + r]), 150.0, 30.0, 90.0, 0.5)
    np.testing.assert_allclose(R, [r, r / 0.5], rtol=1e-3)


def test_perspective_vlos():
    from olddisks.geometry import perspective_vlos
    # at the centre: v_sys; with zero proper motion: v_sys cos(rho); proper motion of 1 mas/yr at 100 kpc gives
    # 474 km/s transverse, seen as 474 sin(rho) along the motion direction (east here)
    v = perspective_vlos(np.array([150.0, 150.0, 151.0]), np.array([0.0, 1.0, 0.0]), 150.0, 0.0, 100.0, 0.0, 0.0, 100.0)
    np.testing.assert_allclose(v, [100, 100 * np.cos(np.radians(1)), 100 * np.cos(np.radians(1))], rtol=1e-9)
    v = perspective_vlos(np.array([151.0]), np.array([0.0]), 150.0, 0.0, 0.0, 1.0, 0.0, 100.0)
    np.testing.assert_allclose(v, [4.740470 * 100 * np.sin(np.radians(1))], rtol=1e-4)
