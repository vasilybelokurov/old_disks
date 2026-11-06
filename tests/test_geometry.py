import numpy as np

from olddisks.geometry import elliptical_radius


def test_major_and_minor_axis():
    # Galaxy at (150, 30) with PA = 90 deg (major axis along East), ellipticity 0.5.
    r = 0.1  # deg
    ra_major = 150 + r / np.cos(np.radians(30))
    R = elliptical_radius(np.array([ra_major, 150.0]), np.array([30.0, 30 + r]), 150.0, 30.0, 90.0, 0.5)
    np.testing.assert_allclose(R, [r, r / 0.5], rtol=1e-3)
