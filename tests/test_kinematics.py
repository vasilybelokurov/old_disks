import numpy as np

from olddisks.kinematics import fit_rotation, mock_velocities


def test_recovers_rotation_and_dispersion():
    rng = np.random.default_rng(1)
    n = 4000
    r, phi = rng.uniform(0, 3, n), rng.uniform(0, 2 * np.pi, n)
    ev = rng.uniform(1, 5, n)
    v = mock_velocities(r, phi, ev, vsys=-130.0, vrot=15.0, theta=1.0, sigma=12.0, rng=rng)
    f = fit_rotation(r, phi, v, ev)
    assert abs(f["vrot"] - 15.0) < 1.0 and abs(f["sigma"] - 12.0) < 0.6
    assert abs(f["vsys"] + 130.0) < 0.8 and abs(np.angle(np.exp(1j * (f["theta"] - 1.0)))) < 0.1
