import numpy as np
import pytest

from olddisks.defit import permutation_test, sample


@pytest.mark.parametrize("galaxy, n_members, n_feh", [("NGC 147", 520, 227), ("NGC 185", 442, 312)])
def test_sample_counts(galaxy, n_members, n_feh):
    """All Geha+2010 members (520, 442), Ho+2015 [Fe/H] attached to the matched ones (227, 312)."""
    s = sample(galaxy)
    assert len(s) == n_members and np.isfinite(s.feh).sum() == n_feh
    assert (s.ev >= 1.5).all() and np.isfinite(s[["x", "R", "v", "imag"]]).all().all()


def _mock(rng, n, slope):
    """Stars on a major-axis strip; rotation amplitude (km/s) linear in [Fe/H] with `slope`; [Fe/H] with a radial gradient."""
    R = rng.uniform(0.1, 3, n)
    x = R * np.sign(rng.uniform(-1, 1, n))
    feh = -1 - 0.2 * R + rng.normal(0, 0.5, n)
    v = -200 + (15 + slope * (feh + 1)) * np.tanh(R) * np.sign(x) + rng.normal(0, 20, n)
    return x, R, v, feh


def test_permutation_null_rate():
    """No dependence on [Fe/H] (but a radial gradient in [Fe/H] and a rising curve): p < 0.05 in about 5 % of mocks."""
    rng = np.random.default_rng(3)
    p = [permutation_test(*_mock(rng, 300, 0.0), vsys=-200, nperm=500, rng=rng)["p_rot"] for _ in range(200)]
    assert 0.02 < np.mean(np.array(p) < 0.05) < 0.09


def test_permutation_detects_strong_trend():
    rng = np.random.default_rng(4)
    p = [permutation_test(*_mock(rng, 300, 30.0), vsys=-200, nperm=500, rng=rng)["p_rot"] for _ in range(20)]
    assert np.mean(np.array(p) < 0.01) > 0.9
