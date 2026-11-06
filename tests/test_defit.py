import numpy as np
import pandas as pd
import pytest

from olddisks.defit import blind_fit, make_mock, permutation_test, sample, statistics, tercile_fit


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
    p = []
    for _ in range(200):
        x, R, v, feh = _mock(rng, 300, 0.0)
        p.append(permutation_test((v + 200) * np.sign(x), R, feh, nperm=500, rng=rng)["p"])
    assert 0.02 < np.mean(np.array(p) < 0.05) < 0.09


def test_permutation_detects_strong_trend():
    rng = np.random.default_rng(4)
    p = []
    for _ in range(20):
        x, R, v, feh = _mock(rng, 300, 30.0)
        p.append(permutation_test((v + 200) * np.sign(x), R, feh, nperm=500, rng=rng)["p"])
    assert np.mean(np.array(p) < 0.01) > 0.9


def _truth_sample(rng, n=4000):
    R = rng.uniform(0.1, 4, n)
    x = R * np.sign(rng.uniform(-1, 1, n)) * 0.95
    return pd.DataFrame(dict(x=x, y=np.zeros(n), R=R, ev=np.full(n, 3.0), feh=-1 + rng.normal(0, 0.6, n),
                             efeh=np.full(n, 0.4)))


def test_blind_fit_recovers_parameters():
    rng = np.random.default_rng(5)
    s = _truth_sample(rng)
    truth = dict(vsys=-200.0, V=15.0, sigma0=20.0, beta=-0.1, f=0.05, R0=1.0, Rstar=2.5)
    s["v"] = make_mock(s, truth, slope=0.0, k=1.0, rng=rng).v
    b = blind_fit(s, R0=1.0, Rstar=2.5)
    assert abs(b["vsys"] + 200) < 1.0 and abs(b["V"] - 15) < 1.5 and abs(b["sigma0"] - 20) < 1.0
    assert abs(b["beta"] + 0.1) < 0.05 and abs(b["f"] - 0.05) < 0.02


def test_tercile_contrast_follows_slope():
    """With exact [Fe/H] (k -> 0) the tercile contrast in V/sigma equals slope x (mean z difference of the terciles)."""
    rng = np.random.default_rng(6)
    s = _truth_sample(rng, 6000)
    truth = dict(vsys=-200.0, V=15.0, sigma0=20.0, beta=0.0, f=0.0, R0=1.0, Rstar=2.5)
    s["efeh"] = 1e-6
    m = make_mock(s, truth, slope=0.5, k=1.0, rng=rng)
    t = tercile_fit(m, blind_fit(m, R0=1.0, Rstar=2.5))
    expect = 0.5 * (t["zmean"][2] - t["zmean"][0])
    assert abs(t["delta"] - expect) < 0.12


def test_statistics_runs_on_real_sample():
    st = statistics(sample("NGC 185"), R0=1.0, Rstar=2.5, nperm=200, rng=np.random.default_rng(0))
    assert set(st) >= {"p_rot", "p_sig", "delta", "blind"}


@pytest.mark.parametrize("galaxy", ["Bootes I", "Draco"])
def test_dsph_sample(galaxy):
    from olddisks.defit import sample_dsph
    s = sample_dsph(galaxy)
    assert len(s) > 30 and (s.R < 6).all()
    # no star twice (1 arcsec), and perspective-corrected velocities centred on v_sys
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    c = SkyCoord(s.ra.to_numpy() * u.deg, s.dec.to_numpy() * u.deg)
    assert (c.match_to_catalog_sky(c, nthneighbor=2)[1].arcsec > 1).all()


def test_mock_slope_sign_independent_of_rotation_sense():
    """A positive slope means more rotational support for metal-rich stars whatever the sense of rotation (sign of V)."""
    rng = np.random.default_rng(8)
    s = _truth_sample(rng, 6000)
    s["efeh"] = 1e-6
    for V in (15.0, -15.0):
        truth = dict(vsys=-200.0, V=V, sigma0=20.0, beta=0.0, f=0.0, R0=1.0, Rstar=2.5)
        m = make_mock(s, truth, slope=0.5, k=1.0, rng=rng)
        t = tercile_fit(m, blind_fit(m, R0=1.0, Rstar=2.5))
        assert abs(t["delta"] - 0.5 * (t["zmean"][2] - t["zmean"][0])) < 0.12


def test_mock_support_never_counter_rotates():
    """For a non-rotator (V = 0) a positive slope gives rotating metal-rich and non-rotating (not counter-rotating)
    metal-poor stars, and the tercile contrast sees it."""
    rng = np.random.default_rng(9)
    s = _truth_sample(rng, 6000)
    s["efeh"] = 1e-6
    truth = dict(vsys=0.0, V=0.0, sigma0=10.0, beta=0.0, f=0.0, R0=1.0, Rstar=2.5)
    m = make_mock(s, truth, slope=1.0, k=1.0, rng=rng)
    t = tercile_fit(m, blind_fit(m, R0=1.0, Rstar=2.5, fit_theta=True))
    assert abs(t["ratio"][0]) < 0.1 and t["delta"] > 0.4
