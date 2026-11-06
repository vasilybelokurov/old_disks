"""Rotation versus [Fe/H] in the M31 dwarf ellipticals NGC 147 and NGC 185 (agreed plan: docs/plans/de_model.md, bottom).

sample(): all Geha+2010 velocity members, Ho+2015 [Fe/H] attached where matched (catalogues.ho2015); x (major-axis coordinate)
    and R (elliptical radius) in R_half along the LVDB photometric axes.
blind_fit(): [Fe/H]-blind maximum likelihood fit to all members,
    v ~ (1 - f) N(vsys + V h(R) x/R, sigma(R)^2 + e_v^2) + f N(vsys, S_OUT^2 + e_v^2),
    sigma(R) = sigma0 exp(beta (R - Rstar)), h(R) = tanh(R/R0) / tanh(Rstar/R0) (R0 = None: h = R/Rstar), so that
    V and sigma0 are the rotation and dispersion at Rstar.
permutation_test(): mean over equal-count radial bins of the Spearman rho of y with [Fe/H]; p from permutations of [Fe/H]
    within the bins (keeps the radial [Fe/H] gradient).
statistics(): T_rot (y = (v - vsys) sign(x) sign(V)) and T_sigma (y = |residual| / sqrt(sigma(R)^2 + e_v^2)) with their p,
    and the tercile contrast of V/sigma at Rstar (tercile_fit).
make_mock(): velocities (and observed [Fe/H]) at the real positions and errors for a slope d(V/sigma)/d[Fe/H].
"""
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import rankdata

from olddisks.catalogues import ho2015, table
from olddisks.geometry import tangent_plane
from olddisks.paths import LVDB

EV_FLOOR = 1.5  # km/s; 41 Geha+2010 rows have e_HRV = 0
S_OUT = 100.0  # km/s, width of the contaminant component
NBINS = 5


def sample(galaxy, dpa=0.0):
    """dpa [deg]: offset added to the LVDB position angle (systematics test)."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    g = pd.read_csv(LVDB).drop_duplicates("name").set_index("name").loc[galaxy]
    t = table("J/ApJ/711/361/veloc")
    t = t[t.NGC == int(galaxy.split()[1])].reset_index(drop=True)
    ho = ho2015()
    ho = ho[ho.galaxy == galaxy].set_index("star_id")
    c = SkyCoord(t.RAJ2000.astype(str), t.DEJ2000.astype(str), unit=(u.hourangle, u.deg))
    xi, eta = tangent_plane(c.ra.deg, c.dec.deg, g.ra, g.dec)
    p = np.radians(g.position_angle + dpa)
    xmaj = (xi * np.sin(p) + eta * np.cos(p)) * 60 / g.rhalf
    xmin = (xi * np.cos(p) - eta * np.sin(p)) * 60 / g.rhalf
    name = t.Name.astype(str)
    return pd.DataFrame(dict(star_id=name, x=xmaj, R=np.hypot(xmaj, xmin / (1 - g.ellipticity)), v=t.HRV.to_numpy(float),
                             ev=np.maximum(t.e_HRV.to_numpy(float), EV_FLOOR), imag=t.Imag.to_numpy(float),
                             feh=name.map(ho.feh).to_numpy(float), efeh=name.map(ho.efeh).to_numpy(float)))


def h(R, R0, Rstar):
    return R / Rstar if R0 is None else np.tanh(R / R0) / np.tanh(Rstar / R0)


def _components(s, q, R0, Rstar):
    """Mean and member variance per star for q = (vsys, V, ln sigma0, beta)."""
    vsys, V, lns, beta = q
    mean = vsys + V * h(s.R.to_numpy(), R0, Rstar) * s.x.to_numpy() / s.R.to_numpy()
    sig = np.exp(lns + beta * (s.R.to_numpy() - Rstar))
    return mean, sig**2 + s.ev.to_numpy() ** 2


def _loglike(s, vsys, mean, var, f):
    v, ev2 = s.v.to_numpy(), s.ev.to_numpy() ** 2
    lin = -0.5 * (v - mean) ** 2 / var - 0.5 * np.log(2 * np.pi * var)
    vo = S_OUT**2 + ev2
    lo = -0.5 * (v - vsys) ** 2 / vo - 0.5 * np.log(2 * np.pi * vo)
    return np.logaddexp(np.log1p(-f) + lin, np.log(f) + lo)


def blind_fit(s, R0, Rstar):
    """Maximum likelihood (vsys, V, sigma0, beta, f) over all stars of s; best of several starts."""
    def nll(p):
        mean, var = _components(s, p[:4], R0, Rstar)
        return -np.sum(_loglike(s, p[0], mean, var, 1 / (1 + np.exp(-p[4]))))
    v0 = np.median(s.v)
    starts = [[v0, V, np.log(np.std(s.v)), 0.0, -3.0] for V in (-30, -10, 10, 30)]
    best = min((minimize(nll, x, method="L-BFGS-B", bounds=[(v0 - 50, v0 + 50), (-100, 100), (0, 5), (-2, 2), (-10, 0)])
                for x in starts), key=lambda r: r.fun)
    vsys, V, lns, beta, lf = best.x
    return dict(vsys=vsys, V=V, sigma0=np.exp(lns), beta=beta, f=1 / (1 + np.exp(-lf)), R0=R0, Rstar=Rstar, nll=best.fun)


def _rbins(R, nbins=NBINS):
    edges = np.quantile(R, np.linspace(0, 1, nbins + 1))
    k = np.clip(np.searchsorted(edges, R, side="right") - 1, 0, nbins - 1)
    return [np.flatnonzero(k == i) for i in range(nbins)]


def _binned_rho(y, z, bins):
    """Mean over bins of the Spearman rho of y and z (z of shape (..., n) for permutations)."""
    out = []
    for b in bins:
        ry = rankdata(y[b])
        rz = rankdata(z[..., b], axis=-1)
        ry, rz = ry - ry.mean(), rz - rz.mean(axis=-1, keepdims=True)
        out.append((rz * ry).sum(-1) / np.sqrt((rz**2).sum(-1) * (ry**2).sum()))
    return np.mean(out, axis=0)


def permutation_test(y, R, feh, nperm=10000, rng=None):
    """Statistic and two-sided permutation p (stars with finite [Fe/H] only)."""
    rng = np.random.default_rng() if rng is None else rng
    ok = np.isfinite(feh)
    y, R, feh = np.asarray(y)[ok], np.asarray(R)[ok], np.asarray(feh)[ok]
    bins = _rbins(R)
    t = _binned_rho(y, feh, bins)
    perm = np.tile(feh, (nperm, 1))
    for b in bins:
        perm[:, b] = rng.permuted(perm[:, b], axis=1)
    null = _binned_rho(y, perm, bins)
    return dict(t=t, p=(np.sum(np.abs(null) >= abs(t)) + 1) / (nperm + 1), null=null)


def tercile_fit(s, b):
    """V/sigma at Rstar for the observed-[Fe/H] terciles (assigned within radial bins); vsys, beta, f, R0 fixed from b."""
    s = s[np.isfinite(s.feh)].reset_index(drop=True)
    terc = np.zeros(len(s), int)
    for idx in _rbins(s.R.to_numpy()):
        terc[idx] = np.floor(3 * (rankdata(s.feh.to_numpy()[idx]) - 0.5) / len(idx)).astype(int)
    ratio, zmean = [], []
    for t in range(3):
        st = s[terc == t]

        def nll(q):
            mean, var = _components(st, (b["vsys"], q[0], q[1], b["beta"]), b["R0"], b["Rstar"])
            return -np.sum(_loglike(st, b["vsys"], mean, var, b["f"]))
        r = min((minimize(nll, [V, np.log(b["sigma0"])], method="L-BFGS-B", bounds=[(-100, 100), (0, 5)])
                 for V in (b["V"], -b["V"], 0.0)), key=lambda r: r.fun)
        ratio.append(np.sign(b["V"]) * r.x[0] / np.exp(r.x[1]))
        zmean.append(st.feh.mean())
    return dict(ratio=ratio, zmean=zmean, delta=ratio[2] - ratio[0])


def statistics(s, R0, Rstar, nperm=10000, rng=None):
    b = blind_fit(s, R0, Rstar)
    mean, var = _components(s, (b["vsys"], b["V"], np.log(b["sigma0"]), b["beta"]), R0, Rstar)
    y_rot = (s.v - b["vsys"]) * np.sign(s.x) * np.sign(b["V"])
    y_sig = np.abs(s.v - mean) / np.sqrt(var)
    rot = permutation_test(y_rot, s.R, s.feh, nperm, rng)
    sig = permutation_test(y_sig, s.R, s.feh, nperm, rng)
    t = tercile_fit(s, b)
    return dict(t_rot=rot["t"], p_rot=rot["p"], t_sig=sig["t"], p_sig=sig["p"], delta=t["delta"], ratio=t["ratio"],
                zmean=t["zmean"], blind=b)


def make_mock(s, truth, slope, k, rng):
    """Mock of sample s: true [Fe/H] z drawn given the observed one (Gaussian shrinkage about the radial trend, with errors
    k e_z); V(z) = V + slope sigma0 (z - mean z), sigma independent of z; a fraction f of contaminants N(vsys, S_OUT);
    returned with mock velocities and mock observed [Fe/H] (NaN kept where the real star has none)."""
    s = s.copy()
    tag = np.isfinite(s.feh).to_numpy()
    R, feh, e = s.R.to_numpy(), s.feh.to_numpy(), k * s.efeh.to_numpy()
    a, c = np.polyfit(R[tag], feh[tag], 1)
    m = a * R + c
    tau2 = max(np.var(feh[tag] - m[tag]) - np.mean(e[tag] ** 2), 0.01 * np.var(feh[tag]))
    z = m + rng.normal(0, np.sqrt(tau2), len(s))
    lam = tau2 / (tau2 + e[tag] ** 2)
    z[tag] = m[tag] + lam * (feh[tag] - m[tag]) + rng.normal(0, np.sqrt(lam) * e[tag])
    V = truth["V"] + slope * truth["sigma0"] * (z - z[tag].mean())
    mean = truth["vsys"] + V * h(R, truth["R0"], truth["Rstar"]) * s.x.to_numpy() / R
    sig = truth["sigma0"] * np.exp(truth["beta"] * (R - truth["Rstar"]))
    v = mean + rng.normal(0, np.sqrt(sig**2 + s.ev.to_numpy() ** 2))
    out = rng.uniform(size=len(s)) < truth["f"]
    v[out] = truth["vsys"] + rng.normal(0, np.sqrt(S_OUT**2 + s.ev.to_numpy()[out] ** 2))
    s["v"] = v
    s.loc[tag, "feh"] = z[tag] + rng.normal(0, e[tag])
    s["ztrue"] = z
    return s
