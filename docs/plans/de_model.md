# Plan B: kinematic model for the M31 dEs, adapted from the LMC model (2026-11-06, draft for adversarial review)

## Summary

- **Model:** the LMC forward model (`../lmc-old-population/src/lmcold/kinfit.py`) reduced to what line-of-sight
  velocities can constrain. Mean velocity v_sys + V_p tanh(R/R0) x/R, a dispersion σ, and a broad contaminant
  component. V_p = V sin i is the projected amplitude, because inclination is not measurable from v_los alone.
- **Metallicity enters as in the LMC (populations in [Fe/H]), but deconvolved.**
  - Primary: V_p and σ linear in the true [Fe/H] z, with z marginalised per star given [Fe/H]_obs, e_[Fe/H] and a
    Gaussian-mixture prior on the metallicity distribution (MDF).
  - Secondary (the LMC-like presentation): two or three bins in observed [Fe/H], with the dilution quantified by mocks.
- **Quantity reported:** V_p/σ at R = 2 R_half for metal-poor (z = MDF 16th percentile) and metal-rich (84th percentile)
  stars, and d(V_p/σ)/d[Fe/H]. Compared with the LMC only after projecting the LMC model to v_los at a matched radius.
- **Inference:** maximum likelihood (JAX gradients, L-BFGS-B) for the mocks; NUTS (numpyro) for the real fits, with flat
  priors in physical bounds, as in the LMC.
- **No real fit before the mocks pass.** Mocks use the real positions, errors, contaminants, binaries and selection.

## 1. What carries over from the LMC model and what does not

| LMC (`kinfit.py`) | M31 dEs | Reason |
|---|---|---|
| Disc plane: inclination, line of nodes, exact 3D projection | Projected model: R = elliptical radius (photometric PA, ellipticity), x = major-axis coordinate, cos φ_disc = x/R | No proper motions, so the inclination is degenerate with V; fit V_p = V sin i. An oblate rotator seen at inclination i has exactly this form when the ellipticity traces i. |
| Rotation v0 tanh(R/R0) | V_p tanh(R/R0); R0 shared by all metallicities | Same functional form. Geha+2010 see flattening only in the outer parts, so R0 needs a broad prior. |
| σ and σ_los | One σ (line of sight), with an optional radial slope dlnσ/dR | Geha+2010: σ "decreases gently with radius". |
| Outlier mixture (f, S_OUT) | Mixture with a broad contaminant: free mean and width, Gaussian prior centred on the M31 halo (NGC 205) or the MW foreground (NGC 147, 185) | NGC 205 sits on M31's halo (v_sys −246 vs M31 −300 km/s). |
| Populations per [Fe/H] bin sharing geometry and v_sys | V_p(z) = V_0 + V_1 (z − z_ref), σ(z) = σ_0 + σ_1 (z − z_ref); geometry, v_sys, R0, contaminant shared | Errors of ~0.5 dex make bins leak (dilution 0.5–0.6 for a step). Deconvolution is needed. |
| Stars with no proper motion use v_los only | Stars with no [Fe/H] (Geha members missing from Ho+2015) use z marginalised over the MDF prior | Handles the Ho+2015 S/N cut (metal-rich stars have lower S/N) without discarding the velocities. |
| Perspective (exact projection) | Not needed | At 650–830 kpc, the perspective gradient from a proper motion of about 0.05 mas/yr is < 0.5 km/s across the field. To be checked against the measured M31-satellite proper motions. |
| Centre and PA fitted | Centre and PA fixed to LVDB; PA free with a prior as a robustness test | The DEIMOS masks follow a strip along the major axis; the minor axis is not sampled. |

## 2. Likelihood per star

For star i with v_i, e_v,i, position (x_i, R_i) and, if available, [Fe/H]_obs,i ± e_z,i:

L_i = (1 − f) ∫ N(v_i | μ(x_i, R_i, z), σ(z, R_i)² + e_v,i² + s_bin²) p(z | z_obs,i) dz  +  f N(v_i | μ_c, σ_c² + e_v,i²)

- μ = v_sys + V_p(z) tanh(R/R0) x/R.
- p(z | z_obs) ∝ N(z_obs | z, e_z²) × MDF(z). The MDF is a 2-component Gaussian mixture fitted to the galaxy's observed
  [Fe/H] with errors deconvolved (extreme deconvolution), fixed during the kinematic fit.
- The integral uses Gauss–Hermite quadrature (20 nodes).
- s_bin is a binary jitter (single-epoch data), with a prior from the binary fraction (default fixed at 1 km/s;
  sensitivity test with 3 km/s).
- Stars without [Fe/H] use p(z) = MDF.

Parameters per galaxy: v_sys, V_0, V_1, R0, σ_0, σ_1, dlnσ/dR, f, μ_c, σ_c (10).

## 3. Mocks and tests (acceptance criteria)

1. **Unit tests.**
   - Likelihood equals the simple Gaussian when e_z → 0 and f = 0.
   - Quadrature equals a brute-force integral to 1e-6.
   - Recovery at N = 5000: all parameters within 3σ.
2. **Mocks at the real positions, velocity errors and [Fe/H] errors** of NGC 147, NGC 185 and NGC 205. True z drawn from
   the deconvolved MDF; observed z = true + error; the Ho+2015 selection is reproduced by dropping [Fe/H] for the same
   fraction of stars as a function of magnitude. 5 % contaminants, 10 % binaries (velocity offsets from a log-normal
   period distribution, about 2–5 km/s rms for red giants). 200 mocks per scenario. Scenarios:
   - null: V_1 = σ_1 = 0;
   - Geha-like: V_0 and σ_0 at the Geha+2010 values, V_1 = 0;
   - trend: V_p/σ changes by 0.3 per dex.

   Acceptance:
   - bias < 0.3 σ_stat for V_1, σ_1 and V_p/σ at 2 R_half;
   - null false-positive rate ≤ 5 %;
   - the posterior σ of d(V_p/σ)/d[Fe/H] is reported, setting the detectable effect size.
3. **Reproduction of Geha+2010** ([2010ApJ...711..361G](https://ui.adsabs.harvard.edu/abs/2010ApJ...711..361G)) using
   all Geha members (520 and 442), ignoring [Fe/H]:
   - NGC 147: V_max 17 ± 2, σ 16 ± 1 km/s;
   - NGC 185: V_max 15 ± 5, σ 24 ± 1 km/s;
   - agreement within 2σ, allowing that their V_max is a profile maximum, not V_p of a tanh.
4. **Binned (LMC-style) fit** with 2–3 observed-[Fe/H] bins, compared with the continuous fit through the mock-calibrated
   dilution.
5. **Robustness:** PA free; dlnσ/dR = 0; Kirby-like velocity zero-point offsets (n/a: single instrument); NGC 205 with and
   without a spatial cut against M31's disc.

## 4. Order of work

1. Install JAX/numpyro.
2. `src/olddisks/defit.py` (likelihood, ML, NUTS) and `tests/test_defit.py`.
3. `scripts/de_mocks.py` → `results/de_mocks.csv`, `figures/de_mocks.png`.
4. `scripts/fit_de_geha.py`: reproduction of Geha+2010 → `figures/de_geha_repro.png`.
5. `scripts/fit_de.py`: real fits → `results/de_fits.csv`, `figures/de_fit_<gal>.png` (data, model and mocks).
6. Only then: the verdict on V/σ versus [Fe/H] for each galaxy.

## Open questions for the review

- Is V_p (= V sin i) comparison across metallicities valid if metal-poor and metal-rich stars have different intrinsic
  shapes (and hence different effective inclinations)?
- Linear V_p(z) versus a two-population (discrete) model: which is better suited to errors of ~0.5 dex?
- Is fixing the MDF before the kinematic fit acceptable, or should it be fitted jointly?
- How to handle the Ho+2015 S/N selection more rigorously than "missing at random given magnitude"?
- Is NGC 205 salvageable (M31 contamination, tides, 11.5 km/s errors)?
