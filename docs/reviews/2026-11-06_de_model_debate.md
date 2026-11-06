# Adversarial review of the dE model plan (2026-11-06)

Critic and proposer (both Opus), critic signs off in round 2. Scratch scripts are in the session scratchpad
(de_critic/, de_proposer/).

## Round 1, critic
1. **Power.** σ(d(V/σ)/d[Fe/H]) = 0.52/dex (NGC 147), 0.21/dex (NGC 185). Minimum detectable slope ≈ 1.0 and 0.45 per dex.
2. **Deconvolved amplitude is set by the e_z scale.** Attenuation λ = 0.30 (NGC 147) and 0.51 (NGC 185); ±10 % in e_z
   changes the amplitude by −30 % to +92 %.
3. **Stars without [Fe/H] are not missing at random:** they are more central (KS p ~ 1e−5), and NGC 185 has a gradient of
   −0.23 dex per R_half.
4. **R0 is unconstrained;** report V/σ inside the data.
5. **NGC 205 is unusable (FATAL).**
6. **The Geha reproduction is not a test of the [Fe/H] part.**
7. **Use a model-free permutation statistic.**
8. **sin i** cancels within a galaxy; not for the LMC comparison.
9. **CaT [Fe/H] near the TRGB** (AGB, intermediate age).
10. Binaries negligible (Δσ ≈ 0.04 km/s); floor e_v.
11. **The model is effectively v versus major-axis distance** (median |x|/R ≈ 0.93).

## Round 1, proposer
- Conceded 1, 3–7, 9–11.
- On 2: the error-scale fit gives k = 0.99 (NGC 185) and k = 0.48 (NGC 147). In NGC 147 the [Fe/H] variance is flat
  across e_z terciles, and the brightest stars have the largest spread: a CaT systematic.
- Two-Gaussian MDF preferred (Δnll 15.4 and 9.8).

## Round 2, critic: sign-off with must-fix items
- Within-R permutation nulls are calibrated at the real positions, errors and [Fe/H]: FPR 0.03–0.05 at p < 0.05.
- **Must-fix:**
  1. Do not fit k; report k = 0.5 and k = 1.
  2. Normalise T_σ.
  3. Multiple testing: pre-register the primary test.
- **Simplify:** drop the joint MDF + k + quadrature model; use mock-calibrated tercile contrasts.
