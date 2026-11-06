# Plan v3: rotational support versus [Fe/H] in Local Group dwarfs (2026-11-06)

Agreed after an adversarial review (critic and proposer, two rounds; record in
`docs/reviews/2026-11-06_plan_debate.md`). Supersedes `archival_rereduction.md` (the dIrr re-reduction is dropped).
Constraints: published programmes only; external datasets are the core.

## Summary

- **Primary sample: M31 dwarf ellipticals NGC 147, NGC 185 and NGC 205** (and And VII).
  - [Fe/H]: Ho+2015 (Ca II triplet, Carrera+2013 calibration, Keck/DEIMOS;
    [2015ApJ...798...77H](https://ui.adsabs.harvard.edu/abs/2015ApJ...798...77H)).
  - Velocities: Geha+2010 ([2010ApJ...711..361G](https://ui.adsabs.harvard.edu/abs/2010ApJ...711..361G)) and Geha+2006.
  - These galaxies rotate (V/σ ≈ 0.6–1) and the spectra reach 5–8 R_eff.
  - **No [Fe/H]-split kinematics found in ADS**, including the 26 papers citing Ho+2015 that mention kinematics.
- **Consistency checks: NGC 6822 and IC 1613, Kirby+2014 stars only.** Must reproduce Belland+2020 and Taibi+2024.
- **Separate test, not in the stack: MW dSphs.** Their velocity gradients may be tidal (Martínez-García+2023,
  [2023MNRAS.518.3083M](https://ui.adsabs.harvard.edu/abs/2023MNRAS.518.3083M)).
- **Method: fit V/σ as a continuous function of true [Fe/H]**, deconvolving each star's [Fe/H] error. No median split:
  with errors of ≥ 0.4 dex, a median split keeps only 0.48–0.58 of a step-like Δ.
- **Stack only systems with clearly detected rotation**, using the fractional change d ln(V/σ)/d[Fe/H].
- **Dropped:** dIrr re-reduction (adds no stars under published-only), Belland+2020 as a source (the same 299 Kirby
  stars), WLM (no table), and any novelty claim for And II (Ho+2012 already did it).

## Sample (counts to be confirmed after cross-matching)

| Galaxy | [Fe/H] source (N) | Velocity source | Status |
|---|---|---|---|
| NGC 185 | Ho+2015 (321) | Geha+2010 VizieR J/ApJ/711/361 (442 rows) | velocities on disk; Ho+2015 Table 2 to fetch |
| NGC 147 | Ho+2015 (230) | Geha+2010 (520 rows) | as above |
| NGC 205 | Ho+2015 (224) | Geha+2006 (to fetch) | inventory has only 31 stars |
| And VII | Ho+2015 (104) | Tollerud+2012 / Kirby+2020 | on disk |
| And II | Ho+2015 (300) | Ho+2012 (not on mirror) | context only, not novel |
| NGC 6822 | Kirby+2013 (274) | Kirby+2014 | on disk |
| IC 1613 | Kirby+2013 (124) | Kirby+2014 | on disk |
| MW dSphs | Tolstoy+2023 (Scl), Walker+2023 (others) | same | on disk; Walker needs a Gaia PM cut |

## Steps, tests and acceptance criteria

1. **Data.**
   - Fetch Ho+2015 Table 2 (MRT or CDS) and the Geha+2006 velocities.
   - Cross-match by position.
   - Record which Geha members lack [Fe/H]: 478 stars were dropped at S/N < 3, and metal-rich stars have the lowest S/N.
   - **Test:** the cross-match reproduces the Ho+2015 Table 1 counts. Compare the radial and magnitude distributions of
     kept and dropped stars (figure).
2. **Model** (`src/olddisks/kinematics.py`, extended).
   - Mixture likelihood: members plus a broad contaminant component; binary/jitter term; repeat epochs combined.
   - Rotation curve that flattens at r_t (prior on r_t); θ free.
   - V and σ linear in true [Fe/H], with [Fe/H] errors marginalised.
   - Stars without [Fe/H] constrain the mean kinematics through the mixture.
   - Perspective-rotation correction for the MW dSphs.
3. **Mocks before any real number.** Mocks use the real positions, velocity errors, [Fe/H] errors and the S/N-driven
   selection, with 5 % contaminants and 10 % binaries. Acceptance:
   - bias < 0.3 σ_stat on V, σ and dV/d[Fe/H];
   - null false-positive rate ≤ 5 %;
   - power for the Belland-like effect size (V/σ 0.3 → 0.5) reported per galaxy;
   - stack σ(μ) ≤ 0.10, otherwise report upper limits.
4. **Reproduction tests.**
   - Geha+2010 V_max and σ for NGC 147 and NGC 185.
   - Belland+2020 quintile trend for NGC 6822.
   - Arroyo-Polonio+2024 Sculptor gradient 4.0 ± 1.5 km/s/deg.
   - Split-by-source nulls (Kirby vs Swan vs BOSS in NGC 6822; Walker vs Tolstoy in Sculptor) consistent with 0 at 2σ.
5. **Real data, then the stack.**
   - Per galaxy: figures of data, model and mocks.
   - Stack the fractional trend over systems with rotation detected at > 3σ.
   - Report the MW dSphs separately, against orbital phase.
6. **LMC comparison.** Only after projecting the LMC model to V_los/σ_los at the same radius in R_half.

## Prior work to cite

Belland+2020 (NGC 6822), Taibi+2024 (IC 1613), Leaman+2012 (WLM), Ho+2012 (And II), Leaman+2017 (age–σ, eight
galaxies), Wheeler+2017 (global v/σ), Arroyo-Polonio+2024 (Sculptor), Martínez-García+2021/2023 (dSph rotation, tides).

## Open risks

- **[Fe/H] errors ≥ 0.4 dex at S/N ≈ 6** dilute any trend; deconvolution recovers it only if the error model is right.
- **M31 tides** on NGC 147 and NGC 205.
- **Fetching tables:** depends on CDS or the journal MRTs.
