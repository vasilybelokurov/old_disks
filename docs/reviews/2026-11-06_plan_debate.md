# Adversarial review of the dIrr plan (2026-11-06)

Two Opus agents with read access to the repository: a critic and a proposer, two rounds. Scratch mocks are in the
session scratchpad (critic/m.py, proposer/p.py); they are not part of the repository.

## Round 1, critic: verdict "abandon as written"
1. **FATAL.** Belland+2020 uses Kirby+2013's 299 NGC 6822 stars (their §2), so it adds no stars.
2. **FATAL.** The [Fe/H]-split rotation is already published for NGC 6822 (Belland+2020), IC 1613 (Taibi+2024 §6.1) and
   WLM (Leaman+2012).
3. **FATAL.** Under published-only, no dIrr reaches the threshold: IC 1613 ≤ 414, NGC 6822 380, WLM 180. Re-reduction adds
   no stars.
4. **MAJOR.** Power for NGC 6822 falls to 0.45 with 5 % contaminants and to 0.16 for a Belland-like effect; about 1000
   stars would be needed.
5. **MAJOR.** The [Fe/H] halves differ in source mix (Kirby, Swan, BOSS), which biases the |V| estimator.
6. **MAJOR.** [Fe/H] is a weak age proxy.
7. **MAJOR.** The LMC comparison is not like for like: 3D versus line of sight, inner radii, sin i, prolate rotation in
   NGC 6822.
8. **Minor/major.** Walker+2023 membership is by velocity window only; Sculptor is already analysed (Arroyo-Polonio+2024).
9. Novelty of a multi-galaxy study unverified.

## Round 1, proposer
- Conceded 1, 2, 4–6, 8.
- Partial on 3 and 7.
- New: Ho+2015 (2015ApJ...798...77H) gives CaT [Fe/H] for 1243 DEIMOS red giants in NGC 147, 185, 205, And II,
  And VII and M32 (Table 1: 230, 321, 224, 300, 104, 64), from published Geha programmes.
- Mixture-fit mocks: NGC 6822 power 0.76 clean, 0.43–0.52 with contaminants, 0.15 for a Belland-like effect.
- Proposed a hierarchical stack of Δ(V/σ) over about 10 galaxies.

## Round 2, critic: verdict "proceed with changes"
- Ho+2015 verified (also checked independently in ADS by the main session).
- [Fe/H] errors ≥ 0.4 dex dilute a median split to 0.48–0.58 of a step Δ (0.69–0.79 of a linear one). Fit the
  continuous dependence with errors deconvolved.
- Ho+2015 selection depends on [Fe/H]: metal-rich stars have the lowest S/N, and 478 stars were dropped at S/N < 3.
  Model it.
- Velocities for NGC 205 (Geha+2006) and And II (Ho+2012) are not in hand.
- Novelty: And II is already split (Ho+2012); Leaman+2017 compares age–σ across eight galaxies. NGC 147/185/205 [Fe/H]-split
  kinematics appear unpublished (26 citing papers checked).
- MW dSph gradients may be tidal: keep them out of the main stack.
- A single absolute Δ across dSphs and dEs is not meaningful: use a fractional trend, and stack only detected rotators.

Outcome: `docs/plans/rotation_vs_feh.md` (plan v3).
