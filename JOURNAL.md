# JOURNAL: rotational support versus metallicity in Local Group dwarfs

Rules: `AGENTS.md`.

---

## 2026-11-06: Session 1, start of the inventory

### Literature checks (ADS)
- WLM: Leaman et al. 2012, ApJ 750, 33 (https://ui.adsabs.harvard.edu/abs/2012ApJ...750...33L). Abstract confirms 180 RGB
  stars, FORS2 + DEIMOS CaT, v_rot = 17 ± 1 km/s, V/σ ~ 1, intrinsically thick (0.39 ≤ q0 ≤ 0.57), σ rising with age.
- Isolated LG dwarfs: Kirby et al. 2014, MNRAS 439, 1015 (https://ui.adsabs.harvard.edu/abs/2014MNRAS.439.1015K). The
  brief's count (862 stars in seven galaxies) has not yet been checked against the paper.
- Galaxy parameters: Local Volume Database v1.1.2, Pace 2025, OJAp 8, 142 (https://ui.adsabs.harvard.edu/abs/2025OJAp....8E.142P).

### Setup
- git repository, package `src/olddisks`, project `.venv` (Python 3.13, numpy 2.5.4, pandas 3.1.1, sqlutilpy 0.30.0, astropy
  8.0.2, pyarrow), xattr `com.dropbox.ignored` set. `pytest`: 1 test (elliptical radius) passes.

### WSDB survey census (`scripts/build_survey_census.py`, figure `figures/survey_census.png`)
Sample: LVDB with confirmed_galaxy, D < 1.5 Mpc, M_V < −7 and a known v_sys, excluding LMC, SMC, Sgr and M33 → 60 galaxies.
Candidate = R_ell < 3 R_half and |v − v_sys| < 3√(σ² + e_v²) (σ = 10 km/s if missing). No PM or CMD cut.

Velocity candidates with pipeline [Fe/H] (n_v_feh / n_v), from `data/survey_census.csv`:

| galaxy | survey | n_v_feh / n_v | median e_v [km/s] |
|---|---|---|---|
| Sextans | DESI DR1 | 516 / 516 | 5.7 |
| Fornax | APOGEE DR17 | 217 / 250 | 0.15 |
| Draco | DESI DR1 | 240 / 240 | 2.2 |
| Sculptor | APOGEE DR17 | 158 / 173 | 0.45 |
| Sextans | APOGEE DR17 | 56 / 80 | 0.8 |
| Ursa Minor | APOGEE DR17 | 53 / 69 | 0.5 |
| Draco | APOGEE DR17 | 51 / 64 | 0.5 |
| Carina | APOGEE DR17 | 54 / 60 | 0.3 |
| NGC 6822 | SDSS-V BOSS | 34 / 34 | 6.2 (contaminated, below) |
| NGC 205 | DESI DR1 | 28 / 28 | 2.3 |
| Sextans | SDSS-V BOSS | 27 / 27 | 8.4 |
| CVn I | DESI DR1 | 25 / 25 | 4.5 |

All other galaxies have < 15 candidates per survey. The isolated dIrrs (WLM, IC 1613, Leo A, Pegasus, Aquarius, Cetus,
Phoenix, Tucana) have ≤ 13 candidates in total. **For them, the published literature catalogues are the only route.**

Seen in the figure:
- Classical dSphs (Sextans, Draco, Fornax, Sculptor, UMi, Carina): narrow peaks at v_sys well above the background.
- NGC 6822 (BOSS): candidates have [Fe/H] ≈ −0.5 to 0 and fill the 3σ window evenly. These are Milky Way disc
  foreground (b ≈ −18°, v_sys = −55 km/s). This sample is not usable as is.
- NGC 205 (DESI): [Fe/H] ≈ −1, positions concentrated on the M31 side; possible M31 contamination, not tested.
- APOGEE DR17 censors the metal-poor end: the minimum [Fe/H] among candidates is −2.47 (ASPCAP grid edge), and 115 of 706
  candidates have a velocity but no [Fe/H]. A rotation-versus-[Fe/H] analysis with APOGEE must model this; it biases the
  metal-poor bin.
- DESI RVS [Fe/H] extends to −4 (5th percentile −3.2) and is broad for Sextans; its precision for faint giants has not
  been checked (compare with literature [Fe/H] on overlapping stars).

### WSDB issue
`sdss5_dr20.allstar` (SDSS-V APOGEE) in WSDB has 389,782 rows, only 258 of them at dec < −20, and none within any dwarf
(including northern Draco and UMi). It looks incomplete; it is dropped from the census. Next: check with S. Koposov or
the SDSS DR20 documentation.

### Next
1. Literature census (subagent running) → `data/literature_census.csv`; check the counts in the papers.
2. Merge the WSDB and literature censuses into `data/inventory.csv`, with Gaia PM + CMD membership for the WSDB samples.
3. Choose a pilot: one classical dSph with the most stars with velocity and [Fe/H] (Sculptor/Fornax/Sextans) and WLM.

## 2026-11-06: Session 1 (continued), external datasets

User: external datasets are the core of the project; WSDB is secondary. Repository public at
https://github.com/vasilybelokurov/old_disks; commit and push after each tested step.

### Literature census (subagent, `data/literature_census.csv`, 60 rows)
Spot-checked in ADS: Belland+2020 (2020ApJ...903...10B), which already measures v and σ in five [Fe/H] bins in NGC 6822 and finds v/σ
rising with [Fe/H], so it is the direct precedent; Walker+2023 (2023ApJS..268...19W; 4492 likely members, 38 systems,
σ[Fe/H] ≈ 0.4 dex); Koposov+2026 DESI DR1 stellar VAC (2026OJAp....955260K); Sestito & Kobayashi 2026 (2026MNRAS.549g1137S),
whose finding that DESI [Fe/H] < −2 is unreliable matches the broad DESI [Fe/H] seen in `figures/survey_census.png`.

### VizieR download (`scripts/fetch_vizier.py`, `data/vizier_manifest.csv`)
CDS Strasbourg timed out (also outside the sandbox; user confirms CDS is down). Used the CfA mirror. A bug was found and fixed:
astroquery ignored `conf.server` and still went to CDS; the server must be passed as `Vizier(vizier_server=...)`.
26 catalogues downloaded. 8 are missing on the mirror, including the key isolated-dwarf ones (WLM Leaman+2012, NGC 6822
Belland+2020, IC 1613 Taibi+2024): see `data/sources.yaml`.

### Archive census (`scripts/build_archive_census.py`, `data/archive_census.csv`, `figures/archive_census.png`)
ESO raw + phase3 and KOA DEIMOS within min(3 R_half, 2°). Validation: recovers DART 71.B-0641 (Hill/Tolstoy) in Sculptor and
KOA C257D (PI Brooks, a Leaman+2012 co-author) in WLM.
- KOA DEIMOS: every M31 satellite and most isolated dwarfs; Draco and UMi 41 h each.
- ESO: classical southern dSphs dominate (Carina 260 h GIRAFFE + 307 h UVES; Sculptor, Fornax, Sextans 90–140 h GIRAFFE).
- Isolated dIrrs with substantial ESO time but no matching public per-star catalogue: IC 1613 (X-shooter 103 h,
  VIMOS 92 h, FORS2 59 h, GIRAFFE 30 h + 4809 phase3 GIRAFFE spectra), NGC 6822 (VIMOS 167 h, FORS2 128 h, UVES 80 h),
  NGC 3109, Sextans A (UVES 106 h), WLM (FORS2 46 h), Phoenix (FORS2 81 h), Tucana (FORS2 66 h).
  These are the candidates for re-reduction if the published catalogues are inadequate.
- ESO phase3 counts are reduced products (one per star for fibres), not telescope time.

### Next
Harmonise the 26 catalogues into one per-star table (galaxy, ra, dec, v, e_v, [Fe/H], e_[Fe/H], member flag, source) and
count members with v + [Fe/H] per galaxy, together with the radial coverage. Get the 8 missing catalogues from the journals (MRT) or from CDS when it is back.

### Harmonised inventory (`src/olddisks/catalogues.py`, `scripts/build_inventory.py`, `data/inventory.csv`)
24 published tables + WSDB surveys → one per-star table. Tests (`tests/test_catalogues.py`, 8 pass) reproduce the published member
counts: Tolstoy+2023 Sculptor 1604 members / 1339 with [Fe/H] (needs `q_[Fe/H] = 0`; the 265 with q = 1 have median S/N 11),
Kacharov+2017 Phoenix 194, Kirby+2014 792 of 862, Kirby+2020 M31 dSphs 256.
Membership: the paper's flag where given; otherwise the velocity window. Union across sources deduplicates at 1".

**Bug found from the figure and fixed (test first):** Ji+2021 table4 is Antlia 2, not Crater 2 (the paper covers both;
2021ApJ...921...32J). `figures/inventory_maps.png` showed 290 Crater II "members" outside 6 R_half. New test: the median
member position of every (source, galaxy) lies within 2 R_half of the LVDB centre; it failed (76 R_half), then passed after the fix.

Members with v + [Fe/H] (union), galaxies with ≥ 40 (`figures/inventory_counts.png`, `figures/inventory_maps.png`):

| galaxy           |   n_v |   n_v_feh |   n_beyond_2rh |   max_r_ell |   med_ev |   med_efeh |    M_V |
|:-----------------|------:|----------:|---------------:|------------:|---------:|-----------:|-------:|
| Sculptor         |  3145 |      2469 |            702 |        7.09 |     1.09 |       0.31 | -10.77 |
| Carina           |  2687 |      2318 |            321 |        6.37 |     1.18 |       0.58 |  -9.43 |
| Draco            |  2038 |      2037 |            602 |       16.48 |     0.89 |       0.4  |  -8.87 |
| Sextans          |  1798 |      1666 |            290 |        7.87 |     1.44 |       0.59 |  -8.72 |
| Ursa Minor       |  1839 |      1618 |            389 |       12.13 |     2.04 |       0.34 |  -8.86 |
| Fornax           |  3782 |      1326 |            567 |        4.83 |     1.26 |       0.37 | -13.39 |
| Leo I            |   528 |       528 |            161 |        4.26 |     1.1  |       0.37 | -11.82 |
| Leo II           |   416 |       416 |            155 |        6.99 |     1.04 |       0.43 |  -9.48 |
| NGC 6822         |   398 |       380 |             15 |        2.97 |     2.5  |       0.11 | -15.68 |
| Antlia II        |   353 |       352 |              2 |        2.5  |     2    |       0.29 |  -9.73 |
| Crater II        |   240 |       240 |             13 |        2.97 |     1.27 |       0.38 |  -8.18 |
| Bootes I         |   301 |       202 |             56 |      318.55 |     1.28 |       0.35 |  -6.02 |
| Phoenix          |   177 |       177 |             41 |        3.64 |     7.94 |       0.16 |  -9.87 |
| Tucana II        |   159 |       159 |              0 |        0.94 |     0.34 |       0.17 |  -3.75 |
| IC 1613          |   154 |       139 |              1 |        2.16 |     3.05 |       0.11 | -15.12 |
| Leo A            |   127 |       113 |             22 |        3.65 |     3.9  |       0.15 | -12.11 |
| Pegasus dIrr     |   105 |        96 |              6 |        2.8  |     2.8  |       0.11 | -12.14 |
| Andromeda VII    |   171 |        87 |             14 |        2.65 |     3.6  |       0.13 | -12.78 |
| Andromeda V      |   115 |        81 |             16 |        3.82 |     2.9  |       0.12 |  -9.3  |
| Aquarius         |    69 |        64 |              6 |        4.02 |     4.38 |       0.13 | -10.47 |
| Cetus            |   156 |        58 |             26 |        5.84 |     5.25 |       0.13 | -11.19 |
| Canes Venatici I |    56 |        56 |             15 |        3    |     1.58 |       0.63 |  -8.73 |
| Tucana           |    67 |        56 |             12 |        5    |     6    |       0.16 |  -9.54 |
| UGC 4879         |    67 |        54 |             26 |        7.35 |     4.1  |       0.14 | -12.23 |
| Sagittarius dIrr |    45 |        43 |             15 |        5.95 |     4.2  |       0.15 | -11.79 |
| Andromeda I      |    61 |        40 |             11 |        3    |     3.55 |       0.13 | -11.35 |

Caveats seen in the figures:
- Walker+2023 (the largest single contributor for Fornax, Sculptor, Carina, Draco, UMi, Sextans, Leo I, Leo II) has no
  membership column in VizieR, so its stars are selected only by the velocity window. Draco, UMi and Sextans show scattered
  "x" points beyond 3 R_half; Bootes I has max R_ell = 318 R_half. Next: add a Gaia proper-motion cut (Walker+2023 carries
  pmRA/pmDE) for window-selected stars.
- [Fe/H] scales are not homogeneous (CaT calibrations, spectral synthesis, Walker+2023 σ ≈ 0.4 dex, APOGEE floor at −2.5,
  DESI unreliable below −2). Offsets between sources must be measured on overlapping stars before any binning in [Fe/H].
- Velocity zero-points between sources not yet compared.
- No public per-star table yet for WLM (2 stars), Sextans A/B, NGC 3109, NGC 147/185 [Fe/H], or the IC 1613 MUSE sample.
  Missing on the mirror: Leaman+2012 (WLM), Belland+2020 (NGC 6822), Taibi+2024 (IC 1613).

### Feasibility for dIrrs and M31 satellites (`scripts/feasibility_mocks.py`, `results/feasibility.csv`, `figures/feasibility.png`)
Comparison: the LMC measurement (lmc-old-population `results/kin_support.csv`) uses APOGEE + BOSS giants (57,557 BOSS and
5,796 APOGEE members) with v_los and spectroscopic [Fe/H], plus Gaia PMs of the same giants: 2,168–28,909 per [Fe/H] bin,
v0/σ_los from 1.84 to 4.64 with ±0.03–0.05. (User correction: the LMC comparison is the APOGEE/BOSS giants, not the RR Lyrae.)
Model: `olddisks.kinematics` (v_sys + V min(r,1) cos(φ−θ), Gaussian σ convolved with errors; recovery test passes at N = 4000).
Mocks at real positions and errors, split at the median [Fe/H], 300 per scenario.
- Power to detect ΔV/σ = 0.5 (0.5 → 1.0) at 95%: Sculptor 1.00 (reference), NGC 6822 0.65, Leo A 0.39, Phoenix 0.37,
  Pegasus 0.35, And V 0.35, IC 1613 0.30, And VII 0.27, Aquarius 0.26, others 0.07–0.20.
- Null floor (true V = 0) of recovered V/σ: 0.20 (NGC 6822) to 0.65 (Tucana); Sculptor 0.06. With 20–90 stars per half,
  a non-rotating population looks like V/σ ≈ 0.3–0.6.
- Per-half scatter of V/σ: Sculptor 0.04–0.05 (LMC-like precision); NGC 6822 0.17; the other dIrrs and M31 satellites 0.2–0.7.
- NGC 6822 and IC 1613 samples sit inside ~1.1 R_half (90th percentile of r_ell), where solid-body rotation is small.
Scaling 1/√N from NGC 6822 (σ_V/σ ≈ 0.17 at 190 per half), 80% power for Δ = 0.5 needs σ(Δ) ≈ 0.2, i.e. about 270 stars
per half (≥ 550 members with [Fe/H]), and coverage to ≥ 2 R_half. Optimistic: perfect membership, homogeneous [Fe/H], known model form.

### Re-reduction plan (`docs/plans/archival_rereduction.md`)
Programme-level audit of ESO data (dp_tech, ins_mode, titles). This corrects the earlier census reading: NGC 6822 VIMOS 118 h is
imaging (Gallart), 25 h is HII-region IFU; IC 1613 VIMOS and GIRAFFE spectra are massive stars; Sextans A UVES is supergiants
or its GC. The old earlier recommendation to re-reduce these was wrong. The real red-giant CaT data are Cole 71.B-0516(B) FORS2 (IC 1613 17 h,
NGC 6822 19 h, 8 GC calibrators; only AAS abstracts), Held 075.D-0838 FORS2 (NGC 6822 18 h; nothing in ADS) and Mayer
072.B-0307 GIRAFFE (IC 1613 RGB 8.7 h; nothing in ADS). Validation set: Tolstoy 69.B-0305(B), the Swan+2016 masks.

### Published-only constraint and adversarial review (user, 2026-11-06)
User: only data from published programmes; asked for an adversarial review by two Opus agents. Result
(`docs/reviews/2026-11-06_plan_debate.md`): the dIrr re-reduction plan was abandoned. Belland+2020 = the same Kirby stars;
the per-galaxy results are already published; no dIrr can reach the threshold. New primary sample: M31 dEs NGC 147, 185, 205
with Ho+2015 CaT [Fe/H] (verified in ADS: 2015ApJ...798...77H, 1243 red giants) plus Geha velocities. The continuous
V/σ([Fe/H]) fit deconvolves [Fe/H] errors (≥ 0.4 dex); MW dSphs are a separate tidal test. Plan v3: `docs/plans/rotation_vs_feh.md`.
Next: fetch Ho+2015 Table 2 and Geha+2006 velocities.

### M31 dEs: data assembled (plan v3, step 1)
- Ho+2015 (J/ApJ/798/77 table2, 1243 rows) and Geha+2006 NGC 205 velocities (J/AJ/131/332 table3, 725 rows) fetched from the
  CfA mirror (`scripts/fetch_vizier.py`).
- `catalogues.ho2015()`: Ho [Fe/H] joined within 1" to Geha velocities, nearest wins. Matched 227/230 (NGC 147), 312/321
  (NGC 185), 221/224 (NGC 205; 3 Ho stars shared a Geha star). New test, written first and seen to fail, now passes (12 tests).
- Samples (`figures/de_data.png`):

  | Galaxy | N | median e_v [km/s] | median e_[Fe/H] [dex] | [Fe/H] spread [dex] | r_ell 90 % [R_half] |
  |---|---|---|---|---|---|
  | NGC 147 | 227 | 3.7 | 0.48 | 0.61 | 2.0 |
  | NGC 185 | 312 | 3.2 | 0.45 | 0.72 | 4.0 |
  | NGC 205 | 221 | 11.5 (single value for all Geha+2006 stars) | 0.54 | 0.76 | 9.4 |

- Seen in the figure:
  - NGC 147 rotates clearly (about ±15–20 km/s over ±2 R_half along PA 34°); metal-poor and metal-rich medians are similar by eye.
  - NGC 185 has a weak gradient.
  - NGC 205 is broad (±100 km/s) with both halves offset negative: possible M31 halo contamination (v_sys −246 km/s),
    11.5 km/s errors, tides. It needs a contaminant mixture, and its zero-point needs checking.
  - All masks lie along the major-axis strip, so the minor axis is barely sampled. Fix θ to the photometric PA, or put a
    prior on it.
- [Fe/H] errors (0.45–0.54 dex) are comparable to the [Fe/H] spread (0.6–0.76), confirming the dilution concern.
- Next: model with a mixture, [Fe/H]-error deconvolution and the Ho selection (stars without [Fe/H]); mocks before any real fit.

### Plan B for the dEs: model from the LMC kinfit, adversarial review, mocks, real fit (2026-11-06)
- **LMC model studied** (`lmc-old-population/src/lmcold/kinfit.py`): disc-plane forward model, v0 tanh(R/R0), mixture with
  broad outliers, populations per [Fe/H] bin, ML + NUTS, v0/σ_los per bin.
- **Adapted** to line-of-sight-only data: V_p h(R) x/R, σ(R), mixture. Draft `docs/plans/de_model.md`.
- **Review** (critic + proposer, Opus, 2 rounds; `docs/reviews/2026-11-06_de_model_debate.md`):
  - NGC 205 dropped.
  - Stars without [Fe/H] are not missing at random.
  - R0 is unconstrained, so quote V/σ inside the data.
  - The [Fe/H] error scale k is not identified (NGC 147 bright stars scatter most).
  - Significance comes from within-R permutation statistics.
  - Amplitude from mock-calibrated tercile contrasts.
  - **Pre-registered** (commit 63fd27e, before the real tests): primary T_rot in NGC 185 at α = 0.01; others at 0.0025;
    headline k = 1.
- **Code:** `src/olddisks/defit.py`, `tests/test_defit.py` (19 tests pass).
- **Mocks** (`scripts/de_mocks.py`, `results/de_mocks_summary.csv`, `figures/de_mocks.png`): 3200 mocks.
  - False-positive rate within α; null p-values flat.
  - Power at k = 1: NGC 185 0.07 / 0.21 / 0.72 at slopes 0.3 / 0.5 / 1.0 dex⁻¹; NGC 147 ≤ 0.07.
  - Dilution C = 0.59 (NGC 185), 0.42 (NGC 147) at k = 1; ≈ 1.0–1.1 at k = 0.5.
- **Real data** (`scripts/fit_de.py`, `results/de_fits.csv`, `figures/de_fit.png`):
  - [Fe/H]-blind fits:
    - NGC 185: V(2.5 R_half) = 8.0 km/s, σ = 23.3 km/s (Geha+2010: σ 24 ± 1).
    - NGC 147: V(1.5 R_half) = 13.9 km/s, σ = 17.0 km/s (Geha: σ 16 ± 1).
  - **NGC 185 (primary):**
    - T_rot = 0.016, p = 0.78: no detection.
    - V/σ by tercile 0.32 / 0.24 / 0.30.
    - d(V/σ)/d[Fe/H] = −0.04 ± 0.30 dex⁻¹ (k = 1), i.e. 95 % range −0.62 to +0.55; −0.02 ± 0.15 at k = 0.5.
    - Stable under the TRGB cut, PA ±10° and the R0 grid.
    - T_σ: p = 0.93.
  - **NGC 147:**
    - T_rot = 0.081, p = 0.22: no detection.
    - V/σ by tercile 0.59 / 0.84 / 1.00.
    - Slope 0.98 ± 0.58 (k = 1) and 0.41 ± 0.24 (k = 0.5): 1.7σ, not significant. PA +10° gives p = 0.06.
- **Caveats:**
  - The TRGB cut assumes Geha I is dereddened (UNVERIFIED).
  - sd taken from the null mocks.
  - The LMC comparison (matched radius, projection) is not yet done; no claim made.

### 11 MW dSphs (pre-registered extension), and fixes (2026-11-06)
**Bugs found and fixed (tests first):**
1. The inventory union did not deduplicate within a source. Walker+2023 rows for the same star from several instruments
   were counted 2–3×, so the earlier dSph counts were inflated.
2. The mock injected the slope with the wrong sign for V < 0.
3. Mock support could go negative (counter-rotation) for non-rotators, and the calibration then averaged to zero. Now
   V/σ(z) = max(0, |V|/σ + slope (z − z̄)).
4. Ji+2021 uses e_[Fe/H] = 99 placeholders; set to NaN.
5. pandas index misalignment of slope_sd.

The NGC 185/147 p-values are unchanged after the fixes. The amplitude calibration changed: NGC 185 −0.06 ± 0.45 and
NGC 147 1.01 ± 0.59 per dex (k = 1).

**Walker+2023 [Fe/H] errors validated on overlaps** (implied k = 0.88–1.02 against Tolstoy+2023, Pace+2020, Tolstoy+2025,
Spencer+2017). Typical members have e = 0.45–0.65 dex, so the [Fe/H] information fraction λ ≈ 0 for most Walker samples:
there is no [Fe/H] leverage. Catalogue offsets reach 0.26 dex (UMi), which justifies using a single source per galaxy.

**Results** (`scripts/fit_dsph.py`, `results/dsph_fits.csv`, `figures/dsph_fit.png`; α = 0.0045):
- **T_rot (velocity gradient versus [Fe/H]): no detection in any galaxy or catalogue.** Smallest p: Fornax (Walker) 0.041,
  Sculptor (Tolstoy) 0.088.
- **T_σ (normalised |residual| versus [Fe/H]): detected, metal-rich colder.**
  - Sculptor: p = 1e−4 in both Walker and Tolstoy; t = −0.17 and −0.19.
  - Fornax: p = 7e−4 (Walker) and 3e−3 (Pace).
  - Marginal: Sextans (Tolstoy) 0.006, Draco 0.037.
  - This reproduces the known chemo-dynamical populations (e.g. Battaglia+2008 Sculptor, 2008ApJ...681L..13B). It is
    a validation of the method, not new.
- **Calibrated slopes** where there is leverage:
  - Sculptor (Tolstoy) +0.82 ± 0.29 dex⁻¹ (2.8σ, tercile V/σ 0.02 → 0.30). The pre-registered T_rot gives p = 0.088, and
    17 runs were made, so this is a hint, not a detection.
  - Others are consistent with zero: UMi (Pace) −1.0 ± 0.7, Sextans (Tolstoy) 0.3 ± 0.36, Antlia II (Ji) 0.06 ± 0.23,
    Leo II (Spencer) −0.3 ± 0.4.
- **Gradients are small** (V/σ ≲ 0.3), except Antlia II (≈0.9, its known tidal gradient, Ji+2021) and UMi.
