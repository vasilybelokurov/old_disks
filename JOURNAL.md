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
