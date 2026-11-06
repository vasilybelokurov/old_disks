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
