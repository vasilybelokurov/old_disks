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
