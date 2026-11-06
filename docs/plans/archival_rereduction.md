# Plan: rebuild per-star catalogues of published programmes for isolated dIrrs (2026-11-06, revised)

**User decision (2026-11-06): only data from published programmes.** Cole 71.B-0516(B), Held 075.D-0838 (NGC 6822 part)
and Mayer 072.B-0307 have no refereed publication and are excluded. The plan is now to recover or rebuild per-star
catalogues for programmes that are published but whose tables we do not have, and to re-reduce published data only
where that is needed for homogeneity (velocity zero-points, CaT calibration).

## Summary

- **Goal.** Maximise the number of red giants with v_los and [Fe/H] in IC 1613, NGC 6822 and WLM from published
  programmes. `results/feasibility.csv` puts the threshold for an 80 % detection of ΔV/σ = 0.5 between [Fe/H] halves at
  about 550 stars (270 per half) reaching ≥ 2 R_half.
- **Missing published catalogues** (absent on the VizieR mirror while CDS is down):
  - Belland+2020, NGC 6822, DEIMOS ([2020ApJ...903...10B](https://ui.adsabs.harvard.edu/abs/2020ApJ...903...10B)): N unknown.
  - Taibi+2024, IC 1613, MUSE 097.B-0373 ([2024A&A...689A..88T](https://ui.adsabs.harvard.edu/abs/2024A%26A...689A..88T)): about 300 red giants per the literature census.
  - Leaman+2012, WLM, FORS2 072.B-0497 + DEIMOS C257D ([2012ApJ...750...33L](https://ui.adsabs.harvard.edu/abs/2012ApJ...750...33L)): 180 red giants.
- **Route A (preferred):** obtain these tables from CDS when it returns, from the journals' machine-readable tables, or
  from the authors.
- **Route B (fallback):** re-reduce the raw data of these published programmes from ESO and KOA.
- **Expected totals (unverified):**
  - IC 1613: 139 + about 300 MUSE (overlap unknown).
  - NGC 6822: 380 + Belland (unknown).
  - WLM: 180.
  - Only NGC 6822 is plausibly near the threshold.
- **Validation before science.** Any re-reduction must reproduce the published values for the same stars.

---
## Original version (before the published-only decision), kept for the record

## Summary

- **Goal.** Raise IC 1613 and NGC 6822 to ≳ 550 red giants with line-of-sight velocity and Ca II triplet (CaT) [Fe/H],
  reaching ≥ 2 R_half. `results/feasibility.csv` puts the threshold for an 80 % detection of ΔV/σ = 0.5 between [Fe/H]
  halves at about 270 stars per half.
- **The audit corrects the census.** Most ESO hours on these galaxies are not red-giant spectroscopy:
  - NGC 6822: VIMOS 118 h is imaging, 25 h is IFU spectroscopy of HII regions.
  - IC 1613: the VIMOS and GIRAFFE spectra are of massive stars.
  - Sextans A: all of the UVES time is supergiants or its globular cluster.
  - **Sextans A is dropped.**
- **The real targets are unpublished FORS2 and GIRAFFE CaT programmes on IC 1613 and NGC 6822**, plus Keck/DEIMOS
  masks that are not in published catalogues.
- **One pipeline for all multi-slit data.** PypeIt (supports VLT/FORS2 MXU and Keck/DEIMOS). GIRAFFE uses the ESO
  `giraf` pipeline.
- **Validation before science.** Re-reduce masks that are already published and reproduce the published velocities and
  [Fe/H] before touching the new ones.

## 1. Data identified (ESO archive, `dbo.raw`, queried 2026-11-06)

| Programme | PI | Instrument | Galaxy: fields (hours) | Published? | Role |
|---|---|---|---|---|---|
| 71.B-0516(B) | Cole | FORS2 MXU, CaT | IC 1613: bar, bar-0, disk (17.4 h); NGC 6822: Halo, S2 (19.2 h); Phoenix (6.7 h); 8 GCs (NGC 104, 288, 362, 2298, 6528, 6681, 6715, 7078) | No refereed paper found (AAS abstracts [2004AAS...205.9311C](https://ui.adsabs.harvard.edu/abs/2004AAS...205.9311C), [2004oee..sympE..62Z](https://ui.adsabs.harvard.edu/abs/2004oee..sympE..62Z)) | **New stars**; GCs = CaT calibrators |
| 075.D-0838(A) | Held | FORS2 MXU | NGC 6822 (6 OBs, 18.1 h); Phoenix (28.3 h); 47 Tuc | NGC 6822 part: none found in ADS. Phoenix probably in Kacharov+2017 (to verify) | **New stars** (NGC 6822) |
| 072.B-0307(A) | Mayer | GIRAFFE MEDUSA | IC 1613_RGB (11 OBs, 8.7 h), IC 1613_AGB (5.4 h) | None found | **New stars** (IC 1613); grating to check (must cover CaT) |
| 69.B-0305(B) | Tolstoy | FORS2 MXU | NGC 6822 f1, f2 (6.9 h); Tucana (10.4 h); GCs | Swan+2016 ([2016MNRAS.456.4315S](https://ui.adsabs.harvard.edu/abs/2016MNRAS.456.4315S), 72 stars, 2 fields) | **Validation** |
| 072.B-0497(B) | Cole | FORS2 MXU | WLM bar fields (26.7 h); GCs | Leaman+2009/2012 (to confirm it is the same data) | Validation and/or WLM rebuild (Leaman catalogue missing on the VizieR mirror) |
| 097.B-0373(A) | Battaglia | MUSE | IC 1613 (3 fields, 9.2 h) | Taibi+2024 ([2024A&A...689A..88T](https://ui.adsabs.harvard.edu/abs/2024A%26A...689A..88T)) | Use the published catalogue; re-extract only if it is unavailable |

Keck/DEIMOS (KOA, `data/archive_census.csv`): NGC 6822 35 h / 11 programmes, IC 1613 19 h / 6, Leo A 29 h / 6,
Pegasus 21 h / 5, Aquarius 25 h / 4, Sag DIG 26 h / 6, plus the M31 satellites. Kirby+2014/2017/2020 and Belland+2020
used part of this. Masks absent from published catalogues are a second tranche (phase 4).

**Star yields are not yet known.** A FORS2 MXU mask holds about 30–60 slits, and GIRAFFE MEDUSA has 132 fibres. Phase 0
counts slits and fibres from the headers. Rough expectations, to be verified:
- IC 1613: +100–150 (Cole) and +100–300 (Mayer), on top of 139 published + about 300 MUSE red giants.
- NGC 6822: +100 (Cole) and +100–200 (Held), on top of 380 published + Belland+2020 (count unknown).

## 2. Phases, tests and acceptance criteria

### Phase 0: audit (≈ 1 day)
- Download headers only. Record per OB: grism/grating (FORS2 1028z or 1400V covering 8450–8700 Å; GIRAFFE LR8/HR21),
  number of slits or fibres, exposure, seeing, and whether raw calibrations exist.
- Cross-match slit/fibre positions with published members (Kirby+2014, Swan+2016, Taibi+2024) to count new stars.
- **Output:** `data/rereduction_targets.csv` (programme, OB, grating, n_slits, n_new_stars).
- **Go / no-go:** proceed for a galaxy only if published + new stars ≥ 550, or ≥ 400 with coverage to 2 R_half.

### Phase 1: validation reduction (≈ 1 week)
- Reduce 69.B-0305(B) (the NGC 6822 f1/f2 masks behind Swan+2016) with PypeIt.
- Velocities:
  - Cross-correlate the CaT region against synthetic templates (or fit with `rvspecfit`).
  - Correct slit mis-centring with the telluric A band or the ~7600 Å O₂ band.
  - Combine exposures.
- [Fe/H]: CaT equivalent widths with a single published calibration (Starkenburg+2010 or Carrera+2013, to be chosen
  and cited after checking), using V − V_HB or absolute magnitude.
- **Tests:**
  1. Recover Swan+2016 velocities: median offset < 1 km/s; scatter of differences consistent with the quoted errors
     (χ²/dof < 2).
  2. Recover Swan+2016 [Fe/H]: median offset < 0.1 dex, scatter < 0.2 dex.
  3. Globular clusters in the same programme reproduce the Harris catalogue [Fe/H] within 0.1 dex and v_sys within 2 km/s.
- **Figure:** `figures/validation_swan2016.png` (Δv and Δ[Fe/H] against S/N).

### Phase 2: new FORS2 masks (≈ 1–2 weeks)
- Reduce Cole 71.B-0516(B) (IC 1613, NGC 6822, and its eight GCs as the CaT calibration set) and Held 075.D-0838
  (NGC 6822) with the validated pipeline.
- **Tests:**
  1. Repeat stars across masks agree within errors.
  2. Stars in common with Kirby+2014 or Taibi+2024 give the velocity zero-point offset and the [Fe/H] offset (on the
     inventory scale). Offsets are applied as measured; none are assumed.
- **Photometry for the CaT calibration:** archival pre-images, or WSDB DELVE/DECaLS for V and I.

### Phase 3: GIRAFFE IC 1613 (≈ 1 week)
- Reduce 072.B-0307(A) with the ESO `giraf` pipeline. Sky-subtract using the sky fibres.
- Same velocity and CaT procedure if the grating covers CaT. If it does not (e.g. a blue HR setting), velocities only.
- **Test:** overlap with the Phase 2 IC 1613 masks and with Taibi+2024.

### Phase 4: KOA DEIMOS second tranche (optional, ≈ 2 weeks)
- Identify DEIMOS masks in NGC 6822, IC 1613, Leo A, Pegasus, Aquarius, Sag DIG and the M31 satellites whose slits are
  not in published catalogues.
- Reduce with PypeIt (DEIMOS is its best-tested instrument).
- **Test:** reproduce Kirby+2014 velocities for one already-published mask (offset < 1 km/s).

### Phase 5: merge and decide
- Add the new catalogues to `src/olddisks/catalogues.py` and `scripts/build_inventory.py`.
- Put [Fe/H] and velocities on one scale using overlap stars.
- Rerun `scripts/feasibility_mocks.py`.
- **Decision:** fit the real data only for galaxies where the mock power is ≥ 0.8. Otherwise report upper limits or
  stack galaxies.

## 3. Risks

- **Old FORS2 calibrations (2002–2005):** arcs and flats may be missing or taken with a different mask. Phase 0 checks this.
- **Slit mis-centring:** produces velocity errors of several km/s, comparable to σ ≈ 10 km/s in IC 1613. Needs the
  telluric correction.
- **Crowding and young contamination:** red supergiants and AGB stars mimic red giants in these star-forming galaxies.
  Select on the CMD; keep only stars below the RGB tip with red-giant colours.
- **CaT calibration at [Fe/H] < −2 and for young stars:** the CaT–[Fe/H] relation assumes old red giants. Age effects
  reach a few tenths of a dex. This matters because the [Fe/H] split is our proxy for age.
- **Ethics:** the data are public, but the Cole and Held programmes were never published. Contact the PIs (A. Cole,
  M. Bellazzini/E. Held) before publishing results from them, and offer co-authorship if they have reduced catalogues.
  They may also simply share them, which would remove Phases 2–3.

## 4. Commands (to be filled as each phase is implemented)

```bash
source .venv/bin/activate
pip install pypeit                          # FORS2 MXU and DEIMOS
python scripts/eso_audit.py                 # Phase 0 -> data/rereduction_targets.csv
python -m pytest tests/test_validation.py   # Phase 1 acceptance tests
```
