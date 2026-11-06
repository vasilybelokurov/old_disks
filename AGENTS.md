# old_disks: project rules

**Question.** How does the rotational support of the stars in dwarf galaxies depend on metallicity (and, where it can be
inferred, age)? Do the oldest, most metal-poor stars rotate, and how does $v/\sigma$ change from metal-poor to metal-rich
populations within one galaxy and across galaxies of different mass and environment?

Related projects (same author, same methods; read them before writing new code):
- `../rrl_lmc_rotation`: LMC RR Lyrae and giants, forward-modelled kinematics, mocks, experiment registry.
- `../lmc-old-population`: release-quality version of the LMC work; data provenance model (`data/sources.yaml`).

## Code rules (user, 2026-11-06)

1. **Write the least code that works.**
   - No speculative options, flags, modes or abstractions; add one only when a second real use exists.
   - No unused functions, no commented-out code, no copies of superseded variants.
   - Shared code lives in `src/olddisks`; scripts import it rather than copy it.
2. **Always reproduce a bug with a failing test first.** Write the test, run it and see it fail, fix the code, see it
   pass, and commit both together.
3. **Show a plot before claiming a result.** A number becomes a result only after a figure showing it (data, model and
   mocks where relevant) has been made and looked at. Cite the figure path next to the number.

## Part 1: spectroscopic inventory

Deliverable: `data/inventory.csv`, one row per (galaxy, dataset), plus `data/inventory.md` describing every column.

Minimum columns:
- galaxy, dataset/paper, instrument, programme ID, public status (catalogue / spectra / proprietary);
- number of members with velocities, number with [Fe/H] (or CaT EW), number with repeat epochs;
- median velocity error and median [Fe/H] error;
- radial coverage in units of $R_\mathrm{half}$, coverage along major and minor axes;
- overlap with photometry deep enough for ages (oldest main-sequence turnoff reached: yes/no);
- URL of the table (ADS/VizieR/archive) and its sha256 once downloaded.

Scope, in priority order:
1. Resolved-star spectroscopy of Local Group dwarfs with per-star velocity **and** metallicity: isolated dIrrs/transition
   dwarfs (e.g. WLM, NGC 6822, IC 1613, Leo A, Pegasus, Aquarius, VV 124, Cetus) and the bright dSphs (e.g. Sculptor,
   Fornax, Carina, Sextans, Draco, Ursa Minor, Leo I/II), plus M31 satellites where samples are large enough.
2. Large surveys covering dwarfs: APOGEE DR17, SDSS-V, DESI MWS, Gaia DR3 RVS, S5, H3 (check WSDB first, `wsdb` skill).
3. Integrated-light IFU data (Fornax3D, GECKOS, MaNGA, SAMI) are **deferred**: they do not give per-star metallicity
   splits. Add them only if the user asks.

The kick-off brief (2026-11-06) quoted sample sizes and results (e.g. WLM, 180 RGB stars, $v_\mathrm{rot}\approx17$ km/s,
arXiv:1202.4474; seven-galaxy DEIMOS sample, 862 stars). Treat every such number as **unverified** until checked against
the paper (ADS/arXiv URL recorded in the inventory). Do not invent tables, programme IDs or counts; mark unknowns as
unknown.

## Part 2: rotational support versus metallicity

Per galaxy, fit members, foreground and kinematics jointly:
- likelihood per star with its own velocity error, repeat epochs combined, binary/jitter term;
- foreground/Milky Way contamination as a mixture component;
- mean velocity field with a rotation amplitude and position angle (or a rotation curve) and a dispersion, each allowed to
  depend on [Fe/H] (metallicity bins first, a continuous model only if bins show a trend);
- correct the line-of-sight velocities for perspective rotation induced by the galaxy's systemic proper motion
  (Gaia) before interpreting a velocity gradient as rotation;
- report $v_\mathrm{rot}/\sigma$ (state the definition and radius used) with full posteriors; keep upper limits;
  compare rotating against non-rotating models.

Before interpreting any real-data fit, run mocks with the observed positions, errors and selection: recover a known
rotation and dispersion per metallicity bin, and show the null (no rotation) case does not produce a detection.

Metallicity is not age. Any statement about ancient (>10 Gyr) stars needs an age argument (CMD + star-formation history);
state it explicitly or restrict the claim to metallicity.

## Data and constants

- `data/sources.yaml` is the single source of truth for inputs: origin and release, URL or query, fetch script, row count,
  sha256, paper reference.
- Galaxy parameters (centre, distance, position angle, ellipticity, $R_\mathrm{half}$, systemic velocity and proper
  motion) come only from the pinned Local Volume Database release listed in `sources.yaml` (path in
  `src/olddisks/paths.py`). Any override must be recorded in `sources.yaml` with its reference URL; no file hard-codes them.
- Large downloads go in `~/data/old_disks`, not in the repository.

## Bookkeeping

- Every step gets a dated entry in `JOURNAL.md`: findings, decisions, failed approaches, limitations, next step. Update it
  without being asked, during long work as well as at the end.
- Figures: `figures/*.png`, made by `scripts/fig_*.py`.
- Python from the project `.venv` (set xattr `com.dropbox.ignored` on it); tests with `pytest`.
- Commit after each tested step.
- Second opinions from Codex (ask-codex skill); save reviews in `docs/reviews/`.
- Citations: real, checked, with ADS/arXiv/DOI URLs; otherwise label UNVERIFIED.
