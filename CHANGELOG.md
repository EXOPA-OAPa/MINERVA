# Changelog

## [Unreleased]

Restructuring of the code used for the paper (`v1.0-paper`) into a modular,
documented package, and fixes of several issues of that version.
**The fixes marked with ⚠️ change the results** with respect to `v1.0-paper`.

### Fixed

- ⚠️ **Eccentricity.** Eccentricity and argument of periastron are now used when
  given (catalog or configuration). In `v1.0-paper` they were always set to zero
  (typo `eccetricity` in the catalog loader, capitalised keys in the single-target
  loader). The time of periastron is computed from the transit time with the exact
  conversion of RadVel; the argument of periastron is in degrees.
- ⚠️ **ΔBIC sign.** ΔBIC = BIC(null) − BIC(full model) is no longer taken in absolute
  value: a planet model worse than the null hypothesis is not counted as evidence.
- ⚠️ **Bootstrap.** Every simulation resamples the stellar parameters from their
  nominal values. In catalog mode, `v1.0-paper` resampled from the previously
  resampled values (random walk). In single-target mode the uncertainty on the
  stellar mass was set equal to the mass.
- ⚠️ **V magnitude fallback.** V magnitude, otherwise Gaia magnitude, otherwise 10, as
  documented. In `v1.0-paper` a valid V could be replaced by 10, and the fallback was
  missing in catalog mode.
- ⚠️ **Stars without a convective envelope** (τ ≈ 0) get no activity signal and no
  activity term in the error bars, and are fitted with a global model. In `v1.0-paper`
  they were flagged with log R'HK = 0, which the code read as extreme activity
  (a ~1 m/s signal was injected and the seasonal fit was used). The check on the
  recovered stellar period is not applied to them.
- ⚠️ **K-dwarf activity relation.** The exponential GAPS/ArMS calibration is not
  extrapolated beyond its range: log R'HK is capped at −3.99 (beyond it, the relation
  reached unphysical amplitudes of thousands of m/s).
- ⚠️ **Stellar age.** Random ages are drawn only up to the end of the main sequence of
  the YaPSI track used for τ (central hydrogen exhausted). In `v1.0-paper` the upper
  limit included the subgiant phase, with τ of hundreds of days for stars with a
  main-sequence T_eff.
- ⚠️ **Failed simulations** count as N_total (all available observations) in the
  "ALL" statistics, instead of 0. When at least half of the simulations fail, the
  median is reported as a lower limit.
- ⚠️ **Pre-whitening.** The loop over the intruder peaks stops at the first accepted
  intruder and otherwise tries the next one. In `v1.0-paper` the `break` was inside
  the debug-figure block, so the results depended on the `DEBUG_PLOT` flag.
- ⚠️ **Reproducibility.** The activity term of the error bars uses the seeded random
  generator of the run: two runs with the same configuration give the same results.
- ⚠️ **Differential rotation.** The relative shear α is drawn in [0.01, 0.08] for
  highly active stars (log R'HK ≥ −4.5) and in [0.10, 0.20] otherwise, as described in
  the paper (Sect. 2.4). `v1.0-paper` used the opposite assignment.
- K_act at T_eff = 3500 K exactly (it was undefined), missing T_eff now reported as an
  error, catalogs with missing optional columns, error messages in the results file,
  per-threshold failures in the summary file, summary header printed once.

### Changed

- Modular package (`src/minerva/`, src layout), one function per relation or analysis step; English
  comments, messages and documentation.
- Configuration: new sections `paths`, `run_options` (previously flags in `main.py`)
  and `analysis` (`n_start`, `step`, `sampling_strategy`, `fap_star`, `fap_planet`,
  `delta_bic_threshold`, `random_seed`, previously fixed in the code; defaults are
  the values of `v1.0-paper`).
- Result and log files are overwritten at every run (they were appended to).
- The GLS of Zechmeister is shipped unchanged in `minerva/_vendor/`; the YaPSI
  tracks are stored gzip-compressed (only those used by the code).
- Debug figures are built only when they are saved or shown (about 2x faster runs
  with the debug PDF enabled).
- Faster runs, same results: YaPSI tracks read once and kept in memory; the RV-series
  figure can be limited or disabled (`rv_series_pdf`); catalog targets can run in parallel
  (`n_workers`); vectorised GLS periodogram (`gls_method`, default `recurrence`, with the
  original implementation available as `reference`).
- Unused code removed (functions never called, disabled code, unused imports;
  `george` is no longer a dependency).
- RadVel is no longer a dependency: the two functions MINERVA used (Keplerian RV and time of
  periastron) are in `minerva/keplerian.py`, with credit. This removes h5py, pytables,
  ultranest and jupyter from the installation, which failed to build on some systems.
  Circular orbits give identical results; eccentric orbits differ by ~1e-11 K (pure NumPy
  Kepler solver instead of RadVel's C solver), with no effect on the results in our tests.

### Added

- Documentation in `docs/`: installation, API reference (architecture, one page per
  module, parameter reference, conventions, command-line tools, tests) and science cases.
- Examples: one folder per science case, each with a complete configuration
  (`examples/NN_<case>/config.yaml`).
- `scripts/run_minerva.py --quick` (5 simulations per target) and `--workers N`.
- Validation checks `tests/check_parallel.py` (parallel = sequential, reproducibility) and
  `tests/check_gls_methods.py` (vectorised GLS vs reference).

- Regression test (`tests/regression/`) and a syntax-tree comparison tool for
  documentation-only changes (`scripts/check_same_code.py`).
- Example catalog with fictitious targets and a synthetic cadence (`examples/`).
