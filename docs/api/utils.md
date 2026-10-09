# Utilities — `utils`, `tracks`, `plotting`, `logs`

[← API reference](README.md)

* [`utils`: input files and configuration](#utils)
* [`tracks`: YaPSI tracks](#tracks)
* [`plotting`: debug figures](#plotting)
* [`logs`](#logs)

---

## utils

| Function | Description |
|---|---|
| `load_config(filename)` | reads the YAML configuration: a flat dictionary with `general_variables`, the active mode and its parameters (`catalog` or `single_target_params`). Exactly one of `read_catalog`, `single_target` must be `True`; a single target needs `pl_mass_Mj` or `pl_mass_Me` |
| `import_file_catalog(namecatalog, InDir=None)` | catalog as a pandas DataFrame; separator `,` or `;` |
| `load_target_params_from_catalog(catalog, i)` | parameters of the i-th target ([catalog columns](parameters.md#catalog-columns)); empty cells become `None`, missing optional columns are allowed; planetary mass and radius converted between Jupiter and Earth units; default uncertainties |
| `load_single_target_params(config)` | the same for `single_target_params` |
| `resolve_magnitude(target_params)` | V magnitude: V, otherwise Gaia, otherwise 10 |
| `import_file_sampling(...)` | epochs of the standard and USP cadences (first column of each file) |
| `import_nsts_noise(...)` | photon-noise table: magnitudes and σ_ph for F, G, K, M stars |
| `import_max_age_from_mass(...)` | maximum stellar age vs mass |
| `import_file_activity(...)` | GAPS activity database (`SP_CLASS`, `logRHK`, `Jitter_ms`, ...) |
| `find_seasons_automatic(time_series, min_points=7, n_cuts=10, min_gap_days=50)` | cuts between observing seasons: the largest gaps longer than `min_gap_days` that leave at least `min_points` per season |
| `define_seasons(time_sampling, xm)` | indices of the observations of each season |
| `print_params(target_params)` | parameters of the target, in the log |

---

## keplerian

From RadVel (Fulton et al. 2018, MIT License), the only functions MINERVA needs:

| Function | Description |
|---|---|
| `rv_drive(t, orbel)` | radial velocity of a Keplerian orbit, `orbel = [per, tp, e, ω (rad), K]` |
| `timetrans_to_timeperi(tc, per, ecc, omega)` | time of periastron from the time of transit |
| `true_anomaly(t, tp, per, e)`, `kepler(M, e)` | true anomaly; Kepler's equation (iterative solver, 1e-12) |

---

## tracks

YaPSI evolutionary tracks (Spada et al. 2017), gzip-compressed in `data/YaPSI_Tracks/`.
For a star of mass M the track with the largest mass ≤ M is used (mixing length
`A1p91804` below 1 M☉, `A1p82126` above). Files are read once and kept in memory.

| Function | Description |
|---|---|
| `track_pro(starmass, age, InDir)` | convective turnover time τ [d] at the model closest in age |
| `ms_lifetime(starmass, InDir, xc_min=1e-3)` | age [yr] at the end of the main sequence: first model with central hydrogen X_c < `xc_min` |
| `fap_levels(fap, time, fbeg, fend)` | GLS power of a false alarm probability (analytic approximation; used by the RV-series figure) |

---

## plotting

Debug figures, enabled in `run_options`; they never change the results. Each function
takes the variables it needs as keyword arguments; `pick(locals(), func)` selects them
from the caller.

| Function | Figure | Option | Saved for |
|---|---|---|---|
| `plot_rv_series` | RV time series, its components and GLS | `rv_series_pdf` | every simulation, or the first N |
| `gls_peaks` | GLS of the current step, FAP thresholds, classified peaks | `export_debug_pdf` | every simulation (`debug_pdf_all_iter`) or the first 3 |
| `prewhitening_attempt` | one pre-whitening attempt: raw vs residual GLS | `debug_plot` | first 3 simulations |
| `reverse_prewhitening` | GLS of the planet-subtracted data | `debug_plot` | first 3 simulations |
| `fit_diagnostics` | data, model, residuals | `debug_fit_plot` | every simulation (`debug_pdf_all_iter`) or the first 3 |
| `phase_folded_fit` | planet and seasonal activity components, phase-folded | `debug_fit_plot` | first 3 simulations |
| `statistics_history` | K, ΔBIC, χ²_red and jitter vs number of observations | `export_debug_pdf` | first 3 simulations |

Debug figures go to `Plots/DEBUG_<target>.pdf`, RV series to `Plots/plot_series_*.pdf`.

---

## logs

`ToggleableLogger(filename, to_file=True, to_screen=False)`: replaces `sys.stdout` during a
run, so that `print` writes to the log file (`Log_*.txt`); progress bars and summaries go
to the terminal. `set_mode(to_file, to_screen)`, `restore_stdout()`.
