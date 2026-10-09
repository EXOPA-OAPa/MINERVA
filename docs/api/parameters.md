# Parameter reference

[← API reference](README.md)

The configuration is a YAML file ([`config/config_example.yaml`](../../config/config_example.yaml),
[examples](../../examples/README.md)). Relative paths are relative to the configuration
file. Unknown keys in `run_options` and `analysis` are rejected.

* [`paths`](#paths) · [`run_options`](#run_options) · [`analysis`](#analysis) ·
  [`general_variables`](#general_variables) · [run mode](#run-mode) · [`catalog`](#catalog) ·
  [`single_target_params`](#single_target_params) · [catalog columns](#catalog-columns)

## paths

| Key | Default | Description |
|---|---|---|
| `data_dir` | `<repository>/data` | reference data: NSTS table, maximum ages, GAPS activity database, YaPSI tracks, default cadences |
| `user_dir` | folder of the configuration | catalog and cadence files; looked up before `data_dir` |
| `output_dir` | `Output/` next to the configuration | results, logs, `Plots/` |

## run_options

| Key | Default | Description |
|---|---|---|
| `n_workers` | 1 | catalog targets in parallel processes; integer or `auto` (all cores). Results identical to a sequential run |
| `rv_series_pdf` | `True` | RV time-series figure: `True` (every simulation), `False`, or N (first N simulations of each target) |
| `bootstrap` | `False` | resample stellar mass, T_eff, radius, age, P_rot, log R'HK within their uncertainties at every simulation |
| `export_debug_pdf` | `True` | write the debug figures to `Plots/DEBUG_<target>.pdf` |
| `debug_plot` | `True` | pre-whitening and reverse pre-whitening figures |
| `debug_fit_plot` | `True` | fit figures |
| `debug_pdf_all_iter` | `True` | GLS and fit figures for every simulation (`False`: first 3) |
| `show_on_screen` | `False` | show the figures interactively (blocks the run) |

The defaults reproduce the version used for the paper; for large runs disable the debug
figures and limit `rv_series_pdf` ([running faster](../../README.md#running-faster)).

## analysis

| Key | Default | Description |
|---|---|---|
| `n_start` | 10 | observations at the first step of the sequential recovery |
| `step` | 5 | observations added at each step |
| `sampling_strategy` | `sequential` | order in which observations are added: `sequential` (time order) or `random` |
| `fap_star` | 0.10 | false alarm probability for stellar (intruder) peaks |
| `fap_planet` | 0.50 | false alarm probability for the planet peak |
| `delta_bic_threshold` | 6.0 | minimum ΔBIC (planet model vs null hypothesis) for a detection |
| `random_seed` | 42 | seed of the random generator; the i-th catalog target uses `random_seed + i` |
| `gls_method` | `recurrence` | GLS implementation: `recurrence` (fastest), `exact`, or `reference` (original code) |

## general_variables

| Key | Unit | Default (example) | Description |
|---|---|---|---|
| `name_sampling` | | `harps.rv` | default cadence for P ≥ 2 d (HD 164922, Benatti et al. 2020) |
| `name_sampling_usp` | | `time_usp.txt` | default cadence for P < 2 d (TOI-1807, Nardiello et al. 2022) |
| `inputfile_nsts` | | `NSTS_noise_new.txt` | photon noise vs V and spectral type |
| `inputfile_mass_age` | | `max_age_from_mass.dat` | maximum age vs mass |
| `inputfile_activity` | | `PLOT_POINTS_DATA_...csv` | GAPS activity database |
| `sig_floor` | m/s | 0.5 | instrumental noise σ_instr (HARPS-N: Cosentino et al. 2014) |
| `sig_pl` | m/s | 2 | unresolved additional planets σ_pl (injected noise only, not in the error bars) |
| `sig_best` | | `[5, 5.5]` | SNR thresholds K/σ_K (5 ↔ 20% mass precision, 3.33 ↔ 30%) |
| `n_iterations` | | 100 | Monte Carlo simulations per target |

## Run mode

| Key | Description |
|---|---|
| `read_catalog` | `True`: run every target of `catalog.namecatalog` |
| `single_target` | `True`: run the target of `single_target_params` |

Exactly one of the two must be `True`.

## catalog

| Key | Description |
|---|---|
| `namecatalog` | catalog file (in `user_dir` or `data_dir`) |
| `real_cadence` | `True`: each target uses the cadence file of its `Cadence FileName` column (default cadences when empty) |

## single_target_params

| Key | Unit | Required | If empty |
|---|---|---|---|
| `targetname` | | yes | |
| `real_cadence` | | | `False` (default cadences) |
| `star_mass`, `err_star_mass` | M☉ | mass yes | error: 10% |
| `star_radius`, `err_star_radius` | R☉ | | error: 10% |
| `star_temp`, `err_star_temp` | K | T_eff yes (≤ 7500 K) | error: 100 K |
| `sp` | | | from T_eff |
| `v_mag`, `gaia_mag` | mag | | V, otherwise Gaia, otherwise 10 |
| `pl_mass_Mj` or `pl_mass_Me` | M_J, M⊕ | one of the two | |
| `pl_radius_Rj`, `pl_radius_Re` | R_J, R⊕ | | |
| `pl_period`, `err_pl_period` | d | period yes | (`err_pl_period` not used) |
| `transit_midpoint` | JD − 2450000 | yes | |
| `eccentricity`, `periastron` | –, deg | | 0, 0 |
| `age`, `err_age` | Gyr | | drawn up to the end of the main sequence; error 10% |
| `vsini_tag` | | | rotation regime drawn at random |
| `Prot`, `err_Prot` | d | | from age and regime; error 10% |
| `log_Rhk`, `err_log_Rhk` | dex | | from the Rossby number; error 0.05 dex |
| `K_act` | m/s | | from log R'HK and T_eff |

Uncertainties are used only with `run_options.bootstrap`.

## Catalog columns

CSV with `,` or `;` as separator; empty cells allowed for optional columns
([`examples/01_catalog/catalog_example.csv`](../../examples/01_catalog/catalog_example.csv)).

| Column | Unit | Equivalent key |
|---|---|---|
| `Planet Name` | | `targetname` |
| `Star Mass [Ms]`, `Star Mass Error [Ms]` | M☉ | `star_mass`, `err_star_mass` |
| `Star Radius [Rs]`, `Star Radius Error [Rs]` | R☉ | `star_radius`, `err_star_radius` |
| `Star Temperature [K]`, `Star Temperature Error [K]` | K | `star_temp`, `err_star_temp` |
| `Star Spectral Type` | | `sp` |
| `Star V Mag`, `Star Gaia Mag` | mag | `v_mag`, `gaia_mag` |
| `Planet Mass [Mj]`, `Planet Mass [Me]` | M_J, M⊕ | `pl_mass_Mj`, `pl_mass_Me` |
| `Planet Radius [Rj]`, `Planet Radius [Re]` | R_J, R⊕ | `pl_radius_Rj`, `pl_radius_Re` |
| `Planet Period [days]`, `Planet Period Error [days]` | d | `pl_period`, `err_pl_period` |
| `Transit Mid Time [JD - 2450000]` | d | `transit_midpoint` |
| `Star Age [Gyr]`, `Star Age Error [Gyr]` | Gyr | `age`, `err_age` |
| `Stellar Rotation Period [d]`, `Star Rotation Period Error [days]` | d | `Prot`, `err_Prot` |
| `logRhk`, `logRhk Error` | dex | `log_Rhk`, `err_log_Rhk` |
| `Stellar Activity RV Amplitude [m/s]` | m/s | `K_act` |
| `Eccentricity`, `Periastron` | –, deg | `eccentricity`, `periastron` |
| `Cadence FileName` (or `Catalog FileName`) | | cadence file of the target |

The rotation regime (`vsini_tag`) cannot be set from a catalog: it is drawn at random.
