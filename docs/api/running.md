# Running MINERVA — `runner`, `context`, `simulation`, `summary`

[← API reference](README.md)

* [`runner.run`](#run)
* [`context.RunContext`](#runcontext)
* [`runner.run_catalog`, `runner.run_single_target`](#catalog-and-single-target-runs)
* [`simulation.run_simulation`](#run_simulation)
* [`simulation.setup_seasons_logic`, `get_indices_at_step`](#cadence-and-steps)
* [`summary.print_execution_summary`](#statistics-over-the-simulations)

---

## run

`minerva.runner.run(config_path, overrides=None)`

Runs MINERVA with a YAML configuration: the catalog (`read_catalog: True`) or the single
target (`single_target: True`). Returns the `RunContext` of the run.

| Argument | Description |
|---|---|
| `config_path` | path of the YAML configuration ([parameter reference](parameters.md)) |
| `overrides` | optional `{"n_iterations": int, "n_workers": int or "auto"}`, replacing the configuration (used by `--quick` and `--workers` of the command line) |

```python
from minerva.runner import run
run("examples/02_single_target/config.yaml", overrides={"n_iterations": 5})
```

---

## RunContext

`minerva.context.RunContext(config_path, overrides=None)`

Everything shared by the simulations of a run: configuration, folders, options,
analysis parameters and reference data. Passed to most functions as `var`.

Main attributes:

| Attribute | Content |
|---|---|
| `config` | flat dictionary of the active configuration (`utils.load_config`) |
| `options` | `run_options`, with defaults ([parameters](parameters.md#run_options)) |
| `analysis` | `analysis` parameters, with defaults ([parameters](parameters.md#analysis)) |
| `Gls` | GLS class chosen by `analysis.gls_method` ([`fast_gls`](recovery.md#fast_gls)) |
| `data_dir`, `user_dir`, `InDir`, `OutDir`, `FigDir` | folders |
| `time_sampling`, `time_sampling_usp` | epochs of the cadences in use (replaced per target with `real_cadence`) |
| `tmp_mag`, `errRV_F`, `errRV_G`, `errRV_K`, `errRV_M` | NSTS photon-noise table |
| `mass_tmp`, `max_age` | maximum stellar age vs mass |
| `df_activity` | GAPS activity database (error-bar term) |
| `seasons` | indices of the observing seasons of the current cadence |
| `numero_iterazioni` | number of simulations per target (`n_iterations`) |
| `sig_floor`, `sig_pl`, `sig_best` | noise terms and SNR thresholds |

**`resolve(filename)`**: full path of an input file. Absolute paths are kept; relative
names are looked up in `user_dir` first, then in `data_dir`.

Unknown keys in `run_options` or `analysis` raise a `ValueError`.

---

## Catalog and single-target runs

**`run_catalog(var)`**: every target of the catalog. For target `i` the random generator
is `numpy.random.default_rng(random_seed + i)`. With `run_options.n_workers > 1` the
targets run in parallel processes (`concurrent.futures`, spawn context): each process
builds its own `RunContext`, and each target writes its own log (`Log_<target>.txt`) and
RV-series PDF; the catalog summary is written at the end in catalog order. Results are
identical to a sequential run.

**`run_single_target(var)`**: the target of `single_target_params`, with generator
`default_rng(random_seed)`.

Both write `results_<target>.txt` (one line per simulation), the statistics of
`print_execution_summary`, the log and, if enabled, the debug and RV-series PDFs.

---

## run_simulation

`minerva.simulation.run_simulation(var, target_params, niter, pdf, debug_pdf, rng)`

One Monte Carlo simulation of a target: synthetic star, RV time series, sequential
recovery (the [flowchart](README.md#architecture)).

| Argument | Description |
|---|---|
| `var` | the `RunContext` |
| `target_params` | parameters of the target (copied: the original is never modified) |
| `niter` | index of the simulation (the first 3 save the debug figures) |
| `pdf`, `debug_pdf` | open `PdfPages` for the RV-series and debug figures, or `None` |
| `rng` | the random generator of the target |

Returns `(log_Rhk, nobs, status)`, with `nobs[th]` the number of observations needed for
the SNR threshold `th` (the total number of observations if never reached) and
`status[th]` ∈ {1 detection, 3 false positive, 0 no detection}; `None` if the star
cannot be simulated (T_eff > 7500 K, missing T_eff, invalid spectral type).

---

## Cadence and steps

**`setup_seasons_logic(var, target_params)`**: loads the cadence of the target (its own
file with `real_cadence`, otherwise the defaults) and splits it into observing seasons
(`var.seasons`).

**`get_indices_at_step(j, n_start, step, N_total, strategy='sequential', seed=42)`**:
indices of the observations used at step `j`: the first `n_start + j·step` epochs in
time order (`sequential`) or a fixed random subset (`random`; nested from one step to
the next). Returned sorted.

---

## Statistics over the simulations

`minerva.summary.print_execution_summary(target_name, log_rhk, nobs_history_dict, status_flag_dict, output_filename, thresholds, is_bootstrapped=False)`

For each SNR threshold: median and 16th/84th percentiles of the number of observations
over all simulations (failures count as the total number of observations; a median
marked `>=` is a lower limit, when at least half of the simulations fail) and over the
successful ones only, and the number of failures (no detection, false positive). Printed
to the log and the screen, appended to the per-target results file; returns
`(statistics, failures)` per threshold.
