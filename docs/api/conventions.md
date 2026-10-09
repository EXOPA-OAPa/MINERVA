# Conventions

[← API reference](README.md)

## Units

| Quantity | Unit |
|---|---|
| times, epochs, transit mid-time | JD − 2450000 [d] |
| periods (planet, rotation) | d |
| radial velocities, K, noise terms, error bars | m/s |
| stellar mass, radius | M☉, R☉ |
| planetary mass, radius | M_J or M⊕, R_J or R⊕ |
| age | Gyr in the inputs, yr inside the code |
| argument of periastron ω | deg (converted to rad internally) |
| chromospheric activity | log R'HK (dex); R5 = 10⁵ R'HK |

## Orbits and signals

* Planet: Keplerian (`keplerian.rv_drive`, from RadVel); the time of periastron is computed from the transit
  mid-time, so for a circular orbit the RV is −K sin(2π(t − T0)/P), zero and decreasing at
  the transit. The fit assumes a circular orbit (paper, Sect. 2.5).
* Activity: sinusoids at P_rot and P_rot(1 + α) (and P_rot/2 for high activity) with
  random phases, one set per observing season for intermediate and high activity.

## Observing seasons

Seasons are separated at the largest gaps of the cadence longer than 50 d that leave at
least 7 epochs per season (`utils.find_seasons_automatic`). They define the seasonal
phases of the injected activity and the seasonal parameters of the simultaneous fit.

## Random numbers

* One `numpy.random.Generator` per target: `default_rng(random_seed + i)` for the i-th
  catalog target, `default_rng(random_seed)` for a single target. Simulations of a target
  share it, in order.
* Inside a simulation the draws follow a fixed order (bootstrap, white noise, age,
  rotation regime, period, activity signal, error bars); the sequential recovery draws
  nothing. Changing the order of the draws changes every result.
* Two runs with the same configuration give identical results; catalog targets are
  independent of each other, so a parallel run gives the results of a sequential run.

## Status codes

For each SNR threshold, at the end of a simulation:

| Code | Meaning | Number of observations stored |
|---|---|---|
| 1 | detection: the fit was physically sound and K/σ_K reached the threshold | at which the threshold was reached |
| 3 | false positive: K/σ_K above the threshold, but the fit was never physically sound | total number of observations |
| 0 | no detection | total number of observations |

## Output files

Written in `output_dir`; overwritten at every run.

**`results_<target>.txt`**: one line per simulation, with the columns

`Target Name, Vmag, Planetary Mass [Mj], Planetary Radius [Rj], Planetary Period [d],
Eccentricity, Periastron [deg], Age [Gyr], Tag, Prot [d], logRHK, Kact [m/s], Kpl [m/s],
std(RV) [m/s], sig_ph [m/s], std(noise_rv) [m/s], Nobs, Niter, Kfit [m/s], Kfit_err [m/s],
SNRfinal, dBIC, jitter_act [m/s], simult_flag, Kstar_fit [m/s], Pstar_fit [d],
status_flag, sig_best, bootstrap_flag`

`Nobs`, `status_flag` and `sig_best` are lists, one value per SNR threshold
(e.g. `"[30, 35]"`). At the end: the statistics of the target
([`print_execution_summary`](running.md#statistics-over-the-simulations)).

**`results_<catalog>.txt`** (catalog runs): one line per target, in catalog order, with
median, 16th and 84th percentiles over all simulations (`ALL`) and over the successful
ones (`SUCC`) for each threshold, the number of simulations and the failures per
threshold.

**`Log_*.txt`**: detailed log (one per catalog, or one per target in parallel runs).

**`Plots/`**: `plot_series_*.pdf` (RV series, `rv_series_pdf`), `DEBUG_<target>.pdf`
(debug figures).
