# Observing effort for a single target

[← Science cases](../../README.md#science-cases) · Example: [`02_single_target`](../../examples/02_single_target/config.yaml)

## Science

Before proposing an RV follow-up, the number of observations needed to measure the mass
of one planet with a given precision, and the risk of not reaching it, can be estimated
from the stellar and planetary parameters and the expected cadence.

## What MINERVA models

The target of `single_target_params` is simulated `n_iterations` times with the same
pipeline as a catalog target. Missing stellar properties are drawn at every simulation;
known ones can be fixed ([known activity](stellar_activity.md#fixing-what-is-known)).

## Configuration

```yaml
read_catalog: False
single_target: True
single_target_params:
  targetname: my_planet
  star_mass: 1.0
  star_temp: 5800
  v_mag: 9.0
  pl_mass_Me: 20
  pl_period: 12.0
  transit_midpoint: 9003.4
run_options:
  debug_plot: True                  # figures of the first 3 simulations
  debug_fit_plot: True
  export_debug_pdf: True
```

## Notes

* The debug PDF (`Plots/DEBUG_<target>.pdf`) shows, for the first simulations, every step:
  periodogram, pre-whitening, fit and the evolution of K, ΔBIC and SNR with the number of
  observations. It is the best way to understand why a simulation fails.
* MINERVA does not estimate planetary masses: for candidates, use a mass–radius relation
  (the paper uses Chen & Kipping 2017, as implemented by Louie et al. 2018).

## Validation

Paper, Sect. 3.1 and Table 2: comparison with the number of observations used for 8
published planets, with the activity parameters derived by MINERVA and fixed to the
literature values.

## See also

[Survey planning](survey_planning.md) · [Stellar activity](stellar_activity.md) · [Signal recovery](signal_recovery.md)
