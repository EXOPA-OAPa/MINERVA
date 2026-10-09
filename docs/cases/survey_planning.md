# Planning a survey on a catalog of targets

[← Science cases](../../README.md#science-cases) · Example: [`01_catalog`](../../examples/01_catalog/config.yaml)

## Science

Large characterisation programmes (Ariel, PLATO follow-up) need the telescope time
required to measure the masses of hundreds or thousands of transiting planets before the
observations start. MINERVA estimates, for every target, the distribution of the number
of RV observations needed to reach a given precision on K, from which the observing
nights follow (paper, Sect. 3.2: 1541 targets of the Ariel Mission Candidate Sample).

## What MINERVA models

Each target of the catalog is simulated `n_iterations` times; every simulation draws the
unknown stellar properties (age, rotation, activity) from their distributions and runs the
sequential recovery. Per target the result is the median and the 16th–84th percentiles of
the number of observations, over all simulations (failures count as the total number of
observations) and over the successful ones, plus the number of failures.

## Configuration

```yaml
read_catalog: True
single_target: False
catalog:
  namecatalog: my_catalog.csv      # columns: docs/api/parameters.md#catalog-columns
  real_cadence: False               # True: a cadence file per target
general_variables:
  sig_best: [3.33]                  # 30% mass precision
  n_iterations: 100
run_options:
  n_workers: auto                   # targets in parallel
  rv_series_pdf: False
  export_debug_pdf: False
```

## Notes

* The **ALL** statistics include failures and are the conservative estimate for a blind
  survey; **SUCC** describes the favourable realisations, comparable with published
  detections (which are biased towards successes).
* Medians are not additive: the total effort of a sample is computed from the
  distributions (e.g. bootstrap over the simulations of every target), not by summing the
  medians.
* Converting observations into nights requires exposure times and overheads per
  magnitude (paper, Table 3).

## Validation

The paper compares MINERVA with the observing effort of 8 published planets (Sect. 3.1,
Table 2). The regression test checks that the pipeline reproduces a reference run of the
example catalog exactly.

## See also

[Single target](single_target.md) · [Cadences](cadences.md) · [Running faster](../../README.md#running-faster)
