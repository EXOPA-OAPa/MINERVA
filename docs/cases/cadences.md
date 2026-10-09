# Observing cadences

[← Science cases](../../README.md#science-cases) · Examples: [`01_catalog`](../../examples/01_catalog/config.yaml), [`06_usp_planet`](../../examples/06_usp_planet/config.yaml)

## Science

The number of observations needed depends strongly on the sampling: weather gaps,
seasonal visibility and scheduling constraints produce aliases and power leakage that an
ideal uniform sampling would not have (paper, Sect. 2.1).

## What MINERVA models

Observations are taken from a real cadence and added in time order (`sampling_strategy:
sequential`), 5 at a time. Two default cadences are included:

| Planets | Cadence | Epochs |
|---|---|---|
| P ≥ 2 d | HARPS-N monitoring of HD 164922 (Benatti et al. 2020) | 301 over six years |
| P < 2 d (USP) | HARPS-N observations of TOI-1807 (Nardiello et al. 2022), ≥ 2 per night | 160 over two years |

Any cadence can be used: a text file with one epoch (JD − 2450000) per line, per target
(`Cadence FileName` with `real_cadence: True`) or as default (`name_sampling`,
`name_sampling_usp`).

## Configuration

```yaml
catalog:
  real_cadence: True        # each target: its Cadence FileName (default cadences if empty)
general_variables:
  name_sampling: harps.rv
  name_sampling_usp: time_usp.txt
analysis:
  n_start: 10
  step: 5
  sampling_strategy: sequential
```

## Notes

* The maximum number of observations of a simulation is the length of the cadence: a
  failed simulation saturates there.
* The default cadences are favourable cases (dense, sustained programmes).

## See also

[Survey planning](survey_planning.md) · [Signal recovery](signal_recovery.md)
