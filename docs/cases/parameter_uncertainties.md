# Uncertainties on the stellar parameters

[← Science cases](../../README.md#science-cases) · Example: [`04_bootstrap`](../../examples/04_bootstrap/config.yaml)

## Science

For a known target, the stellar parameters have uncertainties that propagate into the
activity level and into the observing effort; the risk of under-sampling an active star
is often dominated by them (paper, Sect. 3.1).

## What MINERVA models

With `bootstrap: True`, at every simulation stellar mass, T_eff, radius, age, P_rot and
log R'HK are drawn from Gaussians centred on their values, with their uncertainties
(defaults: 10%, 100 K, 0.05 dex), within physical bounds; every simulation starts from the
nominal values.

## Configuration

```yaml
run_options:
  bootstrap: True
single_target_params:
  star_mass: 0.55
  err_star_mass: 0.03
  log_Rhk: -4.9
  err_log_Rhk: 0.1
```

## Notes

* The required number of observations of successful simulations is usually stable; the
  bootstrap mainly changes the fraction of failures for active or low-SNR targets.

## See also

[Stellar activity](stellar_activity.md) · [`star.bootstrap_parameters`](../api/star.md#bootstrap_parameters)
