# Stellar activity from age and rotation

[← Science cases](../../README.md#science-cases) · Examples: [`02_single_target`](../../examples/02_single_target/config.yaml), [`03_known_activity`](../../examples/03_known_activity/config.yaml)

## Science

Activity-induced RV variations depend on the magnetic activity of the star, tied to its
rotation and age. Without previous RV observations, MINERVA derives the expected activity
from stellar properties through empirical calibrations, treating each target as a member
of a synthetic population of the solar neighbourhood (paper, Sect. 2.3).

## What MINERVA models

Chain of relations ([`minerva.star`](../api/star.md)), each skipped if the quantity is given:

1. **age**: uniform between 10⁸ yr and the end of the main sequence (constant birthrate);
2. **rotation regime**: 25% slow, 50% median, 25% fast;
3. **P_rot**: Gallet & Bouvier (2015), log P_rot = A + B log(age);
4. **τ**: convective turnover time from the YaPSI tracks (Spada et al. 2017); Ro = P_rot/τ;
5. **log R'HK**: from Ro (Mamajek & Hillenbrand 2008) or, for M dwarfs, from P_rot
   (Astudillo-Defru et al. 2017);
6. **K_act**: RV dispersion from log R'HK and T_eff (Santos et al. 2000; for K dwarfs the
   GAPS/ArMS exponential calibration, Eq. 3 of the paper, valid up to log R'HK = −3.99).

Stars without a convective envelope (τ ≈ 0, e.g. hot F stars) are inactive: no activity
signal and no activity term in the error bars.

## Fixing what is known

Any quantity of the chain can be given in the catalog or in `single_target_params`; the
chain then starts from it. With `bootstrap: True` the given values are resampled within
their uncertainties at every simulation.

```yaml
single_target_params:
  age: 0.8          # Gyr
  Prot: 12.5        # d
  log_Rhk: -4.45
  K_act:            # empty: from log R'HK
```

## Notes

* The relations are valid in a statistical sense, not for individual stars: for a known
  target, fix the measured P_rot and log R'HK.
* Ages beyond the end of the main sequence would take τ from the subgiant phase of the
  track, inconsistent with a main-sequence T_eff: random ages are limited to the main
  sequence.

## Validation

Paper, Sect. 3.1: the validation targets were run both with the activity parameters
derived by MINERVA and with the literature values, with consistent results.

## See also

[Activity signal](activity_signal.md) · [Parameter uncertainties](parameter_uncertainties.md) · [`minerva.star`](../api/star.md)
