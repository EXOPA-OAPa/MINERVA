# Eccentric orbits

[← Science cases](../../README.md#science-cases) · Example: [`05_eccentric_orbit`](../../examples/05_eccentric_orbit/config.yaml)

## Science

The eccentricity of a newly discovered transiting planet is usually unknown and assumed
zero when the follow-up is planned. An eccentric orbit spreads the power of the RV signal
into harmonics and weakens the main periodogram peak (paper, Sect. 3.1).

## What MINERVA models

The injected planet is a Keplerian orbit with eccentricity e and argument of periastron ω
(degrees); the recovery always fits a circular orbit, as a blind follow-up would. The
difference between runs with e > 0 and e = 0 is the hidden cost of the circular
assumption.

## Configuration

```yaml
single_target_params:
  eccentricity: 0.5
  periastron: 0.0       # deg: a strongly asymmetric RV curve
```

## Notes

* Compare runs with the same configuration and e = 0 to quantify the additional
  observations and the fraction of failures (paper, Sect. 3.1).
* K_pl includes the factor 1/√(1 − e²).

## See also

[Signal recovery](signal_recovery.md) · [Conventions](../api/conventions.md#orbits-and-signals)
