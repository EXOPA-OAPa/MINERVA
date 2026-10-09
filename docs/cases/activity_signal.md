# The injected activity signal

[← Science cases](../../README.md#science-cases) · Example: [`02_single_target`](../../examples/02_single_target/config.yaml)

## Science

Active regions evolve, appear at different latitudes and live from weeks to months: the
activity signal is quasi-periodic, not a stable sinusoid. Its complexity grows with the
activity level (paper, Sect. 2.4).

## What MINERVA models

Sinusoids at P_rot and P_rot(1 + α), with α the relative shear of differential rotation
(beating), computed by `signals.create_stellar_activity_signal_dynamic_complex`:

| Regime | log R'HK | Signal |
|---|---|---|
| low | < −5 | phase-coherent: K_act at P_rot + 0.1 K_act at P_rot(1+α) |
| intermediate | −5 to −4.5 | per season: new phases, K_s ~ N(K_act, 0.1 K_act) ≥ 0.5 m/s, secondary 0.4 K_s |
| high | > −4.5 | per season: as above, secondary 0.8 K_s, first harmonic P_rot/2 with 0.3 K_s |

α is uniform in [0.01, 0.08] for log R'HK ≥ −4.5 and in [0.10, 0.20] otherwise: the
relative shear decreases with the rotation rate (Barnes et al. 2005; Reinhold et al. 2013).

## Notes

* Seasonal resets reproduce the loss of coherence between observing seasons, which are
  longer apart than the lifetime of active regions.
* The recovery fits one sinusoid at the stellar period with free amplitude and phase per
  season: part of the injected complexity (secondary period, harmonic) is not modelled,
  as in a real analysis without Gaussian processes.

## Validation

The amplitudes and phases follow the paper (Eqs. 4–6); the RV-series figure
(`rv_series_pdf`) shows the injected activity of each simulation.

## See also

[Stellar activity](stellar_activity.md) · [Signal recovery](signal_recovery.md)
