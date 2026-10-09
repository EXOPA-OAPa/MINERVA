# Noise budget and error bars

[← Science cases](../../README.md#science-cases) · Example: any

## Science

The scatter of the RVs and the uncertainties attached to them are different things: the
injected noise is the true scatter of the data, the error bars are what the observer
knows and uses to weight the data (paper, Sect. 2.2).

## What MINERVA models

* **Injected white noise**: σ = √(σ_instr² + σ_ph² + σ_pl²), with σ_instr the instrumental
  noise (`sig_floor`, 0.5 m/s for HARPS-N), σ_ph the photon noise of the star (HARPS-N
  exposure time calculator, by V magnitude and spectral class) and σ_pl the unresolved
  signal of additional planets (`sig_pl`, 2 m/s). Stellar activity is added separately,
  as a correlated signal.
* **Error bars**: δRV = √(σ_instr² + σ_ph² + σ_act,extr²), where σ_act,extr is the
  degradation of the RV extraction in active stars (distorted line profiles), drawn from
  the GAPS database of 250 FGKM stars for the spectral class and the log R'HK of the star,
  with a 10% epoch-to-epoch scatter. σ_pl is not in the error bars.
* The formal uncertainty on K is inflated by √χ²_red when the residuals exceed the error
  bars (never deflated).

## Configuration

```yaml
general_variables:
  sig_floor: 0.5
  sig_pl: 2
```

Other instruments: replace the photon-noise table (`inputfile_nsts`) and `sig_floor`.

## Notes

* The activity term of the error bars and the injected activity signal come from two
  different calibrations (GAPS jitter database; Santos et al. / GAPS-ArMS K_act).

## See also

[`minerva.signals`](../api/signals.md#error-bars) · [Survey planning](survey_planning.md)
