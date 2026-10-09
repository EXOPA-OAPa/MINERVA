# Sequential recovery and detection criteria

[← Science cases](../../README.md#science-cases) · Test: [regression](../api/testing.md#regression-test)

## Science

A real campaign stops when the planet is detected with the required precision. MINERVA
reproduces this with an automated, unsupervised analysis repeated as observations
accumulate, which gives a conservative estimate (no Gaussian processes or other advanced
activity mitigation; paper, Sect. 2.6).

## What MINERVA models

At each step (10 observations, then 5 more at a time):

1. GLS periodogram; peaks above the stellar FAP threshold are planet candidates (near the
   known period) or intruders;
2. intruders are pre-whitened in order of power; one that improves the planet peak is the
   stellar signal; otherwise, if the planet dominates, the planet is subtracted and the
   residuals are searched for the star (reverse pre-whitening);
3. fit with P and T0 fixed: planet + seasonal activity (if a stellar period was found) or
   planet only, and a null hypothesis without planet;
4. detection: ΔBIC > threshold, K > 0, K consistent with the injected one within 3σ,
   stellar period equal to P_rot or a harmonic (simultaneous fits), K/σ_K above each SNR
   threshold.

The simulation stops when the highest SNR threshold is reached.

## Configuration

```yaml
general_variables:
  sig_best: [5, 5.5]          # SNR thresholds
analysis:
  fap_star: 0.10
  fap_planet: 0.50
  delta_bic_threshold: 6.0
```

## Notes

* The check on the injected K and on P_rot uses the truth of the simulation: it rejects
  detections at the wrong amplitude or with the wrong stellar period, which a real analysis
  could not recognise.
* A status 3 (false positive) means that the SNR exceeded the threshold but the fit was
  never physically sound.

## See also

[`minerva` recovery API](../api/recovery.md) · [Cadences](cadences.md)
