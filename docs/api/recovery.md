# Signal recovery — `periodogram`, `fitting`, `detection`, `fast_gls`

[← API reference](README.md)

At every step of the sequential analysis (paper, Sect. 2.6) `run_simulation` calls, in
order: the GLS periodogram, `significant_peaks`, `classify_peaks`, `prewhitening`
(and `reverse_prewhitening` if no intruder is accepted), one of the two fits,
`fit_statistics` and `evaluate_detection`.

* [`periodogram`](#periodogram)
* [`fitting`](#fitting)
* [`detection`](#detection)
* [`fast_gls`](#fast_gls)

---

## periodogram

**`frequency_range(time_f, ind, target_params) → (f_beg, f_end)`**: from 1/baseline of
the current data to their Nyquist frequency 0.5/median(Δt), extended to 1.5 f_planet if
needed.

**`significant_peaks(per, percFAP_star, percFAP_planet, ind)`**: power thresholds of the
two FAP levels (`analysis.fap_star`, `analysis.fap_planet`) and the local maxima above the
stellar threshold, sorted by decreasing power.

**`classify_peaks(t_sel, target_params, sorted_periods)`**: a peak within 1.5 frequency
resolutions (1/baseline) of the known planet frequency is a **planet candidate**; daily
aliases of the planet (|1.0027 c/d − f| ≈ f_planet within 2 resolutions) are ignored; every
other peak is an **intruder** (typically rotation or its harmonics). Also returns whether
the planet is the highest peak.

**`prewhitening(...) → (intruder_accepted, use_simultaneous, P_star, P_intruder)`**: the
intruders are tested in decreasing order of power. Each one is fitted with a sinusoid and
subtracted; in the GLS of the residuals the planet peak is searched in its window:

* removing the intruder **destroys** the planet (power drops by more than 30%, or falls
  below the planet FAP threshold) → rejected, next intruder; unless the intruder is
  dominant (power > 0.40: "blind override"), in which case it is the stellar signal;
* the planet **improves** (higher power, or ≥ 90% of the power and half the distance from
  the known frequency) → accepted: simultaneous fit with P_star = P_intruder;
* otherwise, no benefit → next intruder.

**`reverse_prewhitening(...) → (P_star, use_simultaneous)`**: when no intruder is
accepted and the planet is the highest peak, a preliminary planet model (P and T0 fixed)
is subtracted; if the highest peak of the residuals is above the stellar FAP threshold it
is the stellar period and triggers the simultaneous fit. Otherwise: planet-only fit.

---

## fitting

Fits always use the raw data, never the pre-whitening residuals; P and T0 of the planet
are fixed to the ephemerides. `scipy.optimize.curve_fit` (trust-region reflective with
bounds).

| Function | Description |
|---|---|
| `model_planet_only(t, gamma, K_pl, P_pl, T0_pl)` | circular orbit γ − K sin(2π (t − T0)/P) |
| `seasonal_simultaneous_model(P, T0, season_boundaries, ω_star)` | returns the model: −K sin(2π(t−T0)/P) + Σ_s [γ_s + A_s sin ω★t + B_s cos ω★t] in each season; parameters K, then (γ, A, B) per season |
| `seasonal_star_model(season_boundaries, ω_star)` | the same without the planet (null hypothesis) |
| `activity_seasons(is_active, log_Rhk, t_fit, time_f, ind, all_season_indices)` | one global season for inactive stars and log R'HK < −5; otherwise the observing seasons with at least 3 points |
| `fit_simultaneous(...)` | planet + seasonal activity (K ≥ 0) and star-only null hypothesis; k = 1 + 3 N_seasons, k_null = 3 N_seasons; on failure, planet-only fit |
| `fit_planet_only(t_fit, rv_fit, err_fit, P, T0)` | (γ, K ≥ 0), k = 2; null hypothesis: weighted mean, k_null = 1 |

---

## detection

**`fit_statistics(...)`**: χ²_red = χ²/dof; σ_K = σ_K,formal · max(1, √χ²_red);
SNR = K/σ_K; BIC = χ² + k ln N; ΔBIC = BIC_null − BIC_full (positive when the planet is
favoured).

**`evaluate_detection(..., delta_bic_threshold=6.0)`**: the fit is **physically sound**
when

* ΔBIC > `analysis.delta_bic_threshold` and K > 0;
* |K_fit − K_injected| ≤ 3 σ_K;
* for simultaneous fits of active stars, the stellar period matches P_rot, P_rot/2 or
  P_rot/3 within 15% (at most 2.5 d).

For each SNR threshold, the number of observations at which K/σ_K first exceeds it is
recorded, and reset if the detection is lost. The simulation is complete when the highest
threshold is reached.

**`save_simulation_result(...)`**: status per threshold (1 detection, 3 false positive:
SNR above the threshold but the fit was never physically sound, 0 no detection; failures
saturated to the total number of observations) and one line in `results_<target>.txt`
([conventions](conventions.md#output-files)).

---

## fast_gls

`minerva.fast_gls.FastGls(*args, method="exact" | "recurrence", **kwargs)` — subclass of the
GLS of Zechmeister & Kürster (2009) (`minerva/_vendor/gls.py`, unchanged) with the same
equations and interface. The sums over the data are computed for blocks of frequencies
with matrix products instead of a Python loop:

* `exact`: cosines and sines computed directly; powers equal to the reference within
  ~10⁻¹³ (order of the sums);
* `recurrence`: cosines and sines by trigonometric recurrences, restarted every 1000
  frequencies (as the `fast` option of the reference code); ~10⁻⁹.

`make_gls(method)` returns the class for `analysis.gls_method`
(`reference` = the original implementation).
