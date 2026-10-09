# `minerva.signals` — synthetic RV time series

[← API reference](README.md)

The RV time series is white noise + stellar activity + planet; each epoch has a formal
error bar (paper, Sects. 2.1, 2.2, 2.4, 2.5).

* [`photon_noise`, `white_noise`](#noise)
* [`select_cadence`](#select_cadence)
* [`activity_signal`, `create_stellar_activity_signal_dynamic_complex`](#activity-signal)
* [`add_planetary_signal`](#planet-signal)
* [`error_bars`, `simulate_activity_error`](#error-bars)

---

## Noise

**`photon_noise(var, target_params) → σ_ph [m/s]`**: HARPS-N photon noise from the NSTS
table, at the magnitude of the grid closest to V and for the spectral class of the star.

**`white_noise(var, target_params, sig_ph, rng) → noise`**: Gaussian noise with

  σ = √(σ_instr² + σ_pl² + σ_ph²),

σ_instr = `sig_floor`, σ_pl = `sig_pl` (unresolved additional planets), one value per
epoch of the cadence.

---

## select_cadence

`select_cadence(var, target_params) → epochs`: the USP cadence for P < 2 d, the standard
one otherwise (the target's own file with `real_cadence`).

---

## Activity signal

**`activity_signal(var, target_params, time_f, Prot, K_act, log_Rhk, is_active, rng)`**:
zero for inactive stars, otherwise
`create_stellar_activity_signal_dynamic_complex(Prot, time_f, K_act, log_Rhk, target_params, var, rng)`:

* **Differential rotation.** Secondary period P₂ = P_rot (1 + α), with α uniform in
  [0.01, 0.08] for log R'HK ≥ −4.5 and in [0.10, 0.20] otherwise (relative shear decreases
  with the rotation rate: Barnes et al. 2005; Reinhold et al. 2013).
* **Low activity** (log R'HK < −5): phase-coherent signal over the whole time series,
  K_act at P_rot plus 0.1 K_act at P₂.
* **Intermediate** (−5 ≤ log R'HK ≤ −4.5) and **high activity** (log R'HK > −4.5): in each
  observing season, new random phases and an amplitude K_s ~ N(K_act, 0.1 K_act), at least
  0.5 m/s; secondary component at P₂ with 0.4 K_s (intermediate) or 0.8 K_s (high); for high
  activity, the first harmonic P_rot/2 with 0.3 K_s.

The sinusoids are computed with `keplerian.rv_drive` (circular orbits).

> At exactly log R'HK = −4.5 the star uses the high-activity range of α and the
> intermediate-activity amplitudes.

---

## Planet signal

`add_planetary_signal(target_params, eccentricity, periastron, time_f, noise_rv, RVact, sig_ph) → (RV, err, K_pl, RV_pl)`

Keplerian orbit (`keplerian.rv_drive`, from RadVel, Fulton et al. 2018), with

  K_pl = 28.4 m/s · M★^(−2/3) · (P / 1 yr)^(−1/3) · m_p [M_J] · sin i / √(1 − e²),

sin i = 1 (transiting planets), M★ in M☉ (the planet mass is neglected with respect to
the stellar mass). The time of periastron is computed from the transit time with
`keplerian.timetrans_to_timeperi`; `periastron` (ω) is in degrees. Returns the total RV
(noise + activity + planet), provisional error bars, K_pl and the planet signal.

The planetary mass must be given (catalog or configuration): MINERVA does not estimate it
from the radius.

---

## Error bars

**`error_bars(var, sig_ph, log_Rhk, is_active, time_f, subset_for_simulation, debug_pdf, save_this_plot, show_on_scr, rng) → (err, jitter)`**

Formal uncertainty of each epoch:

  δRV_i = √(σ_instr² + σ_ph² + σ_act,extr²),

where σ_act,extr is the activity-driven RV extraction term (paper, Sect. 2.2), drawn from
the GAPS database of stars of the same spectral class (`subset_for_simulation`) by
`simulate_activity_error`, with a 10% epoch-to-epoch scatter. σ_pl is not included: the
observer does not know it. Inactive stars: √(σ_instr² + σ_ph²) with a 10% scatter.

**`simulate_activity_error(target_rhk, n_obs, df_data, sigma_inst, window=0.15, var_pct=0.20, pdf_obj=None, show_plot=False, rng=None)`**
(called with `window=0.20`, `var_pct=0.10`):

* **local window**: the GAPS stars with log R'HK within ± `window` of the target; a value
  is drawn from their empirical distribution (CDF);
* **broken-line fit**: for log R'HK > −4.5 or beyond the database, a two-slope fit of
  jitter vs log R'HK (`func_broken_line`), and a draw around the prediction with the
  scatter of the residuals.

Returns the error bars and the activity term; optionally a diagnostic figure.
