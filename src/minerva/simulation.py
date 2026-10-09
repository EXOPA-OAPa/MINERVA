"""
One Monte Carlo realisation of MINERVA (paper, Fig. A.1):

1. synthetic star: spectral type, age, rotation, activity level and amplitude (star.py);
2. synthetic RV time series: white noise, activity signal, planet, error bars (signals.py);
3. sequential recovery: starting from n_start observations and adding `step` at a time,
   GLS periodogram and pre-whitening (periodogram.py), planet-only or simultaneous fit
   (fitting.py), detection criteria (detection.py), until the highest SNR threshold is
   reached or the observations are exhausted.
"""
import copy
import sys

import numpy as np

from . import detection, fitting, periodogram, plotting, signals, star
from . import utils as minerva_utils
from .tracks import fap_levels


def get_indices_at_step(j, n_start, step, N_total, strategy='sequential', seed=42):
        """
        Indices of the observations used at step j (n_start + j * step, at most N_total).
        - strategy='sequential': the first observations in time order (realistic campaign);
        - strategy='random': observations drawn at random (fixed seed, so the sets are nested).
        The indices are returned sorted, to keep the time order of the data.
        """

        # 1. Number of observations at this step
        n_current = n_start + (j * step)

        # Never more than the available observations
        if n_current > N_total:
            n_current = N_total

        # 2. Sequential selection
        if strategy == 'sequential':
            all_indices = np.arange(N_total)
            ind = all_indices[:n_current]

        # 3. Random selection
        elif strategy == 'random':
            rng = np.random.default_rng(seed)
            all_indices = rng.permutation(N_total)
            ind = all_indices[:n_current]

        else:
            raise ValueError("Unsupported strategy: choose 'sequential' or 'random'.")

        return np.sort(ind)


def setup_seasons_logic(var, target_params):
    """Load the cadence of the target (real cadence from the catalog if requested)
    and split it into observing seasons."""
    # Cadence files: the target's own file when real_cadence is set, the default ones otherwise
    if var.config['real_cadence']:
        input_s = target_params.get('catalogs', var.inputfile_sampling)
        input_s_usp = target_params.get('catalogs', var.inputfile_sampling_usp)
        # Fallback when the catalog has no cadence file for this target
        if input_s is None: input_s = var.inputfile_sampling
        if input_s_usp is None: input_s_usp = var.inputfile_sampling_usp
    else:
        input_s = var.inputfile_sampling
        input_s_usp = var.inputfile_sampling_usp

    var.time_sampling, var.time_sampling_usp = minerva_utils.import_file_sampling(
        InDir="", inputfile_sampling=var.resolve(input_s), inputfile_sampling_usp=var.resolve(input_s_usp)
    )

    # Observing seasons of the cadence used for this planet (USP cadence for P < 2 d)
    if target_params['pl_period'] >= 2:
        xs = minerva_utils.find_seasons_automatic(var.time_sampling)
        var.seasons = minerva_utils.define_seasons(var.time_sampling, xs)
        prefix = 's'
    else:
        xs_usp = minerva_utils.find_seasons_automatic(var.time_sampling_usp)
        var.seasons = minerva_utils.define_seasons(var.time_sampling_usp, xs_usp)
        prefix = 's_usp'

    # Attributes var.s1, var.s2, ... (or var.s1_usp, ...)
    for ii, s in enumerate(var.seasons):
        attr_name = f"{prefix}{ii+1}" if prefix == 's' else f"s{ii+1}_usp"
        setattr(var, attr_name, s)


def run_simulation(var, target_params, niter, pdf, debug_pdf, rng):
    """One Monte Carlo simulation of a target.

    Returns log_Rhk, the number of observations needed for each SNR threshold (failures
    saturated to the total number of observations) and the status for each threshold
    (1 detection, 3 false positive, 0 no detection); None if the star cannot be simulated.
    """
    # Work on a copy: every simulation starts from the nominal parameters (the bootstrap
    # must not accumulate across simulations)
    target_params = copy.deepcopy(target_params)

    # Run options
    BOOTSTRAP          = var.options["bootstrap"]
    DEBUG_PLOT         = var.options["debug_plot"]
    DEBUG_FIT_PLOT     = var.options["debug_fit_plot"]
    DEBUG_PDF_ALL_ITER = var.options["debug_pdf_all_iter"]

    # =========================================================================
    # 1. Synthetic star
    # =========================================================================
    if BOOTSTRAP == True:
        star.bootstrap_parameters(target_params, rng)

    if not star.check_spectral_type(var, target_params):
        return None

    # Noise: photon noise of the star and white noise of the time series
    sig_ph = signals.photon_noise(var, target_params)
    noise_rv = signals.white_noise(var, target_params, sig_ph, rng)

    age = star.stellar_age(var, target_params, rng)
    tag = star.rotation_class(target_params, rng)
    Prot, log_Prot = star.rotation_period(target_params, age, tag, rng)
    tau, R0 = star.rossby_number(var, target_params, age, Prot)
    log_Rhk, is_active = star.activity_level(target_params, tau, R0, log_Prot, Prot)
    K_act = star.activity_amplitude(target_params, log_Rhk, is_active)

    # =========================================================================
    # 2. Synthetic RV time series
    # =========================================================================
    time_f = signals.select_cadence(var, target_params)
    RVact = signals.activity_signal(var, target_params, time_f, Prot, K_act, log_Rhk, is_active, rng)
    eccentricity = target_params['eccentricity']
    periastron = target_params['periastron']
    RVfinal, err, kpl, RVpl = signals.add_planetary_signal(target_params,eccentricity,periastron,time_f,noise_rv,RVact,sig_ph)
    # GAPS stars of the same spectral class, for the activity term of the error bars
    subset_for_simulation = var.df_activity [var.df_activity ['SP_CLASS'] == target_params['sp'][0]].dropna(subset=['logRHK', 'Jitter_ms'])

    # Debug figures are built only if they will be saved or shown:
    # make_plots     -> saved for the first 3 simulations only (or shown on screen)
    # make_plots_all -> GLS and fit diagnostics, saved for every simulation if debug_pdf_all_iter
    show_on_scr = getattr(var, 'SHOW_ON_SCREEN', False)
    save_this_plot = (debug_pdf is not None) and (niter < 3)
    make_plots     = save_this_plot or show_on_scr
    make_plots_all = ((debug_pdf is not None) or show_on_scr) if DEBUG_PDF_ALL_ITER else make_plots

    err, jitter = signals.error_bars(var, sig_ph, log_Rhk, is_active, time_f, subset_for_simulation,
                                     debug_pdf, save_this_plot, show_on_scr, rng)

    print(f"\n--- INJECTED STELLAR PARAMETERS ---")
    print(f"   Sp. Class : {target_params['sp']}")
    print(f"   Age       : {age/1e9:.3f} Gyr")
    print(f"   Rot. Tag  : {tag}")
    print(f"   P_rot     : {Prot:.2f} d")
    print(f"   log(R'hk) : {log_Rhk:.3f}")
    print(f"   K_act     : {K_act:.2f} m/s")
    print(f"   jitter_act: {jitter:.4f} m/s")
    print(f"--- PLANET PARAMETERS ---")
    print(f"   Period       : {target_params['pl_period']} d")
    print(f"   Eccentricity : {eccentricity}")
    print(f"   Periastron   : {periastron} deg")
    print(f"------------------------------------")

    # RV time series figure (run_options.rv_series_pdf: True, False, or the number of simulations)
    rv_series = var.options["rv_series_pdf"]
    if rv_series is True or (rv_series is not False and rv_series is not None and niter < int(rv_series)):
        plotting.plot_rv_series(
                pdf=pdf,
                time_f=time_f,
                RVfinal=RVfinal,
                err=err,
                RVact=RVact,
                noise_rv=noise_rv,
                target_params=target_params,
                Prot=Prot,
                niter=niter,
                fap_levels_func=fap_levels
            )

    # =========================================================================
    # 3. Sequential recovery
    # =========================================================================
    n_start = int(var.analysis["n_start"])   # observations at the first step
    step = int(var.analysis["step"])         # observations added at each step
    N_total = len(time_f)
    max_iter = (N_total - n_start) // step
    # 'sequential' (realistic campaign) or 'random'
    sampling_strategy = var.analysis["sampling_strategy"]
    ind = np.array([], dtype=int)

    # SNR thresholds, in increasing order (e.g. [5, 5.5])
    if isinstance(var.sig_best, (list, tuple, np.ndarray)):
        thresholds = sorted(list(var.sig_best))
    else:
        thresholds = [float(var.sig_best)]
    # Number of observations at which each threshold was reached (0 = not reached)
    nobs_reached = {th: 0 for th in thresholds}
    # History of the fit statistics, for the debug figure
    history_N = []
    history_dBIC = []
    history_Kpl = []
    history_jitter = []
    history_chi2 = []
    # Default status (no detection), in case the loop ends early
    status_dict = {th: 0 for th in thresholds}
    # Last intruder period tested by the pre-whitening (kept across steps)
    P_intruder = None
    for j in range(max_iter+1):
        ind = get_indices_at_step(j, n_start, step, N_total, strategy=sampling_strategy)
        print(f"\n{'-'*50}\n>>> Step {j}: analysis with {len(ind)} points ({sampling_strategy})")
        t_sel = time_f[ind]
        y_sel = RVfinal[ind]
        err_sel = err[ind]
        # Points added at this step (highlighted in the debug figures)
        if j == 0:
            n_new_points = 0
        else:
            n_new_points = min(step, len(ind) - (n_start + (j-1)*step))
        idx_split = len(ind) - n_new_points
        # FAP levels for stellar and planetary peaks
        percFAP_star = float(var.analysis["fap_star"])
        percFAP_planet = float(var.analysis["fap_planet"])

        # --- GLS periodogram ---
        fbeg, fend = periodogram.frequency_range(time_f, ind, target_params)
        try:
            per = var.Gls((t_sel, y_sel, err_sel), verbose=False, fbeg=fbeg, fend=fend)
        except Exception as e:
            print(f"GLS computation failed: {e}")
            per = None
            break
        if per is not None:
            (freq, power, threshold_star, threshold_planet, sorted_periods,
             sorted_powers) = periodogram.significant_peaks(per, percFAP_star, percFAP_planet, ind)
            if make_plots_all:
                plotting.gls_peaks(**plotting.pick(locals(), plotting.gls_peaks))

            # --- Planet and intruder peaks, pre-whitening ---
            (freq_res, f_planet_known, idx_planet_candidates, idx_intruders, peak_uncert,
             planet_is_highest) = periodogram.classify_peaks(t_sel, target_params, sorted_periods)
            (intruder_accepted, use_simultaneous, P_star_final, P_intruder) = periodogram.prewhitening(
                per, t_sel, y_sel, err_sel, fbeg, fend, sorted_periods, sorted_powers, idx_intruders,
                idx_planet_candidates, f_planet_known, freq_res, peak_uncert, percFAP_star, percFAP_planet,
                threshold_star, threshold_planet, P_intruder, DEBUG_PLOT, make_plots, var, niter, debug_pdf, Prot)
            if not intruder_accepted:
                P_star_final, use_simultaneous = periodogram.reverse_prewhitening(
                    per, t_sel, y_sel, err_sel, fbeg, fend, target_params, planet_is_highest, percFAP_star,
                    freq_res, P_star_final, use_simultaneous, DEBUG_PLOT, make_plots, var, niter, debug_pdf, Prot)

            # --- Fit: always on the raw data, never on the pre-whitening residuals ---
            t_fit = t_sel
            rv_fit = y_sel
            err_fit = err_sel
            # Planet ephemerides, fixed in the fit
            T0_PLANET = target_params['transit_midpoint']
            SNR_final, chisq_red, delta_bic = 0.0, 0.0, 0.0
            is_detected = False
            popt = None
            fit_desc = "NONE"
            y_model_final = None
            K_star_fit = 0.0
            P_PLANET = target_params['pl_period']
            if use_simultaneous:
                # The stellar period must be defined: last intruder tested, or the true Prot
                if P_star_final is None:
                    if P_intruder is not None:
                        P_star_final = P_intruder
                    else:
                        P_star_final = Prot
                # Seasons of the activity model in the current data
                all_season_indices = var.seasons
                (current_seasons_idx, N_seasons, season_boundaries,
                 fit_desc) = fitting.activity_seasons(is_active, log_Rhk, t_fit, time_f, ind, all_season_indices)
                if N_seasons == 0:
                    print("   ⚠️ No season with enough points: falling back to the planet-only fit.")
                    use_simultaneous = False
            if use_simultaneous:
                # Simultaneous fit: planet + seasonal activity at the stellar period
                omega_star = 2 * np.pi / P_star_final
                model_simultaneous_seasonal = fitting.seasonal_simultaneous_model(P_PLANET, T0_PLANET,
                                                                                  season_boundaries, omega_star)
                model_star_only_seasonal = fitting.seasonal_star_model(season_boundaries, omega_star)
                (popt, K_fit, K_err_formal, dof, k_params, bic_null, y_model_final, rv_plot_data,
                 K_star_fit, use_simultaneous) = fitting.fit_simultaneous(
                    t_fit, rv_fit, err_fit, current_seasons_idx, N_seasons,
                    model_simultaneous_seasonal, model_star_only_seasonal, P_PLANET, T0_PLANET, K_star_fit)
            # Planet-only fit (also when the simultaneous fit failed)
            if not use_simultaneous:
                (fit_desc, popt, K_fit, K_err_formal, dof, k_params, bic_null, y_model_final,
                 rv_plot_data) = fitting.fit_planet_only(t_fit, rv_fit, err_fit, P_PLANET, T0_PLANET)
            # Both fits failed: add points and try again
            if popt is None or y_model_final is None:
                print("   ⚠️ Both fits failed to converge. Adding data...")
                continue
            if DEBUG_FIT_PLOT and (popt is not None) and make_plots_all:
                plotting.fit_diagnostics(**plotting.pick(locals(), plotting.fit_diagnostics))

            # --- Fit statistics: chi2, error rescaling, SNR, Delta BIC ---
            (residuals, chisq_red, scale_factor, K_err_final, SNR_final, delta_bic) = detection.fit_statistics(
                rv_fit, y_model_final, err_fit, t_fit, dof, k_params, K_fit, K_err_formal, bic_null)
            sys.stdout.set_mode(to_file=True, to_screen=False)
            print(f"\n--- {fit_desc} ---")
            print(f"   K_planet = {K_fit:.2f} +/- {K_err_final:.2f} m/s")
            print(f"   Chi2_red = {chisq_red:.2f} -> Scale Factor = {scale_factor:.2f}")
            print(f"   SNR = {SNR_final:.2f}")
            print(f"   DeltaBIC (vs Flat) = {delta_bic:.2f}")
            history_N.append(len(ind))
            history_dBIC.append(delta_bic)
            history_Kpl.append(K_fit)
            history_jitter.append(jitter)
            history_chi2.append(chisq_red)

            # --- Detection criteria, for every SNR threshold ---
            (is_detected, SNR_final, rejection_reasons, title_color) = detection.evaluate_detection(
                K_fit, K_err_final, SNR_final, delta_bic, kpl, use_simultaneous, is_active, P_star_final, Prot,
                thresholds, nobs_reached, ind, var.analysis["delta_bic_threshold"])
            sys.stdout.set_mode(to_file=True, to_screen=False)
            if DEBUG_FIT_PLOT and popt is not None and make_plots:
                plotting.phase_folded_fit(**plotting.pick(locals(), plotting.phase_folded_fit))

            # --- End of the simulation: detection reached, or no observations left ---
            is_last_iteration = (j == max_iter)
            if is_detected or is_last_iteration:
                status_dict = detection.save_simulation_result(
                    var, target_params, time_f, thresholds, nobs_reached, SNR_final, use_simultaneous, K_star_fit,
                    P_star_final, BOOTSTRAP, age, tag, Prot, log_Rhk, K_act, kpl, RVfinal, sig_ph, noise_rv, niter,
                    K_fit, K_err_final, delta_bic, jitter, make_plots, debug_pdf,
                    history_N, history_dBIC, history_Kpl, history_jitter, history_chi2)
                if is_detected:
                    print(f"\n SIMULATION {niter} COMPLETED: DETECTION")
                    break
                elif is_last_iteration:
                    print(f"\n SIMULATION {niter} FAILED: all the observations used.")
            else:
                if len(rejection_reasons) > 0:
                    reasons_msg = " | ".join(rejection_reasons)
                    print(f"   ❌ REJECTED: {reasons_msg}")
                print(f"   -> Adding points and trying again... \n")

    # Failures (no detection or false positive) are saturated to the total number of
    # observations: "at least N_total observations are needed"
    N_avail = len(time_f)
    nobs_out = {th: (nobs_reached[th] if status_dict[th] == 1 else N_avail) for th in thresholds}
    return log_Rhk, nobs_out, status_dict
