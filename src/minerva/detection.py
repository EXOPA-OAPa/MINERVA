"""
Fit statistics, detection criteria and per-simulation results (paper, Sect. 2.6).
"""
import numpy as np

from . import plotting


def fit_statistics(rv_fit, y_model_final, err_fit, t_fit, dof, k_params, K_fit, K_err_formal, bic_null):
    """Goodness of fit and significance of the planet.

    The formal uncertainty on K is scaled by max(1, sqrt(chi2_red)); the SNR is K / sigma_K;
    Delta BIC = BIC(null hypothesis) - BIC(full model), positive when the planet is favoured."""
    residuals = rv_fit - y_model_final
    chisq = np.sum((residuals / err_fit)**2)
    chisq_red = chisq / dof if dof > 0 else 1.0

    # Error inflation: the formal uncertainty is never reduced, only increased
    scale_factor = max(1.0, np.sqrt(chisq_red))

    K_err_final = K_err_formal * scale_factor
    SNR_final = K_fit / K_err_final if K_err_final > 0 else 0

    # BIC = chi2 + k ln(N) of the full model
    N_pts = len(t_fit)
    bic_full = chisq + k_params * np.log(N_pts)

    # Signed: positive when the model with the planet is favoured over the null hypothesis
    delta_bic = bic_null - bic_full

    return residuals, chisq_red, scale_factor, K_err_final, SNR_final, delta_bic


def evaluate_detection(K_fit, K_err_final, SNR_final, delta_bic, kpl, use_simultaneous, is_active, P_star_final, Prot,
                       thresholds, nobs_reached, ind, delta_bic_threshold=6.0):
    """Detection criteria (paper, Sect. 2.6).

    The fit is physically sound when Delta BIC is above threshold, K > 0, the fitted K is within
    3 sigma of the injected one and, for simultaneous fits of active stars, the stellar period
    matches Prot, Prot/2 or Prot/3 (15% tolerance, at most 2.5 d). nobs_reached (modified in
    place) records, for each SNR threshold, the number of points at which it was first reached;
    it is reset if the detection is lost. The simulation is complete when the highest threshold
    is reached."""
    THRESH_BIC = float(delta_bic_threshold)  # minimum Delta BIC (analysis.delta_bic_threshold in the config)
    is_valid_K = (K_fit > 0)
    rejection_reasons = []  # reasons for rejecting the detection at this step (for the log)

    # 1. K must be positive
    if K_fit <= 0:
        print(f"   ⚠️ K zero or negative ({K_fit:.2f}). Detection rejected.")
        SNR_final = 0
        is_valid_K = False
    else:
        is_valid_K = True

    # 2. Statistical evidence against the null hypothesis
    if delta_bic <= THRESH_BIC:
        rejection_reasons.append(f"Insufficient statistical evidence H1 vs H0 (DeltaBIC={delta_bic:.2f} <= {THRESH_BIC})")

    # 3. Physical consistency of the model
    is_physically_sound = False
    if (delta_bic > THRESH_BIC) and is_valid_K:
        # 3a. The fitted K must match the injected one within 3 sigma
        diff_k_pl = np.abs(K_fit - kpl)
        is_k_pl_correct = diff_k_pl <= (3.0 * K_err_final)

        if not is_k_pl_correct:
            rejection_reasons.append(f"K_fit inconsistent with the injected value (distance: {diff_k_pl:.2f} m/s > 3 sigma)")

        # 3b. Stellar period. The check makes sense only if a stellar signal was actually injected:
        #     for a star without activity any "intruder" is noise, and there is no true period.
        if use_simultaneous and is_active:
            diff_p1 = np.abs(P_star_final - Prot)
            diff_p2 = np.abs(P_star_final - (Prot / 2.0))
            diff_p3 = np.abs(P_star_final - (Prot / 3.0))

            # Tolerance: 15% of the period, but never more than 2.5 days
            # (avoids accepting differences of several days for very slow rotators)
            max_tol_p1 = min(0.15 * Prot, 2.5)
            max_tol_p2 = min(0.15 * (Prot / 2.0), 2.5)
            max_tol_p3 = min(0.15 * (Prot / 3.0), 2.5)

            is_p1 = diff_p1 <= max_tol_p1
            is_p2 = diff_p2 <= max_tol_p2
            is_p3 = diff_p3 <= max_tol_p3

            is_p_star_correct = is_p1 or is_p2 or is_p3
            is_star_correct = is_p_star_correct

            if is_p1: harm_str = "fundamental"
            elif is_p2: harm_str = "harmonic P/2"
            elif is_p3: harm_str = "harmonic P/3"
            else: harm_str = "wrong"

            print(f"   [Star check] P_fit: {P_star_final:.2f}d (true: {Prot:.2f}d -> found: {harm_str})")

            if not is_p_star_correct:
                rejection_reasons.append(f"Wrong stellar period: found {P_star_final:.2f}d (out of tolerance)")
        else:
            is_star_correct = True

        is_physically_sound = is_k_pl_correct and is_star_correct

    # 4. SNR milestones, updated at every step
    if is_physically_sound:
        title_color = 'olive'

        for th in thresholds:
            if SNR_final >= th:
                # First time above the threshold (or first time after losing it): record the number of points
                if nobs_reached[th] == 0:
                    nobs_reached[th] = len(ind)
                    print(f"   ✅ MILESTONE: SNR >= {th} reached with {len(ind)} points.")
                else:
                    # Already reached: keep the smallest number of points recorded
                    print(f"   ✅ MILESTONE: SNR >= {th} already reached with {nobs_reached[th]}/{len(ind)} points.")
            else:
                # The SNR dropped below the threshold: the detection is lost
                if nobs_reached[th] > 0:
                    print(f"   📉 SIGNAL LOST: SNR dropped to {SNR_final:.2f} (< {th}). Milestone reset.")
                    nobs_reached[th] = 0

        # The sequential loop stops when the highest threshold is reached
        max_th = thresholds[-1]
        if SNR_final >= max_th:
            print(f"   ✅✅✅ FINAL MILESTONE: SNR >= {max_th} reached!")
            is_detected = True
            title_color = 'green'
        else:
            is_detected = False  # keep adding points
            rejection_reasons.append(f"Highest SNR threshold not reached yet (SNR={SNR_final:.2f} < {max_th})")

    else:
        # The model is no longer physically sound: all temporary milestones are lost
        for th in thresholds:
            if nobs_reached[th] > 0:
                print(f"   💥 PHYSICALLY WRONG MODEL: all temporary milestones lost!")
                nobs_reached[th] = 0

        is_detected = False
        title_color = 'red'

    return is_detected, SNR_final, rejection_reasons, title_color


def save_simulation_result(var, target_params, time_f, thresholds, nobs_reached, SNR_final, use_simultaneous, K_star_fit,
                           P_star_final, BOOTSTRAP, age, tag, Prot, log_Rhk, K_act, kpl, RVfinal, sig_ph,
                           noise_rv, niter, K_fit, K_err_final, delta_bic, jitter, make_plots, debug_pdf,
                           history_N, history_dBIC, history_Kpl, history_jitter, history_chi2):
    """End of one simulation: status per SNR threshold and one line in the per-target results file.

    Status: 1 = detection, 3 = false positive (SNR above threshold but fit not physically sound),
    0 = no detection. Failed thresholds are saturated to the total number of observations
    (nobs_reached is modified in place)."""
    max_pts = len(time_f)
    status_dict = {}

    for th in thresholds:
        if nobs_reached[th] > 0 and nobs_reached[th] <= max_pts:
            # Physically validated during the sequential loop: true positive
            status_dict[th] = 1
        else:
            # Never validated: false positive if the SNR was above the threshold, otherwise no detection
            if SNR_final >= th:
                status_dict[th] = 3
            else:
                status_dict[th] = 0

            # Failed threshold: at least all the available observations would be needed
            nobs_reached[th] = max_pts

    simult_flag = 1 if use_simultaneous else 0
    out_K_star = K_star_fit if use_simultaneous else 0.0
    out_P_star = P_star_final if use_simultaneous else 0.0

    # Lists written as "[a, b]" in the results file
    nobs_str_to_save = ", ".join([str(nobs_reached[th]) for th in thresholds])
    status_str_to_save = ", ".join([str(status_dict[th]) for th in thresholds])

    # Evolution of K, Delta BIC, chi2 and jitter with the number of points (first iterations only)
    if make_plots:
        plotting.statistics_history(**plotting.pick(locals(), plotting.statistics_history))

    boot_flag = 1 if BOOTSTRAP == True else 0

    with open(var.outputfilename, 'a') as f:
        f.write(f"{target_params['targetname']}, {target_params['v_mag']}, "
                f"{np.round(target_params['pl_mass_Mj'],4)}, "
                f"{np.round(target_params['pl_radius_Rj'],4) if target_params['pl_radius_Rj'] is not None else 0.0}, "
                f"{target_params['pl_period']}, "
                f"{target_params.get('eccentricity', 0.0)}, {target_params.get('periastron', 0.0)}, " # eccentricity and periastron
                f"{np.round(age/1e9,4)}, {tag}, "
                f"{np.round(Prot,4)}, {np.round(log_Rhk,4)}, {np.round(K_act,4)}, "
                f"{np.round(kpl,4)}, {np.round(np.std(RVfinal),4)}, {sig_ph}, "
                f"{np.round(np.std(noise_rv),4)}, "
                f"\"[{nobs_str_to_save}]\", {niter}, " 
                f"{np.round(K_fit,4)}, {np.round(K_err_final,4)}, {np.round(SNR_final,4)}, " 
                f"{np.round(delta_bic,2)}, {np.round(jitter,4)}, "                    
                f"{simult_flag}, {np.round(out_K_star,4)}, {np.round(out_P_star,4)}, "        
                f"\"[{status_str_to_save}]\", \"{var.sig_best}\", {boot_flag}\n")

    return status_dict
