"""
Periodogram analysis of the sequential recovery (paper, Sect. 2.6): GLS frequency
range, significant peaks, classification of planet and intruder peaks,
pre-whitening of the intruders and reverse pre-whitening of the planet.
"""
import numpy as np
import scipy.signal as signal
from scipy.optimize import curve_fit

from . import plotting


def frequency_range(time_f, ind, target_params):
    """Frequency range of the GLS: from 1/baseline to the Nyquist frequency of the current data
    (extended to 1.5 times the planet frequency if needed)."""
    t_sorted = np.sort(time_f[ind])
    diff = np.diff(t_sorted)
    median_diff = np.median(diff)
    nyquist = 0.5 / median_diff

    # Highest frequency (fend)
    # 1. Default: at least up to 2 c/d, to see the daily aliases
    fend = 2.0

    # 2. Dense sampling (several points per night): use the actual Nyquist frequency
    if median_diff < 0.4:
        fend = nyquist

    # NOTE: this line overrides the two choices above, so fend is always the Nyquist frequency
    # (behaviour of the version used for the paper).
    fend = nyquist

    # 3. Very short periods (USP): extend the range to include the planet
    f_planet = 1.0 / target_params['pl_period']
    if fend < (f_planet * 1.5):
        fend = f_planet * 1.5

    # Lowest frequency (fbeg): 1 / baseline
    baseline = np.max(t_sorted) - np.min(t_sorted)
    fbeg = 1.0 / baseline if baseline > 0 else 1e-5

    return fbeg, fend


def significant_peaks(per, percFAP_star, percFAP_planet, ind):
    """Peaks of the GLS above the FAP threshold for stellar signals, sorted by decreasing power."""
    freq = per.freq
    power = per.power
    # 1. Power thresholds corresponding to the two FAP levels (powerLevel of the GLS)
    threshold_star = per.powerLevel(percFAP_star)
    threshold_planet = per.powerLevel(percFAP_planet)

    # 2. Local maxima above the stellar threshold
    peaks_idx, _ = signal.find_peaks(power, height=threshold_star)

    peak_freqs = freq[peaks_idx]
    peak_powers = power[peaks_idx]
    peak_periods = 1.0 / peak_freqs

    # Sorted by decreasing power
    sort_mask = np.argsort(peak_powers)[::-1]
    sorted_periods = peak_periods[sort_mask]
    sorted_powers = peak_powers[sort_mask]

    # 3. Log
    print(f"\n--- PEAK ANALYSIS (Zechmeister GLS, FAP < {percFAP_star*100}%) with {len(ind)} points ---")
    print(f"Power threshold for FAP {percFAP_star*100}%: {threshold_star:.4f}")

    if len(sorted_periods) > 0:
        print(f"Found {len(sorted_periods)} peaks above the threshold:")
        for i in range(min(5, len(sorted_periods))):
            this_fap = per.FAP(sorted_powers[i])
            print(f"  {i+1}. Period: {sorted_periods[i]:.4f} d (Power: {sorted_powers[i]:.4f}, FAP: {this_fap:.4f})")
    else:
        print(f"  No peak above the {percFAP_star*100}% FAP.")

    return freq, power, threshold_star, threshold_planet, sorted_periods, sorted_powers


def classify_peaks(t_sel, target_params, sorted_periods):
    """Classify the significant peaks as planet candidates (within the frequency resolution of the
    known planet frequency) or intruders (e.g. stellar rotation and its harmonics)."""
    baseline = np.max(t_sel) - np.min(t_sel)
    freq_res = 1.0 / baseline
    f_planet_known = 1.0 / target_params['pl_period']

    # Indices referring to sorted_periods
    idx_planet_candidates = []
    idx_intruders = []

    peak_uncert = 1.5  # planet window: 1.5 times the frequency resolution
    for i, P_val in enumerate(sorted_periods):
        f_val = 1.0 / P_val
        dist = abs(f_val - f_planet_known)

        if dist < (peak_uncert * freq_res):
            # The planet
            idx_planet_candidates.append(i)
        else:
            # Daily alias of the planet (|f - 1 c/d| ~ f_planet)? Slightly wider tolerance,
            # because of the imperfect sampling
            is_alias_minus = abs((1.0027 - f_val) - f_planet_known) < (2.0 * freq_res)
            is_alias_plus  = abs((f_val - 1.0027) - f_planet_known) < (2.0 * freq_res)

            if is_alias_minus or is_alias_plus:
                print(f"   [!] Ignoring the peak at {P_val:.3f}d: it is an alias of the planet.")
            else:
                idx_intruders.append(i)  # a real intruder (possibly the star)

    # Is the planet the highest peak of the periodogram?
    planet_is_highest = (len(idx_planet_candidates) > 0 and idx_planet_candidates[0] == 0)
    if planet_is_highest:
        print("   > The planet is the dominant peak.")

    return freq_res, f_planet_known, idx_planet_candidates, idx_intruders, peak_uncert, planet_is_highest


def prewhitening(per, t_sel, y_sel, err_sel, fbeg, fend, sorted_periods, sorted_powers, idx_intruders,
                 idx_planet_candidates, f_planet_known, freq_res, peak_uncert, percFAP_star, percFAP_planet,
                 threshold_star, threshold_planet, P_intruder, DEBUG_PLOT, make_plots, var, niter, debug_pdf, Prot):
    """Iterative pre-whitening of the intruder peaks, in decreasing order of power.

    Each intruder is fitted and subtracted; it is accepted as the stellar signal if the planet
    peak improves in the residuals (this triggers the simultaneous fit), rejected otherwise.
    P_intruder is passed in and returned so that its last value survives across steps, as
    in the version used for the paper."""
    best_mode = "RAW"
    use_simultaneous = False
    P_star_final = None

    # A. Planet peak in the raw data, before touching them
    search_radius = peak_uncert * freq_res
    f_min = f_planet_known - search_radius
    f_max = f_planet_known + search_radius

    mask_window = (per.freq >= f_min) & (per.freq <= f_max)
    if np.any(mask_window):
        idx_local_max = np.argmax(per.power[mask_window])
        power_planet_raw = per.power[mask_window][idx_local_max]
        freq_planet_raw  = per.freq[mask_window][idx_local_max]
    else:
        idx_closest = np.argmin(np.abs(per.freq - f_planet_known))
        power_planet_raw = per.power[idx_closest]
        freq_planet_raw  = per.freq[idx_closest]

    dist_raw = abs(freq_planet_raw - f_planet_known)

    # One attempt per intruder, in decreasing order of power
    max_attempts = len(idx_intruders)
    intruder_accepted = False

    if max_attempts > 0:
        print(f"   > Pre-whitening test (max {max_attempts} attempts). Raw planet power: {power_planet_raw:.4f}")

    for attempt in range(max_attempts):
        curr_intruder_idx = idx_intruders[attempt]
        P_intruder = sorted_periods[curr_intruder_idx]
        Pow_intruder = sorted_powers[curr_intruder_idx]

        print(f"\n   --- ATTEMPT {attempt+1}: intruder at P={P_intruder:.4f}d (Pow={Pow_intruder:.4f}) ---")

        # B. Fit and subtract the intruder
        try:
            omega_int = 2 * np.pi / P_intruder
            def model_intruder(t, gamma, K, ph):
                return gamma + K * np.sin(omega_int * t + ph)

            guess_K = np.std(y_sel) * 1.4
            p0 = [np.mean(y_sel), guess_K, 0.0]
            popt_i, _ = curve_fit(model_intruder, t_sel, y_sel, sigma=err_sel, p0=p0, maxfev=2000)

            rv_trial_residuals = y_sel - model_intruder(t_sel, *popt_i)

            # C. GLS of the trial residuals
            per_trial = var.Gls((t_sel, rv_trial_residuals, err_sel), fbeg=fbeg, fend=fend, verbose=False)
            new_threshold_planet = per_trial.powerLevel(percFAP_planet)

            # D. Planet peak in the residuals
            mask_res = (per_trial.freq >= f_min) & (per_trial.freq <= f_max)
            if np.any(mask_res):
                idx_local_max_res = np.argmax(per_trial.power[mask_res])
                power_planet_res = per_trial.power[mask_res][idx_local_max_res]
                freq_planet_res  = per_trial.freq[mask_res][idx_local_max_res]
            else:
                idx_closest_res = np.argmin(np.abs(per_trial.freq - f_planet_known))
                power_planet_res = per_trial.power[idx_closest_res]
                freq_planet_res  = per_trial.freq[idx_closest_res]

            dist_res = abs(freq_planet_res - f_planet_known)
            print(f"       Planet: Pow {power_planet_raw:.4f} -> {power_planet_res:.4f}")

            # E. Did removing the intruder destroy the planet?
            #    - planet power drops by more than 30%, or
            #    - the planet falls below the FAP threshold
            cond_power_drop = power_planet_res < (0.70 * power_planet_raw)
            cond_lost = power_planet_res < new_threshold_planet

            # A weak planet cannot produce dominant aliases: a very strong intruder is the stellar
            # activity, and the drop of the planet power only means that the spectral leakage
            # has been cleaned ("blind override").
            is_dominant_activity_peak = (Pow_intruder > 0.40)
            if cond_power_drop or cond_lost:
                if is_dominant_activity_peak:
                    print(f"       🎯 [Blind Override] The peak at {P_intruder:.2f}d is very strong (Pow={Pow_intruder:.4f}).")
                    print(f"       🎯 Not an alias of the planet: it is the dominant activity. The planet drop is only spectral cleaning. Accepting it.")
                    cond_power_drop = False
                    cond_lost = False
                else:
                    print(f"       ❌ WARNING: this minor peak was an alias/harmonic: removing it destroyed the planet.")
                    print(f"       Rejecting this intruder and trying the next one (if any).")
                    continue  # rejected: next intruder

            # The planet survived: did it improve?
            cond1 = (power_planet_res > power_planet_raw)
            cond2 = (power_planet_res > 0.9 * power_planet_raw) and (dist_res < 0.5 * dist_raw)

            if cond1 or cond2 or is_dominant_activity_peak:
                print(f"       ✅ IMPROVEMENT CONFIRMED: using the residuals of intruder {attempt+1}.")
                best_mode = "RESIDUALS"
                P_star_final = P_intruder
                use_simultaneous = True
                intruder_accepted = True
                is_improved = (cond1 or cond2)
                color_res = 'green'
                status_msg = "SUCCESS: Better Power/Accuracy"

            else:
                print(f"   >>> NO BENEFIT. Keeping the raw data.")
                best_mode = "RAW"
                color_res = 'red'
                status_msg = "FAILED: Worsened Signal"
                is_improved = False

            # Debug figure of this attempt
            if DEBUG_PLOT and make_plots:
                plotting.prewhitening_attempt(**plotting.pick(locals(), plotting.prewhitening_attempt))

            if intruder_accepted:
                break  # success: stop the attempts

            else:
                print(f"       ➖ No net benefit (but no destruction). Trying the next intruder.")

        except Exception as e:
            print(f"       ⚠️ Fit of intruder {attempt+1} failed: {e}")

    return intruder_accepted, use_simultaneous, P_star_final, P_intruder


def reverse_prewhitening(per, t_sel, y_sel, err_sel, fbeg, fend, target_params, planet_is_highest, percFAP_star, freq_res,
                         P_star_final, use_simultaneous, DEBUG_PLOT, make_plots, var, niter, debug_pdf, Prot):
    """Reverse pre-whitening: subtract a preliminary planet model and look for a significant stellar
    peak in the residuals; if found, it triggers the simultaneous fit."""
    if planet_is_highest:
        print("   >>> No intruder accepted, but the planet is the dominant peak.")
        print("   >>> REVERSE PRE-WHITENING: subtracting the planet to look for the star...")

        # 1. Preliminary planet-only model (fixed phase)
        T0_temp = target_params['transit_midpoint']
        P_temp = target_params['pl_period']

        def temp_planet_model(t, gamma, K):
            omega = 2 * np.pi / P_temp
            return gamma - K * np.sin(omega * (t - T0_temp))

        try:
            p0_temp = [np.median(y_sel), np.std(y_sel)]
            popt_temp, _ = curve_fit(temp_planet_model, t_sel, y_sel, sigma=err_sel, p0=p0_temp, maxfev=2000)

            # 2. Residuals: raw data - planet model
            rv_no_planet = y_sel - temp_planet_model(t_sel, *popt_temp)

            # 3. GLS of the residuals and its highest peak
            per_rev = var.Gls((t_sel, rv_no_planet, err_sel), fbeg=fbeg, fend=fend, verbose=False)
            thresh_star_rev = per_rev.powerLevel(percFAP_star)

            idx_max_rev = np.argmax(per_rev.power)
            pow_max_rev = per_rev.power[idx_max_rev]
            f_star_rev = per_rev.freq[idx_max_rev]
            P_star_rev = 1.0 / f_star_rev

            # 4. Does the star emerge from the noise?
            if pow_max_rev > thresh_star_rev:
                print(f"       ✅ SUCCESS: star found in the residuals: P_star={P_star_rev:.4f}d (Pow={pow_max_rev:.4f} > FAP threshold)")
                use_simultaneous = True
                P_star_final = P_star_rev
                rev_status = "SUCCESS: Star Detected in Residuals"
                rev_color = "green"

            else:
                print(f"       ❌ FAILED: no significant stellar signal (Pow={pow_max_rev:.4f} < FAP threshold).")
                print("       Planet-only fit on the raw data.")
                use_simultaneous = False
                rev_status = "FAILED: No Star Found in Residuals"
                rev_color = "red"

            # Debug figure of the reverse pre-whitening
            if DEBUG_PLOT and make_plots:
                plotting.reverse_prewhitening(**plotting.pick(locals(), plotting.reverse_prewhitening))

        except Exception as e:
            print(f"       ⚠️ Reverse pre-whitening fit failed: {e}. Planet-only fit on the raw data.")
            use_simultaneous = False

    else:
        print("   >>> No intruder accepted and the planet is not the dominant peak: keeping the raw data.")
        use_simultaneous = False

    return P_star_final, use_simultaneous
