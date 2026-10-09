"""Statistics of the required number of RV observations over the Monte Carlo realisations."""
import sys

import numpy as np


def print_execution_summary(target_name, log_rhk, nobs_history_dict, status_flag_dict, output_filename, thresholds, is_bootstrapped=False):
    """Median and 16th/84th percentiles of the number of RV observations needed, for each SNR
    threshold: over all simulations (failures saturated to the total number of observations)
    and over the successful ones only, plus the number of failures.

    Printed to the log and to the screen, and appended to the per-target results file.
    Returns the statistics and the failures per threshold (used for the catalog summary).
    """
    # Only the simulations actually run (the skipped ones are NaN)
    n_iter = int(np.sum(~np.isnan(np.asarray(status_flag_dict[thresholds[0]], dtype=float))))

    boot_text = "ON" if is_bootstrapped else "OFF"

    summary_msg = (
        "\n" + "*"*80 + "\n" +
        f" SUMMARY: {target_name}\n" +
        "*"*80 + "\n" +
        f" logrhk        = {log_rhk}\n"
        f" Tot. Runs     = {n_iter}\n"
        f" --- RV POINTS DISTRIBUTION (ALL vs SUCCESS ONLY) & FAILURES --- [BOOTSTRAP: {boot_text}]\n"
    )

    results_stats = {}
    fails_stats = {}
    fails_detail = {}   # (total, no detection, false positive) per threshold

    for th in thresholds:
        numeros = np.asarray(nobs_history_dict[th], dtype=float)
        flags = np.asarray(status_flag_dict[th], dtype=float)
        valid = ~np.isnan(flags)
        numeros, flags = numeros[valid], flags[valid]

        # 1. All simulations (failures are saturated to the total number of observations)
        if len(numeros) > 0:
            med_all = np.median(numeros)
            p16_all = np.percentile(numeros, 16)
            p84_all = np.percentile(numeros, 84)
            fail_frac = np.mean(flags != 1)
        else:
            med_all = p16_all = p84_all = fail_frac = np.nan

        # If >= 50% of the simulations fail, the median equals N_total: it is a lower limit
        lim_med = ">=" if fail_frac >= 0.5 else ""
        lim_84  = ">=" if fail_frac >= 0.16 else ""

        # 2. Successful simulations only (flag == 1)
        mask_success = (flags == 1)
        if np.any(mask_success):
            med_succ = np.median(numeros[mask_success])
            p16_succ = np.percentile(numeros[mask_success], 16)
            p84_succ = np.percentile(numeros[mask_success], 84)
        else:
            med_succ, p16_succ, p84_succ = np.nan, np.nan, np.nan  # no success

        # Failures for this threshold: no detection (flag 0) and false positives (flag 3)
        fails_total = len(np.where((flags == 0) | (flags == 3))[0])
        fails_nd = len(np.where(flags == 0)[0])
        fails_fp = len(np.where(flags == 3)[0])

        results_stats[th] = (med_all, p16_all, p84_all, med_succ, p16_succ, p84_succ)
        fails_stats[th] = fails_total
        fails_detail[th] = (fails_total, fails_nd, fails_fp)

        summary_msg += f" SNR >= {th:<4} | ALL: {lim_med}{med_all:.1f} [{p16_all:.0f}-{lim_84}{p84_all:.0f}] | SUCC: {med_succ:.1f} [{p16_succ:.0f}-{p84_succ:.0f}] | Fails: {fails_total} (ND:{fails_nd}, FP:{fails_fp})\n"

    summary_msg += ("*"*80 + "\n")

    # 1. Log and screen
    print(summary_msg)
    sys.__stdout__.write(summary_msg)

    # 2. Per-target results file
    with open(output_filename, 'a') as f:
        f.write('*' * 180 + ' \n')
        file_str = f"{target_name} "
        for th in thresholds:
            m_all, p16_all, p84_all, m_succ, p16_succ, p84_succ = results_stats[th]
            fails_total, fails_nd, fails_fp = fails_detail[th]

            file_str += f"Med_All(SNR_{th})={m_all:.1f} [{p16_all:.0f}-{p84_all:.0f}] Med_Succ(SNR_{th})={m_succ:.1f} [{p16_succ:.0f}-{p84_succ:.0f}] Fails: {fails_total} (ND:{fails_nd}, FP:{fails_fp}) "

        file_str += f"Niter = {n_iter}\n"
        f.write(file_str)

    return results_stats, fails_stats
