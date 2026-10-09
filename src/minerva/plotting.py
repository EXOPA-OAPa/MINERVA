"""
Debug figures of the sequential analysis (enabled with run_options in the config).
They are pure output: they do not change any result.
"""
import inspect

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np

from ._vendor.gls import Gls
from .fitting import model_planet_only
from .tracks import fap_levels



def pick(local_vars, func):
    """Arguments for a debug-figure function, taken from the caller's local variables.

    Only the variables that exist at the time of the call are passed: some of them
    (e.g. the seasonal parameters) exist only for simultaneous fits, and the figure
    code reads them only in that case. Missing ones default to None.
    """
    wanted = inspect.signature(func).parameters
    return {k: v for k, v in local_vars.items() if k in wanted}

def gls_peaks(debug_pdf=None, freq=None, ind=None, percFAP_planet=None, percFAP_star=None, power=None, Prot=None, sorted_periods=None, sorted_powers=None, target_params=None, threshold_planet=None, threshold_star=None, var=None):
    """GLS periodogram of the current data set with the FAP thresholds and the classified peaks."""
    fig = plt.figure(figsize=(10, 5))
    plt.plot(1/freq, power, 'k-', lw=1, label='GLS (Zechmeister)')
    plt.axhline(threshold_star, color='red', linestyle='--', label=f'FAP {percFAP_star*100}%')
    plt.axhline(threshold_planet, color='red', linestyle='-.', label=f'FAP {percFAP_planet*100}%')

    if target_params['Prot'] is not None: plt.axvline(target_params['Prot'], ls='-', lw=2, alpha=0.5, color='red', label='Prot')
    else: plt.axvline(Prot, ls='-', lw=2, alpha=0.5, color='red', label='Prot')

    plt.axvline(target_params['pl_period'], ls='-', lw=2, alpha=0.5, color='blue', label='Pp')  

    if len(sorted_periods) > 0:
        plt.plot(sorted_periods, sorted_powers, 'rx', markersize=8)

    plt.xscale('log')
    plt.xlabel("Period [days]")
    plt.ylabel("Normalised power") # normalised between 0 and 1
    plt.title(f"GLS periodogram (N={len(ind)})")
    plt.legend()
    # Save to the PDF and/or show on screen

    # 1. PDF: only for the first simulations (niter < 3)
    if debug_pdf is not None:
        debug_pdf.savefig(fig)

    # 2. Screen: only if show_on_screen is set in run_options
    if getattr(var, 'SHOW_ON_SCREEN', False): 
        plt.show()

    # 3. Always close the figure, to avoid running out of memory
    plt.close(fig)  


def prewhitening_attempt(best_mode=None, color_res=None, debug_pdf=None, dist_raw=None, dist_res=None, err_sel=None, f_max=None, f_min=None, f_planet_known=None, freq_planet_raw=None, freq_planet_res=None, freq_res=None, idx_intruders=None, idx_planet_candidates=None, is_improved=None, model_intruder=None, niter=None, P_intruder=None, per=None, per_trial=None, percFAP_planet=None, percFAP_star=None, popt_i=None, Prot=None, sorted_periods=None, sorted_powers=None, status_msg=None, t_sel=None, threshold_planet=None, threshold_star=None, var=None, y_sel=None):
    """Diagnostic figure of one pre-whitening attempt (raw vs residual periodogram)."""
    f_rot_known = 1.0 / Prot
    fig, ax = plt.subplots(3, 1, figsize=(10, 14), constrained_layout=False)

    # --- 1. Time domain ---
    t_mod = np.linspace(np.min(t_sel), np.max(t_sel), 500)
    y_mod = model_intruder(t_mod, *popt_i)

    ax[0].errorbar(t_sel, y_sel, yerr=err_sel, fmt='o', color='black', alpha=0.6, label='Raw Data')
    ax[0].plot(t_mod, y_mod, color='orange', lw=2, label=f'Intruder Fit (P={P_intruder:.2f}d)')
    ax[0].set_ylabel("RV [m/s]")
    ax[0].set_title(f"1. Time Domain: Modeling the Intruder ({best_mode} selected)", fontweight='bold')
    ax[0].legend(loc='upper right')
    ax[0].grid(True, alpha=0.3)

    # --- 2. Frequency domain: overview ---
    # GLS of the raw and cleaned data
    ax[1].plot(per.freq, per.power, color='black', alpha=0.5, label='Raw GLS')
    ax[1].plot(per_trial.freq, per_trial.power, color=color_res, lw=2, label='Cleaned GLS')

    # --- FAP level and classified peaks ---
    # 1. FAP level
    ax[1].axhline(threshold_star, color='gray', linestyle='-.', lw=1.5, label=f'FAP {percFAP_star*100}%')
    ax[1].axhline(threshold_planet, color='gray', linestyle='--', lw=1.5, label=f'FAP {percFAP_planet*100}%')

    # Vertical lines: planet and intruder
    ax[1].axvline(f_planet_known, color='blue', linestyle='--', alpha=0.7, label='Planet Freq')
    ax[1].axvline(f_rot_known, color='red', linestyle='--', alpha=0.7, label='Star Freq')
    ax[1].axvline(1.0/P_intruder, color='orange', linestyle='--', alpha=0.7, label='Intruder Freq')

    # 2. Classified peaks
    # Planet candidates (green)
    if len(idx_planet_candidates) > 0:
        f_cand = 1.0 / sorted_periods[idx_planet_candidates]
        p_cand = sorted_powers[idx_planet_candidates]
        ax[1].scatter(f_cand, p_cand, s=80, facecolors='none', edgecolors='green', lw=2, label='Planet Candidates', zorder=5)

    # Intruders (red)
    if len(idx_intruders) > 0:
        f_int = 1.0 / sorted_periods[idx_intruders]
        p_int = sorted_powers[idx_intruders]
        ax[1].scatter(f_int, p_int, s=80, facecolors='none', edgecolors='red', lw=2, label='Intruders', zorder=5)

    # 3. The removed intruder (filled star)
    target_f = 1.0 / P_intruder
    # its power (approximate, for the figure)
    target_p_idx = np.argmin(np.abs(per.freq - target_f))
    target_p = per.power[target_p_idx]
    ax[1].scatter([target_f], [target_p], s=150, marker='*', color='red', label='REMOVED TARGET', zorder=10)


    ax[1].set_ylabel("Power")
    ax[1].set_title(f"2. Frequency Domain: Peak Classification (Threshold={threshold_star:.3f})")
    ax[1].legend(loc='upper right', fontsize='small', ncol=2)
    ax[1].grid(True, alpha=0.3)

    # --- 3. Frequency domain: zoom on the planet ---
    zoom_window = 8 * freq_res
    ax[2].set_xlim(f_planet_known - zoom_window, f_planet_known + zoom_window)

    ax[2].plot(per.freq, per.power, color='black', alpha=0.5, label='Raw')
    ax[2].plot(per_trial.freq, per_trial.power, color=color_res, lw=2, label='Cleaned')

    # FAP level, if inside the zoom
    ax[2].axhline(threshold_star, color='gray', linestyle='-.', lw=1, alpha=0.5, label=f'FAP {percFAP_star*100}%')
    ax[2].axhline(threshold_planet, color='gray', linestyle='--', lw=1, alpha=0.5, label=f'FAP {percFAP_planet*100}%')

    # Expected planet frequency
    ax[2].axvline(f_planet_known, color='blue', linestyle='--', lw=2, label='True Freq')

    # Planet peaks found (raw and cleaned)
    ax[2].axvline(freq_planet_raw, color='black', linestyle=':', alpha=0.5)
    ax[2].axvline(freq_planet_res, color=color_res, linestyle=':', alpha=0.8)

    # Shift of the planet peak
    ax[2].annotate("", xy=(f_planet_known, 0), xytext=(freq_planet_raw, 0),
                   arrowprops=dict(arrowstyle="->", color='black', alpha=0.5))
    if is_improved: # o cond_tradeoff
        ax[2].annotate("Shift", xy=(f_planet_known, 0), xytext=(freq_planet_res, 0),
                       arrowprops=dict(arrowstyle="->", color='green', lw=2))

    ax[2].set_title(f"3. Verdict: {status_msg}\nRaw Dist: {dist_raw:.4f} -> Res Dist: {dist_res:.4f}", 
                    color=color_res, fontweight='bold')

    ax[2].axvspan(f_min, f_max, color='gold', alpha=0.1, label='Target Window')
    ax[2].legend()
    ax[2].grid(True, alpha=0.3)

    # Save to the PDF and/or show on screen

    # 1. PDF: only for the first simulations (niter < 3)
    if debug_pdf is not None and niter < 3:
        debug_pdf.savefig(fig)

    # 2. Screen: only if show_on_screen is set in run_options
    if getattr(var, 'SHOW_ON_SCREEN', False): 
        plt.show()

    # 3. Always close the figure, to avoid running out of memory
    plt.close(fig)


def reverse_prewhitening(debug_pdf=None, err_sel=None, f_star_rev=None, freq_res=None, niter=None, P_star_rev=None, P_temp=None, per=None, per_rev=None, percFAP_star=None, popt_temp=None, pow_max_rev=None, Prot=None, rev_color=None, rev_status=None, t_sel=None, temp_planet_model=None, thresh_star_rev=None, var=None, y_sel=None):
    """Diagnostic figure of the reverse pre-whitening (periodogram of the planet-subtracted data)."""
    fig_rev, ax_rev = plt.subplots(3, 1, figsize=(10, 14), constrained_layout=True)

    # --- 1. Time domain ---
    t_mod_rev = np.linspace(np.min(t_sel), np.max(t_sel), 500)
    y_mod_rev = temp_planet_model(t_mod_rev, *popt_temp)

    ax_rev[0].errorbar(t_sel, y_sel, yerr=err_sel, fmt='o', color='black', alpha=0.6, label='Raw Data')
    ax_rev[0].plot(t_mod_rev, y_mod_rev, color='blue', lw=2, label=f'Temp Planet Fit (P={P_temp:.2f}d)')
    ax_rev[0].set_ylabel("RV [m/s]")
    ax_rev[0].set_title("1. Time Domain: Modeling the Dominant Planet", fontweight='bold')
    ax_rev[0].legend(loc='upper right')
    ax_rev[0].grid(True, alpha=0.3)

    # --- 2. Frequency domain: raw vs residuals ---
    ax_rev[1].plot(per.freq, per.power, color='black', alpha=0.3, label='Raw GLS (Planet Dominates)')
    ax_rev[1].plot(per_rev.freq, per_rev.power, color=rev_color, lw=2, label='Residuals GLS (Planet Removed)')
    ax_rev[1].axhline(thresh_star_rev, color='gray', linestyle='-.', lw=1.5, label=f'FAP {percFAP_star*100}%')

    # Reference lines
    ax_rev[1].axvline(1.0/P_temp, color='blue', linestyle='--', alpha=0.7, label='Planet Freq')

    # True stellar frequency (if known)
    f_rot_known_plot = 1.0 / Prot
    ax_rev[1].axvline(f_rot_known_plot, color='red', linestyle='--', alpha=0.7, label='True Star Freq')

    # Stellar peak found in the residuals
    ax_rev[1].scatter([f_star_rev], [pow_max_rev], s=150, marker='*', color=rev_color, label='Max Residual Peak', zorder=10)

    ax_rev[1].set_ylabel("Power")
    ax_rev[1].set_title(f"2. Frequency Domain: Reverse Pre-Whitening (Max Pow={pow_max_rev:.3f})")
    ax_rev[1].legend(loc='upper right', fontsize='small', ncol=2)
    ax_rev[1].grid(True, alpha=0.3)

    # --- 3. Frequency domain: zoom on the peak found ---
    zoom_win_rev = 15 * freq_res
    ax_rev[2].set_xlim(f_star_rev - zoom_win_rev, f_star_rev + zoom_win_rev)

    ax_rev[2].plot(per_rev.freq, per_rev.power, color=rev_color, lw=2, label='Residuals')
    ax_rev[2].axhline(thresh_star_rev, color='gray', linestyle='-.', lw=1, alpha=0.5, label=f'FAP {percFAP_star*100}%')

    ax_rev[2].axvline(f_rot_known_plot, color='red', linestyle='--', lw=2, label='True Star Freq')

    ax_rev[2].axvline(f_star_rev, color=rev_color, linestyle=':', lw=2, label=f'Found Peak (P={P_star_rev:.2f}d)')

    ax_rev[2].set_title(f"3. Verdict: {rev_status}", color=rev_color, fontweight='bold')
    ax_rev[2].legend(loc='upper right')
    ax_rev[2].grid(True, alpha=0.3)

    # Save to the PDF and/or show on screen
    if debug_pdf is not None and niter < 3:
        debug_pdf.savefig(fig_rev, bbox_inches='tight')

    if getattr(var, 'SHOW_ON_SCREEN', False): 
        plt.show()

    plt.close(fig_rev)


def fit_diagnostics(debug_pdf=None, dof=None, err_fit=None, K_err_formal=None, K_fit=None, model_simultaneous_seasonal=None, P_PLANET=None, popt=None, rv_fit=None, T0_PLANET=None, t_fit=None, use_simultaneous=None, var=None, y_model_final=None):
    """Diagnostic figure of the final fit: data, model and residuals."""
    fig_diag, ax_diag = plt.subplots(figsize=(10, 6))
    t_dense = np.linspace(np.min(t_fit), np.max(t_fit), 1000)

    # Model on a dense time grid
    if use_simultaneous:
        y_dense = model_simultaneous_seasonal(t_dense, *popt)
    else:
        y_dense = model_planet_only(t_dense, *popt, P_PLANET, T0_PLANET)

    # Residuals and reduced chi2
    residuals = rv_fit - y_model_final
    chisq = np.sum((residuals / err_fit)**2)
    chisq_red = chisq / dof if dof > 0 else 1.0

    # Data and fit
    ax_diag.errorbar(t_fit, rv_fit, yerr=err_fit, fmt='o', color='black', ecolor='blue', alpha=0.7, label='RV data')
    ax_diag.plot(t_dense, y_dense, color='blue', lw=2, label='Model fit')

    # Fit statistics
    text_info = f"$\\chi^2_{{\\text{{red}}}}$: {chisq_red:.2f}\nSNR K_pl: {K_fit/K_err_formal:.2f}"
    ax_diag.text(0.02, 0.05, text_info, transform=ax_diag.transAxes, 
                 bbox=dict(facecolor='white', alpha=0.8, edgecolor='blue'))

    # Axes and labels
    ax_diag.set_title("RV model fit", fontsize=11, fontweight='bold')
    ax_diag.set_xlabel("Time [d]")
    ax_diag.set_ylabel("RV [m/s]")
    ax_diag.legend(loc='upper right')
    ax_diag.grid(True, alpha=0.3)

    # Save and/or show
    if debug_pdf is not None: 
        debug_pdf.savefig(fig_diag, bbox_inches='tight')
    if getattr(var, 'SHOW_ON_SCREEN', False): 
        plt.show()
    plt.close(fig_diag)


def phase_folded_fit(chisq_red=None, current_seasons_idx=None, debug_pdf=None, delta_bic=None, err_fit=None, fit_desc=None, idx_split=None, K_err_final=None, K_fit=None, n_new_points=None, N_seasons=None, niter=None, P_PLANET=None, P_star_final=None, popt=None, residuals=None, rv_fit=None, rv_plot_data=None, scale_factor=None, SNR_final=None, T0_PLANET=None, t_fit=None, title_color=None, use_simultaneous=None, var=None):
    """Phase-folded planet and stellar-activity components of the final fit."""

    # 3 panels if the star is fitted too, 2 otherwise
    n_panels = 3 if use_simultaneous else 2
    fig = plt.figure(figsize=(10, 4 * n_panels))

    if use_simultaneous:
        gs = gridspec.GridSpec(3, 1, height_ratios=[3, 3, 1.5], hspace=0.3)
        ax_phase = plt.subplot(gs[0])
        ax_star  = plt.subplot(gs[1])
        ax_res   = plt.subplot(gs[2])
    else:
        gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1.5], hspace=0.25)
        ax_phase = plt.subplot(gs[0])
        ax_res   = plt.subplot(gs[1])

    # Models and phases
    phase_dense = np.linspace(0, 1, 200)

    # Offset (gamma) of the planet model
    if use_simultaneous:
        gamma_plot = 0.0       # the planet is centred on 0 (the offsets are in the stellar model)
    else:
        gamma_plot = popt[0]   # planet-only model: popt[0] is gamma

    # Planet
    t_dense_pl = phase_dense * P_PLANET + T0_PLANET
    y_plot_dense_pl = model_planet_only(t_dense_pl, gamma_plot, K_fit, P_PLANET, T0_PLANET)
    phase_pl = ((t_fit - T0_PLANET) / P_PLANET) % 1.0

    label_data_pl = 'Old Data' 
    label_new_pl  = f'New Data (+{n_new_points})'

    # --- Planet curve ---
    ax_phase.plot(phase_dense, y_plot_dense_pl, color=title_color, lw=3, 
                  label=f'Planet Fit (K={K_fit:.2f} m/s)')

    # Planet uncertainty band
    osc_pl = y_plot_dense_pl - gamma_plot
    rel_err = min(K_err_final / K_fit, 1.0) if K_fit > 0 else 1.0
    ax_phase.fill_between(phase_dense, gamma_plot + osc_pl*(1-rel_err), gamma_plot + osc_pl*(1+rel_err), 
                          color=title_color, alpha=0.2, label=f'Error (x{scale_factor:.1f})')

    ax_res.axhline(0, color='gray', linestyle='--')

    # Data points (simultaneous vs planet-only fit)
    if use_simultaneous:
        # Stellar component of the data
        rv_pl_only_model = -K_fit * np.sin(2*np.pi/P_PLANET * (t_fit - T0_PLANET))
        rv_star_data = rv_fit - rv_pl_only_model
        phase_st = (t_fit % P_star_final) / P_star_final

        omega_star = 2 * np.pi / P_star_final
        t_dense_st = phase_dense * P_star_final 

        import matplotlib.cm as cm
        colors = [cm.tab10(i % 10) for i in range(N_seasons)]

        # One colour per season, on all three panels
        for s, loc_idx in enumerate(current_seasons_idx):
            gamma_s = popt[1 + 3*s]
            A_s     = popt[2 + 3*s]
            B_s     = popt[3 + 3*s]
            K_local = np.sqrt(A_s**2 + B_s**2)
            phi_local = np.arctan2(B_s, A_s)


            # 1. Stellar curve
            y_plot_dense_st = gamma_s + A_s * np.sin(omega_star * t_dense_st) + B_s * np.cos(omega_star * t_dense_st)
            # Legend with the seasonal parameters
            label_season = f'S{s+1} (K={K_local:.1f}, $\\gamma$={gamma_s:.1f}, $\\phi$={phi_local:.2f} rad)'

            ax_star.plot(phase_dense, y_plot_dense_st, color=colors[s], lw=2.5, alpha=0.7, 
                         label=label_season)

            # 2. Old and new points of this season
            idx_old = loc_idx[loc_idx < idx_split]
            idx_new = loc_idx[loc_idx >= idx_split]

            # --- Old points (season colour) ---
            if len(idx_old) > 0:
                # Planet
                ax_phase.errorbar(phase_pl[idx_old], rv_plot_data[idx_old], yerr=err_fit[idx_old], 
                                  fmt='o', color=colors[s], markeredgecolor='black', markeredgewidth=0.5, alpha=0.85, zorder=5)
                # Star
                ax_star.errorbar(phase_st[idx_old], rv_star_data[idx_old], yerr=err_fit[idx_old], 
                                 fmt='o', color=colors[s], markeredgecolor='black', markeredgewidth=0.5, alpha=0.85, zorder=5)
                # Residuals
                ax_res.errorbar(t_fit[idx_old], residuals[idx_old], yerr=err_fit[idx_old], 
                                fmt='o', color=colors[s], markeredgecolor='black', markeredgewidth=0.5, alpha=0.85, zorder=5)

            # --- New points (red crosses) ---
            if len(idx_new) > 0:
                # Planet
                ax_phase.errorbar(phase_pl[idx_new], rv_plot_data[idx_new], yerr=err_fit[idx_new], 
                                  fmt='X', color='red', markeredgecolor='black', markersize=9, alpha=1.0, zorder=10)
                # Star
                ax_star.errorbar(phase_st[idx_new], rv_star_data[idx_new], yerr=err_fit[idx_new], 
                                 fmt='X', color='red', markeredgecolor='black', markersize=9, alpha=1.0, zorder=10)
                # Residuals
                ax_res.errorbar(t_fit[idx_new], residuals[idx_new], yerr=err_fit[idx_new], 
                                fmt='X', color='red', markeredgecolor='black', markersize=9, alpha=1.0, zorder=10)

        # Stellar panel
        ax_star.set_title(f"Stellar Activity Phase Fold ({N_seasons} seasons, P={P_star_final:.2f}d)", color='orange', fontweight='bold', fontsize=11)
        ax_star.set_xlabel("Star Phase")
        ax_star.set_ylabel("RV [m/s]")
        # Smaller font, 2 columns, to fit the legend
        legend_fontsize = 'x-small' if N_seasons > 3 else 'small'
        ax_star.legend(loc='upper right', fontsize=legend_fontsize, ncol=2)
        ax_star.grid(True, alpha=0.3)

    else:
        # Planet-only fit (no seasons)
        # old points black, new points red
        ax_phase.errorbar(phase_pl[:idx_split], rv_plot_data[:idx_split], yerr=err_fit[:idx_split], 
                          fmt='o', color='black', alpha=0.6, zorder=5, label=label_data_pl)
        ax_res.errorbar(t_fit[:idx_split], residuals[:idx_split], yerr=err_fit[:idx_split], 
                        fmt='o', color='black', alpha=0.5)

        if n_new_points > 0:
            ax_phase.errorbar(phase_pl[idx_split:], rv_plot_data[idx_split:], yerr=err_fit[idx_split:], 
                              fmt='X', color='red', markeredgecolor='black', markersize=9, alpha=1.0, zorder=10, label=label_new_pl)
            ax_res.errorbar(t_fit[idx_split:], residuals[idx_split:], yerr=err_fit[idx_split:], 
                            fmt='X', color='red', markeredgecolor='black', markersize=9, alpha=1.0, zorder=10)

    # Titles
    #ax_phase.set_title(f"Fit Mode: {fit_desc}\nSNR={SNR_final:.2f}, $\Delta$BIC={delta_bic:.1f} (Chi2_red={chisq_red:.2f})", 
    ax_phase.set_title(f"Fit Mode: {fit_desc}\nSNR={SNR_final:.2f}, $\\Delta$BIC={delta_bic:.1f} (Chi2_red={chisq_red:.2f})", 
                        color=title_color, fontweight='bold', fontsize=12)
    ax_phase.set_xlabel("Planet Phase (0 = Transit T0)")
    ax_phase.set_ylabel("RV [m/s]")
    ax_phase.legend(loc='upper right')
    ax_phase.grid(True, alpha=0.3)

    ax_res.set_xlabel("Time [d]")
    ax_res.set_ylabel("Resid. [m/s]")
    ax_res.grid(True, alpha=0.3)

    # Save to the PDF and/or show on screen
    if debug_pdf is not None and niter < 3:
        debug_pdf.savefig(fig, bbox_inches='tight')

    if getattr(var, 'SHOW_ON_SCREEN', False): 
        plt.show()

    plt.close(fig)


def statistics_history(debug_pdf=None, history_chi2=None, history_dBIC=None, history_jitter=None, history_Kpl=None, history_N=None, jitter=None, kpl=None, niter=None, var=None):
    """Evolution of K, Delta BIC, reduced chi2 and activity jitter with the number of observations."""
    fig_track, axs = plt.subplots(2, 2, figsize=(12, 9), sharex=True)

    # 1. Delta BIC
    axs[0, 0].plot(history_N, history_dBIC, 'o-', color='purple', lw=2)
    axs[0, 0].axhline(6.0, color='red', linestyle='--', label='Confidence threshold (6.0)')
    axs[0, 0].set_ylabel(r'$\Delta$BIC (vs Ipotesi Nulla)')
    axs[0, 0].set_title('Evolution of the statistical evidence')
    axs[0, 0].grid(True, alpha=0.3); axs[0, 0].legend()

    # 2. Fitted vs injected K
    axs[0, 1].plot(history_N, history_Kpl, 'o-', color='blue', lw=2)
    axs[0, 1].axhline(kpl, color='green', linestyle='-', lw=2, label=f'True K ({kpl:.2f} m/s)')
    axs[0, 1].set_ylabel('Fitted K planet [m/s]')
    axs[0, 1].set_title('Stabilisation of the planet amplitude')
    axs[0, 1].grid(True, alpha=0.3); axs[0, 1].legend()

    # 3. Activity jitter
    axs[1, 0].plot(history_N, history_jitter, 'o-', color='orange', lw=2)
    axs[1, 0].axhline(jitter, color='black', linestyle=':', label=f'True jitter ({jitter:.2f} m/s)')
    axs[1, 0].set_ylabel('Activity jitter [m/s]')
    axs[1, 0].set_title('Evolution of the activity jitter')
    axs[1, 0].grid(True, alpha=0.3); axs[1, 0].legend()

    # 4. Reduced chi2
    axs[1, 1].plot(history_N, history_chi2, 'o-', color='green', lw=2)
    axs[1, 1].axhline(1.0, color='red', linestyle='-', alpha=0.5)
    axs[1, 1].set_ylabel(r'Reduced $\chi^2$')
    axs[1, 1].set_xlabel('Number of RV points')
    axs[1, 1].set_title('Quality of the global fit')
    axs[1, 1].grid(True, alpha=0.3)

    # Layout before drawing
    fig_track.tight_layout()

    # Force the rendering (avoids frozen windows with some backends)
    fig_track.canvas.draw()
    fig_track.canvas.flush_events()

    if debug_pdf is not None and niter < 3:
        debug_pdf.savefig(fig_track, bbox_inches='tight')

    if getattr(var, 'SHOW_ON_SCREEN', False):
        # short pause to refresh the GUI
        plt.pause(0.1) 
        plt.show(block=True) # keep the window open until it is closed

    plt.close(fig_track)


def extract_params_gls(res_gls):
  p_max = 1./res_gls.best['f']
  err_p_max = res_gls.best['e_P']
  maxpow = res_gls.pmax
  k_max = res_gls.best['amp']
  err_k = res_gls.best['e_amp']
  sig_K = k_max/err_k 
  fopt = res_gls.best['f']
  T0_fit = res_gls.best['T0']
  return p_max, err_p_max, maxpow, k_max, err_k, sig_K, fopt, T0_fit


def plot_rv_series(pdf, time_f, RVfinal, err, RVact, noise_rv, target_params, Prot, niter, fap_levels_func):
    fig, ax = plt.subplots(2, 1, gridspec_kw={'height_ratios': [2, 1]})
    
    # RV plots (top panel)
    ax[0].errorbar(time_f, RVfinal, yerr=err, ls='', marker='o', label='final', color='tab:green', markeredgecolor='k')
    ax[0].scatter(time_f, RVact, marker='o', edgecolor='k', label='stellar activity', zorder=30, alpha=0.6)
    ax[0].scatter(time_f, RVfinal - RVact - noise_rv, marker='o', edgecolor='k', label='planet', zorder=30, alpha=0.6)
    ax[0].set_title(f'Simulation {niter}')
    ax[0].legend(ncol=3)
    ax[0].set_xlabel('BJD')
    ax[0].set_ylabel('RV [m/s]')

    # Periodogram plots (bottom panel)
    pl_period = target_params.get('pl_period', 1)
    fend = 0.6 if pl_period > 2 else 2.5
    
    FAP = np.asarray([0.1, 0.01])
    color_FAP = ['red', 'orange']
    
    xx_plot = [1e-5, fend]
    
    signalgls = Gls((time_f, RVfinal, err), verbose=False, fbeg=xx_plot[0], fend=xx_plot[1])
    p_max, err_p_max, maxpow, k_max, err_k, sig_K, fopt, T0_fit = extract_params_gls(signalgls)
    fap_pow = fap_levels(FAP, time_f, 0.0, fend=1)
    freqs = signalgls.freq
    powers = signalgls.power


    ax[1].plot(freqs, powers, lw=0.5)
    for mm, m in enumerate(FAP):
        ax[1].axhline(fap_pow[mm], color=color_FAP[mm], ls='--', label='FAP %s%%'%(m*100)) 
        
    ax[1].axvline(1 / pl_period, c='orange', lw=4, alpha=0.4, label='planet')
    ax[1].axvline(1 / Prot, c='tab:blue', lw=4, alpha=0.4, label='rotation')
    ax[1].set_xlabel('Freq [c/d]')
    ax[1].set_ylabel('Power')
    ax[1].legend()

    pdf.savefig(fig)
    plt.close(fig)
