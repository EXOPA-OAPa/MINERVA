"""
Models and fits of the sequential analysis (paper, Sect. 2.6).

The planet is a circular orbit with period and transit time fixed to the
ephemerides. Stellar activity is a sinusoid at the stellar period, with free
amplitude and phase in each observing season (seasonal model) or one global
sinusoid for quiet or inactive stars. Fits use scipy.optimize.curve_fit (trf).
"""
import numpy as np
from scipy.optimize import curve_fit


# =============================================================================
# Models
# =============================================================================
def model_planet_only(t, gamma, K_pl, P_pl, T0_pl):
    """Circular orbit: gamma - K sin(2 pi (t - T0) / P), zero at the transit and decreasing."""
    omega = 2 * np.pi / P_pl
    return gamma - K_pl * np.sin(omega * (t - T0_pl))


def seasonal_star_model(season_boundaries, omega_star):
    """Star-only seasonal model (null hypothesis of the simultaneous fit).

    Parameters of the returned model: (gamma, A, B) for each season.
    """
    def model_star_only_seasonal(t, *params_star):
        rv_star_null = np.zeros_like(t)
        for s, (t_min_s, t_max_s) in enumerate(season_boundaries):
            gamma_s = params_star[0 + 3*s]
            A_s     = params_star[1 + 3*s]
            B_s     = params_star[2 + 3*s]

            mask = (t >= t_min_s) & (t <= t_max_s)
            if np.any(mask):
                rv_star_null[mask] = gamma_s + A_s * np.sin(omega_star * t[mask]) + B_s * np.cos(omega_star * t[mask])
        return rv_star_null
    return model_star_only_seasonal


def seasonal_simultaneous_model(P_PLANET, T0_PLANET, season_boundaries, omega_star):
    """Planet (fixed P and T0) + seasonal stellar sinusoid: parameters K_pl, then (gamma, A, B) per season."""
    def model_simultaneous_seasonal(t, *params):
        K_pl = params[0]
        rv_pl = -K_pl * np.sin((2 * np.pi / P_PLANET) * (t - T0_PLANET))
        rv_star = np.zeros_like(t)
        for s, (t_min_s, t_max_s) in enumerate(season_boundaries):
            gamma_s = params[1 + 3*s]
            A_s     = params[2 + 3*s]
            B_s     = params[3 + 3*s]

            mask = (t >= t_min_s) & (t <= t_max_s)
            if np.any(mask):
                rv_star[mask] = gamma_s + A_s * np.sin(omega_star * t[mask]) + B_s * np.cos(omega_star * t[mask])
        return rv_pl + rv_star
    return model_simultaneous_seasonal


# =============================================================================
# Fits
# =============================================================================
def activity_seasons(is_active, log_Rhk, t_fit, time_f, ind, all_season_indices):
    """Seasons for the activity model: one global season for quiet or inactive stars,
    otherwise the observing seasons with at least 3 points.

    Returns current_seasons_idx, N_seasons, season_boundaries, fit_desc.
    """
    if (not is_active) or log_Rhk < -5.0:
        # Quiet or inactive star: one global season (phase-coherent signal)
        current_seasons_idx = [np.arange(len(t_fit))]
        N_seasons = 1
        season_boundaries = [(np.min(time_f), np.max(time_f))]
        fit_desc = f"SIMULTANEOUS GLOBAL"
        print(f"   > Quiet star (logRhk = {log_Rhk:.2f}): global fit with a single phase.")
    else:
        # Active star: independent amplitude and phase in each observing season
        current_seasons_idx = []
        season_boundaries = []
        for s_global_idx in all_season_indices:
            local_idx = np.where(np.isin(ind, s_global_idx))[0]
            if len(local_idx) >= 3:
                current_seasons_idx.append(local_idx)
                t_min_s = np.min(time_f[s_global_idx])
                t_max_s = np.max(time_f[s_global_idx])
                season_boundaries.append((t_min_s, t_max_s))

        N_seasons = len(current_seasons_idx)
        fit_desc = f"SIMULTANEOUS SEASONAL ({N_seasons} Seasons)"
        print(f"   > Active star (logRhk = {log_Rhk:.2f}): seasonal fit ({N_seasons} blocks).")

    return current_seasons_idx, N_seasons, season_boundaries, fit_desc


def fit_simultaneous(t_fit, rv_fit, err_fit, current_seasons_idx, N_seasons,
                     model_simultaneous_seasonal, model_star_only_seasonal,
                     P_PLANET, T0_PLANET, K_star_fit):
    """Simultaneous fit planet + seasonal activity, and star-only null hypothesis.

    If the fit fails, popt is None and use_simultaneous is False (the caller then
    falls back to the planet-only fit). K_star_fit keeps the value reached before
    the failure, as in the version used for the paper.
    Returns popt, K_fit, K_err_formal, dof, k_params, bic_null, y_model_final,
    rv_plot_data, K_star_fit, use_simultaneous.
    """
    use_simultaneous = True
    popt = K_fit = K_err_formal = dof = k_params = bic_null = y_model_final = rv_plot_data = None
    guess_rms = np.std(rv_fit)

    # Initial guesses and bounds: K_pl first (K_pl >= 0), then (gamma, A, B) for each season
    p0 = [guess_rms * 0.5]
    lower_bounds = [0.0]
    upper_bounds = [np.inf]

    for s, loc_idx in enumerate(current_seasons_idx):
        gamma_guess = np.median(rv_fit[loc_idx])  # seasonal median as offset guess
        p0.extend([gamma_guess, guess_rms, 0.0])

        lower_bounds.extend([-np.inf, -np.inf, -np.inf])
        upper_bounds.extend([np.inf, np.inf, np.inf])

    bounds = (lower_bounds, upper_bounds)

    try:
        # --- Full model: planet + seasonal activity ---
        popt, pcov = curve_fit(model_simultaneous_seasonal, t_fit, rv_fit,
                               p0=p0, sigma=err_fit, absolute_sigma=True, bounds=bounds, maxfev=10000)

        K_fit = popt[0]
        K_err_formal = np.sqrt(pcov[0, 0])

        y_model_final = model_simultaneous_seasonal(t_fit, *popt)

        # For the figures: planet-only component (centred on 0) and data cleaned of the stellar component
        y_plot_model = model_planet_only(t_fit, 0.0, K_fit, P_PLANET, T0_PLANET)
        rv_star_model = y_model_final - y_plot_model
        rv_plot_data = rv_fit - rv_star_model

        # Mean stellar amplitude over the fitted seasons (for the results file)
        A_vals = popt[2::3]
        B_vals = popt[3::3]
        K_star_fit = np.mean([np.sqrt(a**2 + b**2) for a, b in zip(A_vals, B_vals)])

        # Number of free parameters
        k_params = 1 + 3 * N_seasons
        dof = len(t_fit) - k_params

        # --- Null hypothesis: star only, same seasons (starts from the full-model solution) ---
        p0_null = popt[1:]
        popt_null, _ = curve_fit(model_star_only_seasonal, t_fit, rv_fit, p0=p0_null, sigma=err_fit, absolute_sigma=True, maxfev=10000)
        y_null = model_star_only_seasonal(t_fit, *popt_null)

        chisq_null = np.sum(((rv_fit - y_null) / err_fit)**2)
        k_null = 3 * N_seasons
        bic_null = chisq_null + k_null * np.log(len(t_fit))

    except Exception as e:
        print(f"   ⚠️ Seasonal simultaneous fit failed: {e}. Falling back to the planet-only model.")
        popt = None
        use_simultaneous = False  # the caller runs the planet-only fit in the same step

    return popt, K_fit, K_err_formal, dof, k_params, bic_null, y_model_final, rv_plot_data, K_star_fit, use_simultaneous


def fit_planet_only(t_fit, rv_fit, err_fit, P_PLANET, T0_PLANET):
    """Planet-only fit (gamma, K) and flat null hypothesis (weighted mean).

    If the fit fails, popt is None.
    Returns fit_desc, popt, K_fit, K_err_formal, dof, k_params, bic_null, y_model_final, rv_plot_data.
    """
    popt = K_fit = K_err_formal = dof = k_params = bic_null = y_model_final = rv_plot_data = None
    fit_desc = "PLANET ONLY"
    func_to_fit = lambda t, g, k: model_planet_only(t, g, k, P_PLANET, T0_PLANET)

    p0 = [np.median(rv_fit), np.std(rv_fit)*1.4]   # [gamma, K]
    bounds = ([-np.inf, 0.0], [np.inf, np.inf])     # K >= 0
    try:
        popt, pcov = curve_fit(func_to_fit, t_fit, rv_fit, p0=p0, sigma=err_fit, absolute_sigma=True, bounds=bounds)

        gamma_fit, K_fit = popt
        K_err_formal = np.sqrt(pcov[1, 1])

        k_params = 2  # gamma, K_pl
        dof = len(t_fit) - k_params

        y_model_final = model_planet_only(t_fit, *popt, P_PLANET, T0_PLANET)

        # Planet-only case: the data shown in the figures are the raw data
        rv_plot_data = rv_fit

        # Null hypothesis: flat (weighted mean)
        w = 1.0 / err_fit**2
        gamma_flat = np.sum(rv_fit * w) / np.sum(w)
        chisq_flat = np.sum(((rv_fit - gamma_flat) / err_fit)**2)
        k_null = 1  # gamma
        bic_null = chisq_flat + k_null * np.log(len(t_fit))
    except Exception as e:
        print(f"   ⚠️ Planet-only fit failed: {e}")
        popt = None

    return fit_desc, popt, K_fit, K_err_formal, dof, k_params, bic_null, y_model_final, rv_plot_data
