"""
Synthetic RV signals: photon and white noise, observing cadence, stellar activity
signal and formal error bars (paper, Sects. 2.1, 2.2 and 2.4).
"""
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.optimize import curve_fit
from scipy.stats import norm

from . import keplerian


def photon_noise(var, target_params):
    """Photon-noise RV uncertainty for the star, from the NSTS table (V magnitude, spectral class)."""
    v_mag = target_params['v_mag']
    sp_class = target_params['sp'][0].upper()  # first letter of the spectral type

    # 1. Closest magnitude of the grid: robust to rounding, and it clips to the edges
    #    of the table when v_mag is outside its range.
    idx = np.argmin(np.abs(var.tmp_mag - v_mag))

    # 2. Photon-noise column of the spectral class
    err_arrays = {
        'F': var.errRV_F,
        'G': var.errRV_G,
        'K': var.errRV_K,
        'M': var.errRV_M
    }

    # 3. Photon noise at that magnitude (G column as a fallback for an unknown class)
    sig_ph = err_arrays.get(sp_class, var.errRV_G)[idx]

    return sig_ph


def white_noise(var, target_params, sig_ph, rng):
    """Gaussian white noise: instrumental floor, additional planets and photon noise in quadrature."""
    sigma = np.sqrt(var.sig_floor**2 + var.sig_pl**2 + sig_ph**2)
    # One value per epoch of the cadence used for this planet (USP cadence for P < 2 d)
    if (target_params['pl_period'] >= 2.): time_arr = var.time_sampling
    else: time_arr = var.time_sampling_usp

    unit = rng.normal(0, 1, size=np.size(time_arr))
    noise_rv = unit * sigma

    return noise_rv


def select_cadence(var, target_params):
    """Observing epochs: the USP cadence for P < 2 d, the standard one otherwise."""
    if target_params['pl_period'] >= 2:
        time_f = var.time_sampling
    else:
        time_f = var.time_sampling_usp

    return time_f


def activity_signal(var, target_params, time_f, Prot, K_act, log_Rhk, is_active, rng):
    """Activity-induced RV signal on the observing epochs (zero for inactive stars)."""
    if is_active and K_act > 0:
        RVact = create_stellar_activity_signal_dynamic_complex(Prot, time_f, K_act, log_Rhk, target_params, var, rng)
    else:
        # Star without activity (no convective envelope): no signal
        RVact = np.zeros_like(time_f)

    return RVact


def error_bars(var, sig_ph, log_Rhk, is_active, time_f, subset_for_simulation, debug_pdf, save_this_plot, show_on_scr, rng):
    """Formal RV error bars: instrumental and photon noise plus, for active stars, the activity-driven
    extraction term drawn from the GAPS database (paper, Sect. 2.2).

    The additional-planet term sig_pl is part of the injected noise but not of the error bars:
    an observer does not know about it.
    """
    errore_base = np.sqrt(var.sig_floor**2 + sig_ph**2)
    if is_active:
        total_errors, jitter = simulate_activity_error(
            target_rhk=log_Rhk,          # activity level of the star
            n_obs=len(time_f),           # number of epochs
            df_data=subset_for_simulation,
            sigma_inst=errore_base,      # instrumental + photon noise
            var_pct=0.10,                # 10% epoch-to-epoch scatter of the error bars
            window=0.20,                 # half-width of the log R'HK window in the GAPS database
            pdf_obj=debug_pdf if save_this_plot else None,
            show_plot=show_on_scr,
            rng=rng
        )
    else:
        # Inactive star: instrumental + photon noise only, with a 10% scatter so that
        # the error bars are not all identical
        base_sigma = np.sqrt(var.sig_floor**2 + sig_ph**2)
        sigma_var = base_sigma * 0.10
        # These are error-bar lengths, not noise values: they must be positive
        total_errors = np.abs(rng.normal(base_sigma, sigma_var, size=len(time_f)))
        jitter = 0

    err = total_errors.copy()

    return err, jitter


def create_stellar_activity_signal_dynamic_complex(Prot, time_f, K_act, log_Rhk, target_params, var, rng):
    """
    Activity-induced RV signal (paper, Sect. 2.4): rotational modulation with differential
    rotation (beating between Prot and Prot(1+alpha)) and, for intermediate and high activity,
    new random phases and amplitudes in each observing season (finite lifetime of active
    regions) plus the first harmonic for high activity.
    """
    
    # output signal
    RVact_final = np.zeros_like(time_f)



    season_indices = var.seasons
    
    # Spots at different latitudes rotate at different rates (differential rotation).
    # alpha = relative shear (e.g. 0.1 = 10% difference in period); it decreases with the rotation
    # rate (Barnes et al. 2005; Reinhold et al. 2013), so it is small for active stars (paper, Sect. 2.4).
    if log_Rhk >= -4.5:
        alpha = rng.uniform(0.01, 0.08)   # high activity: small relative shear
    else:
        alpha = rng.uniform(0.10, 0.20)   # low/intermediate activity: solar-like shear

    # secondary period (higher latitudes)
    Prot_2 = Prot * (1.0 + alpha) 

    # --- 1. Low activity: one phase-coherent sinusoid with a weak modulation ---
    if log_Rhk < -5.0:
        # a weak secondary component (10% of K_act), so that the signal
        # is not a mathematically perfect sinusoid
        orbd1 = np.asarray([Prot, time_f[0], 0, 0, K_act])
        rv1 = keplerian.rv_drive(time_f, orbd1)
        
        orbd2 = np.asarray([Prot_2, time_f[0], 0, 0, K_act * 0.1])
        rv2 = keplerian.rv_drive(time_f, orbd2)
        
        RVact_final = rv1 + rv2

    # --- 2. Intermediate and high activity ---
    # Each observing season has its own active regions, born with random phases
    else: 
        # complexity of the signal by activity regime
        if log_Rhk <= -4.5:  # intermediate activity
            ratio_p2 = 0.4  # secondary signal at 40% of the amplitude
            include_harmonic = False
        else:  # high activity
            ratio_p2 = 0.8  # secondary signal at 80% of the amplitude, plus the first harmonic
            include_harmonic = True

        for idx_season in season_indices:
            # skip empty seasons
            if len(idx_season) == 0:
                continue

            # --- A. Finite lifetime of active regions: new phases in each season ---
            # (old spots have decayed, new ones have emerged)
            t0_base = rng.uniform(0, Prot)
            t0_diff = rng.uniform(0, Prot_2)  # independent phase of the secondary component
            
            # seasonal amplitude: N(K_act, 0.1 K_act), at least 0.5 m/s
            k_season = rng.normal(K_act, K_act * 0.1)
            k_season = max(0.5, k_season)

            # --- B. Differential rotation: beating within the season ---
            # signal 1: Prot
            orbd1 = np.asarray([Prot, time_f[0] + t0_base, 0, 0, k_season])
            rv_season = keplerian.rv_drive(time_f[idx_season], orbd1)
            
            # signal 2: Prot(1+alpha): the beating changes amplitude and shape within the season
            orbd2 = np.asarray([Prot_2, time_f[0] + t0_diff, 0, 0, k_season * ratio_p2])
            rv_season += keplerian.rv_drive(time_f[idx_season], orbd2)

            # --- C. First harmonic (high activity only) ---
            # Prot/2: several spot groups on opposite hemispheres
            if include_harmonic:
                orbd3 = np.asarray([Prot / 2.0, time_f[0] + rng.uniform(0, Prot), 0, 0, k_season * 0.3])
                rv_season += keplerian.rv_drive(time_f[idx_season], orbd3)

            # signal of this season
            RVact_final[idx_season] = rv_season

    return RVact_final


def add_planetary_signal(target_params,eccentricity,periastron,time_f,noise_rv,RVact,sig_ph):
  # Planet
  p_yr = (target_params['pl_period']/365.25)
  pl_angle = np.pi/2
 
  kpl = 28.4*(target_params['star_mass'])**(-2/3)*(p_yr)**(-1/3)*target_params['pl_mass_Mj']*np.sin(pl_angle)/(np.sqrt(1-eccentricity**2))
  print('Kp = ', kpl)
  # Time of periastron from the transit time (exact conversion; for e = 0 it is Tc - P/4).
  # The argument of periastron is given in degrees and converted to radians.
  omega_rad = np.deg2rad(periastron)
  t0 = keplerian.timetrans_to_timeperi(target_params['transit_midpoint'], target_params['pl_period'],
                                          eccentricity, omega_rad)
  orbd = np.asarray([target_params['pl_period'], t0, eccentricity, omega_rad, kpl])

  print('planet' , target_params['pl_period'],t0,kpl)
  RVpl = keplerian.rv_drive(time_f, orbd)
  # total signal: noise + activity + planet
  RVfinal = noise_rv + RVact + RVpl
  err = np.zeros_like(time_f)+sig_ph

  return RVfinal, err, kpl, RVpl


def func_broken_line(x, x_break, y_break, m_quiet, m_active):
    """
    Broken line: jitter vs log R'HK.
    x_break:  log R'HK where the slope changes (e.g. -4.8)
    y_break:  jitter at the break
    m_quiet:  slope below the break (expected low/flat)
    m_active: slope above the break (expected steep)
    """
    return np.where(x < x_break, 
                    m_quiet  * (x - x_break) + y_break,  # below the break
                    m_active * (x - x_break) + y_break)  # above the break


def simulate_activity_error(target_rhk, n_obs, df_data, sigma_inst, window=0.15, var_pct=0.20, pdf_obj=None, show_plot=False, rng=None):
    
    # random generator of the caller (reproducible); without it, an unseeded one
    if rng is None:
        rng = np.random.default_rng()
    mode = "CDF"
    base_jitter = 0
    sigma_resid = 0
    popt = None
    residuals = []
    n_points = 0

    max_rhk_data = df_data['logRHK'].max()

    # Very active star, or outside the data: use the broken-line fit (extrapolation)
    ACTIVE_THRESHOLD = -4.5 
    force_fit = (target_rhk > max_rhk_data) or (target_rhk > ACTIVE_THRESHOLD)

    # Branch 1: broken-line fit
    if force_fit:
        mode = "BROKEN"
        print(f"Target ({target_rhk}) in the active/extrapolation regime: using the broken-line fit.")
        
        # all the data of the spectral class, to locate the break
        fit_data = df_data.dropna(subset=['logRHK', 'Jitter_ms'])
        x_data = fit_data['logRHK'].values
        y_data = fit_data['Jitter_ms'].values

        try:
            # parameters: [x_break, y_break, m_quiet, m_active]
            
            
            # 1. ranges of the data, for the initial guesses
            x_min, x_max = np.min(x_data), np.max(x_data)
            y_min, y_max = np.min(y_data), np.max(y_data)
            x_median = np.median(x_data)
            y_median = np.median(y_data)

            # 2. initial guesses
            # x_break: median of the data
            guess_x_break = x_median
            
            # y_break: median of the jitter
            guess_y_break = y_median
            
            # m_quiet: (Y_median - Y_min) / (X_median - X_min)
            # 1e-6 avoids divisions by zero
            guess_m_quiet = (y_median - y_min) / (x_median - x_min + 1e-6)
            
            # m_active: (Y_max - Y_median) / (X_max - X_median)
            guess_m_active = (y_max - y_median) / (x_max - x_median + 1e-6)
            
            # keep the guesses positive
            p0 = [guess_x_break, guess_y_break, max(0.1, guess_m_quiet), max(1.0, guess_m_active)]

            # 3. bounds
            # x_break within the range of the data
            # y_break between 0 and slightly above the observed maximum
            # m_quiet may be slightly negative
            # m_active positive (the jitter grows with activity)
            
            # 5% margin, so that the fit does not stick to the edges
            x_range = x_max - x_min
            bound_x_low = x_min + 0.05 * x_range
            bound_x_high = x_max - 0.05 * x_range

            lower_bounds = [bound_x_low,  0.0,         -2.0,   0.0]
            upper_bounds = [bound_x_high, y_max * 2.0, 50.0, 200.0]
            
            popt, pcov = curve_fit(func_broken_line, x_data, y_data, p0=p0, bounds=(lower_bounds, upper_bounds))

            # best-fit parameters
            x_knee, y_knee, m1, m2 = popt
            
            # 1. predicted mean jitter at the star's log R'HK
            predicted_mean = func_broken_line(target_rhk, *popt)
            
            # 2. scatter of the residuals
            model_values = func_broken_line(x_data, *popt)
            residuals = y_data - model_values
            sigma_resid = np.std(residuals)
            
            # 3. random draw around the prediction
            u_star = rng.normal(0, 1)  
            base_jitter = predicted_mean + (u_star * sigma_resid)
            
            if base_jitter < 0.5: base_jitter = 0.5
                
            print(f"📈 Fit Broken Line:")
            print(f"   Knee at logRhk = {x_knee:.3f} (Jitter={y_knee:.2f} m/s)")
            print(f"   Quiet slope: {m1:.2f} | Active slope: {m2:.2f}")
            print(f"   Predicted mean value = {predicted_mean:.2f} m/s")
            print(f"🎲 Drawn value: {base_jitter:.2f} m/s")
            
        except Exception as e:
            print(f"Fit failed ({e}). Falling back to the maximum jitter.")
            base_jitter = fit_data['Jitter_ms'].max()
            sigma_resid = 0

    # Branch 2: local window in log R'HK (empirical CDF)
    else:
        rhk_min = target_rhk - window; rhk_max = target_rhk + window
        mask = (df_data['logRHK'] >= rhk_min) & (df_data['logRHK'] <= rhk_max)
        local_data = df_data[mask]
        
        if len(local_data) < 5:
            rhk_min -= 0.1; rhk_max += 0.1
            mask = (df_data['logRHK'] >= rhk_min) & (df_data['logRHK'] <= rhk_max)
            local_data = df_data[mask]

        local_jitters = local_data['Jitter_ms'].values
        local_jitters = local_jitters[local_jitters > 0.001]
        
        if len(local_jitters) > 0:
            sorted_jitters = np.sort(local_jitters)
            n_points = len(sorted_jitters)
            u_star = rng.uniform(0, 1)  
            base_jitter = np.quantile(local_jitters, u_star, method='linear')
            n_points = len(sorted_jitters)
        else:
            base_jitter = 0.5
            n_points = 0

    '''
    sigma_var = base_jitter * var_pct
    nightly_jitters = np.abs(np.zeros(n_obs) + base_jitter)
    total_errors = np.sqrt(sigma_inst**2 + nightly_jitters**2)
    '''
    # epoch-to-epoch scatter of the activity term
    sigma_var = base_jitter * var_pct
    jitters_random = rng.normal(base_jitter, sigma_var, n_obs)         
    # a jitter cannot be negative
    nightly_jitters = np.abs(jitters_random)
    total_errors = np.sqrt(sigma_inst**2 + nightly_jitters**2)
    
    
    sns.set_style("whitegrid")
    # the figure is built only if it will be saved or shown
    if pdf_obj is None and not show_plot:
        return total_errors, base_jitter

    fig, axes = plt.subplots(2, 2, figsize=(16, 10))
    ax1, ax2 = axes[0]; ax3, ax4 = axes[1]
    
    # --- Panel 1 ---
    ax1.scatter(df_data['logRHK'], df_data['Jitter_ms'], c='lightgray', alpha=0.5, s=15, label='Dataset')
    ax1.axvline(target_rhk, c='black', ls='-', label='Target')

    if mode == "BROKEN" and popt is not None:
        # broken-line model
        x_plot_min = -5.2 
        x_plot_max = max(target_rhk + 0.2, df_data['logRHK'].max())
        x_model = np.linspace(x_plot_min, x_plot_max, 200)
        y_model = func_broken_line(x_model, *popt)
        
        ax1.plot(x_model, y_model, 'r-', lw=2, label='Broken-line fit')
        ax1.fill_between(x_model, y_model - sigma_resid, y_model + sigma_resid, color='red', alpha=0.1)
        
        # break point
        ax1.scatter([popt[0]], [popt[1]], c='red', marker='o', s=50, zorder=5, label='Knee')
        ax1.set_title(f"1. Broken-line fit (knee @ {popt[0]:.2f})", fontsize=12, fontweight='bold')
    else:
        rhk_min = target_rhk - window; rhk_max = target_rhk + window
        ax1.axvspan(rhk_min, rhk_max, color='red', alpha=0.1)
        ax1.set_title("1. Local window in log R'HK", fontsize=12, fontweight='bold')

    ax1.scatter(target_rhk, base_jitter, c='gold', marker='*', s=300, edgecolors='k', zorder=10, label='Simulation')
    ax1.legend()

    # --- Panel 2 ---
    if mode == "BROKEN":
        sns.histplot(residuals, kde=True, ax=ax2, color='orange', stat='density')
        ax2.set_title(f"2. Fit residuals (sigma={sigma_resid:.2f})")
    elif n_points > 0:
        y_prob = np.arange(1, n_points+1) / n_points
        ax2.step(sorted_jitters, y_prob, where='post', c='blue', alpha=0.4, lw=1)
        ax2.set_title("2. Draw from the local CDF")
    
    # --- Panels 3 and 4 ---
    x_g = np.linspace(base_jitter - 4*sigma_var, base_jitter + 4*sigma_var, 100)
    ax3.plot(x_g, norm.pdf(x_g, base_jitter, sigma_var), 'r-', lw=2)
    ax3.set_title(f"3. Epoch-to-epoch scatter ({var_pct*100:.0f}%)")
    
    ax4.plot(np.arange(1, n_obs + 1), total_errors, 'ko-', lw=1)
    ax4.set_title("4. Total estimated error")

    plt.tight_layout()
    # Save to the PDF and/or show on screen
    if pdf_obj is not None:
        pdf_obj.savefig(fig, bbox_inches='tight')
        
    if show_plot:
        plt.show()
        
    plt.close(fig)  # always close the figure
    
    return total_errors, base_jitter
