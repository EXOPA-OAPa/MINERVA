"""
Stellar properties of a synthetic star: parameter bootstrap, spectral type, age,
rotation, convective turnover time, chromospheric activity and amplitude of the
activity-induced RV signal (paper, Sect. 2.3).

Each function implements one relation: to change a relation, change its function.
"""
import numpy as np

from .tracks import ms_lifetime, track_pro


def bootstrap_parameters(target_params, rng):
    """Resample the stellar parameters within their uncertainties (modifies target_params in place).

    Each parameter is drawn from a Gaussian centred on its catalog value, with the
    catalog uncertainty as sigma, and bounded to physically meaningful values.
    """
    # 1. Stellar mass [Msun]
    if target_params.get('star_mass') is not None and target_params.get('err_star_mass') is not None:
        val_mass = rng.normal(target_params['star_mass'], target_params['err_star_mass'])
        target_params['star_mass'] = max(0.08, val_mass)

    # 2. Effective temperature [K]
    if target_params.get('star_temp') is not None and target_params.get('err_star_temp') is not None:
        val_temp = rng.normal(target_params['star_temp'], target_params['err_star_temp'])
        target_params['star_temp'] = max(2500, val_temp)  # lower limit for M dwarfs

    # 3. Stellar radius [Rsun]
    if target_params.get('star_radius') is not None and target_params.get('err_star_radius') is not None:
        val_rad = rng.normal(target_params['star_radius'], target_params['err_star_radius'])
        target_params['star_radius'] = max(0.1, val_rad)

    # 4. Age [Gyr]
    if target_params.get('age') is not None and target_params.get('err_age') is not None:
        val_age = rng.normal(target_params['age'], target_params['err_age'])
        target_params['age'] = max(0.01, val_age)  # no negative ages

    # 5. Rotation period [d]
    if target_params.get('Prot') is not None and target_params.get('err_Prot') is not None:
        val_prot = rng.normal(target_params['Prot'], target_params['err_Prot'])
        target_params['Prot'] = max(0.5, val_prot)

    # 6. log R'HK (chromospheric activity index)
    if target_params.get('log_Rhk') is not None and target_params.get('err_log_Rhk') is not None:
        val_log_rhk = rng.normal(target_params['log_Rhk'], target_params['err_log_Rhk'])
        # Upper limit: -3.5 is already extreme activity; the relations are not valid beyond it
        target_params['log_Rhk'] = min(-3.5, val_log_rhk)


def check_spectral_type(var, target_params):
    """Validate T_eff and assign or validate the spectral type.

    Returns False (after writing the error to the results file) if the star
    cannot be simulated, True otherwise. Sets target_params["sp"] if missing."""
    temp = target_params.get('star_temp')

    sp = target_params.get('sp')
    # Strip the string: blank or empty strings become empty, i.e. False
    if isinstance(sp, str):
        sp = sp.strip()

    # 1. Temperature limit, valid in every case: only FGKM stars
    if temp is not None and temp > 7500:
        error_msg = f"ERROR: The star is too hot (T={temp}K)! Must be an FGKM star (<= 7500K)."
        print(error_msg)
        with open(var.outputfilename, 'a') as f:
            f.write(error_msg + "\n")
        return False

    # The effective temperature is mandatory
    if temp is None:
        error_msg = f"ERROR: star_temp missing for {target_params.get('targetname')} (mandatory parameter)."
        print(error_msg)
        with open(var.outputfilename, 'a') as f:
            f.write(error_msg + "\n")
        return False

    # 2. Assign the spectral type from T_eff, or validate the one given in input
    if not sp:
        if temp > 6000:
            target_params['sp'] = 'F0'
        elif temp > 5000:
            target_params['sp'] = 'G0'
        elif temp >= 3500:
            target_params['sp'] = 'K0'
        else:
            target_params['sp'] = 'M0'

    else:
        if not isinstance(sp, str) or sp[0].upper() not in ['F', 'G', 'K', 'M']:
            error_msg = f"ERROR: Spectral type '{sp}' is invalid. Must start with F, G, K, or M."
            print(error_msg)
            with open(var.outputfilename, 'a') as f:
                f.write(error_msg + "\n")
            return False

    return True


def stellar_age(var, target_params, rng):
    """Stellar age in years: from the catalog, or drawn uniformly up to the end of the main sequence."""
    if target_params['age'] is None:
        # Upper limit: 10 Gyr below 1 Msun, the maximum age of max_age_from_mass.dat above
        if target_params['star_mass'] < 1:
            age_max = 1e+10
        else:
            okm = np.argmin(np.abs(var.mass_tmp - target_params['star_mass']))
            age_max = var.max_age[okm]
        # The age must not exceed the end of the main sequence of the track used for tau:
        # beyond it the track is a subgiant (lower T_eff, tau of hundreds of days),
        # inconsistent with the catalog T_eff. max_age_from_mass.dat includes the post-MS phase.
        t_ms = ms_lifetime(target_params['star_mass'], var.InDir)
        age_max = min(age_max, t_ms)

        # For very massive stars the main sequence can last less than 1e8 yr: avoid a negative range
        age_min = min(1e+8, 0.1 * age_max)
        age = age_min + (age_max - age_min) * rng.uniform(0, 1)
    else:
        age = target_params['age'] * 1e9
        t_ms = ms_lifetime(target_params['star_mass'], var.InDir)
        if age > t_ms:
            print(f"   [WARN] age {age/1e9:.2f} Gyr beyond the end of the main sequence ({t_ms/1e9:.2f} Gyr): "
                  f"tau from the subgiant phase")

    return age


def rotation_class(target_params, rng):
    """Rotation regime ('slow', 'intermediate', 'fast'): from the catalog, or drawn at random (25/50/25%)."""
    if target_params['vsini_tag'] is None:
        vsini = rng.uniform(0, 20)
        if vsini < 5:    tag = 'slow'
        elif vsini > 15: tag = 'fast'
        else:            tag = 'intermediate'
    else:
        tag = target_params['vsini_tag']

    return tag


def rotation_period(target_params, age, tag, rng):
    """Rotation period [d]: from the catalog, or from age and rotation regime (Gallet & Bouvier 2015).

    log Prot = A + B log(age), with the coefficients of Table 1 of the paper
    (G: 1 Msun, K: 0.8 Msun, M: 0.5 Msun tracks). F stars: log-uniform ranges in T_eff bins.
    """
    if target_params['Prot'] is None:

        # --- F stars ---
        if target_params['sp'][0] == 'F':
            if target_params['star_temp'] <= 6500:
                log_Prot = np.log10(9.) + (np.log10(12.) - np.log10(9.)) * rng.uniform(0, 1)
            elif ((target_params['star_temp'] > 6500) & (target_params['star_temp'] <= 7000)):
                log_Prot = np.log10(5.) + (np.log10(9.) - np.log10(5.)) * rng.uniform(0, 1)
            elif ((target_params['star_temp'] > 7000) & (target_params['star_temp'] <= 7500)):
                log_Prot = np.log10(2.) + (np.log10(5.) - np.log10(2.)) * rng.uniform(0, 1)

        # --- G stars ---
        elif target_params['sp'][0] == 'G':
            if ((age < 1.e9) & (tag == 'slow')):           log_Prot = -1.28 + 0.273 * np.log10(age)
            elif ((age < 1.e9) & (tag == 'intermediate')): log_Prot = -3.12 + 0.477 * np.log10(age)
            elif ((age < 1.e9) & (tag == 'fast')):         log_Prot = -12.99 + 1.574 * np.log10(age)
            else:                                          log_Prot = -3.53 + 0.523 * np.log10(age)

        # --- K stars ---
        elif target_params['sp'][0] == 'K':
            if ((age < 1.e9) & (tag == 'slow')):           log_Prot = -1.03 + 0.234 * np.log10(age)
            elif ((age < 1.e9) & (tag == 'intermediate')): log_Prot = -3.22 + 0.477 * np.log10(age)
            elif ((age < 1.e9) & (tag == 'fast')):         log_Prot = -10.63 + 1.301 * np.log10(age)
            else:                                          log_Prot = -3.10 + 0.465 * np.log10(age)

        # --- M stars ---
        elif target_params['sp'][0] == 'M':
            if ((age < 3.e9) & (tag == 'slow')):           log_Prot = -3.71 + 0.539 * np.log10(age)
            elif ((age < 3.e9) & (tag == 'intermediate')): log_Prot = -6.44 + 0.827 * np.log10(age)
            elif ((age < 3.e9) & (tag == 'fast')):         log_Prot = -11.43 + 1.354 * np.log10(age)
            else:                                          log_Prot = -1.31 + 0.301 * np.log10(age)

        Prot = 10**log_Prot
    else:
        Prot = target_params['Prot']
        log_Prot = np.log10(Prot)

    return Prot, log_Prot


def rossby_number(var, target_params, age, Prot):
    """Convective turnover time from the YaPSI tracks and Rossby number R0 = Prot / tau."""
    tau = track_pro(target_params['star_mass'], age, InDir=var.InDir)
    R0 = Prot / tau
    print('Convective turnover time tau = ', tau)

    return tau, R0


def activity_level(target_params, tau, R0, log_Prot, Prot):
    """Chromospheric activity log R'HK: from the catalog, or from the Rossby number (FGK, Mamajek &
    Hillenbrand 2008) or the rotation period (M, Astudillo-Defru et al. 2017). Stars without a
    convective envelope (tau ~ 0) are flagged as inactive (log R'HK = NaN)."""
    is_active = True
    if target_params['log_Rhk'] is None:
        if (tau != None):
            if (tau < 0.001):
                # No convective envelope -> no spots.
                # NaN rather than 0: 0 would be read as extremely high activity by the comparisons below
                log_Rhk = np.nan   # K_act = 0 is set by activity_amplitude
                is_active = False
            else:
                # FGK stars: log R'HK from the Rossby number (Mamajek & Hillenbrand 2008)
                if R0 < 0.4:
                    log_Rhk = -4.230 - 1.451 * (R0 - 0.233)
                else:
                    log_Rhk = -4.522 - 0.337 * (R0 - 0.814)
                # M stars: log R'HK from the rotation period (Astudillo-Defru et al. 2017)
                if target_params['sp'][0] == 'M':
                    if Prot > 10:
                        log_Rhk = -1.509 * log_Prot - 2.550
                    else:
                        log_Rhk = -4.045
                    print("log R'HK of the M star = ", log_Rhk)
        else:
            raise ValueError("Error: cannot compute log R'HK without the convective turnover time tau")
    else:
        log_Rhk = target_params['log_Rhk']

    return log_Rhk, is_active


def activity_amplitude(target_params, log_Rhk, is_active):
    """Amplitude of the activity-induced RV signal [m/s]: from the catalog, or from log R'HK and T_eff
    (Santos et al. 2000; exponential GAPS/ArMS calibration for K dwarfs, capped at log R'HK = -3.99)."""
    if target_params['K_act'] is None:
        if is_active:
            # Relations of Santos et al. (2000, CORALIE survey IV), with R5 = 1e5 R'HK
            R5 = 1.e5 * 10**log_Rhk
            if   (target_params['star_temp'] > 6000) & (target_params['star_temp'] <= 7500): K_act = 9.2 * R5**(0.75)
            elif (target_params['star_temp'] > 5000) & (target_params['star_temp'] <= 6000): K_act = 7.9 * R5**(0.55)
            # K dwarfs: exponential relation calibrated on GAPS/ArMS K stars (MINERVA paper, Eq. 3).
            # Valid in the range of the data, log R'HK ~ [-5.03, -3.99]: no extrapolation beyond it.
            elif (target_params['star_temp'] > 4000) & (target_params['star_temp'] <= 5000):
                LOGRHK_K_MAX = -3.99
                log_rhk_eff = min(log_Rhk, LOGRHK_K_MAX)
                if log_rhk_eff < log_Rhk:
                    print(f"   [WARN] log R'HK = {log_Rhk:.3f} beyond the calibration range: using {LOGRHK_K_MAX}")
                R5_eff = 1.e5 * 10**log_rhk_eff
                K_act = 1.558 * np.exp(0.605 * R5_eff)
            elif (target_params['star_temp'] >= 3500) & (target_params['star_temp'] <= 4000): K_act = 7.8 * R5**(0.13)
            elif (target_params['star_temp'] < 3500) & (log_Rhk <= -5):     K_act = 3.5
            elif (target_params['star_temp'] < 3500) & (log_Rhk > -5):      K_act = 5.0
            else:
                raise ValueError(f"K_act cannot be computed: T_eff={target_params['star_temp']} outside the expected ranges")

        else:
            K_act = 0
    else:
        K_act = target_params['K_act']

    return K_act
