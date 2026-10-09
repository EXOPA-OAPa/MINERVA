"""Input and configuration: reading the input files and the YAML configuration,
parameters of the targets (catalog or single target), observing seasons of a cadence."""
import numpy as np
import sys 
import yaml
import pandas as pd 
from astropy.constants import M_jup, M_earth, R_jup, R_earth

def title():
  print('                                                               ') 
  print(' _______ _________ _        _______  _______           _______ ') 
  print(r'(       )\__   __/( (    /|(  ____ \(  ____ )|\     /|(  ___  )') 
  print(r'|       |   ) (   |  \  ( || (    \/| (    )|| )   ( || (   ) |')
  print(r'| || || |   | |   |   \ | || (__    | (____)|| |   | || (___) |')
  print(r'| ||_|| |   | |   | (\ \) ||  __)   |     __)( (   ) )|  ___  |')    
  print(r'| |   | |   | |   | | \   || (      | (\ (    \ \_/ / | (   ) |')
  print(r'| )   ( |___) (___| )  \  || (____/\| ) \ \__  \   /  | )   ( |')
  print(r'|/     \|\_______/|/    )_)(_______/|/   \__/   \_/   |/     \|')
  print('                                                               ') 
  print('                                                               ') 
  print('                                                               ') 
  print('                    Welcome to MINERVA!!!                      ')
  sys.stdout.flush()  # show the banner immediately


def import_file_catalog(namecatalog,InDir=None):
  if InDir == None:
    InDir = './utility_files/'
  # separator: ',' or ';'
  catalog = pd.read_csv(InDir+namecatalog,header=0,sep='[,;]')
  return catalog 

def _is_missing(x):
    """True if the value is None, NaN or 0 (0 was used to mean 'missing' in the paper version)."""
    if x is None:
        return True
    try:
        return bool(np.isnan(x)) or x == 0
    except TypeError:
        return False


def resolve_magnitude(target_params):
    """
    V magnitude used for the photon noise: v_mag if available, otherwise gaia_mag,
    otherwise 10 (rule documented in config.yaml). Modifies target_params in place.
    """
    if _is_missing(target_params.get('v_mag')):
        if not _is_missing(target_params.get('gaia_mag')):
            print(f"[WARN] {target_params.get('targetname')}: v_mag missing, using gaia_mag = {target_params['gaia_mag']}")
            target_params['v_mag'] = target_params['gaia_mag']
        else:
            print(f"[WARN] {target_params.get('targetname')}: v_mag and gaia_mag missing, using v_mag = 10")
            target_params['v_mag'] = 10.0
    return target_params


def load_target_params_from_catalog(catalog, i):
    """
    Parameters of the i-th target of the catalog (columns described in README.md).
    Missing planetary masses and radii are converted from the other unit; missing
    uncertainties get default values.
    """
    def safe_get(col, default=np.nan):
        """Value of the column for target i, or the default if the column does not exist."""
        if col in catalog.columns:
            return catalog[col].iloc[i]
        return default

    # Cadence file: column 'Cadence FileName' or, for older catalogs, 'Catalog FileName'
    cadence_file = safe_get('Cadence FileName')
    if cadence_file is None or pd.isna(cadence_file):
        cadence_file = safe_get('Catalog FileName')

    target_params = {
        "targetname": safe_get('Planet Name', None),
        "star_mass": safe_get('Star Mass [Ms]'),
        "err_star_mass": safe_get('Star Mass Error [Ms]', None),
        "star_radius": safe_get('Star Radius [Rs]'),
        "err_star_radius": safe_get('Star Radius Error [Rs]', None),
        "star_temp": safe_get('Star Temperature [K]'),
        "err_star_temp": safe_get('Star Temperature Error [K]', None),
        "sp": safe_get('Star Spectral Type', None),
        "v_mag": safe_get('Star V Mag'),
        "gaia_mag": safe_get('Star Gaia Mag'),
        "pl_mass_Mj": float(safe_get('Planet Mass [Mj]')) if safe_get('Planet Mass [Mj]') not in [None, '', 'NaN'] else np.nan,

        "pl_mass_Me": safe_get('Planet Mass [Me]'),
        "pl_radius_Rj": safe_get('Planet Radius [Rj]'),
        "pl_radius_Re": safe_get('Planet Radius [Re]'),
        "pl_period": safe_get('Planet Period [days]'),
        "err_pl_period": safe_get('Planet Period Error [days]'),
        "transit_midpoint": safe_get('Transit Mid Time [JD - 2450000]'),
        "age": safe_get('Star Age [Gyr]'),
        "err_age": safe_get('Star Age Error [Gyr]', None),
        "vsini_tag": None,  # not read from the catalog
        "Prot": safe_get('Stellar Rotation Period [d]'),
        "err_Prot": safe_get('Star Rotation Period Error [days]'),
        "log_Rhk": safe_get('logRhk'),
        "err_log_Rhk": safe_get('logRhk Error', None),
        "K_act": safe_get('Stellar Activity RV Amplitude [m/s]'),
        "catalogs": cadence_file,
        "eccentricity": safe_get('Eccentricity'),
        "periastron": safe_get('Periastron')
    }
    
    # empty cells (NaN) -> None
    for key, value in target_params.items():
        if isinstance(value, (float, int)) and np.isnan(value):
            target_params[key] = None
    # planetary mass: convert between Jupiter and Earth masses
    if target_params["pl_mass_Mj"] is None and target_params["pl_mass_Me"] is None:
        raise ValueError("Error: the configuration must provide at least 'pl_mass_Mj' or 'pl_mass_Me'.")
    if target_params["pl_mass_Mj"] is None:
        target_params["pl_mass_Mj"] = target_params["pl_mass_Me"] * M_earth.value / M_jup.value
    if target_params["pl_mass_Me"] is None:
        target_params["pl_mass_Me"] = target_params["pl_mass_Mj"] * M_jup.value / M_earth.value

    # planetary radius: convert between Jupiter and Earth radii
    if target_params["pl_radius_Rj"] is None and target_params["pl_radius_Re"] is not None:
        target_params["pl_radius_Rj"] = target_params["pl_radius_Re"] * R_earth.value / R_jup.value
    if target_params["pl_radius_Re"] is None and target_params["pl_radius_Rj"] is not None:
        target_params["pl_radius_Re"] = target_params["pl_radius_Rj"] * R_jup.value / R_earth.value
    if target_params["pl_radius_Rj"] is None and target_params["pl_radius_Re"] is  None:
        target_params["pl_radius_Rj"] = None
        target_params["pl_radius_Re"] = None 

    # Default uncertainties
    if target_params.get('err_pl_period') is None:
        if target_params.get('pl_period') is not None: target_params['err_pl_period'] = target_params['pl_period']/10
    
    if target_params.get('err_star_mass') is None: 
        if target_params.get('star_mass') is not None: target_params['err_star_mass'] = target_params["star_mass"]*0.1

    if target_params.get('err_star_radius') is None: 
        if target_params.get('star_radius') is None: target_params["err_star_radius"] = None
        else: target_params["err_star_radius"] = target_params['star_radius']*0.1
    
    if target_params.get('err_star_temp') is None: 
        target_params['err_star_temp'] = 100
    
    if target_params.get('err_age') is None:
        if target_params.get("age") is not None: target_params['err_age'] = target_params["age"]*0.1
        else: target_params['err_age'] = None 

    if target_params.get('err_Prot') is None:
        if target_params.get('Prot') is None: target_params['err_Prot'] = None
        else: target_params['err_Prot'] = target_params['Prot'] * 0.1

    if target_params.get('err_log_Rhk') is None:
        if target_params.get('log_Rhk') is None: 
            target_params['err_log_Rhk'] = None
        else: 
            target_params['err_log_Rhk'] = 0.05  # typical uncertainty: 0.05 dex

    # Circular orbit unless the eccentricity is given
    if target_params.get('eccentricity') is None:
        target_params['eccentricity'] = 0.0

    if target_params.get('periastron') is None:
        target_params['periastron'] = 0.0

    # V magnitude: V, otherwise Gaia, otherwise 10
    resolve_magnitude(target_params)
    
    return target_params


def load_single_target_params(config):
    """
    Parameters of the single target defined in single_target_params of the configuration.
    Missing planetary masses and radii are converted from the other unit; missing
    uncertainties get default values.
    """
    
    # parameters from config.yaml
    target_params = {
        "targetname": config['targetname'],
        "star_mass": config['star_mass'],
        "err_star_mass": config.get('err_star_mass'), 
        "star_radius": config['star_radius'],
        "err_star_radius": config['err_star_radius'],
        "star_temp": config['star_temp'],
        "err_star_temp": config['err_star_temp'],
        "sp": config['sp'],
        "v_mag": config['v_mag'],
        "gaia_mag": config['gaia_mag'],
        "pl_mass_Mj": config.get('pl_mass_Mj'),
        "pl_mass_Me": config.get('pl_mass_Me'),
        "pl_radius_Rj": config.get('pl_radius_Rj'),
        "pl_radius_Re": config.get('pl_radius_Re'),
        "pl_period": config['pl_period'],
        "err_pl_period": config['err_pl_period'],
        "transit_midpoint": config['transit_midpoint'],
        "age": config.get('age'),
        "err_age": config.get('err_age'),
        "vsini_tag": config.get('vsini_tag'),
        "Prot": config.get('Prot'),
        "err_Prot": config.get('err_Prot'),
        "log_Rhk": config.get('log_Rhk'),
        "err_log_Rhk": config.get('err_log_Rhk'),
        "K_act": config.get('K_act'),
        "eccentricity": config.get('eccentricity'),
        "periastron": config.get('periastron')
    }
    
    # empty cells (NaN) -> None
    for key, value in target_params.items():
        if isinstance(value, (float, int)) and np.isnan(value):
            target_params[key] = None

    # planetary mass: convert between Jupiter and Earth masses
    if target_params["pl_mass_Mj"] is None and target_params["pl_mass_Me"] is None:
        raise ValueError("Error: the configuration must provide at least 'pl_mass_Mj' or 'pl_mass_Me'.")
    if target_params["pl_mass_Mj"] is None:
        target_params["pl_mass_Mj"] = target_params["pl_mass_Me"] * M_earth.value / M_jup.value
    if target_params["pl_mass_Me"] is None:
        target_params["pl_mass_Me"] = target_params["pl_mass_Mj"] * M_jup.value / M_earth.value

    # planetary radius: convert between Jupiter and Earth radii
    if target_params["pl_radius_Rj"] is None and target_params["pl_radius_Re"] is None:
       pass
    if target_params["pl_radius_Rj"] is None and target_params["pl_radius_Re"] is not None:
        target_params["pl_radius_Rj"] = target_params["pl_radius_Re"] * R_earth.value / R_jup.value
    if target_params["pl_radius_Re"] is None and target_params["pl_radius_Rj"] is not None:
        target_params["pl_radius_Re"] = target_params["pl_radius_Rj"] * R_jup.value / R_earth.value

    # Default uncertainties
    if target_params.get('err_pl_period') is None:
        if target_params.get('pl_period') is not None: target_params['err_pl_period'] = target_params['pl_period']/10
    
    if target_params.get('err_star_mass') is None: 
        if target_params.get('star_mass') is not None: target_params['err_star_mass'] = target_params["star_mass"]*0.1

    if target_params.get('err_star_radius') is None: 
        if target_params.get('star_radius') is None: target_params["err_star_radius"] = None
        else: target_params["err_star_radius"] = target_params['star_radius']*0.1
    
    if target_params.get('err_star_temp') is None: 
        target_params['err_star_temp'] = 100
    
    if target_params.get('err_age') is None:
        if target_params.get("age") is not None: target_params['err_age'] = target_params["age"]*0.1
        else: target_params['err_age'] = None 

    if target_params.get('err_Prot') is None:
        if target_params.get('Prot') is None: target_params['err_Prot'] = None
        else: target_params['err_Prot'] = target_params['Prot'] * 0.1

    if target_params.get('err_log_Rhk') is None:
        if target_params.get('log_Rhk') is None: 
            target_params['err_log_Rhk'] = None
        else: 
            target_params['err_log_Rhk'] = 0.05  # typical uncertainty: 0.05 dex


    # Circular orbit unless the eccentricity is given
    if target_params.get('eccentricity') is None:
        target_params['eccentricity'] = 0.0
    if target_params.get('periastron') is None:
        target_params['periastron'] = 0.0
    
     # V magnitude: V, otherwise Gaia, otherwise 10
    resolve_magnitude(target_params)


    
    return target_params

    


def print_params(target_params):
  print('****************************')
  print('         PARAMETERS         ')
  print('****************************')
  print('Target: ' ,target_params['targetname'])
  print('Stellar Mass: ', target_params['star_mass'])     #Unit: Msun    (mandatory)
  print('Stellar Mass Error: ', target_params['err_star_mass']) 
  print('Stellar Radius: ', target_params['star_radius']) #Unit: Rsun    (optional)
  print('Stellar Radius Error: ', target_params['err_star_radius'])
  print('Stellar Temperature: ', target_params['star_temp']) #Unit: Kelvin  (mandatory)
  print('Stellar Temperature Error: ', target_params['err_star_temp'])
  print('Spectral Type:' , target_params['sp'])  #Type: string, e.g. G0 V (at least the first letter)   (mandatory)

  
  print('Vmag: ', target_params['v_mag'])           
  print('Gaia mag: ', target_params['gaia_mag'])                

  print('Planetary Mass [Mj]: ', target_params['pl_mass_Mj'])     #Unit: Mass of Jupipter
  print('Planetary Mass [Me]: ', target_params['pl_mass_Me'])     #Unit: Mass of Earth

  print('Planetary Period [d]: ', target_params['pl_period'])                    #Unit: days (mandatory)
  print('Transit MidTime [JD-2450000]: ', target_params['transit_midpoint'])             #Unit: JD-2450000 (mandatory)

  print('Planetary Radius [Rj]: ', target_params['pl_radius_Rj'])            #Unit: Radius of Jupiter 
  print('Planetary Radius [Re]: ', target_params['pl_radius_Re'])            #Unit: Radius of Earth

  print('Stellar Age [Gyr]: ', target_params['age'])            #Unit: Gyr
  print('Stellar Age Error: ' , target_params['err_age'])
  print('vsini tag: ', target_params['vsini_tag'])      #(Available options: 'slow', 'fast', 'intermediate')

  print('Stellar Rotational Period [d]:', target_params['Prot'])           #Unit: days
  print('Stellar Rotation Period Error: ', target_params['err_Prot'])
  print('log R\'HK: ', target_params['log_Rhk'])        
  print('Stellar Activity RV Amplitude [m/s]: ', target_params['K_act'])  #Unit: m/s
  print('****************************')


def import_file_sampling(InDir=None,inputfile_sampling=None,inputfile_sampling_usp=None):
  if InDir == None:
    InDir = './utility_files/'
  if inputfile_sampling==None:
    inputfile_sampling='harps.rv'
  if inputfile_sampling_usp==None:
    inputfile_sampling_usp='time_usp.txt'
    
  # Harps rv
  '''
  Temporal sampling from a real set of observations. 300 time stamps spread over six seasons.
  This temporal sampling is presented in Benatti et al. 2020 A&A 639, 50.
  To use a different time sampling, please provide a different input file with the preferred Julian days.  
  '''
  time_sampling = np.genfromtxt(InDir + inputfile_sampling,usecols=0)

  # check the array dimension
  if time_sampling.ndim == 1:
    time_sampling = time_sampling
  else:
    time_sampling = time_sampling[:, 0]


  '''
  An alternative temporal sampling useful to recover ultra short period planets. It is considered when the planet period exceed the Nyquist frequency of the default temporal sampling.
  Even this secondary sampling is obtained within GAPS, Nardiello et al. 2022 for the USP planet TOI-1807, however, the timestamps are selected in order to produce a higher Nyquist frequency.
  '''
  time_sampling_usp = np.genfromtxt(InDir + inputfile_sampling_usp,usecols=0)

  return time_sampling, time_sampling_usp

def import_nsts_noise(InDir=None, inputfile_nsts=None):
  if InDir == None:
    InDir = './utility_files/'
  if inputfile_nsts == None:
    inputfile_nsts = 'NSTS_noise_new.txt'

  '''
  Read the file with the RV photon noise as a function of magnitude, exposure time and spectral type.
  Data obtained with NSTS (ETC mode) for airmass = 1.5 and seeing = 1"
  '''
  nsts_noise = np.genfromtxt(InDir + inputfile_nsts)
  tmp_mag, t_exp, errRV_F, errRV_G, errRV_K, errRV_M = nsts_noise[:,0], nsts_noise[:,1], nsts_noise[:,2], nsts_noise[:,3], nsts_noise[:,4], nsts_noise[:,5]

  return tmp_mag, t_exp, errRV_F, errRV_G, errRV_K, errRV_M

def import_max_age_from_mass(InDir=None, inputfile_mass_age=None):
  if InDir == None:
    InDir = './utility_files/'
  if inputfile_mass_age == None:
    inputfile_mass_age = 'max_age_from_mass.dat'

  '''
  Read the file containig the upper limit of the star's age according with its mass value 
  '''
  max_age_from_mass = np.genfromtxt(InDir + inputfile_mass_age)
  mass_tmp, max_age = max_age_from_mass[:,0], max_age_from_mass[:,1]

  return mass_tmp, max_age

def import_file_activity(InDir=None,inputfilename=None):
    if InDir == None:
        InDir = './utility_files/'
    if inputfilename==None:
        inputfilename='PLOT_POINTS_DATA_KG-KP-M-MP-YO-OC-KP-SO.csv'

    '''
    Read the file containing the distribution of the error as function of logR'hk for different star spectral type
    '''
    df_activity = pd.read_csv(InDir+inputfilename)

    return df_activity


def load_config(filename="config.yaml"):
    """
    Read the YAML configuration and return a flat dictionary with the general variables
    and the parameters of the active mode (catalog or single target).

    Parameters:
    filename (str): YAML file to read.

    Returns:
    dict: active configuration variables.
    """
    with open(filename, "r") as file:
        config = yaml.safe_load(file)
    
    # active parameters
    active_config = {}

    if "general_variables" in config:
        active_config.update(config["general_variables"])

    active_config["read_catalog"] = config.get("read_catalog", False)
    active_config["single_target"] = config.get("single_target", False)

    # exactly one mode must be active
    if active_config["read_catalog"] == active_config["single_target"]:
        raise ValueError("Error: in the configuration, exactly one of 'read_catalog' and 'single_target' must be True.")
    
    if config.get("read_catalog", False):
        active_config.update(config.get("catalog", {}))  # catalog mode: catalog file name

    if config.get("single_target", False):
        params = config.get("single_target_params", {})

        # at least one planetary mass is required
        pl_mass_Mj = params.get("pl_mass_Mj")
        pl_mass_Me = params.get("pl_mass_Me")

        if pl_mass_Mj is None and pl_mass_Me is None:
            raise ValueError("Error: the configuration must provide at least 'pl_mass_Mj' or 'pl_mass_Me'.")
        
        # only pl_mass_Mj: compute pl_mass_Me
        if pl_mass_Mj is not None and pl_mass_Me is None:
            pl_mass_Me = pl_mass_Mj * M_jup.value/M_earth.value

        # only pl_mass_Me: compute pl_mass_Mj
        if pl_mass_Me is not None and pl_mass_Mj is None:
            pl_mass_Mj = pl_mass_Me * M_earth.value/M_jup.value

        # add the computed values to the active parameters
        active_config.update(params)
        active_config["pl_mass_Mj"] = pl_mass_Mj
        active_config["pl_mass_Me"] = pl_mass_Me

        # planetary radius (optional)
        pl_radius_Rj = params.get("pl_radius_Rj")
        pl_radius_Re = params.get("pl_radius_Re")

        # only pl_radius_Rj: compute pl_radius_Re
        if pl_radius_Rj is not None and pl_radius_Re is None:
            pl_radius_Re = pl_radius_Rj * R_jup.value/R_earth.value

        # only pl_mass_Me: compute pl_mass_Mj
        if pl_radius_Re is not None and pl_radius_Rj is None:
            pl_radius_Rj = pl_radius_Re * R_earth.value/R_jup.value

        # add the computed values to the active parameters
        active_config.update(params)
        active_config["pl_radius_Rj"] = pl_radius_Rj
        active_config["pl_radius_Re"] = pl_radius_Re


    return active_config


def find_seasons_automatic(time_series, min_points=7, n_cuts=10, min_gap_days=50):
    """
    Split a time series into observing seasons at its n_cuts largest gaps, keeping only
    the cuts that:
    - leave at least `min_points` points in every segment;
    - correspond to a gap longer than `min_gap_days`.

    Args:
        time_series (array): 1D array of times (sorted or not).
        min_points (int): minimum number of points of a segment.
        n_cuts (int): number of largest gaps tried at the start.
        min_gap_days (float): minimum gap for a valid cut [d].

    Returns:
        array: times (xm) at which the series is split.
    """

    x = np.sort(np.array(time_series))
    diff = np.diff(x)
    # indices of the n_cuts largest gaps
    top_diff_indices = np.argsort(diff)[-n_cuts:]
    top_diff_indices = np.sort(top_diff_indices) + 1
    
    xm = []
    cut_gaps = []
    valid_indices = []
    for idx, ind in enumerate(top_diff_indices):
        if ind - 1 >= 0 and ind + 1 < len(x):
            left = x[ind - 1]
            right = x[ind + 1]
            gap = right - left
            if gap >= min_gap_days:
                xm.append((left + right) / 2)
                cut_gaps.append(gap)
                valid_indices.append(ind)

    xm = np.array(xm)
    cut_gaps = np.array(cut_gaps)

    # iteratively remove the cuts that leave too few points in a segment
    while True:
        segments = np.split(x, np.searchsorted(x, xm))
        segment_lengths = [len(seg) for seg in segments]

        invalid_segments = [i for i, l in enumerate(segment_lengths) if l < min_points]
        if not invalid_segments:
            break

        # cuts to remove
        candidate_cuts = []
        for i in invalid_segments:
            if i == 0 or i > len(xm):
                continue
            cut_idx = i - 1
            candidate_cuts.append(cut_idx)

        if not candidate_cuts:
            break

        # remove the cut with the smallest gap
        best_cut = min(candidate_cuts, key=lambda i: cut_gaps[i])
        xm = np.delete(xm, best_cut)
        cut_gaps = np.delete(cut_gaps, best_cut)

    return xm

def define_seasons(time_sampling, xm):
    """
    Split time_sampling into seasons at the times xm.

    Args:
        time_sampling (array): times to split.
        xm (array): times separating the seasons.

    Returns:
        list of arrays: indices of the observations of each season.
    """

    time_sampling = np.array(time_sampling)
    i=0
    while (xm.size==0):
      i+=1
      xm = find_seasons_automatic(time_sampling, min_points=7-i, n_cuts=10, min_gap_days=50-i*5)

    xm = np.sort(np.array(xm))  # sorted
    seasons = []

    # first season: before the first xm
    s0 = np.where(time_sampling <= xm[0])[0]
    seasons.append(s0)

    # intermediate seasons
    for i in range(len(xm) - 1):
        si = np.where((time_sampling > xm[i]) & (time_sampling <= xm[i + 1]))[0]
        seasons.append(si)

    # last season: after the last xm
    s_last = np.where(time_sampling > xm[-1])[0]
    seasons.append(s_last)

    return seasons
