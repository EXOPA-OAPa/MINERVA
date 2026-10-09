"""Stellar evolutionary tracks (YaPSI) and periodogram FAP levels."""
import functools
import glob
import numpy as np
import re


@functools.lru_cache(maxsize=None)
def _track_index(InDir):
    """Sorted list of the YaPSI track files used by MINERVA and their masses (read once per folder)."""
    track = np.sort(glob.glob(InDir+'YaPSI_Tracks/M0p*A1p91804*'))
    track = np.append(track, np.sort(glob.glob(InDir+'YaPSI_Tracks/M1p*A1p82126*')))
    track = np.append(track, np.sort(glob.glob(InDir+'YaPSI_Tracks/M2p*A1p82126*')))
    track = np.append(track, np.sort(glob.glob(InDir+'YaPSI_Tracks/M3p*A1p82126*')))
    track = np.append(track, np.sort(glob.glob(InDir+'YaPSI_Tracks/M4p*A1p82126*')))

    numeri = []
    for file_name in track:
        match = re.search(r'M([0-9]+)p([0-9]+)', file_name)
        if match:
            numeri.append(int(match.group(1)) + int(match.group(2)) / 100)
    return tuple(track), np.asarray(numeri)


@functools.lru_cache(maxsize=None)
def _read_track(path):
    """Content of one track file, read (and decompressed) only once. Read-only array."""
    matrix = np.loadtxt(path, comments='#')
    matrix.flags.writeable = False
    return matrix


def _load_track(starmass, InDir):
    '''
    Load the YaPSI track used for the star: the track with the largest mass <= starmass
    (or the lowest-mass track available). Mixing length A1p91804 below 1 Msun, A1p82126 above.
    Shared by track_pro and ms_lifetime, so that tau and the maximum age come from the same track.
    Track files are read once and kept in memory (they are the same for every simulation).
    '''
    track, numeri = _track_index(InDir)

    if np.size(np.where(numeri < starmass)) == 0:
        return _read_track(track[0])
    liminf = np.where(numeri <= starmass)[0][-1]
    return _read_track(track[liminf])


def track_pro(starmass, age, InDir):
    '''
    Convective turnover time tau [d] from the YaPSI tracks (Spada et al. 2017), at the
    model closest in age. Tracks: http://www.astro.yale.edu/yapsi/
    '''
    matrix = _load_track(starmass, InDir)
    age_trk = matrix[:, 2] * 1e9
    tau_c   = matrix[:, 13]
    return tau_c[np.argmin(abs(age_trk - age))]


def ms_lifetime(starmass, InDir, xc_min=1e-3):
    '''
    Age (in years) at the end of the main sequence: first model of the track with
    central hydrogen X_c < xc_min (column 4 of the YaPSI tracks).
    If the track does not reach the TAMS, the maximum age of the track is returned.
    '''
    matrix = _load_track(starmass, InDir)
    age_trk = matrix[:, 2] * 1e9
    Xc      = matrix[:, 3]
    idx = np.where(Xc < xc_min)[0]
    return age_trk[idx[0]] if idx.size else age_trk.max()


def fap_levels(fap, time, fbeg, fend):
    '''GLS power corresponding to a given false alarm probability (analytic approximation).'''
    NN = np.size(time)
    tbase = np.max(time)-np.min(time)
    M = (fend-fbeg)*tbase
    Prob = 1.-(1-fap)**(1/M)
    return 1-Prob**(2/(NN-3))