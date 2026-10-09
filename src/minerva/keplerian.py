"""
Keplerian radial velocities.

The functions below reproduce those of RadVel (Fulton et al. 2018, PASP 130, 044504;
https://github.com/California-Planet-Search/radvel, MIT License), which MINERVA used as
a dependency: they are the only ones MINERVA needs, and including them avoids the many
dependencies of RadVel (h5py, pytables, ultranest, jupyter, ...).

* Circular orbits (e = 0, every activity sinusoid and most planets): the same NumPy
  formula as RadVel, so the results are identical.
* Eccentric orbits: the Kepler solver is the NumPy version of RadVel (RadVel uses an
  equivalent solver written in C by default); results agree to ~1e-12 m/s.

RadVel is Copyright (c) 2016 California Planet Search and distributed under the MIT License.
"""
import numpy as np


def rv_drive(t, orbel):
    """Radial velocity of a Keplerian orbit.

    Args:
        t (array): times of observation [d]
        orbel (array): [per, tp, e, om, K]: period [d], time of periastron [d],
            eccentricity, argument of periastron of the star [rad], semi-amplitude [m/s]

    Returns:
        array: radial velocity [m/s]
    """
    per, tp, e, om, k = orbel

    # circular orbits
    if e == 0.0:
        m = 2 * np.pi * (((t - tp) / per) - np.floor((t - tp) / per))
        return k * np.cos(m + om)

    if per < 0:
        per = 1e-4
    if e < 0:
        e = 0
    if e > 0.99:
        e = 0.99

    nu = true_anomaly(t, tp, per, e)
    return k * (np.cos(nu + om) + e * np.cos(om))


def true_anomaly(t, tp, per, e):
    """True anomaly at the times t (Murray & Dermott, p. 27)."""
    m = 2 * np.pi * (((t - tp) / per) - np.floor((t - tp) / per))
    eccarr = np.zeros(t.size) + e
    e1 = kepler(m, eccarr)
    n1 = 1.0 + e
    n2 = 1.0 - e
    return 2.0 * np.arctan((n1 / n2)**0.5 * np.tan(e1 / 2.0))


def kepler(Marr, eccarr):
    """Solve Kepler's equation E - e sin E = M for the eccentric anomaly E (iterative,
    third-order corrections, convergence 1e-12)."""
    conv = 1.0e-12
    k = 0.85
    Earr = Marr + np.sign(np.sin(Marr)) * k * eccarr      # first guess
    fiarr = (Earr - eccarr * np.sin(Earr) - Marr)
    convd = np.where(np.abs(fiarr) > conv)[0]
    nd = len(convd)
    while nd > 0:
        ecc = eccarr[convd]
        E = Earr[convd]
        fi = fiarr[convd]
        fip = 1 - ecc * np.cos(E)
        fipp = ecc * np.sin(E)
        fippp = 1 - fip
        d1 = -fi / fip
        d2 = -fi / (fip + d1 * fipp / 2.0)
        d3 = -fi / (fip + d2 * fipp / 2.0 + d2 * d2 * fippp / 6.0)
        E = E + d3
        Earr[convd] = E
        fiarr = (Earr - eccarr * np.sin(Earr) - Marr)
        convd = np.abs(fiarr) > conv
        nd = np.sum(convd)
    if Earr.size > 1:
        return Earr
    return Earr[0]


def timetrans_to_timeperi(tc, per, ecc, omega):
    """Time of periastron passage from the time of transit.

    Args:
        tc (float): time of transit [d]
        per (float): period [d]
        ecc (float): eccentricity
        omega (float): argument of periastron [rad]
    """
    try:
        if ecc >= 1:
            return tc
    except ValueError:
        pass
    f = np.pi / 2 - omega
    ee = 2 * np.arctan(np.tan(f / 2) * np.sqrt((1 - ecc) / (1 + ecc)))   # eccentric anomaly
    return tc - per / (2 * np.pi) * (ee - ecc * np.sin(ee))
