"""
Vectorised Generalised Lomb-Scargle periodogram.

FastGls is the GLS of Zechmeister & Kuerster (2009) shipped in minerva/_vendor/gls.py,
with the same equations (8)-(15) and the same results, but the sums over the data are
computed for a whole block of frequencies at once (matrix-vector products) instead of
one frequency at a time in a Python loop. It is several times faster.

Two methods:
- "exact": cosines and sines computed directly, as in the reference implementation; the sums
  are added up in a different order, so the powers differ in the last digits (~1e-13 relative);
- "recurrence": cos/sin obtained by trigonometric recurrences (the `fast` option of the
  reference code), refreshed every 1000 frequencies; much faster, ~1e-9 relative differences.
"""
import functools

import numpy as np

from ._vendor.gls import Gls as _ReferenceGls

# Frequencies processed together: memory ~ BLOCK x N_data x 8 bytes per array.
# A multiple of the refresh interval of the recurrences.
_BLOCK = 4000
_REFRESH = 1000


class FastGls(_ReferenceGls):
    """Drop-in replacement of Gls with a vectorised periodogram computation."""

    method = "exact"   # "exact" or "recurrence"

    def __init__(self, *args, method=None, **kwargs):
        if method is not None:
            if method not in ("exact", "recurrence"):
                raise ValueError(f"Unknown FastGls method: {method}")
            self.method = method
        super().__init__(*args, **kwargs)

    def _calcPeriodogram(self):
        if self.fast:   # trigonometric recurrences: keep the reference implementation
            return super()._calcPeriodogram()

        if self.e_y is None:
            w = np.ones(self.N)
        else:
            w = 1 / (self.e_y * self.e_y)
        self.wsum = w.sum()
        w /= self.wsum

        self._Y = np.dot(w, self.y)       # Eq. (7)
        wy = self.y - self._Y             # subtract the weighted mean
        self._YY = np.dot(w, wy**2)       # Eq. (10), weighted variance with offset only
        wy *= w                           # attach the weights

        C, S, YC, YS, CC, CS = np.zeros((6, self.nf))
        omega = 2. * np.pi * self.f
        for a in range(0, self.nf, _BLOCK):
            b = min(a + _BLOCK, self.nf)
            if self.method == "recurrence":
                cosx, sinx = self._recurrence_block(omega, a, b)
            else:
                x = omega[a:b, None] * self.th[None, :]
                cosx = np.cos(x)
                sinx = np.sin(x)
            C[a:b] = cosx @ w                 # Eq. (8)
            S[a:b] = sinx @ w                 # Eq. (9)
            YC[a:b] = cosx @ wy               # Eq. (11)
            YS[a:b] = sinx @ wy               # Eq. (12)
            CC[a:b] = (cosx * cosx) @ w       # Eq. (13)
            CS[a:b] = (cosx * sinx) @ w       # Eq. (15)

        SS = 1. - CC
        if not self.ls:
            CC -= C * C            # Eq. (13)
            SS -= S * S            # Eq. (14)
            CS -= C * S            # Eq. (15)
        D = CC*SS - CS*CS          # Eq. (6)

        self._a = (YC*SS-YS*CS) / D
        self._b = (YS*CC-YC*CS) / D
        self._off = -self._a*C - self._b*S

        # power
        self.p = (SS*YC*YC + CC*YS*YS - 2.*CS*YC*YS) / (self._YY*D)   # Eq. (5) in ZK09

    def _recurrence_block(self, omega, a, b):
        """cos and sin for frequencies a..b-1 by trigonometric recurrence: exp(i w_k t) =
        exp(i w_k0 t) * exp(i dw t)^(k-k0), restarted every _REFRESH frequencies."""
        eid = np.exp(2j * np.pi * self.fstep * self.th)
        out = np.empty((b - a, self.N), dtype=complex)
        for r in range(a, b, _REFRESH):
            q = min(r + _REFRESH, b)
            block = np.empty((q - r, self.N), dtype=complex)
            block[0] = np.exp(1j * omega[r] * self.th)
            block[1:] = eid
            out[r - a:q - a] = np.cumprod(block, axis=0)
        return out.real, out.imag


def make_gls(method="reference"):
    """GLS class for the method chosen in the configuration (analysis.gls_method):
    "reference" (the vendored implementation, as in the paper), "exact" or "recurrence"."""
    if method == "reference":
        return _ReferenceGls
    if method in ("exact", "recurrence"):
        return functools.partial(FastGls, method=method)
    raise ValueError(f"Unknown gls_method: {method!r} (use 'reference', 'exact' or 'recurrence')")
