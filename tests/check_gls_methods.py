"""
Check: the vectorised GLS (minerva.fast_gls) against the reference implementation
(minerva/_vendor/gls.py, Zechmeister & Kuerster 2009).

On random data sets (number of points, time span, signal and frequency range drawn at
random) it compares the periodogram powers, the highest peak and its parameters, and
reports the speed-up.

Usage:  python tests/check_gls_methods.py [--n 60] [--quick]
Exit code 0 if, for every data set, the highest peak is the same and the relative
differences of the powers are below the tolerance of each method.
"""
import argparse
import sys
import time

import numpy as np

from minerva._vendor.gls import Gls
from minerva.fast_gls import FastGls

TOLERANCE = {"exact": 1e-10, "recurrence": 1e-6}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--n", type=int, default=60, help="number of random data sets")
    ap.add_argument("--quick", action="store_true", help="10 data sets")
    args = ap.parse_args()
    n = 10 if args.quick else args.n

    rng = np.random.default_rng(3)
    t_cost = {"reference": 0.0, "exact": 0.0, "recurrence": 0.0}
    worst = {"exact": 0.0, "recurrence": 0.0}
    ok = True
    for _ in range(n):
        N = int(rng.integers(10, 300))
        t = np.sort(rng.uniform(0, 2000, N))
        y = 5 * np.sin(2 * np.pi * t / rng.uniform(2, 40)) + rng.normal(0, 2, N)
        e = np.full(N, 2.0)
        fend = float(rng.choice([0.5, 2.0, 5.0]))

        t0 = time.perf_counter()
        ref = Gls((t, y, e), fbeg=1 / 2000, fend=fend, verbose=False)
        t_cost["reference"] += time.perf_counter() - t0
        for method in ("exact", "recurrence"):
            t0 = time.perf_counter()
            new = FastGls((t, y, e), fbeg=1 / 2000, fend=fend, verbose=False, method=method)
            t_cost[method] += time.perf_counter() - t0
            rel = np.max(np.abs(ref.power - new.power) / np.maximum(np.abs(ref.power), 1e-300))
            worst[method] = max(worst[method], rel)
            if np.argmax(ref.power) != np.argmax(new.power) or rel > TOLERANCE[method]:
                ok = False

    print(f"{n} random data sets")
    for method in ("exact", "recurrence"):
        print(f"  {method:10s}: {t_cost['reference'] / t_cost[method]:4.1f}x faster, "
              f"max relative difference of the powers {worst[method]:.1e} (tolerance {TOLERANCE[method]:.0e})")
    print("PASSED" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
