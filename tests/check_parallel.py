"""
Check: a catalog run with parallel workers gives the same results as a sequential run,
and two identical runs give identical results (reproducibility).

Runs the regression configuration three times (sequential, sequential again, parallel)
with few simulations and compares every result file.

Usage:  python tests/check_parallel.py [--workers 4] [--n-iterations 3]
"""
import argparse
import filecmp
import shutil
import sys
import tempfile
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "tests" / "regression" / "config_regression.yaml"


def run(tmp, name, n_workers, n_iterations):
    from minerva.runner import run as minerva_run
    cfg = yaml.safe_load(open(CONFIG))
    cfg["paths"] = {"data_dir": str(ROOT / "data"), "user_dir": str(ROOT / "examples" / "01_catalog"),
                    "output_dir": str(tmp / name)}
    cfg["run_options"]["n_workers"] = n_workers
    cfg["general_variables"]["n_iterations"] = n_iterations
    path = tmp / f"{name}.yaml"
    yaml.safe_dump(cfg, open(path, "w"), sort_keys=False)
    minerva_run(str(path))
    return sorted((tmp / name).glob("results_*.txt"))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--n-iterations", type=int, default=3)
    args = ap.parse_args()

    tmp = Path(tempfile.mkdtemp(prefix="minerva_check_"))
    try:
        a = run(tmp, "sequential", 1, args.n_iterations)
        b = run(tmp, "sequential_again", 1, args.n_iterations)
        c = run(tmp, "parallel", args.workers, args.n_iterations)
        ok = True
        for label, other in [("repeated sequential run", b), (f"parallel run ({args.workers} workers)", c)]:
            same = [filecmp.cmp(x, y, shallow=False) for x, y in zip(a, other)] if len(a) == len(other) else [False]
            print(f"{label:32s}: {sum(same)}/{len(a)} result files identical")
            ok &= all(same)
        print("PASSED" if ok else "FAILED")
        return 0 if ok else 1
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
