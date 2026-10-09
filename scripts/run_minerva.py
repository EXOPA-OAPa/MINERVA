"""
Run MINERVA.

Usage
-----
    python scripts/run_minerva.py <config.yaml> [--quick] [--workers N|auto]

    --quick       5 simulations per target: a fast check of the installation, not science
    --workers N   catalog targets in N parallel processes ("auto": all cores); overrides
                  run_options.n_workers
"""
import argparse

from minerva.runner import run
from minerva.utils import title


def main():
    parser = argparse.ArgumentParser(
        description="MINERVA: number of RV observations needed to measure a planetary mass.")
    parser.add_argument("config", help="path to the YAML configuration file")
    parser.add_argument("--quick", action="store_true",
                        help="5 simulations per target (smoke test, not science)")
    parser.add_argument("-w", "--workers", default=None,
                        help="parallel processes for catalog targets (integer or 'auto')")
    args = parser.parse_args()

    overrides = {}
    if args.quick:
        overrides["n_iterations"] = 5
    if args.workers is not None:
        overrides["n_workers"] = args.workers if args.workers == "auto" else int(args.workers)

    title()
    run(args.config, overrides)


if __name__ == "__main__":
    main()
