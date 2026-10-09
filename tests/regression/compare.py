"""
Regression check for MINERVA.

Compares the result files of a run against the reference outputs stored in
``tests/regression/reference/``. The reference was produced with the example
catalog (``examples/catalog_example.csv``), 10 iterations per target, debug
plots disabled.

Any refactoring must reproduce the reference *exactly* (same text, same digits).
A difference means that some result changed, even if the logic looks the same:
a common cause is a change in the order of random draws from ``rng``.

Usage
-----
    python tests/regression/compare.py path/to/Output

Exit code 0 if every reference file is reproduced exactly, 1 otherwise.
"""

import csv
import io
import sys
from pathlib import Path

REFERENCE_DIR = Path(__file__).parent / "reference"


def _parse(line):
    """Split one CSV-like result line, keeping quoted lists such as "[30, 35]" together."""
    return next(csv.reader(io.StringIO(line), skipinitialspace=True))


def compare_file(ref_path, new_path):
    """Return a list of human-readable differences between two result files."""
    if not new_path.exists():
        return [f"missing file: {new_path.name}"]

    ref_lines = ref_path.read_text().splitlines()
    new_lines = new_path.read_text().splitlines()
    if ref_lines == new_lines:
        return []

    diffs = []
    header = [h.strip() for h in ref_lines[0].split(",")] if ref_lines else []
    if len(ref_lines) != len(new_lines):
        diffs.append(f"different number of lines: {len(ref_lines)} (reference) vs {len(new_lines)}")

    for n, (a, b) in enumerate(zip(ref_lines, new_lines), start=1):
        if a == b:
            continue
        fa, fb = _parse(a), _parse(b)
        changed = [
            f"{header[k] if k < len(header) else f'col {k}'}: {x} -> {y}"
            for k, (x, y) in enumerate(zip(fa, fb)) if x != y
        ]
        diffs.append(f"line {n}: " + ("; ".join(changed) if changed else "formatting differs"))
        if len(diffs) >= 10:
            diffs.append("... (further differences not shown)")
            break
    return diffs


def main(output_dir):
    output_dir = Path(output_dir)
    ref_files = sorted(REFERENCE_DIR.glob("results_*.txt"))
    if not ref_files:
        print(f"No reference files found in {REFERENCE_DIR}")
        return 1

    n_fail = 0
    for ref in ref_files:
        diffs = compare_file(ref, output_dir / ref.name)
        if diffs:
            n_fail += 1
            print(f"[DIFFERENT] {ref.name}")
            for d in diffs:
                print(f"    {d}")
        else:
            print(f"[IDENTICAL] {ref.name}")

    print()
    if n_fail:
        print(f"{n_fail}/{len(ref_files)} files differ from the reference.")
        return 1
    print(f"All {len(ref_files)} files reproduce the reference exactly.")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(2)
    sys.exit(main(sys.argv[1]))
