"""
Prepare the YaPSI tracks for the repository.

Copies, gzip-compressed, only the tracks that MINERVA actually reads
(mixing length A1p91804 for M < 1 Msun, A1p82126 for M >= 1 Msun).
numpy.loadtxt reads .gz files directly, so the code needs no change.

Usage
-----
    python scripts/compress_yapsi_tracks.py <folder with the original .trk files> data/YaPSI_Tracks
"""

import glob
import gzip
import os
import shutil
import sys


def tracks_used_by_minerva(src):
    """Same file selection as the track loader in MINERVA."""
    files = sorted(glob.glob(os.path.join(src, "M0p*A1p91804*.trk")))
    for k in "1234":
        files += sorted(glob.glob(os.path.join(src, f"M{k}p*A1p82126*.trk")))
    return files


def main(src, dst):
    os.makedirs(dst, exist_ok=True)
    files = tracks_used_by_minerva(src)
    for f in files:
        out = os.path.join(dst, os.path.basename(f) + ".gz")
        with open(f, "rb") as a, gzip.open(out, "wb", compresslevel=9) as b:
            shutil.copyfileobj(a, b)
    total = sum(os.path.getsize(os.path.join(dst, os.path.basename(f) + ".gz")) for f in files)
    print(f"{len(files)} tracks written to {dst} ({total / 1e6:.1f} MB)")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(2)
    main(sys.argv[1], sys.argv[2])
