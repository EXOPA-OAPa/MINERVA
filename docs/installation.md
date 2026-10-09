# Installation

[← Back to README](../README.md)

## 1. Python environment

MINERVA requires Python ≥ 3.9. A dedicated virtual environment is recommended:

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
python -m pip install --upgrade pip setuptools
```

or with conda:

```bash
conda create -n minerva python=3.12
conda activate minerva
```

## 2. Install MINERVA

From the repository root, in **editable** mode, so that `import minerva` points to
`src/minerva/` and every change is immediately visible:

```bash
pip install -e .
```

Dependencies (installed automatically): numpy, scipy, pandas, matplotlib, seaborn,
astropy, pyyaml, tqdm: all available as pre-built packages, nothing is compiled. The GLS periodogram of Zechmeister & Kürster (2009) is
included in `src/minerva/_vendor/`.

Check that the package is the one in the repository:

```bash
python -c "import minerva; print(minerva.__file__)"     # .../MINERVA/src/minerva/__init__.py
```

## 3. Reference data

`data/` contains the reference data used by every run: the HARPS-N photon-noise table,
the maximum stellar ages, the GAPS activity database, the default cadences and the
YaPSI evolutionary tracks (gzip-compressed, only the 90 tracks used). They are read from
`data/` by default; another folder can be set with `paths.data_dir`.

To rebuild the compressed tracks from the original YaPSI files
(http://www.astro.yale.edu/yapsi/):

```bash
python scripts/compress_yapsi_tracks.py <folder with the .trk files> data/YaPSI_Tracks
```

## 4. Check the installation

```bash
python scripts/run_minerva.py examples/01_catalog/config.yaml --quick
python scripts/run_minerva.py tests/regression/config_regression.yaml
python tests/regression/compare.py Output_regression     # "All 5 files reproduce the reference exactly."
```

## Troubleshooting

* **`ModuleNotFoundError: minerva`**: install with `pip install -e .` from the repository
  root, or run with `PYTHONPATH=src`.
* **`FileNotFoundError` for a catalog or cadence**: relative paths in the configuration
  are relative to the configuration file; catalog and cadence files are looked up in
  `paths.user_dir`, then in `paths.data_dir`.
* **Regression test not reproduced**: results depend at the last digits on the versions of
  numpy and scipy; regenerate the reference on your machine
  ([tests](api/testing.md#regression-test)).
* **Parallel runs on macOS or Windows**: `n_workers > 1` uses the `spawn` start method;
  run the script from the command line (not from an interactive session that cannot be
  re-imported).
