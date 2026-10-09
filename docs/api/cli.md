# Command-line tools

[← API reference](README.md)

## scripts/run_minerva.py

```bash
python scripts/run_minerva.py <config.yaml> [--quick] [-w N | --workers N]
```

| Option | Description |
|---|---|
| `config.yaml` | the configuration ([parameter reference](parameters.md)) |
| `--quick` | 5 simulations per target: a fast check of the installation (smoke test, not science) |
| `-w, --workers N` | catalog targets in N parallel processes (`auto`: all cores); overrides `run_options.n_workers` |

Outputs in `paths.output_dir` ([conventions](conventions.md#output-files)). The same run
from Python:

```python
from minerva.runner import run
run("examples/01_catalog/config.yaml", overrides={"n_iterations": 5, "n_workers": "auto"})
```

## Examples

`examples/NN_<case>/config.yaml` — see the [examples overview](../../examples/README.md):

```bash
python scripts/run_minerva.py examples/01_catalog/config.yaml
```

Results in `examples/NN_<case>/results/`.

## Utilities

| Script | Description |
|---|---|
| `scripts/compress_yapsi_tracks.py <tracks folder> data/YaPSI_Tracks` | copies, gzip-compressed, only the YaPSI tracks used by MINERVA (90 of 111; ~90 MB) |
| `scripts/make_logo.py [output.svg] [Montserrat font]` | regenerates the vector logo `docs/images/minerva_logo.svg` (needs fontTools and the Montserrat font) |
| `scripts/make_flowchart.py [output.svg]` | regenerates the flowchart `docs/images/minerva_flowchart.svg` (edit the script when the code changes) |
| `scripts/check_same_code.py old.py new.py` | compares the syntax trees of two versions of a file, ignoring comments, docstrings and messages: for documentation-only changes |
| `tests/regression/compare.py <output folder>` | compares a run with the regression reference ([tests](testing.md)) |
