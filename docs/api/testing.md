# Tests and validation

[← API reference](README.md)

Run from the repository root after `pip install -e .`.

| Script | What it checks | Reference |
|---|---|---|
| `tests/regression/` | the complete pipeline on the example catalog (4 targets × 10 simulations): every result file must be reproduced exactly | outputs of the current version (`tests/regression/reference/`) |
| `tests/check_parallel.py` | a parallel catalog run gives the same results as a sequential run; two identical runs give identical results | sequential run |
| `tests/check_gls_methods.py` | vectorised GLS (`exact`, `recurrence`) vs the original implementation on random data sets: same highest peak, relative differences of the powers, speed | `minerva/_vendor/gls.py` |
| `scripts/check_same_code.py` | two versions of a file contain the same code (comments, docstrings and messages excluded) | previous version |

```bash
# regression test
python scripts/run_minerva.py tests/regression/config_regression.yaml
python tests/regression/compare.py Output_regression

# parallel and reproducibility; GLS implementations
python tests/check_parallel.py --workers 4
python tests/check_gls_methods.py
```

The examples are integration tests too:

```bash
for f in examples/*/config.yaml; do python scripts/run_minerva.py "$f" --quick || break; done
```

### Regression test

The targets of the example catalog exercise the main branches of the code: quiet M dwarf,
active K dwarf (GAPS/ArMS activity relation), hot F star without convective envelope (no
activity), eccentric orbit, planet-only and simultaneous fits, detections, false positives
and non-detections. The reference is produced with `tests/regression/config_regression.yaml`.

A change that is not meant to change the results must reproduce the reference exactly.
A change that is meant to change them (a new relation, a bug fix) requires a new
reference: regenerate it and explain why in the commit message.

```bash
cp Output_regression/results_*.txt tests/regression/reference/
```

Results depend at the last digits on the versions of numpy and scipy: generate the
reference on the machine where MINERVA is run.

### Validation against the literature

The scientific validation (number of RV observations of published planets, HARPS-N
cadences) is described in the paper (Sect. 3.1, Table 2); the literature targets and
cadences are not distributed with the code.
