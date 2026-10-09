# MINERVA examples

Each folder contains a complete configuration (`config.yaml`), meant to be copied and
customised; the catalog example also contains its catalog and cadence. All targets are
fictitious. Run from the repository root:

```bash
python scripts/run_minerva.py examples/01_catalog/config.yaml
python scripts/run_minerva.py examples/02_single_target/config.yaml --quick   # 5 simulations
```

Results are written to `results/` inside the example folder
([output files](../docs/api/conventions.md#output-files)).

| Example | Mode | Science case | What it shows |
|---|---|---|---|
| `01_catalog` | catalog, 4 targets | [survey planning](../docs/cases/survey_planning.md) | own cadence per target (`real_cadence`), parallel targets (`n_workers: auto`): quiet M dwarf, active K dwarf, inactive F star, eccentric G dwarf |
| `02_single_target` | single target | [single target](../docs/cases/single_target.md), [stellar activity](../docs/cases/stellar_activity.md) | activity derived from the stellar age; debug figures of the first simulations |
| `03_known_activity` | single target | [stellar activity](../docs/cases/stellar_activity.md#fixing-what-is-known) | age, P_rot and log R'HK fixed to known values; K_act from the K-dwarf relation |
| `04_bootstrap` | single target | [parameter uncertainties](../docs/cases/parameter_uncertainties.md) | stellar parameters resampled within their uncertainties |
| `05_eccentric_orbit` | single target | [eccentric orbits](../docs/cases/eccentric_orbits.md) | e = 0.5 injected, circular fit (compare with `02_single_target`) |
| `06_usp_planet` | single target | [cadences](../docs/cases/cadences.md) | P = 0.55 d: high-cadence USP sampling, 160 epochs |

To use your own targets: copy a folder, put your catalog and cadence files in it (or set
`paths.user_dir`), edit the configuration ([parameter reference](../docs/api/parameters.md)).
For science runs use at least 100 simulations per target (`n_iterations`).
