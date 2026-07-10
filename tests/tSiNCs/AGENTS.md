# Si / Diamond Nanocrystal Vibrations & FTIR

## Purpose

Own the **end-to-end nanocrystal vibration workflow** in FireCore: structure generation (JS/Python), MMFF relax/Hessian/spectrum (Python + C++), QM reference spectra, NPZ fixtures, and C++/Vispy viewers. Canonical home for all nanocrystal **orchestration scripts** and test outputs (`OUT_*`).

## Ownership

- Node CLIs: `nanocrystals.mjs`, ensemble/atlas configs, small symmetric batch generator
- Python test drivers: `run_vib_spectra.py`, `crosscheck_nanocrystal_generators.py`, NPZ bootstrap/export helpers
- FFfit scripts: `test_FFfit.py` (thin CLI wrapper), `test_fffit_hybrid.py` (hybrid fitting tests), `test_parity_py_cpp.py`, `test_parity_graph_cpp.py` — reusable logic lives in `pyBall/FFfit_utils.py` and `pyBall/FFfit_plots.py`
- Committed fixtures: `fixtures/si_1nm_passivation/`, `fixtures/npz_viewer/`, `fixtures/vibration_*`
- Viewer launchers: `run_cpp_mol_browser.sh`, `run_vispy_mol_browser.sh`, `test_cpp_npz_load.sh`
- Generated output trees under this folder (`OUT_small_nc/`, `OUT_nc_ensemble_v2/`, `OUT_nc_atlas/`) — not committed unless promoted to fixtures

## Local Contracts

- **Scripts resolve `REPO` as FireCore root** (`path.resolve(__dirname, '../..')`); **outputs default under `TEST_DIR`** (this folder).
- **Libraries stay outside:** `web/molgui_webgpu/Nanocrystals.js`, `web/common_js/npzIO.js`, `pyBall/nanocrystal_pipeline.py`, `pyBall/FTIR.py`, `pyBall/FFfit_utils.py`, `pyBall/FFfit_plots.py`.
- **Production spectrum:** dense `eigh` + mode histogram via `pyBall/nanocrystal_pipeline.py` (full MMFF Hessian FD at relax geometry). Linearized `03_topology.npz` is built in ensemble runs but **not** used for production spectra yet (Path B deferred).
- **Legacy `scripts/` copies** of nanocrystal CLIs are deprecated; edit and run **`tests/tSiNCs/`** versions only.
- **CIF and MMFF `.dat` files** come from `cpp/common_resources/` (read-only from REPO).

## Work Guidance

### Entrypoints (from repo root)

```bash
node tests/tSiNCs/nanocrystals.mjs ensemble --config tests/tSiNCs/ensemble.example.json
node tests/tSiNCs/make_small_symmetric_nc.mjs
python3 tests/tSiNCs/crosscheck_nanocrystal_generators.py
python3 -m pyBall.nanocrystal_pipeline relax --init-npz …
```

### Subcommands (`nanocrystals.mjs`)

| Subcommand | Role |
|------------|------|
| `generate` | Single crystal mol2/xyz |
| `ensemble` | Batch gen → relax → topology → Hessian → spectrum → accumulate |
| `topology` | MMFFL linearized `03_topology.npz` |
| `audit` | Geometry presets G0–G5 |
| `nonbond` | Collision / nonbond group debug |
| `rings` | Ring detection SVG |

### Documentation triad (keep in sync)

| Doc | Role |
|-----|------|
| [`README.md`](README.md) | Folder index + run commands |
| [`doc/Topics/FTIR_Nanocrystals/README.md`](../../doc/Topics/FTIR_Nanocrystals/README.md) | Topic guides and progress logs |
| [`doc/topical_audit/Nanocrystal_Vibrations.md`](../../doc/topical_audit/Nanocrystal_Vibrations.md) | Cross-implementation audit |

Open work: [`ToDo_Nanocrystal.md`](ToDo_Nanocrystal.md).

## Verification

- Generator parity: `python3 crosscheck_nanocrystal_generators.py` (from this folder)
- Small NC batch: `node make_small_symmetric_nc.mjs` → `OUT_small_nc/index.md`
- NPZ gallery: `fixtures/si_1nm_passivation/` stages 01→05 per crystal
- Classical ladder: `tests/tMMFF/run.sh test_vibration_spectra.py`
- FFfit hybrid tests: `python3 -m pytest tests/tSiNCs/test_fffit_hybrid.py -v` (5 test cases)
- Graph + batch dihedral parity: `PYTHONPATH=.. python3 test_parity_graph_cpp.py` (14 tests)

## Child DOX Index

| Path | Role |
|------|------|
| [`fixtures/si_1nm_passivation/README.md`](fixtures/si_1nm_passivation/README.md) | Nine-crystal NPZ pipeline gallery |
| [`fixtures/npz_viewer/README.md`](fixtures/npz_viewer/README.md) | Minimal viewer smoke fixtures |
| [`fixtures/vibration_benchmarks/README.md`](fixtures/vibration_benchmarks/README.md) | Benchmark NPZ structures |
