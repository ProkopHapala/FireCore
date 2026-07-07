# FireCore Tests and Examples

This directory contains examples and test scripts that demonstrate FireCore functionality. **Each subfolder has its own `README.md`** with entry points, useful scripts, and obsolete files to ignore.

> **Reorganization in progress (planning only):** see [`REORGANIZATION.plan.md`](REORGANIZATION.plan.md). Major split: `tVibrations`, `tSurfaces`, `tRigidBody` (new siblings); FDBM → `tAFM`; shrink `tMMFF`. Not implemented yet.

## Policies (agreed direction)

### Naming: `test_*` is often a misnomer

Most `test_*.py` files are **runnable demos** (plots, scans, exploratory workflows), not pytest/automated tests. Going forward:

- **`demo_*`** — exploration / workflow scripts
- **`check_*`** / **`verify_*`** — validation with explicit pass/fail
- Do not add new `test_*` unless it is a real automated check

### Output hygiene (bonding rule)

**Scripts must not write results next to source code.** All runtime output goes under `out/<script_stem>/` (or `fixtures/` for committed reference data).

- We **avoid `.gitignore` on output trees** that agents need to read for debugging.
- Large ephemeral runs live under `out/` with per-script subdirs; README documents what is safe to delete.
- See [`REORGANIZATION.plan.md`](REORGANIZATION.plan.md) for full layout and current violations (`tMMFF/result_Trajectroy_Opt/`, etc.).

## Quick start

Most directories with C++/GPU libs use `run.sh` (compile + run):

```bash
cd tests/tMMFF && ./run.sh
cd tests/tUFF && python3 test_parity_suite.py
cd tests/Fireball/t01_H2 && ../../../build/fireball.x
```

Python-only dirs: run scripts directly (`tests/dftb/`, `tests/pyFireball/`, etc.).

## Directory index

Each row links to the folder README. **Status** reflects current maintenance level.

### Molecular mechanics and force fields

| Directory | Purpose | Entry | Status |
|-----------|---------|-------|--------|
| [tMMFF/](tMMFF/README.md) | MMFF core demos, assembly (shrinking) | `./run.sh` | **Messy** — split planned; see REORGANIZATION.plan.md |
| *`tVibrations/`* | *planned:* general phonons, Hessians, solvers | — | **Not created** — from tMMFF; not tSiNCs |
| *`tSurfaces/`* | *planned:* GridFF, Ewald2D, folded electrostatics | — | **Not created** — absorbs tEwald2D |
| *`tRigidBody/`* | *planned:* rigid scans & dynamics on surfaces | — | **Not created** — from tMMFF |
| [tMMFFmulti/](tMMFFmulti/README.md) | Multi-replica MMFF on GPU | `./run.sh` | Active |
| [tMMFFsp3/](tMMFFsp3/README.md) | MMFFsp3 profiles, H-bond scans | `./run.sh` | Active |
| [tUFF/](tUFF/README.md) | UFF/MMFF CPU vs OCL vs PyOCL parity | `test_parity_suite.py` | **Active** — primary parity suite |
| [tCUDA/](tCUDA/README.md) | CUDA vs OpenCL MMFF | `./run.sh` | Active (needs CUDA) |
| [tFitFF/](tFitFF/README.md) | Fireball-based FF parameter fitting | `./run.sh` | Active |
| [tFitREQ/](tFitREQ/README.md) | REQH non-bonded fitting (Rvdw, Evdw, Q, H-bond) | `./run.sh` | Active |
| [tEwald2D/](tEwald2D/README.md) | 2D Ewald electrostatics | `./run.sh` | Active — **→ merge into `tSurfaces/`** |
| [tLattice2D/](tLattice2D/README.md) | 2D lattice vector matching | `./run.sh` | Active |
| [tFF2D/](tFF2D/README.md) | 2D force-field topology builder | `./run.sh` | Active |
| [NonBondSampling/](NonBondSampling/README.md) | vdW surface sampling | `sample.py` | Exploratory |

### GUI applications

| Directory | Purpose | Entry | Status |
|-----------|---------|-------|--------|
| [tMolGUIapp/](tMolGUIapp/README.md) | Single-molecule editor | `./run.sh` | Active — `run.sh` is mostly commented recipes |
| [tMolGUIapp_multi/](tMolGUIapp_multi/README.md) | Multi-replica OCL editor | `./run.sh` | Active |
| [tMolGUIapp_QMMM/](tMolGUIapp_QMMM/README.md) | QM/MM GUI | `./run.sh` | Active |
| [tMolGUIapp_QMMM_multi/](tMolGUIapp_QMMM_multi/README.md) | Multi-replica QM/MM | `./run.sh` | Active |

### DFT and quantum methods

| Directory | Purpose | Entry | Status |
|-----------|---------|-------|--------|
| [Fireball/](Fireball/README.md) | Fortran `fireball.x` molecule tests (24 cases) | per-subdir `run.sh` | **Mixed** — path bugs; `t02_CH4` broken |
| [pyFireball/](pyFireball/README.md) | Python FireCore, STM, ribbon scans, OCL parity | `./run.sh` | **Active** — ignore 19 `copy`/`bak` files |
| [dftb/](dftb/README.md) | DFTB+ API and OpenCL waveplot parity | `test_python_api.py` | **Active** — legacy HBsmall batch scripts obsolete |
| [tDFT/](tDFT/README.md) | GPU density projection | `./run.sh` | Active |
| [tDFT_CO/](tDFT_CO/README.md) | CO density projection | `./run.sh` | Active |
| [tDFT_pentacene/](tDFT_pentacene/README.md) | Pentacene density / CO convolution | `./run.sh` | Active |
| [pySCF/](pySCF/README.md) | PySCF experiments | `try_pyscf.py` | Exploratory |
| [pyocl_dft/](pyocl_dft/README.md) | OpenCL DFT density milestones | `test_firecore_data.py` | Active |
| [tPsi4resp/](tPsi4resp/README.md) | Psi4 RESP fitting | `psi4resp.py` | Active (needs conda `p4env`) |

### eFF, Kekulé, QM/MM

| Directory | Purpose | Entry | Status |
|-----------|---------|-------|--------|
| [tEFF/](tEFF/README.md) | eFF/RARFF CPU vs GPU parity | `./run.sh` | Active |
| [tEFFapp/](tEFFapp/README.md) | Standalone `EFFapp.x` | `./run.sh` | Active |
| [tKekule/](tKekule/README.md) | C++ Kekulé optimizer | `./run.sh` | Legacy — use tKekuleExplorer |
| [tKekuleExplorer/](tKekuleExplorer/README.md) | KekuleBackend tests | per-script CLI | **Active** |
| [tQMMM_diacetylene/](tQMMM_diacetylene/README.md) | QM/MM diacetylene | `./run.sh` | Legacy paths |

### AFM, nanocrystals, specialized

| Directory | Purpose | Entry | Status |
|-----------|---------|-------|--------|
| [tAFM/](tAFM/README.md) | AFM + FDBM (consolidating from tMMFF) | `./run.sh` | **Active** — see `pyocl_fdbm/` |
| [tSiNCs/](tSiNCs/README.md) | Si/diamond nanocrystal FTIR hub (specialized) | `run_vib_spectra.py` | **Active** — links to `tVibrations` for general APIs |
| [tXRD/](tXRD/README.md) | XRD Debye scattering | `test_debye_histogram.py` | Active |
| [tAttach/](tAttach/README.md) | Molecular attachment / polymers | `attach_new3.py` | Active |
| [tIsing/](tIsing/README.md) | Hubbard/Ising OpenCL MC | `run_hubbard_cli.py` | Active |
| [tMQCA/](tMQCA/README.md) | Molecular quantum cellular automata | `test_mqca.py` | Active |
| [tLammpsTrj/](tLammpsTrj/README.md) | LAMMPS trajectory analysis | `run.py` | Legacy — needs local traj data |

### Solvers and utilities

| Directory | Purpose | Entry | Status |
|-----------|---------|-------|--------|
| [tQuadrature/](tQuadrature/README.md) | 3D quadrature rules | `./run.sh` | Active |
| [tSchroedinger1D/](tSchroedinger1D/README.md) | 1D Green's function solver | `./run.sh` | Active |
| [tSchroedinger2D/](tSchroedinger2D/README.md) | 2D Green's function solver | `./run.sh` | Active |
| [pyutils/](pyutils/README.md) | Ad-hoc utilities | per-script | Exploratory |
| [blender/](blender/README.md) | Blender molecular rendering | `p3_atoms_pbr.py` | Active |
| [tmp/](tmp/README.md) | Scratch experiments | — | **Ignore** |

## Common patterns

### `run.sh` recipe menus

Many GUI and MMFF `run.sh` files are **mostly commented** example invocations. The active command is usually the last uncommented `python3` or `./Binary` line. Use `./run.sh no` in GUI tests to skip recompile.

### Obsolete file naming

Across `tests/`, files with **`copy`**, **`bak`**, **`legacy`**, or **`old`** in the name are almost always stale duplicates. Prefer the canonical name without suffix. Examples:

- `tests/pyFireball/` — 19 obsolete copies
- `tests/tMMFF/` — `run_tipSpline_scan_bak.py`, `TipSplineOptimizer copy.py`, etc.
- `tests/dftb/` — `example_orbitals copy.py`

### Resource symlinks

Many dirs symlink `data/` or `common_resources/` → `../../cpp/common_resources`. Regenerate if missing:

```bash
ln -s ../../cpp/common_resources data
```

### Prerequisites

- **Fdata:** Fortran tests need external `Fdata_HCNO` (see [fireball-qmd.github.io](https://fireball-qmd.github.io))
- **GPU:** CUDA/OpenCL tests need drivers and hardware
- **Build:** compile from `cpp/` or `fortran/` before running; `run.sh` handles this where present

## Recommended smoke tests

| Goal | Command |
|------|---------|
| MMFF phonons | `cd tests/tMMFF && ./run.sh` |
| UFF/MMFF parity | `cd tests/tUFF && python3 test_parity_suite.py` |
| eFF parity | `cd tests/tEFF && ./run.sh` |
| Fireball Fortran | `cd tests/Fireball/t01_H2 && ../../../build/fireball.x` |
| Fireball Python | `cd tests/pyFireball && ./run.sh` |
| DFTB+ API | `cd tests/dftb && python3 test_python_api.py` |
| Nanocrystal FTIR | `cd tests/tSiNCs && python3 run_vib_spectra.py adamantane` |

## Adding new tests

1. Create `tests/t<Name>/` following naming convention (sibling folder, not nested under tMMFF)
2. Add `run.sh` if compilation is needed
3. Symlink `data` → `../../cpp/common_resources` when appropriate
4. **All outputs → `out/<script_stem>/`** — never beside `.py` files; golden refs → `fixtures/`
5. Use `demo_*` / `check_*` naming — not `test_*` unless automated
6. Write `README.md` with entry point, useful scripts, and what to ignore
7. Update this index and `REORGANIZATION.plan.md` if layout changes
