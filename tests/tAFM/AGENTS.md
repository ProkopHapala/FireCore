# AFM Simulation Tests

## Purpose

AFM (Atomic Force Microscopy) simulation tests: FDBM (Full Density Based Model) force-field pipeline, pyOpenCL density-based path, Morse/LJ reference path, and rigid-body tip dynamics.

## Ownership

- FDBM pipeline tests (DFTB backend, pySCF integration)
- Gradient kernel tests for tip-surface interaction
- PTCDA-on-CaF2 and single-atom validation cases
- `pyocl_fdbm/` subdirectory: modular FDBM pipeline with Pauli fitting and DFTB reference

## Local Contracts

- **Run from this directory** — scripts use relative paths to data and kernels.
- **Use `run.sh` / `make.sh` scripts** when available; never invoke `make` directly.
- **FDBM pipeline:** Two paths exist — Morse/LJ (classical) and density-based (DFTB/pySCF). Prefer density-based for accuracy; use Morse/LJ for speed.
- **Pauli fitting:** `fit_fdbm_pauli.py` fits Pauli repulsion from DFTB density; verify against `diagnostic_forcefield.py`.

## Work Guidance

### Core Tests
- `test_fdbm.py` — basic FDBM force-field test
- `test_gradient_kernel.py` / `test_gradient_simple.py` / `test_gradient_visual.py` — gradient evaluation parity
- `test_ptcda.py` — PTCDA molecule test case
- `test_single_atom.py` — single-atom reference
- `afm_morse_pbc.py` — Morse potential with PBC

### FDBM Pipeline (`pyocl_fdbm/`)
- `run_pyocl_fdbm.py` / `run_pyocl_fdbm_dftb.py` — full FDBM pipeline
- `run_pyocl_fdbm_dftb_pentacene.py` — pentacene-specific test
- `fit_fdbm_pauli.py` — fit Pauli repulsion parameters
- `run_fitted_afm.py` — run with fitted parameters
- `compare_fireball_dftb.py` — compare Fireball vs DFTB reference
- `test_full_pipeline.py` / `test_modular_pipeline.py` — end-to-end pipeline tests

## Verification

- Run `test_fdbm.py` or relevant pipeline script from this directory
- Gradient tests should show smooth, physically meaningful force fields
- Compare FDBM results against DFTB reference where available

## Child DOX Index

- `pyocl_fdbm/` — Modular FDBM pipeline with DFTB backend and Pauli fitting
