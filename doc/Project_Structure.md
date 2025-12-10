# FireCore Project Structure and Development Guide

## Project Overview

FireCore is a comprehensive computational chemistry and physics repository containing multiple subprojects, programs, and test scripts. It implements various numerical methods for:

- Density Functional Theory (DFT) — Fortran core, callable from C++/Python via library
- Classical Molecular Dynamics — C++ and OpenCL
- Force Field Development and Fitting — reactive and non‑reactive
- Quantum‑Classical Hybrid Methods (QM/MM) — DFT + classical FF
- Numerical Method Demonstrations — standalone Python scripts

Important: FireCore is not a single build target but a repository of interconnected subprojects with different compilation and execution workflows.

## Directory Structure and Components

### 1. Core C++ Implementation (`cpp/`)

- standalone visual applications (`cpp/apps/`)
  - SDL interactive molecular editors and simulators
  - GPU‑accelerated variants: `cpp/apps_OCL/` (OpenCL), `cpp/apps_CUDA/` (CUDA)
- libraries
  - `cpp/libs/Molecular/` — core molecular mechanics
  - `cpp/libs_SDL/` — SDL visualization
  - `cpp/libs_OCL/` — OpenCL accelerated libs
- core implementation
  - `cpp/common/` — algorithms, data structures, math
  - `cpp/common_SDL/` — GUI and visualization
  - `cpp/common_resources/` — shared resources (params, molecules)

### 2. Fortran DFT Implementation (`fortran/`, `fortran2/`)

- `fortran/` — current Fireball DFT implementation
- `fortran2/` — reorganized/refactored version (WIP)
- key subdirs in `fortran/`: `MAIN/`, `MODULES/`, `ASSEMBLERS/`, `INTERACTIONS/`, `GRID/`, `NEIGHBORS/`

### 3. Python Interface (`pyBall/`)

- core modules: `FireCore.py` (DFT), `MMFF.py`, `MMFFsp3.py`, `AtomicSystem.py`, `Forces.py`, `Forces_cpp.py`
- specialized: `pyBall/OCL/` (pyOpenCL), `pyBall/DFT/`, `pyBall/GUI/`
- utilities: `atomicUtils.py`, `plotUtils.py`, `buildUtils.py`

### 4. OpenCL Implementation (`pyBall/OCL/`)

- standalone pyOpenCL implementations independent of C++ libs
- direct OpenCL kernels; alternative to C++ OpenCL bindings

### 5. Tests and Examples (`tests/`)

- Best place to see implemented functionality and usage
- categories: `tests/Fireball/`, `tests/tMMFF*/`, `tests/tMolGUIapp*/`, `tests/tFitREQ*/`, `tests/tDFT*/`, `tests/tEFF*/`, `tests/pyFireball/`, `tests/pySCF/`, `tests/dftb/`

### 6. Documentation (`doc/`)

- `doc/Markdown/` — technical docs (GUI, force fields, etc.)
- `doc/DevNotes/` — development notes and TODOs
- `doc/py/` — Python sketches and demonstrations
- `doc/Julia/` — Julia implementations of algorithms
- `doc/Maxima/` — derivations and symbolic computations

### 7. Build System

- `cpp/Build/` — CMake build directory for C++ components
- `Build/` or `build/` — Fortran build output (see `make.sh`)
- `make.sh` — main build script for Fortran components
- `cpp/CMakeLists.txt` — C++ CMake configuration

## Workflows

- C++/Apps: Prefer using `run.sh` scripts in `tests/` to build and run example targets; they set paths and arguments
- Fortran: `make.sh` in project root, then `tests/Fireball/*/run.sh`
- Python: set `PYTHONPATH` to the repo root, use `tests/*/run.sh` to run examples

## See Also

- Global coding rules: `doc/Global_Coding_Rules.md`
