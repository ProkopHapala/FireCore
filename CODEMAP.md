# CODEMAP.md

# FireCore Repository Structure

FireCore is a multi-component computational chemistry and physics repository containing:

- DFT/DFTB electronic structure methods
- Classical molecular mechanics and MD
- GPU/OpenCL/CUDA acceleration
- QM/MM coupling
- Force-field fitting tools
- Visualization and GUI applications
- Numerical method experiments and prototypes

The repository is NOT a single build target. Different subsystems have different workflows and build procedures.

# Main Components

- `/cpp/` : High-performance C/C++ implementation.
   -  `/cpp/common/` : Core numerical algorithms, math, forcefields, molecular mechanics, neighbor lists, grid operations.
   -  `/cpp/libs/` : Reusable compiled libraries.
   -  `/cpp/apps/` : SDL-based standalone visualization and simulation applications.
   -  `/cpp/apps_OCL/` : OpenCL-accelerated applications.
   -  `/cpp/apps_CUDA/` : CUDA-accelerated applications.
   -  `/cpp/common_resources/` : Shared molecular data, forcefield parameters, basis sets, test molecules.
- `/fortran/`: Original Fireball DFT/DFTB implementation.
    - `/fortran/MAIN/` : SCF drivers
    - `/fortran/MODULES/` : global modules and data structures
    - `/fortran/ASSEMBLERS/` : Hamiltonian assembly
    - `/fortran/INTERACTIONS/` : interaction integrals
    - `/fortran/GRID/` : FFT/grid operations
    - `/fortran/NEIGHBORS/` : neighbor construction
- `/fortran2/` contains reorganized/refactored work-in-progress code.

- `/pyBall/` : Python bindings and utilities
    - `FireCore.py` — Fireball interface
    - `MMFF.py` — molecular mechanics
    - `AtomicSystem.py` — structure manipulation
    - `Forces.py` — force interfaces
    - `plotUtils.py` — plotting helpers
    - `buildUtils.py` — rebuild/recompile helpers
- `/pyBall/OCL/` : Standalone pyOpenCL implementations independent from C++ bindings.
- `/tests/` : Primary reference for implemented functionality and expected workflows.
    - `tests/Fireball/`
    - `tests/tMMFF*/`
    - `tests/tMolGUIapp*/`
    - `tests/tFitREQ*/`
    - `tests/tDFT*/`
    - `tests/tEFF*/`
    - `tests/pyFireball/`
- `/doc/` - Technical notes, derivations, experiments, and developer documentation. Includes:
    - Markdown docs
    - development notes
    - Julia experiments
    - Maxima derivations
    - Python prototypes

# Build & Execution

## Critical Rule

Do NOT compile manually using ad-hoc `make` commands.

Use provided:
- `run.sh`
- `make.sh`
- project test scripts

These scripts:
- configure paths,
- rebuild dependencies,
- launch correct binaries,
- manage runtime resources.

# Typical Workflow

## Run molecular mechanics tests

```bash
cd tests/tMMFF
./run.sh