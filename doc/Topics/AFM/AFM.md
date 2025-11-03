# AFM Simulation Package Documentation

This document provides a comprehensive overview of the Atomic Force Microscopy (AFM) simulation package, detailing its structure, workflow, and usage.

## 1. Code Map: Key Files and Modules

The package is a multi-language system combining Python for high-level logic, C++/OpenCL for GPU-accelerated computations, and Fortran for the core DFT calculations.

| Component | Path | Language | Role |
|---|---|---|---|
| **High-Level Jobs** | `/pyBall/DFT/jobs.py` | Python | Orchestrates the entire AFM simulation workflow, from DFT to final image generation. |
| **Mid-Level Logic** | `/pyBall/DFT/high_level.py` | Python | Provides functions for specific tasks like coefficient conversion and convolution. |
| **Python-C++ Bridge**| `/pyBall/DFT/oclfft.py` | Python | `ctypes` wrapper that loads the C++ shared library and exposes its functions to Python. |
| **C++/OpenCL Host** | `/cpp/libs_OCL/OCL_GridFF.cpp`| C++ | Main C++ library file. Manages OpenCL environment and defines the C functions called by Python. |
| **DFT/OCL Interface** | `/cpp/common/OpenCL/OCL_DFT.h` | C++ | Header defining the `OCL_DFT` class, which encapsulates DFT-specific OpenCL tasks like density projection. |
| **Primary Kernels** | `/cpp/common_resources/cl/myprog.cl` | OpenCL C | **Contains the critical density projection kernels** (`projectOrbDenToGrid_texture`, `projectDenmatToGrid`, etc.). |
| **Force Field Kernels**| `/cpp/common_resources/cl/GridFF.cl` | OpenCL C | Contains kernels for generating force fields from potentials (e.g., `make_MorseFF`). |
| **Relaxation Kernels**| `/cpp/common_resources/cl/relax.cl` | OpenCL C | Contains kernels for probe particle relaxation (`relaxStrokesTilted`) and post-processing (`convolveZ`). |
| **DFT Engine** | `/fortran/MAIN/libFireCore.f90`| Fortran | Core Fireball DFT library for self-consistent wavefunction and density calculations. |
| **Tests** | `/tests/tDFT/`, `/tests/tDFT_CO/`, `/tests/tDFT_pentacene/` | Python | Example scripts demonstrating how to use the package. |

---

## 2. Workflow Review

The simulation follows a five-step process to generate AFM frequency shift maps. Each step is detailed below with the full function call chain.

### Step 1: DFT Calculation (Fireball)

-   **Goal:** Calculate self-consistent wavefunctions and electron density matrix for the sample and tip molecules.
-   **Concept:** A quantum chemistry calculation is performed using the Fireball DFT code to determine the electronic structure of the molecules.
-   **Call Chain:**
    1.  **Python Layer:** The process is initiated by calling the function `prepare_fireball()` in the file `/home/prokop/git/FireCore/pyBall/DFT/jobs.py`.
    2.  **C++ Bridge:** This Python function uses a `ctypes` wrapper around the `libFireCore.so` shared library, which is compiled from the Fortran source code in `/home/prokop/git/FireCore/fortran/`.
    3.  **Fortran Engine:** The `init`, `assembleH`, and `solveH` functions within the Fortran library (e.g., in `/home/prokop/git/FireCore/fortran/MAIN/libFireCore.f90`) are executed to run the SCF calculation and produce the molecular orbital coefficients (`wfcoef`).

### Step 2: Density Projection to Grid

-   **Goal:** Project the electron density, derived from the DFT calculation, onto a 3D real-space grid.
-   **Concept:** The total electron density `rho(r)` is constructed by summing the contributions of the atomic basis functions `phi_i` weighted by the density matrix `P_ij`, according to the formula `rho(r) = sum_{i,j} P_{ij} * phi_i(r) * phi_j(r)`.
-   **Call Chain:**
    1.  **Python Layer:** The function `project_dens_GPU()` in `/home/prokop/git/FireCore/pyBall/DFT/jobs.py` is called, passing the wavefunction coefficients.
    2.  **C++ Bridge:** This calls the wrapper function `oclfft.projectAtomsDens()` in `/home/prokop/git/FireCore/pyBall/DFT/oclfft.py`.
    3.  **C++/OCL Host:** The wrapper calls the C-exported function `projectAtomsDens()` located in `/home/prokop/git/FireCore/cpp/libs_OCL/OCL_GridFF.cpp`. This function is a method of the `OCL_DFT` class (defined in `/home/prokop/git/FireCore/cpp/common/OpenCL/OCL_DFT.h`), which manages the OpenCL kernel arguments and execution.
    4.  **OpenCL Kernel:** The C++ host code enqueues the kernel **`projectOrbDenToGrid_texture`** from `/home/prokop/git/FireCore/cpp/common_resources/cl/myprog.cl`. This kernel iterates over each grid point and calculates the density by summing the square of the molecular orbital values. Alternative kernels in the same file support related projections:
        - **`projectDenmatToGrid` / `_simp`** consume a full density matrix and mirror the traversal used on the GPU for better numerical symmetry.
        - **`projectAtomDenToGrid_texture`** accumulates per-atom reference densities (used when constructing difference densities with `projectAtomsDens0`).
        - The shared mechanics include looping over orbitals (`iorb0` to `iorb1`), sampling basis functions from the 2D texture (`imgIn`), and blending results with the `acumCoef` weights.
        - The dominant path (`projectOrbDenToGrid_texture`) squares the molecular orbital amplitude (`wf.x*wf.x`) and adds it to `rho(r)` for that grid point.
        - Typical density projections pass `acumCoef=[0.0, 2.0]` so that each occupied orbital contributes two electrons (see `project_dens_GPU()` in `jobs.py`). Difference-density runs toggle `bDen0diff`, which first records neutral atom densities through `projectAtomsDens0()` and then combines them via `acumCoef=[1.0,-1.0]`.

    Additional GPU projection modes are exposed for specialized workflows:
    - **Density-matrix projection (`project_denmat_GPU`)** in `jobs.py` switches to `oclfft.projectDenmat()`, which drives the `projectDenmatToGrid` kernel and can reuse the precomputed density matrix directly (skipping `convCoefs`).
    - **Reference density accumulation (`projectAtomsDens0`)** uploads atomic templates via `oclfft.projectAtomsDens0()`. This calls `projectAtomDenToGrid_texture` to build neutral or background densities that can later be combined with molecular results through `acumCoef` weights.

### Step 3: Potential & Energy Calculation

-   **Goal:** Calculate the full interaction potential (Electrostatic, Pauli Repulsion, Van der Waals) on the 3D grid.

#### A) Electrostatic Interaction

-   **Concept:** The electrostatic potential `V(r)` is calculated from the charge density `rho(r)` by solving the Poisson equation, `nabla^2 V = -rho/epsilon_0`. This is done efficiently in Fourier space where the equation becomes an element-wise multiplication. The final interaction energy is a convolution of the tip's charge distribution with the sample's potential.
-   **Call Chain (Hartree Potential):**
    1.  **Python Layer:** The `poisson()` function in `/home/prokop/git/FireCore/pyBall/DFT/high_level.py` is called.
    2.  **C++ Bridge:** This calls `oclfft.poisson()` in `/home/prokop/git/FireCore/pyBall/DFT/oclfft.py`.
    3.  **C++/OCL Host:** This calls the `poisson()` function in `/home/prokop/git/FireCore/cpp/libs_OCL/OCL_GridFF.cpp`. This function first performs an FFT on the density grid, then calls the `poissonW` kernel, and finally performs an inverse FFT.
    4.  **OpenCL Kernel:** The **`poissonW`** kernel registered from `/home/prokop/git/FireCore/cpp/common_resources/cl/myprog.cl` performs the element-wise multiplication in Fourier space (a legacy variant also exists in `GridFF.cl`, but the host-side bindings load the version in `myprog.cl`).

#### B) Pauli Repulsion and Van der Waals (VdW) Interaction

There are two primary models for handling these short-range interactions:

**Method 1: Full Analytical Potential (Lennard-Jones or Morse)**

-   **Concept:** Uses a classical analytical potential like Lennard-Jones or Morse, which has both a long-range attractive part (VdW) and a short-range repulsive part (Pauli). `E(r) = E0 * [exp(-2a(r-r0)) - 2*exp(-a(r-r0))]`.
-   **Call Chain:**
    1.  **Python Layer:** A high-level function in `jobs.py` or a test script calls `oclfft.evalLJC_QZs()`.
    2.  **C++ Bridge:** The wrapper in `oclfft.py` calls the C-exported function `evalLJC_QZs()` in `/home/prokop/git/FireCore/cpp/libs_OCL/OCL_GridFF.cpp`.
    3.  **C++/OCL Host:** This C++ function enqueues the OpenCL kernel.
    4.  **OpenCL Kernel:** The **`evalLJC_QZs`** kernel in `/home/prokop/git/FireCore/cpp/common_resources/cl/relax.cl` is executed. It iterates over all grid points and calculates the contribution from each atom based on the supplied Lennard-Jones/Morse parameters, generating grids for both the repulsive (`E_Paul`) and attractive (`E_Lond`) components. (The `make_MorseFF` kernel in `GridFF.cl` provides similar physics but is currently accessed through the PyOpenCL path rather than the C++ exports.)

**Method 2: Hybrid DFT-D (Density-based Pauli + Analytical VdW)**

-   **Concept:** This is a more accurate hybrid method. The Pauli repulsion is calculated from the overlap of the tip and sample electron densities, `E_Pauli ~ integral( (rho_tip * rho_sample)^alpha )`. The attractive VdW part, however, is *still required* and is calculated separately using the attractive part of an analytical potential.
-   **Call Chain (Pauli Repulsion from Density):**
    1.  **Python Layer:** The `convolve()` function in `/home/prokop/git/FireCore/pyBall/DFT/high_level.py` is called with the tip and sample density grids as input.
    2.  **C++ Bridge:** This calls `oclfft.convolve()` in `/home/prokop/git/FireCore/pyBall/DFT/oclfft.py`.
    3.  **C++/OCL Host:** This calls the `convolution()` function in `/home/prokop/git/FireCore/cpp/libs_OCL/OCL_GridFF.cpp`. This function orchestrates a forward FFT on both density grids, calls a kernel to multiply them in Fourier space, and then performs a backward FFT.
    4.  **OpenCL Kernel:** The element-wise multiplication is performed by the **`mul`** kernel in `/home/prokop/git/FireCore/cpp/common_resources/cl/myprog.cl`.
-   **Call Chain (VdW Attraction):**
    *   This is identical to Method 1, but only the attractive part (`E_Lond`) of the output is used. The same **`evalLJC_QZs`** kernel provides the attractive grid, and the host discards the repulsive component when building the total interaction.

### Step 4: Probe Particle Relaxation

-   **Goal:** Find the relaxed 3D position of the AFM tip’s apex atom in the combined potential field from Step 3.
-   **Concept:** The probe particle is modeled as a point mass attached by a spring to the tip anchor. It moves iteratively according to the net force from the pre-calculated potential grid until the force is minimized.
-   **Call Chain:**
    1.  **Python Layer:** A function like `relaxStrokesTilted()` in `/home/prokop/git/FireCore/pyBall/DFT/oclfft.py` is called.
    2.  **C++ Bridge:** This calls the C-exported function `relaxStrokesTilted()` in `/home/prokop/git/FireCore/cpp/libs_OCL/OCL_GridFF.cpp`.
    3.  **C++/OCL Host:** The C++ function sets up the arguments (including the 3D force field texture) and enqueues the relaxation kernel.
    4.  **OpenCL Kernel:** The **`relaxStrokesTilted`** kernel in `/home/prokop/git/FireCore/cpp/common_resources/cl/relax.cl` is executed. It performs the iterative relaxation using the FIRE algorithm (`update_FIRE` function within the kernel) to find the equilibrium position.

### Step 5: AFM Frequency Shift Calculation

-   **Goal:** Convert the calculated forces into a frequency shift (`df`) map, which represents the final AFM image.
-   **Concept:** The frequency shift is calculated from the force-vs-distance curve `F(z)` using the Giessibl formula, which involves a weighted integral (convolution) of the force over the oscillation amplitude of the cantilever.
-   **Call Chain:**
    1.  **Python Layer:** A high-level function would orchestrate this, likely calling a wrapper for the convolution.
    2.  **C++ Bridge:** A `ctypes` wrapper would call a C++ function to perform the convolution.
    3.  **C++/OCL Host:** The C++ function would enqueue the convolution kernel.
    4.  **OpenCL Kernel:** The **`convolveZ`** kernel in `/home/prokop/git/FireCore/cpp/common_resources/cl/relax.cl` is designed for this. It takes the array of forces along Z (`Fin`) and a set of weights (`weighs`) and performs a 1D convolution to produce the final frequency shift values (`Fout`).

---

## 2.1 Supporting Assets and Data Formats

Before running the workflow, ensure the following assets and conventions are in place:

- **Basis sets (`Fdata/basis/`):** `loadWfBasis()` searches the Fireball basis directory for `.wf{1,2}` tiles. The default sampling radius (`RcutSamp`) and per-element radial cutoffs (`Rcuts`) determine how much of each orbital is tabulated. Missing data triggers a runtime error when `projectAtomsDens` tries to sample the textures.
- **OpenCL sources:** `init()` and `initPP()` expect the `cl_src_dir` argument to point at `/cpp/common_resources/cl/`. The shared library loads `myprog.cl`, `GridFF.cl`, and `relax.cl` relative to this directory.
- **FFT buffers:** Density grids are stored as complex single-precision arrays (`float2`), matching the interfaces of `poisson()`, `convolve()`, and `gradient()`. Upload/download helpers automatically convert between double precision NumPy arrays and these GPU buffers (`upload_d`, `download`).
- **Grid descriptors:** Grid shape is defined by `pos0`, `dA`, `dB`, and `dC` vectors (stored as `float4`). Keep these consistent between `setGridShape` and the Python metadata (`setGridShape_dCell`) so that interpolation and relaxation kernels sample the same lattice.
- **Probe/force textures:** `evalLJC_QZs` outputs into `float4` buffers where `.xyz` stores force components and `.w` stores the potential energy. Relaxation and visualization code rely on this layout.

---

### 3. Tutorial: How to Run the Tests

The test scripts in `/tests/` are the best way to understand and run the package. The tests for pentacene and CO are particularly relevant.

### Prerequisites

Use the repository-provided helper scripts so that all dependencies are rebuilt with consistent flags. For AFM workflows this typically means:

```bash
# From the project root directory /home/prokop/git/FireCore/
cd tests/tDFT_pentacene/
./run.sh       # rebuilds required libraries and launches the example
```

The `run.sh` scripts take care of compiling (Fortran, C++, OpenCL) and exporting the necessary environment variables. Avoid invoking `make` manually unless you are working on low-level build system changes.

### Running a Test (Example: Pentacene)

1.  **Navigate to the test directory & read the script:**
    ```bash
    cd /home/prokop/git/FireCore/tests/tDFT_pentacene/
    less run.py
    ```
    The `run.py` script demonstrates how the Python side composes the workflow (density projection → potentials → relaxation → analysis).

2.  **Launch the example:**
    ```bash
    ./run.sh             # preferred entry point; wraps the command below
    # python run.py      # underlying command executed by run.sh
    ```

3.  **Other scenarios:** Similar harnesses exist in `tests/tDFT_CO/`, `tests/tDFT/`, and AFM-specific directories; always use their local `run.sh` wrappers for reproducible builds.

### What the Test Script Does (Conceptual Example)

A typical test script (`run_test.py` or similar) for a full AFM simulation would perform the following actions:

```python
# Conceptual example of a test script

import numpy as np
from pyBall.DFT import jobs, oclfft

# 1. Define the sample and tip geometry
pentacene_atoms, pentacene_types = load_xyz('pentacene.xyz')
co_atoms, co_types = load_xyz('CO.xyz')

# 2. Calculate sample density grid using Fireball + OpenCL
# This runs steps 1 and 2 of the workflow
jobs.projectDens(
    atomType=pentacene_types,
    atomPos=pentacene_atoms,
    ngrid=(128, 64, 192),
    dcell=[0.15, 0.15, 0.15],
    bSCF=True,
    saveName="pentacene_dens"
)

# 3. Calculate tip density grid
jobs.projectDens(
    atomType=co_types,
    atomPos=co_atoms,
    ngrid=(64, 64, 64),
    dcell=[0.15, 0.15, 0.15],
    bSCF=True,
    saveName="co_dens"
)

# 4. Calculate Potentials (Hartree, Pauli, VdW) - Step 3
# Load sample density, calculate potential
oclfft.loadFromBin("pentacene_dens.bin", 0)
oclfft.poisson(iA=0, iOut=1, dcell=[0.15, 0.15, 0.15]) # Hartree potential in buffer 1
oclfft.saveToBin("pentacene_Vhartree.bin", 1)

# ... similar steps for Pauli and VdW potentials ...
# The C++ code might combine these steps into a single call

# 5. Define scan grid and run relaxation - Step 4
scan_dim = (50, 50)
scan_window = ((-5.0, -5.0), (5.0, 5.0))
oclfft.makeStartPointGrid(scan_dim[0], scan_dim[1], p0=[...], da=[...], db=[...])

# Load all potentials into textures/buffers
# ...

# Relax the probe particle at all scan points
oclfft.relaxStrokesTilted(iBuffOut=4, nz=30, dtip=-0.1)

# 6. Calculate frequency shift and save results - Step 5
# The forces are in buffer 4. Now calculate df.
# This might be another kernel call or a CPU-based calculation.
# ...

print("Simulation finished. Output files are in the current directory.")
```

By studying the scripts in the test directories, you can see the exact function calls and parameters used to perform each step of the simulation.

---

## Appendix A: Host Exports and Kernel Mapping

The table below summarizes the main C exports provided by `libOCL_GridFF.so`, the Python wrappers that invoke them, and the OpenCL kernels they ultimately launch. Use it as a quick reference when extending the workflow.

| Python Wrapper (`pyBall/DFT/oclfft.py`) | C Export (`OCL_GridFF.cpp`) | Primary Kernel(s) | Notes |
|---|---|---|---|
| `projectAtoms()` | `projectAtoms` | `projectOrbDenToGrid_texture` (MO accumulation) | Projects a single orbital amplitude; mainly debugging. |
| `projectAtomsDens()` | `projectAtomsDens` | `projectOrbDenToGrid_texture` | Default density pipeline using MO coefficients and `acumCoef`. |
| `projectAtomsDens0()` | `projectAtomsDens0` | `projectAtomDenToGrid_texture` | Builds reference/neutral densities for difference calculations. |
| `projectDenmat()` | `projectDenmat` | `projectDenmatToGrid`, `projectDenmatToGrid_simp` | Consumes density matrix without converting to MO coefficients. |
| `poisson()` | `poisson` | `poissonW` (from `myprog.cl`) | FFT → kernel → inverse FFT on the selected buffer. |
| `convolve()` | `convolve` | `mul` | Performs element-wise multiplication in Fourier space. |
| `gradient()` | `gradient` | `gradient` | Computes ∇V in real space after inverse FFT. |
| `evalLJC_QZs()` | `evalLJC_QZs` | `evalLJC_QZs` | Produces Pauli/VdW grids stored in float4 buffers. |
| `evalLJC_QZs_toImg()` | `evalLJC_QZs_toImg` | `evalLJC_QZs_toImg` | Same physics, writes to a 3D image for texture sampling. |
| `relaxStrokesTilted()` | `relaxStrokesTilted` | `relaxStrokesTilted` | FIRE-based probe relaxation over the scan grid. |
| `getFEinStrokes()` | `getFEinStrokes` | `relaxStrokesTilted` (single step) | Samples forces along predefined strokes without relaxation. |
| `convolveZ()` (planned) | – (kernel consumed directly) | `convolveZ` | Currently launched from C++ utilities; add wrapper if df needs to run from Python. |

When porting functionality to PyOpenCL (see `AFM_migration_plan.md`), this mapping shows which kernels are already reusable and which orchestration steps must be replicated in Python.
