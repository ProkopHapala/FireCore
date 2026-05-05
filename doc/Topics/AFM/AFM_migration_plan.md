AFM Density-Based Simulation: Classical Force Fields and Density Projection Workflows
https://windsurf.com/codemaps/7a2c010f-2075-4a19-893c-1b54ec8df995-fe86ab10a43f3d18

AFM PyOpenCL System: Morse/LJ Path and FDBM Density-Based Path
https://windsurf.com/codemaps/9bb4c2a5-0c38-4943-abe9-254cfdcc75af-fe86ab10a43f3d18


# AFM PyOpenCL Migration Plan (2025 Update)

This plan describes how to replace the `libOCL_GridFF.so` C/C++ host layer with a pure Python implementation that relies on NumPy for data preparation and PyOpenCL for GPU execution. It expands previous drafts by synchronizing with the verified documentation in @/doc/Topics/AFM/AFM.md and enumerating every C/C++ function, class, and data dependency that must be re-created in Python.

## 1. Scope and Goals

- **Goal:** deliver a Python package (working name `pyBall.pyocl_dft`) that exposes the same API currently provided by `pyBall.DFT.oclfft`, while driving the existing OpenCL kernels (`myprog.cl`, `GridFF.cl`, `relax.cl`) directly from Python.
- **Out of scope:** rewriting the OpenCL kernels themselves; re-deriving Fireball DFT routines; changing physics models.
- **Success criteria:**
  1. All Python tests under `tests/tDFT*` and AFM tutorials run using the new package without loading `libOCL_GridFF.so`.
  2. Density grids, potentials, relaxation trajectories, and frequency-shift outputs match C++ results within numerical tolerance (define per dataset in §6).
  3. Documentation and examples default to the Python host while keeping the legacy C++ path as fallback until feature parity is proven.

## 2. Summary of the Current Workflow

The AFM workflow has five computational stages (@/doc/Topics/AFM/AFM.md#37-112):

1. **Fireball DFT** (Fortran) produces wavefunction coefficients `wfcoef` and density matrices.
2. **Density projection** uses `projectAtomsDens`, `projectAtomsDens0`, or `projectDenmat` to populate FFT-aligned grids (`float2` buffers) via `projectOrbDenToGrid_texture`, `projectAtomDenToGrid_texture`, or `projectDenmatToGrid` kernels (@/doc/Topics/AFM/AFM.md#41-54).
3. **Potential assembly** performs FFT-based Hartree/Pauli/VdW operations (`poisson`, `convolve`, `evalLJC_QZs`) (@/doc/Topics/AFM/AFM.md#55-90).
4. **Probe relaxation** runs `relaxStrokesTilted` to locate equilibrium tip positions (@/doc/Topics/AFM/AFM.md#88-100).
5. **Frequency-shift post-processing** applies `convolveZ` weights to relaxed force traces (@/doc/Topics/AFM/AFM.md#103-111).

Each stage is currently orchestrated in C++ via `OCL_DFT` and `OCL_PP` classes, exported in `OCL_GridFF.cpp` and invoked through `pyBall.DFT.oclfft.py` (summarized in Appendix A of AFM.md @/doc/Topics/AFM/AFM.md#214-233). The migration replicates those orchestrations in Python.

## 3. Target Python Package Layout

Create a new directory `pyBall/pyocl_dft/` (parallel to `pyBall/DFT/`) with the following modules:

| Module | Responsibility | C++ Source to Replace |
|---|---|---|
| `context.py` | Device/context management, kernel compilation, buffer registry | `OCL_GridFF.cpp` (init/release), `OCL_DFT::makeMyKernels`, `OCL_PP::makeKernels_PP` |
| `fft.py` | clFFT/pyclfft planning, FFT/IFFT helpers, buffer allocation | `OCL_GridFF.cpp::initFFT`, `OCL_DFT::runFFT`, buffer bookkeeping |
| `io.py` | Loads/saves (`saveToBin`, `loadWfBasis`, `saveToXsf*`) using NumPy | `OCL_GridFF.cpp`, `OCL_DFT::loadWfBasis` |
| `density.py` | Implements `project_dens_GPU`, `project_denmat_GPU`, `projectAtomsDens0` using PyOpenCL | `OCL_DFT::projectAtomsDens`, `projectAtomsDens0`, `projectDenmat`, helper preparation routines |
| `potentials.py` | Wraps `poisson`, `convolve`, `gradient`, `evalLJC_QZs` | `OCL_GridFF.cpp` functions & `OCL_PP::evalLJC_QZs` |
| `relax.py` | Probe-particle relaxation, `getFEinStrokes`, parameter utilities | `OCL_PP::relaxStrokesTilted`, `getFEinStrokes` |
| `df.py` | Implements Giessibl frequency-shift convolution using `convolveZ` | `relax.cl::convolveZ` orchestration (currently C++ utilities) |
| `assets.py` | Discovers `Fdata/basis/`, validates grid descriptors, manages `acumCoef` presets | New |
| `tests/` (package) | Python-level regression harness mirroring `run.sh` scripts | New |

Provide a top-level facade `pyBall/pyocl_dft/__init__.py` exposing functions with the same names/signatures as `pyBall.DFT.oclfft` so existing user scripts import the new module transparently.

## 4. Detailed Migration Tasks

### 4.1 Stage 0 – Preparatory Work

1. **Dependency inventory**
   - Confirm availability of `pyopencl`, `pyclfft` (or alternative FFT bindings). Pin versions in `requirements.txt` or project docs.
   - Audit `Fdata/basis/` contents; document mandatory elements (per AFM.md §2.1 @/doc/Topics/AFM/AFM.md#115-123).
2. **Kernel packaging**
   - Copy `myprog.cl`, `GridFF.cl`, `relax.cl` to a runtime-accessible directory; ensure `context.py` loads them relative to repository root or an environment variable to mimic `cl_src_dir` semantics.
3. **API compatibility layer**
   - Design decorators to match current `ctypes` signatures (e.g., `poisson(ibuff_in:int, ibuff_out:int, dcell)`), raising `NotImplementedError` until each feature is ported.

### 4.2 Stage 1 – Context & Resource Management (C++: `OCL_GridFF.cpp`, `OCL_DFT` core)

| Task | Python Target | Notes |
|---|---|---|
| Context initialisation | `context.py:init()` | Create PyOpenCL context/queue, compile kernels, mirror `OCL_DFT::makeMyKernels()`; manage logging/verbosity. |
| Buffer bookkeeping | `context.BufferTable` class | Store metadata (`shape`, `dtype`, `format`) for each logical buffer; replace global arrays `buffers[]`, `Ns[]`. |
| Upload/download helpers | `context.upload`, `context.download` | Support float/float2/double conversions (cf. `upload_d`) and ensure host/device offsets follow memory rules (check retrieved memory `e22a0816...`). |
| Texture/image handling | `context.create_image3d` etc. | Mirror `newFFTimage`, respecting CL image formats required by `evalLJC_QZs_toImg`. |
| Grid descriptors | `assets.make_grid_descriptor` | Provide `pos0`, `dA/B/C`, ensure compatibility with kernels (per AFM.md supporting assets). |

### 4.3 Stage 2 – FFT Infrastructure (C++: `initFFT`, `runFFT`, `convolution` scaffolding)

1. Implement `fft.PlanCache` to encapsulate clFFT plan creation keyed by grid dimensions.
2. Provide `fft.run_fft(buffer_id, direction)` to enqueue forward/inverse FFTs, matching the existing three-buffer setup (`inputA`, `inputB`, `outputC`).
3. Add helpers to resize FFT buffers when grid dimensions change (mimic `newFFTbuffer` logic).

### 4.4 Stage 3 – Density Projection Pipeline (C++: `OCL_DFT` methods)

| Python Function | Replaces | Key Steps | Reference |
|---|---|---|---|
| `density.project_atoms_dens()` | `OCL_DFT::projectAtomsDens` | 1) prepare `float4` atom arrays & coefficient matrices; 2) set `acumCoef=[0.0, 2.0]` (default) or alternative; 3) call `projectOrbDenToGrid_texture`. | AFM.md Step 2 @/doc/Topics/AFM/AFM.md#41-54 |
| `density.project_atoms_dens0()` | `OCL_DFT::projectAtomsDens0` | Upload neutral atom templates, set `acumCoef=[1.0, -1.0]` for difference densities. | AFM.md Step 2 bullet on difference density |
| `density.project_denmat()` | `OCL_DFT::projectDenmat` | Re-implement `atoms2box` chunking in Python (likely using NumPy slicing) before calling `projectDenmatToGrid`. | AFM.md Step 2 alternative workflow |
| `density.prepare_basis()` | `OCL_DFT::loadWfBasis`, `convCoefs`, `assignAtomDensCoefs` | Parse `.wf1/.wf2`, perform 1D interpolation (NumPy/scipy), reorder orbitals, compute accumulators. | Appendix A, Supporting Assets |

Additionally, mirror supportive structures:

- **Coefficient conversion**: reproduce the `convCoefs` logic for s/p decomposition; include assertions for expected orbital counts to "fail loudly" per user rules.
- **Atom boxing**: implement `atoms2box` in Python for chunk-wise processing; optional but needed for large systems.

### 4.5 Stage 4 – Density-Based Force Fields (FDBM)

### 4.5.1 Physics Model Corrections

**CRITICAL**: The following corrections apply to the density-based AFM (FDBM) model:

1. **Pauli Repulsion** (NOT density summation):
   - Formula: E_pauli(R) = A * ∫ ρ_tip(r+R)^b * ρ_sample(r)^b dr
   - This is a **convolution** operation, not a sum of densities
   - Can be evaluated efficiently via FFT: FFT(ρ_tip^b) * FFT(ρ_sample^b), then inverse FFT
   - Parameters A and b are fitted to reference data

2. **Electrostatics** (NOT density summation):
   - Compute electrostatic potential V from sample density ρ_sample via Poisson equation:
     ∇²V = -4π * ρ_sample
   - Solve via FFT: V = IFFT[ -4π * FFT(ρ_sample) / |k|² ]
   - Then convolve with tip potential (tip charge density, e.g., CO molecule)
   - Force: F = -∫ ρ_tip(r) ∇V(r+R) dr (convolution of tip density with potential gradient)

3. **NEVER sum electron densities directly** - this is physically incorrect

### 4.5.2 Implementation Details

| Python Function | Responsibilities | Implementation Notes |
|---|---|---|
| `potentials.poisson()` | Solve ∇²V = -4πρ via FFT | FFT density, multiply by -4π/|k|², inverse FFT. Handle k=0 singularity (set to 0). |
| `potentials.convolve()` | Convolution via FFT | FFT both inputs, multiply element-wise, inverse FFT. Scaling factors must account for voxel size. |
| `potentials.pauli_convolution()` | Pauli repulsion via FFT | Raise densities to power b, FFT, multiply, inverse FFT, scale by A. |
| `potentials.gradient()` | Compute ∇V or ∇E fields | Two options: (1) Fourier derivative (multiply by ik), (2) Finite differences on grid. |
| `potentials.eval_ljc_qzs()` | Classical LJ + point charges | Legacy classical FF - keep for comparison. |
| `potentials.eval_ljc_qzs_to_img()` | Classical FF to 3D image | Legacy - keep for comparison. |

### 4.5.3 Practical Implementation Issues

#### Grid Alignment for Convolution
- **CO tip grid must have O-atom at grid origin (corner)**, not center
- Otherwise convolution result will be shifted incorrectly
- When generating tip density grid, ensure tip is positioned so O-atom is at (0,0,0) in grid coordinates
- Sample grid origin should be at corner (standard OpenCL convention)

#### Voxel Size and Units
- Voxel size (dA, dB, dC) must be consistent between sample and tip grids
- FFT convolution assumes same voxel spacing - resample if needed
- Units: density in e/Å³, potential in eV, forces in eV/Å
- When using Fourier derivatives, multiply by ik accounts for 1/Å units automatically

#### Force Computation Methods

**Option 1: Fourier Derivative (Preferred)**
- Compute gradient in Fourier space: ∇f = IFFT[ ik * FFT(f) ]
- More numerically stable than finite differences
- Must correctly account for voxel size: k values depend on grid dimensions and physical size
- Implementation: k_x = 2π * n_x / (N_x * dA), etc.
- Force on tip: F = -∇E(R) where E is the convolution result

**Option 2: Trilinear Interpolation with Explicit Forces**
- Store pre-computed force field on grid (F_x, F_y, F_z at each voxel)
- Interpolate forces at tip position using trilinear interpolation
- Avoids slow B-spline fitting (BsplineConv3D in GridFF.cl is expensive)
- Forces computed once, then sampled during relaxation
- Less accurate than B-spline but much faster

**Option 3: B-spline Interpolation (Accurate but Slow)**
- GridFF.cl implements cubic B-spline interpolation (basis(), dbasis() functions)
- Requires BsplineConv3D preprocessing step
- Analytical derivatives are accurate
- Too slow for real-time use in large systems

**Recommendation**: Use Fourier derivative for initial implementation, fallback to trilinear with explicit forces if performance issues arise.

### 4.5.4 Delta-Density for Electrostatics

**CRITICAL**: For electrostatics, we must use delta-density to ensure charge neutrality because Fireball only includes valence electrons.

#### Option 1: Delta-density with Neutral Atoms (RECOMMENDED - Already Implemented)

**Formula**: `delta_rho = rho_SCF - rho_NA`

- `rho_SCF`: Self-consistent density from Fireball (valence electrons only)
- `rho_NA`: Neutral atom density (sum of isolated atomic densities)
- Result: Only the deformation density (bonding effects) contributes to electrostatics

**Implementation Status**:
- **Fortran**: `firecore_dens2points(points, f_den=1.0, f_den0=-1.0)` - works
- **C++ (oclfft.py)**: `projectAtomsDens0(iOut, atypes, apos, acumCoef=[1.0,-1.0])` - works
- **Reference tests**: `tests/tDFT_pentacene/run.py` uses this approach
- **PyOpenCL (Grid.py)**: NOT IMPLEMENTED - only projects SCF density from sparse blocks

**Advantages**:
- Numerically stable (neutral atom density is smooth/delocalized)
- Already implemented in C++ reference system
- Charge neutrality enforced by construction

**Implementation for PyOpenCL**:
```python
# Need to add to GridProjector class:
def project_neutral_atom_density(self, atoms, grid_spec):
    """Project neutral atom (promolecule) density to grid."""
    # Build on-site density matrix from neutral atom occupations
    # Use same projection kernel as SCF but with different coefficients
    # acumCoef = [1.0, -1.0] for delta_rho = rho_SCF - rho_NA
    pass
```

#### Option 2: Delta-density with Nuclear Charges (Alternative)

**Formula**: `delta_rho = rho_SCF - rho_Nuc`

- `rho_Nuc`: Nuclear charge density (Gaussian or B-spline smeared)
- For C: +4 charge (compensates 4 valence electrons)
- For H: +1 charge (compensates 1 valence electron)
- For Si: +4 charge (compensates 4 valence electrons, not +14)

**Implementation Status**:
- NOT IMPLEMENTED anywhere in current codebase
- GridFF has B-spline smearing for point charges but not for nuclear density
- Would require Gaussian smearing of nuclear charges

**Advantages**:
- More physically direct (explicit nuclear charges)
- Can tune smearing width for numerical stability

**Disadvantages**:
- Not implemented - requires new code
- Smearing width parameter adds complexity
- May be less numerically stable than neutral atom approach

**Recommendation**: Implement Option 1 first (neutral atoms), add Option 2 later for comparison.

### 4.5.5 Pauli Repulsion Parameters

**Formula**: `E_pauli(R) = A * ∫ ρ_tip(r+R)^b * ρ_sample(r)^b dr`

#### Current Parameter Values

From code analysis:
- **A (amplitude)**: `A_pauli = 16.0` (default in `pyBall/OCL/AFM.py` and `tests/tAFM/test_fdbm.py`)
- **b (exponent)**: `b = 1.0` (hardcoded in current implementation)
- User mentioned 18.0 as alternative value for A
- Literature suggests b in range [0.8, 1.4]

#### Parameter Determination Strategy

**Step 1: Compare with Classical Force Fields**
Plot line profiles of:
1. Pauli repulsion with A=16.0, b=1.0
2. Lennard-Jones repulsion (C12/r^12 term)
3. Morse repulsion (De^{-α(r-r0)} term)
4. van der Waals attraction (C6/r^6 term)

Compare magnitudes and decay rates to ensure Pauli is in physically reasonable range.

**Step 2: Fit to Reference Data**
If available, fit A and b to:
- Reference AFM force curves
- DFT-calculated tip-sample interaction energies
- Experimental frequency shift data

**Step 3: Parameter Sensitivity Analysis**
Test different values:
- A: [12.0, 16.0, 18.0, 20.0]
- b: [0.8, 1.0, 1.2, 1.4]

Evaluate impact on:
- Force magnitude
- Force decay rate
- AFM image contrast

**Implementation Note**: For initial testing, use A=16.0, b=1.0 as starting point. Make A and b configurable parameters.

### 4.5.6 Step-by-Step Debugging Policy for FDBM Implementation

**CRITICAL**: The FDBM implementation must be validated step-by-step with rich visualization and logging. Each component must be tested independently before composition.

#### Debugging Philosophy

1. **Separate validation for each physical component**
2. **Rich visualization at every step** (2D cuts, line profiles, 3D projections)
3. **Quantitative metrics** (integrals, max/min values, charge neutrality)
4. **Reference comparison** (vs C++ implementation, analytical solutions)
5. **Persistent outputs** (PNG plots, log files, .npy arrays for review)

#### Step 1: Density Projection to Grid

**Goal**: Project SCF density from Fireball sparse matrix to 3D grid using PyOpenCL

**Inputs**:
- Fireball sparse density matrix (from `fc.get_rho_sparse()`)
- Atomic positions and types
- Grid specification (origin, step, ngrid)

**Outputs to validate**:
- `rho_grid`: 3D electron density [e/Å³]
- `rho_na_grid`: Neutral atom density [e/Å³]
- `rho_diff`: Delta density = rho_SCF - rho_NA [e/Å³]

**Validation checks**:
- Total charge: ∫ rho_grid dV ≈ total valence electrons
- Neutral atom charge: ∫ rho_na_grid dV ≈ total valence electrons
- Delta charge: ∫ rho_diff dV ≈ 0 (charge neutrality)
- Single-atom test: H atom should integrate to 1.0 e, C to 4.0 e

**Visualization**:
- XY slice at z=2.0 Å (above molecule plane)
- XZ slice through selected carbon atom
- Line profile above selected carbon atom (z-direction)
- Max projection (for overview)

**Files to save**:
- `rho_grid.npy`, `rho_na_grid.npy`, `rho_diff_grid.npy`
- `step1_rho_slices.png` (XY, XZ, YZ slices)
- `step1_rho_maxproj.png` (max projections)
- `step1_rho_lineprofile.png` (line profile above C atom)
- `step1_log.txt` (quantitative metrics)

**Reference**: Compare with C++ `projectAtomsDens0` from `tests/tDFT_pentacene/run.py`

#### Step 2: Electrostatics via Poisson Equation

**Goal**: Compute electrostatic potential from delta-density using FFT Poisson solver

**Inputs**:
- `rho_diff` from Step 1
- Grid step (voxel size)

**Outputs to validate**:
- `V_ES`: Electrostatic potential [eV]

**Validation checks**:
- Poisson equation: ∇²V = -4π * rho_diff (verify via finite difference)
- Potential range: reasonable magnitude (not diverging)
- Decay behavior: ~1/r for point charge test

**Visualization**:
- XY slice at z=2.0 Å
- XZ slice through selected carbon atom
- Line profile above selected carbon atom
- Compare with analytical Coulomb potential for point charge test

**Files to save**:
- `V_ES.npy`
- `step2_VES_slices.png`
- `step2_VES_lineprofile.png`
- `step2_VES_poisson_check.png` (verify ∇²V vs -4πρ)
- `step2_log.txt`

**Reference**: Compare with C++ `poisson` from `oclfft.py`

#### Step 3: Pauli Repulsion via Density Convolution

**Goal**: Compute Pauli repulsion energy via FFT convolution of tip and sample densities

**Inputs**:
- `rho_grid` from Step 1 (sample density)
- Tip density (CO molecule or Gaussian blob)
- Parameters: A_pauli, b (exponent)

**Outputs to validate**:
- `E_Pauli_field`: Pauli energy field [eV]
- `grads_E_Pauli`: Pauli force field [eV/Å]

**Validation checks**:
- Energy scale: comparable to classical LJ repulsion
- Decay behavior: exponential decay with distance
- Force magnitude: reasonable for AFM imaging (0-10 eV/Å)
- Grid alignment: tip O at origin, convolution shift verified

**Visualization**:
- XY slice at z=2.0 Å
- XZ slice through selected carbon atom
- Line profile above selected carbon atom
- Compare with LJ repulsion profile

**Files to save**:
- `E_Pauli_field.npy`, `grads_E_Pauli.npy`
- `step3_Epauli_slices.png`
- `step3_Epauli_lineprofile.png`
- `step3_Epauli_vs_LJ.png` (comparison with classical repulsion)
- `step3_log.txt`

**Reference**: Compare with existing `test_fdbm.py` implementation

#### Step 4: Electrostatics Tip-Sample Convolution

**Goal**: Compute electrostatic interaction energy via convolution of tip density with sample potential

**Inputs**:
- `V_ES` from Step 2 (sample electrostatic potential)
- Tip delta-density (or full density with charge compensation)
- Tip charge q_CO

**Outputs to validate**:
- `E_ES_field`: Electrostatic energy field [eV]
- `grads_E_ES`: Electrostatic force field [eV/Å]

**Validation checks**:
- Energy scale: reasonable for partial charges (~0.1-1 eV)
- Force magnitude: comparable to Pauli and vdW
- Charge effect: test with different q_CO values

**Visualization**:
- XY slice at z=2.0 Å
- XZ slice through selected carbon atom
- Line profile above selected carbon atom
- Compare with point-charge electrostatics

**Files to save**:
- `E_ES_field.npy`, `grads_E_ES.npy`
- `step4_EES_slices.png`
- `step4_EES_lineprofile.png`
- `step4_EES_vs_pointcharge.png`
- `step4_log.txt`

**Reference**: Compare with C++ convolution from `tests/tDFT_pentacene/run.py`

#### Step 5: Dispersion (Lennard-Jones C6/r^6)

**Goal**: Compute van der Waals dispersion as sum of C6/r^6 terms

**Inputs**:
- Sample atomic positions
- Tip position (or scan grid)
- C6 parameters per atom type
- C6_CO parameter for tip

**Outputs to validate**:
- `E_vdw_field`: Dispersion energy field [eV]
- `grads_E_vdw`: Dispersion force field [eV/Å]

**Validation checks**:
- Energy scale: attractive, ~0.1-10 eV
- Decay behavior: 1/r^6
- Force magnitude: comparable to Pauli and ES
- Sum over atoms vs grid-based convolution

**Visualization**:
- XY slice at z=2.0 Å
- XZ slice through selected carbon atom
- Line profile above selected carbon atom
- Compare with analytical C6/r^6

**Files to save**:
- `E_vdw_field.npy`, `grads_E_vdw.npy`
- `step5_Evdw_slices.png`
- `step5_Evdw_lineprofile.png`
- `step5_Evdw_vs_analytical.png`
- `step5_log.txt`

**Reference**: Compare with existing LJ implementation in `test_fdbm.py`

#### Step 6: Final Composed Calculation

**Goal**: Combine all components and run full AFM simulation

**Inputs**:
- All fields from Steps 3-5
- Scan grid (xy positions, z heights)
- Relaxation parameters

**Outputs to validate**:
- `E_total_field`: Total energy = E_Pauli + E_ES + E_vdw
- `F_total_field`: Total force = F_Pauli + F_ES + F_vdw
- `Fz_relax`: Relaxed force field after PP relaxation
- `df`: Frequency shift = -dFz/dz

**Validation checks**:
- Energy conservation in relaxation
- Force continuity across scan
- Reasonable AFM contrast (ring/bond features)
- Comparison with classical Morse/LJ results

**Visualization**:
- Component traces (Pauli, ES, vdW, total) at selected atoms
- Fz vs height curves
- df maps at different heights
- Comparison with classical FF

**Files to save**:
- `E_total_field.npy`, `F_total_field.npy`
- `Fz_relax.npy`, `df.npy`
- `step6_component_traces.png`
- `step6_Fz_vs_height.png`
- `step6_df_maps.png`
- `step6_vs_classical.png`
- `step6_log.txt`

**Reference**: Compare with existing `test_fdbm.py` full calculation

#### Plotting Conventions

**Coordinate system**:
- Sample molecule at z=0 plane
- XY view: slice at z=2.0 Å (above molecule)
- XZ view: slice through selected carbon atom (y=constant)
- Line profile: z-direction above selected carbon atom

**Color maps**:
- Density: 'magma' or 'inferno' (positive values)
- Delta density: 'bwr' with symmetric normalization (diverging)
- Potential: 'bwr' with symmetric normalization
- Energy fields: 'viridis' or 'plasma'
- Force fields: 'RdBu' for signed quantities

**Line profiles**:
- Plot all components on same axes for comparison
- Include classical FF reference curves
- Mark key distances (e.g., R0=3.2 Å for vdW minimum)

**Quantitative metrics in logs**:
- Min/max values
- Integral over grid
- Charge neutrality check
- Peak positions
- Decay rates

#### Directory Structure for Debug Outputs

```
tests/tAFM/pyocl_fdbm/
├── pentacene.xyz
├── run_pyocl_fdbm.py
├── debug/
│   ├── step1_density/
│   │   ├── rho_grid.npy
│   │   ├── rho_na_grid.npy
│   │   ├── rho_diff_grid.npy
│   │   ├── step1_rho_slices.png
│   │   ├── step1_rho_maxproj.png
│   │   ├── step1_rho_lineprofile.png
│   │   └── step1_log.txt
│   ├── step2_electrostatics/
│   │   ├── V_ES.npy
│   │   ├── step2_VES_slices.png
│   │   ├── step2_VES_lineprofile.png
│   │   ├── step2_VES_poisson_check.png
│   │   └── step2_log.txt
│   ├── step3_pauli/
│   │   ├── E_Pauli_field.npy
│   │   ├── grads_E_Pauli.npy
│   │   ├── step3_Epauli_slices.png
│   │   ├── step3_Epauli_lineprofile.png
│   │   ├── step3_Epauli_vs_LJ.png
│   │   └── step3_log.txt
│   ├── step4_electrostatics_conv/
│   │   ├── E_ES_field.npy
│   │   ├── grads_E_ES.npy
│   │   ├── step4_EES_slices.png
│   │   ├── step4_EES_lineprofile.png
│   │   ├── step4_EES_vs_pointcharge.png
│   │   └── step4_log.txt
│   ├── step5_dispersion/
│   │   ├── E_vdw_field.npy
│   │   ├── grads_E_vdw.npy
│   │   ├── step5_Evdw_slices.png
│   │   ├── step5_Evdw_lineprofile.png
│   │   ├── step5_Evdw_vs_analytical.png
│   │   └── step5_log.txt
│   └── step6_composed/
│       ├── E_total_field.npy
│       ├── F_total_field.npy
│       ├── Fz_relax.npy
│       ├── df.npy
│       ├── step6_component_traces.png
│       ├── step6_Fz_vs_height.png
│       ├── step6_df_maps.png
│       ├── step6_vs_classical.png
│       └── step6_log.txt
└── README.md (summary of validation results)
```

#### Implementation Order

1. **Step 1**: Implement density projection with validation
2. **Step 2**: Implement Poisson solver with validation
3. **Step 3**: Implement Pauli convolution with validation
4. **Step 4**: Implement electrostatics convolution with validation
5. **Step 5**: Implement dispersion with validation
6. **Step 6**: Compose all components and run full simulation

Each step must pass validation before proceeding to the next.

### 4.5.7 Implementation Findings from Step-by-Step Validation

Based on actual implementation in `tests/tAFM/pyocl_fdbm/run_pyocl_fdbm.py`, the following critical issues were discovered and fixed:

#### Issue 1: vdW Divergence at Close Distances

**Problem**: Pure C6/r^6 dispersion diverges when the tip approaches the sample (z < 2 Å), producing forces of -5000 eV/Å that completely dominate all other components.

**Root Cause**: The simple pairwise C6/r^6 model has no repulsive wall. At typical AFM scan ranges (z_end ~ 0.4-2.0 Å), the tip can get arbitrarily close to atoms, causing r^-6 to blow up.

**Fix**: Added regularization `RA2_VDW = 1.5^2` Å² (from existing `test_fdbm.py`):
- Energy: `E = -C6_eff / (r^2 + RA2)^3`
- Force: `F_a = -6 * C6_eff * r_a / (r^2 + RA2)^4`

**Result**: vdW forces now in reasonable range [-0.4, 0.1] eV/Å, comparable to Pauli repulsion.

#### Issue 2: Force Sign Convention Consistency

**Problem**: vdW force had wrong sign (repulsive instead of attractive) because of inconsistent conventions between energy gradient storage and force computation.

**Root Cause**: In Steps 3-4 (Pauli, ES), `grads_E` stores energy gradients (∂E/∂x), and force is computed as `-grads_E` in Step 6. But in Step 5 (vdW), the direct force formula was stored instead of the gradient, leading to double negation.

**Fix**: Store energy gradients (∇E) consistently in all steps, compute force as `-∇E` only in Step 6 during composition.

**Lesson**: All intermediate fields must use the same convention - store either energy gradients or forces, but not mix them.

#### Issue 3: np.gradient Sign with Dec probe_heights

**Problem**: `np.gradient` for z-component gave wrong sign because `probe_heights` is decreasing (e.g., 6.0 → 0.4).

**Root Cause**: `np.gradient(E_vdw)` without spacing computes derivative in index units. When the array index increases as physical z decreases, the returned derivative is the negative of the physical derivative.

**Fix**: Pass actual spacing to `np.gradient`:
```python
dz = probe_heights[1] - probe_heights[0]  # negative value
grads_z = np.gradient(E_vdw, dz, axis=2)
```

**Lesson**: Always pass physical spacing to `np.gradient`, especially when array coordinates decrease with index.

#### Issue 4: Charge Conservation with Grid Projection

**Problem**: `delta_rho` integrates to ~0.23 e instead of 0.0 (target < 0.1).

**Root Cause**: Grid projection with 0.2 Å step and 10 SCF iterations has finite numerical error. Pentacene has 102 valence electrons; integrated charges were 102.15 (SCF) and 101.92 (NA).

**Mitigation**: Using finer grid (0.15 Å) and more SCF iterations (50) improves accuracy. For production, use step ≤ 0.15 Å and nscf ≥ 50.

#### Issue 5: Electrostatics Very Weak

**Problem**: With `q_CO = -0.05`, electrostatic forces are only [-0.008, 0.038] eV/Å - negligible compared to Pauli [0, 2.2] and vdW [-0.36, -0.005].

**Analysis**: The small tip charge produces weak electrostatic interaction. This may be realistic for CO tip, but means AFM contrast is dominated by Pauli + vdW, not ES.

**Recommendation**: For systems with strong charge transfer or polar molecules, increase `q_CO` or use Mulliken-based electrostatics model.

#### Issue 6: Interpolation from Density Grid to Scan Grid

**Problem**: Pauli and ES forces computed on density grid (152×88×96) must be interpolated to scan grid (40×40×16) for composition with vdW.

**Solution**: Use `scipy.ndimage.map_coordinates` with trilinear interpolation at scan positions mapped to density grid fractional coordinates.

**Performance**: Interpolation is fast (< 1s for 40×40×16 scan points), but for larger scans consider computing Pauli/ES directly on scan grid via FFT convolution of subsampled density.

#### Validation Results (Pentacene, 40×40 scan, z=5.0-2.0 Å)

| Component | Fz Range [eV/Å] | Notes |
|---|---|---|
| Pauli | [0.0, 2.24] | Repulsive, dominant at z < 3.2 Å |
| vdW | [-0.36, -0.005] | Attractive, dominant at z > 3.5 Å |
| ES | [-0.008, 0.007] | Negligible with q_CO=-0.05 |
| Total (raw) | [-0.037, 1.88] | Balance point ~3.2-3.5 Å |
| Total (relax) | [-0.037, 1.85] | PP relaxation has small effect at this resolution |
| df | [-3.79, 0.02] | Mostly negative (attractive regime) |

#### Files Created

- `tests/tAFM/pyocl_fdbm/run_pyocl_fdbm.py` - Main implementation script
- `tests/tAFM/pyocl_fdbm/debug/step{1-6}/` - Debug outputs for each step
- `tests/tAFM/pyocl_fdbm/pentacene.xyz` - Input molecule

### 4.6 Stage 5 – Probe Relaxation & Frequency Shift (C++: `OCL_PP`)

1. `relax.relax_strokes_tilted()` – replicate FIRE-based iteration (`update_FIRE` logic lives in kernel) by setting kernel arguments, including scan grid, step sizes, damping parameters.
2. `relax.get_fe_in_strokes()` – capture force/energy along predefined trajectories without relaxation.
3. `df.convolve_z()` – wrap the `convolveZ` kernel to produce frequency-shift values using weight arrays; provide convenience functions to generate Giessibl weights.
4. Provide Python utilities to configure scan grids mirroring `makeStartPointGrid` and `setGridShapePP`.

### 4.7 Stage 6 – I/O and Utilities

- Recreate `saveToBin`, `loadFromBin`, `saveToXsf`, and friends using NumPy and shared helper functions (`io.py`).
- Implement logging/verbosity toggles matching `setVerbosity`, `setErrorCheck`.
- Ensure Python functions reuse the "Supporting Assets" rules from AFM.md (basis availability, buffer formats, etc.).

## 5. Validation & Testing Strategy

1. **Unit tests** for each Python wrapper using small synthetic systems (e.g., single hydrogen atom) to confirm kernel argument wiring.
2. **Integration tests** replicating `tests/tDFT_pentacene/run.py` and related harnesses. Update each directory’s `run.sh` to allow `PYTHON_AFMPATH=pyBall.pyocl_dft` (feature flag) to switch between implementations.
3. **Numerical parity checks**:
   - Compare density grid norms and pointwise differences after `project_atoms_dens`.
   - Compare Hartree potentials (`poisson`), VdW grids (`eval_ljc_qzs`), and relaxed positions.
   - Validate `convolveZ` outputs against saved reference traces.
4. **Performance benchmarking** to ensure PyOpenCL overhead is acceptable; gather timings similar to existing `GridFF_cl::make_MorseFF()` prints.

## 6. Risks and Mitigation

| Risk | Impact | Mitigation |
|---|---|---|
| Missing basis data | Python path fails to load orbital textures | Add pre-flight checks in `assets.py`; provide descriptive error messages referencing @/doc/Topics/AFM/AFM.md#115-123. |
| FFT API differences | Incorrect scaling/orientation | Cross-check with C++ output; wrap clFFT plan creation to mimic stride/order; add regression tests. |
| Memory alignment | Kernel reads wrong offsets | Mirror buffer-offset logic noted in internal memory reminders (host/device offsets) and validate with assertions before kernel launch. |
| Numeric drift | Breaks regression comparisons | Keep computations in `float32` to match GPU; where double precision is required, document conversions. |

## 6.1 Critical Implementation Notes

### Reuse of `OpenCLBase`

All new PyOpenCL solver classes introduced by this migration must inherit from `pyBall.OCL.OpenCLBase.OpenCLBase`. The base class already encapsulates:

- Device selection (`select_device`) and queue creation with logging hooks.
- Kernel compilation (`load_program`) plus header parsing utilities used across existing OpenCL components.
- Buffer bookkeeping (`create_buffer`, `check_buf`, `toGPU`, `fromGPU`).

Adopting it ensures the new Python host remains aligned with current GPU tooling (e.g., atomic MMFF solvers). When implementing modules such as `pyocl_dft/context.py`, wrap the functionality in subclasses that extend `OpenCLBase` rather than rebuilding context logic from scratch.

### FFT Initialisation Patterns

The project already integrates clFFT/pyclfft in several locations:

- `tests/tMMFF/run_test_GridFF_ocl.py`
- `pyBall/tests/ocl_GridFF.py`
- `pyBall/OCL/GridFF.py`
- `pyBall/OCL/clUtils.py`

These files demonstrate correct plan creation, buffer registration, and execution order (FFT → kernel → inverse FFT) using the shared utility functions. The PyOpenCL AFM migration should follow the same conventions—initialise FFT plans via the helper routines in `clUtils.py` (or refactor shared helpers if necessary) and reuse the `GridFF_cl` approach for organising kernels and buffers. Documenting this reuse avoids subtle bugs such as incorrect strides or plan lifetime mismatches.

### Fireball DFT Integration

The Python AFM host continues to obtain molecular orbital data from the Fireball Fortran library via `pyBall.FireCore`. Critical entry points:

- **Initialisation sequence** (`FireCore.initialize`, `FireCore.preinit`, `FireCore.init`): sets verbosity, prepares species tables, and loads atomic coordinates. Example usage is shown in `tests/pyFireball/relax_molecules.py`.
- **SCF driver** (`FireCore.SCF`, `FireCore.assembleH`, `FireCore.solveH`, `FireCore.updateCharges`): run after `initialize` to populate the wavefunction coefficients (`wfcoef`) and density matrix.
- **Coefficient access**:
  - `FireCore.get_wfcoef(ikp, wfcoef)` retrieves the full molecular-orbital coefficient matrix as a NumPy array (shape `norbitals × norbitals`).
  - `FireCore.getPointer_wfcoef` exposes the raw pointer when zero-copy access is required (ensure lifetime management is understood before using).
- **Density exports**:
  - `FireCore.setupGrid(Ecut, g0, ngrid, dCell)` configures the real-space grid used by Fireball’s own projection routines and returns `(ngrid, dCell, lvs)` descriptors that should be forwarded to the OpenCL host.
  - `FireCore.getGridDens(ngrid=..., Cden=1.0, Cden0=0.0)` outputs an electron density grid directly from Fireball if we want to bypass the GPU projection for validation.
  - `FireCore.dens2points(points, f_den=1.0, f_den0=0.0)` evaluates the density at arbitrary coordinates, useful for spot checks.
  - **Grid sizing & ordering:** Fireball’s `setupGrid` will snap the requested lattice to FFT-friendly dimensions (e.g. returning `48×52×52` instead of the nominal `50×50×50`). Always persist the exact `ngrid` reported by the Fortran driver and forward it unchanged to the GPU pipeline. The Fortran buffers are column-major; when exposing them to NumPy reverse the shape (`ngrid[::-1]`) or allocate Fortran-contiguous arrays so that `z` remains the fastest axis when plotting or projecting.
- **Reference scripts**: explore `tests/pyFireball/relax_molecules.py` for end-to-end initialisation and relaxation, and `pyBall/FireCore.py` for the authoritative list of bindings to `fortran/MAIN/libFireCore.f90`.

During the migration, design the PyOpenCL pipeline so that the Fireball SCF step can hand off either (a) the MO coefficient matrix for GPU projection (`projectOrbDenToGrid_texture`) or (b) the already-sampled density grid for comparison. Record grid metadata (`pos0`, `dA/B/C`, `ngrid`) immediately after `setupGrid` so both paths remain consistent.

## 7. Roadmap & Milestones

1. **Milestone A – Fireball Data Extraction**
   - Verify Fireball SCF workflow via `pyBall.FireCore` to obtain wavefunction coefficients (`get_wfcoef`) and real-space density grids (`getGridDens`).
   - Capture grid descriptors from `setupGrid` and persist both the density matrix (2D NumPy array) and reference density (3D NumPy array) to disk for GPU comparison.
   - Provide `tests/pyocl_dft/test_firecore_data.py` CLI script that loads a sample molecule (from `tests/pyFireball`), runs the Fireball calculation, saves the coefficient/density files, and plots diagnostic slices using `matplotlib`.
2. **Milestone B – Skeleton & Context**
   - Implement `context.py`, `fft.py`, stub wrappers.
   - Validate kernel compilation and simple FFT round-trip.
   - Produce a standalone smoke-test script (now superseded by milestone A artefacts).
3. **Milestone C – Density Projection**
   - Port basis loading, coefficient conversion, `project_atoms_dens`.
   - Achieve parity on simple systems; add unit tests.
   - Deliver `tests/pyocl_dft/test_density.py` CLI script that loads Fireball-produced coefficient/density files, runs the Python host projection, saves the complex grid to `.bin`/`.xsf`, and plots a density slice using `matplotlib`.
4. **Milestone D – Potentials**
   - Implement FFT-based Hartree/Pauli/VdW wrappers.
   - Validate against C++ outputs for pentacene/CO tests.
   - Provide `tests/pyocl_dft/test_potentials.py` CLI script that loads density files, performs Hartree/Pauli convolutions, saves resulting grids, and plots radial profiles with `matplotlib`.
5. **Milestone E – Relaxation & df**
   - Port probe relaxation, stroke utilities, `convolveZ`.
   - End-to-end AFM image produced via Python host.
   - Ship `tests/pyocl_dft/test_relaxation.py` CLI script that relaxes a minimal scan grid, stores trajectories, and plots force traces/frequency shifts using `matplotlib`.
6. **Milestone F – Integration & Documentation**
   - Update tutorials, ensure `run.sh` scripts can switch hosts.
   - Collect benchmarks, document known differences.
   - Create `tests/pyocl_dft/test_end_to_end.py` CLI script orchestrating the full pipeline (density → potentials → relaxation → df), writing summary artifacts and generating figure panels for documentation.

## 8. Deliverables Checklist

- [ ] `pyBall/pyocl_dft/` package with modules described in §3.
- [ ] Unit and integration tests passing with Python host enabled.
- [ ] Updated documentation referencing the new Python path (AFM.md already compatible).
- [ ] Migration guide for developers (how to switch, debug, extend).
- [ ] Performance comparison report (Python vs C++ host).

Once all boxes are checked and parity is confirmed, the project can decide whether to deprecate `libOCL_GridFF.so` or keep it as a fallback for legacy workflows.

