# Molecule-Surface Interactions

FireCore specializes in simulations of molecules on rigid substrates (e.g., AFM/STM tips, crystal surfaces). Rather than treating the substrate atomistically, the code precomputes the interaction potential on a 3D grid or fits it to a compact basis. This document covers both approaches.

**Related Windsurf Codemaps:**
- [FireCore Force Field Navigation: From Audit Docs to Implementation](https://windsurf.com/codemaps/d550a435-7c6f-47b1-aeb5-efc9b098564f-fe86ab10a43f3d18) — Interactive trace through GridFF B-spline, Ewald2D, Surface.cl, SurfaceEwald.py, and test parity.
- [Surface Potential Evaluation: GridFF B-spline and XYZ Rigid Kernels](https://windsurf.com/codemaps/2a639fae-c9cb-407a-9d45-7b806c90c749-8796fe608a7d71c1) — GridFF interpolation kernel trace (basis(), fe3d_pbc(), sample3D()).
- [FoldedAtomicFunctions: Surface Potential Basis Fitting System](https://windsurf.com/codemaps/c9fc44a7-57a2-47c5-906f-886fa301ccc7-8796fe608a7d71c1) — FAF basis fitting and evaluation pipeline.
- [Interactive GridFF Scanning: PTCDA-on-CaF2 Constrained Relaxation System](https://windsurf.com/codemaps/99d506e2-223b-4ae7-bb60-8c2498fedfb9-8796fe608a7d71c1) — Real-world GridFF scanning with constrained relaxation.
- [Molecule-Substrate Interaction Energy Scanning: Assembly, GUI Placement, Force Fields & Surface Evaluation](https://windsurf.com/codemaps/38bd3cb6-31c0-45b6-9e09-fda94257999c-8796fe608a7d71c1) — Full molecule-on-surface energy scanning pipeline.
- [Molecule-on-Surface Systems: GridFF, XYZ Scanning, Surface Sampling, and Assembly](https://windsurf.com/codemaps/f8407e23-3a2e-41f1-abcf-9c15f3644c41-8796fe608a7d71c1) — GridFF + rigid-body assembly for surface adsorption.
- [AFM PyOpenCL System: Morse/LJ Path and FDBM Density-Based Path](https://windsurf.com/codemaps/9bb4c2a5-0c38-4943-abe9-254cfdcc75af-8796fe608a7d71c1) — AFM force evaluation using GridFF and FDBM models.
- [AFM Simulation: GPU Rigid Body Dynamics, CPU GridFF Relaxation, and Interactive GUI](https://windsurf.com/codemaps/594f7eaf-c3ab-4139-8f20-d1d2d7f8d401-fe86ab10a43f3d18) — GPU rigid-body AFM with GridFF substrate potential.
- [AFM FDBM Pipeline: DFTB Backend & pySCF Integration Points](https://windsurf.com/codemaps/02d559c9-de47-4058-b07b-3318664b454e-fe86ab10a43f3d18) — DFTB-based force-field parameter derivation for AFM.
- [DFTB Reference Calculation & FDBM AFM Forcefield Comparison System](https://windsurf.com/codemaps/1153fe89-ff29-4d4b-b4a6-e97d8f37047f-fe86ab10a43f3d18) — DFTB reference vs classical surface interaction comparison.
- [Rigid Body Dynamics on Surfaces (pyOpenCL)](https://windsurf.com/codemaps/b5d9c2d2-50f0-4ba7-bc65-60db6e06e423-8796fe608a7d71c1) — Rigid-body dynamics with GridFF surface sampling.
- [Rigid Body Dynamics System for AFM Simulation](https://windsurf.com/codemaps/c9f13e1f-edfa-4702-814f-5036d03ea6c9-fe86ab10a43f3d18) — 6-DOF rigid body AFM tip mechanics.

---

## 1. GridFF (Grid-based Force Field)

### Physics & Purpose

GridFF precomputes the substrate's interaction potential on a regular 3D grid, then evaluates it at atom positions via interpolation. The potential is separated into three independent channels:

1. **Pauli repulsion** — short-range, steeply rising, prevents penetration.
2. **London dispersion** — medium-range, attractive, $~1/r^6$ decay.
3. **Coulomb electrostatics** — long-range, $~1/r$ decay, requires special treatment for convergence.

The total interaction energy for atom $i$ at position $\mathbf{r}_i$ is:

$$E_i = P_i \cdot V_{\text{Pauli}}(\mathbf{r}_i) + L_i \cdot V_{\text{London}}(\mathbf{r}_i) + Q_i \cdot V_{\text{Coulomb}}(\mathbf{r}_i)$$

where $P_i$, $L_i$, $Q_i$ are the atom's PLQ parameters (see `nonbonding_forcefields.md`). This linear scaling allows the same grid to be reused for all atom types.

### Implementation Files

- **`cpp/common/molecular/GridFF.h`** — Main `GridFF` class (inherits from `NBFF`):
  - `FFPaul[]`, `FFLond[]`, `FFelec[]` — precomputed 3D grid arrays (Quat4f or Quat4d).
  - `grid` — `GridShape` object defining origin, spacing, and dimensions.
  - `allocateFFs()` — allocates grid memory.
  - `evalDipole()` — computes dipole moment and center of charge for the substrate.
  - `autoNPBC()` — automatically determines periodic image count from cell dimensions.
  - `Bspline_Pauli[]`, `Bspline_London[]`, `Bspline_Coulomb[]` — B-spline coefficient arrays.
  - `ewald` — pointer to `EwaldGrid` for long-range electrostatics.
- **`cpp/common_resources/cl/GridFF.cl`** — OpenCL interpolation kernels:
  - `basis(float u)` / `dbasis(float u)` — cubic B-spline basis functions and derivatives.
  - `fe1D()` / `fe2d()` / `fe3d_pbc()` — 1D, 2D, and 3D B-spline interpolation with periodic boundary conditions.
  - `sample3D()` — samples the grid at arbitrary points, returning $(f_x, f_y, f_z, E)$.
  - `sample3D_grid()` — resamples the grid onto a different grid (e.g., for visualization).
- **`cpp/libs_OCL/OCL_GridFF.cpp`** — OpenCL wrapper for grid construction and sampling.
- **`pyBall/OCL/GridFF.py`** — PyOpenCL interface:
  - `GridFF_cl` class: compiles `GridFF.cl`, manages buffers, and exposes `sample3D()`, `sample3D_comb()`.
  - `sample3D_comb()` — samples all three channels (Pauli, London, Coulomb) simultaneously.
- **`pyBall/OCL/GridFF_new.py`** — Refactored variant with improved buffer management.
- **`pyBall/OCL/GridFFRelaxedScan.py`** — Routines for relaxed potential energy surface (PES) scans.

### Key Physics

**Grid Construction**:
The grids are populated by summing pairwise potentials from all substrate atoms (with periodic images). For a grid point $\mathbf{g}$:

$$V_{\text{Pauli}}(\mathbf{g}) = \sum_{j \in \text{substrate}} \sum_{\mathbf{s} \in \text{PBC}} \frac{R_j^6}{|\mathbf{g} - \mathbf{r}_j + \mathbf{s}|^{12}}$$

$$V_{\text{London}}(\mathbf{g}) = \sum_{j \in \text{substrate}} \sum_{\mathbf{s} \in \text{PBC}} -2 \frac{R_j^6}{|\mathbf{g} - \mathbf{r}_j + \mathbf{s}|^6}$$

The Coulomb grid requires **Poisson solving** because the $1/r$ tail does not converge with a finite number of PBC images. FireCore supports two approaches:

1. **FFT-based Poisson solver** — computes the potential in reciprocal space via $\tilde{V}(\mathbf{k}) = 4\pi \tilde{\rho}(\mathbf{k}) / k^2$.
2. **Ewald summation** — splits the potential into short-range (real-space) and long-range (reciprocal-space) parts. See `EwaldGrid.h`.

**B-Spline Interpolation**:
GridFF uses **cubic B-splines** (degree 3) for interpolation. The basis functions are:

$$B_0(u) = \frac{1}{6}(1-u)^3$$
$$B_1(u) = \frac{1}{6}(3u^3 - 6u^2 + 4)$$
$$B_2(u) = \frac{1}{6}(-3u^3 + 3u^2 + 3u + 1)$$
$$B_3(u) = \frac{1}{6}u^3$$

with $u = (x - x_i)/\Delta x$ being the fractional coordinate within the cell. The derivatives are:

$$B'_0(u) = -\frac{1}{2}(1-u)^2, \quad B'_1(u) = \frac{1}{2}(3u^2 - 4u)$$
$$B'_2(u) = \frac{1}{2}(-3u^2 + 2u + 1), \quad B'_3(u) = \frac{1}{2}u^2$$

The interpolation is **C² continuous**, meaning both the potential and its gradient (force) vary smoothly across cell boundaries. This is critical for energy conservation in MD.

**Periodic Boundary Conditions in Z**:
The interpolation kernel handles PBC in the x-y plane (substrate periodicity) but clamps in the z-direction (no periodic images above/below the slab). Atoms far above the surface sample the tail of the potential; atoms below the surface receive a large repulsive force.

### Performance Considerations

- **Grid Resolution**: Typical grid spacing is $0.1$–$0.2$ Å. Finer grids improve accuracy but increase memory and construction time quadratically.
- **Memory Layout**: Each grid cell stores a `float4` (3 force components + energy). For a $200 \times 200 \times 100$ grid, this is ~60 MB per channel, or ~180 MB total.
- **GPU Sampling**: The `sample3D` kernel is memory-bandwidth bound. Each atom requires reading 64 neighboring grid points (4×4×4 for tricubic interpolation). Coalescing is improved by processing atoms in spatially contiguous batches.
- **PBC Index Selection**: The kernel precomputes `xqs[4]` and `yqs[4]` in `__local` memory to avoid repeated modulo arithmetic for periodic indices.

### GridFF Modes

`GridFF.h` defines several interpolation modes via the `GridFFmod` enum:

| Mode | Interpolation | Precision | Use Case |
|------|--------------|-----------|----------|
| `Direct` | No interpolation (exact sum) | Double | Reference calculations |
| `LinearFloat` | Trilinear | Float | Fast, low accuracy |
| `LinearDouble` | Trilinear | Double | Medium accuracy |
| `HermiteFloat` | Hermite spline | Float | Smooth gradients |
| `HermiteDouble` | Hermite spline | Double | High accuracy |
| `BsplineFloat` | Cubic B-spline | Float | Production default |
| `BsplineDouble` | Cubic B-spline | Double | Highest accuracy |

The default is `BsplineDouble` (see `GridFF.h:155`).

### Dipole Approximation

For very large substrates, constructing the full Coulomb grid is expensive. `GridFF` offers a dipole approximation (`evalDipole()`, line 52):

1. Compute the substrate's total charge $Q$, dipole moment $\mathbf{p}$, and center of charge $\mathbf{R}_0$.
2. For grid points far from the substrate ($r > r_{\text{dipole}}$), replace the exact sum with the multipole expansion:
   $$V_{\text{Coulomb}}(\mathbf{r}) \approx \frac{Q}{|\mathbf{r} - \mathbf{R}_0|} + \frac{\mathbf{p} \cdot (\mathbf{r} - \mathbf{R}_0)}{|\mathbf{r} - \mathbf{R}_0|^3}$$

This reduces the construction cost from O(N_substrate × N_grid) to O(N_grid) for the far-field region.

### Test Coverage

- `tests/tMMFF/run_test_GridFF.py` — Basic functionality: construct grid, sample at points, compare to direct sum.
- `tests/tMMFF/run_test_GridFF_CaF2.py` — CaF₂(111) surface, validates against known adsorption sites.
- `tests/tMMFF/run_test_GridFF_gauss_smear.py` — Tests Gaussian charge smearing for Coulomb grid convergence.
- `tests/tMMFF/run_test_GridFF_ocl.py` / `run_test_GridFF_ocl_new.py` — OpenCL vs CPU parity.
- `tests/tMMFF/GridFF_CaF2_doc_tutorial.md` — Step-by-step tutorial.

---

## 2. FoldedAtomicFunctions (FAF)

### Physics & Purpose

While GridFF is general and accurate, it requires storing millions of grid points. **FoldedAtomicFunctions** (FAF) provide a compact analytic representation of 2D-periodic substrate potentials using a separable basis:

$$V(x, y, z) = \sum_{n=1}^{N_z} \sum_{m=1}^{N_x} c_{nm} \cos(k_m x) \cdot \phi_n(z)$$

where:
- $\cos(k_m x)$ are plane waves (or sines) capturing the lateral periodicity.
- $\phi_n(z)$ are exponential decay functions (or polynomials) capturing the vertical decay.
- $c_{nm}$ are fitted coefficients.

This reduces storage from O(N_x × N_y × N_z) grid points to O(N_x × N_z) coefficients, typically a few hundred numbers instead of millions.

### Implementation Files

- **`doc/py/FoldedAtomicFunctions/potentials.py`**:
  - `morse_potential(z, D, a, r0)` — 1D Morse potential.
  - `generate_morse_samples()` — creates a matrix of Morse potentials with varying parameters.
  - `cut_2d_potential_profiles()` — extracts 1D slices from a 2D FAF potential.
- **`doc/py/FoldedAtomicFunctions/faf_func.py`** — Minimal functional re-implementation:
  - `makePotentialXZ()` — constructs a 2D periodic potential by summing Morse + Coulomb contributions from substrate atoms with periodic replicas.
  - `scan_Qs()` — fits the potential for several charge values and plots error maps.
  - `scan_basis_a0()` — scans the exponential decay parameter $a_0$ to optimize basis accuracy.
- **`doc/py/FoldedAtomicFunctions/FoldedAtomicFunction.md`** — Design document.
- **`doc/Topics/OnSurfaceAssembly/FoldedSubstratePotential_OpenCL.md`** — Discussion of OpenCL implementation.

### Key Physics

**Basis Functions**:
The standard FAF basis uses:

$$\phi_{nm}(x, z) = \cos\left(\frac{2\pi n x}{L_x}\right) \cdot e^{-a_m z}$$

where $L_x$ is the substrate periodicity and $a_m$ are decay constants (typically logarithmically spaced). For symmetric cells, only cosine terms are needed; for asymmetric cells, sine terms are also included.

**Fitting Procedure**:
1. Compute the reference potential $V_{\text{ref}}(x_i, z_j)$ on a dense grid by direct summation of pairwise potentials.
2. Construct the basis matrix $\Phi$ where $\Phi_{ij, nm} = \phi_{nm}(x_i, z_j)$.
3. Solve the least-squares problem:
   $$\mathbf{c} = \arg\min_{\mathbf{c}} ||\mathbf{V}_{\text{ref}} - \Phi \mathbf{c}||^2$$
   using `numpy.linalg.lstsq` (SVD-based, robust to rank deficiency).

**Error Analysis**:
The fit quality is assessed in the physically relevant region ($z > z_{\text{min}}$, typically 2–3 Å above the surface):

$$\text{MaxErr} = \max_{z > z_{\text{min}}} |V_{\text{ref}} - V_{\text{fit}}|$$
$$\text{AvgErr} = \frac{1}{N} \sum_{z > z_{\text{min}}} |V_{\text{ref}} - V_{\text{fit}}|$$

Errors near the surface (small $z$) are typically larger because the potential varies rapidly there, but this region is physically inaccessible due to Pauli repulsion.

### Performance Considerations

- **Evaluation Cost**: Computing $V(x,z)$ requires $N_x \times N_z$ basis evaluations. For $N_x=8$, $N_z=5$, this is only 40 terms—comparable to evaluating a single cosine, and ~1000× fewer operations than GridFF interpolation.
- **GPU Local Memory**: The coefficient array (40 floats) fits entirely in GPU local memory, eliminating the global memory bandwidth bottleneck that limits GridFF.
- **Pareto Frontier**: `tutorial_folded_basis_pareto.md` discusses the tradeoff between basis size (speed) and accuracy. The optimal basis is typically on the "knee" of the error vs. size curve.

### Charge-Dependent Fitting

For systems with polar substrates, the Coulomb potential depends on the adsorbate's charge. Rather than fitting separate potentials for each charge state, FAF can fit the potential as a linear function of charge:

$$V(x, z; Q) = V_0(x, z) + Q \cdot V_1(x, z)$$

where $V_0$ and $V_1$ are each expanded in the FAF basis. This doubles the coefficient count but allows arbitrary charge states to be evaluated at no additional cost.

### Test Coverage

- `doc/py/FoldedAtomicFunctions/tutorial_folded_basis_pareto.md` — Tutorial with Pareto analysis.
- `doc/py/FoldedAtomicFunctions/optimize_z_basis.md` — Optimization of the z-basis decay parameters.
- `tests/tMMFF/run_test_GridFF.py` — GridFF tests also serve as reference data for FAF fitting.

---

## 3. Surface.cl — Unified Surface Interaction Kernel

### Purpose

`Surface.cl` (`cpp/common_resources/cl/Surface.cl`) is a unified OpenCL kernel file that consolidates multiple surface-interaction evaluation strategies into a single compilation unit. It contains:

1. **Brute-force surface interaction** (`getSurfMorse`) — pairwise Morse/LJ/Coulomb summation against substrate atoms with PBC.
2. **Folded basis evaluation** (`getSurfFolded`, `getSurfFolded_workgroup`, `getSurfFolded_harmonics`) — compact FAF potential on GPU.
3. **2D Ewald electrostatics** (`compute_ewald_coefficients`, `eval_potential_vacuum`, `eval_potential_full`, `eval_potential_brute`) — exact long-range Coulomb for charged slabs.
4. **Macroscopic layer corrections** (`getMacroRectLayers`, `macro_phi_rect_charge`, `macro_phi_rect_dipole`) — multipole corrections for thick substrates.

### Key Kernels

**`getSurfMorse`** (line 334):
- Evaluates molecule-surface interaction by direct pairwise summation over substrate atoms with periodic images.
- Uses `getMorsePLQH()` for Morse + PLQ interactions.
- Includes an optional macroscopic correction (`bMacro`) for long-range electrostatics using precomputed layer multipoles.
- Uses `__local` memory tiling (`LATOMS[32]`, `LCLJS[32]`) to coalesce substrate atom reads.

**`getSurfFolded`** (line 432):
- Evaluates the FoldedAtomicFunctions basis directly on the GPU.
- Reads coefficients and basis parameters into `__local` memory.
- Computes basis values and gradients via `folded_eval_basis()` and `folded_eval_grad()`.
- Supports up to 64 basis functions and 8 atom types.

**`getSurfFolded_workgroup`** (line 512):
- Workgroup-optimized variant that precomputes 1D basis factors (`L_BX`, `L_BY`, `L_BZ`, `L_dBX`, `L_dBY`, `L_dBZ`) into `__local` memory.
- Avoids register spilling by streaming precalculated factors from local memory during the triple loop.
- Uses `native_cos`, `native_sin`, `native_exp` for speed.

**`getSurfFolded_harmonics`** (line 656):
- Variant using separable harmonic (sine/cosine) basis instead of exponential decay.
- Stores 1D parameters in `LBASIS` array: `[Nx params, Ny params, Nz params]`.
- **Status**: Partially implemented (TODO at line 687).

### Helper Functions

| Function | Line | Purpose |
|----------|------|---------|
| `getMorsePLQH()` | 172 | Morse potential with PLQ factorization + H-bond |
| `getCoulomb()` | 189 | Damped Coulomb potential and force |
| `macro_phi_rect_dipole()` | 197 | Rectangular sheet dipole potential (analytic) |
| `macro_phi_rect_charge()` | 227 | Rectangular sheet monopole potential (analytic) |
| `getMacroRectLayers()` | 237 | Layered substrate macroscopic correction |
| `folded_eval_basis()` | 260 | FAF basis function $\cos(k_x x)\cos(k_y y)\exp(-a z)$ |
| `folded_eval_grad()` | 268 | Gradient of FAF basis (force computation) |
| `getR4repulsion()` | 294 | Short-range R⁻⁴ repulsive blob |

### Ewald2D Kernels

**`compute_ewald_coefficients`** (line 716):
- Computes complex coefficients $C_G$ for vacuum evaluation and per-ion weights $w[g,i]$ for full evaluation.
- Formula: $C_G = \frac{2\pi}{A|G|} \sum_i q_i \exp(|G|z_i) \exp(-i\mathbf{G}\cdot\boldsymbol{\rho}_i)$.
- One work item per G-vector; loops over all ions.

**`eval_potential_vacuum`** (line 772):
- Evaluates potential for $z > z_{\max}$ (above all ions).
- Formula: $\phi = \text{Re}\left[ \sum_G C_G \exp(i\mathbf{G}\cdot\boldsymbol{\rho}) \exp(-|G|z) \right]$.
- **Key optimization**: Precomputes $z1_{b1} = \exp(i\mathbf{b}_1\cdot\boldsymbol{\rho})$ and $z1_{b2} = \exp(i\mathbf{b}_2\cdot\boldsymbol{\rho})$, then uses complex multiplication to compute $\exp(ih\mathbf{b}_1\cdot\boldsymbol{\rho})$ by repeated squaring. Reduces $N_G$ trig evaluations to just 2 per point.

**`eval_potential_full`** (line 834):
- Evaluates potential at arbitrary $z$ (inside or outside slab).
- Formula: $\phi = -\frac{2\pi}{A} \sum_i q_i |z-z_i| + \text{Re}\left[ \sum_G \sum_i w_{g,i} \exp(i\mathbf{G}\cdot\boldsymbol{\rho}) \exp(-|G||z-z_i|) \right]$.

**`eval_potential_brute`** (line 911):
- Direct Coulomb sum over periodic images for validation.
- Sums over circular shells of replicas (N_rep shells).
- Slow but exact; used for parity testing against Ewald2D.

---

## 4. SurfaceEwald.py — OpenCL Ewald2D Python Wrapper

### Purpose

`pyBall/OCL/SurfaceEwald.py` provides the `SurfaceEwaldCL` class, a PyOpenCL wrapper around the Ewald2D kernels in `Surface.cl`. It mirrors the Python reference implementation (`pyBall/Ewald2D.py`) but accelerates coefficient computation and potential evaluation on the GPU.

### Key Methods

| Method | Purpose |
|--------|---------|
| `__init__(platform)` | Initialize OpenCL context and compile `Surface.cl` |
| `make_reciprocal_2d(a_vec, b_vec)` | Compute reciprocal lattice vectors and unit cell area |
| `generate_G_vectors(b1, b2, n_harm)` | Generate G-vectors for \|h\|,\|k\| ≤ n_harm |
| `prepare_system(ion_data, a_vec, b_vec, n_harm)` | Upload ions, compute $C_G$ and $w$ on GPU |
| `eval_vacuum(X, Y, z)` | Evaluate potential on XY grid at fixed $z$ (above slab) |
| `eval_full(X, Y, Z)` | Evaluate potential at arbitrary 3D positions |
| `eval_brute(X, Y, Z, N_rep)` | Brute-force Coulomb sum for validation |

### Performance Considerations

- **Fast mode** (`use_fast=True` in `eval_vacuum`): Uses `__local` memory for G-vector caching when $n_{\text{harm}} \leq 16$.
- **Complex multiplication trick**: The kernel precomputes $\exp(i\mathbf{b}_1\cdot\boldsymbol{\rho})$ and $\exp(i\mathbf{b}_2\cdot\boldsymbol{\rho})$ once per point, then raises to powers $h$ and $k$ via repeated complex multiplication. This is ~$N_G/2$× faster than evaluating $\cos(G\cdot\rho)$ for each G-vector independently.
- **Buffer reuse**: $C_G$ and $w$ are computed once per system and reused for multiple evaluation grids.

---

## 5. Surface_utils.py — GridFF Alignment & FDBM Fitting Utilities

### Purpose

`pyBall/OCL/Surface_utils.py` is a high-level glue layer that imports and reuses existing modules (`GridFF.py`, `RigidBodyAFM.py`, `InteractionEnergy.py`) with minimal new code. It provides utilities for:

1. **GridFF I/O** — Loading `.npy` grids with JSON metadata validation.
2. **GridFF alignment verification** — Detecting origin/shift conventions.
3. **Visualization** — Plotting grids with atom overlays.
4. **FDBM (Force-Field Density-Based Model) fitting** — Building feature matrices and mock references for linear fitting.

### Key Functions

| Function | Purpose |
|----------|---------|
| `load_gridff_metadata(grid_path)` | Load JSON metadata (g0, dg, ns, lvec) |
| `load_gridff_array(path)` | Load `.npy` grid with channel validation (3→4 channel expansion) |
| `load_bspline_gridff(grid_path)` | Load grid + metadata, validate shape consistency |
| `init_gridff_sampler_md(grid_path, apos0, nSystems)` | Initialize `MolecularDynamics` for fast batch GridFF sampling |
| `sample_gridff_channels_rigid(md, transforms, PLQH_channels)` | Sample Pauli/London/Coulomb channels for many rigid transforms |
| `fdbm_build_feature_matrix(Es_PLQ, type_ids, ntypes)` | Build linear regression feature matrix from per-type channel sums |
| `fdbm_make_mock_reference(...)` | Generate mock reference energies with controlled noise |
| `compare_electrostatics_methods(...)` | **Comprehensive comparison**: GridFF vs Ewald2D vs Brute Force |
| `load_substrate_xyz_with_lvec(path)` | Load substrate XYZ and extract lattice vectors from comment line |
| `infer_grid_metadata(grid_path, substrate_info)` | Infer grid origin conventions from substrate geometry |
| `plot_gridff_diagnostics(...)` | Diagnostic plotting with atom overlay and molecule positions |

### FDBM (Force-Field Density-Based Model)

The FDBM fitting utilities in `Surface_utils.py` support a linear model:

$$E_{\text{ref}} = \sum_t P_t \sum_{i \in t} V_{\text{Pauli}}(\mathbf{r}_i) + \sum_t L_t \sum_{i \in t} V_{\text{London}}(\mathbf{r}_i) + \sum_i Q_i V_{\text{Coulomb}}(\mathbf{r}_i) + \text{noise}$$

where $P_t$ and $L_t$ are per-type Pauli and London scaling parameters fitted to DFT reference data. The `fdbm_build_feature_matrix()` function constructs the design matrix for this linear least-squares problem.

---

## 6. Test Coverage for Surface Interactions

### Electrostatics Parity Tests

- **`tests/tMMFF/test_electrostatics_comparison.py`** — **Primary validation script**.
  - Compares three methods on NaCl surfaces: GridFF (OpenCL B-spline), Ewald2D (Python), Brute Force (direct sum).
  - Generates 1D line scans (z on Na, z on Cl, z midpoint) and 2D XY/XZ slices.
  - Optionally tests OpenCL Ewald (`--test_opencl`) against Python reference.
  - Reports RMSE and max error; asserts RMSE < 1e-5 eV.
  - Saves diagnostic PNGs and JSON report.

- **`tests/tMMFF/test_gridff_alignment.py`** — GridFF alignment verification.
  - Uses `pyBall.OCL.Surface_utils.run_alignment_verification`.
  - Detects origin conventions (centered XY, z at top atom vs. bottom vs. zero).
  - Validates that sampled forces/energies match direct pairwise sums.

### Folded Basis Tests

- **`tests/tMMFF/test_folded_fit_nacl1x1.py`** — FoldedAtomicFunctions fit on NaCl(001).
  - Uses `SurfaceEwaldCL` to compute Coulomb reference potential.
  - Fits FAF basis to Morse + Coulomb reference.
  - Validates fit quality (RMSE, max error) across parameter scans.

### GridFF Generation & FDBM Tests

- **`tests/tMMFF/gen_gridff_nacl_gpu.py`** — GPU-accelerated GridFF generation for NaCl surfaces.
  - Uses `pyBall.OCL.Surface_utils.load_substrate_xyz_with_lvec` to load substrate.
  - Computes DFT electrostatics via `pySCF` or `DFTB+` (optional).
  - Generates B-spline PLQd grids with metadata JSON.

- **`tests/tMMFF/test_fdbm_fit_dft.py`** — FDBM fit against DFT reference data.
  - Uses `pyBall.OCL.Surface_utils.load_bspline_gridff` to load precomputed grids.
  - Fits per-type Pauli/London parameters to match DFT adsorption energies.

- **`tests/tMMFF/test_fdbm_fit_gridff_mock.py`** — Mock FDBM fitting with synthetic data.
  - Tests the linear fitting pipeline without requiring DFT calculations.

- **`tests/tMMFF/gui_fdbm_fit.py`** — Interactive PyQt5 GUI for FDBM parameter tuning.
  - Live 2×2 plot: potential, energy, force, error vs. parameter sliders.
  - Uses `pyBall.OCL.Surface_utils` as the shared backend.

---

## Comparison: GridFF vs. FAF

| Property | GridFF | FoldedAtomicFunctions |
|----------|--------|----------------------|
| **Storage** | O(N_x × N_y × N_z) grid points (~100 MB) | O(N_x × N_z) coefficients (~1 KB) |
| **Evaluation** | B-spline interpolation (64 reads/point) | Direct basis summation (40 terms) |
| **GPU Memory** | Global (bandwidth-bound) | Local (register-bound, fast) |
| **Accuracy** | Exact up to interpolation order | Approximate (fit residual) |
| **PBC** | Native (x-y periodic, z clamped) | Native (x periodic via cos basis) |
| **Coulomb** | FFT/Ewald required | Fitted empirically (short-range) |
| **Use Case** | Production MD, arbitrary surfaces | Rapid prototyping, repetitive evaluations |

---

## Surface Ewald Summation

For charged substrates, the Coulomb grid requires Ewald summation to converge. FireCore implements this via:

- **`cpp/common/molecular/EwaldGrid.h`** — Ewald grid construction and solving.
- **`pyBall/MMFF.py`** — Python bindings for `setupEwaldGrid()`, `projectAtomsEwaldGrid()`, `EwaldGridSolveLaplace()`.

The Ewald method splits the potential:

$$V(\mathbf{r}) = V_{\text{short}}(\mathbf{r}) + V_{\text{long}}(\mathbf{r})$$

$$V_{\text{short}}(\mathbf{r}) = \sum_{\mathbf{n}} \sum_i \frac{Q_i}{|\mathbf{r} - \mathbf{r}_i + \mathbf{n}|} \text{erfc}\left(\alpha |\mathbf{r} - \mathbf{r}_i + \mathbf{n}|\right)$$

$$V_{\text{long}}(\mathbf{r}) = \text{FFT}^{-1}\left[ \frac{4\pi}{k^2} e^{-k^2/4\alpha^2} \tilde{\rho}(\mathbf{k}) \right]$$

where $\alpha$ is the splitting parameter. The short-range part is evaluated directly with a cutoff; the long-range part is computed in reciprocal space via FFT.

---

## See Also

- [Topical Audit Index](topical_audit.md) — priority ranking, dependency graph, missing topics
- [Forcefields Overview](forcefields_overview.md) — high-level taxonomy of all force field classes
- [Non-Bonding Forcefields](nonbonding_forcefields.md) — NBFF, exclusion schemes, FMM, PME
- [AFM/STM Simulation](afm_stm_simulation.md) — AFM pipeline that uses GridFF for substrate interactions
- [Molecular Topology](molecular_topology.md) — topology and graph representations underlying force field evaluation

---

*Last updated: 2026-06-23*
