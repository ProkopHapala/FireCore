# Intramolecular Force Fields

This document covers force fields that describe bonded geometry: bond lengths, angles, dihedrals, and related constraints. FireCore implements several distinct approaches, from traditional topology-based force fields (UFF, MMFFsp3) to position-based solvers (ProjectiveDynamics, XPBD) and rigid-body approximations.

**Related Windsurf Codemaps:**
- [FireCore Classical Forcefields: MMFFsp3 & UFF (CPU/GPU/Python)](https://windsurf.com/codemaps/53f2fe2c-ac5c-4c0b-b905-af6653adde97-8796fe608a7d71c1) — Architecture of UFF, MMFFsp3, GPU kernels, Python wrappers.
- [MMFF/UFF CPU vs GPU Testing: C++ OpenCL and PyOpenCL Parity Infrastructure](https://windsurf.com/codemaps/8d1b056f-1502-4363-b52d-8257de4be453-8796fe608a7d71c1) — How CPU vs GPU parity tests are structured.
- [XPBD Molecular Dynamics pyOpenCL](https://windsurf.com/codemaps/2e558e51-fdbe-4bd4-8732-7818724d4ced-8796fe608a7d71c1) — Position-based dynamics solver trace.
- [Rigid Body Dynamics on Surfaces (pyOpenCL)](https://windsurf.com/codemaps/b5d9c2d2-50f0-4ba7-bc65-60db6e06e423-8796fe608a7d71c1) — Rigid body quaternion integration and GridFF coupling.
- [Rigid Body Dynamics System for AFM Simulation](https://windsurf.com/codemaps/c9f13e1f-edfa-4702-814f-5036d03ea6c9-fe86ab10a43f3d18) — 6-DOF rigid body with AFM tip mechanics.
- [FitREQ_PN: Hydrogen-Bond Parameter Fitting System](https://windsurf.com/codemaps/d977d597-94b4-42c3-a92a-0cefe34a3e82-8796fe608a7d71c1) — Automated H-bond parameter optimization.
- [FitREQ Interactive GUI: Monte Carlo Optimization & Energy Decomposition Integration](https://windsurf.com/codemaps/e25a0dfc-f9a8-42ab-b8bb-1d959037ca68-fe86ab10a43f3d18) — Interactive parameter tuning GUI.
- [FitREQ Hydrogen Bond Fitting System - GPU-Accelerated Parameter Optimization](https://windsurf.com/codemaps/bf59a960-ac6c-4eea-b828-9bd18c3d44ac-fe86ab10a43f3d18) — GPU-accelerated H-bond parameter search.

---

## 1. UFF (Universal Force Field)

### Physics & Purpose

UFF is a general-purpose molecular mechanics force field designed for the entire periodic table. It expresses the intramolecular energy as a sum of bonded terms:

$$E_{\text{intra}} = \sum_{\text{bonds}} E_{\text{bond}} + \sum_{\text{angles}} E_{\text{angle}} + \sum_{\text{dihedrals}} E_{\text{dihedral}} + \sum_{\text{inversions}} E_{\text{inversion}}$$

Each term is a simple analytical function parameterized by equilibrium values and stiffness constants. The formulation is deliberately generic so that parameters can be assigned automatically from atomic properties (radii, electronegativity, etc.).

### Implementation Files

- **`cpp/common/molecular/UFF.h`** — Main `UFF` class. Holds arrays of bond, angle, dihedral, and inversion data structures. Implements force assembly, neighbor list construction, PBC handling, and total energy evaluation.
- **`cpp/common_resources/cl/relax_multi.cl`** — OpenCL kernel with inline evaluation functions:
  - `evalBond(float3 h, float dl, float k, float3* f)` (line 157) — harmonic bond stretching.
  - `evalAngCos(...)` (line 104) — angle bending via cos formulation (fast but unstable >90°).
  - `evalAngleCosHalf(...)` (line 120) — improved angle bending using cos(θ/2), quasi-harmonic beyond 90°.
  - `evalPiAling(...)` (line 140) — π–π alignment interaction.
- **`pyBall/OCL/UFF.py`** / **`pyBall/OCL/UFFbuilder.py`** — PyOpenCL wrappers and topology builders.
- **`pyBall/MMFF_multi.py`** — Multi-system UFF interface.

### Key Physics

**Bond stretching** (harmonic):
$$E_{\text{bond}} = \frac{1}{2} k (r - r_0)^2$$
Implemented in `evalBond()` as `fr = dl*k; *f = h*fr; return fr*dl*0.5;`.

**Angle bending** — two formulations exist in the code:

1. `evalAngCos()` uses `c = cos θ = ĥ₁·ĥ₂` and force `f = -K(c-c₀)·2·(ĥ₂ - ĥ₁c)`.
   - Fast (~10 FLOPs) but fails for θ > 90° because the energy is not quasi-harmonic.

2. `evalAngleCosHalf()` uses the identity `cos(θ/2) = |ĥ₁+ĥ₂|/2`.
   - More expensive (requires two `sqrt`s) but behaves quasi-harmonically for all angles.
   - Energy: `E = k(1 - cs.x)` where `cs = (cos(θ/2-θ₀/2), sin(θ/2-θ₀/2))`.

### Performance Considerations

- **Force Assembly Scheme**: To avoid atomic writes on the GPU, UFF uses an *assembly* pattern:
  1. Each thread (per node atom) writes forces into a per-node auxiliary buffer.
  2. A second pass assembles these into the global `aforce` array.
  This is critical for parallel performance because atomic adds on global memory are ~10× slower than coalesced writes.
- **Multi-System GPU**: The `relax_multi.cl` kernel uses a 2D NDRange where `global.y = nSystems`. Each system has its own offset in the shared buffers, allowing thousands of molecules to relax in parallel.
- **Neighbor Lists**: Bonds, angles, and dihedrals are precomputed into compact neighbor lists (`neighs[4]`, `bkNeighs[4]` per atom), reducing the topological search to O(1) per term.

### Test Coverage

- `tests/tUFF/test_UFF_multi.py` — CPU vs GPU parity for multi-system evaluation.
- `tests/tUFF/test_UFF_ocl.py` — OpenCL-specific tests.
- `tests/tUFF/run_parity_uff.sh` — automated parity suite.

---

## 2. MMFFsp3_loc (Molecular Mechanics Force Field sp3 Localized)

### Physics & Purpose

MMFFsp3 is *not* the Merck Molecular Force Field. It is a FireCore-specific formulation optimized for sp³-hybridized carbon and first-row organics (H, C, N, O). The key insight is that traditional 4-body dihedral and improper terms are:

1. **Expensive** — require four-atom loops, poorly suited to GPU SIMT execution.
2. **Redundant** — the same physics (orbital alignment) can be captured by 2-body π–π and π–σ terms.

MMFFsp3 therefore **replaces explicit dihedrals** with:

- **π–π alignment** (`k_pp`): penalizes misalignment of p-orbitals on conjugated systems.
- **π–σ alignment** (`k_sp`): enforces orthogonality between π and σ orbital frameworks.

This is motivated by quantum chemistry: in sp-hybridized atoms, the three σ-bonds form a trigonal plane and the remaining p-orbital is perpendicular. The alignment terms are cheap 2-body interactions that preserve this geometry without explicit 4-body torsion terms.

### Implementation Files

- **`cpp/common/molecular/MMFFsp3_loc.h`** — Main `MMFFsp3_loc` class. Localized data layout where each "node" atom stores:
  - `k_pp` — π–π stiffness.
  - `k_sp` — π–σ stiffness.
  - `pipos` — π-orbital orientation vector.
  - `neighs`, `bonds`, `angles` — compact neighbor lists.
- **`cpp/common_resources/cl/relax_multi.cl`** — Same kernel file as UFF; MMFF terms reuse `evalBond()`, `evalAngCos()`, and add `evalPiAling()`.
- **`pyBall/MMFF.py`** — CPython `ctypes` bindings to the C++ library.
- **`pyBall/OCL/MMFF.py`** — Pure PyOpenCL implementation for rapid prototyping.

### Key Physics

**Node vs. Capping Atoms**:
- **Node atoms** (C, N, O with ≥2 neighbors) carry full angular terms and π-orbital orientation.
- **Capping atoms** (H, lone-pair dummies) have simplified evaluation: no angular terms, force is purely *recoil* from the host node atom. This reduces the number of active atoms in the kernel by ~50% for typical organic molecules.

**Recoil Force Assembly**:
For a bond A–B where B is a capping atom, the force on B is not computed independently. Instead, the force on A is computed, and B receives the equal-and-opposite recoil. This avoids divergent threads (capping atoms would otherwise execute empty angle kernels).

### Performance Considerations

- **Localized Data Layout**: All parameters for a node atom (bonds, angles, π-orientation) are stored in a contiguous struct, enabling coalesced GPU memory access.
- **Collision Damping**: A velocity-dependent damping term is added to stabilize the integrator:
  $$\mathbf{f}_{\text{damp}} = -c_{\text{damp}} \cdot \mathbf{v}$$
  This is essential because the π-alignment terms can create high-frequency oscillations.
- **PyOpenCL Overhead**: The Python wrapper (`pyBall/OCL/MMFF.py`) is convenient for debugging but incurs ~10× overhead in topology preparation. For production, the C++ harness is preferred.

### Test Coverage

- `tests/tMMFF/test_diamond_phonon_bands.py` — Phonon band structure of diamond; validates force constants by finite differences.
- `tests/tMMFF/test_MMFF_ocl_parity.py` — CPU vs GPU parity.
- `tests/tMMFF/run_parity_mmff.sh` — automated parity suite.

---

## 3. ProjectiveDynamics

### Physics & Purpose

Projective Dynamics (PD) is a position-based solver for stiff spring networks. Instead of integrating forces explicitly (which requires tiny timesteps for stiff bonds), PD solves a linear system at each timestep:

$$\mathbf{A} \mathbf{x} = \mathbf{b}$$

where:
- $\mathbf{A} = \mathbf{K} + \frac{\mathbf{M}}{h^2} + \mathbf{D}$ — system matrix (stiffness + inertia + damping).
- $\mathbf{b}$ — right-hand side constructed from predicted positions and spring constraints.
- $\mathbf{x}$ — corrected positions.

This is an *implicit* integration scheme: the inertia term $\mathbf{M}/h^2$ on the diagonal both stabilizes the solver and introduces proper dynamics. The timestep can be 10–100× larger than explicit integration.

### Implementation Files

- **`cpp/common/math/ProjectiveDynamics_d.h`** — Main `ProjectiveDynamics_d` class.
  - `makePDMatrix()` — assembles the sparse system matrix.
  - `updateJacobi()`, `updateGaussSeidel()`, `updateCholesky()`, `updateCG()` — iterative and direct solvers.
  - `updateMomentumMethod()` — momentum-accelerated Jacobi (heavy-ball / Nesterov).
- **`pyBall/pyTruss/truss.cl`** — OpenCL kernels:
  - `jacobi_iteration_sparse()` — parallel Jacobi with diagonal preconditioning.
  - `gauss_seidel_iteration_colored()` — Gauss-Seidel with graph coloring for parallelism.
  - `accum_vertex_hessian()` — local assembly of 3×3 blocks per vertex.
- **`doc/py/ProjectiveDynamics/projective_dynamics.py`** — Reference Python implementation using dense NumPy solvers.
- **`web/molgui_web/js/ProjectiveDynamics.js`** — JavaScript port for web-based demos.

### Key Physics

**Jacobi Iteration**:
$$x_i^{(k+1)} = \frac{1}{A_{ii}} \left( b_i - \sum_{j \neq i} A_{ij} x_j^{(k)} \right)$$
Fully parallel but converges slowly (spectral radius dependent).

**Gauss-Seidel**:
Uses updated values immediately: strictly convergent for SPD matrices, but sequential. Graph coloring (3–6 colors for typical molecular graphs) recovers parallelism.

**Momentum-Accelerated Jacobi**:
$$\mathbf{v}^{(k+1)} = \beta \mathbf{v}^{(k)} + \alpha (\mathbf{x}^{(k+1)} - \mathbf{x}^{(k)})$$
$$\mathbf{x}^{(k+1)} = \mathbf{x}^{(k)} + \mathbf{v}^{(k+1)}$$
With `β ≈ 0.9`, this achieves Gauss-Seidel-like convergence rates while remaining fully parallel.

### Performance Considerations

- **Sparse vs. Dense**: The C++ class uses a custom sparse format (CSR-like) for large systems (>10k atoms). The Python reference uses dense matrices and is limited to ~1k atoms.
- **Local Memory**: The OpenCL kernel preloads the diagonal 3×3 block into `__local` memory to avoid repeated global reads.
- **Vertex Block Descent (VBD)**: Each vertex owns a 3×3 block; updates are local and communication-free for Jacobi.

### Test Coverage

- `doc/py/ProjectiveDynamics/example_pd.py` — Basic truss dynamics demo.
- `pyBall/pyTruss/truss_solver_ocl.py` — OpenCL solver tests.

---

## 4. Rigid-Atom Based FF (XPBD_2D / XPDB_AVBD)

### Physics & Purpose

eXtended Position-Based Dynamics (XPBD) treats each atom as a **rigid body with 6 DOFs**: position $\mathbf{r}$ and rotation (quaternion $\mathbf{q}$). Each atom maintains a "wishlist" of where it wants its neighbors to be, expressed in its local frame. These *ports* are completely rigid, but flexibility emerges because:

1. The atom's frame can translate and rotate under neighbor forces.
2. Neighbors are pulled toward the ports by harmonic springs.

This is conceptually similar to "As-Rigid-As-Possible" (ARAP) mesh deformation, adapted to molecular graphs.

### Implementation Files

- **`pyBall/XPBD_2D/XPBD_2D.py`** / **`XPBD_2D.cl`** — 2D simulator using **complex numbers** for rotation:
  - Rotation of a vector $\mathbf{v}$ by angle θ: $\mathbf{v}' = \mathbf{v} \cdot e^{i\theta}$ (in OpenCL, `float2` complex multiplication).
  - Simplifies 2D angular constraints to complex arithmetic.
- **`pyBall/XPDB_AVBD/XPDB.py`** / **`XPDB.cl`** — 3D XPBD with **Angular-Velocity-Based Dynamics**:
  - Uses quaternions for 3D rotation.
  - Integrates angular velocity $\boldsymbol{\omega}$ explicitly, then updates quaternion via:
    $$\dot{\mathbf{q}} = \frac{1}{2} \mathbf{q} \otimes (0, \boldsymbol{\omega})$$
- **`pyBall/XPDB_AVBD/RRsp3.py`** / **`RRsp3.cl`** — Rigid-atom sp³ solver with port-based topology.

### Key Physics

**Port-Based Constraints**:
For atom A with neighbor B, the port is the ideal position of B in A's local frame: $\mathbf{p}_{AB}^{\text{local}}$. During simulation:

1. Transform port to world space: $\mathbf{p}_{AB}^{\text{world}} = \mathbf{R}_A \mathbf{p}_{AB}^{\text{local}} + \mathbf{r}_A$.
2. Spring force pulls B toward port: $\mathbf{f}_B = k (\mathbf{p}_{AB}^{\text{world}} - \mathbf{r}_B)$.
3. Torque on A: $\boldsymbol{\tau}_A = (\mathbf{p}_{AB}^{\text{world}} - \mathbf{r}_A) \times \mathbf{f}_B$.

This naturally encodes bond lengths and angles without explicit angle potentials.

**Analytic Procrustes Problem**:
After a position update, the optimal rotation that maps current port positions to target positions is found by solving the orthogonal Procrustes problem via SVD of the cross-covariance matrix. This is implemented in `Analytic_Procrustes_Problem.md`.

### Performance Considerations

- **Cluster Collision Detection**: Atoms are grouped into spatial clusters with AABB bounding boxes. Clusters far apart are skipped, reducing collision checks from O(N²) to O(N log N).
- **Local Memory**: The OpenCL kernel caches neighbor positions and ports in `__local` arrays for the workgroup.
- **Heavy-Ball Momentum**: Like ProjectiveDynamics, XPBD uses momentum to accelerate convergence of the constraint solver.

### Test Coverage

- Extensive discussion documents in `pyBall/XPBD_2D/` and `pyBall/XPDB_AVBD/`.

---

## 5. Rigid Body FF

### Physics & Purpose

Rigid Body Dynamics treats entire molecules as rigid bodies with only **6 translational/rotational DOFs** (or 7 with quaternion normalization). This is ideal for:

- **Large time steps** — no high-frequency intramolecular vibrations to resolve.
- **Machine learning** — fixed-size descriptor (6 DOFs per molecule, regardless of atom count).
- **Quantum mechanics coupling** — molecular orbitals transform rigidly with the molecule, enabling STM/AFM simulations where electronic structure is precomputed.
- **On-surface assembly** — rigid molecules can be rapidly packed and evaluated for steric clashes.

### Implementation Files

- **`cpp/common/molecular/RigidBodyFF.h`** — `RigidBodyFF` class:
  - `poses` — position + quaternion per body.
  - `updateQuaternions()` — integrates angular velocity to update rotation.
  - `projectAtoms()` — maps body pose to atomic coordinates.
- **`cpp/common_resources/cl/Rigid.cl`** — OpenCL kernels:
  - `rigid_body_dynamics_kernel()` — integrates position and quaternion with optional external forces (GridFF).
  - `quat_mult()`, `make_qrot()`, `qrot_omega()` — quaternion operations with Taylor expansion for small angles to avoid trigonometric functions.
- **`pyBall/OCL/RigidBodyDynamics.py`** — PyOpenCL wrapper:
  - Handles REQ → PLQ conversion for mixed-species interactions.
  - Computes mass properties (total mass, inertia tensor) from atomic masses.
  - Initializes GridFF for surface interactions.
- **`pyBall/OCL/Assembly.py`** / **`cl/Assembly.cl`** — Assembly evaluation:
  - `evaluate_packing_3d()` — scores packing density with clash penalties.
  - Super-fibonacci rotations for uniform orientation sampling.

### Key Physics

**Quaternion Integration**:
For angular velocity $\boldsymbol{\omega}$, the quaternion update is:
$$\mathbf{q}(t+h) = \mathbf{q}(t) \otimes \left( \frac{\boldsymbol{\omega}}{|\boldsymbol{\omega}|} \sin\frac{|\boldsymbol{\omega}|h}{2}, \cos\frac{|\boldsymbol{\omega}|h}{2} \right)$$

For small $|\boldsymbol{\omega}|h$, `Rigid.cl` uses a 6th-order Taylor expansion to avoid `sin`/`cos`:
```c
float2 sc = quat_factors_taylor(r2);  // r2 = |ω|²
return (float4)(omega * sc.x, sc.y);
```
This is ~10× faster than calling trigonometric functions for typical MD timesteps.

**Symmetric Recoil Forces**:
Forces from GridFF (sampled at each atom) are converted to net force and torque on the body:
$$\mathbf{F} = \sum_i \mathbf{f}_i, \quad \boldsymbol{\tau} = \sum_i (\mathbf{r}_i - \mathbf{r}_{\text{cm}}) \times \mathbf{f}_i$$
The lever arms are cached to avoid recomputation.

### Performance Considerations

- **One Workgroup per Body**: The OpenCL kernel uses a single workgroup (≤128 atoms) per rigid body, enabling fast `__local` memory reductions for force/torque summation.
- **GridFF Integration**: Surface forces are evaluated via B-spline interpolation on precomputed grids (see `surface_interactions.md`), avoiding O(N_surface) pairwise summation.
- **REQ → PLQ**: For interactions between different species, the mixing rule uses the geometric mean: $P_{ij} = \sqrt{P_i P_j}$, which requires storing `sqrt(E)` in the PLQ representation.

### Test Coverage

- `doc/Topics/RigidBodyAssembly/RigidBodyAssemblyDiscussion.md` — Comprehensive design discussion.
- `tests/tMMFF/test_assembly.py` — Packing and assembly tests.

---

## 6. Integration: MolWorld_sp3

### Physics & Purpose

`MolWorld_sp3` is the central orchestrator that integrates all force field components. It manages the simulation state, constructs the molecular graph, and routes forces to the appropriate solver.

### Implementation Files

- **`cpp/common/molecular/MolWorld_sp3.h`** — Main `MolWorld_sp3` class:
  - Instantiates `MMFFsp3`, `MMFFsp3_loc`, `UFF`, `ProjectiveDynamics_d`, `RigidBodyFF`, `NBFF`, and `GridFF`.
  - `makeFFs()` — builds force field parameters from atomic coordinates and types.
  - `optimize()` / `relax()` — high-level relaxation routines.
  - Handles PBC, constraints, and surface interactions.

### Key Physics

**Constraint Handling**:
Fixed atoms (anchors) are handled by zeroing their forces after each evaluation. For rigid bodies, anchor constraints are applied to the center-of-mass position and orientation.

**Solver Routing**:
The user selects the active force field at runtime (`bMMFF`, `bUFF`, `bRigid`). `MolWorld_sp3` dispatches to the appropriate `eval()` method and accumulates forces into the global `fapos` array.

### Performance Considerations

- **Lazy Initialization**: Force field objects are instantiated only when first requested (`bMMFF` must be `true` for UFF initialization—see memory note).
- **Buffer Sharing**: `MolWorld_sp3` owns the main `apos`, `fapos`, `REQs` arrays and binds them to sub-force fields via pointer sharing, avoiding data duplication.

---

## Summary Table

| Force Field | DOFs per Atom | Body Terms | Angle Terms | Dihedral Terms | Solver | GPU Support |
|-------------|---------------|------------|-------------|----------------|--------|-------------|
| **UFF** | 3 (x,y,z) | Harmonic | cos(θ), cos(θ/2) | Explicit 4-body | Explicit Euler | OpenCL |
| **MMFFsp3** | 3 + 2 (π-orientation) | Harmonic | cos(θ) | π–π, π–σ alignment | Explicit Euler | OpenCL |
| **ProjectiveDynamics** | 3 | Harmonic (stiff) | — | — | Implicit (Jacobi/GS/Cholesky) | OpenCL |
| **XPBD (2D)** | 2 + 1 (rotation) | Port-based | Implicit in ports | Implicit in ports | Position-based | OpenCL |
| **XPBD (3D)** | 3 + 3 (rotation) | Port-based | Implicit in ports | Implicit in ports | Position-based | OpenCL |
| **RigidBody** | 6 (per molecule) | — | — | — | Explicit Euler + quaternion | OpenCL |

---

## See Also

- [Topical Audit Index](topical_audit.md) — priority ranking, dependency graph, missing topics
- [Forcefields Overview](forcefields_overview.md) — high-level taxonomy of all force field classes
- [Non-Bonding Forcefields](nonbonding_forcefields.md) — NBFF, exclusion schemes, FMM, PME
- [Surface Interactions](surface_interactions.md) — GridFF, FoldedAtomicFunctions, Ewald2D
- [Web Force Fields](forcefields_web_implementation.md) — WebGL/WebGPU shader implementations
- [Molecular Topology](molecular_topology.md) — topology and graph representations underlying FF evaluation
- [Molecular Topology Types](molecular_topology_types.md) — atom type assignment and parameter loading

---

*Last updated: 2026-06-23*
