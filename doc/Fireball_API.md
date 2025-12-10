# FireBall Python API (High-Level Overview)

This document describes what the FireBall Fortran DFT code provides through the Python module `pyBall.FireCore`.

It is **use‑case–oriented**: focused on *what you can do* with the API, not on argument details. At the end there is a concise function list mapping Python functions to their underlying Fortran entry points.

---

## 1. Initialization and Basic Workflow

Typical pattern for using FireBall from Python is:

1. **Prepare a molecule** (atoms + positions), usually via `pyBall.AtomicSystem` or `atomicUtils`.
2. **Initialize FireBall** with those atoms.
3. **Run a calculation** (SCF, forces, relaxation, density projection, Hamiltonian export, …).

Conceptually:

```text
Atomic structure  →  FireCore.initialize(...)  →  run what you need (SCF/forces/density/HS export/relaxation)
```

Key points:

- Initialization ties FireBall’s internal data structures (basis, integrals, neighbor lists, etc.) to your atomic configuration.
- After initialization you can reuse the same FireBall state for many operations: repeated SCF calls, geometry scans, density evaluation, etc.

The tests in `tests/pyFireball` are good examples of complete workflows built on this pattern.

---

## 2. Single‑Point SCF and Forces (Energy/Force Evaluations)

**Use case:** Given a fixed geometry, compute electronic structure self‑consistently and obtain:

- Total energy and its components
- Forces on atoms
- Converged density / density matrix, later used for densities on grids/points

### Typical usage

From `tests/pyFireball/distort_molecule.py` and the H2O density demos:

1. Construct or load a molecule (`AtomicSystem`, XYZ file, etc.).
2. Call `fc.initialize(atomType=..., atomPos=...)` once.
3. Call `fc.evalForce(positions, nmax_scf=...)` whenever you want an SCF energy/force at a given geometry.

This enables:

- **Energy scans** along bond distances, angles, torsions, etc.
- **Interactive exploration** of potential‑energy surfaces by repeatedly perturbing coordinates and re‑evaluating.
- **Preparing a converged SCF density** for later density projection (see §§3–4).

`tests/pyFireball/distort_molecule.py` shows several examples:

- Distance scan `E(R)` between two atoms.
- Combined angle–distance or angle–scale scans, returning 1D/2D energy maps.

---

## 3. Geometry Relaxation (FIRE Optimizer)

**Use case:** Starting from an initial geometry, relax the structure to a local minimum using FireBall forces and an internal FIRE optimizer.

You give FireBall:

- Atomic types and coordinates
- Maximum number of relaxation steps and SCF iterations

FireBall then:

- Performs SCF + force evaluation at each step
- Updates atomic positions with FIRE dynamics
- Stops when forces satisfy a built‑in tolerance or max steps are reached

### Typical usage

- `tests/pyFireball/relax_molecule.py`:
  - CLI tool: `python relax_molecule.py input.xyz output_dir`
  - Loads a single molecule via `atomicUtils.AtomicSystem`.
  - Calls `fc.initialize(...)`, then `fc.relax(...)`.
  - Writes the relaxed geometry as XYZ.

- `tests/pyFireball/relax_molecules.py`:
  - Loops over many XYZ molecules in a directory.
  - For each: initialize → relax → print or save results → `fc.reload()` to reset FireBall between systems.

This covers use cases such as:

- Building small relaxed molecule libraries
- Pre‑relaxing structures before more expensive calculations
- Systematic comparison of equilibrium geometries across a set of molecules

---

## 4. Real‑Space Electron Density and Orbitals

Once SCF has been run (e.g. via `evalForce` or `SCF`), FireBall holds a converged density matrix and wavefunction coefficients. `FireCore` exposes multiple ways to project this information into real space.

### 4.1 Density at Arbitrary Points

**Use case:** Evaluate the electronic density (or combinations of SCF and neutral-atom density) at a set of arbitrary 3D points.

- In Python you pass an array of 3D points.
- FireBall returns a 1D array of density values at those locations.

This is used in:

- `tests/pyFireball/density_along_line.py`:
  - Load H2O, run SCF once via `evalForce`.
  - Build points along an O–H bond.
  - Call `fc.dens2points(points_on_line, f_den=1.0, f_den0=0.0)`.
  - Plot ρ(r) along the bond.

- `tests/pyFireball/plot_h2o_density_plane.py`:
  - Define a 2D plane through the molecule (O at origin, axes aligned with molecular geometry).
  - Build a dense 2D grid of points in that plane.
  - Call `fc.dens2points(points_on_grid, f_den=1.0, f_den0=0.0)`.
  - Reshape results to a 2D array and plot a density heatmap.

This enables:

- Line profiles of density along bonds or arbitrary paths
- 2D slices of density in molecular planes
- Arbitrary probe lines/planes for custom visualization or analysis

### 4.2 Density and Orbitals on 3D Grids

**Use case:** Build a regular 3D grid around the system and project either:

- Individual molecular orbitals
- Total electron density

onto that grid.

Workflow (as used in the example code at bottom of `FireCore.py`):

1. Call `fc.setupGrid(...)` to define or auto‑choose a real‑space grid.
2. Use `fc.getGridMO(iMO, ...)` to get a 3D array for the i‑th MO.
3. Use `fc.getGridDens(...)` to get a 3D array of density.

This is useful for:

- Building volumetric datasets for visualization (e.g. via external tools that read raw 3D arrays)
- Grid‑based analysis of density or orbitals (integration, projections, etc.)

### 4.3 Export to XSF (Orbital/Density Maps)

**Use case:** Write orbitals or density directly to `.xsf` files for visualization in tools such as XCrySDen.

- `fc.orb2xsf(iMO)` writes the specified MO onto a real‑space grid and exports an XSF file.
- `fc.dens2xsf(f_den0)` writes the density (and optionally additional reference density) to an XSF file.

These routines encapsulate both the projection to the grid and the file formatting.

---

## 5. Hamiltonian / Overlap Export and k‑Space Matrices

FireBall internally stores neighbor‑based Hamiltonian and overlap matrix blocks. The Python bindings expose this in a structured way.

### 5.1 Export Sparse H, S and Indexing Data

**Use case:** Inspect or reuse FireBall’s Hamiltonian/overlap at the neighbor block level, including all indexing data needed to reconstruct the full matrices.

High‑level workflow (see `tests/pyFireball/export_HS_sparse.py`):

1. Initialize a system and assemble the Hamiltonian once (`fc.assembleH(...)`).
2. Query dimensions via `dims = fc.get_HS_dims()`.
3. Allocate and fill all sparse data via `sparse_data = fc.get_HS_sparse(dims)`.
4. Use fields such as:
   - `sparse_data.h_mat`, `sparse_data.s_mat` – neighbor‑block matrices
   - `sparse_data.neighn`, `sparse_data.neigh_j`, `sparse_data.neigh_b` – neighbor lists and shells
   - `sparse_data.num_orb`, `sparse_data.lssh`, `sparse_data.mu`, `sparse_data.nu`, `sparse_data.mvalue`, `sparse_data.nzx`, etc. – orbital indexing and species data

In the demo, this is turned into a human‑readable description per atom:

- For each atom and neighbor, code prints the H and S sub‑blocks with simple orbital labels (s, px, py, pz, ...), mapped back to chemical elements via `pyBall.elements`.

This enables:

- Debugging orbital interactions and neighbor connectivity
- Building customized tight‑binding models from FireBall’s data
- Data export for post‑processing or fitting

### 5.2 Dense k‑Space Hamiltonian and Overlap

**Use case:** Obtain full complex Hermitian matrices H(k) and S(k) for a given k‑point.

After initialization and assembly, you can:

- Call `Hk, Sk = fc.get_HS_k(kvec, norbitals)` to get dense k‑space matrices.

`export_HS_sparse.py` demonstrates this at the Γ‑point (`k = (0,0,0)`) and prints the resulting matrices.

Applications include:

- Band‑structure calculations (once you loop over k)
- Interface with other codes that expect dense H(k), S(k)
- Eigenvalue analyses or custom electronic‑structure workflows in Python

---

## 6. Charges, MO Coefficients, and Internal Pointers

Beyond forces and grids, FireBall exposes more detailed electronic data.

### Charges

- `fc.getCharges(...)` allows you to retrieve Mulliken or Löwdin charges (depending on FireBall options) after SCF.
- This is useful for analysis of charge distribution, population analysis, and custom post‑processing.

### MO Coefficients

- `fc.get_wfcoef(...)` retrieves the MO coefficients matrix for a given k‑point.
- `fc.set_wfcoef(...)` writes a new MO coefficient vector back into FireBall’s storage.
- `fc.getPointer_wfcoef(...)` and related pointer routines expose internal Fortran arrays by pointer, for advanced users who want to manage the data at the C level.

These features are important if you want to:

- Export MO coefficients to other codes
- Implement custom projections or transformations of the MO basis in Python
- Perform more advanced electronic‑structure analysis directly on top of FireBall’s wavefunctions

---

## 7. Python API Reference (Functions & Fortran Bindings)

Below is a concise list of all functions defined in `pyBall/FireCore.py`, each mapped to its underlying Fortran entry point in `fortran/MAIN/libFireCore.f90` (where applicable), with a one‑sentence description.

### 7.1 Library and Utility Functions

- `fc.reload()` – (Python helper, reloads the shared `FireCore` library and reapplies ctypes signatures) Reloads the FireBall shared library, useful when running multiple independent calculations in a single Python session.

### 7.2 Core Control / Initialization

- `fc.setVerbosity(verbosity=0, idebugWrite=0)` → `firecore_setVerbosity` – Sets global Fortran verbosity and debug‑output flags, controlling how much diagnostic information FireBall prints.
- `fc.preinit()` → `firecore_preinit` – Performs FireBall’s basic initialization (constants, options, defaults) before any system is defined.
- `fc.set_lvs(lvs)` → `firecore_set_lvs` – Sets the lattice vectors for periodic calculations or clusters, controlling the simulation cell.
- `fc.init(atomTypes, atomPos)` → `firecore_init` – Initializes a specific atomic configuration in FireBall (allocates internal arrays, reads `info.dat`, sets up orbitals, charges, grid structures, etc.) and returns the number of orbitals.
- `fc.initialize(atomType=None, atomPos=None, verbosity=1)` – (Python helper calling `setVerbosity` → `preinit` → `init`) High‑level convenience function that sets verbosity and fully initializes FireBall for a given set of atom types and positions.

### 7.3 SCF, Forces, and Relaxation

- `fc.evalForce(pos, forces=None, nmax_scf=100, Es=None, ixyz=-1)` → `firecore_evalForce` – Runs an SCF cycle (up to `nmax_scf` iterations), computes forces and detailed energy components, and returns them for a given geometry.
- `fc.SCF(positions, iforce=0, nmax_scf=200)` → `firecore_SCF` – Performs a pure SCF loop for the given positions (optionally with forces), updating the density matrix without directly returning forces.
- `fc.relax(pos, forces=None, fixPos=None, nstepf=1000, nmax_scf=100, Es=None)` → `firecore_relax` – Runs a full geometry relaxation using FIRE, repeatedly doing SCF + forces and updating positions until convergence or maximum steps.

### 7.4 Charges and Basic Electronic Data

- `fc.getCharges(charges)` → `firecore_getCharges` – Fills an array with Mulliken or Löwdin charges per atom (depending on FireBall’s `iqout` option) after SCF.
- `fc.getPointer_wfcoef(bbnkre_ptr)` → `firecore_getPointer_wfcoef` – Returns a C pointer to the internal MO coefficient array `bbnkre`, for low‑level external access.
- `fc.get_wfcoef(wfcoef=None, norb=None, ikp=1)` → `firecore_get_wfcoef` – Retrieves the MO coefficient matrix for a given k‑point into a NumPy array.
- `fc.set_wfcoef(wfcoef, iMO=1, ikp=1)` → `firecore_set_wfcoef` – Writes a given MO coefficient vector into FireBall’s internal `bbnkre` array for a chosen orbital and k‑point.

### 7.5 Hamiltonian / Overlap Assembly

- `fc.assembleH(positions, iforce=0, Kscf=1)` → `firecore_assembleH` – Assembles the Hamiltonian and overlap matrices for the current positions and SCF iteration index, preparing data for diagonalization or export.
- `fc.solveH(k_temp=None, ikpoint=1)` → `firecore_solveH` – Diagonalizes the Hamiltonian at a specified k‑point (with given k‑vector) and updates eigenvalues/eigenvectors.
- `fc.updateCharges(sigmatol=1e-6, sigma=None)` → `firecore_updateCharges` – Updates the electronic density and charges using the current density matrix and returns the SCF residual `sigma`.

### 7.6 Real‑Space Grids and Projections

- `fc.setupGrid(Ecut=100, g0=None, ngrid=None, dCell=None)` → `firecore_setupGrid` – Configures or auto‑builds a 3D real‑space grid (cutoff, origin, spacing, lattice vectors) for projecting orbitals and density.
- `fc.getGridMO(iMO, ewfaux=None, ngrid=None)` → `firecore_getGridMO` – Projects a selected molecular orbital onto the current real‑space grid and returns a 3D array of orbital amplitudes.
- `fc.getGridDens(ewfaux=None, ngrid=None, Cden=1.0, Cden0=0.0)` → `firecore_getGridDens` – Projects the SCF and optional reference density onto the current grid and returns a 3D array of density values.
- `fc.orb2xsf(iMO)` → `firecore_orb2xsf` – Projects the chosen MO to the grid and writes it as an XSF file for volumetric visualization.
- `fc.dens2xsf(f_den0=0.0)` → `firecore_dens2xsf` – Projects the density (and optionally additional reference density) to the grid and writes it as an XSF file.

### 7.7 Point‑wise Projections and Atomic Orbitals

- `fc.getpsi(poss, ys=None, in1=1, issh=1, l=0, m=1)` → `firecore_getpsi` – Evaluates radial + angular atomic orbital values (Ylm * radial) for a given shell/type along a set of 3D positions.
- `fc.orb2points(poss, ys=None, iMO=1, ikpoint=1)` → `firecore_orb2points` – Projects a molecular orbital at a specific k‑point onto an arbitrary list of 3D points and returns the values.
- `fc.dens2points(points, f_den=1.0, f_den0=0.0, ewfaux_out=None)` → `firecore_dens2points` – Evaluates SCF and/or reference density at a set of 3D points (used for line scans and 2D planes in the H2O demos).

### 7.8 Hamiltonian/Overlap Export Routines

- `fc.get_HS_dims()` → `firecore_get_HS_dims` – Returns a small `FireballDims` object describing system and matrix dimensions (natoms, norbitals, neighbor count, sizes of indexing arrays, etc.).
- `fc.get_HS_sparse(dims)` → `firecore_get_HS_sparse` – Allocates and fills a `FireballData` object containing sparse neighbor‑block H and S matrices and all related indexing data (orbitals, neighbors, species, etc.).
- `fc.get_HS_k(kpoint_vec, norbitals)` → `firecore_get_HS_k` – Builds dense complex H(k) and S(k) matrices for a specified k‑vector from the internal sparse representation.

### 7.9 High‑Level Convenience

- `fc.run_nonSCF(atomType, atomPos)` – (Python helper combining several calls: `preinit` → `init` → `assembleH` → `solveH` → `updateCharges`) Performs a minimal, mostly non‑iterative run for testing or quick checks, returning the number of orbitals and the SCF residual `sigma`.

---

This document is intentionally high‑level; for parameter details or advanced usage, refer to the source of `pyBall/FireCore.py`, the Fortran modules in `fortran/`, and the concrete examples in `tests/pyFireball`.
