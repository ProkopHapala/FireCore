# Rigid Surface Potential & GridFF: Surface Interaction Evaluation Systems

## Overview

This document reviews the GridFF (Grid Force Field) and rigid surface potential evaluation systems in FireCore. These systems compute molecule-substrate interactions for AFM simulation, surface adsorption studies, and molecular assembly. Multiple implementations exist across different approaches:

- **GridFF**: Precomputed B-spline substrate potentials for fast GPU evaluation
- **XYZ Rigid**: On-the-fly pairwise evaluation using substrate atom coordinates
- **Surface Sampling**: Height map extraction via iso-surface scanning
- **Assembly**: Rigid-body pose generation and clash scoring

**Note**: The FoldedAtomicFunctions codemap (c9fc44a7...) was not available for this review. This document focuses on the 4 provided codemaps.

## Codemaps Analyzed

- Interactive GridFF Scanning: PTCDA-on-CaF2 Constrained Relaxation System https://windsurf.com/codemaps/99d506e2-223b-4ae7-bb60-8c2498fedfb9-8796fe608a7d71c1 
- Surface Potential Evaluation: GridFF B-spline and XYZ Rigid Kernels https://windsurf.com/codemaps/2a639fae-c9cb-407a-9d45-7b806c90c749-8796fe608a7d71c1
- FoldedAtomicFunctions: Surface Potential Basis Fitting System https://windsurf.com/codemaps/c9fc44a7-57a2-47c5-906f-886fa301ccc7-8796fe608a7d71c1
- Molecule-on-Surface Systems: GridFF, XYZ Scanning, Surface Sampling, and Assembly https://windsurf.com/codemaps/f8407e23-3a2e-41f1-abcf-9c15f3644c41-8796fe608a7d71c1
- Molecule-Substrate Interaction Energy Scanning: Assembly, GUI Placement, Force Fields & Surface Evaluation https://windsurf.com/codemaps/38bd3cb6-31c0-45b6-9e09-fda94257999c-8796fe608a7d71c1

## GridFF Subsystem

### Purpose

GridFF precomputes substrate electrostatic and non-electrostatic potentials on a 3D grid, then uses B-spline interpolation for fast GPU force evaluation during MD/relaxation.

### Implementation Locations

**Python Generation:**
- `pyBall/tests/ocl_GridFF_new.py` - GridFF generation pipeline
- `tests/tMMFF/run_test_GridFF_CaF2.py` - Test script for CaF2 substrate

**GPU Evaluation:**
- `pyBall/OCL/MolecularDynamics.py` - GPU initialization and kernel launch
- `cpp/common_resources/cl/relax_multi.cl` - OpenCL kernels for B-spline interpolation

### GridFF Generation Pipeline

**Entry Point:** `test_gridFF_ocl()` in `ocl_GridFF_new.py`

**Steps:**
1. Load substrate XYZ into AtomicSystem
2. Define grid parameters: origin `g0`, shape `ns = (nx, ny, nz)`
3. Generate Pauli/London grids via OpenCL kernel `make_MorseFF` (Morse potential on 3D grid)
4. Generate Coulomb grid via Ewald summation for slab geometry (`makeCoulombEwald_slab`)
5. Apply B-spline basis (`Bspline_basis5`)
6. Save to disk as `Bspline_PLQd.npy` (shape: nz,ny,nx,3 for Pauli/London/Coulomb channels)

**Key Files:**
- `ocl_GridFF_new.py:600` - Main generation function
- `GridFF.cl:201` - Morse potential kernel
- `GridFF.cl:301` - Ewald Coulomb kernel

### GridFF GPU Evaluation

**Initialization:** `MolecularDynamics.initGridFF()`

**Steps:**
1. Store grid metadata: `grid_ns` (shape), `grid_invStep` (1/voxel_size), `grid_p0` (origin), `GFFParams` (r_damp, alpha_morse)
2. Transpose B-spline data to (z,y,x,channels) layout for OpenCL 3D texture
3. Allocate 3D texture on GPU: `cl.Image(READ_ONLY, RGBA, FLOAT)`
4. Upload B-spline PLQH data via `cl.enqueue_copy()`
5. Generate kernel argument list

**Key Files:**
- `MolecularDynamics.py:1343` - GridFF initialization
- `MolecularDynamics.py:1372` - Data transpose for texture
- `MolecularDynamics.py:1376` - GPU texture allocation

### GridFF Force Evaluation Kernel

**Kernel:** `getNonBond_GridFF_Bspline_tex` in `relax_multi.cl:3263`

**Per-atom workflow:**
1. Load atom position from global memory
2. Transform world coordinates to grid fractional coords: `u = (posi - grid_p0) * grid_invStep`
3. Call `fe3d_pbc_comb_tex()` for tricubic B-spline interpolation with PBC
   - Apply PBC to X indices via `choose_inds_pbc_3()`
   - Interpolate YZ slices via `fe2d_comb_tex()` × 4
   - 1D interpolation via `fe1Dcomb_tex()` using `read_imagef(texture)`
4. Scale grid-space forces to physical units: `fg.xyz *= -grid_invStep.xyz`
5. Write result to forces buffer

**Key Files:**
- `relax_multi.cl:3263` - Main GridFF kernel
- `relax_multi.cl:3184` - 3D B-spline interpolation function
- `relax_multi.cl:2944` - PBC index handling
- `relax_multi.cl:3126` - 1D B-spline basis

## XYZ Rigid Subsystem

### Purpose

XYZ rigid evaluation computes molecule-substrate interactions on-the-fly using pairwise potentials (LJ, Coulomb, H-bond) with explicit substrate atom coordinates. No precomputation required.

### Implementation Locations

**Python Interface:**
- `pyBall/OCL/InteractionEnergy.py` - Rigid scanning interface
- `tests/tMMFF/test_interaction_scan.py` - 1D z-scan test

**GPU Kernel:**
- `cpp/common_resources/cl/relax_multi.cl` - `getSurfMorse` kernel

### XYZ Rigid 1D Scan (Z-Approach)

**Entry Point:** `scanner.scan_1D()` in `test_interaction_scan.py`

**Workflow:**
1. Create InteractionScanner instance
2. Loop over z positions in `z_array`
3. For each z:
   - Update probe z position
   - Call `run_getSurfMorse()` to launch GPU kernel
   - Download forces, extract energy from `forces[0].w`
4. Return energy curve E(z)

**Key Files:**
- `test_interaction_scan.py:100` - Scanner creation
- `InteractionEnergy.py:250` - Z-position loop
- `InteractionEnergy.py:260` - Energy evaluation

### XYZ Rigid Batch Evaluation

**Entry Point:** `eval_rigid_getSurfMorse()` in `MolecularDynamics.py`

**Workflow:**
1. Wrap rigid transforms with PBC using surface lattice vectors
2. Chunk loop for batch processing:
   - Upload rigid transforms (apply rotation/translation to molecule template)
   - Zero force buffer
   - Launch `getSurfMorse` kernel (no relaxation, just energy evaluation)
   - Download forces from GPU
   - Convert forces to total energy (sum force.w over atoms, negate)
3. Return energy dictionary

**Key Files:**
- `MolecularDynamics.py:1630` - Batch evaluation entry
- `MolecularDynamics.py:1635` - PBC wrapping
- `MolecularDynamics.py:1596` - Transform upload
- `MolecularDynamics.py:739` - Kernel launch

### getSurfMorse Kernel

**Kernel:** `getSurfMorse` in `relax_multi.cl:4013`

**Workflow:**
1. Loop over substrate atoms
2. Compute LJ/Coulomb/HBond interactions via `getLJQH_SR()`
3. Accumulate forces and energy
4. Write results to buffer

**Key Files:**
- `relax_multi.cl:4013` - Kernel entry
- `relax_multi.cl:4070` - Substrate atom loop
- `relax_multi.cl:4080` - LJ/Coulomb force computation

## Surface Sampling Subsystem

### Purpose

Extract 2D height maps z(x,y) representing iso-surfaces or first minima of molecule-substrate interaction potentials.

### Current Slow Python Path

**File:** `pyBall/SurfaceSampling.py`

**Problem:** Nested Python loops over (x,y) grid points, not vectorized

**Workflow:**
1. `sample_surface_map()` orchestrator
2. Outer loop over y: `for iy, y in enumerate(ys)`
3. Inner loop over x: `for ix, x in enumerate(xs)`
4. For each (x,y):
   - Generate z-scan array
   - Call `backend.eval_component()` (GridFF: trilinear in Python, XYZ: upload/download nz times)
   - Refine first minimum in Python loop
   - Store result
5. Return SurfaceResult

**Key Files:**
- `SurfaceSampling.py:378` - Main orchestrator
- `SurfaceSampling.py:396` - Y loop
- `SurfaceSampling.py:397` - X loop
- `SurfaceSampling.py:152` - GridFF trilinear interpolation (slow)
- `SurfaceSampling.py:311` - Minimum search (slow)

### GPU Surface Height Map Extraction

**Entry Point:** `eval_surface_iso()` in `SurfaceSampling.py`

**GPU Kernel:** `getSurfaceIsoGridFF` in `relax_multi.cl:2960`

**Workflow:**
1. Python: Build GridFF backend
2. Python: Call `eval_surface_iso()`
3. Python: Setup kernel args, launch with `global_size=(nx,ny)`
4. GPU Kernel:
   - Each workitem handles one (x,y) point
   - Loop over z to find iso-surface crossing or first minimum
   - Store height zh and color to output buffer
5. Python: Download results via `enqueue_copy()`

**Key Files:**
- `render_surface_iso.py:150` - Backend building
- `SurfaceSampling.py:400` - GPU evaluation request
- `MolecularDynamics.py:1630` - Kernel launch
- `relax_multi.cl:2960` - Surface iso kernel
- `relax_multi.cl:2975` - Get (x,y) from workitem ID
- `relax_multi.cl:3010` - Z-scan loop
- `relax_multi.cl:3040` - Write surface point

## Assembly Subsystem

### Purpose

Generate thousands of molecular poses (rotations + translations) and score them via GPU clash detection for molecular assembly on surfaces.

### Implementation Locations

**Python:**
- `tests/tMMFF/test_assembly.py` - CLI driver
- `pyBall/OCL/Assembly.py` - GPU interface

**GPU Kernel:**
- `cpp/common_resources/cl/Assembly.cl` - Assembly scoring kernel

### Assembly Pipeline

**Entry Point:** `test_assembly.py` main()

**Workflow:**
1. Parse args, load molecule
2. Generate transform buffer:
   - Build rotation grid (full3d/inplane/tilt) via `generate_rotations()`
   - Build translation grid
   - Apply 6-fold symmetry
   - Tile PBC lattice replicas
3. Z-span prefiltering (host-side)
4. Launch GPU evaluation via `AssemblyOCL.evaluate_packing()`
5. Filter & export results

**Key Files:**
- `test_assembly.py:199` - CLI driver
- `test_assembly.py:291` - Transform generation
- `test_assembly.py:44` - Rotation generation
- `test_assembly.py:355` - Symmetry application
- `test_assembly.py:358` - PBC tiling
- `test_assembly.py:465` - GPU evaluation

### GPU Assembly Scoring

**Kernel:** `evaluate_packing_3d` in `Assembly.cl:46`

**Workflow:**
1. Load reference molecule (mol 0)
2. Loop over replicas (mol 1..N)
3. Apply transforms to replicas
4. Clash check: compare atom-pair distance to sum of radii, accumulate quadratic overlap penalty
5. Workgroup reduction
6. Thread 0 writes total penalty to results buffer

**Key Files:**
- `Assembly.py:22` - GPU evaluation entry
- `Assembly.py:29` - Transform upload
- `Assembly.py:56` - Kernel enqueue
- `Assembly.cl:46` - Kernel entry
- `Assembly.cl:72` - Load ref molecule
- `Assembly.cl:90` - Loop replicas
- `Assembly.cl:126` - Clash detection
- `Assembly.cl:155` - Workgroup reduction
- `Assembly.cl:171` - Write result

### Assembly Export

**Workflow:**
1. Host-side filtering by thresholds (z/clash/min-dist)
2. Sort by composite score
3. Loop over exported configs:
   - Select replica subset
   - Call GPU emit via `AssemblyOCL.emit_configuration()`
   - Launch `emit_configuration_xyz` kernel
4. Write XYZ/mol2 files

**Key Files:**
- `test_assembly.py:483` - Filtering
- `test_assembly.py:495` - Sorting
- `test_assembly.py:596` - Replica selection
- `test_assembly.py:597` - GPU emit
- `Assembly.py:62` - Emit configuration
- `Assembly.cl:21` - Emit kernel
- `Assembly.cl:29` - Thread per atom
- `Assembly.cl:43` - Write transformed atom

## Interactive GUI Systems

### ExplorerVisPy

**Purpose:** Interactive scanning and visualization of molecule-surface interactions

**Implementation:** `pyBall/ExplorerVisPy.py`

**Initialization:**
1. `__init__()` constructor
2. `_init_scanner()` - Create InteractionScanner instances (reference, fast GPU, folded)
3. `_reload_system()`:
   - Load molecule via `scanner.load_molecule_xyz()`
   - Load substrate via `scanner.load_substrate_xyz()`
   - Initialize GPU scanners (MolecularDynamics, fast_scanner.set_surface())
4. `_update_visuals()`, `_eval_current_pose()`
5. `_build_ui()` - Setup VisPy canvas & controls

**Key Files:**
- `ExplorerVisPy.py:51` - Constructor
- `ExplorerVisPy.py:92` - Scanner initialization
- `ExplorerVisPy.py:328` - System reload
- `ExplorerVisPy.py:337` - Molecule loading
- `ExplorerVisPy.py:356` - GPU scanner creation
- `ExplorerVisPy.py:358` - Substrate configuration

**User Interaction:**
- `_on_scan_surface()` callback - User clicks "Scan Surface" button
- Surface sampling via `build_surface_map()`
- Backend evaluation via `backend.eval_surface_iso_*()`
- GPU kernel execution
- Visualization update via `make_surface_mesh()` and `surface_visual.set_data()`

**Key Files:**
- `ExplorerVisPy.py:504` - Surface scan trigger
- `SurfaceSampling.py:400` - Surface map building
- `ExplorerVisPy.py:550` - Mesh generation
- `ExplorerVisPy.py:560` - Visualization update

### Constraint-based Relaxation

**Purpose:** Interactive relaxation with harmonic constraints on anchor atoms

**Workflow (Single MD Step):**
1. `run_cleanForceMMFFf4()` - Zero force buffers
2. `run_getMMFFf4_rot()` - Evaluate bonded MMFF forces
3. `run_getNonBond_GridFF_Bspline_ex2()` - Evaluate GridFF surface forces
4. `run_updateAtomsMMFFf4()` - Integrate with constraints
   - Kernel reads constraint target: `float4 cons = constr[iaa]`
   - Compute spring force: `fe.xyz += fc`
   - Velocity Verlet integration

**Key Files:**
- `MolecularDynamics.py:859` - Zero forces
- `MolecularDynamics.py:839` - Bonded forces
- `MolecularDynamics.py:805` - GridFF forces
- `MolecularDynamics.py:848` - Integration with constraints
- `relax_multi.cl:1089` - Update kernel
- `relax_multi.cl:1208` - Read constraint
- `relax_multi.cl:1213` - Apply constraint force

**Constraint Update Flow:**
1. VisPy canvas with keyboard input
2. Allocate GPU buffers: `constr` (target pos), `constrK` (stiffness)
3. Upload via `toGPU()` → `cl.enqueue_copy()`
4. Kernel execution: read `cons` and `cK`, compute `fc = (cons - pe) * cK`, accumulate

**Key Files:**
- `ExplorerVisPy.py:103` - Keyboard enable
- `MolecularDynamics.py:265` - Constraint buffer
- `MolecularDynamics.py:266` - Stiffness buffer
- `OpenCLBase.py:183` - GPU upload
- `relax_multi.cl:1210` - Read stiffness
- `relax_multi.cl:1212` - Compute force

## MMFF Surface MD

### Purpose

Run molecular dynamics on surfaces with MMFF bonded forces and GridFF surface forces

### Implementation

**Test Script:** `tests/tMMFF/test_ditetraceno_surface.py`

**Workflow:**
1. Create MMFFL builder
2. Build topology:
   - Assign atom types
   - Detect bonds/angles
   - Convert angles → harmonic bonds
   - Add pi/epair dummies
3. Create GPU MD harness (MolecularDynamics)
4. Upload substrate GridFF via `md.initGridFF()`
5. MD loop:
   - Clean forces
   - Evaluate nonbonded forces (Pauli/LJ + Coulomb)
   - Add surface interaction (Hamaker or Morse)
   - Add bonded forces (getMMFFf4)
   - Integrate positions (Velocity Verlet + damping)
   - Download & write trajectory

**Key Files:**
- `test_ditetraceno_surface.py:94` - Entry point
- `test_ditetraceno_surface.py:50` - MMFF builder
- `MMFFL.py:200` - Topology building
- `MMFFL.py:210` - Atom typing
- `MMFFL.py:220` - Bond/angle detection
- `MMFFL.py:250` - Angle conversion
- `MMFFL.py:300` - Pi/epair dummies
- `test_ditetraceno_surface.py:155` - GPU harness
- `test_ditetraceno_surface.py:158` - GridFF upload
- `test_ditetraceno_surface.py:192` - MD loop
- `test_ditetraceno_surface.py:202` - Nonbonded forces
- `test_ditetraceno_surface.py:214` - Surface forces
- `test_ditetraceno_surface.py:220` - Bonded forces
- `test_ditetraceno_surface.py:231` - Integration

## MMFF Topology Construction

### Purpose

Convert molecular XYZ coordinates into GPU-ready MMFF topology with bonds, angles, and parameters

### Implementation

**Entry Point:** `MMFF.toMMFFsp3_loc()` in `pyBall/OCL/MMFF.py`

**Workflow:**
1. Load AtomTypes.dat parameters via `MMparams.read_atom_types()`
2. Assign atom properties via `initAtomProperties()`:
   - Determine pi-orbitals (npi)
   - Determine electron pairs (nep)
   - Classify nodes vs caps
3. Allocate topology arrays via `MMFF.realloc()`:
   - `apos[nvecs, 4]` - positions
   - `neighs[natoms, 4]` - connectivity
   - `REQs[natoms, 4]` - nonbonded parameters
   - `bKs, bLs[nnode, 4]` - bond parameters
4. Build bonds and angles:
   - Detect connectivity
   - Assign bond parameters
   - Generate pi-orbital positions

**Key Files:**
- `MMFF.py:377` - Topology builder entry
- `MMparams.py:126` - Load parameters
- `MMparams.py:155` - Atom type creation
- `MMFF.py:51` - Atom properties
- `MMFF.py:68` - Pi-orbitals
- `MMFF.py:69` - Electron pairs
- `MMFF.py:79` - Node/cap classification
- `MMFF.py:182` - Array allocation
- `MMFF.py:203` - Positions
- `MMFF.py:206` - Connectivity
- `MMFF.py:211` - REQs
- `MMFF.py:213` - Bond parameters

### MMFF System Packing

**Purpose:** Transfer MMFF topology data from CPU to GPU

**Entry Point:** `pack_system(iSys, mmff)` in `MolecularDynamics.py`

**Workflow:**
1. Upload atomic positions
2. Upload bond connectivity
3. Upload bond parameters
4. Upload angle parameters
5. Upload nonbonded REQs
6. Upload pi-orbital stiffness

**Key Files:**
- `MolecularDynamics.py:1639` - Pack entry
- `MolecularDynamics.py:1648` - Positions
- `MolecularDynamics.py:1650` - Connectivity
- `MolecularDynamics.py:1654` - Bond parameters
- `MolecularDynamics.py:1651` - REQs
- `MolecularDynamics.py:1655` - Angle parameters
- `MolecularDynamics.py:1656` - Pi stiffness

## Path Optimization with Substrate

### Purpose

Batch replica optimization with causal tethers and substrate GridFF forces

### Implementation

**File:** `pyBall/OCL/ManipulationPathOpt.py`

**Initialization:**
1. Build base MMFF topology via `mm.toMMFFsp3_loc()`
2. Create MolecularDynamics harness with `md.realloc(nSystems=n_rep)`
3. Pack MMFF to all replicas via `md.pack_system()`
4. Set substrate GridFF via `md.initGridFF()`
5. Initialize replica states:
   - Set tip trajectory
   - Build apos_all with shifts
   - Upload to GPU via `md.toGPU('apos', apos_all)`
6. Set causal tethers:
   - Build sysneighs array
   - Build sysbonds array
   - Upload to GPU
7. Run relaxation loop:
   - Apply tip Morse force via `prog.getTipMorse()`

**Key Files:**
- `ManipulationPathOpt.py:13` - Constructor
- `ManipulationPathOpt.py:25` - MMFF build
- `ManipulationPathOpt.py:37` - GPU harness
- `ManipulationPathOpt.py:42` - Batch allocation
- `ManipulationPathOpt.py:46` - Pack MMFF
- `ManipulationPathOpt.py:73` - GridFF upload
- `ManipulationPathOpt.py:129` - Tip trajectory
- `ManipulationPathOpt.py:144` - Position shifts
- `ManipulationPathOpt.py:155` - GPU upload
- `ManipulationPathOpt.py:162` - Causal neighbors
- `ManipulationPathOpt.py:176` - Causal bonds
- `ManipulationPathOpt.py:208` - Relax loop
- `ManipulationPathOpt.py:185` - Tip Morse

## Probe Particle Configuration

### Purpose

Load probe particle parameters from AtomTypes.dat and convert to PLQH format for GPU kernels

### Implementation

**File:** `pyBall/SurfaceSampling.py`

**Workflow:**
1. `ProbeParticle.__init__(name, charge, ...)`:
   - Lookup atom types via `get_atom_types()`
   - Load element types from ElementTypes.dat
   - Load atom types from AtomTypes.dat
   - Cache in `_ATOM_TYPES` global
   - Extract `at.RvdW → self.R0`
   - Extract `at.EvdW → self.E0`
   - Store charge, alpha, H
2. `ProbeParticle.plqh()` method:
   - Call `getPLQH(R0, E0, a, Q, H)`
   - Return PLQH array for GPU kernels

**Key Files:**
- `SurfaceSampling.py:32` - Probe initialization
- `SurfaceSampling.py:24` - Atom type lookup
- `SurfaceSampling.py:17` - Element types
- `SurfaceSampling.py:18` - Atom types
- `SurfaceSampling.py:27` - Cache
- `SurfaceSampling.py:40` - RvdW extraction
- `SurfaceSampling.py:41` - EvdW extraction
- `SurfaceSampling.py:38` - Charge storage
- `SurfaceSampling.py:48` - PLQH conversion

## Comparison: GridFF vs XYZ Rigid

### GridFF (Precomputed B-spline)

**Advantages:**
- Fast GPU evaluation via tricubic B-spline interpolation
- Suitable for repeated evaluations (MD, relaxation, scanning)
- Handles electrostatics via Ewald summation during generation
- Efficient for large substrates (precompute once, evaluate many times)

**Disadvantages:**
- Requires offline preprocessing step
- Grid resolution affects accuracy
- Memory intensive (3D grid storage)
- Less flexible for dynamic substrate changes

**Use Cases:**
- MD relaxation on surfaces
- Interactive scanning with real-time feedback
- Surface height map extraction
- Batch pose evaluation

### XYZ Rigid (On-the-fly Pairwise)

**Advantages:**
- No preprocessing required
- Exact pairwise evaluation (no grid interpolation error)
- Flexible for dynamic substrates
- Lower memory footprint

**Disadvantages:**
- Slower for large substrates (O(N_substrate * N_molecule) per evaluation)
- GPU kernel complexity for pairwise loops
- Less suitable for real-time interactive use

**Use Cases:**
- Small substrates
- Dynamic substrate geometries
- One-time evaluations
- Reference calculations for validation

## Status and Consolidation Opportunities

### Active Implementations

1. **GridFF B-spline with texture sampling** - Active, used in ExplorerVisPy and surface sampling
2. **XYZ rigid getSurfMorse** - Active, used for 1D scans and batch evaluation
3. **GPU surface height map extraction** - Active, but slow Python path still exists
4. **Assembly scoring** - Active, GPU-accelerated clash detection

### Deprecated/Experimental

- Slow Python surface sampler in `SurfaceSampling.py` - Should be replaced with GPU kernel
- Multiple duplicate surface sampling backends - Need consolidation

### Consolidation Priorities

1. **Replace Python surface sampler** with GPU `getSurfaceIsoGridFF` kernel
2. **Unify GridFF initialization** across different modules (currently duplicated)
3. **Consolidate surface sampling backends** - GridFF, XYZ, and Folded should share common interface
4. **Standardize probe particle configuration** - Currently scattered across modules

## FoldedAtomicFunctions Subsystem

### Purpose

FoldedAtomicFunctions provides an alternative to GridFF by fitting periodic surface potentials to compact basis functions (plane waves × exponentials) for fast evaluation. Unlike GridFF's B-spline interpolation on a 3D grid, folded functions use analytical basis expansion.

### Implementation Locations

**Python Framework:**
- `doc/py/FoldedAtomicFunctions/FoldedAtomicFunctions.py` - Core fitting framework
- `doc/py/FoldedAtomicFunctions/FoldedAtomicFunction.md` - Design document
- `doc/py/FoldedAtomicFunctions/FoldedAtomicFunction_rigorous_test_plan.md` - Test methodology

**GPU Integration:**
- `pyBall/OCL/MolecularDynamics.py` - Folded basis fitting and evaluation
- `cpp/common_resources/cl/relax_multi.cl` - `getSurfFolded` kernel

### Folded Basis Fitting Pipeline

**Entry Point:** `fit_folded_surface_basis()` in `MolecularDynamics.py`

**Steps:**
1. Generate reference potential grid using embedding-corrected Coulomb + Morse
2. Fit periodic P/L/Q components to separable basis:
   - X-direction: Cosine plane waves `cos(2πnx/L)` for n=0,1,...,nu
   - Z-direction: Exponential decay `exp(-αz)` with α=[α1,α2,...,nzbasis]
3. Solve linear least squares with L2 regularization
4. Store coefficients for GPU evaluation

**Key Files:**
- `MolecularDynamics.py:1029` - Folded basis fitting
- `MolecularDynamics.py:1034` - Folded evaluation (2D NDRange)
- `FoldedAtomicFunctions.py:226` - Basis generation
- `FoldedAtomicFunctions.py:334` - Least squares fitting

### Folded Basis Evaluation

**Kernel:** `getSurfFolded` in `relax_multi.cl`

**Workflow:**
1. Load fitted coefficients from global memory
2. Evaluate basis functions at probe position
3. Sum weighted basis functions to get potential
4. Add analytic macro Coulomb correction (dipole/quadrupole)
5. Return total potential

**Key Files:**
- `MolecularDynamics.py:1034` - GPU evaluation entry
- `relax_multi.cl` - Folded evaluation kernel
- `MolecularDynamics.py:568` - Macro correction (surface moments)

### Comparison: Folded vs GridFF

### FoldedAtomicFunctions

**Advantages:**
- Analytical basis functions (no grid interpolation error)
- Compact coefficient storage (few hundred floats vs 3D grid)
- Naturally periodic (plane waves satisfy PBC)
- Flexible basis family (can swap exponential for other functions)
- Component-wise fitting (P/L/Q separate)

**Disadvantages:**
- Requires fitting step (offline preprocessing)
- Basis expressivity limits accuracy for sharp Coulomb fields
- More complex parameter tuning (nu, nv, nzbasis, regularization)
- Coulomb term difficult to fit with smooth exponentials

**Use Cases:**
- Systems where basis functions match potential shape
- When memory footprint is critical
- When analytical derivatives are needed
- Component-wise analysis (P/L/Q separation)

### Electrostatic Continuum Embedding

**Purpose:** Correct finite supercell dipole/quadrupole artifacts for periodic slab electrostatics

**Implementation:**
- `pyBall/OCL/InteractionEnergy.py` - Macro correction calculation
- `pyBall/Ewald2D.py` - Rigorous 2D Ewald solution (alternative approach)
- `cpp/common_resources/cl/relax_multi.cl` - `getMacroRectLayers` kernel

**Key Insight:** The long-range electrostatic problem was solved via Ewald2D, making folded basis fitting unnecessary for Coulomb. Folded basis is now used primarily for short-range P/L terms.

### AFM Rigid Body Dynamics

**Purpose:** GPU-accelerated rigid body simulation for AFM tip-molecule interaction

**Implementation Locations:**
- `pyBall/OCL/RigidBodyDynamics.py` - Python wrapper
- `cpp/common_resources/cl/Rigid.cl` - OpenCL kernels
- `tests/tMMFF/test_rigid_gridff_ptcda_batch.py` - Batch relaxation
- `tests/tMMFF/test_rigid_gridff_surface.py` - Single-body relaxation

**Architecture:**
- One workgroup per rigid body (32 threads)
- Max 128 atoms per body
- Two kernels:
  - `rigid_body_dynamics_kernel`: Anchor support, no GridFF
  - `rigid_body_gridff_kernel`: GridFF support, no anchor (critical gap)

**Critical Gap:** Need to merge anchor support into GridFF kernel for AFM simulation

### Surface Potential Visualization

**Purpose:** Render z(x,y) height maps of molecule-substrate interaction surfaces

**Implementation:**
- `pyBall/SurfaceSampling.py` - Surface sampling interface
- `tests/tMMFF/render_surface_iso.py` - Headless PNG rendering
- `pyBall/ExplorerVisPy.py` - GUI integration

**Modes:**
- GridFF: Precomputed B-spline evaluation
- XYZ Reference: Slow CPU reference
- XYZ Fast GPU: On-the-fly pairwise evaluation
- XYZ Folded GPU: Folded basis evaluation

**Key Issues:**
- Coulomb scaling discrepancy between GridFF and XYZ (factor of 14.39 eV·Å)
- Folded Coulomb accuracy limited by basis expressivity
- Need proper periodic boundary conditions and macro corrections

## Key File Locations Summary

### GridFF Generation
- `pyBall/tests/ocl_GridFF_new.py` - Main generation pipeline
- `tests/tMMFF/run_test_GridFF_CaF2.py` - Test script
- `GridFF.cl` - OpenCL kernels for potential generation

### GridFF Evaluation
- `pyBall/OCL/MolecularDynamics.py` - GPU initialization and kernel launch
- `cpp/common_resources/cl/relax_multi.cl` - B-spline interpolation kernels

### XYZ Rigid
- `pyBall/OCL/InteractionEnergy.py` - Rigid scanning interface
- `tests/tMMFF/test_interaction_scan.py` - 1D scan test
- `cpp/common_resources/cl/relax_multi.cl` - getSurfMorse kernel

### Surface Sampling
- `pyBall/SurfaceSampling.py` - Python sampling layer (slow path)
- `cpp/common_resources/cl/relax_multi.cl` - GPU surface iso kernel
- `tests/tMMFF/render_surface_iso.py` - Headless surface rendering

### Assembly
- `tests/tMMFF/test_assembly.py` - CLI driver
- `pyBall/OCL/Assembly.py` - GPU interface
- `cpp/common_resources/cl/Assembly.cl` - Assembly scoring kernel

### GUI
- `pyBall/ExplorerVisPy.py` - Interactive scanning GUI
- `pyBall/SequencePlacerVisPy.py` - Sequence placement GUI
- `pyBall/MolecularPlacerVisPy.py` - Symmetric placement GUI

### MMFF
- `pyBall/OCL/MMFF.py` - MMFF topology builder
- `pyBall/OCL/MMFFL.py` - Linearized MMFF (angles → bonds)
- `pyBall/OCL/MMparams.py` - Parameter loading

### Path Optimization
- `pyBall/OCL/ManipulationPathOpt.py` - Batch replica optimization with causal tethers

## FDBM/GridFF Fitting

### Purpose
Fit GridFF potential parameters from Full Density Based Model (FDBM) using DFT interaction energy data

### Implementation
- `doc/Topics/OnSurfaceAssembly/GridFF_FDBM_Fitting.md` - Documentation for FDBM fitting methodology and REQ/PLQ convention
- `tests/tMMFF/test_fdbm_fit_gridff_mock.py` - Mock test for FDBM fitting with CLI options for config generation (xy-range, z-range)
- `tests/tMMFF/test_gridff_alignment.py` - GridFF alignment verification test
- `doc/Topics/OnSurfaceAssembly/GridFF_atoms_alignment.md` - Documentation for GridFF atom alignment issues and solutions

## Electrostatics and Surface Potential

### Purpose
Rigorous electrostatic evaluation for periodic slab geometries using Ewald summation

### Implementation
- `pyBall/OCL/SurfaceEwald.py` - Ewald summation for surface electrostatics
- `tests/tMMFF/test_electrostatics_comparison.py` - Test comparing electrostatics implementations
- `cpp/common_resources/cl/Surface.cl` - OpenCL kernels for surface potential evaluation
- `doc/Topics/OnSurfaceAssembly/FoldedSubstratePotential_OpenCL.md` - Documentation for folded substrate potential with OpenCL implementation

## Surface Utilities

### Purpose
Diagnostic plotting and utility functions for surface potential analysis

### Implementation
- `pyBall/OCL/Surface_utils.py` - Surface utility functions including diagnostic plotting with reusable functions (plot_xz_slice, plot_1d_profile)
- `tests/tMMFF/test_folded_fit_nacl1x1.py` - Test for folded potential fitting on NaCl substrate with 4x4x4 basis and Ewald2D Coulomb reference

## Folded Substrate Potential Documentation

### Purpose
Implementation plan and design documentation for folded atomic functions with OpenCL kernels

### Implementation
- `doc/Topics/OnSurfaceAssembly/FoldedSubstratePotential_OpenCL.md` - Detailed implementation plan for folded basis projection (Pauli/vdW/Electrostatics) with separable basis functions, Ewald2D integration, and GPU kernel optimization patterns
