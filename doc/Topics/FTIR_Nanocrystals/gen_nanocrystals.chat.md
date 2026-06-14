---
description: Silicon nanocrystal generator (scripts/gen_nanocrystals.mjs)
---

# Silicon Nanocrystal Generator

Reference for developers and students on how to generate silicon nanocrystals with configurable cutting planes, pruning, capping, and surface bridge operations (collapse/insert) using `scripts/gen_nanocrystals.mjs`.

## What it does

- Builds a replicated Si supercell from CIF (`cpp/common_resources/crystals/Si-sym.cif` by default).
- Applies plane cuts to shape the nanocrystal (e.g., {111} and {100}).
- Recalculates bonds using MM parameters (ElementTypes/AtomTypes/BondTypes/AngleTypes).
- Prunes undercoordinated Si (iteratively).
- Adds caps (H) to dangling bonds.
- Optionally performs surface bridge insertion (adds SiH2 bridges) and collapse (removes SiH2 bridges, bonds the neighbors).
- Outputs per-sample `.mol2` and optional stacked multi-frame `.xyz` for Jmol/visual inspection.

## CLI (key options)

```
node scripts/gen_nanocrystals.mjs [options]
```

- Geometry & planes:
  - `--cif <path>`: CIF file (default Si-sym.cif)
  - `--nx-range a,b` `--ny-range a,b` `--nz-range a,b`: replication ranges (single value allowed)
  - `--centered 0|1`: center replication
  - `--planeTemplates a111,a100,...`: plane families
  - `--planeSymC <float>`: symmetry scaling for templates
  - `--planeMode ang|frac`: plane definition mode
  - `--planeCScale <float>` `--planeCJitter <float>`: scale/jitter plane offsets

- Bonding & pruning:
  - `--bondFactor <float>` `--defaultRcut <float>`: bonding thresholds
  - `--minSiDegree <int>` `--pruneMaxIter <int>`: prune undercoordinated Si iteratively

- Caps:
  - `--caps H|0|none`: add caps (H) or disable

- Surface bridges:
  - `--collapseProb <0..1>`: probability to collapse a selected surface SiH2 bridge (remove bridge, bond neighbors)
  - `--insertProb <0..1>`: probability to insert a SiH2 bridge on each surface Si–Si bond candidate
  - `--collapseAll 0|1`: legacy carbon bridge collapse

- Energetics (simple counting model):
  - `--E_SiH2` `--E_SiH3` `--E_bare` `--E_bridge` `--muH`

- Output:
  - `--outDir <path>`: output directory
  - `--prefix <str>`: file prefix
  - `--samples <int>`: number of samples (or use `--maxFiles`)
  - `--statsCsv <path>`: append stats (counts, E)
  - `--stackedXyzOut <path>`: write all frames into one multi-frame XYZ
  - `--seed <int>`: deterministic RNG

## Current surface-bridge logic (Si)

- Surface filter: `nSi < 4` (targets surface atoms).
- Insertion (surface-only):
  - Enumerate surface Si–Si bonds; with probability `insertProb`, call `insertBridge` (adds SiH2) and promote the new atom to Si (Z=14).
- Collapse:
  - Select surface Si with 2 heavy neighbors and ≥2 H (SiH2-like) and collapse with probability `collapseProb` (remove bridge + H, bond neighbors).
- No bond rebuild after bridge ops (to preserve newly added R1–R2 bonds).

## Usage examples

```bash
# Generate a set of nanocrystals
node scripts/gen_nanocrystals.mjs --planeTemplates a111,a100 --planeSymC 2.0 --samples 5
```

## Linear Response Vibration Spectroscopy (FTIR)

We have implemented a super-fast linear response method for computing vibration spectra of nanocrystals using the MMFF forcefield. This method avoids full diagonalization and instead computes the dynamic response to an external driving force.

### Key Components

- **`cpp/libs/Molecular/MMFF_lib.cpp`**: Implements `getHessian3Nx3N` which calculates the full 3Nx3N dynamical matrix using finite differences of forces.
- **`pyBall/FTIR.py`**: Reusable Python module for:
    - Building the mass-weighted Hessian.
    - Projecting out rigid body (translation/rotation) degrees of freedom using a penalty shift method.
    - Solving the linear response equation $A(\omega)u = f_{ext}$ where $A(\omega) = K - (\omega + i\eta)^2 M$.
    - Exporting vibrational modes to `.xyz` format for Jmol visualization.
- **`tests/tMMFF/test_vibration_spectra.py`**: CLI tool to run the spectroscopy.

### CLI for Vibration Spectroscopy

```bash
cd tests/tMMFF/
python3 test_vibration_spectra.py --input adamantane.mol2 --fmax 10.0 --eta 0.05
```

Options:
- `-i, --input`: Input molecule file (.mol2 or .xyz).
- `--fmin, --fmax`: Frequency range for the scan.
- `--eta`: Damping parameter (controls peak broadening).
- `--save_modes`: Save eigenmodes to `eigenmodes.xyz`.
- `--save_responses`: Save response vectors at resonance to `response_x.xyz`.

### Results and Verification (Adamantane)

Test run on Adamantane ($C_{10}H_{16}$) with a bonding-only forcefield:
- **Spectrum**: The computed absorption spectrum shows sharp Lorentzian peaks.
- **Comparison**: Diagonalization of the mass-weighted Hessian yields exact eigenfrequencies which align perfectly with the peaks in the linear response spectrum.
- **Broadening**: Peak widths are controlled by the artificial damping $\eta$.
- **Visualization**: Displacement vectors saved in `eigenmodes.xyz` and `response_x.xyz` can be animated in Jmol using the `vector on` command, showing characteristic C-H and C-C stretching/bending modes.
- **Excitation**: Using a symmetric stretching excitation (proportional to displacement from COG) ensures that only internal modes are excited, avoiding large artifacts from rigid body motion at $\omega \approx 0$.


Fixed shape, random passivation:
```
node scripts/gen_nanocrystals.mjs \
  --samples 10 \
  --nx-range 2,2 --ny-range 2,2 --nz-range 2,2 \
  --planeTemplates a111,a100 --planeCScale 0.45 --planeCJitter 0 \
  --caps H --minSiDegree 2 \
  --collapseProb 0.3 --insertProb 0.2 \
  --seed 2 \
  --outDir OUT_fixedshape \
  --prefix fixed_si \
  --statsCsv OUT_fixedshape/stats.csv \
  --stackedXyzOut OUT_fixedshape/stacked.xyz
```

Smaller crystals for inspection:
```
node scripts/gen_nanocrystals.mjs \
  --samples 10 \
  --nx-range 1,2 --ny-range 1,2 --nz-range 1,2 \
  --planeTemplates a111,a100 --planeCScale 0.45 --planeCJitter 0.20 \
  --caps H --minSiDegree 2 \
  --collapseProb 0.3 --insertProb 0.2 \
  --seed 2 \
  --outDir OUT_nanocrystals_small10 \
  --prefix small_si \
  --stackedXyzOut OUT_nanocrystals_small10/stacked.xyz
```

## Troubleshooting / gotchas

- **Do not recalc bonds after bridge ops**: recalculating removes the newly formed R1–R2 bond; we removed those recalc calls in the bridge section.
- If inserts stay zero: check the debug log `surface Si-Si candidate bonds: N`. If N>0 but inserts=0, raise `insertProb` or move insertion earlier (before heavy pruning/capping). If N=0, loosen the surface filter or pruning.
- The script promotes inserted atoms to Si (Z=14) after `insertBridge` (which defaults to carbon geometry but works for SiH2 here).

## Development notes / lessons learned

- Root cause of missing R1–R2 bonds was bond recalculation after collapse; fixed by removing post-op recalc.
- Insertion initially failed due to overly tight filters and try/catch masking; simplified to deterministic candidate loop with surface-only filter and direct `insertBridge` calls.
- Debug logging (candidate count) is essential to see if the surface actually exposes insertable Si–Si bonds.
- Keep randomness explicit (`--seed`) when comparing passivation variants on fixed shapes.

---

# Vibration Spectroscopy via Green's Function Probing

Reference for computing vibration spectra (FTIR/Raman) of nanocrystals using linear response methods and Green's function probing, implemented in `doc/py/OrderN_QM/VibrationProbing.py`.

## Overview

The script implements **mechanical Green's function probing** for vibration spectra of spring-mass systems. It uses dipole-driven linear response to compute IR-active vibrational modes without explicit diagonalization of the dynamical matrix.

## Mathematical Framework

### Dynamic Stiffness

For a system with stiffness matrix K and mass matrix M, the frequency-dependent dynamic stiffness is:

```
A(ω) = K - (ω + iη)²M
```

where:
- ω = real driving frequency
- η = small damping parameter (shifts poles off real axis for stability)
- stabilize = small diagonal shift to prevent singularity

### Dipole-Driven Probing

Instead of explicit eigendecomposition, the method probes the system with a homogeneous electric field:

1. **Force**: `f_i = q_i * E` (same direction for all nodes)
2. **Linear Solve**: `A(ω) * U = f` for each frequency ω
3. **Dipole Response**: `Δp(ω) = Σ_i q_i * u_i(ω)`

The induced dipole `Δp(ω)` highlights vibrational modes that couple strongly to the electric field (IR-active modes).

### Key Insight

This avoids O(N³) diagonalization. For each frequency:
- Factorize A(ω) once (Cholesky)
- Solve for dipole RHS
- Scan ω to sample spectrum

## Implementation Details

### Core Functions

- `dynamic_stiffness(K, M, omega, eta, stabilize)`: Constructs A(ω)
- `mechanical_greens_probing(K, M, omegas, eta, direction_vec, charges)`: Main probing loop
- `solve_response(K, M, omega, eta, charges, direction_vec)`: Single frequency solve
- `corner_quadrupole_charges(nx, ny, charge_val)`: Neutral charge pattern (+,+,-,-)

### Outputs

- `omega`: sampled frequencies
- `energy`: proxy amplitude per ω (||U||² average)
- `dipole`: complex dipole response per ω (3-vector)
- `n_probes`: total frequency points sampled

### Demo Usage

```bash
cd doc/py/OrderN_QM
python VibrationProbing.py \
  --nx 3 --ny 3 \
  --k1 5.0 --k2 5.0 --kdiag 5.0 \
  --fmin 0.1 --fmax 6.0 --nfreq 1000 \
  --eta 0.001
```

## Applicability to Silicon Nanocrystals

### Classical Vibration (Force Fields)

For silicon nanocrystals generated by `gen_nanocrystals.mjs`:

1. **Construct Dynamical Matrix**:
   - K = Hessian from MMFF force field (second derivatives of potential)
   - M = diagonal mass matrix (Si: 28.085 u, H: 1.008 u)
   - Use `mmff.getHessian3x3()` or similar for local blocks

2. **Apply VibrationProbing**:
   - Replace the triangular grid K with actual nanocrystal Hessian
   - Use atomic charges (e.g., partial charges from MMFF) for dipole coupling
   - Scan relevant frequency range (0-1000 cm⁻¹ for Si-Si modes)

### Connection to OrderN_QM Methods

Several techniques from `doc/py/OrderN_QM/` can enhance vibration spectroscopy:

1. **BlockCholesky.md**: Tiled Cholesky decomposition
   - Pre-factorize K and M once
   - Efficiently solve shifted systems for multiple ω
   - GPU-friendly for large nanocrystals (>1000 atoms)

2. **Chebyshev Acceleration** (LinearLagebra docs):
   - Accelerate iterative solvers when direct factorization is too expensive
   - Particularly useful for very large systems

3. **Green's Function Contour Integration** (GF.py):
   - Adapt electronic GF methods to mechanical GF
   - Use contour integration in complex ω-plane instead of real ω scan
   - Stochastic probing for massive systems

4. **Iterative Solvers** (NumricalStabilityLinearSolvers.md):
   - Jacobi/Gauss-Seidel with coloring for sparse K
   - Momentum acceleration for faster convergence

### Quantum Vibrational Corrections

For quantum mechanical accuracy (zero-point energy, anharmonicity):

- **Approximate_Overlap.md**: NDDO methods for electronic structure
- Could compute vibrational corrections via perturbation theory
- More advanced: path-integral MD or vibrational SCF

## Current Limitations

- VibrationProbing.py demo uses simple spring-mass truss (triangular grid)
- Needs adaptation for:
  - Actual MMFF force field Hessian
  - Correct atomic charges for dipole coupling
  - 3D geometry (currently 2D planar)
  - Boundary conditions (fixed atoms, free surfaces)

## Next Steps for Integration

1. Extract Hessian from MMFF for generated nanocrystals
2. Implement charge model (e.g., electronegativity equalization)
3. Adapt VibrationProbing for 3D systems with realistic masses
4. Benchmark against explicit diagonalization for small systems
5. Explore Block Cholesky for larger systems (>500 atoms)

---

# USER

OK, but we can construct dynamical matrix also purely on GPU using pyOpenCL althouhg maybe that will be slower than using CPU C++ because topology buding in in python is much slower than on GPU. Also it definitely help if we start from .mol file wich contain bonding topology than with .xyz where the bonds need to be found O(n^2). I'm now thingking what is the fastest way how to build hessian for Si Nanocrystal or carbon (diamond) cryeta. We generate themn using @gen_nanocrystals.mjs 
then we need to load them and assing the topology (bonds, angles) and evaluate forcefield. We can do it either using pyOpenCL or using C++ librarty binding

pyOpenCL is here:
@MMparams.py @MolecularDynamics.py @UFF.py @MMFF.py @UFF.cl  
@UFF.cl @relax_multi_mini.cl 

the C++ version (with python API) is here:
@MMFF_lib.cpp @MMFF.py @MolWorld_sp3.h @UFF.h @MMFFsp3_loc.h @MMFFBuilder.h @MMFFBuilderBase.h @UFFbuilder.h 

please examine the relevant files and make comprehensive description what is which module responsible for and how are all related, and how should we approach if we want to implement high-throughput workflow where nancrstals are generated (in large number for variations of shape, size, surface variations (defects etc.), then are loaded their Hessians are calculated by forcefield, and then we calculated their vibration spectra by response function (scanning energy and evaluating (EI-H)x=b ), eventually we want to do this in parallel for many cystals on GPU using pyOpenCL

Do not implment any code, make just detailed notes what we can reuse, how to integrate it, and what remains to be done

now I would like to add some consideration for the Hessian (Dynamical matrix). We can consider two distinct cases:
1) if we use only bonding forcefield (bonds, nangles, ihedras etc) but not no-bondong (pauli repulsion, vnad der waals, electrostatics) the hessian matrix should be very sparse (only nearest neibhors or next nearest neibors 1-2 and 1-3 interactions). We should explotit that for efficient evaluation. 
2) But the extranl exitation (the b-vectror, the impulse which caouse the response) is typically electric field acting of partical charges (e.g. negative ozygen atoms on the surface =O , -O-, -OH groups have strong FTIR response). To do it rigorously there should by also internal electrostatics. But one can think about splitting it. We aplly field to charges crerating some external force or displacement to compute vectro b (right-hand-side), but then internaly when we solve the matrix (EI-H)x=b we may use H comprising only bonding forcefield (ignoring interal interaction of charges), that would be sligtly inconsistent (ignoring second-order electrostatic effect while the molecule is vibrating) but it would be much faster. We should consider both ways.


---

# SWE-1.6

#/@ CODE MAP High-Throughput Hessian & Vibration Spectra Workflow for Si Nanocrystals
https://windsurf.com/codemaps/772fabe6-9c9e-4376-bc6e-b974c2ebb806-fe86ab10a43f3d18

# Comprehensive Module Analysis for Si Nanocrystal Hessian & Vibration Spectra

## 1. Module Responsibilities

### 1.1 Nanocrystal Generation

**[scripts/gen_nanocrystals.mjs](cci:7://file:///home/prokophapala/git/FireCore/scripts/gen_nanocrystals.mjs:0:0-0:0)** (JavaScript/Node.js)
- **Purpose**: Generate silicon nanocrystals with configurable geometry and surface chemistry
- **Key Functions**:
  - Builds Si supercell from CIF (Si-sym.cif)
  - Applies plane cuts using Miller index templates ({111}, {100}, etc.)
  - Recalculates bonds using MM parameters
  - Prunes undercoordinated Si atoms iteratively
  - Adds H caps to dangling bonds
  - Surface bridge operations (insertion/collapse of SiH2)
  - Outputs .mol2 files (with topology) and .xyz files
- **Output Format**: .mol2 files contain bonding topology, which is critical for Hessian computation

### 1.2 Parameter Loading (Common to Both pyOpenCL and C++)

**[pyBall/OCL/MMparams.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/OCL/MMparams.py:0:0-0:0)** (Python)
- **Purpose**: Load force field parameters from data files
- **Key Classes**:
  - [ElementType](cci:2://file:///home/prokophapala/git/FireCore/pyBall/OCL/MMparams.py:8:0-38:58): Element parameters (Z, valence, Rcov, RvdW, EvdW, Quff, UFF params, QEq params)
  - [AtomType](cci:2://file:///home/prokophapala/git/FireCore/pyBall/MMFF.py:96:0-124:5): Atom type parameters (parent, element, epair, valence, symmetry, Ruff, MMFF force constants)
  - [MMFFparams](cci:2://file:///home/prokophapala/git/FireCore/pyBall/OCL/MMparams.py:277:0-294:51): Container class with parameter dictionaries
- **Data Files**: ElementTypes.dat, AtomTypes.dat, BondTypes.dat, AngleTypes.dat
- **Reuse**: Can be reused directly by both pyOpenCL and C++ approaches

### 1.3 pyOpenCL Force Field Modules

**[pyBall/OCL/UFF.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/OCL/UFF.py:0:0-0:0)** (Python + OpenCL)
- **Purpose**: UFF force field evaluation on GPU
- **Key Class**: [UFF_CL](cci:2://file:///home/prokophapala/git/FireCore/pyBall/OCL/UFF.py:12:0-388:51)
  - Manages OpenCL buffers and kernels
  - Delegates topology building to `UFF_Builder`
  - Runs UFF kernels: evalBondsAndHNeigh_UFF, evalAngles_UFF, evalDihedrals_UFF, evalInversions_UFF
  - Supports multiple systems in parallel (nSystems parameter)
- **Buffers**: apos, fapos, fint, REQs, bonAtoms/Params, angAtoms/Params, dihAtoms/Params, invAtoms/Params, neighs
- **Topology**: Uses atom-to-force mapping (a2f_offsets, a2f_counts, a2f_indices) for efficient force assembly

**[pyBall/OCL/MolecularDynamics.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/OCL/MolecularDynamics.py:0:0-0:0)** (Python + OpenCL)
- **Purpose**: MD simulation using relax_multi_mini.cl kernel
- **Key Class**: [MolecularDynamics](cci:2://file:///home/prokophapala/git/FireCore/pyBall/OCL/MolecularDynamics.py:59:0-1968:21)
  - Supports multiple systems in parallel (nSystems)
  - Manages MMFF force evaluation on GPU
  - Includes non-bonded interactions (optional)
  - Surface interaction support
- **Buffers**: apos, aforce, REQs, neighs, apars, bLs, bKs, Ksp, Kpp, MDparams, constraints
- **Kernels**: getMMFFf4, updateAtomsMMFFf4, runMD, cleanForceMMFFf4

**[pyBall/OCL/MMFF.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/OCL/MMFF.py:0:0-0:0)** (Python - OCL version)
- **Purpose**: MMFF force field topology building for pyOpenCL
- **Key Functions**: Parameter loading, atom property initialization, REQ generation
- **Reuse**: Shares parameter loading with MMparams.py

### 1.4 OpenCL Kernels

**[cpp/common_resources/cl/UFF.cl](cci:7://file:///home/prokophapala/git/FireCore/cpp/common_resources/cl/UFF.cl:0:0-0:0)**
- **Purpose**: UFF force field kernels for GPU
- **Key Kernels**:
  - `clear_fapos_UFF`, `clear_fint_UFF`: Zero force buffers
  - `evalBondsAndHNeigh_UFF`: Bond evaluation + H-neigh vector computation
  - `evalAngles_UFF`: Angle evaluation
  - `evalDihedrals_UFF`: Dihedral evaluation
  - `evalInversions_UFF`: Inversion evaluation
  - `assembleForces_UFF`: Assemble force pieces from fint into fapos
  - `getNonBond`, `getNonBond_ex2`: Non-bonded interactions (LJ + Coulomb)
- **Features**: Debug controls, exclusion lists, PBC support
- **Hessian Relevance**: Currently computes forces only, but kernel structure could be extended for Hessian

**`cpp/common_resources/cl/relax_multi_mini.cl`**
- **Purpose**: MD relaxation and force evaluation kernels
- **Key Kernels**: getMMFFf4 (force evaluation), updateAtomsMMFFf4 (integration), runMD
- **Hessian Relevance**: Focuses on forces, not Hessian computation

### 1.5 C++ Force Field Implementation

**[cpp/common/molecular/UFF.h](cci:7://file:///home/prokophapala/git/FireCore/cpp/common/molecular/UFF.h:0:0-0:0)** (C++ header)
- **Purpose**: UFF force field implementation (inherits from NBFF)
- **Key Data Structures**:
  - `fint`: Force pieces array [ndihedrals*4 + ninversions*4 + nangles*3 + nbonds]
  - `hneigh`: Bond vectors [natoms*4]
  - Topology arrays: bonAtoms, bonParams, angAtoms, angParams, dihAtoms, dihParams, invAtoms, invParams
  - Neighbor indices: neighBs, angNgs, dihNgs, invNgs
  - `Buckets a2f`: Atom-to-force mapping for efficient assembly
- **Methods**: realloc, eval (bond/angle/dihedral/inversion), assembleForces
- **Hessian**: Not implemented, but structure supports extension

**[cpp/common/molecular/MMFFsp3_loc.h](cci:7://file:///home/prokophapala/git/FireCore/cpp/common/molecular/MMFFsp3_loc.h:0:0-0:0)** (C++ header)
- **Purpose**: MMFFsp3 localized force field (inherits from NBFF)
- **Key Data Structures**:
  - `fneigh`, `fneighpi`: Temporary force storage for neighbors
  - `apars`, `bLs`, `bKs`, `Ksp`, `Kpp`: Per-atom parameters
  - `angles`: Angle data [nnode*6]
- **Methods**: eval_atom (per-atom evaluation), assemble_atom
- **Hessian**: Not implemented, but per-atom structure is Hessian-friendly

**[cpp/common/molecular/MMFFBuilder.h](cci:7://file:///home/prokophapala/git/FireCore/cpp/common/molecular/MMFFBuilder.h:0:0-0:0)** (C++ header)
- **Purpose**: Topology builder for molecular systems
- **Key Class**: [Builder](cci:2://file:///home/prokophapala/git/FireCore/cpp/common/molecular/MMFFBuilder.h:92:0-3963:1) (inherits from BuilderBase)
- **Methods**: 
  - [toMMFFsp3_loc](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/MMFF.py:376:4-690:31): Convert AtomicSystem to MMFFsp3_loc representation
  - [makeConfGeom](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/MMFF.py:899:4-1026:108): Build configuration geometry from bond directions
  - [assignBondParams](cci:1://file:///home/prokophapala/git/FireCore/cpp/common/molecular/MMFFBuilder.h:525:4-550:5): Assign bond parameters from tables or UFF
- **Topology**: Handles bonds, angles, pi-orbitals, electron pairs, capping atoms
- **Reuse**: Can be called from Python via MMFF_lib.cpp bindings

**[cpp/common/molecular/MMFFBuilderBase.h](cci:7://file:///home/prokophapala/git/FireCore/cpp/common/molecular/MMFFBuilderBase.h:0:0-0:0)** (C++ header)
- **Purpose**: Base class for topology building
- **Key Structures**: [Atom](cci:2://file:///home/prokophapala/git/FireCore/cpp/common/molecular/MMFFBuilderBase.h:49:0-66:1), [AtomConf](cci:2://file:///home/prokophapala/git/FireCore/cpp/common/molecular/MMFFBuilderBase.h:80:0-191:1), [Bond](cci:2://file:///home/prokophapala/git/FireCore/cpp/common/molecular/MMFFBuilderBase.h:197:0-220:1), [Angle](cci:2://file:///home/prokophapala/git/FireCore/cpp/common/molecular/MMFFBuilderBase.h:226:0-246:1), [Dihedral](cci:2://file:///home/prokophapala/git/FireCore/cpp/common/molecular/MMFFBuilderBase.h:252:0-274:1), [Inversion](cci:2://file:///home/prokophapala/git/FireCore/cpp/common/molecular/MMFFBuilderBase.h:280:0-299:1)
- **Methods**: insertAtoms, insertBonds, addCaps, autoTypes
- **Reuse**: Foundation for both MMFF and UFF builders

**[cpp/common/molecular/UFFbuilder.h](cci:7://file:///home/prokophapala/git/FireCore/cpp/common/molecular/UFFbuilder.h:0:0-0:0)** (C++ header)
- **Purpose**: UFF-specific topology builder
- **Key Method**: [toUFF(UFF& ff, bool bRealloc=true)](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/UFF.py:50:4-67:23)
- **Topology**: Builds UFF-specific topology (bonds, angles, dihedrals, inversions)
- **Reuse**: Can be called from Python via MMFF_lib.cpp bindings

### 1.6 C++ Library and Python Bindings

**[cpp/libs/Molecular/MMFF_lib.cpp](cci:7://file:///home/prokophapala/git/FireCore/cpp/libs/Molecular/MMFF_lib.cpp:0:0-0:0)** (C++)
- **Purpose**: Python bindings for C++ force field library
- **Key Functions**:
  - [init()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/MMFF.py:1009:0-1035:14): Initialize world with parameters
  - [init_buffers()](cci:1://file:///home/prokophapala/git/FireCore/cpp/libs/Molecular/MMFF_lib.cpp:65:0-90:1), [init_buffers_UFF()](cci:1://file:///home/prokophapala/git/FireCore/cpp/libs/Molecular/MMFF_lib.cpp:23:0-63:1): Expose C++ arrays to Python
  - [print_debugs()](cci:1://file:///home/prokophapala/git/FireCore/cpp/libs/Molecular/MMFF_lib.cpp:93:0-106:1): Debug output
  - [makeGridFF()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/MMFF.py:308:0-318:15): Grid force field construction
- **Global Object**: `MolWorld_sp3 W` (contains MMFFsp3_loc or UFF instance)
- **Buffers Exposed**: apos, fapos, REQs, neighs, apars, bLs, bKs, Ksp, Kpp, etc.
- **Hessian**: No Hessian computation exposed yet

**[pyBall/MMFF.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/MMFF.py:0:0-0:0)** (Python bindings)
- **Purpose**: Python interface to C++ MMFF_lib
- **Key Functions**: 
  - [getAtomTypes()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/MMFF.py:131:0-136:17): Get atom type parameters from C++
  - [samplePBCindexes](cci:1://file:///home/prokophapala/git/FireCore/pyBall/MMFF.py:155:0-159:15), `projectBspline1D/2D/3D`: Utility functions
  - [fit_Bspline](cci:1://file:///home/prokophapala/git/FireCore/pyBall/MMFF.py:195:0-205:17), [fit2D_Bspline](cci:1://file:///home/prokophapala/git/FireCore/pyBall/MMFF.py:220:0-226:17), [fit3D_Bspline](cci:1://file:///home/prokophapala/git/FireCore/pyBall/MMFF.py:231:0-237:13): Fitting functions
- **Library Loading**: Loads MMFF_lib.so from cpp/Build
- **Hessian**: No Hessian functions exposed

## 2. Module Relationships

### 2.1 Data Flow

```
gen_nanocrystals.mjs (.mol2 output)
    ↓
Python loading (AtomicSystem)
    ↓
Parameter loading (MMparams.py)
    ↓
Topology building (MMFF.py OCL or C++ MMFFBuilder)
    ↓
Force field evaluation (UFF.py OCL or C++ UFF/MMFFsp3_loc)
    ↓
[MISSING] Hessian computation
    ↓
Vibration spectra (VibrationProbing.py)
```

### 2.2 Parallel Architectures

**pyOpenCL Path**:
- Python script loads .mol2 → builds topology in Python → uploads to GPU → OpenCL kernels evaluate forces
- **Pros**: Pure Python workflow, easier to modify, GPU parallelism built-in
- **Cons**: Topology building in Python is slower than C++, limited to GPU memory

**C++ Path**:
- Python script loads .mol2 → calls C++ via ctypes → C++ builds topology → C++ evaluates forces → results back to Python
- **Pros**: Fast topology building in C++, can use OpenMP for CPU parallelism
- **Cons**: Requires recompilation for changes, Python-C++ interface overhead

### 2.3 Force Field Options

**UFF** (Universal Force Field):
- Generic, parameter-based
- Supports bonds, angles, dihedrals, inversions
- Good for diverse systems (Si, C, organic molecules)
- Implemented in: UFF.h (C++), UFF.py (pyOpenCL), UFF.cl (OpenCL kernel)

**MMFFsp3_loc** (Localized MMFF):
- Specialized for sp3 systems (diamond, silicon)
- Pi-orbital explicit representation
- Simulates dihedrals via pi-pi interactions
- Implemented in: MMFFsp3_loc.h (C++), MolecularDynamics.py (pyOpenCL), relax_multi_mini.cl (OpenCL kernel)

## 3. High-Throughput Workflow Design

### 3.1 Proposed Architecture

```
Stage 1: Generation (CPU, Sequential)
├── gen_nanocrystals.mjs
│   ├── Generate N nanocrystals with variations
│   └── Output: N .mol2 files (with topology)

Stage 2: Preprocessing (CPU, Parallelizable)
├── Load .mol2 files
├── Assign atom types (from ElementTypes/AtomTypes)
├── Generate REQ parameters
└── Output: N AtomicSystem objects with topology

Stage 3: Hessian Computation (GPU or CPU, Parallel)
├── Option A: C++ + OpenMP (CPU)
│   ├── MMFFBuilder::toMMFFsp3_loc (fast topology)
│   ├── MMFFsp3_loc::eval (force evaluation)
│   ├── [NEW] Hessian computation (finite difference or analytical)
│   └── Output: N Hessian matrices (3N×3N sparse)
│
├── Option B: pyOpenCL (GPU)
│   ├── MMFF.py OCL (topology in Python - slower)
│   ├── UFF.py or MolecularDynamics.py (GPU force eval)
│   ├── [NEW] Hessian computation on GPU
│   └── Output: N Hessian matrices
│
└── Option C: Hybrid
    ├── C++ for topology building (fast)
    ├── pyOpenCL for force evaluation (GPU parallel)
    ├── [NEW] Hessian computation on GPU
    └── Output: N Hessian matrices

Stage 4: Vibration Spectra (CPU/GPU, Parallel)
├── VibrationProbing.py
│   ├── Dynamic stiffness: A(ω) = K - (ω + iη)²M
│   ├── Dipole-driven probing: solve A(ω)U = f
│   ├── Scan ω to sample spectrum
│   └── Output: N vibration spectra
│
└── [OPTIMIZATION] Block Cholesky / iterative solvers
    ├── Pre-factorize K once
    ├── Solve shifted systems for multiple ω
    └── GPU-friendly for large systems
```

### 3.2 Performance Considerations

**Topology Building**:
- C++ MMFFBuilder: O(N) with low constant factor (compiled code)
- Python MMFF.py OCL: O(N) with higher constant factor (interpreted)
- **Recommendation**: Use C++ for topology building if possible

**Force Evaluation**:
- pyOpenCL: Excellent for many systems in parallel (batch processing)
- C++ OpenMP: Good for single large system
- **Recommendation**: Use pyOpenCL for batch processing of many nanocrystals

**Hessian Computation**:
- Finite difference: 6N force evaluations per Hessian (expensive)
- Analytical: Requires new kernel implementation
- **Recommendation**: Start with finite difference on GPU, optimize later

**Vibration Spectra**:
- Direct diagonalization: O(N³), not scalable
- Green's function probing: O(N×Nfreq) with linear solves
- **Recommendation**: Use Green's function probing with Block Cholesky

### 3.3 Hessian Sparsity and Electrostatic Considerations

**Case 1: Bonding-Only Forcefield (Sparse Hessian)**
- If we use only bonding terms (bonds, angles, dihedrals, inversions) without non-bonded interactions (Pauli repulsion, van der Waals, electrostatics), the Hessian matrix becomes very sparse
- **Sparsity Pattern**: Only nearest neighbors (1-2 interactions) and next-nearest neighbors (1-3 interactions) contribute
- **Bond terms**: Connect directly bonded atoms (sparse, O(N) non-zeros)
- **Angle terms**: Connect 1-3 atoms (still sparse, O(N) non-zeros)
- **Dihedral/Inversion terms**: Connect 1-4 atoms (sparse, O(N) non-zeros)
- **Exploiting Sparsity**:
  - Use sparse matrix formats (CSR, COO, or specialized block-sparse for 3×3 atomic blocks)
  - Sparse linear solvers (e.g., sparse Cholesky, sparse LU)
  - For Green's function probing: sparse matrix-vector products are very efficient
  - GPU sparse linear algebra (cuSPARSE, pyOpenCL sparse kernels)
- **Advantages**: Much faster computation, lower memory usage, better scalability
- **Limitations**: Less accurate for systems where non-bonded interactions are important (e.g., hydrogen bonding, long-range electrostatics)

**Case 2: Electrostatic Excitation (Rigorous vs Approximate)**

**Rigorous Approach (Full Electrostatics)**:
- Include non-bonded electrostatics in both the Hessian and the excitation
- **Hessian**: Include charge-charge interaction terms (∂²E/∂x∂x from Coulomb interactions)
- **Excitation**: Apply electric field to partial charges to generate b-vector (right-hand-side)
- **Consistency**: Second-order electrostatic effects are included during vibration
- **Implementation**:
  - Need to compute partial charges (QEq, Mulliken, or from MMFF parameters)
  - Include Coulomb interaction in force field evaluation
  - Hessian becomes less sparse (long-range 1/r interactions)
  - Use Ewald summation or PME for periodic systems
- **Pros**: Physically accurate FTIR/Raman spectra
- **Cons**: Computationally expensive, dense Hessian (or long-range sparse), harder to scale

**Approximate Approach (Split Electrostatics)**:
- **Excitation**: Apply electric field to partial charges to generate b-vector (right-hand-side)
  - b = q_i * E (force on atom i from external field)
  - For FTIR: E is the electric field of incident light
  - For Raman: similar but with polarizability derivatives
- **Hessian**: Use only bonding forcefield (ignore internal electrostatic interactions in H)
  - Sparse Hessian (Case 1)
  - Second-order electrostatic effects during vibration are ignored
- **Consistency**: Slightly inconsistent - electrostatic forces drive vibration but electrostatic stiffness is ignored
- **Approximation Validity**:
  - Good for qualitative spectra
  - Reasonable for systems where bonding dominates vibrational modes
  - Less accurate for modes strongly affected by charge redistribution
- **Advantages**: Much faster, sparse Hessian, scalable to large systems
- **Limitations**: Misses frequency shifts from electrostatic coupling, may mispredict intensities for strongly ionic modes

**Recommendation for Si Nanocrystals**:
- **Start with Approximate Approach**: Bonding-only Hessian + electrostatic excitation
  - Si-Si bonds dominate vibrational modes
  - Surface groups (Si-H, Si-O) can be modeled with partial charges for excitation
  - Fast enough for high-throughput screening
- **Validate with Rigorous Approach**: Compare subset of systems with full electrostatics
  - Check if approximation introduces significant errors
  - Refine if needed for publication-quality results
- **Hybrid Strategy**: Use bonding-only Hessian for screening, full electrostatics for final refinement

**Implementation Notes**:
- **Partial Charges**: Use MMFF QEq parameters or compute from electronegativity equalization
- **Surface Groups**: Assign charges to surface terminations (H: +δ, O: -δ, OH: partial)
- **FTIR Intensity**: Proportional to dipole derivative, which depends on charges and their motion
- **Raman Intensity**: Depends on polarizability derivatives (more complex, may need empirical models)

## 8. Simplified Bonding-Only Approach for Fast Vibration Spectroscopy

### 8.1 Goal: Super Fast Linear Response with Minimal Forcefield

For high-throughput screening of silicon nanocrystals, we want the fastest possible workflow:
- **Forcefield**: Only bonds and angles (1-2 and 1-3 interactions)
- **No non-bonded**: Ignore Pauli repulsion, van der Waals, electrostatics in Hessian
- **Sparse Hessian**: Exploit sparsity for efficient linear algebra
- **Linear response**: Use Green's function probing instead of full diagonalization

### 8.2 Existing Hessian Functions

**Python Interface** (`pyBall/MMFF.py`):
```python
# 3x3 Hessian blocks for selected atoms (diagonal blocks only)
def getHessian3x3(inds, dx=1e-4, bDiag=True):
    # Returns: (n, 4, 3) array
    # For each atom: 3x3 eigenvector matrix (first 3 rows) + eigenvalues (4th row)
    # Uses central difference finite difference on C++ side
    # Cost: 3 force evaluations per atom (x, y, z displacements)

# Full 3N×3N Hessian for selected atoms
def getHessian3Nx3N(inds, dx=1e-4):
    # Returns: (3N, 3N) dense matrix
    # Cost: 3N force evaluations (displace each atom in x, y, z)
    # Uses central difference: ∂F_j/∂x_i = (F_j(x+δ) - F_j(x-δ)) / (2δ)
```

**C++ Implementation** (`cpp/libs/Molecular/MMFF_lib.cpp`):
```cpp
void getHessian3x3(int n, int* inds, double* Hess_, double dx, bool bDiag):
    // For each selected atom:
    // 1. Save original position
    // 2. Displace in x: evaluate forces, restore
    // 3. Displace in y: evaluate forces, restore
    // 4. Displace in z: evaluate forces, restore
    // 5. Compute Hessian block via central difference
    // 6. Optionally diagonalize to get eigenvectors/eigenvalues
    // Uses W.eval_no_omp() for force evaluation

void getHessian3Nx3N(int n, int* inds, double* out_hessian, double dx):
    // For each selected atom p and coordinate k:
    // 1. Displace atom p in direction k by +dx
    // 2. Evaluate forces on all selected atoms
    // 3. Displace atom p in direction k by -dx
    // 4. Evaluate forces on all selected atoms
    // 5. Compute ∂F_o/∂x_p via central difference
    // Total: 3N force evaluations
```

### 8.3 Simplified Linear Response from VibrationProbing.py

**Core Algorithm** (`doc/py/OrderN_QM/VibrationProbing.py`):
```python
def mechanical_greens_probing(K, M, omegas, eta=1e-3, direction_vec=None, charges=None):
    """
    Dipole-driven probing of mechanical Green's function.
    
    Key steps:
    1. For each frequency ω in omegas:
       a. Build dynamic stiffness: A(ω) = K - (ω + iη)² M
       b. Build dipole RHS: f_i = q_i * direction_vec (uniform field)
       c. Solve A(ω) U = f (linear system)
       d. Compute induced dipole: Δp = Σ_i q_i U_i
       e. Store |Δp(ω)| as spectral response
    
    2. Peaks in |Δp(ω)| correspond to IR-active vibrational modes
    
    Advantages:
    - O(N × Nfreq) vs O(N³) for full diagonalization
    - Only need to solve linear systems, not compute all eigenvectors
    - Can reuse Cholesky factorization for multiple frequencies
    """
```

**Adaptation for Molecular Systems**:
- Replace simple spring K with actual Hessian from forcefield
- Use actual atomic masses (Si: 28.085 u, H: 1.008 u, C: 12.011 u)
- Use partial charges for dipole coupling (from MMFF parameters)
- Handle 3D geometry (current demo is 2D)

### 8.4 Topology Building Workflow

**MolWorld_sp3::init()** (`cpp/common/molecular/MolWorld_sp3.h`):
```cpp
virtual int init(){
    // 1. Load parameters if not already loaded
    if( params.atypes.size() == 0 ){
        initParams( "common_resources/ElementTypes.dat", 
                   "common_resources/AtomTypes.dat", 
                   "common_resources/BondTypes.dat", 
                   "common_resources/AngleTypes.dat", 
                   "common_resources/DihedralTypes.dat" );
    }
    
    // 2. Load molecule from XYZ file
    if( xyz_name && (bMMFF || bBUFF) ){
        buildMolecule_xyz( xyz_name );  // Uses Builder to build topology
    }
    
    // 3. Build force fields
    if(bMMFF || bBUFF){
        makeFFs();  // Calls makeMMFFs() internally
    }
}
```

**MolWorld_sp3::makeFFs()** → **makeMMFFs()**:
```cpp
void makeMMFFs(){
    // 1. Validate topology (check bonds in neighbors)
    builder.checkBondsInNeighs(true);
    
    // 2. Reorder atoms (non-capping atoms first)
    builder.numberAtoms();
    builder.sortConfAtomsFirst();
    
    // 3. Build UFF or MMFF topology
    if( bBUFF ){
        MM::UFFBuilder uff_builder;
        uff_builder.cloneFrom(builder);
        uff_builder.assignUFFtypes(...);
        uff_builder.assignUFFparams(...);
        uff_builder.toUFF( ffu, true );  // Convert to UFF representation
    }else{
        builder.assignTypes();  // Assign MMFF atom types
        builder.toMMFFsp3_loc( ffl, true, bEpairs, bBUFF );  // Convert to MMFFsp3_loc
        builder.toMMFFsp3( ff, true, bEpairs );
    }
    
    // 4. Setup PBC if needed
    if(bPBC){
        // Set lattice vectors and PBC shifts
        ffu.setLvec( builder.lvec );
        npbc = makePBCshifts( nPBC, builder.lvec );
        ffu.bindShifts(npbc, pbc_shifts);
        ffu.makeNeighCells( npbc, pbc_shifts );
    }
}
```

**Python Interface** (`pyBall/MMFF.py`):
```python
def init( xyz_name, bMMFF=True, bEpairs=False, bBUFF=False, ... ):
    # 1. Set world parameters
    W.xyz_name = xyz_name
    W.bMMFF = bMMFF
    W.bBUFF = bBUFF
    
    # 2. Initialize parameters
    W.initParams( sElementTypes, sAtomTypes, sBondTypes, sAngleTypes, sDihedralTypes )
    
    # 3. Initialize world (calls C++ MolWorld_sp3::init())
    ret = W.init()
    
    # 4. Initialize buffers for Python access
    init_buffers()  # Expose C++ arrays to Python
```

### 8.5 Recommended Workflow for Bonding-Only Hessian

**Step 1: Load molecule and build topology (C++ path)**
```python
from pyBall import MMFF as mmff

# Initialize with bonding-only MMFF (no pi-orbitals, no dihedrals)
mmff.init(
    xyz_name="data/xyz/CH4.xyz",  # or H2.xyz for smallest test
    bMMFF=True,      # Use MMFF (not UFF)
    bEpairs=False,   # No electron pairs (simplifies topology)
    bBUFF=False,     # Not using UFF
)

# Turn off non-bonded interactions
mmff.setSwitches(
    NonBonded=-1,    # Disable non-bonded (LJ, Coulomb)
    MMFF=+1,         # Enable MMFF bonding
    SurfAtoms=-1,    # Disable surface atoms
    GridFF=-1        # Disable grid force field
)
```

**Step 2: Compute Hessian (finite difference)**
```python
# Get all atom indices
natoms = mmff.natoms
all_inds = np.arange(natoms, dtype=np.int32)

# Compute full 3N×3N Hessian
H = mmff.getHessian3Nx3N(all_inds, dx=1e-4)

# H is now a (3N, 3N) symmetric matrix
# For bonding-only forcefield, H should be sparse (1-2 and 1-3 interactions only)
```

**Step 3: Build mass matrix**
```python
# Get atomic masses from atom types
atom_types = mmff.getAtomTypes()
masses = np.array([atom_types[atype].mass for atype in mmff.atypes])

# Build diagonal mass matrix (3N×3N)
M_diag = np.repeat(masses, 3)  # Each atom has 3 DOFs
M = np.diag(M_diag)
```

**Step 4: Linear response vibration spectroscopy**
```python
from doc.py.OrderN_QM.VibrationProbing import mechanical_greens_probing

# Define frequency range for scanning
omegas = np.linspace(0.01, 0.5, 1000)  # Adjust range for your system

# Define charges for dipole coupling (simplified: use partial charges)
# For Si nanocrystals: Si (0), H (+δ), O (-δ), OH (partial)
charges = np.zeros(natoms)
# Assign charges based on atom types or element
for i, atype in enumerate(mmff.atypes):
    element = atom_types[atype].element
    if element == 1:   # H
        charges[i] = +0.1  # Simplified
    elif element == 8: # O
        charges[i] = -0.2  # Simplified

# Compute spectrum via Green's function probing
result = mechanical_greens_probing(
    K=H,              # Hessian from forcefield
    M=M,              # Mass matrix
    omegas=omegas,    # Frequency range
    eta=1e-3,         # Damping parameter
    direction_vec=np.array([1.0, 0.0, 0.0]),  # Electric field direction
    charges=charges,   # Partial charges
    dim=3              # 3D system
)

# Extract dipole response
dipole_response = result["dipole"]  # Complex (Nfreq, 3)
spectrum = np.linalg.norm(dipole_response, axis=1)  # Magnitude per frequency

# Plot spectrum
import matplotlib.pyplot as plt
plt.plot(omegas, spectrum)
plt.xlabel("Frequency (ω)")
plt.ylabel("|Dipole Response|")
plt.title("FTIR Spectrum (Bonding-Only)")
plt.show()
```

### 8.6 Debugging with Small Test Systems

**Smallest Test: H2 molecule** (`tests/pyFireball/relaxed_mols/H2.xyz`):
```
2

H -0.356688 0.000000 0.000000 0.000000
H 0.356688 0.000000 0.000000 0.000000
```
- 2 atoms, 1 bond
- Hessian: 6×6 matrix
- Expected: 1 vibrational mode (stretch), 5 zero modes (3 translations + 2 rotations)

**Small Diamond-like: CH4 molecule** (`tests/pyFireball/relaxed_mols/CH4.xyz`):
```
5

C 0.000000 0.000000 0.000000 0.000000 1.908000 
H 0.630886 0.000000 -0.892211 0.000000 1.487000 
H 0.630886 0.000000 0.892210 0.000000 1.487000 
H -0.630888 -0.892209 -0.000002 0.000000 1.487000 
H -0.630888 0.892209 -0.000004 0.000000 1.487000 
```
- 5 atoms, 4 bonds, 6 angles
- Hessian: 15×15 matrix
- Expected: 9 vibrational modes (3N-6 for non-linear molecule)

**Debugging Workflow**:
```python
# Test with H2 first
mmff.init(xyz_name="tests/pyFireball/relaxed_mols/H2.xyz", bMMFF=True)
mmff.setSwitches(NonBonded=-1, MMFF=+1)

# Compute Hessian
H = mmff.getHessian3Nx3N(np.arange(2), dx=1e-4)
print("Hessian shape:", H.shape)
print("Hessian:\n", H)

# Check sparsity
print("Non-zero elements:", np.count_nonzero(H))
print("Expected non-zeros for H2: ~12 (bond stiffness)")

# Test linear response
# ... (as above)

# Then test with CH4
mmff.init(xyz_name="tests/pyFireball/relaxed_mols/CH4.xyz", bMMFF=True)
# ... (same workflow)
```

### 8.7 Performance Optimizations for Bonding-Only Hessian

**Exploit Sparsity**:
```python
from scipy.sparse import csr_matrix

# Convert dense Hessian to sparse format
# For bonding-only, most elements are zero
H_sparse = csr_matrix(H)

# Use sparse linear algebra
from scipy.sparse.linalg import spsolve
# Note: mechanical_greens_probing currently uses dense numpy
# Would need to adapt for sparse matrices
```

**Analytical Hessian (Future Work)**:
- Extend C++ UFF.h or MMFFsp3_loc.h with analytical Hessian kernels
- Bond Hessian: ∂²E_bond/∂x∂y (analytic formula)
- Angle Hessian: ∂²E_angle/∂x∂y (analytic formula)
- Cost: O(N) vs O(N) for finite difference, but with much smaller constant
- Would require new OpenCL kernels for GPU implementation

**Block Cholesky for Multiple Frequencies**:
```python
# Pre-factorize K (Hessian) once
L = np.linalg.cholesky(K + 1e-6 * np.eye(K.shape[0]))  # Add regularization

# For each frequency, solve shifted system efficiently
for omega in omegas:
    z = omega + 1j * eta
    A_shifted = K - (z * z) * M + 1e-6 * np.eye(K.shape[0])
    # Can use Woodbury formula or update Cholesky factor
    # More efficient than full factorization per frequency
```

### 8.8 Limitations of Bonding-Only Approach

**Missing Physics**:
- No electrostatic coupling between charged groups
- No van der Waals interactions (important for surface-surface interactions)
- No hydrogen bonding (important for OH-terminated surfaces)
- No long-range effects (important for large nanocrystals)

**Accuracy Trade-offs**:
- Frequency shifts from electrostatic coupling are ignored
- Intensity predictions may be less accurate for ionic modes
- Surface modes may be less accurate
- Not suitable for publication-quality spectra without validation

**When to Use**:
- High-throughput screening of many nanocrystals
- Qualitative trends (size dependence, shape dependence)
- Initial exploration before full QM calculations
- Debugging and development of workflow

**When to Avoid**:
- Publication-quality FTIR/Raman spectra
- Systems with strong electrostatic effects (ionic compounds)
- Surface chemistry dominated by hydrogen bonding
- Validation against experimental data

## 4. Integration Strategy

### 4.1 What Can Be Reused

**From gen_nanocrystals.mjs**:
- .mol2 output format (contains topology)
- Nanocrystal generation logic (can be called from Python if needed)

**From pyOpenCL modules**:
- MMparams.py: Parameter loading (reuse directly)
- UFF.py: GPU force evaluation (reuse directly)
- MolecularDynamics.py: MMFF force evaluation on GPU (reuse directly)
- OpenCL kernels: Force evaluation (extend for Hessian)

**From C++ modules**:
- MMFFBuilder.h: Fast topology building (call via Python bindings)
- UFF.h / MMFFsp3_loc.h: Force evaluation (call via Python bindings)
- MMFF_lib.cpp: Python bindings (extend for Hessian)

**From VibrationProbing.py**:
- Green's function probing algorithm (reuse directly)
- Dynamic stiffness construction (adapt for 3D)
- Dipole-driven probing (adapt for Si nanocrystals)

### 4.2 Integration Points

**Option 1: Pure pyOpenCL Path**
```python
# Workflow
for mol2_file in nanocrystal_files:
    # Load .mol2 (topology already present)
    mol = load_mol2(mol2_file)
    
    # Build topology in Python (slower but pure Python)
    mmff = MMFF()
    mmff.toMMFFsp3_loc(mol, atom_types)
    
    # Upload to GPU and evaluate forces
    md = MolecularDynamics()
    md.realloc(mmff, nSystems=1)
    md.pack_system(mmff, iSys=0)
    md.run_eval_step()
    
    # Compute Hessian (finite difference)
    H = compute_hessian_finite_diff(md, mol)
    
    # Compute vibration spectra
    spectrum = vibration_probing(H, mol.masses)
```

**Option 2: C++ + pyOpenCL Hybrid**
```python
# Workflow
for mol2_file in nanocrystal_files:
    # Load .mol2
    mol = load_mol2(mol2_file)
    
    # Build topology in C++ (fast)
    W.init()  # Initialize C++ world
    W.load_mol(mol2_file)  # C++ handles topology
    
    # Get force evaluation function from C++
    # Or extract parameters and use pyOpenCL
    
    # Compute Hessian using C++ or pyOpenCL
    H = compute_hessian_hybrid()
    
    # Compute vibration spectra
    spectrum = vibration_probing(H, mol.masses)
```

**Option 3: Batch GPU Processing**
```python
# Process many nanocrystals in parallel on GPU
all_mols = [load_mol2(f) for f in nanocrystal_files[:batch_size]]
all_mmff = [build_topology(mol) for mol in all_mols]

# Upload all to GPU at once
md = MolecularDynamics()
md.realloc(all_mmff[0], nSystems=batch_size)
for i, mmff in enumerate(all_mmff):
    md.pack_system(mmff, iSys=i)

# Evaluate all in parallel
md.run_eval_step()

# Compute Hessians (one at a time or in batches)
for i in range(batch_size):
    H[i] = compute_hessian_finite_diff(md, iSys=i)
    spectrum[i] = vibration_probing(H[i], all_mols[i].masses)
```

## 5. Remaining Work

### 5.1 Hessian Computation

**Status**: Not implemented in current codebase

**Required Implementation**:
1. **Finite Difference Hessian** (easiest, immediate):
   - For each atom i and coordinate α:
     - Displace atom i by δ in direction α
     - Evaluate forces on all atoms
     - Compute ∂F/∂x via central difference
   - Cost: 6N force evaluations per Hessian
   - Can be parallelized across atoms on GPU

2. **Analytical Hessian** (more efficient, more work):
   - Extend UFF.cl or relax_multi_mini.cl kernels
   - Add Hessian kernels for each term (bond, angle, dihedral, inversion)
   - Accumulate Hessian contributions
   - Cost: O(N) with good constant factor

**Recommendation**: Start with finite difference on GPU, optimize with analytical kernels later

### 5.2 Vibration Spectra Adaptation

**Status**: VibrationProbing.py exists for 2D truss, needs adaptation

**Required Adaptation**:
1. **3D Geometry**: Current demo is 2D triangular grid
   - Adapt for 3D nanocrystals
   - Use actual atomic masses (Si: 28.085 u, H: 1.008 u)

2. **Hessian Integration**: 
   - Replace simple spring K with actual Hessian from force field
   - Use sparse matrix representation for efficiency

3. **Dipole Coupling**:
   - Compute atomic charges (QEq from MMFF parameters)
   - Implement dipole-driven probing for IR spectra

4. **Boundary Conditions**:
   - Handle fixed atoms vs free surfaces
   - Consider periodic boundary conditions if needed

### 5.3 Linear Solver Optimization

**Status**: VibrationProbing.py uses direct Cholesky, not optimized for large systems

**Required Optimization**:
1. **Block Cholesky** (from BlockCholesky.md):
   - Implement tiled Cholesky decomposition
   - Pre-factorize K and M once
   - Efficiently solve shifted systems for multiple ω
   - GPU-friendly for large nanocrystals

2. **Iterative Solvers** (from NumricalStabilityLinearSolvers.md):
   - Jacobi/Gauss-Seidel with coloring for sparse K
   - Chebyshev acceleration for faster convergence
   - Suitable for very large systems (>1000 atoms)

3. **Green's Function Contour Integration** (from GF.py):
   - Adapt electronic GF methods to mechanical GF
   - Use complex ω-plane contour integration
   - Stochastic probing for massive systems

### 5.4 Data Management

**Status**: No infrastructure for batch processing

**Required Infrastructure**:
1. **Batch Loading**: Load many .mol2 files efficiently
2. **Parameter Caching**: Cache atom types and force field parameters
3. **Result Storage**: Save Hessians and spectra in efficient format (HDF5, npy)
4. **Parallel Orchestration**: Manage GPU memory for batch processing

### 5.5 Python-C++ Interface

**Status**: MMFF_lib.cpp exposes basic functions, not Hessian

**Required Extensions**:
1. **Hessian Computation Function**: Add C++ function for Hessian computation
2. **Topology Extraction**: Extract Hessian from C++ to Python efficiently
3. **Batch Processing**: Support multiple systems in C++ interface

## 6. Recommended Implementation Path

### Phase 1: Proof of Concept (1-2 weeks)
- Use existing gen_nanocrystals.mjs to generate small test set
- Implement finite difference Hessian in Python (using pyOpenCL force evaluation)
- Adapt VibrationProbing.py for 3D nanocrystals
- Test on small systems (<100 atoms)

### Phase 2: Optimization (2-3 weeks)
- Implement C++ topology building (call MMFFBuilder from Python)
- Implement batch GPU processing for force evaluations
- Add Block Cholesky for linear solves
- Test on medium systems (100-500 atoms)

### Phase 3: Scaling (3-4 weeks)
- Implement analytical Hessian kernels in OpenCL
- Add iterative solvers for very large systems
- Implement Green's function contour integration
- Test on large systems (500-2000 atoms)

### Phase 4: Production (2-3 weeks)
- Optimize batch processing pipeline
- Add result storage and analysis tools
- Benchmark and optimize performance
- Document workflow

## 7. Key Decision Points

**Topology Building**: C++ (fast) vs Python (slow but easier to modify)
- **Recommendation**: Start with Python for flexibility, move to C++ if bottleneck

**Force Evaluation**: pyOpenCL (GPU batch) vs C++ OpenMP (CPU single)
- **Recommendation**: Use pyOpenCL for batch processing of many nanocrystals

**Hessian Computation**: Finite difference (easy) vs Analytical (efficient)
- **Recommendation**: Start with finite difference, implement analytical kernels later

**Linear Solver**: Direct (Cholesky) vs Iterative (Jacobi/GS)
- **Recommendation**: Use Block Cholesky for medium systems, iterative for large systems

**Vibration Method**: Full diagonalization vs Green's function probing
- **Recommendation**: Green's function probing (scales better)

---


