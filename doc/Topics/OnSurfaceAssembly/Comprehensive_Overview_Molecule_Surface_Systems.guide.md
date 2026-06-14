# Comprehensive Overview: Molecule-on-Surface Systems in FireCore

**Purpose**: This document consolidates all implemented modules, scripts, and visualization programs for generating rigid-body configurations of molecules on surfaces, scanning their interaction energies (rigid and relaxed), and evaluating with different force fields (GridFF, FoldedAtomic Functions, Electrostatic Embedding). It references code maps, key files, and reports to provide a complete picture of what works, what doesn't, and where to find each component.


## 1. High-Level Architecture Overview

The ecosystem consists of several interconnected subsystems:

### 1.1 Core Subsystems
- **GridFF**: Precomputed B-spline substrate potentials (Pauli/London/Coulomb) for fast GPU evaluation.
- **InteractionScanner (XYZ)**: On-the-fly LJ/Coulomb/HBond evaluation with substrate atoms (CPU/GPU).
- **FoldedAtomicFunctions**: Separable basis fitting for surface potentials; experimental.
- **Rigid-Body Assembly**: GPU-accelerated pose generation and clash scoring for multi-molecule arrangements.
- **Surface Sampling**: Height-map extraction (iso-surfaces) via z-search on (x,y) grids.
- **MMFF Force Field**: Molecular mechanics with linearized topology for surface relaxation.
- **GUI Placement Tools**: Interactive VisPy interfaces for manual pose design.
- **Path Optimization**: Batch replica optimization with causal tethers.

### 1.2 Key Directories
- `pyBall/`: Python wrappers, OpenCL interfaces, GUIs, and utilities.
- `cpp/common_resources/cl/`: OpenCL kernels (relax_multi.cl, Assembly.cl, Rigid.cl).
- `tests/tMMFF/`: Test scripts, tutorials, and batch workflows.
- `doc/Topics/OnSurfaceAssembly/`: Reports and discussions.

## 2. Subsystem Deep Dives

### 2.1 GridFF (Precomputed Substrate Potentials)

#### 2.1.1 What It Does
- Generates 3D B-spline grids (Pauli/London/Coulomb) from a substrate XYZ.
- Uploads to GPU as 3D texture; enables fast tricubic interpolation.
- Used for rigid scans, batch evaluation, and surface sampling.

#### 2.1.2 Key Files and Flow
- **Generation**: `pyBall/tests/ocl_GridFF_new.py:test_gridFF_ocl()`
  - Loads substrate XYZ; defines grid with `g0 = (-Lx/2, -Ly/2, z_top)` [Trace 1a].
  - Calls `make_MorseFF()` for Pauli/London grids via OpenCL kernel [Trace 1b].
  - Calls `makeCoulombEwald_slab()` for Coulomb grid via Ewald summation [Trace 1c].
  - Saves `Bspline_PLQd.npy` (nz,ny,nx,3) to disk [Trace 1d].
- **Loading & GPU Setup**: `tests/tMMFF/run_test_GridFF_CaF2.py`
  - Loads `.npy` with `np.load()` [Trace 2a].
  - Calls `gpu.initGridFF()` → `MolecularDynamics.initGridFF()` [Trace 2b].
  - Sets grid metadata (`grid_ns`, `grid_invStep`, `grid_p0`, `GFFParams`) [Trace 2c].
  - Uploads B-spline data to GPU buffer/texture [Trace 2d].
  - Generates kernel arguments for `getNonBond_GridFF_Bspline` [Trace 2e].
- **Evaluation**: `MolecularDynamics.run_getNonBond_GridFF_Bspline()`
  - Launches kernel `getNonBond_GridFF_Bspline` (or texture variant) [Trace 3a].
  - Per-atom kernel transforms world→grid coords, evaluates B-spline via `fe3d_pbc_comb()` [Trace 3d], scales forces [Trace 3f], writes result [Trace 3g].

#### 2.1.3 Status
- **Works**: Generation, loading, GPU force evaluation for single molecules and batches.
- **Known Issues**: Coulomb scaling discrepancy (possible double 14.39 factor) vs. InteractionScanner; parity tests pending (see SurfacePotentialVisualization_discussion.md).
- **References**: `tests/tMMFF/GridFF_CaF2_doc_tutorial.md`, `cpp/common_resources/CaF2_6L_Ni3_rect_nx2_nz1_L2_top/CaF2_report.md`.

### 2.2 InteractionScanner (XYZ Rigid Scanning)

#### 2.2.1 What It Does
- Rigid scans: molecule-substrate interaction curves along z (1D) or XY grids (2D).
- On-the-fly evaluation: LJ/Coulomb/HBond with PBC; optional macro embedding.
- Backends: `reference` (CPU), `fast_gpu` (GPU kernel), `folded_gpu` (experimental, not wired).

#### 2.2.2 Key Files and Flow
- **1D Rigid Scan**: `tests/tMMFF/test_interaction_scan.py:scan_1D()`
  - Calls `InteractionEnergy.scan_1D()` [Trace 4a].
  - Loops over z steps; updates probe position [Trace 4b]; uploads to GPU [Trace 4c].
  - Calls `run_getSurfMorse()` → kernel `getSurfMorse` [Trace 4d].
  - Downloads forces; extracts energy from `force.w` [Trace 4e].
- **Batch Rigid Evaluation**: `pyBall/OCL/MolecularDynamics.eval_rigid_getSurfMorse_components()`
  - Generates transform grid (rotations + translations) [Trace 5a].
  - Wraps transforms with PBC; uploads chunks to GPU [Trace 5c].
  - Launches `getSurfMorse` for each chunk [Trace 5d]; downloads forces [Trace 5e].
  - Returns per-component energies (LJ, Coulomb, HBond).

#### 2.2.3 Status
- **Works**: 1D/2D rigid scans (reference/fast_gpu); batch evaluation; macro embedding for ≤3 layers.
- **Known Issues**: `folded_gpu` crashes (missing folded buffers); parity with GridFF needs audit (Coulomb scaling).
- **References**: `doc/Topics/OnSurfaceAssembly/ElectrostaticContinuumEmbeding_report.md`.

### 2.3 FoldedAtomicFunctions (Basis Fitting for Surface Potentials)

#### 2.3.1 What It Does
- Constructs separable basis functions (plane waves × exponentials) for 2D surface potentials.
- Fits Morse/Coulomb samples via regularized least-squares.
- Intended as a compact alternative to GridFF; experimental.

#### 2.3.2 Key Files and Flow
- **Basis Generation**: `doc/py/FoldedAtomicFunctions/FoldedAtomicFunctions.py:BasisFunctions.generate_plane_wave_exponential_basis()`
  - Creates φ(x,z) = cos(2πnx/L)·exp(-αz) for n=0..3, α=1..4 [Trace 1b].
  - Stores basis grids in list.
- **Potential Calculation**: `PotentialCalculator.calculate_morse_periodic()`
  - Sums Morse over periodic images (±3) on grid [Trace 1c].
- **Fitting**: `PotentialFitter.fit()`
  - Builds matrix A from basis; solves (AᵀA+λI)⁻¹Aᵀb [Trace 2d, 2e].
  - Reconstructs potential via linear combination [Trace 2f, 2g].
- **Exploration**: `svd_basis_zcut_scan_new.py`
  - Generates random Morse curves; fits with cutoff-polynomial basis and repulsive masking [Traces 4a-4e].

#### 2.3.3 Status
- **Prototype**: Basis generation and fitting work; not yet integrated with GPU kernels.
- **Missing**: GPU evaluation kernel; integration with rigid scanning/surface sampling.
- **References**: `doc/py/FoldedAtomicFunctions/` (all files).

### 2.4 Rigid-Body Assembly (Multi-Molecule Pose Generation)

#### 2.4.1 What It Does
- Generates massive numbers of rigid-body poses (rotations + translations) with symmetry and PBC.
- Scores poses via GPU clash detection and min-distance evaluation.
- Exports best configurations as XYZ/mol2.

#### 2.4.2 Key Files and Flow
- **Pose Generation**: `tests/tMMFF/test_assembly.py:custom_generate_transform_buffer()`
  - Builds rotation grid (`generate_rotations()`) [Trace 1b].
  - Builds translation grid; applies 6-fold symmetry [Trace 1c].
  - Tiles PBC replicas [Trace 1d].
- **GPU Scoring**: `pyBall/OCL/Assembly.evaluate_packing()`
  - Packs transforms; launches `evaluate_packing_3d` kernel [Trace 1e].
  - Kernel loops replicas; performs clash check [Trace 1f]; reduces and writes results [Trace 1g].
- **Export**: `AssemblyOCL.emit_configuration()`
  - Launches `emit_configuration_xyz` kernel to apply transforms and output coordinates [Trace 2e, 2f].

#### 2.4.3 Status
- **Works**: Pose generation, GPU scoring, export for large systems (e.g., PTCDA on NaCl).
- **References**: `tests/tMMFF/test_assembly.py`, `pyBall/OCL/Assembly.py`, `cpp/common_resources/cl/Assembly.cl`.

### 2.5 Surface Sampling (Iso-Height Extraction)

#### 2.5.1 What It Does
- Samples (x,y) grid; for each point finds z where a potential crosses a threshold or reaches its first minimum.
- Returns height map z(x,y) and optional color map (e.g., Coulomb at z_node).
- Supports GridFF and XYZ backends.

#### 2.5.2 Key Files and Flow
- **Python Sampler**: `pyBall/SurfaceSampling.sample_surface_map()`
  - Loops over xy grid [Trace 4b]; evaluates potential along z via backend [Trace 4c].
  - For GridFF: `_interp_channel()` does trilinear interpolation in Python (slow) [Trace 4d].
  - Finds z_node via `refine_first_minimum()` (Python loop) [Trace 4e].
  - Samples color at z_node [Trace 4f]; returns SurfaceMapResult.
- **GPU Kernels (Fast Path)**: `MolecularDynamics.eval_surface_iso_gridff()` / `eval_surface_iso_morse()`
  - GridFF: `getSurfaceIsoGridFF` kernel scans z for each (x,y) in parallel.
  - Morse: `getSurfaceIsoSurfMorse` kernel similarly scans z.
  - Return `points_world` (xyz) and `ok_mask`.

#### 2.5.3 Status
- **Partial**: Python path works but slow; GPU kernels implemented and used in GUI/headless.
- **Known Issues**: Python sampler bottleneck; GridFF vs. Morse parity discrepancies; seam artifacts on cropped cells.
- **References**: `pyBall/SurfaceSampling.py`, `pyBall/OCL/MolecularDynamics.py`, `doc/Topics/OnSurfaceAssembly/SurfacePotentialVisualization_*`.

### 2.6 MMFF Force Field (Surface Relaxation)

#### 2.6.1 What It Does
- Builds linearized molecular mechanics topology (bonds, angles, pi/epair dummies).
- Performs GPU-accelerated MD on surfaces with nonbonded + surface forces.
- Supports rigid and relaxed scans.

#### 2.6.2 Key Files and Flow
- **Topology Builder**: `pyBall/OCL/MMFFL.build_topology()`
  - Assigns atom types via table/UFF [Trace 5a].
  - Detects bonds/angles; converts angles to harmonic bonds [Trace 5b].
  - Adds pi/epair dummies; classifies bonds by tag [Trace 5c].
  - Returns topology dict [Trace 5d].
- **Surface MD**: `tests/tMMFF/test_ditetraceno_surface.py:run_ditetraceno()`
  - Initializes MMFF; builds topology [Trace 6a, 6b].
  - Creates GPU MD harness; uploads MMFF [Trace 6c, 6d].
  - MD loop: clean forces; evaluate nonbonded [Trace 6e]; add surface interaction [Trace 6f]; add bonded forces; integrate [Trace 6g].

#### 2.6.3 Status
- **Works**: Surface MD for molecules like ditetraceno; MMFF topology robust.
- **References**: `tests/tMMFF/test_ditetraceno_surface.py`, `pyBall/OCL/MMFFL.py`.

### 2.7 GUI Placement Tools

#### 2.7.1 What It Does
- Interactive VisPy interfaces for placing molecules on substrates.
- Supports sequence placement (SequencePlacerVisPy) and symmetric flower patterns (MolecularPlacerVisPy).

#### 2.7.2 Key Files and Flow
- **SequencePlacerVisPy**: Builds NaCl substrate; places molecules along rows with rotations; updates VisPy markers.
  - Substrate generation: `gen_nacl_step()` [Trace 3a].
  - Placement: `place_sequence()` calls `make_rotmat()` and computes shifts [Trace 3b-3d].
  - Visual update: `mol_markers.set_data()` [Trace 3f].
- **MolecularPlacerVisPy**: N-fold symmetric copies; hex lattice shifts; PBC replication.
  - Flower generation: `_generate_flower()` with rotations and offsets [Trace 4b-4c].
  - PBC shifts: `_get_all_copies()` [Trace 4d].
  - Visual update: concatenates copies and renders [Trace 4e].

#### 2.7.3 Status
- **Works**: Both GUIs functional; real-time placement and visualization.
- **References**: `pyBall/SequencePlacerVisPy.py`, `pyBall/MolecularPlacerVisPy.py`, `pyBall/MolGUI_common.py`.

### 2.8 Path Optimization (Batch Replica with Tethers)

#### 2.8.1 What It Does
- Initializes many replicas of a molecule; applies substrate GridFF forces.
- Sets causal tethers between replicas for path optimization.
- Runs relaxation loops.

#### 2.8.2 Key Files and Flow
- **Initialization**: `pyBall/OCL/ManipulationPathOpt.__init__()`
  - Builds MMFF topology [Trace 7a]; creates MD harness with batch allocation [Trace 7b].
  - Packs MMFF to all replicas [Trace 7b].
- **Substrate**: `set_substrate_gridff()` uploads GridFF [Trace 7c].
- **Replica States**: `init_replica_states()` uploads shifted positions [Trace 7d].
- **Tethers**: `set_causal_tethers()` builds sysneighs/sysbonds arrays [Trace 7e].
- **Relax Loop**: `run_relax()` calls `_run_tip_morse()` kernel [Trace 7f].

#### 2.8.3 Status
- **Works**: Batch initialization and tether setup; relaxation loops functional.
- **References**: `pyBall/OCL/ManipulationPathOpt.py`.

## 3. Cross-Subsystem Workflows

### 3.1 Rigid Scanning Workflows
1. **GridFF Path**: Generate GridFF → Load → `eval_surface_iso_gridff()` → height map.
2. **XYZ Path**: `InteractionScanner` (reference/fast_gpu) → 1D/2D scans → energy curves.
3. **Batch Path**: Generate transforms → `eval_rigid_getSurfMorse_components()` → per-config energies.

### 3.2 Surface MD Workflows
1. Build MMFF topology → upload to GPU → add substrate forces (GridFF or flat) → MD loop.

### 3.3 Assembly + Scanning
1. Generate poses via Assembly → export best → scan each pose with InteractionScanner or GridFF.

### 3.4 GUI + Scanning
1. Place molecules via VisPy GUI → extract pose → run rigid scan → visualize energy.

## 4. Status Summary (What Works vs. What Doesn't)

| Subsystem | Works | Issues/Gaps | References |
|-----------|-------|-------------|------------|
| GridFF | Generation, loading, GPU force evaluation, batch | Coulomb scaling vs. InteractionScanner; parity tests needed | `tests/tMMFF/run_test_GridFF_CaF2.py`, `doc/Topics/OnSurfaceAssembly/SurfacePotentialVisualization_*` |
| InteractionScanner | 1D/2D rigid scans (CPU/GPU), macro embedding, batch | `folded_gpu` crashes; parity with GridFF | `tests/tMMFF/test_interaction_scan.py`, `pyBall/OCL/InteractionEnergy.py` |
| FoldedAtomicFunctions | Basis generation, fitting (Python) | No GPU kernel; not integrated with scanning | `doc/py/FoldedAtomicFunctions/` |
| Rigid-Body Assembly | Pose generation, GPU scoring, export | None | `tests/tMMFF/test_assembly.py`, `pyBall/OCL/Assembly.py` |
| Surface Sampling | Python sampler (slow), GPU kernels (GridFF/Morse) | Python bottleneck; GridFF vs. Morse parity; seam artifacts | `pyBall/SurfaceSampling.py`, `pyBall/OCL/MolecularDynamics.py` |
| MMFF | Topology building, surface MD | None | `tests/tMMFF/test_ditetraceno_surface.py`, `pyBall/OCL/MMFFL.py` |
| GUIs | SequencePlacer, MolecularPlacer | None | `pyBall/SequencePlacerVisPy.py`, `pyBall/MolecularPlacerVisPy.py` |
| Path Optimization | Batch replica setup, tethers, relaxation | None | `pyBall/OCL/ManipulationPathOpt.py` |

## 5. Key Scripts and How to Run Them

### 5.1 GridFF Generation and Testing
```bash
# Generate GridFF for CaF2
cd pyBall/tests
python ocl_GridFF_new.py --substrate ../../cpp/common_resources/Substrates/generated_rect/CaF2_6L_Ni3_rect_nx2_nz1_L2_top.xyz --sub_types Ca=Ca+2,F=F- --save_grid

# Test GridFF loading and vertical scan
cd ../../tests/tMMFF
python run_test_GridFF_CaF2.py
```

### 5.2 Rigid Scans (InteractionScanner)
```bash
# 1D scan for H2O on NaCl
cd tests/tMMFF
python test_interaction_scan.py --mol H2O_O.xyz --sub NaCl_8x8_L3.xyz --scan_type 1D --probe H --charge 0.2

# 2D scan with fast_gpu
python test_interaction_scan.py --mol H2O_O.xyz --sub NaCl_8x8_L3.xyz --scan_type 2D --backend fast_gpu
```

### 5.3 Batch Rigid Evaluation
```bash
# Batch evaluation for PTCDA on NaCl
cd tests/tMMFF
python test_rigid_gridff_ptcda_batch.py --mol PTCDA.xyz --sub NaCl_8x8_L3.xyz --nrot 36 --nshift 5
```

### 5.4 Surface Sampling (Headless)
```bash
# Generate surface height map (GridFF)
cd tests/tMMFF
python render_surface_iso.py --source gridff --gridff ../../cpp/common_resources/CaF2_6L_Ni3_rect_nx2_nz1_L2_top/Bspline_PLQd.npy --probe H --charge 0.2 --mode first_minimum --selector nonel --color coulomb --nx 81 --ny 71 --nz 160

# Generate surface height map (XYZ fast_gpu)
python render_surface_iso.py --source xyz --sub ../../cpp/common_resources/xyz/NaCl_1x1_L3.xyz --xyz_backend fast_gpu --probe H --charge 0.2 --mode first_minimum --selector nonel --color coulomb --nx 81 --ny 81 --nz 160
```

### 5.5 GUI Placement
```bash
# Sequence placer
python pyBall/ExplorerVisPy.py --mol H2O_O.xyz --sub NaCl_1x1_L3.xyz

# Molecular placer (symmetric)
python pyBall/MolecularPlacerVisPy.py --mol PTCDA.xyz --sub NaCl_8x8_L3.xyz
```

### 5.6 Assembly
```bash
cd tests/tMMFF
python test_assembly.py --mol PTCDA.xyz --sub NaCl_8x8_L3.xyz --nrot 36 --nshift 5 --export_max 100
```

## 6. Force Field Comparisons and Parity Issues

### 6.1 GridFF vs. InteractionScanner (XYZ)
- **Discrepancy**: Coulomb channel stronger in GridFF; suspected double application of 14.39 eV·Å factor.
- **Action**: Audit generation (`ocl_GridFF_new.py`) and sampling (`getNonBond_GridFF_Bspline`) for Coulomb scaling; rebuild GridFF if needed.
- **Reference**: `doc/Topics/OnSurfaceAssembly/SurfacePotentialVisualization_discussion.md` (Current discrepancies section).

### 6.2 Folded vs. GridFF vs. XYZ
- **Folded**: Not yet integrated; basis fitting works but lacks GPU kernel.
- **GridFF**: Fast but needs parity validation.
- **XYZ**: Baseline physics; slower but validated.

### 6.3 Electrostatic Embedding
- **Macro Embedding**: Enabled automatically for ≤3 distinct z-layers; improves Coulomb accuracy for thin slabs.
- **Continuum Embedding**: Documented in `ElectrostaticContinuumEmbeding_report.md`; not yet integrated into scanning workflows.

## 7. Recommendations and Next Steps

1. **Parity Validation**: Implement quantitative Δz/ΔE maps between GridFF and XYZ for NaCl and CaF2; audit Coulomb scaling.
2. **Folded Integration**: Develop GPU kernel for folded basis evaluation; integrate with rigid scanning.
3. **Surface Sampling Optimization**: Replace Python loops with GPU batch evaluation; fix seam artifacts.
4. **Macro Embedding Extension**: Enable macro for thicker slabs via layer grouping.
5. **Unified Interface**: Create a high-level API that selects backend (GridFF/XYZ/Folded) and returns standardized energy/force outputs.

## 8. Document References

- **Discussions/Reports**: 
  - `doc/Topics/OnSurfaceAssembly/SurfacePotentialVisualization_discussion.md` (detailed reproducibility and issues)
  - `doc/Topics/OnSurfaceAssembly/SurfacePotentialVisualization_report.md`
  - `doc/Topics/OnSurfaceAssembly/ElectrostaticContinuumEmbeding_report.md`
  - `doc/Topics/RigidBodyAssembly/RigidGridFF_H2O_NaCl_report.md`
  - `tests/tMMFF/GridFF_CaF2_doc_tutorial.md`
  - `cpp/common_resources/CaF2_6L_Ni3_rect_nx2_nz1_L2_top/CaF2_report.md`
- **Code Maps**: Refer to provided codemaps for detailed flow traces.

## 9. VisPy GUI Programs

### 9.1 Overview of Visual Programs

The FireCore project includes several VisPy-based GUI programs for molecular visualization and manipulation on surfaces. Each serves a distinct purpose and uses different backend modules.

### 9.2 ExplorerVisPy.py

**Purpose**: Interactive molecule-substrate interaction explorer with real-time energy scanning and surface visualization.

**Key Features**:
- Real-time molecule pose manipulation (translation/rotation)
- Surface potential scanning with multiple backends (Reference, Fast GPU, GridFF, Folded GPU)
- Interactive 3D visualization of molecules, substrates, and surface isopotentials
- Energy landscape plotting and analysis

**Backend Modules Imported**:
```python
from pyBall.SurfaceSampling import ProbeParticle, build_surface_map
from pyBall.VispyUtils import make_surface_mesh
from pyBall.MolGUI_common import ELEMENT_COLORS, ELEMENT_SIZES, colors_for_enames, sizes_for_enames, find_bonds
```

**Force Field Backends**:
- **XYZ Reference**: Uses `InteractionScanner` (CPU) - slow but baseline physics
- **XYZ Fast GPU**: Uses `MolecularDynamics.eval_surface_iso_morse()` → OpenCL kernel `getSurfaceIsoSurfMorse`
- **XYZ Folded GPU**: Intended for folded kernels (currently crashes)
- **GridFF**: Uses `MolecularDynamics.eval_surface_iso_gridff()` → OpenCL kernel `getSurfaceIsoGridFF`

**Unique Capabilities**:
- Multi-backend surface scanning with headless rendering support
- Real-time energy evaluation during molecule manipulation
- Surface mesh generation and visualization
- Probe particle configuration (charge, polarizability)

**Usage**:
```bash
python ExplorerVisPy.py                              # default H2O on NaCl
python ExplorerVisPy.py --mol PTCDA.xyz --sub NaCl_8x8_L3.xyz --type_map C=C_R
```

### 9.3 VispyUtils.py

**Purpose**: Reusable VisPy widget library for molecular visualization with orthographic top-down view and picking/dragging capabilities.

**Key Features**:
- Generic MD/FF viewer for atoms, bonds, forces, and various overlays
- 2D/3D picking and dragging with atom constraints
- Group-based visualization and coloring
- Extensive overlay support (bboxes, links, ports, debug vectors)

**Backend Modules Imported**:
```python
import vispy
from vispy import scene
from vispy.scene import visuals
from vispy.color import Colormap
```

**Force Field Backends**: None (pure visualization library)

**Unique Capabilities**:
- Pseudo-2D drag mode for top-down manipulation
- Comprehensive overlay system for debugging (bboxes, links, ports, force vectors)
- Group-based rendering with masks
- Camera controls (zoom, pan, rotation)

**Core Classes**:
- `AtomScene`: Main visualization widget with picking/dragging
- Supports multiple visual layers: atoms, bonds, forces, bboxes, links, labels

### 9.4 SequencePlacerVisPy.py

**Purpose**: Interactive placement of molecular sequences on NaCl step-edge substrates with instant visual feedback.

**Key Features**:
- Build ionic substrates with configurable steps
- Place sequences of molecules with rotations and spacing
- Real-time visualization of placed molecules and bonds
- Export to XYZ format

**Backend Modules Imported**:
```python
from pyBall import SubstrateBuilder as SB
from pyBall.MolGUI_common import (
    ELEMENT_COLORS, ELEMENT_SIZES, colors_for_enames, sizes_for_enames,
    rotmat_x, rotmat_y, rotmat_z, make_rotmat, load_molecule_xyz, save_xyz,
    find_bonds, replicate_bonds, place_sequence, load_default_molecules, DEFAULT_MOLS,
)
```

**Force Field Backends**: None (geometric placement only)

**Unique Capabilities**:
- Step-edge substrate generation
- Sequential molecule placement with custom sequences
- Bond replication across periodic boundaries
- Headless mode support for batch processing

**Usage**:
```bash
python SequencePlacerVisPy.py                    # interactive GUI
python SequencePlacerVisPy.py --headless ...     # headless mode
```

### 9.5 MolecularPlacerVisPy.py

**Purpose**: Interactive molecular placement tool with symmetry operations, N-fold rotation, hexagonal lattice placement, and dual-flower patterns.

**Key Features**:
- Molecular rotation and symmetry operations
- Hexagonal lattice placement with triangular coordinates
- Dual-flower patterns with mirroring options
- Periodic boundary condition replication

**Backend Modules Imported**:
```python
from pyBall.MolGUI_common import (
    ELEMENT_COLORS, ELEMENT_SIZES, colors_for_enames, sizes_for_enames,
    rotmat_x, rotmat_y, rotmat_z, make_rotmat, flip_matrix_x, flip_matrix_y,
    load_molecule_xyz, save_xyz, find_bonds, HexLattice,
)
```

**Force Field Backends**: None (geometric placement only)

**Unique Capabilities**:
- N-fold symmetry generation
- Hexagonal lattice integration
- Dual-flower pattern creation
- PBC replication with bond generation

**Usage**:
```bash
python MolecularPlacerVisPy.py [xyz_file]
```

### 9.6 test_RRsp3_vispy.py

**Purpose**: Debug and visualization tool for RRsp3 (Rigid-Rigid Spring System) molecular dynamics with collision detection and port-based interactions.

**Key Features**:
- RRsp3 molecular dynamics simulation with real-time visualization
- Collision detection with configurable radii
- Port-based inter-molecular connections
- Group-based topology management
- Interactive picking/dragging with constraints

**Backend Modules Imported**:
```python
from RRsp3 import RRsp3, build_neighs_bk_from_bonds, make_bk_slots_clustered, make_exclusions_1st_2nd
from XPTB_utils import pack_molecules_contiguous, make_h2o_geometry, masses_from_elems, load_xyz, perturb_state
from MolGUI_common import colors_for_enames, sizes_for_enames
from VispyUtils import AtomScene
```

**Force Field Backends**:
- **RRsp3**: Custom rigid-body dynamics with spring connections
- **Collision Detection**: Geometric sphere-based collision forces
- **Port Forces**: Directional spring forces between designated ports

**Unique Capabilities**:
- Rigid-body molecular dynamics visualization
- Real-time collision detection and response
- Port-based molecular assembly
- Group-based topology with padding
- Interactive manipulation during simulation

**Key Components**:
- Molecule packing with group-based padding
- Neighbor list and exclusion management
- Port local geometry and spring constants
- Bounding box and topology visualization

### 9.7 Comparison Summary

| Program | Primary Purpose | Force Field Support | Unique Features | Backend Dependencies |
|---------|------------------|-------------------|-----------------|---------------------|
| ExplorerVisPy | Molecule-substrate interaction scanning | GridFF, Morse, Folded | Real-time energy evaluation, surface visualization | SurfaceSampling, MolecularDynamics |
| VispyUtils | Generic molecular visualization library | None (pure viz) | Picking/dragging, extensive overlays | vispy only |
| SequencePlacer | Sequential placement on step-edges | None (geometric) | Step-edge substrates, sequences | SubstrateBuilder, MolGUI_common |
| MolecularPlacer | Symmetry-based molecular placement | None (geometric) | N-fold symmetry, hex lattice, dual-flower | MolGUI_common |
| test_RRsp3_vispy | RRsp3 dynamics debugging | RRsp3, collision, ports | Rigid-body MD, real-time manipulation | RRsp3, XPTB_utils, VispyUtils |

**Common Patterns**:
- All use PyQt5 + VisPy for 3D visualization
- Import from `MolGUI_common` for basic molecular utilities
- Support both interactive and headless modes where applicable
- Use similar UI patterns with controls panel and 3D canvas

**Backend Integration Levels**:
- **Full Physics**: ExplorerVisPy (multiple force field backends)
- **Custom Physics**: test_RRsp3_vispy (specialized RRsp3 system)
- **Geometric Only**: SequencePlacer, MolecularPlacer (placement tools)
- **Visualization Only**: VispyUtils (reusable library)

**End of Overview**
