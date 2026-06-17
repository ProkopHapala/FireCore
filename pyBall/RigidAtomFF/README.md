# RigidAtomFF - Rigid Atom Force Field Implementations

This directory contains implementations of Position-Based Dynamics (PBD) and Extended Position-Based Dynamics (XPBD) for molecular simulations with rigid body constraints.

## Directory Structure

```
RigidAtomFF/
├── RRsp3/           # Active implementation (cluster-sorted rigid ports + collisions)
├── XPDB_new/        # Experimental implementation (force-based + position-based)
├── XPDB_legacy/     # Deprecated simple PBD implementation
├── XPBD_2D/         # 2D Position-Based Dynamics implementation
└── shared/          # Shared utilities and documentation
```

## Subdirectories

### RRsp3/ (Active Development)

**Primary implementation** - Cluster-sorted rigid body dynamics with "ports" (ARAP-style constraints) and collision detection.

**Core Files:**
- `RRsp3.cl` - OpenCL kernels (1311 lines)
  - `update_bboxes_rigid` - Per-cluster AABB for broad-phase
  - `build_local_topology_rigid` - Ghost atom finding and local index remapping
  - `compute_collision_cluster_rigid` - Cluster-local collisions in LDS
  - `compute_ports_cluster_rigid` - Node-only ARAP/ports with multiple rotation solvers
  - `apply_corrections_rigid_ports` - Gathers collision + node + recoil deltas

- `RRsp3.py` - Python harness (624 lines)
  - Buffer allocation and management
  - State upload/download (pos, quat, fixmask, radius)
  - Topology upload (neighs, exclusions, ports, stiffness, bkSlots)
  - Solver execution (step_cluster, compute_cluster_deltas, apply_cluster_corrections)
  - Momentum management

**Rotation Solver Variants:**
- `compute_ports_cluster_rigid_orig` - One-pass XPBD, lightest but under-drives rotation
- `compute_ports_cluster_rigid_substep_optimized` - Newton substeps, cheap but not momentum-conserving
- `compute_ports_cluster_rigid_shapematch` - Polar decomposition, accurate but heavy
- `compute_ports_cluster_rigid_eigen` - Davenport q-method, lighter than shapematch

**Test Scripts:**
- `test_RRsp3_convergence.py` - Convergence testing with parameter sweeps
- `test_RRsp3_debug.py` - Kernel debug log parsing
- `test_RRsp3_momentum.py` - Momentum conservation validation
- `test_RRsp3_smoke.py` - Basic smoke tests
- `test_RRsp3_vispy.py` - Interactive Vispy GUI

**Documentation:**
- `RRsp3.progress.md` - Development history and session summaries
- `RRsp3_momentum_design.md` - Momentum conservation design
- Various `.chat.md` files - Development discussions

**Key Features:**
- Cluster-sorted layout (64 atoms per workgroup)
- Nodes first, then caps, then padding within each cluster
- Ghost atom halo for inter-cluster interactions
- Recoil buffers (dpos_neigh) for momentum conservation
- Jacobi iteration with heavy-ball momentum acceleration
- 1-2 and 1-3 neighbor exclusion for collisions

### XPDB_new/ (Experimental)

**Experimental implementation** - Larger version with both force-based and position-based dynamics paths.

**Core Files:**
- `XPDB_new.cl` - OpenCL kernels (3584 lines)
- `XPDB_new.py` - Python harness (1303 lines)

**Test Scripts:**
- `test_TiledJacobi_molecules.py` - Tiled Jacobi solver testing
- `test_rigid_XPBD_molecules.py` - Rigid XPBD testing on molecules
- `test_XPDB_new_dump_headless.py` - Headless testing
- `test_XPDB_new_relax_benchmark.py` - Relaxation benchmarking
- `test_procrustes_ch4.py` - Procrustes rotation testing on CH4
- `ffdebug_xpdb_rigid_viz.py` - Debug visualization for rigid XPDB
- `ffdebug_xpdb_tiled_viz.py` - Debug visualization for tiled XPDB
- `test_viz_rigid_XPBD.py` - Rigid XPBD visualization testing

**Status:** Superseded by RRsp3 for cluster-sorted PBD. Contains experimental features not yet migrated.

### XPDB_legacy/ (Deprecated)

**Older implementation** - Simple PBD without rigid bodies, quaternions, or clustering.

**Core Files:**
- `XPDB.cl` - OpenCL kernels (1130 lines)
- `XPDB.py` - Python harness (734 lines)

**Test Scripts:**
- `test_XPDB.py` - Basic XPDB testing
- `test_Jacobi_Convergence.py` - Jacobi convergence analysis
- `test_Tiled_Jacobi.py` - Tiled Jacobi testing
- `test_TiledJacobi_new.py` - New tiled Jacobi variant
- `test_XPDB_relax_benchmark.py` - XPDB relaxation benchmarking

**Status:** Historical reference only. Not actively maintained.

### XPBD_2D/ (2D Implementation)

**2D Position-Based Dynamics** - Specialized implementation for 2D simulations.

**Core Files:**
- `XPBD_2D.cl` - OpenCL kernels for 2D PBD
- `XPBD_2D.py` - Python harness for 2D PBD
- `XPBD_2D_utils.py` - Utility functions for 2D simulations

**Test Scripts:**
- `test_XPBD_2D.py` - Main 2D PBD testing
- `test_viz_XPBD_2D.py` - 2D visualization testing
- `ffdebug_xpbd2d_viz.py` - Debug visualization for 2D XPBD

**Documentation:**
- `XPBD_2D.md` - 2D implementation documentation
- `XPBD_2D_progress.chat.md` - Development progress

**Status:** Specialized for 2D simulations, separate from 3D implementations.

### shared/ (Shared Utilities)

**Utilities and documentation** used across implementations.

**Python Utilities:**
- `XPTB_utils.py` - Molecule packing, trajectory I/O, permutation utilities
  - `pack_molecules_contiguous()` - Pack molecules into cluster-sorted layout
  - `load_xyz()` - Load XYZ files
  - `masses_from_elems()` - Get atomic masses
  - `write_xyz_with_ports()` - Write XYZ with port information
  - `write_pdb_trajectory()` - Write PDB trajectories

**Visualization:**
- `generate_viz.py` - Visualization generation for scan outputs (implementation-agnostic)

**Documentation:**
- `Analytic_Procrustes_doc.md` - Rotation optimization theory
- Various `.chat.md` files - Historical discussions

## Import Path Updates

After reorganization, Python imports need to be updated:

**Old imports:**
```python
from XPDB_AVBD.RRsp3 import RRsp3
from XPDB_AVBD.XPTB_utils import pack_molecules_contiguous
```

**New imports:**
```python
from pyBall.RigidAtomFF.RRsp3.RRsp3 import RRsp3
from pyBall.RigidAtomFF.shared.XPTB_utils import pack_molecules_contiguous
```

Or add to sys.path:
```python
import sys
sys.path.append('pyBall/RigidAtomFF/RRsp3')
sys.path.append('pyBall/RigidAtomFF/shared')
```

## Running Tests

### RRsp3 Tests
```bash
cd pyBall/RigidAtomFF/RRsp3
python test_RRsp3_convergence.py --scan
python test_RRsp3_debug.py
python test_RRsp3_momentum.py
python test_RRsp3_vispy.py
```

### XPDB_new Tests
```bash
cd pyBall/RigidAtomFF/XPDB_new
python test_TiledJacobi_molecules.py
python test_rigid_XPBD_molecules.py
```

### XPDB_legacy Tests
```bash
cd pyBall/RigidAtomFF/XPDB_legacy
python test_XPDB.py
python test_Jacobi_Convergence.py
```

### XPBD_2D Tests
```bash
cd pyBall/RigidAtomFF/XPBD_2D
python test_XPBD_2D.py
python test_viz_XPBD_2D.py
```

## Key Concepts

### Cluster-Sorted Layout
- Atoms packed into groups of 64 (GROUP_SIZE)
- Within each cluster: nodes first, then caps, then padding
- Enables efficient local memory access and reduced divergence

### Ghost Atoms / Halo Approach
- External atoms within interaction margin stored as ghosts
- Global indices remapped to local (0-63 cluster, 64+ ghosts)
- Enables pure local memory access in collision/port kernels

### Recoil Buffers
- `dpos_neigh` stores Newton's 3rd law recoil forces
- Gather approach avoids atomic writes on GPU
- Critical for momentum conservation in Jacobi scheme

### Rotation Solvers
- **XPBD**: Position-based dynamics with impulse-based updates
- **Shapematch**: Polar decomposition for optimal rotation (massless)
- **Eigen**: Davenport q-method via quaternion power iteration (massless)
- **Substep**: Newton substeps for aggressive rotational tightening

### Momentum Conservation
- Linear: ΔP = Σ m·dx = 0
- Angular: ΔL = Σ r×(m·dx) + Σ I·dθ = 0
- Verified by `test_RRsp3_momentum.py`

## Known Issues

### Slow Convergence
- Current: ~73-140 iterations to reach 1e-2 error
- Expected: 3-5 iterations (game-physics style)
- Root cause: Jacobi diffusion on linear chains (O(N²) information propagation)
- Proposed solution: Inner local iterations (Block Gauss-Seidel) in LDS

### Collision Mapping
- Need to verify ghost atom building and local index remapping
- Exclusion lists (1-2, 1-3 neighbors) may need validation

## Development Notes

### File Chronology
1. **XPDB.cl/XPDB.py** - Original simple PBD (deprecated)
2. **XPDB_new.cl/XPDB_new.py** - Experimental version with force-based + position-based
3. **RRsp3.cl/RRsp3.py** - Consolidated cluster-sorted PBD (active)

### Migration Path
XPDB_new.cl → RRsp3.cl (extracted cluster PBD path only, dropped force-based path)

### Subproblem Independence
- **Collisions (Subproblem A)**: Can be tested independently with `ENABLE_PORT=0`
- **Ports/Rotations (Subproblem B)**: Can be tested independently with `ENABLE_COLL=0`
- Shared: `build_local_topology_rigid`, `apply_corrections_rigid_ports`

## References

- Projective Dynamics papers
- Augmented Vertex Block Descent (AVBD)
- Position-Based Dynamics (PBD)
- Quaternion-based rotation optimization
