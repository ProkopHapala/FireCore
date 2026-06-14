# GPU Flat Surface Potential (Hamaker + Morse) — PyOpenCL Implementation

## Background

The objective was to implement and test GPU flat surface potentials (Hamaker LJ9-3 and Morse) in MMFF simulations using PyOpenCL. This involved:

1. Implementing the potential in OpenCL kernels
2. Updating Python bindings in MolecularDynamics.py
3. Creating test scripts for validation
4. Running relaxation tests on small (H2O) and large (DiTetracenoHelice, 624 atoms) molecules
5. Using PBC with nPBC=(1,1,0) for the large molecule

## Relevant Files and Functions

### OpenCL Kernels

**`cpp/common_resources/cl/relax_multi.cl`**
- `getSurfFlat` kernel (lines ~1328-1403): Main GPU kernel for flat surface potential evaluation
- Helper functions:
  - `combineREQ(float4 a, float4 b)`: Combines surface and atom REQ parameters
  - `getHamakerLJ93(float3 dp, float3 n, float3* f, float4 REQij)`: Hamaker LJ9-3 potential
  - `getMorseSurface(float3 dp, float3 n, float3* f, float4 REQij, float K)`: Morse surface potential

**`cpp/common_resources/cl/relax_multi_mini.cl`**
- Same `getSurfFlat` kernel and helpers (lines ~611-684) — this file is smaller but contains the same kernel

**`cpp/common/OpenCL/OCL_MM.h`**
- `setup_getSurfFlat` method (lines ~186-187): Registers the kernel in the kernel factory
- Kernel registration in `makeKernels_MM` (lines ~793-819)

### Python Bindings

**`pyBall/OCL/MolecularDynamics.py`**
- `run_getSurfFlat()` method (lines ~403-): Executes the surface force kernel
- `update_pbc_shifts(lvec, nPBC)` (lines ~584-615): Generates PBC shift vectors
- `pack_system(iSys, mmff)` (lines ~617-): Uploads MMFF data to GPU, including lattice vectors and PBC shifts
- `run_cleanForceMMFFf4()` (lines ~392-401): Clears forces, now also zeros full fneigh buffer
- `allocate_cl_buffers()` (lines ~80-160): Buffer allocation — fixed fneigh size and buffer zeroing

**`pyBall/OCL/MMFF.py`**
- `realloc(nnode, ncap, ntors=0, nPBC=(0,0,0))` (lines ~154-): Allocates MMFF arrays, now passes nPBC through
- `toMMFFsp3_loc(mol, atom_types, ...)` (lines ~281-): Converts AtomicSystem to MMFF, now reads nPBC before realloc
- Fixed `neighCell` initialization to use correct zero-shift index

**`pyBall/atomicUtils.py`**
- `loadMol2(fname, ...)` (lines ~1290-): Now parses `#lvs` comment lines in mol2 files

### Test Scripts

**`tests/tMMFF/test_flat_surface_gpu.py`**
- Single-atom scan test: places one O atom at varying heights above surface
- Compares GPU energy/force to analytical Hamaker and Morse formulas
- Validates kernel correctness

**`tests/tMMFF/test_h2o_surface_relax.py`**
- H2O relaxation on flat surface
- Tests both Hamaker (mode=1) and Morse (mode=2)
- Molecule drops from z=4.5 to z≈3.25 (near R_eff), demonstrating surface binding

**`tests/tMMFF/test_ditetraceno_surface.py`**
- DiTetracenoHelice (624 atoms) relaxation on flat surface
- Uses nPBC=(1,1,0), reads lattice from mol2 `#lvs` comment
- Saves trajectory every 100 steps
- Uses `force_node_all` mode for stability

## Problems Encountered

### 1. Python IndentationError in MolecularDynamics.py

**Problem**: Indentation bug at lines 685-687 in PBC shift upload code.

**Solution**: Fixed indentation to properly handle the conditional block.

### 2. Missing getSurfFlat Kernel in relax_multi.cl

**Problem**: `getSurfFlat` was only in `relax_multi_mini.cl`, but `MolecularDynamics.py` loads `relax_multi.cl`.

**Solution**: Added `getSurfFlat` kernel and helper functions to `relax_multi.cl`.

### 3. Duplicate pack_system in MolecularDynamics.py

**Problem**: Two definitions of `pack_system` — Python uses the last one, but the first was dead code causing confusion.

**Solution**: Commented out the first (dead) definition.

### 4. Kernel Argument Generation Errors

**Problem**: `setup_kernels` failed when buffers were missing for non-essential kernels.

**Solution**: Wrapped kernel argument generation in try/except to handle missing buffers gracefully.

### 5. #lvs Comment Not Parsed in mol2 Files

**Problem**: `loadMol2` only handled `@LVS` but DiTetracenoHelice mol2 uses `#lvs`.

**Solution**: Added parsing for `#lvs`/`#LVS` prefix in `atomicUtils.py`.

### 6. fneigh Buffer 4x Too Small (CRITICAL)

**Problem**: `fneigh` buffer was allocated with `float_size` (4 bytes) instead of `float4_size` (16 bytes). Kernel accesses `nnode*2*4` float4 entries, but buffer was only `nnode*2*4*4` bytes → buffer overrun, NaN.

**Solution**: Changed to use `float4_size = 4 * float_size` in `MolecularDynamics.py:124-125`.

### 7. Uninitialized GPU Buffers (CRITICAL)

**Problem**: `avel`, `cvf`, `aforce_old`, `fneigh`, `constr`, `constrK`, `TDrives` were never zeroed after allocation. Large GPU allocations don't guarantee zeroed memory → NaN in velocity updates.

**Solution**: Added buffer zeroing after allocation in `allocate_cl_buffers()` (lines ~154-160).

### 8. nPBC Not Passed Through to realloc (CRITICAL)

**Problem**: Setting `mm.nPBC = (1,1,0)` before `toMMFFsp3_loc` had no effect because `realloc` wasn't receiving it.

**Solution**: Modified `toMMFFsp3_loc` to read `getattr(self, 'nPBC', None) or (0,0,0)` and pass to `realloc`.

### 9. neighCell Zero-Shift Index Wrong (CRITICAL — Root Cause of DiTetracenoHelice NaN)

**Problem**: With nPBC=(1,1,0), PBC shifts array has 9 entries. The zero shift (0,0,0) is at index 4 (center of 3x3 grid), NOT index 0. But `neighCell` was initialized to all zeros, causing every neighbor lookup to apply shift index 0 = (-a-b) ≈ (-49, -28, 0) Å displacement → NaN.

**Solution**: Compute correct zero-shift index: `zero_shift_idx = nPBC[2]*(2*nPBC[1]+1)*(2*nPBC[0]+1) + nPBC[1]*(2*nPBC[0]+1) + nPBC[0]`

### 10. force_node_all Required for Large Aromatics

**Problem**: Using `reorder_nodes_first=True` with `capping_atoms={'H'}` (default) caused NaN for 624-atom DiTetracenoHelice. The `relax_MMFF_pyocl.py` works because it uses `force_node_all` mode.

**Solution**: Use `capping_atoms=set()` and `reorder_nodes_first=False` for large aromatic molecules.

## Solutions Implemented

### Kernel Implementation
- Added `getSurfFlat` kernel to `relax_multi.cl` and `relax_multi_mini.cl`
- Kernel evaluates Hamaker LJ9-3 or Morse potential based on `surf_param.y` mode
- Combines surface and atom REQ parameters via `combineREQ`
- Accumulates forces and energy in `fapos` buffer

### Python Infrastructure
- Added `run_getSurfFlat()` to MolecularDynamics.py
- Fixed buffer allocation, zeroing, and indexing
- Fixed PBC shift generation and neighbor cell indices
- Added mol2 lattice vector parsing

### Test Infrastructure
- Created scan test for single-atom validation
- Created H2O relaxation test showing molecule settling on surface
- Created DiTetracenoHelice test with PBC, demonstrating stable 5000-step relaxation

## Key Takeaways

1. **PBC neighbor indexing is subtle**: The zero-shift index depends on nPBC grid size, not always 0.

2. **GPU buffer allocation matters**: Use correct element sizes (float4 vs float), and zero-initialize.

3. **force_node_all for large aromatics**: Default capping behavior works for small molecules but fails for 600+ atom systems.

4. **Test incrementally**: Start without surface, then add surface, then add PBC to isolate issues.

5. **Float32 precision limits**: Hamaker LJ9-3 involves (R/z)^9 which can exceed float32 precision; relaxed error threshold to 5e-5.

## How to Restart This Work

### Prerequisites
```bash
cd /home/prokop/git/FireCore/tests/tMMFF
export PYTHONPATH=/home/prokop/git/FireCore:$PYTHONPATH
```

### Run Single-Atom Scan Test
```bash
python3 test_flat_surface_gpu.py
```
This validates the GPU kernel against analytical formulas for both Hamaker and Morse.

### Run H2O Relaxation
```bash
python3 test_h2o_surface_relax.py
```
Molecule should drop from z=4.5 to z≈3.25 (near R_eff) over ~2000 steps.

### Run DiTetracenoHelice Relaxation
```bash
python3 test_ditetraceno_surface.py --conf conf1 --mode 1 --nsteps 5000 --stride 100
python3 test_ditetraceno_surface.py --conf conf2 --mode 1 --nsteps 5000 --stride 100
```
Uses nPBC=(1,1,0), reads lattice from mol2 `#lvs`, saves trajectory every 100 steps.

### Key Parameters
- `surf_z`: Surface plane position (default 0.0)
- `surf_REQ`: (R_eff, epsilon, 0, 0) — effective radius and Hamaker constant
- `surf_param`: (K, mode, 0, 0) — K is Morse width, mode=1 for Hamaker, mode=2 for Morse
- `dt`: Timestep (use 0.001 for large molecules, 0.01 for H2O)
- `damp`: Velocity damping factor (0.1-0.15 works well)

### Important Code Locations
- Kernel source: `cpp/common_resources/cl/relax_multi.cl` (search for `getSurfFlat`)
- Python bindings: `pyBall/OCL/MolecularDynamics.py` (search for `run_getSurfFlat`)
- MMFF setup: `pyBall/OCL/MMFF.py` (search for `toMMFFsp3_loc`)
- Buffer allocation: `pyBall/OCL/MolecularDynamics.py` (search for `allocate_cl_buffers`)
- PBC shift generation: `pyBall/OCL/MolecularDynamics.py` (search for `update_pbc_shifts`)

## Trajectory Files Generated

- `trj_h2o_hamaker_surf.xyz` — H2O on Hamaker surface, full trajectory
- `trj_h2o_morse_surf.xyz` — H2O on Morse surface, full trajectory
- `trj_ditetraceno_conf1_hamaker.xyz` — DiTetracenoHelice conf1, every 100 steps
- `trj_ditetraceno_conf2_hamaker.xyz` — DiTetracenoHelice conf2, every 100 steps

## Future Work

1. **Morse surface test**: Currently only Hamaker tested for DiTetracenoHelice. Add Morse mode.

2. **PBC convergence**: Test with larger nPBC (e.g., (2,2,0)) to verify periodic images work correctly.

3. **Force parity**: Compare GPU forces to CPU implementation for the same geometry.

4. **Larger systems**: Test with even larger molecules or multiple molecules in the unit cell.

5. **Performance**: Benchmark GPU performance vs CPU for large molecules.

6. **Surface parameter fitting**: Use this infrastructure to fit surface parameters to experimental data.