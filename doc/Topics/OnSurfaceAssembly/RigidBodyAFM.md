
See following code maps:
/@AFM Simulation: GPU Rigid Body Dynamics, CPU GridFF Relaxation, and Interactive GUI
https://windsurf.com/codemaps/594f7eaf-c3ab-4139-8f20-d1d2d7f8d401-fe86ab10a43f3d18
/@Rigid Body Dynamics System for AFM Simulation
https://windsurf.com/codemaps/c9f13e1f-edfa-4702-814f-5036d03ea6c9-fe86ab10a43f3d18


# USER

I'm looking what rigid body dynamics implementations on GPU I have, I see these

@Comprehensive_Overview_Molecule_Surface_Systems.md @RigidBodyDynamics.py @Rigid.cl 

and some test scripts
@run_rigid_body.py @test_rigid_gridff_ptcda_batch.py @test_rigid_gridff_surface.py @RigidGridFF_H2O_NaCl_report.md 


please read all the relevant files and tell me how they are related.

What I would like to do in the end is to make AFM simulation program where a molecule (e.g. PTCDA) is attached to a tip, and simulated as rigid body (it has one atom attached by harmonic spring to anchor point represeting metalic tip apex, otherwise it moves in the potential described by GridFF). Figure out what is missing to do it and which from exsiting modules or test scripts is the best starting point.

Notice I wat a very fast simulation. I want to run relaxation of many replicas in parallel and each molecule should be represented by single workgrou (like 16-32 threads)

plese do not program anything yet, just analyse the problem, exsiting code, and write notes

---

# SWE-1.6

I'll read all the relevant files to analyze the existing rigid body dynamics implementations and their relationships.



I'll analyze the existing rigid body dynamics implementations and their relationships to help you understand what's available for your AFM simulation.

## File Relationships

**Core Implementation Chain:**
- `Rigid.cl` (OpenCL kernel) ← [RigidBodyDynamics.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/OCL/RigidBodyDynamics.py:0:0-0:0) (Python wrapper) ← Test scripts
- The kernel is at `cpp/common_resources/cl/Rigid.cl` (gitignored, but referenced by RigidBodyDynamics.py)

**Test Scripts Hierarchy:**
1. [run_rigid_body.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/OCL/run_rigid_body.py:0:0-0:0) - Basic single-body test (no GridFF, just anchors + E-field)
2. [test_rigid_gridff_surface.py](cci:7://file:///home/prokophapala/git/FireCore/tests/tMMFF/test_rigid_gridff_surface.py:0:0-0:0) - Single-body GridFF relaxation with detailed diagnostics
3. [test_rigid_gridff_ptcda_batch.py](cci:7://file:///home/prokophapala/git/FireCore/tests/tMMFF/test_rigid_gridff_ptcda_batch.py:0:0-0:0) - Multi-replica production runner with clustering

## Existing Components Analysis

### RigidBodyDynamics.py (Core Wrapper)

**Architecture:**
- Workgroup size: 32 threads (hardcoded, DEFAULT_WORKGROUP_SIZE)
- Each rigid body = 1 workgroup (exactly what you want)
- Max atoms per body: 128 (configurable)
- Two kernels:
  - `rigid_body_dynamics_kernel`: Basic dynamics with anchors + E-field
  - `rigid_body_gridff_kernel`: GridFF surface forces + body dynamics

**Key Features:**
- **Anchor support**: `anchors` buffer with stiffness in w-component (line 262, 329-338)
  - anchors[:,3] = -1.0 means no anchor
  - anchors[:,3] = stiffness value enables harmonic spring to anchor point
  - This is exactly what you need for tip attachment!
- **GridFF integration**: [init_gridff()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/RigidBodyDynamics.py:345:4-364:80) loads B-spline grid, [run_gridff()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/RigidBodyDynamics.py:382:4-392:27) evaluates surface forces
- **Mass separation**: `mass_trans` (linear) vs `mass_rot` (rotational) for tuning
- **Batch support**: Can simulate multiple replicas in parallel via `n_bodies`
- **State management**: [reset_pose()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/RigidBodyDynamics.py:556:4-575:9) for multi-stage workflows

**What's Missing for AFM:**
- The anchor implementation exists but may need verification for tip-like behavior
- No explicit "tip apex" concept - anchors are just fixed points in space
- Need to check if anchor force is computed correctly in the kernel

### Test Scripts Comparison

| Script | Purpose | GridFF | Anchors | Multi-replica | Best for AFM? |
|--------|---------|--------|---------|---------------|---------------|
| run_rigid_body.py | Basic test | No | Yes (indices + stiffness) | No | Simple anchor demo |
| test_rigid_gridff_surface.py | Single-body GridFF | Yes | No | No | GridFF validation |
| test_rigid_gridff_ptcda_batch.py | Production batch | Yes | No | Yes (1000s) | **Best starting point** |

## What You Need for AFM Simulation

**Requirements vs Existing:**
1. ✅ **Single workgroup per molecule** - Already implemented (32 threads per body)
2. ✅ **GridFF surface forces** - Already implemented via `rigid_body_gridff_kernel`
3. ✅ **Harmonic spring to anchor** - Already implemented via anchors buffer
4. ✅ **Many replicas in parallel** - Already implemented (test_rigid_gridff_ptcda_batch.py does 1000s)
5. ❓ **Tip apex anchor point** - Need to implement (anchors are currently just fixed positions)
6. ❓ **Tip movement control** - Need to add (anchor positions are currently static)

## Recommended Starting Point

**Best base: [test_rigid_gridff_ptcda_batch.py](cci:7://file:///home/prokophapala/git/FireCore/tests/tMMFF/test_rigid_gridff_ptcda_batch.py:0:0-0:0)**

**Reasons:**
- Already has GridFF + multi-replica infrastructure
- Has convergence criteria and stability diagnostics
- Has restart capability (useful for tip scanning)
- Has clustering for analyzing results
- Well-tested production code

**Modifications needed:**
1. Add anchor support to the batch runner (currently only run_rigid_body.py has it)
2. Make anchor positions controllable (for tip scanning)
3. Add tip apex geometry (single atom attachment point)
4. Potentially add tip trajectory control (raster scan paths)

## Implementation Notes

**Anchor Implementation (from run_rigid_body.py lines 106-114):**
```python
def _make_anchor_array(apos, indices, stiffness):
    anchors = np.zeros((apos.shape[0], 4), dtype=np.float32)
    anchors[:, :3] = apos.astype(np.float32)
    anchors[:, 3] = -1.0  # -1 means no anchor
    if stiffness <= 0.0 or not indices:
        return anchors
    anchors[indices, :3] = apos[indices].astype(np.float32)
    anchors[indices, 3] = np.float32(stiffness)  # stiffness value
    return anchors
```

This shows anchors are per-atom, with the anchor point set to the initial atom position. For AFM, you'd want the anchor point to be the tip apex position (controlled externally), not the initial atom position.

**Missing from RigidBodyDynamics.py:**
- The class method [from_xyz_and_grid()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/RigidBodyDynamics.py:491:4-554:18) doesn't accept anchor parameters
- Need to add anchor setup to the batch workflow
- Need to verify the kernel computes anchor forces correctly (can't see Rigid.cl due to gitignore)

## Next Steps

1. **Verify kernel behavior**: Check if Rigid.cl computes anchor forces as harmonic springs to anchor points
2. **Add anchor support to batch runner**: Merge anchor logic from run_rigid_body.py into test_rigid_gridff_ptcda_batch.py
3. **Implement tip control**: Add functions to update anchor positions during simulation (for scanning)
4. **Test with simple case**: Single PTCDA with one anchor atom, verify tip-like behavior

Would you like me to proceed with examining the kernel implementation or start modifying the batch runner to add anchor support?

---

# USER

sorry I have gave you wrong link this is the proper path
@Rigid.cl 

plese read it, and then write comprehenbsive desing document with detailed description of the status, existing entry points, relevant function, relevant arryas and variables, etc. And plan how to implement what we want.

---



# Comprehensive Design Document: AFM Rigid Body Simulation

## 1. Current Status

### Existing Implementation
**File**: [cpp/common_resources/cl/Rigid.cl](cci:7://file:///home/prokophapala/git/FireCore/cpp/common_resources/cl/Rigid.cl:0:0-0:0)

**Two kernels implemented:**
1. **`rigid_body_dynamics_kernel`** (lines 195-310)
   - Basic rigid body dynamics
   - Supports harmonic spring anchors
   - Supports uniform electric field
   - Fixed damping (0.90)
   - NO GridFF surface forces

2. **`rigid_body_gridff_kernel`** (lines 312-443)
   - Rigid body dynamics with GridFF surface forces
   - Configurable damping (lin_damp, ang_damp)
   - Configurable force/torque scaling
   - NO anchor support (critical gap)

### Critical Gap for AFM
The GridFF kernel lacks anchor support, which is essential for tip-molecule attachment. The basic kernel has anchors but no surface forces.

## 2. Entry Points and Function Signatures

### Kernel 1: rigid_body_dynamics_kernel
```c
__kernel void rigid_body_dynamics_kernel(
    __global const int*      mols,           // Molecule atom offsets
    __global       float4*   poss,           // Body positions (xyz, mass in w)
    __global       float4*   qrots,          // Body orientations (quaternions)
    __global       float4*   vposs,          // Linear velocities
    __global       float4*   vrots,          // Angular velocities
    __global const cl_Mat3*  I_body_inv,     // Inverse inertia tensors
    __global const float4*   apos_body,      // Atom positions in body frame
    __global       float4*   apos_world,     // Atom positions in world frame (output)
    __global const float4*   anchors,        // Anchor points (xyz, stiffness in w)
    const int   natoms,                      // Total atoms across all bodies
    const int   niter,                       // Integration steps per launch
    const float dt,                          // Timestep
    const float3  Efield                     // Electric field vector
)
```

### Kernel 2: rigid_body_gridff_kernel
```c
__kernel void rigid_body_gridff_kernel(
    __global const int*      mols,
    __global       float4*   poss,
    __global       float4*   qrots,
    __global       float4*   vposs,
    __global       float4*   vrots,
    __global const cl_Mat3*  I_body_inv,
    __global const float4*   apos_body,
    __global       float4*   apos_world,
    __global const float4*   atom_PLQ,       // Atom PLQ coefficients
    __global const float4*   BsplinePLQ,     // GridFF B-spline data
    __global       float4*   atom_force,     // Per-atom forces (output)
    __global       float4*   body_force,     // Per-body forces (output)
    __global       float4*   body_torque,    // Per-body torques (output)
    const int4               grid_ns,        // Grid dimensions
    const float4             grid_invStep,   // Inverse grid spacing
    const float4             grid_p0,        // Grid origin
    const float              dt,
    const float4             md_params,      // (lin_damp, ang_damp, force_scale, torque_scale)
    const int                niter
)
```

## 3. Key Data Structures

### cl_Mat3 (lines 1-5)
```c
typedef struct __attribute__ ((packed)){
    float4 a; // Row 0
    float4 b; // Row 1
    float4 c; // Row 2
} cl_Mat3;
```
Packed 3x3 matrix stored as three float4 vectors (unused w component).

### Constants (lines 191-193)
```c
#define WORKGROUP_SIZE     32      // Threads per workgroup (1 per body)
#define MAX_ATOMS_PER_BODY 128     // Maximum atoms per rigid body
#define ATOMS_PER_THREAD   4       // Atoms processed per thread
```

## 4. Relevant Functions

### Quaternion Operations
- **`quat_mult`** (line 14): Quaternion multiplication
- **`make_qrot`** (line 53): Create rotation quaternion from angular velocity vector
- **`make_qrot_taylor`** (line 47): Taylor series approximation (faster for small rotations)
- **`quat_to_a/b/c`** (lines 80-82): Extract rotation matrix rows from quaternion

### Matrix Operations
- **`rotate_vec_by_matrix`** (line 70): Rotate vector by 3x3 matrix
- **`rotate_vec_by_matrix_T`** (line 78): Rotate by transpose
- **`mat3_dot`** (line 84): Matrix-vector multiplication
- **`mat3_dot_T`** (line 88): Transpose matrix-vector multiplication

### B-spline GridFF Evaluation
- **`basis`** (line 116): Cubic B-spline basis functions
- **`dbasis`** (line 128): B-spline derivatives
- **`fe1Dcomb`** (line 139): 1D B-spline evaluation
- **`fe2d_comb`** (line 144): 2D B-spline evaluation
- **`fe3d_pbc_comb`** (line 156): 3D B-spline with PBC - core GridFF force evaluator

### PBC Helpers
- **`modulo`** (line 92): Modulo with positive result
- **`make_inds_pbc`** (line 98): Generate PBC index offsets
- **`choose_inds_pbc`** (line 108): Select PBC indices

## 5. Anchor Implementation (Existing in Kernel 1 Only)

**Location**: Lines 250-255 in `rigid_body_dynamics_kernel`

```c
float4 anchor   = anchors[ia0+atom_idx];
if(anchor.w > 0.0f){
    float3 d  = p_world.xyz - anchor.xyz;
    float3 fa = d * -anchor.w;
    f.xyz    += fa;
}
```

**Physics**: Harmonic spring F = -k * (r - r_anchor)
- `anchor.xyz`: Anchor point position in world coordinates
- `anchor.w`: Spring stiffness (k)
- `anchor.w < 0`: Disabled (no anchor)
- Force applied when atom deviates from anchor point

**Critical Issue**: This anchor logic is NOT present in `rigid_body_gridff_kernel`

## 6. Workgroup Parallelization Strategy

**Each workgroup = 1 rigid body** (32 threads)

**Atom distribution**:
- Each thread processes 4 atoms (ATOMS_PER_THREAD)
- Total atoms per body = 32 * 4 = 128 (MAX_ATOMS_PER_BODY)
- Loop: `for (int i=0; i<<ATOMS_PER_THREAD; i++)`

**Force reduction** (lines 262-274):
- Thread-local forces stored in `Lforce[WORKGROUP_SIZE]`
- Tree reduction: stride halving (WORKGROUP_SIZE/4, then full tree)
- Final sum in thread 0

**Torque reduction** (lines 392-397 in GridFF kernel):
- Parallel tree reduction with stride halving
- More efficient than the manual reduction in basic kernel

## 7. Integration with Python Wrapper

**File**: [pyBall/OCL/RigidBodyDynamics.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/OCL/RigidBodyDynamics.py:0:0-0:0)

**Current mapping**:
- Kernel headers defined in `self.kernelheaders` (lines 182-221)
- Buffer allocation in [realloc()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/RigidBodyDynamics.py:241:4-278:31) (lines 242-279)
- GridFF setup in [init_gridff()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/RigidBodyDynamics.py:345:4-364:80) (lines 346-365)
- Basic kernel: [run()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/RigidBodyDynamics.py:366:4-380:27) (lines 367-381)
- GridFF kernel: [run_gridff()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/RigidBodyDynamics.py:382:4-392:27) (lines 383-393)

**Missing**: Anchor support in GridFF kernel signature and argument generation

## 8. Implementation Plan for AFM Simulation

### Phase 1: Add Anchor Support to GridFF Kernel

**Step 1.1: Modify kernel signature**
```c
__kernel void rigid_body_gridff_kernel(
    // ... existing parameters ...
    __global const float4*   anchors,        // ADD THIS
    // ... rest ...
)
```

**Step 1.2: Add anchor force calculation**
Insert after line 378 (after GridFF force calculation):
```c
float4 anchor = anchors[ia];
if(anchor.w > 0.0f){
    float3 d = p_world.xyz - anchor.xyz;
    float3 fa = d * -anchor.w;
    f.xyz += fa;
}
```

**Step 1.3: Update Python wrapper**
- Add `anchors` buffer to [realloc()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/RigidBodyDynamics.py:241:4-278:31) (already exists, line 262)
- Add anchors to `rigid_body_gridff_kernel` header in `kernelheaders`
- Update `generate_kernel_args()` to include anchors for GridFF kernel
- Ensure [init_gridff()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/RigidBodyDynamics.py:345:4-364:80) doesn't overwrite anchors buffer

### Phase 2: Add Anchor Support to Batch Runner

**Step 2.1: Modify [test_rigid_gridff_ptcda_batch.py](cci:7://file:///home/prokophapala/git/FireCore/tests/tMMFF/test_rigid_gridff_ptcda_batch.py:0:0-0:0)**
- Add `--anchor-indices` parameter (atom indices to anchor)
- Add `--anchor-k` parameter (stiffness)
- Add `--anchor-pos` parameter (optional custom anchor positions)
- Implement [_make_anchor_array()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/run_rigid_body.py:105:0-113:18) function (copy from run_rigid_body.py)
- Pass anchors to [RigidBodyDynamics.from_xyz_and_grid()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/RigidBodyDynamics.py:491:4-554:18) (need to extend this method)

**Step 2.2: Extend [RigidBodyDynamics.from_xyz_and_grid()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/RigidBodyDynamics.py:491:4-554:18)**
- Add `anchor_indices=None`, `anchor_k=0.0`, `anchor_positions=None` parameters
- Compute anchor positions if not provided (use initial atom positions)
- Call [upload_state()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/RigidBodyDynamics.py:280:4-343:27) with anchors parameter

### Phase 3: Implement Tip Control

**Step 3.1: Add anchor position update method**
```python
def update_anchors(self, new_anchor_positions):
    # Upload new anchor positions to GPU
    self.toGPU('anchors', new_anchor_positions)
```

**Step 3.2: Implement tip scanning trajectory**
- Add functions to generate raster scan paths
- Add loop to update anchor positions between relaxation steps
- Support different scan patterns (raster, spiral, custom)

### Phase 4: Multi-Replica Tip Scanning

**Step 4.1: Parallel tip positions**
- Each replica can have different tip position
- Support tip position arrays (one per body)
- Enable parallel AFM simulation with different tip heights

**Step 4.2: Analysis tools**
- Track force on anchor (tip-molecule interaction)
- Record energy vs tip position
- Generate force-distance curves

## 9. Recommended Starting Point

**Best approach**: Start with single-body AFM test

1. **Create new test script**: `test_afm_single.py` based on [test_rigid_gridff_surface.py](cci:7://file:///home/prokophapala/git/FireCore/tests/tMMFF/test_rigid_gridff_surface.py:0:0-0:0)
2. **Use modified kernel**: After Phase 1 (anchor support in GridFF)
3. **Simple case**: PTCDA with 1 anchor atom, fixed tip position
4. **Validate**: Verify anchor force is computed correctly
5. **Extend**: Add tip movement, then multi-replica

**Why this approach**:
- Minimal changes initially
- Easy to debug
- Can verify physics before scaling up
- Batch runner can be adapted later

## 10. Open Questions

1. **Anchor point semantics**: Should anchor position be in world coordinates or relative to something?
2. **Tip stiffness**: What's realistic range for AFM tip-molecule spring constant?
3. **Tip geometry**: Single atom vs multiple atoms attached to tip?
4. **Scan speed**: How fast can tip move without causing instability?
5. **Force measurement**: Need to track anchor force separately for AFM signal?


---

# USER

plese investigate also @FoldedAtomicFunction_rigorous_test_plan.md @ExplorerVisPy.py @SurfacePotentialVisualization_discussion.md @ExplorerVisPy_GridFF_Interactive_Plan.md @SequencePlacerVisPy.py @VispyUtils.py 


For better testing and debugging we need to do two things
1) solver module should be separate, in RigidBodyDynamics.py as backend, but the frontend should be ither CLI or GUI, both lightweight with most functionality generall calleg from the shared backedn module
2) we should do two tests
i) relaxed scan, where user drag molecule along some trajectroy, and we store .xyz with final relaxed geometries of the molecule along the scan (optionally we should plot also the substrate to se where are substrate atoms, but we always strugle to alling them properly relative to each other). This is in its nature sequeantial test which will use only single workgroup (single system replica)
ii) then when this first relaxed scan  test (i) give satisfactroy results, we will run parallel test which will relax tip over 2x array of initial geometris and anchor point (simulating different position of the tip as it scan the image), and this will produce AFM image which is basically force on the anchor from the relaxed molecule. Notice that there is always force balance sum of force from the substrate and from the tip is equal to zero at relaxed geometry, so basically force between the molecule and sample is the same as force betwen the tip and sample = force n the achor (up to Force convergence tolerance).  

Perhaps to solve this problem it may be usefully how we did optimization of relaxation trajectroy on CaF2, but that was using CPU and full molecular classical frocefield MMFF/UFF not rigid body.


@run_test_GridFF_CaF2.py 

@RigidGridFF_H2O_NaCl_report.md 

@run_ptcda_caf2_relaxed_scan.py 

@run_tipSpline_scan.py @GridFF_RelaxedScan_cpp_notes.md 


Sorry I a bit fogrot what was each script doijng, plese read it all and make some order from this mess, what is doing what how are each of these files related etc, how they are related or usefull for our problem...


---

# SWE-1.6

# AFM Simulation: File Analysis and Architecture Design

## **📋 OBJECTIVE**

Analyze existing FireCore files related to visualization, relaxation scanning, and rigid body dynamics to design a modular architecture for AFM simulation. The goal is to separate the solver backend from CLI/GUI frontends and plan two test modes: sequential relaxed scan and parallel AFM image generation.

## **🔍 FILE CATEGORIZATION AND MAPPING**

### **Category 1: GPU Rigid Body Dynamics (Core Backend)**

#### **RigidBodyDynamics.py** (`pyBall/OCL/RigidBodyDynamics.py`)
- **Role**: Python wrapper for OpenCL rigid body dynamics
- **Key Features**:
  - GPU buffer management (positions, quaternions, velocities, forces)
  - Kernel loading and execution
  - GridFF initialization and batch processing
  - Anchor support infrastructure (anchors buffer exists in `realloc()`)
- **Critical Gap**: `rigid_body_gridff_kernel` lacks anchor support
- **Relevance**: This is the core solver backend for AFM simulation

#### **Rigid.cl** (`cpp/common_resources/cl/Rigid.cl` or `tests/tUFF/data/cl/Rigid.cl`)
- **Role**: OpenCL kernels for rigid body dynamics
- **Two Kernels**:
  1. `rigid_body_dynamics_kernel`: Has anchor support, no GridFF
  2. `rigid_body_gridff_kernel`: Has GridFF, no anchor support
- **Critical Gap**: Need to merge anchor support into GridFF kernel
- **Workgroup Strategy**: One rigid body per workgroup (32 threads)
- **Max Atoms**: 128 atoms per body

#### **test_rigid_gridff_ptcda_batch.py** (`tests/tMMFF/test_rigid_gridff_ptcda_batch.py`)
- **Role**: Production-level multi-replica rigid GridFF relaxation
- **Key Features**:
  - Batch processing of many replicas in parallel
  - Clustering of distinct minima
  - Stability diagnostics and convergence criteria
- **Relevance**: Template for parallel AFM image generation test

#### **test_rigid_gridff_surface.py** (`tests/tMMFF/test_rigid_gridff_surface.py`)
- **Role**: Single-molecule rigid GridFF relaxation with plotting
- **Key Features**:
  - Z-scans and detailed trajectory plotting
  - Convergence visualization
- **Relevance**: Template for sequential relaxed scan test

#### **run_rigid_body.py** (`pyBall/OCL/run_rigid_body.py`)
- **Role**: Basic single-body test with anchor usage
- **Key Features**:
  - Demonstrates anchor array construction (`_make_anchor_array()`)
  - Shows anchor force calculation in kernel
- **Relevance**: Reference for anchor implementation

### **Category 2: CPU GridFF Relaxation (Reference Backend)**

#### **GridFFRelaxedScan.py** (`pyBall/OCL/GridFFRelaxedScan.py`)
- **Role**: CPU-based GridFF relaxation backend (full molecular forcefield, not rigid body)
- **Key Features**:
  - Harmonic constraint on anchor atom
  - Linear anchor path scanning
  - GridFF visualization (XY/XZ slices)
  - Warm-start vs cold-start relaxation
  - Interactive step relaxation (`step_relaxation()`)
- **Architecture Pattern**:
  - `prepare()`: Initialize molecule, substrate, GridFF, MMFF, MD
  - `set_constraint()`: Set harmonic tether target
  - `relax_to_constraint()`: Run relaxation loop
  - `run_path_scan()`: Scan along trajectory
- **Relevance**: **Key reference for backend separation pattern** - this is exactly the modular design we want for GPU rigid body

#### **run_ptcda_caf2_relaxed_scan.py** (`tests/tMMFF/run_ptcda_caf2_relaxed_scan.py`)
- **Role**: Headless CLI frontend for GridFFRelaxedScan
- **Key Features**:
  - Command-line argument parsing
  - Linear anchor path generation
  - GridFF plotting
  - Output to .xyz and .npz
- **Relevance**: Template for CLI frontend design

#### **run_tipSpline_scan.py** (`tests/tMMFF/run_tipSpline_scan.py`)
- **Role**: C++ MMFF backend with spline-based manipulation
- **Key Features**:
  - Spline constraint paths
  - Simulated annealing optimization
  - Force penalty and bounding box constraints
- **Relevance**: Shows advanced constraint patterns (spline paths, optimization)

#### **run_test_GridFF_CaF2.py** (`tests/tMMFF/run_test_GridFF_CaF2.py`)
- **Role**: GridFF generation and validation
- **Key Features**:
  - Gaussian smearing scans
  - Vertical z-scans with test particles
  - Component-wise potential plotting (Pauli, London, Coulomb)
- **Relevance**: For GridFF data generation and validation

#### **GridFF_RelaxedScan_cpp_notes.md** (`doc/Topics/OnSurfaceAssembly/GridFF_RelaxedScan_cpp_notes.md`)
- **Role**: Design document for C++ high-performance parallel scanning
- **Key Features**:
  - Multi-system framework documentation
  - OpenCL driver interface
  - Constraint system in OpenCL kernels
- **Relevance**: Understanding multi-system GPU architecture

### **Category 3: Visualization and GUI Frontends**

#### **ExplorerVisPy.py** (`pyBall/ExplorerVisPy.py`)
- **Role**: VisPy+PyQt5 GUI for interactive molecule-substrate scanning
- **Key Features**:
  - Multiple scanning backends (Reference, Fast GPU, Folded GPU)
  - Interactive GridFF mode with relaxation
  - Keyboard control for anchor movement (arrow keys, page up/down)
  - Real-time force visualization
- **Architecture Pattern**:
  - `_ensure_relaxed_scanner()`: Initialize backend
  - `_interactive_tick()`: Step relaxation
  - `keyPressEvent()`: User control
- **Relevance**: Template for GUI frontend design

#### **VispyUtils.py** (`pyBall/VispyUtils.py`)
- **Role**: Reusable VisPy widget for atom visualization
- **Key Features**:
  - `AtomScene` class: orthographic top-down view
  - Atom positions, colors, sizes, bonds
  - Picking and dragging
  - Force visualization
- **Relevance**: Shared visualization component for GUI

#### **SequencePlacerVisPy.py** (`pyBall/SequencePlacerVisPy.py`)
- **Role**: GUI for placing sequences on NaCl step edges
- **Key Features**:
  - Interactive visual feedback
  - Headless operation support
- **Relevance**: GUI pattern reference

#### **SurfacePotentialVisualization_discussion.md** (`doc/Topics/OnSurfaceAssembly/SurfacePotentialVisualization_discussion.md`)
- **Role**: Design discussion for surface rendering
- **Key Features**:
  - Iso-surface extraction
  - Test particle evaluation
  - Headless PNG rendering before GUI integration
- **Relevance**: Visualization design principles

#### **ExplorerVisPy_GridFF_Interactive_Plan.md** (`doc/Topics/OnSurfaceAssembly/ExplorerVisPy_GridFF_Interactive_Plan.md`)
- **Role**: Implementation plan for interactive GridFF scanning
- **Key Features**:
  - GridFF infrastructure details
  - Morse-derived PLQ coefficients
  - Coulomb reference subtraction
  - Problems encountered and solutions
- **Relevance**: GridFF implementation details

### **Category 4: Testing and Methodology**

#### **FoldedAtomicFunction_rigorous_test_plan.md** (`doc/py/FoldedAtomicFunctions/FoldedAtomicFunction_rigorous_test_plan.md`)
- **Role**: Test plan for folded atomic functions
- **Key Features**:
  - Component-wise outputs (Pauli, London, Coulomb)
  - Parallel 1D scans using OpenCL
  - Rigorous output generation (plots + numeric text)
- **Relevance**: Testing methodology and component separation

## **🏗️ PROPOSED MODULAR ARCHITECTURE**

### **Backend: RigidBodyAFM.py** (New file in `pyBall/OCL/`)

Inspired by `GridFFRelaxedScan.py` pattern:

```python
class RigidBodyAFM:
    """
    GPU-accelerated rigid body dynamics for AFM simulation.
    Molecule attached to tip via harmonic spring, moving in GridFF potential.
    """
    
    def __init__(self, mol_path, gridff_path, ...):
        # Initialize RigidBodyDynamics backend
        # Load molecule and GridFF
        # Set up anchor parameters
        
    def prepare(self, anchor_z_above, lateral_shift, ...):
        # Initialize GPU buffers
        # Set initial pose
        # Load GridFF
        # Set anchor indices and stiffness
        
    def set_anchor_target(self, target_pos):
        # Update anchor position
        # Upload to GPU
        
    def relax_to_constraint(self, target_pos, niter_max, Fconv):
        # Run rigid body relaxation with GridFF and anchor
        # Return relaxed pose, forces, convergence info
        
    def run_path_scan(self, anchor_path, warm_start=True):
        # Scan along trajectory
        # Save relaxed geometries to .xyz
        # Return scan results
        
    def run_parallel_relaxation(self, anchor_positions):
        # Batch relaxation of many replicas
        # Each replica has different anchor position
        # Return force map (AFM image)
        
    def step_relaxation(self, target_pos, nsub):
        # Single relaxation step for interactive use
        # Return current state
```

**Key Design Principles**:
1. **Backend-only**: No GUI, no CLI, just simulation logic
2. **Reuse RigidBodyDynamics**: Extend existing GPU backend
3. **Follow GridFFRelaxedScan pattern**: Same API structure for consistency
4. **Anchor support**: Add to `rigid_body_gridff_kernel` (critical gap)

---

### **Frontend 1: CLI Script** (New file in `tests/tMMFF/`)

**File**: `run_rigid_afm_scan.py`

Inspired by `run_ptcda_caf2_relaxed_scan.py`:

```python
# Command-line arguments:
# --mol, --sub, --gridff
# --anchor_k, --anchor_z_above
# --x0, --x1, --y0, --z0, --nscan (for linear path)
# --dt, --damp, --Fconv, --nstep_max
# --mode: scan | parallel
# --out_xyz, --out_npz

def main():
    # Parse arguments
    # Initialize RigidBodyAFM backend
    # Prepare system
    if args.mode == 'scan':
        # Generate linear anchor path
        # Run path scan
        # Save .xyz and .npz
    elif args.mode == 'parallel':
        # Generate 2D grid of anchor positions
        # Run parallel relaxation
        # Save force map as AFM image
```

**Key Design Principles**:
1. **Lightweight**: Just argument parsing and backend calls
2. **Two modes**: Sequential scan and parallel AFM image
3. **Output**: .xyz for geometries, .npz for data, PNG for images

---

### **Frontend 2: GUI** (Modify `ExplorerVisPy.py`)

**Add new scan type**: "Rigid AFM"

```python
# In ExplorerVisPy:
class ExplorerVisPy:
    def _ensure_rigid_afm_scanner(self):
        # Initialize RigidBodyAFM backend
        # Set up anchor control
        
    def _interactive_tick_rigid_afm(self):
        # Call backend.step_relaxation()
        # Update visualization
        
    def keyPressEvent(self, event):
        # Add arrow key control for anchor movement
```

**Key Design Principles**:
1. **Reuse existing GUI framework**: Add as new scan type
2. **Interactive control**: Keyboard anchor movement
3. **Real-time visualization**: Show forces and relaxed geometry
4. **Headless option**: Allow batch processing without GUI

---

## **🧪 TEST PLANS**

### **Test 1: Sequential Relaxed Scan**

**Objective**: User drags molecule along trajectory, stores relaxed geometries.

**Script**: `tests/tMMFF/test_rigid_afm_scan.py`

```python
def test_sequential_relaxed_scan():
    # Initialize RigidBodyAFM with PTCDA on CaF2
    # Set anchor stiffness K = 2000 eV/Å²
    # Generate linear path: x from 2 to 12 Å, y = 2 Å, z = 6 Å
    # Run path scan with warm_start=True
    # Save to .xyz with comments: target, anchor_real, opposite, fmax
    # Plot opposite oxygen trajectory
    # Verify convergence: fmax < 1e-4 eV/Å
```

**Expected Output**:
- `ptcda_caf2_rigid_scan.xyz`: 102 frames (101 scan points + initial)
- `ptcda_caf2_rigid_scan.npz`: Full data (positions, forces, anchor paths)
- `ptcda_caf2_rigid_scan_opposite.png`: Trajectory plot

**Validation**:
- Compare with CPU GridFFRelaxedScan results
- Check energy and force convergence
- Verify anchor force calculation

---

### **Test 2: Parallel AFM Image Generation**

**Objective**: Relax many replicas with different tip positions, produce force map.

**Script**: `tests/tMMFF/test_rigid_afm_parallel.py`

```python
def test_parallel_afm_image():
    # Initialize RigidBodyAFM with PTCDA on CaF2
    # Generate 2D grid of anchor positions:
    #   x: 0 to 20 Å, 50 points
    #   y: 0 to 15 Å, 40 points
    #   Total: 2000 replicas
    # Set anchor stiffness K = 2000 eV/Å²
    # Run parallel relaxation (batch processing)
    # Extract force on anchor for each replica
    # Save force map as 2D image (AFM image)
    # Save relaxed geometries for selected minima
```

**Expected Output**:
- `ptcda_caf2_afm_image.png`: Force map (50x40 pixels)
- `ptcda_caf2_afm_data.npz`: Full force and position data
- `ptcda_caf2_minima.xyz`: Relaxed geometries of distinct minima

**Validation**:
- Verify force map shows expected features
- Check convergence across all replicas
- Compare with sequential scan at selected points

## **🔧 IMPLEMENTATION ROADMAP**

### **Phase 1: Backend Extension (Critical)**

1. **Add anchor support to `rigid_body_gridff_kernel`**:
   - Copy anchor force calculation from `rigid_body_dynamics_kernel`
   - Add anchors parameter to kernel signature
   - Test with single molecule

2. **Create `RigidBodyAFM.py`**:
   - Extend `RigidBodyDynamics` class
   - Implement `prepare()`, `set_anchor_target()`, `relax_to_constraint()`
   - Follow `GridFFRelaxedScan` API pattern
   - Test with single relaxation step

3. **Add path scanning to backend**:
   - Implement `run_path_scan()`
   - Support warm_start and cold_start
   - Save to .xyz and .npz

### **Phase 2: Sequential Test**

4. **Create `test_rigid_afm_scan.py`**:
   - Implement sequential relaxed scan
   - Compare with CPU GridFFRelaxedScan
   - Validate convergence and forces

### **Phase 3: Parallel Test**

5. **Add batch processing to backend**:
   - Implement `run_parallel_relaxation()`
   - Use existing batch infrastructure from `test_rigid_gridff_ptcda_batch.py`
   - Extract anchor forces

6. **Create `test_rigid_afm_parallel.py`**:
   - Implement 2D anchor grid generation
   - Run parallel relaxation
   - Generate AFM force map image

### **Phase 4: Frontends**

7. **Create CLI script `run_rigid_afm_scan.py`**:
   - Add argument parsing
   - Support both scan and parallel modes
   - Output formatting

8. **Add GUI support to `ExplorerVisPy.py`**:
   - Add "Rigid AFM" scan type
   - Implement interactive anchor control
   - Real-time visualization

---

## **📊 FILE RELATIONSHIP DIAGRAM**

```
┌─────────────────────────────────────────────────────────────────┐
│                    BACKEND (Solver Layer)                        │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  RigidBodyDynamics.py  ──►  Rigid.cl (kernels)                 │
│       │                              │                          │
│       │                              ├── rigid_body_dynamics_kernel
│       │                              │    (has anchor, no GridFF)
│       │                              │                          │
│       │                              └── rigid_body_gridff_kernel
│       │                                   (has GridFF, NO anchor)
│       │                                   ⚠️ CRITICAL GAP         │
│       │                                                          │
│  ┌─────▼─────┐  ┌──────────────────────────────────────────────┐ │
│  │  NEW:     │  │  Reference: GridFFRelaxedScan.py              │ │
│  │RigidBodyAFM│  │  (CPU backend, same API pattern)            │ │
│  │  .py      │  └──────────────────────────────────────────────┘ │
│  └─────┬─────┘                                                    │
│         │                                                         │
│         │  Methods to implement:                                 │
│         │  - prepare()                                           │
│         │  - set_anchor_target()                                 │
│         │  - relax_to_constraint()                               │
│         │  - run_path_scan()                                     │
│         │  - run_parallel_relaxation()                           │
│         │                                                         │
└─────────┼─────────────────────────────────────────────────────────┘
          │
          ▼
┌─────────────────────────────────────────────────────────────────┐
│                    FRONTEND (Interface Layer)                    │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  ┌─────────────────┐  ┌──────────────────────────────────────┐  │
│  │  NEW: CLI       │  │  Modify: ExplorerVisPy.py          │  │
│  │  run_rigid_     │  │  (Add "Rigid AFM" scan type)       │  │
│  │  afm_scan.py    │  │                                    │  │
│  └────────┬────────┘  └──────────────┬─────────────────────┘  │
│           │                          │                         │
│           │                          │                         │
│           ▼                          ▼                         │
│  ┌─────────────────┐  ┌──────────────────────────────────────┐  │
│  │  Tests:         │  │  Visualization: VispyUtils.py       │  │
│  │  - test_rigid_  │  │  (AtomScene widget)                 │  │
│  │    afm_scan.py  │  │                                    │  │
│  │  - test_rigid_  │  │                                    │  │
│  │    afm_parallel │  │                                    │  │
│  │    .py          │  │                                    │  │
│  └─────────────────┘  └──────────────────────────────────────┘  │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

---

## **🎯 KEY INSIGHTS FROM EXISTING CODE**

### **From GridFFRelaxedScan.py (CPU Backend Pattern)**

1. **Separation of Concerns**: Backend handles simulation, frontend handles I/O
2. **Constraint System**: Harmonic tether on specific atom
3. **Path Scanning**: Linear trajectory with warm/cold start options
4. **Interactive Mode**: `step_relaxation()` for GUI integration
5. **Visualization**: GridFF slicing and plotting

### **From ExplorerVisPy.py (GUI Pattern)**

1. **Backend Initialization**: `_ensure_*_scanner()` lazy initialization
2. **Interactive Loop**: `_interactive_tick()` per-frame updates
3. **User Control**: Keyboard events for parameter adjustment
4. **Real-time Feedback**: Force and position display

### **From test_rigid_gridff_ptcda_batch.py (Parallel Pattern)**

1. **Batch Processing**: Many replicas in parallel
2. **Clustering**: Identify distinct minima
3. **Stability Diagnostics**: Convergence checking
4. **Data Management**: Large-scale result handling

### **From run_tipSpline_scan.py (Advanced Constraints)**

1. **Spline Paths**: Flexible trajectory definition
2. **Optimization**: Simulated annealing for path refinement
3. **Force Penalties**: Constrain maximum forces
4. **Bounding Box**: Keep molecule in reliable region

## **⚠️ CRITICAL GAPS AND SOLUTIONS**

### **Gap 1: Anchor Support in GridFF Kernel**

**Problem**: `rigid_body_gridff_kernel` has GridFF but no anchor support.

**Solution**: 
- Copy anchor force calculation from `rigid_body_dynamics_kernel` (lines with `anchor.w > 0.0f`)
- Add anchors parameter to kernel signature
- Update `RigidBodyDynamics.run_gridff()` to pass anchors buffer

**Reference**: `run_rigid_body.py` shows anchor array construction pattern.

### **Gap 2: Modular Backend Separation**

**Problem**: RigidBodyDynamics.py mixes low-level buffer management with high-level simulation logic.

**Solution**:
- Create `RigidBodyAFM.py` as high-level wrapper
- Keep `RigidBodyDynamics.py` as low-level GPU interface
- Follow `GridFFRelaxedScan.py` pattern for API design

### **Gap 3: Batch Anchor Processing**

**Problem**: Existing batch processing doesn't handle per-replica anchor positions.

**Solution**:
- Extend batch infrastructure to support anchor arrays
- Each replica has its own anchor position
- Extract per-replica anchor forces for AFM image

## **📝 SUMMARY**

### **File Categories**:
1. **GPU Rigid Body Dynamics**: RigidBodyDynamics.py, Rigid.cl, test scripts
2. **CPU GridFF Relaxation**: GridFFRelaxedScan.py, run_ptcda_caf2_relaxed_scan.py, run_tipSpline_scan.py
3. **Visualization/GUI**: ExplorerVisPy.py, VispyUtils.py, SequencePlacerVisPy.py
4. **Testing/Methodology**: FoldedAtomicFunction test plan, GridFF generation scripts

### **Architecture**:
- **Backend**: New `RigidBodyAFM.py` extending `RigidBodyDynamics.py`
- **CLI**: New `run_rigid_afm_scan.py` for batch processing
- **GUI**: Modify `ExplorerVisPy.py` to add "Rigid AFM" scan type

### **Tests**:
1. **Sequential**: Linear anchor path scan, save relaxed geometries
2. **Parallel**: 2D anchor grid, generate AFM force map image

### **Critical Path**:
1. Add anchor support to `rigid_body_gridff_kernel`
2. Create `RigidBodyAFM.py` with GridFFRelaxedScan API pattern
3. Implement sequential test
4. Implement parallel test
5. Add CLI and GUI frontends

This architecture maximizes code reuse, follows existing patterns, and provides a clear path from backend to frontends.

## Progress Notes (2026-04-26)

### Achievements
1.  **Parallel Grid Sampler**: Implemented `sample_gridff_single_atom` for ultra-fast grid verification. It evaluates forces/energies for 1D/2D scans in parallel on GPU.
2.  **Diagnostic Tools**: Added `plot_gridff_diagnostics` to `RigidBodyAFM.py` to visualize XY/XZ potential slices with substrate atom overlays.
3.  **Extended Scans**: Simulation script `run_rigid_afm_scan.py` now supports long scans and automated projection plotting (COM and corner atoms).

### Identified Issues with CaF2 Grid
1.  **Non-Periodicity**: `CaF2_6L_Ni3_rect` grid has a large discontinuity at the X/Y boundaries (~12.6 in Pauli, ~33.2 in London). The last voxel is effectively empty/zero, creating a "strip" artifact.
2.  **Alignment Shift**: The GridFF was generated centered in XY ($p_0 = -0.5 \times L_s$), while the substrate atoms are in the original coordinate system. This causes a half-box shift if not accounted for.
3.  **Electrostatic Sign**: Above the F-terminated surface ($Z > 4.87\text{\AA}$), the potential is negative everywhere, repelling $Q=-0.2$ atoms. This is physically consistent for an anion-terminated surface but explains the "inverse effect" observed.
4.  **Isotropy**: Suspected anisotropic scaling if lattice vectors and voxel counts are mismatched.

### Next Steps
- Switch to **NaCl (100)** surface to verify if these artifacts are specific to the CaF2 grid generation.
- Correct $p_0$ handling in `RigidBodyDynamics.from_xyz_and_grid` once generation parameters are standardized.

### Final Implementation (2026-04-26)

1.  **Consolidated CLI**: `run_rigid_afm_scan.py` now integrates diagnostics, single-atom tests, and full molecule simulations (both 1D trajectory and 2D parallel images).
2.  **Parallel Image Simulation**: The `--scan2d` mode relaxes thousands of molecules simultaneously on the GPU (one per pixel), enabling rapid generation of high-resolution AFM images.
3.  **Host-Side Sync**: `RigidBodyDynamics.py` now correctly manages `self.anchors` on the host, ensuring that multiple tip updates in a loop are properly uploaded to the GPU.
4.  **Lattice Autodetect**: All tools now automatically detect lattice vectors from substrate `.xyz` files to ensure proper scaling, periodicity, and aspect ratios in plots.

### **How to Use**

All simulation and diagnostic tasks are now consolidated into a single script:
`tests/tMMFF/run_rigid_afm_scan.py`

#### **1. Full All-in-One Run**
Run diagnostics, single-atom alignment tests, 1D trajectory, and 2D parallel images in one command:
```bash
python3 run_rigid_afm_scan.py --sub NaCl_1x1_L3 --mol xyz/PTCDA.xyz --all --outdir results_nacl
```

#### **2. Standalone Tasks**
- **Grid Diagnostics**: `python3 run_rigid_afm_scan.py --sub NaCl_1x1_L3 --diag`
- **One-Atom Alignment Test**: `python3 run_rigid_afm_scan.py --sub NaCl_1x1_L3 --test_grid --test_z 5.0`
- **1D Relaxed Trajectory**: `python3 run_rigid_afm_scan.py --sub NaCl_1x1_L3 --scan1d`
- **2D Parallel AFM Image**: `python3 run_rigid_afm_scan.py --sub NaCl_1x1_L3 --scan2d`

### **Results Summary (NaCl vs CaF2)**
- **NaCl (100)**: Verified perfect periodicity ($4\text{\AA}$ unit cell) and correct alignment between GridFF wells and Na/Cl ion positions.
- **CaF2**: Identified significant non-periodicity and XY centering artifacts in the current `Bspline_PLQd.npy` grid. Future work should regenerate this grid using the standardized centered origin ($p_0 = -0.5 \times L_s$) and ensure periodic padding.

### **Advanced Visualization Features (2026-04-26)**

The 2D parallel AFM scan now generates comprehensive 11-panel visualizations in `scan_molecule_2d_full.png`:

**Panel Layout**:
- **Row 1**: Force Z, Energy, Convergence map
- **Row 2**: COM Z shift, COM XY shift (HSV), COM X shift, COM Y shift
- **Row 3**: Rotation angle (from quaternion), Tilt axis X, Tilt axis Y, Tilt axis Z

**COM XY Shift (HSV Plot)**:
- **Hue**: Direction of XY shift (angle in XY plane)
- **Brightness**: Magnitude of shift
- **Purpose**: Visualize both direction and magnitude of lateral displacement in a single plot

**Quaternion Tilt Visualization**:
- **Rotation Angle**: $2 \times \arccos(|q_w|) \times 180/\pi$ degrees
- **Tilt Axis**: Normalized rotation axis from quaternion $(q_x, q_y, q_z)$ components
- **Purpose**: Characterize molecular reorientation during relaxation

**Diagnostic Output**:
- Anchor displacement after relaxation (mean and max in Å)
- Convergence rate (number of converged replicas / total)

### **2D Scan Parameters**

Customizable scan area and resolution:
```bash
--nx 160 --ny 160                    # Number of pixels (160×160 = 25,600 replicas)
--x0_2d 0.0 --x1_2d 16.0            # X range in Å
--y0_2d 0.0 --y1_2d 16.0            # Y range in Å
```

**Example**: Scan 16×16 Å area with 0.1 Å resolution:
```bash
python3 run_rigid_afm_scan.py --scan2d --nx 160 --ny 160 \
    --x0_2d 0.0 --x1_2d 16.0 --y0_2d 0.0 --y1_2d 16.0 \
    --sub NaCl_1x1_L2
```

**Note**: Large scans (e.g., 160×160 = 25,600 replicas) require significant GPU memory and computation time. Start with smaller tests (e.g., 40×40) to verify convergence.

### **Debug Mode**

Enable diagnostic plots with 2D scan using `--debug` flag:
```bash
python3 run_rigid_afm_scan.py --scan2d --debug \
    --nx 160 --ny 160 --x0_2d 0.0 --x1_2d 16.0 \
    --y0_2d 0.0 --y1_2d 16.0 --sub NaCl_1x1_L2
```

The `--debug` flag enables:
- GridFF diagnostics (`diag_raw_grids.png`)
- Single-atom 2D alignment test (`test_single_atom_2d.png`)
- 2D relaxed AFM scan with full visualization

### **Convergence Optimization**

Current convergence rate is low (238/25600 = 0.9% for default parameters). To improve:

1. **Relax convergence criteria**:
   ```bash
   --fconv 1e-2 --tconv 1e-2    # Less strict (default: 1e-3)
   ```

2. **Increase relaxation steps**:
   ```bash
   --nsteps 2000                 # More steps (default: 1000)
   ```

3. **Adjust anchor stiffness**:
   ```bash
   --kanchor 50.0                # Softer spring (default: 20.0)
   ```

**Diagnostic check**: The script outputs anchor displacement after relaxation. Mean displacement should be small (<0.01 Å) if the anchor constraint is working correctly.