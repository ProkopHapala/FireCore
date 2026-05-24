# AFM Simulation: File Analysis and Architecture Design

## **📋 OBJECTIVE**

Analyze existing FireCore files related to visualization, relaxation scanning, and rigid body dynamics to design a modular architecture for AFM simulation. The goal is to separate the solver backend from CLI/GUI frontends and plan two test modes: sequential relaxed scan and parallel AFM image generation.

---

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

---

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

---

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

---

### **Category 4: Testing and Methodology**

#### **FoldedAtomicFunction_rigorous_test_plan.md** (`doc/py/FoldedAtomicFunctions/FoldedAtomicFunction_rigorous_test_plan.md`)
- **Role**: Test plan for folded atomic functions
- **Key Features**:
  - Component-wise outputs (Pauli, London, Coulomb)
  - Parallel 1D scans using OpenCL
  - Rigorous output generation (plots + numeric text)
- **Relevance**: Testing methodology and component separation

---

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

---

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

---

## **⚠️ CRITICAL GAPS AND SOLUTIONS**

### **Gap 1: Anchor Support in GridFF Kernel**

**Problem**: `rigid_body_gridff_kernel` has GridFF but no anchor support.

**Solution**: 
- Copy anchor force calculation from `rigid_body_dynamics_kernel` (lines with `anchor.w > 0.0f`)
- Add anchors parameter to kernel signature
- Update `RigidBodyDynamics.run_gridff()` to pass anchors buffer

**Reference**: `run_rigid_body.py` shows anchor array construction pattern.

---

### **Gap 2: Modular Backend Separation**

**Problem**: RigidBodyDynamics.py mixes low-level buffer management with high-level simulation logic.

**Solution**:
- Create `RigidBodyAFM.py` as high-level wrapper
- Keep `RigidBodyDynamics.py` as low-level GPU interface
- Follow `GridFFRelaxedScan.py` pattern for API design

---

### **Gap 3: Batch Anchor Processing**

**Problem**: Existing batch processing doesn't handle per-replica anchor positions.

**Solution**:
- Extend batch infrastructure to support anchor arrays
- Each replica has its own anchor position
- Extract per-replica anchor forces for AFM image

---

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
