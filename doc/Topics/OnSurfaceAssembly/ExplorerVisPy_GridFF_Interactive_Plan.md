# Interactive GridFF Scanning Implementation Plan


## **🎯 OBJECTIVE**

Implement high-performance parallel relaxed scanning using C++ OpenCL framework with spline-based constraint paths for 1000+ parallel system replicas.

---

## **🔍 C++ SYSTEM ARCHITECTURE ANALYSIS**

### **📋 Core Multi-System Framework**

#### **MolWorld_sp3_multi.h - Multi-Replica Container**
- **Primary Role**: Manages `nSystems` parallel molecular replicas
- **Key Architecture**:
  ```cpp
  int nSystems;                    // Number of parallel system replicas
  OCL_MM ocl;                      // OpenCL driver for multi-system execution
  std::vector<MMFFsp3_loc> ffls;   // Per-system force field instances
  FIRE fire[MAX_SYSTEMS];         // Per-system optimization engines
  ```
- **Memory Layout**: All systems stored in contiguous GPU buffers
  - `atoms[nSystems * nvecs]` - positions and velocities
  - `constr[nSystems * natoms]` - constraint targets and stiffness
  - `REQs[nSystems * natoms]` - atomic parameters

#### **OCL_MM.h - High-Performance OpenCL Driver**
- **Key Advantage**: Native C++ OpenCL (no Python overhead)
- **Multi-System Support**:
  ```cpp
  int initAtomsForces(int nSystems_, int nAtoms_, int nnode_, int npbc_);
  // Allocates: atoms[nSystems*nvecs], constr[nSystems*natoms], etc.
  ```
- **GPU Kernel Execution**: 2D work groups `(natoms, nSystems)`
- **Constraint System**: Already implemented in `relax_multi.cl`

### **⚡ EXISTING CONSTRAINT SYSTEM**

#### **OpenCL Kernel: relax_multi.cl - updateAtomsMMFFf4()**
```cpp
// Lines 1207-1216: Constraint force calculation
float4 cons = constr[iaa];  // (x,y,z,K) - constraint target and stiffness
if(cons.w>0.f){             // Active constraint if stiffness > 0
    float4 cK = constrK[iaa];
    const float3 fc = (cons.xyz - pe.xyz)*cK.xyz;
    fe.xyz += fc;           // Add harmonic constraint force
}
```

#### **Constraint Buffer Structure**
- **constr[iaa]**: `{x_target, y_target, z_target, stiffness_K}`
- **constrK[iaa]**: `{kx, ky, kz, ?}` - individual stiffness components
- **Per-System Indexing**: `iaa = isys * natoms + iatom`

#### **Multi-System Parallel Execution**
```cpp
const int iS = get_global_id(1);  // System index (0 to nSystems-1)
const int iG = get_global_id(0);  // Atom index (0 to natoms-1)
const int iaa = iG + iS*natoms;   // Global atom index
```

---

## **🎯 PROPOSED SPLINE-BASED CONSTRAINT SYSTEM**

### **📍 Concept: Spline-Controlled Constraint Paths**

#### **Core Innovation**
- **Spline Definition**: Control points define smooth 3D paths
- **Per-Replica Control**: Each system replica follows different spline parameters
- **Parallel Execution**: 1000+ systems run simultaneously with different constraint trajectories

#### **Mathematical Framework**
```cpp
// Spline types available
#include "Bspline.h"           // Cubic B-spline interpolation
#include "spline_hermite.h"    // Hermite cubic spline
```

#### **Implementation Architecture**
```cpp
struct SplineConstraint {
    Vec3d* control_points;     // Control points for spline
    int    n_control_points;   // Number of control points
    double t_parameter;        // Current parameter along spline (0 to 1)
    int    spline_type;        // BSPLINE or HERMITE
    double dt_step;            // Parameter increment per relaxation step
};
```

### **🔄 Multi-Replica Spline System**

#### **Per-Replica Spline Control**
```cpp
class SplineConstraintManager {
    int nSystems;
    std::vector<SplineConstraint> splines;  // [nSystems]
    
    void updateConstraints(int isys, double t) {
        // Evaluate spline at parameter t
        Vec3d pos = evaluateSpline(splines[isys], t);
        
        // Update GPU constraint buffer
        int iaa = isys * natoms + anchor_atom;
        constr[iaa] = {pos.x, pos.y, pos.z, stiffness_K};
    }
};
```

#### **Parallel Relaxation Loop**
```cpp
for(int istep=0; istep<nsteps; istep++){
    // Update all spline constraints
    for(int isys=0; isys<nSystems; isys++){
        double t = istep * dt_step;  // Different t per system
        updateConstraints(isys, t);
    }
    
    // Upload updated constraints to GPU
    uploadConstraints();
    
    // Run parallel relaxation step
    ocl.run_updateAtomsMMFFf4();  // All nSystems simultaneously
    
    // Check convergence for each system
    checkConvergence();
}
```

---

## **🚀 PERFORMANCE ADVANTAGES**

### **✅ Massive Parallelization**
- **GPU Utilization**: 1000+ systems exploit full GPU capability
- **Native C++ Driver**: 10-100x faster than pyOpenCL
- **Contiguous Memory**: Optimal GPU memory access patterns

### **✅ Flexible Trajectory Design**
- **Complex Paths**: B-splines and Hermite splines enable sophisticated trajectories
- **Per-Replica Variation**: Each system can explore different parameter space
- **Smooth Interpolation**: Continuous constraint motion avoids discontinuities

### **✅ Scientific Applications**
- **Potential Energy Surface Mapping**: Systematically explore molecular configurations
- **Reaction Path Sampling**: Multiple parallel reaction coordinates
- **Optimization Landscape**: Global optimization with diverse starting points

---

## **📋 IMPLEMENTATION PLAN**

### **Phase 1: Spline Infrastructure**
1. **Spline Evaluation Classes**
   - Implement B-spline evaluation from `Bspline.h`
   - Implement Hermite spline from `spline_hermite.h`
   - Create unified `SplineConstraint` interface

2. **Integration with MolWorld_sp3_multi**
   - Add spline constraint manager to multi-system framework
   - Extend constraint buffer updates for spline evaluation
   - Implement per-replica spline parameter storage

### **Phase 2: Parallel Execution**
1. **Constraint Update System**
   - Efficient GPU constraint buffer updates
   - Batch updates to minimize CPU-GPU transfers
   - Parameter synchronization across systems

2. **Relaxation Loop Optimization**
   - Parallel convergence checking
   - Adaptive step sizes per system
   - Early termination for converged systems

### **Phase 3: Interface Development**
1. **Command-Line Interface**
   - Spline definition files (control points)
   - Multi-system configuration parameters
   - Output management for parallel results

2. **Python Integration**
   - Extend `MMFFmulti_lib.cpp` for spline access
   - Interface to `MMFF_multi.py` for workflow integration
   - Visualization tools for parallel trajectories

---

## **🔬 TECHNICAL SPECIFICATIONS**

### **Performance Targets**
- **System Count**: 1000+ parallel replicas
- **Speedup**: 10-100x over Python implementation
- **Memory**: Efficient GPU buffer utilization
- **Scalability**: Linear performance scaling to GPU limits

### **Constraint System Requirements**
- **Spline Types**: B-spline, Hermite cubic spline
- **Control Points**: User-defined 3D trajectories
- **Per-Replica Parameters**: Independent spline parameters per system
- **Update Frequency**: Configurable constraint update intervals

### **Integration Points**
- **Existing Kernels**: `relax_multi.cl` constraint system already implemented
- **Memory Layout**: Compatible with current multi-system buffer structure
- **OpenCL Driver**: `OCL_MM.h` provides native C++ performance
- **Force Field**: `MMFFsp3_loc` handles molecular topology

---

## **📊 EXPECTED OUTCOMES**

### **Scientific Impact**
- **High-Throughput Screening**: 1000x parallel molecular configuration exploration
- **Complex Trajectory Mapping**: Sophisticated constraint paths beyond linear scans
- **Energy Landscape Analysis**: Comprehensive potential surface mapping

### **Performance Impact**
- **Order of Magnitude Speedup**: Native C++ + massive parallelization
- **GPU Utilization**: Full exploitation of modern GPU capabilities
- **Scalability**: Linear scaling to available GPU resources

### **Technical Advancement**
- **Framework Extension**: Leverages existing robust C++ molecular dynamics
- **Minimal Disruption**: Uses established constraint system and OpenCL infrastructure
- **Future-Proof**: Extensible to other force fields and optimization methods

---

## **Goals**

### Primary Objectives
1. **Interactive PTCDA-on-CaF2 Scanning**: Real-time manipulation of PTCDA molecule on CaF2 substrate using GridFF potentials
2. **Harmonic Constraint Control**: Drag corner oxygen atom with keyboard arrows while molecule relaxes on surface
3. **GridFF Force Evaluation**: Use precomputed B-spline substrate potentials (Pauli/London/Coulomb) instead of on-the-fly XYZ evaluation
4. **Real-time Visualization**: Update atom positions and energy display during interactive relaxation
5. **Debug XZ Slice Visualization**: Display cross-sections of GridFF potentials with adjustable PLQ coefficients

### Secondary Objectives
- Batch scanning workflow for offline XYZ generation
- Energy landscape analysis and path optimization
- Integration with existing ExplorerVisPy framework

## **Existing Functionality Analysis**

### **Core Module Architecture**

#### **@[pyBall/OCL/MMparams.py] - Parameter Database**
**Responsibility**: Load and manage atomic force field parameters
- **Classes**: `ElementType`, `AtomType`
- **Key Functions**: `read_atom_types()`, `read_element_types()`
- **Data Sources**: `cpp/common_resources/AtomTypes.dat`, `ElementTypes.dat`
- **Connection**: Provides REQ parameters to MMFF topology builder
- **Integration Point**: Used by MMFF class for force field parameter assignment

#### **@[pyBall/OCL/MMFF.py] - Force Field Topology Builder**
**Responsibility**: Convert molecular structure to MMFF representation
- **Key Functions**:
  - `toMMFFsp3_loc()`: Convert AtomicSystem to MMFF data structures
  - `build_topology()`: Generate bonds, angles, pi-orbitals from connectivity
  - `realloc()`: Allocate arrays for atoms, nodes, constraints
- **Data Structures**:
  - `apos`: Atomic positions (natoms + nnode + pi-orbitals)
  - `neighs`: Bond connectivity and neighbor lists
  - `REQs`: Nonbonded parameters (RvdW, EvdW, charge)
  - `bKs`, `bLs`: Bond force constants and equilibrium lengths
- **Connection**: Packs topology data into MolecularDynamics GPU buffers
- **Integration Point**: Creates force field data for GridFF scanner

#### **@[pyBall/OCL/MolecularDynamics.py] - GPU MD Engine**
**Responsibility**: PyOpenCL wrapper for OpenCL MD kernels
- **Key Functions**:
  - `pack_system()`: Upload MMFF topology to GPU buffers
  - `initGridFF()`: Load B-spline GridFF data to GPU texture/buffers
  - `run_updateAtomsMMFFf4()`: Velocity Verlet integration with constraints
  - `run_getNonBond_GridFF_Bspline()`: GridFF B-spline force evaluation
  - `run_cleanForceMMFFf4()`: Zero force buffers
  - `download_results()`: Retrieve positions and forces from GPU
- **Constraint System**: 
  - Buffers: `constr` (target positions), `constrK` (spring constants)
  - Kernel Support: `updateAtomsMMFFf4` already handles constraints
- **Connection**: Used by ExplorerVisPy's `fast_scanner` for XYZ-based scanning
- **Integration Point**: Will be extended for GridFF interactive mode

#### **@[cpp/common_resources/cl/relax_multi.cl] - OpenCL Kernels**
**Responsibility**: GPU implementation of force evaluation and integration
- **Key Kernels**:
  - `updateAtomsMMFFf4()`: Integration with harmonic constraints (lines 1097-1098)
  - `getNonBond_GridFF_Bspline()`: Tricubic B-spline GridFF interpolation
  - `getMMFFf4()`: Bonded force evaluation (bonds, angles, pi-orbitals)
  - `cleanForceMMFFf4()`: Initialize force buffers to zero
- **Constraint Implementation**: 
  - Format: `constr[i] = (x_target, y_target, z_target, K_flag)`
  - Format: `constrK[i] = (kx, ky, kz, padding)`
  - Force: `F_constraint = -K * (r_current - r_target)`
- **Connection**: Called by MolecularDynamics Python wrapper
- **Integration Point**: No modifications needed - already supports required functionality

### **Current ExplorerVisPy Framework**

#### **Existing Scanner Architecture**
- **Reference Scanner**: `InteractionScanner` - CPU-based evaluation
- **Fast GPU Scanner**: `MolecularDynamics` - XYZ-based surface forces
- **Folded GPU Scanner**: `MolecularDynamics` - Experimental folded basis
- **Missing**: GridFF-specific scanner instance

#### **Current Relaxation Support**
- **UI Controls**: 
  - `cb_relax`: Enable constrained relaxation
  - `sp_springk`, `sp_rdt`, `sp_rsteps`: Relaxation parameters
- **Scanner Properties**: `relax_dt`, `relax_nsteps`, `spring_k`
- **Backend Support**: `_run_scan_*()` functions support `relax=True`
- **Missing**: Real-time interactive relaxation loop

#### **Current Visualization Framework**
- **3D Canvas**: VisPy scene with molecule, substrate, surface mesh
- **Energy Display**: Real-time energy component display
- **Scan Visualization**: Trajectory points and paths
- **Missing**: XZ slice debugging visualization

## **New Functionality Implementation**

### **Phase 1: GridFF Backend Integration**

#### **GridFF Scanner Instance**
- **Location**: `ExplorerVisPy.__init__()`
- **Purpose**: Dedicated MolecularDynamics instance for GridFF
- **Connection**: Parallel to existing `fast_scanner` and `folded_scanner`
- **Data Flow**: Load GridFF → Build MMFF → Pack to GPU → Setup constraints

#### **GridFF Data Loading**
- **Location**: `ExplorerVisPy._init_gridff_scanner()`
- **Purpose**: Load CaF2 B-spline data and initialize GPU
- **Data Sources**: 
  - GridFF: `cpp/common_resources/CaF2_6L_Ni3_rect_nx2_nz1_L2_top/Bspline_PLQd.npy`
  - Molecule: User-sified XYZ file
  - Parameters: `AtomTypes.dat`, `ElementTypes.dat`
- **Connection**: Uses `MolecularDynamics.initGridFF()` and `MMFF.build_topology()`

#### **Constraint System Setup**
- **Location**: `ExplorerVisPy._setup_gridff_constraints()`
- **Purpose**: Configure harmonic constraint on corner oxygen atom
- **Target Atom**: Find oxygen with minimum (x+y) coordinates (corner)
- **Parameters**: Spring constant K=50.0 eV/Å² (adjustable via UI)
- **Connection**: Uploads to GPU `constr`/`constrK` buffers

### **Phase 2: Interactive Control System**

#### **Interactive Mode Toggle**
- **Location**: `ExplorerVisPy._on_gridff_interactive()`
- **Purpose**: Start/stop real-time relaxation loop
- **UI Control**: New button "Start/Stop GridFF Interactive"
- **State Management**: `gridff_active` flag, timer management
- **Connection**: Controls relaxation timer and keyboard input

#### **Real-time Relaxation Loop**
- **Location**: `ExplorerVisPy._gridff_relax_step()`
- **Purpose**: Single MD step with GridFF forces and constraints
- **Force Evaluation Sequence**:
  1. `run_cleanForceMMFFf4()` - Zero forces
  2. `run_getMMFFf4()` - Bonded MMFF forces
  3. `run_getNonBond_GridFF_Bspline()` - GridFF surface forces
  4. `run_updateAtomsMMFFf4()` - Integration with constraints
- **Update Rate**: 50ms timer (20 FPS) for smooth interaction
- **Connection**: Updates visualization and energy display

#### **Keyboard Control System**
- **Location**: `ExplorerVisPy.keyPressEvent()`
- **Purpose**: Move constraint target with arrow keys
- **Controls**:
  - Arrow keys: Move anchor in XY plane (±0.2 Å per press)
  - PageUp/PageDown: Move anchor in Z direction
  - Active only during `gridff_active` mode
- **Connection**: Updates GPU constraint buffers in real-time

### **Phase 3: Visualization Enhancements**

#### **Real-time Position Updates**
- **Location**: `ExplorerVisPy._update_gridff_visuals()`
- **Purpose**: Download relaxed positions and update 3D view
- **Data Flow**: GPU → `download_results()` → `mol_apos` → VisPy markers
- **Energy Display**: Evaluate current GridFF energy components
- **Connection**: Integrates with existing `_update_mol_visual()` framework

#### **XZ Slice Debug Visualization**
- **Location**: New UI group "XZ Slice Debug"
- **Purpose**: Display cross-section of GridFF potentials
- **Controls**:
  - PLQ coefficient sliders (Pauli, London, Coulomb)
  - Y-position selector for slice location
  - Plot button to generate matplotlib slice
- **Data Processing**: Extract XZ slice + combine PLQ components
- **Connection**: Uses existing matplotlib panel in ExplorerVisPy

#### **Constraint Visualization**
- **Location**: Enhanced 3D scene rendering
- **Purpose**: Show constrained atom and target position
- **Implementation**: Different color/marker for anchor atom
- **Connection**: Extends existing VisPy marker system

### **Phase 4: Batch Scanning Workflow**

#### **Offline Batch Script**
- **Location**: `tests/tMMFF/ptcda_caf2_gridff_scan.py`
- **Purpose**: Generate XYZ configurations for relaxed scan paths
- **Workflow**:
  1. Load PTCDA + CaF2 GridFF
  2. Define scan grid over surface
  3. For each point: Set constraint → Relax → Save XYZ
  4. Generate trajectory file and energy log
- **Connection**: Reuses GridFF scanner from interactive mode

#### **Path Optimization**
- **Location**: Extension of batch script
- **Purpose**: Find minimum energy paths using constrained relaxation
- **Methods**: Gradient descent, simulated annealing, genetic algorithms
- **Output**: Optimized trajectory and energy profile
- **Connection**: Uses same constraint system as interactive mode

## **Module Connection Map**

### **Data Flow Architecture**
```
AtomTypes.dat → MMparams.py → MMFF.py → MolecularDynamics.py → relax_multi.cl
     ↑                                                            ↓
ExplorerVisPy ←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←
```

### **Function Call Hierarchy**
```
ExplorerVisPy._init_gridff_scanner()
  ↓
MolecularDynamics.initGridFF()
  ↓
MolecularDynamics.pack_system()
  ↓
ExplorerVisPy._setup_gridff_constraints()
  ↓
ExplorerVisPy._gridff_relax_step()
  ↓
MolecularDynamics.run_*() kernels
  ↓
ExplorerVisPy._update_gridff_visuals()
```

### **Buffer Management**
- **MMFF Topology**: `mmff.*` → GPU buffers via `pack_system()`
- **GridFF Data**: B-spline arrays → GPU texture/buffers via `initGridFF()`
- **Constraints**: `constr`/`constrK` → GPU via `toGPU()`
- **Positions**: GPU → CPU via `download_results()`

## **Implementation Dependencies**

### **Required Data Files**
1. **GridFF Data**: `cpp/common_resources/CaF2_6L_Ni3_rect_nx2_nz1_L2_top/Bspline_PLQd.npy`
2. **Molecule Files**: PTCDA.xyz (or user-specified)
3. **Parameter Files**: `AtomTypes.dat`, `ElementTypes.dat`
4. **Kernel Files**: `cpp/common_resources/cl/relax_multi.cl`

### **Parameter Requirements**
- **GridFF**: Grid shape, origin, spacing, B-spline coefficients
- **MMFF**: Bond connectivity, atom types, force constants
- **Constraints**: Anchor atom index, spring constant, target position
- **Visualization**: Colors, markers, plot ranges

### **Performance Considerations**
- **GPU Memory**: GridFF texture + MMFF topology + constraint buffers
- **Update Rate**: 20 FPS balance between responsiveness and stability
- **Kernel Launch**: Minimize CPU-GPU transfers in relaxation loop
- **Visualization**: Efficient VisPy marker updates

## **Testing Strategy**

### **Unit Tests**
1. **GridFF Loading**: Verify B-spline data loads correctly
2. **MMFF Topology**: Check PTCDA bond/angle generation
3. **Constraint Forces**: Test harmonic constraint implementation
4. **Energy Evaluation**: Validate GridFF force components

### **Integration Tests**
1. **Interactive Loop**: Test real-time relaxation stability
2. **Keyboard Control**: Verify constraint position updates
3. **Visualization**: Check position updates in 3D view
4. **Energy Conservation**: Monitor total energy during relaxation

### **Validation Tests**
1. **Force Field Parity**: Compare GridFF vs XYZ reference
2. **Constraint Response**: Verify anchor follows keyboard input
3. **Relaxation Convergence**: Test system reaches equilibrium
4. **XZ Slice Accuracy**: Validate cross-section visualization

## **Timeline and Milestones**

### **Phase 1: Core GridFF Integration (Week 1)**
- [ ] GridFF scanner instance creation
- [ ] B-spline data loading
- [ ] MMFF topology building
- [ ] Basic constraint setup

### **Phase 2: Interactive Controls (Week 2)**
- [ ] Timer-based relaxation loop
- [ ] Keyboard event handling
- [ ] Real-time visualization updates
- [ ] Energy display integration

### **Phase 3: Visualization Enhancements (Week 3)**
- [ ] XZ slice debugging tools
- [ ] Constraint visualization
- [ ] UI improvements and controls
- [ ] Performance optimization

### **Phase 4: Advanced Features (Week 4)**
- [ ] Batch scanning workflow
- [ ] Path optimization algorithms
- [ ] Documentation and tutorials
- [ ] Testing and validation

## **Success Criteria**

### **Functional Requirements**
- [ ] Real-time interactive manipulation of PTCDA on CaF2
- [ ] Harmonic constraint following keyboard control
- [ ] Stable relaxation loop with GridFF forces
- [ ] Accurate energy component evaluation
- [ ] Smooth 3D visualization updates

### **Performance Requirements**
- [ ] 20 FPS update rate during interaction
- [ ] <100ms response time to keyboard input
- [ ] Stable energy conservation during relaxation
- [ ] Memory usage < 2GB for typical systems

### **Usability Requirements**
- [ ] Intuitive keyboard controls
- [ ] Clear energy and position feedback
- [ ] Robust error handling and recovery
- [ ] Comprehensive documentation

## **Future Extensions**

### **Advanced Interaction Modes**
- Multi-point constraints for complex manipulation
- Path following with automatic constraint updates
- Force feedback for haptic devices
- VR/AR integration for immersive manipulation

### **Enhanced Analysis Tools**
- Energy landscape mapping and visualization
- Reaction pathway discovery
- Transition state search with constraints
- Machine learning guided manipulation

### **Performance Optimizations**
- Multi-GPU support for larger systems
- Adaptive time stepping for stability
- Cached force evaluations for repeated positions
- Compressed GridFF representations
