I want to implement UFF on GPU using OpenCL and test it. Making simple test script which check forces on atom calculated using CPU version of UFF and GPU/OpenCL version is integral part of this task. We need to be able to check difference of total energy and forces from each method and see where are discrepancies. It should be also to possible turn of and off different components (bond-strech, angles, torsion, improper dihedrals) separately so that we can check each component separately.

Your job is to review all relevant files, and note where is what (relevant functions), note also what is missing and what shoud be implemented, and based on that make step-by-step plan of implementation, noting exactly whare makes sense to implement what. But do not write the actuall code yet in this phase.
You should write all these notes in a single markdonw file UFF_ocl_plan.md

In second step you shoud write the test scrip in python which include MMFF_multi.py and depending on some value of internal switches can switch it to use UFF rather than MMFF on either CPU or GPU, and turn on/off the individual components. In this cript we do not run dynamics, nor initialize surface, and the non-covalent interaction should be off. We only evaluate force for given geometry of simple molecule (like H20.mol).


These are relevant parts where you shoud start   

* /home/prokophapala/git/FireCore/cpp/common/molecular/UFF.h   : C++ implementation of UFF solvers
* /home/prokophapala/git/FireCore/tests/tDFT/data/cl/UFF.cl    : OpenCL kernels for UFF solver 
* /home/prokophapala/git/FireCore/pyBall/OCL/UFF.py            : python pyOpenCL interface to OpenCL UFF solver (not yet tested)
* /home/prokophapala/git/FireCore/cpp/common/OpenCL/OCL_MM.h   : C++ interface to OpenCL version of MMFFsp3_loc forcefield (this is alternative of UFF which we implemented earlier, and now we want to be able to switch to UFF by some runtime flag)
   * /home/prokophapala/git/FireCore/cpp/common/OpenCL/OCL.h   : OpenCL integration utilities in C/C++ used in OCL_MM.h (perhaps part of it can be replaced by official OpenCL C++ interface)
* cpp/common/molecular/MMFFsp3_loc.h                           : CPU implementation of MMFFsp3 forcefield (alternative to UFF). I show it because we want to create OCL_UFF.h which has the same relation to UFF.h as the relation between OCL_MM.h and MMFFsp3_loc.h 
* /home/prokophapala/git/FireCore/cpp/common/molecular/MolWorld_sp3.h        : Top level class for Moleculer-Dynamics simualtions which includes many different methods incluging UFF and MMFFsp3_loc and allows to switch between them 
* /home/prokophapala/git/FireCore/cpp/common/molecular/MolWorld_sp3_multi.h  : GPU accelerated version of MolWorld_sp3.h, which currently offers GGPU acceleration only for MMFFsp3 (by OCL_MM.h) now we want to add  OCL_UFF.h and be able to switch between them similarly as we can switch in MolWorld_sp3.h between CPU versions of UFF and MMFFsp3 forcefields.

# GUI app integration
* /home/prokophapala/git/FireCore/cpp/apps/MolecularEditor/MolGUIapp.cpp              ; graphical GUI application which integrates  MolWorld_sp3.h  with GUI (   MolGUI.h )
* /home/prokophapala/git/FireCore/cpp/apps_OCL/MolecularEditorOCL/MolGUIapp_multi.cpp : graphical GUI app which integrates GPU accelerated MolWorld_sp3_multi.h  with  MolGUI.h 
   
# python / library integration
*   /home/prokophapala/git/FireCore/cpp/libs/Molecular/MMFF_lib.cpp   : includes MolWorld_sp3.h  and exposes its function in `extern "C"` block so they can be easily called externally (bypasing C++ name mangling)
*   /home/prokophapala/git/FireCore/cpp/libs_OCL/MMFFmulti_lib.cpp    : includes  MolWorld_sp3_multi.h  and exposes its function in `extern "C"` block so they can be easily called externally (bypasing C++ name mangling)
*   /home/prokophapala/git/FireCore/pyBall/MMFF_multi.py              : python wraper which loads MMFFmulti_lib.cpp  using ctypes
*   /home/prokophapala/git/FireCore/pyBall/OCL/MMFF.py                : python wraper which loads MMFF_lib.cpp   using ctypes


---

I want to implement Universal Force Field (UFF) GPU acceleration using OpenCL and create comprehensive testing infrastructure to validate the implementation. This task has two main phases:

## Phase 1: Analysis and Planning
Review the existing FireCore codebase to understand the current UFF and MMFF implementations, identify what components exist, what's missing, and create a detailed implementation plan. Specifically:

1. **Analyze existing implementations:**
   - C++ UFF solver: `cpp/common/molecular/UFF.h`
   - OpenCL UFF kernels: `tests/tDFT/data/cl/UFF.cl`
   - Python OpenCL interface: `pyBall/OCL/UFF.py` (untested)
   - Reference MMFF implementations for comparison:
     - CPU: `cpp/common/molecular/MMFFsp3_loc.h`
     - OpenCL: `cpp/common/OpenCL/OCL_MM.h`
     - Integration: `cpp/common/molecular/MolWorld_sp3.h` and `cpp/common/molecular/MolWorld_sp3_multi.h`

2. **Document integration points:**
   - Library interfaces: `cpp/libs/Molecular/MMFF_lib.cpp` and `cpp/libs_OCL/MMFFmulti_lib.cpp`
   - Python wrappers: `pyBall/MMFF_multi.py` and `pyBall/OCL/MMFF.py`
   - GUI applications: `cpp/apps/MolecularEditor/MolGUIapp.cpp` and `cpp/apps_OCL/MolecularEditorOCL/MolGUIapp_multi.cpp`

3. **Create implementation plan:**
   - Identify missing components needed to create `OCL_UFF.h` (analogous to `OCL_MM.h`)
   - Plan integration into `MolWorld_sp3_multi.h` for runtime switching between UFF and MMFF
   - Document required changes to library interfaces and Python wrappers
   - Note dependencies on `cpp/common/OpenCL/OCL.h` utilities

**Deliverable:** Write all findings and the step-by-step implementation plan in `UFF_ocl_plan.md`

## Phase 2: Test Script Implementation
Create a Python test script that validates UFF forces and energies between CPU and GPU implementations with the following requirements:

1. **Core functionality:**
   - Use `MMFF_multi.py` as base but switch to UFF instead of MMFF via runtime flags
   - Support both CPU and GPU/OpenCL execution modes
   - Load simple molecules (e.g., `H2O.mol`) for testing
   - Calculate and compare total energy and per-atom forces between CPU and GPU

2. **Component isolation:**
   - Allow enabling/disabling individual UFF components:
     - Bond stretching
     - Angle bending  
     - Proper dihedrals (torsions)
     - Improper dihedrals (inversions)
   - This enables testing each component separately to isolate discrepancies

3. **Constraints:**
   - No molecular dynamics simulation (static geometry evaluation only)
   - No surface initialization
   - Non-covalent interactions disabled
   - Focus purely on covalent UFF terms

4. **Output:**
   - Report energy differences between CPU/GPU
   - Report force differences per atom and component
   - Highlight discrepancies above specified thresholds
   - Support verbose debugging output

The test script should serve as both validation tool and reference implementation for integrating UFF into the broader FireCore ecosystem.

---

# UFF OpenCL Implementation Plan

## Current State Analysis

### Existing UFF Components

#### 1. C++ UFF Implementation (`cpp/common/molecular/UFF.h`)
**Status: Complete CPU implementation**
- Core UFF force field evaluation with bonds, angles, dihedrals, inversions
- OpenMP parallelization support
- Non-bonded interactions (LJ + Coulomb) with PBC support
- Force assembly and optimization loops
- Key methods:
  - `evalBonds()`, `evalAngles()`, `evalDihedrals()`, `evalInversions()`
  - `assembleAtomForce()` for force accumulation
  - `run()` main optimization loop with configurable components

#### 2. OpenCL UFF Kernels (`tests/tDFT/data/cl/UFF.cl`)
**Status: Needs verification - location unclear from retrieval**
- GPU kernels for UFF force evaluation
- **TODO: Verify existence and completeness**

#### 3. Python OpenCL Interface (`pyBall/OCL/UFF.py`)
**Status: Partially implemented, untested**
- Inherits from `OpenCLBase` 
- Has buffer allocation (`realloc_buffers()`)
- Has parameter upload (`upload_topology_params()`)
- Has UFF parameter assignment (`toUFF()`, `assignUFFParameters()`)
- **Missing:**
  - Kernel execution methods
  - Force/energy retrieval
  - Component enable/disable flags
  - Integration with main evaluation loop

### Reference MMFF Implementation (for comparison)

#### 1. CPU MMFF (`cpp/common/molecular/MMFFsp3_loc.h`)
- Similar structure to UFF.h
- Bond, angle, torsion evaluation
- Force assembly patterns

#### 2. OpenCL MMFF (`cpp/common/OpenCL/OCL_MM.h`)
**Status: Complete reference implementation**
- Shows integration pattern with `OCL.h` utilities
- Buffer management and kernel execution
- Force/energy retrieval methods
- **Key pattern to replicate for UFF**

#### 3. Integration Layer (`cpp/common/molecular/MolWorld_sp3_multi.h`)
**Status: MMFF integrated, UFF missing**
- Runtime switching between CPU/GPU modes
- **Missing: UFF integration**

### Integration Infrastructure

#### 1. Library Interfaces
- `cpp/libs/Molecular/MMFF_lib.cpp` - CPU version
- `cpp/libs_OCL/MMFFmulti_lib.cpp` - GPU version
- **Missing: UFF equivalents**

#### 2. Python Wrappers
- `pyBall/MMFF_multi.py` - wraps MMFFmulti_lib
- `pyBall/OCL/MMFF.py` - direct OpenCL interface
- **Missing: UFF equivalents**

### Parameter Assignment Workflow (`MMFFBuilder`)

The `MMFFBuilder` class (`cpp/common/molecular/MMFFBuilder.h`) is central to the force field setup process in `FireCore`. It acts as a versatile tool for constructing molecular topology (bonds, angles, etc.) and assigning the corresponding interaction parameters, regardless of the target force field (UFF or MMFF). This process is orchestrated within the `MolWorld_sp3::makeMMFFs()` method.

The workflow is as follows:

1.  **Parameter Loading:** The `MMFFparams` class (`cpp/common/molecular/MMFFparams.h`) first loads the raw force field parameters from `.dat` files (e.g., `ElementTypes.dat`, `AtomTypes.dat`, `BondTypes.dat`). These files contain tabulated values for atom properties, bond lengths, stiffness constants, etc.

2.  **Topology and Type Assignment:** The `MMFFBuilder` instance (`builder`) uses this parameter database to:
    *   Build the molecular graph (finding bonds from atomic positions).
    *   Assign initial atom types.
    *   Perform advanced type assignment based on the chemical environment. This is where the logic diverges for UFF and MMFF:
        *   **For UFF:** The `builder.assignUFFtypes()` and `builder.assignUFFparams()` methods are called. These functions implement the rules from the UFF paper to dynamically determine atom types (e.g., `C_3`, `O_R`) and calculate interaction parameters (bond lengths, angles, force constants) based on elemental properties, hybridization, and bond orders.
        *   **For MMFF:** The `builder.assignTypes()` method is called, which typically relies on more specific, pre-tabulated parameters found in the `.dat` files.

3.  **Data Transfer to Force Field Object:** Once the `builder` has fully constructed the system with all its parameters, the data is transferred to a dedicated, optimized force field object.
    *   `builder.toUFF(ffu, ...)`: Populates a `UFF` object (`ffu`) with the calculated parameters.
    *   `builder.toMMFFsp3_loc(ffl, ...)`: Populates an `MMFFsp3_loc` object (`ffl`).

#### How this applies to the GPU Implementation:

This established workflow is critical for the GPU implementation. The `MolWorld_sp3_multi.h` class will leverage this exact process to prepare a *template* CPU force field object.

-   A single `UFF` object will be created and populated using `builder.toUFF()`.
-   The `MolWorld_sp3_multi::pack_system()` function will then read the parameters from this template `UFF` object and copy them into the large, contiguous host arrays for each of the `nSystems` replicas.
-   This ensures that the parameter generation logic is centralized and consistent between the CPU and GPU paths. The only new requirement for the GPU implementation is the logic to correctly pack the UFF parameters (bond params, angle params, etc.) into the multi-system arrays before the bulk upload to the GPU.

### Multi-Replica Simulation on GPU

The `MolWorld_sp3_multi.h` class is designed to accelerate simulations by running multiple replicas of the same molecular system in parallel on the GPU. This is particularly effective for small systems where a single simulation would not fully utilize the GPU's parallel processing capabilities.

The core strategy involves:
1.  **Data Replication on the Host:** Instead of creating multiple full C++ force field objects in memory, large contiguous arrays are allocated on the host (CPU memory) to hold the data for all replicas. These arrays are sized `nSystems * per_system_size`. Key examples from `MolWorld_sp3_multi.h` include `atoms`, `aforces`, `REQs`, and `MMpars`.
2.  **Packing Parameters:** A function, `pack_system(isys, ff, ...)` is used to copy the parameters from a single template force field object (e.g., `ffl` of type `MMFFsp3_loc` or `ffu` of type `UFF`) into the appropriate slice of the large host arrays for a given replica index `isys`. This allows for efficient replication of the system parameters. It also provides a point to introduce variations between replicas, such as different initial positions or velocities. The `pack_systems()` method orchestrates this for all replicas. See the `Parameter Assignment Workflow` section for details on how the template force field object is populated.
3.  **Bulk GPU Upload:** Once the host-side arrays are fully packed with data for all replicas, they are uploaded to the corresponding GPU buffers in a single, efficient bulk transfer using the `upload()` method. This method copies data from the host arrays (e.g., `atoms`) to the GPU buffers managed by the `OCL_MM` class (e.g., `ibuff_atoms`).
4.  **Parallel Kernel Execution:** The OpenCL kernels (e.g., in `OCL_MM.h`) are written to work on these large, multi-system buffers. Each work-group on the GPU is typically assigned to a single replica, identified by `get_group_id(1)`. This allows all replicas to be processed in parallel.

#### How this applies to UFF implementation:

This multi-replica architecture will be directly applicable to the new `OCL_UFF.h` implementation.

-   **`OCL_UFF.h` Buffers:** The buffers declared in `OCL_UFF.h` (e.g., `buf_bonds`, `buf_bond_params`) will need to be allocated on the GPU to hold data for `nSystems`.
-   **`MolWorld_sp3_multi.h` Integration:**
    -   A new `OCL_UFF* uff_ocl` member will be added.
    -   The `realloc()` function will be updated to also initialize the UFF-specific buffers in `uff_ocl` for `nSystems`.
    -   A new packing function, analogous to `pack_system`, will be needed to copy UFF parameters from a template `UFF` object into the multi-system host arrays.
    -   The `upload()` function will be extended to transfer these UFF-specific host arrays to the GPU.
-   **Python Wrappers:** The Python test script and wrappers will need to manage the creation and configuration of these multiple replicas, although the C++ library will handle the low-level details of packing and uploading.

This approach ensures that the UFF implementation can leverage the existing high-throughput simulation framework in `FireCore`.

## Missing Components Analysis

### Critical Missing Components

1. **`cpp/common/OpenCL/OCL_UFF.h`**
   - OpenCL wrapper for UFF (analogous to OCL_MM.h)
   - Buffer management for UFF-specific data structures
   - Kernel execution interface
   - Force/energy retrieval

2. **UFF Integration in `MolWorld_sp3_multi.h`**
   - Runtime switching between UFF and MMFF
   - UFF-specific initialization
   - Parameter passing interface

3. **Library Interface for UFF**
   - `cpp/libs_OCL/UFFmulti_lib.cpp`
   - Extern "C" functions for Python integration

4. **Python Wrapper for UFF**
   - `pyBall/UFF_multi.py`
   - ctypes interface to UFFmulti_lib

5. **Component Control Interface**
   - Flags to enable/disable UFF components
   - Both in C++ and Python interfaces

### OpenCL Kernel Verification Needed

1. **Verify `tests/tDFT/data/cl/UFF.cl` exists and is complete**
2. **Check kernel signatures match expected interface**
3. **Validate against CPU implementation**

## Implementation Plan

### Phase 1: Core OpenCL UFF Integration

#### Step 1.1: Create `OCL_UFF.h`
**Location: `cpp/common/OpenCL/OCL_UFF.h`**
**Pattern: Follow `OCL_MM.h` structure**

```cpp
class OCL_UFF : public OCL {
    // Buffer declarations for UFF-specific data
    cl_mem buf_bonds, buf_angles, buf_dihedrals, buf_inversions;
    cl_mem buf_bond_params, buf_angle_params, buf_dihedral_params, buf_inversion_params;
    cl_mem buf_forces, buf_energies;
    
    // Component enable/disable flags
    bool bBonds, bAngles, bDihedrals, bInversions;
    
    // Methods
    void init_UFF();
    void upload_topology();
    void upload_parameters();
    void eval_UFF(bool bClearForce=true);
    void download_forces();
    double get_total_energy();
};
```

#### Step 1.2: Verify/Complete OpenCL Kernels
**Location: Verify `tests/tDFT/data/cl/UFF.cl` or create in `cpp/common_resources/cl/UFF.cl`**
- Ensure kernels for bonds, angles, dihedrals, inversions
- Add component enable/disable support
- Validate kernel signatures

#### Step 1.3: Complete `pyBall/OCL/UFF.py`
**Missing methods to implement:**
```python
def run_eval_step(self, bClearForce=True):
    """Execute UFF evaluation on GPU"""
    
def get_forces(self, iSys=0):
    """Download forces from GPU"""
    
def get_energies(self):
    """Download energies from GPU"""
    
def set_component_flags(self, bBonds=True, bAngles=True, bDihedrals=True, bInversions=True):
    """Enable/disable UFF components"""
```

### Phase 2: Integration Layer

#### Step 2.1: Integrate UFF into `MolWorld_sp3_multi.h`
**Location: `cpp/common/molecular/MolWorld_sp3_multi.h`**
```cpp
class MolWorld_sp3_multi {
    OCL_UFF* uff_ocl;  // Add UFF OpenCL instance
    bool bUseUFF;      // Runtime switch between UFF/MMFF
    
    void init_UFF_OCL();
    void switch_to_UFF();
    void switch_to_MMFF();
};
```

#### Step 2.2: Create UFF Library Interface
**Location: `cpp/libs_OCL/UFFmulti_lib.cpp`**
```cpp
extern "C" {
    void init_UFF_multi(const char* xyz_file);
    void set_UFF_components(int bBonds, int bAngles, int bDihedrals, int bInversions);
    void eval_UFF_forces();
    double get_UFF_energy();
    void get_UFF_forces_array(double* forces);
}
```

#### Step 2.3: Create Python UFF Wrapper
**Location: `pyBall/UFF_multi.py`**
```python
class UFF_multi:
    def __init__(self):
        self.lib = ctypes.CDLL("./libs_OCL/UFFmulti_lib.so")
        
    def init(self, xyz_file):
        """Initialize UFF system"""
        
    def set_components(self, bBonds=True, bAngles=True, bDihedrals=True, bInversions=True):
        """Enable/disable UFF components"""
        
    def eval_forces(self):
        """Evaluate forces and return energy"""
        
    def get_forces(self):
        """Get force array"""
```

### Phase 3: Testing Infrastructure

#### Step 3.1: Create Component Testing Interface
**Add to both C++ and Python:**
- Individual component evaluation methods
- Energy/force breakdown by component
- Comparison utilities

#### Step 3.2: Validation Against CPU Implementation
- Force comparison utilities
- Energy comparison utilities
- Threshold-based discrepancy detection

## Implementation Priority

### High Priority (Core Functionality)
1. **Verify/complete OpenCL kernels** - Foundation for everything
2. **Complete `pyBall/OCL/UFF.py`** - Direct testing interface
3. **Create basic test script** - Immediate validation capability

### Medium Priority (Integration)
4. **Create `OCL_UFF.h`** - C++ integration layer
5. **Integrate into `MolWorld_sp3_multi.h`** - Runtime switching
6. **Create library interface** - Python integration

### Low Priority (Polish)
7. **GUI integration** - User interface
8. **Advanced testing features** - Comprehensive validation

## Testing Strategy

### Phase 1: Direct OpenCL Testing
- Use `pyBall/OCL/UFF.py` directly
- Compare against CPU UFF implementation
- Test individual components

### Phase 2: Integration Testing  
- Test through `UFF_multi.py` wrapper
- Validate runtime switching
- Performance benchmarking

### Phase 3: Application Testing
- GUI application testing
- Complex molecule validation
- Stress testing

## Risk Assessment

### High Risk
- **OpenCL kernel completeness** - May need significant development
- **Parameter compatibility** - UFF vs MMFF parameter differences

### Medium Risk
- **Integration complexity** - Multiple layers of abstraction
- **Performance optimization** - GPU efficiency tuning

### Low Risk
- **Python wrapper creation** - Well-established pattern
- **Testing infrastructure** - Clear requirements

## Success Criteria

1. **Functional parity** - GPU UFF matches CPU UFF results within tolerance
2. **Component isolation** - Individual UFF terms can be tested separately  
3. **Integration completeness** - Runtime switching between UFF/MMFF works
4. **Performance improvement** - GPU implementation faster than CPU for relevant system sizes
5. **Maintainability** - Code follows established FireCore patterns