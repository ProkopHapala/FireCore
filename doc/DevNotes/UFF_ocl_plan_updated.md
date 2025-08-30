# UFF OpenCL Implementation - Updated Technical Plan

## Refined Development Strategy: A Two-Path Approach

The project will proceed along two parallel but complementary paths. This strategy allows for rapid, isolated testing of the low-level OpenCL kernels while systematically building the high-performance, integrated C++ solution.

### Path A: Direct PyOpenCL Kernel Validation
*   **Goal**: To quickly verify the correctness of each OpenCL kernel (`UFF.cl`) in isolation.
*   **Tools**: This path uses the direct Python-to-OpenCL wrapper `pyBall/OCL/UFF.py` and the test script `tests/tUFF/test_UFF_ocl.py`.
*   **Advantage**: This approach bypasses the C++ compilation cycle, allowing for fast iteration on kernel logic. It is the ideal environment for debugging the physics of each UFF component (bonds, angles, etc.) by comparing GPU results directly against the CPU "ground truth".

### Path B: High-Performance C++ Integration
*   **Goal**: To integrate the validated kernels into the `MolWorld_sp3_multi` framework for high-throughput, multi-replica simulations.
*   **Tools**: This path involves implementing `OCL_UFF.h`, integrating it into `MolWorld_sp3_multi.h`, and exposing it to Python via `MMFFmulti_lib.cpp`.
*   **Advantage**: This is the final, production-ready implementation that will deliver the highest performance and be fully integrated with the rest of the FireCore simulation ecosystem.

The plan below is structured to follow these two paths, ensuring that the kernels are validated (Path A) before being integrated into the more complex C++ environment (Path B).

## Critical Implementation Insights

### 1. Kernel Development Strategy

**Priority Order:**
1. **`evalBondsAndHNeigh_UFF`** - Foundation kernel that:
   - Calculates bond forces
   - Pre-computes normalized bond vectors (`hneigh`)
   - These `hneigh` vectors are consumed by angle/dihedral kernels
   - **Critical**: This kernel determines the efficiency of all subsequent kernels

2. **`evalAngles_UFF`** - Depends on `hneigh` from bonds kernel
3. **`evalDihedrals_UFF`** - Also depends on `hneigh`
4. **`evalInversions_UFF`** - Improper dihedrals
5. **`assembleForces_UFF`** - Final force assembly from `fint` buffer

### 2. Buffer Management Strategy

**Key Technical Details:**
- All buffers must be sized for `nSystems * nAtoms/nBonds/etc.`
- Use buffer index mapping: `iSys * nAtoms + iAtom`
- Intermediate force buffer (`fint`) prevents race conditions
- Atom-to-force mapping (`a2f_*`) enables efficient force assembly

### 3. Integration with MMFFBuilder Workflow

**Critical Path:**
```cpp
// In MolWorld_sp3_multi::realloc()
1. MMFFBuilder creates template UFF object
2. Allocate host arrays for nSystems
3. Call pack_uff_system() for each replica
4. Upload to GPU via uff_ocl->upload()
```

**New Function Required:**
```cpp
void pack_uff_system(int isys, const UFF& template_ff) {
    // Copy bond parameters: template_ff.bond_* → host_bond_params[isys*nBonds + ...]
    // Copy angle parameters: template_ff.angle_* → host_angle_params[isys*nAngles + ...]
    // Copy dihedral parameters: template_ff.dihedral_* → host_dihedral_params[isys*nDihedrals + ...]
    // Copy inversion parameters: template_ff.inversion_* → host_inversion_params[isys*nInversions + ...]
}
```

## Implementation Roadmap Refinement

### Phase 1: Core Kernel Development (High Priority)

#### Step 1.1: Verify/Create UFF.cl Kernels
**Location: `cpp/common_resources/cl/UFF.cl`**

**Critical Requirements:**
- Each kernel must handle multi-system indexing
- Use intermediate force buffer strategy
- Support component enable/disable flags
- Match CPU UFF.h calculation exactly

**Kernel Signatures:**
```opencl
__kernel void evalBondsAndHNeigh_UFF(
    __global float4* apos,           // Atom positions [nSystems*nAtoms]
    __global int2* bonds,            // Bond topology [nSystems*nBonds] 
    __global float4* bond_params,    // Bond parameters [nSystems*nBonds]
    __global float4* fint,           // Intermediate forces [nSystems*nForceComponents]
    __global float4* hneigh,         // Normalized bond vectors [nSystems*nBonds]
    int nSystems, int nAtoms, int nBonds,
    int bBonds                       // Component enable flag
);
```

#### Step 1.2: Direct Kernel Testing
**Use `pyBall/OCL/UFF.py` for isolated testing:**
- Load simple molecules (H2O, CH4, C2H6)
- Test each kernel individually
- Compare intermediate results against CPU calculations
- Validate `hneigh` vector calculations

### Phase 2: C++ Integration Layer

#### Step 2.1: Complete OCL_UFF.h Implementation
**Based on existing progress, focus on:**

```cpp
class OCL_UFF : public OCL {
    // Multi-system buffer management
    void realloc(int nSystems, int nAtoms, int nBonds, int nAngles, int nDihedrals, int nInversions);
    
    // Parameter upload for multi-system
    void upload_topology(const std::vector<int2>& bonds, const std::vector<int3>& angles, ...);
    void upload_parameters(const std::vector<float4>& bond_params, ...);
    
    // Component control
    void set_component_flags(bool bBonds, bool bAngles, bool bDihedrals, bool bInversions);
    
    // Evaluation
    void eval_UFF(bool bClearForce=true);
    
    // Results retrieval
    void download_forces(std::vector<float4>& forces);
    double get_total_energy();
};
```

#### Step 2.2: Integration into MolWorld_sp3_multi.h
**Key additions:**
```cpp
class MolWorld_sp3_multi {
    OCL_UFF* uff_ocl;
    bool bUFF;  // Runtime switch
    
    // Host arrays for multi-system UFF data
    std::vector<float4> host_bond_params;
    std::vector<float4> host_angle_params;
    std::vector<float4> host_dihedral_params;
    std::vector<float4> host_inversion_params;
    
    // New methods
    void pack_uff_system(int isys, const UFF& template_ff);
    void init_UFF_mode();
};
```

### Phase 3: Python Interface Enhancement

#### Step 3.1: Complete pyBall/OCL/UFF.py
**Missing critical methods:**
```python
def run_eval_step(self, bClearForce=True):
    """Execute UFF evaluation on GPU"""
    
def get_forces(self, iSys=0):
    """Download forces for specific system"""
    
def set_component_flags(self, bBonds=True, bAngles=True, bDihedrals=True, bInversions=True):
    """Control UFF components"""
```

#### Step 3.2: Create UFF_multi.py Wrapper
**Following MMFF_multi.py pattern:**
```python
class UFF_multi:
    def set_forcefield(self, ff_name="UFF"):
        """Switch to UFF mode"""
        
    def set_uff_switches(self, bBonds=True, bAngles=True, bDihedrals=True, bInversions=True):
        """Control UFF components"""
```

## Critical Technical Considerations

### 1. Parameter Compatibility
- UFF and MMFF have different parameter structures
- Ensure `MMFFBuilder.toUFF()` correctly populates all required fields
- Validate parameter units and conventions

### 2. Force Assembly Correctness
- The `fint` buffer strategy is crucial for correctness
- Each force component must write to unique buffer locations
- Final assembly must correctly sum all contributions

### 3. Performance Optimization
- Pre-computed `hneigh` vectors are essential for angle/dihedral efficiency
- Buffer layout affects memory coalescing
- Work-group size optimization for different kernel types

### 4. Debugging Strategy
- Start with single-system testing using `pyBall/OCL/UFF.py`
- Validate each kernel individually before integration
- Use CPU UFF.h as ground truth for all comparisons

## Risk Mitigation

### High Risk: Kernel Correctness
- **Mitigation**: Extensive single-kernel testing with known molecules
- **Validation**: Compare every intermediate result against CPU

### Medium Risk: Multi-System Integration
- **Mitigation**: Start with single system, then scale to multiple
- **Validation**: Ensure each replica produces identical results

### Low Risk: Python Interface
- **Mitigation**: Follow established MMFF patterns exactly
- **Validation**: Direct comparison with existing MMFF_multi.py workflow

## Success Metrics

1. **Kernel Validation**: Each kernel produces bit-identical results to CPU (within floating-point precision)
2. **Component Isolation**: Individual UFF terms can be enabled/disabled and tested separately
3. **Multi-System Scaling**: Performance scales linearly with number of systems
4. **Integration Completeness**: Runtime switching between UFF/MMFF works seamlessly
5. **Maintainability**: Code follows FireCore patterns and is well-documented