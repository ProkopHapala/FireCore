# Plan: Verifying eFF Force Calculations (CPU vs. GPU)

## 1. Overview

The primary objective is to validate that the pyOpenCL (GPU) implementation of the Electron Force Field (eFF) produces results identical to the reference C++ (CPU) implementation.

This will be accomplished by developing a test script that executes a single-point force calculation on a simple molecule (e.g., H₂O) with both implementations. The script will compare the resulting force vectors on all atoms and electrons to ensure they match exactly.

## 2. Key Files for Implementations

### C++ (CPU) Implementation
- **Core Logic**: `cpp/common/molecular/eFF.h`
- **C++ Interface**: `cpp/libs/Molecular/eFF_lib.cpp`
- **Python Wrapper**: `pyBall/eFF.py`

### pyOpenCL (GPU) Implementation
- **Python Host Code**: `pyBall/OCL/eFF_ocl.py`
- **GPU Kernel**: `cpp/common_resources/cl/eFF.cl`

## 3. Verification Strategy

The verification will proceed in two phases. Phase 1 will compare the total forces. If discrepancies are found, Phase 2 will enable a detailed, pairwise comparison of each force component.

### Phase 1: Total Force Comparison

A Python script (`test_ocl_vs_cpu.py`) will be created to perform the following:

1.  **Define Common Geometry**: A static H₂O molecule, including electron positions and sizes, will be hardcoded in an XYZ-like format.
2.  **Load into Both Models**: The geometry will be loaded into the C++ model via `pyBall/eFF.py` and into the GPU model via `pyBall/OCL/eFF_ocl.py`.
3.  **Verify Settings**: The script will print the active settings (e.g., `KRSrho` parameters, active force components) for both implementations to ensure they are identical.
4.  **Calculate Total Forces**: A single force evaluation will be triggered in both models for the static geometry.
5.  **Compare Output**: The total force on each particle will be printed in a compact, single-line format for easy comparison.

#### Proposed Output Format (Phase 1)

first we should check the ouput forces on particles (atoms, electrons)
```
--- CPU
CPU atom 0 ( Fx, Fy, Fz )
CPU atom 1 ( Fx, Fy, Fz )
CPU elec 0 ( Fx, Fy, Fz, Fw )

--- GPU
GPU atom 0 ( Fx, Fy, Fz )
GPU atom 1 ( Fx, Fy, Fz )
GPU elec 2 ( Fx, Fy, Fz, Fw )
...

Then if we find discrepancies, we need to check the individual
```
--- Pairwise Force Debugging ---
CPU EE(0,1) Coul: ( Fx, Fy, Fz ) | Fsi, Fsj
CPU AE(0,1) Paul: ( Fx, Fy, Fz ) | Fsi, Fsj
CPU AA(0,2) Coul: ( Fx, Fy, Fz )

GPU EE(0,1) Coul: ( Fx, Fy, Fz ) | Fsi, Fsj
GPU AE(0,1) Paul: ( Fx, Fy, Fz ) | Fsi, Fsj
GPU AA(0,2) Coul: ( Fx, Fy, Fz )
...

### Phase 2: Pairwise Force Component Debugging

If the total forces from Phase 1 do not match, we will activate detailed logging to inspect the individual pairwise interactions that contribute to the total force.

1.  **Enable Debug Prints**: We will uncomment and enhance existing `printf` statements within the core evaluation loops of both implementations.
    *   **C++ (`eFF.h`)**: Add `printf` inside the loops of `evalEE()`, `evalAE()`, etc., guarded by a `verbosity` flag.
    *   **GPU (`eFF.cl`)**: Utilize the `if (get_global_id(0) == idDBG)` guard to print forces for each interaction type (`Ion-Ion`, `Elec-Elec`, etc.) directly from the kernel.
2.  **Run Comparison Script Again**: The script will be re-run, now capturing the detailed pairwise output.
3.  **Analyze Component Forces**: This allows for a granular comparison of each component (Coulomb, Pauli) for each particle pair `(i, j)`, making it possible to pinpoint the exact source of any discrepancy.

#### Proposed Output Format (Phase 2)


```

This two-phase approach provides a clear and efficient path to validating the GPU implementation, starting with a high-level check and escalating to a detailed, component-wise analysis only if required.
