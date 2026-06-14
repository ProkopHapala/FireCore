## UFF OpenCL Implementation Plan

### 1. High-Level Goal

The primary objective is to create a GPU-accelerated version of the UFF force field using OpenCL, fully integrated into the `FireCore` multi-replica simulation framework. This will allow for high-throughput simulations of molecular systems using UFF, with a testing mechanism to ensure correctness by comparing results against the existing CPU implementation. The implementation must support toggling individual energy terms (bonds, angles, etc.) for detailed validation.

### 2. System Architecture

The GPU-accelerated UFF will be integrated into the existing `MolWorld_sp3_multi` infrastructure, which is designed for running many replicas of a system in parallel. The data and execution flow will be as follows:

1.  **Python Layer (`MMFF_multi.py`)**: The user interacts with the system, loading a molecule and selecting the "UFF" force field.
2.  **C++ Library (`MMFFmulti_lib.cpp`)**: The Python call is forwarded to the C++ library wrapper.
3.  **Core C++ Orchestrator (`MolWorld_sp3_multi.h`)**:
    *   Receives the request to use UFF.
    *   Uses `MMFFBuilder` to create a single *template* instance of the CPU `UFF` force field object. This centralizes the complex parameter assignment logic.
    *   Allocates large, contiguous host-side arrays to hold the parameters and coordinates for all `nSystems` replicas.
    *   Calls a new `pack_uff_system()` function to copy parameters from the template `UFF` object into the correct slice of the host arrays for each replica.
    *   Calls `upload()` to perform a bulk transfer of these host arrays to the GPU buffers managed by `OCL_UFF`.
4.  **OpenCL Wrapper (`OCL_UFF.h`)**:
    *   Manages all UFF-specific GPU buffers (`cl_mem` objects).
    *   Manages and enqueues the OpenCL kernels from `UFF.cl`.
5.  **OpenCL Kernels (`UFF.cl`)**:
    *   Execute in parallel on the GPU, with each work-group typically processing one system replica.
    *   Kernels are designed to avoid race conditions by writing intermediate forces to a temporary buffer (`fint`), which are then summed in a final assembly step.

### 3. Component Analysis

#### Existing Components
- **`cpp/common/molecular/UFF.h`**: A complete and correct CPU implementation of the UFF force field. This will serve as the "ground truth" for validation.
- **`cpp/common/molecular/MMFFBuilder.h`**: The parameterization engine. Its `assignUFFtypes()` and `toUFF()` methods are essential for generating the template force field.
- **`cpp/common/molecular/MolWorld_sp3_multi.h`**: The multi-replica orchestrator. It already contains the logic for packing and running MMFF simulations on the GPU, which we will adapt for UFF.
- **`pyBall/OCL/UFF.py`**: A direct Python-to-OpenCL interface. While not used for the final high-performance path, it is an invaluable tool for interactively testing and debugging the `UFF.cl` kernels directly from Python.

#### Missing or Incomplete Components
- **`cpp/common/OpenCL/OCL_UFF.h`**: **This is the primary C++ component to be created.** It will be the OpenCL wrapper for UFF, analogous to `OCL_MM.h`. It will manage all UFF-specific GPU buffers and kernel execution.
- **`cpp/common_resources/cl/UFF.cl`**: The OpenCL kernel file. Kernels for bonds, angles, dihedrals, and inversions must be implemented and validated. A key feature will be the use of an intermediate force buffer to prevent race conditions.
- **Integration code in `MolWorld_sp3_multi.h`**: Logic to handle the UFF case, including calling the UFF-specific packing and upload functions.
- **Updates to `MMFFmulti_lib.cpp` and `pyBall/MMFF_multi.py`**: Exposing the new functionality (selecting UFF, toggling components) to the Python level.

### 4. Detailed Implementation Plan

This plan proceeds from the lowest level (OpenCL kernels) up to the Python interface.

#### Step 1: Implement and Test OpenCL Kernels (`UFF.cl`)
1.  **Create `UFF.cl`**: Place this file in `cpp/common_resources/cl/`.
2.  **Implement Kernels**:
    *   `evalBondsAndHNeigh_UFF`: This kernel should calculate bond forces and also compute and store normalized bond vectors (`hneigh`). This pre-calculation is crucial for the efficiency of subsequent kernels.
    *   `evalAngles_UFF`, `evalDihedrals_UFF`, `evalInversions_UFF`: These kernels will consume the pre-computed `hneigh` vectors to calculate angular forces.
    *   **Force Strategy**: Each of the above kernels should write their force contributions to a large, non-overlapping temporary buffer (`fint`). This avoids race conditions and the need for atomic operations.
    *   `assembleForces_UFF`: A final kernel that reads from the `fint` buffer and assembles the total force on each atom into the final force buffer (`fapos`). This requires an atom-to-force-piece mapping (`a2f_*` buffers).
3.  **Direct Kernel Testing**: Use the `pyBall/OCL/UFF.py` script to test each kernel individually.
    *   Load a simple molecule (e.g., water, ethane).
    *   Upload the parameters and positions.
    *   Run one kernel at a time and download the intermediate forces from `fint`.
    *   Compare these forces against a manual calculation or the CPU version to ensure correctness.

#### Step 2: Implement the C++ OpenCL Wrapper (`OCL_UFF.h`)
1.  **Create `OCL_UFF.h`**: This class will inherit from `OCL` and manage all UFF-related GPU resources.
2.  **Define Buffers**: Declare all necessary `cl_mem` buffer handles as integer indices, following the design of `OCL_MM.h`. This includes buffers for topology, parameters, intermediate forces (`ibuff_fint`), and the atom-to-force map (`ibuff_a2f_*`).
3.  **Implement `realloc()`**: This method will allocate all GPU buffers based on the system size (`nAtoms`, `nBonds`, etc.) and the number of replicas (`nSystems`).
4.  **Implement `makeKernels()` and `setup_kernels()`**: These methods will build the OpenCL program from `UFF.cl` and set the arguments for each kernel task.
5.  **Implement `eval()`**: This method will enqueue the kernels in the correct order: clear forces -> bonds -> angles/dihedrals/inversions -> assemble.

#### Step 3: Integrate into the Multi-Replica Orchestrator (`MolWorld_sp3_multi.h`)
1.  **Add UFF Members**: Add an `OCL_UFF* uff_ocl` pointer and a `bool bUFF` flag to the class.
2.  **Update `realloc()`**: When `bUFF` is true, call `uff_ocl->realloc()` to allocate the UFF-specific GPU buffers.
3.  **Create Host Arrays**: In the `realloc()` method, allocate host-side `std::vector`s for all UFF parameters (e.g., `host_bon_params`, `host_ang_params`). These will be sized for `nSystems`.
4.  **Implement `pack_uff_system(int isys, const UFF& ff)`**: This is a new, critical function. It will read parameters from the template `UFF` object `ff` and write them into the correct slice (for replica `isys`) of the host-side arrays.
5.  **Update `upload()`**: When `bUFF` is true, extend this function to upload the UFF-specific host arrays to the corresponding GPU buffers managed by `uff_ocl`.
6.  **Update `eval()`**: When `bUFF` is true, call `uff_ocl->eval()` instead of `mmff_ocl->eval()`.

#### Step 4: Expose Functionality to Python
1.  **Update `MMFFmulti_lib.cpp`**:
    *   Add an `extern "C"` function `setForcefield(const char* ffname)` that sets the `bUFF` flag in `MolWorld_sp3_multi`.
    *   Add `setUFFSwitches(bool bonds, bool angles, ...)` to control the component flags in the `OCL_UFF` instance.
2.  **Update `pyBall/MMFF_multi.py`**:
    *   Add Python methods that call the new C-level functions (e.g., `set_forcefield("UFF")`, `set_uff_switches(...)`).

#### Step 5: Create the Final Test Script
1.  **Create `tests/tUFF/test_GPU.py`**: This script will be the final validation tool.
2.  **Script Logic**:
    *   Initialize `MMFF_multi` and select the UFF force field.
    *   Load a molecule (e.g., `H2O.xyz`).
    *   Run a CPU-only evaluation using `MolWorld_sp3` to get reference forces and energy.
    *   Run a GPU evaluation using `MolWorld_sp3_multi`.
    *   Download the forces and energy from the GPU.
    *   Compare the CPU and GPU results and report any discrepancies.
    *   Implement loops to test each component (bonds, angles, etc.) in isolation by calling `set_uff_switches`.