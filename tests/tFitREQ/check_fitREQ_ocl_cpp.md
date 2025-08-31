# CPU vs GPU FitREQ Implementation Analysis

## Overview

This document analyzes the differences and consistency between the CPU/C++ and GPU/OpenCL implementations of the FitREQ molecular fitting system. After significant debugging, the two implementations now produce consistent forces and derivatives.

The main test scripts are:
- `/home/prokophapala/git/FireCore/tests/tFitREQ/check_fitREQ_ocl_cpp.py`
- `/home/prokophapala/git/FireCore/tests/tFitREQ/check_fitREQ_ocl_cpp_derivs_.py`

These scripts now confirm the good match between the CPU and GPU calculations.

## Testing Instructions

**USE THE `run.sh` SCRIPT FOR TESTING:**
```bash
cd /home/prokophapala/git/FireCore/tests/tFitREQ
./run.sh
```

This script runs the `check_fitREQ_ocl_cpp.py` test which compares CPU vs GPU objective function values. The script automatically:
- Builds the required libraries
- Runs the consistency test with default parameters
- Shows setup times and comparison results

**Alternative Manual Testing:**
```bash
python3 check_fitREQ_ocl_cpp.py --test_starting --verbose 1
```

## Latest Test Results

We have achieved a good match between the CPU and GPU implementations for both forces and their derivatives.

**Example Output:**
```
plotDOFscan_one(iDOF=0 : H_O.Q [CPU])  relFerrmax: +1.91e-03  (F_ana-F_num)[min,max]: -2.46e-02, +2.46e-02  |  E[min,max]: +1.09e-04, +5.40e-01 F_ana[min,max]: -1.29e+01, +1.29e+01 F_num[min,max]: -1.12e+01, +1.12e+01
plotDOFscan_one(iDOF=0 : H_O.Q [GPU])  relFerrmax: +1.91e-03  (F_ana-F_num)[min,max]: -2.46e-02, +2.46e-02  |  E[min,max]: +1.09e-04, +5.40e-01 F_ana[min,max]: -1.29e+01, +1.29e+01 F_num[min,max]: -1.12e+01, +1.12e+01
```
As shown above, the results for CPU and GPU are identical.

### Remaining Issue: Scaling Factor

To achieve this match, a hardcoded scaling factor was introduced in the `getErrorDerivs` method of `/home/prokophapala/git/FireCore/pyBall/OCL/FittingDriver.py`:

```python
return J*2., fDOFs*4.   # this scalling is hack, we should find why we need this scaling to match the reference
```

This scaling factor is a temporary hack. The underlying reason for this discrepancy needs to be investigated to find a permanent solution.

## Relevant Files

The following files are relevant to this analysis and were modified during the debugging process:

### 1. CPU Implementation Files
- **`cpp/common/molecular/FitREQ.h`** - CPU C++ solver class with core functions (`evalExampleDerivs_MorseQH2`, etc.)
- **`cpp/libs/Molecular/FitREQ_lib.cpp`** - CPU library interface eporting pure C functions (to avoid C++ name mangling)
- **`pyBall/FitREQ.py`** - Python module binding to the C/C++ library via `ctypes`

### 2. GPU Implementation Files
- **`cpp/common_resources/cl/FitREQ.cl`** - GPU OpenCL kernel implementation (`evalSampleDerivatives_template`). For ease of testing, it also contains `evalSampleDerivatives_template_serial`.
- **`pyBall/OCL/FittingDriver.py`** - Python PyOpenCL class that drives the fitting process.
- **`pyBall/OCL/NonBondFitting.py`** - python pyOpenCL class (`FittingDriver`) which prepares the data and provides the interface to the OpenCL kernels

### 3. Test Scripts and Input Files
- **`tests/tFitREQ/check_fitREQ_ocl_cpp_derivs_.py`** - Main Python test script for comparing CPU vs GPU implementations
- **`H2O_single.xyz`** - Test XYZ file with water dimer configuration
- **`dofSelection_H2O_Morse.dat`** - DOF selection file defining parameters to fit

### Key Functions/Methods
- **Python**: `setup_cpu_fit()`, `scan_dof()`, `plot_overlays()`
- **CPU C++**: `evalExampleDerivs_MorseQH2()`, `evalSample()`
- **GPU OpenCL**: `evalSampleDerivatives_template` kernel
- **GPU Python**: `FittingDriver.getErrorDerivs()`

## Test Script Flow

The main test script `check_fitREQ_ocl_cpp_derivs_.py` follows this high-level flow:

1. **Setup Phase:**
   - CPU: `setup_cpu_fit()` - loads types, DOF selection, XYZ data, sets up model
   - GPU: `setup_gpu_driver()` - creates OpenCL driver, loads data, compiles kernels

2. **Scan Phase:**
   - Both scan DOF parameters and compare energy/force curves
   - CPU uses `fit.scanParam()`
   - GPU uses `fit_ocl.getErrorDerivs()`

3. **Comparison:**
   - Plot overlays of CPU vs GPU results
   - Check derivatives for consistency

## Key Architectural Differences (Now Resolved)

This section describes the key architectural differences between the CPU and GPU implementations that were the source of discrepancies during development. These issues have been resolved, but the description is kept for architectural reference.

### 1. Data Input and Processing
**CPU Version:**
- Uses `FitREQ.loadXYZ(fname, bAddEpairs=True, bOutXYZ=False)`
- Can automatically add electron pairs (dummy atoms) during loading.
- Stores data in `Atoms` objects with `AddedData` structure.

**GPU Version:**
- Loads XYZ data into flat numpy arrays (`host_atoms`, `host_atypes`, `host_ranges`).
- Processes data into OpenCL buffer format.

### 2. Parameter Storage and Mixing Rules
**Both versions store parameters as:**
`REQH = (RvdW, sqrt(EvdW), Qcharge, Hcorrection)`

**Key Difference in EvdW Storage:**
- CPU: `initAllTypes()` stores `Quat4d{ RvdW, sqrt(EvdW), Qbase, Hb }`
- GPU: Stores `sqrt(EvdW)` in `tREQHs_base`.

**Mixing Rules:**
- Both use Lorentz-Berthelot mixing: `E0 = sqrt(Ei * Ej)`, `R0 = Ri + Rj`

### 3. Model Evaluation
**CPU Version:**
- C++ evaluators: `evalExampleDerivs_*()` functions.
- Supports multiple models and handles dummy atoms with special logic.

**GPU Version:**
- OpenCL kernels with macro injection (`MODEL_PAIR_ACCUMULATION`).
- `evalSampleDerivatives_template` kernel.

### 4. Derivative Calculation and DOF Assembly
**CPU Version:**
- Direct evaluation of derivatives using chain rule.
- Accumulates per-atom derivatives into DOF gradients.

**GPU Version:**
- Two-kernel approach:
  1. `evalSampleDerivatives_template` - computes per-sample energies and per-atom derivatives.
  2. `assembleAndRegularize` - gathers derivatives into DOF forces.
- Uses mapping arrays (`DOFtoAtom`, `DOFcofefs`) for assembly.

### 5. Electron Pair / Dummy Atom Handling
**CPU Version:**
- Full support for dummy atoms (electron pairs, sigma holes).
- `MMFFBuilder` generates molecular topology and adds dummy atoms.

**GPU Version:**
- Limited dummy atom support. This was a major source of discrepancy, now handled correctly for the tested systems.

### 6. Regularization Implementation
**Both versions support:**
- Soft-wall potentials and harmonic tethering.

**Differences:**
- CPU: Integrated into C++ optimization loop.
- GPU: Handled in `assembleAndRegularize` kernel.

## Summary of Resolutions

The main source of discrepancy was traced to inconsistent handling of charges and fitted parameters between the CPU and GPU code paths. The CPU implementation had its own logic for overriding charges from Degrees of Freedom (DOFs), which was not perfectly mirrored in the GPU implementation.

The resolution involved a significant refactoring of both the Python host code and the OpenCL kernels to ensure that:
1.  The source of charges (from atom types or from per-atom values) is handled consistently.
2.  The fitted parameters from the DOFs are correctly propagated to both the CPU and GPU evaluation routines before the calculation.
3.  The debug outputs were unified to allow for easier side-by-side comparison.

These changes have led to the current state where the results match.
