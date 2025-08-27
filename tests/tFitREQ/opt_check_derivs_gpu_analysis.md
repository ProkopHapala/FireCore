# CPU vs GPU FitREQ Implementation Analysis

## Overview

This document analyzes the differences between the CPU/C++ implementation and GPU/OpenCL implementation of the FitREQ molecular fitting system. The analysis focuses on identifying why the two implementations might produce different energies and fitting errors.

## Relevant Files

The following files are relevant to this analysis and were modified during the debugging process:

### 1. CPU Implementation Files
- **`cpp/common/molecular/FitREQ.h`** - CPU C++ solver class with core functions (`evalExampleDerivs_MorseQH2`, etc.)
- **`cpp/libs/Molecular/FitREQ_lib.cpp`** - CPU library interface eporting pure C functions (to avoid C++ name mangling)
- **`pyBall/FitREQ.py`** - Python module binding to the C/C++ library via `ctypes`

### 2. GPU Implementation Files  
- **`cpp/common_resources/cl/FitREQ.cl`** - GPU OpenCL kernel implementation (`evalSampleDerivatives_template`)
- **`pyBall/OCL/NonBondFitting.py`** - python pyOpenCL class (`FittingDriver`) which prepates the data and provides the interface to the OpenCL kernels 

### 3. Test Scripts and Input Files
- **`tests/tFitREQ/opt_check_derivs_gpu.py`** - Main Python test script for comparing CPU vs GPU implementations
- **`H2O_single.xyz`** - Test XYZ file with water dimer configuration
- **`dofSelection_H2O_Morse.dat`** - DOF selection file defining parameters to fit

### Key Functions/Methods
- **Python**: `setup_cpu_fit()`, `scan_dof()`, `plot_overlays()`
- **CPU C++**: `evalExampleDerivs_MorseQH2()`, `evalSample()`
- **GPU OpenCL**: `evalSampleDerivatives_template` kernel
- **GPU Python**: `FittingDriver.getErrorDerivs()`

## Test Script Flow

The main test script `opt_check_derivs_gpu.py` follows this high-level flow:

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

## Key Differences Identified

### 1. Data Input and Processing

**CPU Version:**
- Uses `FitREQ.loadXYZ(fname, bAddEpairs=True, bOutXYZ=False)` 
- Automatically adds electron pairs (dummy atoms) during loading
- Stores data in `Atoms` objects with `AddedData` structure containing:
  - Bonds and directions for dummy atoms
  - Host atom indices for dummy atoms
  - Fragment separation index (`HBn0`)

**GPU Version:**
- Loads XYZ data into flat numpy arrays (`host_atoms`, `host_atypes`, `host_ranges`)
- No automatic dummy atom generation during loading
- Processes data into OpenCL buffer format

**Impact:** The CPU version includes electron pair corrections by default, while the GPU version may not be handling dummy atoms consistently.

### 2. Parameter Storage and Mixing Rules

**Both versions store parameters as:**
```
REQH = (RvdW, sqrt(EvdW), Qcharge, Hcorrection)
```

**Key Difference in EvdW Storage:**
- CPU: `initAllTypes()` stores `Quat4d{ RvdW, sqrt(EvdW), Qbase, Hb }`
- GPU: Stores `sqrt(EvdW)` in `tREQHs_base` but may have different scaling

**Mixing Rules:**
- Both use Lorentz-Berthelot mixing: `E0 = sqrt(Ei * Ej)`, `R0 = Ri + Rj`
- GPU implementation: Uses macro-injected code from `Forces.cl`

### 3. Model Evaluation

**CPU Version:**
- C++ evaluators: `evalExampleDerivs_*()` functions
- Supports multiple models: LJQH2, LJr8QH2, MorseQ
- Handles dummy atoms with special logic for electron pairs
- Uses different evaluators for corrected vs uncorrected calculations

**GPU Version:**
- OpenCL kernels with macro injection
- `evalSampleDerivatives_template` kernel
- Model-specific code injected via `MODEL_PAIR_ACCUMULATION` macro
- Currently limited to specific models (MorseQ, LJQH2, etc.)

**Potential Issue:** The GPU version may not fully implement the same physics as the CPU version, especially regarding dummy atom handling.

### 4. Derivative Calculation and DOF Assembly

**CPU Version:**
- Direct evaluation of derivatives using chain rule
- Accumulates per-atom derivatives into DOF gradients
- Handles electron pair charge conservation via `acumHostDerivs()`

**GPU Version:**
- Two-kernel approach:
  1. `evalSampleDerivatives_template` - computes per-sample energies and per-atom derivatives
  2. `assembleAndRegularize` - gathers derivatives into DOF forces
- Supports dual-pass evaluation for both fragments
- Uses mapping arrays (`DOFtoAtom`, `DOFcofefs`) for assembly

**Key GPU Feature:** Fragment-1 restriction - derivatives are only computed for the first fragment in single-pass mode, requiring dual-pass for complete coverage.

### 5. Electron Pair / Dummy Atom Handling

**CPU Version:**
- Full support for dummy atoms (electron pairs, sigma holes)
- `MMFFBuilder` generates molecular topology and adds dummy atoms
- Separate charge handling for dummy atoms
- Short-range corrections for electron pairs

**GPU Version:**
- Limited dummy atom support (marked as "not yet implemented" in documentation)
- No automatic dummy atom generation
- May not handle charge conservation between host and dummy atoms

**Major Issue:** This is likely the primary source of energy/fitting differences. The CPU version includes electron pair corrections by default (`bAddEpairs=True`), while the GPU version may not.

### 6. Regularization Implementation

**Both versions support:**
- Soft-wall potentials with `xlo`/`xhi` bounds
- Harmonic tethering to initial values (`K0`)
- Different stiffness parameters (`Klo`, `Khi`)

**Differences:**
- CPU: Integrated into C++ optimization loop
- GPU: Handled in `assembleAndRegularize` kernel
- Host-side regularization energy calculation may differ

## Recommended Debugging Steps

### 1. Verify Dummy Atom Handling
- Check if GPU version is properly loading/processing electron pairs
- Compare the number of atoms processed by CPU vs GPU
- Ensure charge conservation between host and dummy atoms

### 2. Compare Parameter Values
- Verify that `tREQHs_base` contains the same values in both implementations
- Check if `sqrt(EvdW)` scaling is applied consistently
- Compare mixing rule implementations

### 3. Energy Evaluation Consistency
- Run single-sample energy calculations and compare
- Check if the same model equations are being used
- Verify fragment separation logic (`HBn0`, `n0`)

### 4. Derivative Chain Rule
- Ensure the derivative calculations follow the same mathematical formulation
- Check if chain rule for `sqrt(EvdW)` is implemented identically
- Verify force sign conventions

## Input Files Used

- `H2O_single.xyz`: Single water dimer configuration with reference energy
- `dofSelection_H2O_MorseSR_plusE.dat`: DOF definitions for O_3 and H_O parameters

## Conclusion

The most likely source of discrepancies is the **dummy atom/electron pair handling**. The CPU version includes these corrections by default, while the GPU version appears to lack full implementation of this feature. This would cause significant differences in both energies and derivatives.

Secondary sources of differences could include:
- Subtle differences in the model evaluation code
- Different handling of fragment separation
- Implementation differences in the derivative chain rule

The debugging should start by ensuring both versions process the same molecular structures with identical parameters and corrections.
