
# FitREQ Molecular Fitting Library Documentation

This document provides a comprehensive overview of the FitREQ system, a library designed for fitting non-covalent interaction parameters between rigid molecular fragments. It covers the underlying physical principles, implementation details, handling of dummy atoms and degrees of freedom, and practical guidance for users and developers.

## 1. Introduction to FitREQ

FitREQ is a specialized tool for optimizing parameters that describe non-covalent interactions, particularly Van der Waals (vdW) and electrostatic (Coulomb) forces, between molecular fragments. The primary objective is to minimize the root-mean-square error (RMSE) between calculated interaction energies (from the model) and reference interaction energies (e.g., from high-level quantum mechanical calculations) over a training set of molecular samples.

The system is built upon a C++ core for performance, with a C-compatible interface exposed for easy integration, and Python bindings for user-friendly interaction and data preparation.

## 2. Physical Principles and Fitting Objective

The core idea is to refine atomic parameters that govern non-covalent interactions. Each atom type is associated with a set of parameters, typically a 4-vector `REQH` = (RvdW, EvdW, Qcharge, Hcorrection), representing:
*   `RvdW`: Van der Waals radius
*   `EvdW`: Van der Waals well depth
*   `Qcharge`: Partial atomic charge
*   `Hcorrection`: A correction term, often related to hydrogen bonding or other specific interactions.

The fitting process involves:
1.  **Defining a Model:** The interaction energy between two molecular fragments is calculated based on the current `REQH` parameters of their constituent atoms. This typically involves summing pairwise vdW and Coulomb interactions.
2.  **Reference Data:** A training set consists of multiple molecular configurations (samples), each with a known reference interaction energy (e.g., from DFT or other accurate methods).
3.  **Error Minimization:** The system calculates the difference between the model-predicted energy and the reference energy for each sample. The objective function is the RMSE of these energy differences.
4.  **Optimization:** Gradient-based optimization algorithms (e.g., FIRE algorithm, as suggested by the `OptRandomWalk` in `FitREQ_lib.cpp`) are used to iteratively adjust the `REQH` parameters to minimize the RMSE. The system computes variational derivatives of the error with respect to the `REQH` parameters to guide the optimization.

### Hydrogen Bond Corrections

We model hydrogen bonds by two methods
1. we decrease repulsive of non-covalent interatomic pairwise potential ( i.e. Pauli repulsion, modeled by $A/r_{ij}^{12}$ in Lennard-Jones or by $A \exp(-2 a r_{ij})$ in Morse and Buckingham.
   * This is done by multiplying amplitude of repulsion $A$ by some constant, $A_H = A(1 + H_{ij})$ where $H_{ij} = H_i H_j$
   * This correction is in effect only for atoms which have opposite Hbond correction charges  $H_i, H_j$, which means $H_{ij}<0$.
2. we add electron pairs as dummy atoms at certain distance (~0.5A) from host Hbond-acceptor (e.g. electronegative O,N atoms). We transfer some charge (0.1-0.6e) from host atom to electron pair. 
   * The electron pair is typically closer to Hbond-donor 
   * This is supposed to simulate the angular dependence of hydrogen bonds
   
#### Only for selected atoms

* Fit REQ only for atoms which have non-zero Hb correction ( like `if(abs(Hbond)<1e-300) continue;` )
   * Or we may make a flag `bool* bHBs;   if(bHBs[ia]) continue`
   * Or we may make Hb atoms selection `int* iHbs;  for(int ia in iHbs){ ... };`
   * Or we may make new temp array which lists only hydrogen-bond corrected atoms
* For all of this we need do non-Hbonded reference
   * we compute reference Hbond energy first (using non-corrected LJQ or MorseQ) and store it in temporary array `atoms->userData->Emodel0`
   * to evaluate error we simply substract `dE = Eref - Emodel0 - Ecor` where Emodel0 is pre-calculated and Ecor is calculated just for selected Hbonded atoms.
 
#### Electron pairs - short-range 
   * Electron pairs have certain charge on itself ( which is subtracted from host atom)
   * However they can have also some short range function (e.g. $\exp(-a r_{ij})$ ) which is attractive for all Hbond-donors (electron depleted Hydrogen atoms). 
      * we can check this by `Hb[ja]>0` 


## 3. System Architecture and Class Design

The FitREQ system is structured across three main components:

### 3.1. C++ Core (`cpp/common/molecular/FitREQ.h`)

The `FitREQ` class is the central C++ component. It extends `MolecularDatabase` and manages the entire fitting process.

*   **Parameters (`REQH`):** Stores `Quat4d` parameters for each atom type. These are the degrees of freedom (DOFs) that are optimized.
*   **Degrees of Freedom (DOFs):**
    *   DOFs are defined *per atom type*, not per individual atom. This allows for a more generalized and transferable parameter set.
    *   The system maintains arrays for `DOFs`, `fDOFs` (forces/gradients), `vDOFs` (velocities), and `fDOFbounds` (bounds for parameters).
    *   A `dofSelection.dat` file (loaded via `loadDOFSelection`) specifies which `REQH` components for which atom types are active DOFs and thus subject to optimization. This provides great flexibility to fix certain parameters or optimize only a subset.
*   **Training Samples:** Stored as a `std::vector<Atoms*>` (pointers to `Atoms` objects). Each `Atoms` object represents a molecular configuration from the training set.
*   **`AddedData` Struct:** This is a crucial extension to the `Atoms` object for `FitREQ`. It stores additional information necessary for the fitting, particularly for handling dummy atoms and fragment separation:
    *   `bonds`, `directions`: Information related to bonding and orientation of dummy atoms.
    *   `host_atom_indices`: Links dummy atoms back to their parent "real" atoms.
    *   `HBn0`: Fragment separation index. This integer (`n0`) indicates the split point between two molecular fragments within a single `Atoms` object. This is essential for calculating inter-fragment non-covalent interactions.
*   **Optimization Control:**
    *   `bAddEpairs`: Flag to control the addition of electron pairs (dummy atoms).
    *   `bRegularization`: Enables regularization terms to prevent overfitting.
    *   `bCheckRepulsion`: Checks for overly repulsive samples.
    *   `bParallel`: Enables parallel computation (e.g., OpenMP) for gradient accumulation.
    *   `bFilterSamples`: Filters out problematic samples.
*   **Gradient Accumulation:** Supports parallel accumulation and reduction of gradients across multiple CPU cores.

### 3.2. C-Interface Library (`cpp/libs/Molecular/FitREQ_lib.cpp`)

This file provides a C-compatible API to the `FitREQ` C++ class, enabling interaction from other languages (like Python via `ctypes`). It uses a global `FitREQ` instance (`W`) and a global `OptRandomWalk` optimizer.

Key functions exposed:
*   **Configuration:** `setVerbosity`, `setGlobalParams`, `setup`, `setFilter`.
*   **Data Loading:** `loadTypes` (for `ElementTypes.dat`, `AtomTypes.dat`), `loadDOFSelection`, `loadWeights`, `loadXYZ`. The `loadXYZ` function is critical as it can optionally add electron pairs (`bAddEpairs` flag) and handle output XYZ files.
*   **Optimization:** `run` initiates the fitting process with configurable algorithm and parallelization options.
*   **Evaluation:** Functions to `getEs` (energies) and `getFs` (forces/gradients).
*   **Parameter Scanning:** `scanParam`, `scanParam2D` for exploring the energy landscape.
*   **Buffer Access:** Exposes internal data buffers (like `DOFs`, `typeREQs`) for direct memory access from Python, allowing efficient data transfer and manipulation.

### 3.3. Python Bindings (`pyBall/FitREQ.py`)

This Python module uses `ctypes` to create a user-friendly interface to the C-interface library. It wraps the C functions, handles argument type conversions, and integrates with NumPy for efficient array operations.

*   **`ctypes` Integration:** Defines `argtypes` and `restype` for C functions, ensuring correct data passing.
*   **NumPy Bridge:** Provides helper functions (`_np_as`, `getBuff`, `getIBuff`) to convert between NumPy arrays and C pointers, allowing Python to directly access and modify the internal C++ data buffers.
*   **Global Buffer Access (`getBuffs`):** After calling `getBuffs()`, key internal arrays (e.g., `DOFs`, `fDOFs`, `typeREQs`) become accessible as global NumPy arrays in Python, enabling interactive analysis and manipulation.
*   **Utility Functions:**
    *   `loadTypes`, `loadDOFSelection`, `loadWeights`, `loadXYZ`: Python wrappers for the C loading functions.
    *   `EnergyFromXYZ`: Parses energy data from XYZ files.
    *   `find_all_dirs`, `extract_fragment_names`, `combine_fragments`, `concatenate_xyz_files`: These are crucial for preparing complex training datasets, especially when dealing with multiple molecular fragments and organizing them into directories. They help automate the process of combining individual XYZ files into a single large training file for `loadXYZ`.
*   **Workflow:** The typical Python workflow involves loading types and DOF selections, preparing and loading training XYZ samples (potentially using the concatenation utilities), configuring the system, running the optimization, and then evaluating results.

## 4. Handling of Dummy Atoms and Molecular Topology

A significant aspect of FitREQ, particularly for non-covalent interactions, is the treatment of "dummy atoms." These are not physical atoms present in the initial `.xyz` files but represent important electronic features like free electron pairs (lone pairs) or sigma holes.

The `MMFFBuilder` (Molecular Mechanics Force Field Builder) library is used to manage molecular topology and integrate these dummy atoms.

### 4.1. `MMFFBuilder.h` and `MMFFBuilder.md`

The `MM::Builder` class is responsible for:
*   **Molecular Graph Construction:** Building the connectivity (bonds, angles, dihedrals) from atomic coordinates.
*   **Atomic Typing and Hybridization:** Determining the chemical environment of each atom (e.g., sp3, sp2, sp hybridization) which is crucial for assigning correct force field parameters and for identifying positions for dummy atoms.
*   **AtomConf Struct:** This struct within `MMFFBuilder` explicitly tracks npi (number of pi-orbitals) and ne (number of electron pairs) for each atom. This is where dummy atoms are conceptually represented in the topology.
*   **Adding Dummy Atoms:** `MMFFBuilder` provides methods to programmatically add these dummy atoms (electron pairs, capping hydrogens) to the molecular structure. These are not explicitly in the input `.xyz` files but are generated based on the inferred bonding environment.
    *   `addEpair()`: Method in `AtomConf` to add an electron pair.
    *   `addH()`: Method in `AtomConf` to add a capping hydrogen.
    *   `bAddEpairs` flag in `FitREQ` and `loadXYZ` in `FitREQ_lib.cpp` and `pyBall/FitREQ.py` indicates that this process is integrated into the sample loading.
*   **Fragment Separation (`splitByBond`):** The `MM::splitByBond` function is used to logically divide a molecule into two fragments by "breaking" a specified bond. This is essential for `FitREQ` to calculate inter-fragment interactions. The `HBn0` (fragment separation index) in `FitREQ`'s `AddedData` likely originates from this mechanism.

### 4.2. Importance for FitREQ

*   **Accurate Electrostatics:** Electron pairs and sigma holes contribute significantly to the electrostatic potential around molecules. By including them as "dummy atoms" with their own `REQH` parameters, `FitREQ` can achieve a more accurate representation of non-covalent interactions, especially hydrogen bonding and halogen bonding.
*   **Topology-Dependent Parameters:** The ability of `MMFFBuilder` to infer hybridization and bonding environment allows `FitREQ` to use atom-type-specific parameters that reflect the local chemical context.
*   **Automated Preparation:** The integration means that users don't need to manually define dummy atoms; they are generated automatically during the data loading process based on the molecular topology.

## 5. Degrees of Freedom and Parameter Selection

As mentioned, DOFs in FitREQ are per atom type. This design choice offers significant advantages:
*   **Generalizability:** Parameters optimized for a specific atom type (e.g., sp3 carbon, carbonyl oxygen) can be applied across many different molecules containing that atom type, promoting transferability of the force field.
*   **Reduced Complexity:** Instead of optimizing parameters for every single atom in a large system, the problem is reduced to optimizing a smaller set of parameters for a limited number of atom types.
*   **Flexible Selection:** The `dofSelection.dat` file allows users to precisely control which `REQH` components (R, E, Q, H) for which atom types are allowed to vary during optimization. This is crucial for:
    *   **Fixing known parameters:** If some parameters are well-established, they can be kept constant.
    *   **Targeted optimization:** Focusing the optimization on specific interaction types (e.g., only charges, or only vdW parameters).
    *   **Debugging:** Isolating the effect of certain parameters.

## 6. Input File Formats and User Workflow

The typical workflow for using FitREQ involves several input files:

### 6.1. `ElementTypes.dat` and `AtomTypes.dat`

These files define the mapping from element names and atom type names to internal integer IDs. `AtomTypes.dat` also specifies initial `REQH` parameters for each atom type.

### 6.2. `dofSelection.dat`

This file is critical for defining the degrees of freedom. It specifies which components of the `REQH` vector for each atom type are active DOFs. For example:

Here, for `C_sp3`, R, E, and Q are active, while H is fixed. For `O_carbonyl`, R is fixed, and E, Q, H are active.

### 6.3. Training Samples (`.xyz` files)

Training data consists of multiple `.xyz` files, each representing a molecular configuration. Each `.xyz` file should contain:
*   Atomic coordinates and element types.
*   A comment line (the second line in `.xyz` format) that typically includes the reference energy for that configuration.
*   For multi-fragment systems, the `concatenate_xyz_files` utility in `pyBall/FitREQ.py` can combine individual fragment `.xyz` files into a single training file, adding metadata to the comment lines to indicate fragment origins.

### 6.4. `weights.dat` (Optional)

Allows assigning different weights to individual training samples, giving more importance to certain configurations during the RMSE calculation.

### 6.5. Typical User Workflow (Python)

1.  **Prepare Input Files:**
    *   Ensure `ElementTypes.dat`, `AtomTypes.dat` are correctly defined.
    *   Create `dofSelection.dat` to specify active DOFs.
    *   Prepare training `.xyz` files. If fragments are in separate directories, use `pyBall.FitREQ.concatenate_xyz_files` to combine them.
2.  **Initialize and Load Data:**
    ```python
    import pyBall.FitREQ as fr
    fr.loadTypes()
    fr.loadDOFSelection()
    fr.loadXYZ("path/to/combined_training_data.xyz", bAddEpairs=True) # bAddEpairs is important for dummy atoms
    fr.getBuffs() # Access internal buffers as NumPy arrays
    ```
3.  **Configure Optimization:**
    ```python
    fr.setGlobalParams(dt=0.1, damp=0.5, ...)
    fr.setup(bRegularization=True, ...)
    # Access and modify fr.DOFs, fr.typeREQs directly if needed
    ```
4.  **Run Optimization:**
    ```python
    fr.run(nstep=1000, iparallel=1) # Run for 1000 steps, parallelized
    ```
5.  **Evaluate Results:**
    ```python
    energies = fr.getEs()
    # Analyze energies, compare with reference, plot convergence etc.
    ```
6.  **Save Parameters:** Extract optimized parameters from `fr.typeREQs` or `fr.DOFs`.

## 7. Developer Notes

### 7.1. Code Structure and Dependencies

*   **`FitREQ.h`:** Core logic, data structures, optimization loop. Depends on `MMFFBuilder` for molecular topology.
*   **`FitREQ_lib.cpp`:** C-interface, global instances, initialization, and direct function calls to `FitREQ` methods.
*   **`pyBall/FitREQ.py`:** Python `ctypes` wrapper, NumPy integration, and data preparation utilities.
*   **`MMFFBuilder`:** Provides fundamental molecular topology services, including handling of dummy atoms and fragment splitting.

### 7.2. Dummy Atom Implementation Details

*   The `ne` (electron pairs) and `nH` (capping hydrogens) fields in `MM::AtomConf` are key. When `bAddEpairs` is true during `loadXYZ`, `MMFFBuilder` analyzes the bonding environment of atoms and adds these "virtual" neighbors.
*   These dummy atoms are then assigned their own `REQH` parameters (or use default ones) and are included in the energy calculations, even though they don't correspond to explicit atoms in the original `.xyz` file.
*   The `AddedData` struct in `FitREQ` is used to store the specific information about these added features (e.g., their relative positions, host atom).

### 7.3. Parallelization

*   `FitREQ` supports parallel gradient accumulation (e.g., via OpenMP) to speed up the optimization process, especially for large training sets. The `iparallel` flag in `run()` controls this.
*   Care is taken to minimize synchronization overhead in parallel regions.

### 7.4. Future Development Considerations

*   The current OpenCL-based implementation (`FitREQ_ocl.md`) does not yet handle dummy atoms. The detailed documentation of dummy atom handling in this C++/Python version is crucial for guiding future re-implementation in the OpenCL version to ensure feature parity and performance.
*   The modular design (C++ core, C-interface, Python bindings) facilitates independent development and testing of each layer.

## 8. Practical Examples and Use Cases

FitREQ is particularly useful for:
*   **Developing custom force fields:** Tailoring non-covalent parameters for specific chemical systems or interaction types (e.g., metal-ligand interactions, specific functional groups).
*   **Refining existing force fields:** Improving the accuracy of vdW and electrostatic terms in established force fields for better agreement with high-level reference data.
*   **Studying intermolecular interactions:** Gaining insights into the nature of non-covalent forces by analyzing the optimized parameters.
*   **Teaching and Research:** Providing a flexible platform for students and researchers to explore force field development and molecular simulation.

This documentation aims to serve as a comprehensive guide for both users seeking to apply FitREQ to their problems and developers looking to extend or maintain the library.