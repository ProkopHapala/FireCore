
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

*   The `ne` (electron pairs) and `nH` (capping hydrogens) fields in `MM::AtomConf` are key. When `bAddEpairs` is true during `loadXYZ`, `MMFFBuilder` analyzes the bonding environment of atoms and adds these "virtual" neighbors. The `loadXYZ` function calls `addAndReorderEpairs` to manage the creation and reordering of these dummy atoms.
*   These dummy atoms are then assigned their own `REQH` parameters (or use default ones) and are included in the energy calculations, even though they don't correspond to explicit atoms in the original `.xyz` file.
*   The `AddedData` struct in `FitREQ` is used to store the specific information about these added features (e.g., their relative positions, host atom). Their positions are initialized using `initFittedAdata`.
*   In the end, a database is created in the form of a list of `Atoms` objects (from `Atoms.h`), with an added pointer to `void* userData = 0;`. The instances of the `Atoms` class representing each training sample are then stored in the `samples<Atoms*>` vector.
*   The `FitREQ` class extends the `MolecularDatabase` class (from `MolecularDatabase.h`) to efficiently work with this `userData`.

#### 7.2.1 Mapping dummy atoms to their hosts (critical)

Location of code:

* `MMFFBuilder::buildBondsAndEpairs()` in `cpp/common/molecular/MMFFBuilder.h`
* `MMFFBuilder::listEpairBonds(Vec2i*& bs, Vec3d*& dirs)` in `cpp/common/molecular/MMFFBuilder.h`
* `FitREQ::addAndReorderEpairs(Atoms*& atoms)` in `cpp/common/molecular/FitREQ.h`

Flow and data structures (final, correct behavior):

* __Build & detect__: After auto-bonding and adding dummy atoms, `buildBondsAndEpairs()` finalizes the molecular graph. Dummy atoms (electron pairs) are capping atoms (usually `iconf == -1`).
* __List attachments (`bs`, `dirs`)__:
  * `listEpairBonds(bs, dirs)` produces arrays of length `nep` with entries
    `bs[k] = { ia_host, je_old }` and `dirs[k]` the attachment direction (unit vector from host to e‑pair).
* __Reordering and indexing policy__ (`addAndReorderEpairs`):
  * Final atom ordering inside one sample is enforced as
    `[mol1 core (0..n0bak-1)] [mol1 epairs (n0bak..n0bak+nE0-1)] [mol2 core] [mol2 epairs]`.
  * `nE0` = number of epairs whose host is in fragment 1 (`ia_host < n0bak`).
  * We copy e‑pairs explicitly into the target slots; only those final slots are marked as epairs.
  * After placing mol1 epairs we set `atoms->n0 = n0bak + nE0`. Mol2 core is copied starting at `atoms->n0`, ensuring it never overwrites mol1 epairs.
* __Updated `bs` (new indices)__:
  * For mol1 e‑pairs, we remap to the final epair index `j = n0bak + iE0` and keep the host as the original index in fragment 1:
    `bs_new[idx0].x = ia_host;  bs_new[idx0].y = j`.
  * For mol2 e‑pairs, host indices must be shifted by `nE0` because mol2 core is moved after mol1 epairs:
    `bs_new[idx1].x = ia_host + nE0;  bs_new[idx1].y = j_final`.
* __Dense host array__:
  * We fill `AddedData.host` of size `atoms->natoms` with `-1` and set
    `host[ bs_new[i].y ] = bs_new[i].x`.
  * Core (non‑dummy) atoms retain `-1`.

Invariants (checked with debug prints):

* If an index `i` is an epair slot, `params->atomTypeNames[ atypes[i] ][0] == 'E'`.
* Every e‑pair index present in `bs` has a defined host (`host[i] >= 0`).
* No core atom has a host (`host[i] == -1`).

This mapping and ordering are used downstream for forces/energy attribution and are mirrored in the extended XYZ export.

#### 7.2.2 Extended XYZ export with dummy-atom annotations

When `FitREQ::loadXYZ(fname, bAddEpairs=true, bOutXYZ=true)` is used, FitREQ writes an auxiliary file `<fname>_Epairs.xyz` containing the processed system with dummy atoms. The format has been extended to make debugging and data exchange explicit:

Header:

```
<natoms>
# n0 <split_index> E_tot <energy> <rest_of_original_comment>
```

Body (one line per atom):

```
<atype_name>  <x>  <y>  <z>  <q>  <host>
```

Where:

* __atype_name__: Full atom type name from `MMFFparams::atomTypeNames[atype]` (not just element symbol), so you can differentiate dummy types (e.g., lone-pair type, sigma-hole type) from real atoms.
* __x, y, z__: Cartesian coordinates (Angstrom).
* __q__: Charge of the atom as stored in `Atoms::charge[i]` (0 if the charge array is not allocated).
* __host__: Index of the host atom for dummy atoms, taken from `AddedData.host[i]`. For real atoms, this is `-1`.

Notes:

* The `host` column encodes attachment for both electron pairs and sigma holes (any capping dummy) because the detection traverses all capping neighbors.
* The original `.xyz` second-line comment is preserved after the energy and `n0` metadata so that fragment/source info remains intact.
* If desired, a consumer can easily reconstruct the set of dummy atoms via `{ i | host[i] >= 0 }` and the mapping to hosts via `host[i]`.

#### 7.2.3 Sigma holes (E_h): typing and export

Recent updates introduce a dedicated sigma‑hole dummy type `E_h` and ensure it propagates through building and export:

* __Type definition__: Add `E_h` to `AtomTypes.dat` alongside `E` with appropriate parameters (see `tests/tFitREQ/data/AtomTypes.dat`).
* __Binding__: `MMFFBuilder::bindParams()` tries to resolve `E_h` via `params->atomTypeDict`. If not found, it logs a note and falls back to generic `E`.
* __Lazy init__: `MMFFBuilder::addSigmaHole()` now lazily initializes `E_h` from `params->atomTypeDict` if `bindParams()` was not called beforehand, guaranteeing `E_h` usage when available.
* __Export fidelity__: `FitREQ::writeSampleToFile()` preserves full dummy atom names. The optional one‑letter truncation is applied only to core atoms (`host<0`), so sigma holes appear explicitly as `E_h` in `<fname>_Epairs.xyz`.

This makes sigma holes distinguishable end‑to‑end and simplifies downstream analysis and visualization.

#### 7.2.4 Geometry robustness for electron pairs

To avoid overlapping dummy atoms and ill‑defined directions:

* __Pi‑geometry guard__: `MMFFBuilder::addEpairsByPi()` only uses pi‑based geometry when a valid pi system and direction exist; otherwise it falls back to non‑pi geometry.
* __Zero‑direction warning__: `MMFFBuilder::addEpair()` detects `|h|≈0` and emits a warning: dummy would be placed on host. This surfaces problematic cases early and prevents silent failures. A future improvement can add a final placement fallback when this occurs.

#### 7.2.5 Debug logging defaults

Verbose tracing remains available but is muted by default to reduce noise in normal runs. The following prints are commented out (search for `printf` in the listed functions to toggle):

* `MMFFBuilder::addEpair()` detailed placement line
* `MMFFBuilder::addEpairsByPi()` per‑atom vector dump
* `MMFFBuilder::addSigmaHole()` placement line (lazy `E_h` discovery remains noted in code, also commented by default)
* `MMFFBuilder::listEpairBonds()` per‑bond dump
* `FitREQ::addAndReorderEpairs()` remapping prints (sanity warnings preserved)
* `AddedData::fill_host()` per‑entry dump

Keep the sanity checks enabled (e.g., type checks for epair slots) to fail fast on inconsistencies.

### 7.3. Parallelization

*   `FitREQ` supports parallel gradient accumulation (e.g., via OpenMP) to speed up the optimization process, especially for large training sets. The `iparallel` flag in `run()` controls this.
*   Care is taken to minimize synchronization overhead in parallel regions.

### 7.4. Future Development Considerations

*   The current OpenCL-based implementation (`FitREQ_ocl.md`) does not yet handle dummy atoms. The detailed documentation of dummy atom handling in this C++/Python version is crucial for guiding future re-implementation in the OpenCL version to ensure feature parity and performance.
*   The modular design (C++ core, C-interface, Python bindings) facilitates independent development and testing of each layer.

### 7.5. DOF testing workflow and Python utilities

* __Script location__: `tests/tFitREQ/opt_check_derivs.py` drives DOF scans/derivative checks.
* __Typical flow__ (see the script):
  * `fit.loadTypes()` → `fit.loadDOFSelection(dof_fname)` → `fit.loadXYZ(...)` → `fit.setup(...)` → `fit.setWeights(...)` → `fit.getBuffs()`.
  * For scans, the script iterates DOFs and calls `fit.plotDOFscan_one(...)` (or `fit.plotDOFscans(...)`).
* __Reload DOFs before each scan__ (scan independence): the script explicitly calls `fit.loadDOFSelection(dof_fname)` inside the per‑DOF loop to restore all DOFs to their initial selection values before every scan. This avoids subtle state coupling between scans even though `plotDOFscan_one()` backs up and restores the scanned DOF locally.
  * Reference: `opt_check_derivs.py` lines where `fit.loadDOFSelection(dof_fname)` is called right before `fit.plotDOFscan_one(...)`.
* __Python helpers__ in `pyBall/FitREQ.py`:
  * `loadDOFnames(fname, comps="REQH", return_specs=False)`: parses DOF selection and, if requested, returns per‑DOF specs (min/max, regularization). Useful to set scan ranges/labels.
  * `plotDOFscan_one(iDOF, ...)`: scans one DOF, plots energy and optionally force; internally backs up/restores the single DOF value.
  * `plotDOFscans(iDOFs, xs, DOFnames, ...)`: multi‑DOF overlay using `plotDOFscan_one()`.

### 7.6. Mixing rules and internal sqrt(E) convention

FitREQ uses Lorentz–Berthelot‑style mixing, however, we store Ei=sqrt(Eii) to avoid extra sqrt at runtime.

* __Internal storage__ (`FitREQ::initAllTypes()` in `cpp/common/molecular/FitREQ.h`):
  * `typeREQs[i] = ( R_i, E_i=sqrt(E_ii), Q_i, H_i )` where `R_i` is vdW radius, `E_ii` is vdW well-depth, `Q_i = base charge`, `H_i = H‑bond flag/strength`.
  * Code reference: `initAllTypes()` constructs `Quat4d{ RvdW, sqrt(EvdW), Qbase, Hb }`.
* __Pairwise mixing__ in evaluators (`evalExampleDerivs_*` in `FitREQ.h`):
  * Radii: `R_ij = R_i + R_j` (arithmetic sum).
  * At runtime, we compute `E_ij = E_i * E_j` without an extra sqrt ( because we store sqrt of well-depth E_i=sqrt(E_ii) not well depth E_ii ).
  * Charges: electrostatics use `Qij = q_i q_j` 
  * H‑correction: `H_ij = h_i h_j` is applied only if `H_ij<0` (via a sign gate `sH`); for epairs, vdW is replaced by a short‑range term.
* __Why sqrt(E)__: Pre‑storing `sqrt(E)` makes Lorentz–Berthelot `epsilon_ij = sqrt( epsilon_i epsilon_j )` a single multiply per pair (`E0 = y_i y_j`). This reduces cost and simplifies derivatives (e.g., `∂E/∂y_i ∝ y_j`). The appendix formulas use `E_0 = ε_i ε_j` with `ε_i ≡ sqrt(E_i)` to match the implementation.
* __Electron pairs (SR replacement)__: If either partner is a dummy epair (`host[idx]>=0`), vdW is replaced by a short‑range attraction `getSR2()` (Gaussian). The SR width `w` is taken from the epair partner’s radius component (`R`), and the derivative `∂E/∂w` is accumulated to that epair’s `R` DOF.


## Appendix: REQH energy models and variational derivatives

This appendix summarizes the exact energy expressions and variational derivatives implemented in `cpp/common/molecular/FitREQ.h` for the example evaluators `evalExampleDerivs_*`. It is intended to enable reimplementation in other languages.

Notation (per atom i):

$$\mathrm{REQH}_i=(R_i,\,\varepsilon_i,\,q_i,\,h_i),\quad E_0=\varepsilon_i\varepsilon_j$$

For a pair i–j:

$$R_0=R_i+R_j,\quad H_{ij}=\begin{cases}h_i h_j,& h_i h_j<0\\0,& \text{otherwise}\end{cases},\quad r=\lVert\mathbf r_j-\mathbf r_i\rVert,\quad u=\frac{R_0}{r}$$

Coulomb: $k=\mathrm{COULOMB\_CONST}$ and

$$E_{el}=\frac{k\,q_i q_j}{r}$$

Unless noted, $E_{ij}=E_{el}+E_{vdW}$.

1) LJQH2 (12–6, H gates the 12-term):

$$E_{vdW}=E_0\left[ (1+H_{ij})\,u^{12}-2\,u^{6}\right]$$

Derivatives (intermediates in code):

$$\frac{\partial E}{\partial E_0}=u^{6}\big[(1+H_{ij})u^{6}-2\big],\quad \frac{\partial E}{\partial R_0}=12\,\frac{E_0}{R_0}\,u^{6}\big[(1+H_{ij})u^{6}-1\big],\quad \frac{\partial E}{\partial H}=-E_0 u^{12}$$

$$\frac{\partial E}{\partial q_i}=\frac{k q_j}{r},\quad \frac{\partial E}{\partial q_j}=\frac{k q_i}{r}$$

Per-atom accumulation (`evalExampleDerivs_LJQH2_SR`):

$$\begin{aligned}
\frac{\partial E}{\partial R_i}&=-\frac{\partial E}{\partial R_0},\quad &\frac{\partial E}{\partial \varepsilon_i}&=-\tfrac{1}{2}\,\varepsilon_j\,\frac{\partial E}{\partial E_0},\\
\frac{\partial E}{\partial q_i}&=-\frac{k q_j}{r},\quad &\frac{\partial E}{\partial h_i}&=\;h_j\,s_H\,\frac{\partial E}{\partial H}
\end{aligned}$$

Symmetric contributions are added for atom j; see `evalExampleDerivs_LJQH2_SR()` for exact accumulation.

2) LJr8QH2 (8–6, H gates the r^8 term):

$$E_{vdW}=E_0\left[3(1+H_{ij})u^{8}-4u^{6}\right]$$

Let $u^2\equiv u_2$ and $u_2'=(1+H_{ij})u_2$. Then

$$\frac{\partial E}{\partial E_0}=u^{6}(3u_2'-4),\quad \frac{\partial E}{\partial R_0}=24\,\frac{E_0}{R_0}u^{6}(u_2'-1),\quad \frac{\partial E}{\partial H}=-3E_0 u^{8}$$

Coulomb and per-atom accumulation as in LJQH2.

3) MorseQ_SR (Morse, H gates the $e^{-2\alpha\Delta}$ term):

Let $\Delta=r-R_0$, $e=\exp(-\alpha\Delta)$, $e^2=e^2$.

$$E_{vdW}=E_0\big[(1+H_{ij})e^{2}-2e\big]$$

$$\frac{\partial E}{\partial E_0}=(1+H_{ij})e^{2}-2e,\quad \frac{\partial E}{\partial R_0}=2\alpha E_0\big[(1+H_{ij})e^{2}-e\big],\quad \frac{\partial E}{\partial H}=-E_0 e^{2}$$

Coulomb and per-atom accumulation as in LJQH2.

4) Short-range SR terms for electron pairs (replace vdW if i or j is an epair):

If `host[idx] >= 0` for either partner, vdW is replaced by SR:

- Gaussian form used in `*_SR` (code `getSR2`), with $w$ taken from the epair partner’s radius:

  $$E_{SR}=H_{ij}\,e^{-(r/w)^2},\quad \frac{\partial E}{\partial H}=e^{-(r/w)^2},\quad \frac{\partial E}{\partial w}=2H_{ij}\,e^{-(r/w)^2}\,\frac{r^{2}}{w^{3}}$$

- Exponential form in older `_uncorr` (code `getSR`):

  $$E_{SR}=H_{ij}\,e^{-r/w}$$

Accumulation: if i is an epair, `w = R_i` and the R-derivative is added to i; if j is an epair, `w = R_j` and the R-derivative is added to j. The code flips signs to match its negative-gradient convention (see `evalExampleDerivs_*_SR()`).

5) Electron-pair charge conservation (host balancing):

Charges for an epair are taken from its host atom, enforcing `q_host + q_epair = const`. After per-atom accumulation into `fREQs`, `acumHostDerivs()` subtracts the host’s charge derivative from the epair DOF:

- For each epair–host link `(host=h, epair=e)`: the host’s `fREQs[h].z` is subtracted into the DOF index of the epair type, ensuring the constraint in gradient space.

5b) Mapping per-atom derivatives to DOFs:

Per-atom gradients `dEdREQs[i] = (∂E/∂R_i, ∂E/∂ε_i, ∂E/∂q_i, ∂E/∂h_i)` are reduced to the global DOF vector via the type map `typToREQ` in `acumDerivs()`:

* For atom i of type t, each component k∈{x,y,z,w} is added into DOF index `typToREQ[t].k` (if non-negative).
* After this, `acumHostDerivs()` applies the charge-balance correction described above to epair/host links.

5c) Sign convention in code vs formulas here:

The code accumulates negative energy derivatives into `fREQs` (e.g., `fREQi.x += -dE_dR0`). All formulas in this appendix are written as plain partial derivatives `∂E/∂(·)`. When comparing to code, account for the minus sign used to form a descent direction.

6) corr vs uncorr evaluators:

- `_uncorr` (e.g., `evalExampleDerivs_LJQH2_uncorr`): computes energies only for pairs where neither atom is fitted; skips derivatives and fitted atoms entirely. Used to keep a static background energy.
- `_corr` / `_SR` (e.g., `evalExampleDerivs_LJQH2_corr`, `*_SR`): computes energies and accumulates variational derivatives for selected/fitted atoms; includes SR treatment for epairs and gates H via `s_H`. The flag `bCheckAddJ` prevents double counting when both partners are fitted by zeroing `Eij` on the second pass.

Implementation references: `getSR`/`getSR2`, `evalExampleDerivs_LJQH2_SR`, `evalExampleDerivs_LJr8QH2_SR`, `evalExampleDerivs_MorseQ_SR`, and host balancing in `acumHostDerivs`.

## 8. Practical Examples and Use Cases

FitREQ is particularly useful for:
*   **Developing custom force fields:** Tailoring non-covalent parameters for specific chemical systems or interaction types (e.g., metal-ligand interactions, specific functional groups).
*   **Refining existing force fields:** Improving the accuracy of vdW and electrostatic terms in established force fields for better agreement with high-level reference data.
*   **Studying intermolecular interactions:** Gaining insights into the nature of non-covalent forces by analyzing the optimized parameters.
*   **Teaching and Research:** Providing a flexible platform for students and researchers to explore force field development and molecular simulation.

This documentation aims to serve as a comprehensive guide for both users seeking to apply FitREQ to their problems and developers looking to extend or maintain the library.

---

## Development / ToDo

### Exporting Pre-processed Training Samples with Dummy Atoms

**Goal:** Create a function to export the `Atoms` objects (training samples) from `FitREQ`, including the dummy atoms added during `loadXYZ` and `addAndReorderEpairs`, into a standard `.xyz` file format. This will facilitate testing and comparison with the OpenCL implementation.

**Location:** The new function should ideally be added to `FitREQ.h` as a method of the `FitREQ` class, as it directly interacts with the internal `samples` vector and `AddedData` structures. A corresponding C-interface function can be added to `FitREQ_lib.cpp` for Python binding.

**Function Signature (proposed for `FitREQ.h`):**

```cpp
// In FitREQ.h
void exportSampleToXYZ( int iSample, const char* filename );
```

**Detailed Steps:**

1.  **Understand `Atoms` and `AddedData` Structure:**
    *   Recall that each training sample is an `Atoms` object within the `samples` vector.
    *   Each `Atoms` object has a `void* userData` pointer, which `FitREQ` uses to store an instance of `AddedData`.
    *   The `AddedData` struct contains information about the dummy atoms (e.g., `nAdded`, `addedAtomType`, `addedAtomPos`, `addedAtomCharge`).
    *   The `Atoms` object itself contains the original atom positions (`pos`), types (`atypes`), and charges (`qs`).

2.  **Accessing Sample Data:**
    *   Inside `exportSampleToXYZ`, retrieve the `iSample`-th `Atoms` object from the `samples` vector.
    *   Cast the `userData` pointer back to `AddedData*` to access dummy atom information.

3.  **Constructing Output Data:**
    *   For the specified `iSample`:
        *   Get the number of original atoms (`samples[iSample]->natoms`).
        *   Get the number of added dummy atoms (`((AddedData*)samples[iSample]->userData)->nAdded`).
        *   The total number of atoms to write will be `samples[iSample]->natoms + ((AddedData*)samples[iSample]->userData)->nAdded`.
        *   Iterate through the original atoms: write their type, x, y, z coordinates, and charge.
        *   Iterate through the dummy atoms: write their type (e.g., "X" or "EP" for electron pair, "H" for capping hydrogen, based on `addedAtomType`), x, y, z coordinates, and charge. We need to decide on a naming convention for dummy atom types (e.g., `EP` for electron pairs, `CH` for capping hydrogens).

4.  **File Writing:**
    *   Open the `filename` in write mode.
    *   Write the total number of atoms (original + dummy) on the first line.
    *   Write a comment line (e.g., including sample index, energy if available, and any other relevant metadata from the `FitREQ` sample). Try to replicate the comment line format from `input_example.xyz`: `# n0 3 Etot 3.46679898796992282901 x0 01.40 z -90 H2O-D1_H2O-A1`.
    *   For each atom (original and dummy): write `AtomType X Y Z Charge` on a new line, separated by spaces. Ensure proper floating-point precision for coordinates and charges.

5.  **Error Handling:**
    *   Check if the `filename` can be opened for writing.
    *   Check if `iSample` is within valid bounds.
    *   Handle cases where `userData` might be null or not of type `AddedData*`.

**Testing:**
*   Create a simple test case in Python that:
    1.  Loads an XYZ file into `FitREQ`.
    2.  Calls the new `exportSampleToXYZ` function.
    3.  Compares the generated XYZ file with a reference file or visually inspects it.

**Integration with `FitREQ_lib.cpp` and `pyBall/FitREQ.py`:**
*   Add a C-interface function in `FitREQ_lib.cpp` (e.g., `_exportSampleToXYZ`) that calls the `FitREQ::exportSampleToXYZ` method.
*   Add a corresponding Python wrapper function in `pyBall/FitREQ.py` that calls the C-interface function.


---

# DEBUGGING (Resolved):

We observed incorrect epair ordering and bogus host assignments in `tests/tFitREQ/input_example_1.xyz_Epairs.xyz` caused by overwriting of mol1 e‑pair slots when copying mol2 core.

__Root cause__
* `addAndReorderEpairs()` incremented `atoms->n0` while placing mol1 e‑pairs, then copied mol2 core starting at the evolving `n0`, which overwrote some epair slots. `isep` flags and `bs` became inconsistent with `atypes`.

__Fix implemented__ (in `cpp/common/molecular/FitREQ.h`):
* Initialize `isep` to zeros and only mark final epair slots.
* Do not change `n0` during epair placement; after placing mol1 e‑pairs, set `atoms->n0 = n0bak + nE0` once.
* Copy mol2 core starting at this final `atoms->n0`.
* Keep `bs` consistent: mol2 host indices are shifted by `nE0`.
* Add invariant checks warning if any `isep[i]==1` is non‑`E`.

__Current guarantees__
* Final layout: `[mol1 core][mol1 E][mol2 core][mol2 E]`.
* Every `E` line has a valid host (`>=0`); core atoms print host `-1`.
* `AddedData.host` is filled from `bs` and is authoritative for export/debug.