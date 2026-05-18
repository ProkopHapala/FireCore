# Documentation of OpenCL-accelerated Non-Covalent Interaction Fitting Program (`NonBondFitting.py`, `FitREQ.cl`)

## Part 1: Conceptual Overview

### 1.1 Introduction

This program is a high-performance tool designed to optimize the parameters of a classical force field for describing non-covalent interactions (NCIs). NCIs, such as hydrogen bonds and van der Waals forces, are crucial for accurately modeling molecular systems in chemistry and biology.

The core purpose of this software is to adjust a set of parameters in our force field model so that the energies it calculates for various molecular arrangements match the "ground truth" energies obtained from high-accuracy (but computationally expensive) quantum mechanical (QM) calculations.

The program is built using Python for overall control and OpenCL for massively parallel computation on Graphics Processing Units (GPUs), enabling the rapid processing of thousands of training examples required for robust parameter fitting.

### 1.2 Core Physical and Mathematical Principles

#### 1.2.1 The Objective Function: Minimizing Squared Error

The fundamental goal of the fitting process is to minimize the difference between our model's predictions and the reference data. We quantify this difference using a **Root Mean Square Error (RMSE)**-based objective function, which we aim to minimize. The function `L` is defined as the weighted sum of squared errors over all training samples:

`L = Σᵢ wᵢ * (E_model,i - E_ref,i)²`

Where:
*   `i` is the index of a training sample (a specific geometry of two interacting molecular fragments).
*   `E_ref,i` is the reference interaction energy for sample `i` from a QM calculation.
*   `E_model,i` is the interaction energy for sample `i` calculated by our force field model using the current set of parameters.
*   `wᵢ` is a weight factor for sample `i`, allowing us to prioritize the fitting of certain interactions.

The optimization process is essentially a search for the set of parameters that makes the value of `L` as close to zero as possible.

#### 1.2.2 The Force Field Model: REQH Parameters

Our model describes the interaction energy between two atoms as a sum of a Lennard-Jones (LJ) potential and a Coulomb (electrostatic) potential. The behavior of each atom type is defined by four primary parameters, which are the **Degrees of Freedom (DOFs)** in our fitting process:

*   **R (RvdW):** The van der Waals radius of the atom type. It primarily influences the position of the potential energy minimum.
*   **E (EvdW):** The van der Waals well-depth. It determines the strength of the attractive LJ interaction.
*   **Q (Charge):** The partial atomic charge. It governs the strength of the electrostatic interaction.
*   **H (H-bond):** A correction term applied to the LJ potential, typically used to fine-tune the behavior of hydrogen bonds.

The interaction energy `E_pair` between two atoms, `a` and `b`, is given by:

`E_pair = E_LJ + E_el`

`E_LJ = E₀ * [ (1 + H₂) * (R₀/r)¹² - 2 * (R₀/r)⁶ ]`
`E_el = k * (Qₐ * Qᵦ) / r`

Where `r` is the distance between the atoms, and `R₀`, `E₀`, `H₂`, `Qₐ`, and `Qᵦ` are derived from the individual atomic REQH parameters using standard mixing rules (e.g., `R₀ = Rₐ + Rᵦ`, `E₀ = sqrt(Eₐ * Eᵦ)`).

#### 1.2.3 The Optimization Process: Gradient-Based Minimization

To efficiently find the minimum of our objective function `L`, we use a gradient-based optimizer (like FIRE or L-BFGS). These methods require knowing the "force" on each parameter—that is, the negative derivative (gradient) of `L` with respect to that parameter. The force tells us in which direction to adjust the parameter to lower the error.

We use the **chain rule** of calculus to find this force. For a single parameter `p` (e.g., the `RvdW` of an oxygen atom):

`F_p = -∂L/∂p = - Σᵢ [ ∂L/∂E_model,i * ∂E_model,i/∂p ]`

Breaking this down:

1.  **`∂L/∂E_model,i`**: This is the derivative of the error term with respect to the model energy. It is straightforward to calculate: `2 * wᵢ * (E_model,i - E_ref,i)`. This term tells us how much a change in the model energy affects the total error. It is large when our model is very wrong.

2.  **`∂E_model,i/∂p`**: This is the derivative of the model energy with respect to the parameter `p`. This term is calculated by our OpenCL kernels by summing the analytic derivatives of the LJ and Coulomb equations.

The core of our OpenCL program is to efficiently compute `∂E_model,i/∂p` for all atoms and all samples, multiply it by the error term `-2 * wᵢ * (E_model,i - E_ref,i)`, and then assemble the final forces for each DOF.

---

## Part 2: User Guide

This section describes how to prepare the necessary input files and run the program.

### 2.1 Overview of Input Files

To run a fitting, you must provide three text files:

1.  **`AtomTypes.dat`**: A library of default chemical and physical parameters for all atom types known to the system.
2.  **`input.xyz`**: A file containing the geometries, charges, and reference energies for all the training samples.
3.  **`dofSelection.dat`**: The main control file where you specify which REQH parameters you want to optimize and how to constrain them.

### 2.2 Input File Formats

#### 2.2.1 `AtomTypes.dat`: The Base Parameter Library

This file acts as a database of default parameters. The program uses these values as a starting point before applying any overrides from `dofSelection.dat`.

*   **Purpose:** To provide a robust set of default `RvdW` and `EvdW` values for all atom types.
*   **Format:** A space-separated text file. Lines starting with `#` are ignored. The program primarily uses columns 1, 10, and 11.

| Column | Name      | Description                                     | Used? |
| :----: | :-------- | :---------------------------------------------- | :---: |
| 1      | **name**  | **The unique name of the atom type (e.g., `O_3`, `H_OH`)** | **Yes** |
| 2-9    | ...       | Other chemical properties (not used by this program) | No    |
| 10     | **RvdW**  | **The default van der Waals radius (Å)**            | **Yes** |
| 11     | **EvdW**  | **The default van der Waals well-depth (kcal/mol)** | **Yes** |
| 12+    | ...       | Other parameters (not used)                     | No    |

**Example:**
```
#name  ...      RvdW        EvdW
C_3    ...      1.9255    0.0045532
O_3    ...      1.7500    0.0026018
H_OH   ...      1.4870    0.0006808
```

#### 2.2.2 `input.xyz`: The Training Data

This file contains all the reference molecular geometries and their corresponding energies that the model will be trained on. It is a concatenation of many standard XYZ format blocks.

*   **Purpose:** To provide the ground truth data for the fitting.
*   **Format:** Each sample entry has the following structure:
    1.  **Line 1:** The total number of atoms (`na`) in the sample (both fragments combined).
    2.  **Line 2 (Comment):** A line starting with `#` containing critical information.
        *   `# n0 [int] ... Etot [float] ...`
        *   `n0`: The number of atoms in the *first* fragment. This is used to split the system into two interacting molecules.
        *   `Etot`: The reference QM interaction energy for this geometry.
    3.  **Lines 3 to `na`+2:** The atom data lines.
        *   `atom_type   x   y   z   charge`
            *   `atom_type`: The name of the atom type, which must correspond to names in `AtomTypes.dat` and `dofSelection.dat`.
            *   `x, y, z`: Cartesian coordinates of the atom.
            *   `charge`: The initial partial charge of the atom. This value is used by the `Q` component of the REQH model. The charge parameter Q can come from two sources: the initial partial charge from the XYZ file (`atoms.w`) or from the type-based parameter matrix (`tREQHs[:,2]`). The choice is controlled by a runtime flag `use_type_charges` in the `FittingDriver` constructor. If `use_type_charges=True`, charges are taken from `tREQHs[:,2]`; otherwise from `atoms.w`. This allows switching between per-atom and type-based charges without recompilation.

**Example (`input.xyz` containing two samples):**
```
6
# n0 3 Etot 3.4667 ...
O_3 -0.968   0.555   0.000  -0.679
H_O -0.000   0.000   0.000   0.339
H_O -1.182   0.944   0.000   0.339
O_3  0.857  -1.399   0.000  -0.679
H_O -0.756  -2.004   0.000   0.339
H_O  0.756  -2.004   0.000   0.339
8
# n0 5 Etot -0.0374 ...
O_2 -2.573   1.559   0.000  -0.394
C_2 -1.410   1.280   0.000   0.488
... (6 more atoms) ...
```

#### 2.2.3 `dofSelection.dat`: The Optimization Control File

This is the most important file for controlling the fitting process. It specifies exactly which parameters to optimize and defines their behavior during the optimization.

*   **Purpose:** To select DOFs, set their initial values, and define regularization potentials to constrain them.
*   **Format:** A space-separated text file. Lines starting with `#` are ignored.

| Column | Name      | Description                                                                 |
| :----: | :-------- | :-------------------------------------------------------------------------- |
| 1      | **typename**| The name of the atom type to which this DOF belongs (e.g., `O_3`).            |
| 2      | **comp**    | The component of the REQH vector to optimize: `0`=R, `1`=E, `2`=Q, `3`=H.  |
| 3      | **Min**     | The hard lower bound for this parameter. The value will be clamped.        |
| 4      | **Max**     | The hard upper bound for this parameter.                                   |
| 5      | **xlo**     | The lower boundary of the soft-wall potential.                              |
| 6      | **xhi**     | The upper boundary of the soft-wall potential.                              |
| 7      | **Klo**     | Stiffness of the lower soft-wall (harmonic potential).                      |
| 8      | **Khi**     | Stiffness of the upper soft-wall.                                           |
| 9      | **K0**      | Stiffness of a harmonic potential tethering the parameter to its `xstart`. |
| 10     | **xstart**  | The **initial value** for this parameter at the start of the optimization. |
| 11     | **invMass** | (Legacy) Inverse mass for MD-based optimizers. Used by FIRE.               |

**Example:**
```
# typename comp   Min    Max      xlo    xhi      Klo   Khi   K0     xstart   invMass
O_3        1    0.001   0.01    0.002   0.008    1.0   1.0   0.1    0.0026   1.0
H_OH       0    1.2     1.8     1.4     1.6      1.0   1.0   0.0    1.487    1.0
```
This example defines two DOFs:
1.  The `EvdW` parameter (`comp=1`) for atom type `O_3` will be optimized, starting from a value of `0.0026`.
2.  The `RvdW` parameter (`comp=0`) for atom type `H_OH` will be optimized, starting from `1.487`.

---

## Part 3: Developer and Implementation Guide

This section describes the internal architecture of the program, focusing on the OpenCL implementation for performance.

### 3.1 High-Performance Architecture

The fitting problem is "pleasantly parallel," as the error and derivative for each of the thousands of training samples can be calculated independently. This makes it a perfect candidate for GPU acceleration.

The core design principle is the **Work-Group-Per-Sample Strategy**. Instead of launching one GPU thread for every atom in the entire dataset, we launch one **work-group** (e.g., 32 threads) for each **training sample**. All threads in a work-group cooperate on the same sample, which has two major advantages:
1.  **Data Locality:** All data for a sample can be loaded and processed, maximizing cache hits.
2.  **Efficient Reductions:** As shown below, this allows for the summation of values (like energy) across a sample to be done using ultra-fast `__local` memory and parallel reduction algorithms, completely avoiding slow global atomic operations.

### 3.2 OpenCL Kernel Breakdown (current implementation)

We use three kernels that separate concerns cleanly and allow a single-kernel path for derivative testing/optimization.

1)  **`evalSampleDerivatives_template`** (per-sample evaluation + error scaling)

   - __Role__: For each sample (one work-group), compute the sample energy Emol and the per-atom REQH derivatives of Emol. Compute the error scale `LdE = W*(Emol - Eref)`. Write:
     - Scaled per-atom derivatives into `dEdREQs`.
     - Per-sample objective contribution into `Jmols[iS] = 0.5*(Emol - Eref)*LdE`.

   - __Inputs__:
     - `ranges [nSamples] (int4)`: `(i0, j0, ni, nj)` indices into `atoms`/`atypes` for fragment-1 and fragment-2 of sample iS.
     - `tREQHs [nTypes] (float4)`: REQH parameters; column 1 stores `sqrt(E)`.
     - `atypes [nAtomTot] (int)`: type index per atom.
     - `ieps   [nAtomTot] (int2)`: pair indices to subtract charges (electron-pair handling). If `>=0`, subtract corresponding `tREQHs[iep].z` from `Q`.
     - `atoms  [nAtomTot] (float4)`: `(x,y,z,q)` positions + charges; `q` is used as the actual atomic charge.
     - `ErefW  [nSamples] (float2)`: `(Eref, W)` per sample.
     - `globParams (float4)`: currently `{alphaMorse,0,0,0}`.

   - __Outputs__:
     - `dEdREQs [nAtomTot] (float4)`: for each atom `ia` in the sample, stores the vector `∂Emol/∂REQH(ia)` scaled by `LdE`. In code: `dEdREQs[i0 + iL] = fREQi * LdE` for lanes `iL<ni` (fragment-1).  See notes below on coverage of both fragments.
     - `Jmols   [nSamples] (float)`: per-sample objective term `0.5 * W * (Emol - Eref)^2`.

   - __Work sizes__: one work-group per sample; typical `local_size = 32`.

   - __Model injection__: uses a macro marker `//<<<MODEL_PAIR_ACCUMULATION` where model-specific pair energy and its REQH-derivatives are accumulated.

   - __Important semantics__:
     - `dEdREQs` written by this kernel are already scaled by the error factor `LdE`. This makes assembly a pure linear gather/reduction.
     - Intended coverage is per-atom for the whole sample (both fragments). The present implementation writes for lanes `iL<ni` (fragment-1 atoms). If DOFs involve atoms from fragment-2, either (a) also populate those entries, or (b) ensure `DOFtoAtom` excludes fragment-2 atoms for the corresponding DOFs. Mismatch here can lead to near-zero assembled derivatives.

2)  **`evalSampleEnergy_template`** (energy-only)

   - __Role__: For each sample (one work-group), compute the sample energy `Emol` only. Used by visualization or energy-only scoring.

   - __Inputs__: same structural inputs as above; uses model code injected at `//<<<MODEL_PAIR_ENERGY`.

   - __Output__:
     - `Emols [nSamples] (float)`: per-sample energies.

3)  **`assembleAndRegularize`** (gather per-atom derivatives into DOF forces; add regularization)

   - __Role__: Each work-group handles one DOF. It reduces over the atoms contributing to that DOF using the mapping arrays and produces the total force on the DOF. Then it adds regularization forces based on per-DOF constraints and current DOF values.

   - __Inputs__:
     - `nDOFs (int)`
     - `DOFnis [nDOFs] (int2)`: for DOF k, `DOFnis[k] = (i0, ni)` indexes its block in `DOFtoAtom`/`DOFcofefs`.
     - `DOFtoAtom [nInds] (int)`: flattened list of global atom indices `ia` that contribute to each DOF; blocks are laid out back-to-back.
     - `DOFcofefs [nInds] (float4)`: coefficients to linearly map an atom's REQH derivative vector into the scalar contribution for this DOF. The kernel uses `dot(DOFcofefs[j], dEdREQs[ia])`.
     - `dEdREQs   [nAtomTot] (float4)`: scaled per-atom derivatives written by `evalSampleDerivatives_template`.
     - `DOFs      [nDOFs] (float)`: current DOF values (needed for regularization only).
     - `regParams [nDOFs] (float8)`: `{min,max, xlo,xhi, Klo,Khi, K0,x0}` per DOF.

   - __Output__:
     - `fDOFs [nDOFs] (float)`: total forces per DOF (physical assembly + regularization).

   - __Work sizes__: one work-group per DOF; typical `local_size = NLOC_assembleDOFderivatives` (e.g., 128).

#### 3.2.1 Objective composition and buffers

- __Per-sample objective__: `Jmols[iS] = 0.5 * W * (Emol - Eref)^2` written by `evalSampleDerivatives_template`.
- __Global objective__: host sums `J_phys = Σ_s Jmols[s]`. Regularization energy `J_reg` is computed on host (clamped to `[min,max]` then soft-walls and tether) and the total objective used in optimization is `J = J_phys - J_reg` to match the sign convention of the kernel that subtracts regularization gradients.
- __Scaled per-atom derivatives__: `dEdREQs[ia]` holds `LdE * ∂Emol/∂REQH(ia)`; assembly uses a pure sum of linear projections.

#### 3.2.2 DOF indexing: `DOFnis`, `DOFtoAtom`, `DOFcofefs` (must be consistent)

- __Global atom indices__: The atom indices stored in `DOFtoAtom` are global indices into `atoms`/`atypes`/`dEdREQs` (not sample-local). They refer to the same `ia` indices used by the kernels.

- __Block structure__:
  - For DOF k: `DOFnis[k].x = i0` (start in flattened arrays), `DOFnis[k].y = ni` (number of contributing atoms).
  - The contributing atoms for k are `ia_j = DOFtoAtom[i0 + j]` for `j=0..ni-1`.
  - The mapping coefficients for these atoms are `cof_j = DOFcofefs[i0 + j]`.

- __Assembly formula__ (what the kernel computes):
  - Physical part per DOF k: `f_phys[k] = Σ_{j=0..ni-1} dot( DOFcofefs[i0+j], dEdREQs[ DOFtoAtom[i0+j] ] )`.
  - Then regularization contributions are added inside the kernel based on `DOFs[k]` and `regParams[k]`.

- __Construction rule-of-thumb__:
  - For a DOF that tweaks component `c ∈ {R,E,Q,H}` of a given atom type T, include all atoms of type T in `DOFtoAtom` and set `DOFcofefs` to select the `c`-component, optionally premultiplied by the chain-rule factors coming from mixing rules (e.g., derivatives of `R0`, `E0`, `Q` with respect to the underlying per-atom REQH values). Keep signs consistent with model definitions.

- __Consistency requirement (important to avoid zero forces)__:
  - For every `ia` present in `DOFtoAtom`, the producer kernel must have written a valid `dEdREQs[ia]` for the current launch. If the derivative kernel writes only fragment-1 atoms for a sample, do not include fragment-2 atoms in `DOFtoAtom` (or extend the kernel to also populate those entries). Otherwise, assembled forces for such DOFs can be spuriously zero.

### 3.3 The Python Driver (`FittingDriver`)

The Python class acts as the high-level orchestrator for the entire process. Its key responsibilities include:

*   **File I/O:** Loading and parsing the three input files (`AtomTypes.dat`, `input.xyz`, `dofSelection.dat`).
*   **Data Preparation:**
    *   Building the initial parameter matrix (`tREQHs_base`) by combining default values from `AtomTypes.dat` with `xstart` values from `dofSelection.dat`.
    *   Creating the crucial mapping arrays (`DOFnis`, `DOFtoAtom`, `DOFcofefs`) that tell the `assembleAndRegularize` kernel how to gather scattered derivative contributions for each DOF.
*   **Buffer Management:** Allocating OpenCL memory buffers on the GPU and managing data transfers.
*   **Kernel Execution:** Launching the OpenCL kernels with the correct global and local work sizes.
*   **Optimization Loop / API:**
    - `getError(dofs)` runs the energy-only kernel to return `J_phys` built from `Emols` and `ErefW`. No regularization.
    - `getErrorDerivs(dofs)` runs `evalSampleDerivatives_template` and `assembleAndRegularize`, reads `Jmols` and `fDOFs`, and returns `(J_total, forces)` with `J_total = Σ Jmols - J_reg`.
    - `get_forces(dofs)` is available for force-only use; for derivative testing and optimization, prefer `getErrorDerivs()` to keep the single-kernel semantics.

#### 3.3.1 Template compilation and model injection

- The derivative and energy kernels support model injection via markers:
  - `//<<<MODEL_PAIR_ACCUMULATION` (derivatives + energy accumulation with REQH partials)
  - `//<<<MODEL_PAIR_ENERGY` (energy-only accumulation)
- The driver compiles with `FittingDriver.compile_with_model(macros={...})`, setting `use_template=True` so `set_kernel_args()` binds the templated kernels and their extra arguments (e.g., `Jmols`).

### 3.4 Debugging notes: zero derivatives and dual-pass workaround

This section records a recurring pitfall that led to near-zero assembled derivatives and the practical fix adopted.

- __[symptom]__ Final DOF forces were ~0 even when energies and per-sample objectives were non-zero.
- __[root cause]__ The derivative kernel `evalSampleDerivatives_template` writes `dEdREQs` only for fragment-1 lanes (`iL < ni`), i.e. indices `[i0, i0+ni)`. If the DOF assembly mapping (`DOFtoAtom`) includes fragment-2 atoms (`[j0, j0+nj)`), those entries remain zero and cancel the reduction in `assembleAndRegularize`.
- __[first fix]__ On the host, restrict `DOFtoAtom` to fragment-1 atoms only so the assembler only gathers indices that are written by the derivative kernel. Implemented in `FittingDriver.prepare_host_data()` by masking `host_ranges` fragment-1 spans. This immediately restores non-zero assembled forces for fragment-1–only DOFs.
- __[workaround to cover both fragments]__ To also obtain derivatives for fragment-2 without complicating the kernel, run the derivative kernel twice with swapped fragment ranges:
  - Build a swapped ranges array per-sample: `(i0, j0, ni, nj) -> (j0, i0, nj, ni)`.
  - Launch pass-1 with original `ranges` to populate `[i0, i0+ni)`.
  - Launch pass-2 with `ranges_swapped` to populate `[j0, j0+nj)`.
  - `Jmols` is a per-sample scalar and is simply re-written; reading after pass-2 yields the correct objective. `dEdREQs` becomes fully populated across both fragments after the two passes.

#### Driver support (host implementation)

- __[helpers]__ `FittingDriver` adds:
  - `self._ensure_swapped_ranges()`: constructs `host_ranges_swapped` and uploads `ranges_swapped_buff` (size in bytes) once.
  - `self._set_eval_ranges(rng_buf)`: rebinds only the `ranges` argument of the active eval kernel (templated or non-templated signatures handled).
- __[dual-pass API]__
  - `getErrorDerivs(..., dual_pass=False)` and `get_forces(..., dual_pass=False)` accept `dual_pass=True` to execute the second derivative pass with swapped ranges. Assembly remains unchanged and sums both fragments.
- __[assembly invariants]__ `assembleAndRegularize` is unchanged. It always reduces `dot(DOFcofefs[j], dEdREQs[ia])` over `DOFtoAtom`. Consistency rule still applies: every `ia` in `DOFtoAtom` must be written by the derivative producer kernels for the launches performed. With `dual_pass=True`, both fragments' `dEdREQs` are produced within the same host call, enabling broader `DOFtoAtom` coverage if desired.

### 3.5 Recommended debugging strategy (single H2O)

Use a minimal dataset to validate kernels, argument binding, DOF mapping, and the EvdW chain rule.

- __[data]__ `tests/tFitREQ/H2O_single.xyz` (1 sample, 3 atoms per fragment)
- __[DOFs]__ `tests/tFitREQ/dofSelection_H2O_MorseSR.dat` (R,Q,H). To exercise EvdW chain rule (E stored as sqrt(E)), optionally add comp=1 lines:

```text
E_HO 1  1e-6  1.0   0.00  0.00  0 0 0  0.0007  1.0
E_O3 1  1e-6  1.0   0.00  0.00  0 0 0  0.0026  1.0
```

- __[compile templated kernels]__ Use `FittingDriver.compile_with_model(...)` with `MODEL_MorseQ_PAIR` and `HBOND_GATE` define to enable the single-kernel derivative path (`use_template=True`).

- __[single-pass vs dual-pass]__ Keep `DOFtoAtom` restricted to fragment-1 (already done in `FittingDriver.prepare_host_data()`), which matches `evalSampleDerivatives_template` writing only frag-1. Use `dual_pass=True` only when you deliberately include fragment-2 atoms in `DOFtoAtom`.

- __[run scan and compare]__ Prefer the small H2O case for readable output and fast iteration:

```bash
python3 tests/tFitREQ/test_derivatives_ocl.py \
  --xyz H2O_single.xyz \
  --dof dofSelection_H2O_MorseSR.dat \
  --model MODEL_MorseQ_PAIR \
  --points 41 --eps 1e-6 --tol 1e-6 \
  --regularize 0 --present_only
```

If you added the optional comp=1 (E) DOFs, run the same command pointing to your modified DOF file.

- __[inspect per-atom derivatives]__ Call `FittingDriver.dump_dEdREQs(sample=0)` after a derivative evaluation to print a few `dEdREQs` rows. Expect non-zero entries for the frag-1 atoms only. For the E-component, the macro fix applies the chain rule via `REQj.y/REQi.y` scaling inside the model accumulation.

- __[expected checks]__
  - Analytic vs numeric dJ/dx agree: `max|F_ana - F_num| <= tol` (e.g., 1e-6).
  - `Jmols` is finite; `dEdREQs` not all zeros; E-component varies when scanning comp=1 DOFs.
  - If zeros occur: verify `DOFtoAtom` excludes frag-2 for single-pass, templated kernels are bound (`use_template=True`), `globParams` passed only for templated path, and `HBOND_GATE` is set as intended.

#### Recommendations

- __[simple path]__ Keep `DOFtoAtom` restricted to fragment-1 for stability during development (single-pass). Enable `dual_pass=True` only when you deliberately include fragment-2 atoms in `DOFtoAtom`.
- __[future work]__ Consider extending the derivative kernel to write both fragments in one pass (e.g., second write for `iL < nj` using `j0` base) to avoid the extra launch. Ensure local-memory usage and barriers remain correct.

### 3.6 MorseQ-only: Zero-derivatives issue — findings and focused debug plan

* __[scope]__ Strictly use `MODEL_MorseQ_PAIR` and focus on why analytic derivatives are zero. All other forcefields are out-of-scope until this is solved.

* __[what MorseQ macro does]__ In `cpp/common_resources/cl/Forces.cl` the block `//>>>macro MODEL_MorseQ_PAIR` (lines ~93–116) accumulates per-atom REQH partials and energy. For EvdW (component=1), it applies the chain rule for `E0 = sqrt(Ei*Ej)` directly in the pair accumulation:
  - `fREQi.y += dE_dE0 * 0.5f * (REQj.y / REQi.y);`  // derivative w.r.t. `REQi.y = sqrt(Ei)`
  - `tREQHs` stores EvdW as `sqrt(E)`; host updates use `update_tREQHs_from_dofs()`.

* __[where zeros can come from]__
  - __A. LdE = 0__: `evalSampleDerivatives_template` (in `cpp/common_resources/cl/FitREQ.cl`) writes `dEdREQs = fREQi * LdE`, where `LdE = W * (Emol - Eref)`. If `Emol ≈ Eref`, all analytic derivatives become zero even if `fREQi` is non-zero.
  - __B. Mapping mismatch__: Assembler gathers only indices present in `DOFtoAtom`. If it includes atoms the derivative kernel didn’t write (e.g., fragment-2 in single-pass), assembled DOF gradients are ~0. We already restrict to fragment-1 in `pyBall/OCL/NonBondFitting.py: FittingDriver.prepare_host_data()`.
  - __C. Missing EvdW DOFs__: DOF file lacks `comp=1` for types present in the dataset’s fragment-1 (e.g., `O_3`, `H_O` in `tests/tFitREQ/H2O_single.xyz`). Then no EvdW scan happens, or forces assemble to zero.

* __[minimal reproducible setup]__
  - __Data__: `tests/tFitREQ/H2O_single.xyz`.
  - __DOFs__: `tests/tFitREQ/dofSelection_H2O_MorseSR_plusE.dat` (duplicates MorseSR and adds `comp=1` for `O_3` and `H_O`).
  - __Run (MorseQ only, EvdW-only, no regularization)__:
    ```bash
    python3 tests/tFitREQ/test_derivatives_ocl.py \
      --xyz tests/tFitREQ/H2O_single.xyz \
      --dof tests/tFitREQ/dofSelection_H2O_MorseSR_plusE.dat \
      --model MODEL_MorseQ_PAIR \
      --only_comp E --present_only \
      --points 41 --eps 1e-6 --tol 1e-6 \
      --regularize 0 --hb_gate 1
    ```

* __[diagnostics to separate causes]__
  - __Check LdE__: Kernel prints (when `iDBG==0`) show `Emol`, `Eref`, `LdE`. If `LdE≈0`, analytic derivatives will be zero by construction; numeric `dJ/dx` should also be near zero over small scans.
  - __Inspect per-atom derivatives__: After `fit.getErrorDerivs(x)`, dump scaled per-atom derivatives for sample 0:
    - `FittingDriver.dump_dEdREQs(sample=0)` prints `dEdREQs[ia] = fREQi * LdE`. If `LdE≠0`, non-zero `dEdREQs[:,1]` confirms EvdW path is active. If `LdE=0`, divide by printed `LdE` (mentally) to infer `fREQi` sign/magnitude, or add a temporary print of raw `fREQi` before scaling.
  - __Verify mapping__: Ensure `DOFtoAtom` includes fragment-1 atoms of the scanned types (`O_3`, `H_O`). This is enforced by the frag-1 mask in `prepare_host_data()`; if you switch to `dual_pass`, you may include fragment-2 atoms as well.
  - __HBOND gating__: `HBOND_GATE` impacts the H-term only; it should not zero EvdW derivatives. For sanity, you can re-run with `--hb_gate 0` to exclude gating effects.

* __[success criteria]__
  - `Jmols` finite; during scans `Es(x)` varies (not identically zero).
  - Over inner scan points: `max|F_ana - F_num| <= tol` (e.g., 1e-6) for EvdW DOFs.
  - `dEdREQs[:,1]` non-zero for fragment-1 atoms when `LdE≠0`.

* __[if still zero]__
  - Confirm the DOF file actually includes `comp=1` for types present in the sample’s fragment-1 (use the provided `..._plusE.dat`).
  - Temporarily adjust `Eref` or pick scan points where `Emol−Eref` is not vanishing, to ensure `LdE` is not suppressing derivatives.
  - Add a debug print in `evalSampleDerivatives_template` to print raw `fREQi` (before scaling) for lanes `iL<ni` to confirm pair accumulation is producing EvdW partials.

#### Additional troubleshooting: EvdW derivative scaling and validation

- __EvdW stored as sqrt(E): chain rule fix__
  - Parameters store `EvdW` as `sqrt(E)` in both kernels and Python DOFs.
  - Correct accumulation for component `E` is:
    ```c
    // in model macros (e.g., MODEL_MorseQ_PAIR) and non-templated kernel
    fREQi.y += dE_dE0 * REQj.y;  // REQj.y = sqrt(Ej), accumulates ∂E/∂sqrt(Ei)
    ```
  - Fixed in `cpp/common_resources/cl/Forces.cl` (templated macros) and `cpp/common_resources/cl/FitREQ.cl` (non-templated), removing an erroneous `0.5f` factor and division by `REQi.y`.
  - Objective consistently uses `0.5 * W * (Emol - Eref)^2` and gradients multiply by `LdE = W*(Emol - Eref)`.

- __Derivative test: robust diagnostics__
  - `tests/tFitREQ/test_derivatives_ocl.py` now prints min/max and abs-max for analytic and numeric derivatives and accepts `--signal_eps`.
  - A DOF comparison is judged only when both analytic and numeric abs-max exceed `signal_eps`; otherwise it is marked `SKIP(weak)` to avoid false "matches" from near-zero derivatives.

- __Exact test command (MorseQ, E-only)__
  - Path semantics: the script joins non-absolute `--xyz/--dof` with its own directory `tests/tFitREQ/`, so pass filenames relative to that directory (or use absolute paths).
  - Example:
    ```bash
    python3 tests/tFitREQ/test_derivatives_ocl.py \
      --xyz H2O_single.xyz \
      --dof dofSelection_H2O_MorseSR_plusE.dat \
      --model MODEL_MorseQ_PAIR \
      --only_comp E \
      --points 5 --eps 1e-6 --tol 1e-6 \
      --regularize 0 --verbose 2 --signal_eps 1e-12 | tee tests/tFitREQ/OUT-nonBond
    ```

- __Common pitfalls checklist__
  - `LdE ≈ 0`: If `Emol ≈ Eref`, analytic derivatives scale to zero. Expect numeric dJ/dx to be small as well over tiny scans.
  - DOF mapping includes fragment-2 atoms while derivative kernel writes only fragment-1 (single-pass). Use frag-1-only mapping (default) or `dual_pass=True`.
  - `HBOND_GATE`: Only gates H-term; does not zero EvdW. For sanity, try `--hb_gate 0`.

### 3.5.1 Identical-type pairs (ti==tj): factor-of-two in dE/dsqrt(E)

When parameters store EvdW as `sqrt(E)` and the mixing rule is `E0 = Ei * Ej` with `Ei = sqrt(Eii)`, then for identical types (`ti == tj`) we have `Ei == Ej` and `E0 = Ei^2`. The derivative with respect to the stored variable `Ei = sqrt(E)` therefore doubles compared to the generic case:

- Generic: `∂E/∂Ei = (∂E/∂E0) * (∂E0/∂Ei) = dE_dE0 * Ej`
- Identical types (`ti==tj`): `Ej == Ei` and `E0 = Ei^2` ⇒ `∂E0/∂Ei = 2*Ei` ⇒ `∂E/∂Ei` is larger by a factor of 2

Implementation in OpenCL kernels:

- Non-templated `evalSampleDerivatives` in `cpp/common_resources/cl/FitREQ.cl`:
  - Caches `tj` per local tile in `LTYPES[32]`.
  - Accumulates E-component as `fREQi.y += (ti==tj ? 2.f : 1.f) * dE_dE0 * REQj.y;`
  - Heavy printf debugging is guarded by `if(iG==iDBG)` to avoid excessive output.

- Templated `evalSampleDerivatives_template` in `cpp/common_resources/cl/FitREQ.cl`:
  - Leaves injected model code unchanged for portability.
  - Before the injected block, stores `float y0_pair = fREQi.y;` and after it, post-scales only this pair's incremental E-component: `dy = fREQi.y - y0_pair; if(ti==tj) fREQi.y = y0_pair + 2.f*dy;`
  - Uses a local `LTYPES[32]` buffer to access `tj` without extra global reads.

Notes:

- CPU reference (`cpp/common/molecular/FitREQ.h`) is intentionally left unchanged as a stable baseline; the ti==tj doubling is an OpenCL-side correction only.
- The ti==tj doubling is independent of whether you run single-pass or dual-pass over fragments; dual-pass remains the recommended way to populate `dEdREQs` for both fragments when needed.

This plan keeps the investigation strictly within MorseQ, isolates the LdE scaling from the pairwise derivatives, and verifies host-side DOF mapping so we can eliminate the zero-derivatives failure mode.

NOTE: but is this really a problem? Because anyway we accumulate this derivative 2x (?) if we accumulate derivatives from both atoms (both fragments) ????

# Part 4: Electron Pair & Sigma Hole Addition (`add_epairs.py`)

## 4.1 Purpose

The script `tests/tFitREQ/add_epairs.py` adds dummy atoms representing **electron pairs**
(lone pairs on O, N, F) and **sigma holes** (positive regions on H, F, Cl, Br) to
multi-frame XYZ dimer files. This is a preprocessing step: the resulting XYZ files with
explicit `E` atoms can be loaded by the fitting pipeline (`loadXYZ(fname, bAddEpairs=False)`)
without needing the C++ builder to add them on-the-fly.

The script is separate from the main fitting pipeline and works purely in Python using
`AtomicSystem` from `pyBall`.

## 4.2 What Gets Added

The script reads `AtomTypes.dat` to decide per atom type:

| Type | Electron pairs (nepair) | Sigma-hole epair_name |
|------|------------------------|----------------------|
| O_3 (water oxygen) | 2 | — |
| O_2 (carbonyl oxygen) | 2 | — |
| N_3 (ammonia N) | 1 | — |
| N_2 (pyridine N) | 1 | — |
| N_1 (HCN nitrile N) | 1 | — |
| F_ (fluorine) | 1 | — |
| H_O (hydroxyl H) | 0 | E_H_O |
| H_N (amine H) | 0 | E_H_N |
| H_F (HF H) | 0 | E_H_F |
| H_C1 (alkyne H) | 0 | E_H_C1 |
| C, H_a, etc. | 0 | none |

**Geometry rules** (implemented in `AtomicSystem.make_epair_geom`):

| Hybridization | Example | npi | nσ | Geometry |
|---------------|---------|-----|-----|----------|
| sp³ (nb=3) | NH₃ | 0 | 3 | 1 epair opposite the 3-bond average |
| sp³ (nb=2) | H₂O | 0 | 2 | 2 epairs in tetrahedral positions |
| sp² (nb=2) | =N− | 1 | 2 | 1 epair along bisector |
| sp² (nb=1) | =O | 1 | 1 | 2 epairs in pi-plane |
| sp (nb=1) | ≡N | 2 | 1 | 1 epair opposite the bond |

Sigma holes are placed along the bond direction, outward from the atom (opposite its
single neighbor).

## 4.3 Usage

```
python tests/tFitREQ/add_epairs.py -i <input> -o <output> [options]
```

### Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `-i, --input` | required | Input .xyz file or directory of .xyz files |
| `-o, --output` | `output` | Output file or directory |
| `--mode` | `epairs` | What to add: `epairs`, `sigma`, or `both` |
| `--lepair` | `1.0` | Distance (Å) of epair dummy atoms from host |
| `--sigma-dist` | `0.5` | Distance (Å) of sigma-hole dummy atoms from host |
| `--atypes` | auto | Path to `AtomTypes.dat` |
| `--simple-names` | off | Output simple element names (C, N, O, H, E) instead of full type names |
| `-v, --verbose` | off | Show per-frame details |

### Examples

```bash
# Single file, electron pairs only:
python tests/tFitREQ/add_epairs.py -i scan.xyz -o scan_ep.xyz

# Single file, both epairs and sigma-holes, custom distances:
python tests/tFitREQ/add_epairs.py -i scan.xyz -o scan_full.xyz \
    --mode both --lepair 1.2 --sigma-dist 0.6

# Process all .xyz files in a directory, output to another directory:
python tests/tFitREQ/add_epairs.py \
    -i /path/to/confs/wb97m-split/ \
    -o /path/to/output/ \
    --mode both --lepair 1.0

# With simple names for Jmol visualization:
python tests/tFitREQ/add_epairs.py -i scan.xyz -o scan_ep.xyz \
    --mode both --simple-names

# Verbose output (shows each file and frame detail):
python tests/tFitREQ/add_epairs.py -i input_dir/ -o out_dir/ -v

# Re-run a sub-set that was already processed (skips _bak and _Epairs files):
python tests/tFitREQ/add_epairs.py -i input_dir/ -o output_dir/ --mode both
```

## 4.4 Output Format

Each frame in the output follows the same multi-frame XYZ format as the input:

```
<natoms_after_addition>
# n0 <n0_updated> Etot <energy> x0 <dist> z <angle> <molname>
<atom_type_1> <x1> <y1> <z1> <charge1>
...
<E_dummy> <x> <y> <z> <charge=0>
...
```

Key differences from input:

- **n0 is updated** to `n0_original + n_added_to_molA`
- Added dummy atoms are listed **after** all their host molecule's real atoms, in the order:
  `[molA_real] [molA_epairs] [molA_sigma] [molB_real] [molB_epairs] [molB_sigma]`
- Atom type names for real atoms are preserved (e.g., `N_3`, `C_2`, `H_O`)
- Added dummies are labeled `E` (atomic number 200)
- Charges on real atoms are preserved; dummies get `0.0`

## 4.5 Implementation Details

### 4.5.1 File: `tests/tFitREQ/add_epairs.py`

**Build-fragment loop** (for each of molecule A and B in each frame):

1. Create an `AtomicSystem` from the fragment's atoms
2. Compute bonds via distance-based detection (`findBondsNP`)
3. For each atom with `nepair > 0` in `AtomTypes.dat`: call `make_epair_geom(i, npi, nb, distance)`
4. For each atom with `nepair == 0` but an explicit `epair_name`: place a sigma-hole dummy
   opposite its single bond
5. Restore original type names for real atoms
6. Concatenate both fragments and update `n0`

### 4.5.2 File: `pyBall/AtomicSystem.py` (methods used)

| Method | Lines | Role |
|--------|-------|------|
| `add_electron_pairs(distance)` | 951-964 | Iterates atoms, dispatches to `make_epair_geom` |
| `make_epair_geom(i, npi, nb, distance)` | 980-1025 | Places 1-2 dummies per atom depending on hybridization |
| `place_electron_pair(i, dir, distance)` | 1027-1052 | Appends dummy atom to arrays |
| `neighs()` | 171-175 | Builds neighbor list from bonds |

### 4.5.3 Distance conventions

The `--lepair` default is **1.0 Å** (matching the C++ `Lepairs` default in `FitREQ.h`
line 490). The `--sigma-dist` default is **0.5 Å** (matching the C++ sigma-hole placement
convention).

### 4.5.4 Dependency on `AtomTypes.dat`

The script reads columns 5 (nepair) and 4 (epair_name) to decide what each atom type gets:

```
#name parent element epair  nv ne npi sym  Ruff   RvdW    EvdW   Qbase  Hb ...
O_3    O      O      E_O_3  2  2   0   0   0.658  1.7500  0.0026 -0.1  -0.9  ...
H_O    H_     H      E_H_O  1  0   0   0   0.354  1.4430  0.0019  0.1   0.9  ...
```

- `column 5` (nv/valence): total bond order
- `column 6` (nepair): number of lone pairs
- `column 4` (epair_name): for atoms with nepair=0, if this is not `*`, the atom is a
  sigma-hole donor and gets 1 dummy opposite its bond

The auto-search path looks for `AtomTypes.dat` in:
1. `tests/tFitREQ_PN/data/`
2. `tests/tFitREQ/data/`
3. `tests/tFitREQ/` (same directory as the script)

## 4.6 Practical Notes

- The script is **fast**: ~260k frames in under 2 minutes on a modern CPU.
- Output files can be large: each frame grows by ~5–15 atoms depending on molecule size.
- The `--simple-names` option is useful for visualizing in Jmol, Avogadro, or VMD where
  type names like `H_O` are not recognized.
- If you get an error about `ELEMENT_DICT` lookup failing for `E` atoms, you are trying
  to process a file that already contains dummy atoms. Use a fresh input file, or the
  script will skip files with `_bak` or `_Epairs` in the name.
- For molecules with no N or O atoms (e.g., HF dimers), `--mode epairs` adds nothing;
  use `--mode sigma` or `--mode both` to add sigma holes instead.


# Part 5: Reference Energy Map Plotting (`plot_ref.py`)

## 5.1 Purpose

The script `tests/tFitREQ/plot_ref.py` produces publication-style multi-panel figures
for 2D angle–distance scan data. Given one or more multi-frame XYZ files it can
generate a grid of:

- **Row 1** — 2D imshow of the reference energy as a function of angle (x) and
  distance (y), with a white `+` at the global minimum and an optional black
  R_min(angle) curve overlaid.
- **Row 2** — Line plot of E_min(angle) (black solid) and the angular energy
  slice at the distance of the global minimum (red dashed).
- **Row 3** — Ball-and-stick geometry of the dimer at the global-minimum frame,
  coloured by element (dummy `E` atoms in purple, others by
  `elements.ELEMENT_DICT[e][8]`), labelled with 1-based indices (primed `'` for
  the acceptor fragment when n0 is known).

## 5.2 Files and Dependencies

| File | Role |
|------|------|
| `tests/tFitREQ/plot_ref.py` | The main script (CLI entry point) |
| `pyBall/FitREQutils.py` | Pure-Python utilities extracted from `FitREQ.py` + `split_scan_imshow_new.py` |
| `pyBall/elements.py` | Element colour / radius look-up tables |
| `pyBall/atomicUtils.py` | `findBondsNP()`, `load_xyz_movie()`, etc. |

The script requires **no C library** and **no OpenCL** — it works entirely in
Python + numpy + matplotlib.

## 5.3 Usage

```
python tests/tFitREQ/plot_ref.py -i <input>... [options]
```

### Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `-i, --input` | required | One or more XYZ file paths (or names with `--dir`) |
| `--dir` | — | Directory to resolve partial names in `--input` |
| `--glob` | — | Glob pattern to select files inside `--dir` (e.g. `"*HCN*"`) |
| `-o, --output` | — | Save figure to this PNG path (otherwise `plt.show()`) |
| `--kcal` | off | Convert energies to kcal/mol |
| `--sym` | off | Symmetric colour scale `vmin = min(E), vmax = -vmin` |
| `--min-lines` | off | Overlay R_min curve + show energy panel + molecule geometry |
| `--ncols` | auto | Number of columns in the subplot grid |

### Examples

```bash
# Simple 2D maps (one row per system):
python tests/tFitREQ/plot_ref.py \\
    -i HCN-D1_CH2NH-A1-y H2O-D1_NH3-A1-y \\
    --dir /path/to/scans/ --kcal

# Full three-row layout:
python tests/tFitREQ/plot_ref.py -i scan.xyz --min-lines --kcal

# 8 systems, 4 columns:
python tests/tFitREQ/plot_ref.py -i system1.xyz ... system8.xyz \\
    --min-lines --kcal --ncols 4 -o out.png

# Using a glob:
python tests/tFitREQ/plot_ref.py --dir /path/to/scans/ --glob "*H2O*"
```

## 5.4 How It Works

### 5.4.1 Pipeline

For each XYZ file, `plot_ref.py` calls `_build_frame_grid()` which replicates
the processing inside `compute_panel_data()` but additionally builds a
**frame-index grid**:

1. **Parse** the XYZ via `read_scan_atomicutils()` (or `parse_xyz_blocks()` as
   fallback) → `(Es, Ps)` or `(Es, Ts, Ps)`.
2. **Geometry** → `compute_ra_vec(Ps)` → `r` (distance per frame), `a` (angle
   per frame).
3. **Rows** → `detect_rows_by_r(r)` → list of `(start, end)` frame ranges, one
   per unique angle, sorted by appearance in the file.
4. **Frame-index grid** — a `(ny, nx)` integer array where cell `[iy, ix]`
   holds the original frame index `frame_idx[iy, ix]`, and `-1` for NaN-padded
   cells.
5. **Reshape** → `reshape_to_grid(Es, r, a, rows)` → `(Vraw, R, A, rv)`.
6. **Shift** → `compute_shift_from_grid(Vraw)` → `V = Vraw - shift`.

The resulting `V` (shifted energy grid) and `frame_idx` are then used for all
three rows.

### 5.4.2 Row 1 — 2D imshow

`plot_imshow()` from `FitREQutils.py` is called with `bSym=False`, then the
returned `im` object has its colour limits overridden:

```python
vmin_ref = np.nanmin(V) * conv
if vmin_ref < 0:
    im.set_clim(vmin_ref, -vmin_ref)       # symmetric around the well
else:
    im.set_clim(-vabs, vabs)               # all-positive → symmetric max|V|
```

The **global minimum** is marked with a white `+` symbol at the pixel centre
(`y0 + (ix + 0.5) * (y1 - y0) / nx`), matching the imshow `extent`.

The **R_min(angle) curve** overlaid on the image is drawn by converting the
pixel column index of each row's minimum into a data y-coordinate via the same
extent formula, then plotting `x = A[sorted], y = pixel_center_y`.

### 5.4.3 Row 2 — Energy curves

`extract_min_curves(A, rv, V.T)` returns `(rmin, emin)` — the per-angle minimum
distance and energy. The global minimum indices `(iy_g, ix_g)` give the
fixed-r slice: `V[:, ix_g]` plotted against `A`, with a legend entry "E @
r=xx.xx".

### 5.4.4 Row 3 — Molecule geometry

The frame index at the global minimum `frame_idx[iy_g, ix_g]` is used to
extract atom type names and positions from `Ts_raw[frame]` and `Ps_raw[frame]`.
The function `plot_molecule()` draws:

- **Scatter markers** coloured by element symbol (`E` → `'#9932CC'` purple,
  others → `elements.ELEMENT_DICT[e][8]` hex colour).
- **Atom labels** as 1-based indices; if `n0` is known, the acceptor fragment
  indices are primed (`1', 2', …`).
- **n0 separator** — a vertical dotted line at the midpoint between the two
  fragments, with "Donor" / "Acceptor" annotations.

## 5.5 Colour Scaling

The default auto-scale (`--sym` or implied by `--min-lines`) is:

    vmin = min(E_map)              # most negative (deepest well)
    vmax = -vmin                   # symmetric

If every data point is repulsive (`min > 0`), `vmax = max|V|` is used instead,
avoiding an inverted colour range.

This differs from matplotlib's auto-scale which would centre on the data mean.
The symmetric scaling ensures attractive wells are clearly visible even when
the repulsive wall is much larger.

## 5.6 Known Pitfalls and Future Work

### 5.6.1 Frame-index mismatch for certain files

Files such as `H2O-A1_HCN-D1-y.xyz` and `H2O-A1_HCN-D1-z.xyz` produce a
mismatch between the shift-subtracted grid `V` and the frame-index grid built
from the raw energy array `Es`. The likely cause is that `compute_panel_data()`
internally replaces `Es` with header-derived energies `Eh` when they have the
same size and finite values (`if Eh.size == Es.size and
np.any(np.isfinite(Eh)): Es = Eh`).  The header energies come from
`parse_headers_ra()`, while the frame indexing uses the raw geometry parsing.
If the ordering or count of header-derived values differs from the raw parse,
the frame index no longer points to the correct configuration.

**Workaround**: these files still render the 2D map and energy curves
correctly; only the molecule geometry in row 3 may show a wrong or empty
frame.  A proper fix would require `compute_panel_data()` to also return (or
accept) a frame-index grid so that the energy replacement does not desync the
indices.

### 5.6.2 All-zero angles (single-angle scans)

Some scan files contain data at only one angle (e.g. `CH2O-A1_HCN-D1-y.xyz`
where all `A` values are `0`).  In that case:

- The 2D imshow still displays correctly (the extent pads to avoid zero
  width).
- The R_min overlay is skipped (all angles identical).
- The molecule geometry panel still shows the correct frame.

### 5.6.3 E (dummy) atom bond detection

`findBondsNP()` from `pyBall.atomicUtils` calls `getAtomRadiusNP()` which
indexes into `elements.ELEMENTS` by atomic number.  Dummy `E` atoms have
Z=200, which is beyond the 113-element table.  The script avoids this by not
drawing bonds — only scatter markers are used for the geometry panel.

If bond lines are ever needed again, they must be computed while skipping `E`
atoms or by treating each fragment separately with the host–epair connectivity
defined by `n0`.

### 5.6.4 Speed

Parsing a 703-frame XYZ file and building the three-row figure takes ~1–2
seconds per system. For 8 systems at 4 columns this yields a figure in under
15 seconds.

### 5.6.5 `pyBall/FitREQutils.py` maintenance

All pure-Python functions from `pyBall/FitREQ.py` and
`tests/tFitREQ/split_scan_imshow_new.py` have been consolidated into
`pyBall/FitREQutils.py`.  The file `split_scan_imshow_new.py` now imports from
`FitREQutils` and only keeps its own `if __name__ == '__main__':` CLI block.
When adding new features to any shared parsing/plotting function, edit
`FitREQutils.py` directly.


# Part 6: Frame Reordering (`reorder_frames.py`)

## 6.1 Purpose

The script `tests/tFitREQ/reorder_frames.py` reorders frames in multi-frame XYZ
files based on geometric parameters extracted from comment headers. This is
useful when scan data is generated with frames in a non-sequential order (e.g.,
angle=0 appears at the start instead of in the middle of the -90...+90 sequence).

The script sorts frames by:
1. **Angle** (y or z field) - ascending order
2. **Distance** (x0 field) - ascending order

This ensures that subsequent analysis and plotting tools receive data in the
expected sequential order.

## 6.2 Usage

```
python tests/tFitREQ/reorder_frames.py -i <input> [options]
```

### Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `-i, --input` | required | Input .xyz file or directory |
| `-o, --output` | overwrite input | Output .xyz file path |
| `--dir` | off | Process all .xyz files in input directory |

### Examples

```bash
# Reorder single file (overwrites input):
python tests/tFitREQ/reorder_frames.py -i scan.xyz

# Reorder single file to new output:
python tests/tFitREQ/reorder_frames.py -i scan.xyz -o scan_reordered.xyz

# Reorder all .xyz files in a directory (in-place):
python tests/tFitREQ/reorder_frames.py -i /path/to/scans/ --dir
```

## 6.3 Implementation

The script uses `reorder_frames_by_angle()` from `pyBall/FitREQutils.py` which:

1. Parses the XYZ file to extract:
   - Number of atoms (natoms)
   - Comment line containing: `n0`, `Etot`, `x0` (distance), and `y` or `z` (angle)
   - Atom lines for each frame
2. Extracts geometric parameters from comments using regex:
   - `energy_re`: matches `Etot` or `E=` followed by a float
   - `reR`: matches `x0` followed by a float
   - `reAy`/`reAz`: matches `y` or `z` followed by a float
3. Stores each frame as a dict with: `natoms`, `comment`, `atoms`, `angle`, `distance`
4. Sorts frames by `(angle, distance)` tuple
5. Writes reordered frames back to the output file

The parsing correctly handles XYZ files where the `natoms` line appears **before**
the comment line (standard format).

## 6.4 Integration with `add_epairs.py`

Typical workflow:

```bash
# 1. Generate electron-paired files:
python tests/tFitREQ/add_epairs.py -i scans/ -o scans_ep/ --mode both

# 2. Reorder frames for proper sequential ordering:
python tests/tFitREQ/reorder_frames.py -i scans_ep/ --dir

# 3. Plot reordered data:
python tests/tFitREQ/plot_ref.py --dir scans_ep/ --min-lines --kcal
```

The reordering step is necessary because some scan generators produce frames with
angle=0 at the beginning of the file, which breaks the expected -90...0...+90
sequence used by plotting tools.

## 6.5 Notes

- The script preserves all frame data exactly, only reordering the sequence
- Comment lines are preserved with their original content (including molecule names)
- The sorting is stable: frames with identical (angle, distance) maintain their
  original relative order
- For directories, the script processes files in alphabetical order
- The script reports the number of frames reordered for each file
