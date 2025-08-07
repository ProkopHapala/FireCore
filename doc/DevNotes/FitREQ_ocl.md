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
            *   `charge`: The initial partial charge of the atom. This value is used by the `Q` component of the REQH model.

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

### 3.2 OpenCL Kernel Breakdown

The program relies on two main, highly optimized OpenCL kernels:

1.  **`evalAndScalePerSample`**: This is the primary workhorse kernel. It performs a three-stage process in a single pass for each sample:
    *   **Evaluate:** Each thread in the work-group calculates the per-atom energy contribution and the raw (unscaled) parameter derivatives (`∂E_model/∂p`) for its assigned atom.
    *   **Reduce:** The threads perform a **parallel reduction** using `__local` memory to efficiently sum all per-atom energy contributions into a single, total `E_model` for the sample. This is vastly superior to a serial loop or global atomics.
    *   **Scale:** The first thread in the work-group calculates the final error scaling factor `dEw = -2 * w * (E_model - E_ref)`. This single float value is broadcast to all other threads via `__local` memory. Each thread then multiplies its private, raw derivatives by this common `dEw` factor and writes the final, scaled derivative to global memory.

2.  **`assembleAndRegularize`**: This kernel takes the final scaled derivatives and calculates the total force on each DOF.
    *   **Assemble:** Each work-group is assigned a single DOF. It performs another parallel reduction, looping through all atoms that are affected by its assigned DOF and summing their contributions to get the total physical force on that DOF.
    *   **Regularize:** After the physical force is calculated, the first thread of the work-group reads the current value of the DOF and applies the regularization forces based on the soft/hard limits and tethering potentials defined in `dofSelection.dat`. The final result (physical force + regularization force) is written to global memory.

### 3.3 The Python Driver (`FittingDriver`)

The Python class acts as the high-level orchestrator for the entire process. Its key responsibilities include:

*   **File I/O:** Loading and parsing the three input files (`AtomTypes.dat`, `input.xyz`, `dofSelection.dat`).
*   **Data Preparation:**
    *   Building the initial parameter matrix (`tREQHs_base`) by combining default values from `AtomTypes.dat` with `xstart` values from `dofSelection.dat`.
    *   Creating the crucial mapping arrays (`DOFnis`, `DOFtoAtom`, `DOFcofefs`) that tell the `assembleAndRegularize` kernel how to gather scattered derivative contributions for each DOF.
*   **Buffer Management:** Allocating OpenCL memory buffers on the GPU and managing data transfers.
*   **Kernel Execution:** Launching the OpenCL kernels with the correct global and local work sizes.
*   **Optimization Loop:** Providing a `get_forces(dofs)` method that encapsulates the entire GPU-side calculation. This method can be plugged into any external, Python-based optimizer. The default implementation uses the **FIRE (Fast Inertial Relaxation Engine)** algorithm, a robust and efficient choice for force-based energy minimization.