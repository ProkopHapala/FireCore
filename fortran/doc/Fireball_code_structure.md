## Fireball Code Structure Documentation

This document outlines the structure of the Fireball DFTB Fortran code, based on the provided file listing and analysis.

### Overall Philosophy of Fireball (as inferred)

Fireball is a Density Functional Tight Binding (DFTB) code, characterized by:

1.  **Numerical Atomic Orbitals (NAOs):** It utilizes pre-calculated, spherically symmetric numerical atomic orbitals as its basis set. These are typically generated from atomic DFT solutions and are not analytical functions like Gaussians.
2.  **Pre-calculated Integrals (`Fdata`):**
    *   Two-center integrals (overlap $S_{ij}$, kinetic energy $T_{ij}$, and parts of the potential) between NAOs on different atoms are pre-calculated and stored in tables (e.g., in a directory named `Fdata`) as a function of interatomic distance.
    *   Three-center integrals, crucial for more accurate XC functionals and charge-dependent effects (like in McWEDA), are also pre-calculated and tabulated, depending on the geometry of the three atoms.
    *   One-center terms, especially for exchange-correlation, are also pre-calculated.
3.  **Tight-Binding Framework:** The Hamiltonian matrix elements are constructed using these pre-calculated integrals, often involving interpolation and rotation of these tabulated values to match the actual molecular geometry.
4.  **Self-Consistent Field (SCF):** The electronic structure is solved iteratively:
    *   An initial guess for atomic charges (or electron density distribution) is made.
    *   A Hamiltonian is constructed based on these charges.
    *   The Schrödinger-like equation $Hc = eSc$ is solved to get new orbitals and energies.
    *   New charges/density are calculated from the new orbitals.
    *   The new and old charges are mixed, and the process repeats until self-consistency.
5.  **Approximations for XC and Self-Consistency:**
    *   **Harris Functional:** The **Harris Functional** (also known as Harris-Foulkes functional) is a non-self-consistent (non-SCF) approach that approximates the total energy of a system using a superposition of neutral atomic charge densities. It avoids the iterative self-consistency loop, making it computationally very efficient, especially for obtaining initial energy estimates or for systems where a full SCF treatment is not strictly necessary. While not fully self-consistent, it provides a robust and often surprisingly accurate first approximation.
    *   **Sankey-Niklewski (SNXC):** The **Sankey-Niklewski (SNXC)** method represents an early self-consistent field (SCF) scheme within the tight-binding framework. It was a foundational step towards incorporating charge self-consistency into FIREBALL, allowing for a more accurate description of charge transfer effects compared to non-SCF methods. However, it had limitations, particularly in its treatment of exchange-correlation (XC) interactions, which led to the development of more sophisticated approximations.
    *   **McWEDA (Multi-center Weighted Exchange-correlation Density Approximation):** The **Multi-center Weighted Exchange-correlation Density Approximation (McWEDA)** is a significant advancement in FIREBALL for accurately calculating exchange-correlation (XC) matrix elements. It addresses deficiencies found in earlier approximations (like SNXC) by employing a multi-center approach and weighted densities, which allows for a more precise description of the XC potential and energy. McWEDA is crucial for improving the accuracy of the self-consistent field (SCF) calculations within FIREBALL, particularly for systems with complex electronic structures.
    *   **Kohn-Sham Grid:** A **Kohn-Sham Grid** refers to a real-space numerical grid commonly employed in conventional Density Functional Theory (DFT) codes to perform numerical integrations, particularly for exchange-correlation terms. While FIREBALL's core methodology prioritizes efficiency by largely avoiding a real-space grid for its primary charge density and potential calculations (relying on local-orbital basis sets and tabulated integrals), the framework does allow for the use of a real-space grid for specific parts of the XC problem. This can enable the use of more standard and potentially more accurate DFT functionals, albeit at a higher computational cost compared to the highly optimized, grid-less approximations like McWEDA.

### Directory Overview

A brief overview of the main directories and their purpose:

*   `MAIN`         : This is the heart of the program, containing the main `fireball.f90` program, the `libFireCore.f90` library interface, the SCF loop drivers, and the core solver routines.
*   `MODULES`      : Defines global data structures, variables, and constants used throughout the code. These modules manage the state of the calculation.
*   `ALLOCATIONS`  : Contains routines for allocating and reallocating memory for major arrays (Hamiltonian, overlap, density, neighbor lists).
*   `ASSEMBLERS`   : Houses routines responsible for constructing (assembling) the Hamiltonian matrix, overlap matrix, and various energy components from pre-calculated integrals and current electronic/atomic configurations.
*   `DASSEMBLERS`  : Contains routines for calculating atomic forces, which are essentially derivatives of the energy terms assembled in the `ASSEMBLERS` directory.
*   `GRID`         : Includes functionalities related to real-space grid operations, such as projecting densities or orbitals onto a grid, FFTs, and potentially parts of a Kohn-Sham DFT solver on a grid.
*   `INITIALIZERS` : Provides routines for initializing various aspects of the calculation, such as basic constants, atomic information, neighbor lists, and initial charge states.
*   `INTERACTIONS` : Contains the core routines for calculating the actual one-center, two-center, and three-center interactions (integrals) based on atomic positions and orbital types.
*   `INTERPOLATERS`: Provides functions for interpolating values from the pre-calculated integral tables (e.g., using splines) to get values at specific interatomic distances/angles.
*   `MATH`         : Includes basic mathematical utility functions like factorial or cross products.
*   `NEIGHBORS`    : Contains algorithms for finding and managing neighbor lists for atoms, which is crucial for efficient calculation of short-ranged interactions.
*   `READFILES`    : Handles the parsing of input files (e.g., options, atomic coordinates) and reading of the `Fdata` integral tables.
*   `ROTATIONS`    : Provides routines for performing rotations of orbitals or integral components to align with the global coordinate system of the molecule.

### Analysis of Fortran Modules (`/home/prokop/git/FireCore/fortran/MODULES/`)

These modules define global variables and data structures.

*   **`options.f90`**:
    *   **Purpose:** Defines global variables that act as switches and parameters controlling the calculation's behavior, typically read from an input file.
    *   **Key Variables (examples):** `itheory`, `itheory_xc` (select theory level and XC functional), `iforce` (enable/disable forces), `icluster` (cluster vs. periodic), `iqout` (charge population analysis type), `ifixcharge` (fix atomic charges), `idebugWrite`, `verbosity` (debug output levels), `iwrtewf`, `igrid` (wavefunction/grid output).

*   **`interactions.f90`**:
    *   **Purpose:** Defines parameters and allocatable arrays for atomic orbitals and interaction matrix elements. Central to Hamiltonian and overlap matrix construction.
    *   **Key Variables (examples):** `norbitals`, `numorb_max` (orbital counts), `ME2c_max`, `ME3c_max` (max matrix elements for dimensioning), `index_max*`, `lssh`, `mu`, `nu`, `mvalue`, `nssh`, `degelec` (orbital indexing and mapping), `iatyp`, `imass` (atom types and masses), `h_mat`, `s_mat` (Hamiltonian and overlap matrices), `vna`, `vnl`, `vxc`, `vxc_1c` (potential components: nuclear-attraction, non-local pseudopotential, XC).

*   **`integrals.f90`**:
    *   **Purpose:** Declares arrays to hold pre-calculated integral values read from `Fdata` files.
    *   **Key Variables (examples):** `fdataLocation` (path to `Fdata`), `icon3c` (3-center integral indexing), `numx3c_bcna`, `hx_bcna`, `bcna_01` (parameters and arrays for tabulated 3-center "neutral atom" integrals), similar sets for `xc3c` (3-center XC) and `den3` (3-center density) integrals, `ind2c`, `numz2c`, `xintegral_2c`, `splineint_2c` (parameters and arrays for 2-center integrals and splines), `exc1c_0`, `xcnu1c` (1-center XC energies/potentials).

*   **`charges.f90`**:
    *   **Purpose:** Manages atomic charges and related quantities.
    *   **Key Variables (examples):** `Q`, `Qin`, `Qout` (atomic charges: current, input to SCF, output from SCF), `qtot`, `ztot` (total charge, total valence electrons).

*   **`density.f90`**:
    *   **Purpose:** Holds density matrix, wavefunctions, and related electronic structure data.
    *   **Key Variables (examples):** `rho` (density matrix), `cape` (energy-weighted density matrix), `bbnkre`, `bbnkim` (wavefunction coefficients, real and imaginary), `eigen_k` (eigenvalues), `foccupy` (occupation numbers), `efermi` (Fermi energy).

*   **`molecule.f90`**:
    *   **Purpose:** Likely a high-level container for the system's geometry and atomic properties.
    *   **Expected Key Variables:** `natoms` (number of atoms), `ratom` (atomic coordinates), `vatom` (velocities), `fatom` (forces), `avec` (lattice vectors).

*   **Other Notable Modules:**
    *   `constants_fireball.f90`: Defines physical and mathematical constants.
    *   `dimensions.f90`: Defines array dimensions and limits.
    *   `energy.f90`: Stores components of the total energy.
    *   `forces.f90`: Stores atomic forces.
    *   `grid.f90`: Variables related to real-space grid calculations.
    *   `kpoints.f90`: Manages k-points for periodic calculations.
    *   `neighbor_map.f90`: Data structures for neighbor lists.
    *   `timing.f90`: Variables for timing different parts of the code.
    *   `workmat.f90`: Allocatable work arrays used in various calculations (e.g., diagonalization).

### Analysis of Main Entry Points (`/home/prokop/git/FireCore/fortran/MAIN/`)

*   **`fireball.f90` (Standalone Program Driver):**
    *   **Overall Flow:**
        1.  **Initialization:**
            *   Calls `initbasics()` (from `INITIALIZERS`).
            *   Calls `readdata_mcweda()` (from `READFILES`) to read inputs and `Fdata`.
            *   Initializes wavefunctions, FIRE optimizer parameters.
        2.  **Outer Loop (Geometry Steps / MD Time Steps):**
            *   Sets `iforce = 1` for force calculation.
            *   Resets SCF convergence flags.
            *   **Inner Loop (SCF Cycle):** Iterates until charge self-consistency.
                1.  `assemble_mcweda()` (from `ASSEMBLERS`): Constructs Hamiltonian (`h_mat`) and overlap (`s_mat`) matrices using current charges and tabulated integrals.
                2.  `solveH(ikpoint, k_temp)`: Solves Hc = eSc for the current k-point.
                3.  `denmat()`: Calculates density matrix, occupation numbers, Fermi energy, and output charges (`Qout`).
                4.  Calculates charge convergence criterion (`sigma`).
                5.  Checks for SCF convergence.
                6.  `mixer()`: Mixes input (`Qin`) and output (`Qout`) charges for the next SCF iteration.
            *   **Post-SCF (within the outer loop step):**
                1.  `getenergy_mcweda()` (from `ASSEMBLERS`): Calculates total energy.
                2.  `getforces_mcweda()` (from `DASSEMBLERS`): Calculates atomic forces.
                3.  `move_ions_FIRE()` (from `FIRE` module, likely in `MAIN` or `INITIALIZERS`): If optimizing geometry, moves ions.
                4.  Checks for geometry convergence.
        3.  **Finalization:** Optional outputs, timing summary.

*   **`libFireCore.f90` (Library Interface):**
    *   **Purpose:** Exposes Fireball functionalities as C-bindable subroutines for use from Python, C++, etc.
    *   **Key Subroutines (mapping to `fireball.f90` logic):**
        *   **Initialization:**
            *   `firecore_preinit()`: Basic setup.
            *   `firecore_initdir()`: Full initialization from disk-based inputs.
            *   `firecore_init(natoms_, atomTypes, atomsPos)`: Initializes with geometry/types passed as arguments. Performs detailed setup similar to `fireball.f90` start.
            *   `firecore_set_lvs()`: Sets lattice vectors.
        *   **SCF Cycle Control (Granular):**
            *   `firecore_assembleH()`: Corresponds to `assemble_mcweda()`.
            *   `firecore_solveH()`: Corresponds to `solveH()`.
            *   `firecore_updateCharges()`: Corresponds to `denmat()` followed by `mixer()` and charge convergence check.
        *   **SCF Cycle Control (Combined):**
            *   `firecore_SCF(nmax_scf, positions_, iforce_)`: Encapsulates the entire SCF loop.
        *   **Force/Energy Evaluation:**
            *   `firecore_evalForce()`: Performs SCF, then calls `getenergy_mcweda()` and `getforces_mcweda()`.
        *   **Relaxation:**
            *   `firecore_relax()`: Implements a geometry relaxation loop.
        *   **Data Access:**
            *   `firecore_getPointer_*`: Provides C pointers to internal Fortran arrays (positions, forces, wavefunction coefficients, charges).
            *   `firecore_getCharges()`: Retrieves calculated atomic charges.
            *   `firecore_get_wfcoef()`, `firecore_set_wfcoef()`: Get/set wavefunction coefficients.
        *   **Grid Operations:**
            *   `firecore_setupGrid()`, `firecore_getGridMO()`, `firecore_getGridDens()`, `firecore_orb2xsf()`, `firecore_dens2xsf()`, `firecore_orb2points()`, `firecore_dens2points()`, `firecore_getpsi()`.

### Core SCF Subroutines (from `/home/prokop/git/FireCore/fortran/MAIN/`)

*   **`anderson2.f90` (called by `mixer`):**
    *   Implements Anderson mixing scheme to accelerate SCF convergence by extrapolating input charges based on a history of previous charge differences.

*   **`denmat.f90`:**
    *   Calculates the density matrix `rho(mu,nu,ineigh,iatom)` and energy-weighted density matrix `cape`.
    *   Sums contributions from occupied states: $\sum_j f_i c^*_{j,\mu} c_{j,\nu} \exp(i\mathbf{k}\cdot\mathbf{R})$.
    *   Calls `fermie()` to get occupation numbers `foccupy` and Fermi energy `efermi`.
    *   Calculates output charges `Qout` (Lowdin or Mulliken) based on `iqout`.

*   **`fermie.f90`:**
    *   Calculates the Fermi energy (`efermi`) by ensuring the sum of occupation numbers (`foccupy`) equals the total number of electrons (`ztot`), typically using a bisection method.
    *   Occupation numbers are determined by the Fermi-Dirac distribution.

*   **`ktransform.f90` (called by `solveH`):**
    *   Transforms the real-space, block-sparse Hamiltonian (`h_mat`, `vnl`) and overlap (`s_mat`) matrices to k-space for a given `kpoint`.
    *   $H_k(j\mu, j\nu) = \sum_l \exp(i\mathbf{k}\cdot\mathbf{R}_l) H_{real}(j\mu, j\nu)$.

*   **`mixer.f90`:**
    *   Orchestrates the charge mixing step of the SCF cycle.
    *   Calculates $\sigma = \sqrt{\sum (Q_{in} - Q_{out})^2}$.
    *   If not converged, calls a mixing algorithm (e.g., `anderson2`) to update `Qin`.
    *   Renormalizes `Qin` for charge conservation.

*   **`solveH.f90`:**
    *   Core routine for solving the electronic structure for a fixed potential.
    *   1.  Calls `ktransform()` to get S(k) and H(k).
    *   2.  Calls `sqrtS()` to compute S(k)<sup>-1/2</sup>.
    *   3.  Transforms Hamiltonian to an orthogonal basis (Lowdin transform): $H_k^{ortho} = S^{-1/2} H_k S^{-1/2}$.
    *   4.  Diagonalizes `Hk_ortho` (e.g., using LAPACK's `zheevd`) to get eigenvalues (`eigen_k`) and eigenvectors in the orthogonal basis.
    *   5.  Transforms eigenvectors back to the original non-orthogonal atomic orbital basis (`bbnkre`, `bbnkim`).

*   **`sqrtS.f90`:**
    *   Diagonalizes the k-space overlap matrix S(k).
    *   Identifies and handles linear dependencies (small eigenvalues).
    *   Constructs $S(k)^{-1/2}$ using the eigenvalues and eigenvectors of S(k).

### Mapping the Overall Code Execution Flow

1.  **Setup Phase (e.g., `firecore_init` or start of `fireball.f90`):**
    *   **Input Reading (`READFILES` directory):** Parse `options.input`, atomic coordinates, `Fdata` paths. Key routines: `readparam`, `readdata`.
    *   **Basic Initialization (`INITIALIZERS` directory):** Set up constants, species information (`readinfo`), orbital mappings (`make_munu*`), neighbor search parameters. Key routines: `initbasics`, `initconstants`, `initamat`.
    *   **Data Loading (`READFILES` directory):** Load integral tables from `Fdata` into arrays defined in the `integrals` module. Key routines: `readdata_2c`, `readdata_3c`, `read_1c`.
    *   **Neighbor Finding (`NEIGHBORS` directory):** Determine interacting atom pairs/triplets. Key routines: `neighbors`, `neighborsPP`.
    *   **Memory Allocation (`ALLOCATIONS` directory):** Allocate large arrays for H, S, rho, charges, wavefunctions (defined in `MODULES`). Key routines: `allocate_h`, `allocate_s`, `allocate_rho`.
    *   Initial guess for `Qin` (input charges, in `charges` module), often from `initcharges` (in `INITIALIZERS`).

2.  **SCF Cycle (Main loop in `fireball.f90` or controlled by `firecore_SCF`/`firecore_evalForce`):**
    *   **`assemble_mcweda` (from `ASSEMBLERS` directory):**
        *   Iterates through atom pairs/triplets.
        *   Retrieves pre-calculated integrals from `integrals` module.
        *   Uses `INTERPOLATERS` (`interpolate_2d`, `interpolate_1d`) for values at current geometry.
        *   Uses `ROTATIONS` (`rotate`, `makeDmat`) to transform integrals to global coordinates.
        *   Computes Hamiltonian contributions (kinetic, pseudopotential, XC, electrostatic based on `Qin`).
        *   Sums contributions into `h_mat` and `s_mat` (in `interactions` module).
        *   Involves specific assemblers like `assemble_2c`, `assemble_3c`, `assemble_olsxc_on`, `buildh`.
        *   `get_ewald` (from `INTERACTIONS`) for long-range electrostatics if applicable.
    *   **`solveH` (in `MAIN` directory):**
        *   Calls `ktransform` to get S(k), H(k).
        *   Calls `sqrtS` to get $S(k)^{-1/2}$.
        *   Diagonalizes to get `eigen_k` and wavefunction coefficients (`bbnkre`/`bbnkim` in `density` module).
    *   **`denmat` (in `MAIN` directory):**
        *   Calls `fermie` to get occupations `foccupy`.
        *   Calculates `rho` (density matrix in `density` module) and `Qout` (output charges in `charges` module).
    *   **`mixer` (in `MAIN` directory):**
        *   Compares `Qin` and `Qout`.
        *   If not converged, calls `anderson2` to update `Qin`.

3.  **Post-SCF (after convergence):**
    *   **`getenergy_mcweda` (from `ASSEMBLERS` directory):** Calculates total energy components (stored in `energy` module).
    *   **`getforces_mcweda` (from `DASSEMBLERS` directory):** Calculates atomic forces (`ftot` in `forces` module). Involves routines like `Dassemble_2c`, `Dassemble_3c`.
    *   **Optional:** Geometry optimization (e.g., `move_ions_FIRE` using routines related to the `FIRE` module) or Molecular Dynamics.
    *   **Optional:** Output to grid/files (using `GRID` module routines like `project_orb`, `writeout_xsf`).

This structure highlights a modular design where different tasks (initialization, assembly, solving, I/O) are handled by specific sets of routines and modules, facilitating the complex workflow of a DFTB calculation.


---

## Fireball Module Variable Descriptions

This document details key global variables within specified Fireball Fortran modules, their purpose, and primary usage context.

### `/home/prokop/git/FireCore/fortran/MODULES/options.f90`

This module contains global variables that control the behavior of a Fireball calculation. They are typically set by reading an `options.input` file or through the library interface.

*   **`iparam_file`**: If `1`, parameters are read from a file. If `0`, parameters are set by default or via API calls. Used during initialization in `readparam.f90` or `firecore_preinit`.
*   **`itheory`**: Selects the main electronic structure theory level.
    *   `0`: Harris functional (non-SCF or simplified SCF).
    *   `1`: DOGS (Density-Optimized Generalized Slater-Koster) SCF method.
    *   `2`: Extended Hubbard model.
    *   Primarily used in `assemble_mcweda` and `getenergy_mcweda` to dispatch to theory-specific routines.
*   **`itheory_xc`**: Selects the specific exchange-correlation (XC) functional approximation.
    *   `0`: Horsfield XC.
    *   `1`: Generalized Sankey-Niklewski (GSN) XC.
    *   `2`: McWEDA XC.
    *   Used in XC assembly routines like `assemble_olsxc_on`, `assemble_snxc_on` within `assemble_mcweda`.
*   **`igsn`**: Enables the Generalized Sankey-Niklewski XC theory. Checked in XC assembly parts of `assemble_mcweda`.
*   **`iharris`**: Enables the Harris functional theory. Checked in `assemble_mcweda` and energy/force calculation routines.
*   **`idogs`**: Enables the DOGS SCF method. Checked in `assemble_mcweda` for DOGS-specific terms (`vca`, `vxc_ca`).
*   **`imcweda`**: Enables the McWEDA SCF theory (often default). Checked in `assemble_mcweda` for McWEDA/OLSXC terms.
*   **`iks`**: Enables Kohn-Sham DFT calculations on a real-space grid. Used in `GRID` module routines.
*   **`igrid`**: Enables grid projection functionalities (e.g., for plotting densities/orbitals). Used by `GRID` module routines and `libFireCore.f90` grid functions.
*   **`iwrtewf`**: If greater than 0, enables writing out wavefunctions, typically to grid files. Checked in `fireball.f90` post-SCF and `libFireCore.f90` grid functions.
*   **`idipole`**: Enables inclusion of long-range dipole corrections to electrostatics. Affects `assemble_lr` and related force calculations.
*   **`ivec_2c`**: If `1`, use vectorized versions of interpolation routines for 2-center integrals (e.g., `interpolate_1d_vec`). Checked in `INTERACTIONS/doscentros*.f90`.
*   **`ivec_3c`**: If `1`, use vectorized versions of interpolation routines for 3-center integrals (e.g., `interpolate_2d_vec`). Checked in `INTERACTIONS/trescentros*.f90`.
*   **`i2dlin`**: If `1`, use bi-linear interpolation for 3-center integrals instead of bi-cubic. Checked in `INTERPOLATERS/interpolate_2d.f90`.
*   **`iforce`**: If `1`, calculate atomic forces. Checked in main SCF loops (`fireball.f90`, `libFireCore.f90`) to call force calculation routines like `getforces_mcweda`.
*   **`icluster`**: If `1`, perform a gas-phase (isolated molecule/cluster) calculation. If `0`, assume periodic boundary conditions (PBC). Affects neighbor finding, k-point usage, and Ewald sums.
*   **`iimage`**: Controls how often periodic images are considered or updated in PBC calculations. Used in `NEIGHBORS` and geometry updates.
*   **`iqout`**: Specifies the type of charge population analysis.
    *   `0` or `1`: Lowdin population analysis.
    *   `2`: Mulliken population analysis.
    *   Used in `MAIN/denmat.f90` to calculate `Qout`.
*   **`ifixcharge`**: If `1`, atomic charges (`Qin`) are kept fixed during the SCF cycle (disables charge self-consistency). Checked in `MAIN/mixer.f90`.
*   **`ifixneigh`**: If `1`, neighbor lists are kept fixed throughout the calculation (e.g., after the first step). Checked before calling neighbor list routines.
*   **`timing_verbosity`**: Controls the level of detail for timing information output. Used in `MAIN/fireball.f90` and potentially other main loops.
*   **`idebugWrite`**: General flag for enabling detailed debug output to files or console. Checked throughout the code.
*   **`verbosity`**: Overall verbosity level for program output (e.g., `0` for minimal, `10` for maximum). Checked in many print statements.
*   **`ntpr`**: Frequency for printing information during a run (e.g., print every `ntpr` MD steps or SCF iterations).
*   **`restartxyz`**: If `1`, the simulation attempts to start from a `restart.xyz` file. Checked during initialization.
*   **`inputxyz`**: If `1`, the input coordinate file is assumed to be in `.xyz` format. Checked in `READFILES/readdata*.f90`.
*   **`itestrange`**: Enables diagnostic tests related to interaction ranges or cutoffs.
*   **`iconstraints`**: Used to apply specific constraints during calculations (e.g., fixing atom positions, cell shape).
*   **`ioff2c`**: Diagnostic flags to selectively turn off different types of 2-center integral contributions. Checked in `ASSEMBLERS/assemble_2c*.f90`.
*   **`ioff3c`**: Diagnostic flags to selectively turn off different types of 3-center integral contributions. Checked in `ASSEMBLERS/assemble_3c*.f90`.
*   **`testrange`**: A value used in conjunction with `itestrange` for range-dependent diagnostics.
*   **`iwrtxyz`**: Controls the writing of XYZ coordinate files during a simulation (e.g., every N steps). Checked in geometry update routines.

### `/home/prokop/git/FireCore/fortran/MODULES/interactions.f90`

This module defines arrays and parameters related to atomic orbitals, their interactions, and the construction of Hamiltonian and overlap matrices.

*   **`wrtout`**: If `.TRUE.`, enables extensive write output for debugging purposes, often for matrix elements or intermediate interaction values.
*   **`smt_elect`**: Smoothing parameter used in Ewald summation or other long-range electrostatic calculations to partition real and reciprocal space contributions.
*   **`smt_vnl`**: Smoothing parameter used for non-local pseudopotential (VNL) terms, potentially for range separation or damping.
*   **`norbitals`**: Total number of atomic basis orbitals in the system. Defines the primary dimension of H and S matrices. Set in `INITIALIZERS/initamat.f90`.
*   **`norbitals_new`**: Actual number of orbitals used after handling linear dependencies (e.g., by removing orbitals with very small eigenvalues in the overlap matrix). Modified in `MAIN/sqrtS.f90`.
*   **`nbands`**: Total number of electronic bands (eigenstates) solved for. Typically equal to `norbitals_new`.
*   **`nsh_max`**: Maximum number of shells (s, p, d, etc.) allowed on a single atom type. Used for dimensioning some Fdata reading arrays.
*   **`ideriv_max`**: Maximum derivative order of integrals read from `Fdata` files (e.g., for force calculations). Set in `READFILES/readheader_*.f90`.
*   **`isorpmax`**: Maximum shell index (related to s, p, d character) encountered in the `Fdata` files for general integrals. Set in `READFILES/readheader_*.f90`.
*   **`isorpmax_xc`**: Maximum shell index specifically for exchange-correlation integrals from `Fdata`.
*   **`numorb_max`**: Maximum number of atomic orbitals on any single atom type in the basis set. Used for dimensioning fixed-size parts of `h_mat`, `s_mat`, etc. Calculated in `INITIALIZERS/initamat.f90`.
*   **`ME2c_max`**: Maximum number of unique two-center matrix element types (e.g., ss-sigma, sp-sigma, pp-sigma, pp-pi) between any pair of atom types. Used for dimensioning arrays holding 2c integral tables. Calculated in `INITIALIZERS/make_munu.f90`.
*   **`ME2cPP_max`**: Similar to `ME2c_max`, but specifically for two-center pseudopotential matrix elements.
*   **`ME2cDipY_max`**, **`ME2cDipX_max`**: Max unique 2c matrix elements for Y and X dipole components, if dipole interactions are included.
*   **`ME3c_max`**: Maximum number of unique three-center matrix element types. Used for dimensioning arrays holding 3c integral tables. Calculated in `INITIALIZERS/make_munu.f90`.
*   **`MES_max`**: Maximum number of matrix elements when using a spherical density approximation (e.g., in OLSXC). Calculated in `INITIALIZERS/make_munuS.f90`.
*   **`index_max2c`**: Maps `(species_i, species_j)` to the maximum index/count of 2-center integral types for that pair. Used in `ASSEMBLERS` to loop over relevant integrals.
*   **`index_max3c`**: Maps `(species_i, species_j)` to the maximum index/count of 3-center integral types involving these two species and a third.
*   **`index_max2cDipY`**, **`index_max2cDipX`**: Similar to `index_max2c` for Y and X dipole integrals.
*   **`lssh`**: Stores the angular momentum quantum number (s=0, p=1, d=2) for each shell of each atomic species `(shell_index, species_index)`. Used extensively for integral selection and rotation.
*   **`mu`**, **`nu`**: Map `(species_index, shell_index, m_value_sub_index)` to a local orbital index within that species. `mu` and `nu` are used to index the two orbitals involved in a matrix element. Crucial for constructing `h_mat` and `s_mat`.
*   **`mvalue`**: Stores the actual magnetic quantum number (m-value) for each orbital defined by `(species_index, shell_index, m_value_sub_index)`. Used in rotations.
*   **`nssh`**: Stores the number of shells for each atomic species.
*   **`nssh_tot`**: Total number of shells summed over all species (used for some Fdata indexing).
*   **`num_orb`**: Stores the total number of orbitals for each atomic species.
*   **`num_orb_sh`**: Stores the number of orbitals per shell type (s, p, d). E.g., `num_orb_sh(1)=1` (s), `num_orb_sh(2)=3` (p).
*   **`muDipY`**, **`nuDipY`**, **`muDipX`**, **`nuDipX`**: Orbital index mappings similar to `mu`/`nu` but for dipole integrals.
*   **`getmssh`**, **`getlssh`**, **`getissh`**, **`getiatom`**: These provide a flat mapping from a global orbital index (1 to `norbitals`) to its m-value, l-value, shell index, and parent atom index, respectively. Populated in `INITIALIZERS/get_info_orbital.f90`.
*   **`degelec`**: Stores the starting global orbital index for each atom in the system. `degelec(iatom)` gives the index of the first orbital belonging to `iatom`. Essential for addressing blocks in `h_mat`, `s_mat`, and wavefunctions.
*   **`mu2shell`**: Maps a local orbital index `mu` (within a species) to its shell index. Used for extended Hubbard interactions.
*   **`index_maxPP`**, **`lsshPP`**, **`muPP`**, **`nuPP`**, **`nsshPP`**, **`num_orbPP`**: Arrays analogous to their non-PP counterparts, but specifically for Kleinmann-Bylander pseudopotential projectors.
*   **`index_maxS`**, **`muS`**, **`nuS`**, **`mvalueS`**: Arrays analogous to their non-S counterparts, but for the basis functions used in spherical density approximations (e.g., OLSXC).
*   **`iatyp`**: Stores the species index (type) for each atom `iatom` in the system. Links `ratom(iatom)` to `Fdata` for that species.
*   **`imass`**: Stores the atomic mass for each *species type*. Used in MD or dynamics.
*   **`cl_PP`**: Stores the coefficients for Kleinmann-Bylander pseudopotential projectors.
*   **`h_mat`**: The Hamiltonian matrix elements, stored in a block-sparse format: `(orbital_index_on_atom_i, orbital_index_on_atom_j, atom_index_k, neighbor_index_l_of_k)`. Assembled in `ASSEMBLERS`.
*   **`s_mat`**: The overlap matrix elements, with the same block-sparse structure as `h_mat`. Assembled in `ASSEMBLERS`.
*   **`sm_mat`**: Potentially a smoothed version of the overlap matrix, used for Ewald sums or other long-range corrections.
*   **`sVNL`**: Matrix elements of the overlap with non-local pseudopotential projectors $\langle\phi_{\mu} | V_{nl} | \phi_{\nu}\rangle$. Assembled in `ASSEMBLERS/assemble_sVNL.f90`.
*   **`t_mat`**: Kinetic energy matrix elements. Assembled in `ASSEMBLERS`.
*   **`vna`**: Nuclear attraction potential matrix elements (typically two-center). Assembled in `ASSEMBLERS`.
*   **`vnl`**: Non-local pseudopotential matrix elements. Assembled in `ASSEMBLERS`.
*   **`vnl2c`**, **`vnl3c`**: Two-center and three-center components of the non-local pseudopotential, if separated.
*   **`vxc`**: Exchange-correlation potential matrix elements, often dominated by three-center contributions in McWEDA/OLSXC. Assembled in `ASSEMBLERS`.
*   **`vxc_1c`**: One-center (on-site) exchange-correlation potential matrix elements. Assembled in `ASSEMBLERS`.
*   **`dip`**: Dipole matrix elements, used for DOGS theory or explicit long-range dipole interactions.
*   **`ewald`**: The Ewald sum matrix `(atom_i, atom_j)` representing long-range electrostatic interactions between atoms. Calculated in `INTERACTIONS/get_ewald.f90`.
*   **`ewaldlr`**, **`ewaldsr`**: Long-range and short-range parts of Ewald sum contributions to `h_mat`, if orbital-resolved.
*   **`vca`**: Charge-dependent potential correction matrix elements, specific to DOGS theory. Assembled in `ASSEMBLERS`.
*   **`vxc_ca`**: Charge-dependent exchange-correlation correction matrix elements, specific to DOGS theory. Assembled in `ASSEMBLERS`.
*   **`dipcm`**, **`dipc`**: Arrays for dipole matrix elements if `idipole` is active, potentially in condensed or detailed forms.
*   **`xc_overtol`**: A tolerance value used in XC calculations, possibly for integral cutoffs or convergence.
*   **`wavefxn`**, **`napot`**: Store filenames for wavefunctions and neutral atom potentials when using grid-based functionalities (`IF_DEF_GRID`).

### `/home/prokop/git/FireCore/fortran/MODULES/integrals.f90`

This module declares arrays that store the pre-calculated integral values read from the `Fdata` files. These are the raw numerical tables used for interpolation.

*   **`fdataLocation`**: Stores the path to the directory containing the `Fdata` files. Set in `READFILES/findFdata.f90`.
*   **`nsup`**, **`nsu`**: Related to skipping certain atom species during `Fdata` reading (likely a legacy feature).
*   **`icon3c`**: An index mapping `(species_i, species_j, species_k)` to a specific table index for three-center integrals. Used to locate the correct 3c table in `Fdata`.
*   **`numx3c_bcna`**, **`numy3c_bcna`**: Store the number of grid points in the x and y dimensions for tabulated three-center "bare core neutral atom" (bcna) integrals for each species combination. Read from `Fdata` headers.
*   **`hx_bcna`**, **`hy_bcna`**: Store the grid spacing in x and y for 3c bcna integrals.
*   **`x3cmax_bcna`**, **`y3cmax_bcna`**: Store the maximum range in x and y for 3c bcna integral tables.
*   **`bcna_01`** to **`bcna_05`**: Store the tabulated values of 3c bcna integrals. The dimensions are typically `(grid_x_index, grid_y_index, ME_type_index, derivative_type_index, species_combination_index)`. The `_01` to `_05` usually correspond to different theta angle slices of the 3D integral. Read in `READFILES/read_3c.f90`.
*   **`numx3c_xc3c`**, **`numy3c_xc3c`**, **`hx_xc3c`**, **`hy_xc3c`**, **`x3cmax_xc3c`**, **`y3cmax_xc3c`**: Analogous to `_bcna` counterparts, but for three-center exchange-correlation (`xc3c`) integrals.
*   **`xc3c_01`** to **`xc3c_05`**: Store tabulated 3c xc3c integrals.
*   **`numx3c_den3`**, **`numy3c_den3`**, **`hx_den3`**, **`hy_den3`**, **`x3cmax_den3`**, **`y3cmax_den3`**: Analogous, but for three-center density (`den3`) integrals, primarily used in the SNXC method.
*   **`den3_01`** to **`den3_05`**: Store tabulated 3c den3 integrals (for SNXC).
*   **`den3S_01`** to **`den3S_05`**: Store tabulated 3c integrals related to spherical density components (used in OLSXC).
*   **`bcna_vec`**, **`xc3c_vec`**, **`den3_vec`**, **`den3S_vec`**: Potentially alternative, vectorized, or restructured storage for the corresponding `_0x` integral tables.
*   **`bcna_t`**, **`xc3c_t`**, **`den3_t`**, **`den3S_t`**: Arrays of derived type `t_data3c`. This type likely encapsulates the grid parameters and data for a specific 3c integral table, providing a more structured way to handle them.
*   **`ind2c`**: A fixed-size mapping table. `ind2c(interaction_type, derivative_type)` gives an index or offset into the 2-center integral tables. `interaction_type` could be overlap, kinetic, Vna, Vnl, etc. `derivative_type` is 0 for value, 1 for 1st derivative, etc.
*   **`numz2c`**: Stores the number of grid points (along the interatomic distance z) for two-center integrals for each `(species_i, species_j, integral_type_index_from_ind2c)`. Read from `Fdata` headers.
*   **`xintegral_2c`**: Stores the tabulated values of 2-center integrals: `(grid_z_index, ME_type_index, derivative_type_index, species_i, species_j)`. Read in `READFILES/read_2c.f90`.
*   **`splineint_2c`**: Stores the spline coefficients for 2-center integrals, allowing for smooth interpolation: `(spline_coeff_index, grid_z_index, ME_type_index, derivative_type_index, species_i, species_j)`. Calculated in `INTERPOLATERS/buildspline_1d.f90`.
*   **`z2cmax`**: Stores the maximum interatomic distance (range) for which 2-center integrals are tabulated for `(species_i, species_j, integral_type_index)`.
*   **`exc1c_0`**: Stores the base one-center exchange-correlation energy for `(species, shell_type)`. Read in `READFILES/read_1c.f90`.
*   **`exc1c`**: Stores more detailed one-center XC energy terms, potentially dependent on charge state or derivatives: `(species, shell_type_1, shell_type_2, charge_index_or_derivative_type)`.
*   **`xcnu1c`**: Stores one-center XC potential terms: `(species, shell_type, charge_index_or_derivative_type)`.
*   **`xcnu1cs`**: One-center XC potential terms specifically for spherical density approximations (OLSXC).
*   **`exc1c0`**, **`nuxc1c`**, **`dexc1c`**, **`d2exc1c`**, **`dnuxc1c`**, **`d2nuxc1c`**: Arrays for alternative or derivative forms of one-center XC energies and potentials, likely used for more advanced XC treatments or force calculations.

### `/home/prokop/git/FireCore/fortran/MODULES/forces.f90`

This module contains arrays for storing atomic forces and derivatives of interaction terms needed for force calculations.

*   **`dewald`**: Stores derivatives of Ewald sum terms with respect to atomic positions `(atom_i, atom_j, cartesian_component_xyz)`. Used in `DASSEMBLERS` to calculate Ewald forces.
*   **`fewald`**: Stores the Ewald forces acting on each atom `(atom_i, cartesian_component_xyz)`.
*   **`flrew`**: Stores the long-range component of Ewald forces on atoms.
*   **`dipp`**, **`dippcm`**, **`dippc`**: Store derivatives of various dipole matrix elements, used if dipole interactions or corrections are active.
*   **`sp_mat`**: Stores derivatives of the overlap matrix elements ($dS/dR$), often called Pulay terms. Structure mirrors `s_mat`. Calculated in `DASSEMBLERS`.
*   **`spm_mat`**: Derivatives of the smoothed overlap matrix.
*   **`spVNL`**: Derivatives of the `sVNL` matrix elements.
*   **`tp_mat`**: Derivatives of the kinetic energy matrix elements ($dT/dR$).
*   **`ftot`**: Stores the total calculated forces on each atom `(cartesian_component_xyz, atom_i)`. This is the primary output for forces.
*   **`ftotold`**: Stores forces from the previous ionic step, used in some geometry optimizers (e.g., FIRE).
*   **`ftotnew`**: Stores newly calculated forces in the current ionic step before mixing or use by an optimizer.
*   **`dusr`**: Stores the derivative of the short-range repulsive potential ($U_{SR}$) with respect to atomic positions. Contributes to `ftot`.
*   **`dxcv`**: Stores the derivative of the XC correction term ($\delta U_{xc}$ from Harris functional) with respect to atomic positions. Contributes to `ftot`.
*   **`fro`**: Stores force contributions arising from the density matrix and overlap derivatives (Pulay forces related to `sp_mat * rho`).
*   **`ft`**: Stores force contributions from the kinetic energy term (related to `tp_mat`).
*   **`deltaF`**: Represents the change in forces or the forces used directly by a geometry optimization algorithm.
*   **`fana`**, **`fanl`**, **`faxc`**: Store force contributions from the nuclear attraction ($V_{na}$), non-local ($V_{nl}$), and exchange-correlation ($V_{xc}$, often 3-center) potential terms, respectively. These are typically Hellmann-Feynman type terms.
*   **`fotna`**, **`fotnl`**, **`fotxc`**: Store "on-top" force contributions, which are derivatives of the potential matrix elements themselves (e.g., $d(V_{na})/dR$).
*   **`f3naa`**, **`f3nab`**, **`f3nac`**: Store components of 3-center nuclear attraction forces, distributed among the three atoms (a, b, c) involved.
*   **`f3nla`**, **`f3nlb`**, **`f3nlc`**: Similar to `f3naa` etc., but for 3-center non-local forces.
*   **`f3xca`**, **`f3xcb`**, **`f3xcc`**: Similar, for 3-center exchange-correlation forces.
*   **`dxcdcc`**: Derivative of the XC double counting correction term, specifically for OLSXC forces.
*   **`arhop_off`**, **`arhopij_off`**, **`rhop_off`**, **`rhopij_off`**, **`rhop_on`**, **`arhop_on`**: Derivatives of various on-site and off-site (average) density terms used in calculating OLSXC forces.
*   **`faca`**, **`faxc_ca`**: Force contributions from charge-dependent potential ($V_{ca}$) and XC ($V_{xc\_ca}$) terms in DOGS theory.
*   **`fotca`**, **`fotxc_ca`**: "On-top" force contributions from $V_{ca}$ and $V_{xc\_ca}$ in DOGS theory.
*   **`f3caa`**, **`f3cab`**, **`f3cac`**: Components of 3-center $V_{ca}$ forces (DOGS).
*   **`f3xca_ca`**, **`f3xcb_ca`**, **`f3xcc_ca`**: Components of 3-center $V_{xc\_ca}$ forces (DOGS).

### `/home/prokop/git/FireCore/fortran/MODULES/density.f90`

This module holds arrays related to the electronic density, wavefunctions, and derived quantities like the Fermi energy and occupation numbers.

*   **`bbnkre`**: Stores the real part of the wavefunction coefficients (eigenvectors) in the non-orthogonal atomic orbital (NAO) basis: `(orbital_index, band_index, kpoint_index)`. Calculated in `MAIN/solveH.f90`.
*   **`bbnkim`**: Stores the imaginary part of the wavefunction coefficients in the NAO basis.
*   **`blowim`**, **`blowre`**: Store the imaginary and real parts of wavefunction coefficients in the Lowdin (orthogonalized, $S^{1/2}$ transformed) basis. Used for Lowdin population analysis in `MAIN/denmat.f90`.
*   **`eigen_k`**: Stores the eigenvalues (electronic energy levels) for each band and k-point: `(band_index, kpoint_index)`. Calculated in `MAIN/solveH.f90`.
*   **`cape`**: The energy-weighted density matrix (EWDW), $P^E_{\mu\nu} = \sum_i f_i \epsilon_i c^*_{i\mu} c_{i\nu}$. Used for calculating the band structure energy component. Structure mirrors `h_mat`. Calculated in `MAIN/denmat.f90`.
*   **`rho`**: The density matrix, $P_{\mu\nu} = \sum_i f_i c^*_{i\mu} c_{i\nu}$. Structure mirrors `h_mat`. Calculated in `MAIN/denmat.f90`.
*   **`rhoPP`**: A density matrix related to pseudopotential projectors, if a specific formulation requires it.
*   **`rhoA`**: (IF_DEF_GRID). Stores atomic density contributions, e.g., a sum of neutral atom densities, possibly on a grid or per atom. Used in Harris functional or as initial guess. Populated in `GRID/project_dens0.f90`.
*   **`rho_old`**: (IF_DEF_GRID). Stores the density matrix from the previous SCF iteration. Used in density mixing schemes.
*   **`rho_on`**, **`rhoi_on`**: On-site total density and individual atomic density components used in the OLSXC average density calculation.
*   **`rho_off`**, **`rhoij_off`**: Off-site total density and pair (atom i, atom j) density components for OLSXC.
*   **`arho_on`**, **`arhoi_on`**: On-site "average rho" quantities derived from `rho_on` and `rhoi_on` for OLSXC.
*   **`arho_off`**, **`arhoij_off`**: Off-site "average rho" quantities for OLSXC.

### `/home/prokop/git/FireCore/fortran/MODULES/energy.f90`

This module defines variables for storing various components of the total energy of the system.

*   **`atomic_energy`**: Sum of the energies of isolated neutral atoms. Used as a reference to calculate binding or cohesive energy. Calculated once at initialization.
*   **`etot`**: The total energy of the system for the current configuration. This is the primary energy output. Sum of various components calculated in `ASSEMBLERS/getenergy_mcweda.f90`.
*   **`ebs`**: The band structure energy, typically calculated as $\sum_i f_i \epsilon_i$, where $f_i$ are occupation numbers and $\epsilon_i$ are eigenvalues. Calculated in `ASSEMBLERS/getenergy_mcweda.f90`.
*   **`etotnew`**: Stores the newly calculated total energy in the current step/iteration.
*   **`etotold`**: Stores the total energy from the previous step/iteration. Used for convergence checks or in optimizers.
*   **`etotper`**: Total energy per atom (`etot / natoms`).
*   **`etotxc_1c`**: Contribution to the total energy from one-center (on-site) exchange-correlation terms. Calculated in `ASSEMBLERS/getenergy_mcweda.f90`.
*   **`getot`**: Often represents the Gibbs free energy if entropic terms are included, or can be another variant of the total energy.
*   **`getot_initial`**: Initial value of `getot`.
*   **`getotper`**: `getot` per atom.
*   **`deltaE`**: The change in total energy between successive iterations or steps (e.g., `etotnew - etotold`). Used for SCF or geometry convergence.
*   **`deltaFmax`**: The maximum force component on any atom. Used as a convergence criterion in geometry optimization. Updated in `FIRE/move_ions_FIRE.f90`.
*   **`uiiuee`**: Represents the classical electrostatic energy, including ion-ion repulsion and the electron-electron Hartree energy double-counting term ($1/2 \iint \rho(r)\rho(r')/|r-r'| drdr' - \int V_H[\rho]\rho(r)dr$). Calculated in `ASSEMBLERS/getenergy_mcweda.f90`.
*   **`uxcdcc`**: The exchange-correlation double-counting correction energy. This term corrects for using $\mu_{xc}[\rho_{in}]$ in the Hamiltonian construction while the XC energy functional is $E_{xc}[\rho_{in}]$. Specifically, $\int (E_{xc}[\rho_{in}]/\rho_{in} - \mu_{xc}[\rho_{in}]) \rho_{in} dr$. Calculated in `ASSEMBLERS/getenergy_mcweda.f90`.
*   **`uxcdcc_sn`**: The XC double-counting correction specific to the Sankey-Niklewski XC scheme.
*   **`uxcdcc_ols`**: The XC double-counting correction specific to the OLSXC scheme.
*   **`uxcdcc_hf`**: The XC double-counting correction specific to the Horsfield XC scheme.



## Fireball Default Parameters from `set_default_params`

This section outlines the default parameters set in the `set_default_params` subroutine within `/home/prokop/git/FireCore/fortran/READFILES/readparam.f90` and their implications for a Fireball calculation if not overridden by an input file.

### Electronic Theory and XC
*   **`nkpoints = 1`**: By default, calculations are performed at the Gamma point (k=0,0,0) only. Suitable for molecules/large supercells; for periodic solids, change via `input.kpts` or API.
*   **`ivec_2c = 0`**, **`ivec_3c = 0`**: Vectorized 2-center and 3-center integral interpolation routines are disabled by default; scalar versions used.
*   **`i2dlin = 0`**: Bi-cubic interpolation (more accurate, potentially slower) is used for 3-center integrals by default, not bi-linear.
*   **`iharris = 0`**: The Harris functional (non-SCF or simplified SCF) is **not** the default theory.
*   **`idogs = 1`**: The DOGS (Density-Optimized Generalized Slater-Koster) SCF method is **enabled** by default (implies `itheory=1`).
*   **`imcweda = 1`**: The McWEDA XC scheme is **enabled** by default (implies `itheory_xc=2`).
*   **`iks = 0`**: Kohn-Sham DFT calculations on a real-space grid are disabled by default.
*   **`igrid = 0`**: Grid projection functionalities (e.g., for plotting densities/orbitals) are disabled by default.
*   **`iwrtewf = 0`**: Wavefunctions will not be written out to grid files by default.
*   **`igsn = 0`**: The Generalized Sankey-Niklewski (GSN) XC scheme is disabled by default.
*   **`iqout = 1`**: Lowdin population analysis is used by default to determine atomic charges.
*   **`qstate = 0.0d0`**: The system is assumed to be neutral by default.

### System and SCF Control
*   **`icluster = 1`**: Calculations are for isolated molecules/clusters (gas-phase) by default; PBC off.
*   **`ifixcharge = 0`**: Atomic charges are allowed to change and achieve self-consistency during SCF.
*   **`ifixneigh = 0`**: Neighbor lists will be updated as needed.
*   **`iimage = 0`**: Atoms are not re-imaged into the central cell (relevant for PBC, but default is `icluster=1`).
*   **`iconstraints(1)=0, iconstraints(2)=1, iconstraints(3)=1, iconstraints(4)=1`**: Default constraints: Center of mass NOT constrained; Total momentum IS constrained; Kinetic energy IS scaled to initial temperature; Total angular momentum IS constrained.
*   **`basisfile = 'input.bas'`**: Atomic coordinates/species info read from `input.bas`.
*   **`lvsfile = 'input.lvs'`**: If PBC used (`icluster=0`), lattice vectors read from `input.lvs`.
*   **`kptpreference = 'input.kpts'`**: If PBC used, k-point info read from `input.kpts`.
*   **`fdataLocation = 'Fdata'`**: Pre-calculated integral tables expected in `./Fdata/`.
*   **`nstepf = 100`**: Default number of MD/geometry optimization steps is 100.
*   **`force_tol = 1.0E-4`**: Default force convergence criterion for geometry optimization is $10^{-4}$ (Ha/Bohr or eV/Ang).
*   **`dt = 0.5`**: Default MD timestep is 0.5 fs.
*   **`max_scf_iterations = 200`**: Maximum of 200 SCF iterations.
*   **`bmix = 0.04d0`**: Default linear charge mixing parameter is 0.04 (small, suitable for Anderson mixing).
*   **`sigmatol = 1.0E-5`**: SCF converged when RMS charge difference (`sigma`) < $10^{-5}$.
*   **`tempfe = 100.0d0`**: Fermi temperature for smearing occupations is 100 K.
*   **`rescal = 1.0d0`**: No rescaling of positions, lattice vectors, or k-points by default.
*   **`xyz2line = 2`**: Second line of `answer.xyz` likely contains energy, temperature, and time.

### Debugging and I/O
*   **`verbosity = 0`**: Minimal console output.
*   **`idebugwrite = 0`**: Detailed debug outputs disabled.
*   **`timing_verbosity = 0`**: Minimal timing information printed.
*   **`ntpr = 1`**: Information printed every 1 step (MD/geometry optimization).
*   **`restartxyz = 0`**: Simulation will not attempt to start from `restart.xyz`.
*   **`inputxyz = 0`**: Primary input coordinate file is `input.bas`, not generic `.xyz`.
*   **`iwrtxyz = 0`**: XYZ coordinate files (e.g., `answer.xyz`) not written out by default during simulation.

**Summary of Default Behavior:**

By default, Fireball is set up to run a **DOGS SCF calculation with McWEDA XC** for an **isolated molecule/cluster** at the **Gamma point**. It will attempt to achieve charge self-consistency with a tolerance of $10^{-5}$ using Anderson mixing (implied by small `bmix`) and a Fermi smearing of 100K. Forces will be calculated, and a geometry optimization or MD run will proceed for 100 steps with a force tolerance of $10^{-4}$ Ha/Bohr (or equivalent). Output will be minimal. This setup is a reasonable starting point for many molecular systems. For periodic solids, `icluster`, `lvsfile`, and `kptpreference` would need to be adjusted.
