## Hamiltonian Assembly in Fireball (`assemble_mcweda.f90`)

The `assemble_mcweda.f90` subroutine is the central routine for constructing the Hamiltonian (`h_mat`) and Overlap (`s_mat`) matrices in Fireball for the McWEDA/DOGS theories. It orchestrates calls to various subroutines that calculate different interaction terms based on atomic positions, orbital types, pre-calculated integral tables, and the current self-consistent charge distribution.

### Key Hamiltonian/Overlap Matrix Components (Variables in module `interactions.f90` ):
The final Hamiltonian $H$ is constructed as a sum of several terms:
$H = T + V_{na} + V_{nl} + V_{xc} + V_{xc\_1c} + V_{ca} + V_{xc\_ca} + V_{Ewald}$

*   $T$ (**`t_mat`**): Kinetic energy matrix elements, $T_{\mu\nu} = \langle\phi_{\mu} | \hat{T} | \phi_{\nu}\rangle$.
*   $V_{na}$ (**`vna`**): Nuclear attraction potential matrix (interaction with atomic cores), $V_{\mu\nu}^{na} = \langle\phi_{\mu} | \sum_I V_{core,I} | \phi_{\nu}\rangle$.
*   $V_{nl}$ (**`vnl`**): Non-local pseudopotential matrix, $V_{\mu\nu}^{nl} = \langle\phi_{\mu} | \hat{V}_{nl} | \phi_{\nu}\rangle$.
*   $V_{xc}$ (**`vxc`**): Multi-center (off-site) exchange-correlation potential matrix, $V_{\mu\nu}^{xc,mc}$.
*   $V_{xc\_1c}$ (**`vxc_1c`**): One-center (on-site) exchange-correlation potential matrix, $V_{\mu\nu}^{xc,1c}$.
*   $V_{ca}$ (**`vca`**): Charge-dependent potential correction (DOGS theory), $V_{\mu\nu}^{ca}$.
*   $V_{xc\_ca}$ (**`vxc_ca`**): Charge-dependent exchange-correlation correction (DOGS theory), $V_{\mu\nu}^{xc,ca}$.
*   $V_{Ewald}$ (contributions from **`ewaldlr`**, **`ewaldsr`**): Orbital-resolved long-range electrostatic contributions derived from the Ewald sum.

**Other relevant terms and their connection to the Hamiltonian:**
*   **`s_mat`**: Overlap matrix $S_{\mu\nu} = \langle\phi_{\mu} | \phi_{\nu}\rangle$. This is not part of $H$ itself but is essential for solving the generalized eigenvalue problem $Hc = eSc$.
*   **`sVNL`**: Overlap with non-local pseudopotential projectors. These are intermediate quantities used in the calculation of the `vnl` matrix elements.
*   **`ewald`**: Atom-pair Ewald sum matrix. This matrix of classical electrostatic interactions between atomic point charges is used to derive the orbital-resolved `ewaldlr` and `ewaldsr` contributions that are added to `h_mat`.
*   **`h_mat`**: This variable stores the final, total Hamiltonian matrix, which is the sum of all the above potential and kinetic energy terms.

### Common Acronyms and Terms:

*   **XC**: Exchange-Correlation.
*   **NA**: Neutral Atom (referring to properties derived from isolated, neutral atoms).
*   **T**: Kinetic Energy.
*   **Vna**: Nuclear Attraction potential (interaction with atomic cores).
*   **Vnl**: Non-Local pseudopotential.
*   **SCF**: Self-Consistent Field.
*   **Kscf**: Index of the current SCF iteration.
*   **McWEDA**: Multi-center Weighted Exchange-correlation Density Approximation.
*   **OLSXC**: Orthogonalized Linear combination of Slater-type orbitals for Exchange and Correlation (a key part of McWEDA for 3-center XC).
*   **DOGS**: Density-Optimized Generalized Slater-Koster (an SCF theory level).
*   **SNXC**: Sankey-Niklewski Exchange-Correlation (an older XC scheme).
*   **PP**: PseudoPotential.
*   **`Qin`**: Input atomic charges for the current SCF iteration (from `charges` module).
*   **`Fdata`**: Directory/files containing pre-calculated integral tables.

The assembly process can be divided into steps performed only once at the beginning of the SCF loop (`Kscf = 1`) and steps performed in every SCF iteration.

### Operations Performed ONCE per SCF Loop (`Kscf = 1`)

These steps primarily depend on the system's geometry and the pre-calculated `Fdata` integrals. They establish the non-charge-dependent parts of the Hamiltonian.

1.  **Neighbor List Setup:**
    *   If `ifixneigh` is 0 (default), neighbor lists are recalculated.
    *   Calls `reallocate_neigh`: Allocates memory for neighbor list arrays.
    *   Calls `neighbors`: Finds atom pairs within the cutoff radius.
    *   Calls `neighborsPP`: Finds atom pairs for pseudopotential interactions.
    *   Calls `num_neigh_tot`: Calculates the total count of neighbors.
    *   Calls `backnay`: Sets up reverse mapping for neighbor lists.
    *   Calls `neighbors_pairs`: Finds pairs of neighbors (relevant for 3-center terms).
    *   Calls `common_neighbors`: Finds atoms that are common neighbors to a pair (for 3-center terms).
    *   Calls `common_neighborsPP`: Finds common neighbors for pseudopotential 3-center terms.
    *   Calls `check_neighbors`: Verifies the consistency of neighbor lists.
    *   *Variables Populated*: Internal neighbor list arrays (not in provided modules).
    *   *Dependency*: Geometry only.

2.  **Initial Ewald Sum (if `itheory` is 1 or 2):**
    *   Calls `get_ewald`: Calculates the Ewald sum for long-range electrostatic interactions. This is typically based on initial point charges (e.g., neutral atoms or the initial guess for `Qin`).
    *   *Variables Populated*: `ewald`.
    *   *Dependency*: Geometry and initial charges (`Qin` at Kscf=1).

3.  **Base 2-Center Assembly:**
    *   Calls `assemble_sVNL`: Assembles overlap matrix elements with non-local pseudopotential projectors.
        *   *Variables Populated*: `sVNL`.
        *   *Equation/Purpose*: Calculates $\langle\phi_{\mu} | V_{nl} | \phi_{\nu}\rangle$ for the overlap matrix contribution from non-local projectors.
    *   Calls `assemble_2c`: Assembles two-center kinetic energy and nuclear attraction potential matrix elements.
        *   *Variables Populated*: `t_mat`, `vna`.
        *   *Equation/Purpose*: Calculates $\langle\phi_{\mu} | \hat{T} | \phi_{\nu}\rangle$ and $\langle\phi_{\mu} | V_{na} | \phi_{\nu}\rangle$ (2-center parts).
    *   Calls `assemble_2c_PP`: Assembles two-center pseudopotential matrix elements.
        *   *Variables Populated*: `vnl`.
        *   *Equation/Purpose*: Calculates $\langle\phi_{\mu} | V_{nl} | \phi_{\nu}\rangle$ (2-center part).
    *   *Dependency*: Geometry, Fdata tables.

4.  **Base 3-Center Assembly:**
    *   Calls `assemble_3c`: Assembles three-center contributions to the Hamiltonian (e.g., parts of kinetic, nuclear attraction, or Vnl).
        *   *Variables Populated*: Adds to `t_mat`, `vna`, `vnl`.
        *   *Equation/Purpose*: Calculates $\langle\phi_{\mu} | \hat{T} | \phi_{\nu}\rangle$, $\langle\phi_{\mu} | V_{na} | \phi_{\nu}\rangle$, and $\langle\phi_{\mu} | V_{nl} | \phi_{\nu}\rangle$ (3-center parts).
    *   Calls `assemble_3c_PP`: Assembles three-center pseudopotential matrix elements.
        *   *Variables Populated*: Adds to `vnl`.
        *   *Equation/Purpose*: Calculates $\langle\phi_{\mu} | V_{nl} | \phi_{\nu}\rangle$ (3-center part).
    *   *Dependency*: Geometry, Fdata tables.

### Operations Performed EVERY SCF Step

These steps depend on the current self-consistent charge distribution (`Qin`) and are recalculated in each iteration of the SCF loop to achieve self-consistency.

1.  **Self Interaction Setup:**
    *   Sets up `neigh_self` and `neighPP_self`: Indices for the "self" neighbor (atom interacting with itself in the central cell).
    *   *Variables Populated*: `neigh_self`, `neighPP_self` (likely in `neighbor_map`).
    *   *Dependency*: Geometry (indices are fixed after Kscf=1 if `ifixneigh=1`).

2.  **Assemble 1-Center XC (if `itheory_xc` is 2 - McWEDA):**
    *   Calls `assemble_olsxc_1c`: Assembles the one-center (on-site) exchange-correlation potential matrix elements using the OLSXC method.
    *   *Variables Populated*: `vxc_1c`.
    *   *Equation/Purpose*: Calculates $\langle\phi_{i}^m | V_{xc}[\rho] | \phi_{i}^n \rangle$ (1-center part, OLSXC) based on the current density $\rho$.
    *   *Dependency*: Self-consistent charges (`Qin`).

3.  **Assemble SNXC or OLSXC (if `itheory_xc` is 1 or 2):**
    *   Assuming default `itheory_xc = 2` (McWEDA):
        *   Calls `average_ca_rho` (if `itheory` is 1 - DOGS): Calculates quantities related to the average density for charge-dependent terms in DOGS/OLSXC.
            *   *Variables Populated*: `rho_on`, `rhoi_on`, `rho_off`, `rhoij_off`, `arho_on`, `arhoi_on`, `arho_off`, `arhoij_off`.
            *   *Purpose*: Prepares density-related terms for OLSXC/DOGS assembly.
            *   *Dependency*: Self-consistent charges (`Qin`).
        *   Calls `assemble_olsxc_on`: Assembles the on-site (one-center) OLSXC exchange-correlation potential matrix elements and calculates the correction energy.
            *   *Variables Populated*: Adds to `vxc_1c`, populates `uxcdcc_ols`.
            *   *Equation/Purpose*: Calculates on-site OLSXC potential matrix elements and the corresponding double-counting correction energy based on the current self-consistent density.
            *   *Dependency*: Self-consistent charges (`Qin`).
        *   Calls `assemble_olsxc_off`: Assembles the off-site (two-center and three-center) OLSXC exchange-correlation potential matrix elements.
            *   *Variables Populated*: `vxc`.
            *   *Equation/Purpose*: Calculates off-site OLSXC potential matrix elements based on the current self-consistent density.
            *   *Dependency*: Self-consistent charges (`Qin`).

4.  **Assemble 2-Center DOGS (if `itheory` is 1):**
    *   Calls `assemble_ca_2c`: Assembles two-center charge-dependent potential corrections specific to DOGS theory.
    *   *Variables Populated*: `vca`, `vxc_ca`.
    *   *Equation/Purpose*: Calculates $\langle\phi_{\mu} | V_{ca}[\Delta Q] | \phi_{\nu}\rangle$ and $\langle\phi_{\mu} | V_{xc\_ca}[\Delta Q] | \phi_{\nu}\rangle$ (2-center parts) based on the current charge deviations $\Delta Q$ (derived from `Qin`).
    *   *Dependency*: Self-consistent charges (`Qin`).

5.  **Assemble 3-Center DOGS (if `itheory` is 1):**
    *   Calls `assemble_ca_3c`: Assembles three-center charge-dependent potential corrections specific to DOGS theory.
    *   *Variables Populated*: Adds to `vca`, `vxc_ca`.
    *   *Equation/Purpose*: Calculates $\langle\phi_{\mu} | V_{ca}[\Delta Q] | \phi_{\nu}\rangle$ and $\langle\phi_{\mu} | V_{xc\_ca}[\Delta Q] | \phi_{\nu}\rangle$ (3-center parts) based on the current charge deviations $\Delta Q$.
    *   *Dependency*: Self-consistent charges (`Qin`).

6.  **Assemble Long-Range Electrostatics (if `itheory` is 1 - DOGS, or if `idipole` is active):**
    *   Calls `assemble_lr`: Assembles long-range electrostatic contributions. This uses the `ewald` matrix (calculated at Kscf=1 based on initial charges) and applies it with the *current* self-consistent charges `Qin` to generate orbital-resolved `ewaldlr` and `ewaldsr` contributions to the Hamiltonian. If `idipole` is active, dipole corrections are also added here.
    *   *Variables Populated*: `ewaldlr`, `ewaldsr`.
    *   *Dependency*: Geometry, Self-consistent charges (`Qin`).

7.  **Build Final Matrices:**
    *   Calls `buildh`: Combines all the assembled matrix element pieces into the final Hamiltonian and Overlap matrices.
    *   *Variables Populated*: `h_mat`, `s_mat`.
    *   *Equation/Purpose*: Constructs $H = T + V_{na} + V_{nl} + V_{xc} + V_{xc\_1c} + V_{ca} + V_{xc\_ca} + V_{Ewald}$ and $S$ (Overlap matrix) by summing the contributions stored in the various `v*` and `t_mat`, `s_mat` arrays.
    *   *Dependency*: All previously assembled matrix element components.

**Summary of Dependencies:**

*   **Geometry-dependent (Kscf=1 only, non-SCF part):** Neighbor lists, base 2-center integrals (contributing to T, Vna, Vnl), base 3-center integrals (contributing to T, Vna, Vnl), initial Ewald sum matrix (`ewald`).
*   **Self-consistent charge (`Qin`)-dependent (Every SCF step):** 1-center XC (OLSXC), off-site multi-center XC (OLSXC), 2-center charge-dependent potentials (DOGS), 3-center charge-dependent potentials (DOGS), application of long-range electrostatics (using current `Qin` with pre-calculated `ewald` matrix).

The *base* 3-center integrals (parts of T, Vna, Vnl from `assemble_3c` and `assemble_3c_PP`) are indeed calculated at `Kscf=1` and are geometry-dependent, not SCF-charge dependent. The *self-consistent* 3-center character comes primarily from the OLSXC terms (`assemble_olsxc_off`) and the 3-center DOGS terms (`assemble_ca_3c`), which are recalculated in every SCF step using the current charges.

The Hartree potential is implicitly handled: the Ewald sum (`ewald`, `ewaldlr`, `ewaldsr`) accounts for long-range classical electrostatics between the self-consistent atomic charges (`Qin`). Short-range electron-electron interactions are effectively included within the density-dependent XC potentials (McWEDA/OLSXC).
