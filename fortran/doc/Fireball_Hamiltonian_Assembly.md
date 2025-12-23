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

### Subroutines Called in `assemble_mcweda.f90`

- **`assemble_2c`**: Assembles two-center kinetic energy and nuclear attraction potential matrix elements from pre-calculated integral tables.
- **`assemble_ca_2c`**: Assembles two-center charge-dependent potential corrections specific to DOGS theory based on current self-consistent charges.
- **`assemble_3c`**: Assembles three-center contributions to the Hamiltonian, including kinetic energy, nuclear attraction, and non-local pseudopotential terms.
- **`assemble_3c_PP`**: Assembles three-center pseudopotential matrix elements by contracting projector overlaps with coefficients.
- **`assemble_olsxc_off`**: Assembles off-site (two-center and three-center) exchange-correlation potential matrix elements using the OLSXC approximation.
- **`assemble_ca_3c`**: Assembles three-center charge-dependent potential corrections specific to DOGS theory based on current self-consistent charges.
- **`assemble_lr`**: Assembles long-range electrostatic contributions to the Hamiltonian from the Ewald sum using current atomic charges.


### Computational Complexity (O(n) Estimates)

Based on for-loop analysis in each subroutine, with scaling in terms of:  
- **N_a**: number of atoms.  
- **N**: total molecular orbitals/basis functions (sum of `num_orb(iatom)`).  
- **Assumptions**: Dense systems, `neighn(iatom)` ~ O(N_a), `neigh_comn(ialp)` ~ O(N_a^2), `num_orb`, `nssh` constants; `doscentros`/`trescentros` O(1).

- **`assemble_2c`**: O(N_a^2) or O(N^2) – loops over atoms and neighbors, with O(num_orb^2) work per pair.
- **`assemble_ca_2c`**: O(N_a^2) or O(N^2) – similar to `assemble_2c`, plus shell loops for charge corrections.
- **`assemble_3c`**: O(N_a^3) – loops over atoms and common neighbor pairs, with O(num_orb^2) per triplet.
- **`assemble_3c_PP`**: O(N_a^3) – similar to `assemble_3c`, with contraction over PP orbitals.
- **`assemble_olsxc_off`**: O(N_a^2) or O(N^2) – loops over atoms/neighbors, with shell/orbital matrix operations per pair.
- **`assemble_ca_3c`**: O(N_a^3) – like `assemble_3c`, with shell loops for charge-dependent terms.
- **`assemble_lr`**: O(N_a^2) or O(N^2) – applies precomputed Ewald matrix over orbital pairs per neighbor pair.


## Assmelbe McWEDA pseudocodes

Here is the mathematical documentation for the Fireball Hamiltonian assembly subroutines.

### Nomenclature

*   **Indices:**
    *   $I, J, K$: Atomic indices (Nuclear centers).
    *   $\mu, \nu$: Atomic orbital indices (basis functions).
    *   $n$: Projector index (for pseudopotentials).
*   **Coordinates:**
    *   $\mathbf{R}_I$: Position vector of atom $I$.
    *   $\mathbf{r}$: Electronic coordinate.
    *   $\mathbf{r}_{ij} = \mathbf{R}_J - \mathbf{R}_I$: Bond vector.
*   **Functions:**
    *   $\chi_{I\mu}(\mathbf{r}) \equiv \chi_\mu(\mathbf{r} - \mathbf{R}_I)$: Basis function (Fireball numerical orbital) centered at $I$.
    *   $\rho(\mathbf{r})$: Electronic density.
    *   $\bar{\rho}$: Average/Spherically symmetrized density.
*   **Potentials:**
    *   $V_{NA}$: Neutral Atom potential (short-range).
    *   $V_{CA}$: Charged Atom potential (long-range deviation from neutral).
    *   $V_{XC}$: Exchange-Correlation potential.
    *   $\hat{V}_{NL}$: Non-Local Pseudopotential operator.
*   **Matrices:**
    *   $S_{\mu\nu} = \langle \chi_\mu | \chi_\nu \rangle$: Overlap matrix.
    *   $H_{\mu\nu} = \langle \chi_\mu | \hat{H} | \chi_\nu \rangle$: Hamiltonian matrix.

---

### 1. `assemble_3c`
**Purpose:** Assembles the **Neutral Atom (NA)** 3-center interactions.
In Fireball, the potential is a sum of atomic potentials: $V(\mathbf{r}) = \sum_K V_{NA}^K(\mathbf{r} - \mathbf{R}_K)$.
This routine calculates terms where the potential center $K$ is distinct from the basis centers $I$ and $J$ (or at least treated as a distinct "third body" geometric entity).

**Mathematical Formulation:**
The goal is to compute:
$$ H_{I\mu, J\nu}^{3c} = \sum_K \langle \chi_{I\mu} | V_{NA}^K | \chi_{J\nu} \rangle $$

**Pseudocode / Math-Logic:**

*   **Loop** over potential centers $K$ (Atoms):
    *   $\mathbf{R}_{NA} = \mathbf{R}_K$
    *   **Loop** over neighbors $I$ of $K$:
        *   **Loop** over neighbors $J$ of $K$:
            *   *(Note: This creates a triplet $I, J, K$. $I$ and $J$ are the basis centers, $K$ provides the potential)*
            
            *   **Geometry Definition:**
                *   $\mathbf{y} = \mathbf{R}_J - \mathbf{R}_I$ (Bond vector between basis centers)
                *   $\mathbf{x} = \mathbf{R}_K - \frac{\mathbf{R}_I + \mathbf{R}_J}{2}$ (Vector from bond center to potential)
                *   $\cos\theta = (\hat{\mathbf{y}} \cdot \hat{\mathbf{x}})$
            
            *   **Interpolation (`trescentros`):**
                *   Retrieve precalculated 3-center integral table $T(x, y, \theta)$ for species tuple $(Z_I, Z_J, Z_K)$.
                *   $v_{\mu\nu}^{local} = \text{SplineInterp}(T, x, y, \theta)$
            
            *   **Accumulation:**
                *   $H_{I\mu, J\nu} \leftarrow H_{I\mu, J\nu} + v_{\mu\nu}^{local}$
                *   $H_{J\nu, I\mu} \leftarrow H_{J\nu, I\mu} + v_{\mu\nu}^{local}$ (Symmetrization)

---

### 2. `assemble_ca_3c`
**Purpose:** Assembles the **Charged Atom (CA)** 3-center interactions.
This accounts for the fact that atoms in the solid are not neutral. The potential is corrected by a term proportional to the charge transfer $\delta Q_K$. It uses an Ewald-like summation where short-range terms are exact (interpolated) and long-range terms are approximated as multipoles, blended by a switching function.

**Mathematical Formulation:**
$$ H_{I\mu, J\nu}^{CA} = \sum_K \delta Q_K \left[ \underbrace{S_{smooth} \cdot \langle \chi_{I\mu} | V_{short}^K | \chi_{J\nu} \rangle}_{\text{Exact 3C}} + \underbrace{(1 - S_{smooth}) \cdot \langle \chi_{I\mu} | V_{long}^K | \chi_{J\nu} \rangle}_{\text{Approximate (Multipole)}} \right] $$

**Pseudocode / Math-Logic:**

*   **Loop** over potential centers $K$:
    *   Calculate charge deviation: $\delta Q_K = \sum_{\text{shell } s} (Q_K(s) - Q_K^{neutral}(s))$
    *   **Loop** over neighbors $I, J$ of $K$:
        *   **Geometry:**
            *   Calculate distances $R_{IK} = |\mathbf{R}_K - \mathbf{R}_I|$, $R_{JK} = |\mathbf{R}_K - \mathbf{R}_J|$
            *   Calculate 3C coords $x, y, \theta$ (same as `assemble_3c`).
        
        *   **Smoothing Factor (`smoother`):**
            *   $f_I = \text{Smooth}(R_{IK}, R_{cut})$, $f_J = \text{Smooth}(R_{JK}, R_{cut})$
            *   $S_{smooth} = f_I \cdot f_J$
        
        *   **Long-Range Term (Multipole):**
            *   Approximate potential using overlaps $S_{\mu\nu}$ and dipoles $D_{\mu\nu}$.
            *   $V_{mono} = \frac{1}{2} S_{\mu\nu}^{IJ} (\frac{1}{R_{IK}} + \frac{1}{R_{JK}}) + \text{DipoleCorrections}$
            *   $E_{waldSR} \leftarrow E_{waldSR} + (1 - S_{smooth}) \cdot V_{mono} \cdot \delta Q_K$
        
        *   **Short-Range Term (Exact):**
            *   $v_{exact} = \text{SplineInterp}(T_{CA}, x, y, \theta)$
        
        *   **Accumulation:**
            *   $V_{CA}^{matrix} \leftarrow V_{CA}^{matrix} + \delta Q_K \cdot [ S_{smooth} \cdot v_{exact} + (1-S_{smooth}) \cdot V_{mono} ]$

---

### 3. `assemble_3c_PP`
**Purpose:** Assembles the **Non-Local Pseudopotential** interactions.
Unlike NA/CA potentials which are multiplicative functions $V(r)$, the Pseudopotential is a projection operator sum.

**Mathematical Formulation:**
$$ \hat{V}_{NL} = \sum_K \sum_{n} \epsilon_{K,n} | \alpha_{K,n} \rangle \langle \alpha_{K,n} | $$
Matrix element:
$$ H_{I\mu, J\nu}^{PP} = \sum_K \sum_{n} \epsilon_{K,n} \underbrace{\langle \chi_{I\mu} | \alpha_{K,n} \rangle}_{S^{VNL}_{I\mu, K n}} \underbrace{\langle \alpha_{K,n} | \chi_{J\nu} \rangle}_{S^{VNL}_{J\nu, K n}} $$

**Pseudocode / Math-Logic:**

*   **Loop** over projection centers $K$:
    *   Load PP coefficients $\epsilon_{K}$ (`cl_value`).
    *   **Loop** over neighbor atoms $I$ and $J$ (where $\chi_I$ and $\chi_J$ overlap with $\alpha_K$):
        *   **Contraction:**
            *   Initialize $V_{local} = 0$
            *   **Loop** over projector orbitals $n$ on $K$:
                *   Fetch precomputed Overlap: $A = \langle \chi_{I\mu} | \alpha_{K,n} \rangle$ (`sVNL`)
                *   Fetch precomputed Overlap: $B = \langle \chi_{J\nu} | \alpha_{K,n} \rangle$
                *   $V_{local} \leftarrow V_{local} + \epsilon_{K,n} \cdot A \cdot B$
        *   **Accumulate:**
            *   $V_{NL}(\mu, \nu) \leftarrow V_{NL}(\mu, \nu) + V_{local}$

---

### 4. `build_olsxc_off`
**Purpose:** Computes the Exchange-Correlation (XC) matrix elements for a specific pair $(I, J)$ using the **Ortega-Lewis-Sankey (OLS)** approximation. This is a first-order expansion around a superposition of spherical atomic densities.

**Mathematical Formulation:**
*   Reference Density: $\bar{\rho} = \rho_I^{atom} + \rho_J^{atom}$ (at the bond center).
*   Potential Expansion: $V_{XC}(\rho) \approx V_{XC}(\bar{\rho}) + \frac{dV_{XC}}{d\rho}(\bar{\rho}) \cdot (\rho - \bar{\rho})$.
*   Matrix Element:
    $$ \langle \chi_\mu | V_{XC} | \chi_\nu \rangle \approx V_{XC}(\bar{\rho}) S_{\mu\nu} + V'_{XC}(\bar{\rho}) \left[ \langle \chi_\mu | \rho | \chi_\nu \rangle - \bar{\rho} S_{\mu\nu} \right] $$

**Pseudocode / Math-Logic:**

*   **Inputs:**
    *   $S_{\mu\nu}$: Overlap matrix element.
    *   $\rho_{\mu\nu}^{matrix}$: The actual density matrix element $\langle \chi_\mu | \hat{\rho} | \chi_\nu \rangle$.
    *   $\bar{\rho}_{avg}$: Spherically averaged atomic density values.
*   **Compute Potentials (`cepal`):**
    *   $\mu_{xc} = V_{XC}(\bar{\rho}_{avg})$
    *   $\mu'_{xc} = \frac{dV_{XC}}{d\rho}(\bar{\rho}_{avg})$
*   **Construct Matrix Element:**
    *   Term 1 (0th order): $val = \mu_{xc} \cdot S_{\mu\nu}$
    *   Term 2 (1st order correction): $val \leftarrow val + \mu'_{xc} \cdot (\rho_{\mu\nu}^{matrix} - \bar{\rho}_{avg} \cdot S_{\mu\nu})$
*   **Correction for Non-Linearity (McWeda/OLS terms):**
    *   Subtract double counting terms involving single-center densities (equations in code are essentially $V_{XC}(\rho_{tot}) - V_{XC}(\rho_I) - V_{XC}(\rho_J)$ logic).

---

### 5. `assemble_olsxc_off`
**Purpose:** The driver loop that calculates the "Off-Diagonal" (2-center) XC terms. It gathers the geometric overlaps and calls `build_olsxc_off`.

**Pseudocode / Math-Logic:**

*   **Loop** over Atom $I$:
    *   **Loop** over Neighbor $J$:
        *   **Geometry:** Calculate $R_{IJ}$.
        *   **Get Atomic Density Overlaps (`doscentros`):**
            *   Compute $\langle \chi_I | \rho_{atom} | \chi_J \rangle$.
            *   This populates temporary arrays `denmx`, `den1x` representing the reference density overlaps.
        *   **Get Basis Overlap:**
            *   Retrieve $S_{\mu\nu}$ (`s_mat`).
        *   **Calculate Matrix Element:**
            *   Call `build_olsxc_off(denmx, den1x, S)` to get $V_{XC}^{\mu\nu}$.
        *   **Accumulate:**
            *   $V_{XC}^{global}(\mu, \nu) \leftarrow V_{XC}^{global}(\mu, \nu) + V_{XC}^{\mu\nu}$.

---

### 6. `assemble_mcweda`
**Purpose:** The **Master Driver**. It orchestrates the entire Hamiltonian construction by summing all physical contributions.
$$ H_{total} = T + V_{NA}^{2c} + V_{NA}^{3c} + V_{NL}^{2c} + V_{NL}^{3c} + V_{CA}^{2c} + V_{CA}^{3c} + V_{XC} $$

**Pseudocode / Math-Logic:**

1.  **Initialization:**
    *   Update Neighbor Lists (if atoms moved).
    *   Calculate Ewald Energy (Electrostatics).

2.  **1-Center Terms (On-site):**
    *   $H \leftarrow \text{AtomicLevels}$.
    *   Add On-site XC contributions (`assemble_olsxc_1c`).

3.  **2-Center Terms (Pairwise):**
    *   Call `assemble_sVNL()`: Precompute $\langle \chi | \alpha \rangle$.
    *   Call `assemble_2c()`: Computes $S, T, V_{NA}^{2c}$.
    *   Call `assemble_2c_PP()`: Computes $V_{NL}^{2c}$ (where $I,J$ are neighbors of PP center $K$ but $K=I$ or $K=J$).
    *   Call `assemble_snxc_off` / `assemble_olsxc_off`: Computes $V_{XC}^{2c}$.

4.  **3-Center Terms (Triplets):**
    *   Call `assemble_3c()`: Computes $V_{NA}^{3c}$ (Geometric triplets).
    *   Call `assemble_3c_PP()`: Computes $V_{NL}^{3c}$ (Projection triplets).
    *   Call `assemble_ca_3c()`: Computes $V_{CA}^{3c}$ (Charged atom corrections).

5.  **Long Range:**
    *   Call `assemble_lr()`: Ewald summation for widely separated atoms.

6.  **Final Summation:**
    *   $H_{final} = \sum (\text{All Matrix Arrays})$.
    *   Symmetrize $H$ if necessary.


Here is the mathematical documentation and pseudo-code analysis for the provided Fortran files.

### 1. `trescentros.f` (The 3-Center Integrator)
**Purpose:** Calculates the matrix element $\langle \chi_{I\mu} | V_K | \chi_{J\nu} \rangle$ where the potential center $K$ is geometrically distinct from the bond center of $I-J$.
It maps the geometry $(x, y, \cos\theta)$ to the integral value using a 2D spline interpolation combined with a Legendre polynomial expansion for the angular dependence.

**Mathematical Formulation:**
*   **Coordinates:**
    *   $y = |\mathbf{R}_J - \mathbf{R}_I|$ (Bond length)
    *   $\mathbf{R}_{mid} = (\mathbf{R}_I + \mathbf{R}_J) / 2$
    *   $x = |\mathbf{R}_K - \mathbf{R}_{mid}|$ (Distance from bond center to potential)
    *   $\cos\theta = \widehat{(\mathbf{R}_J - \mathbf{R}_I)} \cdot \widehat{(\mathbf{R}_K - \mathbf{R}_{mid})}$
*   **Expansion:**
    $$ M_{\mu\nu}^{local}(x, y, \theta) = \sum_{l=0}^{4} C_{l}^{\mu\nu}(x, y) P_l(\cos\theta) $$
    Where $P_l$ are Legendre polynomials.
*   **Rotation:**
    The result $M^{local}$ is in the molecular frame (z-axis along bond). It must be rotated to the crystal frame using the Slater-Koster rotation matrix $\mathcal{R}$ derived from `epsilon`.
    $$ M_{\mu\nu}^{crystal} = \mathcal{R} \cdot M_{\mu\nu}^{local} \cdot \mathcal{R}^T $$

**Math-Pseudocode:**
*   **Input:** Geometry ($x, y, \cos\theta$), Species, Shells.
*   **Interpolate coefficients:**
    *   For $l = 0 \dots 4$:
        *   $Q_l = \text{Spline2D}(\text{Table}_{species, isorp, l}, x, y)$
*   **Legendre Summation:**
    *   $P_0 = 1$
    *   $P_1 = \cos\theta$
    *   $P_2 = \frac{1}{2}(3\cos^2\theta - 1)$
    *   ... (up to $P_4$)
    *   $Val_{\text{molecular}} = \sum_l Q_l \cdot P_l$
    *   *Correction:* If $m$-value requires sine dependence: $Val \leftarrow Val \cdot \sin\theta$
*   **Basis Recovery:**
    *   Map the 1D list of values to the 2D shell pairs $(\mu, \nu)$.
*   **Rotation:**
    *   $Val_{\text{crystal}} = \text{Rotate}(Val_{\text{molecular}}, \epsilon)$

---

### 2. `doscentros.f` (The 2-Center Integrator)
**Purpose:** Calculates 2-center integrals $\langle \chi_{I\mu} | \hat{O} | \chi_{J\nu} \rangle$. This includes Overlaps ($S$), Kinetic Energy ($T$), and 2-center Potential terms ($V_{NA}^{2c}$).

**Mathematical Formulation:**
$$ M_{\mu\nu}(\mathbf{R}_{IJ}) = \hat{D}(\mathbf{R}_{IJ}) \left[ \text{Spline1D}(|\mathbf{R}_{IJ}|) \right] $$
Where $\hat{D}$ is the Slater-Koster rotation operator.

**Math-Pseudocode:**
*   **Input:** Interaction Type (S, T, V, etc.), Atoms $I, J$.
*   **Interpolate:**
    *   $R = |\mathbf{R}_J - \mathbf{R}_I|$
    *   $v_{slater}(R) = \text{Spline1D}(\text{Table}_{type}, R)$
    *   Returns array of components: $\sigma, \pi, \delta, \dots$
*   **Recover & Rotate:**
    *   Construct diagonal matrix $M_{mol}$ from $\sigma, \pi, \delta$.
    *   $M_{cryst} = \text{Rotate}_{SK}(M_{mol}, \epsilon)$
*   **Forces (if `iforce=1`):**
    *   Compute $\frac{d M}{d R}$ via spline derivative.
    *   Apply chain rule for direction vectors using `eta` vector from `epsilon`.

---

### 3. `epsilon.f` (Coordinate System)
**Purpose:** Constructs the Local-to-Global rotation matrix (Metric Tensor).

**Mathematical Formulation:**
Constructs an orthonormal basis $\{ \hat{\mathbf{x}}', \hat{\mathbf{y}}', \hat{\mathbf{z}}' \}$ where $\hat{\mathbf{z}}'$ points along the bond.
*   $\hat{\mathbf{z}}' = \frac{\mathbf{r}_{ij}}{|\mathbf{r}_{ij}|}$
*   $\hat{\mathbf{y}}' = \frac{\hat{\mathbf{z}}' \times \mathbf{r}_{ref}}{|\dots|}$ (with handling for collinear singularities)
*   $\hat{\mathbf{x}}' = \hat{\mathbf{y}}' \times \hat{\mathbf{z}}'$

---

### 4. `assemble_sVNL.f90` (VNL Overlap Assembler)
**Purpose:** Pre-calculates the overlaps between atomic basis functions $\chi$ and the Non-Local Pseudopotential projectors $\alpha_{Kn}$. These are needed for the $V_{NL}$ contraction.

**Mathematical Formulation:**
$$ S^{VNL}_{I\mu, K n} = \langle \chi_{I\mu} | \alpha_{K n} \rangle $$
This is structurally identical to an orbital overlap $S_{\mu\nu}$ but uses specific VNL tables.

**Math-Pseudocode:**
*   **Loop** $I$ (Atom):
    *   **Loop** $K$ (Projector Atom nearby):
        *   **Calculation:**
            *   Call `doscentros(interaction=5)` (Type 5 = Non-Local)
            *   This returns the vector of overlaps between all orbitals $\mu$ on $I$ and all projectors $n$ on $K$.
        *   **Storage:**
            *   Store in `sVNL(mu, n, neighbor_index, atom_index)`.
            *   Note: This table is not symmetric in indices (Basis vs Projector).

---

### 5. `average_rho.f90` (Density Averaging)
**Purpose:** Computes the spherically averaged density $\bar{\rho}$ required for the OLS-XC functional.
It computes: $\bar{\rho} = \frac{\int \rho(\mathbf{r}) \chi_i \chi_j d\mathbf{r}}{\int \chi_i \chi_j d\mathbf{r}}$.

**Mathematical Formulation:**
$$ \bar{\rho}_{IJ}^{\mu\nu} = \frac{\langle \chi_{I\mu} | \hat{\rho} | \chi_{J\nu} \rangle}{S_{IJ}^{\mu\nu}} $$
The numerator includes contributions from on-site density, neighbor densities (2-center), and other neighbors (3-center).

**Math-Pseudocode:**
*   **Initialize:** $\bar{\rho} = 0$.
*   **1. On-Site Part ($I=J$):**
    *   $\text{Num} = \langle \chi_I | \rho_I + \sum \rho_{neigh} | \chi_I \rangle$
    *   $\text{Denom} = S_{II} = 1$
    *   $\bar{\rho}_{on} \leftarrow \text{Num} / \text{Denom}$
*   **2. Off-Site Part ($I \ne J$):**
    *   **Loop** Neighbors $(I, J)$:
        *   **3-Center Contribution:**
            *   Loop $K$ (Common Neighbor):
            *   $val_{3c} = \text{trescentros}(\text{interaction}=3)$ (Density triplet)
            *   $\langle \chi_I | \rho_K | \chi_J \rangle \leftarrow val_{3c} \cdot Q_K$
        *   **2-Center Contribution:**
            *   Call `doscentros(interaction=15)`: $\langle \chi_I | \rho_I | \chi_J \rangle$
            *   Call `doscentros(interaction=16)`: $\langle \chi_I | \rho_J | \chi_J \rangle$
        *   **Normalization:**
            *   $S_{IJ} = \text{doscentrosS}(I, J)$ (Overlap)
            *   Avoid singularity: If $|S_{IJ}| < \epsilon$, set $S_{IJ} = \epsilon$.
            *   $\bar{\rho}_{off} = \frac{\text{Sum of 2C and 3C density terms}}{S_{IJ}}$

---

### 6. `assemble_ca_2c.f90` (Charged Atom 2-Center)
**Purpose:** Assembles 2-center terms derived from charge transfer $\delta Q$.
It includes a Short-Range (exact) term and a Long-Range (Multipole) term, blended by a smoother.

**Mathematical Formulation:**
$$ H_{I\mu, J\nu}^{CA-2c} = \sum_{K \in \{I,J\}} \delta Q_K \left[ f(r) \cdot V_{exact}^{2c} + (1-f(r)) \cdot V_{multipole} \right] $$

**Math-Pseudocode:**
*   **Loop** Pairs $(I, J)$:
    *   $\delta Q_I = Q_I - Q_I^{neutral}$
    *   $\delta Q_J = Q_J - Q_J^{neutral}$
    *   **Dipoles:**
        *   $D_{\mu\nu} = \text{doscentros(interaction=9)}$ (Gets dipole integrals)
    *   **Ewald/Multipole Term (Long Range):**
        *   $V_{mono} \approx \frac{S_{IJ}}{R_{IJ}}$ (plus dipole corrections $D/R^2$)
        *   $W_{long} = V_{mono} \cdot (\delta Q_I + \delta Q_J)$
        *   Calculate Smoother $f(R_{IJ}) = \text{smoother}(R_{IJ}, R_{cut})$.
        *   $E_{wald} \leftarrow (1 - f) \cdot W_{long}$
    *   **VNA Term (Short Range):**
        *   *Ontop Left:* $\langle \chi_I | V_I | \chi_J \rangle$ via `doscentros(interaction=2)`.
        *   *Ontop Right:* $\langle \chi_I | V_J | \chi_J \rangle$ via `doscentros(interaction=3)`.
        *   *Accumulate:* $V_{CA} \leftarrow f \cdot (\delta Q_I V_{left} + \delta Q_J V_{right})$.

---

### 7. `smoother.f`
**Purpose:** A polynomial blending function.

**Mathematical Formulation:**
$$ f(r) = \begin{cases} 1 & r < R_{start} \\ 0 & r > R_{end} \\ P(\frac{r - R_{start}}{R_{end} - R_{start}}) & \text{otherwise} \end{cases} $$
Where $P(x)$ is typically $(1-x^n)^m$ (Old Method) or a polynomial $1 + ax^2 + bx^3 + cx^4$ (New Method).

**Math-Pseudocode:**
*   $x = (r - R_{start}) / (R_{end} - R_{start})$
*   `stn` = Evaluate Polynomial($x$)
*   `dstn` = Evaluate Derivative($x$)


---

## How to implentn on GPU

Based on the analysis of the mathematical structure and the computational constraints, here is the proposed OpenCL implementation strategy.

### Core Philosophy
1.  **Gather, Don't Scatter:** Invert the Fortran loops. Instead of iterating over atoms $K$ and adding to neighbors $(I,J)$, we iterate over output pairs $(I,J)$ and gather contributions from neighbors $K$. This eliminates global atomic writes.
2.  **Interaction Lists:** We will not use $N_{atom} \times N_{atom}$ matrices. We will use a **Sorted Interaction List** of pairs $(I, J)$ representing the non-zero block-sparse structure of the Hamiltonian.
3.  **Kernel Specialization:** Split kernels by *mathematical operation type* (Spline Interpolation vs. Matrix Contraction vs. Grid Sampling) rather than by physical term ($T$ vs $V$).

---

### Phase 0: Data Structure Preparation (CPU/GPU Helper)

Before physics kernels run, we must organize the geometry.

**Structure:** `InteractionList`
*   A flat array of structs: `{atom_i, atom_j, species_pair_id, output_index}`.
*   **Sorted** by `species_pair_id` (Crucial for L1 cache hits on spline tables).

**Structure:** `CommonNeighborList (CNL)`
*   For every pair $(I, J)$ in the Interaction List, we need a list of atoms $K$ that are within cutoff of *both*.
*   **Format:** A Compact Row Storage (CRS) style list.
    *   `CNL_Offsets`: Index where the neighbor list for pair $(I,J)$ starts.
    *   `CNL_Data`: The indices of atoms $K$.

---

### Phase 1: The "Universal 2-Center" Kernel
**Math:** 1D Spline Interpolation + Slater-Koster Rotation.
**Fortran Equivalents:** `assemble_2c`, `assemble_sVNL`, `doscentros`.

This kernel computes all "simple" pairwise terms in one go to minimize memory reads of geometry $\mathbf{R}_{IJ}$.

**Inputs:** `ratoms`, `InteractionList`, Spline Tables (S, T, Vna, Vnl).
**Outputs:**
1.  **`Buffer_S`**: Overlap Matrix (Needed for Density).
2.  **`Buffer_H_2c`**: Sum of $T + V_{NA}^{2c} + V_{XC}^{0}$.
3.  **`Buffer_sVNL`**: The Basis-Projector overlaps $\langle \chi_I | \Psi_K \rangle$. *Note: This requires a slightly different interaction list (Pairs $I, K$), but can likely be fused or run as a second pass of the same kernel template.*

**Optimization:**
*   Fuse $S, T, V$ interpolations.
*   Compute distance $R_{IJ}$ and rotation vectors $\hat{x}, \hat{y}, \hat{z}$ **once** per thread and apply to all terms.

---

### Phase 2: Density Assembly (The Dependency)
**Math:** Weighted Summation (Gather).
**Fortran Equivalents:** `average_rho`, `build_olsxc_off`.

$V_{XC}$ depends on the average density $\bar{\rho}$. We cannot compute the final XC potential until we have $\bar{\rho}$.

**Kernel: `compute_avg_rho`**
*   **Work Item:** Pair $(I, J)$.
*   **Logic:**
    1.  Read `Buffer_S` (Overlap $S_{IJ}$).
    2.  Gather On-site densities $\rho_I, \rho_J$.
    3.  Gather 3-Center densities (using `CNL`): Loop $K$, interpolate Density-Triplet tables.
    4.  Compute $\bar{\rho}_{IJ} = \frac{\langle \rho \rangle}{S_{IJ}}$.
*   **Output:** `Buffer_Rho_Avg` (Scalar per orbital pair).

---

### Phase 3: The "Heavy" 3-Center Gather Kernel
**Math:** 3D Interpolation (2D Grid + Legendre) + Geometric Invariants.
**Fortran Equivalents:** `assemble_3c`, `assemble_ca_3c`, `trescentros`.

This is the most expensive kernel. It handles Neutral Atom ($V_{NA}^{3c}$) and Charged Atom ($V_{CA}^{3c}$) terms.

**Kernel: `assemble_3c_universal`**
*   **Work Item:** Pair $(I, J)$.
*   **Inputs:** `InteractionList`, `CNL`, `Fdata_3C` (Textures/Buffers), `Charges_Q`.
*   **Logic:**
    1.  Initialize `acc_NA = 0`, `acc_CA = 0`.
    2.  **Loop over $K$** in `CNL(I, J)`:
        *   Fetch $\mathbf{R}_K$, Charges $Q_K$.
        *   Compute 3C Invariants: $x, y, \cos\theta$.
        *   **Interpolate:**
            *   $v_{NA} = \text{Sample}(Table_{NA}, x, y, \theta)$
            *   $v_{CA\_short} = \text{Sample}(Table_{CA}, x, y, \theta)$
        *   **Ewald/Smoother:**
            *   Compute `smoother(r)` function.
            *   Compute Multipole approximation $v_{long}$.
            *   Blend: $v_{CA} = \text{smooth} \cdot v_{CA\_short} + (1-\text{smooth}) \cdot v_{long}$.
        *   Accumulate: `acc_NA += v_NA`, `acc_CA += v_CA * delta_Q`.
    3.  Write result to `Buffer_H_3c`.

**Optimization:**
*   This kernel is register-heavy. Do not include XC or PP logic here.
*   The "Grid Sampling" math is identical for NA and CA; only the tables and charge pre-factors differ.

---

### Phase 4: Pseudopotential Contraction (VNL)
**Math:** Sparse Matrix Contraction (No Grids).
**Fortran Equivalents:** `assemble_3c_PP`.

This is mathematically distinct from Phase 3. It's linear algebra, not interpolation.

**Kernel: `contract_pp`**
*   **Work Item:** Pair $(I, J)$.
*   **Inputs:** `Buffer_sVNL` (Computed in Phase 1), `PP_Coeffs`.
*   **Logic:**
    1.  **Loop over $K$** (Projector Atoms nearby):
        *   *Note: This loop might need a slightly different neighbor list than CNL, specifically atoms $K$ where $S_{VNL}(I,K)$ is non-zero.*
    2.  Load Matrix Block $A = \text{sVNL}[I, K]$.
    3.  Load Matrix Block $B = \text{sVNL}[J, K]$.
    4.  Load Coeffs $C_K$.
    5.  Multiply: `acc += A * C * Transpose(B)`.
*   **Output:** Accumulate into `Buffer_H_PP`.

---

### Phase 5: XC & Final Assembly
**Math:** Non-linear Functional Evaluation + Summation.
**Fortran Equivalents:** `assemble_olsxc`, `buildh`.

**Kernel: `assemble_final`**
*   **Work Item:** Pair $(I, J)$.
*   **Inputs:** All intermediate buffers (`H_2c`, `H_3c`, `H_PP`, `Rho_Avg`), `Buffer_S`.
*   **Logic:**
    1.  Load $\bar{\rho}$ from `Buffer_Rho_Avg`.
    2.  **Compute OLS-XC:**
        *   Calculate $\mu_{xc}(\bar{\rho})$ and derivatives.
        *   Compute $V_{XC} = \mu_{xc} \cdot S_{IJ} + \dots$ (OLS formula).
    3.  **Summation:**
        *   $H_{total} = H_{2c} + H_{3c} + H_{PP} + V_{XC}$.
    4.  **Write:** Store final Hamiltonian block.

---

### Summary of Kernel Split

| Kernel Name | Mathematical Operation | Inputs | Dependencies |
| :--- | :--- | :--- | :--- |
| **`kernel_2c_gen`** | 1D Spline + Rotation | Geom | None |
| **`kernel_rho_gather`** | 3C Gather + 1D Spline | Geom, `Buffer_S` | `kernel_2c_gen` |
| **`kernel_3c_gather`** | 3D Interpolation + Ewald | Geom, Charges | Charges updated |
| **`kernel_pp_contract`** | Sparse MatMul ($A \cdot C \cdot B^T$) | `Buffer_sVNL` | `kernel_2c_gen` |
| **`kernel_finalize`** | Functional ($V(\rho)$) + Sum | All Buffers | All above |

### Why this strategy?
1.  **Decouples Memory Patterns:** Spline interpolation (random access), 3D Grid interpolation (heavy calculation), and Matrix Contraction (linear memory access) are separated into different kernels to maximize occupancy for each specific workload.
2.  **Solves Dependencies:** The XC density dependency is handled by an explicit intermediate buffer.
3.  **No Global Atomics:** Every kernel parallelizes over the *output* element $(I,J)$.
4.  **Load Balancing:** Sorting the Interaction List by species ensures that within a warp/workgroup, all threads are likely accessing the same tables and executing the same branches.