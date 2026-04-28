

**The Harris Functional**
The Harris functional is an approximate energy functional used for calculating ground-state electronic properties. It is based on the Kohn-Sham implementation of DFT and was developed to provide simplified but accurate computational schemes for complex electronic systems, particularly the interaction between weakly interacting fragments. A key feature is that **the energy functional is defined entirely in terms of an input charge density**, often taken as a sum of overlapped frozen fragment densities. The error in the total energy estimate is only second order in the errors in this input density.

Within the context of the Sankey-Niklewski method, the approximate total energy functional based on the Harris approach is given as :
$$E_{tot}^{approx} = E_{BS} + U_{SR} + \delta U_{xc}$$
where:
*   **EBS** is the **band-structure energy** , calculated as twice the sum of the occupied eigenvalues (assuming spin degeneracy) :
    $$E_{BS} = 2 \sum \epsilon_i$$
    The eigenvalues ε_i are obtained by solving a **non-self-consistent single-particle Schrödinger-like equation** with a Hamiltonian (h) constructed using the input density (n₀), typically a sum of neutral atomic densities :
    $$\hat{h}[n_0] \psi_i = \epsilon_i \psi_i$$
    The Hamiltonian h includes the kinetic energy, the ionic potential, the neutral-atom Hartree potential (VN₀), and the exchange-correlation potential (μxc) evaluated at the input density n₀ :
    $$\hat{h} = \hat{T} + V_{ions} + V_{N_0}[n_0] + \mu_{xc}(n_0)$$
*   **USR** is the **short-range repulsive potential** . It represents the repulsive interactions between atoms, including the screened nuclear repulsion and the repulsive contributions from the overlapping neutral-atom charge densities . It can be written as a sum of one-body terms and a two-body central potential between neutral atoms :
    $$U_{SR} = \sum_i U_1(R_i) + \frac{1}{2} \sum_{i,j} V_{SR}(|R_i-R_j|)$$
*   **δUxc** is the **exchange-correlation correction term** . It corrects for the fact that the exchange-correlation potential μxc(n₀) is used in the Hamiltonian h, while the total energy involves the exchange-correlation energy functional εxc(n₀). It is given by the integral over the input density n₀ :
    $$\delta U_{xc} = \int n_0(r)[\epsilon_{xc}(n_0) - \mu_{xc}(n_0)]dr$$


**The Sankey-Niklewski (Original Method)**

The original Sankey-Niklewski (SN) method is an **ab initio multicenter tight-binding model** developed for molecular-dynamics simulations and other applications in covalent systems . It is founded on **Density Functional Theory within the Local Density Approximation (LDA) and the pseudopotential scheme** . It is a first-principles method that **does not require any empirical fitting** . It utilizes the approximate Harris energy functional, taking the input density (n₀) as a **simple sum of neutral-atom spherical atomic densities** .

The total energy of the system is approximated by the Harris-like functional mentioned above . Atomic forces are determined by taking the derivative of this total energy with respect to atomic positions , and the method **includes Pulay corrections exactly** .

The core of the method involves solving a one-electron eigenvalue equation in a basis of atomic-like orbitals :
$$\sum_\nu (h)_{\mu\nu} a_{i\nu} = \epsilon_i \sum_\nu (S)_{\mu\nu} a_{i\nu}$$ (Generalized eigenvalue problem det|h - εS| = 0 )
where ε_i are the electronic eigenvalues and a_iν are the expansion coefficients of the electronic eigenstates ψ_i in terms of atomic-like orbitals ϕ_ν :
$$\psi_i = \sum_\nu a_{i\nu} \phi_\nu$$

The **Hamiltonian matrix elements $(h)_{\mu\nu}$** and **overlap matrix elements $(S)_{\mu\nu}$** are calculated between these atomic-like orbitals :
$$\langle\phi_\mu | \hat{h} | \phi_\nu\rangle$$
$$\langle\phi_\mu | \phi_\nu\rangle$$

A key innovation of the SN method is the use of **slightly excited pseudo-atomic-orbitals (PAO's) called "fireballs"** as the basis set. These orbitals are defined by the boundary condition that they **vanish outside and at a predetermined radius (rᶜ)**. This gives the Hamiltonian and overlap matrix elements a **short range**, meaning they are exactly zero beyond a certain distance between the atoms.

These matrix elements are **calculated in real space** using a **multicenter approach** (including one-, two-, and three-center integrals). Four-center integrals are avoided. The values of these integrals for different geometries are **pre-calculated and stored in one- and two-dimensional data tables**, which are then used via interpolation during molecular dynamics simulations.

The original SN method calculated the exchange-correlation matrix elements using an approximation based on an **"average density" (n̄)** defined for each matrix element:
$$\bar{n} = \langle\phi_\mu|n|\phi_\nu\rangle / S_{\mu\nu}$$
where n is the electron density and S_μν is the overlap integral. The exchange-correlation matrix element was approximated based on the value of the exchange-correlation energy functional at this average density.

**The McWEDA Method**

The Multi-center Weighted Exchange-correlation Density Approximations (McWEDA) method is an **improved approach for calculating the exchange-correlation (XC) contributions** within first-principles tight-binding methods like FIREBALL. It was developed to address certain deficiencies in the original Sankey-Niklewski and the Horsfield approximations for XC terms, particularly in describing on-site terms and computational efficiency.

McWEDA treats **on-site and off-site XC matrix elements separately**. The core idea is to write each matrix element as a **dominant main contribution** (one-center for on-site, two-center for off-site) that is calculated exactly, **plus a correction term**. The correction term is then approximated using a **Generalized Sankey-Niklewski (GSN) approach**.

The GSN approach generalizes the original SN average density definition by using new **positive-definite weighting functions (w)** based on the radial part of the orbitals. The GSN average density ρ̄mn is defined as:
$$\bar{\rho}_{mn} = \langle w_m|\rho|w_n\rangle / \langle w_m|w_n\rangle$$
where ρ is the electron density and wm, wn are the weighting functions associated with orbitals ϕm, ϕn. This definition avoids issues encountered with the original SN average density, such as being undefined for zero overlap.

The GSN formula approximates an XC matrix element based on the XC potential and its derivative evaluated at this average density:
$$\langle\phi_m|V_{xc}[\rho]|\phi_n\rangle \approx V_{xc}[\bar{\rho}_{mn}] \langle\phi_m|\phi_n\rangle + V_{xc}'[\bar{\rho}_{mn}][\langle\phi_m|\rho r|\phi_n\rangle - \bar{\rho}_{mn}\langle\phi_m|\phi_n\rangle]$$

The explicit McWEDA equations for calculating the XC potential matrix elements (<ϕm|Vxc[ρ]|ϕn>) are as follows:

*   **On-site matrix elements (<ϕ_iᵐ|Vxc[ρ]|ϕ_iⁿ>):**
    $$\langle\phi_i^m|V_{xc}[\rho]|\phi_i^n\rangle \approx \langle\phi_i^m|V_{xc}[\rho_i]|\phi_i^n\rangle + V_{xc}[\bar{\rho}_{mn}]\langle\phi_i^m|\phi_i^n\rangle + V_{xc}'[\bar{\rho}_{mn}][\langle\phi_i^m|\rho r|\phi_i^n\rangle - \bar{\rho}_{mn}\langle\phi_i^m|\phi_i^n\rangle] - V_{xc}[\bar{\rho}_i]\langle\phi_i^m|\phi_i^n\rangle - V_{xc}'[\bar{\rho}_i][\langle\phi_i^m|\rho_i r|\phi_i^n\rangle - \bar{\rho}_i\langle\phi_i^m|\phi_i^n\rangle]$$
    where ρ_i is the density of atom i, ρ̄mn is the GSN average density involving the full system density (ρ), and ρ̄_i is the GSN average density involving only the atomic density (ρ_i). The first term (<ϕ_iᵐ|Vxc[ρ_i]|ϕ_iⁿ>) is the exact one-center contribution.

*   **Off-site matrix elements (<ϕ_iᵐ|Vxc[ρ]|ϕ_jⁿ>):**
    $$\langle\phi_i^m|V_{xc}[\rho]|\phi_j^n\rangle \approx \langle\phi_i^m|V_{xc}[\rho_i+\rho_j]|\phi_j^n\rangle + V_{xc}[\bar{\rho}_{mn}]\langle\phi_i^m|\phi_j^n\rangle + V_{xc}'[\bar{\rho}_{mn}][\langle\phi_i^m|\rho r|\phi_j^n\rangle - \bar{\rho}_{mn}\langle\phi_i^m|\phi_j^n\rangle] - V_{xc}[\bar{\rho}_{ij}]\langle\phi_i^m|\phi_j^n\rangle - V_{xc}'[\bar{\rho}_{ij}][\langle\phi_i^m|\rho_{ij}r|\phi_j^n\rangle - \bar{\rho}_{ij}\langle\phi_i^m|\phi_j^n\rangle]$$
    where ρ_i+ρ_j is the sum of densities of atoms i and j, ρ̄mn is the GSN average density involving the full system density (ρ), and ρ̄_ij is the GSN average density involving the sum of atomic densities (ρ_i+ρ_j). The first term (<ϕ_iᵐ|Vxc[ρ_i+ρ_j]|ϕ_jⁿ>) is the exact two-center contribution.

In McWEDA, the GSN approximation is applied *only* to calculate the correction terms, not the dominant one- or two-center contributions.

**General Main Equation and Operation of FIREBALL**

FIREBALL is a **real-space local-pseudoatomic-orbital molecular dynamics implementation of DFT cast in a tight-binding-like form**. Its fundamental goal is to provide a good balance between accuracy and computational efficiency, particularly for complex and large systems where other DFT methods might be computationally prohibitive.

At its core, FIREBALL is based on a **self-consistent extension of the Harris-Foulkes functional**. While standard Kohn-Sham DFT typically achieves self-consistency by iteratively adjusting the output electron density until it matches the input density on a real-space grid, FIREBALL achieves self-consistency on the **orbital occupation numbers**.

The total energy functional in FIREBALL is based on the Harris-Foulkes form (similar to Eq. 1 in Source 167 and Source 395, and Eq. 1 in Source 221, with slightly varying explicit terms across sources):
$$E_{tot} \approx E_{BS} + E_{dc} + E_{ion-ion}$$
where Etot is the total energy, EBS is the band structure energy, Edc is the double counting energy, and Eion-ion is the ion-ion interaction energy. The energy is calculated based on an **input charge density** ρ(r).

The **input charge density ρ(r)** is represented as a **sum of confined spherical atomic-like densities**:
$$\rho(r) = \sum_i \rho_i(r) = \sum_{\mu=(i,l,m)} n_\mu |\phi_\mu(r)|^2$$
Here, ρ_i(r) is the density centered at atom i, ϕ_μ(r) are the local atomic-like orbitals (fireballs), and n_μ are the **orbital occupation numbers**. These occupation numbers determine the electron distribution of the input density.

The band structure energy (EBS) comes from solving the one-electron Schrödinger equation:
$$\hat{H}[\rho(r)] \psi_n = \epsilon_n \psi_n$$
where ψ_n are the eigenstates, ε_n are the eigenvalues, and **Ĥ** is the Hamiltonian, which **depends on the input density ρ(r)**:
$$\hat{H} = \hat{T} + V_{ions} + V_{Ha}[\rho] + V_{xc}[\rho]$$

**How FIREBALL Works and Differs from Other DFT Codes:**

1.  **Basis Set:** FIREBALL uses **localized numerical pseudo-atomic orbitals called "fireballs"**. These orbitals are explicitly constructed to vanish beyond a certain cutoff radius. This is a significant difference from plane-wave DFT codes (which use a delocalized basis tiling the entire space) or Gaussian basis set codes (which use analytical functions like Gaussians). This localization makes the Hamiltonian and overlap matrices **sparse**, which is computationally advantageous for large systems.
2.  **Real-Space and Multicenter Integrals:** Calculations are performed predominantly in **real space**. Key Hamiltonian and overlap matrix elements are evaluated using a **multicenter expansion approach** (typically up to three centers exactly). These integrals are **pre-calculated and stored in data tables** (1D or 2D) and interpolated during the simulation, rather than being calculated on-the-fly. This contrasts with plane-wave methods that work in reciprocal space or other methods that calculate integrals on the fly.
3.  **Input Density Representation:** The input electron density is represented as a **sum of localized, typically spherical, atomic-like densities**. This simple form allows for efficient calculation of Hartree terms without requiring four-center integrals.
4.  **Self-Consistency:** Instead of converging the electron density on a grid, FIREBALL converges the **orbital occupation numbers**. The output occupation numbers (n_μ⁰ᵘᵗ), derived from projecting the self-consistent eigenvectors onto orthogonalized atomic orbitals (Lowdin orbitals), are matched with the input occupation numbers (n_μ):
    $$n_\mu^{out} = 2 \sum_n |\langle\psi_n|\chi_\mu\rangle|^2 f_n$$
    where ψ_n are the occupied eigenvectors, χ_μ are orthogonalized atomic-like orbitals, and f_n are occupation factors. The input occupation numbers (n_μ) in ρ(r) (Eq. 7) are then updated based on n_μ⁰ᵘᵗ until convergence.
5.  **Exchange-Correlation:** FIREBALL employs the **McWEDA method** for calculating exchange-correlation contributions, which improves upon previous approximations while maintaining computational efficiency.
6.  **Forces:** Atomic forces are calculated exactly from the derivative of the total energy with respect to atomic positions. Pulay corrections are naturally included .

These features, particularly the localized basis, real-space multicenter approach with pre-calculated tables, and the orbital occupation self-consistency, contribute to FIREBALL's computational efficiency, making it suitable for molecular dynamics simulations and studies of larger, complex systems that are challenging for more computationally intensive ab initio methods.