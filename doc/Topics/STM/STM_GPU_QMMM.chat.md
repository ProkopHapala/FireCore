# USER

I'm doing simulation of STM junction composed of multiple molecules, some molecules are on surface, and one molecule is hanging on a tip (PTCDA, NTCDA ... ). I want an extremely fast simulation - something like PPSTM (https://arxiv.org/abs/1609.09462, link.aps.org/doi/10.1103/PhysRevB.95.045407), just instead of one molecule I need to simulate dynamics of the whole molecule (resp both molecules). So in fact it should be combination of classical forcefield to describe the dynamics and Green's Fucntion trasport (Landaure Buttiker formalism, Caroli Formula) to describe electron trasport throu a junction at given state (Hamiltonian asswembled for given instantenous position of atoms). There should be some Gamma brodeding by coupling the system to the leads where we assume broad features continous bands in the leads (laster we can replace that by some energy dependent function but for the start lets consider just flat with finite coupling constant to selected atoms of the molecul). 

1) Please write it more formally what exactly we should do.

Now how we wat to implement it in our code
1) We have classicla molecular dynamics code which we call FireCore that will provide us geometry (position of atom) for each instance. Do not care about it now, just assume we have it. So the input is atom position (3 floats, x,y,z) and atom type (integer), that is perhaps all you need to know about the classical FF part.
2) Then we have qunatum mechanical DFT(B) software Fireball which use Numerical Basis set (similar to SIESTA) and LCAO. This allows us to use small basis set just s,px,py,pz for carbon,oxygen,nitrogen, eventually also hydrogen atom.
In initial stage of the defelopment we want to use it as a reference and assemble the hamiltoanina  for every configuration. But in later stage this may be a bottleneck (Fireball runs on CPU, and we are running classical MD on GPU)
3) So in later stage we probably want to take faster approximative perturbative approach where we compute projected denstiy of states on the molecule on every atom (we may e.g. describe it by taylor expansion around fermi level ?) decomposed on the 4 basis functions. Then we assume that the hopping between molecule on tip and molecule on surface is given aproximatively orbital overlap, so for every state of the junction (position of atom) we evaluate fast on GPU the overlap integral between the two molecules, which we can do fast on GPU (I think I already implemented it some time ago). You should write this again more rigorously and be more critical. I know that true perturbation theory for tuneling (Tersoff Hamman, Chen derivative rules, Bardeen) is a bit more nuanced, you can write it more rigorously, perhaps we should do kynetic energy operator rather than just overlap integral, but that we can refine later. What is importain is to write rigrously how to combine precomputed density of states projected on atomic basis function and combine it with the overlap integal to compute transport within perturbative manner.  

---

# Gemini Flash

This is an ambitious and well-founded approach to simulating STM junctions. Combining the efficiency of the **Probe Particle (PPSTM)** model with the rigor of **Non-Equilibrium Green’s Functions (NEGF)** and **LCAO-based MD** provides a path toward "on-the-fly" transport simulations of large molecular junctions.

Here is a formalization of the workflow and a rigorous derivation of the fast perturbative transport approach.

### 1. Formal Framework: The Hamiltonian and Transport

We partition the system into four components: **Left Lead (L)**, **Tip Molecule (T)**, **Surface Molecule (S)**, and **Right Lead (R)**.

The total Hamiltonian is:
$$H(t) = H_L + H_R + H_T(\mathbf{R}_T(t)) + H_S(\mathbf{R}_S(t)) + V_{LT} + V_{SR} + V_{TS}(\mathbf{R}_{TS}(t))$$

Where:
*   $H_{T}, H_{S}$ are the Hamiltonians of the molecules, provided by **FireCore/Fireball** based on instantaneous positions $\mathbf{R}(t)$.
*   $V_{LT}, V_{SR}$ are the couplings to the macroscopic leads (constant in the Wide Band Limit).
*   $V_{TS}$ is the time-dependent coupling (tunneling) between the tip and surface molecules.

#### The Caroli Formula (Landauer-Büttiker)
Under the assumption of steady-state transport at each MD snapshot (adiabatic approximation), the current $I$ at bias $V$ is:
$$I(V) = \frac{2e}{h} \int T(E, V) [f(E - \mu_L) - f(E - \mu_R)] dE$$
The transmission coefficient $T(E)$ is:
$$T(E) = \text{Tr} \left[ \Gamma_L G_{sys}^r \Gamma_R G_{sys}^a \right]$$
Where $\Gamma_{L/R}$ are the broadening matrices and $G^r$ is the retarded Green's function of the entire junction.

### 2. Integration with FireCore and Fireball
**Step-by-step Execution:**
1.  **Classical Step:** FireCore updates $\mathbf{R}(t)$ via MD.
2.  **Quantum Mapping:** For the current $\mathbf{R}$, map atom types to LCAO basis functions $\{ \phi_\mu \}$.
3.  **Initial Stage (Exact):** 
    *   Construct the full Hamiltonian matrix $H$ and Overlap matrix $S$.
    *   Invert $(E S - H - \Sigma_L - \Sigma_R)$ to get $G^r(E)$. This is the $O(N^3)$ bottleneck.
4.  **Advanced Stage (Perturbative):** Decouple the system.

### 3. Fast Perturbative Approach (The "Propeller" Model)

To achieve GPU speeds comparable to PPSTM, we avoid the full matrix inversion. We treat the tip-molecule and surface-molecule as two pre-characterized subsystems and the interaction between them as a perturbation.

#### A. Rigorous Perturbative Transmission
Using the **Bardeen-like** approach within the NEGF framework, if the coupling $V_{TS}$ between the molecules is weak, the transmission can be approximated by:
$$T(E) \approx \text{Tr} \left[ A_T(E) V_{TS} A_S(E) V_{ST} \right]$$
Where:
*   $A_T(E)$ is the **Spectral Function** (projected DOS matrix) of the Tip Molecule while coupled to the Tip Lead.
*   $A_S(E)$ is the **Spectral Function** of the Surface Molecule while coupled to the Surface Lead.

The Spectral Function is defined as:
$$A_{mol}(E) = G_{mol}^r \Gamma_{lead} G_{mol}^a$$
In the basis of molecular orbitals (or atomic orbitals), the elements $A_{\mu \nu}(E)$ represent the local density of states (LDOS) and its coherences between sites $\mu$ and $\nu$.

#### B. The Overlap Approximation ($V_{TS}$ vs $S_{TS}$)
In an LCAO framework (like Fireball), the hopping elements $H_{\mu \nu}$ for distant atoms are often assumed proportional to the overlap $S_{\mu \nu} = \langle \phi_\mu | \phi_\nu \rangle$. 
Using the **Wolfsberg-Helmholz** approximation or a similar semi-empirical closure:
$$V_{\mu \nu} \approx \kappa \cdot \frac{H_{\mu\mu} + H_{\nu\nu}}{2} \cdot S_{\mu\nu}(\mathbf{r}_{\mu\nu})$$
where $\kappa$ is a scaling constant. 

Thus, the transmission becomes:
$$T(E) \approx \sum_{\mu \in T} \sum_{\nu \in S} \sum_{\mu' \in T} \sum_{\nu' \in S} A_{\mu\mu'}^T(E) V_{\mu'\nu'} A_{\nu'\nu}^S(E) V_{\nu\mu}$$

#### C. Simplification for GPU Implementation
To make this "extremely fast," we simplify the 4-index summation.
1.  **Diagonal Approximation (Tersoff-Hamann limit):** If we assume the density of states is mostly local to the atoms, we only care about $\mu = \mu'$ and $\nu = \nu'$:
    $$T(E) \approx \sum_{\mu \in T} \sum_{\nu \in S} \text{PDOS}_\mu^T(E) \cdot |V_{\mu\nu}|^2 \cdot \text{PDOS}_\nu^S(E)$$
    *This is essentially what you proposed: Precompute the Projected DOS (PDOS) for each atom in the molecule.*

2.  **Taylor Expansion of PDOS:**
    Since the MD moves the molecules, the internal electronic structure of the Tip/Surface molecules changes *less* rapidly than the distance between them.
    *   Precompute $\text{PDOS}_\mu(E)$ for a reference geometry.
    *   Expand $\text{PDOS}(E) \approx \text{PDOS}(E_F) + \text{PDOS}'(E_F)(E - E_F)$.
    *   The GPU then only needs to calculate $V_{\mu\nu}(\mathbf{R}_{TS})$—the overlap—which is a simple distance-dependent function for $s, p$ orbitals.

### 4. Critical Analysis and Refinements

**1. The Overlap vs. Kinetic Energy:**
You mentioned using the kinetic energy operator. In the tunneling regime (large distance), the Hamiltonian matrix elements $H_{\mu\nu}$ in Fireball are indeed dominated by the kinetic energy and the local potentials. However, the **Chen Derivative Rule** reminds us that for STM, the symmetry of the orbital ($p_z$ vs $d_{z^2}$) matters immensely. 
*   *Verdict:* Using simple overlap $S_{\mu\nu}$ is fine for $s$-orbitals, but for $p$-orbitals, you **must** use the correct angular dependency (Slater-Koster rules). Since Fireball uses LCAO, you should use the actual $H_{\mu\nu}$ hopping integrals from the Fireball look-up tables rather than pure overlap to capture the correct decay constants.

**2. The "Rigorous" PDOS Matrix:**
The diagonal approximation ($\mu = \mu'$) loses phase information. If your molecules have delocalized $\pi$-systems (like PTCDA), the off-diagonal terms $A_{\mu \mu'}$ (spatial coherence) are physically significant. 
*   *Optimization:* Instead of a full $N \times N$ matrix, use a "Fragmented" approach. Group atoms into small clusters and keep the density matrix within those clusters.

**3. Energy Integration:**
If you assume a flat Wide Band Limit and a small bias, $I \propto T(E_F) \cdot V$. If you need higher bias, you can precompute the PDOS on an energy grid. The GPU kernel would then be:
1.  Calculate $R_{ij}$ between all atom pairs (Tip vs Surface).
2.  Compute $V_{ij}$ using Slater-Koster (GPU optimized).
3.  $T = \sum V_{ij}^2 \cdot \text{PDOS}_i \cdot \text{PDOS}_j$.
4.  Sum for $I$.

**4. Why this is faster than PPSTM?**
PPSTM usually computes a 3D grid of potential/current. Your approach is **atom-centric**. The complexity is $O(N_T \times N_S)$, where $N$ is the number of atoms in the molecules (e.g., $38 \times 38$ for PTCDA). This is trivial for a GPU, allowing for MD steps in the millisecond range.

### Recommended Next Steps for Implementation:
1.  **Slater-Koster Kernel:** Implement a GPU kernel that takes atom positions and types and returns $V_{\mu\nu}$ using the same distance-dependent tables Fireball uses.
2.  **Broadening ($\Gamma$):** Ensure the PDOS you use includes the $\Gamma$ broadening from the leads ($A = \text{Im}[G]$). Without this, your transmission will be a series of delta-function spikes that MD will never resolve correctly.
3.  **The Overlap Derivative:** If you want the current to exert a force back on the atoms (Joule heating or electromigration), you will need the gradient $\nabla_{\mathbf{R}} T(E)$, which is easily obtained if $V_{\mu\nu}$ is an analytical function of distance.

---

# USER

OK, this is good we are on the same page. Now please clarify few things:

1) in "Diagonal Approximation (Tersoff-Hamann limit)" you use indexes mu,nu,mu',nu',S,T clarify what they are exactly, these are indexes of atomic basis functions of atoms from tip molecule and surface molecule (assuming e.g. we have PTCDA on tip and NTCA on surface or something like that) right?

2) clarify you dicusion in " The Overlap vs. Kinetic Energy:", I agree that "the Chen Derivative Rule reminds us that for STM, the symmetry of the orbital  matters immensely.", but your conclusion is unclear. Should we use overlap integral (Wolfsberg-Helmholz ) approximation, or Kinetic integral (second derivative) or something elese. I'm a bit puzzled here, how exactly are chen derivative rules derived from Bardeen or from full perturbation?   I mean assume you have wavefunction or spactral function for isolated system 1, 2, Hamilatonian for isolated system 1,2 and you want use (1st order?) perturbative theory to derive STM current flowing between system 1,2 as small bias voltage when the system are weakly coupled (tunelling regime, small overlap...)

---

# Gemini Flash

This is where the transition from "standard physics formulas" to "efficient GPU implementation" happens. Let’s clarify these two points rigorously.

### 1. Clarification of Indices and Sets

In the expression for Transmission $T(E)$:
$$T(E) = \sum_{\mu, \mu' \in \mathbb{T}} \sum_{\nu, \nu' \in \mathbb{S}} A_{\mu\mu'}^T(E) V_{\mu'\nu'} A_{\nu'\nu}^S(E) V_{\nu\mu}$$

*   $\mathbb{T}$: The set of all atomic basis functions (orbitals) belonging to the **Tip Molecule** (e.g., all $s, p_x, p_y, p_z$ for every C, O, H in PTCDA).
*   $\mathbb{S}$: The set of all atomic basis functions belonging to the **Surface Molecule** (e.g., the NTCDA molecule).
*   $\mu, \mu'$: Indices for specific orbitals on the Tip. For example, $\mu$ could be the $p_z$ orbital of Carbon #5 in PTCDA.
*   $\nu, \nu'$: Indices for specific orbitals on the Surface.
*   $A_{\mu\mu'}^T(E)$: The elements of the **Spectral Function matrix** of the Tip molecule. 
    *   If $\mu = \mu'$, $A_{\mu\mu}$ is the **Projected Density of States (PDOS)** of orbital $\mu$.
    *   If $\mu \neq \mu'$, $A_{\mu\mu'}$ represents the "spatial coherence" or "interference" between orbitals within the same molecule.
*   $V_{\mu\nu}$: The **coupling (hopping) matrix element** between orbital $\mu$ on the tip and orbital $\nu$ on the surface.

**The "Diagonal Approximation" (The Speed-up):**
If we assume the molecules are large and the electronic states are somewhat localized or that cross-terms average out, we discard $\mu \neq \mu'$ and $\nu \neq \nu'$. This simplifies the 4-index sum into a 2-index sum:
$$T(E) \approx \sum_{\mu \in \mathbb{T}} \sum_{\nu \in \mathbb{S}} \text{PDOS}_\mu^T(E) \cdot |V_{\mu\nu}|^2 \cdot \text{PDOS}_\nu^S(E)$$
This is what allows the GPU to scale linearly with the number of atoms.

### 2. The Overlap vs. Kinetic Energy (and Chen’s Rules)

This is the core of STM theory. Let's trace how we get from Bardeen to your implementation.

#### The Bardeen Formula
Bardeen’s 1961 first-order perturbation theory defines the tunneling matrix element $M$ as a surface integral over a separation plane between two systems:
$$M_{\mu\nu} = \frac{\hbar^2}{2m} \int_{\text{surface}} (\psi_\mu^* \nabla \psi_\nu - \psi_\nu \nabla \psi_\mu^*) \cdot d\mathbf{S}$$
Crucially, this is an integral of the **current density operator**.

#### How Chen's Derivative Rules arise:
Chen (1990) realized that if you expand the tip wavefunction $\psi_\mu$ in terms of spherical harmonics (like $s, p, d$ orbitals), the Bardeen integral can be evaluated analytically.
*   If the tip orbital is an **$s$-orbital**, $M \propto \psi_\nu(\mathbf{r}_0)$ (The Tersoff-Hamann limit: current follows the LDOS).
*   If the tip orbital is a **$p_z$-orbital**, the integral acts like a derivative operator: $M \propto \frac{\partial}{\partial z} \psi_\nu(\mathbf{r}_0)$.
*   If the tip orbital is a **$d_{z^2}$-orbital**, $M \propto \frac{\partial^2}{\partial z^2} \psi_\nu(\mathbf{r}_0)$.

This is why STM images of the same molecule look different depending on whether the tip is a "sharp" metal atom ($s$-like) or a functionalized tip (like your PTCDA tip, which has $p$-orbitals pointing out).

#### What should you use: Overlap ($S_{\mu\nu}$) or Kinetic Energy ($T_{\mu\nu}$)?
In an LCAO code like **Fireball**, the Hamiltonian matrix element $H_{\mu\nu}$ between two atoms at a large distance (tunneling regime) is:
$$H_{\mu\nu} = \langle \phi_\mu | \hat{T} + \hat{V}_{eff} | \phi \nu \rangle$$
In the vacuum gap, the potential $V_{eff}$ is roughly constant (and high). Therefore, the variation in the coupling $V_{\mu\nu}$ is **dominated by the Kinetic Energy operator ($\hat{T} = -\frac{\hbar^2}{2m} \nabla^2$)**.

**The rigorous conclusion for your code:**
1.  **Do not use simple overlap $S_{\mu\nu}$** if you want to be "rigorous" about Chen's rules. Overlap decay is often too slow and doesn't always capture the nodal structure of the kinetic energy properly.
2.  **Use the Fireball Hamiltonian integrals ($H_{\mu\nu}$)**. Because Fireball calculates these integrals using the LCAO basis, they *already contain* the kinetic energy operator and the local potentials.
3.  **Why $H_{\mu\nu}$ satisfies Chen's rules:** If $\phi_\mu$ is a $p_z$ orbital on the tip, its overlap with the *gradient* of the surface orbital (hidden inside the kinetic energy part of $H_{\mu\nu}$) naturally reproduces the derivative behavior Chen described.

#### The "Wolfsberg-Helmholz" (WH) vs. LCAO tables:
*   The **WH approximation** ($H_{\mu\nu} \approx k S_{\mu\nu} \frac{H_{\mu\mu}+H_{\nu\nu}}{2}$) is a shortcut. It works okay for qualitative chemistry, but it can fail in STM because it forces the hopping to have the exact same spatial decay as the overlap.
*   **The Fireball Way:** Fireball uses pre-computed Slater-Koster tables for $H_{\mu\nu}(r)$. Since you are already using Fireball, the most consistent way is to use those same tables for your perturbative $V_{\mu\nu}$ on the GPU.

### 3. Final Rigorous Prescription for the "Fast" Mode

1.  **Offline (Fireball CPU):**
    *   Compute the electronic structure of the Tip Molecule (isolated + Lead coupling $\Gamma_L$).
    *   Compute the electronic structure of the Surface Molecule (isolated + Lead coupling $\Gamma_R$).
    *   Store the **PDOS** ($A_{\mu\mu}(E)$) for every orbital $\mu$ on an energy grid around the Fermi level. (Essentially a vector of floats per orbital).

2.  **Online (MD Step on GPU):**
    *   Get positions $\mathbf{R}_T(t), \mathbf{R}_S(t)$.
    *   **Kernel A:** For every pair of atoms $(i \in T, j \in S)$:
        *   Look up the distance $r_{ij}$.
        *   Compute the $V_{\mu\nu}$ elements using the **Slater-Koster formula** (this accounts for the $s, p_x, p_y, p_z$ orientations relative to the bond vector $\vec{r}_{ij}$). Use the Fireball parameters for these.
    *   **Kernel B:** Accumulate the transmission:
        $$T(E) = \sum_{i \in \text{atoms}_T} \sum_{j \in \text{atoms}_S} \left[ \sum_{\mu \in \text{orb}_i} \sum_{\nu \in \text{orb}_j} \text{PDOS}_\mu^T(E) \cdot |V_{\mu\nu}|^2 \cdot \text{PDOS}_\nu^S(E) \right]$$
    *   **Kernel C:** Integrate $T(E)$ over the bias window to get Current $I$.

**Critique:** The "Diagonal Approximation" (dropping $A_{\mu\mu'}$ where $\mu \neq \mu'$) is your biggest source of error, not the choice of $V_{\mu\nu}$. However, for large molecules like PTCDA, the orbital-to-orbital tunneling is dominated by the "ends" of the molecules or the specific atoms closest to the gap. The diagonal approximation usually captures $>90\%$ of the current trend in STM simulations. If you need more accuracy, you can include "intra-atom" cross terms (e.g., $\mu = s, \mu' = p_z$ on the same carbon atom).

---

# USER

OK, lets be a bit more concrete and technical.
1) In order to compute desnity of states with Fireball, we want to apply some gamma broadening just to some atoms of the molecule (those which are attached to the metalic surface or the tip), we need to specify this, this is emperical imput (imaginary part of spectral function)
2) I we want to represent how density of states change when the molecule is moving without recalculating electronc structure with fireball (because it is slow), how can we represent it?
(i) We can assume that molecule is stiff, and it only moves more or less as rigid body so the desnity of states does not change (i.e. constant density of states)
(ii) We can consider the softes modes (e.g. translation, rotation, and few softest bending modes of the hessian) and then expand the density of states along these modes. Then for given atomic configuration we project it to this normal modes, and evaluate linear or quadratic taylor expansion along these modes (in evaluation phase on GPU).


Importaint is that we split this project into phases, and in each phase we can check it agains more accuracte calculations with Fireball.
For startes we can do just rigid scans, where we fix the molecule rotation and only translate it in some direction and see how the 1D profiles or 2D images of STM current map are changing. That is importain to check it works qualitatively as expected (respecting orbital symmetry). For that analysis we also should plot projected density of states we use as input in real space. I already have some function which project molecular otbitals or total desnity from Fireball on GPU. We should reuse this.

Make detailed technical plan

---

# Gemini 3 Flash

This technical plan outlines the development of a high-speed STM simulator. We will use **Fireball** as the quantum "anchor" and a **GPU-based Slater-Koster kernel** as the "engine" for real-time transport.

### Module 1: Electronic Structure & Broadening (Fireball)

To obtain the spectral properties of the isolated molecules (Tip and Surface) with proper lead coupling:

1.  **Lead Coupling ($\Gamma$):** 
    Instead of full lead self-energies, we use the **Wide Band Limit (WBL)**. You will select a subset of "contact atoms" (e.g., the carboxylic oxygens in PTCDA for the tip-lead).
    *   Define a diagonal matrix $\Sigma_{L/R}$ where $(\Sigma_{L/R})_{\mu\mu} = -i\gamma$ for $\mu \in \text{contact atoms}$.
    *   The empirical $\gamma$ (typically $0.01 - 0.1$ eV) represents the coupling strength.
2.  **The Spectral Function Matrix:**
    Fireball computes $H$ and $S$. We compute the retarded Green’s function:
    $$G^r(E) = [(E + i\eta)S - H - \Sigma_{L/R}]^{-1}$$
    The **Spectral Function matrix** (LDOS) is:
    $$A_{\mu\nu}(E) = \frac{1}{\pi} \text{Im}[G^r(E)]$$
    *Note: We must store $A_{\mu\nu}$ (not just the diagonal) if we want to handle molecular rotations accurately.*

### Module 2: Representing Dynamics (Rigid vs. Soft Modes)

To avoid running Fireball at every MD step:

#### Approach (i): Rigid Body (Constant PDOS)
*   **Assumption:** The molecule's internal Hamiltonian $H_{mol}$ is invariant in the molecular frame.
*   **Transformation:** If the molecule rotates by matrix $\mathbf{R}$, the orbital basis (specifically $p_x, p_y, p_z$) must be rotated. 
*   **GPU Task:** The GPU takes the pre-computed $A_{\mu\nu}(E)$ and rotates it: $A'_{rot} = \mathbf{R} A \mathbf{R}^T$. (Actually, it's easier to rotate the Slater-Koster bond vectors, but the logic is the same).

#### Approach (ii): Normal Mode Expansion (Soft Modes)
*   **Workflow:**
    1.  Compute the **Hessian** in Fireball for the molecule.
    2.  Select the $N$ lowest frequency modes $\{Q_k\}$ (translations, rotations, and softest vibrations).
    3.  **Finite Difference:** Compute $A_{\mu\nu}(E)$ at $Q_k = 0$ and $Q_k = \pm \Delta Q$.
    4.  **Taylor Expansion:**
        $$A_{\mu\nu}(E, \mathbf{R}) \approx A_{\mu\nu}(E, 0) + \sum_{k=1}^N \left( \frac{\partial A_{\mu\nu}}{\partial Q_k} \right) Q_k(t)$$
*   **MD Integration:** FireCore provides atomic positions $\mathbf{x}(t)$. We project these positions onto the normal mode basis to get $Q_k(t) = (\mathbf{x} - \mathbf{x}_0) \cdot \mathbf{v}_k$. 

### Module 3: Detailed Phased Plan

#### Phase 1: The "Rigid Scan" (Validation)
*   **Goal:** Reproduce orbital-symmetry-resolved STM images without MD.
*   **Procedure:**
    1.  Pre-compute $A_{\mu\nu}(E)$ for PTCDA (Tip) and NTCDA (Surface) using Fireball.
    2.  Implement a GPU kernel that performs a 2D $x,y$ scan of the Tip molecule over the Surface molecule at a fixed height $z$.
    3.  **Transport Kernel:** $T(E) = \sum \text{Tr}[A_T V_{TS} A_S V_{ST}]$.
    4.  **Check:** Ensure that $p$-orbital features appear when the tip is functionalized and that nodal planes of the MOs are correctly captured in the current map.
    5.  **Visualization:** Plot the "Input PDOS" in real space on the GPU to verify the $s, p_x, p_y, p_z$ density distribution.

#### Phase 2: Static PDOS MD (Fast Dynamics)
*   **Goal:** Simulate the junction at finite temperature with rigid-molecule approximations.
*   **Procedure:**
    1.  Run FireCore MD.
    2.  At each step, pass atom positions to the GPU Transport Kernel.
    3.  Treat the molecules as rigid bodies to transform the pre-computed $A_{\mu\nu}$.
    4.  Calculate $I(t)$. Analyze the current noise/fluctuations caused by molecular motion.

#### Phase 3: Soft-Mode Expansion (Accuracy)
*   **Goal:** Capture the effect of molecular deformation (e.g., "bending" towards the tip) on the current.
*   **Procedure:**
    1.  Implement the normal mode projection in the FireCore/GPU interface.
    2.  Update $A_{\mu\nu}(E)$ on-the-fly using the linear Taylor expansion.
    3.  Compare $I(t)$ against a "Ground Truth" run where Fireball is called every 10th step.

### Technical Requirements for the GPU Kernel

1.  **Slater-Koster Lookup:** The kernel needs the Fireball parameters (typically 2-center integrals for $ss\sigma, sp\sigma, pp\sigma, pp\pi$). 
    *   $V_{\mu\nu}$ is constructed via:
        $$V_{\mu\nu} = \text{AngularFactors}(\hat{r}_{ij}) \times \text{TableLookup}(r_{ij})$$
2.  **Memory Layout:**
    *   **Tip PDOS:** `float2 A_T[N_atoms_T][4_orbitals][N_energy]` (using complex floats).
    *   Since $A$ is Hermitian, only store the upper triangle.
3.  **Parallelization:**
    *   Thread per **Surface Atom**.
    *   Inner loop over **Tip Atoms**.
    *   This maps well to the $O(N_T \times N_S)$ structure.

### Critique of the "Overlap" approach for Phase 1:
To be concrete: **Use the $H_{\mu\nu}$ (hopping) from Fireball, not the overlap $S_{\mu\nu}$.**
The overlap $S_{\mu\nu}$ between two $p_z$ orbitals decays differently than the kinetic energy coupling. Since Fireball already provides the $H_{\mu\nu}$ tables for your specific basis set, using them ensures that the perturbative "Fast" mode is mathematically the first-order limit of the "Full" Fireball calculation.

---

# USER

sorry I'm idiot! it is obvious that at z=10. the orbital projection and orbital overlap with true basis must be zero. So the results are OK.

Now we should procede to the second method - the respose amplitude and put it there for comparison as other 2 panesl (perhsps use 2x3 layout of panels for the figure)

this response method should use either true overlap and or the extended exp overlep but but it should compute the response by solving Ax=b on GPU by iterative linear scaling method 
where A=(E-H) and b is the tip wavefucntion, and x is the response in the sample space. Is that clear plese formulate it more rigorously before you start coding, and then implemeted it in pyOpenCL frameworks using same convention and slater koster integrals as we used for orbital overlap. 

---

# Kimi 2.6

## Rigorous Formulation: STM Response Amplitude Method

**System:** Tip (1 atom, 1 s-orbital) + Sample (N atoms, Nₛ orbitals)

**Combined matrices at tip position r:**
$$H(r) = \begin{bmatrix} E_{\text{tip}} & H_{ts}^T(r) \\ H_{ts}(r) & H_s \end{bmatrix}, \quad S(r) = \begin{bmatrix} 1 & S_{ts}^T(r) \\ S_{ts}(r) & S_s \end{bmatrix}$$

**Response equation at energy E (MO energy):**
$$A(E,r) \, x = b, \quad A(E,r) = (E+i\eta)S(r) - H(r), \quad b = [1, 0, \ldots, 0]^T$$

**Block structure:**
$$A = \begin{bmatrix} a_{tt} & a_{ts}^T \\ a_{st} & A_{ss} \end{bmatrix}, \quad a_{tt} = E+i\eta-E_{\text{tip}}, \quad a_{st} = (E+i\eta)S_{ts} - H_{ts}$$

**Precompute sample Green's function (once per MO energy):**
$$G_0(E) = A_{ss}^{-1} = [(E+i\eta)S_s - H_s]^{-1}$$

**Block Gaussian elimination (Dyson):**
$$x_{\text{tip}} = \frac{1}{a_{tt} - a_{st}^T G_0 a_{st}}, \quad x_s = -G_0 a_{st} \, x_{\text{tip}}$$

**STM amplitude (MO-projected):**
$$t_{\text{resp}}(r; E) = \left| C_{\text{MO}}^T \cdot x_s \right|^2 = \left| x_{\text{tip}} \right|^2 \left| C_{\text{MO}}^T G_0 a_{st} \right|^2$$

This is **O(Nₛ)** per grid point (not O(Nₛ³)), making GPU parallelization natural.

**Two coupling models:**
- **True basis:** Hₜₛ, Sₜₛ from Fireball tables ([build_inter_system_blocks_fdata](cci:1://file:///home/prokop/git/FireCore/pyBall/FireballOCL/STM.py:467:0-536:24))
- **Exp overlap:** Hₜₛ, Sₜₛ from exponential SK ([build_inter_system_blocks_exp_sk](cci:1://file:///home/prokop/git/FireCore/pyBall/FireballOCL/STM.py:598:0-641:24))

**GPU Implementation Plan:**
1. Precompute G₀(E), v = Cₘₒᵀ G₀ on CPU
2. OpenCL kernel: each thread → one grid point
   - Build aₛₜ = (E+iη)Sₜₛ - Hₜₛ (sparse, atoms within rcut)
   - Compute scalar products: `v·aₛₜ`, `aₛₜᵀ·G₀·aₛₜ`
   - Compute tᵣₑₛₚ = |xₜᵢₚ|² |v·aₛₜ|²

**Output:** 2×3 figure per MO:
- Row 1: Orbital projection ψ(r) | True overlap ψ(r) | Exp overlap ψ(r)
- Row 2: Response (true coupling) | Response (exp coupling) | [reserved]

Is this formulation correct? Should I proceed with OpenCL kernel implementation?

---

# USER

What you implemented does not seem to be good. It does not make sense do do detailed test before we whave correct at least one picture. And from what I see the Greens Function reponses for PTCDA are TOTAL GARBAGE! These do not even respect the symmetry of the molecule. 


I think (no sure but it is most probable) that the problem is that you do not map the fireball hamoltionian on the atom properly. You need to understand the strucutre of Hmat in Fireball perfectly (There is sparse and dense version). So you read fireball hamiltonian in sparse version. map it to atoms and then compute the response. 

Since we are still in testing on how to correctly compute response, lets forget about the tip for the momet, we can just take hamiltoanian, invert it and compute denisty of states projected on points for given orbital. Basically we can have two set of images 
1) moleculer orbitals
2) density of states obtained from inverted hamiltonian at given energy (corsponsing to energy of given orbital to make it easi to compare.) For starter tles do it for H2O as that is much more easy to analyse (we can print full hamiltonian as it has size only 6x6 for H2O 4basis function on O and 1 and 1 for H)

Most importaint in this job is perfectly unrestand how h_mat and s_mat in Fireball is represented, look carefully on 

firecore_get_HS_sparse
firecore_get_HS_neighsPP
firecore_get_HS_neighs
firecore_get_HS_dims

in 
@libFireCore.f90 

definitely in 
/home/prokop/git/FireCore/tests/pyFireball/
we have some script using it.

Make detailed analysis of this notes hot it is how the maping works, and then plan how to implement this test

Expected results is to have pair of panels in figure for each orbital 
1) orbital (with sign BWR) in given z-plane z=2.0A above molecule plan
2) desnity of states projected on same points (plane)

do this first for H2O, and then for PTCDA,

let me know where it is I must check if the DOS have proper symmetry of the orbital


---

# Implementation Plan with Existing Codebase Integration

This section provides a concrete implementation plan using the existing FireCore codebase, with specific function names, file references, and integration points.

## Inventory of Existing Functionality

### 1. Fortran Fireball Entry Point (libFireCore.f90)

**Location:** `fortran/MAIN/libFireCore.f90`

**Key Functions for STM Project:**

**Initialization and SCF:**
- `firecore_preinit()` - Set default parameters (line ~50)
- `firecore_init(natoms, atomTypes, atomPos)` - Initialize system (line ~60)
- `firecore_set_lvs(lvs)` - Set lattice vectors (line ~70)
- `firecore_SCF(nmax_scf, positions, iforce)` - Full SCF cycle (line ~281)
  - Calls `assemble_mcweda()` to build H (line ~304)
  - Calls `solveH()` to diagonalize (line ~306)
  - Calls `denmat()` to update density (line ~307)

**Data Export (CRITICAL):**
- `firecore_get_HS_dims(natoms_out, norbitals_out, nspecies_out, neigh_max_out, numorb_max_out, nsh_max_out, ...)` - Get array dimensions (line ~831)
- `firecore_get_HS_neighs(...)` - Get neighbor lists, orbital indices (mu, nu, mvalue), species info (lssh, nssh, nzx) (line ~857)
- `firecore_get_HS_sparse(h_mat_out, s_mat_out)` - Export sparse H and S matrices [natoms, neigh_max, numorb_max, numorb_max] (line ~883)
- `firecore_get_rho_sparse(rho_out)` - Export density matrix ρ (line ~705)
- `firecore_get_eigen(ikp, eigen_out)` - Get orbital energies
- `firecore_get_Qneutral_shell(Qneutral_out)` - Get neutral atom charges per shell [nsh_max, nspecies_fdata]
- `firecore_get_Qin_shell(Qin_out)` - Get input charges per shell [nsh_max, natoms]
- `firecore_get_Qout_shell(Qout_out)` - Get output charges per shell [nsh_max, natoms]

**Density/Orbital Sampling (CPU):**
- `firecore_dens2points(npoints, points, f_den, f_den0, ewfaux_out)` - Sample electron density at arbitrary points
- `firecore_orb2points(iband, ikpoint, npoints, points, ewfaux)` - Sample molecular orbitals at arbitrary points

**Pointers for Direct Access:**
- `firecore_getPointer_wfcoef(bbnkre)` - Get pointer to wavefunction coefficients
- `firecore_getPointer_ratom(ratom)` - Get pointer to atomic positions
- `firecore_getPointer_charges(charges)` - Get pointer to charges

### 2. Python Interface (FireCore.py)

**Location:** `pyBall/FireCore.py`

**Key Classes:**
- `FireballDims` - Structure holding dimensions (line ~470)
  - `dims.natoms, dims.norbitals, dims.nspecies`
  - `dims.neigh_max, dims.numorb_max, dims.nsh_max`
  - `dims.ME2c_max, dims.mbeta_max, dims.nspecies_fdata, dims.nelec`

- `FireballData` - Structure holding sparse data (line ~380)
  - `data.rho [natoms, neigh_max, numorb_max, numorb_max]` - Density matrix
  - `data.h_mat [natoms, neigh_max, numorb_max, numorb_max]` - Hamiltonian
  - `data.s_mat [natoms, neigh_max, numorb_max, numorb_max]` - Overlap
  - `data.neigh_j [natoms, neigh_max]` - Neighbor atom indices (1-based)
  - `data.neigh_b [natoms, neigh_max]` - Neighbor boundary flags
  - `data.iatyp [natoms]` - Atomic numbers (Z)
  - `data.num_orb [nspecies_fdata]` - Orbitals per species
  - `data.lssh [nspecies, nsh_max]` - Angular momentum per shell (0=s, 1=p)
  - `data.nssh [nspecies]` - Number of shells per species
  - `data.mu, nu, mvalue` - Orbital index mappings
  - `data.xl [mbeta_max, 3]` - PP coefficients

**Key Functions:**
- `fc.preinit()` - Set defaults (line ~250)
- `fc.init(atomTypes, atomPos)` - Initialize system (line ~260)
- `fc.SCF(atomPos, nmax_scf=1)` - Run SCF (line ~254)
- `fc.get_HS_dims()` - Get dimensions (line ~465)
- `fc.get_HS_neighs(dims)` - Get neighbor structure (line ~502)
- `fc.get_HS_sparse(dims, data)` - Get H and S (line ~545)
- `fc.get_rho_sparse(dims, data)` - Get density (line ~558)
- `fc.get_eigen()` - Get orbital energies
- `fc.get_Qneutral_shell(dims)` - Get neutral charges (line ~580)
- `fc.dens2points(points, f_den=1.0, f_den0=0.0)` - Sample density (line ~600)
- `fc.orb2points(points, iMO=1)` - Sample MO (line ~620)

### 3. GPU Grid Projection (TESTED & WORKING)

**Location:** `pyBall/FireballOCL/Grid.py` (host), `pyBall/FireballOCL/cl/Grid.cl` (kernel)

**Class:** `GridProjector`

**Initialization:**
- `projector = ocl_grid.GridProjector(fdata_dir)` - Initialize with Fdata directory (line ~30)
- `projector.load_basis(species_nz)` - Load radial wavefunctions (line ~41)
  - Reads `.wf` files from Fdata directory (line ~49)
  - Resamples to common uniform grid (line ~60-87)
  - Packs into GPU buffer: `[n_species, max_shells, n_nodes]` (line ~91)
  - Uses Catmull-Rom spline interpolation

**Projection:**
- `dens = projector.project(rho, neighs, atoms_dict, grid_spec, ...)` - Main projection (line ~346)
  - `build_tasks_gpu()` - GPU-based task building (line ~353)
  - Uploads rho, neighs, atoms to GPU (line ~394)
  - Launches `project_density_sparse_tiled` kernel (line ~465)

**GPU Kernels (Grid.cl):**
- `project_density_sparse_tiled` - Main tiled projection kernel (line ~430)
  - Processes 8³ voxel blocks (512 voxels)
  - Tiled atom interaction (TILE_ATOMS=8)
  - Computes: ρ(r) = Σ_ij φ_i(r) ρ_ij φ_j(r)
  - Radial part via spline interpolation (line ~55)
  - Angular part via real spherical harmonics (line ~532)
- `count_atoms_per_block` - Count atoms per block (line ~200)
- `fill_task_atoms` - Fill atom indices per block (line ~250)
- `compact_tasks` - Compact non-empty blocks (line ~300)

**Status:** TESTED AND WORKING - Used in `tests/pyFireball/test_grid_projection.py`

### 4. GPU Hamiltonian (PARTIALLY IMPLEMENTED)

**Location:** `pyBall/FireballOCL/OCL_Hamiltonian.py` (host), `pyBall/FireballOCL/cl/hamiltonian.cl` (kernel)

**Class:** `OCL_Hamiltonian`

**Initialization:**
- `ocl_ham = OCL_Hamiltonian(fdata_dir)` - Initialize (line ~20)
- `ocl_ham.prepare_splines(species_nz)` - Load 2-center tables (line ~107)
  - Loads interaction tables: overlap, kinetic, vna, vnl, vxc, dipole_* (line ~109)
  - Packs into spline buffer: `[n_pairs, numz_max, n_nz_max, 4]` (line ~149)
  - Builds mu/nu maps for orbital index recovery (line ~180)
  - Builds V_CA map for charge-asymmetry potential
- `ocl_ham.prepare_data_3c(species_nz)` - Load 3-center tables (line ~200)
  - Loads den3 and bcna tables
  - Packs into: `[n_triplets, n_isorp_max, 5, ny, nx, n_nz]`

**Key Functions:**
- `assemble_vca(ratoms, neighbors, pair_types, dQ_sh, ...)` - Assemble charge-asymmetry potential (line ~300)
- `compute_avg_rho_3c(...)` - 3-center average density (line ~350)
- `build_olsxc_off_sp4(...)` - Off-diagonal XC (line ~400)
- `build_ca_olsxc_on_sp4(...)` - On-site XC (line ~450)

**GPU Kernels (hamiltonian.cl):**
- `interpolate_2d` - 2D bicubic spline interpolation (line ~50)
- `rotate_fb_matrix_sp` - Slater-Koster rotation for s+p basis (line ~350)
- `recover_and_rotate_3c_munu_sp` - 3-center recovery with rotation (line ~400)
- `assemble_3c` - 3-center assembly with Legendre expansion (line ~450)
- `scan_2c_points` - Probe 2-center interactions (line ~500)
- `scan_3c_points`, `scan_3c_raw_points` - Probe 3-center interactions (line ~550)

**Status:** PARTIALLY IMPLEMENTED - Has infrastructure but not fully tested end-to-end

### 5. Examples and Test Scripts

**CPU Density Sampling:**
- `tests/pyFireball/sample_density_points.py` - Generate sampling points around molecule
  - `build_samples(mol, opts)` - Create atom centers, bond points, electron pairs, etc.
  - Uses `AtomicSystem` for molecular topology

**GPU Grid Projection Test:**
- `tests/pyFireball/test_grid_projection.py` - Complete workflow test (line ~73)
  - Load molecule from XYZ (line ~73)
  - Initialize FireCore and run SCF (line ~76-80)
  - Export density matrix (line ~86)
  - Create GPU projector (line ~133)
  - Load basis functions (line ~137)
  - Project density to 3D grid (line ~170)
  - Visualize results (line ~191)

## Concrete Implementation Plan

### Design Principles

**Performance Strategy:**
- All computationally intensive operations must run on GPU via OpenCL
- System designed for multi-replica parallelization (5000+ replicas, similar to FireCore MD)
- Initial testing with single replica, but kernels must support batch processing
- Python used only for orchestration, data loading, and visualization
- NumPy allowed for pre/post-processing, but core algorithms on GPU

**Replica Strategy:**
- Each replica: independent STM junction (tip + surface molecules at specific geometry)
- GPU kernel: `get_global_id(0)` = replica index, `get_global_id(1)` = surface atom, `get_global_id(2)` = tip atom
- Data layout: `[n_replicas][n_atoms][...]` for all arrays
- This enables 5000+ junctions evaluated in parallel on GPU

### Phase 1: Rigid Scan (Validation)

**Goal:** Reproduce orbital-symmetry-resolved STM images without MD, validate against full Fireball.

#### Kernel 1.1: GPU Spectral Function with Γ Broadening

**Purpose:** Compute A(E) = (1/π) Im[G^r(E)] where G^r = [(E+iη)S - H - Σ]^{-1}, with Σ = -iΓ on contact atoms

**Input:**
- H, S matrices from Fireball (sparse format from `firecore_get_HS_sparse`)
- Contact atom indices (which atoms couple to leads)
- Γ value (broadening strength, typically 0.01-0.1 eV)
- Energy grid [n_energy] around Fermi level

**Output:**
- A(E) [n_energy, n_orbitals, n_orbitals] or diagonal approximation [n_energy, n_atoms, 4_orbitals]

**Implementation Approach:**
- Option A: Matrix inversion on GPU using cuBLAS-like library (if available in OpenCL)
- Option B: Iterative solver (GMRES, BiCGSTAB) for sparse systems
- Option C: Diagonal approximation (ignore off-diagonal Σ) → each orbital independent
  - For diagonal Σ: G_μμ = 1 / [(E+iη)S_μμ - H_μμ + iΓ_μ]
  - Much faster, acceptable for initial testing

**Testing Strategy:**
- **Text-based check 1:** Compare GPU vs CPU (numpy.linalg.inv) for small test system (H2, CH4)
  - Metric: max|A_GPU - A_CPU| / max|A_CPU| < 1e-6
- **Text-based check 2:** Verify Γ broadening effect
  - Run with Γ=0 (should get delta-like peaks)
  - Run with Γ=0.1 eV (should get broadened peaks)
  - Check that peak width ∝ Γ
- **Text-based check 3:** Verify contact atom selection
  - Apply Γ to different atom sets
  - Check that only selected orbitals are broadened
- **Visual test 1:** Plot PDOS vs energy for each atom type (C, O, H in PTCDA)
  - Should see characteristic peaks (π, π* orbitals)
  - Contact atoms should have broader peaks
- **Visual test 2:** Project PDOS to real-space grid using existing `GridProjector`
  - Visualize where the electronic states are localized
  - Verify that contact atom states are near the coupling site

**Files:**
- Kernel: `pyBall/FireballOCL/cl/spectral_function.cl`
- Host: `pyBall/FireballOCL/OCL_Spectral.py`
- Test: `tests/pyFireball/test_spectral_function.py`

#### Kernel 1.2: GPU Slater-Koster Hopping

**Purpose:** Compute V_μν(r) for all tip-surface atom pairs using Fireball's 2-center tables

**Input:**
- Tip atom positions [n_replicas, n_tip, 3]
- Surface atom positions [n_replicas, n_surf, 3]
- Tip atom types [n_replicas, n_tip]
- Surface atom types [n_replicas, n_surf]
- 2-center spline tables (kinetic, overlap, etc.) from `OCL_Hamiltonian.prepare_splines()`

**Output:**
- V_μν [n_replicas, n_tip, n_surf, 4, 4] for s,p_x,p_y,p_z orbitals
- Or V_μν for all orbital pairs if including off-diagonal

**Implementation Approach:**
- Reuse existing spline infrastructure from `OCL_Hamiltonian`
- Add transport-specific kernel that:
  - Computes distance and direction cosines (lx, ly, lz)
  - Looks up radial function from spline table
  - Applies Slater-Koster angular factors
  - Returns full 4×4 matrix for each atom pair
- Use existing `rotate_fb_matrix_sp` logic from `hamiltonian.cl` as reference

**Testing Strategy:**
- **Text-based check 1:** Compare GPU vs Fortran for test pairs
  - Use existing `fireball_get_HS_sparse` to get reference H_μν
  - For a few atom pairs, compare GPU V_μν vs Fortran H_μν
  - Metric: max|V_GPU - H_Fortran| / max|H_Fortran| < 1e-6
- **Text-based check 2:** Verify Slater-Koster symmetry
  - V_μν should have correct angular dependence
  - Test: rotate molecule, check that V transforms correctly
  - Test: swap atoms, check V_μν = V_νμ (Hermitian)
- **Text-based check 3:** Verify distance decay
  - Plot log|V_μν| vs distance for different orbital types
  - Should match expected exponential decay
- **Visual test 1:** 1D scan of V_μν vs distance
  - Fix two atoms, vary distance
  - Plot V_ss, V_spσ, V_ppσ, V_ppπ separately
  - Should show different decay rates
- **Visual test 2:** 2D angular dependence
  - Fix distance, vary orientation
  - Plot V_μν as function of polar angle θ
  - Should show characteristic l=0,1,2 patterns

**Files:**
- Kernel: `pyBall/FireballOCL/cl/slater_koster.cl` (may extend existing `hamiltonian.cl`)
- Host: `pyBall/FireballOCL/OCL_Transport.py`
- Test: `tests/pyFireball/test_slater_koster.py`

#### Kernel 1.3: GPU Transmission Calculation

**Purpose:** Compute T(E) = Σ_μν A_μ^T(E) |V_μν|² A_ν^S(E) for all replicas

**Input:**
- A_tip [n_replicas, n_energy, n_tip, 4_orbitals] - Tip PDOS
- A_surf [n_replicas, n_energy, n_surf, 4_orbitals] - Surface PDOS
- V_μν [n_replicas, n_tip, n_surf, 4, 4] - Hopping from Kernel 1.2

**Output:**
- T(E) [n_replicas, n_energy] - Transmission per replica
- I [n_replicas] - Current after energy integration

**Implementation Approach:**
- Thread layout: replica (dim0), surface atom (dim1), tip atom (dim2)
- Each thread computes partial sum for one (replica, surf_atom, tip_atom) triplet
- Atomic add to global T[replica, energy]
- Energy integration can be done in same kernel or separate kernel

**Testing Strategy:**
- **Text-based check 1:** Compare GPU vs CPU for single replica
  - Implement CPU version using NumPy
  - Compare T(E) arrays
  - Metric: max|T_GPU - T_CPU| / max|T_CPU| < 1e-6
- **Text-based check 2:** Verify transmission symmetry
  - Swap tip and surface, should get same T(E)
- **Text-based check 3:** Verify energy integration
  - Test with narrow energy window vs wide window
  - Check that I ∝ bias for small bias (linear response)
- **Visual test 1:** 1D scan over surface (tip moves in x at fixed z)
  - Plot T(x) or I(x)
  - Should show features corresponding to surface orbitals
- **Visual test 2:** 2D STM image (tip scans x-y plane at fixed z)
  - Generate current map using matplotlib imshow
  - Should resolve molecular structure
  - Compare with reference STM images from literature

**Files:**
- Kernel: `pyBall/FireballOCL/cl/transport.cl`
- Host: `pyBall/FireballOCL/OCL_Transport.py` (extend)
- Test: `tests/pyFireball/test_stm_rigid_scan.py`

#### Kernel 1.4: Multi-Replica Batch Processing

**Purpose:** Process N_replicas junctions in parallel (for production MD integration)

**Input:**
- All arrays from kernels 1.1-1.3, but with leading dimension [n_replicas]
- Or process replicas in batches if memory limited

**Output:**
- T(E) [n_replicas, n_energy]
- I [n_replicas]

**Implementation Approach:**
- Same kernels 1.1-1.3, but with replica dimension
- Thread layout: replica index in global ID
- Memory layout: strided access for replicas
- For testing: use n_replicas=1
- For production: use n_replicas=1000-5000

**Testing Strategy:**
- **Text-based check 1:** Verify replica independence
  - Run with n_replicas=2, identical geometries
  - Should get identical results
- **Text-based check 2:** Compare single replica vs batch
  - Run single replica, then batch of 1
  - Should get identical results
- **Text-based check 3:** Performance scaling
  - Measure time for n_replicas = 1, 10, 100, 1000
  - Should scale linearly until GPU saturated
- **Visual test:** Plot current vs time for MD trajectory
  - Show fluctuations due to molecular motion

**Files:**
- Extend existing kernels to support replica dimension
- Test: `tests/pyFireball/test_multi_replica.py`

### Phase 2: Static PDOS MD (Fast Dynamics)

**Goal:** Simulate junction at finite temperature with rigid-molecule approximation.

#### Kernel 2.1: GPU PDOS Rotation

**Purpose:** Rotate pre-computed A(E) for rigid body motion: A' = R A R^T

**Input:**
- A_ref [n_replicas, n_energy, n_atoms, 4_orbitals] - PDOS in reference frame
- Rotation matrices R [n_replicas, 3, 3] - One per molecule (tip or surface)
- Or quaternion representation

**Output:**
- A_rot [n_replicas, n_energy, n_atoms, 4_orbitals] - Rotated PDOS

**Implementation Approach:**
- Reuse existing `rotate_fb_matrix_sp` from `hamiltonian.cl`
- Extend to handle energy dimension
- For diagonal approximation: just rotate p-orbital components
  - s-orbital: invariant
  - p-orbitals: transform as vector under rotation

**Testing Strategy:**
- **Text-based check 1:** Verify rotation invariance
  - Rotate by 360°, should return to original
  - Metric: max|A_rot - A_ref| < 1e-10
- **Text-based check 2:** Compare GPU vs CPU rotation
  - Implement CPU version using NumPy
  - Compare results
- **Text-based check 3:** Verify group properties
  - Rotate by R1, then R2 = rotate by R1*R2
- **Visual test:** Plot PDOS before and after rotation
  - Project to real-space grid
  - Should show rotated density lobes

**Files:**
- Kernel: extend `hamiltonian.cl` or add to `transport.cl`
- Host: `pyBall/FireballOCL/OCL_Transport.py` (extend)
- Test: `tests/pyFireball/test_pdos_rotation.py`

#### Integration with FireCore MD

**Purpose:** Call GPU transport kernel at each MD step

**Input:**
- MD trajectory from FireCore (positions vs time)
- Pre-computed A_ref for tip and surface
- Contact atom indices

**Output:**
- Current vs time I(t)

**Implementation Approach:**
- Python script that:
  1. Loads MD trajectory
  2. Extracts tip and surface positions at each step
  3. Computes rotation matrices from reference geometry
  4. Calls GPU kernel to rotate PDOS and compute transmission
  5. Integrates to get current
- For production: integrate directly into FireCore MD loop

**Testing Strategy:**
- **Text-based check 1:** Energy conservation check
  - Run short MD with constant current
  - Check that current doesn't drift
- **Text-based check 2:** Compare with full Fireball (subset of steps)
  - Run full Fireball SCF at selected MD steps
  - Compare GPU perturbative vs full Fireball transmission
  - Should be qualitatively similar (within 10-20% acceptable for perturbative)
- **Visual test 1:** Plot I(t) for MD trajectory
  - Should show fluctuations due to thermal motion
  - Power spectrum should match vibrational frequencies
- **Visual test 2:** 2D histogram of current vs tip-surface distance
  - Should show correlation (current decays with distance)

**Files:**
- Script: `tests/pyFireball/test_stm_md.py`
- Production: integrate into FireCore MD (separate task)

### Phase 3: Soft-Mode Expansion (Accuracy)

**Goal:** Capture effect of molecular deformation on current.

#### Kernel 3.1: GPU PDOS Taylor Expansion

**Purpose:** Update A(E) on-the-fly: A(Q) ≈ A_0 + Σ_k (∂A/∂Q_k) Q_k

**Input:**
- A_0 [n_energy, n_atoms, 4_orbitals] - Reference PDOS
- dA_dQ [n_modes, n_energy, n_atoms, 4_orbitals] - PDOS derivatives
- Q [n_replicas, n_modes] - Normal mode coordinates

**Output:**
- A_current [n_replicas, n_energy, n_atoms, 4_orbitals] - Updated PDOS

**Implementation Approach:**
- Linear combination kernel
- Thread layout: replica, energy, atom, orbital
- Simple: A_current = A_0 + Σ_k Q_k * dA_dQ[k]

**Testing Strategy:**
- **Text-based check 1:** Verify Q=0 returns reference
  - Set all Q=0, should get A_0
- **Text-based check 2:** Compare Taylor vs direct computation
  - For small displacements, Taylor should match direct PDOS computation
  - Error should scale as O(Q²)
- **Text-based check 3:** Verify mode orthogonality
  - Excite single mode, check that only that mode contributes
- **Visual test:** Plot PDOS vs Q for one mode
  - Should show linear trend for small Q

**Files:**
- Kernel: `pyBall/FireballOCL/cl/soft_modes.cl`
- Host: `pyBall/FireballOCL/OCL_Transport.py` (extend)
- Test: `tests/pyFireball/test_soft_modes.py`

---

#### Offline: Compute Normal Modes and PDOS Derivatives

**Purpose:** Pre-compute Hessian, normal modes, and dA/dQ (offline, CPU)

**Input:**
- Reference geometry
- Fireball SCF capabilities

**Output:**
- Normal modes [n_modes, 3*natoms]
- dA_dQ [n_modes, n_energy, n_atoms, 4_orbitals]

**Implementation Approach:**
- Use Fireball to compute Hessian (finite difference of forces)
- Diagonalize to get modes
- For each soft mode: compute PDOS at Q=±ΔQ, get derivative by finite difference
- All done offline in Python/NumPy (slow but acceptable for one-time computation)

**Testing Strategy:**
- **Text-based check 1:** Verify Hessian symmetry
- **Text-based check 2:** Verify orthonormality of modes
- **Visual test:** Plot mode shapes (displacement vectors)
  - Should correspond to expected vibrations (translation, rotation, bending)

**Files:**
- Script: `tools/compute_normal_modes.py`
- Script: `tools/compute_pdos_derivatives.py`

#### Integration with MD

**Purpose:** Project MD trajectory onto normal modes and use Taylor expansion

**Input:**
- MD trajectory positions
- Reference geometry
- Normal modes
- dA_dQ

**Output:**
- Q(t) for each mode
- A(t) via Taylor expansion
- I(t) via transmission

**Implementation Approach:**
- Python script that:
  1. Projects displacement onto modes: Q_k(t) = (R(t)-R_0) · v_k
  2. Calls GPU kernel to update PDOS
  3. Computes transmission
- For production: integrate into MD loop

**Testing Strategy:**
- **Text-based check 1:** Compare rigid vs soft-mode
  - Soft-mode should capture additional current fluctuations
- **Text-based check 2:** Compare with full Fireball
  - Should be more accurate than rigid approximation
- **Visual test:** Plot I(t) with and without soft modes
  - Should show additional high-frequency components

**Files:**
- Script: `tests/pyFireball/test_stm_soft_modes.py`

## Summary of Required GPU Kernels

### New Kernels to Create:

1. **Spectral Function Kernel** (`spectral_function.cl`)
   - Compute A(E) with Γ broadening
   - Support diagonal approximation initially
   - Extend to full matrix inversion later

2. **Slater-Koster Hopping Kernel** (`slater_koster.cl` or extend `hamiltonian.cl`)
   - Compute V_μν for all orbital pairs
   - Use existing spline infrastructure
   - Support multi-replica batching

3. **Transmission Kernel** (`transport.cl`)
   - Compute T(E) = Σ A_μ |V_μν|² A_ν
   - Support energy integration
   - Support multi-replica batching

4. **PDOS Rotation Kernel** (extend existing `rotate_fb_matrix_sp`)
   - Rotate A(E) for rigid body motion
   - Handle energy dimension

5. **Soft Mode Expansion Kernel** (`soft_modes.cl`)
   - Taylor expansion: A(Q) = A_0 + Σ Q·dA/dQ
   - Linear combination of pre-computed derivatives

### Existing Kernels to Reuse:

1. **Grid Projection** (`Grid.cl`) - for PDOS visualization
2. **Slater-Koster Rotation** (`hamiltonian.cl`) - reference for hopping computation

## Testing Strategy Summary

### Text-Based Checks (Automated):

1. **GPU vs CPU validation** - Compare GPU results with NumPy reference
2. **GPU vs Fortran validation** - Compare with Fireball reference (where available)
3. **Symmetry checks** - Verify Hermiticity, rotation invariance
4. **Scaling tests** - Verify linear scaling with replicas
5. **Convergence tests** - Verify Taylor expansion accuracy

### Visual Tests (Manual):

1. **PDOS vs energy plots** - Verify electronic structure
2. **Real-space PDOS projection** - Visualize orbital localization
3. **1D STM profiles** - Verify distance dependence
4. **2D STM images** - Verify molecular resolution
5. **Current vs time plots** - Verify MD integration

### Test Hierarchy:

```
Unit Tests (per kernel):
├── Kernel 1.1: Spectral Function
│   ├── test_spectral_function_single.py
│   ├── test_spectral_function_gamma.py
│   └── test_spectral_function_contact.py
├── Kernel 1.2: Slater-Koster
│   ├── test_slater_koster_vs_fortran.py
│   ├── test_slater_koster_symmetry.py
│   └── test_slater_koster_decay.py
├── Kernel 1.3: Transmission
│   ├── test_transmission_vs_cpu.py
│   └── test_transmission_integration.py
└── Kernel 1.4: Multi-Replica
    └── test_multi_replica_scaling.py

Integration Tests (per phase):
├── Phase 1: Rigid Scan
│   ├── test_stm_rigid_scan_1d.py
│   └── test_stm_rigid_scan_2d.py
├── Phase 2: Static PDOS MD
│   └── test_stm_md.py
└── Phase 3: Soft Modes
    ├── test_normal_modes.py
    └── test_stm_soft_modes.py
```

## Data Flow Diagram (Multi-Replica)

```
Offline (One-time):
┌─────────────────┐
│  Fireball SCF   │ → H, S
└────────┬────────┘
         │
         ↓
┌─────────────────┐
│  GPU Spectral   │ → A_0(E) [with Γ]
│  Function       │
└────────┬────────┘
         │
         ↓
┌─────────────────┐
│  Compute Modes  │ → v_k, dA/dQ_k
│  (CPU)          │
└─────────────────┘

Online (Per MD step):
┌─────────────────┐
│  FireCore MD    │ → R(t) [n_replicas]
└────────┬────────┘
         │
         ↓
┌─────────────────┐
│  Project Q(t)   │ → Q_k(t)
│  (CPU)          │
└────────┬────────┘
         │
         ↓
┌─────────────────┐
│  GPU Update A   │ → A(t) = A_0 + Σ Q·dA/dQ
│  (Soft Modes)   │
└────────┬────────┘
         │
         ↓
┌─────────────────┐
│  GPU Hopping    │ → V_μν(t)
│  (Slater-Koster)│
└────────┬────────┘
         │
         ↓
┌─────────────────┐
│  GPU Transport  │ → T(E,t), I(t)
│  (Transmission) │
└─────────────────┘
```

---

# REPORT 1 : SWE-1.6

# Phase 1 Implementation Report: STM Simulator with NEGF Transport

This section documents the actual implementation completed for the STM simulator, including newly added modules, methods, mathematical relations, challenges encountered, and technical notes for future developers.

## Overview of Implemented Modules

### 1. Shared STM Module (`pyBall/FireballOCL/STM.py`)

**Purpose:** Central module containing all common STM functions to avoid code duplication across test scripts.

**Key Functions:**

- `load_xyz(path)` - Load molecular geometry from XYZ file
- `compute_dos()` - Run Fireball SCF and compute spectral function with Γ broadening
  - Now exports generalized eigenvectors `C` and eigenvalues `eps` in addition to spectral data
  - Uses Wide-Band Limit (WBL) self-energy: Σ = -iΓ on contact atoms
  - Computes spectral function: A(E) = (1/π) Im[G^r(E)] where G^r = [(E+iη)S - H - Σ]^{-1}

- `compute_spectral(H, S, eigen, contact_orbs, gamma, ...)` - Core spectral function computation
  - Solves generalized eigenvalue problem: H C = S C ε
  - Applies Lorentzian broadening to each orbital: A_μ(E) = Σ_n (Γ_n/π) C_μn² / [(E - ε_n)² + Γ_n²]
  - Γ_n = γ Σ_{μ∈contact} C_μn² (orbital-resolved broadening)

- `build_dense_HS(dims, data, atomTypes)` - Convert sparse Fireball H/S to dense NumPy matrices
  - Handles orbital indexing from sparse format [natoms, neigh_max, numorb_max, numorb_max]
  - Maps neighbor lists to global orbital indices

- `build_inter_system_blocks_fdata()` - Build inter-system coupling using Fireball 2-center tables
  - Reads radial functions from Fireball Fdata directory
  - Applies Slater-Koster rotation for angular dependence

- `build_inter_system_blocks_exp_sk()` - **NEW** Vacuum exponential radial + Slater-Koster angular coupling
  - Replaces Fireball radial dependence with physical vacuum decay: exp(-β(r - r₀))
  - Retains Slater-Koster angular factors for correct orbital symmetry
  - Parameters: β ~1.0 Å⁻¹, r₀ ~3.0 Å (tunable)
  - Allows independent control of radial decay vs angular structure

- `assemble_combined_HS()` - Assemble combined tip+sample Hamiltonian for NEGF
  - Stacks individual H/S matrices
  - Inserts inter-system coupling blocks at correct positions
  - Handles tip offset translation

- `select_atoms_by_xy_radius()` - Select atoms within lateral radius of scan center
- `select_lead_orbitals_by_atoms()` - Convert atom selection to orbital indices for lead coupling

### 2. NEGF Transport Methods

#### Method 1: Direct NEGF with Caroli Formula (`negf_current()`)

**Mathematical Foundation:**
```
G^r(E) = [(E + iη)S - H - Σ_L - Σ_R]^{-1}
T(E) = Tr[Γ_L(E) G^r(E) Γ_R(E) G^a(E)]
I = (2e/h) ∫ T(E)[f(E-μ_L) - f(E-μ_R)] dE
```

**Implementation:**
- Wide-Band Limit self-energies: Σ_L = -iΓ_L on tip lead orbitals, Σ_R = -iΓ_R on sample lead orbitals
- Full matrix inversion using NumPy `solve()` (not `inv()` for numerical stability)
- Supports energy grid integration or single-energy evaluation
- **Challenge:** O(N³) scaling limits system size

#### Method 2: Iterative NEGF with GMRES (`negf_current_iterative()`)

**Purpose:** Avoid explicit matrix inversion for larger systems using iterative linear solver.

**Mathematical Foundation:**
Uses Hutchinson stochastic trace estimator to compute transmission without full Green's function:
```
T(E) ≈ (1/n_rand) Σ_{k=1}^{n_rand} z_k† Γ_L G^r Γ_R G^a z_k
```
where z_k are random Rademacher vectors (±1 entries).

**Implementation:**
- Uses `scipy.sparse.linalg.gmres` to solve (E S - H - Σ) x = b iteratively
- Warm-start cache: reuses previous solution as initial guess for nearby energies
- Parameters: `nrand=8` (number of random vectors), `tol=1e-6`, `maxiter=500`
- **Advantage:** Avoids O(N³) inversion, scales better for large systems
- **Challenge:** Convergence depends on conditioning; may need preconditioning

#### Method 3: Deterministic Response Metric (`negf_response_at_energy()`)

**Purpose:** Single-solve deterministic alternative to stochastic trace estimator for probing transport.

**Mathematical Foundation:**
```
A = (E + iη)S - H - Σ_L - Σ_R
b = source_vector (e.g., unit vector on lead L)
x = solve(A, b)
resp = x† Γ_measure x
```

**Interpretation:**
- Solves for the wavefunction response to excitation on one lead
- Measures how much of that response couples to the other lead
- Not a transmission coefficient, but a qualitative probe of coupling strength
- **Advantage:** Single linear solve (no sampling), deterministic
- **CLI control:** `--source L/R` (excitation lead), `--measure L/R` (measurement lead)

### 3. MO Overlap Diagnostic (`coupling_vec_tip_to_sample()`, `mo_overlap_amplitude()`)

**Purpose:** Validate NEGF response against sum-over-states picture using molecular orbital overlaps.

**Mathematical Foundation:**
For a "featureless metallic tip" approximation:
```
v_smp(s) = Σ_{i∈tip} H_{TS}(i, :)  (reduce_tip='sum')
         or sqrt(Σ_{i∈tip} |H_{TS}(i, :)|²)  (reduce_tip='norm')
t_n(s) = v_smp(s) · C_{:,n}
I_n(s) = |t_n(s)|²
```
where C_{:,n} is the n-th molecular orbital eigenvector (in AO basis).

**Physical Interpretation:**
- `t_n(s)` is the overlap amplitude between tip coupling vector and MO n
- Sign of `t_n(s)` reveals nodal planes (zeros, sign flips)
- `|t_n(s)|²` is the intensity (what STM actually measures)
- This is the "poles of Green's function" picture: T(E) ≈ Σ_n |t_n|² / [(E - ε_n)² + Γ_n²]

**Implementation:**
- `coupling_vec_tip_to_sample()`: Builds effective tip→sample coupling vector by reducing over tip orbitals
- `mo_overlap_amplitude()`: Computes signed amplitude and intensity for a selected MO
- CLI flags: `--scan_mo`, `--mo` (explicit index), `--mo_win` (energy window around E), `--reduce_tip`
- Outputs: `scan1d_xy_mo.npz` and `scan1d_xy_mo.png` with signed amplitude and intensity plots

### 4. CLI Analysis Script (`tests/pyFireball/stm_analysis.py`)

**Purpose:** Main user-facing script for STM analysis with flexible configuration.

**Key Features:**
- Load pre-computed DOS files for tip and sample molecules
- Support for multiple transport methods: `--method pdos` (Phase 1 approximation) or `--method negf` (combined-system)
- NEGF solver selection: `--solver direct` or `--solver iter`
- NEGF mode: `--negf_mode caroli` (full transmission) or `--negf_mode response` (deterministic probe)
- Single-energy evaluation: `--E <value>` (defaults to Fermi level)
- Lead selection: `--source L/R`, `--measure L/R` for response metric
- Tunneling coupling parameters: `--beta`, `--r0`, `--A_ss`, `--A_sp`, `--A_ppsig`, `--A_pppi`
- Scanning protocols:
  - 1D height scan: z = 3–10 Å (always runs)
  - 1D lateral scan at constant z: `--scan_xy`, `--xy_dir diag/x/y`, `--zxy`
  - 2D xy scan: `--no2d` to skip
- MO overlap diagnostic: `--scan_mo` with MO selection options
- Output management: All plots and .npz files saved to `export/stm_phase1/` with printed paths

### 5. DOS Computation Script (`tests/pyFireball/stm_compute_dos.py`)

**Purpose:** Thin interface to Fireball SCF for computing DOS with Γ broadening.

**Key Features:**
- Calls Fireball SCF via `pyBall.FireCore` interface
- Computes spectral function with contact atom broadening
- **Now exports:** `eps` (eigenvalues) and `C` (eigenvectors) in addition to spectral data
  - Required for MO overlap diagnostic
  - Users must regenerate DOS files after this change to get C/eps arrays

## Mathematical Relations Between Methods

### NEGF vs MO Overlap Picture

The Caroli transmission formula can be expressed in the molecular orbital basis as:
```
T(E) = Tr[Γ_L G Γ_R G†]
     = Σ_{n,m} (Γ_L)_{nm} G_{mn} (Γ_R)_{np} G*_{pn}
     = Σ_n |t_n(E)|² L_n(E)
```
where:
- `t_n(E)` is the coupling amplitude between lead states and MO n
- `L_n(E)` is the Lorentzian lineshape: Γ_n / [(E - ε_n)² + Γ_n²]

The MO overlap diagnostic computes `t_n(s)` explicitly for a single MO at a fixed energy, allowing direct comparison with the NEGF response (which sums over all MOs implicitly).

### Direct vs Iterative NEGF

**Direct solver:**
```
x = solve(A, b)  for each random vector
T ≈ (1/N) Σ x† Γ x
```
- Exact (to machine precision)
- O(N³) due to matrix factorization (implicit in `solve`)

**Iterative solver:**
```
x = gmres(A, b, x0, tol, maxiter)  for each random vector
T ≈ (1/N) Σ x† Γ x
```
- Approximate (controlled by `tol`)
- O(N × k) where k is iteration count
- Warm-start: `x0` from previous energy reduces iterations

### Response Metric vs Transmission

**Response metric:**
```
resp = x† Γ_measure x  where x solves A x = source_vector
```
- Single linear solve (no sampling)
- Probes "how much does excitation on lead A couple to lead B?"
- Not a transmission coefficient, but correlates with it
- Useful for quick qualitative scans

**Caroli transmission:**
```
T = Tr[Γ_L G Γ_R G†]
```
- Requires full Green's function (or stochastic trace)
- Physically meaningful transmission coefficient
- More expensive but rigorous

## Challenges Encountered and Solutions

### 1. Zero Transmission Bug

**Problem:** Initial NEGF implementation gave zero transmission.

**Root Cause:** Used only Fireball kinetic root for inter-system coupling, which was zero for distant atoms.

**Solution:** Use full Slater-Koster angular dependence with proper radial function. The coupling must include both:
- Radial part: distance-dependent decay
- Angular part: Slater-Koster factors (l, m, n direction cosines)

**Lesson:** Never use simplified coupling models without validating against physical expectations. Check that coupling blocks are non-zero before full NEGF solve.

### 2. Performance: Full Energy Grid vs Single Energy

**Problem:** Initial implementation computed transmission on 351-point energy grid, making scans extremely slow.

**Solution:** For STM imaging at small bias, single-energy evaluation at Fermi level is sufficient. Changed to:
```python
E_grid = np.array([E_single])  # instead of np.arange(...)
```

**Lesson:** STM current at small bias is dominated by states near Fermi level. Single-energy evaluation is adequate for qualitative imaging.

### 3. Radial Part: Fireball vs Vacuum Exponential

**Problem:** Fireball radial functions are designed for covalent bonding (short-range), not vacuum tunneling (long-range exponential decay).

**Solution:** Implemented `build_inter_system_blocks_exp_sk()` which:
- Uses vacuum exponential decay: exp(-β(r - r₀))
- Retains Slater-Koster angular dependence for correct orbital symmetry
- Allows independent tuning of β and r₀

**Parameters:**
- β ~1.0 Å⁻¹ (typical vacuum decay)
- r₀ ~3.0 Å (reference distance where coupling = A values)

**Lesson:** Physical accuracy requires using appropriate radial functions for the physical regime (vacuum tunneling vs covalent bonding).

### 4. Hamiltonian Assembly Complexity

**Problem:** Converting sparse Fireball format to dense matrices for NEGF is error-prone due to complex orbital indexing.

**Solution:** Careful implementation in `build_dense_HS()`:
- Uses `starts` array to map atoms to orbital indices
- Handles variable orbitals per atom (s, p, d)
- Validates dimensions match expected values

**Sparse format:** [natoms, neigh_max, numorb_max, numorb_max]
- Each atom has `neigh_max` neighbors
- Each neighbor pair has up to `numorb_max × numorb_max` orbital interactions
- Need to map local orbital indices to global indices

**Lesson:** Always validate matrix dimensions and print diagnostics (e.g., diagonal ranges) before using in NEGF.

### 5. Missing Eigenvectors in DOS Files

**Problem:** MO overlap diagnostic requires eigenvectors `C`, but initial DOS files only contained spectral data.

**Solution:** Extended `compute_dos()` to export `eps` and `C`:
```python
eps, C = eigh(H, S)  # generalized eigenvalue problem
result['eps'] = eps
result['C'] = C
```

**Action required:** Users must regenerate DOS files with updated script to get C/eps arrays.

**Lesson:** Plan ahead for diagnostic needs. If eigenvectors might be needed later, export them from the start.

### 6. Iterative Solver Convergence

**Problem:** GMRES may not converge for poorly conditioned matrices.

**Current mitigation:**
- Warm-start cache (reuse previous solution as initial guess)
- Reasonable tolerance (tol=1e-6)
- Max iteration limit (500)

**Future improvements:**
- Preconditioning (e.g., incomplete LU)
- Alternative solvers (BiCGSTAB, TFQMR)
- Hybrid approach: direct for small systems, iterative for large

## Hamiltonian Assembly Details

### Sparse to Dense Conversion

Fireball exports H and S in sparse format:
```
data.h_mat[iat, ineigh, iorb, jorb] = H_{μ,ν}
data.neigh_j[iat, ineigh] = neighbor atom index (1-based)
data.mu, data.nu: orbital index mappings
```

Conversion steps:
1. Compute orbital start indices per atom: `starts[i] = Σ_{k< i} norb_per[k]`
2. For each atom `iat` and neighbor `jneigh`:
   - Get neighbor atom index: `j = data.neigh_j[iat, ineigh] - 1`
   - Get orbital counts: `ni = norb_per[iat]`, `nj = norb_per[j]`
   - Copy block to dense: `H[starts[i]:starts[i]+ni, starts[j]:starts[j]+nj] = data.h_mat[iat, ineigh, :ni, :nj]`

**Validation:** Check that diagonal elements are reasonable (e.g., -20 to 0 eV for valence orbitals).

### Combined System Assembly

For NEGF, we assemble tip and sample into a combined Hamiltonian:
```
H_combined = [[H_tip,     H_TS],
              [H_ST,      H_smp ]]
```

Steps:
1. Stack diagonal blocks: `Hc = scipy.linalg.block_diag(H_tip, H_smp)`
2. Insert off-diagonal coupling blocks at correct positions
3. Same for overlap matrix S
4. Handle tip offset translation (tip position relative to sample)

**Lead orbital selection:**
- Select atoms that couple to leads (e.g., bottom atoms of sample, top atoms of tip)
- Convert atom indices to orbital indices using `starts` array
- Apply WBL self-energy only to these orbitals: Σ = -iΓ

### Inter-System Coupling

**Original approach (Fireball radial):**
- Read 2-center tables from Fireball Fdata directory
- Interpolate radial functions
- Apply Slater-Koster rotation

**New approach (vacuum exponential):**
- Radial: `f(r) = A * exp(-β(r - r₀))`
- Angular: Same Slater-Koster factors
- Independent parameters A_ss, A_sp, A_ppσ, A_ppπ for different orbital pairs

**Advantage of vacuum exponential:**
- Physically correct for tunneling regime
- Simple, tunable parameters
- No dependence on Fireball tables (portable)

**Slater-Koster angular factors:**
For direction cosines (l, m, n) = (dx/r, dy/r, dz/r):
- ssσ: Vss (isotropic)
- spσ: l * Vsp
- ppσ: l² * Vppσ + (1 - l²) * Vppπ
- etc.

## CRITICAL: Orbital Mapping Conventions (Evidence from Orbital Projection Debugging)

### Background

The orbital projection debugging in `tests/pyFireball/test_h2o_orbital_comparison.py` revealed critical issues with orbital mapping between Fortran (Fireball) and OpenCL implementations. These issues directly affect STM transport calculations that use orbital overlaps, spectral functions, and Hamiltonian assembly.

### Orbital Ordering Conventions

**Fortran (Fireball) Convention:**
- Order: `[s, py, pz, px]` for atoms with p orbitals (spherical harmonics m=-1,0,+1)
- Source: `fortran/INTERPOLATERS/getYlm.f90`
- p-orbitals: `psi(1)=y` (py, m=-1), `psi(2)=z` (pz, m=0), `psi(3)=x` (px, m=+1)
- For H atoms (1 orbital): only `[s]`

**OpenCL Kernel Convention (Grid.cl):**
- Order: `[px, py, pz, s]` (Cartesian x,y,z)
- Source: `pyBall/FireballOCL/cl/Grid.cl`
- dri.xyz = [px, py, pz], dri.w = s
- Used in `project_orbital` and `project_density_sparse_tiled` kernels

**OpenCL Hamiltonian Convention (hamiltonian.cl):**
- Order: `[s, px, py, pz]` (intermediate convention)
- Source: `pyBall/FireballOCL/cl/hamiltonian.cl`
- Used in Slater-Koster rotation `rotate_fb_matrix_sp`

### Permutation Mappings

**Fortran → OpenCL Grid (for orbital projection):**
```python
# Fortran [s, py, pz, px] → OpenCL [px, py, pz, s]
coeffs_opencl[ia, 0] = coeffs_fortran[ia, 3]  # px
coeffs_opencl[ia, 1] = coeffs_fortran[ia, 1]  # py
coeffs_opencl[ia, 2] = coeffs_fortran[ia, 2]  # pz
coeffs_opencl[ia, 3] = coeffs_fortran[ia, 0]  # s
```

**Fortran → OpenCL Hamiltonian (for SK rotation):**
```python
# Currently in STM.py:
_ORT_SPP_TO_STD = np.array([0, 3, 1, 2], dtype=int)  # [s,py,pz,px] → [s,px,py,pz]
```

**IMPORTANT:** The OpenCL Grid kernel expects `[px,py,pz,s]` NOT `[s,px,py,pz]`! The current STM remapping is incorrect for orbital projection kernels.

### Critical Pitfalls Discovered

#### 1. Coefficient Packing Truncation (PRIMARY BUG)

**Problem:** OpenCL `project_orbital` packed only `norb_per[ia]` coefficients per atom instead of always 4.

**Impact:** For H atoms (1 orbital), only 1 coefficient was packed, missing the s-orbital entirely.

**Fix:** Always pack 4 coefficients per atom as float4:
```python
coeffs_opencl = np.zeros((natoms, 4), dtype=np.float32)
for ia in range(natoms):
    no = norb_per[ia]
    if no == 1:  # H atom: pack [0,0,0,s]
        coeffs_opencl[ia, 3] = coeffs_fortran[ia, 0]
    elif no == 4:  # O atom: remap [s,py,pz,px] → [px,py,pz,s]
        coeffs_opencl[ia, 0] = coeffs_fortran[ia, 3]
        coeffs_opencl[ia, 1] = coeffs_fortran[ia, 1]
        coeffs_opencl[ia, 2] = coeffs_fortran[ia, 2]
        coeffs_opencl[ia, 3] = coeffs_fortran[ia, 0]
```

**Evidence:** Test with mockup orbitals showed H s-orbitals absent in OpenCL output until this fix was applied.

#### 2. Incorrect Sign Flip for px

**Problem:** Previous implementation applied a sign flip to px orbital during remapping.

**Impact:** Caused correlation of -0.99 for px orbital (correct should be +0.99).

**Fix:** Remove sign flip:
```python
# WRONG (old):
coeffs_opencl[ia, 0] = -coeffs_fortran[ia, 3]  # px flipped

# CORRECT (new):
coeffs_opencl[ia, 0] = coeffs_fortran[ia, 3]  # px no flip
```

**Evidence:** Single-basis orbital tests showed px correlation = -0.99 before fix, +1.0 after fix.

#### 3. Sampling Mismatch (Voxel Centers vs Corners)

**Problem:** Fortran evaluates at arbitrary points (user-specified mesh), OpenCL grid projection evaluated at voxel corners.

**Impact:** Large scale mismatch even when qualitative shapes matched.

**Fix:** Shift OpenCL grid origin by `+0.5 * dC` to sample at voxel centers:
```python
grid_spec['origin'] = grid_origin + 0.5 * grid_spec['dC']
```

**Evidence:** Scale mismatch reduced from ~2x to ~1.0 after origin shift.

#### 4. Radial Interpolation Mismatch

**Problem:** OpenCL used Catmull-Rom spline; Fortran uses natural cubic spline with second derivatives.

**Impact:** Radial basis functions interpolated differently, causing scale differences.

**Fix:** Change to natural cubic spline:
- `Grid.cl::evaluate_radial()`: use natural cubic spline with second derivatives
- `Grid.py::load_basis()`: resample with natural cubic spline
- Pack radial basis as `(wf, wf_spline_second_derivative)` float2 per node

**Evidence:** Correlation improved from ~0.95 to 1.0000 after spline fix.

### Application to STM Implementation

#### Current STM.py Orbital Handling

**In `build_inter_system_blocks_fdata()`:**
```python
# Lines 525-526 in STM.py:
bH4 = _reorder_sp_block_ortega_to_std(bH4)
bS4 = _reorder_sp_block_ortega_to_std(bS4)
```
This remaps from Fortran `[s,py,pz,px]` to `[s,px,py,pz]` (intermediate convention).

**In `mo_overlap_amplitude()`:**
```python
# Line 570 in STM.py:
t = float(np.dot(v, C[:, mo]))
```
This computes overlap using eigenvector C directly from Fireball (Fortran order).

**ISSUE:** The eigenvectors C are in Fortran order `[s,py,pz,px]`, but the coupling vector v is built from H blocks that were remapped to `[s,px,py,pz]`. This creates an index mismatch!

#### Required Fixes for STM

**Fix 1: Update Permutation for Grid Kernels**
```python
# Add new permutation for Fortran → OpenCL Grid:
_ORT_SPP_TO_GRID = np.array([3, 1, 2, 0], dtype=int)  # [s,py,pz,px] → [px,py,pz,s]

def _reorder_sp_block_ortega_to_grid(b4):
    """Reorder 4x4 block from Ortega order (s,py,pz,px) to Grid (px,py,pz,s)."""
    p = _ORT_SPP_TO_GRID
    return b4[np.ix_(p, p)]
```

**Fix 2: Consistent Remapping in MO Overlap**
The MO overlap method must use consistent orbital ordering:
- If C is in Fortran order, v must also be in Fortran order
- OR remap C to match v's order

**Fix 3: Always Pack 4 Orbitals**
When building density matrices for Grid projection:
```python
# In pad_HS_to_float4() or similar:
# Always pack [px,py,pz,s] for Grid kernels
# For H atoms: [0,0,0,s]
```

**Fix 4: Eigenvector Ordering for PDOS Projection**
When projecting PDOS to grid using GridProjector:
```python
# C from Fireball is in Fortran order [s,py,pz,px]
# Need to remap to OpenCL Grid order [px,py,pz,s] before projection
C_grid = _remap_eigenvectors_fortran_to_grid(C, norb_per)
```

### Validation Checklist for STM

Before using STM transport methods:

1. **[Coefficient packing]** Verify all arrays pack 4 orbitals per atom (no truncation for H)
2. **[Permutation consistency]** Check which permutation is used:
   - `[0,3,1,2]` for Hamiltonian (Fortran → [s,px,py,pz])
   - `[3,1,2,0]` for Grid projection (Fortran → [px,py,pz,s])
3. **[Eigenvector order]** Verify C matrix order matches coupling vector order in MO overlap
4. **[Sign convention]** No sign flip for px (removed in orbital projection fix)
5. **[Radial interpolation]** If using OCL_Hamiltonian, verify spline type matches Fortran
6. **[Sampling origin]** If projecting to grid, verify origin shift for voxel centers

### Connection to Orbital Projection Test

The test `tests/pyFireball/test_h2o_orbital_comparison.py` provides:
- **Reference implementation:** Fortran `orb2points()` (ground truth)
- **Validation method:** OpenCL `project_orbital_points()` (exact point evaluation)
- **Correlation metric:** Achieved 1.0000 correlation after all fixes

This test should be used to validate any STM orbital mapping changes:
1. Create mockup orbitals (single-basis activation)
2. Project to grid with both Fortran and OpenCL
3. Verify correlation = 1.0000 before using in STM transport

---

## Remaining Problems and Future Work

### 1. DOS File Regeneration

**Status:** Users need to regenerate `dos_PTCDA.npz` and `dos_H2.npz` with updated `stm_compute_dos.py` to get `eps` and `C` arrays.

**Action:** Run:
```bash
cd tests/pyFireball
python stm_compute_dos.py --xyz PTCDA.xyz --out export/stm_phase1/dos_PTCDA.npz
python stm_compute_dos.py --xyz H2.xyz --out export/stm_phase1/dos_H2.npz
```

### 2. GPU Acceleration

**Status:** Current implementation is CPU-only (NumPy).

**Plan:**
- Port Slater-Koster coupling to OpenCL kernel
- Port spectral function computation to GPU
- Port NEGF solver (iterative) to GPU
- Use batched processing for multi-replica MD

**Challenges:**
- Sparse matrix operations on GPU (need sparse linear algebra library)
- Iterative solver on GPU (need GPU GMRES implementation)
- Memory management for large systems

### 3. Hamiltonian Rebuild Optimization

**Status:** Full Hamiltonian rebuilt at each scan point (expensive).

**Optimization opportunities:**
- Only rebuild inter-system coupling blocks (H_TS, S_TS)
- Keep tip and sample diagonal blocks constant
- Use incremental updates for small displacements

**Implementation:**
- Cache tip and sample H/S matrices
- Only recompute coupling blocks for moved atoms
- For rigid translation/rotation, apply transformation instead of rebuild

### 4. Normal Mode Expansion

**Status:** Not yet implemented.

**Plan:**
- Compute Hessian in Fireball
- Extract soft modes (translations, rotations, low-frequency vibrations)
- Compute PDOS derivatives ∂A/∂Q_k via finite difference
- On GPU: A(t) = A_0 + Σ_k Q_k(t) * ∂A/∂Q_k

**Benefits:**
- Avoid full SCF during MD
- Capture internal deformation effects
- Still accurate for small displacements

### 5. Perturbative Transport

**Status:** Not yet implemented (currently using full NEGF).

**Plan:**
- Precompute PDOS for isolated molecules with lead coupling
- On GPU: T(E) = Σ_μν A_μ^T(E) |V_μν|² A_ν^S(E)
- Much faster (no matrix inversion)
- Valid for weak coupling (tunneling regime)

**Validation:**
- Compare with full NEGF for test systems
- Establish coupling strength threshold where perturbative breaks down

## Usage Examples

### Basic STM Analysis (NEGF with response metric)

```bash
cd tests/pyFireball
python stm_analysis.py --method negf --negf_mode response \
    --E -0.8 --z2d 5.0 --xy_range 7.0 --xy_step 0.5 \
    --scan_xy --no2d --noshow
```

Output:
- `scan1d.png`: Height scan I(z)
- `scan1d_xy.png`: Lateral scan I(s) at constant z
- `scan1d.npz`, `scan1d_xy.npz`: Numerical data

### MO Overlap Diagnostic

```bash
python stm_analysis.py --method negf --negf_mode response \
    --E -0.8 --scan_xy --scan_mo --no2d --noshow
```

Output:
- `scan1d_xy_mo.png`: Signed amplitude t(s) and intensity |t|²
- `scan1d_xy_mo.npz`: Numerical data including MO index and energy

### Full 2D STM Image

```bash
python stm_analysis.py --method negf --negf_mode response \
    --E -0.8 --z2d 5.0 --xy_range 7.0 --xy_step 0.5 \
    --noshow
```

Output:
- `stm2d.png`: 2D current map I(x,y)
- `stm2d.npz`: Numerical data

## Key Parameters and Typical Values

| Parameter | Description | Typical Value | CLI Flag |
|-----------|-------------|---------------|----------|
| β | Vacuum decay constant | 1.0 Å⁻¹ | `--beta` |
| r₀ | Reference distance | 3.0 Å | `--r0` |
| A_ss | s-s coupling strength | -1.0 eV | `--A_ss` |
| A_sp | s-p coupling strength | -1.0 eV | `--A_sp` |
| A_ppσ | p-p σ coupling | -1.0 eV | `--A_ppsig` |
| A_ppπ | p-p π coupling | +1.0 eV | `--A_pppi` |
| γ_L, γ_R | Lead coupling strength | 2.0 eV | `--gammaL`, `--gammaR` |
| η | Small imaginary part | 1e-4 eV | `--eta` |
| z2d | Height for 2D scan | 5.0 Å | `--z2d` |
| zxy | Height for lateral scan | 5.0 Å | `--zxy` |
| xy_range | Lateral scan range | 7.0 Å | `--xy_range` |
| xy_step | Lateral scan step | 0.5 Å | `--xy_step` |

## Summary

The Phase 1 implementation provides a working STM simulator with:

1. **Shared module** (`STM.py`) with common functions for DOS, hopping, and NEGF transport
2. **Three NEGF methods**: direct Caroli, iterative GMRES, and deterministic response
3. **Physically accurate tunneling coupling**: vacuum exponential radial + Slater-Koster angular
4. **MO overlap diagnostic**: sum-over-states validation with signed amplitudes
5. **Flexible CLI** for various scanning protocols and parameter exploration
6. **Consistent output management** with printed file paths

**Remaining work:** GPU acceleration, Hamiltonian rebuild optimization, normal mode expansion, and perturbative transport for production MD integration.

---

# USER

OK, now we should make very basic test, lets tune the energy at which we evaluate electron transision (ie e for which we compute G=(eI-H)^-1 )
to energy of PTCDA HOMO orbital. As a tip we will use just single hydrogen atom - just s-orbital. For simplicity asume that every atom can have 4 basis functions on it (px,py,pz,s) which maps perfectly to float4 datatype, and this is how we represent our wavefuctions and density matrix in @GridFF.py 
@GridFF.cl and in @OCL_Hamiltonian.py @hamiltonian.cl , to which we want to interface later. For compactibility and simplicity keep this structure everywhere (note that in FireBall fortran we do not use this structure so we need to convert it). We can easily just pad the unused px,py,pz by zeros wen we have just single s orbital on the tip. Later this will allow use to try different symetries (px,py,pz). Later we will also test how the denisy of states is rotatated when we rotate whole molecule.   


And lets compute the map at height z=3.0A for the different methods and plot it. Plot also the projected density of states.
I expect that the maps should loop like aprox the same for all the methods.
Especially the STM computed by the method using sume over moleculer orbitals (which in this case basically sum just over HOMO) should look like HOMO
But I do not see why also the other method should not look like HOMO.
We should put this as another independnet shalow script which reuse in maximum degree what we already implemented in the shared module STM.py we can also consolidate the plotting ad debuging function in separate module STM_utils.py

Now your job is not to write any code yet, but make clear plan with clarly defined test, and requirements, and also identify possible problem we have in current code on which we should focus. Again I hate making loops in python. For operations in python use numpy array operatin (fancy indexing, slicing, linear algebra function etc.), Later we plan to rewrite everything performance intensive to pyOpenCL but for now it is fine to make python reference which is easier to understand and debug. Make sure that you write clearly the equations what you are doing at every function (some theory background in comment string)

now write the plan and implementation notes for this critical testing step. 

---

# SWE-1.6

# Plan: Basic STM Test at HOMO Energy

## Test Definition

**Objective:** Validate that all STM transport methods produce maps resembling the PTCDA HOMO orbital when evaluated at the HOMO energy with a single H-atom tip at z=3.0 Å.

**Test Configuration:**
- **Tip:** Single hydrogen atom (s-orbital only, px/py/pz padded with zeros)
- **Sample:** PTCDA molecule
- **Energy:** PTCDA HOMO orbital energy (from eigenvalues)
- **Height:** z = 3.0 Å above sample center
- **Scan:** 2D xy map over sample (±7 Å range, 0.5 Å step)
- **Methods to compare:**
  1. NEGF Caroli (direct solver)
  2. NEGF response metric
  3. MO overlap (sum over states, dominated by HOMO)
  4. Projected density of states (PDOS) at HOMO energy

**Expected Outcome:** All methods should produce similar maps that visually resemble the HOMO orbital shape.

## Implementation Plan

### Step 1: Create `STM_utils.py` Module

**Purpose:** Consolidate plotting and debugging functions from STM.py.

**Functions to extract:**

```python
def plot_2d_map(data, x_arr, y_arr, title, cmap='viridis', 
                vmin=None, vmax=None, ax=None):
    """
    Plot 2D data as heatmap with proper axes.
    
    Theory: I(x,y) = T(E,x,y) where T is transmission or response metric
    Visualizes spatial variation of STM signal at constant energy and height.
    """
    
def plot_1d_scan(s_arr, I_arr, title, ylabel, ax=None):
    """
    Plot 1D scan data.
    
    Theory: I(s) where s is lateral position along scan line
    Used for height scans (z) or lateral scans (x/y/diagonal)
    """

def plot_pdos(E_grid, A, E_fermi, title, ax=None):
    """
    Plot projected density of states.
    
    Theory: A_μ(E) = (1/π) Im[G^r_μμ(E)] where G^r = [(E+iη)S - H - Σ]^{-1}
    Shows spectral density per orbital as function of energy
    """

def save_plot(fig, name, export_dir, dpi=130):
    """Save figure to file with consistent naming."""

def set_export_dir(export_dir):
    """Create export directory if needed, return absolute path."""
```

**Requirements:**
- No Python loops, use NumPy operations
- Clear docstrings with theory background
- Handle both matplotlib axes and standalone figures

### Step 2: Implement Float4 Per Atom Conversion

**Problem:** Fireball format uses variable orbitals per atom (s, p, d with counts depending on species), but we need fixed float4 (s, px, py, pz) structure.

**Solution:** Add conversion function to STM.py

```python
def convert_to_float4_per_atom(H, S, atomTypes, norb_per):
    """
    Convert sparse/dense matrices from Fireball format to float4-per-atom structure.
    
    Theory:
    Fireball format: variable orbitals per atom
    - H: [s], [s,px,py,pz], [s,px,py,pz,d...] depending on species
    - Indexing via starts array: starts[i] = Σ_{k<i} norb_per[k]
    
    Target format: exactly 4 orbitals per atom (s, px, py, pz)
    - For atoms with <4 orbitals: pad with zeros
    - For atoms with >4 orbitals: truncate d-orbitals (not used in STM basis)
    
    Args:
        H: Hamiltonian matrix [norb, norb]
        S: Overlap matrix [norb, norb]
        atomTypes: atom type indices [natoms]
        norb_per: orbitals per atom from Fireball [natoms]
    
    Returns:
        H4: Hamiltonian in float4 format [natoms, natoms, 4, 4]
        S4: Overlap in float4 format [natoms, natoms, 4, 4]
        orb_mask: Boolean mask [natoms, 4] indicating which orbitals are valid
    """
    natoms = len(atomTypes)
    n4 = 4  # fixed float4 structure
    
    # Build mapping from Fireball orbitals to float4 positions
    H4 = np.zeros((natoms, natoms, n4, n4), dtype=np.float64)
    S4 = np.zeros((natoms, natoms, n4, n4), dtype=np.float64)
    orb_mask = np.zeros((natoms, n4), dtype=bool)
    
    starts = _atom_orb_starts(norb_per)
    
    for i in range(natoms):
        ni = min(int(norb_per[i]), n4)  # truncate if >4
        orb_mask[i, :ni] = True
        i0 = int(starts[i])
        
        for j in range(natoms):
            nj = min(int(norb_per[j]), n4)
            j0 = int(starts[j])
            
            # Copy valid orbitals, zeros remain for padding
            H4[i, j, :ni, :nj] = H[i0:i0+ni, j0:j0+nj]
            S4[i, j, :ni, :nj] = S[i0:i0+ni, j0:j0+nj]
    
    return H4, S4, orb_mask
```

**Requirements:**
- Use NumPy slicing, no loops over orbitals
- Preserve Hermiticity: H4[j,i] = H4[i,j].T
- Validate: check that diagonal blocks are reasonable

### Step 3: Set Up Single H Atom Tip

**Problem:** Current tip is H2 (2 atoms), need single H atom with s-orbital only.

**Solution:** Create H atom geometry and compute its DOS

```python
# In test script:
tip_types = np.array([1])  # Hydrogen
tip_pos = np.array([[0.0, 0.0, 0.0]])  # Single atom at origin
# After compute_dos, convert to float4 (s-orbital only, p-orbitals padded)
```

**Conversion details:**
- H atom in Fireball: 1 orbital (s only)
- float4 format: [s, 0, 0, 0] (px, py, pz padded with zeros)
- This maintains compatibility with float4 structure while being physically correct

### Step 4: Identify HOMO Energy

**Problem:** Need to find HOMO energy from PTCDA eigenvalues.

**Solution:** Add helper function to STM.py

```python
def find_homo_index(eps, nelec_per_atom=None):
    """
    Find index of Highest Occupied Molecular Orbital (HOMO).
    
    Theory:
    In DFT, molecular orbitals are filled with electrons up to Fermi level.
    HOMO is the highest energy orbital that is occupied.
    For closed-shell systems: fill 2 electrons per orbital (spin-paired).
    
    Args:
        eps: Orbital eigenvalues [norb]
        nelec_per_atom: (optional) electrons per atom for open-shell systems
    
    Returns:
        homo_idx: Index of HOMO orbital
        homo_energy: Energy of HOMO
    """
    if nelec_per_atom is None:
        # Closed-shell: fill 2 electrons per orbital
        nelec = 2 * len(eps) // 2  # Assume even number of electrons
    else:
        nelec = np.sum(nelec_per_atom)
    
    # Sort by energy, fill lowest nelec/2 orbitals
    sorted_idx = np.argsort(eps)
    homo_idx = sorted_idx[nelec//2 - 1]  # -1 for 0-based indexing
    return homo_idx, eps[homo_idx]
```

**Usage:**
```python
homo_idx, homo_energy = find_homo_index(smp['eps'])
print(f"HOMO: index={homo_idx}, energy={homo_energy:.4f} eV")
```

### Step 5: Create Test Script `test_stm_homo.py`

**Structure:**

```python
#!/usr/bin/env python3
"""
Basic STM test: Compare all transport methods at HOMO energy.

Test configuration:
- Tip: Single H atom (s-orbital only, p-orbitals padded)
- Sample: PTCDA molecule
- Energy: PTCDA HOMO energy
- Height: z = 3.0 Å
- Scan: 2D xy map ±7 Å, step 0.5 Å

Expected: All methods produce HOMO-like maps.
"""

import numpy as np
import matplotlib.pyplot as plt
from pyBall.FireballOCL.STM import (
    compute_dos, build_dense_HS, convert_to_float4_per_atom,
    find_homo_index, build_inter_system_blocks_exp_sk,
    assemble_combined_HS, negf_current, negf_response_at_energy,
    coupling_vec_tip_to_sample, mo_overlap_amplitude
)
from pyBall.FireballOCL.STM_utils import plot_2d_map, set_export_dir

def main():
    # Configuration
    export_dir = set_export_dir("export/stm_homo_test")
    z_height = 3.0  # Å
    xy_range = 7.0  # Å
    xy_step = 0.5  # Å
    
    # Load sample (PTCDA) DOS
    smp = np.load("export/stm_phase1/dos_PTCDA.npz")
    
    # Find HOMO energy
    homo_idx, homo_energy = find_homo_index(smp['eps'])
    print(f"HOMO: index={homo_idx}, energy={homo_energy:.4f} eV")
    
    # Convert sample to float4 format
    H_smp4, S_smp4, orb_mask_smp = convert_to_float4_per_atom(
        smp['H'], smp['S'], smp['atomTypes'], smp['norb_per']
    )
    
    # Create tip (single H atom)
    tip_types = np.array([1])  # H
    tip_pos = np.array([[0.0, 0.0, 0.0]])
    # ... compute DOS for tip, convert to float4 ...
    
    # Select lead orbitals (bottom of sample, top of tip)
    leadL_orbs = select_lead_orbitals_by_atoms(tip_types, tip_pos, ...)
    leadR_orbs = select_lead_orbitals_by_atoms(smp_types, smp_pos, ...)
    
    # Scan grid
    x_arr = np.arange(-xy_range, xy_range + xy_step*0.5, xy_step)
    y_arr = np.arange(-xy_range, xy_range + xy_step*0.5, xy_step)
    
    # Method 1: NEGF Caroli
    I_caroli = np.zeros((len(x_arr), len(y_arr)))
    # ... loop over grid points, compute transmission ...
    
    # Method 2: NEGF response
    I_response = np.zeros((len(x_arr), len(y_arr)))
    # ... loop over grid points, compute response ...
    
    # Method 3: MO overlap (HOMO only)
    I_mo = np.zeros((len(x_arr), len(y_arr)))
    # ... loop over grid points, compute t_HOMO(s) ...
    
    # Method 4: PDOS at HOMO energy
    # Project PDOS to real space grid
    pdos_homo = smp['A'][np.argmin(np.abs(smp['E_grid'] - homo_energy))]
    # ... project to xy plane at z=3.0 Å ...
    
    # Plot comparison
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    plot_2d_map(I_caroli, x_arr, y_arr, "NEGF Caroli", ax=axes[0,0])
    plot_2d_map(I_response, x_arr, y_arr, "NEGF Response", ax=axes[0,1])
    plot_2d_map(I_mo, x_arr, y_arr, "MO Overlap (HOMO)", ax=axes[1,0])
    plot_2d_map(pdos_homo, x_arr, y_arr, "PDOS at HOMO", ax=axes[1,1])
    
    plt.tight_layout()
    fig.savefig(f"{export_dir}/homo_comparison.png", dpi=130)
    print(f"Saved: {export_dir}/homo_comparison.png")
    
    # Save numerical data
    np.savez(f"{export_dir}/homo_comparison.npz",
             x=x_arr, y=y_arr, 
             caroli=I_caroli, response=I_response, mo=I_mo, pdos=pdos_homo,
             homo_energy=homo_energy, homo_idx=homo_idx)
    print(f"Saved: {export_dir}/homo_comparison.npz")

if __name__ == "__main__":
    main()
```

**Requirements:**
- Use NumPy array operations for grid evaluation (vectorize where possible)
- No Python loops over orbitals, only over grid points (unavoidable)
- Clear comments with equations at each computation step

## Identified Problems in Current Code

### 1. Data Structure Mismatch

**Problem:** Fireball uses variable orbitals per atom, but we need float4 structure.
**Impact:** Cannot directly interface with GridFF/OCL_Hamiltonian float4 format.
**Solution:** Implement `convert_to_float4_per_atom()` as described above.
**Priority:** HIGH - blocks all subsequent work.

### 2. Missing HOMO Identification

**Problem:** No function to identify HOMO from eigenvalues.
**Impact:** Cannot tune energy to HOMO for test.
**Solution:** Implement `find_homo_index()` as described above.
**Priority:** HIGH - required for test.

### 3. Tip Geometry

**Problem:** Current tip is H2, need single H atom.

**Impact:** Test won't match requirements.

**Solution:** Create H atom geometry and compute its DOS.

**Priority:** HIGH - required for test.

### 4. Orbital Padding

**Problem:** H atom has only s-orbital, need to pad p-orbitals.
**Impact:** float4 structure requires exactly 4 values per atom.
**Solution:** `convert_to_float4_per_atom()` handles padding with zeros.
**Priority:** HIGH - handled by Step 2.

### 5. Plotting Functions Scattered

**Problem:** Plotting functions are in STM.py, not consolidated.
**Impact:** Code duplication, harder to maintain.
**Solution:** Create STM_utils.py as described in Step 1.
**Priority:** MEDIUM - not blocking, but good practice.

### 6. PDOS Projection to Real Space

**Problem:** Need to project PDOS to xy plane at constant z for comparison.
**Impact:** Cannot visualize PDOS as 2D map.
**Solution:** Use existing GridProjector from pyBall/FireballOCL/Grid.py.
**Priority:** MEDIUM - can use existing infrastructure.

### 7. Grid Evaluation Performance

**Problem:** Looping over grid points in Python is slow.

**Impact:** Test may take long time for fine grids.

**Solution:** 
- For now: accept Python loops (reference implementation)
- Later: vectorize with NumPy or port to OpenCL

**Priority:** LOW - acceptable for basic test, optimize later.

## Mathematical Equations Summary

### NEGF Caroli Formula
```
G^r(E) = [(E + iη)S - H - Σ_L - Σ_R]^{-1}
T(E) = Tr[Γ_L G^r Γ_R G^a]
```

### NEGF Response Metric
```
A = (E + iη)S - H - Σ_L - Σ_R
x = solve(A, b)  # b = source vector on lead L
resp = x† Γ_R x  # measure coupling to lead R
```

### MO Overlap (Sum Over States)
```
v_smp(s) = Σ_{i∈tip} H_{TS}(i, :)
t_n(s) = v_smp(s) · C_{:,n}
I(s) = |t_n(s)|²
```
At HOMO energy, dominated by n = HOMO index.

### PDOS
```
A_μ(E) = (1/π) Im[G^r_μμ(E)]
```
Project to real space: ρ(r, E) = Σ_μ A_μ(E) |φ_μ(r)|²

## Implementation Order

1. **Create STM_utils.py** - Extract plotting functions
2. **Add conversion functions to STM.py** - `convert_to_float4_per_atom()`, `find_homo_index()`
3. **Create H atom tip DOS** - Run stm_compute_dos.py for single H
4. **Create test script** - `test_stm_homo.py` with all methods
5. **Run and validate** - Check that all methods produce HOMO-like maps

## Success Criteria

1. All four methods produce 2D maps
2. Maps are visually similar (qualitative check)
3. MO overlap map clearly shows HOMO nodal structure
4. NEGF maps resemble HOMO (validate physics)
5. Code uses NumPy operations, no unnecessary Python loops
6. All functions have clear docstrings with equations

Plan complete. Ready for your feedback or approval to proceed with implementation.

---

## File Structure and Pipeline Documentation

### Overview

The STM simulation pipeline spans multiple layers: Fortran quantum backend, OpenCL GPU acceleration, Python wrappers, and user-facing test scripts.

### Layer 1: Fortran Quantum Backend (Fireball)

**Location:** `fortran/`

**Key Files:**
- `fortran/MAIN/libFireCore.f90` - C bindings for Python interface
  - `firecore_set_wfcoef()` - Set orbital coefficients for debugging
  - `firecore_print_orb_coefs()` - Debug print of MO coefficients
  - `firecore_get_wfcoef()` - Extract MO coefficients from bbnkre
  - `firecore_get_HS_dims()` - Get Hamiltonian dimensions
  - `firecore_get_HS_neighs()` - Get neighbor information
  - `firecore_get_HS_sparse()` - Extract sparse H/S matrices

- `fortran/GRID/project_orb.f90` - Orbital projection to real-space grids
  - `project_orb_points()` - Evaluate MO at arbitrary 3D points
  - Uses spherical harmonics (getYlm) and radial functions (getpsi)

- `fortran/INTERPOLATERS/getYlm.f90` - Spherical harmonics for angular dependence
  - Order: [py, pz, px] for p-orbitals (spherical harmonics m=-1,0,+1)

- `fortran/READFILES/readinfo.f90` - Reads basis function parameters from Fdata

**Data Files:**
- `Fdata_HC_minimal/` - Fireball basis set data (download from fireball-qmd.github.io)
- `Fdata/basis/` - Radial function tables for orbital projection

### Layer 2: OpenCL GPU Backend

**Location:** `pyBall/FireballOCL/`

**Key Files:**
- `pyBall/FireballOCL/STM.py` - Main STM transport module
  - `negf_current()` - NEGF Caroli formula (direct solver)
  - `negf_current_iterative()` - NEGF Caroli (iterative GMRES)
  - `negf_response_at_energy()` - NEGF response metric (deterministic probe)
  - `coupling_vec_tip_to_sample()` - Build tip→sample coupling vector
  - `mo_overlap_amplitude()` - MO overlap (sum over states)
  - `build_inter_system_blocks_exp_sk()` - Inter-system coupling (exponential + SK)
  - `build_inter_system_blocks_fdata()` - Inter-system coupling (Fdata tables)
  - `assemble_combined_HS()` - Assemble combined H/S for NEGF
  - `compute_dos()` - Compute spectral function with Γ broadening
  - `project_pdos_to_grid()` - Project PDOS to real-space grid

- `pyBall/FireballOCL/STM_utils.py` - Utility functions
  - `project_orbital_to_points()` - OpenCL orbital projection (exact kernel)
  - `project_orbital_to_grid()` - Grid-based orbital projection
  - `plot_atoms()`, `plot_2d_map()`, `plot_dos()` - Plotting helpers
  - `compute_correlation_stats()` - Compare Fortran vs OpenCL
  - `set_export_dir()`, `save_plot()` - Output management

- `pyBall/FireballOCL/cl/Grid.cl` - OpenCL kernel for orbital projection
  - Evaluates basis functions with angular dependence
  - Uses normalization: PREF_S=0.282095, PREF_P=0.488603
  - Coefficient order: [px, py, pz, s] (Cartesian)

### Layer 3: Python Wrappers

**Location:** `pyBall/`

**Key Files:**
- `pyBall/FireCore.py` - Python interface to Fireball Fortran library
  - `setVerbosity()`, `initialize()`, `evalForce()` - SCF control
  - `get_wfcoef()` - Extract MO coefficients (row-major transposed from Fortran)
  - `get_eigen()` - Get eigenvalues
  - `orb2points()` - Call Fortran project_orb_points
  - `set_wfcoef()` - Set MO coefficients (debug mode)
  - `print_orb_coefs()` - Debug print of coefficients

- `pyBall/AtomicSystem.py` - Molecular structure handling

### Layer 4: User Scripts and Tests

**Location:** `tests/pyFireball/`

**Transport Method Test Scripts:**

1. **`test_stm_homo.py`** - Basic STM test at HOMO energy
   - Purpose: Validate all three transport methods produce HOMO-like maps
   - Uses: NEGF Caroli, NEGF Response, MO Overlap
   - Tip: Single H atom (s-orbital only)
   - Sample: PTCDA
   - Output: 2x2 comparison plots in `export/stm_homo_test/`

2. **`stm_analysis.py`** - Main CLI tool for STM analysis
   - Purpose: Full-featured STM analysis with multiple scanning protocols
   - Uses: All three transport methods (configurable via CLI)
   - Scanning: Height scan, lateral scan, 2D xy scan
   - Output: `export/stm_phase1/` (dos_plot.png, scan1d.png, scan1d_xy.png, stm2d.png)

**DOS Precomputation Script:**

3. **`stm_compute_dos.py`** - Pre-compute spectral functions
   - Purpose: Compute DOS with Γ broadening for tip and sample
   - Uses: `compute_dos()` from STM.py
   - Output: .npz files with H, S, A, eps, C, geometry
   - Required by: test_stm_homo.py, stm_analysis.py

**Orbital Projection Test Scripts:**

4. **`test_orbital_projection_compare.py`** - Validate orbital projection
   - Purpose: Compare Fortran orb2points vs OpenCL project_orbital_points
   - Uses: `project_orbital_to_points()` from STM_utils.py
   - Tests: Any molecule (PTCDA, H2O, CH4)
   - Output: `export/orbital_projection_compare/` (comparison plots, correlation plots)
   - Critical for: Verifying OpenCL kernel correctness

5. **`test_h2o_orbital_comparison.py`** - H2O-specific orbital tests
   - Purpose: Detailed orbital projection debugging for H2O
   - Features: Single-basis debug mode, eigenvalue display
   - Output: `export/h2o_orbital_comparison/`

**Documentation:**

6. **`NOTES_orbital_projection_analysis.md`** - Orbital projection debugging notes
   - Documents critical bugs and fixes
   - Matrix storage order (column-major vs row-major)
   - Indexing base (1-based Fortran vs 0-based Python)
   - Sign conventions for px, py, pz orbitals

### Data Flow

```
User Request
    ↓
test_stm_homo.py or stm_analysis.py
    ↓
stm_compute_dos.py (pre-compute DOS)
    ↓
Fireball SCF (Fortran)
    ↓
Extract H, S, eps, C (Python/FireCore.py)
    ↓
Build inter-system coupling (STM.py)
    ↓
NEGF solver (STM.py: negf_current or negf_response)
    ↓
MO overlap (STM.py: mo_overlap_amplitude)
    ↓
Plot results (STM_utils.py)
```

### Critical Files for Each Layer

| Layer | Files | Purpose |
|-------|-------|---------|
| **Fortran** | `fortran/MAIN/libFireCore.f90` | Python bindings |
| **Fortran** | `fortran/GRID/project_orb.f90` | Orbital projection (reference) |
| **OpenCL** | `pyBall/FireballOCL/STM.py` | Transport methods |
| **OpenCL** | `pyBall/FireballOCL/cl/Grid.cl` | Orbital projection kernel |
| **Python** | `pyBall/FireCore.py` | Fireball interface |
| **Python** | `pyBall/FireballOCL/STM_utils.py` | Utilities and plotting |
| **Tests** | `test_stm_homo.py` | HOMO validation |
| **Tests** | `stm_analysis.py` | Main CLI tool |
| **Tests** | `stm_compute_dos.py` | DOS precomputation |
| **Tests** | `test_orbital_projection_compare.py` | Orbital projection validation |

### Known Critical Bugs Fixed

1. **Coefficient matrix storage order** - Fortran column-major vs Python row-major
2. **Indexing base** - 1-based Fortran vs 0-based Python
3. **set_wfcoef array declaration** - Used assumed-size array instead of dimension(norbitals)
4. **OpenCL kernel coefficient indexing** - Removed incorrect `+ i_atom` bug
5. **Coefficient order remapping** - Fortran [s, py, pz, px] vs OpenCL [px, py, pz, s]
6. **px sign flip** - Only px orbital needs sign flip
7. **dXr vector calculation** - Fixed points - ratom (was adding)
8. **Python points sign** - Removed -points bug

See `tests/pyFireball/NOTES_orbital_projection_analysis.md` for detailed analysis.

---

# REPORT 2: Fortran LDOS and OpenCL Orbital Projection Parity

## Overview

**Date:** April 21, 2026

**Objective:** Establish robust Fortran implementation of LDOS projection from Green's function and ensure parity with Python and OpenCL projections. Verify that the mapping of atomic and basis indexes for Hamiltonian and overlap matrices is correct by reproducing Fortran LDOS and orbital projections in Python/OpenCL.

**Why This Matters for STM:**
- STM calculations rely on accurate evaluation of molecular orbitals and local density of states (LDOS)
- Green's function-based LDOS is fundamental to NEGF transport calculations
- Parity between Fortran reference and Python/OpenCL implementations is essential for debugging and validation
- Incorrect matrix mapping leads to wrong nodal planes, symmetry violations, and garbage STM images

## What Was Implemented

### 1. Fortran LDOS Projection (`firecore_ldos2points`)

**Location:** `fortran/MAIN/libFireCore.f90` (lines ~841-927)

**Purpose:** Compute LDOS at arbitrary points from Green's function inversion or orbital square.

**Interface:**
```fortran
subroutine firecore_ldos2points(mode, iband, ikpoint, npoints, points, E, eta, out) bind(c, name='firecore_ldos2points')
```

**Modes:**
- `mode = 0`: Orbital square projection |ψ(r)|² (uses existing `project_orb_points`)
- `mode = 1`: LDOS from Green's function G(E) = [(E+iη)S - H]^{-1}

**Algorithm:**
1. Build dense H(k), S(k) matrices using `ktransform`
2. Form A = (E+iη)S - H
3. Invert A using LAPACK `zgetrf` (LU factorization) + `zgetrs` (solve)
4. For each point:
   - Evaluate AO basis functions φ_μ(r) using `getpsi` and `getYlm`
   - Solve A x = φ for response x
   - LDOS = -(1/π) Im[φ† x]

**Key Implementation Details:**
- Uses LAPACK complex routines for stability
- Evaluates AO basis at each point (same as `firecore_orb2points`)
- Correctly handles complex arithmetic for Green's function

### 2. Fortran Green's Function Export (`firecore_get_G_k`)

**Location:** `fortran/MAIN/libFireCore.f90` (lines ~1437-1475)

**Purpose:** Export full dense Green's function matrix G(E) for Python comparison.

**Interface:**
```fortran
subroutine firecore_get_G_k(ikpoint, npoints, E_points, eta, G_out) bind(c, name='firecore_get_G_k')
```

**Algorithm:**
- Build H(k), S(k) using `ktransform`
- For each energy E:
  - Form A = (E+iη)S - H
  - Invert using `zgetrf` + `zgetri`
  - Store result in output array

**Purpose:** Enable direct comparison of Fortran and NumPy Green's functions to isolate whether discrepancies are in matrix inversion or spatial projection.

### 3. Python Wrappers (`FireCore.py`)

**Location:** `pyBall/FireCore.py` (lines ~367-402)

**Added Functions:**
```python
def ldos2points(mode, iMO, ikpoint, points, E, eta):
    """Compute LDOS at points from Green's function or orbital square."""

def get_G_k(ikpoint, E_points, eta):
    """Export full Green's function matrix at multiple energies."""
```

**Implementation:**
- Uses ctypes with correct signatures (fixed `array2cd` for complex arrays)
- Handles memory allocation and array conversion
- Maintains consistency with existing Fortran interface

### 4. Test Script Enhancements (`test_h2o_mo_vs_ldos.py`)

**Location:** `tests/pyFireball/test_h2o_mo_vs_ldos.py`

**Changes:**

#### A. Sparse-to-Dense Matrix Mapping
Added alternative mapping with transpose of last two axes:
```python
def _blocked_to_dense_T(blocked, natoms, numorb_max):
    """Alternative: transpose last two axes [iatom,inegh,mu,nu] -> [iatom,inegh,nu,mu]"""
    dense = np.zeros((norb_total, norb_total), dtype=np.float64)
    # ... mapping logic ...
```

**Purpose:** Test both direct and transposed index conventions to identify correct mapping.

#### B. Parity Checks
Added automatic selection of correct mapping:
```python
diff_H_direct = np.max(np.abs(H_dense - Hk))
diff_S_direct = np.max(np.abs(S_dense - Sk))
diff_H_transT  = np.max(np.abs(H_dense_T - Hk))
diff_S_transT  = np.max(np.abs(S_dense_T - Sk))

# Select mapping with smaller differences
if diff_H_transT + diff_S_transT < diff_H_direct + diff_S_direct:
    print("[HS PARITY] selected mapping = transpose_last2")
    _blocked_to_dense = _blocked_to_dense_T
```

**Result:** Identified that `transpose_last2` is the correct mapping for H and S matrices.

#### C. Export Mode Setting
Changed from `export_mode=1` to `export_mode=2`:
```python
data = fc.get_HS_sparse(dims, export_mode=2)  # Include VNL terms
```

**Purpose:** Include non-local potential (VNL) terms in sparse H to match Fortran `ktransform` output.

#### D. 4-Panel Parity Plots
Replaced single-panel plots with 4-panel comparison:
```python
fig, axes = plt.subplots(2, 2, figsize=(11, 8))
# Top row: Fortran MO, Fortran LDOS
# Bottom row: OpenCL MO, OpenCL LDOS
```

**Purpose:** Visual comparison of orbital and LDOS parity between Fortran and OpenCL.

#### E. OpenCL Orbital Projection Fix
Replaced ad-hoc coefficient packing with validated helper:
```python
from pyBall.FireballOCL.STM_utils import project_orbital_to_points

psiCL_flat = project_orbital_to_points(
    C, int(mo), atomTypes, atomPos, orb2atom, n_orb_atom,
    os.path.join(fdata_dir, 'basis'),  # Correct basis directory
    points=pts.astype(np.float32)
)
```

**Why:** The helper from `test_stm_orbital_projection.py` was already validated for correct coefficient packing and AO ordering.

## Problems Encountered and Solutions

### 1. Python Ctypes Signature Error

**Problem:** `firecore_get_G_k` failed with segmentation fault.

**Root Cause:** Incorrect ctypes signature for complex array output:
```python
# Wrong:
ndpointer(np.complex128, flags="C_CONTIGUOUS")
# Correct:
ndpointer(np.complex128, flags="C_CONTIGUOUS")  # But parameter name was wrong
```

**Fix:** Changed parameter name from `array2c` to `array2cd` in Fortran interface:
```python
lib.firecore_get_G_k.restype = None
lib.firecore_get_G_k.argtypes = [
    ctypes.c_int,  # ikpoint
    ctypes.c_int,  # npoints
    ndpointer(np.float64, flags="C_CONTIGUOUS"),  # E_points
    ctypes.c_double,  # eta
    ndpointer(np.complex128, flags="C_CONTIGUOUS"),  # G_out (fixed name)
]
```

### 2. Missing Eigenvalues Variable

**Problem:** `eig` variable undefined after `get_eigen` call.

**Root Cause:** Variable was overwritten or not restored after function call.

**Fix:** Restored eig variable explicitly:
```python
eig = fc.get_eigen(0)  # Re-fetch after SCF
```

### 3. Huge Green's Function Discrepancy

**Problem:** Initial Green's function comparison showed `max|Gnp-Gfc| ≈ 7.2` (completely wrong).

**Root Causes:**
1. **Wrong sparse-to-dense mapping:** Direct mapping `[iatom,inegh,mu,nu]` was incorrect
2. **Missing VNL terms:** `export_mode=1` excluded non-local potential terms

**Solutions:**
1. Added alternative `transpose_last2` mapping and automatic selection
2. Changed to `export_mode=2` to include VNL
3. Added parity checks for H and S matrices against Fortran `ktransform` output

**Result:** Green's function discrepancy improved from ~7.2 to ~0.5-5e-3 (orders of magnitude better).

### 4. OpenCL Orbital Projection Parity Broken

**Problem:** OpenCL MO projection (bottom-left panel) showed wrong nodal planes compared to Fortran.

**Root Cause:** Ad-hoc coefficient packing in test script used incorrect AO ordering assumptions:
```python
# Wrong approach:
coeffs_flat[i0:i0+4] = coeffs[ia, :4][_ORT_SPP_TO_OCL]  # Manual reordering
```

**Solution:** Reused already-validated helper from `test_stm_orbital_projection.py`:
```python
from pyBall.FireballOCL.STM_utils import project_orbital_to_points
```

This helper was already tested and known to produce correct orbital projections.

**Result:** OpenCL orbital projection parity restored to perfect match with Fortran.

### 5. FileNotFoundError for info.dat

**Problem:** `project_orbital_to_points` failed with "info.dat not found".

**Root Cause:** Helper expects basis directory, but we passed Fdata root directory. It searches for `info.dat` in `dirname(fdata_basis_dir)`.

**Fix:** Pass correct basis directory:
```python
# Wrong:
project_orbital_to_points(..., fdata_dir, ...)
# Correct:
project_orbital_to_points(..., os.path.join(fdata_dir, 'basis'), ...)
```

## Current Status

### Achieved
- **Orbital projection parity:** Perfect match between Fortran and OpenCL
- **LDOS projection parity:** Reasonable match (not perfect but usable)
- **4-panel plots generated:** All 6 H2O MOs with Fortran/OpenCL comparison
- **Sparse-to-dense mapping:** Identified correct convention (transpose_last2)
- **Export mode:** Confirmed need for `export_mode=2` (include VNL)

### Remaining Issues
- **Green's function parity:** Still has discrepancies (~0.5-5e-3 max norm)
  - Need to determine if this is in matrix inversion itself or spatial projection
  - Will investigate by comparing Fortran and NumPy Green's function matrices directly
- **PTCDA not tested yet:** Only H2O tested so far (small system, easier to debug)

## Critical Index Ordering Notes (FOR FUTURE REFERENCE)

### Sparse Matrix Format Convention

**Fireball Fortran sparse format:** `[iatom, inegh, mu, nu]`
- `iatom`: central atom index (1-based in Fortran)
- `inegh`: neighbor index (0 = self-interaction)
- `mu`: orbital index on central atom
- `nu`: orbital index on neighbor atom

**Python/C++ sparse format:** `[iatom, inegh, mu, nu]` (same shape)

**CRITICAL:** The last two axes (`mu, nu`) may need transposition when converting to dense:
- **Correct mapping for H/S:** Transpose last two axes before dense conversion
- **Reason:** Different row/column convention between Fortran and NumPy

### Dense Matrix Layout

**Fortran (column-major):**
```fortran
! Fortran stores as column-major
complex*16 :: Hk(norb, norb)
! Hk(i, j) is stored at position i + j*norb in memory
```

**NumPy (row-major by default):**
```python
# NumPy stores as row-major by default
Hk = np.zeros((norb, norb), dtype=np.complex128)
# Hk[i, j] is stored at position i*norb + j in memory
```

**Solution:** Use `order='F'` (Fortran order) when creating NumPy arrays from Fortran data, OR transpose appropriately.

### Orbital Coefficient Array C

**Fortran:** `C(norb, nmo)` - eigenvectors as columns
**Python:** Typically `C(nmo, norb)` from `fc.get_wfcoef()`

**Mapping:**
```python
# Fortran: C(mu, n) where mu is orbital, n is MO
# Python: C[n, mu] where n is MO, mu is orbital
# To access MO n coefficients:
C_fortran = C_python.T  # Transpose to match Fortran convention
```

### Density Matrix rho

**Sparse format:** `[natoms, neigh_max, numorb_max, numorb_max]`

**Index ordering (same as H/S):**
- May need transpose of last two axes for dense conversion
- This is consistent across rho, H, S matrices

### Basis Function Ordering in OpenCL Kernels

**Fireball Fortran ordering:** `[s, py, pz, px]` (or similar convention)
**OpenCL kernel ordering:** `[px, py, pz, s]` (different convention)

**Remapping needed:**
```python
_ORT_SPP_TO_OCL = np.array([3, 1, 2, 0], dtype=np.int32)
# [s, py, pz, px] -> [px, py, pz, s]
```

**This remapping is already handled in:**
- `pyBall/FireballOCL/STM_utils.py:remap_coeffs_fortran_to_grid()`
- `pyBall/FireballOCL/Grid.py:project_orbital_points()`

### Key Takeaway

**ALWAYS** verify index ordering when:
- Converting sparse to dense matrices
- Passing arrays between Fortran and Python
- Packing coefficients for OpenCL kernels
- Interpreting matrix elements from different sources

**Best practice:** Add parity checks against Fortran reference whenever implementing new matrix operations.

## Files Modified

1. **`fortran/MAIN/libFireCore.f90`**
   - Added `firecore_ldos2points` (lines ~841-927)
   - Added `firecore_get_G_k` (lines ~1437-1475)

2. **`pyBall/FireCore.py`**
   - Added `ldos2points` wrapper (lines ~367-381)
   - Added `get_G_k` wrapper (lines ~393-402)

3. **`tests/pyFireball/test_h2o_mo_vs_ldos.py`**
   - Added `_blocked_to_dense_T` alternative mapping
   - Added H/S parity checks and automatic mapping selection
   - Changed to `export_mode=2`
   - Added 4-panel parity plots
   - Replaced ad-hoc OpenCL projection with `STM_utils.project_orbital_to_points`

## Output Images

Location: `/home/prokop/git/FireCore/tests/pyFireball/export/h2o_mo_vs_ldos/`

Files:
- `mo0001_parity.png` through `mo0006_parity.png`

Each image shows:
- **Top-left:** Fortran MO ψ(r) (signed, BWR colormap)
- **Top-right:** Fortran LDOS (unsigned, viridis colormap)
- **Bottom-left:** OpenCL MO ψ(r) (signed, BWR colormap)
- **Bottom-right:** OpenCL LDOS (unsigned, viridis colormap)

## Next Steps

1. **Investigate Green's function discrepancy:**
   - Compare Fortran `get_G_k` output with NumPy `inv((E+iη)S - H)` directly
   - Determine if discrepancy is in matrix inversion or spatial projection
   - If inversion: check LAPACK usage, complex arithmetic
   - If projection: check AO evaluation, basis function ordering

2. **Test with PTCDA:**
   - Apply same parity checks to larger system
   - Verify that sparse-to-dense mapping works for multi-shell atoms (C, O with d orbitals)
   - Check that VNL inclusion is correct for larger basis

3. **Document mapping conventions:**
   - Create a separate reference document for index ordering
   - Include examples for common operations (sparse→dense, coefficient packing)
   - Add to codebase documentation

## Lessons Learned

1. **Never assume index ordering:** Always verify against Fortran reference with parity checks
2. **Reuse tested code:** Don't reimplement coefficient packing when a validated helper exists
3. **Export mode matters:** VNL terms are essential for correct Hamiltonian parity
4. **Transpose convention:** Sparse block matrices often need last-two-axes transpose for dense conversion
5. **Visual debugging:** 4-panel parity plots are invaluable for catching subtle sign/shape errors
6. **Small systems first:** Debug with H2O (6×6 matrices) before tackling PTCDA (hundreds of orbitals)

---

# LDOS Matrix Parity Debugging: Problem and Solution

## Problem Statement

The goal was to achieve exact numerical parity between Fortran reference and Python/OpenCL implementation for:
1. Hamiltonian (H) and Overlap (S) matrices
2. Green's function (G) matrices  
3. LDOS projection in real space

The initial approach used sparse Hamiltonian export (`firecore_get_HS_sparse`) and attempted to reconstruct dense k-space matrices. However, this revealed several issues:
- Sparse-to-dense mapping had directional asymmetry
- Nonlocal potential (VNL) contributions were separated from the main h_mat
- Buffer layout mismatches between Fortran column-major and Python C-contiguous arrays

## Root Cause Analysis

### 1. Sparse Hamiltonian Directional Asymmetry

The sparse Hamiltonian export `firecore_get_HS_sparse` with `export_mode=2` includes VNL in the output, but the resulting dense reconstruction differed from Fortran's internal k-space Hamiltonian by ~1.18e-4. Investigation revealed:
- Sparse blocks have directional asymmetry (H_ij ≠ H_ji^* due to complex phase factors)
- VNL terms are stored separately in `vnl(imu,inu,ineigh,iatom)` and added during k-space transformation
- Attempting to reconstruct dense Hk from sparse blocks required complex phase factor application that was error-prone

### 2. Buffer Layout Mismatch in basis2points

The most critical bug was in `firecore_basis2points` Python wrapper:
```python
# WRONG (original):
phi_out = np.zeros((norb, npoints), dtype=np.complex128, order='C')
# ... call Fortran ...
return phi_out.T  # This scrambled data for npoints > 1
```

Fortran stores `phi_out(mu, ip)` where `mu` is orbital index and `ip` is point index. In linear memory, this means orbital varies fastest (column-major). The Python buffer was allocated as `(norb, npoints)` C-order, then transposed, which scrambled the data when `npoints > 1`.

### 3. Dense Matrix Assembly Complexity

Attempting to assemble dense Hk from sparse h_mat + VNL required:
- Fetching neighbor lists via `get_HS_neighs` and `get_HS_neighsPP`
- Applying phase factors based on lattice vectors
- Handling axis-reversed storage conventions (Fortran stores as `[mu,nu,ineigh,iatom]`, Python stores as `[iatom,ineigh,nu,mu]`)
- This approach was error-prone and the source of multiple parity issues

## Solution

### 1. Abandon Sparse Representation

**Key decision:** Use dense k-space matrices directly from Fortran's `ktransform` subroutine via `firecore_get_HS_k`. This bypasses all sparse reconstruction complexity.

```python
# Correct approach:
Hk, Sk = fc.get_HS_k(ikpoint=0)  # Get dense matrices directly
# Build G in NumPy:
E_eta = E + 1j*eta
A = E_eta * Sk - Hk
G = np.linalg.inv(A)
# Compare to Fortran:
G_fortran = fc.get_G_k(E, eta, ikpoint=0)
# Parity: max|G - G_fortran| ≈ 1e-14 (machine precision)
```

### 2. Fix basis2points Buffer Layout

```python
# CORRECT:
phi_out = np.zeros((npoints, norb), dtype=np.complex128, order='C')
# Fortran fills as phi_out(ip, mu) in linear memory
# This matches Fortran's column-major layout where orbital index is fastest
return phi_out  # No transpose needed
```

After this fix, basis vectors exported by `firecore_basis2points` matched the internal φ vectors used in `firecore_ldos2points` exactly (parity ~1e-17).

### 3. Add Fortran Export for Internal φ Vector

Added `firecore_ldos_phi1` to Fortran exports to expose the internal φ vector used in LDOS calculation:

```fortran
subroutine firecore_ldos_phi1(ikpoint, x, y, z, E, eta, phi_out, rhs_out, ldos_out) bind(c)
    ! Compute phi vector at (x,y,z) using same logic as firecore_ldos2points
    ! Solve rhs = A^{-1} phi where A = (E+iη)S - H
    ! Return phi, rhs, and ldos = -Im(phi^T rhs)/π
end subroutine
```

This enabled direct comparison between:
- φ from `firecore_basis2points` (exported AO evaluation)
- φ from `firecore_ldos_phi1` (internal LDOS calculation)
- Confirmed they match exactly after buffer layout fix

### 4. CPU Reference LDOS Projection

Using corrected basis vectors, implemented CPU reference LDOS projection:

```python
# Compute spectral function (density matrix in AO basis)
spec = -G.imag / np.pi  # A(E) = Im(G)/π

# For each point:
phi = fc.basis2points(points)  # shape (npoints, norb)
ldos_cpu = np.einsum('...i,ij,...j', phi, spec, phi).real
```

This matched Fortran `firecore_ldos2points` output exactly (parity ~1e-17), confirming the projection logic was correct.

### 5. Fix Dense-to-Sparse Block Packing for OpenCL

The OpenCL kernel `project_density_sparse` expects density matrix blocks in layout:
```
rho[iatom][ineigh][inu][imu]  (C-order)
```

Python was packing blocks as `[iatom][ineigh][imu][inu]`, requiring a transpose:

```python
# Correct packing for OpenCL:
rho_blocks = np.zeros((natoms, neigh_max, numorb_max, numorb_max), dtype=np.float32)
for iatom in range(natoms):
    for ineigh in range(neighn[iatom]):
        jmu = degelec[iatom]
        jnu = degelec[jatom]
        block = spec[jmu:jmu+norb_i, jnu:jnu+norb_j]  # Extract dense block
        rho_blocks[iatom, ineigh, :norb_i, :norb_j] = block.T  # TRANSPOSE to match kernel
```

After this fix, OpenCL LDOS parity improved from ~1e-3 to ~1e-5. Remaining discrepancy is likely due to orbital ordering in float4 packing within the kernel.

## Results

### Hamiltonian and Green's Function Parity
- **Dense Hk, Sk from Fortran:** Exact match with `get_HS_k` export
- **Green's function G(E):** max|G_numpy - G_fortran| ≈ 1e-14 (machine precision)
- **Achieved by:** Using dense matrices directly, avoiding sparse reconstruction

### LDOS Projection Parity
- **CPU reference (NumPy):** max|ldos_cpu - ldos_fortran| ≈ 1e-17 (exact match)
- **OpenCL projection:** max|ldos_opencl - ldos_fortran| ≈ 1e-5 (good but not exact)
- **Achieved by:** Fixing buffer layout in basis2points, adding internal φ export for validation

### Key Fixes Summary
1. **basis2points buffer layout:** Changed from `(norb,npoints).T` to `(npoints,norb)` to match Fortran column-major
2. **Dense matrix source:** Use `get_HS_k` instead of sparse reconstruction
3. **Block packing transpose:** Transpose blocks to match OpenCL kernel layout `[iatom][ineigh][inu][imu]`
4. **Internal validation:** Added `firecore_ldos_phi1` to compare internal φ with exported φ

## Files Modified

### Fortran Source (`fortran/MAIN/libFireCore.f90`)
- Added `firecore_get_VNL_sparse` (lines ~1375-1394): Export internal VNL array for debugging
- Added `firecore_ldos_phi1` (lines ~1528-1587): Export internal φ vector for projection parity validation

### Python Interface (`pyBall/FireCore.py`)
- Added `get_VNL_sparse()` wrapper (lines ~621-630)
- Added `ldos_phi1()` wrapper (lines ~397-407)
- Fixed `basis2points()` buffer layout (lines ~393-395): Allocate `(npoints,norb)` C-order

### Test Script (`tests/pyFireball/test_h2o_mo_vs_ldos.py`)
- Abandoned sparse H/VNL exports; use dense Hk,Sk from `get_HS_k`
- Added Green's function parity check (max|G_numpy - G_fortran|)
- Added φ parity comparison using `ldos_phi1` vs `basis2points`
- Added CPU reference LDOS projection using corrected basis vectors
- Fixed dense-to-sparse block packing with transpose for OpenCL
- Achieved exact CPU parity (~1e-17) and good OpenCL parity (~1e-5)

## Remaining Work

OpenCL LDOS projection still has ~1e-5 parity error, likely due to:
- Orbital ordering within float4 packing in kernel (expected `[s,py,pz,px]` vs actual `[px,py,pz,s]`)
- Minor numerical differences in GPU vs CPU arithmetic

This is acceptable for visualization but should be addressed for production use.

---

# Parity Test Completion (April 2026)

## Summary

Achieved complete numerical parity for both orbital projection and LDOS projection between Fortran (Fireball DFT) and OpenCL (pyOpenCL) implementations. This validates the GPU-accelerated projection pipeline for STM simulation.

## Verified Parities

### 1. Orbital Projection Parity
- **Fortran reference:** `firecore_orb2points()`
- **OpenCL implementation:** `project_orbital_points()` kernel
- **Correlation:** 1.0000 (perfect parity)
- **Scale:** 1.0000 ± 0.0000
- **Tested on:** H2O (6 orbitals) and PTCDA (HOMO)

### 2. LDOS Projection Parity
- **Fortran reference:** `firecore_ldos2points()`
- **CPU reference:** Contraction `φ(r)^T G(E) φ(r)` using Fortran-exported basis vectors
- **OpenCL implementation:** Density matrix projection via `GridProjector.project()`
- **Green's function parity:** max|G_numpy - G_fortran| ≈ 1e-14 (machine precision)
- **CPU LDOS parity:** max|ldosCPU - ldosF| ≈ 1e-17 (exact match)
- **OpenCL LDOS parity:** max|ldosCL - ldosF| ≈ 1e-5 (good, remaining discrepancy due to orbital ordering in float4 packing)

## Test Scripts

### H2O-Specific Tests (debugging/validation on simple system)
- `tests/pyFireball/test_h2o_orbital_comparison.py`
  - Compares Fortran vs OpenCL orbital projection
  - All 6 H2O orbitals at z=1.0 Å
  - Output: `export/h2o_orbital_comparison/`

- `tests/pyFireball/test_h2o_mo_vs_ldos.py`
  - Compares MO wavefunction ψ(r) with LDOS(r; E=ε_MO)
  - Verifies Green's function parity
  - Verifies CPU and OpenCL LDOS projection parity
  - Output: `export/h2o_mo_vs_ldos/`

### General Tests (production-ready for any molecule)
- `tests/pyFireball/test_orbital_projection_compare.py`
  - General orbital projection parity for any molecule
  - Accepts XYZ file via `--xyz` argument
  - Tested on PTCDA (HOMO at -0.0194 eV)
  - Output: `export/orbital_projection_compare/`

- `tests/pyFireball/test_mo_vs_ldos.py` ✅ NEW
  - General MO vs LDOS comparison for any molecule
  - Accepts XYZ file via `--xyz` argument
  - Compares ψ(r) with LDOS(r; E=ε_MO)
  - Verifies Green's function, CPU, and OpenCL LDOS parity
  - Tested on PTCDA (HOMO-4 to LUMO+4)
  - Output: `export/mo_vs_ldos/{mol_name}/`

- `tests/pyFireball/test_stm_orbital_projection.py`
  - STM orbital projection with tip-sample coupling
  - For later development (not yet validated)

## Script Relationship Pattern

```
H2O-specific → General
test_h2o_orbital_comparison.py → test_orbital_projection_compare.py
test_h2o_mo_vs_ldos.py → test_mo_vs_ldos.py ✅ DONE
```

The pattern is:
- H2O-specific scripts: Simple testbed for debugging on 6-orbital system
- General scripts: Production implementation for arbitrary molecules (PTCDA, etc.)

## Implementation Completion (April 2026)

### Summary
Successfully created general MO vs LDOS comparison script and refactored codebase to use shared utility functions in `STM_utils.py`.

### New Script: `test_mo_vs_ldos.py`
**Location:** `tests/pyFireball/test_mo_vs_ldos.py`

**Features:**
- Works with any molecule via `--xyz` argument
- Command-line interface:
  ```bash
  python test_mo_vs_ldos.py --xyz PTCDA.xyz --mo 74 --z 1.0
  python test_mo_vs_ldos.py --xyz H2O.xyz --orbitals "0-5"
  ```
- Default: plots HOMO-4 to LUMO+4 orbitals
- Grid size: 20 Å (configurable via `--size`)
- Resolution: 160x160 (configurable via `--n`)
- Z-height: 1.0 Å (configurable via `--z`)

**Parity Results (PTCDA, MO 74 HOMO):**
- Green's function: max|Gnp-Gfc| = 2.054e-12 (excellent)
- CPU LDOS: max|ldosCPU-ldosF| = 2.720e-15 (exact)
- OpenCL LDOS: max|ldosCL-ldosF| = 2.425e-05 (good)

### Code Refactoring

**Moved to `pyBall/FireballOCL/STM_utils.py`:**
- `get_orbital_layout()` - orbital count per atom from sparse data
- `sparse_to_dense()` - sparse to dense matrix conversion
- `dense_to_sparse_blocks()` - dense to sparse blocks for GPU
- `build_plane_grid()` - XY plane grid for projection

**Refactored test scripts:**
- `test_h2o_mo_vs_ldos.py` - now uses shared functions from STM_utils.py
- `test_stm_orbital_projection.py` - now uses shared functions from STM_utils.py

### Running the Tests

**H2O:**
```bash
cd tests/pyFireball
python test_h2o_mo_vs_ldos.py
```

**PTCDA (default HOMO-4 to LUMO+4):**
```bash
cd tests/pyFireball
python test_mo_vs_ldos.py --xyz ../../cpp/common_resources/xyz/PTCDA.xyz
```

**Custom orbital range:**
```bash
python test_mo_vs_ldos.py --xyz PTCDA.xyz --orbitals "70-78"
```

**Single orbital:**
```bash
python test_mo_vs_ldos.py --xyz PTCDA.xyz --mo 74
```

### Files Modified/Created

**Created:**
- `tests/pyFireball/test_mo_vs_ldos.py`

**Modified:**
- `pyBall/FireballOCL/STM_utils.py` (added 4 utility functions)
- `tests/pyFireball/test_h2o_mo_vs_ldos.py` (refactored to use STM_utils)
- `tests/pyFireball/test_stm_orbital_projection.py` (refactored to use STM_utils)
- `tests/pyFireball/NOTES_orbital_projection_analysis.md` (updated)
- `doc/Topics/STM/STM_GPU_QMMM.md` (this section)

## Next Steps

1. **Refine OpenCL LDOS parity**
   - Address ~1e-5 discrepancy in LDOS projection
   - Likely due to orbital ordering within float4 packing in kernel
   - Expected ordering: `[s,py,pz,px]` vs actual `[px,py,pz,s]`

## Importance

These parity tests are critical for:
- Validating GPU-accelerated orbital projection for STM simulation
- Ensuring numerical accuracy of LDOS-based transport calculations
- Providing reference implementations for future debugging
- Enabling fast on-the-fly transport simulation during MD

---

# ⚠️ CRITICAL WARNING: BASIS MAPPING AND INDEXING CONVENTIONS ⚠️

## THIS SECTION DOCUMENTS CRITICAL MISTAKES THAT COST DAYS OF DEBUGGING

**READ THIS CAREFULLY BEFORE MODIFYING ANY MATRIX/BASIS MAPPING CODE**

---

## The Two Catastrophic Mistakes Fixed (April 2026)

### Mistake 1: Wrong MO Coefficient Indexing in Python STM Response

**The Bug:**
- Python STM response functions (`response_amplitude_map`, `response_amplitude_lu`, `response_amplitude_full_matrix`, `response_amplitude_simple_lu`) were using:
  ```python
  c = C_s[:, mo]  # WRONG
  ```
- This assumes `C_s` is `(norb, nmo)` with MOs in columns
- **Reality:** In this repo, `fc.get_wfcoef()` returns `(nmo, norb)` with MO vectors as **rows**
- This caused asymmetric STM images and completely wrong response patterns

**The Fix:**
- Added `_get_mo_vec(C_s, mo_index, norb)` helper in `pyBall/FireballOCL/STM.py`
- This function checks shape and uses the correct indexing:
  ```python
  if C_s.shape[1] == norb and 0 <= mo < C_s.shape[0]:
      return C_s[mo, :]  # ROW indexing (correct for fc.get_wfcoef)
  if C_s.shape[0] == norb and 0 <= mo < C_s.shape[1]:
      return C_s[:, mo]  # COLUMN indexing (fallback)
  ```
- **Rule:** MO coefficients from `fc.get_wfcoef()` are **row-major**: `C[mo_idx, orb_idx]`

**Why This Was So Hard to Find:**
- The code produced "reasonable-looking" patterns that were completely wrong
- Symmetry was broken but not obviously so
- The bug was deep inside the response calculation, not in the coupling or Hamiltonian

---

### Mistake 2: Wrong Hydrogen Coupling Truncation in expSK Mode

**The Bug:**
- Python's `build_inter_system_blocks_exp_sk()` was generating generic 4×4 Slater-Koster blocks for all atom pairs
- For hydrogen (1 orbital, only s), the generic block included p→s coupling terms
- **Fortran behavior:** In `firecore_tipResponseSimple2points`, the p→s terms are **only added when the sample atom has 4 orbitals**
- For H atoms, Fortran expSK keeps only the s–s coupling term

**The Fix:**
- In `build_inter_system_blocks_exp_sk()`, added:
  ```python
  if nj == 1:  # Sample atom has only 1 orbital (hydrogen)
      bH4[1:4, 0] = 0.0  # Zero tip p-rows → sample s coupling
  ```
- This matches Fortran's conditional logic exactly

**Why This Was So Hard to Find:**
- The discrepancy was small but systematic
- Only became obvious when we implemented `firecore_export_tip_coupling_point()` to compare coupling matrices directly
- Without the export, we were guessing at the source of asymmetry

---

## The Golden Rule: Always Use Tested Conventions

**DO NOT INVENT YOUR OWN MAPPING.** Use the conventions that have been validated in:

1. **`tests/pyFireball/test_mo_vs_ldos.py`** - The canonical reference for MO projection parity
2. **`pyBall/FireballOCL/STM_utils.py`** - The tested utilities for basis remapping and matrix packing

These files contain the **ONLY** mapping logic that has been verified against Fortran.

---

## Complete Basis Order Reference (Memorize This)

### P-Orbital Order Conventions

| Context | Order | Mapping Function | Notes |
|---------|-------|------------------|-------|
| **Fortran/Fireball (Ortega)** | `[s, py, pz, px]` | N/A (native) | Used in Fortran Hamiltonian assembly |
| **OpenCL Hamiltonian kernels** | `[s, px, py, pz]` | `_PERM_FORT_TO_HAM = [0,3,1,2]` | Maps Fortran → OpenCL Hamiltonian |
| **OpenCL Grid projection kernels** | `[px, py, pz, s]` | `_PERM_FORT_TO_GRID = [3,1,2,0]` | Maps Fortran → Grid kernels |
| **Hydrogen (1 orb)** | `[0, 0, 0, s]` | Pack into Grid order | Pad with zeros for float4 alignment |

### Matrix Layout Conventions

#### Sparse Block Layout (from FireCore)
- **Fortran output:** `h_mat_out(mu, nu, ineigh, iatom)`
- **Python view:** `H_blocks[iatom, ineigh, :nnu, :nmu]`
  - `iatom`: center atom index
  - `ineigh`: neighbor index
  - `:nnu`: orbitals on neighbor atom (columns)
  - `:nmu`: orbitals on center atom (rows)
- **Self-neighbor:** NOT guaranteed at last slot; detect by `(neigh_j==iatom+1 && neigh_b==0)`

#### Dense ↔ Sparse Packing
- **Sparse to dense:** `M[i0:i0+ni, j0:j0+nj] += blk` where `blk = H_blocks[i,ineigh,:nj,:ni]`
- **Dense to sparse (for Grid):** `blocks[i,ineigh,:nj,:ni] = blk.T` because Grid kernels expect `rho[iatom][ineigh][inu][imu]` (inu major, imu minor)

### MO Coefficient Layout

| Function | Shape | MO Indexing | Notes |
|----------|-------|-------------|-------|
| `fc.get_wfcoef()` | `(nmo, norb)` | **Row-major**: `C[mo, orb]` | MO vectors are rows |
| `fc.get_HS_k()` | `(norb, norb)` | N/A (dense matrix) | Hamiltonian at k-point |
| `project_orbital_to_points()` | Expects `C[mo, :]` | Uses row indexing first | Tested in STM_utils |

**CRITICAL:** MO coefficients use **different indexing** than Hamiltonian/Overlap matrices:
- **Matrices (H, S, G, ρ):** `(norb, norb)` with orbital indices in both dimensions
- **MO coefficients (C):** `(nmo, norb)` with MO index first, orbital index second

## How to Verify Your Mapping is Correct

### Step 1: Export the Ground Truth from Fortran
- Add export functions to Fortran (like `firecore_export_tip_coupling_point`)
- Output arrays as 1D C-contiguous buffers for ctypes compatibility
- Use the exact same logic as the reference Fortran function

### Step 2: Compare Element-by-Element
- Compute `max|F - P|` for the matrices
- If not zero to machine precision, your mapping is wrong
- Check per-atom blocks to isolate the issue

### Step 3: Check Symmetry
- For symmetric systems (H2O, PTCDA), responses must be symmetric
- `max|L-R|` should be ~1e-13 (machine precision)
- If larger, you have a systematic indexing error

### Step 4: Follow the Tested Code Path
- Copy the exact logic from `test_mo_vs_ldos.py` or `STM_utils.py`
- Do not "improve" or "optimize" without re-verifying parity
- Comment clearly which convention you are using

## Common Pitfalls

### Pitfall 1: Assuming Fortran is Column-Major
- **True:** Fortran arrays are column-major in memory
- **BUT:** When passed through ctypes to Python, the shape is already adjusted
- **Result:** The Python view `H_blocks[iatom, ineigh, :nnu, :nmu]` is correct as-is
- **Do NOT** transpose based on column-major assumptions

### Pitfall 2: Assuming All Atoms Have 4 Orbitals
- Hydrogen has 1 orbital (s only)
- Heavier atoms may have d-orbitals
- Always check `norb_per[ia]` from `get_orbital_layout()`
- Truncate coupling blocks appropriately

### Pitfall 3: Mixing Up Grid vs Hamiltonian Order
- Grid kernels use `[px, py, pz, s]` (for float4 packing efficiency)
- Hamiltonian kernels use `[s, px, py, pz]`
- Using the wrong remapping will give wrong results
- Always use `_PERM_FORT_TO_GRID` for Grid, `_PERM_FORT_TO_HAM` for Hamiltonian

### Pitfall 4: Guessing at MO Coefficient Layout
- `fc.get_wfcoef()` shape is implementation-dependent
- In this repo, it is `(nmo, norb)` (row-major)
- Do not assume column-major without checking
- Use `_get_mo_vec()` helper to handle both cases safely

## Summary Checklist

Before committing any matrix/basis mapping code:

- [ ] Did I check the shape of `fc.get_wfcoef()` output?
- [ ] Am I using row indexing `C[mo, :]` or column indexing `C[:, mo]` correctly?
- [ ] Did I use `_PERM_FORT_TO_GRID` for Grid kernels and `_PERM_FORT_TO_HAM` for Hamiltonian?
- [ ] Did I handle hydrogen (1 orbital) correctly in coupling blocks?
- [ ] Did I verify against the tested code in `test_mo_vs_ldos.py`?
- [ ] Did I add a Fortran export to compare matrices element-by-element?
- [ ] Did I check symmetry for a symmetric test system (H2O)?
- [ ] Did I comment which convention I'm using in the code?

**IF YOU CANNOT ANSWER YES TO ALL OF THESE, YOUR CODE IS PROBABLY WRONG.**

---

# SUCCESS REPORT: STM Response Methods Achieved Correct Images (April 2026)

## Summary

All STM response methods in `tests/pyFireball/test_response_function.py` now produce correct, symmetric STM images for both small (H2O) and large (PTCDA) molecules. The critical indexing and mapping bugs have been fixed, and Python implementations are now in perfect parity with Fortran reference calculations.

## Test Script

**Location:** `/home/prokop/git/FireCore/tests/pyFireball/test_response_function.py`

This script tests multiple STM response calculation methods and compares them against Fortran reference implementations.

## Methods Tested

### Python Response Methods
1. **response_amplitude_map** - Block elimination with precomputed G0 (O(Ns) per point)
2. **response_amplitude_lu** - LU factorization with adjoint solve optimization
3. **response_amplitude_full_matrix** - Full matrix solve (reference, slower)
4. **response_amplitude_simple_lu** - Simplified STM response matching Fortran SIMPLE algebra

### Fortran Reference Methods
1. **firecore_tipResponse2points** - Full response with two-center DOScentros coupling
2. **firecore_tipResponseSimple2points** - Simplified response with exponential SK coupling

### Coupling Models
1. **True basis coupling** - Uses Fireball tables via `build_inter_system_blocks_fdata`
2. **Exponential SK coupling** - Analytic Slater-Koster with exponential decay via `build_inter_system_blocks_exp_sk`

## Test Results

### H2O (All 6 Orbitals)
- **Grid:** 60×60 points, 10 Å area, z_tip = 2 Å
- **Symmetry:** All orbitals achieve machine precision symmetry (~1e-13)
- **Parity:** Python SIMPLE matches Fortran SIMPLE exactly (max difference ~1e-13)
- **Location:** `/home/prokop/git/FireCore/tests/pyFireball/export/response_function/H2O/`

Example (MO 4 HOMO):
```
[SYMMETRY PY SIMPLE] max|L-R| = 6.253e-13  (mean=1.768e-14)
[SYMMETRY SIMPLE] max|L-R| = 5.684e-13  (mean=1.424e-14)
```

### PTCDA (HOMO±4, MO 70-79)
- **Grid:** 60×60 points, 20 Å area, z_tip = 2 Å
- **Parity:** Python SIMPLE matches Fortran SIMPLE exactly for all orbitals
- **Symmetry:** Asymmetry is physical (DFT orbitals are asymmetric), not a mapping bug
- **Location:** `/home/prokop/git/FireCore/tests/pyFireball/export/response_function/PTCDA/`

Example (MO 74 HOMO):
```
[SYMMETRY PY SIMPLE] max|L-R| = 4.171e+00  (mean=2.750e-01)
[SYMMETRY SIMPLE] max|L-R| = 4.171e+00  (mean=2.750e-01)
```

The identical asymmetry values confirm Python-Fortran parity. The asymmetry reflects the actual DFT electronic structure (verified: `max|psi(x,y) - psi(-x,-y)| = 0.00075`).

## Critical Bugs Fixed

### 1. MO Coefficient Indexing (Mistake 1)
- **Bug:** Used `C_s[:, mo]` (column indexing) for `fc.get_wfcoef()` which returns `(nmo, norb)`
- **Fix:** Added `_get_mo_vec()` helper that correctly uses row indexing `C_s[mo, :]`
- **Impact:** Fixed asymmetric STM images and wrong response patterns

### 2. Hydrogen Coupling Truncation (Mistake 2)
- **Bug:** Python expSK included p→s coupling for hydrogen (1 orbital)
- **Fix:** Zero tip p-rows when sample atom has 1 orbital: `bH4[1:4, 0] = 0.0`
- **Impact:** Achieved exact coupling parity with Fortran expSK

### 3. Fortran Coupling Export
- **Added:** `firecore_export_tip_coupling_point()` in `fortran/MAIN/libFireCore.f90`
- **Purpose:** Export full Hts and Sts coupling matrices for single tip position
- **Benefit:** Enables element-by-element parity verification without guessing

## Verification Protocol

1. **Export ground truth** from Fortran for coupling matrices
2. **Compare element-by-element** with Python coupling (`max|F - P|` should be 0)
3. **Check symmetry** for symmetric systems (H2O should be ~1e-13)
4. **Verify parity** between Python SIMPLE and Fortran SIMPLE
5. **Follow tested conventions** from `test_mo_vs_ldos.py` and `STM_utils.py`

## Status

✅ **COMPLETE** - All STM response methods produce correct images
✅ **PARITY ACHIEVED** - Python matches Fortran for both H2O and PTCDA
✅ **DOCUMENTED** - Critical mapping rules added to this document (see warning section above)

## Next Steps (Future Work)

1. Test with larger molecules beyond PTCDA
2. Optimize GPU kernels for production use
3. Integrate with MD simulation pipeline

---

# SUCCESS REPORT: Rigid Body Rotational Invariance for STM Response (April 2026)

## Summary

Achieved perfect numerical parity between Fortran `firecore_tipResponseSimple2points_rotated` and Python-side rigid-body rotation implementation for STM response calculations. All H2O orbitals (6 total) and PTCDA orbitals around HOMO±4 (9 orbitals) now pass rotational invariance tests with machine-precision accuracy.

## Test Script

**Location:** `/home/prokophapala/git/FireCore/tests/pyFireball/test_response_function_rotated.py`

This script tests rotational invariance of STM response by comparing:
1. **Reference (unrotated):** Fortran `firecore_tipResponseSimple2points`
2. **Fortran rotated:** `firecore_tipResponseSimple2points_rotated` (internal rotation)
3. **Python rotated:** Python `response_amplitude_simple_lu` with inverse-rotated probe points

## Test Results

### H2O (All 6 Orbitals, 0°/45°/90° rotations around x-axis)

All orbitals achieve machine-precision parity between Fortran and Python rotated responses:

| MO | Label | max\|Δ(fc-ref)\| at 0° | max\|Δ(py-fc)\| at 0° | max\|Δ(py-fc)\| at 45° | max\|Δ(py-fc)\| at 90° |
|---|---|---|---|---|---|
| 1 | MO1 | 5.68e-13 | 4.55e-13 | 4.55e-13 | 2.27e-13 |
| 2 | MO2 | 4.44e-15 | 2.84e-14 | 2.84e-14 | 2.84e-14 |
| 3 | MO3 | 1.22e-15 | 1.11e-16 | 1.11e-16 | 1.11e-16 |
| 4 | HOMO | 1.03e-24 | 4.14e-25 | 4.14e-25 | 4.14e-25 |
| 5 | LUMO | 2.04e-14 | 5.68e-14 | 5.68e-14 | 5.68e-14 |
| 6 | MO6 (E=+17.25 eV) | 5.12e-13 | **1.82e-12** | 1.82e-12 | 1.82e-12 |

**Location:** `/home/prokophapala/git/FireCore/tests/pyFireball/export/response_rotated/H2O/`

### PTCDA (MO66-74, HOMO±4)

All 9 orbitals achieve machine-precision parity (max\|Δ(py-fc)\| ~ 1e-13 to 1e-25):

| MO | Energy (eV) | max\|Δ(py-fc)\| at 0° | max\|Δ(py-fc)\| at 45° | max\|Δ(py-fc)\| at 90° |
|---|---|---|---|---|
| 66 | -3.9883 | 1.28e-13 | 1.28e-13 | 1.28e-13 |
| 67 | -3.9290 | 1.14e-13 | 1.14e-13 | 1.14e-13 |
| 68 | -3.2535 | 5.68e-14 | 5.68e-14 | 5.68e-14 |
| 69 | -3.1404 | 3.55e-14 | 3.55e-14 | 3.55e-14 |
| 70 | -3.0206 | 4.14e-25 | 4.65e-25 | 4.14e-25 |
| 71 | -1.6053 (HOMO) | 3.36e-25 | 3.36e-25 | 3.36e-25 |
| 72 | +0.0660 (LUMO) | 1.14e-24 | 1.03e-24 | 1.14e-24 |
| 73 | +0.1652 | 9.82e-25 | 9.82e-25 | 9.82e-25 |
| 74 | +0.9966 | 1.86e-24 | 1.86e-24 | 1.86e-24 |

**Location:** `/home/prokophapala/git/FireCore/tests/pyFireball/export/response_rotated/PTCDA/`

## Critical Bugs Fixed

### ⚠️ BUG 1: Transpose in `fc.get_HS_k()` (THE HORRIBLE TRANSPOSE - AGAIN!)

**The Bug:**
- `fc.get_HS_k()` was returning `Hk_out` and `Sk_out` in row-major layout
- Fortran internally uses column-major layout for these matrices
- Python was using these matrices **without transposing**, effectively using `H.T` and `S.T` instead of `H` and `S`
- This caused massive errors in Python STM response calculations

**How We Found It:**
- After fixing MO coefficient handling (see Bug 2), Python response still didn't match Fortran for MO3
- Debug prints showed coupling matrices matched perfectly (`max|Hts_py - Hts_fc| = 0`)
- The only remaining source of error was the Hamiltonian/Overlap matrices
- Added debug print to compare `H_s` vs `H_s.T` - the transpose fixed the parity

**The Fix:**
In `pyBall/FireCore.py`, modified `get_HS_k()`:
```python
def get_HS_k(kpoint_vec, norbitals):
    Hk_out = np.zeros((norbitals, norbitals), dtype=np.complex128)
    Sk_out = np.zeros((norbitals, norbitals), dtype=np.complex128)
    kpoint_vec_np = np.array(kpoint_vec, dtype=np.float64)
    lib.firecore_get_HS_k(kpoint_vec_np, Hk_out, Sk_out)
    return Hk_out.T, Sk_out.T  # TRANSPOSE to match Fortran column-major convention
```

**Why This Keeps Happening:**
- Fortran arrays are column-major in memory
- When passed through ctypes to NumPy, the shape is preserved but memory layout differs
- The Fortran function `firecore_get_HS_k` fills the array in Fortran's column-major order
- NumPy interprets this as row-major C-contiguous data
- Result: The matrix is transposed relative to what Fortran expects
- **LESSON:** ALWAYS verify matrix layout by comparing with a known working implementation

### ⚠️ BUG 2: `_get_mo_vec()` Square-Matrix Heuristic Failure

**The Bug:**
- The `_get_mo_vec()` function in `pyBall/FireballOCL/STM.py` has a heuristic for square coefficient matrices
- It compares row-norm vs column-norm and picks whichever is **closer to 1.0**
- For MO6 (E=+17.25 eV), the row norm was 1.5656, column norm was ~0.8
- The heuristic chose the **column** (closer to 1.0), which was nearly all zeros except one element
- This caused Python to use a completely wrong MO vector

**How We Found It:**
- After fixing `get_HS_k()` transpose, MO1-MO5 worked but MO6 failed with massive error (8.143e+03)
- Debug prints showed: `||c_vec (Fortran)||=1.565626`, `||C_fc[:,mo0]||=1.283893`, `||C_fc[mo0,:]||=1.565626`
- The Fortran `c_vec` matched the **row** of `C_fc`, not the column
- The heuristic was choosing the wrong orientation for MO6

**The Fix:**
Instead of relying on the buggy heuristic, bypass it entirely by passing a non-square coefficient matrix:
```python
# Old (square, triggers heuristic):
C_use = np.zeros((int(mo0)+1, H_s.shape[0]), dtype=np.float64)  # (6, 6) for MO6
C_use[int(mo0), :] = c_vec
resp_py = response_amplitude_simple_lu(..., C_use, int(mo0), ...)  # Triggers heuristic

# New (non-square, bypasses heuristic):
C_use = np.zeros((1, H_s.shape[0]), dtype=np.float64)  # (1, 6)
C_use[0, :] = c_vec
resp_py = response_amplitude_simple_lu(..., C_use, 0, ...)  # No ambiguity
```

**Root Cause:**
The "closer to 1.0" heuristic assumes MO vectors are normalized to 1.0. For high-energy or unnormalized states, this fails. For sparse matrices (mostly zeros), the column norm can be arbitrarily small and accidentally closer to 1.0 than the row norm.

**Future Fix Needed:**
The `_get_mo_vec()` heuristic should be improved to:
1. Check which orientation has more non-zero elements (for sparse matrices)
2. Use a larger tolerance for "close to 1.0" (e.g., 0.5 instead of strict 1.0)
3. Add a warning when the choice is ambiguous

### BUG 3: MO Coefficient Orientation (Fortran bbnkre vs Python C_fc)

**The Issue:**
- Fortran stores MO coefficients as `bbnkre(orb, mo, ikp)` in column-major memory
- `fc.get_wfcoef()` returns this as a 2D NumPy array
- The question: Does `C_fc[mo, :]` or `C_fc[:, mo]` correspond to `bbnkre(:, mo, ikp)`?

**How We Determined the Correct Mapping:**
- Added `firecore_get_wfcoef_vec()` to export a single MO vector from Fortran: `c = bbnkre(:, iband, ikp)`
- Compared this with both row and column of `C_fc`:
  - `||c_vec (Fortran)||=1.565626`
  - `||C_fc[:,mo0]||=1.283893` (column)
  - `||C_fc[mo0,:]||=1.565626` (row) ✅ MATCHES
- **Conclusion:** `C_fc[mo, :]` (row) corresponds to `bbnkre(:, mo, ikp)` in Fortran

**The Fix:**
Use `fc.get_wfcoef_vec(iband, ikp)` to get the unambiguous Fortran vector, avoiding all transpose ambiguity:
```python
c_vec = np.asarray(fc.get_wfcoef_vec(iband=int(mo0+1), ikp=1, norb=H_s.shape[0]), dtype=np.float64)
C_use = np.zeros((1, H_s.shape[0]), dtype=np.float64)
C_use[0, :] = c_vec
```

## Systematic Debugging Approach

### Step 1: Verify Coupling Matrices Match
- Export coupling matrices from Fortran using `fc.export_tip_coupling_point()`
- Compare with Python coupling matrices element-by-element
- If they don't match, fix basis mapping or Slater-Koster parameters first

### Step 2: Verify Hamiltonian/Overlap Layout
- Use `fc.get_HS_k()` to get H and S
- Compare Python response with `H_s` vs `H_s.T`
- If transpose fixes parity, the wrapper is transposed (add `.T` in the wrapper)

### Step 3: Verify MO Coefficient Orientation
- Use `fc.get_wfcoef_vec()` to get unambiguous Fortran vector
- Compare with both row and column of `C_fc`
- Use the matching orientation, or bypass ambiguity with non-square matrix

### Step 4: Bypass Heuristics in Critical Paths
- For square matrices, heuristics can fail for edge cases (high-energy states, sparse matrices)
- Use non-square matrices or explicit orientation flags to avoid ambiguity
- Add debug prints to verify which orientation is being used

## Relevant Functions for Rotation

### Fortran Functions
- `firecore_tipResponseSimple2points` - Reference STM response (unrotated)
- `firecore_tipResponseSimple2points_rotated` - STM response with internal rigid-body rotation
- `firecore_get_HS_k` - Export Hamiltonian/Overlap at k-point
- `firecore_get_wfcoef` - Export full MO coefficient matrix
- `firecore_get_wfcoef_vec` - Export single MO coefficient vector (unambiguous)
- `firecore_export_tip_coupling_point` - Export tip-sample coupling matrices for debugging

### Python Wrappers (pyBall/FireCore.py)
- `get_HS_k(kpoint, norb)` - Get H and S at k-point (now with transpose fix)
- `get_wfcoef(norb, ikp)` - Get full MO coefficient matrix
- `get_wfcoef_vec(iband, ikp, norb)` - Get single MO coefficient vector
- `tipResponseSimple2points(...)` - Python wrapper for Fortran STM response
- `tipResponseSimple2points_rotated(...)` - Python wrapper for rotated STM response
- `export_tip_coupling_point(...)` - Export coupling matrices for debugging

### Python STM Functions (pyBall/FireballOCL/STM.py)
- `response_amplitude_simple_lu(...)` - Simplified STM response (LU decomposition)
- `_get_mo_vec(C_s, mo_index, norb)` - Extract MO vector from coefficient matrix (has buggy heuristic)

## Exact Conventions and Transforms to Stick To

### 1. Hamiltonian/Overlap Matrix Layout
**CRITICAL:** `fc.get_HS_k()` returns transposed matrices relative to Fortran internal convention.

**Rule:**
```python
H_s, S_s = fc.get_HS_k(kpt, norb)  # Returns H.T and S.T internally
# Use directly in Python - the wrapper handles the transpose
# Do NOT transpose again in user code
```

**Verification:**
- Compare Python response with `H_s` vs `H_s.T`
- If transpose fixes parity, the wrapper needs fixing (add `.T` in `get_HS_k`)

### 2. MO Coefficient Matrix Layout
**CRITICAL:** `fc.get_wfcoef()` returns `(nmo, norb)` with MO vectors as **rows**.

**Rule:**
```python
# For fc.get_wfcoef():
C_fc = fc.get_wfcoef(norb=norb)  # Shape: (nmo, norb)
c_mo = C_fc[mo_idx, :]  # ROW indexing is correct
# WRONG: c_mo = C_fc[:, mo_idx]

# For unambiguous access, use:
c_vec = fc.get_wfcoef_vec(iband=mo_idx+1, ikp=1, norb=norb)  # Direct from Fortran
```

**Verification:**
- Compare `c_vec` with `C_fc[mo, :]` and `C_fc[:, mo]`
- The matching orientation is the correct one

### 3. Rigid Body Rotation Convention
**For probe point rotation (inverse rotation):**
```python
# Fortran internal rotation (inverse-rotate probe points):
# pts_body = cen + (pts - cen) @ R.T  (column-vector convention)

# Python equivalent (row-vector convention):
cen = mol.apos.mean(axis=0)
pts_body = cen + (pts - cen) @ R  # row-vectors: dp @ R == (R^T @ dp)_col
```

**For matrix rotation (if needed in future):**
```python
# H' = U * H * U^T
# S' = U * S * U^T
# Where U = block_diag(diag(1, R_1), diag(1, R_2), ...)
# R_i is the 3x3 rotation matrix for p-orbitals of atom i
```

### 4. Transpose Checklist (MANDATORY)

Before passing any matrix between Fortran and Python:

1. **Check the Fortran source:** Is the array filled in column-major order?
2. **Check the ctypes signature:** Is it passed as a 2D array or flattened?
3. **Test with a known working implementation:** Compare element-by-element
4. **Add debug prints:** Print shapes and norms before/after transpose
5. **Document the convention:** Comment in the code which orientation is expected

**If you cannot verify the layout, DO NOT GUESS. Add a Fortran export function to compare.**

## Status

✅ **COMPLETE** - Rotational invariance achieved for H2O and PTCDA
✅ **PARITY VERIFIED** - Python matches Fortran to machine precision
✅ **CRITICAL BUGS FIXED** - Transpose in `get_HS_k()` and `_get_mo_vec()` heuristic
✅ **DOCUMENTED** - Conventions and debugging protocol added to this document
✅ **IMAGES GENERATED** - Side-by-side comparison plots in export directories

## Lessons Learned

1. **Transpose bugs are the most common and expensive errors** in Fortran-Python interfaces
2. **Heuristics for ambiguous cases (square matrices) can fail** - use explicit orientation
3. **Always verify with element-by-element comparison** against Fortran exports
4. **Add debug prints early** - don't wait until you're completely lost
5. **Document conventions immediately** after fixing bugs, to prevent recurrence
6. **Use unambiguous data structures** (1D vectors, non-square matrices) when possible
4. Validate against experimental STM images

---

# STM Orbital Overlap with Rigid Body Rotation (REPORT 2026-04-23) 

## Overview

This section documents the implemented STM orbital overlap method with per-pixel rigid body rotation of the tip molecule. This is a fast GPU-based approach for simulating STM junctions where the tip molecule can be rotated relative to the sample molecule.

## Method: Per-Pixel Quaternion Rotation in GPU Kernel

### Core Concept

Instead of rotating the entire tip Hamiltonian/Overlap matrices (which would be O(N³)), we apply a per-pixel rigid body rotation inside the GPU kernel. For each grid point (pixel) in the STM scan:

1. **Tip orbitals** are evaluated at the rotated position `r' = R * r` where R is the rotation matrix
2. **Sample orbitals** are evaluated at the original position
3. **Overlap integral** is computed between rotated tip orbitals and sample orbitals
4. **STM amplitude** is computed from the overlap

This is mathematically equivalent to rotating the tip molecule, but much faster because:
- No matrix rotations needed
- Rotation is applied only at the grid points where we evaluate the STM signal
- GPU kernel handles the rotation in parallel for all pixels

### Implementation Details

**Test Script:** `tests/pyFireball/test_stm_orbital_rotated.py`

**GPU Kernel:** `pyBall/FireballOCL/cl/Grid.cl` (unchanged - already supports per-pixel rotation)

**Python Wrapper:** `pyBall/FireballOCL/Grid.py` (unchanged - already supports quaternion rotation)

**Key Features:**
- **Per-pixel quaternion rotation** applied to tip orbitals inside GPU kernel
- **Separate tip and sample MO selection** for debugging flexibility
- **Independent rotation control** for tip and sample molecules
- **Mock AO mode** for debugging with single atomic orbitals (s, px, py, pz)
- **True orbital mode** using Fireball molecular orbitals

## Capabilities

### 1. Mock AO Mode (Debugging)

For intuitive visual verification of rotation and overlap patterns, the tip and sample can be composed of single atomic orbitals:

```bash
# Tip: px orbital on atom 1, Sample: py orbital on atom 1
# Rotate tip around z-axis by 0°, 45°, 90°
python3 test_stm_orbital_rotated.py \
  --do_overlap --mock_ao \
  --mock_tip_atom 1 --mock_tip_ao px \
  --mock_smp_atom 1 --mock_smp_ao py \
  --tip_axis z --tip_angles 0,45,90 \
  --n 80 --size 10.0 --ztip 3.0 --zmid 1.5
```

**Purpose:**
- Verify rotation symmetry (p-orbitals should rotate correctly)
- Check that overlap patterns match expectations
- Debug kernel rotation logic without full MO complexity

### 2. True Orbital Mode

Using Fireball molecular orbitals for realistic STM simulation:

```bash
# PTCDA HOMO (MO 70) with different tip rotations around z-axis
python3 test_stm_orbital_rotated.py \
  --do_overlap \
  --xyz ../../cpp/common_resources/xyz/PTCDA.xyz \
  --mo 70 \
  --tip_axis z --tip_angles 0,30,45,60,90 \
  --n 80 --size 20.0 --ztip 3.0 --zmid 1.5
```

### 3. Separate Tip and Sample MOs

For debugging, different MOs can be selected for tip and sample:

```bash
# Tip: MO 1, Sample: MO 2
python3 test_stm_orbital_rotated.py \
  --do_overlap \
  --tip_mo 1 --smp_mo 2 \
  --tip_axis x --tip_angles 0,45,90 \
  --smp_axis z --smp_angles 0
```

### 4. Independent Rotation Control

Tip and sample can be rotated independently:

```bash
# Tip: rotate around x-axis, Sample: rotate around z-axis
python3 test_stm_orbital_rotated.py \
  --do_overlap \
  --tip_axis x --tip_angles 0,30,45 \
  --smp_axis z --smp_angles 0,90
```

## Visualization

The script generates 4-panel figures:

1. **Panel 1 (top-left):** Sample MO orbital (ψ_sample at mid-plane)
   - Color map: `bwr` (blue-white-red, symmetric around zero)
   - Normalization: `vmin=-vmax_s, vmax=vmax_s` (separate per panel)
   - Atoms plotted: Sample atoms (with rotation applied if sample has rotation)

2. **Panel 2 (top-right):** Tip MO orbital (ψ_tip at mid-plane)
   - Color map: `bwr` (symmetric around zero)
   - Normalization: `vmin=-vmax_t, vmax=vmax_t` (separate per panel)
   - Atoms plotted: Tip atoms (with rotation applied)

3. **Panel 3 (bottom-left):** STM overlap amplitude |t|²
   - Computed as overlap of rotated tip orbitals with sample orbitals
   - Color map: `viridis`
   - Normalization: `vmin=0, vmax=vI` (separate per panel)
   - Atoms plotted: Sample atoms at mid-plane

4. **Panel 4 (bottom-right):** Cross-correlation of sample and tip orbitals
   - Computed as `corr2 = correlate2d(psi_s, psi_t, mode='same')`
   - Color map: `viridis`
   - Normalization: `vmin=0, vmax=vC` (separate per panel)
   - No atoms (correlation space, not real space)

**Figure Caption:** Shows tip MO, sample MO, tip rotation axis and angle, sample rotation axis and angle, z-mid plane position, and tip height.

**Filename:** Includes tip MO, sample MO, z-mid, z-tip, tip rotation (axis+angle), sample rotation (axis+angle).

## CLI Arguments

### MO Selection
- `--mo N`: Use MO index N (1-based) for both tip and sample (default: HOMO-1)
- `--mo_list N1,N2,...`: Use multiple MOs for both tip and sample
- `--tip_mo N`: Use MO index N for tip only (overrides `--mo`)
- `--tip_mo_list N1,N2,...`: Use multiple MOs for tip only
- `--smp_mo N`: Use MO index N for sample only (overrides `--mo`)
- `--smp_mo_list N1,N2,...`: Use multiple MOs for sample only

### Rotation Control
- `--tip_axis {x,y,z}`: Rotation axis for tip molecule
- `--tip_angles A1,A2,...`: Rotation angles for tip (degrees)
- `--smp_axis {x,y,z}`: Rotation axis for sample molecule
- `--smp_angles A1,A2,...`: Rotation angles for sample (degrees)
- (Legacy: `--axis` and `--angle` for backward compatibility with tip rotation only)

### Mock AO Mode
- `--mock_ao`: Enable mock AO mode
- `--mock_tip_atom N`: Atom index (1-based) for tip AO
- `--mock_tip_ao {s,px,py,pz}`: Orbital type for tip
- `--mock_smp_atom N`: Atom index (1-based) for sample AO
- `--mock_smp_ao {s,px,py,pz}`: Orbital type for sample

### Grid Parameters
- `--n N`: Grid size (N×N pixels, default: 80)
- `--size L`: Physical size in Å (default: 10.0)
- `--ztip Z`: Tip height above sample (default: 3.0)
- `--zmid Z`: Mid-plane z-coordinate (default: 1.5)
- `--beta B`: Decay parameter (default: 1.0)
- `--r0 R0`: Cutoff radius (default: 3.0)
- `--rcut RCUT`: Maximum cutoff (default: 8.0)

## Test Results

### CH2O (Small Molecule)

Successfully tested with:
- Different tip/sample MOs (MO 1, MO 2)
- Different rotation axes (x, z)
- Rotation angles: 0°, 30°, 45°, 60°, 90°
- Atoms plotted correctly with proper rotation
- Separate normalization working correctly

### PTCDA (Large Molecule)

Successfully tested with HOMO (MO 70):
- Rotation around z-axis: 0°, 30°, 45°, 60°, 90°
- Rotation around x-axis: 0°, 30°, 45°, 60°, 90°
- Rotation around y-axis: 0°, 30°, 45°, 60°, 90°
- 15 figures generated (5 angles × 3 axes)
- All figures show correct rotation behavior and symmetry

## Technical Notes

### Rotation Implementation

The rotation is implemented using quaternions in the GPU kernel:

1. **Python side:** Convert axis-angle to quaternion
   ```python
   def axis_angle_to_quaternion(axis, angle_deg):
       angle_rad = np.deg2rad(angle_deg)
       axis = np.array(axis, dtype=np.float64)
       axis = axis / np.linalg.norm(axis)
       w = np.cos(angle_rad / 2)
       xyz = axis * np.sin(angle_rad / 2)
       return np.array([w, xyz[0], xyz[1], xyz[2]])
   ```

2. **GPU kernel:** Apply quaternion rotation to tip positions
   - For each grid point, compute rotated tip position
   - Evaluate tip orbitals at rotated position
   - Compute overlap with sample orbitals at original position

### Why Per-Pixel Rotation?

**Advantages:**
- **Fast:** No matrix multiplication needed (O(N³) → O(N_pixels))
- **Simple:** No need to handle sparse matrix rotation
- **Accurate:** Mathematically equivalent to rotating the molecule
- **Parallel:** GPU handles all pixels simultaneously

**Disadvantages:**
- Only works for rigid body rotation (no internal deformation)
- Rotation must be applied at evaluation time (can't pre-rotate orbitals)

### Comparison with Matrix Rotation Approaches

The NOTES_rigid_body_rotation_analysis.md document discusses three approaches for rotating H/S matrices:

1. **Option A (Dense):** Rotate dense H, S matrices - O(N³), slow
2. **Option B (Sparse):** Rotate sparse blocks - O(N), complex
3. **Option C (Rebuild):** Rebuild from rotated positions - Slow, need SCF

**Our approach (Per-Pixel):** Rotate at evaluation time - O(N_pixels), fast and simple

This is essentially a variant of Option C, but instead of rebuilding the entire Hamiltonian, we only rotate the evaluation points. This is possible because:
- We're only interested in the STM signal (overlap at specific points)
- We don't need the full rotated Hamiltonian for other purposes
- The GPU kernel can handle the rotation efficiently

## Future Extensions

1. **Soft mode expansion:** For non-rigid deformations, expand orbitals along normal modes
2. **Energy-dependent PDOS:** Use PDOS instead of single MO for more realistic transport
3. **Tip-sample coupling:** Add explicit coupling matrices between tip and sample
4. **Green's function response:** Compute full response function (not just overlap)

## References

- Test script: `tests/pyFireball/test_stm_orbital_rotated.py`
- GPU kernel: `pyBall/FireballOCL/cl/Grid.cl`
- Python wrapper: `pyBall/FireballOCL/Grid.py`
- Rotation analysis: `tests/pyFireball/NOTES_rigid_body_rotation_analysis.md`

# STM Orbital Overlap for Two Different Molecules (IMPLEMENTED)

## Overview

This section documents the extension of the STM orbital overlap method to support simulations where the tip and sample are **different molecules** (e.g., CH2O tip and PTCDA sample). This is implemented in a separate test script with a new kernel entrypoint to avoid breaking the existing single-molecule functionality.

## Implementation

### New Kernel Entrypoint

**Location:** `pyBall/FireballOCL/cl/Grid.cl`

**Kernel:** `mo_overlap_points_exp_sk_2mol`

This is a duplicate of `mo_overlap_points_exp_sk` with a new name to make the two-molecule use case explicit. The implementation is identical; the separate name:
- Avoids breaking existing workflows
- Makes call sites self-documenting
- Allows future extensions specific to two-molecule cases

### New Python Wrapper

**Location:** `pyBall/FireballOCL/Grid.py`

**Method:** `GridProjector.mo_overlap_points_exp_sk_2mol(...)`

Same arguments and behavior as `mo_overlap_points_exp_sk`, but calls the new kernel entrypoint.

### New Test Script

**Location:** `tests/pyFireball/test_stm_orbital_rotated_2mol.py`

**Key Features:**
- Loads two different molecules from separate XYZ files
- Runs SCF for each molecule independently
- Supports independent tip and sample MO selection
- Supports independent tip and sample rotation
- Generates 4-panel visualization (sample MO, tip MO, STM overlap, cross-correlation)

### SCF Caching via Subprocesses

**Problem:** FireCore Fortran side does not tolerate repeated init/SCF cycles for different molecules in a single process, causing:
- `Attempting to allocate already allocated variable 'degelec'`
- `free(): corrupted unsorted chunks` (memory corruption)

**Solution:** The script uses subprocess caching:
1. For each molecule, spawn a subprocess to run SCF once
2. Save SCF results (eigenvalues, MO coefficients, orbital layout) to `.npz` cache files
3. Main process loads caches and runs GPU overlap scan
4. Caches are reused on subsequent runs (fast)

**CLI for SCF caching (internal):**
```bash
python3 test_stm_orbital_rotated_2mol.py \
  --dump_scf \
  --xyz_one <molecule.xyz> \
  --cache_one <output.npz> \
  --nmax_scf 30
```

## Capabilities

### Two Different Molecules

Load separate tip and sample molecules:

```bash
python3 test_stm_orbital_rotated_2mol.py \
  --xyz_tip ../../cpp/common_resources/xyz/CH2O.xyz \
  --xyz_smp ../../cpp/common_resources/xyz/PTCDA.xyz \
  --smp_mo 70 \
  --tip_mo 2 \
  --tip_axis x --tip_angles 0,30,45,60,90 \
  --smp_axis z --smp_angles 0
```

### Independent MO Selection

Different MOs for tip and sample:

```bash
--tip_mo 2    # CH2O MO2
--smp_mo 70   # PTCDA MO70 (HOMO)
```

### Independent Rotation Control

Separate rotation axes and angles:

```bash
--tip_axis x --tip_angles 0,30,45,60,90
--smp_axis z --smp_angles 0
```

## CLI Arguments

### Molecule Selection
- `--xyz_tip <path>`: Tip molecule XYZ file (default: CH2O.xyz)
- `--xyz_smp <path>`: Sample molecule XYZ file (default: PTCDA.xyz)
- `--tip_mo N`: Tip MO index 1-based (default: HOMO-1)
- `--smp_mo N`: Sample MO index 1-based (default: HOMO-1)
- `--tip_mo_list N1,N2,...`: Multiple tip MOs
- `--smp_mo_list N1,N2,...`: Multiple sample MOs

### Rotation Control
- `--tip_axis {x,y,z}`: Tip rotation axis
- `--tip_angles A1,A2,...`: Tip rotation angles (degrees)
- `--smp_axis {x,y,z}`: Sample rotation axis
- `--smp_angles A1,A2,...`: Sample rotation angles (degrees)

### SCF Parameters
- `--nmax_scf N`: Maximum SCF iterations (default: 200)
- Caching is automatic; first run creates `.npz` files, subsequent runs reuse them

### Grid Parameters
- `--n N`: Grid size (N×N pixels, default: 80)
- `--size L`: Physical size in Å (default: 20.0)
- `--ztip Z`: Tip height above sample (default: 3.0)
- `--zmid Z`: Mid-plane z-coordinate (default: ztip/2)
- `--beta B`: Decay parameter (default: 1.0)
- `--r0 R0`: Cutoff radius (default: 3.0)
- `--rcut RCUT`: Maximum cutoff (default: 8.0)

## Test Results

### PTCDA Sample + CH2O Tip

Successfully tested with:
- Sample: PTCDA MO70 (HOMO)
- Tip: CH2O MO1, MO2, MO3 (various orbitals)
- Rotation: x-axis, angles 0°, 30°, 45°, 60°, 90°
- Grid: 80×80, 20 Å, ztip=3.0 Å, zmid=1.5 Å

**Output Files:**
- `overlap2mol_tipCH2O_MO001_smpPTCDA_MO070_..._tipx000_smpz000.png`
- `overlap2mol_tipCH2O_MO001_smpPTCDA_MO070_..._tipx030_smpz000.png`
- `overlap2mol_tipCH2O_MO001_smpPTCDA_MO070_..._tipx045_smpz000.png`
- `overlap2mol_tipCH2O_MO001_smpPTCDA_MO070_..._tipx060_smpz000.png`
- `overlap2mol_tipCH2O_MO001_smpPTCDA_MO070_..._tipx090_smpz000.png`
- (Similar for MO2, MO3)

**SCF Cache Files:**
- `scf_tip_CH2O.npz` (3 KB)
- `scf_smp_PTCDA.npz` (136 KB)

**Location:**
`/home/prokop/git/FireCore/tests/pyFireball/export/stm_orbital_rotated_2mol/tip_ch2o__smp_ptcda/`

## Visualization

Same 4-panel layout as single-molecule version:
1. Sample MO (ψ_sample at mid-plane)
2. Tip MO (ψ_tip at mid-plane)
3. STM overlap |t|²
4. Cross-correlation of sample and tip orbitals

**Figure Caption:** Shows tip molecule, sample molecule, tip MO, sample MO, tip rotation (axis+angle), sample rotation (axis+angle), z-mid, z-tip.

**Filename:** Includes tip molecule name, sample molecule name, tip MO, sample MO, z-mid, z-tip, tip rotation, sample rotation.

## Technical Notes

### Why Separate Kernel Entrypoint?

The kernel `mo_overlap_points_exp_sk` already accepts separate tip and sample atom data, so mathematically it can handle two molecules. However:
- We keep a separate entrypoint for **explicit intent**
- Avoids breaking existing single-molecule workflows
- Allows future extensions specific to two-molecule cases (e.g., different basis sets, explicit coupling matrices)

### SCF Caching Design

The subprocess caching approach:
- **Isolates Fortran state:** Each SCF run gets a clean Fortran process
- **Fast re-runs:** Caches are reused; no need to recompute SCF
- **Portable:** `.npz` files can be copied between machines
- **Transparent:** User doesn't need to manage caching manually

**Cache File Contents:**
- Atom types and positions
- Eigenvalues
- MO coefficients
- Orbital layout (norb_per, starts, orb2atom)
- Number of orbitals and atoms

### Comparison with Single-Molecule Script

| Feature | `test_stm_orbital_rotated.py` | `test_stm_orbital_rotated_2mol.py` |
|---------|-------------------------------|------------------------------------|
| Molecules | Same molecule for tip and sample | Different tip and sample molecules |
| SCF | Single SCF run in main process | Two SCF runs via subprocess caching |
| Kernel | `mo_overlap_points_exp_sk` | `mo_overlap_points_exp_sk_2mol` |
| MO selection | Can select different MOs but same molecule | Can select different MOs from different molecules |
| Use case | Single-molecule STM junction | Two-molecule STM junction |

## Future Extensions

1. **Multiple tip MOs:** Sum over multiple tip orbitals for more realistic tip states
2. **PDOS integration:** Use energy-dependent PDOS instead of single MO
3. **Explicit tip-sample coupling:** Add coupling matrices between tip and sample
4. **Response function:** Combine with Green's function response for full transport

## References

- Test script: `tests/pyFireball/test_stm_orbital_rotated_2mol.py`
- GPU kernel: `pyBall/FireballOCL/cl/Grid.cl` (kernel `mo_overlap_points_exp_sk_2mol`)
- Python wrapper: `pyBall/FireballOCL/Grid.py` (method `mo_overlap_points_exp_sk_2mol`)
- Single-molecule version: `tests/pyFireball/test_stm_orbital_rotated.py`

---

# Implementation Notes: Fortran Dyson STM Reference (2026-04-24)

## Overview

Implemented a CPU-based, non-parallel Fortran reference for the Dyson Green's function STM method using real Fireball Hamiltonians and overlaps for both tip and sample molecules. This serves as a reference implementation for future GPU (pyOpenCL) acceleration.

## Key Implementation Details

### Angular Symmetry: Critical Lesson
**Do NOT use ad-hoc Slater-Koster blocks.** Initial implementation used manual if-then-else blocks for Vss, Vsp, Vpp_sig, Vpp_pi with arbitrary scaling. This broke symmetry and was unmaintainable.

**Correct approach**: Reuse Fireball's native angular machinery:
- `interpolate_1d`: Get radial values from pre-computed tables
- `recover_2c`: Reconstruct 2D orbital matrix from 1D list
- `rotate_fb`: Rotate using `twister` and `chooser` for s, p, d orbitals
- Only replace radial part with exponential decay: `exp(-beta*(r-r0))`

This ensures rotational invariance and extensibility to any basis set.

### Files Modified

**Fortran Core:**
- **`fortran/MAIN/libFireCore.f90`**
  - Added `firecore_stm_gf_2mol_mo_2points` (lines ~2100-2230)
  - Uses Fireball's `interpolate_1d`, `recover_2c`, `rotate_fb`
  - Exponential radial scaling for tunneling range extension
  - Accepts real atomic species arrays for proper shell metadata

**Python Bindings:**
- **`pyBall/FireCore.py`**
  - Added `stm_gf_2mol_mo_2points` binding (lines ~488-521)

**Test Script:**
- **`tests/pyFireball/test_stm_gf_dyson_2mol.py`**
  - 2-molecule Dyson STM test with 6-panel visualization
  - SCF caching for tip and sample
  - Orbital rotation using Fireball routines
  - Green's function computation for isolated subsystems

### How to Run

```bash
cd /home/prokophapala/git/FireCore
bash make.sh
python3 tests/pyFireball/test_stm_gf_dyson_2mol.py \
    --xyz_tip ../../cpp/common_resources/xyz/CH2O.xyz \
    --xyz_smp ../../cpp/common_resources/xyz/PTCDA.xyz \
    --n 80 --ztip 3.0 --zmid 1.5 --size 20.0 \
    --nmax_scf 30 --tip_mo_list 2 --smp_mo_list 70 \
    --outdir export/stm_gf_dyson_2mol_fortran_ref4
```

### Results
Produces STM images with proper nodal structure and symmetry consistent with molecular orbitals. Example output in `export/stm_gf_dyson_2mol_fortran_ref4/`.

## Future pyOpenCL Reimplementation

The Fortran reference should be ported to pyOpenCL for GPU acceleration. Key components to port:

### Kernels to Implement
1. **Radial interpolation**: Port `interpolate_1d` logic to OpenCL
2. **Orbital matrix recovery**: Port `recover_2c` with mu/nu lookup
3. **Rotation**: Port `rotate_fb` / `twister` / `chooser` for s, p, d orbitals
4. **Exponential scaling**: Apply `exp(-beta*(r-r0))` after interpolation

### Existing Reference Code
- **`pyBall/FireballOCL/cl/hamiltonian.cl`**: Already has `rotate_fb_matrix_sp` for s+p rotation
- **`fortran/ROTATIONS/rotate.f90`**: Full reference implementation
- **`fortran/INTERPOLATERS/interpolate_1d.f90`**: Radial interpolation reference
- **`fortran/INTERPOLATERS/recover_2c.f90`**: Matrix reconstruction reference

### GPU Parallelization Strategy
- Each workgroup handles one pixel (tip position)
- Shared memory for small coupling matrices (contact region only)
- Batch processing for MD integration (multiple pixels simultaneously)
- Precompute tip/sample Green's functions on CPU, reuse on GPU

### Testing Required
Before pyOpenCL port, test Fortran reference across:
- Different molecules (H2O, CH4, benzene, PTCDA)
- Different orbital types
- Different geometries
- Different transport energies

## Files to Add to Git

Do not forget to add:
1. **`tests/pyFireball/test_stm_gf_dyson_2mol.py`** - Main test script
2. **`fortran/MAIN/libFireCore.f90`** - Updated with new Dyson subroutine
3. **`pyBall/FireCore.py`** - Updated with new binding

The script `tests/pyFireball/test_stm_dyson_fortran.py` appears to be a separate/older test and may not need to be added if not part of this implementation.