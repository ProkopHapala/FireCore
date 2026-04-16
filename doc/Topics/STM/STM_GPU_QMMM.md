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

---

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

---

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

---

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

---

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

---

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

---

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

---

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

---

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

---

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

---

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

---

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

---

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

---

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

---

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