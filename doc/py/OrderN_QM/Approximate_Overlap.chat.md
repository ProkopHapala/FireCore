https://aistudio.google.com/app/prompts?state=%7B%22ids%22:%5B%2213BfBgHtNbkazBS-RZIc0wpSbfLnUqEKW%22%5D,%22action%22:%22open%22,%22userId%22:%22100958146796876347936%22,%22resourceKeys%22:%7B%7D%7D&usp=sharing

## USER

Could you explain be better semiempirical methods based on neglect of differential overlap? (https://en.wikipedia.org/wiki/NDDO and similar) I realized that those methods would hugely help me to define fast linear scalling methods, as overlap matrix is often source of numerical  complications. 

https://en.wikipedia.org/wiki/Semi-empirical_quantum_chemistry_method

what about semiemperical quantum chemistry methods? How they differ?



---

You have hit on a crucial insight. The handling of the **Overlap Matrix ($S$)** and the **Two-Electron Integrals (ERIs)** are indeed the two biggest bottlenecks in Hartree-Fock and DFT calculations.

By utilizing approximations based on the **Neglect of Differential Overlap (NDO)**, you effectively assume the basis set is orthonormal ($S \approx I$). This simplifies the computationally expensive generalized eigenvalue problem ($FC = SCE$) into a standard eigenvalue problem ($FC = CE$). More importantly, it drastically reduces the number of electron repulsion integrals from $\sim N^4$ to $\sim N^2$.

Here is a detailed breakdown of the hierarchy of these methods, how they differ, and why NDDO is the "sweet spot" for fast calculations.

---

### 1. The Core Concept: Zero Differential Overlap (ZDO)

In standard Quantum Chemistry, the "Differential Overlap" is the product of two atomic orbitals, $\phi_\mu$ and $\phi_\nu$, at the same position in space for the same electron:
$$ d\Omega = \phi_\mu(r_1) \phi_\nu(r_1) $$

If $\mu$ and $\nu$ are on different atoms, this product is small. The **ZDO approximation** sets this product to zero unless $\mu = \nu$.

**The Impact on Integrals:**
A standard Two-Electron Integral looks like this:
$$ (\mu \nu | \lambda \sigma) = \iint \phi_\mu^*(1)\phi_\nu(1) \frac{1}{r_{12}} \phi_\lambda^*(2)\phi_\sigma(2) \, d\tau_1 \, d\tau_2 $$

If we apply strict ZDO:
1.  **3-Center and 4-Center Integrals vanish.** (These are the most numerous and expensive).
2.  We are left only with **1-Center** and **2-Center** integrals.
3.  The Overlap Matrix $S_{\mu\nu}$ becomes the Kronecker delta $\delta_{\mu\nu}$ (Identity matrix).

---

### 2. The Hierarchy: CNDO vs. INDO vs. NDDO

The difference between these methods lies in **how strict** they are about neglecting the overlap.

#### A. CNDO (Complete Neglect of Differential Overlap)
*   **The Approximation:** The most severe. It assumes that atomic orbitals are spherical spheres of charge. It ignores the directionality of p-orbitals and d-orbitals inside the electron repulsion integrals.
*   **The Rule:** $(\mu \nu | \lambda \sigma)$ is non-zero *only* if $\mu=\nu$ and $\lambda=\sigma$.
*   **Pros:** Extremely fast.
*   **Cons:** It cannot predict bond lengths correctly or distinguish between singlet and triplet states well because it destroys the "Exchange Interaction" (which depends on overlap).
*   **Verdict:** Obsolete.

#### B. INDO (Intermediate Neglect of Differential Overlap)
*   **The Approximation:** Recognizes that electrons on the *same* atom do interact differently depending on which orbital they are in.
*   **The Rule:** Keeps all **1-center** exchange integrals (where $\mu, \nu, \lambda, \sigma$ are all on the same atom).
*   **Pros:** Can calculate spectroscopic states (splitting between singlet and triplet).
*   **Cons:** Still treats atoms as "spheres" when they interact with *other* atoms.

#### C. NDDO (Neglect of Diatomic Differential Overlap)
This is the method you are interested in. It is the foundation of modern semi-empirical chemistry (MNDO, AM1, PM3, PM6, PM7).

*   **The Approximation:** It neglects overlap only when the orbitals $\mu$ and $\nu$ belong to **different atoms**.
*   **The Rule:** It keeps all integrals where the charge cloud $\phi_\mu \phi_\nu$ is on Atom A, and $\phi_\lambda \phi_\sigma$ is on Atom B.
*   **Physical Meaning:** This preserves the multipole moments (dipoles, quadrupoles) of the charge distribution. It knows that a $p_x$ orbital points along the x-axis and will repel another electron differently than a $p_y$ orbital.
*   **Why it works for you:** It maintains the directional nature of chemical bonding (essential for organic chemistry) while still killing the $N^4$ scaling.

---

### 3. Comparison of Modern NDDO-based Methods

Since NDDO provides the mathematical framework, the actual "methods" differ in how they parameterize the remaining integrals against experimental data.

| Method | Full Name | Characteristics |
| :--- | :--- | :--- |
| **MNDO** | Modified NDDO | The "grandfather" of modern codes. Good for simple organics. Fails badly at hydrogen bonding and steric crowding. |
| **AM1** | Austin Model 1 | Modified core-core repulsion functions to fix Hydrogen bonding. Much better for biochemistry, but still has issues with phosphorus/sulfur. |
| **PM3** | Parametric Method 3 | Released by Stewart. Purely mathematical optimization of parameters (less physical intuition, more data fitting). Often better than AM1 for transition states. |
| **PM6 / PM7** | Parametric Method 6/7 | The modern standard. Trained on massive databases (crystal structures, enthalpies). PM7 includes corrections for dispersion (van der Waals) and hydrogen bonding. |
| **RM1** | Recife Model 1 | A re-parameterization of AM1 to fix its errors while keeping the same math. |

---

### 4. Beyond NDDO: The Linear Scaling Perspective

You mentioned that the Overlap Matrix is a source of numerical complications. Here is how these methods solve that for your linear scaling goal:

**1. The Orthogonalization Trick**
In standard DFT, as the system gets large, the overlap matrix $S$ becomes nearly singular (linear dependencies in the basis set), causing numerical instability when inverting or diagonalizing.
NDDO methods perform calculations in a pre-orthogonalized basis. The parameters (experimental data) "absorb" the effects of the non-orthogonality.
*   **Result:** You generally solve $FC = CE$ (Standard Eigenvalue) instead of $FC = SCE$.

**2. Sparse Matrices**
Because 3- and 4-center integrals are zero, and 2-center integrals decay rapidly with distance, the Hamiltonian matrix in NDDO methods is naturally **sparse**.
*   This makes them perfect candidates for **"Divide and Conquer"** or **"Mosaic"** linear scaling algorithms.

**3. The Modern Competitor: DFTB (Density Functional Tight Binding)**
If you are writing code today, you must look at **DFTB**.
*   **Relation to NDDO:** It is mathematically very similar to NDDO but derived from Density Functional Theory rather than Hartree-Fock.
*   **Advantage:** It is generally more robust for materials science and metals than PM6/PM7.
*   **SCC-DFTB:** This version includes Self-Consistent Charge (iterative solving), which fixes charge transfer errors in standard tight-binding.

### Summary Advice for Your Implementation

If your goal is to define a fast, linear-scaling method:

1.  **Stick to NDDO formalism:** It offers the best balance of speed (neglects distant overlap) and accuracy (keeps local directional overlap).
2.  **Look into OMx (Orthogonalization Model x):** This is a newer class of methods that explicitly includes corrections for the fact that the basis is *not* actually orthogonal, offering higher accuracy than PM6 without much extra cost.
3.  **Numerical Stability:** By forcing $S=I$, you eliminate the Cholesky decomposition step usually required to orthogonalize the basis, which is a major speedup and stability boost for large $N$.

The current state-of-the-art for general purpose semi-empirical code is usually **PM7** (in the MOPAC software) or **GFN2-xTB** (a Tight-Binding method by Grimme, which is incredibly popular right now for being fast and robust across the periodic table).

---

## USER

yes, definitely I'm mor interested in SCC-DFTB, but in any case I need more explanation of mathematical background how to get-rid of the overlap matrix. To what basis the system hamiltionan and schrodedinger equation correspond (how can I recover molecular orbitals and density in the real space? is it even rigorously possible)

rather then rely on extensive fittoing of expoeriment or higher levels of theory, I would prefer to use some systematic mathematical method how to get from standard LCAO basis set to SCC-DFT or NDDO with diagonal overlap.

Also I do not undestand clearly what it means "differential overlap" I was thinking we neglect density matrix completely (?)



---

This is a deep question that touches the very heart of why tight-binding and semi-empirical methods work. You are asking about the **rigorous link** between full DFT and these fast approximations.

Let’s dismantle your doubts systematically.

---

### 1. Clarification: Differential Overlap vs. Density Matrix

You wrote: *"I was thinking we neglect density matrix completely (?)"*

**No.** If you neglect the Density Matrix ($P_{\mu\nu}$), you have no electrons and no chemistry.

**Differential Overlap** refers to the physical product of two basis functions in real space:
$$ d\Omega(\vec{r}) = \chi_\mu(\vec{r}) \cdot \chi_\nu(\vec{r}) $$

In the exact theory, electron repulsion integrals require us to calculate the interaction between "cloud" $d\Omega_1$ and "cloud" $d\Omega_2$.
$$ (\mu\nu|\lambda\sigma) = \int \int \frac{\chi_\mu(1)\chi_\nu(1) \cdot \chi_\lambda(2)\chi_\sigma(2)}{r_{12}} d\tau_1 d\tau_2 $$

**The ZDO (Zero Differential Overlap) Approximation** says:
"If $\mu \neq \nu$, the product $\chi_\mu(\vec{r}) \cdot \chi_\nu(\vec{r})$ is everywhere zero."

Mathematically, this imposes:
$$ \chi_\mu \chi_\nu = \delta_{\mu\nu} |\chi_\mu|^2 $$
This kills 3-center and 4-center integrals, but the **Density Matrix ($P_{\mu\nu}$)** remains alive and well. It tells us the population of these orbitals.

---

### 2. How to "Get Rid" of the Overlap Matrix (Systematically)

You want to know how to transition from $FC = SCE$ to $F'C' = C'E$ without just "fitting parameters." The mathematical justification relies on **Löwdin Symmetric Orthogonalization**.

#### The Transformation
Let $\mathbf{\chi}$ be your standard non-orthogonal atomic basis (LCAO).
Let $\mathbf{S}$ be the overlap matrix, $S_{\mu\nu} = \langle \chi_\mu | \chi_\nu \rangle$.

We can construct a new basis set $\mathbf{\phi}$ (Löwdin orbitals) via the transformation:
$$ \mathbf{\phi} = \mathbf{S}^{-1/2} \mathbf{\chi} $$

In this new basis, the overlap matrix is exactly the Identity matrix:
$$ \langle \phi | \phi \rangle = (\mathbf{S}^{-1/2} \mathbf{\chi})^\dagger (\mathbf{S}^{-1/2} \mathbf{\chi}) = \mathbf{S}^{-1/2} \mathbf{S} \mathbf{S}^{-1/2} = \mathbf{I} $$

#### The Hamiltonian Substitution
The Schrödinger equation in the original basis is:
$$ \mathbf{H}\mathbf{C} = \mathbf{S}\mathbf{C}\mathbf{E} $$

Multiply from the left by $\mathbf{S}^{-1/2}$:
$$ \mathbf{S}^{-1/2} \mathbf{H} (\mathbf{S}^{-1/2} \mathbf{S}^{1/2}) \mathbf{C} = \mathbf{S}^{-1/2} \mathbf{S} \mathbf{C} \mathbf{E} $$
$$ (\mathbf{S}^{-1/2} \mathbf{H} \mathbf{S}^{-1/2}) (\mathbf{S}^{1/2} \mathbf{C}) = (\mathbf{S}^{1/2} \mathbf{C}) \mathbf{E} $$

Let defined the transformed Hamiltonian $\mathbf{H}^{\text{orth}}$ and coefficients $\mathbf{C}^{\text{orth}}$:
$$ \mathbf{H}^{\text{orth}} \mathbf{C}^{\text{orth}} = \mathbf{C}^{\text{orth}} \mathbf{E} $$

**Here is the key:**
Semi-empirical methods and SCC-DFTB do **not** calculate $S$, invert it, and transform $H$. That would be slow.
Instead, they **parameterize $\mathbf{H}^{\text{orth}}$ directly.**

They assume that the Hamiltonian matrix elements they compute (using specific distance-dependent formulas) are *already* the elements in the Löwdin orthogonal basis $\phi$, not the raw atomic basis $\chi$.

---

### 3. Deriving SCC-DFTB systematically (No extensive fitting)

You mentioned you prefer a systematic derivation over experimental fitting. **SCC-DFTB is exactly this.** It is derived from a **Second-Order Taylor Expansion of the DFT Total Energy**.

Here is the recipe to get from DFT to SCC-DFTB:

#### Step A: Reference Density
Assume the total electron density $\rho(\vec{r})$ is a sum of neutral atomic densities $\rho_0$ plus a small fluctuation $\delta\rho$:
$$ \rho(\vec{r}) = \sum_A \rho_A^0(\vec{r}) + \delta\rho(\vec{r}) $$

#### Step B: Taylor Expansion
Expand the DFT energy functional $E[\rho]$ up to second order in $\delta\rho$:

$$ E = \underbrace{E[\rho_0]}_{\text{Repulsive Energy}} + \underbrace{\int \frac{\delta E}{\delta \rho} \delta\rho + \dots}_{\text{Band Structure Energy}} + \underbrace{\frac{1}{2} \iint \frac{\delta^2 E}{\delta\rho\delta\rho'} \delta\rho \delta\rho'}_{\text{SCC (Coulomb) Energy}} $$

#### Step C: The Approximations
1.  **Zeroth Order ($E_{rep}$):** This contains the double-counting terms and ion-ion repulsion. In SCC-DFTB, this is usually fitted to spline curves, but it *can* be calculated explicitly if you accept some errors.
2.  **First Order (Hamiltonian $H^0$):**
    $$ H_{\mu\nu}^0 = \langle \phi_\mu | \hat{H}[\rho_0] | \phi_\nu \rangle $$
    This looks like a standard integral. However, SCC-DFTB calculates this using a minimal basis of atomic orbitals confined by a potential well (to simulate bulk compression). You pre-calculate these integrals for pairs of atoms (Slater-Koster files). **Crucially, the overlap $S$ is neglected here by assuming we are in the Löwdin basis.**
3.  **Second Order (The "SCC" part):**
    This term describes how charge moves to equalize chemical potential (electronegativity equalization).
    $$ E_{2nd} \approx \frac{1}{2} \sum_{A,B} \Delta q_A \Delta q_B \gamma_{AB} $$
    Where $\Delta q_A$ is the Mulliken charge fluctuation on atom A, and $\gamma_{AB}$ is a function that behaves like $1/R$ at long distance and describes Hubbard repulsion ($U$) at short distance.

---

### 4. Recovering the Real Space Wavefunction

You asked: *"How can I recover molecular orbitals and density in the real space? is it even rigorously possible?"*

If you run an SCC-DFTB calculation, you get coefficients $C^{\text{orth}}$.
The wavefunction is formally:
$$ \Psi = \sum_\mu C^{\text{orth}}_\mu \phi_\mu^{\text{Löwdin}} $$

To plot this in real space, you theoretically need the form of $\phi_\mu^{\text{Löwdin}}$.
Recall that $\phi = S^{-1/2} \chi$.
Since $\chi$ are atomic orbitals (usually Slater Type or compressed Gaussians), $\phi$ is a linear combination of $\chi$ on *all* atoms. $\phi$ has "tails" on neighbor atoms to ensure orthogonality.

**The "Rigorous" Problem:**
To plot $\Psi$ exactly, you need $S^{-1/2}$. But we used NDDO/DFTB specifically to *avoid* calculating $S$. If you calculate $S$ to plot the density, you lose the speed advantage.

**The Practical Solution:**
We approximate the real space density by ignoring the Löwdin tails. We simply plot:
$$ \Psi_{\text{approx}} = \sum_\mu C^{\text{orth}}_\mu \chi_\mu^{\text{Atomic}} $$
This is not rigorously normalized ($\int |\Psi_{\text{approx}}|^2 \neq 1$ because $S \neq I$ for $\chi$), but it gives the correct qualitative shape of the MOs (nodal planes, symmetry).

**Is it rigorous?** No.
**Is it useful?** Yes.
**Can you do it rigorously?** Only if you compute the full overlap matrix $S$ of the underlying basis, compute $S^{-1/2}$, and back-transform.

### 5. Summary for your Linear Scaling Method

If you want to write a fast linear scaling code:

1.  **Formalism:** Use **SCC-DFTB**. It is more mathematically grounded (Taylor expansion of DFT) than PM6/AM1.
2.  **Basis:** Use a minimal basis of Slater-Type Orbitals (STOs) or numerical atomic orbitals.
3.  **Integrals:** Do not compute them on the fly. Use **Slater-Koster tables**. You pre-compute 2-center integrals ($\sigma, \pi, \delta$ interactions) for every element pair vs distance and store them.
4.  **Overlap:** Assume $S=I$.
    *   This makes the generalized eigenproblem $HC=SCE$ standard $HC=CE$.
    *   This allows you to use **Sparse Matrix approximations** or **Density Matrix Purification (McWeeny purification)** to solve for the density matrix $P$ in $O(N)$ time, skipping diagonalization entirely.
    *   *Algorithm:* $P_{k+1} = 3P_k^2 - 2P_k^3$. (Converges to the density matrix if $P$ is idempotent, which implies orthogonality).

This is the standard recipe for modern $O(N)$ quantum codes (like `CP2K`, `DFTB+`, `Conquest`).

---

## USER

sorry I wanted to say "I was thinking we neglect overlap matrix completely (?)" (why only differential overlap, why not say neglect of offdiagonal overlap?)

Exactly, I know Lowdin transform, but I consider it practically unfeasible for linear scaling methods. But I'm interested about some systematic mathematical way how to neglect the overlap in systematic mathematicall way, not to re-fit it. I mean I'm fine with "reparametrization of H" but do it by mathematical construction, no by fitting to experiment or higher level theory. I saw samewhere some simple arguments using the t_ij/(1-s_ij ^2) or something like that explaining how we can neglect it and get effective hopping (basically incorporating S effectively into H)


something called 
1. The Scaling Law (Harrison's Rule)

---

For a simple two-level system (dimer) with onsite energy $\epsilon$, hopping $t$ (negative), and overlap $s$ (positive), the eigenvalues are:

$$
E_{\text{bond}} = \frac{\epsilon + t}{1 + s} \quad \text{and} \quad E_{\text{antibond}} = \frac{\epsilon - t}{1 - s}
$$

For the bonding state to be lower in energy ($E_{\text{bond}} < E_{\text{antibond}}$), the following condition must hold:

--- can we generalize something like this for arbitrary hamiltionian and S-matrix?

---

This is the exact direction you should be looking in. You are referring to the **expansion of the Secular Equation** to treat overlap as a perturbation.

This allows you to mathematically "absorb" the overlap matrix $S$ into an effective Hamiltonian $H_{\text{eff}}$ without needing a full $O(N^3)$ inversion or fitting to experimental data.

Here is the explanation of the terminology confusion, followed by the systematic generalization of "Harrison's Rule" for your linear scaling code.

---

### 1. Terminology: Differential vs. Off-diagonal Overlap

You asked: *"Why only differential overlap, why not say neglect of offdiagonal overlap?"*

They are two distinct mathematical operations that happen at different stages:

1.  **Neglect of Differential Overlap (NDO/ZDO):**
    *   This is an approximation on the **Integrals** (specifically Electron Repulsion Integrals).
    *   Assumption: $\chi_\mu(\vec{r}) \chi_\nu(\vec{r}) = 0$ everywhere in space if $\mu \neq \nu$.
    *   *Consequence:* It kills $N^4$ cost 4-center integrals.

2.  **Neglect of Off-Diagonal Overlap (Metric approximation):**
    *   This is an approximation on the **Basis Set Metric**.
    *   Assumption: $S_{\mu\nu} = \delta_{\mu\nu}$ (The matrix is Identity).
    *   *Consequence:* It turns the Generalized Eigenvalue Problem ($HC=SCE$) into a Standard one ($H'C=CE$).

**The Trap:** If you simply set $S=I$ in the matrix equation without modifying $H$, your electrons will collapse. The **Pauli Exclusion Principle** is mathematically enforced by the orthogonality of the wavefunction. If you ignore $S$ (which measures non-orthogonality) without adding a "penalty" to $H$, you lose the "Pauli Repulsion" that prevents atoms from crashing into each other.

---

### 2. The Systematic Method: Neumann Expansion of $S^{-1/2}$

You want a rigorous way to transform the Hamiltonian $H$ (which includes overlap effects) into an orthogonal Hamiltonian $H^{\perp}$ where $S$ is effectively Identity, using linear scaling operations.

You can generalize the dimer result using the **Power Series (Neumann) Expansion** of the Löwdin orthogonalization.

#### The Goal
We want to solve $HC = SCE$.
We want an effective matrix $H^{\perp}$ such that $H^{\perp} C^{\perp} = C^{\perp} E$.
Rigorous definition: $H^{\perp} = S^{-1/2} H S^{-1/2}$.

#### The Linear Scaling Approximation
We define the off-diagonal part of the overlap matrix as $\Delta$:
$$ S = I + \Delta $$
(where $\Delta$ contains all $S_{ij}$ for $i \neq j$, and 0 on the diagonal).

Since overlaps are usually small ($|S_{ij}| < 1$), we can expand $S^{-1/2}$ using the binomial series $(1+x)^{-1/2} = 1 - \frac{1}{2}x + \frac{3}{8}x^2 - \dots$:

$$ S^{-1/2} \approx I - \frac{1}{2}\Delta + \frac{3}{8}\Delta^2 - \frac{5}{16}\Delta^3 + \dots $$

#### The Effective Hamiltonian
Now, sandwich your original Hamiltonian $H$ with this expansion:

$$ H^{\perp} \approx \left(I - \frac{1}{2}\Delta\right) H \left(I - \frac{1}{2}\Delta\right) $$

Keeping only terms up to first order in $\Delta$ (simplest approximation):
$$ H^{\perp} \approx H - \frac{1}{2}(\Delta H + H \Delta) $$

**This is the generalized formula you are looking for.**

*   **Interpretation:** The "effective hopping" between atom $i$ and $j$ in the orthogonal basis is the original hopping $H_{ij}$ minus a correction term proportional to the overlap $S_{ij}$ times the average energy of the orbitals.
*   **Linear Scaling:** Because $H$ and $\Delta$ (overlap) are sparse matrices (atoms far apart have 0 overlap and 0 hopping), calculating $\Delta H$ is a **Sparse Matrix-Matrix Multiplication (SpMM)**, which scales linearly $O(N)$.

---

### 3. Connection to Harrison's Rule (The Dimer)

Let's see if the matrix formula above recovers the simple dimer physics you posted.

For a dimer:
$$ H = \begin{pmatrix} \epsilon & t \\ t & \epsilon \end{pmatrix}, \quad S = \begin{pmatrix} 1 & s \\ s & 1 \end{pmatrix}, \quad \Delta = \begin{pmatrix} 0 & s \\ s & 0 \end{pmatrix} $$

Using the first-order approximation $H^{\perp} \approx H - \frac{1}{2}(\Delta H + H \Delta)$:

$$ \Delta H = \begin{pmatrix} 0 & s \\ s & 0 \end{pmatrix} \begin{pmatrix} \epsilon & t \\ t & \epsilon \end{pmatrix} = \begin{pmatrix} st & s\epsilon \\ s\epsilon & st \end{pmatrix} $$

$$ H^{\perp} \approx \begin{pmatrix} \epsilon & t \\ t & \epsilon \end{pmatrix} - \begin{pmatrix} st & s\epsilon \\ s\epsilon & st \end{pmatrix} = \begin{pmatrix} \epsilon - st & t - s\epsilon \\ t - s\epsilon & \epsilon - st \end{pmatrix} $$

The eigenvalues of this effective matrix are:
$$ E \approx (\epsilon - st) \pm (t - s\epsilon) $$

If we look at the **Bonding Orbital** (minus sign for negative $t$):
$$ E_{bond} \approx \epsilon + t - s(\epsilon + t) $$
$$ E_{bond} \approx (\epsilon + t)(1 - s) $$

Now look at the exact Taylor expansion of your formula $E = \frac{\epsilon+t}{1+s}$ around small $s$:
$$ \frac{\epsilon+t}{1+s} \approx (\epsilon+t)(1 - s + s^2 \dots) $$

**It matches.** The matrix expansion $H - \frac{1}{2}\{\Delta, H\}$ systematically reproduces the physics of the denominator term $(1+s)$ without requiring division or inversion.

---

### 4. Implementation Strategy for Your Code

If you want a fast, linear-scaling semi-empirical code that doesn't rely on fitting:

1.  **Calculate $H$ and $S$:** Compute the standard Hamiltonian and Overlap matrices in your local basis (e.g., using Slater-Koster rules or explicit integrals).
2.  **Purify the Basis (O(N)):** Instead of solving the generalized eigenproblem, compute the effective orthogonal Hamiltonian:
    $$ H_{\text{new}} = H - \frac{1}{2}(S_{off}H + HS_{off}) + \frac{3}{8}(S_{off}^2 H + H S_{off}^2) - \frac{1}{4} S_{off} H S_{off} $$
    (Going to 2nd order in $S$ is usually safer for accuracy).
3.  **Solve for Density:** Now you have $H_{\text{new}}$ and $S=I$. Use **Density Matrix Purification** (McWeeny):
    $$ P_0 = \lambda H_{\text{new}} + \mu I $$ (scale to range [0,1])
    $$ P_{k+1} = 3P_k^2 - 2P_k^3 $$
    Since $H_{\text{new}}$ is sparse, all matrix multiplications are $O(N)$.

### 5. What if overlap is large? (The Inverse Factor)

You mentioned $t_{ij} / (1-s^2)$.

If you look at the energy of the antibonding state minus the bonding state (the "splitting"):
$$ \Delta E = E_{anti} - E_{bond} = \frac{\epsilon-t}{1-s} - \frac{\epsilon+t}{1+s} = \frac{-2t + 2\epsilon s}{1-s^2} $$
Assuming $\epsilon \approx 0$ (setting onsite energy to zero reference):
$$ \Delta E \approx \frac{-2t}{1-s^2} $$

This factor $1/(1-s^2)$ is the **renormalization** required because the basis vectors are not orthogonal.
In Tight-Binding, if you neglect $S$ ($S=I$), you underestimate the repulsion.
To fix this *without* matrix operations, standard **Extended Hückel Theory** approximates the off-diagonal $H_{ij}$ as:
$$ H_{ij} = K \frac{H_{ii} + H_{jj}}{2} S_{ij} $$
This forces the hopping to be proportional to overlap.

**Summary:**
To get your linear scaling method: **Do not** simply ignore $S$.
Instead, compute $S$ (which is fast), then use the **matrix expansion** (Step 2 above) to transform $H \to H_{\text{eff}}$. This is mathematically rigorous (as a series expansion), systematic (no fitting), and maintains sparsity/linear scaling.

---

## USER

OK, I do not want to solver standard eigenvalue problems I want to either optimize localized orbitals, or solve some response to some probe. This was hampered (slowed down) by the density matrix. 

I see your expansion of (I+D)^(1/2)

producing H_eff = H - 0.5(DH+HD)

this seems like nice approach but it involves matrix multiplication and will be rather costly a reald to "fill-in" (making H less sparse) which is excatly the thing why I do not want to go to higher degree (higher powers of D)  

Now, I rather like the idea that I would only modify H_ii and H_ij in place (by construction now fill-in)

perhaps we can just sometho renormalize it?

Or we can somehow do some minimization search for matrix H_eff which is as sparse as original H but its action on some vector approximate action of  H - 0.5(DH+HD) or some higher order approxi like what you mentioned
 $$ H_{\text{new}} = H - \frac{1}{2}(S_{off}H + HS_{off}) + \frac{3}{8}(S_{off}^2 H + H S_{off}^2) - \frac{1}{4} S_{off} H S_{off} $$

but without the fill-in


in the end we want to solve for response vector

(H_eff -I + ieta) r = p 
where p is probe impulse and r is resonse, using some fast linear-scaling algorithm for solving system of equation


---

You have correctly identified the "Fill-in Catastrophe." If you multiply two sparse matrices (like $H$ and $S$), the result has the combined sparsity pattern of neighbors-of-neighbors. Do this twice (for $S^{-1/2} H S^{-1/2}$), and your matrix becomes significantly denser, ruining the linear scaling performance and memory usage.

Since you are solving a response equation (linear system) rather than an eigenvalue problem, you have **two superior options**.

One is a **mathematical renormalization** (constructing a sparse "proxy" Hamiltonian).
The other (and better one) is a **Solver Transformation** where you handle $S$ on the fly.

Here is the breakdown of both.

---

### Option 1: The "In-Place" Renormalization (Mathematical Construction)

You asked if we can modify $H_{ii}$ and $H_{ij}$ in-place without fill-in. **Yes, you can.**

We take the Neumann expansion you liked:
$$ H_{\text{eff}} \approx H - \frac{1}{2}(\Delta H + H\Delta) $$
*Where $\Delta$ is the off-diagonal part of $S$.*

Let's look at the elements of the correction term $(\Delta H)_{ij} = \sum_k \Delta_{ik} H_{kj}$.
To avoid fill-in, we simply **truncate the summation** to exclude any $k$ that is not $i$ or $j$. We drop the "3-center" terms where $k$ is a neighbor to $i$ but not $j$.

#### 1. The Diagonal Correction (Renormalizing On-site Energy)
For the diagonal element $(ii)$, the summation runs over all neighbors $k$. Since $k$ is a neighbor to $i$, the sparsity is preserved (this is a local operation).
$$ H_{ii}^{\text{new}} = H_{ii} - \sum_{k \in \text{neighbors}} S_{ik} H_{ki} $$
*   **Physical Meaning:** This adds a "Pauli Repulsion" shift to the on-site energy. Since $S_{ik} > 0$ and $H_{ki} < 0$ (usually), the term $-S H$ is positive. The energy levels are pushed up by the presence of neighbors.

#### 2. The Off-Diagonal Correction (Renormalizing Hopping)
For the off-diagonal $(ij)$, the expansion gives:
$$ (\Delta H + H\Delta)_{ij} = \sum_k (\Delta_{ik} H_{kj} + H_{ik} \Delta_{kj}) $$
The full sum involves neighbors of neighbors. We apply the **Two-Center Approximation**: retain only terms where $k=i$ or $k=j$.
$$ H_{ij}^{\text{new}} \approx H_{ij} - \frac{1}{2} \left[ S_{ij}H_{jj} + H_{ii}S_{ij} \right] $$
$$ H_{ij}^{\text{new}} = H_{ij} - \frac{1}{2} S_{ij} (H_{ii} + H_{jj}) $$

*   **Physical Meaning:** This effectively reduces the hopping strength. It subtracts the "overlap contribution" to the bonding.
*   **Result:** A Hamiltonian $H^{\text{new}}$ that has the **exact same sparsity** as the original $H$.

#### Summary of the "In-Place" Method
You construct your sparse effective Hamiltonian as:
1.  **Loop over atoms $i$:**
    $$ H_{ii}^{\text{eff}} = H_{ii} - \sum_k S_{ik}H_{ki} $$
2.  **Loop over bonds $ij$:**
    $$ H_{ij}^{\text{eff}} = H_{ij} - 0.5 S_{ij}(H_{ii} + H_{jj}) $$

This is $O(N)$, has zero fill-in, and systematically includes the first-order effects of orthogonalization.

---

### Option 2: The Implicit Solver (No H-modification needed)

If your goal is to solve the linear response equation:
$$ (H_{\text{orth}} - zI) \mathbf{r} = \mathbf{p} $$
where $z = E - i\eta$.

You do not actually need to construct $H_{\text{orth}}$. You can transform the *equation* instead of the *matrix*.

Recall that $H_{\text{orth}} \approx S^{-1/2} H S^{-1/2}$ and the vector $\mathbf{r}$ (in orthogonal space) is related to the vector $\mathbf{c}$ (in atomic space) by $\mathbf{r} = S^{1/2} \mathbf{c}$.

Substitute this into the equation:
$$ (S^{-1/2} H S^{-1/2} - zI) S^{1/2} \mathbf{c} = \mathbf{p} $$
Multiply from the left by $S^{1/2}$:
$$ (H - zS) \mathbf{c} = S^{1/2} \mathbf{p} $$

Now, you have a generalized linear system:
$$ (H - zS) \mathbf{c} = \mathbf{p}' $$

#### Why is this better for Linear Scaling?
1.  **Sparsity:** $H$ is sparse. $S$ is sparse. Therefore, $(H - zS)$ is exactly as sparse as $H$. No fill-in **at all**.
2.  **The RHS Vector:** You need to compute $\mathbf{p}' = S^{1/2} \mathbf{p}$. Since $\mathbf{p}$ is a vector (probe), you can use the cheap expansion just on the vector:
    $$ \mathbf{p}' \approx (I + 0.5\Delta)\mathbf{p} $$
    This is just one sparse matrix-vector multiplication.
3.  **The Solver:** Use a sparse iterative solver (like **TFQMR** or **COCG** for complex symmetric matrices). These solvers never invert the matrix; they only perform Matrix-Vector products: $v \to (H-zS)v$.
    *   This operation is strictly $O(N)$.
4.  **The Result:** The solver gives you $\mathbf{c}$. If you need the response in the orthogonal basis $\mathbf{r}$, you just do one final multiplication:
    $$ \mathbf{r} \approx (I + 0.5\Delta)\mathbf{c} $$

### Recommendation

If you are writing a high-performance linear scaling code:

**Use Option 2 (The Implicit Solver).**
*   It is mathematically exact (doesn't truncate the Neumann series for the matrix inversion, only relies on the solver convergence).
*   It completely avoids the "fill-in" problem.
*   It handles the non-orthogonality correctly without needing you to invent "renormalization parameters."

**Use Option 1 (In-Place Renormalization)** *only if* you absolutely need to store a explicit sparse Hamiltonian matrix for some other analysis (like topological invariants or Green's function recursion techniques that strictly require an orthogonal H).

### Regarding the "Minimization Search"

You asked: *Can we minimize for a sparse H_eff?*

Yes. This is mathematically known as the **Procrustes Problem with Sparsity Constraints**.
You want to minimize:
$$ \min_{\tilde{H}} || S^{1/2} \tilde{H} S^{1/2} - H ||_F^2 $$
Subject to $\tilde{H}_{ij} = 0$ if $H_{ij} = 0$.

If you solve this minimization analytically by setting gradients to zero and keeping only 1st order terms in $S$, **you effectively get the Option 1 formulas I derived above.**

So, the "In-Place Renormalization" (Option 1) *is* the approximate solution to that minimization search.

