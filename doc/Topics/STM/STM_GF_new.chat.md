https://aistudio.google.com/app/prompts?state=%7B%22ids%22:%5B%221u6ceIskV52PNfeEJXMqyAWBGy-xFgKTW%22%5D,%22action%22:%22open%22,%22userId%22:%22100958146796876347936%22,%22resourceKeys%22:%7B%7D%7D&usp=sharing



---

# USER

OK, this is rather good, now we should try the oteher method based on green-functions.

lets not write any code, rather just discuss. 

lets have some Gamma which is applied selectively to atoms of the tip as a source of electron, and other gamma applied selectively to atoms of sample which is drain of electrons (or vice versa), for the sample molecule and tip molecule we take the hamiltionain in the dense form Hk from fireball. And fill the tip-sample hopping by those exponential overlaps, something similar we already did in @test_response_function.py 
@test_response_function_rotated.py 

so the path of electron is like

1) it jump from out drain to state of the tip
2) it propagate thrgou the tip molecule
3) then it hop by the exponential hopping to the sample
4) then it propagare through the sample 
5) then it hopp to the drain gamma on sample

this progrss should be possible to describe by greensfunction formalism efficiently right

Consider we do everuthing in the basisset pf 4 basis per atom s,px,py,pz with Fireball basiset

Perhaps we aready have something like that noted in @STM_GPU_QMMM.md right?

do not write any code, just equations, physics/math notes, analyse problem do not write any code

https://windsurf.com/codemaps/9fa40c64-e78c-42f2-9573-574936c8040d-fe86ab10a43f3d18

---

# Gemini-3-flash

This Green's function (GF) approach is a powerful alternative to the response-function/perturbation theory because it naturally handles the "all-orders" propagation and the effects of the leads (source/drain) through self-energies.

Here is the mathematical and physical analysis of the transport problem as you described it.

### 1. System Partitioning
We divide the Hilbert space into two main blocks: the **Tip (T)** and the **Sample (S)**.
The total Hamiltonian $H$ is:
$$H = \begin{pmatrix} H_{TT} & H_{TS} \\ H_{ST} & H_{SS} \end{pmatrix}$$

*   $H_{TT}$: Dense Hamiltonian of the tip molecule (Fireball).
*   $H_{SS}$: Dense Hamiltonian of the sample molecule (Fireball).
*   $H_{TS}(R)$: The coupling matrix, which depends on the relative tip-sample position $R$. This is populated by the exponential overlaps of the basis functions (s, px, py, pz).

### 2. Modeling Source and Drain (Self-Energies)
Instead of including explicit leads, we use the **Wide Band Limit (WBL)** or specific Broadening matrices $\Gamma$ applied to "selective atoms."

Let $\Gamma_L$ be the source (Lead on Tip) and $\Gamma_R$ be the drain (Lead on Sample).
We define the **Self-Energies** $\Sigma$ as:
$$\Sigma_L = -i \frac{\Gamma_L}{2}, \quad \Sigma_R = -i \frac{\Gamma_R}{2}$$

$\Gamma_L$ is a matrix that is non-zero only for the orbitals of the "top" atoms of the tip. Similarly, $\Gamma_R$ is non-zero only for the "bottom" atoms of the sample (or whichever atoms are connected to the substrate).

### 3. The Retarded Green's Function
The retarded Green's function for the coupled system at energy $E$ (usually Fermi level $E_F$) is:
$$G^R(E) = [ (E + i\eta)I - H - \Sigma_L - \Sigma_R ]^{-1}$$

To make this efficient, we define the "isolated" Green's functions (including their respective leads but not each other):
$$g_{T}(E) = [ (E + i\eta)I - H_{TT} - \Sigma_L ]^{-1}$$
$$g_{S}(E) = [ (E + i\eta)I - H_{SS} - \Sigma_R ]^{-1}$$

Since $H_{TT}$ and $H_{SS}$ are constant during a scan, **$g_T$ and $g_S$ can be pre-calculated and stored.**

### 4. The Path Analysis (Dyson Equation)
The electron path you described (Source $\to$ Tip $\to$ Tunneling $\to$ Sample $\to$ Drain) is captured by the off-diagonal block of the full Green's function $G_{TS}$.

Using the Dyson equation, the full Green's function of the combined system can be expressed in terms of the isolated $g_T, g_S$ and the coupling $H_{TS}$:
$$G_{TS} = g_T H_{TS} G_{SS}$$
$$G_{SS} = [ I - g_S H_{ST} g_T H_{TS} ]^{-1} g_S$$

In the **Weak Coupling Regime** (typical for STM tunneling), we can use the first-order approximation:
$$G_{TS} \approx g_T H_{TS} g_S$$

### 5. Transmission and Current
The Transmission function $T(E)$, which is proportional to the conductance, is given by the Fisher-Lee formula:
$$T(E) = \text{Tr} [ \Gamma_L G_{TS} \Gamma_R G_{ST} ]$$

Substituting the weak-coupling approximation:
$$T(E) \approx \text{Tr} [ \Gamma_L (g_T H_{TS} g_S) \Gamma_R (g_S^\dagger H_{ST} g_T^\dagger) ]$$

### 6. Computational Efficiency for Scanning
This is where the GF method becomes superior for the GPU implementation mentioned in `@STM_GPU_QMMM.md`:

1.  **Pre-computation:** Invert $H_{TT}$ and $H_{SS}$ once to get $g_T$ and $g_S$.
2.  **Projection:** $\Gamma_L$ and $\Gamma_R$ are very sparse (only a few atoms). We only need the sub-blocks of $g_T$ and $g_S$ that correspond to:
    *   The atoms connected to the leads.
    *   The atoms involved in the tunneling (the "interface" atoms).
3.  **The Scan Loop:** As the tip moves, only $H_{TS}$ changes. $H_{TS}$ is also sparse because it only involves atoms within a certain cutoff distance.
4.  **Matrix Contraction:** The calculation of $T(E)$ becomes a sequence of small matrix multiplications:
    $$T(E) = \text{Tr} [ \tilde{g}_T H_{TS} \tilde{g}_S H_{ST} ]$$
    where $\tilde{g}$ are the "effective" propagators filtered by the $\Gamma$ matrices.

### 7. Physical Interpretation of the Path
*   **$g_T \Gamma_L g_T^\dagger$:** Represents the Spectral Function (density of states) of the tip as "seen" from the tunneling interface, including the broadening from the source wire.
*   **$H_{TS}$:** The "jump" across the vacuum gap.
*   **$g_S \Gamma_R g_S^\dagger$:** Represents the Spectral Function of the sample as "seen" from the tunneling interface, including the broadening from the substrate.

### Summary of the "Equation Flow"
To calculate the current at a point $(x, y, z)$:
1.  Construct $H_{TS}$ based on atomic distances and Fireball orbital overlaps.
2.  Compute $M = g_T H_{TS} g_S$.
3.  Compute $T = \text{Tr} [ \Gamma_L M \Gamma_R M^\dagger ]$.
4.  Current $I \propto \int T(E) [f_T(E) - f_S(E)] dE$. (Usually simplified to $T(E_F) \cdot V_{bias}$).

This avoids the $O(N^3)$ cost of a full diagonalization at every pixel, reducing it to a sparse matrix multiplication $O(N_{active}^2)$, which is perfect for GPU acceleration.

---

# Kimi 2.6

I'll try to access the codemap you shared to get the full context, then we'll dive into the physics.
The codemap didn’t fully load for me, but from your description this is the classic **partitioned non-equilibrium Green’s function (NEGF) problem for STM/molecular junctions**, and yes — it maps exactly to the 5-step process you described. The key is that the tip–sample coupling is localized and (usually) weak, so you never need to invert the full $(N_{\rm tip}+N_{\rm sample})$ Hamiltonian. Below is the full physical/mathematical breakdown.

### 1. Hamiltonian and Basis Structure

You have a Fireball basis with **4 functions per atom**: $s, p_x, p_y, p_z$.  
Label tip basis states by Greek indices $\mu,\nu \in {\cal T}$ and sample states by $\alpha,\beta \in {\cal S}$.

Because Fireball uses numerical atomic orbitals, the full problem is **non-orthogonal**:

$$
S = \begin{pmatrix} S_{tt} & S_{ts} \\ S_{st} & S_{ss} \end{pmatrix}, \qquad
H = \begin{pmatrix} H_{tt} & H_{ts} \\ H_{st} & H_{ss} \end{pmatrix}
$$

* $H_{tt}, S_{tt}$ — dense Fireball blocks for the isolated tip molecule.
* $H_{ss}, S_{ss}$ — dense Fireball blocks for the isolated sample molecule.
* $H_{ts}, S_{ts}$ — the **exponential overlap/hopping** you already compute in `test_response_function.py` (Slater–Koster direction cosines $\times$ exponential radial decay). For STM distances this is small.

The source and drain are not macroscopic leads here; they are **local broadenings** (self-energies) attached to selected atoms:

$$
\Sigma_L(E) \;\; \text{non-zero only on tip atoms (source)}
$$
$$
\Sigma_R(E) \;\; \text{non-zero only on sample atoms (drain)}
$$

In the simplest wide-band limit you just add a constant imaginary part to those atomic diagonal blocks:
$\Sigma_{L,R} \approx \mp i\Gamma_{L,R}/2$.

### 2. Uncoupled Propagators (Steps 2 & 4)

Before coupling tip and sample, define the **retarded Green’s functions** of the two isolated, damped subsystems:

$$
g_t(E) = \big[(E+i0^+)S_{tt} - H_{tt} - \Sigma_L(E)\big]^{-1}
$$

$$
g_s(E) = \big[(E+i0^+)S_{ss} - H_{ss} - \Sigma_R(E)\big]^{-1}
$$

These describe steps 2 and 4 of your list: propagation through the tip and through the sample, **already including** the source/drain coupling at their endpoints.

### 3. The Coupling Matrix (Step 3)

In the non-orthogonal basis the energy-dependent coupling block is:

$$
M_{ts}(E) \equiv (E+i0^+)S_{ts} - H_{ts}
$$

*If* the tip–sample distance is large enough that $S_{ts} \approx 0$, this reduces to $M_{ts} \approx -H_{ts} \equiv V$, the exponential hopping matrix.  
If you are in the **contact regime** ($< 3$ Å), keep the full $M_{ts}(E)$ because the $ES_{ts}$ term matters.

### 4. Dyson Partitioning — The Efficient Resummation

You do **not** invert the full matrix. Instead, use the block-Dyson identity. The exact off-diagonal Green’s function (the amplitude to go from sample to tip) is:

$$
G_{ts}(E) = \big[I - g_t(E)\,M_{ts}(E)\,g_s(E)\,M_{st}(E)\big]^{-1}\, g_t(E)\,M_{ts}(E)\,g_s(E)
$$

Equivalently, one often writes the “contact” self-energy that the sample feels from the tip:

$$
\Sigma_{\rm contact}(E) = M_{st}(E)\,g_t(E)\,M_{ts}(E)
$$

then

$$
G_{ss}(E) = \big[g_s^{-1}(E) - \Sigma_{\rm contact}(E)\big]^{-1}
$$

and $G_{ts} = g_t M_{ts} G_{ss}$.

**Why this is efficient:**  
$M_{ts}$ is sparse and short-ranged. Only atoms near the tip apex and sample surface have non-zero entries. If that “contact” region contains $N_c$ basis functions (e.g. 4 bases $\times$ 10 atoms = 40), the matrix to invert is $N_c \times N_c$ — **tiny** compared to inverting the full $H$.

### 5. Transmission — Your 5 Steps in One Formula (Caroli Formula)

The transmission probability for an electron at energy $E$ is:

$$
T(E) = {\rm Tr}\Big[\,\Gamma_L(E)\,G^R_{ts}(E)\,\Gamma_R(E)\,G^A_{st}(E)\,\Big]
$$

with the broadening matrices:

$$
\Gamma_{L,R}(E) = i\big[\Sigma^R_{L,R}(E) - (\Sigma^R_{L,R})^\dagger(E)\big]
$$

Because $\Gamma_L$ lives only on the tip and $\Gamma_R$ only on the sample, the trace picks out exactly the $G_{ts}$ and $G_{st}$ blocks — the tunneling amplitude.

### 6. Weak-Coupling (Bardeen / Transfer-Hamiltonian) Limit

For STM, $M_{ts}$ is exponentially small. Expanding $G_{ts}$ to first order:

$$
G_{ts}(E) \;\approx\; g_t(E)\,M_{ts}(E)\,g_s(E)
$$

The transmission becomes:

$$
T(E) \approx {\rm Tr}\Big[\, \underbrace{(g_t \Gamma_L g_t^\dagger)}_{A_t(E)}\, M_{ts}\, \underbrace{(g_s \Gamma_R g_s^\dagger)}_{A_s(E)}\, M_{st} \,\Big]
$$

Now your 5 steps are manifest:

| Step | Physical object | Mathematical expression |
|------|----------------|------------------------|
| 1. Injection into tip | Source broadening | $\Gamma_L$ |
| 2. Propagate through tip | Tip spectral function | $A_t = g_t \Gamma_L g_t^\dagger$ |
| 3. Hop to sample | Exponential coupling | $M_{ts}(E)$ |
| 4. Propagate through sample | Sample spectral function | $A_s = g_s \Gamma_R g_s^\dagger$ |
| 5. Extract at drain | Drain broadening | $\Gamma_R$ |

So $T(E)$ is the spectral overlap of tip and sample, weighted by the exponential coupling:

$$
T(E) \approx {\rm Tr}\big[\,A_t(E)\,M_{ts}(E)\,A_s(E)\,M_{st}(E)\,\big]
$$

This is the **Bardeen–Todorov–Pendry** STM formula, just written in your 4-basis Fireball notation.

### 7. What About Non-Orthogonality?

Fireball’s $S \neq I$ creates two subtleties:

1. **Within each block** ($S_{tt}, S_{ss}$): You can Löwdin-orthogonalize each molecule independently before coupling, turning $g_t, g_s$ into standard orthogonal propagators. This is exact for the isolated blocks.

2. **Cross-overlap $S_{ts}$**: In the tunneling limit you can set $S_{ts}=0$ and $M_{ts}=-H_{ts}$. If you are in contact, keep $S_{ts}$ inside $M_{ts}(E)$ as shown above. The partitioned formulas remain identical; you just use the full $M_{ts}$.

### 8. Current and Bias

If you apply a bias $V_b$ between tip and sample, shift the Hamiltonians:

$$
H_{tt} \to H_{tt} + \frac{eV_b}{2}, \qquad H_{ss} \to H_{ss} - \frac{eV_b}{2}
$$

(or whatever potential division your geometry implies). The current is:

$$
I = \frac{2e}{h}\int dE\, T(E)\,\big[f_L(E) - f_R(E)\big]
$$

where $f_{L,R}$ are Fermi functions of the source and drain reservoirs.  
At finite bias you technically need the **non-equilibrium** density matrix (Keldysh formalism):

$$
G^< = G^R \Sigma^< G^A, \qquad \Sigma^< = \Gamma_L f_L + \Gamma_R f_R
$$

but the current still reduces to the Caroli formula above for non-interacting electrons.

### 9. Relation to `STM_GPU_QMMM.md`

I couldn’t read the full file, but from the title (*STM GPU QMMM*) it likely deals with:

* **QM/MM partitioning**: Tip = QM (Fireball), Sample = QM or MM.
* **GPU density projection**: Computing $\rho({\bf r})$ from the density matrix.

The Green’s function framework above gives you exactly the **non-equilibrium density matrix** you need for that projection:

$$
\rho = \frac{1}{2\pi i}\int dE\, G^<(E)
$$

Once you have $G^R$ and $G^<$ from the partitioned formulas, you can project the electron density onto a real-space grid GPU-side without ever storing the full $(N_t+N_s)^2$ matrix in the coupled basis — you only need the contact-region $G_{ts}$ and the local blocks $G_{tt}, G_{ss}$.

### Summary

Your physical picture is **exactly** the partitioned NEGF prescription:

1. Build $g_t(E)$ and $g_s(E)$ from Fireball $H_{k}$ + local Gammas.
2. Build $M_{ts}(E)$ from exponential Fireball overlaps/hoppings.
3. Compute $G_{ts}(E)$ via the Dyson resummation (small contact-matrix inverse).
4. Obtain $T(E)$ from the Caroli trace.
5. Integrate for current, or use $G^<$ for real-space density.

The exponential coupling makes step 3 cheap, and the partitioning means you never pay the $O((N_t+N_s)^3)$ price — you pay $O(N_t^3) + O(N_s^3)$ once per energy, plus $O(N_c^3)$ for the contact resummation. That is precisely why Green’s functions are the natural language for this problem.


---

# USER

OK, lets now thinkg about numerical efficiency. You say it is fast because we can comput it by sparse matrix multiplication. But it seems to me that G of tip and of sample is not sparse at all, and also when we do matrix multiplication they are densifying. 

I was rather thingking aboutpropagating vectro (representing wavefunction, but not eigenfunction, rather some impulse response) over the system by either just multiplying it by the matrixes, or by solving linear system.

I mean instead of getting G matrix explicitly as G(E) = 1/(EI-H) and then applying to vectro 
response = G(E)*source

we can always just solve the linear equation
 (EI-H) response = source
A=(EI-H)
b=source
x=response

Ax=b

but I'm now thinking how to apply it on this situation - basically how to replace proper Caroli formula with L,R lead broadening (gamma) by solving linear system with modified hamilatioan

The H_tot = H_sample + H_tip + H_hopping (tip-sample) and hermitean conjugate

H_hopping is made form exponential overlaps between tip and sample basisfuctions

I can compute big G_tot(E) = (E I-H_tot) right ?

---

# Gemini-3-flash

You are right to point out that explicitly calculating and storing the Green's function matrix $G = (EI - H)^{-1}$ is inefficient because, while $H$ is sparse (Fireball cutoffs), $G$ is dense.

Using a **Linear System Solver ($Ax = b$)** approach is much more efficient, especially for the "Factorize Once, Solve Many" scenario of an STM scan.

Here is the physical and mathematical analysis of how to map the Caroli/Green's function problem onto a linear system of vector propagations.

### 1. The Effective Non-Hermitian Hamiltonian
To account for the electron source (Tip lead) and drain (Sample lead), we modify the Hamiltonian by adding the self-energies $\Sigma$. 
Let $A$ be the operator we want to invert:
$$A(E) = (E + i\eta)I - H_{tot} - \Sigma_L - \Sigma_R$$
Where:
*   $H_{tot} = H_{T} + H_{S} + H_{TS}(R)$
*   $\Sigma_L = -i\Gamma_L/2$ (Source on Tip)
*   $\Sigma_R = -i\Gamma_R/2$ (Drain on Sample)

Because of $\Sigma$, $A$ is **complex and non-Hermitian**.

### 2. Decomposing the Transmission (The "Source" Vector)
The Caroli formula for transmission is $T = \text{Tr}[\Gamma_L G \Gamma_R G^\dagger]$.
If we assume the source $\Gamma_L$ is localized on a few "top" atoms of the tip, we can decompose $\Gamma_L$ into a sum of outer products (eigen-vectors of the lead coupling):
$$\Gamma_L = \sum_k \gamma_{L,k} |s_k\rangle \langle s_k|$$
Where $|s_k\rangle$ are the orbital basis vectors of the tip-source atoms.

The total transmission then becomes a sum over independent "source" excitations:
$$T = \sum_k \gamma_{L,k} \left( \langle s_k | G^\dagger \Gamma_R G | s_k \rangle \right)$$

### 3. The "Impulse Response" Vector Propagation
For a single source vector $|s_k\rangle$, we define the **response vector** $|\psi_k\rangle$:
$$|\psi_k\rangle = G |s_k\rangle$$
This is exactly the solution to the linear system:
$$A(E) |\psi_k\rangle = |s_k\rangle$$

Once we have the response vector $|\psi_k\rangle$, the transmission for that source is simply the "drainage" on the sample side:
$$T_k = \langle \psi_k | \Gamma_R | \psi_k \rangle$$

### 4. Efficient Scanning (The Dyson/Woodbury Strategy)
In an STM scan, $H_T$, $H_S$, $\Sigma_L$, and $\Sigma_R$ are **constant**. Only $H_{TS}(R)$ changes.
Solving the full system $A |\psi\rangle = |s\rangle$ at every pixel $(x,y,z)$ is expensive if the sample is large. We can use the fact that $H_{TS}$ is a **low-rank coupling** (only the tip apex and the nearby surface atoms are involved).

#### Step A: Pre-calculate Tip Response
Since the source is on the tip, and if we assume the tip-sample coupling is a perturbation (weak coupling):
1.  Solve for the tip-only response: $(E I - H_T - \Sigma_L) |\phi_T\rangle = |s\rangle$.
2.  This $|\phi_T\rangle$ is the "active" wavefunction at the tip apex.

#### Step B: Solve for Sample Response (The Scan Loop)
For each tip position $R$:
1.  Define the "Injection Vector" into the sample: $|b_S(R)\rangle = H_{ST}(R) |\phi_T\rangle$.
    *   This vector is very sparse (only non-zero on sample atoms near the tip).
2.  Solve the linear system for the sample response:
    $$((E + i\eta)I - H_S - \Sigma_R) |\psi_S\rangle = |b_S(R)\rangle$$
3.  Calculate $T = \langle \psi_S | \Gamma_R | \psi_S \rangle$.

### 5. Why this is numerically superior
1.  **Sparsity:** $H_S$ is the Fireball Hamiltonian, which is highly sparse. $A = (EI - H_S - \Sigma_R)$ is a sparse matrix.
2.  **Factorize Once:** We can perform an **LU decomposition** (or LDLT) of the sample matrix $A$ **only once** at the beginning of the scan.
3.  **Back-substitution:** For every pixel $(x,y,z)$, we only perform a "forward/backward substitution" (solving $Ly=b$ and $Ux=y$). This is $O(N^2)$ for dense but for sparse matrices with good ordering (like nested dissection), it is extremely fast—approaching $O(N)$.
4.  **No G-matrix:** We never compute $G$. We only compute the vector $|\psi_S\rangle$.

### 6. Summary of the Path in Vector Form
1.  **Source Injection:** $|s\rangle$ (Source on top of tip).
2.  **Tip Propagation:** $|\phi_T\rangle = A_T^{-1} |s\rangle$ (Pre-computed).
3.  **Vacuum Tunneling:** $|b_S\rangle = H_{ST} |\phi_T\rangle$ (Sparse MVM, computed every pixel).
4.  **Sample Propagation:** $|\psi_S\rangle = A_S^{-1} |b_S\rangle$ (Sparse back-substitution, computed every pixel).
5.  **Drain Collection:** $T = \sum_{j \in \text{drain}} \Gamma_{R,j} |\psi_{S,j}|^2$ (Vector dot product).

This approach handles the full chemistry of the tip and sample, allows for multiple orbitals (s, px, py, pz), and perfectly captures the "path" of the electron while maintaining the computational efficiency needed for 3D scanning. 

Would you like to analyze how the **Self-Energies ($\Gamma$)** should be structured for the Fireball basis specifically (e.g., how to "attach" the lead to the s-orbital vs p-orbitals of a boundary atom)?

---

# Kimi 2.6

You are absolutely right to call that out. The statement “sparse matrix multiplication makes it fast” was sloppy. Here is the honest assessment:

- $g_t(E) = [(E+i\eta)S_{tt}-H_{tt}-\Sigma_L]^{-1}$ is the **inverse of a sparse matrix**, therefore it is **dense**.
- In the partitioned Dyson formula $G_{ts} \approx g_t M_{ts} g_s$, even if the hopping $M_{ts}$ is sparse, multiplying it by dense $g_t$ and $g_s$ immediately densifies everything.
- So the partitioned approach is efficient **only** if the contact region is tiny enough that the dense matrices are small. For a full molecular cluster it is not a sparse algorithm.

Your intuition about **solving linear systems** is the correct and much more powerful route. Below is the full mapping from the Caroli formula to the “impulse-response / wavefunction-solver” picture.

### 1. The Full Sparse System

Build the **total** Hamiltonian and overlap in the full 4-basis-per-atom space:

$$
H_{\rm tot} = \begin{pmatrix} H_{tt} & H_{ts} \\ H_{st} & H_{ss} \end{pmatrix}, \qquad
S = \begin{pmatrix} S_{tt} & S_{ts} \\ S_{st} & S_{ss} \end{pmatrix}
$$

$H_{ts}$ contains your exponential Fireball overlaps; it is sparse and extremely short-ranged (only a handful of tip atoms see a handful of sample atoms). $S_{ts} \approx 0$ in the tunneling regime, but keep it if you are in contact.

Attach the source and drain as **local imaginary potentials** (wide-band limit of the self-energy):

$$
\Sigma_L = -\frac{i}{2}\Gamma_L, \qquad \Sigma_R = -\frac{i}{2}\Gamma_R
$$

where $\Gamma_L$ and $\Gamma_R$ are diagonal matrices (or small dense blocks) living only on the selected tip atoms and sample atoms, respectively.

The retarded Green’s function of the full coupled system is:

$$
G^R(E) = \big[(E+i\eta)S - H_{\rm tot} + \tfrac{i}{2}\Gamma_L + \tfrac{i}{2}\Gamma_R\big]^{-1}
$$

Define the **sparse linear-system matrix**:

$$
A(E) \equiv (E+i\eta)S - H_{\rm tot} + \tfrac{i}{2}\Gamma_L + \tfrac{i}{2}\Gamma_R
$$

$A(E)$ is large but **sparse**. We never compute $A^{-1}$ explicitly.

### 2. From Caroli to Linear Solves

Assume for clarity that $\Gamma_L$ and $\Gamma_R$ are diagonal (one broadening value per basis function on the source/drain atoms). The Caroli transmission is:

$$
T(E) = \operatorname{Tr}\big[\Gamma_L \, G^R(E) \, \Gamma_R \, G^A(E)\big]
$$

Because the Gammas are diagonal, this trace collapses to a double sum over the **source basis** ($a \in L$) and **drain basis** ($b \in R$):

$$
T(E) = \sum_{a \in L}\sum_{b \in R} \Gamma_{L,a}\,\Gamma_{R,b}\,|G^R_{ab}(E)|^2
\tag{1}
$$

Now observe:

$$
G^R_{ab} = \langle a | G^R | b \rangle = \big[G^R \, e_b\big]_a
$$

where $e_b$ is the unit vector localized on drain basis function $b$. Therefore, if we solve the linear system:

$$
A(E)\,\psi_b = e_b \qquad \Longrightarrow \qquad \psi_b = G^R(E)\,e_b
\tag{2}
$$

the $a$-th component of the solution is exactly the Green’s function element we need:

$$
(\psi_b)_a = G^R_{ab}(E)
$$

Insert this into Eq. (1):

$$
\boxed{
T(E) = \sum_{b \in R} \Gamma_{R,b} \left[ \sum_{a \in L} \Gamma_{L,a} \, |(\psi_b)_a|^2 \right]
}
\tag{3}
$$

### 3. Physical Interpretation: Impulse Response

Equation (2) is precisely the **impulse response** you described:

1. **Injection (Step 1):** Place a point source $e_b$ on one basis function $b$ of the sample drain.
2. **Propagation (Steps 2–4):** Solve $A\psi_b = e_b$. The matrix $A$ contains the full Hermitian $H_{\rm tot}$ plus the non-Hermitian “absorption” terms $+i\Gamma/2$ on both source and drain atoms. Because of those imaginary parts, $A$ describes **outgoing boundary conditions**: the electron can propagate through the tip, hop across $H_{ts}$, propagate through the sample, and leak out at the drains. The multiple scattering between tip and sample is resummed to all orders automatically by the linear solve.
3. **Extraction (Step 5):** Look at the amplitude $|(\psi_b)_a|^2$ on each tip basis $a$, weight it by the source strength $\Gamma_{L,a}$, and sum.

Repeat for every basis function $b$ in the drain region, weight by $\Gamma_{R,b}$, and sum. That is the full transmission at energy $E$.

If the drain region contains, say, 10 atoms $\times$ 4 basis = 40 functions, you need **40 sparse linear solves** with the **same matrix** $A(E)$ and 40 different right-hand sides. With a sparse direct solver (MUMPS, PARDISO, SuperLU), you pay for one sparse LU factorization of $A(E)$, then each solve is a cheap back-substitution.

### 4. The Modified Hamiltonian

The matrix $A(E)$ can be viewed as the Schrödinger equation with a **complex, non-Hermitian Hamiltonian**:

$$
\mathcal{H}_{\rm eff}(E) = H_{\rm tot} - \frac{i}{2}\Gamma_L - \frac{i}{2}\Gamma_R
$$

so that:

$$
A(E) = (E+i\eta)S - \mathcal{H}_{\rm eff}
$$

The imaginary terms $-i\Gamma/2$ act as **optical potentials** or **sinks**. They absorb probability current at the source and drain atoms, which is exactly what allows a steady-state current to flow: the electron is injected by the right-hand side source vector, propagates through the molecule, and is absorbed by the $\Gamma$ terms. This is the mathematical implementation of your 5-step physical picture.

### 5. What About “Propagating by Matrix Multiplication”?

You mentioned iterating $\psi \to H\psi$ instead of solving $A\psi = b$. This corresponds to expanding the inverse in a power series:

$$
G^R = \frac{1}{E - \mathcal{H}_{\rm eff}} \approx \frac{1}{E}\sum_{n=0}^{N} \left(\frac{\mathcal{H}_{\rm eff}}{E}\right)^n
$$

This is essentially the **Chebyshev** or **Lanczos** expansion of the Green’s function. It works, but:

- Near a resonance (when $E$ is close to an eigenvalue of $H_{\rm tot}$), the series converges terribly slowly.
- You need $N \sim 1/\eta$ iterations, where $\eta$ is your energy broadening. For sharp molecular levels ($\eta \sim$ meV), this is thousands of iterations.
- Direct sparse solve targets the exact energy $E$ in one shot (up to the linear solver tolerance) and is usually faster for moderate system sizes ($N \lesssim 10^5$ basis functions).

The matrix-multiplication picture becomes advantageous only if:
- The system is so large that sparse direct factorization is impossible, **and**
- You only need a spectral average (then stochastic trace estimation + Chebyshev expansion is the way to go).

For a Fireball STM cluster with a few thousand atoms ($\sim 10^4$ basis functions), the sparse direct solve is the natural choice.

### 6. Non-Orthogonality Note

Because Fireball uses numerical orbitals, do not use $EI - H$. The correct operator is:

$$
A(E) = (E+i\eta)S - H_{\rm tot} + \tfrac{i}{2}\Gamma_L + \tfrac{i}{2}\Gamma_R
$$

If you are in the tunneling regime ($> 3$ Å), $S_{ts} \approx 0$ and $S$ is block-diagonal. You can then Löwdin-orthogonalize each block independently, which turns $A$ into a standard $EI - H_{\rm eff}$ form in the orthogonalized basis. If you are in contact, keep the full $S$; modern sparse direct solvers handle the generalized form $(E S - H)\psi = b$ without trouble.

### Summary

| Quantity | Explicit (bad) | Linear Solve (good) |
|----------|---------------|---------------------|
| Green’s function | $G^R = A^{-1}$ (dense, $N^3$) | Never form $A^{-1}$ |
| Transmission | Tr$[\Gamma_L G^R \Gamma_R G^A]$ with dense matrices | Eq. (3): weighted sum of $|\psi_b|^2$ |
| Propagation | Dense matrix products $g_t M_{ts} g_s$ | Sparse solve $A\psi_b = e_b$ |
| Cost | $O(N^3)$ | One sparse factorization + $N_R$ back-substitutions |

Your instinct is correct: **do not compute Green’s function matrices explicitly**. Build the sparse $A(E)$, factor it, and shoot impulse vectors $e_b$ from the drain. The responses $\psi_b$ contain everything needed for the Caroli transmission.


---

# USER

No the task is to solve all operations on GPU in parallel, each workgroup should be solving one pixel, workgroup can have like 16-32 threads and it should use shared memory (prelod things in shared memory for efficiecy), we want to solve all pixels (thousands pixel of the image) in parallel each by one workgroup. All performance intensive operations should be done on GPU. 
Te linears solver we should therefoer implement direcly in OpenCL kernel (perhaps as separate function) operating with this shared memory. Not sure what is good compromise between simplicity and efficiecy, Not sure what type or matrixes we need to solve there (are they Hermitian, symmetric?.... can we use CG, or Jacobi, iteration, Gauss elimination, GMRESS), plese focus on this aspect.

We expect the size of system to be molecule composed of ~50-100 atom like PTCDA so with 4 atoms on each this is like 200-400 dimension of the matrix. We want to solve this matrix using wogroup of reasonable size (max 128 threads - in a way to be efficient). We wat to paralelize this over pixel to explot paralleism with 10k threads on GPU, why we need to do this? We want to solve each pixel in miliseconds or mircroscetods, generating hundesrtsns or  thoushads images per second during moleculer dynamics of the molecule (molecular dynamics we do by forcefield by other module, this is not the task now, focus on how to quicky evaluate the electron trnsmission resp STM current, whewre the linear solver is probably performance bottleneck)

plese explicitly metion all relevant files and modules for implementing this like what is in

there are some detailes about the code
https://windsurf.com/codemaps/f398c2cf-5ff8-4d75-a398-c83e788e27b4-fe86ab10a43f3d18


I can also probide you with the mors important file Grid.cl with the existing kernels

---

# Gemini-3-flash

This is a classic high-performance computing problem: **Many-small-dense-systems**.

Since your matrix size is $N \approx 200-400$ (for a molecule like PTCDA) and you want to solve it for thousands of pixels in parallel, the strategy depends heavily on whether you need the **full interaction** (multiple reflections between tip and sample) or if the **single-pass path** (Source $\to$ Tip $\to$ Sample $\to$ Drain) is sufficient.

Here is the analysis of the solver strategy for `Grid.cl` and your GPU implementation.

### 1. The Matrix Dimension vs. GPU Architecture
*   **The Problem:** A $400 \times 400$ complex matrix (`float2`) takes about **1.28 MB** of memory.
*   **The Constraint:** Modern GPU **Shared Memory (`__local`)** is typically limited to **32 KB - 64 KB** per workgroup. 
*   **Conclusion:** You **cannot** fit the full $400 \times 400$ matrix into shared memory. You must either:
    1.  Use **Global Memory** for the matrix (slower but large).
    2.  Use the **Dyson/Inverse-MVM** approach (Pre-calculate the large static parts, solve a tiny coupling system in shared memory).

### 2. Recommended Strategy: The "Small Dyson" Solver
The most efficient way to get "exact" Green's functions including all tip-sample interactions without solving a $400 \times 400$ system every pixel is to pre-invert the static parts.

1.  **Pre-computation (Once per frame/MD step):**
    *   Invert $A_T = (E I - H_{TT} - \Sigma_L)$ and $A_S = (E I - H_{SS} - \Sigma_R)$ on the GPU or CPU. 
    *   These are now dense $G_T$ and $G_S$ matrices.
2.  **Per-Pixel Kernel (Workgroup tasks):**
    *   Identify the **Active Orbitals** ($k_{tip}, k_{smp}$): These are the only orbitals where the hopping $H_{TS}(R)$ is non-zero (e.g., the tip apex and the few surface atoms under it). Usually $k \approx 10-20$.
    *   **Shared Memory:** Load only the sub-blocks of $G_T$ and $G_S$ corresponding to these active orbitals into `__local` memory.
    *   **Solve the Coupling:** The full $G$ is given by the Dyson equation. For the active subspace, you solve a $k \times k$ system (where $k \approx 20$). 
    *   **Matrix-Vector Multiplication (MVM):** Once the coupling is solved, the propagation through the rest of the molecule is just a dense MVM.

### 3. If you want a Full Linear Solver in the Kernel
If you insist on solving the full $(A_{TT} + A_{SS} + H_{TS})x = b$ per pixel, here is how to implement it:

#### The Algorithm: Parallel Gauss-Jordan (Direct)
For $N=400$, iterative solvers (CG, GMRES) might struggle with convergence if the energy $E$ is near an eigenvalue (resonance). A direct solver is more predictable.

*   **Type:** Gauss-Jordan elimination with **Partial Pivoting**.
*   **Parallelization:** 
    *   Each Workgroup = 1 Pixel.
    *   Each Thread = Handles 4-8 rows of the matrix.
*   **Step-by-step in `Grid.cl`:**
    1.  **Pivot Search:** Threads collaborate to find the max element in the current column (using a local reduction).
    2.  **Row Swap:** Threads swap their assigned rows.
    3.  **Division:** One thread (or all) divides the pivot row.
    4.  **Elimination:** Every thread updates its assigned rows in parallel: $Row_i = Row_i - factor \times PivotRow$.
    5.  **Barrier:** `barrier(CLK_LOCAL_MEM_FENCE)` between every column elimination.

#### Complexity Check:
*   $N^3$ operations. For $N=400$, that's $64 \times 10^6$ ops.
*   With 10,000 pixels, you are looking at $0.6$ Tera-ops. 
*   A modern GPU (e.g., RTX 3080/4090) can handle this in **milliseconds** if you saturate the cores.

### 4. Relevant Files and Modules

To implement this, the following files in your repository are critical:

1.  **`Grid.cl`**: This is where the OpenCL kernel lives. You would add a new kernel `kernel_stm_transport_gf`.
2.  **`STM_GPU_QMMM.md`**: Update this with the Green's Function math (Self-energies $\Sigma$ and the Caroli formula).
3.  **`test_response_function_rotated.py`**: This currently calculates overlaps. It should be modified to prepare the $H_{TT}$ and $H_{SS}$ matrices and send them to the GPU.
4.  **`pyopencl` wrapper**: You need a Python side that manages the `cl.Buffer` for the dense Hamiltonians.

### 5. Mathematical Setup for the Kernel

The matrix you are solving in the kernel is:
$$ \mathbf{A} = \begin{pmatrix} E-H_{TT}+i\Gamma_L/2 & -H_{TS}(R) \\ -H_{ST}(R) & E-H_{SS}+i\Gamma_R/2 \end{pmatrix} $$

**The "Path" Vector Propagation ($Ax=b$):**
1.  **Source:** $b$ is a vector with $1.0$ at the "source" orbital of the tip.
2.  **Solve:** $A \psi = b$ (Using the Parallel Gauss-Jordan above).
3.  **Current:** $I = \sum_{j \in \text{sample-drain}} \Gamma_{R,j} |\psi_j|^2$.

### Summary Analysis
*   **Simplicity:** **Pre-inverted MVM** is much simpler and likely 10x faster. It assumes the tip-sample coupling doesn't drastically change the electronic structure of the leads.
*   **Physics:** **Full Gauss-Jordan** per pixel is the most "correct" as it allows the tip and sample to polarize each other to all orders.
*   **Bottleneck:** Memory bandwidth. The matrix $A$ will likely have to stay in **Global Memory** (L2 cache will help). You should use `float2` (complex) and ensure coalesced memory access where threads read consecutive columns.

**Would you like me to outline the specific Workgroup-Parallel Gauss-Jordan logic (how threads divide the 400 rows) for the `Grid.cl`?**

---

# Kimi 2.6

You are right to push back on the “sparse matrix multiply” hand-waving. For $N\sim 400$ the explicit inverse $G(E)$ is dense, and any dense matrix–vector pipeline on a GPU workgroup will choke on memory, not arithmetic. The correct move is **not** to invert the full $400\times 400$ matrix per pixel, but to use the **Woodbury identity** to shrink the linear solve to the size of the **tip** (or the contact patch), which is $\sim 20$–$60$ orbitals. That fits in shared memory, can be Gauss-Jordan-ed by $32$ threads in microseconds, and leaves the heavy $O(N_s^3)$ work to a **single CPU precomputation** that is reused for all pixels.

Below is the full physics-to-GPU mapping, the shared-memory solver design, and the exact files/modules you will touch.

## 1. The Physics: Why the Full System is Overkill

Recall the partitioned retarded Green’s function. The total Hamiltonian is

$$
\mathcal{H}_{\rm tot} = \begin{pmatrix} H_{tt} & H_{ts} \\ H_{st} & H_{ss} \end{pmatrix},
\qquad
S = \begin{pmatrix} S_{tt} & S_{ts} \\ S_{st} & S_{ss} \end{pmatrix}
$$

with local optical potentials (source/drain) on selected atoms:

$$
\Sigma_L = -\frac{i}{2}\Gamma_L,\qquad \Sigma_R = -\frac{i}{2}\Gamma_R
$$

The full retarded propagator is

$$
G^R(E) = \big[(E+i\eta)S - \mathcal{H}_{\rm tot} + \tfrac{i}{2}\Gamma_L + \tfrac{i}{2}\Gamma_R\big]^{-1}
\tag{1}
$$

**Key observation:** $H_{ts}$ is non-zero only for a handful of tip atoms near the surface. If the tip has $n_t$ orbitals and the “contact” sample atoms (those within hopping range) have $n_c$ orbitals, then $H_{ts}$ is an $n_t\times n_c$ **tall-skinny** sparse block. Everything else is precomputable.

### 1.1 Woodbury Reduction

Define the uncoupled (block-diagonal) Green’s functions precomputed on the CPU:

$$
g_t = \big[(E+i\eta)S_{tt}-H_{tt}-\Sigma_L\big]^{-1},\qquad
g_s = \big[(E+i\eta)S_{ss}-H_{ss}-\Sigma_R\big]^{-1}
$$

The Dyson equation for the sample block is

$$
G_s = g_s + g_s\,H_{st}\,g_t\,H_{ts}\,G_s
\;\;\Longrightarrow\;\;
G_s = \big[I - g_s\,H_{st}\,g_t\,H_{ts}\big]^{-1} g_s
\tag{2}
$$

The matrix to invert is $N_s\times N_s$ ($\sim 320\times 320$). But the perturbation
$\Sigma_{\rm contact}\equiv H_{st}g_t H_{ts}$ has **rank $n_t$** (the tip size). Apply the Woodbury identity:

$$
\big[I - g_s\,H_{st}\,g_t\,H_{ts}\big]^{-1}
= I + g_s\,H_{st}\,\mathcal{M}^{-1}\,H_{ts}
\tag{3}
$$

where the **small** matrix is

$$
\boxed{
\mathcal{M} \equiv g_t^{-1} - H_{ts}\,g_s(c,c)\,H_{st}
}
\qquad\text{(size } n_t\times n_t\text{)}
\tag{4}
$$

- $g_t^{-1}=(E+i\eta)S_{tt}-H_{tt}-\Sigma_L$ is known **once per energy** on the CPU.
- $g_s(c,c)$ is the $n_c\times n_c$ contact-contact sub-block of $g_s$, also precomputed.
- $H_{ts}$ is built **per pixel** on the GPU from exponential Slater–Koster overlaps.

For a tip of $5$–$10$ atoms ($n_t=20$–$40$), $\mathcal{M}$ is tiny. Inverting it is the only linear solve you do per pixel.

### 1.2 Transmission from the Reduced System

The Caroli formula is

$$
T(E)=\operatorname{Tr}\!\big[\Gamma_L\,G_{ts}\,\Gamma_R\,G_{st}\big]
\tag{5}
$$

with $G_{ts}=g_t H_{ts} G_s$. Using (2)–(3) and the cyclic trace:

$$
T(E)=\sum_{d\in{\rm drain}}\Gamma_{R,d}\,
\Big[ G_s(c,d)^{\dagger}\,B\,G_s(c,d) \Big]_{dd}
\tag{6}
$$

where

$$
B \equiv H_{st}\,g_t^{\dagger}\,\Gamma_L\,g_t\,H_{ts}
\qquad\text{(}n_c\times n_c\text{, Hermitian)}
\tag{7}
$$

and the contact-to-drain block of the full sample propagator is

$$
G_s(c,d)=g_s(c,d)+g_s(c,c)\,H_{st}\,\mathcal{M}^{-1}\,H_{ts}\,g_s(c,d)
\tag{8}
$$

**All large matrices** ($g_s$, $g_t$, $g_s(c,d)$, $g_s(c,c)$) are **precomputed on CPU** and passed as constant/global buffers.  
**All per-pixel work** reduces to:

1. Build $H_{ts}$ (sparse SK/exponential, $n_t\times n_c$).
2. Matmul: $Z = H_{ts}\,g_s(c,d)$  ($n_t\times n_d$).
3. Matmul: $Q = H_{ts}\,g_s(c,c)\,H_{st}$  ($n_t\times n_t$).
4. Form $\mathcal{M}=g_t^{-1}-Q$.
5. **Invert $\mathcal{M}$ in shared memory** ($n_t\times n_t$).
6. Matmul: $W=\mathcal{M}^{-1}Z$.
7. Matmul: $Y=g_s(c,c)\,H_{st}\,W$.
8. Assemble $G_s(c,d)=g_s(c,d)+Y$.
9. Quadratic forms for each drain atom $d$ using $B$.

Total per-pixel arithmetic: $\sim 0.5$–$2$ million complex ops. At GPU throughput this is **microseconds**, not milliseconds.

## 2. GPU Architecture: One Workgroup = One Pixel

### 2.1 Workgroup & Memory Design

| Resource | Size / Role |
|----------|-------------|
| **Workgroup size** | $T=32$ threads (or $16$–$64$; $32$ is the sweet spot for occupancy and warp alignment) |
| **Shared memory** | $\lesssim 48$–$64$ KB. All small dense blocks live here. |
| **Global/Constant** | $g_t$, $g_t^{-1}$, $g_s$ blocks, $\Gamma$ vectors, geometry tables. Read-only, broadcast to all workgroups. |

**Shared memory budget** (example: $n_t=40$, $n_c=48$, $n_d=40$):

| Buffer | Size | Purpose |
|--------|------|---------|
| `l_H_ts` | $40\times 48\times 8=15.4$ KB | Tip-contact hopping |
| `l_M` | $40\times 40\times 8=12.8$ KB | Woodbury matrix $\mathcal{M}$ |
| `l_Minv` | $40\times 40\times 8=12.8$ KB | Inverse (or reuse `l_M` in-place) |
| `l_Z` | $40\times 40\times 8=12.8$ KB | $H_{ts}g_s(c,d)$ |
| `l_W` | $40\times 40\times 8=12.8$ KB | $\mathcal{M}^{-1}Z$ |
| `l_pivot` | $40\times 4=0.16$ KB | Pivot indices |
| **Total** | **$\sim 67$ KB** | Fits in $64$–$96$ KB shared mem on modern AMD/NVIDIA. If you only have $48$ KB, drop `l_W` and stream it through registers, or reduce $n_t$ to $32$. |

For smaller tips ($n_t\le 32$) the total drops under $45$ KB, which is safe on every GPU.

### 2.2 Why Not Iterative Solvers?

You asked about CG, BiCGSTAB, GMRES, Jacobi.

- **CG**: Requires Hermitian positive-definite. $\mathcal{M}$ is **complex symmetric**, not HPD. Breaks.
- **BiCGSTAB / GMRES / COCG**: Work for complex symmetric, but for $n_t\le 60$ the setup cost (dot products, reductions) and convergence iterations ($\sim 20$–$50$) make them **slower** than a direct dense solve.
- **Jacobi / Gauss-Seidel**: Only useful if strongly diagonally dominant. $\mathcal{M}$ has off-diagonal blocks from $H_{ts}g_s H_{st}$ that are comparable to the diagonal. Unreliable.

**Verdict:** For $n_t\le 64$, use **dense Gauss-Jordan elimination with partial pivoting** in shared memory. It is branch-simple, needs no convergence tolerance, and produces the exact inverse in $O(n_t^3)$ with $n_t^3\sim 10^5$ ops.

### 2.3 The Shared-Memory Complex Solver

Here is the algorithm implemented as an OpenCL device function. It operates on a complex matrix stored as `float2` (`.x`=real, `.y`=imag) in `__local` memory.

```c
// Complex helpers
inline float2 cmul(float2 a, float2 b){
    return (float2)(a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x);
}
inline float2 cdiv(float2 a, float2 b){
    float d = b.x*b.x + b.y*b.y;
    return (float2)((a.x*b.x + a.y*b.y)/d, (a.y*b.x - a.x*b.y)/d);
}
inline float cabs2(float2 a){ return a.x*a.x + a.y*a.y; }

// In-place Gauss-Jordan with partial pivoting
// A is __local float2 A[n*n], pivot is __local int pivot[n]
// T = get_local_size(0), tid = get_local_id(0)
void gauss_jordan_inverse(__local float2* A, __local int* pivot, int n, int tid, int T){
    for(int k=0; k<n; k++){
        // --- 1. Find pivot row (thread 0 only) ---
        if(tid==0){
            float maxv = 0.0f;
            int piv = k;
            for(int i=k; i<n; i++){
                float v = cabs2(A[i*n + k]);
                if(v > maxv){ maxv = v; piv = i; }
            }
            pivot[k] = piv;
        }
        barrier(CLK_LOCAL_MEM_FENCE);

        // --- 2. Swap rows k and pivot[k] ---
        int piv = pivot[k];
        for(int j=tid; j<n; j+=T){
            float2 tmp = A[k*n + j];
            A[k*n + j] = A[piv*n + j];
            A[piv*n + j] = tmp;
        }
        barrier(CLK_LOCAL_MEM_FENCE);

        // --- 3. Normalize pivot row ---
        if(tid==0){
            float2 pivval = A[k*n + k];
            float2 invpiv = cdiv((float2)(1.0f,0.0f), pivval);
            A[k*n + k] = invpiv;
            for(int j=k+1; j<n; j++){
                A[k*n + j] = cmul(A[k*n + j], invpiv);
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);

        // --- 4. Eliminate column k from all other rows ---
        for(int i=tid; i<n; i+=T){
            if(i != k){
                float2 factor = A[i*n + k];
                A[i*n + k] = (float2)(0.0f, 0.0f);
                for(int j=k+1; j<n; j++){
                    float2 prod = cmul(factor, A[k*n + j]);
                    A[i*n + j] = (float2)(A[i*n + j].x - prod.x, A[i*n + j].y - prod.y);
                }
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    // A now contains the inverse (with row permutations from pivoting).
    // Because we only used row swaps, the inverse is correct up to a row permutation.
    // For our use (M^{-1} * Z), row permutation does not matter because we
    // apply the inverse immediately to a matrix/vector on the right.
}
```

**Why Gauss-Jordan and not LU?**
- You need $\mathcal{M}^{-1}$ explicitly to form $\mathcal{M}^{-1}Z$ (step 6).
- Gauss-Jordan gives the inverse directly; LU would require $n_d$ forward/back substitutions.
- For $n_t\le 60$ and $n_d\sim n_t$, the cost is identical, but Gauss-Jordan has simpler control flow in OpenCL (no triangular solves).

**Parallel efficiency:**
- With $T=32$ and $n_t=40$, each thread handles $1$–$2$ rows in the elimination step.
- The inner loop over columns ($j$) is sequential per thread, but $j\le 40$ is tiny.
- Barriers are per-pivot ($n_t$ barriers total). For $n_t=40$, negligible overhead.

## 3. Per-Pixel Kernel Flow

Here is the logical flow of the new kernel. It lives in **Grid.cl** and reuses your existing `evaluate_radial`, `sk_contract_sp`, and quaternion rotation utilities.

```c
__kernel void stm_transmission_negf(
    // --- Precomputed CPU buffers (constant or __global) ---
    __global const float2* g0_t_inv,      // [nt*nt]   (E+iη)S_t - H_t - Σ_L
    __global const float2* g0_s_cc,       // [nc*nc]   contact-contact block of g_s
    __global const float2* g0_s_cd,       // [nc*nd]   contact-drain block of g_s
    __global const float*  gamma_L,       // [nt]      diagonal Γ_L
    __global const float*  gamma_R,       // [nd]      diagonal Γ_R on drain atoms
    
    // --- Geometry (per scan point) ---
    __global const float4* tip_centers,   // [n_points] tip center (x,y,z,0)
    __global const float4* tip_quat,      // [n_points] rotation quaternion
    __global const float4* tip_rel,       // [ntip_atoms] relative positions
    __global const float4* smp_pos,       // [nsmp_atoms] sample positions
    
    // --- Parameters ---
    const float beta, const float r0, const float rcut,
    const int ntip_atoms, const int nsmp_atoms,
    const int nt, const int nc, const int nd,
    
    // --- Output ---
    __global float* out_T                 // [n_points]
){
    const int ip  = get_group_id(0);      // one pixel per workgroup
    const int tid = get_local_id(0);
    const int T   = get_local_size(0);
    
    // --- Shared memory declarations ---
    __local float2 l_H_ts[NT_MAX*NC_MAX];
    __local float2 l_M[NT_MAX*NT_MAX];
    __local float2 l_Z[NT_MAX*ND_MAX];
    __local float2 l_W[NT_MAX*ND_MAX];
    __local int    l_piv[NT_MAX];
    __local int    l_contact_atoms[NC_MAX]; // which sample atoms are in contact
    __local int    l_drain_atoms[ND_MAX];   // which sample atoms are drains
    
    // 1. Determine contact atoms for this pixel (tip position + rcut)
    //    Each thread checks a subset of sample atoms.
    //    Build l_contact_atoms[] and l_H_ts[] using existing SK/exponential code.
    
    // 2. Build H_ts in shared memory
    //    For each tip atom ia, sample atom ja in contact:
    //      r = |p_tip(ia) - p_smp(ja)|
    //      if r < rcut: fill 4x4 SK block in l_H_ts
    
    // 3. Compute Z = H_ts * g0_s_cd   (nt x nd)
    //    Collaborative matmul: each thread computes a subset of (i,j) elements.
    
    // 4. Compute M = g0_t_inv - H_ts * g0_s_cc * H_st
    //    First P = H_ts * g0_s_cc  (nt x nc)
    //    Then Q = P * H_st         (nt x nt)   [H_st = transpose(H_ts) because real basis]
    //    Then M = g0_t_inv - Q     (nt x nt)
    //    Store in l_M.
    
    // 5. Invert M in shared memory
    gauss_jordan_inverse(l_M, l_piv, nt, tid, T);
    // l_M now holds M^{-1}
    
    // 6. Compute W = M^{-1} * Z   (nt x nd)
    
    // 7. Compute Y = g0_s_cc * H_st * W   (nc x nd)
    //    Then G_s_cd = g0_s_cd + Y
    
    // 8. Precompute B = H_st * g0_t^H * Gamma_L * g0_t * H_ts  (nc x nc)
    //    This can be done as:
    //      U = g0_t * H_ts   (nt x nc)   [matmul]
    //      B = U^H * diag(gamma_L) * U   (nc x nc)
    
    // 9. For each drain atom d:
    //      T_d = gamma_R[d] * ( G_s_cd(:,d)^H * B * G_s_cd(:,d) )
    //    Sum to local float T_sum.
    
    // 10. Write out_T[ip] = T_sum
}
```

**Key implementation notes:**

- **Matrix multiplies** for small dense blocks ($\le 64$) are done by having each thread compute one or a few output elements. No fancy tiling is needed; the arithmetic intensity is high enough because everything is in shared memory.
- **Complex conjugate transpose** (`^H`) appears only in step 8 and 9. Since Fireball uses real basis functions, $H_{ts}$ is real; the only complex conjugation is on $g_0$.
- **Contact atom list** (`l_contact_atoms`) is built dynamically per pixel because the tip moves. With $10$ tip atoms and $80$ sample atoms, the all-pairs distance check is $800$ comparisons — trivial for $32$ threads.

## 4. CPU Precomputation (Do This Once per MD Step / Energy)

Before launching the GPU kernel, the host code must:

1. **Assemble $H_{tt}, H_{ss}, S_{tt}, S_{ss}$** from Fireball.
2. **Attach $\Sigma_L, \Sigma_R$** as local imaginary potentials on selected source/drain atoms.
3. **LU-factorize** (or directly invert) the tip and sample Green’s functions:
   $$
   g_t(E) = [(E+i\eta)S_{tt}-H_{tt}-\Sigma_L]^{-1}
   $$
   $$
   g_s(E) = [(E+i\eta)S_{ss}-H_{ss}-\Sigma_R]^{-1}
   $$
   Use LAPACK `zgetrf`/`zgetri` or your own dense solver on CPU. These are $40^3$ and $320^3$ — milliseconds total.
4. **Extract sub-blocks**:
   - $g_s(c,c)$: contact-contact block (depends on which sample atoms are within `rcut` of the tip’s *range*; if the tip scans laterally, the contact set changes. You can either extract a generous superset that covers all possible scan positions, or re-extract per pixel on GPU from the full $g_s$ buffer.)
   - $g_s(c,d)$: contact-to-drain block.
5. **Upload to GPU** as `__constant` or `__global` buffers.

If you are integrating over a bias window with $N_E$ energy points, upload $g_t$ and $g_s$ for all energies. The kernel can loop over energies internally (the geometry $H_{ts}$ is energy-independent), accumulating $T(E)$ and then the current.

## 5. Files & Modules to Implement

| File / Module | What to do |
|---------------|------------|
| **`Grid.cl`** | **Add** `stm_transmission_negf` kernel. **Add** `gauss_jordan_inverse` device function. **Reuse** `sk_contract_sp`, `quat_rotate3`, `evaluate_radial` from existing kernels. |
| **`Grid.cl` (existing)** | `response_amplitude_exp` is the **single-orbital limit** ($n_t=1$, $n_c=1$) of the above formalism. Keep it as a fast path for s-orbital tips; the new kernel supersedes it for multi-orbital tips. |
| **CPU wrapper (Python/C)** | Precompute $g_t, g_s$ via dense complex inversion (NumPy/LAPACK). Extract $g_s(c,c), g_s(c,d)$. Pack $\Gamma_L, \Gamma_R$ diagonal vectors. Launch `clEnqueueNDRangeKernel` with 1D grid = number of pixels, local size = 32. |
| **Shared-memory sizing header** | Define `NT_MAX`, `NC_MAX`, `ND_MAX` (e.g. 64, 64, 64). If you exceed these on a specific molecule, the kernel returns 0 or falls back to a global-memory path. |
| **Tip geometry builder** | Reuse the tip quaternion/position logic from `mo_overlap_points_exp_sk`. The tip is a rigid body; only the center and quaternion change per pixel. |

## 6. Performance Estimate

For a PTCDA-like system ($\sim 80$ sample atoms, $5$–$10$ tip atoms):

| Step | Cost | Time (GPU) |
|------|------|------------|
| Build $H_{ts}$ | $\sim 200$ SK evaluations | $\sim 1\,\mu$s |
| Matmul $Z, Q$ | $\sim 200$K complex ops | $\sim 2\,\mu$s |
| Gauss-Jordan $40\times 40$ | $\sim 60$K complex ops | $\sim 3\,\mu$s |
| Matmul $W, Y$, quadratic forms | $\sim 300$K complex ops | $\sim 3\,\mu$s |
| **Total per pixel** | | **$\sim 10\,\mu$s** |

With $10\,000$ pixels and a GPU that runs $1\,000$–$2\,000$ concurrent workgroups, the wall-clock time is **$\sim 50$–$100\,\mu$s per image** (all pixels in parallel). This easily reaches **thousands of images per second** for MD, assuming the CPU can feed geometries fast enough.

## 7. Summary

- **Do not** build or invert the full $400\times 400$ $G(E)$ per pixel.
- **Do not** use iterative solvers (CG, BiCGSTAB) for a $40\times 40$ matrix; dense Gauss-Jordan in shared memory is faster and exact.
- **Do** use the **Woodbury identity** to reduce the per-pixel solve to the tip dimension $n_t$.
- **Do** precompute the large sample/tip Green’s functions $g_s, g_t$ on CPU once per energy.
- **Do** implement the solver as a collaborative OpenCL workgroup (32 threads) operating entirely in `__local` memory, one workgroup per STM pixel.

The physics is still the full NEGF Caroli formula; the numerics just exploit the fact that the tip–sample coupling is low-rank.

---

# USER

plese sketch the solution - I mean the optimized OpenCL code respecting the constraisn, and efficiently using local memory and other resources, which reaspoable tiling / preloading strategy..., write the code, but also analyse possibilities and pros and cones, use provided information about carefully considering memory layout and representation of basis, with spater costre rotation and other things which you can learn from 

---

# Gemini-3.1-pro

This is a fantastic high-performance computing challenge. Calculating the Green's function transport for a molecule with 200–400 basis functions across thousands of pixels directly on the GPU requires a careful strategy.

### 1. Analysis of the Constraints & The "Dyson Subspace" Trick

**The Memory Constraint:** 
A dense $400 \times 400$ complex matrix (`float2`) takes about **1.28 MB**. The absolute maximum shared memory (`__local`) per workgroup on most GPUs is **32 KB to 64 KB**.
*Conclusion:* You **cannot** load the full matrices into shared memory, and running a $400 \times 400$ Gauss-Jordan elimination in global memory for every single pixel would be far too slow.

**The Physical Optimization (The "Active Subspace"):**
The hopping matrix $H_{TS}(R)$ between the tip and the sample is mostly zeros! It only has non-zero entries for the atoms physically close to each other (e.g., the tip apex and the $\approx 1-4$ sample atoms directly beneath it). 
If we use a cutoff radius $R_{cut}$, the number of "active" orbitals $N_{act}$ is very small. For example, 4 tip atoms and 4 sample atoms = $16 \times 16$ active orbitals. 

Instead of solving the full $400 \times 400$ system $Ax=b$, we use the **Dyson Equation projection**:
1.  **Host/Global:** Precalculate the full isolated Green's functions $G_T$ (Tip) and $G_S$ (Sample). Store them in `__global` memory.
2.  **Kernel/Local:** Identify the small active subsets.
3.  **Kernel/Local:** Extract only the small active blocks of $G_T$ and $G_S$ into `__local` memory (e.g., $16 \times 16$ complex = **2 KB**, easily fitting in shared memory!).
4.  **Kernel/Local:** Solve the multiple-scattering transmission matrix strictly inside this active subspace using a local **Gauss-Jordan solver**.

### 2. The Algorithm per Workgroup (Pixel)

Each workgroup handles **1 pixel (1 tip position)**. The threads (e.g., 64 threads) collaborate:
1.  **Distance Filter:** Threads loop over atoms, measure distances, and build a list of active tip and sample atom indices.
2.  **Preload (Tiling):** Threads collaboratively copy the required sub-blocks of $G_S$ and $G_T$ from global to local memory.
3.  **Build $V_{TS}$:** Threads compute the $4 \times 4$ Slater-Koster overlap blocks for the active pairs.
4.  **Matrix Multiplications:** Compute $W = I - G_S V_{TS}^\dagger G_T V_{TS}$ in local memory.
5.  **Linear Solver:** Run an in-place parallel Gauss-Jordan elimination on $W$ to solve $W x_{S} = b_S$.
6.  **Current Integration:** Calculate the trace/dot product representing the transmission.

### 3. Sketch of the OpenCL Solution (`Grid.cl`)

Here is the structurally optimized kernel demonstrating this strategy.

```c
// Complex arithmetic helpers
inline float2 c_add(float2 a, float2 b) { return (float2)(a.x+b.x, a.y+b.y); }
inline float2 c_sub(float2 a, float2 b) { return (float2)(a.x-b.x, a.y-b.y); }
inline float2 c_mul(float2 a, float2 b) { return (float2)(a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x); }
inline float2 c_div(float2 a, float2 b) {
    float den = b.x*b.x + b.y*b.y + 1e-30f;
    return (float2)((a.x*b.x + a.y*b.y)/den, (a.y*b.x - a.x*b.y)/den);
}

// Define the maximum active subspace size. 
// 32 orbitals = 8 atoms (4 orbitals per atom: px, py, pz, s)
// Fits perfectly in shared memory (32x32 complex matrix = 8 KB)
#define MAX_ACT_ORB 32 

__kernel void solve_stm_dyson_wg(
    const int n_pixels,
    __global const float4* tip_centers,
    __global const float4* tip_pos_rel,
    __global const float4* smp_pos,
    const int ntip_atoms,
    const int nsmp_atoms,
    // Precalculated full Green's functions from Host (Global Memory)
    __global const float2* GT_global, // [4*ntip_atoms * 4*ntip_atoms]
    __global const float2* GS_global, // [4*nsmp_atoms * 4*nsmp_atoms]
    // Incident wave vector injected from the source lead into the tip
    __global const float2* uT_source, // [4*ntip_atoms]
    // Hopping parameters
    const float beta, const float r0, const float rcut,
    // Output
    __global float* out_current
) {
    // 1 Workgroup = 1 Pixel
    const int pixel_id = get_group_id(0);
    const int t_idx    = get_local_id(0); // Thread ID (e.g., 0 to 63)
    const int threads  = get_local_size(0);
    
    if (pixel_id >= n_pixels) return;

    // --- SHARED MEMORY ALLOCATIONS ---
    __local int active_T_atoms[8]; // Max 8 active tip atoms
    __local int active_S_atoms[8]; // Max 8 active sample atoms
    __local int num_act_T, num_act_S;

    // Local Matrices (32x32 floats2 = 8KB each)
    __local float2 GS_loc[MAX_ACT_ORB][MAX_ACT_ORB];
    __local float2 GT_loc[MAX_ACT_ORB][MAX_ACT_ORB];
    __local float2 V_ts[MAX_ACT_ORB][MAX_ACT_ORB]; // Hopping
    __local float2 W[MAX_ACT_ORB][MAX_ACT_ORB]; // The Dyson Matrix to invert
    
    // Local Vectors
    __local float2 uT_loc[MAX_ACT_ORB];
    __local float2 bS_loc[MAX_ACT_ORB]; // Right-hand side

    if (t_idx == 0) { num_act_T = 0; num_act_S = 0; }
    barrier(CLK_LOCAL_MEM_FENCE);

    // 1. DYNAMICALLY IDENTIFY ACTIVE ATOMS (Distance < rcut)
    // Threads cooperatively scan atoms
    const float3 cen = tip_centers[pixel_id].xyz;
    const float rcut2 = rcut * rcut;

    // (Pseudo-code for brevity: Check distances, atomic_inc num_act, save to active_T_atoms)
    // Assume we found N_T active tip orbitals and N_S active sample orbitals.
    int N_T = min(num_act_T * 4, MAX_ACT_ORB);
    int N_S = min(num_act_S * 4, MAX_ACT_ORB);

    // 2. PRELOAD ACTIVE BLOCKS FROM GLOBAL TO SHARED MEMORY
    // Threads loop through the 2D local arrays and fetch from global memory
    for(int i = t_idx; i < N_S * N_S; i += threads) {
        int r = i / N_S; int c = i % N_S;
        int glob_r = active_S_atoms[r/4]*4 + (r%4);
        int glob_c = active_S_atoms[c/4]*4 + (c%4);
        GS_loc[r][c] = GS_global[glob_r * (4*nsmp_atoms) + glob_c];
    }
    // ... similarly preload GT_loc and uT_loc ...
    barrier(CLK_LOCAL_MEM_FENCE);

    // 3. COMPUTE SLATER-KOSTER HOPPING MATRIX V_ts
    // Threads distribute the atom pairs and calculate exponential overlaps
    for(int i = t_idx; i < (N_T/4) * (N_S/4); i += threads) {
        int tip_a = i / (N_S/4);
        int smp_a = i % (N_S/4);
        // ... evaluate f = exp(-beta*(r-r0)) and SK directional cosines ...
        // ... fill the 4x4 block V_ts[tip_a*4 + orbT][smp_a*4 + orbS] ...
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    // 4. MATRIX MULTIPLICATIONS: Construct W = I - GS * V_ts^H * GT * V_ts
    // M1 = GT * V_ts
    // M2 = V_ts^H * M1
    // M3 = GS * M2
    // W  = I - M3
    // (This is highly parallelized. Each thread computes a few elements of the resulting matrices)
    for(int i = t_idx; i < N_S; i += threads) {
        for(int j = 0; j < N_S; j++) {
            W[i][j] = (i == j) ? (float2)(1.0f, 0.0f) : (float2)(0.0f, 0.0f); // Init I
            // Subtract M3 (Pseudo-code)
            // W[i][j] = c_sub( W[i][j], M3[i][j] ); 
        }
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    // 5. IN-PLACE PARALLEL GAUSS-JORDAN SOLVER (W * xS = bS)
    // We solve the linear system for the active sample subspace
    for (int k = 0; k < N_S; k++) {
        // Thread 0 finds pivot (partial pivoting optional for N<32, but good for stability)
        if (t_idx == 0) {
            float2 pivot = W[k][k];
            // Normalize pivot row
            for(int j=k; j < N_S; j++) W[k][j] = c_div(W[k][j], pivot);
            bS_loc[k] = c_div(bS_loc[k], pivot);
        }
        barrier(CLK_LOCAL_MEM_FENCE);

        // All threads eliminate the other rows
        for (int i = t_idx; i < N_S; i += threads) {
            if (i != k) {
                float2 factor = W[i][k];
                for (int j = k; j < N_S; j++) {
                    W[i][j] = c_sub(W[i][j], c_mul(factor, W[k][j]));
                }
                bS_loc[i] = c_sub(bS_loc[i], c_mul(factor, bS_loc[k]));
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    // Now bS_loc contains the exact response wave x_S.

    // 6. CALCULATE FINAL TRANSMISSION CURRENT
    // I = x_S^H * Gamma_R * x_S
    float current = 0.0f;
    if (t_idx == 0) {
        for(int i=0; i<N_S; i++) {
            // Assuming Gamma_R is a simple density-of-states weighting for the drain
            current += (bS_loc[i].x*bS_loc[i].x + bS_loc[i].y*bS_loc[i].y);
        }
        out_current[pixel_id] = current;
    }
}
```

### 4. Pros and Cons of this Strategy

#### Pros:
1.  **Exact Physics:** By directly solving the Dyson equation $W = (I - G_S H_{ST} G_T H_{TS})$, we include **infinite-order multiple scattering** (electrons bouncing back and forth between tip and sample). This perfectly matches the Non-Equilibrium Green's Function (NEGF) / Caroli formula.
2.  **Memory Optimal:** Extracts a small subspace (e.g. $32 \times 32$). The entire calculation fits in the GPU's `__local` L1 shared memory, preventing slow accesses to the huge 1.28MB dense matrices during the solver loop.
3.  **High Throughput:** A $32 \times 32$ Gauss-Jordan elimination distributed over 64 threads takes less than $\sim 500$ clock cycles. You can process thousands of pixels in microseconds.
4.  **No `atomic` Bottlenecks:** Unlike grid-projection density, the entire solving process in this kernel uses deterministic barrier synchronization, meaning zero thread collision.

#### Cons:
1.  **Cutoff Radius Dependency:** If the tip gets extremely close or if $R_{cut}$ is too large, the active subspace might exceed `MAX_ACT_ORB = 32`. 
    *   *Fix:* You would need to add a safe fallback inside the kernel. If `num_act_S > MAX_ACT_ORB`, you process the closest atoms and ignore the furthest, or allocate dynamic shared memory.
2.  **Host Precomputation:** The CPU (or another GPU kernel) must pre-invert $A_T$ and $A_S$ to supply the dense $G_T$ and $G_S$ matrices to global memory before the scan begins.

### 5. Alignment with your Python/Fireball conventions

*   **Basis ordering:** You mentioned in `project_orbital` that Fireball uses `(s, py, pz, px)` but OpenCL uses `(px, py, pz, s)`. The $4 \times 4$ Slater-Koster generation step in the kernel (Step 3) will map natively to `float4`, meaning $G_S$ and $G_T$ matrices generated on the Python side **must** be re-ordered into the `px, py, pz, s` structure before passing to `cl.Buffer`.
*   **Rotations:** If the tip rotates (`tip_quat`), you simply apply `quat_rotate3` to the displacement vector $\vec{d}$ before passing it to the exponential/SK parameter calculations, just as you already do in `mo_overlap_points_exp_sk`.

---

# DeepSeek-V4

# Best Method for GPU‑Accelerated Green’s Function STM Transport

After reviewing all the provided discussions and the Fortran prototype results, the optimal approach for simulating the STM current of a rigid‑tip / fixed‑sample junction using an **OpenCL GPU** is a **workgroup‑per‑pixel Woodbury‑Dyson solver** with a **shared‑memory Gauss‑Jordan elimination**. This method:

- Reduces the per‑pixel problem from a large (\(N \approx 400\)) sparse or dense solve to the inversion of a tiny matrix of size **n_tip_active** (typically 20–60) that fits entirely in fast `__local` memory.
- Uses **one workgroup per pixel** to maximise occupancy; thousands of pixels are solved concurrently.
- Achieves **microsecond** per‑pixel times, enabling hundreds to thousands of full images per second during molecular dynamics.
- Preserves exact symmetry by re‑using Fireball’s angular machinery (`rotate_fb`) inside the GPU kernel for the tip–sample hopping.

The rest of this note details the physical framework, the numerical reduction, the GPU kernel design, memory layout considerations, and how to integrate the Fortran‑validated angular routines.

## 1. Physical System and Green’s Function Transmission

We have two molecules described by Fireball Hamiltonian and overlap matrices \(H_{TT}, S_{TT}\) (tip) and \(H_{SS}, S_{SS}\) (sample). The tip can rotate; the sample is fixed. The source and drain are modelled by constant imaginary self‑energies

\[
\Sigma_L = -\frac{i}{2}\Gamma_L \quad (\text{on selected atoms of the tip}),\qquad 
\Sigma_R = -\frac{i}{2}\Gamma_R \quad (\text{on selected atoms of the sample}),
\]

applied only at the atoms that are physically in contact with the macroscopic leads. The full retarded Green’s function at energy \(E\) is

\[
G^R(E) = \big[ (E+i\eta) S - H_{\text{tot}} + \tfrac{i}{2}\Gamma_L + \tfrac{i}{2}\Gamma_R \big]^{-1},
\qquad 
H_{\text{tot}} = \begin{pmatrix} H_{TT} & H_{TS}(R) \\ H_{ST}(R) & H_{SS} \end{pmatrix},
\quad 
S = \begin{pmatrix} S_{TT} & S_{TS}(R) \\ S_{ST}(R) & S_{SS} \end{pmatrix}.
\]

The tip–sample coupling \(H_{TS}(R)\) is built from the exponential overlaps of Fireball orbitals: \(t_{\mu\nu} = t^{(0)}_{\mu\nu}(r_0) \cdot \exp[-\beta(r-r_0)]\), using the correct Slater–Koster angular structure.

The Landauer–Caroli transmission at energy \(E\) (bias‑window integrated later) is

\[
T(E) = \operatorname{Tr}\big[ \Gamma_L \, G^R_{\,ts}(E) \, \Gamma_R \, G^A_{\,st}(E) \big].
\]

## 2. The Woodbury‑Dyson Reduction

Inverting the full \(N_{\text{tot}} \times N_{\text{tot}}\) matrix per pixel is infeasible on a workgroup’s shared memory. Instead we pre‑compute the **isolated** tip and sample Green’s functions (once per energy) on the CPU:

\[
g_T(E) = \big[ (E+i\eta) S_{TT} - H_{TT} - \Sigma_L \big]^{-1}, \qquad
g_S(E) = \big[ (E+i\eta) S_{SS} - H_{SS} - \Sigma_R \big]^{-1}.
\]

These are dense but small: \(O(200^3)\) for a medium molecule. After pre‑computation, the exact multiple scattering between tip and sample is expressed by the Dyson equation. The key insight is that \(H_{TS}\) is **sparse and short‑ranged**: only a few tip and sample atoms are within the cutoff \(r_{\text{cut}}\). The effective coupling matrix has **rank \(n_t\)** (the number of tip orbitals inside the contact region). Using the Woodbury identity, the transmission can be expressed entirely in terms of small sub‑blocks.

Define:

- \(n_t\) – number of active tip orbitals (tip basis functions of atoms with any sample atom within \(r_{\text{cut}}\))
- \(n_c\) – number of active sample orbitals (sample basis functions that are in range)
- \(n_d\) – number of drain orbitals (sample orbitals attached to \(\Gamma_R\))

The crucial reduction is:

\[
\boxed{
T(E) = \sum_{d=1}^{n_d} \Gamma_{R,d} \; \Big\| \Big( G_s(c,d) + G_s(c,c)\, H_{st} \, \mathcal{M}^{-1} \, H_{ts} \, G_s(c,d) \Big) \, \mathbf{e}_d \Big\|^2_{\Gamma_L}
}
\]

where

- \(G_s(c,c)\) is the \(n_c\times n_c\) block of \(g_S\) corresponding to the contact atoms,
- \(G_s(c,d)\) is the \(n_c\times n_d\) block,
- \(H_{ts}\) is the \(n_t\times n_c\) hopping matrix (real, from Fireball),
- \(H_{st}=H_{ts}^T\),
- \(\mathcal{M} = g_T^{-1} - H_{ts} \, G_s(c,c) \, H_{st}\) is an **\(n_t\times n_t\)** complex matrix,
- \(\| \cdot \|^2_{\Gamma_L}\) denotes the quadratic form \(\mathbf{x}^\dagger \, \Gamma_L^{(t)} \, \mathbf{x}\) where \(\Gamma_L^{(t)}\) is the tip‑side broadening in the active subspace.

All manipulations involve only matrices of size \(\le n_c\), with the only true linear solve being the inversion of \(\mathcal{M}\).

For typical STM settings with a 4‑atom apex tip and 5–10 sample atoms in the contact region, \(n_t \le 40\) and \(n_c \le 40\). The inversion of \(\mathcal{M}\) is the computational heart of one pixel.

## 3. GPU Kernel Design (Workgroup = 1 Pixel)

### 3.1 Workgroup and Memory Layout

- **Workgroup size**: 32 or 64 threads. 32 is often optimal for occupancy and warp/wavefront alignment.
- **Shared memory (\(__local\))**: Holds all small, pixel‑dependent matrices:
  - `H_ts` – \(n_t \times n_c\) real, could be packed as `float` (or `float2` if complex needed later). Since the basis is real, the hopping is real; we can store as `float` and cast when needed.
  - \(g_T^{-1}\) or the tip block – but we actually precompute \(g_T\) and its inverse is available on the CPU, so we can send the \(n_t\times n_t\) matrix \(g_T^{-1}\) directly to the kernel. Alternatively, send the full precomputed sub‑block of \(g_T\) (inverse of tip active block).
  - \(G_s(c,c)\) – \(n_c\times n_c\) complex.
  - \(G_s(c,d)\) – \(n_c\times n_d\) complex.
  - \(\mathcal{M}\) – \(n_t\times n_t\) complex (in‑place Gauss‑Jordan, so reused).
  - Temporary buffers for matrix multiplications.

With \(n_t, n_c \le 40\), the shared memory footprint is well under 48 KB, which is available on virtually all GPUs.

- **Global / Constant memory**: Pre‑computed \(g_T\) (full), \(g_S\) (full), \(\Gamma\) vectors, geometry data, and by‑species basis information for rotation.

### 3.2 Per‑Pixel Workflow

1. **Distance filtering**: Threads cooperatively find all tip–sample atom pairs with distance \(< r_{\text{cut}}\). Count active tip and sample orbitals. The list of active atoms determines the local index mapping into the global orbital numbering.

2. **Preload sub‑blocks**:
   - From the global buffers of \(g_T\) and \(g_S\), each thread copies the required \(n_t\times n_t\), \(n_c\times n_c\), \(n_c\times n_d\) blocks into shared memory. Because \(n_c, n_d\) are small, we can load in a loop over rows/cols assigned to threads (no sophisticated tiling needed).

3. **Build \(H_{ts}\) in shared memory**:
   - For each active tip atom–sample atom pair, compute the displacement vector \(\mathbf{d}\) (rotated into the tip’s molecular frame using the quaternion). 
   - Use Fireball’s radial tables to obtain the reference coupling at \(r_0\), then multiply by \(\exp[-\beta(r-r_0)]\).
   - Apply `rotate_fb` (or equivalent `rotate_fb_matrix_sp` from `hamiltonian.cl`) to construct the full \(4\times 4\) angular block in the laboratory frame.
   - This step is the most physics‑intensive, but for ~10–20 atom pairs it costs only a few hundred Slater–Koster evaluations, executed in parallel by the workgroup.

4. **Matrix multiplications** (all in shared memory):
   - Compute \(\tilde{G}_s = G_s(c,c) H_{st}\) (\(n_c \times n_t\))
   - Compute \(Q = H_{ts} \tilde{G}_s\) (\(n_t \times n_t\))
   - Form \(\mathcal{M} = g_T^{-1} - Q\).
   - Compute \(Z = H_{ts} G_s(c,d)\) (\(n_t \times n_d\)).

   Each matrix multiply is parallelised by assigning each thread a small subset of output elements (e.g., column‑wise or block‑wise). Since everything is in shared memory, no global memory accesses occur during the multiplication loops.

5. **Solve \(\mathcal{M} \, W = Z\)**: compute \(\mathcal{M}^{-1} Z\).
   - A dense Gauss‑Jordan elimination with partial pivoting is applied in‑place on the augmented matrix \([\mathcal{M} | Z]\).
   - The algorithm works directly on the complex values stored as `float2`.
   - Only the pivot row operations need to be done sequentially (by thread 0 or a subgroup), but elimination of the remaining rows is fully parallelised across threads.

6. **Assemble the contact‑to‑drain propagator**:
   \[
   G_cd = G_s(c,d) + G_s(c,c) H_{st} W.
   \]
   (This is just another small matrix multiply.)

7. **Compute transmission**:
   - For each drain orbital \(d\), extract the column vector \(\mathbf{x}_d\) from \(G_cd\).
   - Compute the inner product \(\mathbf{x}_d^\dagger \, \Gamma_L^{(t)} \, \mathbf{x}_d\). Since \(\Gamma_L^{(t)}\) is diagonal (or a small sparse block) this reduces to a sum.
   - Multiply by \(\Gamma_{R,d}\) and accumulate into `out_T[pixel_id]`.

## 4. Key Implementation Details

### 4.1 Basis Ordering and `float4`

Fireball’s internal ordering is usually `s, py, pz, px`, but our OpenCL kernels already use `float4` to represent `(px, py, pz, s)` (see `project_orbital`). The Green’s function matrices must be reordered on the CPU **before** uploading to the GPU, or the rotation routine must be adapted. The latter is easier: keep all global matrices in the `(s, py, pz, px)` order used by Fireball, and implement a small conversion in the Slater–Koster block assembly (just permute indices). That preserves compatibility with existing Fireball rotation routines that expect `s, py, pz, px`.

### 4.2 Rotation Machinery

The Fortran prototype demonstrated that ad‑hoc Slater–Koster blocks break rotational symmetry. Therefore, the GPU kernel must reuse Fireball’s angular functions. In the existing codebase, the file `pyBall/FireballOCL/cl/hamiltonian.cl` already contains `rotate_fb_matrix_sp`, which applies the rotation matrix to a \(4\times 4\) block given the Euler angles or rotation matrix of the tip atom. We should invoke that function after building the distance‑dependent scaling, exactly as done in the Fortran subroutine.

**Tip rotation per pixel**: The tip is a rigid body. The rotation quaternion is passed to the kernel. For each active tip atom, we compute the rotation matrix once (or use pre‑computed per‑atom lab positions already rotated by the host), and then pass the rotation matrix to `rotate_fb_matrix_sp`. This can be done in a loop parallelised over atom pairs.

### 4.3 The Gauss‑Jordan Kernel Function

A robust implementation of in‑place complex Gauss‑Jordan with partial pivoting is provided in the earlier Gemini‑3.1‑pro response. The essential features are:

- Column‑wise outer loop; the pivot row is found by a reduction over threads (or by thread 0 in sequential regions).
- Row swap and normalisation executed by thread 0 (or a single thread) after a barrier.
- Elimination of other rows is fully data‑parallel: each thread handles a subset of rows.
- The function directly overwrites the input matrix with its inverse (or a row‑permuted inverse). Because we only need \(\mathcal{M}^{-1} Z\), we place \(Z\) as additional columns and let Gauss‑Jordan solve the linear system simultaneously.

For \(n_t \le 40\), the whole solve requires fewer than \(n_t^3 \approx 64\,000\) complex operations and completes in a few microseconds.

### 4.4 Pre‑computation on the CPU

Once per MD step (or per energy point), the host code must:

- Run Fireball SCF for tip and sample molecules to obtain \(H_{TT}, S_{TT}\) and \(H_{SS}, S_{SS}\).
- Add \(\Sigma_L, \Sigma_R\) as imaginary diagonal matrices on the drain/source atoms.
- Invert the dense complex matrices using LAPACK (`zgetrf` + `zgetri`). Tip matrix \(\approx 40\times 40\) is trivial; sample \(\approx 320\times 320\) takes a few milliseconds.
- Reorder the basis if necessary, extract the full \(g_T, g_S\) into `__global` buffers.

For multi‑energy integration, all energy points can be uploaded; the kernel loops over energies internally (the geometry‑dependent hopping is energy‑independent), accumulating the current.

## 5. Performance Estimation

Assuming:

- \(n_t = 32\), \(n_c = 40\), \(n_d = 40\).
- 32 threads per workgroup.
- GPU capable of running ~2000 concurrent workgroups.

Operations per pixel:
- Build \(H_{ts}\) – ~500 Slater‑Koster evaluations, ~2000 FLOPs.
- Matrix multiplies – ~\(10^5\) complex FLOPs.
- Gauss‑Jordan on \(32\times 32\) – ~\(6\times 10^4\) FLOPs.
Total arithmetic ~\(2\times 10^5\) FLOPs. Memory accesses are almost all local. At 2‑5 TFLOPS, a single pixel takes below \(100\) ns pure arithmetic; the actual wall time is dominated by synchronisation and loading from global to local, giving ~1–5 µs per pixel.

With \(10\,000\) pixels, the kernel launch covers all workgroups. Even a modest GPU can sustain thousands of such workgroups in flight. The limiting factor becomes the CPU throughput for geometry updates, but the kernel itself can produce >2000 images per second for typical scan sizes.

## 6. Code Integration Path

New OpenCL kernel: add to `Grid.cl` (or a new `transport.cl`). Reuse:

- `quat_rotate3` from existing geometry utilities.
- `rotate_fb_matrix_sp` from `hamiltonian.cl`.
- Radial evaluation (probably a small `evaluate_radial_exp` function for the exponential tail, plus table lookups for reference values).

CPU side (Python):

- Extend `pyBall/FireballOCL/STM.py` to call the new kernel.
- Use the shared memory size definitions `MAX_ACT_ORB = 64` etc.
- Upload the pre‑computed Green’s function blocks, gamma vectors, and per‑species basis data (needed for rotation).

Testing should start with the Fortran reference data: compare the GPU kernel output against the validated Fortran results for a few scan points to ensure correctness of the Slater–Koster rotation and the matrix logic.

## 7. Summary of Advantages

- **Exact NEGF** with infinite multiple scattering (no perturbative weak‑coupling approximation).
- **Memory optimal**: tiny shared‑memory direct solve avoids the \(O(N^3)\) global‑memory inversion.
- **Massively parallel**: one pixel per workgroup, leveraging GPU’s many‑core architecture.
- **Symmetry‑preserving**: reuses Fireball’s own angular machinery, avoiding spurious artifacts.
- **Millisecond‑speed MD compatible**: thousands of images per second feasible.

This design synthesises the best ideas from the earlier discussions and the practical Fortran experience, and is ready for implementation in pyOpenCL.


---

# Implementation Notes: Fortran Dyson STM Reference (2026-04-24)

## Overview

Implemented a CPU-based, non-parallel Fortran reference for the Dyson Green's function STM method using real Fireball Hamiltonians and overlaps for both tip and sample molecules. The implementation uses Fireball's native angular functions (via `rotate_fb`) and replaces only the radial part of two-center integrals with an exponential decay function to extend the range for tunneling coupling.

## Files Modified

### Fortran Core
- **`fortran/MAIN/libFireCore.f90`**
  - Added `firecore_stm_gf_2mol_mo_2points` subroutine (lines ~2100-2230)
  - Implements 2-molecule weak-coupling Dyson scan with:
    - Real tip and sample Hamiltonians/overlaps from Fireball
    - Rotated tip molecule orbitals via MO coefficient vectors
    - Per-pixel tunneling coupling using Fireball's two-center machinery
    - Uses `interpolate_1d`, `recover_2c`, `rotate_fb` for proper angular structure
    - Exponential radial scaling: `exp(-beta*(r-r0))` applied after interpolating at reference distance `r0`
  - Interface accepts real atomic species arrays (`tip_atypes`, `smp_atypes`) for proper shell metadata access

### Python Bindings
- **`pyBall/FireCore.py`**
  - Added Python ctypes binding `stm_gf_2mol_mo_2points` (lines ~488-521)
  - Handles input validation and array conversions
  - Matches new Fortran signature with species arrays

### Test Script
- **`tests/pyFireball/test_stm_gf_dyson_2mol.py`**
  - Main test script for 2-molecule Dyson STM
  - Runs SCF for tip and sample molecules (cached to avoid recomputation)
  - Rotates tip orbitals using Fireball's rotation routines
  - Computes Green's functions for isolated subsystems
  - Calls Fortran Dyson routine for scan over tip positions
  - Generates 6-panel matplotlib visualization:
    1. Tip orbital (real-space projection)
    2. Sample orbital (real-space projection)
    3. Overlap STM (perturbative)
    4. Fortran GF current (Dyson)
    5. Sample orbital (rotated to lab frame)
    6. Tip orbital (rotated to lab frame)
  - Removed log10 scale from GF current panel to show raw structure

## How to Run

```bash
cd /home/prokophapala/git/FireCore
bash make.sh  # Rebuild Fortran library
python3 tests/pyFireball/test_stm_gf_dyson_2mol.py \
    --xyz_tip ../../cpp/common_resources/xyz/CH2O.xyz \
    --xyz_smp ../../cpp/common_resources/xyz/PTCDA.xyz \
    --n 80 \
    --ztip 3.0 \
    --zmid 1.5 \
    --size 20.0 \
    --nmax_scf 30 \
    --tip_mo_list 2 \
    --smp_mo_list 70 \
    --outdir export/stm_gf_dyson_2mol_fortran_ref4
```

## Results

The corrected implementation produces STM images that qualitatively resemble molecular orbitals with reasonable symmetry. The generated images show:

- Proper nodal structure consistent with the molecular orbitals
- Symmetry preservation consistent with tip and sample orbital symmetry
- Current patterns that follow the expected lobe structure

Example output:
- `tests/pyFireball/export/stm_gf_dyson_2mol_fortran_ref4/tip_ch2o__smp_ptcda/fortran_dyson2mol_tipCH2O_MO002_smpPTCDA_MO070_zmid1.50_ztip3.00_tipz000_smpz000.png`

## Angular Symmetry Implementation Lessons

### Critical Issue: Ad-Hoc Slater-Koster Blocks
Initial implementation used manual if-then-else blocks to construct orbital coupling matrices (Vss, Vsp, Vpp_sig, Vpp_pi) with arbitrary scaling. This approach:

- **Broke symmetry**: Manual angular handling did not preserve proper rotational invariance
- **Was fragile**: Required special cases for each orbital combination
- **Was unmaintainable**: Difficult to extend to d-orbitals or more complex basis sets

### Correct Approach: Reuse Fireball's Angular Machinery
The corrected implementation uses Fireball's canonical two-center integral assembly path:

1. **`interpolate_1d`**: Get radial values from pre-computed tables at reference distance `r0`
2. **Exponential rescaling**: Apply `exp(-beta*(r-r0))` to extend range for tunneling
3. **`recover_2c`**: Reconstruct 2D orbital matrix from 1D list using mu/nu lookup tables
4. **`rotate_fb`**: Rotate matrix from molecular to crystal coordinates using `twister` and `chooser`

This ensures:
- **Correct angular dependence**: Fireball's rotation machinery handles s, p, d orbitals properly
- **Rotational invariance**: Symmetry is preserved through proper rotation matrices
- **Extensibility**: Same code works for any basis set supported by Fireball

### Initialization Requirements
The corrected path requires Fireball to be initialized in the active process:

- `readinfo()`: Read species information from info.dat
- `make_munu()`: Build orbital index mappings
- `readdata_mcweda()`: Load two-center integral tables

Without this initialization, arrays like `index_max2c` are unallocated, causing runtime errors.

## Future Work

### Rigorous Testing Needed
The current implementation should be tested across:
- Different molecules (H2O, CH4, benzene, larger organics)
- Different orbital types (s, p, d if available)
- Different tip-sample geometries
- Different transport energies (HOMO, LUMO, mid-gap)

### pyOpenCL Reimplementation
The Fortran reference should be reimplemented in pyOpenCL for GPU acceleration. Key considerations:

- **Reuse angular logic**: Port `rotate_fb` and `twister` to OpenCL kernels
- **Parallel scan**: Each workgroup handles one pixel
- **Shared memory**: Use for small coupling matrices (contact region only)
- **Batch processing**: Process multiple pixels simultaneously for MD integration

### Relevant Files for pyOpenCL Port
- **`pyBall/FireballOCL/STM.py`**: Existing STM utilities (overlap calculations)
- **`pyBall/FireballOCL/cl/hamiltonian.cl`**: Contains `rotate_fb_matrix_sp` for Slater-Koster rotation
- **`fortran/ROTATIONS/rotate.f90`**: Reference implementation of `rotate_fb`
- **`fortran/INTERPOLATERS/interpolate_1d.f90`**: Reference for radial interpolation
- **`fortran/INTERPOLATERS/recover_2c.f90`**: Reference for orbital matrix reconstruction

## Files to Add to Git

Do not forget to add the following new/modified files to git:

1. **`tests/pyFireball/test_stm_gf_dyson_2mol.py`** - Main test script for 2-molecule Dyson STM
2. **`fortran/MAIN/libFireCore.f90`** - Updated with new `firecore_stm_gf_dyson_2mol_mo_2points` subroutine
3. **`pyBall/FireCore.py`** - Updated with new Python binding

The test script `tests/pyFireball/test_stm_dyson_fortran.py` appears to be a different/older test and may not need to be added if it's not part of this implementation.

---

# Deep Analysis: Fortran Green's Function STM Implementation

## 1. Fortran Implementation Structure

**Location:** `fortran/MAIN/libFireCore.f90:2102-2217`

**Subroutine:** `firecore_stm_gf_2mol_mo_2points`

### Algorithm Overview

The Fortran implementation uses the **weak-coupling approximation** for the Dyson equation:

```
G_ts ≈ G_tip * M_ts * G_sample
```

Where:
- `G_tip`, `G_sample`: Pre-computed Green's functions (passed in as `GT_global`, `GS_global`)
- `M_ts`: Tip-sample coupling matrix (computed per point using Slater-Koster)

### Per-Point Algorithm

For each scan point `ip`:

```
1. Initialize M_ts = 0 (tip_norb × smp_norb)

2. Build M_ts (tip-sample coupling):
   For each tip orbital it and sample orbital is:
     a. Get atom indices: ia = tip_orb2atom(it), ja = smp_orb2atom(is)
     b. Compute positions: pt = cen + tip_pos_rel(ia), ps = smp_pos(ja)
     c. Compute distance: d = ps - pt, r = |d|, r2 = r²
     d. Apply cutoff: skip if r² > rcut² or r² < 1e-16
     e. Compute exponential: f = exp(-beta*(r-r0))
     f. Compute direction unit vector: sighat = d/r
     g. Call epsilon(r2v, sighat, eps) → 3×3 rotation matrix
     h. Retrieve Slater-Koster parameters via interpolate_1d
     i. Expand to orbital block via recover_2c → 4×4 hmol, smol
     j. Rotate orbital block via rotate_fb → hrot, srot
     k. Extract element: M_ts(it,is) = -hrot[orb_it, orb_is]

3. Green's function propagation (weak coupling):
   tmpS = GS_global * c_smp           # Sample response
   tmpT = M_ts * tmpS                 # Tip-sample coupling
   tmpTG = GT_global * tmpT           # Tip response
   amp = c_tip^H * tmpTG             # Overlap with tip MO
   out_current = |amp|²              # Current
```

## 2. Slater-Koster Rotation Mechanism

### Key Fortran Subroutines

**`recover_2c(in1, in2, hlist, hbox)`** ([fortran/INTERPOLATERS/recover_2c.f90](cci:7://file:///home/prokop/git/FireCore/fortran/INTERPOLATERS/recover_2c.f90:0:0-0:0)):
- **Input:** `hlist` (1D array of Slater-Koster parameters for atom pair in1-in2)
- **Output:** `hbox` (numorb_max × numorb_max orbital block matrix)
- **Purpose:** Expands 1D parameter list to 4×4 orbital block (s, px, py, pz)
- **Mapping:** Uses `mu(index,in1,in2)` and `nu(index,in1,in2)` to map list indices to orbital positions

**`rotate_fb(in1, in2, eps, mmatrix, xmatrix)`** ([fortran/ROTATIONS/rotate.f90](cci:7://file:///home/prokop/git/FireCore/fortran/ROTATIONS/rotate.f90:0:0-0:0)):
- **Input:** `mmatrix` (orbital block in molecular frame), `eps` (3×3 rotation matrix from direction cosines)
- **Output:** `xmatrix` (rotated orbital block in lab frame)
- **Purpose:** Rotates orbital block from molecular frame to lab frame
- **Algorithm:**
  ```fortran
  call twister(eps, dmat, pmat)  # Build d-orbital (5×5) and p-orbital (3×3) rotation matrices
  
  For each shell issh in atom1:
    l1 = angular momentum (s=0, p=1, d=2)
    call chooser(l1, dmat, pmat, left)  # Get rotation matrix for this shell
    
    For each shell jssh in atom2:
      l2 = angular momentum
      call chooser(l2, dmat, pmat, right) # Get rotation matrix for this shell
      
      # Rotate the block: xmatrix = left @ mmatrix @ right^T
      For each orbital m1 in shell1, m2 in shell2:
        For each k1 in shell1, k2 in shell2:
          xmatrix(n1+k1,n2+k2) += left(k1,m1) * mmatrix(n1+m1,n2+m2) * right(k2,m2)
  ```

**`epsilon(r2v, sighat, eps)`**:
- **Input:** `r2v` (position vector), `sighat` (direction unit vector)
- **Output:** `eps` (3×3 rotation matrix)
- **Purpose:** Computes rotation matrix from direction cosines (bond direction to lab frame)

### Orbital Conventions (Ortega/Fireball)

From [rotate.f90](cci:7://file:///home/prokop/git/FireCore/fortran/ROTATIONS/rotate.f90:0:0-0:0) comments:
```
P-shell:   py   pz   px
           1    2    3

D-shell:   xy    yz   z^2  xz   x^2-y^2
           1     2   3     4      5
```

This is the **Fortran/Fireball order**: `[s, py, pz, px]` for p-orbitals.

## 3. Comparison with Python Rotation ([test_response_function_rotated.py](cci:7://file:///home/prokop/git/FireCore/tests/pyFireball/test_response_function_rotated.py:0:0-0:0))

### Python Approach: Dense Matrix Rotation

**`rotate_dense_hs(H, S, norb_per_atom, R)`** (`pyBall/FireballOCL/STM.py:1178`):
```python
# Build full U = block_diag(U_0, U_1, ...)
U = np.zeros((norb, norb), dtype=np.float64)
for ia in range(natoms):
    ni = int(norb_per_atom[ia])
    Ui = build_atom_rotation_matrix(ni, R)
    U[i0:i0+ni, i0:i0+ni] = Ui
    
# Rotate: H' = U @ H @ U^T, S' = U @ S @ U^T
H_rot = U @ H @ U.T
S_rot = U @ S @ U.T
```

**Key difference:**
- **Python:** Rotates the **entire dense Hamiltonian** as a rigid body
- **Fortran (GF):** Does NOT rotate the Hamiltonian; instead rotates the **coupling matrix M_ts** per point using Slater-Koster direction cosines

### Why Fortran Uses Per-Point Slater-Koster Rotation

The Fortran GF implementation:
1. **Does NOT rotate the tip Hamiltonian** - it uses pre-computed `GT_global` (tip Green's function)
2. **Rotates tip positions** via `tip_pos_rel` (passed in from Python)
3. **Rotates MO coefficients** via `c_tip` (passed in from Python)
4. **Computes M_ts per point** using Slater-Koster rotation based on the **instantaneous direction** between each tip orbital and sample orbital

This is fundamentally different from the Python `rotate_dense_hs` approach:
- Python: "Rotate the molecule, then compute coupling"
- Fortran GF: "Keep molecule fixed, rotate the coupling computation based on relative positions"

## 4. GPU Reimplementation Strategy

### Current Fortran Limitations

1. **Per-point Slater-Koster computation** - expensive for many points
2. **Sequential point processing** - no parallelism across points
3. **CPU-only** - no GPU acceleration
4. **Weak coupling only** - no full Dyson resummation

### GPU Implementation Options

**Option A: Per-Pixel Slater-Koster on GPU (Current Fortran approach on GPU)**
- Each workgroup = 1 pixel
- Compute M_ts in parallel using Slater-Koster rotation
- Apply weak-coupling formula
- **Pros:** Matches Fortran exactly, simple
- **Cons:** Still O(N_orb²) per pixel, expensive for large molecules

**Option B: Pre-computed Green's Functions + Small Dyson (Design discussion in STM_GF_new.md)**
- Pre-compute `G_tip` and `G_sample` (as Fortran does)
- For each pixel:
  - Compute only **active orbitals** (those within rcut)
  - Build small coupling matrix M_active (k×k, where k ≈ 10-20)
  - Solve small Dyson equation: `(I - G_tip_active * M_active * G_sample_active)^(-1)`
  - Extract transmission
- **Pros:** Much faster (k×k vs N_orb×N_orb), exact (not weak coupling)
- **Cons:** More complex, requires active orbital detection

**Option C: Linear System Solver (Impulse Response)**
- Build sparse matrix `A = (E+iη)S - H + iΓ/2`
- Factorize once (LU decomposition)
- For each pixel: solve `A ψ = b` where `b` is source vector
- **Pros:** Exact, handles strong coupling
- **Cons:** Requires sparse LU on GPU (complex), Fortran doesn't do this

### Recommended GPU Approach: Option B with Per-Pixel Slater-Koster

**Kernel structure:**
```opencl
__kernel void stm_gf_dyson_2mol(
    // Pre-computed Green's functions (dense)
    __global const float2* GT_global,  // (tip_norb, tip_norb)
    __global const float2* GS_global,  // (smp_norb, smp_norb)
    
    // MO coefficients
    __global const float2* c_tip,
    __global const float2* c_smp,
    
    // Atomic data
    __global const float4* tip_pos_rel,
    __global const float4* smp_pos,
    __global const int* tip_atypes,
    __global const int* smp_atypes,
    __global const int* tip_orb2atom,
    __global const int* smp_orb2atom,
    
    // Slater-Koster lookup tables
    __global const float* sk_h_table,  // (nZ, nZ, nSK_params)
    __global const float* sk_s_table,
    
    // Scan points
    __global const float4* points,
    
    // Parameters
    float beta, r0, rcut,
    
    // Output
    __global float* out_current
)
```

**Per-workgroup algorithm:**
1. Load Green's functions into `__local` memory (if small enough)
2. For each point (assigned to workgroup):
   - Compute M_ts using Slater-Koster (parallel over orbital pairs)
   - Apply weak coupling: `tmpS = GS * c_smp`, `tmpT = M_ts * tmpS`, `tmpTG = GT * tmpT`
   - Compute amplitude: `amp = c_tip^H * tmpTG`
   - Output `|amp|²`

### Slater-Koster Rotation on GPU

**Challenge:** `rotate_fb` uses nested loops over shells and orbitals

**GPU implementation:**
- Pre-compute rotation matrices for p and d orbitals from direction cosines
- Use shared memory for small orbital blocks (4×4)
- Parallelize over orbital pairs within workgroup

**Orbital order mapping:**
- Fortran: `[s, py, pz, px]`
- GPU kernel: Need to match this order exactly for parity

## 5. Comparison with [test_response_function_rotated.py](cci:7://file:///home/prokop/git/FireCore/tests/pyFireball/test_response_function_rotated.py:0:0-0:0)

### What [test_response_function_rotated.py](cci:7://file:///home/prokop/git/FireCore/tests/pyFireball/test_response_function_rotated.py:0:0-0:0) Does

1. **Rotates dense Hamiltonian:** `H_rot = U @ H @ U^T`
2. **Rotates MO coefficients:** `C_rot = C @ U^T`
3. **Rotates atomic positions:** `apos_rot = R @ apos`
4. **Rotates scan grid:** `grid_rot = R @ grid`
5. **Tests rotational invariance:** Unrotated system + rotated grid ≈ rotated system + unrotated grid

### What Fortran GF Does

1. **Does NOT rotate Hamiltonian** - uses pre-computed `GT_global`, `GS_global`
2. **Rotates MO coefficients** - Python passes in rotated `c_tip`
3. **Rotates atomic positions** - Python passes in rotated `tip_pos_rel`
4. **Does NOT rotate scan grid** - uses fixed lab-frame points
5. **Computes M_ts per point** - Slater-Koster rotation based on instantaneous direction

### Key Difference

**Python rotation test:**
- Validates that rotating the **entire system** (H, S, C, positions) gives the same physics

**Fortran GF:**
- Relies on Python to rotate positions and MO coefficients
- Does NOT rotate the Hamiltonian (uses pre-computed G)
- Computes coupling on-the-fly using Slater-Koster based on relative positions

### For Proper Tip Rotation in Fortran GF

To properly rotate the tip molecule in the Fortran GF approach, you need to:

1. **Rotate tip positions in Python:**
   ```python
   tip_pos_rot = rotate_points(tip.apos, R_tip)
   tip_pos_rot = tip_pos_rot - tip_pos_rot.mean(axis=0)  # Center at origin
   ```

2. **Rotate tip MO coefficients in Python:**
   ```python
   C_tip_rot = rotate_mo_coefficients(C_tip0.T, tip_scf['norb_per'], R_tip).T
   ```

3. **Pass rotated data to Fortran:**
   - `tip_pos_rel`: Rotated tip positions (relative to tip center)
   - `c_tip`: Rotated MO coefficients
   - `GT_global`: **Unrotated** tip Green's function (or recompute if needed)

**Important:** The Fortran GF does NOT rotate the Hamiltonian itself. It assumes the Green's functions are pre-computed and valid for the rotated geometry. If you rotate the tip significantly, you may need to recompute `GT_global` with the rotated positions.

## 6. Summary

| Aspect | Fortran GF Implementation | Python Rotation Test |
|--------|-------------------------|---------------------|
| Hamiltonian rotation | None (uses pre-computed G) | Rotates dense H via U @ H @ U^T |
| MO coefficients | Rotated in Python, passed to Fortran | Rotated in Python |
| Positions | Rotated in Python, passed to Fortran | Rotated in Python |
| Coupling computation | Per-point Slater-Koster rotation | Per-point Slater-Koster rotation |
| Green's function | Pre-computed, weak coupling | Pre-computed, weak coupling |
| Rotation philosophy | "Rotate positions, recompute coupling" | "Rotate entire system, test invariance" |

**For GPU implementation:** Use the Fortran approach (per-point Slater-Koster with pre-computed G), but parallelize across pixels and implement the Slater-Koster rotation efficiently in OpenCL.

---

# pyOpenCL GPU Implementation of `firecore_stm_gf_2mol_mo_2points` (2026-04-24)

## Overview

A pyOpenCL GPU implementation of the 2-molecule Green's Function STM scan. The GPU kernel computes the weak-coupling chain `amp = c_tip^H · GT · M_ts · GS · c_smp` for thousands of pixels in parallel using a simple work-item-per-pixel kernel with simplified exponential Slater-Koster hopping.

## Files Modified

### `pyBall/FireballOCL/cl/Grid.cl` (+95 lines)

- **New kernel `stm_gf_dyson_2mol_mo_scan`** (line 1366)
  - Work-item-per-pixel (single thread, no workgroup parallelism needed)
  - Pre-computed vectors `u_T = c_tip^H · GT` and `v_S = GS · c_smp` are passed in
  - Kernel loops over all orbital pairs accumulating `amp = Σ u_T[it] · V(it,is) · v_S[is]`
  - Simplified exponential SK: `Vss=Vsp=Vps=f`, `Vpp_sig=f`, `Vpp_pi=0.2f` where `f = exp(-beta*(r-r0))`
  - Orbital convention: OpenCL Grid order `[px,py,pz,s]` (same as all existing Grid.cl kernels)
  - Output: `|amp|²` as float per pixel

### `pyBall/FireballOCL/Grid.py` (+145 lines)

- **New method `stm_gf_dyson_2mol_mo_scan`** (line 868)
  - Accepts GT/GS/c_tip/c_smp in **Fortran convention**
  - Pre-computes on CPU: `v_S = GS @ c_smp`, `u_T = c_tip^H @ GT`
  - Remaps from Fortran order `[s,py,pz,px]` to OpenCL Grid order `[px,py,pz,s]` via permutation `[3,1,2,0]`
  - Handles H atoms: 1 orbital in Fortran → padded to 4 in OCL (px=py=pz=0, s in slot 3)
  - Builds OCL-padded `orb2atom` arrays for kernel
  - Pack float4/complex buffers and launches kernel with `(n_pixels,)` global size

### `tests/pyFireball/test_stm_gf_dyson_2mol_ocl.py` (new file)

Test script that:
- Runs SCF for tip (CH₂O) and sample (PTCDA) molecules
- Pre-computes Green's functions at transport energy `E = 0.5*(E_tip + E_smp)`
- Rotates tip geometry and MO coefficients
- Computes three STM images for each MO pair:
  1. **Fortran GF** reference via `fc.stm_gf_2mol_mo_2points`
  2. **GPU GF** via `gp.stm_gf_dyson_2mol_mo_scan` (simplified SK)
  3. **Overlap STM** via `gp.mo_overlap_points_exp_sk_2mol` (reference)
- 6-panel figure: sample MO, tip MO, overlap STM, Fortran GF, GPU GF, log-log scatter plot
- Each intensity panel uses its own vmin/vmax (not log scale)

Example run:
```bash
cd tests/pyFireball
python test_stm_gf_dyson_2mol_ocl.py --n 80 --nmax_scf 30 --size 20 --ztip 3 --zmid 1.5 --rcut 8 --beta 1.5 --tip_mo_list 4,5 --smp_mo_list 72
```

## Critical Convention Handling

### Orbital Order Remapping: Fortran → OpenCL Grid

```
Fortran per-atom:   [s,  py, pz, px]   (Ortega/Fireball spherical harmonics)
OpenCL Grid order:  [px, py, pz, s]    (Cartesian, used in all Grid.cl kernels)
Permutation:        [3,  1,  2,   0]   (Fortran[idx] → OCL[perm[idx]])
```

H atoms (1 orbital in Fortran): mapped to `[0, 0, 0, s]` in OCL (only slot 3 populated).

### Green's Function Matrix Convention

The Green's function matrices `GT_global` and `GS_global` come from `np.linalg.inv(z*S - H)` where H and S are in Python row-major storage (= Fortran column-major transpose). Since H and S are symmetric, the transpose is identity and the convention is correct.

### Column-Major vs Row-Major

Fortran writes matrices column-major; Python reads them row-major. For symmetric matrices this is invisible. The Fortran function receives pointers to Python's row-major arrays and interprets them as column-major. Since both `GT_global` and `GS_global` are complex-symmetric (inverse of symmetric matrices), the storage ordering introduces no error.

## Performance

| Grid Size | Fortran (CPU) | GPU (RTX 3090) | Speedup |
|-----------|---------------|----------------|---------|
| 80×80     | 2.0 s         | 3.0 ms         | ~670×   |
| 100×100   | 3.1 s         | 3.0 ms         | ~1000×  |

GPU time is nearly constant (dominated by kernel launch overhead, ~3 ms) regardless of grid size.

## Spatial Correlation with Fortran Reference

| Height | Linear corr | log-space corr | GPU↔Overlap |
|--------|-------------|----------------|-------------|
| z=2.0 | 0.30 | 0.47 | 0.56 |
| z=3.0 | 0.33 | 0.43 | 0.60 |
| z=4.0 | 0.36 | 0.50 | 0.83 |
| z=5.0 | 0.33 | 0.50 | 0.85 |

### Findings

1. **Correlation ~0.3–0.5**: GPU GF captures the same major orbital lobes and nodal structure as Fortran, but differs in fine spatial details (relative peak heights between atom sites).
2. **GPU↔Overlap rises with height**: At larger z, the GF propagator becomes less important and GF & simple overlap converge. At close z (≤2.5Å), GF intra-molecular propagation matters more.
3. **Ratio F/GPU varies with z** (0.34–0.77) because the simplified SK parameters differ from the species-dependent Fireball 2c integrals used by Fortran.
4. **CPU replica confirms kernel correctness**: A CPU-side implementation of the same algorithm (same SK, same vectors) exactly matches the GPU output.

### Root Cause of Discrepancy with Fortran

The Fortran function `firecore_stm_gf_2mol_mo_2points` builds `M_ts` using Fireball's actual **2-center integral tables** (`interpolate_1d(13,...)` at r0, rotated via `recover_2c`/`rotate_fb`). The GPU kernel uses **simplified uniform exponential SK** (`Vss=Vsp=Vps=Vpp_sig=1, Vpp_pi=0.2`). The species-dependent hopping ratios (C-C vs C-O vs C-H) differ from 1:1:1:1, causing different spatial intensity maps.

## How to Match Fortran Exactly

To reproduce Fortran results exactly, the GPU kernel would need to:

1. **Precompute per-species-pair reference SK blocks** on CPU using Fireball's `interpolate_1d` at `r0` (callable from Python via `FireCore` bindings)
2. **Upload reference blocks** to GPU as a small lookup table `[n_species][n_species][16]` (4×4 matrix elements per pair)
3. **In the kernel**, look up species pair, scale by `exp(-beta*(r-r0))`, apply SK rotation via direction cosines

This preserves the ~1000× speedup while matching the Fortran reference exactly.

## Images Location

```
tests/pyFireball/export/stm_gf_dyson_2mol_ocl/tip_ch2o__smp_ptcda/
```

Each PNG: 6 panels (2×3) — sample MO, tip MO, overlap STM, Fortran GF, GPU GF, scatter plot. All intensity panels use their own vmin/vmax (no log-scale).

---

# TODO: Achieve Exact Fortran Parity in GPU Implementation

## Critical Missing Piece: Direction-Dependent Slater-Koster Rotation

The GPU kernel currently uses a **simplified isotropic exponential hopping** (`Vss=Vsp=Vps=Vpp_sig=1, Vpp_pi=0.2`) which completely omits the **direction-dependent rotation of p and d orbitals**. This is the dominant source of the ~0.3–0.5 correlation with the Fortran reference. Without direction rotation, the STM images are physically wrong for any molecule with p-orbitals (which is almost all molecules).

### What Fortran Does (Correct)

For each tip orbital–sample orbital pair at each scan point:

1. **Interpolate radial integral tables**: `interpolate_1d(13, 0, in1, in3, index, 0, r0, href, sref)` → species-dependent 1D parameter list `hlist_ref`
2. **Expand to orbital block**: `recover_2c(in1, in3, hlist_ref, hmol)` → `numorb_max × numorb_max` block in molecular frame
3. **Compute direction cosines**: `epsilon(r2v, sighat, eps)` → 3×3 rotation matrix from bond direction to lab frame
4. **Rotate orbital block**: `rotate_fb(in1, in3, eps, hmol, hrot)` → block in lab frame
   - Internally calls `twister(eps, dmat, pmat)` to build p-orbital (3×3) and d-orbital (5×5) rotation matrices
   - Shell-wise contraction: `xmatrix = left @ mmatrix @ right^T`
5. **Scale by distance decay**: `hrot * exp(-beta*(r-r0))`
6. **Extract element**: `M_ts[it,is] = -hrot[orb_it, orb_is]`

### What GPU Currently Does (Wrong)

1. **No species dependence**: all atom pairs use the same hopping strength
2. **No direction rotation**: isotropic decay only, no `epsilon` / `twister` / `rotate_fb`
3. **Orbital order remapping in kernel**: `[s,py,pz,px] → [px,py,pz,s]` — adds confusion, should be done in Python wrapper

### Required Fixes

#### 1. Implement `epsilon` in OpenCL
- Input: direction unit vector `sighat = d / r`
- Output: 3×3 rotation matrix `eps`
- Reference: Fortran `epsilon` subroutine (find exact location in `fortran/`)
- This is the bond-direction → lab-frame transformation

#### 2. Implement `twister` in OpenCL
- Input: `eps` (3×3 rotation matrix)
- Output: `dmat` (5×5 d-orbital rotation), `pmat` (3×3 p-orbital rotation)
- Reference: `fortran/ROTATIONS/rotate.f90:113` → `call twister(eps, dmat, pmat)`

#### 3. Implement `rotate_fb` in OpenCL
- Input: species pair `(in1, in2)`, `eps`, `mmol` (4×4 block in molecular frame)
- Output: `hrot` (4×4 block in lab frame)
- Algorithm:
  - Call `twister(eps, dmat, pmat)`
  - For each shell pair `(issh, jssh)`:
    - `l1 = lssh(issh, in1)`, `l2 = lssh(jssh, in2)`
    - `left = chooser(l1, dmat, pmat)` (identity for s, pmat for p, dmat for d)
    - `right = chooser(l2, dmat, pmat)`
    - Contract: `hrot += left @ mmol @ right^T`
- Reference: `fortran/ROTATIONS/rotate.f90:74-146`

#### 4. Pre-compute SK Reference Blocks on CPU
- Use Fireball's `interpolate_1d(13, 0, in1, in3, index, 0, r0, ...)` via Python bindings
- Expand via `recover_2c` logic (can be done in Python — just mapping `mu(index,in1,in2)` / `nu(index,in1,in2)`)
- Store as lookup table: `SK_ref[in1][in3][4][4]` (float16 per species pair)
- Upload to GPU as `__constant` or `__global` buffer

#### 5. Update Kernel Algorithm
```opencl
// Per pixel, per orbital pair (it, is):
d = ps - pt;
r = length(d);
if (r < 1e-6f) { sighat = (float3)(0,0,1); }
else { sighat = d / r; }

// 1. Look up reference block for species pair
hmol = SK_ref[tip_atype][smp_atype];  // 4×4, already in molecular frame

// 2. Rotate by direction cosines
eps = epsilon(sighat);
twister(eps, dmat, pmat);
hrot = rotate_fb(tip_atype, smp_atype, eps, hmol);

// 3. Scale by distance decay
f = exp(-beta * (r - r0));
hrot = hrot * f;

// 4. Extract element
M_ts[it,is] = -hrot[orb_it, orb_is];
```

#### 6. Remove Orbital Order Remapping from Kernel
- The GPU kernel should use **Fortran orbital order** directly: `[s, py, pz, px]`
- Any remapping should happen in the Python wrapper before buffer upload, not inside the kernel
- This eliminates a source of indexing bugs

#### 7. H Atom Handling
- H has only 1 orbital (s) in Fortran
- In the GPU, we should either:
  - Keep 1 orbital per H and handle indexing carefully, OR
  - Pad to 4 but ensure only slot 0 (s) is populated and the other 3 are never accessed
- Current padding `[0,0,0,s]` with remapping `[3,1,2,0]` is fragile; better to match Fortran order exactly

### Files to Modify

1. **`pyBall/FireballOCL/cl/Grid.cl`**
   - Add `epsilon`, `twister`, `chooser`, `rotate_fb` functions
   - Rewrite `stm_gf_dyson_2mol_mo_scan` kernel to use full SK rotation
   - Remove simplified isotropic hopping

2. **`pyBall/FireballOCL/Grid.py`**
   - Add Python method to pre-compute `SK_ref` blocks via Fireball `interpolate_1d`
   - Upload `SK_ref` buffer to GPU
   - Fix orbital order: pass data in Fortran order, no kernel-side remapping

3. **`tests/pyFireball/test_stm_gf_dyson_2mol_ocl.py`**
   - Pass `SK_ref` to GPU wrapper
   - Verify parity: correlation with Fortran should jump from ~0.3–0.5 to >0.99
   - Run parity sweep across multiple tip/sample MOs, heights, and rotations

### Expected Result

After implementing direction-dependent SK rotation with species-dependent hopping tables:
- **Spatial correlation** with Fortran reference: **>0.99** (currently ~0.3–0.5)
- **Performance**: ~1000× speedup preserved (3 ms vs 2 s for 80×80 grid)
- **Physical correctness**: p-orbital anisotropy correctly captured, nodal planes and lobe orientations match Fortran exactly