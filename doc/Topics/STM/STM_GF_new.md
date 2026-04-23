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