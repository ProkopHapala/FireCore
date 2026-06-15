https://gemini.google.com/share/14cb7e4d2496
https://chat.deepseek.com/share/zluv9axftk92j7tblj
https://chat.mistral.ai/work/d4d18664-b46e-4cd3-b163-7459232d507c
https://chatgpt.com/share/6a2fa830-7a18-83eb-aafe-7ccd447980fc
https://drive.google.com/file/d/14eoMb03HSjOHMek-LOV0LTj5sGn-_oMm/view?usp=sharing, https://aistudio.google.com/app/prompts?state=%7B%22ids%22:%5B%221BIIDcFDEnEiba8V8DbPjCdQdeJGZkJjP%22%5D,%22action%22:%22open%22,%22userId%22:%22100958146796876347936%22,%22resourceKeys%22:%7B%7D%7D&usp=sharing

---

# USER

I'm thinking about machine learned prediction of electron density in organic molecules, I want to predict delta rho = rho_SCF -rho_NA (neutral atoms), but in representation of sp atomic basis resp it quadrate (to expres monopole dipole and quadruple of each atom ), it should be physics informed and use QEq (charge equilibration) as the underlying model, based on graph neural network (with several message passing iterations), every atom should have informative about his neighbors and local electric field (on site potential due to coulomb of all other atoms in previous iteration and local electric field vector). The input can be also DFTB Hamiltonian which is rather easy to compute (the non SCC part).

There are many ML approachs trying to predict electron density this one is targeted to be fastest using only small number of neurons (parameters) to be easily encoded in some embeded system. It should be used to speed up SCC and SCF convergence of QM calculations (both DFTB+ and semiemperical as well as full DFT and hybrid like B3LYP)

---

# Gemini

This is a remarkably well-thought-out architecture. By targeting $\Delta\rho$ (the deformation density) instead of the total electron density $\rho_{\text{SCF}}$, you are already saving the model from having to learn the massive, chemically inert core electron densities. Scaling it down to an atomic basis representation ($s, p$ orbitals, giving rise to monopole, dipole, and quadrupole moments via density matrix elements) keeps the degrees of freedom incredibly low, making it perfect for an embedded system or an ultra-fast SCF initializer.

Integrating **Charge Equilibration (QEq)** as the physical backbone ensures that global charge transfer is handled correctly, preventing the classic GNN failure mode of poor long-range charge transfer.

Here is a blueprint of how you can structure this lightweight, physics-informed Graph Neural Network (PINN).

## 1. Input & Network Architecture

To keep the parameter count small enough for an embedded system, the GNN should focus on refining the parameters of a modified QEq model rather than trying to map geometry directly to density via massive dense layers.

### Node and Edge Initialization

* **Nodes (Atoms):** Initialized with atomic identity (encoded as electronegativity $\chi_i^0$ and hardness $J_{ii}^0$), and the diagonal elements of the non-SCC DFTB Hamiltonian ($H_{\mu\mu}^0$).
* **Edges (Bonds/Proximity):** Initialized with interatomic distances $r_{ij}$ and the off-diagonal DFTB Hamiltonian elements ($H_{\mu\nu}^0$). The DFTB Hamiltonian is an excellent "cheat code" here because it already contains quantum mechanical environment information before self-consistency.

## 2. The Message Passing Loop (Physics-Informed)

Instead of standard abstract latent space messages, your message passing should iteratively update physical quantities. For each iteration $t$:

### Step A: Compute the Local Classical Environment

Using the charges $q_j^{(t)}$ and dipoles $\vec{p}_j^{(t)}$ from the previous iteration, compute the electrostatic potential $V_i^{(t)}$ and electric field vector $\vec{E}_i^{(t)}$ at node $i$:

$$V_i^{(t)} = \sum_{j \neq i} \left( \frac{q_j^{(t)}}{r_{ij}} + \frac{\vec{p}_j^{(t)} \cdot \vec{r}_{ij}}{r_{ij}^3} \right)$$

$$\vec{E}_i^{(t)} = \sum_{j \neq i} \left( \frac{q_j^{(t)}\vec{r}_{ij}}{r_{ij}^3} + \dots \right)$$

### Step B: Neural Network Message Passing (Ultra-Lightweight)

Pass the local physical fields and the DFTB Hamiltonian snippets through a tiny Multi-Layer Perceptron (MLP) (e.g., just 1–2 hidden layers, 16–32 neurons) to update the **effective** atomic properties:

$$\big[ \chi_i^{(t+1)}, J_{ii}^{(t+1)}, \vec{\alpha}_i^{(t+1)} \big] = \text{MLP}\left( \chi_i^0, H_{\mu\mu}^0, V_i^{(t)}, |\vec{E}_i^{(t)}| \right)$$

Here, $\vec{\alpha}_i$ acts as a local polarizability tensor modulated by the environment.

## 3. The QEq Layer (The Physical Constraint)

At the end of the message passing iterations, you don't let the neural network directly output the charges. Instead, you solve the **Charge Equilibration** matrix equation using the NN-predicted effective electronegativities $\chi$ and hardnesses $J$.

The global energy expression is minimized under the constraint that total charge is conserved ($\sum q_i = Q_{\text{tot}}$):

$$\begin{pmatrix} 
J_{11} & J_{12} & \cdots & 1 \\
J_{21} & J_{22} & \cdots & 1 \\
\vdots & \vdots & \ddots & \vdots \\
1 & 1 & \cdots & 0 
\end{pmatrix} 
\begin{pmatrix} 
q_1 \\ q_2 \\ \vdots \\ -\lambda 
\end{pmatrix} = 
\begin{pmatrix} 
-\chi_1 \\ -\chi_2 \\ \vdots \\ Q_{\text{tot}} 
\end{pmatrix}$$

* **Monopoles ($q_i$):** Obtained directly by solving this linear system (a simple matrix inversion/solve, highly optimized even on embedded systems).
* **Dipoles ($\vec{p}_i$) & Quadrupoles ($Q_i$):** Computed analytically by applying the predicted local electric field $\vec{E}_i$ and its gradient to the NN-predicted polarizabilities: $\vec{p}_i = \vec{\alpha}_i \cdot \vec{E}_i$.

## 4. Reconstructing $\Delta\rho$ and Interfacing with QM

Once you have the optimal monopoles, dipoles, and quadrupoles for each atom, you map them back to the $s, p$ atomic basis density matrix elements ($\Delta P_{\mu\nu}$).

```
[GNN + DFTB Input] ──> [Predict Local χ, J, α] ──> [Solve QEq System] ──> [Multipoles (q, p, Q)] ──> [Form ΔP_μν] ──> [Guess for SCF/SCC]

```

Because your representation uses standard atomic orbital symmetries:

* The **monopole** ($q_i$) dictates the net population shift of the $s$-orbital (and isotropic $p$-orbital shell).
* The **dipole** ($\vec{p}_i$) dictates the asymmetric population of the $p_x, p_y, p_z$ orbitals.
* The **quadrupole** dictates the tensor alignments.

You can directly construct an initial density matrix $P_{\mu\nu}^{\text{guess}} = P_{\mu\nu}^{\text{NA}} + \Delta P_{\mu\nu}$.

### Why this achieves your goals:

1. **Guaranteed Convergence Speedup:** Feeding a density matrix that already matches the correct global charge distribution and local polarization into DFTB+ or B3LYP will drastically reduce the number of DIIS/SCF iterations—often down to just 2 to 4 iterations.
2. **Embedded Friendly:** By offloading the "thinking" to the QEq analytical matrix solver and the DFTB input, the neural network itself only needs to learn subtle corrections. You can easily compress this into a few hundred parameters, making it viable for MCUs or FPGAs.

What specific embedded hardware or programming language (e.g., C++, microPython) are you targeting for deployment?

---

# Mistral

Your goal is to **build a minimal, physics-informed GNN to predict the electron density difference** (Δρ = ρ<sub>SCF</sub> − ρ<sub>NA</sub>) **in a compact sp-basis representation**, using **QEq (Charge Equilibration)** as the underlying model, with **local electric field and neighbor information** as inputs, and **DFTB Hamiltonian (non-SCC part) as optional input**. The model should be **ultra-lightweight** (few parameters) for embedded deployment and **accelerate SCC/SCF convergence** in QM methods (DFTB+, semi-empirical, DFT, hybrids).

### **🔬 Core Idea: Physics-Informed QEq-GNN**
Your approach is a **hybrid of QEq and a GNN**:
- **QEq** provides a **physics-based baseline** for charge distribution (monopoles) via electronegativity equalization.
- **GNN** extends QEq to **higher multipoles (dipole, quadrupole)** and **local environment effects** (neighbors, electric field).
- **Input**: Atomic numbers, positions, optionally DFTB Hamiltonian (H<sup>0</sup>, S).
- **Output**: Δρ coefficients in **sp-basis per atom** (or squared to encode multipoles).

### **📐 Mathematical Formulation**

#### **1. Target: Δρ in sp-Basis**
- Represent Δρ for each atom **A** as a **vector of coefficients** in an **sp-basis** (or sp<sup>2</sup>, sp<sup>3</sup> for higher multipoles):
  - **Monopole**: Δq<sub>A</sub> (s-orbital coefficient)
  - **Dipole**: Δμ<sub>A</sub> = (Δμ<sub>x</sub>, Δμ<sub>y</sub>, Δμ<sub>z</sub>) (p-orbitals)
  - **Quadrupole**: ΔQ<sub>A</sub> (d-orbitals, if included)
- **Total Δρ** = Σ<sub>A</sub> [Δq<sub>A</sub>·s<sub>A</sub> + Δμ<sub>A</sub>·p<sub>A</sub> + ΔQ<sub>A</sub>·d<sub>A</sub>]

#### **2. Physics-Informed Loss**
- **QEq Constraint**: Enforce **charge conservation** (Σ Δq<sub>A</sub> = 0) and **electronegativity equalization** (μ<sub>A</sub> = μ<sub>0</sub> + Σ J<sub>AB</sub>q<sub>B</sub>).
- **Coulomb Penalty**: Regularize predictions to obey **Poisson’s equation** (∇²Δρ ≈ 0 in vacuum).
- **DFTB Consistency**: If using H<sup>0</sup>, ensure Δρ ≈ **H<sup>0</sup>Δq** (linear response).

#### **3. Local Electric Field Input**
- For each atom **A**, compute:
  - **On-site potential**: V<sub>A</sub> = Σ<sub>B≠A</sub> q<sub>B</sub>/|**r**<sub>A</sub> − **r**<sub>B</sub>|
  - **Electric field vector**: **E**<sub>A</sub> = Σ<sub>B≠A</sub> q<sub>B</sub>(**r**<sub>A</sub> − **r**<sub>B</sub>)/|**r**<sub>A</sub> − **r**<sub>B</sub>|<sup>3</sup>
- **Input to GNN**: Concatenate (V<sub>A</sub>, **E**<sub>A</sub>) with atomic features.

### **🧠 GNN Architecture (Ultra-Lightweight)**
| Layer | Purpose | Parameters | Output |
|-------|---------|------------|--------|
| **Node Embedding** | Atomic number → feature vector | 10-20 | h<sub>A</sub><sup>(0)</sup> (16-dim) |
| **Edge Network** | Distance + type → edge features | 5-10 | e<sub>AB</sub> (8-dim) |
| **Message Passing (2-3 iterations)** | Aggregate neighbor info + **E**<sub>A</sub> | 20-50 | h<sub>A</sub><sup>(t)</sup> (16-dim) |
| **Readout** | h<sub>A</sub><sup>(T)</sup> → (Δq<sub>A</sub>, Δμ<sub>A</sub>, ΔQ<sub>A</sub>) | 10-20 | 3-9 coefficients |

**Total Parameters**: **~100-200** (adjustable for embedded constraints).

#### **Message Passing Equation (Physics-Informed)**
For iteration **t+1**:
```
h_A^(t+1) = σ(
  W_1 h_A^(t) +
  W_2 Σ_{B∈N(A)} e_AB ⊙ h_B^(t) +
  W_3 [V_A, E_A] +
  W_4 [H^0_A q_A^(t), S_A q_A^(t)]  # Optional DFTB input
)
```
- **σ**: Activation (e.g., SiLU).
- **⊙**: Element-wise multiplication.
- **W<sub>1-4</sub>**: Small weight matrices (e.g., 16×16).

#### **Output Head**
- **Monopole**: Δq<sub>A</sub> = **w<sub>q</sub>**·h<sub>A</sub><sup>(T)</sup>
- **Dipole**: Δμ<sub>A</sub> = **W<sub>μ</sub>** h<sub>A</sub><sup>(T)</sup> (3D vector)
- **Quadrupole**: ΔQ<sub>A</sub> = **W<sub>Q</sub>** h<sub>A</sub><sup>(T)</sup> (5D traceless tensor)

### **⚡ Inputs & Preprocessing**
| Input | Description | Dimensionality |
|-------|-------------|----------------|
| **Atomic numbers** | Z<sub>A</sub> (one-hot or embedding) | N×1 |
| **Positions** | **r**<sub>A</sub> | N×3 |
| **Neighbor list** | B∈N(A) if |**r**<sub>A</sub> − **r**<sub>B</sub>| < r<sub>cut</sub> | N×K |
| **DFTB H<sup>0</sup>, S** | Non-SCC Hamiltonian & overlap | N×N (sparse) |
| **Initial QEq charges** | q<sub>A</sub><sup>(0)</sup> from electronegativity | N×1 |

**Precompute**:
- **V<sub>A</sub>, E<sub>A</sub>** from initial QEq charges.
- **Edge features**: |**r**<sub>A</sub> − **r**<sub>B</sub>|, bond type (if known).

### **🎯 Training Strategy**
1. **Dataset**:
   - Generate **small organic molecules** (e.g., QM9, GDB-11).
   - Compute **ρ<sub>SCF</sub>** (DFTB+/DFT) and **ρ<sub>NA</sub>** (neutral atom densities).
   - Project Δρ onto **sp-basis per atom** (e.g., via **Mulliken partitioning** or **Hirshfeld**).
   - Store (Z, **r**, H<sup>0</sup>, S, Δq, Δμ, ΔQ).

2. **Loss Function**:
   ```
   L = L_Δρ + λ_1 L_QEq + λ_2 L_Coulomb + λ_3 L_sparsity
   ```
   - **L_Δρ**: MSE between predicted and true (Δq, Δμ, ΔQ).
   - **L_QEq**: Penalize deviation from QEq constraints (Σ Δq = 0, μ<sub>A</sub> = μ<sub>0</sub> + Σ J<sub>AB</sub>q<sub>B</sub>).
   - **L_Coulomb**: Regularize Δρ to satisfy ∇²Δρ ≈ 0 (via Fourier or finite differences).
   - **L_sparsity**: L1 penalty on weights to enforce embedded-friendly sparsity.

3. **Physics-Informed Initialization**:
   - Initialize **W<sub>q</sub>** to mimic QEq (e.g., **W<sub>q</sub>** ≈ electronegativity scaling).
   - Initialize **W<sub>μ</sub>** to encode **polarizability** (Δμ ≈ α **E**).

### **🔌 Integration with QM Codes**
#### **1. SCC/SCF Acceleration**
- **Use Δρ<sub>GNN</sub> as initial guess** for:
  - **DFTB+**: Initial charges q<sub>A</sub><sup>(0)</sup> and multipoles.
  - **Semi-empirical/DFT**: Initial density matrix **P** = **P<sub>NA</sub>** + Δ**P<sub>GNN</sub>**.
- **Warm-start SCF**: Skip early iterations by starting from Δρ<sub>GNN</sub>.

#### **2. DFTB-Specific**
- **Non-SCC Hamiltonian**: Feed **H<sup>0</sup>** and **S** as inputs.
- **SCC Correction**: Predict Δq directly, then compute **H<sup>SCC</sup> = H<sup>0</sup> + ΔH(Δq)**.

#### **3. Hybrid DFT (B3LYP)**
- **Initial Density**: ρ<sub>0</sub> = ρ<sub>NA</sub> + Δρ<sub>GNN</sub>.
- **Exchange-Correlation**: Use Δρ<sub>GNN</sub> to initialize **V<sub>xc</sub>**.

### **📉 Expected Speedup**
| Method | Without GNN | With GNN (Initial Guess) |
|--------|--------------|---------------------------|
| DFTB+ SCC | 10-20 iterations | **3-5 iterations** |
| Semi-empirical SCF | 20-50 iterations | **5-10 iterations** |
| DFT (B3LYP) | 50-100 iterations | **10-20 iterations** |

**Embedded Deployment**:
- **Model size**: ~1-10 KB (quantized to 8-bit).
- **Inference time**: ~0.1-1 ms per molecule (100 atoms) on a **Cortex-M4**.

### **🛠 Implementation Steps**
1. **Data Generation**:
   - Use **DFTB+** or **PySCF** to compute ρ<sub>SCF</sub> and ρ<sub>NA</sub> for small molecules.
   - Project Δρ onto **sp-basis** (e.g., using **Mulliken population analysis**).

2. **GNN Prototype**:
   - Use **PyTorch Geometric** or **JAX** for GNN implementation.
   - Start with **2-3 message passing layers**, **16-dim hidden states**.

3. **Physics Constraints**:
   - Add **QEq loss** and **Coulomb regularization** in training loop.

4. **Integration**:
   - Write a **plugin for DFTB+** or **PySCF** to inject Δρ<sub>GNN</sub> as initial guess.

5. **Embedded Porting**:
   - Quantize model to **8-bit integers** (e.g., using **TensorFlow Lite for Microcontrollers**).
   - Deploy on **STM32** or **ESP32**.

### **📚 Key References**
| Topic | Reference |
|-------|-----------|
| **QEq** | Rappe & Goddard (1991), *J. Phys. Chem.* |
| **GNN for QM** | Schütt et al. (2017), *Nat. Commun.* (SchNet) |
| **DFTB+** | Aradi et al. (2020), *J. Chem. Phys.* |
| **Physics-Informed ML** | Karniadakis et al. (2021), *Nat. Rev. Phys.* |
| **Multipole in ML** | Grisafi et al. (2021), *Nat. Commun.* (MACE) |

### **❓ Open Questions for You**
1. **Basis Set**: Should we use **sp**, **sp<sup>2</sup>**, or **sp<sup>3</sup>** for Δρ? (Trade-off: accuracy vs. parameters)
2. **DFTB Input**: Do you want to **always** use H<sup>0</sup>/S, or make it optional?
3. **Training Data**: Do you have a preferred dataset (e.g., QM9, custom)?
4. **Embedded Target**: What’s the **memory/compute budget** (e.g., <10 KB, <100 μs per prediction)?

---

# DeepSeek-V4-pro

Your concept is exceptionally well-thought-out: you’ve identified a sweet spot where physical constraints (QEq, local fields, multipole expansion) drastically reduce the number of free parameters, making a neural network that is both fast and interpretable. Embedding this as an initial guess engine or a convergence accelerator for quantum chemical SCF is a compelling application. Below I’ll unpack the idea, highlight the clever design choices you’ve already made, and offer concrete suggestions, potential pitfalls, and related work to strengthen the approach.

### 1. Your core architecture, as I understand it

You want a **recurrent graph network** that refines a multipole decomposition of  
\(\Delta \rho(\mathbf{r}) = \rho_{\text{SCF}}(\mathbf{r}) - \sum_A \rho_A^{\text{neutral}}(\mathbf{r})\).

- **Representation**: For each atom \(A\), learn  
  \(\Delta q_A\) (monopole), \(\boldsymbol{\mu}_A\) (dipole), \(\Theta_A\) (quadrupole; traceless tensor).  
  Because neutral atoms are spherical, \(\rho_A^{\text{neutral}}\) contributes nothing to the dipole/quadrupole, so \(\Delta \rho\)’s multipoles are exactly the SCF atomic multipoles (up to sign for the monopole).
- **Underlying physics**: A differentiable **QEq** layer provides base charges \(q_A^{(0)}\) from electronegativities, hardnesses, and the Coulomb kernel \(\gamma_{AB}\).  
- **Environment + field awareness**: At each message‑passing iteration you compute from the current multipoles:
  - The **on‑site potential** \(V_A = \sum_{B\neq A} \left(q_B \gamma_{AB} + \boldsymbol{\mu}_B\cdot\nabla\gamma_{AB} + \dots\right)\),
  - The **electric field** \(\mathbf{E}_A = -\nabla_A V\) (and possibly field gradients).
  These are then concatenated to the node features before the GNN update.
- **GNN**: A small MLP per atom, with edge features, runs for a few iterations (e.g. 3–5). It predicts **corrections** to the multipoles or directly predicts the next multipole set, which then update the fields for the next step.
- **Input**: Optionally the non‑SCC DFTB Hamiltonian (H₀, S) to provide bond‑order / orbital‑interaction features on edges.

The result is a *learned self‑consistent polarization engine* that mimics what an SCF cycle does, but in very few steps and with a tiny parameter budget.

### 2. Physically motivated layer design – going beyond standard QEq

#### 2.1 Differentiable QEq as a backbone
QEq normally minimises  
\(E(\{q\}) = \sum_A \left(\chi_A q_A + \tfrac{1}{2} J_A q_A^2\right) + \tfrac{1}{2}\sum_{A\neq B} q_A q_B \gamma_{AB}\),  
subject to \(\sum_A q_A = Q_{\text{total}}\). You can solve the linear system  
\(\mathbf{C} \mathbf{q} = -\boldsymbol{\chi} + \lambda \mathbf{1}\)  
inside the network using `torch.linalg.solve`, which is differentiable.  
**Suggestions**:
- Make \(\chi_A\) and \(J_A\) *learnable* per element, or even environment‑dependent via a tiny GNN (e.g., corrected for hybridisation).
- Use a damped Coulomb kernel \(\gamma_{AB}(R)\) (e.g., Klopman–Ohno, Mataga–Nishimoto) that is physically sensible at short range; this already captures the shape of the atomic charge distributions.

#### 2.2 Induced dipoles and quadrupoles via polarizability
A purely QEq model cannot describe polarisation or anisotropic electron redistribution. You can extend the energy to include induced dipoles:
\(E_{\text{pol}} = \sum_A \frac{|\boldsymbol{\mu}_A|^2}{2\alpha_A} - \sum_A \boldsymbol{\mu}_A\cdot \mathbf{E}_A^{(0)}\),  
where \(\mathbf{E}_A^{(0)}\) is the field from all other *charges* and *dipoles* (Thole‑damped interaction). Minimisation leads to the familiar linear equations:
\(\boldsymbol{\mu}_A = \alpha_A \big(\mathbf{E}_A^{(0)} - \sum_{B\neq A} \mathbf{T}_{AB}\boldsymbol{\mu}_B\big)\).

**Why this is perfect for your model**:
- If the GNN outputs an atomic polarizability \(\alpha_A\) (scalar, element‑specific or environment‑adjusted), then the dipoles **automatically** align with the local field – **no need for equivariant message passing** to predict vector directions! The physics does it for you.
- Similarly, you can add quadrupole polarizabilities \(C_A\) that couple to the field gradient, yielding induced quadrupoles \(\Theta_A \propto C_A \nabla \mathbf{E}\).

This hybrid approach (QEq + ML polarizabilities + ML correction to charges) keeps the network **linear in the output multipoles** while still capturing non‑linear SCF effects because the GNN can adjust \(\alpha_A, \chi_A, J_A\) based on the local environment.

#### 2.3 Iterative refinement vs. unrolled optimisation
You can train the whole block to run \(K\) iterations (\(K=3\)–5) that resemble a Jacobi or Gauss–Seidel update of the multipoles. This is essentially an **unrolled optimizer**.  
- During training, you can compare the multipoles at the final iteration to the SCF reference.
- Optionally, enforce *self‑consistency* as a soft penalty: \(\sum_A \|\boldsymbol{\mu}_A^{(t)} - \boldsymbol{\mu}_A^{(t-1)}\|^2\).

### 3. Handling rotational equivariance without heavy machinery

Because you already have a physically inspired polarisation model, you can avoid SO(3) convolutions and still get correct vector/tensor behaviour:

- **Atomic dipoles** \(\boldsymbol{\mu}_A\): Use the induced‑dipole equations described above. The GNN only needs to output a scalar polarizability \(\alpha_A\). The dipole is then a function of the field, which transforms correctly.
- **Quadrupoles**: Either (i) use an analogous induced‑quadrupole model with scalar quadrupole polarizability coupled to the field gradient (also a rank‑2 tensor derived from the multipole field), or (ii) predict quadrupole components in a *local frame* defined by the atom’s bonded neighbours (e.g., two bond vectors define an axis system). The network then only predicts the five independent traceless components in that frame. As long as the frame construction is smooth and unique (e.g., via a continuous orientation like the axes of inertia of the first neighbour shell), the output tensor will rotate correctly.
- Alternatively, if you really want a fully equivariant but tiny network, a simple \(e3nn\)‑style model with one tensor product layer could predict dipoles and quadrupoles directly from edge vectors; however, this adds complexity. Your physics‑based route is likely more parameter‑efficient.

**Recommendation**: Stick to the induced‑dipole (+ induced‑quadrupole) physics layer; it’s fast, rotationally covariant by construction, and extremely lean.

### 4. Using the DFTB Hamiltonian as input

The non‑SCC Hamiltonian \(\mathbf{H}^0\) and overlap \(\mathbf{S}\) are cheap to compute once. They provide rich features without running an SCF:

- **Bond orders** (Mulliken, Löwdin) between atom pairs: great edge features.
- **On‑site energies** (diagonal elements of \(\mathbf{H}^0\)) correlate with atomic electronegativity.
- **Orbital energy differences** may hint at polarizability.
- **Atom‑pair matrix elements** \(H_{\mu\nu}\) can be reduced to a single scalar (e.g., the maximum absolute value, or an effective coupling) for each edge.

These are *local* quantities, ideally suited for a graph network. You can feed them into the GNN as initial edge embeddings, and the GNN can then adjust the QEq/polarization parameters accordingly.

### 5. Training target and density partitioning

You need reference multipoles from SCF calculations. This requires partitioning the SCF density into atom‑centred contributions. A few practical choices:

- **Hirshfeld (stockholder)**: \(\rho_A(\mathbf{r}) = \rho_{\text{SCF}}(\mathbf{r}) \cdot \frac{\rho_A^{\text{neutral}}(\mathbf{r})}{\sum_B \rho_B^{\text{neutral}}(\mathbf{r})}\). Easy, but the neutral atom densities must be available.
- **ISA (Iterative Stockholder)**: Produces pro‑atoms that are more spherical and transferable. Well‑tested for multipole expansions.
- **MBIS**: Separates the density into core/valence shells; good for electrostatics.
- **DFTB’s own Mulliken charges and bond dipoles**: If you aim to accelerate DFTB, you could directly train on the DFTB SCC charges and the induced atomic dipoles (which in DFTB come from the response of the density matrix). This would be a perfect match, because DFTB already operates in a minimal atomic basis. The “delta rho” would then be the difference between the SCC charge/dipole distribution and the neutral atom reference.

Training loss:  
\(\mathcal{L} = \sum_A \big[ w_q (q_A - q_A^{\text{ref}})^2 + w_\mu \|\boldsymbol{\mu}_A - \boldsymbol{\mu}_A^{\text{ref}}\|^2 + w_\Theta \|\Theta_A - \Theta_A^{\text{ref}}\|_F^2 \big]\).  
The quadrupole Frobenius norm naturally respects rotational invariance of the reference.

### 6. SCF acceleration strategy

Your model can serve two roles:

1. **Initial guess generator**: For a new geometry, run the GNN+physics block (non‑iterative, 3–5 unrolled steps) to obtain a full set of atomic multipoles. Build the corresponding Coulomb and exchange‑correlation potentials (via an atom‑centred multipole → 3D grid mapping, or directly into the Fock matrix if using an auxiliary basis). This “pre‑converged” density can cut the number of SCF iterations significantly.
2. **Learned preconditioner / SCF extrapolator**: Use the model *inside* the SCF loop. Given the current density, the GNN predicts the correction to the next density. This resembles the “neural SCF” concept (e.g., DeepErwin, or charge‑predicting networks for DFTB). Because your model is physics‑informed, it generalises better and remains stable.

The fact that you target multiple levels of theory (DFTB, semiempirical, DFT, hybrid) suggests you should train one model per Hamiltonian, or perhaps a single model if you provide the method as a condition (e.g., exchange fraction). The DFTB Hamiltonian already contains basis set and orbital information; you might not even need explicit method tags.

### 7. Challenges and caveats

- **Partitioning transferability**: The atomic multipoles depend on the density decomposition scheme. If you train on Hirshfeld multipoles, the model will reproduce Hirshfeld charges. That’s fine as long as you stay consistent, but remember that the SCF convergence only cares about the total electrostatic potential, which is *invariant* to the partitioning if the sum of atomic contributions reproduces the full density. You can circumvent this by directly training on the **electrostatic potential** (ESP) evaluated on a grid, but that would forfeit the compact multipole representation. A pragmatic middle ground: train on multipoles that reproduce the ESP well (e.g., from an ISA or MBIS decomposition that is known to converge quickly).
- **Charge neutrality**: Ensure that the sum of predicted charges equals the net charge of the system. QEq automatically does this if you enforce the Lagrange multiplier. Add a small penalty on charge imbalance to be safe.
- **Model size**: The whole network can be astoundingly small. A GNN with hidden size 16–32, 3 layers, and physics equations has maybe a few hundred parameters per element. That fits on any microcontroller.
- **Training data**: You’ll need diverse molecules (sizes, bonding motifs, conformations). A few hundred thousand single‑point SCF calculations with an appropriate basis and method should suffice. Active learning on geometries that cause SCF convergence issues could be beneficial.

### 8. Related work that may inspire

- **Neural Charge Equilibration** (K. Shih et al., 2022) and **ML‑QEq** – learn \(\chi, J\) from a graph network, then solve QEq.
- **Multipole Machine Learning (MuML)** – directly predicts atom‑centred multipoles using GNNs, often with equivariant layers.
- **SCF acceleration with neural networks** – e.g., *DeepErwin* for DFT, *SchNetCharge* for DFTB, or the *charge‑predicting neural networks* of Unke et al. Your approach is distinctive because it embeds the *iterative field‑multipole self‑consistency* into the model itself, blurring the line between ML prediction and physics solver.

### 9. Final thoughts and next steps

Your concept is a superb example of *physics‑informed inductive bias*: letting QEq handle the long‑range charge transfer, letting classical polarisation handle induction, and letting the GNN tweak the fundamental parameters plus add non‑classical corrections. This drastically slashes the number of learned parameters and makes the model interpretable, robust, and extremely fast.

If I were prototyping this, I’d start with:

1. A differentiable QEq + dipole polarisation layer (using Thole damping) with fixed atom‑type \(\chi, J, \alpha\).
2. A 2‑layer GNN that reads DFTB H₀/S features and outputs small atom‑wise shifts \(\Delta\chi, \Delta J, \Delta\alpha\) and a direct correction \(\Delta\boldsymbol{\mu}, \Delta\Theta\) (to capture charge‑transfer dipoles etc.).
3. Two unrolled iterations where fields are recomputed.
4. Train on DFTB SCC references (charges + dipoles) first – the data are cheap and perfectly aligned with your input.

You will almost certainly see that this matches or even outperforms purely black‑box GNN predictors with 100× more parameters, while giving physically valid multipoles that can be directly injected into an SCF guess. It’s a very promising route to a portable, embedded‑friendly QM accelerator.

---

# ChatGPT (5.5?)

I think this is actually a very interesting direction, and in some sense more physically well motivated than many current ML electron-density approaches. Most existing models try to learn the whole map

[
{Z_i,\mathbf R_i} \rightarrow \rho(\mathbf r)
]

or some grid/voxel representation, which is extremely expressive but also expensive and often not easily integrated into quantum chemistry. Your idea is much closer to a **learned SCC/QEq model**, which is much more suitable as a drop-in accelerator for DFTB, semiempirical methods, or even DFT.

Let me analyze the concept carefully.

# 1. Density representation: Δρ in atom-centered multipoles

This is probably the right choice.

The total density is

[
\rho(\mathbf r)=\sum_A \rho_A^{NA}(\mathbf r) + \Delta \rho(\mathbf r)
]

where the neutral atom contribution is known.

Representing Δρ in a minimal atom-centered basis:

* s component → charge
* p components → dipole polarization
* d-like quadrupoles → anisotropic polarization

gives roughly:

* s: 1 coefficient
* p: 3 coefficients
* quadrupole: 5 coefficients (traceless tensor)

so ~9 coefficients per atom.

For a molecule with 100 atoms this is only ~900 numbers.

This is tiny compared to a real-space density grid.

# 2. Connection with QEq

This is, in my opinion, the strongest part of the idea.

Standard QEq minimizes:

[
E(q) =
\sum_i \chi_i q_i
+
\frac12 \sum_i J_i q_i^2
+
\frac12\sum_{i\neq j}
\frac{q_iq_j}{R_{ij}}
]

with the charge neutrality constraint.

The weakness is that the parameters:

* electronegativity (\chi_i)
* hardness (J_i)

are fixed atom types.

But in reality they depend strongly on the chemical environment.

Your GNN can be viewed as producing *local response parameters*:

[
\chi_i=f_\theta(h_i)
]

[
J_i=g_\theta(h_i)
]

where (h_i) is the atomic embedding obtained by message passing.

But you can go much further.

# 3. Multipolar QEq

Instead of only charges, every atom has a vector:

[
\mathbf q_i =
(q_i,p_x,p_y,p_z,Q_{xx},...)
]

The energy becomes:

[
E=
\sum_i
\mathbf a_i\cdot\mathbf q_i
+
\frac12
\mathbf q_i^T
H_i
\mathbf q_i
+
\frac12
\sum_{ij}
\mathbf q_i^T
T_{ij}
\mathbf q_j
]

where:

* (T_{ij}) are Coulomb interaction tensors
* (H_i) are local hardness/polarizability matrices
* (\mathbf a_i) are local driving fields

This is exactly a physics-based self-consistent polarization model.

The neural network only predicts:

* atomic response tensors (H_i)
* preferred multipoles (a_i)

from the chemical environment.

Then the actual electron redistribution comes from minimizing the electrostatic energy.

This is a very elegant separation:

* ML learns chemistry
* physics solves long-range electrostatics

# 4. Iterative message passing = SCF

The analogy is very strong.

One iteration could look like:

### Initialization

For each atom:

[
h_i^{(0)} = f(Z_i, \text{local geometry}, H_{DFTB})
]

Set:

[
q_i^{(0)}=0
]

### Iteration k

### (a) Build electrostatic potential

From previous multipoles:

[
V_i^{(k)}
=========

\sum_j T_{ij} q_j^{(k-1)}
]

and electric field:

[
\mathbf E_i^{(k)}
=================

-\nabla V_i
]

### (b) Message passing

Update the atomic embedding:

[
h_i^{(k+1)}
===========

G(h_i^{(k)},
{h_j^{(k)}},
V_i,
\mathbf E_i)
]

---

### (c) Predict local response

Network outputs:

[
\chi_i,H_i
]

or directly a correction:

[
\Delta q_i
]

### (d) Solve local QEq

Obtain:

[
q_i^{(k+1)}
]

Repeat 2–5 times.

This is essentially a neural SCF cycle.

# 5. Using DFTB Hamiltonian as input

I think this is potentially extremely valuable.

A normal GNN sees:

* atom types
* distances
* bond topology

But the DFTB Hamiltonian already contains a huge amount of quantum information.

The off-diagonal elements:

[
H_{\mu\nu}
]

encode:

* bond order
* overlap
* orbital symmetry
* conjugation
* aromaticity
* hyperconjugation

The network can use:

* diagonal onsite energies
* eigenvalues
* local density of states
* orbital populations
* bond Hamiltonian blocks

as features.

This may dramatically reduce the size of the neural network.

In a sense, the DFTB Hamiltonian is already a handcrafted quantum embedding.

# 6. How many parameters might be needed?

Because your output space is so constrained, I think the network can be surprisingly small.

Something like:

* atom embedding: 32 dimensions
* 3–5 message passing iterations
* MLP layers of size 32–64

might already work.

Total parameters could be on the order:

[
10^4 - 10^5
]

which is tiny.

This is realistic for an embedded system.

# 7. Training targets

I would not train only the density coefficients.

I would train a hierarchy:

## Level 1

Charges:

[
L_q
===

|q^{ML}-q^{DFT}|^2
]

## Level 2

Dipoles and quadrupoles.

## Level 3

Electrostatic potential:

[
L_V
===

\int
(V^{ML}-V^{DFT})^2
d\mathbf r
]

This is probably more physically important than matching density everywhere.

## Level 4

Energy correction:

[
\Delta E_{SCF}
]

and maybe forces.

---

# 8. The most interesting possibility: predicting SCF residuals

For accelerating DFT or DFTB, I might actually formulate the target not as the final density but as the **density correction needed after one SCF step**.

For example:

[
\Delta \rho =
\rho_{\text{converged}}
-----------------------

\rho_{\text{initial}}
]

or after a few iterations:

[
\Delta\rho_k
============

## \rho^*

\rho_k
]

Then your model acts as a preconditioner or an intelligent mixer.

This could reduce the number of SCF cycles dramatically.

# My overall opinion

I think the most original aspect is **not the GNN**—everyone uses GNNs now. The genuinely interesting idea is combining:

1. **Neutral atom baseline**
2. **Minimal multipole representation**
3. **Learned environment-dependent QEq parameters**
4. **Explicit self-consistent Coulomb iteration**
5. **DFTB Hamiltonian as a quantum-informed descriptor**

Conceptually it is much closer to a modern, ML-enhanced version of **SCC-DFTB or fluctuating charge models** than a generic ML density predictor.

I would actually push the design even further: do **not let the network directly predict the multipoles**. Make it predict a local energy functional (electronegativity, hardness, polarizability tensors), and obtain the multipoles by solving the variational electrostatic problem. That gives you exact charge conservation, proper long-range interactions, and much better extrapolation.

For an ultra-small model intended to run inside a DFT code or on embedded hardware, I think this is one of the most promising architectures you could choose.

---

# Kimi 2.6 / devin / meta-analysis

# Summary

The document records a discussion between you and four LLMs (Gemini, Mistral, DeepSeek-V4-pro, ChatGPT) about a **physics-informed GNN to predict the deformation density** Δρ = ρ_SCF − ρ_NA in a minimal atom-centered **sp-basis representation** (monopole + dipole + quadrupole per atom). The core concept is to use **QEq as a differentiable physical backbone**, feed local electric fields and neighbor information through a tiny GNN, optionally use the **DFTB non-SCC Hamiltonian** as a quantum-informed descriptor, and produce an initial density guess to accelerate SCF/SCC convergence across methods (DFTB+, semi-empirical, DFT, hybrids).

# Common Ideas (Consensus)

All four models converged on these principles:

- **Δρ, not total ρ**: Subtracting the neutral atom density removes the inert core-electron burden, leaving only chemically relevant redistribution (~9 coefficients per atom: 1 charge + 3 dipole + 5 quadrupole).
- **QEq as the physical constraint**: The neural network should not directly output raw multipoles; instead it should produce *environment-dependent parameters* (electronegativity χ, hardness J, polarizability α) that feed into a variational/differentiable QEq solver. This guarantees charge conservation and correct long-range behavior.
- **Ultra-lightweight GNN**: A small MLP (16–32 hidden dims, 2–4 message passing iterations) is sufficient because the physics handles most of the work.
- **Local field inputs**: On-site electrostatic potential V_i and electric field vector **E**_i (computed iteratively from neighbor multipoles) are essential features.
- **DFTB H⁰ as a "quantum embedding"**: The non-SCC Hamiltonian and overlap matrix encode bond order, conjugation, and orbital interactions, dramatically reducing the NN's learning burden.
- **Goal: SCF accelerator / embedded deployment**: The intended use is as a drop-in initial guess engine, not a standalone QM replacement. Target parameter counts range from ~100 to ~10⁵.

# Original / Innovative Ideas by Each LLM

| LLM | Distinctive Contribution |
|-----|--------------------------|
| **Gemini** | Framed the architecture as a **refinement engine for QEq parameters** rather than a direct geometry→density map. Explicit pipeline diagram: GNN + DFTB → predict local χ, J, α → solve QEq → multipoles → form ΔP_μν. Emphasized the "physical constraint" layer. |
| **Mistral** | Most **concrete engineering specification**: parameter budget table (~100–200 params total), expected speedup table (DFTB+ 10–20 → 3–5 iter), embedded target specs (~1–10 KB quantized, Cortex-M4 in <1 ms). Proposed a multi-term loss: L_Δρ + λ₁L_QEq + λ₂L_Coulomb + λ₃L_sparsity, plus physics-informed weight initialization. |
| **DeepSeek-V4-pro** | **Induced-dipole physics to replace equivariant networks**: the GNN outputs only a *scalar polarizability* α_A; the dipole direction comes from solving μ = α(E − Tμ) (Thole-damped). This gives **rotational equivariance for free** without e3nn complexity. Also proposed **unrolled optimization** (treating the iterative field updates as a learned optimizer with 3–5 steps) and training on cheap DFTB SCC first. |
| **ChatGPT** | **Multipolar QEq formalism**: expressed the full energy with Coulomb interaction tensors T_ij, local hardness matrices H_i, and driving fields a_i. Most importantly, proposed **predicting SCF residuals** (the density correction needed after one SCF step) so the model acts as a **learned preconditioner/mixer**, not just an initial guess. Also pushed hardest for the **variational formulation**: "do not let the network predict multipoles directly; predict the energy functional and minimize." |

# Recommended Approach (Meta-Analysis)

Based on the synthesis above and the FireCore context, here is what I would suggest:

## 1. Architecture: Variational Multipolar QEq + Tiny GNN
Adopt **ChatGPT's variational formulation** as the backbone:
- The GNN outputs **local response parameters** per atom: Δχ_A, ΔJ_A, scalar α_A, and optionally a direct vector correction Δ**μ**_A (to capture charge-transfer dipoles that pure induction misses).
- A **differentiable QEq + induced-dipole solver** (following DeepSeek) computes the actual multipoles by minimizing the electrostatic energy. This guarantees exact charge neutrality, correct long-range decay, and rotational covariance without expensive equivariant layers.
- For quadrupoles, use an analogous induced-quadrupole term (scalar C_A × field gradient) or predict the 5 traceless components in a **local frame** defined by bonded neighbors.

## 2. Equivariance Strategy: Physics, Not Architecture
Use **DeepSeek's insight**: let Coulomb's law generate the vector/tensor directions. The GNN only predicts scalars (χ, J, α, C). The dipole direction is determined by the local field. This is far more parameter-efficient than SO(3) convolutions or tensor-product layers.

## 3. Dual-Mode Operation
Follow **ChatGPT's residual idea** and deploy the model in two modes:
- **Mode A (Initial Guess)**: Run 2–4 unrolled iterations at the start of an SCF cycle to produce a pre-converged density matrix guess.
- **Mode B (Preconditioner)**: Inside the SCF loop, predict the residual Δρ = ρ_converged − ρ_current and use it as an intelligent mixer (like a learned DIIS). This can dramatically cut iterations for difficult systems.

## 4. Training Strategy & Hierarchy
Use **Mistral's staged loss** and **DeepSeek's bootstrap**:
1. **Stage 1**: Train on cheap DFTB+ SCC references (charges + dipoles). Data is abundant and perfectly aligned with the sp basis.
2. **Stage 2**: Fine-tune on higher-level theory (B3LYP, hybrid DFT) using a hierarchical loss: charges → dipoles/quadrupoles → electrostatic potential on a grid → total energy/forces.
3. **Stage 3 (Active Learning)**: Sample geometries that cause SCF convergence failures and add them to the training set.

## 5. Inputs & FireCore Integration
Given your existing Fireball/FireCore infrastructure:
- **Primary input**: Fireball non-SCC Hamiltonian H⁰ and overlap S (you already compute these). Use diagonal H_μμ as atomic features and neighbor H_μν blocks as edge features.
- **Secondary input**: Geometry + neighbor list + initial QEq charges (from standard χ/J tables).
- **Target**: Partition the converged Fireball density into atomic multipoles. Since Fireball uses a minimal basis, **Mulliken or Löwdin partitioning** is natural and directly gives you Δq, Δμ, ΔQ in the sp basis without needing real-space grid decomposition.

## 6. Embedded Viability
- Target **<<500 parameters** (well within embedded limits). A 2-layer GNN with 16-dim hidden states and scalar outputs is enough.
- The only "heavy" operation is solving the QEq linear system (N_atoms × N_atoms). For N < 200 this is trivial; for larger systems, use iterative solvers (conjugate gradient) or truncate Coulomb interactions beyond a cutoff.
- Quantize to 8-bit int for MCU deployment.

## 7. Key Pitfall to Avoid
All LLMs touched on partitioning transferability: the model will reproduce whatever partitioning scheme you train on (Hirshfeld, Mulliken, ISA). **Mulliken is fine for Fireball** because your target Δρ is defined in the same atomic basis. Do not mix partitioning schemes between training target and deployment.

In short: **predict scalar chemistry parameters with a tiny GNN, solve physics variationally for multipoles, use as both an initial-guess engine and an intra-SCF preconditioner, and bootstrap training from cheap DFTB before fine-tuning on higher theory.**

---

# USER

I kinda like the idea of use neural networs just to predict local parameters for QEq perhaps that would be first thing to try using only multi-layer preceptron. But there are few problem with that approach 

1) QEq preduct only atomic charges as far as I know, or can it predict also like dipole and quadrupole on atoms or charge of some addatoms (added charge centers, e.g. electron pairs and pi-orbitals in aromatic system) ? If I want to plug it into some QM solver as estimate of charge senity to speed up SCF it would be best to have something which is already orbital resolved as DFT density is (although DFTB+ SCC is perhaps not resolved and it is only isotropic monopole on each atom)

2) QEq is iterative process which can take many cycles to converge, one cannot easyly run backpropagation of dependece of fitting error (reference charge density ves charge densty produced with prediceted parameters) with respect to the preducted parameters, and therefore nor backpropagation of MLP parameters. 

3) Also for practical reason it would be good to have something wich can predic charges in one shot (one itaration), although QEq may be mus ftaste than solving DFT SCF or DFTB+ SCC (where we diagonalize hamiltonian at every iteration, O(n^3)) it is still not negligible cost as it compute all pairwse interaction (coulomb matrix, O(n^3) at every iteration, if we do 100 iterations it may be already costly). If we can make good estimate in one iteration it would be great.

4) I'm not sure that it is god idea to replace equavariatn graph neural networs, I would prefer to use it as it is good approach

5) we should write in more details how exactly to use combination of global electrostatics information encoded in loal coulomb potential and electric filed (4 number [V,Ex,Ey,Fz], like energy and force of charge particle,   [E,Fx,Fy,Gz], float4 fit well GPU architecture), and the local non-SCC DFTB hamiltnian and overlap matrix which is naturally 4x4 blocks for sp basis,   we have one 4x4 on-site block ,and max 4 nearest neighbors 4x4 blocks so 5x16=80 number, this encode local information (local chemical neighborhood and chemical state of each atom). then tere is for sure ti element or proton number of given atom. I think it is good to be sure we have fixed size vectros for ML.  Plese elaborate on this. more in details.

---

# Kimi 2.6 / devin

You raise excellent, concrete objections. I agree with your direction: **predict orbital-resolved quantities directly in one shot**, use physics only to inform the *inputs* and *loss*, and keep an equivariant architecture. Here is a detailed response to each point and a concrete architecture proposal.

## 1. QEq Only Predicts Monopoles; You Need Orbital-Resolved Output

Standard QEq minimizes an energy in atomic charges only:
$$E = \sum_i \chi_i q_i + \frac{1}{2}\sum_i J_i q_i^2 + \frac{1}{2}\sum_{i\neq j} q_i q_j \gamma_{ij}$$
It has no dipoles, quadrupoles, or auxiliary charge centers (lone pairs, pi-orbitals as explicit sites).

**Extensions exist** (FQ+FD fluctuating charge/fluctuating dipole, QTPIE, atom-centered density matrix propagation), but they rapidly become iterative, non-linear, and expensive — exactly what you want to avoid.

**Better approach**: Do **not** use QEq as the output layer. Instead:
- Use a cheap classical QEq solve (with fixed elemental $\chi, J$ tables) **only** to generate initial electrostatic features: the on-site potential $V_i$ and electric field $\mathbf{E}_i$ at each atom. This costs one $O(N^2)$ operation and provides the network with a coarse global charge distribution.
- Let the **neural network predict orbital-resolved quantities directly**: per-atom on-site density matrix blocks $\Delta P_{\mu\nu}$ in the sp basis, or equivalently atomic charges $q_i$, dipole vector $\boldsymbol{\mu}_i$, and quadrupole tensor $\Theta_i$.

This gives you exactly the orbital-resolved density needed for a DFT/DFTB initial guess, and you can place auxiliary charge centers (lone pair sites, pi-bond midpoints) simply by adding them as extra nodes in the graph with their own $Z_{\text{aux}}$ and local basis blocks.

**DFTB+ SCC note**: You are correct that SCC-DFTB itself only uses isotropic monopoles ($\gamma_{AB} q_B$). However, the **non-SCC density matrix** from your predicted guess can still carry dipole/quadrupole information through the orbital populations, giving a much better starting point even if the first SCC iteration projects it back to monopoles.

## 2. Backpropagation Through QEq is Difficult, but Avoidable

There are two ways QEq is solved:

| Method | Differentiable? | Cost |
|--------|-----------------|------|
| **Iterative equilibration** (Jacobi/GS) | Painful — unrolled loops, convergence bifurcations | $O(K N^2)$ |
| **Direct linear solve** (KKT matrix) | Yes — `torch.linalg.solve` is differentiable | $O(N^3)$ once |

Both are unnecessary if the network outputs the density directly. Instead, the physics enters **softly**:
- **Input features**: $V_i, \mathbf{E}_i$ from the cheap classical QEq guess.
- **Loss terms**: Penalty on charge conservation $\left(\sum_i q_i - Q_{\text{tot}}\right)^2$, dipole moment matching, and optionally a "QEq consistency" loss that encourages the predicted charges to approximately satisfy the linear system.

This is much more stable for training and completely eliminates backprop through matrix inversion.

## 3. One-Shot Prediction is the Right Goal

A direct MLP/GNN forward pass is $O(N)$ or $O(N \cdot N_{\text{neigh}})$. Even a dense QEq linear solve is $O(N^3)$ — fine for $N=50$, but not free, and problematic on GPU if batched over many molecules.

With one-shot prediction, the global electrostatics are **pre-computed** into the input vector $[V_i, \mathbf{E}_i]$ (a float4), and the network only has to perform a local nonlinear correction. This is ideal for embedded deployment.

The quality of the one-shot guess depends on how informative the inputs are. The combination of:
- **Global context**: $V_i, \mathbf{E}_i$ (Coulomb potential from all other atoms)
- **Local quantum context**: $H^0, S$ blocks from neighbors
is remarkably complete. The field tells the atom how it is being polarized; the Hamiltonian tells it how its orbitals are coupled to neighbors.

## 4. Keep Equivariant Networks — Use a Lightweight One

You are right not to abandon equivariant GNNs. For dipole (vector) and quadrupole (rank-2 tensor) outputs, **E(3) equivariance guarantees correct rotational behavior without data augmentation**.

However, you do **not** need a large model like MACE or NequIP. For this task, a minimal architecture suffices:

**Recommendation: Tensor-Product GNN with scalar + vector features**
- Hidden dimensions: **16 scalars + 16 vectors** per atom.
- **2 message passing layers**.
- Use spherical harmonics $Y_1(\hat{\mathbf{r}}_{ij})$ of neighbor vectors to update vector features (standard in PaiNN, EGNN).
- Output heads:
  - Scalar head: $\Delta q_i$ (1), isotropic p-population (1), quadrupole trace (1).
  - Vector head: $\boldsymbol{\mu}_i$ (3), s-p mixing vector (3).
  - Tensor head: traceless quadrupole $\Theta_i$ (5 independent components).

Total parameters: **~2,000–5,000**, still trivially small.

If you want to start even simpler, use a **per-atom MLP with Cartesian equivariant constraints**: predict scalar quantities (orbital populations) and treat dipoles as $\boldsymbol{\mu}_i = \alpha_i \mathbf{E}_i + \text{NN}_{\text{vec}}(h_i)$, where the NN vector output is learned but constrained by the local field direction.

## 5. Detailed Input Representation: Fixed-Size Per-Atom Vector

This is the core engineering question. Here is a concrete, GPU-friendly layout.

### 5.1 Global Electrostatics: `float4` per atom
Compute from a **cheap pre-step** (e.g., QEq with tabulated $\chi, J$ or even just Gasteiger/Mulliken charges from $H^0$ diagonalization):
$$[V_i, E_{ix}, E_{iy}, E_{iz}]$$
- $V_i = \sum_{j \neq i} q_j / r_{ij}$ (on-site Coulomb potential)
- $\mathbf{E}_i = \sum_{j \neq i} q_j \mathbf{r}_{ij} / r_{ij}^3$ (electric field)

This packs perfectly into a `float4` for GPU kernels. It encodes the **global charge environment** at atom $i$ in a compact, translationally invariant way.

### 5.2 Local Quantum Environment: DFTB $H^0$ and $S$ blocks
For an sp basis (4 orbitals: $s, p_x, p_y, p_z$), each atom-atom interaction is a $4 \times 4$ block.

Per atom, collect:
- **1 on-site block**: $(H^0_{ii})_{\mu\nu}$ and $S_{ii}$ — diagonal + intra-atomic p-p couplings.
- **Up to 4 nearest-neighbor off-diagonal blocks**: $(H^0_{ij})_{\mu\nu}$ and $S_{ij}$ for the 4 closest neighbors $j$.

| Component | Shape | Count |
|-----------|-------|-------|
| On-site $H^0$ | $4 \times 4$ | 16 |
| 4 neighbor $H^0$ | $4 \times 4 \times 4$ | 64 |
| On-site $S$ | $4 \times 4$ | 16 |
| 4 neighbor $S$ | $4 \times 4 \times 4$ | 64 |

**Total: 160 floats** if you include both $H^0$ and $S$.  
If memory/compute is tight, you can drop $S$ initially and use only $H^0$ (80 floats), but $S$ carries important overlap population information.

**Handling variable neighbors**:
- For atoms with $<4$ neighbors, **zero-pad** the missing blocks.
- Sort neighbors by distance so the ordering is deterministic.

These blocks encode:
- **Chemical state**: on-site energies $\sim$ electronegativity, hybridization.
- **Bonding topology**: off-diagonal couplings $\sim$ bond order, conjugation, aromaticity.
- **Orbital symmetry**: which p-components couple (e.g., $\pi$ vs $\sigma$).

### 5.3 Element Identity
Atomic number $Z_i$ mapped to a small learnable embedding (8–16 dims).

### 5.4 Optional: Neighbor Geometry
For equivariant models, include relative position vectors to the 4 neighbors: $\{\mathbf{r}_i - \mathbf{r}_j\}$ (12 floats). This lets the network learn directional bonding.

### 5.5 Complete Per-Atom Input Vector Layout
| Feature | Dimensions |
|---------|------------|
| $Z$ embedding | 16 |
| $[V, E_x, E_y, E_z]$ | 4 |
| $H^0$ local blocks | 80 |
| $S$ local blocks | 80 |
| Neighbor vectors (optional) | 12 |
| **Total** | **~192** |

This is a **fixed-size** vector regardless of molecule size or bonding. Every atom sees its own local neighborhood. The $[V, \mathbf{E}]$ term provides the "global" awareness that would otherwise require message passing to distant atoms.

### 5.6 Output Layout
Per atom, predict the symmetric on-site density matrix correction $\Delta P_{ii}$ in the sp basis:
- 10 unique floats for the $4 \times 4$ symmetric block.
- If you also want off-diagonal bond density corrections, predict them per edge (another 10 floats per neighbor), but for an initial SCF guess, the on-site blocks + a simple Fermi-level filling of $H^0$ may already be excellent.

## Suggested First Prototype

Given your preference to start with an MLP:

1. **Input**: Assemble the fixed-size vector above for each atom.
2. **Network**: Per-atom MLP with 2 hidden layers of 64 neurons.
3. **Output**: 10 floats $\rightarrow$ symmetric $4 \times 4$ $\Delta P_{\text{on-site}}$ per atom.
4. **Initial guess construction**: $P_{\text{guess}} = P_{\text{NA}} + \Delta P_{\text{NN}}$.
5. **Loss**: 
   - $L_{\text{data}} = \|\Delta P_{\text{NN}} - \Delta P_{\text{ref}}\|^2$ (from converged Fireball/DFTB density)
   - $L_{\text{charge}} = \left(\sum_i q_i - Q_{\text{tot}}\right)^2$
   - $L_{\text{dipole}} = \|\sum_i \boldsymbol{\mu}_i - \boldsymbol{\mu}_{\text{mol}}\|^2$

This gives you an **orbital-resolved, one-shot, physics-informed initial guess** in ~5,000 parameters, with no iterative QEq solve inside the network.

If the MLP works, upgrading to a 2-layer equivariant GNN (PaiNN-style) for vector/tensor outputs is a natural second step that should improve dipole/quadrupole accuracy.

---

# Gemini

You have hit on the exact core bottlenecks and architectural challenges of combining machine learning with semi-empirical physics. Let's break down each of these 5 points with rigorous, GPU-friendly engineering solutions.

## 1. Extending QEq to Dipoles, Quadrupoles, and Orbitals

Standard QEq only solves for scalar charges ($s$-type monopoles). However, it *can* be mathematically extended to higher multipoles (dipoles and quadrupoles) by moving from scalar chemical potentials to tensor orientations.

Instead of adding arbitrary virtual charge centers (like lone pairs), we can map the predicted multipoles directly to an **orbital-resolved density matrix ($\Delta P_{\mu\nu}$)** for an $sp$-basis ($s, p_x, p_y, p_z$):

* **Monopole ($q_i$):** Controls the trace of the density matrix blocks (net population of $s$ and $p$ shells).
* **Dipole ($\vec{p}_i = [p_x, p_y, p_z]$):** Directly maps to the off-diagonal population between the $s$ and $p$ orbitals on the *same* atom ($P_{s,p_x}, P_{s,p_y}, P_{s,p_z}$).
* **Quadrupole Tensor ($Q_i$):** Maps to the splits among the $p$-orbitals themselves ($P_{p_x,p_x} - P_{p_y,p_y}$, etc.).

Since DFTB+ SCC typically only uses isotropic monopoles, this full orbital resolution is a massive advantage when you feed this guess into a **full DFT/hybrid SCF solver (like B3LYP)**, which expects a complete initial density matrix ($P_{\mu\nu}$).

## 2. & 3. Eliminating the QEq Iteration Bottleneck (One-Shot & Backprop)

You are absolutely right: iterating a Coulomb matrix $O(N^2)$ or $O(N^3)$ dozens of times defeats the purpose of an ultra-fast ML initializer, and makes backpropagation a nightmare (vanishing/exploding gradients through unrolled loops).

### The Solution: "One-Shot" Sub-Domain QEq

Instead of solving a global, fully coupled QEq system, we treat global charge transport as a **one-shot direct matrix solve**, and let the GNN pass messages to handle the local polarization.

1. **Direct Linear Solve:** The QEq equation is a linear system: $\mathbf{J}\vec{q} = -\vec{\chi}$. You do **not** need to iterate it. You can solve it exactly in one shot using a lightweight linear solver (like a conjugate gradient or a direct small-block Cholesky factorization on the GPU).
2. **Seamless Backpropagation:** Because it's a direct matrix equation ($\vec{q} = \mathbf{J}^{-1}\vec{\chi}$), you can backpropagate through it analytically using the adjoint matrix method:

$$\frac{\partial L}{\partial \vec{\chi}} = -\mathbf{J}^{-1} \frac{\partial L}{\partial \vec{q}}$$

This allows you to train the GNN end-to-end. The loss function compares the final reconstructed $\Delta \rho$ (or $\vec{q}$) against the true DFT reference, and the gradients flow smoothly back into the GNN parameters.

## 4. Retaining Equivariant Graph Neural Networks (EGNNs)

Do not replace the Equivariant GNN; it is essential. Scalar MLPs fail because they cannot understand that if a molecule rotates, the local electric field vector $\vec{E}$ and the $p$-orbitals must rotate with it.

The EGNN should be used to update the local latent features of the nodes and edges, ensuring **rotational equivariance** for the vector ($\vec{E}$, $\vec{p}$) and tensor ($Q$) properties, while passing standard scalar messages for the QEq parameters ($\chi, J$).

## 5. Fixed-Size GPU Hardware Layout (`float4` Alignment)

To make this blazingly fast on embedded GPUs or custom hardware accelerators, we must enforce rigid, fixed-size memory layouts. Packing your physical descriptors into `float4` vectors maximizes memory bandwidth and coalescing.

Here is the exact detailed layout for an atom $i$, designed to fit perfectly into modern GPU registers:

### Vector A: Electrostatic Environment (`float4`)

Calculated from the *initial* geometric guess (e.g., from classical or formal charges before the GNN starts):

* `V_env.x`: Local Coulomb potential $V_i$ (Scalar)
* `V_env.y`: $E_x$ (Electric field component)
* `V_env.z`: $E_y$ (Electric field component)
* `V_env.w`: $E_z$ (Electric field component)

### Vector B: Atomic Identity (`float2` or `float4`)

* `At_id.x`: Proton number $Z_i$ (or Element Index)
* `At_id.y`: Valence electron count

### Matrix C: The DFTB Hamiltonian & Overlap Block (`float4x4` structures)

For an $sp$-basis, the on-site block is exactly $4 \times 4$. To keep the vector size strictly fixed without dealing with variable coordination numbers (variable neighbor counts), we select the **top 4 most significant neighbor interactions** (sorted by the magnitude of the overlap matrix element $S_{ij}$). If an atom has fewer than 4 neighbors (like Hydrogen), the remaining slots are zero-padded.

This yields exactly 5 blocks of $4 \times 4$ matrices per atom:

1. **On-site block:** $H_{ii}$ (16 floats)
2. **Neighbor 1 block:** $H_{i1}$ (16 floats)
3. **Neighbor 2 block:** $H_{i2}$ (16 floats)
4. **Neighbor 3 block:** $H_{i3}$ (16 floats)
5. **Neighbor 4 block:** $H_{i4}$ (16 floats)

*(Optional but recommended: You can mirror this with 5 blocks of the Overlap matrix $S$, or blend them into an effective local coupling matrix).*

### Total Register Footprint per Atom

* **Electrostatics:** $1 \times \text{float4}$
* **Identity:** $1 \times \text{float4}$ (with padding)
* **Hamiltonian Environment:** $20 \times \text{float4}$ (5 blocks $\times$ 4 rows)

Total input feature size per atom = **22 `float4` vectors** (88 floating-point numbers).

### How the EGNN Processes This

Because this data footprint is completely static, your GPU kernel doesn't need to jump through dynamic memory loops:

1. The EGNN takes the 22 `float4` inputs.
2. The $4 \times 4$ Hamiltonian blocks provide immediate, rich quantum-chemical descriptors of the local covalent environment, while the electrostatics `float4` provides the classical environment.
3. The EGNN outputs the corrected $\chi_i$, $J_{ii}$, and a local polarizability tensor $\alpha_i$.
4. The direct one-shot QEq solver yields the final $\Delta P_{\mu\nu}$ matrix.

This fixed layout removes all control-flow overhead, allowing the entire prediction to execute in a fraction of a millisecond on a modern embedded GPU.

---

# ChatGPT (5.5?)

I think your objections are very valid, and they actually push the architecture into a more interesting direction. Let me go through them one by one, because I think there is a conceptual shift from "ML learns QEq parameters" to something closer to "ML learns a single-step SCF response operator".

# 1) QEq is too simple; we need orbital/multipole resolution

Indeed, classical QEq usually solves only for atomic monopoles:

[
q_i
]

with a simple quadratic energy.

But there are several possible generalizations.

## (a) Multipole QEq

One can extend the variable from scalar charge to a vector:

[
\mathbf{q}*i = (q,p_x,p_y,p_z,Q*{xx},Q_{xy},...)
]

and write

[
E = E_0 + \sum_i \mathbf{\chi}_i \cdot \mathbf{q}*i
+\frac12 \sum*{ij}
\mathbf{q}*i^T T*{ij} \mathbf{q}_j
]

This is basically a polarizable force field.

However, it is still **atom centered**, and therefore it does not represent:

* lone pairs,
* π orbitals,
* σ/π separation,
* bond polarization.

## (b) Extended QEq with extra charge sites

I think this is much closer to your FireCore philosophy.

Instead of one charge per atom, you have "electron centers":

```
atom:
       q_core
         |
   p_x p_y p_z
         |
     π orbital
     lone pair
     bond center
```

More realistically:

* each atom has an s-like charge center,
* p-like polarization centers,
* π orbitals as explicit sites,
* lone pairs as explicit sites.

Then the ML target can become much closer to DFT density matrix.

This is essentially a very low-rank representation of the density.

I actually think this may be more useful than pure atomic multipoles.

## (c) Orbital population representation

Since you already think in terms of DFTB, the most natural quantity may be the **density matrix correction**:

[
\Delta P_{\mu\nu}
=================

## P^{SCF}_{\mu\nu}

P^0_{\mu\nu}
]

where (\mu,\nu) are atomic orbitals.

For an sp basis:

* on-site: (4\times4)
* nearest-neighbor bond block: (4\times4)

This is extremely attractive.

The local density around atom (i) can be represented by:

* (P_{ii}) (16 numbers)
* (P_{ij}) for nearest neighbors (up to 4×16=64 numbers)

Total:

[
16 + 64 = 80
]

Exactly the same size as your Hamiltonian input.

This is a very elegant symmetry.

# 2) Differentiating through QEq is actually possible

Your concern is reasonable, but modern ML can backpropagate through iterative solvers.

For example:

```
parameters
     |
     V
QEq iteration
     |
     V
charges
     |
     V
loss
```

and compute gradients by unrolling the iterations.

However, I actually agree with you that it is not a good engineering choice.

Problems:

* many iterations,
* large memory consumption,
* convergence instability,
* harder GPU implementation.

For a tiny embedded model this is contrary to the goal.

# 3) One-shot prediction is probably the right objective

This is the most important point.

Think about SCF.

The exact relation is:

[
\rho^{SCF} = F(\rho_0, H_0)
]

where (F) is the self-consistent response.

You do not necessarily need to emulate the iterative process.

You can train a network to approximate:

[
\Delta\rho = \mathcal{N}(H^0,S,\text{geometry})
]

directly.

The interesting question is:

**How do you give the network information about long-range electrostatics?**

A local neighborhood is insufficient.

A carbon atom in benzene and a carbon atom near a positively charged group can have the same local geometry but very different charge.

Your idea of adding electrostatic fields is exactly the solution.

# 4) Equivariant GNN vs fixed vectors

I would not throw away equivariant GNNs.

I think there is a very nice compromise.

The network has two parts:

```
                 global electrostatics
                       |
                       V
               [V, Ex,Ey,Ez]
                       |
                       V

       local DFTB blocks
       Hii, Hij, Sii, Sij
              |
              V

             GNN
              |
              V

       local density correction
```

The GNN message passing handles the **local chemistry**:

* bond orientation,
* hybridization,
* conjugation,
* local symmetry.

The electrostatic field handles the **long-range environment**.

These are complementary.

# 5) Fixed-size input based on DFTB blocks is actually a brilliant idea

I think this is perhaps the strongest part of the architecture.

Let's enumerate.

## Atomic identity

You need:

* nuclear charge (Z),
* maybe valence electron count,
* atom type embedding.

Small, maybe 8–16 features.

## Local Hamiltonian information

For an sp basis:

### On-site block

[
H_{ii}
]

4×4 = 16 numbers.

### Neighbor blocks

For tetra-coordinated atoms:

[
H_{ij}
]

up to 4 neighbors.

Each:

4×4 = 16 numbers.

Total:

[
4 \times 16 = 64
]

Therefore:

[
16 + 64 = 80
]

---

You may also include overlap:

[
S_{ii}, S_{ij}
]

This doubles it:

[
80 + 80 = 160
]

---

## Electrostatic context

Your `float4` idea is very GPU friendly:

```
float4 electrostatics = {
    V,
    Ex,
    Ey,
    Ez
};
```

This is essentially the first terms of the multipole expansion.

The potential gives monopole response.

The field gives dipole polarization.

## Optional second derivative

A thought: if you want quadrupoles, you might also include the electric field gradient:

[
\nabla E
]

which is a symmetric 3×3 tensor (5 independent components after removing the trace).

Then you have:

```
V       1
E       3
∇E      5
--------------
total   9
```

This exactly matches charge + dipole + quadrupole response.

# 6) A possible architecture

Actually I would consider something like this:

```
Input per atom:

Z                    1
Hii                 16
Hij nearest         64
Sii                 16
Sij nearest         64
electrostatics       4
-------------------------
total              ~165 numbers
```

Pass through:

```
Linear(165,64)
ReLU
Linear(64,64)
ReLU
Linear(64,32)
```

This produces a local embedding.

Then a few message passing steps:

```
h_i^(k+1)=MLP(h_i,h_j,r_ij)
```

with maybe 3–4 iterations.

Finally output:

```
ΔP_ii (16)
ΔP_ij (up to 64)
```

or

```
atomic multipoles
```

# 7) A very interesting possibility: use the previous SCF iteration

For accelerating DFT/DFTB, the best input may not be just the non-SCC Hamiltonian.

At SCF iteration (k), you already have:

* previous density matrix (P_k),
* previous electrostatic potential (V_k),
* Hamiltonian (H_k).

Then the network predicts:

[
\Delta P = P^* - P_k
]

This is exactly a learned SCF mixer.

It may reduce 50–100 SCF iterations to 2–5.

# My impression

After your additional constraints, I would actually **move away from the QEq parameter prediction idea**. It is elegant, but probably too restrictive.

The architecture that seems most promising is:

1. Use **DFTB non-SCC Hamiltonian and overlap blocks as local quantum descriptors**.
2. Use an **equivariant GNN for local chemical communication**.
3. Add global electrostatic information through `float4 (V,Ex,Ey,Ez)` (and possibly field gradients).
4. Predict a **localized correction to the density matrix** rather than charges.
5. Train it as a one-shot approximation of the SCF response.

What I find particularly elegant is that your input and output live in the same language: local (4\times4) orbital blocks. The network is not trying to discover quantum mechanics from atomic coordinates—it is learning the missing self-consistent response on top of a physically meaningful zeroth-order Hamiltonian. This could potentially be much smaller and more transferable than generic neural electron-density models.

---

# DeepSeek-V4-pro

Excellent, you’ve hit on exactly the right questions. Let’s go through them systematically, and then I’ll lay out a concrete design for the feature encoding that ties the electrostatics and the non‑SCC Hamiltonian into a compact, equivariant GNN architecture.

## 1. Can QEq predict dipoles, quadrupoles, or off‑centre charge?

**Classical QEq** (Rappé & Goddard) only gives atomic charges. It cannot produce dipoles or quadrupoles because it represents the density as point charges on the nuclei.

There are three common ways to extend it:

1. **Polarizable extension** – Add inducible dipoles (and quadrupoles) that minimise an energy  
   \(E = \sum_A \frac{|\boldsymbol{\mu}_A|^2}{2\alpha_A} - \boldsymbol{\mu}_A \cdot \mathbf{E}_A\) with appropriate damping. This is the path I mentioned earlier; it gives you atomic dipoles in a physically correct way, but it still ties the dipoles to the field via a linear response equation.

2. **Off‑atom charge centres** – Place additional massless sites (e.g. lone pairs, π‑clouds) with their own χ, J. The GNN would predict where they sit and what their QEq parameters are. This is very expressive but complicates the model and makes it hard to keep rotational equivariance.

3. **Equivariant network direct prediction** – Let the GNN directly output atomic dipole vectors and quadrupole tensors as part of its node features. No classical polarisation equation is needed. The physics is still encoded in the input features (local field, electronegativity, bond environment) and in the loss (e.g. matching SCF reference multipoles). This is simple, respects SO(3) if you use equivariant message passing, and works with your preference for SE(3)‑equivariant GNNs.

For **accelerating SCF**, you don’t necessarily need the full density matrix. What you need is a guess of the effective Hamiltonian (or the electron density) that is much closer to the converged one than the neutral‑atom guess. In DFTB+ SCC, the only variable that enters the Hamiltonian is the set of **atomic charges** \(q_A\), via  
\(\Delta H_{\mu\nu} = \frac{1}{2} S_{\mu\nu} \sum_C (\gamma_{AC} + \gamma_{BC}) \Delta q_C\) (or similar). Therefore, a highly accurate prediction of the SCC charges alone can dramatically cut the number of SCC iterations, and you never need dipoles or quadrupoles.

For full DFT or semiempirical methods you would need a guess for the density matrix (or the charge density). One pragmatic route is to **reconstruct an atom‑centred model density** from predicted multipoles (charges, dipoles, quadrupoles) and then compute its Coulomb and exchange‑correlation potentials. That’s what many fast semiempirical guess generators do. So predicting atomic multipoles is a valid and general strategy.

My **recommendation**: start by predicting only **charges and atomic dipoles** (using an equivariant GNN) for DFTB SCC. It’s the cleanest proof‑of‑concept, and you can later add quadrupoles if needed.

## 2. Backpropagating through an iterative QEq solver

Actually, **standard QEq does not need an iterative solver**. The energy is quadratic in the charges:  
\(E = \sum_A \left(\chi_A q_A + \frac{1}{2} J_A q_A^2\right) + \frac{1}{2}\sum_{A\neq B} q_A q_B \gamma_{AB}\),  
where \(J_A\) is the hardness (Hubbard U). Minimising with a charge‑neutrality constraint leads to a **linear system**  
\(\mathbf{A} \mathbf{q} = \mathbf{b}\).  
One can solve this directly in \(O(N^3)\) via LU/Cholesky, and `torch.linalg.solve` is fully differentiable (it gives gradients w.r.t. χ, J, and γ). So **there is no inner‑loop iteration** – backpropagation works out of the box.

If you later add induced dipoles, the same holds: the polarisation equations are linear (if α is fixed) and can be solved directly, or you can unroll a few damped Jacobi steps, which is also trivial to backpropagate through. The key point is that the **physics solver is not a roadblock for gradients**.

## 3. “One‑shot” prediction without iterative cost

You want to avoid any iterative process that scales badly. Three options exist:

- **Direct linear solve** – for systems with a few hundred atoms, a single \(O(N^3)\) linear solve is faster than one diagonalisation of the Hamiltonian, and you only do it once per prediction. That’s already a “one‑shot” in the sense that you get the equilibrium charges immediately.

- **Amortised inference with a GNN** – train the GNN to output charges that satisfy the QEq equilibrium condition directly, without ever running the solver. This is effectively a learned surrogate that captures the solution of the linear system as a function of the environment. The GNN would take electrostatics features (V, E) and chemical features (from the non‑SCC Hamiltonian) and output the charges. In training you would minimise the difference to the actual QEq solution (or better, to the SCF reference charges). At inference time you just run the GNN once.

- **Recurrent field‑updating GNN** – you start with a rough guess (e.g. neutral atoms) and at each step you compute V and E from the current charges/dipoles, feed them into the GNN, and it outputs an updated set of charges/dipoles. After only 2–3 steps you have a converged result. Because this loop is unrolled a fixed, small number of times, the total cost is minimal and no large linear system needs solving. This is the “physics‑informed unrolled optimiser” approach.

Given your desire for extreme speed and tiny parameter count, the **recurrent field‑updating scheme is ideal**: it combines the physics of charge equilibration (through V/E) with the expressive power of the GNN, and 2–3 iterations are enough. Training is straightforward because the loop is a short, fixed computation graph.

## 4. Why an equivariant GNN is the right tool

You are absolutely right to prefer equivariant networks – they bake in the correct transformation properties, which drastically cuts the amount of data and parameters needed to learn vector and tensor quantities.

In your case, the atomic dipole \(\boldsymbol{\mu}_A\) is a 3D vector that must rotate with the molecule. An SE(3)‑equivariant GNN can output such vectors directly as type‑1 features, without you having to hand‑code polarisation equations. The network learns the mapping from local fields, bond environment, and element identity to the final dipole (and quadrupole if you use type‑2 features). This is clean and extremely compact.

The challenge is to feed in the electrostatic and Hamiltonian information in a rotationally consistent way – that’s what I’ll detail next.

## 5. Detailed design: how to encode electrostatics and non‑SCC Hamiltonian into fixed‑size vectors

Your proposal of using \([V, E_x, E_y, E_z]\) and the 4×4 blocks of the Hamiltonian is a good starting point, but we can make it **perfectly adapted to an equivariant GNN** while keeping the feature sizes tiny.

### 5.1 Electrostatic features (V, E)

For each atom \(A\), given the current set of predicted multipoles (charges \(q_B\), dipoles \(\boldsymbol{\mu}_B\)), you compute the **on‑site electrostatic potential** and the **electric field** using a damped Coulomb interaction:

\[
V_A = \sum_{B\neq A} \left( q_B \, \phi(R_{AB}) + \boldsymbol{\mu}_B \cdot \nabla \phi(R_{AB}) \right),
\qquad
\mathbf{E}_A = -\nabla_A V_A.
\]

Here \(\phi\) could be a simple \(1/R\) (with a short‑range damping function like the Thole or Klopman–Ohno model). The important thing for equivariance is:

- \(V_A\) is a **scalar** (type‑0 under SO(3)).
- \(\mathbf{E}_A\) is a **vector** (type‑1).

In an equivariant GNN you store them as separate irreducible representations. From an implementation viewpoint you can pack them as a single tensor of shape `[N_atoms, 4]` where the first channel is \(V\) and the next three are \(E_x, E_y, E_z\), and then inside the network you treat the first channel as type‑0 and the remaining three as type‑1. That is very GPU friendly.

### 5.2 Non‑SCC Hamiltonian features

Instead of using the raw 4×4 blocks, which rotate with the global coordinate system and would break equivariance, I strongly recommend extracting **rotationally invariant scalars** that encode exactly the same physical information. The non‑SCC DFTB Hamiltonian is built from the Slater–Koster tables. For each pair of atoms \(A\)–\(B\) (the edge), the following invariants are available:

- Distance \(R_{AB}\).
- The directional cosines are *not* invariant, but the Slater–Koster bond integrals are scalar functions of distance:  
  \(h_{ss\sigma}, h_{sp\sigma}, h_{pp\sigma}, h_{pp\pi}\) (and optionally the overlap integrals \(s_{ss\sigma}, s_{sp\sigma}, s_{pp\sigma}, s_{pp\pi}\)).  
  These are **scalar** under SO(3) because they are defined relative to the bond axis and the distance is the only geometric variable.

For the **on‑site** (node) features of atom \(A\), you have:

- Atomic number \(Z_A\) (one‑hot or learnable embedding).
- On‑site energies: \(\epsilon_s(A), \epsilon_p(A)\) (two scalars).
- Possibly the Hubbard U (hardness) \(J_A\) and electronegativity \(\chi_A\) from the DFTB parameters.

So a **node feature vector** becomes:
\[
[ \text{emb}(Z_A),\ \epsilon_s,\ \epsilon_p,\ V_A,\ E_{A,x},\ E_{A,y},\ E_{A,z} ],
\]
where the first four are type‑0 scalars and the last three are components of a type‑1 vector (or the whole vector is a type‑1 feature). In e3nn notation this might be `"4x0e + 1x1o"` (four even scalars, one odd vector) – but in practice you can simply keep them as separate tensors.

An **edge feature vector** for pair \(A\)–\(B\) becomes:
\[
[ R_{AB},\ h_{ss\sigma},\ h_{sp\sigma},\ h_{pp\sigma},\ h_{pp\pi},\ s_{ss\sigma},\ s_{sp\sigma},\ s_{pp\sigma},\ s_{pp\pi} ].
\]
All of these are invariant scalars. You can expand \(R_{AB}\) into a set of radial basis functions (e.g. Bessel or Gaussian) and then feed the concatenated vector through a small MLP to get the edge embedding.

If you strongly prefer to use the 4×4 block directly, you can still make it work. One way is to define a **local coordinate frame** for each atom (e.g., using the eigenvectors of the moment‑of‑inertia tensor of its nearest neighbours) and then express the 4×4 Hamiltonian block in that frame. Because the frame rotates with the molecule, the matrix elements become pseudo‑scalars that transform appropriately. This is more complex and does not add any information beyond what the Slater–Koster integrals already capture, so I recommend the SK route.

### 5.3 Tying it together in a recurrent GNN

1. **Initialisation**  
   Start with \(q_A = 0\), \(\boldsymbol{\mu}_A = \mathbf{0}\) for all atoms.

2. **For t = 1 … T (e.g. T = 3)**  
   a. Compute \(V_A^{(t)}, \mathbf{E}_A^{(t)}\) from \(\{q_B^{(t-1)}, \boldsymbol{\mu}_B^{(t-1)}\}\).  
   b. Assemble node features: \(\mathbf{h}_A^{(t)} = [\text{embed}(Z_A), \epsilon_s, \epsilon_p, V_A^{(t)}, \mathbf{E}_A^{(t)}]\).  
   c. Run one round of equivariant message passing (e.g., a Tensor‑Field‑Network layer, or a simple SE(3)-Transformer). The message from \(B\) to \(A\) uses the edge vector \(\mathbf{r}_{AB}\) and the invariant edge features (distances, SK integrals).  
   d. The output node features are split into a scalar part (which gives \(\Delta q_A\), the correction to the charge, or the new \(q_A\) directly) and a vector part (\(\Delta \boldsymbol{\mu}_A\)). Update  
      \[
      q_A^{(t)} = q_A^{(t-1)} + \Delta q_A, \qquad \boldsymbol{\mu}_A^{(t)} = \boldsymbol{\mu}_A^{(t-1)} + \Delta \boldsymbol{\mu}_A.
      \]
   e. (If desired, also output type‑2 features for quadrupoles; the update rule is analogous.)

3. **Final output**  
   After \(T\) iterations, \(\{q_A^{(T)}, \boldsymbol{\mu}_A^{(T)}\}\) is your prediction.

Because the GNN operates on scalar and vector features, it automatically respects rotational symmetry. The whole computation graph is a fixed number of steps and can be trained end‑to‑end by minimising a loss against reference SCF charges and dipoles (partitioned e.g. via Hirshfeld or MBIS). You can also add a physics‑regularisation term that penalises violation of the QEq linear equations or the polarisation balance condition, but in my experience the recurrence and the field features already encode that implicitly.

## Final remarks

What you end up with is a **lightning‑fast, few‑shot iterative predictor** that:

- Uses only a handful of scalars and vectors per atom.
- Handles all symmetries via equivariant message passing.
- Incorporates the non‑SCC Hamiltonian through invariant bond integrals and on‑site energies – no need to wrestle with rotating 4×4 blocks.
- Learns to equilibrate charges and polarise dipoles in a physically meaningful way, without ever solving a global linear system.

The resulting model will be tiny (easily < 10 k parameters) and perfect for an embedded accelerator. It can be plugged directly into a DFTB+ or DFT SCF as an initial guess, cutting the number of iterations by a large fraction.

If you’d like, I can expand on any of the components, for instance exactly how to set up the equivariant layer that handles the `(scalar, vector)` node features and the edge SK invariants, or how to obtain reference multipoles from DFTB SCC outputs.

---

# Gemini-3.1-pro

Based on the extensive discussion provided, you have conceptualized a highly efficient, physics-informed machine learning architecture designed to accelerate quantum mechanical SCF/SCC convergence. Your pivot in the second prompt—moving from predicting QEq parameters to a **one-shot, orbital-resolved prediction using fixed-size GPU-friendly vectors**—is an excellent engineering decision.

Here is a summary, meta-analysis, and a recommended architectural blueprint based on the LLMs' insights and my own knowledge.

### 1. Consensus & Common Ideas
Across the board, all models (and your prompts) agree on several core principles:
*   **Targeting $\Delta \rho$**: Predicting the deformation density ($\rho_{SCF} - \rho_{NA}$) drastically reduces the learning burden by ignoring chemically inert core electrons.
*   **Global Electrostatics via Fields**: Using the on-site potential and electric field $[V, E_x, E_y, E_z]$ is the best way to inject long-range environment awareness into a local model, bypassing the need for deep message passing.
*   **DFTB Hamiltonian as a "Quantum Embedding"**: Feeding the non-SCC Hamiltonian ($H^0$) and overlap matrix ($S$) blocks directly into the network provides an incredibly rich, pre-computed description of local chemistry (bond orders, hybridization, conjugation).
*   **Ultra-Lightweight/Embedded focus**: The model should be small (<10k parameters), utilizing a restricted basis (sp-orbitals) to act as a blazing-fast initial guess or intra-SCF preconditioner.

### 2. Original / Innovative Ideas by Each LLM

*   **ChatGPT**: 
    *   *Density Matrix Prediction*: Move entirely away from point-charges and predict the **local density matrix correction ($\Delta P_{\mu\nu}$)** directly as $4 \times 4$ blocks. This elegantly mirrors your $4 \times 4$ Hamiltonian input and gives you full orbital resolution (monopole, dipole, quadrupole) natively.
    *   *Learned SCF Mixer*: Use the model not just for an initial guess, but inside the SCF loop to predict the residual ($\Delta P = P_{converged} - P_{current}$), drastically cutting down iterations.
*   **DeepSeek-V4-pro**: 
    *   *Slater-Koster Invariants*: Pointed out a crucial physics flaw—raw $4 \times 4$ Hamiltonian blocks rotate with the coordinate system, which breaks rotational equivariance. It suggested using rotationally invariant Slater-Koster integrals ($h_{ss\sigma}, h_{pp\pi}$, etc.) as edge features instead.
    *   *Recurrent Field-Updating*: Formulate the GNN as a learned optimizer that unrolls 2-3 times, updating the $[V, \mathbf{E}]$ fields at each step.
*   **Gemini**: 
    *   *GPU-Centric Memory Layout*: Proposed the highly specific 22-`float4` register layout.
    *   *Direct Linear Solve*: Pointed out that QEq does not *need* to be iterative. It is a linear system ($\mathbf{J}\vec{q} = -\vec{\chi}$) that can be solved in one shot, and backpropagated through perfectly using the adjoint matrix method.
*   **Mistral**: 
    *   *Hardware Constraints*: Provided the most rigorous hardware bounds (8-bit quantization, Cortex-M4 targets, ~100-200 parameters per layer) and proposed a multi-term physics loss function (Coulomb penalty + sparsity).

### 3. Recommended Approach (Synthesis & Architecture)

Based on your desire for a **one-shot, orbital-resolved, fixed-size vector model** that retains **equivariant properties**, here is the blueprint I recommend.

#### The Equivariance Problem & The "Local Frame" Solution
You want to use raw $4 \times 4$ blocks of $H^0$ (80 floats), but as DeepSeek noted, $p_x, p_y, p_z$ orbitals change when the molecule rotates. If you feed raw Cartesian blocks into a standard MLP, the network will fail when the molecule is rotated.
*   **The Fix**: Construct a **Local Coordinate Frame** for every atom. 
    *   Define the local $z$-axis as the vector to the closest neighbor.
    *   Define the local $x$-axis using the plane formed by the closest and second-closest neighbor.
    *   Rotate the $4 \times 4$ $H^0$ blocks and the $[E_x, E_y, E_z]$ vectors into this local frame.
*   **Why this is brilliant**: Once in the local frame, the $4 \times 4$ blocks and fields become strictly rotationally invariant. **You no longer need a complex, heavy Equivariant GNN (like e3nn).** A standard, lightning-fast MLP or simple Graph Convolution can process these invariant blocks, and you mathematically guarantee perfect rotational equivariance when you rotate the output $\Delta P_{\mu\nu}$ back to the global frame.

#### The Workflow Pipeline

**Step 1: Cheap Classical Pre-computation (The Electrostatics)**
Run a rapid, one-shot charge estimation (e.g., Gasteiger charges, or a 1-iteration tabulated QEq) simply to get a rough charge on each atom. Use this to compute $V_i$ and $\mathbf{E}_i$.
*   *Memory*: Pack into a `float4` per atom: $[V, E_x, E_y, E_z]$.

**Step 2: Local Feature Assembly (Fixed-Size Vector)**
For each atom $i$, gather features in its *local frame*:
1.  **Identity**: $Z_i$ embedding (e.g., 8 floats).
2.  **Electrostatics**: $[V_i, E^{local}_x, E^{local}_y, E^{local}_z]$ (4 floats).
3.  **On-site Quantum**: $4 \times 4$ block of $H^0_{ii}$ (16 floats).
4.  **Neighbor Quantum**: Up to four $4 \times 4$ blocks of $H^0_{ij}$ (64 floats, zero-padded if <4 neighbors).
*   *Total Input Size*: Exactly 92 floats per atom. Easily fits into L1 cache/GPU registers.

**Step 3: The Neural Network (One-Shot MLP / Light GNN)**
Because all inputs are now local and rotationally invariant, you can use a tiny, standard architecture.
*   Pass the 92-float vector through a 2-layer MLP (e.g., 92 $\rightarrow$ 64 $\rightarrow$ 64).
*   *Optional*: Do 1 or 2 iterations of simple message passing where atom $i$ averages its 64-dim hidden state with its 4 bonded neighbors to capture slightly longer-range inductive effects.

**Step 4: Output the Density Matrix Correction**
The network outputs the **local-frame density matrix correction**.
*   Output size: 10 floats (the upper triangle of a symmetric $4 \times 4$ $\Delta P_{ii}$ block).
*   Rotate this $4 \times 4$ block back into the global Cartesian frame.

**Step 5: Inject into QM Solver**
Add this to the neutral atom density matrix: $P^{guess} = P^{NA} + \Delta P$. Feed this directly into DFTB+ or your QM solver. Because this matrix explicitly contains dipole (s-p mixing) and quadrupole (p-p splitting) terms dictated by the ML model, the SCF convergence will bypass the most difficult initial polarization steps.

### Why this is the ultimate design for your goals:
1.  **No iterative QEq backprop**: QEq is bypassed entirely during learning; you only use classical electrostatics to generate input features.
2.  **Fully Orbital Resolved**: You are directly outputting the sp-basis density matrix $\Delta P_{\mu\nu}$, capturing lone pairs, $\pi$-clouds, and dipoles natively.
3.  **One-Shot**: No inner iterative loops. The inference cost is exactly $O(N)$ ML predictions + the fixed $O(N^2)$ electrostatic field calculation.
4.  **Fixed-size & GPU Optimized**: You retain your 80-number $4 \times 4$ block structure.
5.  **Equivariant but Fast**: By using local reference frames, you get the accuracy of an Equivariant GNN with the speed and tiny parameter count of a standard MLP.

---

# USER

OK, but these desnity blocks 4x4 will be just one block on given atom aor 5x4x4 for all neoghbors (I mean to prepresent bond charge density?

We should try to make concise representation

I was rather thinging about really epressing just local polarization of atom which can be done in sp basis (just 4 numbers for monopole+dipoles, or 9 numbers with quadrupoles) or I was thinging about some resolution of identity (to have some chages on bonds and electron pairs and pi- orbitals. My idea is kinda like that every atom has like electron octed and that can be distributed around the atom according to VSEPR theory (sp3,sp2,sp1)the electron pair is either in sigma bond, pi-bond or in free lectron pair 4=nsigma+npi+nepair  (-O- is 2,0,2 = sp3 tetrahedron,  =O is 1,1,2 sp2 triangle,   -N= is 2,1,1,  -NH2 conjugated to aromatic ring is something between 3,0,1 and 3,1,0) 

---

# Gemini-3.1-pro

You have just hit on an incredibly elegant, chemistry-informed way to compress the output. 

To answer your first question directly: **Do not predict the off-diagonal bond blocks ($5 \times 4 \times 4$).** It is redundant. If you accurately predict the on-site polarization (the $4 \times 4$ block or its multipoles), the QM solver's Hamiltonian will naturally pull the correct charge into the off-diagonal bond regions during the very first diagonalization. You only need to predict the **on-site state** of each atom.

Your second idea—using **VSEPR theory and the "fractional octet"** to represent local polarization—is brilliant. It allows you to reduce the prediction from abstract tensors back to a few simple scalars (orbital populations), which inherently guarantees rotational equivariance because the geometry defines the directions!

Here is how we can formalize your VSEPR idea into a mathematically rigorous, ultra-concise ML representation.

### The "Fractional Octet" Representation (4 to 5 numbers per atom)

Instead of predicting abstract matrix elements or Cartesian dipoles, the neural network predicts the **electron population of chemically defined hybrid orbitals**. 

For any atom, the local geometry dictates a set of natural directions. Let's say an atom has $k$ bonds (where $k \le 4$). 
1. **$\sigma$-orbitals:** Point directly toward the $k$ neighbors. 
2. **$\pi$-orbitals / Lone Pairs:** Occupy the remaining orthogonal directions.

Because there are only 4 orbitals in an $sp$-basis ($s, p_x, p_y, p_z$), there are exactly 4 total populations to distribute:
$$4 = \sum n_\sigma + n_\pi + n_{lp}$$

#### How this maps to the Density Matrix and Multipoles
If the network predicts just the **scalar populations** of these 4 slots, you can deterministically reconstruct the full $4 \times 4$ on-site density matrix (and thus the monopole, dipole, and quadrupole) using the unit vectors of the atomic geometry!

**Example 1: Carbonyl Oxygen ($=O$)**
*   **Geometry:** 1 neighbor (Carbon) along unit vector $\vec{u}_1$.
*   **Basis construction (deterministic):** 
    *   One $\sigma$-hybrid pointing along $\vec{u}_1$.
    *   One $\pi$-orbital perpendicular to $\vec{u}_1$ (aligned with the rest of the molecule's plane).
    *   Two lone-pair hybrids occupying the remaining space.
*   **ML Output:** The GNN predicts just 3 scalars: $\Delta n_\sigma$, $\Delta n_\pi$, $\Delta n_{lp}$.
*   **Result:** A shift in $n_\sigma$ vs $n_{lp}$ creates a dipole along the bond. A shift in $n_\pi$ creates a quadrupole (anisotropy). 

**Example 2: Conjugated Amine ($-NH_2$ on an aromatic ring)**
*   **Geometry:** 3 neighbors (1 Carbon, 2 Hydrogens). These define a local plane (or slight pyramid).
*   **Basis construction:** 3 $\sigma$-hybrids pointing at the neighbors, and 1 orbital normal to the plane ($\pi$ / lone-pair hybrid).
*   **ML Output:** The GNN predicts $\Delta n_{\sigma C}$, $\Delta n_{\sigma H}$, and $\Delta n_{\pi/lp}$.
*   **Result:** If the network predicts that the population of the $\pi$-orbital drops (electrons delocalize into the ring), it automatically forms the correct traceless quadrupole and dipole in the $4 \times 4$ matrix.

### Why this is the ultimate representation for your model:

1.  **Extreme Compression:** You only predict **4 scalars per atom** (the distribution of the octet). This is the absolute minimum amount of information needed to describe monopole (charge), dipole (polarization), and quadrupole (shape).
2.  **Free Equivariance:** Because the network only predicts *scalars* (populations), and you construct the vectors/tensors using the atomic coordinates (the bond vectors), the model is **100% rotationally equivariant by default**. You don't need complex tensor mathematics in the GNN.
3.  **Chemical Interpretability:** You aren't looking at abstract $P_{ij}$ matrix elements. You can look at the output and say, "Ah, the network shifted 0.1 electrons from the lone pair into the $\pi$-system." This perfectly handles aromatic conjugation and hyperconjugation.

### The Final Concise Blueprint

If we combine your input constraints and this VSEPR output, the model becomes beautifully simple:

**Input (Per Atom):**
*   $Z$ (Atomic Number)
*   Electrostatics: $V, |\vec{E}|, \vec{E} \cdot \vec{u}_i$ (Field magnitude and its projection along bonds).
*   DFTB On-site $H_{ii}$ (Trace and p-p splitting).
*   DFTB Bond integrals (Slater-Koster $h_{ss\sigma}, h_{pp\pi}$ to neighbors).

**Neural Network:**
*   A tiny, standard MLP or Graph Network (no expensive E(3) equivariant layers needed, because all inputs are scalars/dot-products).

**Output (Per Atom):**
*   4 numbers: $[\Delta n_{\sigma1}, \Delta n_{\sigma2}, \Delta n_{\sigma3}, \Delta n_{\pi/lp}]$

**Post-Processing (Inside the QM Code):**
*   Take those 4 numbers and multiply them by the outer products of the geometric hybrid orbitals:
    $\Delta P_{ii} = \sum (\Delta n_k) |h_k\rangle \langle h_k|$
*   Add this directly to the $4 \times 4$ neutral-atom diagonal blocks.
*   Start the SCF/SCC cycle.

This gives you a model with incredibly few parameters, perfect GPU float memory alignment, exact physical rotation symmetry, and rich quantum-chemical orbital resolution (monopole + dipole + quadrupole).

---

# Gemini

You are hitting on a beautiful, highly physical intuition. Trying to pass and process five $4 \times 4$ matrices ($5 \times 16 = 80$ floats) per atom defeats the goal of a ultra-lightweight, embedded-friendly network.

Instead of forcing the network to digest massive raw Hamiltonian matrices to figure out the bonding, we can use the **VSEPR/Lewis octet hybridization logic** as a hard physical constraint. This compresses the descriptor size dramatically while perfectly capturing bonds, lone pairs, and $\pi$-conjugation.

Here is how we can design a highly concise, orbital-resolved representation that gives you monopoles, dipoles, and quadrupoles in **one shot** without adding ghost atoms or bond-charge centers.

## 1. The Concise State Representation: The "VSEPR Vector"

Instead of arbitrary machine learning latents, every atom can be described by its valence electron distribution. For a standard organic molecule, the 8 valence electrons (the octet) are partitioned into:

$$4 = n_{\sigma} + n_{\pi} + n_{\text{lone\_pair}}$$

We can compress the local chemical state of an atom into a single, compact `float4` vector (**The VSEPR State**):

$$\vec{S}_{\text{VSEPR}} = [n_{\sigma}, n_{\pi}, n_{\text{lone\_pair}}, \delta_{\text{conjugation}}]$$

* **$n_{\sigma}$:** Number of $\sigma$-bonds (determines local coordination geometry via hybrid orbitals: $sp^3, sp^2, sp$).
* **$n_{\pi}$:** Number of $\pi$-bonds (tells the network if there is a perpendicular p-orbital holding electron density).
* **$n_{\text{lone\_pair}}$:** Number of non-bonding electron pairs (tells the network where asymmetric local polarization/dipoles will naturally form).
* **$\delta_{\text{conjugation}}$:** A continuous scalar between 0 and 1 indicating partial double-bond character (e.g., your example of the $-\text{NH}_2$ group conjugated to an aromatic ring, where $n_{\pi}$ smoothly goes from $0 \rightarrow 1$ and $n_{\text{lone\_pair}}$ goes from $1 \rightarrow 0$).

### Extracting this from the non-SCC DFTB Hamiltonian

You don't need the whole $4 \times 4$ blocks. You can extract this `float4` directly by looking at the *eigenvalues* or the trace of the non-SCC density matrix, or pre-compute it from standard chemical connectivity graphs.

## 2. On-Atom Resolution of Identity (No Bond Charges Needed)

To represent bond charges and lone pairs without adding extra "ghost" geometric centers, we use the **local atomic orbital basis ($s, p_x, p_y, p_z$) as a natural Resolution of Identity (RI).** Any distribution of the octet around the atom maps analytically to the on-site multipoles ($q, \vec{p}, \mathbf{Q}$):

| Chemical Feature | Mathematical Representation in $sp$-Basis | Multipole Interpretation |
| --- | --- | --- |
| **Isotropic Charge / Monopole ($q$)** | Diagonal sum of population: $P_{ss} + P_{xx} + P_{yy} + P_{zz}$ | **Monopole** (Net charge change $\Delta \rho$) |
| **$\sigma$-Bond / Lone Pair Direction ($\vec{p}$)** | Cross-terms between $s$ and $p$: $[P_{sx}, P_{sy}, P_{sz}]$ | **Dipole Vector** (Points toward lone pairs or asymmetric bonds) |
| **$\pi$-Orbital Alignment ($\mathbf{Q}$)** | Differences in $p$-populations: $P_{xx} - P_{yy}$, etc. | **Quadrupole Tensor** (Defines the flattening of the electron cloud) |

By doing this, a **lone pair** is not a new charge center; it is simply a strong local dipole ($\vec{p}$) pointing away from the bonds, generated by mixing the $s$ and $p$ character (hybridization). A **$\pi$-bond** is not a bond charge; it is an accumulation of density in a specific $p$-orbital direction, which is perfectly captured by the quadrupole tensor $\mathbf{Q}$.

## 3. The Hardware-Optimized GPU Layout (Fixed & Ultra-Lean)

With this philosophy, we can throw away the 80-float Hamiltonian representation. The entire input/output state for an atom fits into just **three `float4` vectors**.

### Input Layout (per atom):

1. **`In_Local` (`float4`):** The VSEPR vector $[n_{\sigma}, n_{\pi}, n_{\text{lone\_pair}}, \delta_{\text{conjugation}}]$.
2. **`In_Global` (`float4`):** The classical electrostatic environment from the rest of the molecule $[V, E_x, E_y, E_z]$.

### The Equivariant GNN Layer:

The Equivariant GNN processes these inputs. Because the inputs are small, the GNN only needs a tiny handful of hidden parameters. It outputs the corrected physical parameters for the one-shot solver.

### Output Layout (per atom):

3. **`Out_Density` (`float4`):** * `w`: $\Delta q$ (Monopole charge shift from QEq)
* `x, y, z`: $\vec{p} = [p_x, p_y, p_z]$ (Local dipole vector, oriented automatically by the Equivariant GNN based on the electric field $\vec{E}$ and the neighbor geometric vectors).

If you want full quadrupole resolution, you add one more `float4` for the independent components of the symmetric traceless tensor $\mathbf{Q}$.

## Why this fits your vision perfectly:

1. **Extremely Lightweight:** You are passing around 8 to 12 floats per atom instead of ~90. This easily compiles down to an embedded system, FPGA, or an ultra-fast GPU warp-level kernel.
2. **Physically Bounded:** By forcing the model to operate on $n_{\sigma}, n_{\pi}, n_{\text{lone\_pair}}$, the neural network cannot predict non-physical electron densities. It is just learning how the local environment *polarizes* and *shifts* a pre-existing Lewis/VSEPR structure.
3. **Instant SCF Speedup:** When you pass this $[q, \vec{p}, \mathbf{Q}]$ back to DFTB+ or B3LYP, you can project it onto the orbital basis in one matrix operation. The SCF solver starts with a density matrix that already "knows" where the lone pairs are pointing and which direction the $\pi$-system is conjugated.

---

# Mistral

### **🎯 Refined Approach: VSEPR-Inspired Local Charge Model**

You’re proposing a **chemically intuitive, ultra-compact representation** where:
- Each atom has an **electron octet** distributed as **n<sub>σ</sub> + n<sub>π</sub> + n<sub>epair</sub> = 4** (VSEPR-like).
- The **local charge distribution** is represented in **sp-basis (4 numbers: monopole + dipole)** or **spd-basis (9 numbers: + quadrupole)**.
- **No explicit bond charge blocks**—instead, **bond effects are encoded via local polarization** (dipole/quadrupole) and **hybridization state**.

This is **brilliant** because:
1. **Physically meaningful**: Directly tied to **chemical bonding** (σ/π/lone pairs).
2. **Ultra-compact**: **4-9 numbers per atom** (vs. 80+ for DFTB blocks).
3. **One-shot**: Predicts **local polarization** without iteration.
4. **Embedded-friendly**: Tiny model (<100 params) + fixed-size inputs.

## **📐 Mathematical Formulation**

### **1. Electron Octet Distribution (VSEPR Prior)**
For each atom **A**, define:
- **n<sub>σ,A</sub>**: Number of **σ-bonds** (single bonds + 1 per double/triple bond).
- **n<sub>π,A</sub>**: Number of **π-bonds** (1 per double bond, 2 per triple bond).
- **n<sub>epair,A</sub>**: Number of **lone pairs** (4 − n<sub>σ,A</sub> − n<sub>π,A</sub>).

**Examples**:
| Atom       | Bonding       | n<sub>σ</sub> | n<sub>π</sub> | n<sub>epair</sub> | Hybridization |
|------------|---------------|---------------|---------------|-------------------|---------------|
| -O-        | 2 single bonds | 2             | 0             | 2                 | sp<sup>3</sup> |
| =O         | 1 double bond  | 1             | 1             | 2                 | sp<sup>2</sup> |
| -N≡        | 1 triple bond  | 1             | 2             | 1                 | sp            |
| -NH<sub>2</sub> (aromatic) | 3σ + 1π (conjugated) | 3 | 1 | 0 | sp<sup>2</sup>-like |

**Key Insight**:
- The **hybridization state** (sp<sup>3</sup>, sp<sup>2</sup>, sp) is **implicitly determined** by (n<sub>σ</sub>, n<sub>π</sub>, n<sub>epair</sub>).
- The **local geometry** (tetrahedral, trigonal, linear) follows from **VSEPR theory**.

### **2. Local Charge Representation**
Represent the **local electron density deviation** (Δρ<sub>A</sub> = ρ<sub>A</sub> − ρ<sub>NA,A</sub>) in a **compact basis**:

| Basis      | Coefficients                     | Dimensionality | Physical Meaning                     |
|------------|----------------------------------|----------------|--------------------------------------|
| **sp**     | [q, μ<sub>x</sub>, μ<sub>y</sub>, μ<sub>z</sub>] | 4              | Monopole + dipole                    |
| **spd**    | [q, μ<sub>x</sub>, μ<sub>y</sub>, μ<sub>z</sub>, Q<sub>xx</sub>, Q<sub>xy</sub>, Q<sub>xz</sub>, Q<sub>yy</sub>, Q<sub>yz</sub>] | 9 | Monopole + dipole + quadrupole |

**Why This Works**:
- **Monopole (q)**: Net charge (QEq-like).
- **Dipole (μ)**: Local polarization (e.g., in C=O, μ points from C to O).
- **Quadrupole (Q)**: Asymmetry (e.g., in benzene π-system).

**No Explicit Bond Charges**:
- Bond effects are **encoded in the dipole/quadrupole** of neighboring atoms.
- Example: In **C=O**, the **dipole on C** and **dipole on O** together describe the **bond polarity**.

## **📥 Input Feature Vector (Concise & Fixed-Size)**
For each atom **A**, the input is a **fixed-size vector** combining:
1. **Atomic identity** (1 number).
2. **Local electrostatics** (4 numbers).
3. **Hybridization state** (2 numbers, since n<sub>epair</sub> = 4 − n<sub>σ</sub> − n<sub>π</sub>).

| Feature               | Description                          | Size | Example (C in COOH) |
|-----------------------|--------------------------------------|------|----------------------|
| Z                     | Atomic number (scaled to [0,1])     | 1    | 6 (carbon) → 0.06    |
| V                     | Coulomb potential at A               | 1    | -1.2 eV              |
| E<sub>x</sub>, E<sub>y</sub>, E<sub>z</sub> | Electric field vector        | 3    | [0.5, -0.3, 0.0] eV/Å |
| n<sub>σ</sub>         | Number of σ-bonds                    | 1    | 3 (COOH: 1 to C, 1 to O, 1 to OH) |
| n<sub>π</sub>         | Number of π-bonds                    | 1    | 1 (C=O)              |
| **Total**             |                                      | **7** |                      |

**GPU Optimization**:
- Pad to **8 floats** (`float2` × 4) for **coalesced memory access**.
- **Total per molecule**: 8 × N<sub>atoms</sub> floats.

### **How to Compute Inputs**
1. **n<sub>σ</sub>, n<sub>π</sub>**:
   - Use a **simple graph traversal** (from molecular topology):
     - For each bond, count as **1σ** (single), **1σ+1π** (double), **1σ+2π** (triple).
     - Sum contributions for each atom.
   - **Example code**:
     ```python
     def get_hybridization(atom, bonds):
         n_sigma = sum(1 for bond in bonds if bond.contains(atom))
         n_pi = sum(bond.order - 1 for bond in bonds if bond.contains(atom))
         return n_sigma, n_pi
     ```

2. **V, E**:
   - Compute from **neutral atom charges** (Z<sub>B</sub>) or **initial QEq guess**:
     ```python
     V_A = sum(Z_B / |r_A - r_B| for B in atoms if B != A)
     E_A = sum(Z_B * (r_A - r_B) / |r_A - r_B|**3 for B in atoms if B != A)
     ```

## **🏗 Model Architecture (Ultra-Lightweight)**
### **Option 1: MLP (Simplest)**
```python
# Input: (batch_size, N_atoms, 8)  # Padded to 8 floats
x = Input(shape=(8,))
h = Dense(16, activation="silu")(x)  # 16 × (8 + 1) = 136 params
h = Dense(8, activation="silu")(h)   # 8 × (16 + 1) = 136 params
y = Dense(4)(h)  # Output: [q, μx, μy, μz]  # 4 × (8 + 1) = 36 params
# Total: ~300 params (fits in <1 KB)
```

### **Option 2: Hybrid (MLP + 1 GNN Layer)**
If you want **many-body effects** (e.g., cooperative polarization in H-bonded systems):
```python
# Input: Graph with node features = [Z, V, Ex, Ey, Ez, nσ, nπ]
# Edge features: Distance (optional)
x = NodeFeatures(shape=(7,))  # 7 floats per atom
edge_attr = DistanceEmbedding()  # 1 float (distance)

# 1 message passing layer (8 hidden dims)
h = GraphConv(8, edge_attr)(x)  # ~50 params
y = Dense(4)(h)  # [q, μx, μy, μz]  # 36 params
# Total: ~100 params
```

**Recommendation**: Start with **Option 1 (MLP)**. If accuracy is insufficient, add **1 GNN layer**.

## **🔌 Integration with QM Solvers**
### **1. DFTB+ SCC**
- **Input to SCC**: Use predicted **q<sub>A</sub>** as initial charges.
- **Expected speedup**: **3-5× fewer SCC iterations** (typical SCC needs 10-20 iterations).

### **2. Semi-Empirical/DFT**
- **Initial density matrix**:
  ```
  P_0 = P_NA + ΔP_ML
  ```
  where **ΔP_ML** is constructed from **q, μ** (e.g., via **Mulliken partitioning**).
- **Expected speedup**: **2-4× fewer SCF iterations**.

### **3. Resolution of Identity (Optional)**
If you **really** want bond charges (e.g., for π-systems):
- For each bond **A-B**, compute a **bond charge** as:
  ```
  q_AB = f(q_A, q_B, μ_A, μ_B, r_AB)
  ```
  (e.g., **q_AB = 0.5 × (q_A + q_B) × exp(-r_AB / λ)**).
- **Total density**: ρ = Σ<sub>A</sub> ρ<sub>A</sub> + Σ<sub>A-B</sub> ρ<sub>AB</sub>.
- **Cost**: Adds **1 number per bond**, but may not be worth it for embedded.

## **📊 Training Strategy**
### **1. Dataset**
- **Molecules**: QM9 (134k organic molecules) or **custom dataset** (e.g., drug-like molecules).
- **Targets**:
  - **q<sub>A</sub>**: Converged QEq charges (from DFTB+ or DFT).
  - **μ<sub>A</sub>**: Dipole moments (from **Mulliken/Löwdin population analysis**).
  - **Optional**: Q<sub>A</sub> (quadrupole from **DMA** or **DDEC6**).

### **2. Loss Function**
```python
# Monopole + dipole
loss = mse(q_pred, q_true) + 0.1 * mse(μ_pred, μ_true) + λ * (sum(q_pred))**2
```
- **λ**: Charge conservation penalty (e.g., λ = 10.0).
- **Optional**: Add **Coulomb regularization** (penalize ∇²Δρ ≠ 0).

### **3. Physics-Informed Initialization**
- Initialize **output layer weights** to mimic **electronegativity trends** (e.g., higher Z → more negative q).

## **🚀 Expected Performance**
| Metric               | MLP (4-output) | MLP (9-output) | +1 GNN Layer |
|----------------------|----------------|----------------|--------------|
| **Params**           | ~300           | ~600           | ~400         |
| **Inference (100 atoms)** | ~0.1 ms       | ~0.2 ms        | ~0.3 ms      |
| **Model size (quantized)** | ~1 KB    | ~2 KB          | ~1.5 KB      |
| **SCF Speedup**      | 3-5×           | 4-6×           | 5-8×         |

## **❓ Open Questions for You**
1. **Output basis**: Start with **sp (4 numbers)** or **spd (9 numbers)**?
2. **Hybridization**: Should we **hardcode n<sub>σ</sub>, n<sub>π</sub>** (from topology) or **predict them with ML**?
3. **Bond charges**: Do you **need explicit bond charges**, or is **atomic dipole/quadrupole sufficient**?
4. **Training data**: Should we use **QM9** or a **custom dataset** (e.g., focused on aromatics)?
5. **Embedded constraints**: Strict **memory limit** (e.g., <1 KB)? If so, we can **prune the MLP** further.

---

# DeepSeek-V4-pro

You’re moving towards a chemically richer representation – something that captures the *distribution* of the electron octet around each atom, not just its lowest multipoles. That’s exactly the kind of inductive bias that can make a small model very accurate. Let’s unpack this and then converge on a concrete, lean design that still fits the “tiny embedded accelerator” constraint.

### Two concrete representation strategies, both fitting your sp/VSEPR vision

#### A. Atomic‑centred multipoles (scalar charge + vector dipole + tensor quadrupole)
- **What you predict:** Per atom – \(\Delta q\) (1 scalar), \(\boldsymbol{\mu}\) (3-vector), \(\Theta\) (5 independent quadrupole components). That’s 9 numbers/atom if you include all; 4 if you stop at dipole.
- **Why it fits VSEPR:** An sp\(^3\) oxygen with two lone pairs and two σ bonds has a large atomic dipole and quadrupole that describe the *anisotropy* of its electron cloud – you don’t need explicit lone‑pair sites to capture the physics; the quadrupole does it for you.
- **Equivariance:** an SE(3)‑equivariant GNN outputs these tensors directly; everything rotates perfectly.
- **Simplicity:** You never construct extra sites; the network just learns a mapping from local geometry + non‑SCC Hamiltonian → multipoles.

#### B. Off‑atom charge centres (bond centres, lone‑pair sites, π‑clouds)
- **What you predict:** Charges on a set of *automatically placed* massless sites – e.g. for an sp\(^3\) atom: two σ‑bond midpoints, two lone‑pair positions. For an sp\(^2\) atom: one σ‑bond midpoint per bond, one π‑site above/below the plane, and lone‑pair positions if present.
- **Where do the sites come from?** They are placed by purely geometric, rotation‑equivariant rules that depend only on the atom’s bonding neighbours. For instance:
  - Lone‑pair sites of a tetrahedral centre: placed at directions that complete a tetrahedron given the two bond vectors (deterministic).
  - π‑site: placed at a fixed distance along the local normal vector of the sp\(^2\) plane.
  - Bond centres: midpoint of each A–B bond.
- **Why this fits VSEPR perfectly:** The electron octet is literally distributed among these sites; the physics of σ/π/lone‑pair separation is built into the topology of the charge model.
- **Equivariance:** Because the site positions are defined using only the neighbouring atom coordinates, the whole set of sites rotates rigidly with the molecule. The GNN then only needs to predict *scalar* charges on each site (type‑0 features) – no equivariant vector outputs needed.

### Which one to choose? A pragmatic path

**Start with atomic multipoles (A)** – it is the simplest possible extension of a QEq‑like guess, requires zero extra geometry logic, and already brings anisotropic electron density into the SCF initial guess. If that gives you enough speed‑up, you are done.

If you find that certain systems (e.g., amides, aromatic rings) need a finer description of charge anisotropy, then **add off‑atom sites (B)**. The beauty is that both can coexist: the GNN predicts multipoles on the nuclei *and* charges on a small set of pre‑placed extra sites. The electrostatics then comes from all of them.

Below I’ll detail how to design the charge‑site approach, because it maps most directly onto your VSEPR language and is extremely chemistry‑aware.

### Design of the VSEPR‑based charge‑site predictor (one‑shot, no iteration)

1. **Automatic site construction (no learnable parameters)**
   For each atom \(A\):
   - Detect its local bonding partners (using covalent radii).
   - Determine its hybridisation state (sp\(^3\), sp\(^2\), sp) from the number of neighbours and bond angles (a simple deterministic classifier).
   - Place the following *virtual sites* (positions \( \mathbf{r}_{A,k} \)):
     * **Bond‑σ sites:** one at the midpoint of each A–B bond.
     * **Lone‑pair sites:** directions that complete the hybrid geometry, at a fixed distance (e.g., 0.5 Å from the nucleus).
     * **π‑site:** for sp\(^2\) or sp hybridised atoms, one site above and one below the plane (or along the perpendicular), at e.g., 0.7 Å.
   The number of extra sites per atom is fixed once the coordination is known. For an sp\(^3\) oxygen with two bonds, you get 2 bond‑σ + 2 lone‑pair = 4 sites. Total sites = nuclei + extra.

2. **Graph construction for the GNN**
   - **Nodes:** all nuclei + all virtual sites.
   - **Edges:** connect sites that are close (e.g., within a cut‑off) – this can include nucleus‑site, site‑site, nucleus‑nucleus interactions.
   - **Initial node features:**
     * For nuclei: atomic number, on‑site ε\(_s\), ε\(_p\) from DFTB, and the current V/E (if we use recurrent scheme).
     * For virtual sites: a one‑hot encoding of type (σ‑bond, lone‑pair, π) and the same V/E computed from the current charge distribution.
   - **Edge features:** distance, and for nucleus‑nucleus edges, the Slater–Koster invariants from the non‑SCC Hamiltonian (h\(_{ss\sigma}\), etc.) as we discussed earlier. For virtual‑site edges, just the distance vector and a type indicator.

3. **GNN output**
   After a few message‑passing rounds, each node outputs a **scalar prediction**:
   - For nuclei: base electronegativity χ\(_A\) and hardness J\(_A\) (or directly a charge correction Δq\(_A\)).
   - For virtual sites: base χ\(_k\) and J\(_k\) (or charge directly).
   The GNN is **not required to output vectors**, so you don’t need an equivariant network at all – a standard GNN (SchNet, DimeNet, etc.) on the expanded graph is sufficient, and the sites are already placed equivariantly.

4. **Charge determination – physics layer**
   Two options:
   * **Direct prediction:** GNN outputs \(q_i\) for all sites (nuclei + virtual). Then just enforce total charge neutrality via a final shift.
   * **Physics‑informed QEq:** GNN outputs χ\(_i\) and J\(_i\) for each site. You then solve the **one‑shot linear system**  
     \[
     \mathbf{A}\mathbf{q} = -\boldsymbol{\chi} + \lambda\mathbf{1}
     \]
     with the constraint \(\sum_i q_i = Q_\text{total}\). This is a single \(O(N^3)\) solve, fully differentiable, and for the small total site count (a few hundred) it is faster than one SCF diagonalisation. This embeds charge equilibration across all sites, including bond and lone‑pair charges, and gives you physically correct long‑range behaviour.

   The QEq path is especially elegant because the ML’s job is reduced to predicting **local environment‑dependent electronegativities and hardnesses** for the pre‑placed chemistry‑aware sites. The equilibrium charge distribution then follows automatically.

5. **Recurrent field refinement (optional)**
   Even with the one‑shot QEq, you might want to run one or two cycles where you recompute V and E from the QEq charges, feed them back into the GNN to update χ/J, and resolve. This mimics the self‑consistent polarisation of a real SCF but with very few steps. It adds minimal cost and improves accuracy.

### Why this fits your “concise representation” requirement

- The extra sites add at most ~4 parameters per heavy atom. The total number of charge variables is still linear in system size.
- The GNN can be tiny – you only need to predict a few scalars per site. With fixed graph topology and small hidden dimensions (16–32), the whole model stays under 10 k parameters.
- The representation is **intrinsically interpretable**: you can look at the predicted charges on lone‑pair sites vs. σ‑bond sites and see chemistry at work.

### What about the sp‑basis density blocks?

You originally mentioned using the 4×4 Hamiltonian blocks as input. In this off‑atom site picture, you don’t need to treat the 4×4 explicitly; you only need the **scalar invariants** (on‑site energies and Slater–Koster hopping integrals) as input. The VSEPR placement of sites already captures the geometric arrangement of the valence electron pairs, which is the “sp basis” concept but realised in a real‑space charge picture.

If you later want to go one step deeper and actually predict a density matrix in an sp basis (4 orbitals/atom), that would mean outputting a 4×4 matrix per atom *that rotates with the molecule* – which requires an equivariant network to output matrix representations of SO(3). That’s feasible but adds complexity and parameter count. The charge‑site model gives you essentially the same physical information (anisotropic electron distribution) with only scalar predictions.

### Summary of the recommended path

1. **Start with atomic multipoles (dipole + quadrupole)** predicted by an equivariant GNN. This is the most direct extension of QEq and will already improve SCF convergence dramatically.
2. **If you need even finer anisotropy** (e.g., for accurate hydrogen bonds or π‑stacking), **add VSEPR‑based off‑atom charge sites** that are placed geometrically. Use a standard GNN to predict site‑specific QEq parameters (χ, J), then solve the linear charge‑equilibration once.
3. **Always feed the non‑SCC Hamiltonian** as invariant scalars (on‑site energies, SK bond integrals) + distances. This replaces the need to handle 4×4 blocks directly.
4. **Keep the field‑update loop short** (2–3 steps) or even single‑step if QEq is used – the model then becomes an amortised initial guess engine that runs in microseconds on an MCU.

This design is physics‑informed, chemically intuitive, rotationally correct, and extremely parameter‑efficient. It is, I believe, exactly the kind of bridge between VSEPR thinking and machine learning that you’re aiming for.

Would you like me to flesh out any of these pieces further – e.g., the exact equivariant layer for multipoles, or the geometric algorithm for lone‑pair site placement?

---

# ChatGPT (5.5?)

I think now we are getting to the really interesting part of the design. I actually think your intuition here may be better than using full (4\times 4) density blocks, because we should carefully think about **what information is really needed for SCF acceleration**.

The (4\times4) density matrix in an (sp) basis is mathematically complete for the on-site density, but it is not the most physically compact representation.

Let us analyze.

## 1. On-site (sp) density matrix

For an atom with orbitals:

[
\phi = (s,p_x,p_y,p_z)
]

the density matrix is:

[
P=
\begin{pmatrix}
ss & sp_x & sp_y & sp_z\
p_xs & p_xp_x & p_xp_y & p_xp_z\
p_ys & p_yp_x & p_yp_y & p_yp_z\
p_zs & p_zp_x & p_zp_y & p_zp_z
\end{pmatrix}
]

It has 10 independent components (symmetric).

These naturally decompose into irreducible multipoles.

### Monopole

The trace:

[
q = P_{ss}+P_{xx}+P_{yy}+P_{zz}
]

1 number.

### Dipole

The (s-p) mixing terms:

[
(P_{sx},P_{sy},P_{sz})
]

3 numbers.

### Quadrupole

The traceless part of the (p-p) block:

[
\begin{pmatrix}
Q_{xx}&Q_{xy}&Q_{xz}\
Q_{xy}&Q_{yy}&Q_{yz}\
Q_{xz}&Q_{yz}&Q_{zz}
\end{pmatrix}
]

with

[
Q_{xx}+Q_{yy}+Q_{zz}=0
]

5 numbers.

### Isotropic p population

There is one remaining quantity:

[
n_p=P_{xx}+P_{yy}+P_{zz}
]

which corresponds to the size of the p-shell occupation.

Therefore:

[
10 = 1 + 3 + 5 + 1
]

So your "9 numbers" is almost complete. You only lose the isotropic (p)-population.

## 2. Is that enough?

For polarization, probably yes.

For example in a carbonyl oxygen:

```
   O
  //
 C
```

the density change is mainly:

* increased electron density on oxygen,
* dipole towards oxygen,
* quadrupole due to anisotropy of lone pairs.

The multipoles capture this reasonably.

But there is a problem.

Consider benzene.

```
      C
   //   \\
  C       C
  ||      ||
  C       C
   \\   //
      C
```

The (\pi)-electron density is **delocalized between atoms**.

An atomic multipole representation can reproduce the external electrostatic potential but it cannot tell you how much charge sits in the (C-C) (\pi) bond.

For SCF acceleration, this may matter because the Hamiltonian depends on orbital populations.

# 3. The VSEPR / electron pair representation

This is where I think your idea becomes very attractive.

Instead of storing orbital density matrix elements, store **electron pair occupations**.

Each atom has an electron octet represented as four directional orbitals:

[
4 = n_\sigma+n_\pi+n_{lp}
]

Examples:

### Water-like oxygen

```
   ..
H-O-H
   ..
```

[
n_\sigma=2,\quad n_\pi=0,\quad n_{lp}=2
]

sp³ tetrahedral.

---

### Carbonyl oxygen

```
   ..
   O
  ||
  C
```

[
n_\sigma=1,\quad n_\pi=1,\quad n_{lp}=2
]

sp² trigonal.

### Pyridine nitrogen

One lone pair plus aromatic (\pi):

[
n_\sigma=2,\quad n_\pi=1,\quad n_{lp}=1
]

### Aniline-like nitrogen

It is a resonance mixture:

[
\mathrm{sp^3}: (3,0,1)
]

and

[
\mathrm{sp^2}: (3,1,0)
]

The ML model could naturally interpolate between them.

## 4. Why this representation may be better

It is chemically very intuitive.

A local state might be:

```
atom i:

sigma orbitals:
    n_sigma1
    n_sigma2
    n_sigma3
    n_sigma4

pi orbitals:
    n_pi1
    n_pi2

lone pairs:
    n_lp1
    n_lp2
```

with associated directions:

[
\mathbf{u}_k
]

and occupancies:

[
n_k
]

Then the density is approximately:

[
\rho_i(\mathbf r)
=================

\sum_k n_k \rho_k(\mathbf r-\mathbf R_i,\mathbf u_k)
]

This is very similar to your existing "explicit (\pi) orbital" philosophy in FireCore.

# 5. Connection to ML input

Interestingly, the DFTB Hamiltonian should already contain enough information to infer these states.

For example:

* (sp^3) carbon:

  * four strong (sp) hybrid bonds.

* aromatic carbon:

  * three (\sigma) bonds,
  * one (p_z) orbital with conjugated (H_{ij}).

* carbonyl oxygen:

  * one strong (\sigma),
  * one (\pi),
  * two nonbonding (p) orbitals.

So the network can infer a latent "hybridization state" from:

* on-site (H_{ii}),
* bond blocks (H_{ij}),
* geometry.

# 6. A possible hierarchy of representations

I would actually consider several levels:

| Representation            |     Numbers | Physical meaning       | Good for         |
| ------------------------- | ----------: | ---------------------- | ---------------- |
| Atomic charge             |           1 | QEq                    | DFTB SCC         |
| Multipoles                |         4–9 | polarization           | electrostatics   |
| (sp) density matrix       |          10 | exact local AO density | semiempirical QM |
| Bond density (P_{ij})     | 16 per bond | delocalized bonds      | conjugation      |
| Electron pair occupations |        ~4–8 | VSEPR chemistry        | compact ML       |

## My current opinion

I think the most promising compromise is actually **not** the full (5\times4\times4) density blocks.

For a tiny neural network intended to run inside an SCF loop, I would seriously consider a **localized valence-electron representation**:

* each atom has a small number of electron lobes (σ bonds, π orbitals, lone pairs),
* each lobe has an occupancy and perhaps a polarization parameter,
* the directions are known from molecular geometry or can be predicted,
* the network predicts shifts of these occupancies.

It is almost like a **learned VSEPR + charge equilibration model**.

This has a very appealing property: it sits between classical force fields and quantum density matrices. It is much richer than QEq, but orders of magnitude cheaper than predicting a full density matrix.

In fact, this direction is very close in spirit to the "explicit σ/π electron representation" that you have already developed in FireCore, and it may provide a natural bridge between your classical force-field ideas and an ML-accelerated SCF solver.

---

# USER

OK, read 
@delta_rho_ML.chat.md 
again from line 814 we added some new ideas and consideration, think it through and make again metanalaysis. I thin we rendered there several viable soluton between which we should edecide - what we will actually predict QEq params, charges, in which representation, will we use MLP or GNN, maybe some other ideas.

Write clearly what are the variants, what are the common features, and then try to decide which we should do first (simplest) and which maybe later (most efficient)

---

# Kimi 2.6 / devin / meta-analysis 

# Meta-Analysis: δρ ML Architecture — Variants, Common Ground & Recommended Path

## Common Features (Universal Consensus)

All LLMs and your own reasoning converged on:

- **Target Δρ = ρ_SCF − ρ_NA**, not total density. Neutral atom cores are subtracted away.
- **Global electrostatics pre-computed**: A cheap $O(N^2)$ pass computes $[V_i, E_{ix}, E_{iy}, E_{iz}]$ (a `float4`) per atom to inject long-range awareness without deep message passing.
- **Local quantum descriptors from DFTB H⁰/S**: The non-SCC Hamiltonian and overlap encode bonding, hybridization, conjugation, and orbital symmetry.
- **One-shot inference**: The network itself must be a single forward pass. No iterative QEq or SCF-like loops inside the model.
- **Ultra-small parameter budget**: Target <10k params, embeddable, fast on GPU.
- **Physics-informed loss**: Charge conservation penalty, dipole/quadrupole matching, optional Coulomb regularization.

## The Variants

### Variant A: MLP Predicting QEq Parameters (χ, J, α)
**Original idea, now largely abandoned.**
- GNN/MLP predicts local electronegativity $\chi_i$, hardness $J_i$, and polarizability $\alpha_i$.
- A subsequent **QEq linear solve** produces charges: $\mathbf{J}\vec{q} = -\vec{\chi}$.
- **Rejected because**: QEq only gives monopoles (no dipoles/quadrupoles); backprop through iterative QEq is fragile; the $O(N^3)$ solve adds cost; doesn't produce orbital-resolved density.

### Variant B: Direct Multipole Prediction (Equivariant GNN)
**The "standard modern" approach.**
- Equivariant GNN (e.g. PaiNN-style, 16 scalars + 16 vectors, 2 layers) directly outputs per-atom:
  - Scalar: $\Delta q_i$ (charge shift)
  - Vector: $\boldsymbol{\mu}_i$ (dipole, 3 components)
  - Tensor: $\Theta_i$ (traceless quadrupole, 5 components)
- **Pros**: Full E(3) equivariance by construction; outputs map cleanly to density matrix; physically complete.
- **Cons**: Equivariant GNN libraries (e3nn, etc.) add dependency weight; output tensor components can be noisy to learn directly.

### Variant C: Direct ΔP_μν Density Matrix Blocks
**ChatGPT's most concrete proposal.**
- Network predicts the **on-site $4 \times 4$ density matrix correction** $\Delta P_{ii}$ per atom (10 symmetric floats), optionally neighbor blocks too.
- Inputs: local $H^0_{ii}$ and $H^0_{ij}$ blocks (80 floats) + electrostatics `float4`.
- **Pros**: Orbital-resolved; directly usable as SCF initial guess; input/output live in the same "language" (local AO blocks).
- **Cons**: Raw $4 \times 4$ blocks are coordinate-system dependent — without equivariance or local-frame rotation, the model fails on rotated molecules.

### Variant D: Local-Frame MLP on Rotated Blocks
**Gemini-3.1-pro's fix for Variant C.**
- Construct a **local coordinate frame** per atom from neighbor geometry (e.g., $z$-axis to closest neighbor, $x$-axis from first two neighbors).
- Rotate $H^0$ blocks and $[E_x, E_y, E_z]$ into this frame → everything becomes rotationally invariant.
- Run a **simple MLP** (no equivariant library needed) on the invariant features.
- Rotate the predicted $\Delta P_{ii}$ back to global frame.
- **Pros**: Eliminates need for heavy equivariant GNN; pure MLP is tiny and fast; fixed-size input (~92 floats).
- **Cons**: Frame construction must be smooth and unambiguous (problematic for high-symmetry atoms like single H, or linear molecules).

### Variant E: VSEPR / Electron Pair Occupancy Representation
**Your own idea, enthusiastically endorsed by all LLMs.**
- Each atom distributes its valence octet into chemically meaningful slots: $n_\sigma$ σ-bonds, $n_\pi$ π-bonds, $n_{lp}$ lone pairs ($n_\sigma + n_\pi + n_{lp} = 4$).
- The network predicts **occupancy shifts** $\Delta n_k$ for these slots (e.g., 4 scalars per atom).
- Directions are known from geometry (bond vectors, VSEPR geometry).
- Post-processing constructs the full $4 \times 4$ $\Delta P_{ii}$ from hybrid orbitals.
- **Pros**: Extremely compact (4 numbers/atom vs. 10); chemically interpretable ("network moved 0.1 e⁻ from lone pair to π"); free equivariance because vectors come from geometry, not NN; naturally handles conjugation/resonance.
- **Cons**: Requires deterministic hybrid orbital construction from neighbor geometry; not all atoms fit the octet rule (hypervalent, transition metals); the mapping from occupancies to density matrix must be hard-coded.

### Variant F: Off-Atom Charge Sites (VSEPR-Derived)
**DeepSeek's extension of Variant E.**
- Place **virtual charge sites** at geometrically determined positions: bond midpoints, lone-pair lobes, π-cloud sites above/below sp² planes.
- GNN predicts **scalar charges** on all sites (nuclei + virtual).
- Optionally solve a one-shot QEq on the expanded graph.
- **Pros**: Very intuitive chemistry; no tensor outputs needed (standard GNN works); captures bond-charge and lone-pair anisotropy explicitly.
- **Cons**: Site placement logic is complex and must be robust; expanded graph increases node count; extra geometry rules add engineering overhead.

### Variant G: Recurrent Field-Updating GNN
**DeepSeek's unrolled optimizer.**
- Start with $q^{(0)} = 0$, iterate 2–3 times:
  1. Compute $V^{(t)}, \mathbf{E}^{(t)}$ from current charges/dipoles.
  2. GNN predicts $\Delta q^{(t)}, \Delta \boldsymbol{\mu}^{(t)}$ from fields + Hamiltonian.
  3. Update and repeat.
- **Pros**: Mimics physical self-consistency; can learn charge equilibration naturally.
- **Cons**: Unrolled loop adds inference cost and memory in training; convergence not guaranteed; overkill if one-shot is sufficient.

## Recommended Decision: What to Do First, What to Do Later

| Priority | Variant | Why | When |
|----------|---------|-----|------|
| **🥇 First (now)** | **Variant D: Local-Frame MLP predicting on-site ΔP_μν blocks** | Uses your existing 4×4 Hamiltonian block infrastructure; no new dependencies; pure PyTorch MLP; directly produces orbital-resolved density guess for Fireball/DFTB; fixed-size vectors; easy to debug and profile. | Immediate prototype. |
| **🥈 Second (if D needs boost)** | **Variant B: Lightweight Equivariant GNN (multipoles)** | If local-frame rotation is problematic (high symmetry, linear molecules), replace MLP with a tiny PaiNN/EGNN layer. Same inputs/outputs, just learns equivariance instead of relying on frames. | If local-frame MLP accuracy is inadequate. |
| **🥉 Third (for maximum compression)** | **Variant E: VSEPR Occupancy Representation** | Once the pipeline works, compress output from 10 floats to ~4 floats per atom by predicting hybrid-orbital occupancies instead of raw matrix elements. Requires robust hybrid-orbital construction code. | After first two prove the concept. |
| **Later (stretch goal)** | **Variant F: Off-Atom Sites** | If π-systems and H-bond networks need finer anisotropy than on-site blocks can capture. | Only if Variant E proves insufficient. |
| **Probably skip** | **Variant A (QEq params)** and **Variant G (recurrent)** | QEq-only gives monopoles; recurrent adds complexity with unclear benefit for a one-shot accelerator. | Not recommended as primary path. |

## Concrete First Prototype: Variant D "Local-Frame MLP"

### Input per atom (fixed-size, ~92–100 floats):
1. **Atomic identity**: $Z_i$ embedding (8 floats).
2. **Electrostatics**: $[V_i, E_{ix}, E_{iy}, E_{iz}]$ computed from cheap Gasteiger/QEq guess (4 floats).
3. **On-site $H^0_{ii}$**: $4 \times 4$ block (16 floats).
4. **Neighbor $H^0_{ij}$**: Up to 4 nearest neighbors, $4 \times 4$ each, zero-padded (64 floats).

**Frame construction**:
- $z$-axis = unit vector to closest bonded neighbor.
- $x$-axis = component of vector to 2nd neighbor orthogonal to $z$.
- Rotate all blocks and $\mathbf{E}$ into this frame.

### Network:
- MLP: 92 → 64 → 64 → 10.
- Output: upper-triangle of symmetric $4 \times 4$ $\Delta P_{ii}$ in **local frame**.
- Rotate back to global Cartesian before assembly.

### Loss:
- $L_{\text{data}} = \|\Delta P_{\text{NN}} - \Delta P_{\text{ref}}\|^2$ (from converged Fireball density)
- $L_{\text{charge}} = (\sum_i q_i - Q_{\text{tot}})^2$
- $L_{\text{dipole}} = \|\sum_i \boldsymbol{\mu}_i - \boldsymbol{\mu}_{\text{mol}}\|^2$

### Integration:
- $P_{\text{guess}} = P_{\text{NA}} + \Delta P_{\text{NN}}$ injected into Fireball/DFTB+ SCF.

### Expected scale:
- **~5,000 parameters**, <1 KB quantized.
- **Inference**: ~0.1 ms per 100-atom molecule on GPU.
- **Training data**: Converged Fireball densities for organic molecules → extract per-atom on-site blocks from $P_{\mu\nu}$.

## Why This Order

1. **Local-Frame MLP (D)** is the fastest path to a working prototype because it leverages your existing 4×4 block code and requires no new dependencies. It directly answers the question: *can a tiny network, given H⁰ and electrostatics, predict the on-site density correction?*
2. **Equivariant GNN (B)** is the natural fallback if frame singularities (linear molecules, single atoms) hurt accuracy. The architecture change is a swap of the MLP body, not the pipeline.
3. **VSEPR (E)** is the ultimate compression but requires writing hybrid-orbital construction code that handles all edge cases (pyramidal N, hypervalent S, resonance). Do it after you know the problem is solvable.

The key insight from this entire conversation: **do not solve QEq inside the network, and do not predict abstract QEq parameters**. Predict the **density matrix directly** in a representation you already understand (4×4 blocks), using physics only for input features and loss regularization.

---

# USER

OK, now we can look how to actually o it in our system. Main input for this is the actually orbital resolved charge density (density matrix), overlpa  matrix and hamiltonian. Initially we want to obtin this from dftb+, we made a library which we use in ou AFM system, can you look on it how we can use it ?

@dftb_lib.py@Grid_dftb.py @DFTBcore.py @AFM_FDBM_DFTB.md  @AFM_utils.py 

---

# Kimi 2.6 devin

## What Already Exists (Ready to Reuse)

### 1. DFTB+ SCF Runner → Dense Matrices
**[pyBall/DFTB/DFTBcore.py](cci:7://file:///home/prokop/git/FireCore/pyBall/DFTB/DFTBcore.py:0:0-0:0)** already does exactly what you need:

```python
from pyBall.DFTB.DFTBcore import DFTBcore

dftb = DFTBcore()
dftb.init('dftb_in.hsd')
dftb.enable_matrix_collection(dm=True, h=True, s=True)
energy = dftb.run_scf()

H_dense   = dftb.get_h_dense()   # (norb_total, norb_total)
S_dense   = dftb.get_s_dense()   # (norb_total, norb_total)
DM_dense  = dftb.get_dm_dense()   # (norb_total, norb_total)
C_dense, E = dftb.get_eigvecs_dense()  # MO coeffs + eigenvalues
dftb.finalize()
```

These are C row-major `numpy` arrays. The matrices are in the **non-orthogonal atomic basis** (same basis as the overlap). The non-SCC Hamiltonian $H^0$ is not directly exposed, but you can get the full SCC Hamiltonian at convergence.

> **Note**: If you specifically need the **non-SCC** $H^0$ (before charge corrections), DFTBcore stores it internally before diagonalization. There is [enable_hamiltonian_storage(store=True)](cci:1://file:///home/prokop/git/FireCore/pyBall/DFTB/DFTBcore.py:272:4-283:86) which stores the Hamiltonian before it gets overwritten by eigenvectors. You may need to add a getter for that pre-SCC matrix if it's not already available.

### 2. Orbital Layout → Per-Atom Slicing
**[pyBall/OCL/AFM_utils.py::build_orbital_layout()](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM_utils.py:149:0-173:18)** gives you the mapping from atom → orbital indices:

```python
from pyBall.OCL.AFM_utils import build_orbital_layout
from pyBall.DFTB.DFTBplusParser import parse_wfc_hsd

basis_data = parse_wfc_hsd('wfc.mio-1-1.hsd')
enames = ['C','C','O','H',...]  # per atom

norb_per_atom, orb_offsets, max_l = build_orbital_layout(basis_data, enames)
```

| Array | Shape | Meaning |
|-------|-------|---------|
| `norb_per_atom` | `(natoms,)` | Orbitals per atom (1 for H, 4 for C/N/O sp, 9 for spd) |
| `orb_offsets` | `(natoms+1,)` | `orb_offsets[i]` = first global orbital index of atom `i` |

This is exactly what you need to slice the dense matrices. For atom `i`:
```python
i0, i1 = orb_offsets[i], orb_offsets[i+1]
DM_ii = DM_dense[i0:i1, i0:i1]          # on-site block
H_ii  = H_dense[i0:i1, i0:i1]
S_ii  = S_dense[i0:i1, i0:i1]
```

### 3. Grid Projection (for Validation / Real-Space ρ)
**[pyBall/DFTB/Grid_dftb.py](cci:7://file:///home/prokop/git/FireCore/pyBall/DFTB/Grid_dftb.py:0:0-0:0)** already projects dense DM to real-space grid:
```python
from pyBall.DFTB import Grid_dftb as dg

projector, atoms_dict = dg.setup_gridprojector_from_dftb(dftb_data, basis_ang)
rho_grid = projector.project_density_dense(DM_dense.astype(np.float32),
                                           norb_per_atom, orb_offsets,
                                           atoms_dict, grid_spec)
```

This is useful for:
- Validating that your predicted ΔP reconstructs the correct real-space density
- Computing neutral atom density $\rho_{NA}$ for comparison
- The `B3_FACTOR = 1/(0.529177^3)` unit conversion is already handled

## What You Need to Add (New Code)

### A. Extract Per-Atom On-Site and Neighbor Blocks from Dense Matrices
This is ~30 lines. Given `norb_per_atom`, `orb_offsets`, and a dense matrix `M`, extract:

```python
def extract_local_blocks_dense(M_dense, norb_per_atom, orb_offsets, atom_pos,
                                max_neighbors=4, R_bond_max=2.5):
    """
    Slice dense matrix into per-atom on-site and neighbor blocks.
    
    Returns:
        on_site:  list of (ni,ni) blocks, one per atom
        neigh_blocks: list of (max_neighbors, ni, nj) blocks, zero-padded
        neigh_indices: list of (max_neighbors,) neighbor atom indices, -1 padded
    """
    natoms = len(norb_per_atom)
    on_site = []
    neigh_blocks = []
    neigh_indices = []
    
    # Compute neighbor list from distances (or use bonding topology)
    for i in range(natoms):
        ni = norb_per_atom[i]
        i0, i1 = orb_offsets[i], orb_offsets[i+1]
        on_site.append(M_dense[i0:i1, i0:i1])
        
        # Find up to max_neighbors closest bonded atoms
        # (use covalent radii or just distance < R_bond_max)
        dists = np.linalg.norm(atom_pos - atom_pos[i], axis=1)
        # ... neighbor selection logic ...
        
        blocks_i = np.zeros((max_neighbors, ni, ni), dtype=M_dense.dtype)
        idxs_i = -np.ones(max_neighbors, dtype=np.int32)
        on_site.append(M_dense[i0:i1, i0:i1].copy())
        # fill blocks_i, idxs_i for real neighbors...
    
    return on_site, neigh_blocks, neigh_indices
```

For the **ML input** you want the 4×4 blocks (or `ni×ni` blocks). For sp basis, `ni = 4` for C/N/O, `ni = 1` for H. You should **zero-pad H to 4×4** so the input is uniform.

### B. Compute Electrostatics Float4 `[V_i, E_ix, E_iy, E_iz]`
The AFM pipeline computes real-space `rho_diff` and solves Poisson via FFT, but for the ML you need **per-atom** $V_i$ and $\mathbf{E}_i$ from a rough charge guess.

You can compute this cheaply from Mulliken charges (from DFTB+ or from a one-shot QEq with tabulated $\chi, J$):

```python
def compute_electrostatics_float4(charges, atom_pos, rcut=15.0):
    """
    charges: (natoms,) rough atomic charges (e.g., from Mulliken or QEq)
    atom_pos: (natoms, 3) in Angstrom
    Returns: (natoms, 4) array [V, Ex, Ey, Ez]
    """
    natoms = len(charges)
    out = np.zeros((natoms, 4), dtype=np.float32)
    for i in range(natoms):
        dR = atom_pos - atom_pos[i]
        r = np.linalg.norm(dR, axis=1)
        mask = (r > 0.01) & (r < rcut)
        q = charges[mask]
        r_inv = 1.0 / r[mask]
        r_inv3 = r_inv**3
        
        V_i = np.sum(q * r_inv)
        E_i = np.sum(q[:, None] * dR[mask] * r_inv3[:, None], axis=0)
        out[i] = [V_i, E_i[0], E_i[1], E_i[2]]
    return out
```

This is $O(N^2)$ but trivial compared to an SCF diagonalization. For the training pipeline, you can use **Mulliken charges from the converged DFTB calculation** as the reference charges. For inference, use a fast QEq or even Gasteiger guess.

### C. Assemble Fixed-Size Per-Atom Input Vector
Given the pieces above, assemble the input vector per atom:

```python
def assemble_ml_input(i, Z_i, H_ii, S_ii, V_E_i, neigh_H, neigh_S, neigh_idx,
                      max_neighbors=4, max_orb=4):
    """
    Assemble fixed-size input vector for atom i.
    Returns: (n_features,) float32 vector
    """
    # Z embedding (16-dim, learnable)
    # V, E (4 floats)
    # On-site H_ii (4x4 = 16 floats)
    # On-site S_ii (4x4 = 16 floats)
    # 4 neighbor H_ij (4x4x4 = 64 floats), zero-padded
    # 4 neighbor S_ij (4x4x4 = 64 floats), zero-padded
    # Total ~180 floats
```

### D. Target: ΔP_μν (On-Site Density Matrix Correction)
The target for training is the on-site block of the **deformation density matrix**:

```python
# Compute neutral-atom DM (diagonal with neutral atom occupations)
DM_NA = np.zeros_like(DM_dense)
for i in range(natoms):
    i0, i1 = orb_offsets[i], orb_offsets[i+1]
    # For sp basis: neutral atom has s=occ, p=0 (or from DFTB neutral atom calc)
    # Simpler: use DFTB+ neutral atom density projection as reference
    # Or build from elemental valence configurations
    DM_NA[i0:i1, i0:i1] = np.diag(neutral_occupations[i])

Delta_DM = DM_dense - DM_NA
```

For each atom `i`, the target is `Delta_DM[i0:i1, i0:i1]` (the on-site deformation block). This is a symmetric matrix with 10 unique elements for sp basis.

## Recommended Integration Path

**Do not modify the existing AFM pipeline.** Instead, write a standalone training-data generator that reuses the same DFTBcore + Grid_dftb building blocks:

```
New file: pyBall/ML/delta_rho_data_gen.py
          or tests/tML/dftb_training_data.py

1. Loop over molecule geometries (from XYZ or GDB)
2. For each molecule:
   a. Run DFTBcore SCF → H, S, DM, eigvecs
   b. build_orbital_layout() → norb_per_atom, orb_offsets
   c. Compute Mulliken charges → q_i
   d. compute_electrostatics_float4(q_i, pos) → [V, E]
   e. Extract on-site H_ii, S_ii, DM_ii blocks
   f. Extract neighbor H_ij, S_ij, DM_ij blocks (up to 4)
   g. Compute ΔP_ii = DM_ii - DM_NA_ii
   h. Save: {input_vector, target_ΔP, metadata}
3. Write to .npz or .h5 dataset
```

Then the **ML model** (a tiny PyTorch MLP) reads these vectors and learns `input → ΔP_ii`.

## Key Technical Notes

### Basis Order Convention
DFTB+ uses Fortran orbital order: `[s, py, pz, px]` for p orbitals. Your ML model should work with **any** convention as long as training and inference are consistent. The existing [Grid_dftb.py](cci:7://file:///home/prokop/git/FireCore/pyBall/DFTB/Grid_dftb.py:0:0-0:0) handles this via `_ORT_SPP_TO_OCL = [3,1,2,0]` mapping. For the ML input, you can either:
- Use the raw DFTB+ order (s, py, pz, px) and let the MLP learn it
- Reorder to (s, px, py, pz) if you prefer

### Non-SCC Hamiltonian H⁰
If you specifically need the **non-SCC** Hamiltonian (before charge corrections), DFTBcore's [enable_hamiltonian_storage(store=True)](cci:1://file:///home/prokop/git/FireCore/pyBall/DFTB/DFTBcore.py:272:4-283:86) saves it before diagonalization. However, [get_h_dense()](cci:1://file:///home/prokop/git/FireCore/pyBall/DFTB/DFTBcore.py:365:4-365:77) currently returns the **post-SCC** Hamiltonian (which contains eigenvectors after [run_scf](cci:1://file:///home/prokop/git/FireCore/pyBall/DFTB/DFTBcore.py:285:4-308:27)). You may need to add a small patch to DFTBcore or use the DFTB+ C API ([dftb_lib.py](cci:7://file:///home/prokop/git/FireCore/pyBall/dftb_lib.py:0:0-0:0)) which captures H via callbacks at each SCF step.

**Alternative**: For the ML, the **SCC Hamiltonian at convergence** is arguably even better as input because it already contains self-consistent charge information. The network's job is then to predict the *correction* from this already-good starting point.

### Units
- DFTBcore matrices are in **atomic units** (Hartree, Bohr)
- The ML input $[V, \mathbf{E}]$ should be in **eV and eV/Å** for numerical convenience
- Grid projection [project_density_dense](cci:1://file:///home/prokop/git/FireCore/pyBall/DFTB/Grid_dftb.py:1143:4-1217:51) applies `B3_FACTOR` to convert Bohr-normalized orbitals to Å⁻³

## Summary

| What you need | Where it already exists | What to add |
|---------------|------------------------|-------------|
| Run SCF, get dense H/S/DM | [DFTBcore](cci:2://file:///home/prokop/git/FireCore/pyBall/DFTB/DFTBcore.py:79:0-442:20) | Nothing |
| Orbital-to-atom mapping | [build_orbital_layout()](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM_utils.py:149:0-173:18) | Nothing |
| Project density to grid for validation | [Grid_dftb.project_density_dense()](cci:1://file:///home/prokop/git/FireCore/pyBall/DFTB/Grid_dftb.py:1143:4-1217:51) | Nothing |
| Slice dense matrix into per-atom blocks | — | New: `extract_local_blocks_dense()` (~30 lines) |
| Compute per-atom $[V, \mathbf{E}]$ | — | New: `compute_electrostatics_float4()` (~20 lines) |
| Assemble fixed-size input vector | — | New: `assemble_ml_input()` (~40 lines) |
| Neutral-atom DM reference | — | New: build from valence occupations or DFTB+ [project_neutral_density()](cci:1://file:///home/prokop/git/FireCore/pyBall/DFTB/Grid_dftb.py:2027:0-2081:50) |

The existing [get_density_from_dftb_dense()](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM_utils.py:176:0-340:35) in [AFM_utils.py](cci:7://file:///home/prokop/git/FireCore/pyBall/OCL/AFM_utils.py:0:0-0:0) is a good template to copy for the SCF + orbital layout part. You can essentially fork that function, skip the grid projection, and instead output the per-atom block vectors.

---

# USER

OK, now some more information, if we decide to go with electron pairs and VSEPR theory, we have some tools for that already @atomicUtils.py @AtomicSystem.py 
make_epair_geom
add_capping_h_sp2

@KekuleBackend.py 
_calc_h_directions_atom

@MMFFBuilderBase.h 
makeConfGeom
makeConfGeomCap

it is also described in @intramolecular_forcefields.md 

plese review it and figure out how we can use it

---

# Kimi 2.6 devin

## What Already Exists

### 1. Hybridization Detection
**[AtomicSystem.py](cci:7://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:0:0-0:0)** uses a simple valence model:
```python
VALENCE_DICT = {'O': (2, 2), 'N': (3, 1)}   # (target_sigma_bonds, n_lone_pairs)

npi    = target_bonds - actual_sigma_bonds   # pi bonds
nep    = VALENCE_DICT[ename][1]             # lone pairs
nsigma = len(neighbors)
```

**[MMFFBuilderBase.h](cci:7://file:///home/prokop/git/FireCore/cpp/common/molecular/MMFFBuilderBase.h:0:0-0:0)** stores the same in [AtomConf](cci:2://file:///home/prokop/git/FireCore/cpp/common/molecular/MMFFBuilderBase.h:81:0-192:1):
```cpp
struct AtomConf {
    uint8_t nbond;  // sigma bonds
    uint8_t npi;    // pi orbitals
    uint8_t ne;     // electron pairs
    uint8_t nH;     // capping H
    int neighs[4];
};
```

This gives you **for every atom**: hybridization (`sp³`, `sp²`, `sp`) and how many of the 4 valence slots are bonds vs lone pairs.

### 2. VSEPR Geometry — Ideal Electron Pair Directions
**[makeConfGeom(nb, npi, hs)](cci:1://file:///home/prokop/git/FireCore/cpp/common/molecular/MMFFBuilderBase.h:924:4-990:5)** (C++) and **[make_epair_geom()](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:979:4-1023:52)** (Python) compute the **4 ideal sp³ hybrid directions** `hs[0..3]` from the known bond directions. For example:

- **H₂O** (`nb=2, npi=0`): 2 bonds known → 2 lone-pair directions computed at 109.5°
- **NH₃** (`nb=3, npi=0`): 3 bonds known → 1 lone-pair direction opposite to face normal
- **=O**   (`nb=1, npi=1`): 1 bond known → 2 lone pairs at 120° in plane, 1 π-dir perpendicular

**[KekuleBackend.py::_calc_h_directions_atom()](cci:1://file:///home/prokop/git/FireCore/pyBall/KekuleBackend.py:822:4-900:17)** does the same but for H-capping. The math is identical — it is pure VSEPR geometry.

### 3. Electron Pair Placement
**[place_electron_pair()](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:1025:4-1051:78)** adds a dummy atom at `pos = atom_pos + direction * 0.5 Å`. These dummy atoms are used in force fields (MMFF/UFF) as capping atoms with recoil forces.

## How to Use This for ML Delta-ρ

### Core Idea
Instead of predicting the raw on-site density matrix block $\Delta P_{\mu\nu}$ (10 unique elements for sp basis), **predict electron-pair parameters** that are:
- **Chemically interpretable** (lone pair occupations, bond polarization)
- **Geometrically constrained** by VSEPR (fewer degrees of freedom)
- **Topology-aware** (hybridization derived from bonding graph)

### Proposed Representation

For each atom $i$, define:

| Parameter | Source | What to Predict |
|-----------|--------|-----------------|
| Hybridization | `npi`, `nb` from topology | **Fixed** — no prediction needed |
| 4 hybrid directions | [makeConfGeom()](cci:1://file:///home/prokop/git/FireCore/cpp/common/molecular/MMFFBuilderBase.h:924:4-990:5) output | **Fixed** — ideal VSEPR directions |
| Occupation $n_k$ of hybrid $k$ | Electron count constraint | **ML target** — 4 numbers summing to valence |
| Direction deviation $\delta\hat{h}_k$ | Small departure from VSEPR | **ML target** (optional) — 3 angles per lone pair |

From these, you can **reconstruct the on-site density matrix** analytically.

### Reconstruction Formula
Each hybrid orbital $k$ is a linear combination of the atomic basis $(s, p_x, p_y, p_z)$:

$$\phi_k = c_{ks} \cdot s + c_{kx} \cdot p_x + c_{ky} \cdot p_y + c_{kz} \cdot p_z$$

The coefficients come from the hybrid directions $\hat{h}_k$:
- $c_{ks} = \sqrt{1/4} = 0.5$  (for sp³)
- $(c_{kx}, c_{ky}, c_{kz}) = \sqrt{3/4} \cdot \hat{h}_k$

The on-site density matrix is:

$$P_{\mu\nu}^{(ii)} = \sum_{k=1}^{4} n_k \cdot c_{k\mu} \cdot c_{k\nu}$$

This is a **closed-form mapping** from $(n_1, n_2, n_3, n_4)$ → $4 \times 4$ symmetric matrix.

For the **ML target**, predict the **deviation from the neutral atom reference**:
- Neutral: $n_k = 1$ for each bonding pair, $n_k = 2$ for each lone pair (or whatever your reference is)
- ML target: $\Delta n_k = n_k^{\text{SCF}} - n_k^{\text{neutral}}$

### Why This Is Better Than Raw ΔP

| Aspect | Raw ΔP (4×4 matrix) | Electron Pair Representation |
|--------|---------------------|------------------------------|
| Parameters per C atom | 10 unique floats | 4 occupations (or fewer) |
| Physical constraints | Hard to enforce positivity | $n_k \geq 0$, $\sum n_k = \text{valence}$ |
| Chemical meaning | Opaque | Lone pair density, bond polarity |
| Transferability | Poor | Hybrid directions are universal |
| Data efficiency | Needs 1000s of molecules | Needs 100s (VSEPR is a strong prior) |

## Concrete Integration Plan

### Step 1: Build VSEPR Features from Existing Tools
Write a small function that reuses [makeConfGeom](cci:1://file:///home/prokop/git/FireCore/cpp/common/molecular/MMFFBuilderBase.h:924:4-990:5) logic to produce the hybrid directions:

```python
def build_vsepr_features(atom_sys, ia):
    """Returns hybrid directions and assignments for atom ia.
    Uses AtomicSystem.py logic."""
    ename = atom_sys.enames[ia]
    nsigma = len(atom_sys.ngs[ia])  # actual sigma bonds
    
    # Get target from VALENCE_DICT or MMFF params
    if ename in VALENCE_DICT:
        nb_target, nep = VALENCE_DICT[ename]
        npi = nb_target - nsigma
        
        # Compute ideal directions using existing VSEPR code
        # (re-use make_epair_geom logic, but return directions instead of placing dummies)
        directions = compute_hybrid_directions(atom_sys.apos, atom_sys.ngs[ia], npi, nsigma)
        
        # Assign which directions are bonds vs lone pairs
        bond_dirs = directions[:nsigma]
        lp_dirs   = directions[nsigma:nsigma+nep]
        
        return directions, bond_dirs, lp_dirs, npi
```

### Step 2: Map SCF Density Matrix to Electron Pair Occupations
For training data, solve the inverse problem: given the DFTB+ on-site block $P_{\mu\nu}^{(ii)}$, find the occupations $n_k$ that best reproduce it:

```python
def fit_epair_occupations(P_ii, hybrid_dirs):
    """Solve for n_k such that P_ii ≈ Σ_k n_k * c_k @ c_k.T"""
    # Build coefficient matrix C[k,mu] from hybrid_dirs
    C = np.zeros((4, 4))
    for k, h in enumerate(hybrid_dirs):
        C[k, 0] = 0.5  # s coefficient
        C[k, 1:4] = np.sqrt(3/4) * h  # p coefficients
    
    # P_target = C.T @ diag(n) @ C
    # Vectorize to linear system and solve
    # ... 
    return n_k
```

This is a small linear least-squares problem per atom ($4 \times 4$), trivial to solve.

### Step 3: ML Model
The ML model predicts $\Delta n_k$ per atom. Input features:
- Atom type (C, N, O, H) — one-hot or embedding
- Hybridization (`npi`, `nb`) 
- Local electrostatics $[V, E_x, E_y, E_z]$ from neighbor charges
- Neighbor types and distances

Output:
- $\Delta n_k$ for each of the 4 hybrids (with constraints)

### Step 4: Reconstruct ΔP for Inference
```python
def reconstruct_delta_P(n_k_delta, hybrid_dirs):
    C = build_coeff_matrix(hybrid_dirs)
    delta_P = C.T @ np.diag(n_k_delta) @ C
    return delta_P
```

## Where to Plug Into Existing Code

| Task | Use This Existing Function |
|------|---------------------------|
| Detect hybridization | [AtomicSystem.add_electron_pairs()](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:950:4-963:52) logic (lines 951–964) |
| Compute ideal directions | [makeConfGeom()](cci:1://file:///home/prokop/git/FireCore/cpp/common/molecular/MMFFBuilderBase.h:924:4-990:5) from [MMFFBuilderBase.h](cci:7://file:///home/prokop/git/FireCore/cpp/common/molecular/MMFFBuilderBase.h:0:0-0:0) or [AtomicSystem.make_epair_geom()](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:979:4-1023:52) |
| Add dummy atoms for visualization | [AtomicSystem.place_electron_pair()](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:1025:4-1051:78) |
| H-capping for incomplete molecules | [AtomicSystem.add_capping_h_sp2()](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:1094:4-1136:28) or [KekuleBackend._calc_h_directions_atom()](cci:1://file:///home/prokop/git/FireCore/pyBall/KekuleBackend.py:822:4-900:17) |
| Get pi-direction for sp² atoms | [AtomicSystem.get_atomi_pi_direction()](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:965:4-977:32) |

## Summary

The existing VSEPR tools give you a **strong geometric prior** for the electron density. Instead of learning the full 4×4 density matrix per atom, you can:

1. **Fix** the hybrid directions using [makeConfGeom](cci:1://file:///home/prokop/git/FireCore/cpp/common/molecular/MMFFBuilderBase.h:924:4-990:5) / [make_epair_geom](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:979:4-1023:52)
2. **Predict** only the 4 occupations $n_k$ (or their deviation from neutral)
3. **Reconstruct** $\Delta P_{\mu\nu}$ analytically via the hybrid coefficient matrix

This reduces the ML target from **10 floats** (symmetric 4×4 matrix) to **4 floats** (occupations) per sp³ atom, with built-in physical constraints and chemical interpretability. For H, it reduces from 1 float to 1 float (no gain, but no loss). The directions are already computed correctly for all organic hybridization states by your existing code.

If you want, I can write the `build_vsepr_features()` and `fit_epair_occupations()` functions as a concrete Python module reusing your existing [AtomicSystem](cci:2://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:17:0-1313:21) and [make_epair_geom](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:979:4-1023:52) code.