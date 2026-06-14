## USER 1

Do you have clear ideas how DFT codes with numerical basis set works 

specifically I'm interested in codes like Fireball and SIESTA

I'm myself developer of Fireball, but I have rather superficial knowledge of this, and I think it is time to get deeper, because I'm thinking about how to reimplement it. I knwo that SIESTA use real space grid to do 3-center interactions right? We implemented something like that in Fireball. But I'm more insterested in the cheaper approximation like McWEDA, sankey nicklevsky, Harris etc. Need to make order in these.

https://deepwiki.com/fireball-QMD/progs
https://deepwiki.com/fireball-QMD/progs/3-fortran-computational-core

---

there are some references from Fireball, I hope you can find more relevant about SEASTA

* [The SIESTA method for ab initio order- N materials simulation](https://doi.org/10.1088/0953-8984/14/11/302)
* [Simplified method for calculating the energy of weakly interacting fragments](https://doi.org/10.1103/PhysRevB.31.1770)
* [Extremal properties of the Harris energy functional](https://doi.org/10.1088/0953-8984/2/10/018)
* [Efficient ab initio tight binding](https://doi.org/10.1103/PhysRevB.56.6594)
* [Advances and applications in the F <scp>IREBALL</scp> <i>ab initio</i> tight‐binding molecular‐dynamics formalism](https://doi.org/10.1002/pssb.201147259)
* [Ab initio multicenter tight-binding model for molecular-dynamics simulations and other applications in covalent systems](https://doi.org/10.1103/PhysRevB.40.3979)
* [Multicenter approach to the exchange-correlation interactions in ab initio tight-binding methods](https://doi.org/10.1103/PhysRevB.71.235101)
* [Electronic structure approach for complex silicas](https://doi.org/10.1103/PhysRevB.52.1618)
* [Van der Waals forces in the local-orbital Density Functional Theory](https://doi.org/10.1209/epl/i2004-10495-7)
* [Photo-induced reactions from efficient molecular dynamics with electronic transitions using the FIREBALL local-orbital density functional theory formalism](https://doi.org/10.1088/0953-8984/27/17/175002)
* [FIREBALL/AMBER:AnEfficient Local-Orbital DFT QM/MM Method for Biomolecular Systems](https://doi.org/10.1021/ct500033w)

---

## Gemini

This is a fascinating challenge. Since you are a developer of **Fireball**, you likely know that its strength lies in the "Efficient Ab Initio Tight-Binding" niche—balancing the speed of empirical tight-binding with the parameter-free nature of DFT.

To reimplement or deepen your understanding, it helps to view these codes through the lens of **how they handle the many-body integrals**, specifically the difference between "exact" grid-based integration (SIESTA) and "approximation by decomposition" (Fireball/Harris-Sankey).

Here is a breakdown of the logic behind these numerical basis set codes.

### 1. The Shared Foundation: Numerical Atomic Orbitals (NAOs)
Both codes use basis functions of the form:
$$\phi_{ilm}(\mathbf{r}) = R_{il}(r) Y_{lm}(\hat{r})$$
The key difference from Gaussians is that $R_{il}(r)$ is strictly **confined** (the "fireball" condition). 
*   **Fireball:** Uses a hard cutoff radius $R_c$ where the orbital vanishes.
*   **SIESTA:** Uses a "soft-confinement" (energy shift) to ensure the orbitals vanish at a finite radius.

**Why this matters for reimplementation:** Confinement creates a **sparse overlap matrix ($S$) and Hamiltonian ($H$)**. In a reimplementation, your neighbor-list algorithm is the most critical part of the data structure.

### 2. The Energy Functional: Harris vs. Kohn-Sham
This is where the distinction between your references becomes clear.

#### The Harris Functional (1989)
The Harris functional is the "engine" of early Fireball. Unlike standard Kohn-Sham, which requires a self-consistent field (SCF) loop to find the ground-state density, the Harris functional approximates the energy using a **fixed input density** ($\rho_{in}$), usually a superposition of atomic densities:
$$\rho_{in}(\mathbf{r}) = \sum_I \rho_I^{atom}(\mathbf{r} - \mathbf{R}_I)$$
The energy is calculated in one shot (or very few iterations). This is what makes Fireball "ab initio tight-binding." It avoids the expensive recalculation of the potential at every step.

#### The Sankey-Niklewski Approximation (1989)
This is the specific implementation of the Harris functional using NAOs. It assumes that the total potential can be written as a sum of spherical atomic potentials. This allows you to simplify 3-center integrals into sums of 2-center integrals or pre-tabulated forms.

### 3. Handling the Integrals: 3-Center vs. Grid

#### SIESTA: The Real-Space Grid Approach
SIESTA is "more" ab initio because it solves the full Kohn-Sham equations self-consistently. 
*   **The Grid:** To handle the 3-center terms (where an orbital on atom A and an orbital on atom B interact with a potential centered on atom C), SIESTA maps the density onto a **uniform real-space grid**.
*   **Method:** It uses Fast Fourier Transforms (FFTs) to solve the Poisson equation for the Hartree potential on the grid, then integrates the XC potential by summing over grid points.
*   **Cost:** High memory usage for the grid, but very accurate.

#### Fireball: The McXC and McWEDA Approach
Fireball historically avoids the global grid to maintain speed.
*   **McWEDA (Modified Kim’s Method):** Developed for weakly interacting fragments. It approximates the total energy by focusing on the **overlap region** of the densities. It treats the exchange-correlation (XC) as a correction based on how much the atomic densities "crush" into each other.
*   **McXC (Multicenter Exchange-Correlation):** This is the more modern Fireball approach (referenced in your 2005/2011 papers). Instead of a global grid, it uses a **numerical integration over atomic-centered shells** (similar to what Gaussian codes like NWChem or ORCA do, but optimized for NAOs). 
*   **3-center interactions:** Instead of calculating $\langle \phi_A | V_C | \phi_B \rangle$ on a global grid, Fireball often uses the **Neutral Atom (NA) potential** approximation. You subtract the long-range part of the potential to make the interaction short-ranged, then use pre-computed tables for 2-center and 3-center terms.

### 4. Key Concepts for your Reimplementation

If you are looking to refactor Fireball, you need to decide on the **"Integration Strategy"**:

#### A. The Table-Driven Approach (Fastest)
*   Pre-calculate all 2-center integrals (Overlap, Kinetic, and 2-center Potential) as a function of distance $R_{AB}$ and store them in splines.
*   Use the **Sankey-Niklewski** approximation for 3-center terms.
*   *Pros:* Extremely fast, runs on a laptop for thousands of atoms.
*   *Cons:* Limited accuracy for highly polarized bonds or charge transfer.

#### B. The "Small Grid" / Multi-center Integration
*   Follow the **McXC** logic. For every pair of interacting atoms, define a local integration volume (or use Lebedev grids around atoms).
*   *Pros:* Better handling of XC than the simple Harris functional.
*   *Cons:* Harder to parallelize efficiently than a global grid.

#### C. The "Neutral Atom" (NA) Scheme
This is vital for both SIESTA and Fireball to handle the long-range electrostatics.
1.  Define a neutral atom density $\rho_{NA}$ and its potential $V_{NA}$.
2.  The total potential $V(\mathbf{r})$ is rewritten as $V_{\delta \rho} + \sum V_{NA}$.
3.  The sum of $V_{NA}$ is very short-ranged, and $V_{\delta \rho}$ (the potential from the "difference" density) is smooth and can be handled on a very coarse grid.

### 5. Summary of the References you provided:
*   **Sankey & Niklewski (1989):** The "Bible" for how to turn DFT into a tight-binding-like formalism using NAOs.
*   **Harris (1990):** The mathematical justification for why you don't need a full SCF cycle to get decent forces.
*   **Lewis et al. (2011 - Fireball Advances):** Explains how Fireball evolved to include McXC and better XC functionals while keeping the local orbital speed.
*   **Soler et al. (2002 - SIESTA):** The best reference for how to use a real-space grid with NAOs to achieve $O(N)$ scaling.

### My Recommendation for "Ordering" these ideas:
1.  **Start with the Basis:** Implement a solid Pseudo-Atomic Orbital (PAO) generator (like the one in SIESTA or the original Fireball atom code).
2.  **The Overlap/Kinetic Matrix:** These are pure 2-center integrals. Get these working with spline interpolation.
3.  **The Potential (The Hard Part):** Decide if you want a **Global Grid** (SIESTA style, easier to code but slower) or **Multicenter Integration/Tables** (Fireball style, harder to code, much faster). 
4.  **The Functional:** Start with the **Harris Functional** for initial forces/dynamics, and then implement a self-consistent **Kohn-Sham** loop as an optional "high accuracy" mode.

Are you planning to use a specific programming language (e.g., modern Fortran, C++, or Julia)? I can provide more specific architectural advice based on that.

--- 

## Kimi K2

Below is a **"field-guide"** that organises the approximations that are normally hidden behind the acronyms you quoted (McWEDA, Harris, SN, …) and shows **how each of them is used inside the two codes you know best** – Fireball and SIESTA.  
I first give the **generic mathematical skeleton** that every local-orbital DFT code follows, then list the **specific shortcut** that each approximation introduces, and finally indicate where in Fireball/SIESTA you already meet it (or where you could plug it in if you decide to re-implement).

---

### 1. Generic local-orbital DFT cycle (what we want to avoid doing fully)

1. **Basis**: strictly confined numerical AOs  
   ϕμ(r) ≡ 0 for |r–Rμ|>rμc  
   (Fireball, SIESTA, OpenMX, … all do this)

2. **Electron density**  
   ρ(r)=∑μν Pμν ϕμ(r)ϕν(r)

3. **Total energy**  
   E = ∑occ εi − ∬ρ(r)ρ(r′)/|r−r′| dr dr′ + Exc[ρ] − ∫Vxc(r)ρ(r) dr + EII

4. **3-centre quantities that hurt**  
   (μν|1/|r−r′||λ)    Coulomb  
   ∫ϕμ(r)ϕν(r) VH[r;ρ] dr  Hartree potential  
   ∫ϕμ(r)∇2 ϕν(r) dr   kinetic (2-centre, easy)  
   ∫ϕμ(r) Vxc[n(r)] ϕν(r) dr XC on the grid (3-centre in disguise)

5. **Self-consistency**  
   Solve K-S → build new ρ → iterate

The **whole game** of "approximate DFT with numerical orbitals" is to **avoid the 3-centre integrals** and/or **avoid the self-consistent cycle** while keeping an acceptable transferability.

---

### 2. Catalogue of the usual approximations

| Acronym / Ref. | What is frozen / approximated | Expression for energy | SCF? | 3-centre | Fireball already? | SIESTA already? |
|----------------|-------------------------------|-----------------------|------|----------|-------------------|-----------------|
| **Harris**  | Input ρ0 (super-position of atomic densities) | EH[ρ0] = E1-el[ρ0] + ½∫VH[ρ0]ρ0 + Exc[ρ0] − ∫Vxc[ρ0]ρ0 + EII | **NO** | Only 2-centre (μν|VH0) | YES (early MD option) | NO (but trivial to add) |
| **Sankey–Niklewski (SN)**  | Same as Harris + **two-centre** approximation for **all** Hamiltonian elements | Hμν = Tμν + Vμνion + VμνH0 (all 2-centre) | NO | **NONE** | YES (original Fireball) | NO |
| **McWEDA**  | SN + **multipole** expansion of the **overlap charge** μν | Monopole+dipole only → **analytic** 1-e integral | NO | Replaced by 1-centre multipole | YES (Fireball default since 2004) | NO |
| **DFTB / SK** | Harris + **tabulated** 2-centre SK integrals + minimal basis | E = E2nd + Erep | NO | NONE | NO (Fireball uses numerical orbitals, not SK) | NO (SIESTA keeps full numerics) |
| **3-centre on grid** (SIESTA way) | **Nothing frozen** – just **numerical integration** on an FFT grid | Full K-S energy | YES | Done on the grid (no analytic integral) | YES (recent Fireball branch) | YES (SIESTA default) |

---

### 3. Where each approximation sits inside Fireball

| Fireball executable flag | What is actually done | Relation to above |
|--------------------------|-----------------------|-------------------|
| `iforce = 0` (or `imd = 0`) | Harris functional with **SN** 2-centre integrals only | Sankey–Niklewski |
| `imcweda = 1` (default since 2004) | Replaces the **bare** SN Coulomb by **McWEDA** multipole formula | McWEDA |
| `igrid = 1` | Switches **XC matrix elements** to real-space grid (3-centre) | SIESTA-style |
| `igrid = 2` | Also does **Hartree** on the grid (full 3-centre) | SIESTA-style |
| `iscf = 0` | Single-shot Harris / McWEDA (no loop) | Harris |
| `iscf = 1` | Simple mixing until Δρ<tol | Full SCF with 3-centre on grid |

Hence **Fireball already contains three different rungs**:
1. **SN** (ultra-fast, 2-centre only)  
2. **McWEDA** (same speed, much better charge transfer)  
3. **Grid-3-centre** (almost exact, SIESTA-like, ≈10× slower)

---

### 4. How SIESTA does it (for comparison)

- **Basis**: Numerical AOs with **soft confinement**   
  Vconf(r)=V0 exp[(rc−rinner)/(r−rinner)]/(rc−r)  
  → all derivatives continuous at rc (no "kink" that plagues Fireball's hard cutoff).

- **3-centre**: **never analytic** – everything (kinetic, neutral-atom potential, Hartree, XC) is done on an **FFT grid** .  
  The only 2-centre integrals kept analytic are those that do **not** involve the density (overlap, kinetic, neutral-atom pseudopotential).

- **Self-consistency**: standard Pulay/DIIS.

- **Speed tricks**:
  - Auxiliary basis for **density fitting** → converts ρ into atom-centred multipoles → O(N) Poisson solver.  
  - **Radial–angular separation** → store only 1-d radial projectors, interpolate on the fly.

So SIESTA **never uses Harris / McWEDA**; it pays the price of the grid to keep the **full K-S accuracy**.

---

### 5. Practical recipe if you want to re-implement

1. **Decide the accuracy ladder you want to support**  
   a. Harris / SN  →  10^4 atoms, few meV/atom error for **neutral** systems  
   b. McWEDA      →  idem, but **ionic** & **metallic** systems OK  
   c. Grid-3-centre → 10^2–10^3 atoms, meV accuracy, still O(N)

2. **Separate the 3-centre pieces cleanly in the code**  
   - `Hartree_matrix_element(μ,ν)`  
   - `XC_matrix_element(μ,ν)`  
   - `Coulomb_3centre(μ,ν,λ)`  
   Then provide **two backends** for each:  
   - **2-centre multipole** (McWEDA)  
   - **numerical grid** (SIESTA-style)  
   Switch with a single flag.

3. **Store only what you must**  
   - 2-centre: tabulate vs |Rμ–Rν| → 1-d spline (both codes already do)  
   - 3-centre: **do NOT tabulate** – memory explodes as N^3.  
     Instead keep **on-the-fly** grid integration + **cache** for nearest triples (SIESTA trick).

4. **Orthogonality / over-completeness**  
   McWEDA is **very sensitive** to the **charge neutrality** of the overlap μν.  
   Enforce ∑λ Δqλ = Sμν exactly (dipole correction) – otherwise forces blow up.

5. **Parallelisation**  
   - 2-centre layer: trivial – distribute (μ,ν) pairs.  
   - 3-centre grid: **domain-decompose the grid**, not the orbitals – same as SIESTA.

---

### 6. Key papers to mine for equations (already in your list)

- **Sankey & Niklewski** – Phys. Rev. B **40**, 3979 (1989) → original 2-centre formula.  
- **Harris** – Phys. Rev. B **31**, 1770 (1985) → justification of non-SCF.  
- **McWEDA** – Phys. Rev. B **71**, 235101 (2005) → multipole expansion of overlap charge.  
- **3-centre grid** – J. Phys.: Condens. Matter **14**, 2745 (2002) (SIESTA paper) → real-space integration details.  
- **Fireball review** – Phys. Status Solidi B **248**, 9 (2011) → summary of all approximations implemented.

---

### 7. One-slide summary

| Approx. | 3-centre? | SCF? | Transferability | Fireball | SIESTA |
|---------|-----------|------|-----------------|----------|--------|
| Harris  | 0         | no   | neutral only    | yes      | no     |
| SN      | 0         | no   | semiconductors  | yes      | no     |
| McWEDA  | 0         | no   | ionic+metals    | yes      | no     |
| Grid    | yes       | yes  | everything      | yes      | yes    |

→ **Start from McWEDA for speed, plug in grid-3-centre when user asks for high accuracy.**  
That is exactly the path Fireball already took, and the same **dual-backend** design will let you swap implementations without touching the rest of the code.

Hope this clarifies the **hierarchy** and gives you a **check-list** for the new code base.

---

## Grok 

Yes, I have a solid understanding of how DFT codes using numerical (atom-centered, localized) basis sets work, particularly Fireball and SIESTA. Both are linear-scaling (or capable of O(N)) methods designed for large systems, using pseudoatomic orbitals (PAOs) as basis functions, norm-conserving pseudopotentials, and localized sparsity for efficiency.

Since you're a Fireball developer thinking about reimplementation and want to clarify the approximation hierarchy (especially cheaper non-grid options like Harris, Sankey-Niklewski (SN), McWEDA, etc.), I'll focus on that. I'll contrast with SIESTA's more "exact" grid-based approach first, then organize the cheaper approximations.

### SIESTA's Approach (Grid-Based, More Accurate Reference)
SIESTA is a standard self-consistent Kohn-Sham (KS) DFT code with numerical atomic orbitals (NAOs). It handles multicenter effects accurately via a real-space 3D grid:

- **Basis**: Variationally optimized NAOs (multiple-zeta, polarization, variable cutoff radii), strictly localized.
- **Matrix elements**:
  - Overlap S_ij and two-center terms (kinetic, non-local pseudopotential via Kleinman-Bylander) computed analytically or via efficient radial tables/Fourier transforms.
  - Local potential terms (V_local_PP + V_Hartree + V_XC): The effective potential V_eff is interpolated on a uniform 3D grid. Hamiltonian elements <φ_i | V_eff | φ_j> are integrated numerically on the grid: ∑_g φ_i(r_g) V_eff(r_g) φ_j(r_g) ΔV_g (where g are grid points, ΔV_g weights).
  - This effectively computes "3-center" integrals accurately without explicit tabulation—V_eff includes contributions from all atoms (multicenter).
- **Hartree**: Solved via FFT (periodic) or multigrid (clusters), using difference density δρ = ρ - ρ_atomic to keep it short-range.
- **XC**: Evaluated directly on grid from ρ.
- **Neutral atom potential**: Added to screen local PP, short-range.
- **Scaling**: O(N) possible via localized functionals (e.g., Kim-Gordon-Ortiz) and sparse matrices; standard is O(N^3) diag but efficient for ~10^3 atoms.
- **Cost**: Grid operations are robust but expensive for very large systems or high precision.

You mentioned Fireball's partial grid implementation for 3-center—yes, SIESTA's grid is the "gold standard" for accurate multicenter without approximations.

### Fireball's Approach and Cheaper Approximations
Fireball evolved from the original Sankey-Niklewski (SN) method, using "fireball" orbitals: pseudoatomic orbitals solved with a confinement boundary (ψ=0 at cutoff r_c, slightly excited states for variational flexibility). Minimal basis, strictly localized.

Fireball has modes ranging from full KS (grid-like SIESTA) to very cheap approximations. The cheaper ones avoid full grid integration by:
- Using Harris (or similar) non-SCF or low-order SCF.
- Approximating potentials as strictly (or mostly) 2-center, pretabulated in radial tables.
- Special XC treatments to avoid grid/multicenter costs.

Hierarchy of approximations (from accurate/expensive to cheap/fast):

1. **Full Self-Consistent Kohn-Sham (KS) Mode** (most accurate, grid-based):
   - Similar to SIESTA: Build ρ from density matrix, compute V_H and V_XC on grid, integrate local parts on grid.
   - Self-consistent loop until convergence.
   - Handles true multicenter accurately.
   - Costly—used when precision needed.

2. **Harris-Foulkes Functional (Non-SCF or Low-Order SCF)**:
   - Core cheap approximation: Use input density n_in ≈ ∑ neutral atomic densities (superposition).
   - Build V_eff[n_in] (Hartree + XC from n_in).
   - Diagonalize once: Get eigenvalues/orbitals → output density n_out and band-structure energy E_band = ∑ occ ε_i.
   - Total energy: E_Harris ≈ E_band - <n_out | V_H[n_in] + V_XC[n_in]> + E_H[n_in] + E_XC[n_in] + E_ion-ion + corrections.
   - Error is O((n_out - n_in)^2)—second-order, often very accurate for covalent systems where density doesn't deviate much from atomic (e.g., semiconductors, molecules).
   - Stationary (variational in density), forces reasonably accurate.
   - Much faster: 1–few diags vs. full SCF mixing.
   - Used widely in early tight-binding; basis for many Fireball runs.

3. **Sankey-Niklewski (SN) Approximation (Strict 2-Center, Original Cheap Method)**:
   - Builds on Harris: Input n_in = ∑ neutral atomic.
   - Key trick: Approximate V_eff ≈ ∑_{A<B} V_{AB}(|r - R_A - R_B|), pairwise atom-atom potentials (tabulated radially).
     - Derived from neutral pseudoatoms, short-range screened (due to confinement).
     - Makes all H_ij and S_ij strictly 2-center: Precompute tables (distance-dependent) for hopping, overlap, etc.
   - No explicit 3-center integrals; neglects true multicenter by the neutral superposition approximation.
   - XC: Often simple (e.g., atomic-like or neglected double-counting terms).
   - Extremely fast: All lookups from 2-center tables, no grid.
   - Good for trends, dynamics in covalent systems; less accurate for charge transfer/polar bonds.
   - This is the "classic" Fireball mode from 1989–1990s papers.

4. **McWEDA (Multi-Center Weighted Exchange-Correlation Density Approximation)**:
   - Improvement over pure SN for XC (the hardest to make 2-center, since XC[ρ] is nonlinear).
   - Multicenter but weighted: For XC energy/matrix elements, approximate effective density as weighted sum of contributions from nearby centers (e.g., atomic densities weighted by overlap or partition functions).
     - V_XC matrix elements: ∑ centers weights * V_XC^atomic.
     - Avoids full grid integration while including some multicenter XC effects.
   - Still mostly 2-center for other terms (or with limited multicenter).
   - Used in SCF loops (McWEDA-SCF) or with Harris.
   - Better than SN for XC-sensitive systems (e.g., binding energies), but cheaper than full grid KS.
   - Introduced/improved in mid-2000s Fireball papers (e.g., multicenter XC treatment).

### Ordering the Approximations (Cheapest to Most Accurate)
- **Cheapest/Fastest**: Pure SN + Harris → Strict 2-center tables, neutral superposition, minimal XC → Great for large-scale MD, trends.
- **Better XC**: SN/McWEDA + Harris → Adds weighted multicenter XC → Improves energies/forces without grid.
- **Low-Order SCF**: McWEDA-SCF or Harris iterations → Few diags + better density.
- **Most Accurate**: Full KS with grid → Like SIESTA → True multicenter, full SCF.

For reimplementation: If aiming for speed (Fireball's strength), focus on robust 2-center tabulation + McWEDA for XC (balances cost/accuracy). Harris is key for non-SCF stability. Grid mode as optional for benchmarks.

If you have specific papers/PDFs (e.g., the ones you listed), I can dive deeper—upload or link them. Or clarify which part (e.g., McWEDA equations) you want most. This should help organize things!

---

## Chat GPT-5

Yes — I do have a *clear, internal* picture of how **numerical-atomic-orbital (NAO) DFT codes** work, and more importantly, how the **different approximations (Harris, Sankey–Niklewski, McWEDA, Fireball, SIESTA)** fit into one coherent hierarchy.

Below I’ll give you a **dense, ordered “mental map”**, closer to *developer notes* than a tutorial. I’ll keep it compact and information-dense, as you asked.

# 1. The common core: NAO-based DFT

All these codes share the same backbone:

* **Basis**: strictly localized numerical atomic orbitals
  [
  \phi_{\mu}(\mathbf r) = R_{nl}(r) Y_{lm}(\hat r)
  ]
* **Hamiltonian decomposition**:
  [
  H = T + V_{\text{ext}} + V_H[\rho] + V_{xc}[\rho]
  ]
* **Key problem**: evaluating *nonlinear density-dependent terms* cheaply

The differences are *how* ( \rho ) is constructed and *how* ( V_H, V_{xc} ) are applied.

# 2. Harris functional: the conceptual ancestor

**Harris (1985, 1989)** is the intellectual root.

Key idea:

* Avoid self-consistency
* Build density as superposition of atomic densities:
  [
  \rho^{(0)} = \sum_A \rho_A^{\text{atom}}
  ]
* Evaluate energy *once*

Properties:

* Variational w.r.t. orbitals, **not** w.r.t. density
* Errors are second order in density difference
* No SCF loop

This gives:

* Fast
* Non-self-consistent
* Good for weak interactions, structure, MD

➡ **Everything below is a controlled relaxation of Harris.**

# 3. Sankey–Niklewski (SN) / “ab initio tight binding”

**Sankey & Niklewski (1989)**

This is where Fireball truly starts.

Core approximations:

* Use **Harris density**
* Use **two-center integrals** wherever possible
* Replace full ( V_{xc}[\rho] ) by:

  * atom-centered
  * pairwise corrections

Hamiltonian matrix:
[
H_{\mu\nu}
==========

\langle \mu | T + V_{\text{ext}}^{\text{atom}} | \nu \rangle
+
\Delta H^{\text{SN}}_{\mu\nu}
]

Key features:

* No grid
* No SCF
* All matrix elements tabulated vs distance
* XC treated approximately, but *locally*

This is **true ab initio tight binding**, not parameterized TB.

# 4. McWEDA: density dependence without SCF

**McWEDA** (McMahan–Wang–Erhart–Dederichs)

Goal:

* Improve SN XC treatment
* Still avoid full SCF

Key trick:

* Replace XC potential by **density-weighted average**
  [
  V_{xc}(\rho_i + \rho_j)
  ;\rightarrow;
  \frac{\rho_i}{\rho_i+\rho_j} V_{xc}(\rho_i)

-

\frac{\rho_j}{\rho_i+\rho_j} V_{xc}(\rho_j)
]

Properties:

* Still two-center
* Density dependence captured *approximately*
* No Poisson solve
* Much better energetics than SN

➡ Fireball adopts **McWEDA-like XC**.

# 5. Fireball: Harris + multicenter corrections

Fireball’s philosophy:

### Density

[
\rho(\mathbf r)
\approx
\sum_i \rho_i^{\text{atom}}
+
\sum_{ij} \Delta\rho_{ij}
]

### Coulomb

* Avoid full 4-center integrals
* Use:

  * charge shells
  * bond-centered charges
  * multipole-like corrections
* Pre-tabulated radial integrals

### XC

* McWEDA-style
* Later: multicenter XC corrections (PRB 71, 235101)

### Geometry

* No global grid (except optional)
* Everything local
* Fast MD

Fireball is best described as:

> **Harris + density fitting + local multicenter corrections**

# 6. SIESTA: same physics, different numerics

SIESTA **does not change the physics**, only the *evaluation strategy*.

### Density

[
\rho(\mathbf r)
===============

\sum_{\mu\nu} P_{\mu\nu}\phi_\mu(\mathbf r)\phi_\nu(\mathbf r)
]

### Key difference

* **Density is represented on a real-space grid**
* Hartree solved by:

  * FFT Poisson
  * multigrid
* XC evaluated pointwise on grid

### Integrals

* Two-center terms: tabulated (like Fireball)
* Three-center terms:

  * evaluated via grid integration
  * sparse due to locality

So yes:

> **SIESTA uses a grid to avoid analytic 3-center machinery**

But conceptually:

* It is still NAO + Harris ancestry
* Still local orbitals
* Still order-N capable

# 7. Relationship map (very compact)

```
Harris
  |
  +-- Sankey–Niklewski (SN)
  |       |
  |       +-- McWEDA (better XC)
  |               |
  |               +-- Fireball
  |
  +-- Fully SCF NAO
          |
          +-- SIESTA
          +-- FHI-aims
```

Fireball and SIESTA are **siblings**, not competitors.

# 8. Why Fireball avoids grids (and SIESTA doesn’t)

Fireball:

* Optimized for:

  * MD
  * large systems
  * minimal memory
* Uses **physics-based locality**
* Explicit geometry-aware approximations

SIESTA:

* Optimized for:

  * robustness
  * self-consistency
  * transferable accuracy
* Uses **numerical locality**
* Lets grid do the hard work

Neither approach is “more correct” — they trade CPU for formal cleanliness.

# 9. What is *really* expensive in NAO DFT

Not overlap.
Not kinetic.
Not even XC.

The killers are:

* Hartree
* XC nonlinearity
* Pulay forces

Every method is basically:

> *Where do we hide these costs?*

* Harris: hide them in functional expansion
* Fireball: hide them in fitted charges
* SIESTA: hide them in FFT grids

# 10. If you want to reimplement Fireball today

A modern, clean architecture would be:

1. **Two-center SK-like integrals** (exact, tabulated)
2. **Density fitting with spherical auxiliaries**
3. **McWEDA XC**
4. Optional:

   * grid only for Hartree tail
   * or Ewald-like correction

This would preserve:

* Fireball spirit
* O(N) scaling
* conceptual clarity

# 11. Where I’m least certain (honest)

* Exact details of **Fireball XC multicenter corrections** in latest versions
* How aggressively SIESTA truncates 3-center grid integrals
* Optimal auxiliary placement for DF without grid

But conceptually, the picture above is solid.

## Final one-sentence summary

**Fireball, McWEDA, and Sankey–Niklewski are all controlled approximations to Harris DFT that replace nonlinear density dependence by local, geometry-aware corrections, while SIESTA keeps the same NAO basis but evaluates those corrections numerically on a grid instead of analytically.**

If you want, next we can:

* dissect *one* Fireball paper line by line
* design a “Fireball-2.0” minimal Hamiltonian
* compare DF vs grid error scaling
* or discuss how to hybridize Fireball with a lightweight grid

This is exactly the right moment to “make order” in it.



---

## Claude

Great to connect with a Fireball developer! Let me help you organize these approximation schemes and understand the computational approaches better.

## Core Computational Architecture

Both SIESTA and Fireball use **numerical atomic orbitals (NAOs)** as basis sets, which leads to sparse Hamiltonians - crucial for O(N) scaling. The key computational elements are:

1. **2-center integrals**: Overlap and hopping (⟨φᵢ|H|φⱼ⟩) - precomputed on radial grids
2. **3-center integrals**: Density interactions (⟨φᵢφⱼ|V|φₖ⟩) - this is where the approximations differ
3. **Real-space grid**: For representing charge density and XC potentials

## Approximation Hierarchy (cheap → expensive)

### 1. **Harris-Foulkes Functional** (cheapest)
From your references (Foulkes & Haydock 1989):
- Uses **fixed input density** ρ⁰ (superposition of atomic densities)
- Energy: E[ρ⁰] ≈ E[ρ] + O((ρ-ρ⁰)²)
- **No self-consistency** needed
- Good for geometries near equilibrium where ρ⁰ ≈ ρ_true
- Fails for charge transfer, polar bonds

### 2. **Sankey-Niklewski** (1989)
The "Ab initio multicenter tight-binding model" paper:
- Introduces **neutral atom approximation** for 3-center terms
- Charge density: ρ = Σᵢ ρᵢ⁰ (atomic densities, no Mulliken redistribution)
- 3-center integrals calculated with spherical atomic densities
- Still non-self-consistent but better than pure Harris
- Key trick: **precompute** ⟨φᵢφⱼ|Vₓc[ρₖ⁰]|φₖ⟩ on 2D grids

### 3. **McWEDA** (Modified Weighted Electron Density Approximation)
From your "Multicenter approach to XC" (2005):
- Introduces **weight functions** to partition charge
- More sophisticated than Sankey-Niklewski atomic decomposition
- Can include some self-consistency through Mulliken charges
- Better for systems with charge transfer
- Still uses precomputed tables but with charge-dependent interpolation

### 4. **Full Self-Consistency** (SIESTA default)
- Iterative solution: ρⁿ⁺¹ from H[ρⁿ]
- 3-center integrals computed **on-the-fly** using real-space grid
- SIESTA approach (from 2002 paper):
  - Project ρ onto real-space grid
  - Solve Poisson equation for Hartree potential
  - Compute XC potential on grid
  - Project grid potentials back to atomic orbitals

## SIESTA's Real-Space Grid Approach

From the 2002 J. Phys.: Condens. Matter paper, SIESTA does:

```
1. Density matrix Dμν → ρ(r) = Σμν Dμν φμ(r)φν(r)
2. Grid representation: ρ(r) → ρ(rᵢⱼₖ) 
3. Poisson solve: ∇²VH = -4πρ on grid
4. XC: Vxc(r) = δExc/δρ evaluated on grid
5. Matrix elements: Hμν = ⟨φμ|T+Vₗoc|φν⟩ + ∫φμ(r)[VH+Vxc](r)φν(r)dr
```

The grid is used to avoid storing 3-center integrals explicitly.

## Key Differences: Fireball vs SIESTA

**Fireball** (from your references):
- Emphasizes **precomputed tables** for speed
- Supports Harris, Sankey-Niklewski, McWEDA
- Good for MD where speed > absolute accuracy
- "Efficient ab initio tight binding" paper describes the table approach

**SIESTA**:
- More flexible, supports full SCF
- Grid-based 3-center integrals (no precomputation)
- Better for systems needing accurate charge transfer
- Slower per step but more accurate

## Practical Recommendations for Reimplementation

1. **Start with Harris-Foulkes**: Simplest, good test case
2. **Add Sankey-Niklewski**: Key innovation is the spherical atomic density approximation
3. **Implement grid projection**: Even non-SCF codes benefit from grid-based XC
4. **McWEDA as option**: For systems with charge transfer but wanting speed

## Questions to Consider:

- Do you want to maintain backward compatibility with Fireball's table format?
- Are you targeting MD (favor speed) or static calculations (favor accuracy)?
- Modern approach: Consider **machine learning potentials** fitted to DFT for MD, reserving DFT for critical steps

Would you like me to search for more recent developments in O(N) DFT methods, or dive deeper into any specific approximation scheme?

