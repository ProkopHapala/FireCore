# Non-Bonding Force Fields

This document covers pairwise and many-body non-covalent interactions: van der Waals, electrostatics, and short-range repulsion. FireCore implements these across several abstraction levels, from brute-force O(N²) pairwise loops to approximate O(N log N) tile-based Fast Multipole Methods.

**Related Windsurf Codemaps:**
- [FireCore Classical Forcefields: MMFFsp3 & UFF (CPU/GPU/Python)](https://windsurf.com/codemaps/53f2fe2c-ac5c-4c0b-b905-af6653adde97-8796fe608a7d71c1) — NBFF base class and LJ/Morse/Coulomb evaluation architecture.
- [MMFF/UFF CPU vs GPU Testing](https://windsurf.com/codemaps/8d1b056f-1502-4363-b52d-8257de4be453-8796fe608a7d71c1) — CPU vs GPU parity for non-bonded terms (getLJQH, getMorsePLQH).
- [FitREQ_PN: Hydrogen-Bond Parameter Fitting System](https://windsurf.com/codemaps/d977d597-94b4-42c3-a92a-0cefe34a3e82-8796fe608a7d71c1) — H-bond parameter fitting and REQ optimization.
- [DFTB Reference Calculation & FDBM AFM Forcefield Comparison System](https://windsurf.com/codemaps/1153fe89-ff29-4d4b-b4a6-e97d8f37047f-fe86ab10a43f3d18) — Comparing DFTB reference to classical non-bonded potentials.
- [AFM FDBM Pipeline: DFTB Backend & pySCF Integration Points](https://windsurf.com/codemaps/02d559c9-de47-4058-b07b-3318664b454e-fe86ab10a43f3d18) — DFTB-derived force-field parameter pipeline.

---

## 1. NBFF (Non-Bonded Force Field)

### Physics & Purpose

`NBFF` is the base class for all non-bonded interactions in FireCore. It evaluates the pairwise energy between atoms using a combination of:

1. **Lennard-Jones (12-6)** — dispersion and Pauli repulsion.
2. **Morse potential** — chemically-motivated bond dissociation with exponential repulsion and attraction.
3. **Coulomb electrostatics** — damped to avoid the r=0 singularity.
4. **Hydrogen-bond correction** — pseudo-charge interactions to enhance directional H-bonding.

The total non-bonded energy for a pair (i,j) is:

$$E_{ij} = E_{\text{LJ}}(r_{ij}) + E_{\text{Morse}}(r_{ij}) + E_{\text{Coulomb}}(r_{ij}) + E_{\text{HB}}(r_{ij})$$

### Implementation Files

- **`cpp/common/molecular/NBFF.h`** — Main `NBFF` class.
  - Member arrays: `REQs[natoms]` (R₀, ε₀, Q, H), `neighs[natoms]` (bonded neighbors), `excl[natoms*EXCL_MAX]` (exclusion list for 1-2, 1-3 neighbors).
  - `evalPLQs()` / `makePLQs()` — converts REQs to PLQ (Pauli, London, Charge) factorized form for grid evaluation.
  - `evalPointCoulPBC()` — brute-force Coulomb summation with periodic images.
  - `fitAABB()` / `pointBBs` — bounding-box collision detection for short-range acceleration.
- **`cpp/common/molecular/NBFF_SR.h`** — Short-range variant. Replaces LJ/Coulomb with a purely repulsive R⁻⁴ potential (`repulsion_R4()` in `Forces.h`) plus optional H-bond terms. Used for fast steric clash detection in rigid-body assembly.
- **`cpp/common/molecular/NBFF_old.h`** — Legacy implementation, superseded by `NBFF.h`.
- **`cpp/common_resources/cl/relax_multi.cl`** — OpenCL kernel with the core pairwise functions:
  - `getLJQH()` (line 164) — combined LJ + Coulomb + H-bond evaluation.
  - `getMorseQH()` (line 179) — Morse + Coulomb + H-bond.
  - `getMorsePLQH()` (line 196) — Morse with PLQ-factorized parameters.

### Key Physics

**Lennard-Jones (12-6)**:
```cpp
inline float4 getLJQH(float3 dp, float4 REQ, float R2damp) {
    float r2 = dot(dp,dp);
    float ir2 = 1.f/r2;
    float u2 = REQ.x*REQ.x*ir2;   // (R0/r)²
    float u6 = u2*u2*u2;           // (R0/r)⁶
    float vdW = u6*REQ.y;          // E0*(R0/r)⁶
    float E = (u6 - 2.f)*vdW;      // E0*[(R0/r)¹² - 2*(R0/r)⁶]
    float fr = -12.f*(u6 - 1.f)*vdW*ir2;
    // ... plus Coulomb term ...
}
```
The force is derived analytically from the energy: $\mathbf{f} = -\nabla E = -\frac{dE}{dr} \frac{\mathbf{r}}{r}$.

**Damped Coulomb**:
To avoid the singularity at $r=0$, the Coulomb potential is damped:
$$E_C = \frac{Q_i Q_j}{\sqrt{r^2 + R_{\text{damp}}^2}}$$
With $R_{\text{damp}} = 1.0$ Å by default. This is equivalent to a Gaussian charge distribution of width $R_{\text{damp}}$.

**Morse Potential**:
$$E_M = D_e \left[ \left(1 - e^{-a(r-r_0)}\right)^2 - 1 \right]$$
The parameter $a$ controls the width of the well. In FireCore, the Morse parameters are derived from the LJ parameters: $D_e = \varepsilon$, $a = \alpha_{\text{Morse}} / R_0$ with $\alpha_{\text{Morse}} = 1.5$ by default.

### Performance Considerations

- **OpenMP Parallelization**: The CPU path uses `#pragma omp parallel for` over atom pairs.
- **SIMD (AVX)**: Optional AVX vectorization for the inner loop (see `apos_simd`, `REQs_simd` in `NBFF.h`).
- **PBC Shifts**: Precomputed `shifts[npbc]` array stores all periodic image offsets. For $n_{\text{PBC}} = (2,2,0)$, this yields $5 \times 5 \times 1 = 25$ images. The shifts are reused across all pair evaluations.
- **AABB Acceleration**: For short-range interactions, axis-aligned bounding boxes (AABBs) are constructed per atom or per molecule. The `Buckets` spatial hash (`pointBBs`) skips distant pairs, reducing complexity from O(N²) to O(N) for sparse systems.

### Exclusion Schemes

Bonded neighbors (1-2: directly bonded; 1-3: angle neighbors) must be excluded from non-bonded evaluation to avoid double-counting. FireCore supports two approaches:

1. **Explicit Exclusion List** (`excl[]`):
   - Each atom stores up to `EXCL_MAX=16` excluded neighbor indices.
   - During pair evaluation, the inner loop checks `if (j in excl[i]) skip;`.
   - Simple and exact, but adds branching overhead.

2. **Subtraction Method**:
   - Evaluate *all* pairs with the full non-bonded potential.
   - Subtract the non-bonded contribution for bonded pairs.
   - Used in some UFF variants where the bonded and non-bonded potentials share the same parameterization.
   - Advantage: no branching in the inner loop. Disadvantage: requires computing the full potential for bonded pairs, then negating it.

The choice is controlled at compile time or runtime via `bSubtractBonded` flags.

### Test Coverage

- `tests/tUFF/test_UFF_multi.py` — tests that non-bonded exclusions are handled correctly in multi-system mode.
- `tests/tMMFF/test_diamond_phonon_bands.py` — phonon bands are sensitive to non-bonded cutoff; validates that bulk moduli are reasonable.

---

## 2. Fast Multipole Method (FMM)

### Physics & Purpose

For systems with many atoms (>10⁴), brute-force O(N²) electrostatics becomes prohibitive. The Fast Multipole Method (FMM) approximates long-range interactions by grouping distant atoms into clusters and representing each cluster by its multipole moments (charge, dipole, quadrupole, ...).

FireCore implements a **tile-based FMM** optimized for GPU:
- Single layer (no hierarchical tree) to preserve cache locality.
- Clusters are fixed-size "tiles" that map directly to OpenCL workgroups.
- Interactions between nearby tiles use brute-force; distant tiles use multipole expansion.
- A smooth switching function blends the two regimes to avoid force discontinuities.

### Implementation Files

- **`cpp/common_resources/cl/FMM.cl`** — OpenCL kernel with multipole force functions:
  - `calc_force_MM()` (line 90) — monopole-monopole ($1/R^2$).
  - `calc_force_MD()` (line 104) / `calc_force_DM()` (line 119) — monopole-dipole ($1/R^3$).
  - `calc_force_DD()` (line 134) — dipole-dipole ($1/R^4$).
  - `calc_force_MQ()` (line 155) / `calc_force_QM()` (line 170) — monopole-quadrupole ($1/R^4$).
  - `calc_force_DQ()` (line 185) — dipole-quadrupole ($1/R^5$).
- **`cpp/common/math/Multipoles.h`** — C++ multipole math utilities:
  - `project()` — projects a charge distribution onto multipole moments up to quadrupole order.
  - `Emultipole()` / `EFmultipole()` — evaluates energy and force from multipole moments.
  - `center()` — computes the center of charge (dipole-minimizing origin).
- **`doc/FMM/FMM.md`** — Detailed mathematical derivation of energy blending, force distribution, and switching functions.

### Key Physics

**Multipole Expansion**:
For a cluster of charges $\{q_i\}$ at positions $\{\mathbf{r}_i\}$ relative to the cluster center $\mathbf{R}$, the potential at a distant point $\mathbf{r}$ is:

$$\Phi(\mathbf{r}) = \frac{Q}{|\mathbf{r}-\mathbf{R}|} + \frac{\mathbf{p} \cdot (\mathbf{r}-\mathbf{R})}{|\mathbf{r}-\mathbf{R}|^3} + \frac{1}{2} \sum_{\alpha\beta} Q_{\alpha\beta} \frac{(r_\alpha-R_\alpha)(r_\beta-R_\beta)}{|\mathbf{r}-\mathbf{R}|^5} + \dots$$

where:
- $Q = \sum_i q_i$ — total charge (monopole).
- $\mathbf{p} = \sum_i q_i (\mathbf{r}_i - \mathbf{R})$ — dipole moment.
- $Q_{\alpha\beta} = \sum_i q_i (3(r_{i\alpha}-R_\alpha)(r_{i\beta}-R_\beta) - \delta_{\alpha\beta}|\mathbf{r}_i-\mathbf{R}|^2)$ — quadrupole tensor.

**Energy Blending**:
To avoid force discontinuities at the transition radius $R_{\text{min}}$, the energy is blended:

$$E_{\text{total}} = (1 - S(R)) E_{\text{exact}} + S(R) E_{\text{approx}}$$

where $S(R)$ is the **smootherstep** function:

$$S(x) = 6x^5 - 15x^4 + 10x^3, \quad x = \frac{R - R_{\text{min}}}{R_{\text{max}} - R_{\text{min}}}$$

This ensures $S$, $S'$, and $S''$ are all continuous at the boundaries, preventing energy drift in NVE/NVT ensembles.

**Force on Individual Atoms**:
The force from the approximate multipole interaction is distributed to atoms according to their mass and charge:

$$\mathbf{F}_k = S(R) \left[ \frac{m_k}{M_A} \mathbf{F}_{\text{inter}} + q_k \mathbf{E}_{\text{local}} \right] + (1-S(R)) \mathbf{F}_k^{\text{exact}} + (E_{\text{approx}} - E_{\text{exact}}) \frac{dS}{dR} \frac{m_k}{M_A} \hat{\mathbf{n}}$$

The three terms are:
1. **Multipole force** — proportional to mass and charge.
2. **Exact force** — standard Coulomb, scaled by $(1-S)$.
3. **Switching correction** — ensures energy conservation during the transition.

### Performance Considerations

- **Tile Size**: Typical tile size is 32–64 atoms. Too small → high overhead. Too large → multipole approximation is inaccurate.
- **Register Pressure**: Storing a full quadrupole tensor (9 floats) per tile increases register usage. The current implementation stores packed `float16` moments, which is efficient but limits the maximum expansion order.
- **Atomic Adds**: Force accumulation uses `atomic_add_float3()` (via `cl_khr_fp_atomic` extension). If unavailable, a compare-and-swap loop is used as fallback.
- **Single Precision**: The GPU kernel uses `float` throughout. For high-accuracy electrostatics, the CPU reference in `Multipoles.h` uses `double`.

### Status

- **Experimental**. The tile-based FMM is functional for test cases but not yet integrated into the main production MD loop. The primary bottleneck is tuning $R_{\text{min}}$/$R_{\text{max}}$ for different systems.

---

## 3. REQ → PLQ Conversion

### Physics & Purpose

For grid-based evaluation (see `surface_interactions.md`), the standard Lennard-Jones parameters $R_0$ (van der Waals radius) and $\varepsilon$ (well depth) are converted to a factorized **PLQ** form:

| Symbol | Physical Meaning | Formula |
|--------|------------------|---------|
| $P$ | Pauli repulsion strength | $P = \sqrt{R_0 \cdot \varepsilon}$ |
| $L$ | London dispersion strength | $L = R_0^6 \cdot \varepsilon$ |
| $Q$ | Charge | $Q$ (unchanged) |

This factorization allows the same grid to be precomputed once and reused for different atom types by simple scaling:

$$E_{\text{Pauli}}(r) = P \cdot \frac{R_0^6}{r^{12}}, \quad E_{\text{London}}(r) = L \cdot \frac{1}{r^6}$$

For mixed interactions between species $i$ and $j$, the combination rules become:

$$P_{ij} = \sqrt{P_i P_j}, \quad L_{ij} = \sqrt{L_i L_j}$$

which is naturally compatible with the geometric-mean convention used in most force fields.

### Implementation

- `NBFF::evalPLQs()` / `NBFF::makePLQs()` in `NBFF.h:179-194`.
- `RigidBodyDynamics.py` performs the same conversion for mixed-species rigid-body simulations.

### Performance Impact

- **Grid Reuse**: Without PLQ factorization, a separate grid would be needed for every atom-type pair (O(N_types²) grids). With PLQ, only 3 grids (Pauli, London, Coulomb) are needed, regardless of atom types.
- **Memory Bandwidth**: Grid sampling is memory-bandwidth bound. Reducing the number of grid channels from O(N_types²) to 3 significantly improves throughput.

---

## Summary Table

| Method | Complexity | Range | GPU | Status |
|--------|-----------|-------|-----|--------|
| **NBFF (brute)** | O(N²) | All | OpenCL | Production |
| **NBFF_SR (AABB)** | O(N) | Short | CPU (OpenMP) | Production |
| **FMM (tile)** | O(N log N) | Long | OpenCL | Experimental |
| **Ewald/PME** | O(N log N) | Long | OpenCL (partial) | Experimental |

---

## See Also

- [Topical Audit Index](topical_audit.md) — priority ranking, dependency graph, missing topics
- [Forcefields Overview](forcefields_overview.md) — high-level taxonomy of all force field classes
- [Intramolecular Forcefields](intramolecular_forcefields.md) — UFF, MMFFsp3, ProjectiveDynamics, XPBD, RigidBody
- [Surface Interactions](surface_interactions.md) — GridFF, FoldedAtomicFunctions, Ewald2D
- [Web Force Fields](forcefields_web_implementation.md) — WebGL/WebGPU shader implementations

---

*Last updated: 2026-06-23*
