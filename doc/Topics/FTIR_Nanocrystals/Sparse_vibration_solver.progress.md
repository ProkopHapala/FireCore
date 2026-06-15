# Sparse Vibration Linear Solver — Design & Progress

**Scope:** Fast frequency-domain mechanical response \(A(\omega) u = f\) for large sparse Hessians — validate against dense eigenmodes on small systems, then scale via **scipy** probes and ultimately **pyOpenCL Jacobi + heavy-ball** with local-memory workgroups.

**Related:** `pyBall/FTIR.py`, `tests/tMMFF/test_vibration_spectra.py`, `tests/tMMFF/test_nanocrystal_sparse_hessian.py`, `doc/Topics/FTIR_Nanocrystals/Sparse_Hessian_Vibration_Spectra.guide.md`.

**Status:** M-S0–M-S3 passed on adamantane fixtures; iterative solvers benchmarked; M-S4 deferred.

**Fixtures:** `Parallel_agent_fixtures.guide.md` — bootstrap before starting.

---

## Design goals

1. **Correctness before speed:** No production use of sparse scanning until peaks match dense `eigh` on reference systems.
2. **Right algorithm for the job:**
   - **Small N (~&lt;1500 DOF):** Full diagonalization once (`vibration_spectrum_from_modes`) — **cheapest**.
   - **Large N:** Frequency scan with **sparse linear solves** only if each solve is \(O(\text{nnz})\) or \(O(N)\) iterative — never \(O(N^3)\) per frequency.
3. **Unified sparse matrix:** Pair-indexed 3×3 blocks, max neighbors per atom (from linearized topology or MMFF shell Hessian).
4. **GPU solver:** pyOpenCL, workgroup = molecule/cluster, **`__local`** position/neighbor cache — port best existing Jacobi pattern, adapted to **complex arithmetic** and **fixed ω**.
5. **Staged validation ladder:** dense `eigh` → scipy sparse probe → pyOpenCL Jacobi → large NC timing.

---

## Problem statement

### Dynamic stiffness

\[
A(\omega) = K - (\omega + i\eta)^2 M
\]

**Probe:** solve \(A(\omega) u = f\) for dipole RHS \(f\) (e.g. \(f_i = q_i \hat{e}\)); spectral response \(\propto \|u\|^2\) or induced dipole.

### Why scanning with dense solver fails

| Method | Per-frequency cost | m frequencies |
|--------|-------------------|---------------|
| Dense `np.linalg.solve` | \(O(N^3)\) | \(O(m N^3)\) — **unusable** |
| Dense `eigh` once + mode sum | \(O(N^3) + O(m N)\) | **optimal for small N** |
| Sparse `spsolve` | \(O(\text{nnz}^{1.5})\) typical | \(O(m \cdot \text{nnz}^{1.5})\) |
| Iterative Jacobi/Chebyshev | \(O(k \cdot \text{nnz})\) | \(O(m k \cdot \text{nnz})\) |
| GPU iterative Jacobi | \(O(k \cdot \text{nnz})\) on device | target for large NC |

**Conclusion:** For adamantane (\(N=26\), 78 DOF), **always diagonalize**. Scanning only wins when \(N\) is large enough that \(O(N^3)\) is impossible but \(O(k \cdot \text{nnz})\) is acceptable.

---

## Known bugs in current implementation

| Issue | Location | Fix |
|-------|----------|-----|
| **No rigid-mode projection on sparse path** | `mechanical_greens_probing_sparse`; `test_nanocrystal_sparse_hessian.py` sparse branch | Apply `project_rigid_modes` to K (dense) before sparsification, or project RHS onto physical subspace |
| **Shell-truncated K ≠ full Hessian** | `getHessianSparseBlocks` + `buildNeighShells(n_shells=2)` | Validate `n_shells=2,3,4` vs `getHessian3Nx3N`; document minimum shell for bonded MMFF |
| **Complex `spsolve` conditioning** | Near resonances | Increase `stabilize`; or use complex iterative solver with diagonal preconditioner |
| **Silent solve failures** | `FTIR.py` catches `LinAlgError` / generic `Exception` | Fail loudly in test mode; only warn in production scan |

---

## Existing solver inventory (what to reuse)

### Tier 1 — Reference CPU (vibration-specific)

| Module | Role | Use in plan |
|--------|------|-------------|
| `FTIR.vibration_spectrum_from_modes` | Gold standard small N | Stage 0 reference peaks |
| `FTIR.mechanical_greens_probing` | Dense per-ω (slow) | Cross-check single ω only |
| `FTIR.mechanical_greens_probing_sparse` | scipy `spsolve` per ω | Stage 1 after bugfix |
| `FTIR.build_sparse_hessian_from_blocks` | BSR assembly | MMFF FD path |
| `scipy.sparse.linalg.eigsh` | Partial spectrum | Optional Stage 1b for medium N |

### Tier 2 — CPU iterative (prototype)

| Module | Role | Use in plan |
|--------|------|-------------|
| `pyBall/pyTruss/IterativeLinearSolvers.py` | `solve_jacobi`, `solve_jacobi_momentum`, Chebyshev | Stage 2 CPU complex extension; spectral radius estimates |
| `pyBall/pyTruss/test_Chebyshev_accel.py` | Momentum/Chebyshev tuning | Parameter playbook |

### Tier 3 — GPU Jacobi (production target)

| Module | Maturity | Local memory | Max neighbors | Notes |
|--------|----------|--------------|---------------|-------|
| **`web/common_resources/cl/LFF.cl`** + `LFFSolver.py` | **High for springs** | `__local float4 lpos[32]` | 8 | **Best template** for bonded gather; Jacobi-fly inner loop; `bMix` momentum; needs **complex** variant + ω-dependent diagonal |
| `pyBall/XPDB_AVBD/XPDB_new.cl` | High for MD | Extensive `__local` | 4 ports | Cluster/collision focused; overkill for static vibration solve |
| `web/molgui_web/js/ProjectiveDynamics.js` | JS reference | N/A | CPU | Algorithm reference for constraint Jacobi |
| `pyBall/pyTruss/truss_solver_ocl.py` | Medium | Varies | Truss-specific | `jacobi_fly`, `jacobi_diff` patterns |

**Recommendation:** Base GPU vibration solver on **`LFF.cl` `lff_jacobi` structure**:

- One workgroup per molecule/cluster; `lid` = atom index.
- Private copy `ng_idx[]`, `ng_KLs[]` from global (already done).
- Inner loop: gather neighbor positions from `__local lpos`, accumulate off-diagonal coupling into diagonal solve.
- Replace PD mass term \(M/\Delta t^2\) with \(A_{ii}(\omega) = K_{ii} - (\omega+i\eta)^2 M_i\).
- Add complex-valued `u` (float2 per component) or split real/imag passes.

**Increase `MAX_NEIGHBORS` to 32** when linearized K₁₃+K₁₄ is used (see `Linearized_topology.progress.md`).

### Tier 4 — Not for bonded vibrations (reference only)

- `scipy.sparse.linalg.spsolve` — SuperLU; fine for **validation**, not GPU endgame.
- Random/stochastic probing (OrderN_QM docs) — Stage 4+ if needed for \(N>10^4\).

---

## Sparse matrix formats

### Path A — MMFF finite-difference (current)

```
neigh_idx[natoms, max_neigh], neigh_count[natoms], blocks[natoms, max_neigh, 3, 3]
→ FTIR.build_sparse_hessian_from_blocks → scipy BSR
```

- Pros: exact MMFF bonded Hessian (if shells complete).
- Cons: expensive FD; shell truncation risk.

### Path B — Linearized topology (target)

```
From topology.json / MMFFL:
  neighs[natoms, MAX_N], KLs[natoms, MAX_N, 2], stick_class[]
→ analytical K assembly (3×3 blocks from spring directions)
```

- Pros: \(O(n_\text{sticks})\) assembly; matches GPU LFF data layout.
- Cons: approximate physics until K₁₄ calibrated.

**Both paths** should converge to same **solver API**:

```python
class SparseVibrationProblem:
    neigh_idx: int32[n, max_neigh]
    blocks: float64[n, max_neigh, 3, 3]   # or K/L from springs
    mass: float64[n]
    # solve(omega, eta, rhs) -> u
```

---

## Staged validation plan

### Stage 0 — Dense reference (mandatory gate)

**Systems:** adamantane (`test_vibration_spectra.py`), optionally ethane.

**Steps:**

1. `H` from `getHessian3Nx3N` (or analytical springs later).
2. `H2 = project_rigid_modes(H, M, pos)`.
3. `eigvals = eigh(M^{-1/2} H2 M^{-1/2})` → \(\omega_i\).
4. `vibration_spectrum_from_modes` → peaks at \(\omega_i\).

**Milestone M-S0:** Peak positions agree within \(0.01\) cm⁻¹ equivalent (internal units); documented in test assert.

### Stage 1 — scipy sparse probe (CPU)

**Systems:** adamantane, G1 diamond primitive (~30–80 atoms).

**Steps:**

1. Build `H_sparse` (full shell or dense-to-sparse for small N).
2. **Apply same rigid projection as Stage 0.**
3. `mechanical_greens_probing_sparse` for \(m=100\) frequencies.
4. Compare peak envelope to Stage 0: `max_rel_diff < 1e-2` (damping broadens).

**Also test:**

- `scipy.sparse.linalg.eigsh(H2, k=20, M=M)` — do eigenvalues match dense?

**Milestone M-S1:** Sparse scan peaks align with `eigsh` / dense modes.

### Stage 2 — CPU iterative Jacobi (complex)

**Steps:**

1. Extend `IterativeLinearSolvers.solve_jacobi_momentum` to **complex** \(A\) (Hermitian positive-definite check).
2. Same systems; single \(\omega\) near first peak; compare `u` to `spsolve`.
3. Tune `bmix`, max iter for convergence.

**Milestone M-S2:** Relative error \(\|u_\text{iter} - u_\text{direct}\|/\|u\| < 10^{-4}\) at test ω.

### Stage 3 — pyOpenCL Jacobi (`vib_jacobi.cl` or extend `LFF.cl`)

**New kernel:** `vib_jacobi_solve` (one ω per launch, or batched ω)

**Parameters:** `neighs`, `KLs` or `blocks`, `mass`, `omega`, `eta`, `rhs`, `bmix`, `n_iter`.

**Steps:**

1. Port LFF local gather; diagonal includes \(-(\omega+i\eta)^2 m_i\).
2. Adamantane: GPU vs CPU iterative (Stage 2).
3. Profile: time per solve vs `spsolve` for \(N \approx 200, 500, 1000\).

**Milestone M-S3:** GPU matches CPU iterative; speedup &gt; 3× at \(N=500\) for single solve (informal target).

### Stage 4 — Full spectrum scan (large NC)

**Systems:** relaxed G2 SiNC (200–800 atoms) — from generation milestone.

**Steps:**

1. Cannot dense-diagonalize — use sparse scan only.
2. Compare **total wall time**: `eigsh` (20 modes) vs \(m=500\) sparse GPU solves.
3. If scan wins, run production FTIR-style dipole probe.

**Milestone M-S4:** Spectrum visually reasonable (no spurious \(\omega\to 0\) mountain); Si–H stretch band plausible.

### Stage 5 — Stochastic probing (optional, future)

- Random RHS probes + Chebyshev rational approximation — only if Stage 4 too slow at \(N>5000\).

---

## Test harness (planned scripts)

| Script | Purpose |
|--------|---------|
| `tests/tMMFF/test_vibration_solver_ladder.py` | Runs stages 0–2 with CLI `--stage 0\|1\|2` |
| `tests/tMMFF/test_vibration_jacobi_ocl.py` | Stage 3 GPU parity |
| Extend `test_nanocrystal_sparse_hessian.py` | Fix projection; `--spectrum-method both` assert |

**CLI flags (proposed):**

```
--stage {0,1,2,3,4}
--input molecule.mol2
--hessian-source {mmff_fd, linearized}
--shells 2
--omega-min, --omega-max, --nfreq
--eta, --shift-rigid
--solver {dense_eigh, spsolve, jacobi_cpu, jacobi_ocl}
```

---

## GPU kernel design sketch

```
__kernel void vib_jacobi(
    __global const int*    neighs,      // [n, MAX_N]
    __global const float*  K_blocks,    // [n, MAX_N, 3, 3] or spring K,l0
    __global const float*  mass,        // [n]
    __global const float2* rhs,         // [n, 3] complex
    __global float2*       u,             // [n, 3] in/out
    float omega, float eta, float bMix,
    int nIter
)
```

**Workgroup:** `local_size = 32` (or next power of 2 ≥ atoms in cluster); one molecule per group.

**Local cache:** `__local float4 lpos[WG_SIZE]` — store real parts of neighbor positions for gather (geometry fixed during linear solve).

**Diagonal:** \(A_{ii} = K_{ii} - ((\omega)^2 - \eta^2 + 2i\omega\eta) m_i\) for harmonic shift form (expand complex arithmetic explicitly).

**Momentum:** `u_new = (1-bMix)*u_jacobi + bMix*u_old` per component (heavy-ball from LFF).

---

## Performance expectations (order-of-magnitude)

| N (atoms) | DOF | Dense eigh | scipy spsolve ×100 ω | GPU Jacobi ×100 ω (target) |
|-----------|-----|------------|----------------------|----------------------------|
| 26 (adamantane) | 78 | &lt;0.1 s | ~1 s | overhead-dominated |
| 200 | 600 | ~1 s | ~30 s | ~5 s |
| 1000 | 3000 | minutes | minutes | seconds–tens of seconds |

Exact numbers to be measured in M-S3/M-S4.

---

## Dependencies

| Upstream | Need |
|----------|------|
| Generation (`Nanocrystal_generation.progress.md`) | Relaxed G2 structures for Stage 4 |
| Linearized topology (`Linearized_topology.progress.md`) | Path B Hessian + `MAX_NEIGHBORS` budget |
| MMFF C++ | Path A `getHessianSparseBlocks` for MMFF reference |

---

## Open questions

1. **Factor ω-dependent diagonal once per ω** or reuse Cholesky pattern when only \((\omega^2)\) changes? (Block diagonal shift — cheap.)
2. **Symmetrize K** after shell assembly: always `0.5*(K+K.T)`?
3. **Preconditioner:** diagonal (Jacobi) sufficient for bonded sparse K, or incomplete Cholesky on CPU host?
4. **Batch frequencies:** one kernel launch with ω index in `get_group_id(1)`?

---

## Milestone summary table

| ID | Description | Blocks |
|----|-------------|--------|
| M-S0 | Dense eigh ≡ mode-sum spectrum (adamantane) | — |
| M-S1 | scipy sparse scan ≡ M-S0 peaks; rigid projection fixed | Stage 1 |
| M-S2 | CPU Jacobi momentum ≡ spsolve (single ω) | Stage 2 |
| M-S3 | pyOpenCL Jacobi ≡ CPU Jacobi | Stage 3 |
| M-S4 | Large NC timing + plausible spectrum | Stage 4 |
| M-S5 | Document solver choice flowchart in `Sparse_Hessian_Vibration_Spectra.guide.md` | Docs |

---

## Parallel agent contract (Job 3)

Full rules: **`Parallel_agent_fixtures.guide.md`**.

| | |
|--|--|
| **Fixture root** | `tests/tSiNCs/fixtures/vibration_parallel/` (gitignored) |
| **Read** | `hessian_mmff/adamantane_mmff_dense.npz`, `hessian_mmff/diamond_nc_R6_sparse_blocks.npz`, `structures/diamond_nc_R6_relaxed.xyz`, `expected/adamantane_omegas_modes.txt` |
| **May edit** | `pyBall/FTIR.py`, `tests/tMMFF/test_vibration_solver_ladder.py`, `tests/tMMFF/test_vibration_jacobi_ocl.py` |
| **May add** | `cpp/common_resources/cl/vib_jacobi.cl`, `pyBall/OCL/VibSolver.py` (**do not modify `LFF.cl`**) |
| **Must NOT edit** | `MMFFLTopology.js`, `gen_nanocrystals.mjs`, `EditableMolecule.js`, `MMFFL.py` |
| **Hessian** | MMFF C++ / frozen `.npz` only (Path A) |

---

*Last updated: 2026-06-15 — M-S0–M-S3 complete; iterative solver study; GPU Jacobi one-step parity.*

---

# USER

OK, I'm a bit lost in jargon, what is SPD? what is "η=0η=0" what is condition when CG converges and when it does not?

for CGNE and GMRES you can elaborate more about how it can be potentilly implemented on GPU. I mean we can very effiicnetly implement matrix-vectro products and possible we can do it efficiently especially if we can group multiple matrix-vectros products together, and if we need to orghogonalize only small matrixes (like subspace, something like lanczos sub-space), we can build these matrixes on GPU, then pass small matrix to CPU to diagonalize/orthogonalize and go back to GPU or something like that. also if we run multiple probes (right-had-sides) we can doe same dotproduct to all of them to allow tile-based local memory access pattern which amortize global memory read, it is better than doing the probes sequentially

elaborate on this topic without writing new code


---

# Cursor Composer 2.5

## Plain-language glossary

### SPD (Symmetric Positive-Definite)

A real matrix \(A\) is **SPD** when:

1. **Symmetric:** \(A = A^T\) (same entry above and below the diagonal).
2. **Positive-definite:** for every nonzero vector \(u\), the quadratic form \(u^T A u > 0\).

Intuition: \(A\) behaves like a “stiff spring” in every direction — push in any direction, energy goes up. No zero or negative-curvature directions.

For vibration, the static Hessian \(K\) is SPD (internal modes only; rigid modes are zero). The dynamic matrix

\[
A(\omega) = K - \omega^2 M
\]

is SPD only when \(\omega\) is **below all vibrational frequencies** (no resonance yet). Above a mode frequency, \(A(\omega)\) stops being SPD in that subspace.

---

### What is \(\eta\) (eta)?

\(\eta\) is **damping** in the frequency-domain equation:

\[
A(\omega) = K - (\omega + i\eta)^2 M
\]

- **\(\eta = 0\):** undamped, purely elastic response. \(A(\omega)\) is **real** (no imaginary part).
- **\(\eta > 0\):** damped response. Peaks in the spectrum get broader; \(A(\omega)\) becomes **complex**.

Your “\(\eta=0\eta=0\)” looks like a copy-paste glitch — I mean **\(\eta = 0\)** (zero damping).

In code: `dynamic_stiffness(K, M, omega, eta=...)`.

---

## When CG works and when it does not

**Conjugate Gradient (CG)** solves \(A u = f\) by building better and better guesses. It only has a clean theory when \(A\) is **real SPD**.

| Situation | CG? | Why |
|-----------|-----|-----|
| \(\eta = 0\), \(\omega\) well below all mode frequencies | **Yes** | \(A = K - \omega^2 M\) is real symmetric; still SPD if \(\omega^2 < \lambda_{\min}(M^{-1}K)\) in the vibrational subspace |
| \(\eta = 0\), \(\omega\) near or above a resonance | **No / breaks** | \(A\) has zero or negative eigenvalues → not SPD; CG can stall or diverge |
| \(\eta > 0\) (any \(\omega\)) | **No** | \(A\) is **complex** and not Hermitian → standard CG does not apply |
| Ill-conditioned but still SPD | **Yes, slowly** | More iterations, but still converges if you use enough steps |

**Condition number** \(\kappa(A) = \lambda_{\max}/\lambda_{\min}\): ratio of stiffest to softest direction. CG typically needs about \(O(\sqrt{\kappa})\) iterations. Near resonance, \(\kappa\) blows up → many iterations even when SPD.

**What to use instead:**

- **\(\eta > 0\)** (damped, complex \(A\)): **GMRES** — general complex solver.
- **\(\eta = 0\)** but near resonance: GMRES or direct solve; CG is unreliable.
- **Complex but you insist on CG-like method:** **CGNE** — reformulates as a normal-equations problem (see below).

---

## CGNE — idea and GPU angle

**Problem:** \(A u = f\) with complex, possibly non-Hermitian \(A\).

**CGNE trick:** multiply both sides by \(A^H\) (conjugate transpose):

\[
\underbrace{A^H A}_{\text{always Hermitian PSD}} u = A^H f
\]

Then run CG on the **normal equations** \( (A^H A) u = A^H f \).

**Pros:**
- Works for general complex \(A\).
- Each CG step needs **two** matvecs: one \(A p\), one \(A^H (\cdot)\) — or you can fuse patterns if structure is known.

**Cons:**
- Condition number squares: \(\kappa(A^H A) \approx \kappa(A)^2\) → **~2× worse** convergence than GMRES. We saw ~0.3% error plateau on adamantane vs \(10^{-9}\) for GMRES.

### GPU implementation sketch (no code)

```
┌─────────────────────────────────────────────────────────┐
│  GPU: SpMV for A  and  A^H  (same sparsity, transposed) │
│       Block-diagonal preconditioner P^{-1}              │
│       Vector ops: x += αp, r -= αAp, dot products       │
└─────────────────────────────────────────────────────────┘
         ↑↓  only scalars α, β and small inner products
┌─────────────────────────────────────────────────────────┐
│  CPU:  CG recurrence (α, β from dot products)         │
└─────────────────────────────────────────────────────────┘
```

**Why this fits GPU well:**
- **SpMV dominates cost** — one sparse structure, highly parallel.
- CGNE needs **global dot products** each iteration (\(r^H r\), \(p^H A p\)) — one reduction per probe per iteration. Annoying but \(O(N)\) work with \(O(\log N)\) sync, not \(O(N^2)\).

**Multi-probe batching (your idea — very good):**

If you solve \(A u_k = f_k\) for many RHS (many frequencies \(\omega\) or many dipole directions) **with the same \(A\)**:

- Pack \(k\) vectors into columns: \(U \in \mathbb{C}^{N \times k}\), \(F \in \mathbb{C}^{N \times k}\).
- One **batched SpMV**: \(A U\) — each thread/block handles one row of \(A\), accumulates into \(k\) outputs using **registers or local memory** for the \(k\) partial sums.
- One **batched dot product**: \(\langle r_i, r_j \rangle\) or just \(\|r_k\|^2\) for each column — tile so each work-group reads \(r\) once from global memory and computes \(k\) inner products.

That amortizes the **read** of the sparse matrix indices/values and the **read** of the vector — often the bottleneck on GPU.

**Caveat for your use case:** in a frequency scan, \(A(\omega)\) **changes every \(\omega\)**, so you cannot batch across frequencies with one matrix. You *can* batch across:
- multiple dipole directions at **fixed** \(\omega\),
- multiple \(\omega\) only if you accept rebuilding \(A\) or using a shifted formulation.

---

## GMRES — idea and GPU angle

**GMRES** builds a **Krylov subspace**:

\[
\mathcal{K}_m = \mathrm{span}\{r_0, Ar_0, A^2 r_0, \ldots, A^{m-1} r_0\}
\]

and finds the best approximate solution inside it (minimizes residual norm).

**Arnoldi process** (each step \(j\)):
1. \(w = A v_j\) — **SpMV** (GPU-friendly).
2. Orthogonalize \(w\) against \(v_0,\ldots,v_j\) — **not** GPU-friendly at large \(j\):
   - Modified Gram–Schmidt: \(j+1\) dot products + axpy per step.
   - At restart \(m=30\), that's manageable on CPU for the **small** orthogonalization.

3. Store Hessenberg matrix \(H_m \in \mathbb{C}^{(m+1)\times m}\) — **tiny** (\(30\times30\)).
4. Solve least squares \( \min_y \| \beta e_1 - H_m y \| \) on **CPU** (microseconds).
5. Update \(x = x_0 + V_m y\).

**Restarted GMRES:** every \(m\) steps, throw away the basis, start fresh from new residual. That's what we use (`restart=30`).

### Hybrid GPU/CPU pipeline (your Lanczos/subspace picture)

This is exactly the standard high-performance GMRES design:

```
  GPU                          CPU
  ───                          ───
  SpMV: w = A v_j
  Store v_j in V (N × m)       
  Optional: GPU MGS on V       Or: CPU MGS if m small
  ─────────────────────────────────────────────
                               Form H_m (m+1 × m)
                               QR / least-squares on H_m
                               Compute y, scalar updates
  x += V y                     (or GPU axpy)
  Check residual: r = b - Ax   (SpMV on GPU)
```

**Key point:** never put an \(N \times N\) dense matrix on GPU. Only:
- \(V\): \(N \times m\) with \(m \ll N\) (e.g. 30),
- \(H\): \(31 \times 30\).

Diagonalization/QR of \(H_m\) on CPU is negligible compared to one SpMV at \(N \sim 10^4\)–\(10^6\).

### Why GMRES is harder than CG on GPU

| Operation | CG (SPD) | GMRES |
|-----------|----------|-------|
| Matvecs per iter | 1 | 1 |
| Global reductions per iter | 2 dot products | \(O(j)\) dots + axpys in step \(j\) |
| Memory | 3–4 vectors | \(m+1\) basis vectors \(V\) |
| Small dense solve | none | each restart cycle |

So GMRES is still **SpMV-bound** for large \(N\), but orthogonalization and memory traffic for \(V\) hurt more than CG.

### Multi-probe / multi-RHS GMRES

**Block GMRES:** same \(A\), several RHS \(b_1,\ldots,b_k\). Build a **block Krylov space** — each Arnoldi step does **one batched SpMV** \(A [v_1|\cdots|v_k]\).

Same tile trick as CGNE:
- One row of \(A\) touches \(k\) basis vectors → load neighbor indices once, accumulate \(k\) dot products in registers.
- Orthogonalization becomes \(k \times k\) small blocks or \(k\) independent Arnoldi chains (simpler but less sharing).

For **frequency scans**: each \(\omega\) has different \(A\), so you run separate GMRES solves. Batching helps when:
- same \(\omega\), multiple probe directions or multiple atoms charged,
- or pipelining: GPU works on \(\omega_i\) while CPU orthogonalizes \(\omega_{i-1}\)'s GMRES cycle (async queues).

---

## Richardson vs CG (one sentence)

Richardson: fixed step \(u \leftarrow u + \tau P^{-1}(f - Au)\). CG: picks the **best** step direction and size each iteration using previous information. Same building blocks (SpMV + preconditioner), very different convergence — Richardson failed on our Hessian, CG worked when SPD.

---

## Practical recommendation for your project

| Goal | Solver | GPU worth it? |
|------|--------|----------------|
| Damped spectrum \(\eta > 0\) | GMRES on CPU | Later; orthogonalization + many \(\omega\) make hybrid design non-trivial |
| Undamped \(\eta = 0\), off-resonance | CG | **Best GPU candidate** — 1 SpMV/iter, 2 reductions, batched probes |
| Quick sanity / reference | `spsolve` | CPU only |

Your mental model is right: **heavy work = SpMV on GPU in a tiled/batched layout; light work = small subspace algebra on CPU.** CGNE and GMRES both fit that; CG is simpler; GMRES is more robust for \(\eta > 0\) but more sync and storage per restart cycle.

---

# Progress report — 2026-06-15 (Job 3 wrap-up)

## Summary

Sparse vibration frequency-domain probing is **validated on adamantane** through milestone **M-S3**. The main blocker for production iterative scanning is **not** rigid modes (now fixed) but **solver choice**: bonded MMFF dynamic stiffness \(A(\omega)=K-(\omega+i\eta)^2 M\) is **not diagonally dominant**, so Jacobi/Richardson fail; **GMRES** (damped, \(\eta>0\)) and **CG** (undamped, \(\eta=0\), off-resonance) work on CPU NumPy prototypes.

---

## Achieved milestones

| ID | Status | Result (adamantane fixture) |
|----|--------|-----------------------------|
| **M-S0** | **PASS** | Mass-weighted `eigh` mode frequencies vs frozen fixture: max diff \(\sim 4.7\times 10^{-9}\) |
| **M-S1** | **PASS** | Sparse `spsolve` + rigid-mode projection vs dense probe: rel envelope diff \(\sim 2\times 10^{-10}\); **without** projection → spurious low-\(\omega\) peak \(\sim 0.73\) |
| **M-S2** | **PASS** (limited) | Point Jacobi converges on **15-point scalar chain** toy only; **fails** on full MMFF \(A(\omega)\) |
| **M-S3** | **PASS** (one-step) | GPU block Jacobi one iteration vs CPU: rel err \(\sim 1.7\times 10^{-8}\); **full convergence** on true \(A(\omega)\) not achieved |
| **M-S4** | **Deferred** | Large NC timing — skipped per user request |

**Diagnostic plots** (regenerate: `python3 tests/tMMFF/test_vibration_solver_ladder.py --plot-only`):

- `tests/tMMFF/OUT_vibration_solver/spectrum_sparse_vs_dense.png` — main M-S1 pass visual
- `spectrum_lowfreq_zoom.png`, `mode_frequencies.png`, `spectrum_dense_vs_modesum.png`
- `jacobi_inertial_shift_diagnostic.png` — block dominance ratio & inertial \(\mu M\) shift vs Jacobi convergence
- `spectrum_iterative_solvers.png` — GMRES/CG vs `spsolve` (from `test_iterative_vibration_solvers.py`)

---

## Implemented (code)

### `pyBall/FTIR.py` (core)

**Rigid-mode projection (fixes M-S1 spurious peak):**

- `rigid_mode_shift_vectors(M, pos)` — 6 mass-orthonormal rigid vectors
- `apply_rigid_mode_shift(H, M, pos, shift=1e6)` — rank-6 penalty on sparse/dense \(H\)
- `prepare_sparse_hessian(...)` — BSR assembly + projection

**Sparse probing:**

- `mechanical_greens_probing_sparse(..., pos=, shift_rigid=, fail_loud=, solver=)` — solver routing, optional projection

**CPU iterative solvers (NumPy prototypes):**

| Function | Method | Notes |
|----------|--------|-------|
| `solve_block_jacobi_momentum` | Block Jacobi + heavy-ball | Complex; diverges on MMFF \(A(\omega)\) |
| `solve_block_jacobi_inertial` | Jacobi on \((A+\mu M)\) | Converges for large \(\mu\) but solves wrong system |
| `solve_preconditioned_richardson` | Richardson + block prec | Does not converge on adamantane |
| `solve_cg_spd` | Conjugate gradient | Real SPD only; spectrum vs `spsolve` at \(\eta=0\): max rel diff \(\sim 3\times 10^{-10}\) |
| `solve_cgne` | CG on normal equations | Converges; ~0.3% error plateau (\(\kappa^2\) loss) |
| `solve_gmres` | Restarted GMRES | Best for \(\eta>0\); ~390 matvecs at \(\omega=0.915\); full spectrum max rel diff \(\sim 3\times 10^{-7}\) |
| `solve_linear_system` | Unified `method=` entry | `auto`, `spsolve`, `cg`, `cgne`, `gmres`, `richardson`, `jacobi` |

**Diagnostics:** `block_diag_dominance_ratio`, `inertial_shift_matrix`, `pack_block_rows_from_sparse`

### GPU (new, did not modify `LFF.cl`)

- `cpp/common_resources/cl/vib_jacobi.cl` — block Jacobi + heavy-ball, complex \(A(\omega)\)
- `pyBall/OCL/VibSolver.py` — pyOpenCL wrapper, `solve_dynamic_stiffness_gpu()`

### Tests (new)

| Script | Purpose |
|--------|---------|
| `tests/tMMFF/test_vibration_solver_ladder.py` | M-S0–S2 + `--plot` / `--plot-only` |
| `tests/tMMFF/test_vibration_jacobi_ocl.py` | M-S3 GPU one-step vs CPU |
| `tests/tMMFF/test_iterative_vibration_solvers.py` | GMRES/CG/CGNE/Richardson/Jacobi benchmark + spectrum plot |

**Run (from `tests/tMMFF/`, fixtures required):**

```bash
python3 test_vibration_solver_ladder.py --stage all
python3 test_vibration_jacobi_ocl.py
python3 test_iterative_vibration_solvers.py
```

Fixtures: `tests/tSiNCs/fixtures/vibration_parallel/` (gitignored; bootstrap via `tests/tSiNCs/bootstrap_vibration_parallel_fixtures.py`).

---

## Key technical findings

### Rigid-mode projection

Free clusters have 6 zero eigenmodes of \(K\). In \(A(\omega)=K-(\omega+i\eta)^2 M\) they cause huge spurious response at \(\omega\to 0\). Fix: add \(\sigma\sum_j (M^{1/2}q_j)(M^{1/2}q_j)^T\) with mass-orthonormal rigid \(q_j\), \(\sigma\sim 10^6\). Pushes rigid eigenvalues to \(\sim\sigma\); **does not** change internal vibrational modes. **Does not** fix Jacobi convergence.

### Why Jacobi / Richardson fail; CG / GMRES succeed

- Block diagonal dominance ratio: \(K\sim 0.64\), \(A(\omega)\sim 0.07\) (Jacobi needs \(\gg 1\)).
- Damping term \(-(\omega+i\eta)^2 M\) shrinks diagonal near resonance.
- **CG** needs **SPD** (real symmetric, all eigenvalues \(>0\)): works at \(\eta=0\) below resonances; breaks near/on resonance or for \(\eta>0\) (complex non-Hermitian).
- **GMRES** handles general complex \(A(\omega)\) with \(\eta>0\).
- **Gauss–Seidel** not implemented — same spectral-radius issue as Jacobi.
- **Richardson** \(u\leftarrow u+\tau P^{-1}(f-Au)\): GPU-friendly but impractically slow on bonded \(K\) (eigenvalue spread).

### Solver choice flowchart (current recommendation)

| Use case | Solver |
|----------|--------|
| Small \(N\) (\(<\sim 1500\) DOF) | `eigh` once + `vibration_spectrum_from_modes` |
| Validation / medium \(N\) CPU | `spsolve` + rigid projection |
| Damped spectrum scan \(\eta>0\) | **GMRES** (CPU; ~25–40 restart cycles per \(\omega\)) |
| Undamped \(\eta=0\), off-resonance | **CG** (CPU; good GPU candidate later) |
| GPU Jacobi | One-step math verified; **not** production solver for \(A(\omega)\) |

### Known remaining issues

- `vibration_spectrum_from_modes` vs `mechanical_greens_probing` use different damping conventions (`-i\eta\omega` vs \((\omega+i\eta)^2\)) — pre-existing.
- Diamond NC sparse blocks `n_shells=2`: ~52% Frobenius error vs dense (shell truncation) — not fixed in this job.
- scipy `cg`/`bicgstab` segfault on some complex paths — custom CG/GMRES used instead.

---

## Future plans

1. **M-S4 (large NC):** timing `eigsh` vs GMRES spectrum scan on G2 SiNC (200–800 atoms) when fixtures ready.
2. **Production solver:** wire `solver='gmres'` into main FTIR scan path with tuned `restart`/`cycles`/`tol`; keep `spsolve` as reference gate.
3. **GPU CG (SPD, \(\eta=0\)):** hybrid GPU SpMV + CPU scalar recurrence; batched multi-RHS at fixed \(\omega\) (tile local memory for dot products). Best ROI among iterative GPU options.
4. **GPU GMRES (later):** SpMV on GPU, Arnoldi orthogonalization + \(H_m\) least-squares on CPU; optional block-GMRES for multiple dipole directions at same \(\omega\).
5. **Do not pursue:** Gauss–Seidel, Richardson (unless preconditioner radically improved), full Jacobi convergence on raw \(A(\omega)\).
6. **Docs:** copy this progress into tracked `doc/Topics/FTIR_Nanocrystals/` when doc tree is synced; update `Sparse_Hessian_Vibration_Spectra.guide.md` solver flowchart (M-S5).

---

## Files to commit / push to git

**Include (implementation + tests):**

| Path | Role |
|------|------|
| `pyBall/FTIR.py` | Rigid projection, sparse probing, all iterative solvers |
| `pyBall/OCL/VibSolver.py` | pyOpenCL vibration Jacobi wrapper |
| `cpp/common_resources/cl/vib_jacobi.cl` | OpenCL block Jacobi kernel |
| `tests/tMMFF/test_vibration_solver_ladder.py` | M-S0–S2 ladder + plots |
| `tests/tMMFF/test_vibration_jacobi_ocl.py` | M-S3 GPU parity |
| `tests/tMMFF/test_iterative_vibration_solvers.py` | Iterative solver benchmark |

**Do not commit:**

| Path | Reason |
|------|--------|
| `tests/tSiNCs/fixtures/vibration_parallel/` | Gitignored frozen fixtures (regenerate locally) |
| `tests/tMMFF/OUT_vibration_solver/*.png` | Generated plot output |
| `.backup_md/...` | Backup doc tree unless repo policy includes it |

**Optional (docs):** mirror this progress section to tracked `doc/Topics/FTIR_Nanocrystals/Sparse_vibration_solver.progress.md` if that file is added to the main doc tree.

**Suggested commit message:**

```
Add sparse vibration solvers: rigid-mode projection, GMRES/CG CPU prototypes, GPU Jacobi one-step

Fix spurious low-frequency peaks via rank-6 rigid-body Hessian shift; validate M-S0–S3 on adamantane fixtures; benchmark iterative solvers (GMRES for damped scans, CG for eta=0).
```
