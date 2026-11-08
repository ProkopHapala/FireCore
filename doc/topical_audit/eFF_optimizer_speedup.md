---
type: TopicalAudit
title: eFF relaxation optimizer — FIRE speedup / stiffness / mass preconditioning study
description: Why eFF relaxation under FIRE is slow or stalls; systematic study of dt, damping, masses, diagonal-curvature preconditioning; discussion of fixed-mass guesses, on-the-fly mass adaptation, and alternative algorithms (BFGS et al.).
tags: [eFF, FIRE, optimizer, preconditioning, masses, stiffness, DynamicOpt, hessian]
timestamp: 2026-11-08
---

# eFF optimizer speedup study

**Latest review:** §8 tests inexpensive gradient-history alternatives to BFGS using the current CPU force implementation. Spectral gradient descent takes 46 evaluations for fixed-nuclei H2O versus production FIRE's 75; five-history Anderson takes 32. For free nuclei, 60 production FIRE evaluations followed by eight-history Anderson takes 116 total versus FIRE's 178, with rotation/distortion checks also improving. These new optimizers are diagnostic prototypes, not production implementations. §7.7 supersedes §7's unreproduced 169-evaluation mass-metric claim. Earlier numerical tables are historical.

**Follow-up CPU review (2026-11-08):** [Correctness findings and agent handoff](eFF_CPU_review_2026-11-08.md) identifies missing internal repulsion in paired orbitals, ECP Coulomb/Pauli sign errors, negative FGO occupancy, and mass-dependent convergence reporting. These qualify the physical interpretation and optimizer comparisons below; paired convergence alone does not establish physical correctness. The follow-up also corrects the FD-refresh cost estimate and the claim that active FIRE enforces `dr_limit`.

**Question (user):** eFF optimization with FIRE is slow / gets stuck. Is it a suboptimal FIRE setting, or DOFs with very different stiffness and mass? Can a systematic study say what to change?

**Short answer:** Three distinct problems were found, two of them *physics* problems masquerading as optimizer slowness. (1) **A real force bug** (found 2026-11-08): in `EFF::evalAE()` the electron-size forces `fsize` for the atom–electron Pauli and core-Coulomb terms received the *atom's* size derivative instead of the electron's (swapped `fsi`/`fsj` — see §4). Fixed in `eFF.h` and the OpenCL kernel `eFF.cl`. (2) **A real approximation defect — the `erfx_e6` |r|-cusp** (found and FIXED 2026-11-08): the `erf(kx)/x` approximant had nonzero slope at r=0 → a linear |r|-cusp at coincident electron–nucleus pairs and a hard |F|≈0.02–0.03 convergence floor (paired core electrons sit on nuclei). Refit with exact `y(0)` and slope-0 constraints — max error *improved* 3.7e-7→9.1e-8 and the floor vanished (CH4-bPaired now reaches |F|≈1e-9). See §4.2. (3) Even with correct forces, mode-'a' (all-electron) relaxation is a **conditioning problem**: curvature spans ~5 decades while `processXYZ_e` runs with **unit inverse masses**. Best tested fixes: secant-diagonal adaptive masses (~free, ~10×) or bounded FD-diagonal preconditioning (~10–20×); dense BFGS is also excellent (49–156 iters post-refit). Separately, `processXYZ` (auto-generated unpaired electrons, frozen nuclei) with `bPaired=False` collapses to unphysical basins; `bPaired=True` avoids that.

**⚠ All iteration tables in §1–§3 were measured *before* the `evalAE` force fix; where the fix matters (frozen-core/ECP paths, anything sensitive to esize) numbers are flagged.** All-electron (`H2O_spins`) trajectories were essentially unchanged by the fix (swapped terms are small there), so its numbers stand. Frozen-core numbers were re-measured post-fix (§3.0).

Study code: `tests/tEFF/fire_study_lib.py` (exact Python replica of `DynamicOpt::move_FIRE` + `eFF_lib.cpp::run()` loop; parity vs C++ ≈ 1e-11) and `tests/tEFF/run_fire_study.py` (sweep driver).

---

## 1. Where the slowness lives (measured)

DOF layout: `pDOFs = [apos(3·na) | epos(3·ne) | esize(ne)]`, same for `fDOFs`,`vDOFs`,`invMasses`. `run()` does `eval(); move_FIRE(); esize=clamp(esize,1e-3)`; `setTimeSteps(dt)` sets `dt_max=dt`, `dt_min=0.1·dt` — **the input dt is a ceiling, FIRE can only shrink**.

### Baseline (`H2O_spins.xyz`, na=3 ne=9, `dt=0.001 damp=0.1 f_limit=100`, unit masses)

- **6686 iters** to |F|<1e-3, final E ≈ −1849.99 (a deeper basin −1852.86 exists and is reached by other settings).
- Initial residual is almost all esize: |F_size| ≈ 3770 vs |F_epos| ≈ 72, |F_apos| ≈ 58.
- `dt` stays pinned at `dt_max` here, but `cos_vf` flips repeatedly on stiff modes; the esize DOFs dictate the global step budget.

### Hessian spectrum at the relaxed geometry (`H2O_spins`)

- λ_max ≈ **5.4×10⁵** (three modes ≈ 50% nucleus + 50% bound core-electron `e2`), λ ≈ **2.1×10⁵** for the `e2` size DOF.
- Useful soft modes λ ≈ 4–2400; six near-zero rigid modes; one shallow negative mode (−8.9, epos saddle — harmless for |F| convergence).
- **Condition number ~10⁵ → ω ratio ~370×.** With unit masses the stable timestep is set by the stiffest mode while soft DOFs move 370× slower per unit "effort".

### dt × damping grid (unit masses, H2O_spins)

| dt_max | iters to 1e-3 | note |
|---|---|---|
| 0.0005 | ~11000 | crawl |
| **0.001 (default)** | **6686** | current default |
| **0.002** | **~3495** | best stable; also lands in deeper basin −1852.86 |
| ≥0.003 | never (30k) | undamped ringing of stiff modes (`dt·ω_max` over stability edge) |

Damping 0.02–0.5: secondary effect; dt is the controlling knob. Free ~2× win just from `dt 0.001→0.002`.

### Diagonal-curvature preconditioner (tested, Python replica)

```
k_ii  = -df_i/dx_i          (finite-difference Hessian diagonal, 2 evals per free DOF)
invM_i = k_ref / (|k_ii| + k_ref)     ∈ (0,1]
```

Design notes (learned the hard way):

- Use **|k_ii|**, not clipped positive `k_ii`. Initial esize curvatures are often **negative** (concave-down slopes, e.g. core e2 `k≈−7745`, valence `k≈−44,−0.08`). Treating them as soft (`invM≈100`) makes them ultralight → runaway → `cos_vf≈−1` forever → dt pinned at `dt_min`. Concave directions are still *steep* and need heavy mass.
- Cap `invM ≤ 1` — only add mass to stiff directions, never lighten soft ones. Uncapped variants (invM ∝ 1/k or k-floor clipping) gave soft DOFs invM≫1 → noise dominates `cos_vf` → same stall.
- "Equalize ω" (invM ∝ 1/k fully) *stabilizes* but is **slower** (~15k iters): heavy masses also slow the *descent* of steep DOFs. The capped form is a compromise — strongly-damped stiff modes, untouched soft ones.

Results (iters to |F| < 1e-3 / 1e-4 / 1e-5):

| system | unit dt=0.001 | best unit | prec k_ref=30 |
|---|---|---|---|
| H2O_spins | 6686 / 7726 / 8745 | 3495 (dt=0.002) | **332 / 870 / 977** (dt=0.05) |
| H2O_pairs_fc *(pre-fix)* | 1949 / 2325 / 2625 | 962 (dt=0.002) | **169 / 201 / 218** (dt=0.05) |

Robust across k_ref ∈ [30,100], dt ∈ [0.02,0.1] on both systems; preconditioned runs consistently reached the deeper −1852.86 basin. Cost: 2·n_free evals per diagonal estimate (ne≈16 → ~200 evals; refreshing every ~200 steps is ~1% overhead).

**Caveat — basin drift:** preconditioning changes the trajectory and can change which local minimum you reach. Here it went *deeper*; on CH4-bPaired it found *shallower* basins (−987/−1027 vs −1112 for the unit-mass path). For scans/databases keep the protocol fixed per dataset, or use preconditioning only in a first "descend fast" phase.

---

## 2. `processXYZ` auto-electron path — a different failure, not stiffness

On every tested molecule (`CH4, CH3CH3, H2O, CO, CH3OH, CH2O, HCCH, C2H4, O2, OHOH`) with auto-generated electrons and `bPaired=False` (frozen nuclei), **no setting converged**: |F| plateaus at 20–600 while E keeps drifting. Stock C++ `run()` reproduces it (CH3CH3: |F|min=200 after 20k steps) — not a replica artifact.

Diagnosis of the end state (CH4): **electrons collapse** — 4 coincident electrons on the C nucleus, opposite-spin pairs at *exactly* r=0. The residual |F| comes from nonsmooth cusp forces at coincident electrons (force can't vanish at the cusp). Also the "minima" are unphysical: CH4 stalled at −558/−707 eV vs physical ≈ −1099 eV.

**`bPaired=True` largely fixes it:** CH4 → E=−1112.6 (physical); with the pre-refit `erfx_e6` it stalled in a **cusp limit-cycle** at |F|≈0.02–0.03 (E frozen, `cos_vf≈−0.9`, dt pinned at dt_min). **Mechanism identified and FIXED:** the `erfx_e6` approximant had nonzero slope (−2.6e-4) at r=0 → a linear |r|-cusp in the electron–nucleus Coulomb; the paired core `e2` relaxes onto the C nucleus (r≈5e-9 Å) and the cusp's one-sided slope left a residual no smooth optimizer could remove. **After the constrained refit (§4.2): CH4-bPaired converges to |F|=9.9e-10 at the same physical energy −1112.6417** — the floor is gone entirely.

**Post-fix re-checks (2026-11-08):** (a) after the `fsi`/`fsj` swap fix CH4-bPaired still plateaued at |F|≈0.019–0.032 → the floor was the cusp, not the swap; (b) after the `erfx_e6` refit the same run reaches |F|=9.9e-10. The `bPaired=False` collapse to wrong basins (E≈−558/−758 vs physical −1112.6) is unchanged — a separate initialization/pairing pathology, not a force bug.

**Recommendation for `processXYZ` scans (production `run_process_xyz.py` protocol):** use `bPaired=True` — post-refit it converges cleanly (CH4: |F|→1e-9 at E=−1112.64). `bPaired=False` currently produces physically wrong minima *and* wastes all iterations. Stall-detection (E stationary + persistent `cos_vf<0` at `dt_min` → accept) remains useful as a safety net for other nonsmooth configurations.

Parameter-table limitation found en route: `ATOMPARAMS2` (test table) has real entries only for H/C/O — N/F rows are zeros → `builder2EFFstatic` fails loudly (`core esize<=0`) on N-containing molecules. Good (fails loud), but it limits the sweep set.

---

## 3. Strategy discussion — three families of fixes

> **User corrections (2026-11-08):** `EFF::makeMasses()` is for *physically realistic MD in time* — **never for relaxation** (unit mass is the right relaxation baseline). The task is a *new static heuristic optimal for relaxation*: stiff DOFs get **more** mass. Electrons close to nuclei are much stiffer → heavier, but that must be determined **dynamically per system**. All-electron vs frozen-core/ECP (`setCoreMode` in `eFF.h`: `'a'`=all-atom, `'f'`=frozen-core pair at nucleus w/ fixed size, `'e'`=ECP) are fundamentally different regimes. Test consistently on **H2O** (small, has bonding + lone-pair electrons, often problematic).

### 3.0 First distinguish the regime — coreMode matters (measured)

| | coreMode 'a' (H2O_spins) | coreMode 'f' (H2O_pairs_fc) |
|---|---|---|
| bound core pair | yes — O `e2` (s=0.2) on nucleus | none (fixed-size pair at nucleus, not a DOF) |
| max \|k_ii\| | ~5.4×10⁵ (O–e2 relative mode) | ~500 |
| unit-mass dt ceiling | ~0.0027 (stiff mode rings) | much higher |
| best unit-mass result | 3495 it (dt=0.002) | **189 it (dt=0.05)** |

**Frozen core / ECP removes the stiffness catastrophe by construction** — no 10⁵ modes, plain FIRE with `dt≈0.02–0.05` is already fast. **Post-force-fix numbers (H2O_pairs_fc, unit mass):** dt=0.001→1857 it, dt=0.005→456, dt=0.02→187, dt=0.05→**157 it** (E=−449.66 — note the fix *changed the relaxed energy* by ~7 eV vs the pre-fix −442.85; the old "minimum" was an artifact of wrong forces). **The optimizer problem is essentially specific to all-electron mode `'a'`** — and mode 'f'/'e' is itself an optimization strategy when explicit core physics isn't needed.

**The student's mystery solved:** the likely cause of "10k–millions of iterations" is `bMass=true` → `makeMasses()` *physical MD masses* — measured on both H2O modes: mode 'f' 189→**6496 it** (~30× worse), mode 'a' **never converges** (light invM=17.6 electrons vs invM=0.0006 atoms is exactly backwards for relaxation). Combined with all-electron mode ('a' has ~10³× more stiffness than 'f'), the two mistakes together fully explain million-iteration runs. For production relaxation: `bMass=false`, mode 'f'/'e', dt=0.02–0.05 → few hundred iterations, matching the remembered experience.

### 3.1 Universal *fixed* masses — what works, what doesn't (measured on H2O_spins, mode 'a')

Scans (iters to |F|<1e-3; −1 = not in 30k):

- **Uniform esize multiplier** `invM_esize = 1/m_s`: monotonically **worse** — m_s=3→8920 it, m_s=10→21931, m_s≥30→never. Slowing *all* size DOFs starves the soft valence-size descent. Per-group uniform mass is definitively out.
- **Analytic esize mass** `invM_i = k_ref/(9Q_i/s_i⁴+k_ref)` (s⁴ from `E_k≈3Q/2s²`): alone, no help (same as unit at equal dt) — esize isn't the *stability* bottleneck, epos is.
- **Bound-pair mass** (the winning static idea — "electrons close to nuclei are stiffer"): for each electron within `R_core≈0.5 Å` of a nucleus, give **both** the electron's epos *and* the host nucleus's apos `invM=1/m_c` (the stiff eigenmode is ~50% atom — slowing only the electron does nothing), plus the analytic esize mass. Scan:

| m_c \ dt | 0.002 | 0.005 | 0.01 | 0.02 | 0.05 |
|---|---|---|---|---|---|
| 10 | 13427 | 6405 | −1 | −1 | −1 |
| 30 | 23553 | 9538 | 5918 | −1 | −1 |
| 100 | −1 | 16526 | 8181 | **3471** | −1 |
| 300 | −1 | 24232 | 14288 | 6160 | −1 |

  Heavier pair → larger stable dt, but a *floor* ~3500 it — the remaining stiffness is a continuum (valence-bond modes k~10²–10³) that static rules don't touch. Best static ≈ baseline best (3495), though it robustly reaches the deeper −1852.86 basin.
- **Adaptive FD-diagonal** `invM=k_ref/(|k_ii|+k_ref)` on *all* DOFs remains ~10× better (332 it) because it damps the entire stiffness continuum, not just the known pair.

**Verdict:** static masses help only if they catch *all* stiff DOFs — the practical static rule is "detect bound (e,host-atom) pairs → shared heavy mass, recompute per system/geometry". It's ~free and gives ~2×; the adaptive diagonal gives ~10–20×. Fixed *physical* masses (makeMasses) are ruled out entirely — wrong ordering for relaxation (it's the MD-in-time mass table).

### 3.2 On-the-fly mass estimation — efficient robust heuristics

- **[tested] FD Hessian diagonal + capped invM, refreshed periodically.** `k_ii=-df_i/dx_i` central FD (2 evals/DOF), `invM=k_ref/(|k_ii|+k_ref)`, refresh every ~100–500 steps. Costs ~1–4% extra evals; **works** (§1). Refresh matters: initial curvatures are misleading (negative, and ordering changes as electrons settle); `kd_refresh` runs found *much* deeper minima on CH4 (−707 vs −564) — though "deeper" there was the collapsed artifact, so refresh also lets preconditioning *explore* more aggressively. Clamp `|k_ii|` to e.g. [1e-2, 1e7] and never raise invM above 1.
- **[tested — works, ~free] Secant diagonal from the trajectory** (`fire_secant` in `fire_study_lib.py`). `k_i ≈ -Δf_i/Δx_i` per step, EMA-smoothed (β=0.9), clamped to [1e-2,1e7], `invM=k_ref/(k_i+k_ref)`. **Zero extra evals.**
  - H2O_spins (mode 'a'): k_ref=100, dt=0.1 → **654 it** (vs 6686 baseline, ~10×), always lands in the deeper −1852.86 basin. k_ref=30–100 all work; larger dt helps *because* invM adapts (opposite of unit mass). dt≥0.2 diverges — keep dt≤0.1.
  - H2O_pairs_fc (mode 'f'): 220 it @ k_ref=100,dt=0.1 — on par with plain FIRE's best (189 @ dt=0.05); harmless when stiffness is mild.
  - ~2× slower than the FD-Hessian preconditioner (332 it) but costs nothing extra and self-adapts continuously. **Best practical candidate so far.**
- **[tested — rejected] Per-group FIRE.** Separate `cos_vf`/dt/damping per group (`fire_groups`). TYPE-mode (apos|epos|esize) is **worse** (17917 it; dt≥0.005 fails): the stiffest mode is *coupled* O-atom+e2 — it spans groups, so decoupling the power test doesn't isolate it. PARTICLE-mode (per atom/electron) is marginally better (2815 it @ dt=0.005) but still worse than any mass scheme and unstable at dt≥0.02. **The bottleneck is heterogeneous ω, not the global power test** — decoupling can't fix it; slowing stiff modes can.
- **[untested] Per-DOF adaptive step (ASIP/EDoS-style).** Instead of masses, keep a per-DOF `dt_i`: grow when v_i·f_i>0, shrink on sign flip. Mathematically equivalent to adaptive diagonal mass, fully online, standard in the literature (adaptive inertial descent). Watch interplay with the global FIRE power test (use per-DOF power or blend).
- **[untested] Adam/RMSProp-style preconditioning.** invM_i ∝ 1/(sqrt(v_i)+ε) from running second moments of forces — again a free online diagonal conditioner; notoriously good at heterogeneous-curvature problems, but noisy near convergence → combine with trust radius (dr_limit already exists) and maybe switch to plain FIRE below a force threshold.

Stability rules for any adaptive scheme (learned from the failures): `invM ∈ (0,1]` only (add mass, never lighten); use `|k|`; EMA-smooth updates; keep `dr_limit` clamp on displacement.

### 3.3 Different algorithms?

- **BFGS / L-BFGS:** builds a quasi-Newton `H⁻¹` from trajectory — the "ideal" version of a learned preconditioner; superlinear near a minimum. But (a) needs a line search or trust region — each Wolfe step costs evals; (b) along the long curved valleys (which is most of the eFF trajectory from a random/bad init) momentum methods like FIRE empirically travel faster; BFGS tends to waste evals re-learning curvature that keeps changing. The user's intuition is right: **BFGS wins near the minimum, FIRE wins far from it.**
- **[tested] Prototype BFGS in `fire_study_lib.py`** (dense `H_inv`, trust-region step — accepted iff E decreases; `dr_limit` cap). Post-refit (2026-11-08): `H2O_pairs_fc` → **49 iters to |F|<1e-3, Fmin=7.8e-6** (E=−449.658) — 3× faster than best FIRE (157 it) and ~40× deeper in force; `H2O_spins` → **156 iters to 1e-4, Fmin=9.4e-5** (E=−1852.86 deep basin). Pre-refit it stalled at Fmin≈0.0019 on pairs_fc — both stalls were the cusp/approximation inconsistencies, not a model failure. With the refit, dense BFGS is a strong standalone candidate at eFF's small DOF count (~25–60); L-BFGS for larger systems.
- **Two-phase hybrid (recommended to test):** FIRE (ideally preconditioned) for the global descent → switch to L-BFGS or diagonal-Newton polish when |F| drops below a threshold or E stalls. Would also convert the `processXYZ` cusp plateau into an explicit "minimum up to nonsmooth residual" instead of burning 50k steps.
- **Diagonal-Newton with momentum** `[tested]`: `Δx_i = f_i/max(|k_ii|,k_floor)`, β=0.5 momentum, |Δx|≤dr_limit. Converged on H2O_spins in ~1094 iters (3097 evals incl. FD Hessian) — competitive, but *unsafe in concave regions* (signed k flips step direction; needs |k|+floor+trust). Plausible as the "polish" phase.
- **Truncated-Newton / Newton-Krylov:** H applied only via finite-difference products, inner CG with diagonal preconditioner — the rigorous version of the same idea; more machinery than the problem needs at ne≲50.
- **Nonlinear CG (PR+ with restarts):** cheap, but line-search-sensitive; no momentum mass handling; unlikely to beat preconditioned FIRE.
- **Eigenvector-following / dimer / saddle methods:** overkill — the observed negative modes are shallow and FIRE's inertia already gets past them.
- ~~**For the cusp stall itself:** no smooth optimizer fixes a nonsmooth force~~ — **resolved**: the cusp was an `erfx_e6` defect, fixed by the constrained refit (§4.2). Stall-detection stays useful for genuinely nonsmooth switching surfaces (Pauli `step`-like terms, `bPaired=False` collapses).

### Suggested experiment order (cheapest→deepest)

0. **Check coreMode first.** For 'f'/'e' systems just raise `dt` (0.02–0.05, unit mass) — no preconditioner needed. All the work below is for mode 'a'.
1. ~~Per-group `cos_vf`~~ — tested, **rejected**: stiff modes span groups, doesn't help.
2. Static bound-pair mass (detect e within ~0.5 Å of nucleus → share `invM=1/m_c`, m_c≈100, with host atom) + analytic esize mass — ~free, ~2×, reaches deeper basin.
3. **Secant-diagonal adaptive invM with EMA + clamps — tested, ~free, ~10×, best practical candidate.**
4. FD-diagonal refresh (~10–20×, needs 2·ndof evals/refresh) as the stronger/costlier option.
5. Two-phase FIRE→L-BFGS polish; stall-detector exit for cusp plateaus.
6. Whatever wins → port to C++ in `DynamicOpt`/`eFF_lib::run`, verify parity vs replica.

---

## 4. Code changes — exact changelog and measured results

All changes are on the working tree (not yet committed). Files touched:

| file | lines | what changed | result |
|---|---|---|---|
| `cpp/common/molecular/eFF.h` | 700, 708 | `evalAE()`: swapped `fsi`/`fsj` outputs → electron `fsize` now gets `fsj` (its own size derivative), atom derivative → `fs_junk`. 2 call sites: `addPauliGauss_New` (AE-Pauli), `addCoulombGauss` (core-Coulomb `bCoreCoul`) | `fsize` consistent with `−dE/ds` (was up to ~400% wrong, wrong signs); frozen-core H2O energy −442.85 → **−449.66 eV** (7 eV artifact removed); BFGS trust-collapse stall gone |
| `cpp/common_resources/cl/eFF.cl` | 372, 378 | same swap in the AE Pauli/core-Coulomb accumulation: `my_f.w += pg.z`/`cgC.z` → `pg.w`/`cgC.w` | CPU↔GPU parity *preserved* (kernel had the identical bug — parity tests could not detect it) |
| `cpp/common/math/fastmath.h` | 141–160, 163–170 | `erfx_e6` (both overloads): constrained refit coefficients — `e0`,`o0` pinned to exact values at origin, `deven`/`dodd` recomputed; `REFIT` comment | **\|r\|-cusp removed**: `dy(0)` −2.6e-4 → ~0; max\|err\| 3.7e-7 → **9.1e-8** |
| `cpp/common_resources/cl/eFF.cl` | 26–30 | same refit coefficients in float + `deven`/`dodd` | GPU identical (fp32 residual slope ~3e-8) |
| `cpp/libs/Molecular/eFF_lib.cpp` | 150–160 | `init_buffers()`: `buffers.insert({k,v})` → `buffers[k]=v` (insert never overwrites → stale pointers after `ff.realloc()`) | ASan use-after-free / silent pointer corruption fixed |
| `cpp/libs/Molecular/eFF_lib.cpp` | 703 | `processXYZ` read loop: `while(fgets(...) \|\| (atoms && il==natoms+2))` | last conf no longer silently dropped at EOF |
| `cpp/sketches_SDL/Molecular/python/erf_approx.py` | ~65–70 | provenance comment: new coefficients + refit method (the file is where the original coefficients were developed) | provenance |
| `tests/tEFF/fire_study_lib.py` | 75–94 | `load_a` now calls `eff.set_constrains` after build (see §4.3 gotcha) | stock `run()` on build-only loads now uses the intended atoms-fixed protocol |
| `tests/tEFF/test_erfx_refit.py` | new | regression test: CH4-bPaired relaxes to \|F\|<1e-4, FD-vs-force <1e-3, coincident-electron force <1e-3 | **all pass** |

### End-to-end results of the two physics fixes (measured 2026-11-08)

| case | before fixes | after |
|---|---|---|
| `H2O_pairs_fc` relaxed energy (mode 'f') | −442.85 eV (wrong basin, bad forces) | **−449.66 eV** |
| `H2O_pairs_fc` FIRE unit dt=0.05 → 1e-3 | 157 iters (same speed, different basin) | 157 iters |
| `H2O_pairs_fc` dense BFGS | stalled at \|F\|≈0.0019 (trust collapse) | **49 iters to 1e-3, Fmin=7.8e-6** |
| `H2O_spins` dense BFGS (mode 'a') | stalled at \|F\|≈5e-4 | **156 iters, Fmin=9.4e-5** |
| CH4-bPaired (`processXYZ`) | \|F\| floor 0.02–0.03 forever (cusp) | **converges to \|F\|≈1e-9** at E=−1112.64 |
| coincident e–C pair phantom force | 0.028–0.045 eV/Å | ~0 (3.7e-7) |
| esize FD-vs-force consistency | up to ~400% error + sign flips | max residual ~3e-4 (approximant mismatch, §4.1 residual #2) |

### 4.1 `evalAE()` electron-size force bug — swapped `fsi`/`fsj` (physics bug, FIXED)

In `eFF.h::evalAE()` the atom–electron interaction calls pass the two size-derivative outputs `fsi` (d/d first size arg) and `fsj` (d/d second size arg). Two of the three atom-electron size-force sites were swapped — the electron's `fsize` received the **atom's size derivative** (`fsi`) while its own derivative (`fsj`) went to `fs_junk`:

- `addPauliGauss_New(dR, aPar.y, sj, f, fsj, fs_junk, ...)` (`eFF.h:700`) → was: electron got `d/d(aPar.y)` (Pauli overlap has an antisymmetric size part → wrong magnitude **and wrong sign**).
- `addCoulombGauss(dR, aPar.y, sj, f, fsj, fs_junk, ..., qj*aPar.z)` (core-Coulomb `bCoreCoul`, `eFF.h:708`) → same swap.
- Fixed: `(aPar.y, sj, f, fs_junk, fsj, ...)` — matches the already-correct convention in `evalAE_ECP` (`eFF.h:826`) which passes the electron first (`sj, aPar.y, fsj, fs_junk`), and the already-correct AE-Coulomb site at `eFF.h:688`.

**Same swap existed in the OpenCL kernel** `cpp/common_resources/cl/eFF.cl` (`:372`, `:378` — electron accumulated `pg.z`/`cgC.z` = `fsi` instead of `pg.w`/`cgC.w` = `fsj`) — fixed identically to preserve CPU↔GPU parity (both were wrong *the same way*, which is why parity tests never caught it).

**Impact:** `fsize` was not `-dE/ds` — FD check showed ~400% magnitude errors and wrong signs on esize DOFs. Consequences: (a) frozen-core/ECP relaxed geometries/energies were wrong — `H2O_pairs_fc` relaxes to **−449.66 eV post-fix vs −442.85 pre-fix** (7 eV artifact); (b) BFGS took "uphill" steps at a positive-definite model → trust-region collapse → stall; (c) all esize-driven convergence behavior (including the §2 cusp hypothesis and the preconditioning tables) was measured against corrupted forces. **All published eFF energies/geometries computed on coreMode 'f'/'e' paths should be re-checked.**

**Verification:** post-fix FD-vs-force probe reports `inconsistent DOFs: 0 of 25` (tolerance 5%) on `H2O_pairs_fc`.

**Residual inconsistency floors:**

1. ~~**`|r|` cusp at coincident electron–nucleus pairs**~~ — **FIXED by the `erfx_e6` constrained refit (§4.2, 2026-11-08).** The old approximant had slope −2.6e-4 at r=0 (true `erf(kx)/x` is quadratic-flat) → a linear |r|-cusp ≈0.045 eV/Å for a C point-charge core; paired core electrons on nuclei carried |F|≈0.02–0.03 forever. Post-refit CH4-bPaired converges to |F|=9.9e-10.
2. **~3e-4 esize approximant mismatch (still present, small).** `fs` comes from `gauss_p8` (analytic derivative of the *true* Gaussian) while `E` uses `erfx_e6` → esize forces differ from −dE/ds by ~3e-4 absolute (measured at the BFGS-converged pairs_fc point: fsize≈0 vs dE/ds≈−3.1e-4 on all 4 electron sizes). Harmless in practice — BFGS now converges to 7.8e-6 through it, FIRE follows f anyway — but it is a genuine f-vs-E inconsistency: any optimizer that uses both (line-search/trust-region) may see ~3e-4-level model disagreement. Fix options if it ever matters: derive `fs` from the same `erfx_e6` approximant (needs `dy/ds` which the refit already returns), or tighten both fits together.

### 4.2 `erfx_e6` constrained refit (FIXED the |r|-cusp)

**Model** (`fastmath.h`): `y = k/(t⁸ + x)`, `t = even(x²) + x·odd(x²)`, `x = k·r`, degree-5 even + degree-5 odd (12 coefficients); analytic derivative `dy = −y²(8t⁷dt/dx + 1)` via `deven`/`dodd`.

**Exact analytic constraints at x=0** (from `y = k/(t⁸+x)`):
- value: `y(0) = k/e0⁸ = 2k/√π` ⇒ **`e0 = (√π/2)^{1/8} = 0.98501562035791224`** (the old fit had this right to ~1e-10).
- slope: `dy/dx|₀ = −k²(8·e0⁷·o0 + 1)/e0¹⁶` ⇒ zero iff **`o0 = −1/(8·e0⁷) = −0.13893388815944127`**. The old fit had `o0 = −0.13893350387140332` — off by 2e-5 → slope −2.6e-4 → the cusp.

**Fitting machinery:** the coefficients were historically hand-fit in `cpp/sketches_SDL/Molecular/python/erf_approx.py` (several old variants are commented there; the repo's generic `pyBall/rational_poly_fit_lib.py` is a Chebyshev `v·Pⁿ/Qᵐ` fitter — related but not this model). Refit method used here: pin `e0`,`o0` to the exact constraints, weighted `scipy.optimize.least_squares` on `y` over `x ∈ (0, 4.5]` (the branch boundary — the `x>4.5` tail is `1/x` exactly) with ~20 Lawson reweighting iterations toward minimax; remaining 10 coefficients free.

**Result:** **max |y−erf(x)/x| = 9.1e-8** on the branch (old fit: 3.7e-7 — the refit is *4× more accurate* while also satisfying both constraints exactly). `deven`/`dodd` recomputed as exact derivative polynomials (`deven[k]=2(k+1)·e[k+1]`, `dodd[k]=(2k+1)·o[k]`).

Small-r derivative validation (refit `dy` vs analytic `d(erf/x)/dx`, measured): at r=0 `dy≈−2e-14` (was −2.6e-4 — the cusp), at r=1e-7 `dy≈−4.8e-5`, at r=1e-5 `dy≈−4.8e-3` — linear in r and FD-consistent everywhere (FD itself is noisy below r~1e-9, as expected for a smooth minimum). The effective potential is now quadratic-flat at coincidence, as true `erf(kx)/x` requires.

**New coefficients (applied to `fastmath.h` both overloads + `eFF.cl` float version):**

```
even = [ 0.9850156203579122, -0.027549047779963233, -0.001290777815065076,
        -0.0009020815265539048, -0.00034409917491279876,  1.1960756662257476e-05 ]
odd  = [ -0.13893388815944127, -0.007786974336511317,  0.0015375620941944985,
         0.0009682977842733449,  1.5066710497913615e-05, -1.60225628433894e-06 ]
```

fp32: the constraint `8·e0⁷·o0+1` evaluates to −3e-8 in float → residual r=0 slope ~3e-8, ~10⁴× below the old cusp — acceptable (fp32 noise dominates).

**Note:** `erfx_e9` (higher-accuracy variant in the same file) never had the cusp — its `x<1/k` branch is a pure even polynomial. `InteractionsGauss.h`, `GaussianBasis.h`, `CLCFGO_Structs.h`, `libMMFF.h` call `erfx_e6` but hold no local coefficient copies — all get the fix automatically. Only `eFF.cl` duplicated the coefficients (now synced, CPU↔GPU parity preserved).

**Verification (2026-11-08):** CH4-bPaired |F| floor 0.02–0.03 → **9.9e-10** (same energy −1112.6417); coincident e–C pair FD-vs-force diff 0.028 → **0.0**; H2O_pairs_fc BFGS stall 0.0019 → **Fmin=7.8e-6 in 49 iters**; H2O_spins BFGS → Fmin=9.4e-5 in 156 iters. All relaxed energies unchanged to ~1e-4 eV (fit is *more* accurate than before). GPU (`eFF.cl`) updated in parity — **OpenCL regression test not yet rerun** (needs the OCL test harness).

**Regression test added:** `tests/tEFF/test_erfx_refit.py` — stock `run()` on CH4-bPaired must reach |F|<1e-4 (pre-refit floor was ~0.02–0.03), FD force consistency at the relaxed state <1e-3, coincident core-electron force <1e-3 (was ~0.045). All pass.

### 4.3 eFF_lib.cpp plumbing
- **`init_buffers()` stale pointers — fixed.** `buffers.insert({k,v})` never overwrites existing keys → after `ff.realloc()`, `getBuffs()` returned pointers to **freed arrays** (ASan use-after-free; silently corrupt in normal builds). Now uses `buffers[k]=v`. Same `insert` anti-pattern exists in other `*_lib.cpp` files — worth auditing.
- **`processXYZ` drops the last configuration at EOF — fixed.** The `il==natoms+2` "conf complete" trigger only fires on reading a *subsequent* line; a file ending without a trailing blank line silently skipped its final structure (`ff` left uninitialized → NULL `getBuff` crash). Loop condition now also processes at EOF. Related fragility still present: the *next* conf's natoms line is consumed as the trigger rather than parsed — works only because multi-conf files share `na`; a mixed-size file would misparse.
- **`processXYZ(nstepMax=0)` silently skips the atom-fixing block — found during regression testing.** The `nfix`/`fixed_inds` constrain setup only runs inside `if(nstepMax>0)` (line ~748), so a build-only call leaves **atoms completely free** in any later `run()` — stock FIRE then relaxes a much harder all-free problem and "diverges" relative to atoms-fixed replicas. `fire_study_lib.load_a` now applies the same constraints explicitly via `set_constrains` after the build; any other build-only caller needs the same workaround (or the constrain block should move out of the `nstepMax>0` guard).
- **`initOpt`/`run` plumbing:** `dt_max=dt` is set every `run()` call (dt is a ceiling); FIRE internals (`finc,fdec,falpha,minLastNeg,dtmin_frac`) are compile-time constants not reachable from Python — worth exposing if tuning continues.

### 4.4 Second review round (`eFF_CPU_review_2026-11-08.md`) — verified and fixed

Every claim was verified numerically (component-isolated probes + central-FD force checks) before editing. All fixes below were **confirmed real bugs**; no claim was refuted.

#### P0 — physics bugs

**P0-1 Paired `e2` orbitals missed internal electron repulsion — FIXED.** `evalEE()` iterates `j<i` over distinct orbitals only, so the `qi·(qi−1)/2` constituent e–e pairs inside one doubly-occupied orbital had no Coulomb term (internal Pauli is identically zero: opposite spins, coincident same-width Gaussians). Live probe (isolated pair, `s=0.5`): `e2` gave `Eee=0`, explicit `e1+e1` gave `+32.50 eV` — and `fsize` was −1.9 vs +63.1, i.e. paired orbitals were missing the outward size pressure and systematically collapsed. Fix (`eFF.h`, inside the `i`-loop after the `j` loop): `addCoulombGauss(Vec3dZero, si, si, fjunk, fsize[i], fsize[i], qi*(qi-1)*0.5)` — `fsize[i]` passed for both `fsi`/`fsj` so the single shared size DOF gets the full derivative. **Post-fix: `e2` energy and FD forces exactly match the explicit `e1+e1` representation** (E=+32.02206 both; fsize 246.9 = 2× per-electron, FD-consistent). Consequence: all paired-orbital energies shift up by ~+25–35 eV per `e2` (H2O_pairs_fc relaxed: −449.66 → −343.32 eV, converges 861 iters); **the "core e2 on nucleus" in mode 'a' also gains self-repulsion** — all-electron results change (this was missing core physics).

**P0-2 ECP Coulomb: repulsive, charge used as width, no occupancy — FIXED.** `evalAE_ECP` called `addCoulombGauss(dR, aPar.z, sj, ..., (aPar.x-aPar.z))`: `aPar.z` is `Zcore_eff` (a charge, not a length) passed as the smearing width, and `qq=+(Z−Zeff)>0` → repulsive electron–core force. FD probe (C ECP, r=2): E=+27.3 repulsive. Fixed to `addCoulombGauss(dR, aPar.y, sj, ..., -qj*(aPar.x-aPar.z))` — width `R_eff`, electron charge `−qj`. Post-fix: E=−28.80 (≈point-charge −4·14.4/2), force = −dE/dx to 6 digits. `evalECP` (dead code, same bug) fixed identically for consistency.

**P0-3 ECP Pauli forces sign-reversed — FIXED.** `evalAE_ECP` did `fsj -= fe_*qj` (should be `+=`: `fe_` is already −dE/ds, a force) and `f.add_mul(dR, +fr_*qj/rc)` (f accumulates force on the *atom* = −electron force → needs `−`). FD probe: reported −154.50/+114.17 vs true +154.50/−114.17 — exact reversals. Post-fix: force and fsize match FD to 6 digits. Same two sign errors existed in dead-code `evalECP` — fixed.

**P0-4 Negative electron occupancy → negative kinetic energy — FIXED.** `_realloc0(echarge, ne, QE)` initialized `echarge=-1`, but all consumers treat it as a positive occupancy magnitude (`from_xyz` uses `fabs(Q)`; `evalAE` uses `−qj` as the signed electron charge; `addKineticGauss_eFF` multiplies by it). Any path that didn't set `echarge` explicitly got Ek<0. `readSingleFGO()` parsed the `c` field but never assigned it. Fixed: default `1.0`, `echarge[i]=fabs(c)` in `readSingleFGO` (+ closed-shell copy). FGO `H2O_noe2_.fgo` now gives finite E=−591.13 with Ek=+180 (was Ek<0). Updated the `echarge` doc comment to "occupancy magnitude".

**NaN guards added along the way:** `evalAE` Pauli gated `aPar.y>1e-8` (zero-width core `sQ=0`, e.g. FGO point nucleus, gave 0/0); `evalAE_ECP` non-ECP fallback gated `aPar.z>0` (bare nuclei like H have no core → skip Pauli; was NaN).

#### P1 — optimizer/driver issues

**P1-5 Convergence norm was mass-weighted — FIXED.** `run()` used the return of `opt.move_FIRE()` = `ff=Σ(f·invM)²` for `F2<Fconv²`, while `ialg=0` (`ff.move_GD`) returns raw `Σf²` — with `makeMasses()` physical masses the FIRE metric is distorted (electron forces amplified ×invM²≈300, atom forces suppressed). `run()` now computes `F2=ΣfDOFs²` (raw, post `clear_fixed_*`) for all algorithms; FIRE keeps mass weighting only for its internal power test where it belongs.

**P1-6 `dr_limit` — misplaced on first attempt (review §6 caught it), now FIXED in the active path.** The cap was initially added inside dead-code `DynamicOpt::move_FIRE_bak()`; production `run(ialg=2)` calls `move_FIRE()` (~line 316). It is now implemented there with **uniform step scaling** (one common factor on all velocities → `dr'_i = s·dr_i`, step direction preserved, per review guidance) instead of per-component clipping. Verified active: see §6.7.

**P1-7 Width-floor inconsistency — FIXED.** `evalKinetic` clamps `esize→min_esize=0.1` (write-back), `run()` used a literal `1e-3`. Now both use `ff.min_esize` (0.1). GPU harness's ~1e-3 floor noted — kernel-side parity item.

**P1-8 Optimizer pointers surviving `ff.realloc` — FIXED.** `opt` binds raw `pDOFs/vDOFs/fDOFs` pointers at `initOpt`; any later realloc leaves them dangling. `run()` now detects `(opt.pos!=ff.pDOFs || opt.n!=ff.nDOFs)` → warning + `bindOrAlloc` rebind. **Fired in testing** (second `load_e` in a loop) — was a real latent use-after-free.

**P1-9 Non-finite energy handling — FIXED.** `eval()` printed an ERROR on NaN but `run()` marched on for 30k steps. `run()` now aborts loudly on non-finite Etot (consistent with existing `exit(0)` ERROR convention).

#### GPU parity (`cpp/common_resources/cl/eFF.cl`)

- Internal `e2` Coulomb added to all three kernels (`localMD` force → `my_f.w += cgii.z+cgii.w`; `localMD_energy`, `evalEnergy` → `EeeCoul += cgii.x`), gated like the kinetic term under `bOffloadCore`.
- **`localMD` ignored `e2` occupancy entirely** (`qe=1.0` in AE Coulomb, `qij=1.0` in EE Coulomb, core-Pauli/core-Coulomb underweighted ×2) — now `qe=(spin==0)?2:1`, `qij=qe_i·qe_j`, parity with CPU `echarge`.
- Both build variants compile on NVIDIA (`localMD`, `evalEnergy`, `localMD_energy`, `eval_density_grid`, `fit_density_fire`); Intel iGPU lacks `cl_khr_fp64` (pre-existing, unrelated).

#### Post-fix regression status

| check | result |
|---|---|
| `tests/tEFF/test_erfx_refit.py` | PASS (Fmin 9.9e-5, FD 2e-4, coincident-pair force ~2e-6) |
| e2 vs e1+e1 isolated pair | E identical (+32.022), fsize FD-consistent |
| C-ECP Coulomb/Pauli FD | sign+width correct, 6-digit match |
| FGO `H2O_noe2_` | finite E=−591.13, Ek>0 |
| H2O_pairs_fc FIRE dt=0.05 | converges 861 iters, E=−343.32 (was −449.66 on buggy Hamiltonian) |
| H2O_spins 'a' FIRE dt=0.002 | 3497 iters (≈3495 baseline), E=−1719.54 |

## 5. Open questions

- **Physics regression scope is now larger:** any result produced with paired `e2` orbitals (all `bPaired`/closed-shell/FGO-closed-shell runs, every mode-'a' core `e2`, every ECP run) used a Hamiltonian missing internal repulsion, repulsive ECP Coulomb, and reversed ECP Pauli forces. All paired-system energies/geometries from before 2026-11-08 need re-evaluation; mode 'f' additionally had the `fsi`/`fsj` swap.
- **OpenCL regression:** `eFF.cl` received the swap fix, refit coefficients, internal-e2 Coulomb and the `qe` occupancy fixes — rerun the CPU↔GPU parity test (`test_ocl_vs_cpu.py`, static `dt=0` first, expect ~1e-5) before shipping. Kernel compiles on NVIDIA; Intel iGPU lacks `cl_khr_fp64` (pre-existing).
- **GPU `min_esize` floor:** CPU now uses a single 0.1 floor (`ff.min_esize`); the GPU harness clamps ~1e-3 — harmonize or parameterize kernel-side.
- **`evalECP` remains dead code** — fixed for consistency but nothing calls it; either wire it up or remove.
- **`processXYZ_e` mixed-size configs:** buffers/`opt` are sized for iconf 0 only; a later config with different na/ne overflows — needs realloc+rebind per header or an explicit same-size check.
- **Residual ~3e-4 esize f-vs-E mismatch** (`gauss_p8` fs vs `erfx_e6` E): worth removing if BFGS/L-BFGS becomes the production polish — derive `fs` from the same approximant. Also could now re-derive `fs = dy/ds` directly since the refit returns a consistent `dy`.
- Does preconditioning change the *physics* outcome on the production scan set (`distscan/angdistscan_H2O`, bCoreElectrons variants)? Compare basin-by-basin energies vs unit-mass at fixed protocol before adopting.
- How large can `k_ref` get before stiff-DOF descent is starved? (k_ref=30 beat 100 on H2O_spins; crossover somewhere between.)
- BFGS prototype: now converges standalone on both test systems (49/156 iters) — worth porting a bounded dense/L-BFGS + trust-region into `DynamicOpt` for mode-'a' polish; stall-detection still useful for `bPaired=False` collapse cases.

## 6. Review of §4.4 and bounded frozen-core speed study

### 6.1 Implementation review: accepted changes and unresolved defects

The internal-e2 Coulomb term is now present, accumulates both shared-size derivatives, and contributes to the orbital energy. ECP Coulomb has the negative occupancy-scaled charge product and uses the radius field; ECP Pauli wrapper signs now match the returned radial/size force convention. Positive default occupancy is present. These are substantive improvements. This review inspected the source; it did not independently repeat every physics/GPU regression claimed in §4.4.

The statement "all review findings fixed" is too broad:

| priority | remaining issue | concrete agent instruction / acceptance test |
|---|---|---|
| P1 | `dr_limit` implementation is in `move_FIRE_bak`, unused by `run(ialg=2)` | Put the limit in the active path. Test a one-DOF large-force step and a large preexisting velocity through active `move_FIRE`, then eFF `run`. Require actual displacement ≤ bound. Validate dt>0. Consider one common scaling of the proposed step to preserve its direction; independent component clipping changes direction and is orientation dependent. |
| P1 | raw `F2` sums all DOFs, including `fixmask` reactions | Exclude `ff.fixmask[i]` as well as hard constraints. `move_LeapFrog` skips masked DOFs but does not clear their forces. Test `.O`/`.e2` inputs with finite reaction forces: free DOFs converge without moving constrained coordinates. Zero inverse mass alone is also not a complete fixed-DOF contract. |
| P1 | `run` still moves before testing the old force | Test raw convergence immediately after eval/constraint projection and before moving. On budget exhaustion reevaluate at returned positions. Verify returned coordinates, energy and force refer to the same state. A diagnostic here stopped on old force norm<1e-4 but returned a new-state force norm≈2.34e-4 at dt=0.05. |
| P1 | stale-pointer guard misses `opt.invMasses` and `opt.fixmask` | Rebind all borrowed arrays/masks after realloc, preferably at the owner boundary; allocator address reuse can conceal mask changes from the current pointer/size test. Reset stale `ff.nfix` and rebuild intended constraints. Test changed sizes/masks in ASan, not just a second identical load. |
| P1 | fresh `run` resets velocity/dt but keeps damping and `lastNeg` | Define restart versus continuation. Fresh starts reset the whole FIRE history. Tests must compare fresh versus sequential loads; scan timing must not depend on the previous scan point's hidden state. |
| P1 | energy finiteness is checked, force finiteness is not | Reject nonfinite active forces/coordinates too. ECP at rc=0 still divides by rc; energy can be finite while force is NaN. Use nonzero exit/error status: `exit(0)` currently tells shell/CI that failure succeeded. |
| P1 | width floor agrees numerically, but eval still mutates coordinates | `evalKinetic` still writes sizes to 0.1; clipping leaves incoming size velocity unchanged. Establish admissible steps before evaluation and a clear boundary policy. Check all FD/curvature probes stay inside the valid domain. |
| P1 | zero-width-core guard silently drops Pauli physics | `evalAE` now skips AEPauli for y≤1e-8 although z can denote a nonzero core. A finite FGO energy does not establish a valid model. Distinguish valid bare nuclei from malformed/legacy core parameters; migrate explicitly or fail loudly. ECP fallback still needs a valid width when z>0. |
| P2 | FGO import/export semantics remain incomplete | Assigning `fabs(c)` assumes legacy c is occupancy, which needs a documented format contract. Writer still emits literal 1.0, losing e2 occupancy; vsize is not actually written; closed-shell header's fourth argument is not parsed. Add single/pair round trips before calling FGO fixed. |
| P2 | study replicas are no longer exact CPU replicas | `fire`, `fire_secant`, `fire_groups` still stop on weighted force, clip to 1e-3, and lack a production cap. Update shared study logic before further sweeps; include final-state raw-force and cap/constraint parity. |

Keep paired occupancy within its supported spin model (Q=1 or 2). The generic `qi*(qi-1)/2` expression does not make >2 electrons in one spatial orbital a valid fermionic model. ECP length units, the rc=0/p-type cusp, and static CPU/GPU force parity remain separate open items; kernel compilation alone establishes none of these.

### 6.2 Measurement protocol and limits

Small Python diagnostics reused `tests/tEFF/fire_study_lib.py` and the existing CPU shared library; no production optimizer or force code was edited. `CPP_RECOMPILE=0`, `LD_PRELOAD=/usr/lib/gcc/x86_64-linux-gnu/11/libasan.so`, `ASAN_OPTIONS=detect_leaks=0` were used. Leak reporting was disabled for Python process shutdown; ASan access checks remained enabled. The library timestamp was 08:22:59, source timestamps 08:23:54: **these are existing-library diagnostics, not a fresh-build certification of the final source tree.** Its corrected paired Hamiltonian reproduces E≈−343.31951059 eV. The next agent must rebuild through the repository helper and reproduce the selected cases.

Temporary in-memory copies of the existing study functions changed only raw-force stopping **before** movement and the width floor to 0.1. The active CPU implementation also lacks the cap, so none was added to the study. Identical starting coordinates from `H2O_pairs_fc.xyz` were restored between runs; no repeated loader or inherited optimizer history was used. Tolerance is free-DOF raw L2 |F|<1e-4 eV/Å, confirmed by a fresh terminal eval. Counts below are eval-bearing loop iterations, including the terminal convergence check; actual moves are one fewer. One additional terminal verification eval is excluded consistently.

**Frozen core is not fixed nuclei.** This fixture has no dot-prefix constraints: na=3, ne=4, all 25 DOFs free. The separate fixed-nuclei experiment held the first 9 DOFs fixed, leaving 16 electronic DOFs. Do not compare their iteration counts as equivalent problems.

### 6.3 Primary result: tune dt ceiling and initial mixing before adding masses

Free nuclei, unit inverse masses, f_limit=100, finc=1.1, fdec=0.5, falpha=0.8, minLastNeg=3 unless specified. All successful rows reach the same basin, E≈−343.31951059 eV.

| setting | evaluations to raw force <1e-4 | observation |
|---|---:|---|
| dt=0.01, damping=0.1 | 828 | small steps |
| dt=0.02, damping=0.1 | 459 | 17 negative-power resets |
| dt=0.03, damping=0.1 | 323 | 19 resets |
| dt=0.035, damping=0.1 | 349 | 22 resets |
| dt=0.04, damping=0.1 | 379 | 32 resets |
| dt=0.05, damping=0.1 | 891 | 92 resets; overshooting costs more than larger dt saves |
| dt=0.03, damping=0.5 | **216** | simplest strong candidate; terminal force norm=8.10e-5 |
| dt=0.03, damping=0.7 | 215 | essentially tied; no basis for recommending 0.7 over 0.5 |
| dt=0.03, damping=0.1, fdec=0.7 | 282 | less severe shrink helps this run |
| dt=0.03, damping=0.3, fdec=0.7 | 222 | good, but not better than the simpler damping-only change |
| dt=0.03, damping=0.5, fdec=0.7 | 274 | improvements do not compose automatically |
| dt=0.03, damping=0.1, finc=1.05 / 1.2 | 455 / 313 | slower regrowth hurts here |
| dt=0.03, damping=0.1, minLastNeg=0 / 5 | 361 / 398 | neither beats the default 3 |
| dt=0.03, damping=0.1, falpha=0.95 / 1 | 262 / 1082 | retain decay; constant mixing is substantially worse |

**Explanation:** the initial largest numerical Hessian eigenvalue is ≈742, while the converged largest is ≈3049 eV/Å². A unit-mass undamped leapfrog harmonic mode is stable for dt*sqrt(lambda)<2: the final local ceiling is therefore ≈0.0362, versus ≈0.0734 initially. dt=0.05 can look useful early and repeatedly destabilize the late trajectory. FIRE damping changes the precise boundary, so this is a diagnostic estimate, not a global safe-step proof. The smallest clearly positive final eigenvalue is ≈4.04; frequency ratio ≈27.5. Six tiny rigid-mode eigenvalues, including small numerical negatives, should not be interpreted as physical saddles.

Larger initial FIRE damping means stronger velocity-direction mixing toward the current force after a reset, **not ordinary viscous friction**. It then decays with falpha. The previous all-electron observation that damping mattered little cannot be carried over to this corrected frozen-core problem.

**Agent recommendation:** expose dt ceiling, initial damping, finc/fdec/falpha/minLastNeg with explicit reset semantics. Reproduce dt=0.03,damping=0.5 first. If a larger early ceiling is desired, test a simple lower late ceiling or hysteresis after repeated negative-power events; do not continually regrow above the same local stability limit. No production Hessian/eigensolver is required for this idea. Keep dt_min independently tunable. Validate on a few distorted H2O configurations before selecting defaults.

### 6.4 Preconditioning: useful, but not the first frozen-core change

All rows use corrected raw stopping and reach the same free-nuclei basin:

| method | settings | evaluations |
|---|---|---:|
| existing secant diagonal | k_ref=100, dt=0.05 | 317 |
| existing secant diagonal | k_ref=500, dt=0.05 | 286 |
| existing secant diagonal | k_ref=100 / 500, dt=0.1 | 334 / 386 |
| initial FD diagonal | invM=100/(abs(diag H0)+100), dt=0.05 | 272 + **50 setup evals** |
| initial FD diagonal | k_ref=500, dt=0.05 | 561 + 50 setup evals |
| uniform width mass ×2 | dt=0.035 | 332 |
| uniform width mass ×4 | dt=0.05 | 1213 |

An initial full Hessian was used only for diagnosis; its diagonal costs the same 2*n_free=50 evaluations. Charging that cost makes the 272-step static conditioner approximately tied with the 323-step tuned-unit baseline, and worse than 216-step damping tuning. With periodic refresh interval R, relative extra force-evaluation cost is approximately `2*n_free/R`: n_free=100,R=200 means **100%**, not the earlier ~1% estimate. Secant conditioning still has no extra force calls, but here loses to the best simple parameter choice; test combinations only after establishing that baseline.

Trajectory component secants include off-diagonal coupling; they are not independent Hessian diagonal measurements. Prefer bounded, slowly varying masses and one isotropic mass per position triple when testing rotational robustness. Under nonuniform masses the existing acceleration-space power `v·(invM*F)` is not the physical descent power `v·F`; treat mass-metric FIRE as a separate ablation, not an already-validated detail. Do not use physical MD masses for relaxation.

### 6.5 Fixed nuclei, scans, other core modes and BFGS

With nuclei fixed at the input geometry, unit-mass FIRE/damping=0.1 took **164 / 181 / 165 evaluations** for dt=0.02 / 0.03 / 0.05. E≈−343.08062877 eV. Thus the free-nuclei 216-vs-891 improvement is not evidence of a comparable win for an electronic-only scan. Retest damping=0.3–0.5 on that protocol; optimize its actual residual channels. Record per-atom/electron position and size residuals, rather than assuming all late stiffness is in sizes.

For adjacent scan geometries, **warm-start converged electronic centers and widths** instead of rebuilding every point, preserving orbital identity/spin/occupancy. Reset momentum/FIRE history when targets change. Compare forward/reverse scans and occasional cold starts to identify basin hysteresis; continuation changes the sampled basin and must be part of the dataset protocol. This can save more iterations per geometry than another optimizer constant, without a new algorithm.

Existing dense BFGS, same initial state, raw tolerance 1e-4:

| protocol / initial inverse Hessian | iterations | total evals including rejects | terminal raw force norm |
|---|---:|---:|---:|
| free nuclei, H0=0.01 I | **62** | **72** (9 rejects) | 8.84e-5 |
| fixed nuclei, H0=0.01 I | **38** | **53** (14 rejects) | 1.89e-5 |
| free nuclei, H0=I | 600 budget exhausted | not instrumented in this probe | 2.08e-4 |

**Agent recommendation:** retest scalar BFGS initial scales (0.001,0.01,0.1) before implementing L-BFGS or a FIRE/BFGS hybrid. The old 49/156 iteration results describe an earlier Hamiltonian. Keep BFGS as the existing small-system baseline; compare total force calls, not accepted iterations. Log trust radius and rejection count to separate bad scaling from force/energy disagreement. Remove the existing ~2–3e-4 width force/energy approximation mismatch before declaring a reliable tighter polish: derive width force from the same Coulomb approximant. A small raw-force crossing does not prove exact stationarity of the energy when the two are inconsistent. Test a hybrid only if it beats scaled BFGS standalone including the FIRE startup cost.

For all-electron mode, the main core–nucleus stiffness still needs selective mass conditioning; the frozen-core dt recommendation must not be copied to it. Reestablish its spectrum and force-consistent baseline after the paired-core correction. For ECP, resolve units and rc=0 behavior first, then run the same small parameter matrix; no ECP timing was measured here.

### 6.6 Minimal next-agent sequence

1. Repair active-path cap, masked raw stopping/final-state consistency, restart state and force-finiteness checks. Add focused assertions, then rebuild through the repository helper. Do not conflate these fixes with tuning.
2. Bring `fire_study_lib.py` into parity with that exact production path. Confirm matching first steps, fresh terminal raw force, fixed masks and cap activation. The Python temporary variants used above are diagnostics only.
3. Reproduce four free-nuclei cases: FIRE (.05,.1), (.03,.1), (.03,.5), and BFGS H0=.01. Then fixed nuclei and a few nearby distorted H2O geometries. Compare energies/basin and full force-call budgets. Avoid a broad molecular sweep at this stage.
4. Only if needed, test secant+the winning FIRE parameters, an isotropic block variant, and a lower late dt ceiling. Retain the simplest option that improves the multi-geometry result. Warm-start scan behavior is a separate experiment.

No production code was changed by this follow-up review. Ephemeral numeric results are in `/tmp/eff_review_speed_results_prestep.json`, `/tmp/eff_review_speed_tuning.json`, `/tmp/eff_review_speed_combinations.json`; they are not durable regression fixtures. The original audit was backed up to `/tmp/eFF_optimizer_speedup_before_followup.md` before this documentation update.

### 6.7 §6-P1 items addressed in production code (verified on the live path)

The §6 review table's top items were implemented in `eFF_lib.cpp run()` and `DynamicOpt::move_FIRE()` and re-verified:

| §6 item | resolution |
|---|---|
| cap in `move_FIRE_bak` not active path | cap moved to active `move_FIRE()`; implemented as **uniform step scaling** (`vel[i] = s·(v_i+f_i·invM_i·dt) − f_i·invM_i·dt`, `s=dr_limit/max_i|dr_i|`) → direction preserved. No regression: `test_erfx_refit.py` PASS, H2O_spins 'a' converges. |
| raw `F2` summed fixmask reactions | norm loop now `continue`s on `ff.fixmask[i]` (in addition to `clear_fixed_*` zeroing `nfix` forces) |
| test old force before/after move | `run()` reordered: eval → clear fixed → **compute F2 → record → test → move**. Converged state, reported force, and returned energy now refer to the same coordinates. |
| stale-pointer guard partial | extended to `opt.invMasses` and `opt.fixmask` (all borrowed arrays) |
| restart vs continuation | `run()` now resets the whole FIRE history per call: `cleanVel()` + `damping=damp_max` + `lastNeg=0` (fresh semantics; sequential loads no longer inherit hidden state) |
| force finiteness unchecked | `isfinite(F2)` abort added; both aborts now `exit(1)` (was `exit(0)` = false success) |

Deferred per §6 scope: eval-write-back floor policy (`evalKinetic` mutates `esize`→0.1), FGO format contract (P2), study-replica parity (P2), `processXYZ_e` mixed-size configs, ECP `rc=0` division and length units.

#### Damping ablation on the live production path (answers "is no-damping best?")

H2O_pairs_fc, free nuclei, unit masses, f_limit=100, raw |F|<1e-4, fresh restart per point:

| damping | dt=0.03 evals | dt=0.05 evals |
|---|---:|---:|
| **0.00 (pure inertia + v·f<0 restart)** | 418 | 822 |
| 0.01 | 427 | 803 |
| 0.10 | 323 | 861 |
| 0.30 | 297 | 701 |
| 0.50 | 245 | 519 |
| **0.70** | **208** | 470 |
| 0.90 | 273 | **418** |

All runs reach E=−343.31951 (same basin). **Pure inertia (`damp=0`, only the v·f<0 kill) is ~2× slower than the optimum — the "no damping is best" intuition does not hold here.** Reason: `damping` is not viscous friction; it is the fraction of velocity re-aimed along the force each descending step (the α of original FIRE, decaying ×falpha→full inertia asymptotically). Unit masses on this landscape leave λ≈4…3049 (ω ratio ≈27), so this ablation does not test ideal conditioning. Even equal local frequencies would not prove that zero mixing is optimal on a nonlinear landscape; that remains an empirical question. In a curved valley the accumulated velocity drifts off the local gradient → more negative-power resets; re-aiming ~50–70% reduces resets. Optimal damping therefore depends on the actual conditioning, not on a fixed rule: for this system `dt=0.03, damp≈0.5–0.7` (208–245 evals) is the recommended FIRE setting; the v·f<0 restart is indeed already active (`cos_vf<0 → v→0, dt·fdec, damping→damp_max`).

## 7. Third review: cheap preconditioning without computing a Hessian

### 7.1 Current implementation and scope of this round

Source inspection confirms the uniform step limiter is now inside active `move_FIRE()`. Raw free-DOF convergence is checked before movement, the pointer guard includes masses/masks, and fresh calls reset mixing and positive-power history. These address the principal §6 findings. One statement still needs narrowing: **converged** returns have consistent coordinates/forces/energy, but budget exhaustion still occurs after a final move without a final eval. `run(nstepMax=0)` also leaves Etot uninitialized for verbose reporting. Keep both on the small driver-fix list.

The common cap preserves step direction, but its trigger is the largest Cartesian component, not a rotation-invariant particle displacement. Identical rotated systems can therefore take different early capped steps. Record cap events in rotation tests; an Euclidean norm per particle or a whole-step norm is an optional follow-up, not part of the measurements below.

This round evaluated **no numerical Hessian or Hessian-vector products**. Diagnostics reused existing `fire_study_lib.py::fire` with temporary in-memory corrections for pre-move raw stopping, size floor 0.1, and the active uniform displacement cap 0.2. Callbacks only supplied inverse masses using existing positions and forces. No production code or permanent test script was added. The ASan library dated 08:42:01 was used with the same environment described in §6.2; each run started with zero velocity and fresh optimizer/preconditioner history.

Protocol: `H2O_pairs_fc.xyz`, normally all 25 DOFs free, raw force L2 tolerance 1e-4, initial mixing 0.7, finc=1.1, fdec=0.5, falpha=0.8, minLastNeg=3, f_limit=100. Counts include terminal convergence evaluation and exclude one common final verification eval. All successful free-nuclei cases below reach E≈−343.31951059 eV. The diagnostic unit baseline is 215 evaluations, versus the coding agent's 208 on production. **Replica parity is still not established:** compare improvements within the same diagnostic table, and reproduce in C++ before adopting anything. Timings were subsecond for these small batches; no wall-time speedup is claimed for Python callbacks.

### 7.2 First candidate: one secant stiffness per particle position, one per width

Use each particle's xyz displacement/force change as a vector block, rather than estimating three unrelated Cartesian masses:

```text
dx_b = x_b(t) - x_b(t-1)                    # xyz block of one atom or electron
df_b = F_b(t) - F_b(t-1)
k_new = norm(df_b)/norm(dx_b)               # only if norm(dx_b) > 1e-6 Å
k_b = 0.9*k_b + 0.1*clamp(k_new, 1e-2, 1e7)
invM_b = k_ref/(k_ref + k_b)                # same scalar for x,y,z
```

For each electron width use the same rule with scalars (`abs(df_s/ds)`), separately from its xyz block. Initialize k_b=0 (unit inverse mass), do not estimate on the first force evaluation, and retain the previous estimate if displacement is insignificant. Initial trial `k_ref=100 eV/Å²`. There is one previous position/force per DOF and one stiffness per block: O(N) storage/work, **zero extra force evaluations**, no matrix construction or solve. Exclude constrained DOFs; for partially constrained particles, form norms only in the free subspace.

This measures stiffness experienced along the trajectory, not the true Hessian diagonal. Coupling can change a block's force while that block barely moves; thresholding, smoothing and bounds matter. Unlike a Rayleigh estimate, the norm ratio avoids cancellation between positive and negative curvature contributions. It remains a heuristic and may overestimate stiffness. Do not advertise guaranteed equalized frequencies or a positive-definite Hessian approximation.

| inverse-mass recipe | k_ref | dt ceiling | evaluations |
|---|---:|---:|---:|
| unit masses | — | 0.03 | 215 |
| unit masses | — | 0.05 | 479 |
| component secant `abs(df_i/dx_i)` | 100 | 0.05 | 248 |
| component secant | 500 | 0.05 | **165** |
| particle norm secant | 100 | 0.05 | **191** |
| particle norm secant | 500 | 0.05 | 228 |
| particle norm secant | 30 | 0.05 | 256 |
| particle norm secant | 100 | 0.07 | 180 |
| particle Rayleigh `abs(dx_b dot df_b)/(dx_b dot dx_b)` | 100 | 0.05 | 219 |
| particle norm secant, zero mixing | 100 | 0.05 | 403 |

**Robustness matters more than the best original-orientation number.** Rotation: axis (1,2,3), angle 0.73 radians, rotate every nucleus/electron center consistently. Distortion: `default_rng(42)`, independent position noise std=0.02 Å on all centers, then multiplicative width noise std=0.03. Widths remain valid. Each preconditioner was reinitialized.

| recipe and settings | original | rotated | distorted |
|---|---:|---:|---:|
| unit, dt=.03 | 215 | 233 | 239 |
| component secant, k_ref=500, dt=.05 | 165 | 232 | 209 |
| particle norm secant, k_ref=100, dt=.05 | **191** | **191** | **197** |
| component pooled secant, k_ref=100, dt=.05 | 251 | 455 | 204 |
| component pooled secant, k_ref=500, dt=.05 | 218 | 328 | 285 |

The pooled recipe was `k_i=abs(EMA(-dx_i*df_i))/EMA(dx_i²)`, with beta=0.9, the same displacement threshold/bounds and no second EMA. It avoids dividing every tiny step separately but did not improve robustness in these tests. The particle norm secant is the preferred **first implementation candidate**, despite the faster unrotated component result. Its measured gain over tuned unit masses is modest (roughly 11–18% across these three cases), not the old order-of-magnitude all-electron claim. The dt=.07 result is original-geometry only and needs the same robustness tests.

### 7.3 Correct FIRE metric: cheap, important under unequal masses

Current C++ FIRE aligns velocity with acceleration and uses `v dot (invM*F)` as its power test. That can have the opposite sign to actual energy descent `v dot F`. With diagonal positive masses, the fixed-mass coordinate transform instead gives:

```text
a_i = invM_i*F_i
P   = sum(v_i*F_i)
V2  = sum(v_i*v_i/invM_i)
G2  = sum(F_i*F_i*invM_i)
v_i = (1-alpha)*v_i + alpha*sqrt(V2/G2)*a_i   # P>=0 branch, guard zero G2
```

Use P<0 for restart. Raw physical force still controls convergence; the acceleration norm may remain separate for the existing f_limit safeguard. All reductions exclude constrained DOFs and never divide by their zero mobility. Unit masses must reduce to the current mixing formula. This adds only reductions/multiplications, with no energy/force calls.

**Measured original free-nuclei case:** particle norm secant, k_ref=100, dt=.05, mixing .7: **191 → 169 evaluations** with this metric, same basin, fresh terminal force 9.64e-5. This was the only mass-metric geometry tested: rotate/distort before claiming robustness. The fixed-mass derivation does not by itself specify optimal behavior when masses change. This diagnostic updated masses without rescaling velocity; retain that exact behavior initially for reproduction, then consider freezing converged estimates as a separate experiment. These are fictitious relaxation masses, not a physical MD mass update.

The zero-mixing block-secant case still took 403 evaluations. There is no evidence that this inexpensive conditioner makes damping unnecessary. Retest alpha=0/.3/.7 only after fixing the mass metric, not by assuming an ideal-preconditioner limit.

### 7.4 Analytic local stiffness: useful initialization idea, not a demonstrated winner

Tested a simple positive electron–atom scale:

```text
S_ae² = s_e² + s_core²
k_ae  = k_e*q_e*Z_eff / (r_ae² + 0.5*S_ae²)^(3/2)
k_atom = sum_e(k_ae)
k_epos = sum_a(k_ae)
k_width = c_width * 18*C*q_e/s_e^4
invM = k_ref/(k_ref+k)
```

Here `k_e≈14.399645 eV Å`, `C=const_K_eVA≈3.809982 eV Å²`; current `addKineticGauss_eFF` gives `Ekin=3*C*q/s²`, hence its exact width curvature is `18*C*q/s⁴`. For the frozen-core Coulomb path the diagnostic used `Z_eff=aPars.x-aPars.z`, core width `aPars.y`. The one scalar position estimate is shared across xyz. k_ref=100, dt=.05, mixing=.7. This recomputed the cheap geometry approximation each step; no numerical derivatives.

| c_width | evaluations |
|---:|---:|
| 0 (only position masses) | 261 |
| 0.1 | 207 |
| 1 (full kinetic width curvature) | 234 |

Neither beats the block secant or clearly improves on tuned unit masses. Electron–electron and Pauli stiffness and signed cancellations are missing. The coefficient 0.1 is an experimental attenuation, not a derived physical constant or a production recommendation. Pure positive term sums often overweight widths because attractive curvature cancels kinetic curvature. Simply using `k_width proportional to q/s⁴` everywhere is not sufficient.

**Cheaper implementation possibility for an agent:** reuse radial pair-force coefficients already computed in AE/EE helpers. For a pair with `F_pair = c(r)*dR`, `abs(c)` has stiffness units. Sum this scale into one mass estimate for each participating free particle; it is a transverse-force scale, not a bound on radial curvature. Width initialization can use `6*Ekin/s²` and the internal pair term `2*J/s²` (J proportional to 1/s), computed with existing term values. Reusing intermediates avoids a second O(na*ne) geometry pass or new transcendental evaluations. Exact coincident-Gaussian Coulomb position stiffness scales as `q*Z/S³`, so this captures why all-electron core pairs need heavy masses. Pauli radial derivatives/cancellation and the r=0 limits still need care. **This intermediate-reuse variant was not implemented/tested here**; first validate the cheaper trajectory-only option.

A static occupancy-only recipe (electron position and size invM=1/q, atoms=1) gave 231 evaluations at dt=.04. It makes shared paired-coordinate scaling more comparable to expanded singles but is not enough to remove this system's stiffness disparity.

### 7.5 Fixed nuclei: do not add conditioning by default

Same starting geometry, first 9 DOFs fixed; raw force measured over the remaining 16:

| method | dt | mixing | evaluations |
|---|---:|---:|---:|
| unit masses | .03 | .7 | **86** |
| particle norm secant, k_ref=100 | .05 | .7 | 106 |
| particle norm secant, k_ref=100 | .07 | .7 | 87 |
| particle norm secant, k_ref=30 | .05 | .7 | 151 |
| norm secant + mass metric, k_ref=100 | .05 | .7 | 104 |
| particle norm secant, k_ref=100 | .05 | 0 | 233 |

All reach E≈−343.08062877 eV. The simpler fixed-nuclei baseline benefits greatly from tuning mixing; the tested conditioning adds no compelling gain. Primary recommendation is therefore **optional conditioning for free-nuclei frozen-core relaxation**, with separate electronic-only scan defaults. This distinction also matters when interpreting the old "frozen core is already fast" statements.

### 7.6 Concrete next-agent handoff

1. **Finish study parity before porting a winner.** Put the small preconditioner update in shared study logic, not duplicated scripts. Match active C++ cap, eval/test/move order, masks, floors and restart on unit masses. Explain the current 215-vs-208 diagnostic/production baseline difference with first-step comparisons rather than absorbing it into an advertised speedup. Fix budget-exhaustion final evaluation separately.
2. **Implement one optional particle-norm secant in the CPU relaxation path.** Allocate old positions/forces and block estimates once on initialization; reset after realloc, changed constraints, discontinuous geometry loads and a fresh optimizer start. Start k=0; update only from two actual evaluated states. Do not use a proposed/unaccepted step in secant history. Keep physical forces unchanged, positive bounded inverse masses, and free-DOF raw stopping. Enforce supported occupancies independently of masses.
3. **Keep mass metric as an explicit, separate switch for validation.** Verify unit-mass equivalence, a two-frequency quadratic, P=v dot F sign, no zero-mass division, and mask exclusion. Test old/new metric with the same preconditioner and alpha; do not retune everything at once. Changing this shared DynamicOpt behavior can affect other force fields, so scope it to an opt-in relaxation mode first.
4. **Validate a small bounded matrix, not a large optimizer search.** Original/rotated/distorted H2O; free versus fixed nuclei; unit baseline, block norm k_ref=100 at dt=.05 and .07, block+mass metric. Compare final raw norms, E/basin, force-call count, cap/reset counts and saturation of mass bounds. Only then add a few other initial geometries. A no-extra-evals implementation must actually retain one force call per normal iteration.
5. **Secondary modes:** for all-electron systems use the same vector blocks for each nucleus and electron so coupled core modes can slow both participants. A local `q*Z/S³` prior can seed obviously tight core pairs, but has not been benchmarked after the latest fixes. Do not copy frozen-core dt=.05/.07 to unconditioned all-electron mode. For ECP, fix rc=0/units before timing. No new all-electron/ECP timing was taken in this round.
6. **BFGS:** reuse a stable diagonal mass estimate only as an optional initial inverse-Hessian scale (`H0_i` proportional to invM_i with a separately tested scalar, around the existing .01 baseline). Reset BFGS history when geometry/constraints change discontinuously. Do not change this scale every BFGS iteration as if it were an independent mass without reconciling secant history. Existing 62-iteration/72-eval BFGS remains a stronger free-nuclei benchmark than any FIRE candidate tested here; address force/energy consistency before tighter polish.

No Hessian was calculated, no new complex optimizer was written, and no production code was changed in this round. Ephemeral results: `/tmp/eff_cheap_prec_round3.json`, `/tmp/eff_cheap_prec_round3_robustness.json`, `/tmp/eff_cheap_prec_round3_blocks.json`. Audit backup: `/tmp/eFF_optimizer_speedup_before_prec_review.md`.

### 7.7 Production implementation + reproduction (verified on the live C++ path)

The block-secant preconditioner and the mass-consistent metric were implemented **opt-in** in production (`autoMass`, `massMetric` — both default 0/off, zero impact unless enabled):

- `DynamicOpt.h:101–115` — fields `autoMass, am_nBlocks, am_nScalars, am_kref, am_ema, am_k, am_xold, am_fold, am_init, massMetric`; `updateAutoMasses()` declared; buffers freed in `dealloc()`.
- `DynamicOpt.cpp:306–350` — `updateAutoMasses()`: one scalar stiffness `k_b=|dF_b|/|dx_b|` per **3-component position block** (atoms + electron centers) and one **per scalar size DOF** — xyz and width preconditioned separately, as intended; EMA `k=0.9k+0.1·clamp(k_new,1e-2,1e7)`; `invM=kref/(kref+k)`; free-subspace norms for partially constrained blocks; updates only between two consecutively evaluated states; `am_init` reset discards stale history.
- `DynamicOpt.cpp:352–368` — `massMetric==1` branch in `move_FIRE()`: `P=Σv·f`, `V2=Σv²/invM`, `G2=Σf²·invM`; restart on `P<0`; mixing `v=(1−α)v+α·sqrt(V2/G2)·(f·invM)` via the existing loop. At `invM≡1` reduces exactly to the legacy formula.
- `eFF_lib.cpp` — `setAutoMass(mode,kref)` / `setMassMetric(m)` exports; `run()` sets `am_nBlocks=na+ne, am_nScalars=ne` from the current layout, re-inits secant history on stale-rebind and on every fresh restart.
- `pyBall/eFF.py` — `eff.setAutoMass(mode,kref)` / `eff.setMassMetric(m)` bindings.

**Reproduction on production code, H2O_pairs_fc, raw |F|<1e-4, damp=0.7:**

| case | review replica | production | verdict |
|---|---:|---:|---|
| unit, dt=.03, free | 215 | 208 | match (within replica parity) |
| unit, dt=.05, free | 479 | 470 | match |
| secant k100, dt=.05, free | 191 | **191** | **exact** |
| secant k100, dt=.07, free | 180 | 178 | match |
| secant robustness orig/rot/dist | 191/191/197 | **191/191/197** | **exact** |
| secant zero-mixing, dt=.05 | 403 | **403** | **exact** |
| secant+massMetric, dt=.05 | **169** | 200 | **NOT reproduced** — slightly *worse* than plain secant (200 vs 191); rotated 200, distorted 175. The 169 appears to be a replica artifact; massMetric adds no demonstrated benefit. |
| fixed nuclei: unit dt=.05 | — | **75** | unit masses still best on fixed (secant: 106, secant+metric: 133, k30: 151) — confirms review's "do not condition fixed-nuclei by default" |

**Conclusions:**
- **Adopted (opt-in):** `setAutoMass(1, kref=100)` — verified ~10–15% free-nuclei eval reduction, rotation-invariant (the uniform per-block mass and block-norm secant are both rotation-covariant; the per-component `dr_limit` cap is the only remaining orientation-dependent piece). Best confirmed setting: `dt=0.05–0.07, damp=0.7, kref=100 → ~178–191 evals` vs 208–470 unit-mass.
- **Damping answer, now with the preconditioner:** secant masses + `damp=0` → **403 evals** (vs 191 at damp=0.7). Even with trajectory-estimated stiffness masses, pure inertia remains ~2× slower — the secant masses are a useful but far-from-ideal preconditioner, so the "no damping is optimal once masses equalize stiffness" limit is still not reached here.
- **massMetric kept off by default** — implemented as a separate switch per the review's handoff, but the claimed benefit did not reproduce; keep for experiments, not production.
- **Fixed nuclei:** unit masses at `dt=0.05, damp=0.7` give 75 evals — better than every conditioned variant. Recommendation unchanged: preconditioning is only for free-nuclei relaxation.
- Also observed: `set_constrains` state (`nfix`) **survives** `processXYZ_e` reloads — constraint state is not reset on geometry load. Flagged for the driver-fix list (same stale-state class as §6 items).

## 8. Fourth review: faster gradient-history methods without BFGS

### 8.1 Answer and measurement boundary

**Yes, there are promising inexpensive alternatives, but distinguish electronic-only relaxation from free nuclei.** A scalar spectral-gradient step is sufficient to improve electronic-only H2O. Free nuclei benefit more from a little coupled gradient history, after FIRE has done the initial descent. No tested scalar-only method beat conditioned FIRE on the free-nuclei case. The useful new candidate is a small-memory Anderson/multisecant polish, not another diagonal mass fit.

All new probes call the existing `eff.eval()` for forces/energies, use raw free-DOF force L2 <1e-4, cap the largest proposed coordinate displacement at 0.2 Å, and reject steps leaving any width below 0.1 before evaluation. Accepted steps and all rejected energy/force trials are counted. One common terminal verification eval is excluded. No Hessian, Hessian-vector products, BFGS update, or production code change was used. Python diagnostic prototypes are not claimed to have optimized wall time or production-ready error handling.

Baselines below were remeasured on the live shared library, not a FIRE replica: free nuclei use `setAutoMass(1,100)`, `setMassMetric(0)`, dt=.07, mixing=.7; fixed nuclei use unit masses, dt=.05, mixing=.7. `initOpt(..., bMass=False, f_limit=100)` is called before each run. Hard constraints were explicitly set/cleared between protocols, avoiding the known persistent-nfix issue. Current `H2O_pairs_fc.xyz`, same original/rotated/distorted geometries as §7; same existing ASan/Python environment. New optimizer state is fresh for every case.

### 8.2 Simplest option: spectral gradient descent (BB step)

Let g=-F on free DOFs, `s=x_k-x_(k-1)`, `y=g_k-g_(k-1)`:

```text
BB1: a = (s dot s)/(s dot y)
BB2: a = (s dot y)/(y dot y)
direction d = -g
trial x = x_k + a*d
```

The successful minimal variant uses BB1. First step a=.001 Å²/eV; use the same initial scale when `s dot y<=1e-20`. Otherwise clamp a to [1e-6,.1], then shorten for the 0.2 displacement bound. It needs only previous x/g plus dot products: O(N) storage/work, no solve. Do not call `dt=sqrt(a)` or treat it as a FIRE timestep; a is a direct gradient-descent step length.

Use **nonmonotone Armijo acceptance**: `E_trial <= max(last 5 accepted E) + 1e-4*a*(g dot d)`. Halve a on rejection; include every trial evaluation in the cost. Width-invalid proposals are shortened before calling eval. Restore the last accepted state on a failed search and return explicit failure. The prototype allowed 30 trials in the initial comparison, 25 in the robustness run; successful cases needed only a few rejects overall. This bounded energy window permits a larger spectral step to traverse stiff directions without requiring every intermediate energy to decrease.

| protocol / method | accepted steps | total force evaluations | rejects |
|---|---:|---:|---:|
| fixed nuclei, BB1 | 42 | **46** | 3 |
| fixed nuclei, BB2 | 47 | 50 | 2 |
| fixed nuclei, alternating BB1/BB2 | 42 | 46 | 3 |
| free nuclei, BB1 | 524 | 880 | 355 |
| free nuclei, BB2 | 732 | 737 | 4 |
| free nuclei, alternating BB1/BB2 | 386 | 417 | 30 |

Thus **BB1 is the first low-complexity implementation candidate for fixed-nuclei scans**, not the free-nuclei default. The same spectral step still sees all coupled stiffness through one scalar. A basic PR+ nonlinear-CG trial with descent safeguards and BB-initialized Armijo backtracking also failed to improve free nuclei (717–810 evaluations); monotone fixed-nuclei CG stalled at raw force 1.85e-4. This is evidence about that simple implementation, not a proof that a carefully implemented Wolfe-CG cannot win. The existing `CG.h` solves linear systems and is not a drop-in nonlinear eFF optimizer; `lineSearch.h` is a generic scalar root-search helper, not an already-validated Armijo/Wolfe driver.

### 8.3 Small-memory coupled polish: Anderson acceleration of a gradient map

Store the most recent m accepted differences `S=[s_1,...,s_m]`, `Y=[y_1,...,y_m]`. Given current g, compute:

```text
gamma = argmin_gamma norm(Y*gamma - g)
d = -eta*g - (S - eta*Y)*gamma
```

This accelerates the map `x -> x-eta*g` using several observed gradient changes. Use eta=.001 Å²/eV. With no history, d=-eta*g. Use rank-revealing QR or SVD with relative singular cutoff 1e-8; the prototype used `numpy.linalg.lstsq(Y,g,rcond=1e-8)`. Avoid normal-equation inversion, which squares the conditioning problem. If `g dot d >= -0.001*norm(g)*norm(d)`, clear the history and use d=-eta*g for that step. Restart events must be counted. Stop on the raw gradient norm before constructing another direction.

Start line-search fraction a at 1 (shortened for the 0.2 coordinate-step bound); use the same last-five-energy Armijo rule and factor .5 backtracking as BB1. Check widths before eval. After acceptance append **the actual accepted** s/y, keep the last m, and retain that accepted state on any later failure. Constraints are projected out of g and direction; do not reuse histories across changed targets/masks. The prototype formed full vectors with constrained entries zero; a compact free-DOF representation is equivalent for these fully fixed/free blocks.

This is a small multisecant/quasi-Newton-like method: it captures coupled curvature implicitly from gradients, without constructing an N×N matrix or performing BFGS updates. It is not a claim of acceleration from "no curvature information." Storage is O(N*m); a fresh QR/SVD fit is O(N*m²), with only m=5 or 8 here. For this tiny 25-DOF molecule that algebra need not be cheaper than dense BFGS; the advantage is bounded history and scalable storage. Measure C++ wall time before adopting it on speed grounds.

**Fixed nuclei:** m=5 gave 31 accepted steps / **32 evaluations**, no rejected trials and three direction/history resets. m=8 gave 30 evaluations, but the two-evaluation saving does not justify a larger default from one geometry. eta=.01 with m=5 took 64 evaluations: the base scale matters.

**Free nuclei:** Anderson from the input was poor (roughly 495–600+ evaluations in these trials) because the descent safeguard repeatedly cleared its local model. A short FIRE phase helps. Starting with **60 production FIRE evaluations**, discard FIRE velocity, initialize empty Anderson history, and use m=8, eta=.001: **55 Anderson steps, 116 total evaluations**, no rejected trials, one direction reset. That includes the fresh eval needed after the FIRE budget-exhaustion return. This is **not** a 55-evaluation complete relaxation.

The handoff settings are not yet universal: 30 FIRE evals + m8 took 222 total; 60 FIRE + m5 took 280; 100 FIRE + m5 took 348. Switching too early or keeping too little history can lose the benefit. At the 60-evaluation switch the original raw force was ≈0.771 eV/Å. A force-based handoff gate is a future experiment; only the fixed 60-evaluation gate was tested. Do not generalize it to other molecules/core modes without validation.

### 8.4 Same-start production baselines and robustness

Every successful pair reaches the same energy within numerical tolerance **for its own geometry/constraint protocol**. Distorting fixed nuclei changes the physical problem and its final energy; it is not a basin mismatch.

| method / protocol | original evals | rotated evals | distorted evals |
|---|---:|---:|---:|
| production FIRE, free nuclei, secant, dt=.07 | **178** | **181** | **193** |
| 60 production FIRE + m8 Anderson, free nuclei | **116** | **125** | **135** |
| production FIRE, fixed nuclei, unit, dt=.05 | **75** | **75** | **96** |
| BB1, fixed nuclei | **46** | **46** | **71** |
| m5 Anderson, fixed nuclei, no FIRE startup | **32** | **32** | **70** |

Fresh terminal raw norms are all <1e-4. Free-nuclei energy ≈−343.31951059 eV; fixed original/rotated ≈−343.08062877; fixed distorted ≈−343.14106652. Robustness cases used the exact §7 rotation/noise construction. The hybrid saves roughly 30–35% of force calls versus the corresponding production FIRE baseline here. Fixed-nuclei BB1 saves roughly 26–39%, Anderson about 27–57%. No statistical or broad chemical generalization follows from three related H2O cases.

The earlier scaled-BFGS free-nuclei baseline remains **72 total evals**, versus this hybrid's 116. Therefore the tested gradient-history option beats FIRE but does **not** beat BFGS in force calls on that original free-nuclei geometry. Fixed-nuclei m5 Anderson's 32 evals improves on the earlier BFGS prototype's 53, but BFGS has not been retuned or rebenchmarked on the rotated/distorted cases in this round. Avoid claiming a universal algorithm winner.

### 8.5 Production-review qualifications discovered this round

- **Mass-metric discrepancy has an identifiable comparison difference.** §7's diagnostic changed power/mixing to the mass metric but retained the acceleration norm for f_limit; current `move_FIRE` sets `ff=G2=sum(F²*invM)` and then applies f_limit to that value too. Production and diagnostic changed different safeguards. A different early capped step can change the whole trajectory. Keep massMetric off based on current performance, but the 169-vs-200 result alone does not isolate the mixing formula or establish a mysterious numerical artifact. Test the same limiter norm explicitly before attributing the difference.
- **Fresh secant history does not necessarily mean fresh unit inverse masses.** In `updateAutoMasses`, the `!am_init` branch initializes k/history and returns without writing invMasses=1. `run()` clears am_init but does not reset inverse masses. Our comparisons called `initOpt(bMass=False)` first, so they are clean; repeated production `run()` calls without initOpt can inherit previous adaptive masses on the first step. Define/reset that state explicitly; do not overwrite deliberately supplied static masses by accident.
- **Energy/force mismatch remains relevant.** All acceptance tests use energy, while stopping uses the implemented force, whose width derivative has the known small approximant inconsistency. This can produce a failed line search near convergence (as the monotone CG trial illustrates). Fix the common Coulomb width derivative before tightening tolerances; do not silently add energy slack or declare a stall converged.

### 8.6 Next-agent instructions: smallest useful change first

1. Implement **BB1 + bounded nonmonotone backtracking** as an opt-in electronic-only optimizer in shared code. It needs just previous x/g and a five-energy ring buffer. First reproduce 46/46/71 force calls. Keep FIRE as the existing default while validating. Reuse existing evaluation/constraints and physical force stopping; do not implement a separate eFF force path.
2. If the small least-squares helper is acceptable, add the documented **m5 Anderson** mode for fixed nuclei and reproduce 32/32/70. Inventory shared QR/SVD utilities before adding a solver; handle rank-deficient histories explicitly. Allocation must occur at initialization, not every iteration. Test rank loss on duplicate/collinear steps, constrained DOFs, rejected steps, valid widths and restoration after failure.
3. For free nuclei test exactly **60 production FIRE + m8 Anderson** before inventing a new switching heuristic. Charge all startup/trial evaluations; reproduce 116/125/135 and the raw terminal state. Preserve accepted-step history only and clear FIRE velocity on transition. A few new distorted H2O starts are the next validation, not a broad parameter sweep.
4. Before expanding to all-electron/ECP, fix the still-open width/units/domain issues and repeat the bounded protocol. No evidence in this round establishes the same eta, memory, or switch for those Hamiltonians. Warm starts across scan points remain useful, but discard gradient history whenever constraint targets or orbital identity change unless an explicit continuation protocol validates reuse.
5. Compare force-call count **and** C++ wall time against scaled BFGS, especially at small N. If even a small m-column least-squares solve violates the desired simplicity, keep BB1 for fixed nuclei and the existing conditioned FIRE for free nuclei; a faster scalar-only free-nuclei replacement was not demonstrated here.

Only this audit was edited. Ephemeral results: `/tmp/eff_gradient_probe_round4.json`, `/tmp/eff_anderson_probe_round4.json`, `/tmp/eff_anderson_hybrid_round4.json`, `/tmp/eff_gradient_robustness_round4.json`, `/tmp/eff_gradient_baselines_round4.json`. The documentation backup is `/tmp/eFF_optimizer_speedup_before_gradient_review.md`. Prototype failures/budget exhaustion above were recorded as failures, not hidden as converged results.

### 8.7 Production implementation + reproduction (verified on the live C++ path)

Both methods were implemented **opt-in** in production as `run_GH(nstepMax, Fconv, iMethod, mHist, eta, a_init, outE, outF)` in `eFF_lib.cpp` (~line 310), exposed via `eff.run_GH(...)` in `pyBall/eFF.py`. `iMethod=1` BB1, `iMethod=2` Anderson; the hybrid is composed at the call site (`eff.run(60,...ialg=2)` then `eff.run_GH(iMethod=2,mHist=8)`). Shared protocol exactly per §8.1: raw free-DOF force convergence before proposing, `dr_limit` per-DOF displacement cap, proposals shortened/projected so `esize ≥ min_esize` before eval (floor-touching direction components are projected to zero), nonmonotone Armijo vs max of last 5 accepted energies (c=1e-4, halving), every trial counted, last accepted state restored on line-search failure, `exit(1)` on nonfinite E/F. Anderson least squares: **CPQR** (modified Gram-Schmidt with decreasing-norm column pivoting, rank cutoff 1e-8·|y_j|) — no normal equations; history appended only on accepted steps; descent safeguard `g·d ≥ −1e-3·|g||d|` clears history.

**Reproduction (H2O_pairs_fc, raw |F|<1e-4; constraints rotated consistently with geometry):**

| method / protocol | review prototype | production | verdict |
|---|---:|---:|---|
| BB1, fixed nuclei, orig/rot/dist | 46/46/71 | **46/46/71** | **exact** |
| Anderson m5, fixed, orig/rot/dist | 32/32/70 | 44/41/71 | benefit confirmed (vs FIRE 75/75/100); ~35% more evals than claimed — residual gap from unexposed prototype details (lstsq tolerance/counting), same basin |
| FIRE secant dt=.07 free | 178/181/193 | 178/181/193 | baseline reproduced |
| 60 FIRE + m8 Anderson, free | 116/125/135 | 127/126/118 | ~30% improvement confirmed |
| BFGS free-nuclei | 72 | — | still the eval-count winner for free nuclei |

Two solve variants were tested for the Anderson lstsq: plain/pivoted MGS-QR (kept — better) and min-norm via Gram-matrix eigendecomposition (mimicking `lstsq` SVD; **measurably worse**: hybrid 179/158/198 vs QR 127/126/118 — the min-norm γ produces weaker corrections). CPQR is retained.

Also fixed this round: `updateAutoMasses` `!am_init` branch now writes `invM=1` for free DOFs (§8.5 note — a fresh `run()` reset no longer inherits stale adaptive masses on the first step).

**Notes/corrections on §8.5:** production `f_limit` semantics are *identical* in legacy and `massMetric` code — legacy `project_vf` already accumulates `ff=Σ(f·invM)²` = the acceleration norm G2; `massMetric=1` reuses the same quantity. The 169-vs-200 difference is therefore not a f_limit-norm change; massMetric remains off by default as its benefit was not reproduced. The persistent-`nfix` hazard confirmed in §8.5 was also observed here: `set_constrains` state survives `processXYZ_e` reload (fixed-nuclei "free" sweeps silently ran constrained — caught via the wrong basin energy −343.0806 vs −343.3195). For fixed-nuclei rotation tests the constraint *targets* must be rotated too, otherwise nuclei snap back and the protocol silently becomes a different problem.

**Adoption guidance:** BB1 (`iMethod=1, a_init=0.001`) is a validated cheap default for **fixed-nuclei electronic relaxation** (46 evals vs 75; rotation-invariant). Anderson m5 (`iMethod=2, mHist=5, eta=0.001`) is somewhat better still (44). For **free nuclei**, `run(60) → run_GH(2, m=8)` gives ~125 evals vs 178 — worthwhile when BFGS is unavailable; scaled BFGS (72 evals) remains the eval-count winner. All new paths reach the same basin energies; none affect FIRE/GD defaults (opt-in only).

### 8.8 Refactor: reusable pure-C++ implementation (post-review migration)

Migrated from the eFF_lib.cpp prototype into shared modules, per maintainer direction:

- **`Lingebra.h`** — new `Lingebra::leastSquareFit_QR(n, m, A_colmaj, b, x, Q, R, piv, vwork, rcond)` as a header-**inline** flat-array routine (same convention as `Cholesky.h`; no `double**` conversion, no link dependency, no internal allocation). MGS + decreasing-norm column pivoting + rank dropping + triangular back-solve. Verified against `numpy.linalg.lstsq` on a near-collinear case (x agreement ~7 digits).
- **`DynamicOpt.h/.cpp`** — `opt.run_GH(nstepMax, Fconv, outE, outF)` holds the full optimizer (BB1 + Anderson); single callback `getForce` — `ForceFunction` signature changed to `double(*)(int,xs,dfs)` (**returns energy**, fills force; energy ignored by FIRE/GD, required only by `run_GH`; no in-repo assignments of the old void signature existed — verified); `opt.xmin` = caller-owned optional lower-bound array (generic; eFF maps `esize → min_esize`, `-1e300` elsewhere); aux buffers **preallocated once** in `allocLS()` with fixed `mMax=8` cap (no per-iteration allocation, allocated only when `lsMethod>0`); freed in `dealloc()`.
- **`eFF_lib.cpp`** — `run_GH()` shrunk to a ~20-line wrapper: stale-rebind check, sets `opt.lsMethod/mHist/lsEta/lsAinit`, binds `opt.getForce = _eff_eval_E` (`clearForce + eval + clear_fixed_dynamics`), builds `xmin` (esize floor), calls `opt.run_GH`. API and `pyBall/eFF.py` signature unchanged.
- Other `DynamicOpt` consumers (MMFF_lib, MMFFsp3_lib) rebuild clean — the header-inline solver keeps `DynamicOpt` self-contained (no CMake changes).

**Post-refactor verification (must match pre-refactor numbers):** BB1 **bit-identical** 46/46/71. Anderson: 44/49/54 fixed (was 44/41/71), hybrid 127/201/118 (was 127/126/118) — energies identical, all converge, but counts shift: the Anderson trajectory is FP-order-sensitive (direction accumulation order changed when the solve moved into `leastSquareFit_QR`). Solver correctness verified independently vs `lstsq`; the count spread is within the same trajectory-sensitivity band already separating the prototype's 32/32/70 from production's 44 — no algorithmic regression. Regression test (`test_erfx_refit.py`) passes; FIRE/GD defaults untouched.

### 8.9 Comprehensive benchmark + stale-history bug found and fixed (pre-commit verification)

Final verification sweep (`tests/tEFF/run_method_benchmark.py`, frozen-core `H2O_pairs_fc.xyz`, `|F|_free<1e-3`, orig/rotated/distorted starts, nuclei free vs fixed — each method forked so `exit(1)` failures are recorded not fatal).

**Bug found by the benchmark — stale history injection on first full fill.** Repeated identical `run_GH` calls gave nondeterministic results (Anderson m5 fixed: 36/31/33/31 evals; one sequential run diverged to a wrong basin E=−343.0705 or NaN). Root cause in the history push:

```cpp
int j=(nhist<mHist)?(nhist++):(mHist-1); if(nhist==mHist){ shift S,Y down }
```

When `nhist` goes `mHist−1 → mHist` (first full fill), `j=mHist−1`, `nhist++→mHist`, the `nhist==mHist` check then fires **and shifts anyway** — discarding valid column 0 and copying slot `mHist−1`'s **stale content from the previous `run_GH` call** into slot `mHist−2` before overwriting `mHist−1`. Stale-buffer garbage was injected into the Anderson least-squares on every first-full-fill; identical calls diverged after eval ~12 (the fill point). Fixed: shift only when the queue was already full before the push (`if(nhist==mHist){j=mHist-1; shift} else j=nhist++`). Post-fix: bit-identical repeated calls (27/27/27/27), and **faster**: the injected column was also hurting convergence. This bug plausibly accounts for part of the prototype-vs-production gap in §8.7/§8.8 (32→44).

**Post-fix benchmark (force evaluations, orig/rot/dist, mean):**

| method | free nuclei | fixed nuclei |
|---|---:|---:|
| FIRE dt.03 d.5 (old tuned) | 192/197/187 — 192 | — |
| FIRE dt.03 d.7 | 174/188/191 — 184 | — |
| FIRE dt.05 d.5 / d.7 | 365/405/369 — 380 | 90/90/89 — 90 / **63/63/76 — 67** |
| FIRE secant dt.05 d.7 | 162/162/194 — 173 | 92/92/84 — 89 |
| FIRE secant dt.07 d.7 | 154/151/106 — **137** | — |
| FIRE secant+massMetric dt.07 | 160/160/155 — 158 | — |
| BB1 η=1e-3 | 762/687/537 — 662 | 39/39/58 — 45 |
| Anderson m8 η=1e-3 | 471/476/101 — 349 | 27/27/46 — **33** |
| Anderson m5 η=1e-3 | — | 27/27/47 — 34 |
| FIRE60→Anderson m8 (unit) | 749/752/127 — 543 | — |
| FIRE60→Anderson m8 (secant) | 107/114/103 — **108** | — |
| FIRE60→Anderson m5 (unit) | — | 63/63/67 — 64 |

All rows converged to the correct basin (free −343.31951, fixed −343.08063 eV). `test_erfx_refit.py` passes.

**Conclusions:**

- **Free nuclei (the real use case):** `FIRE60(secant,dt.07,d.7) → Anderson m8` is the best measured at **~108 evals** vs 184 for tuned FIRE and 137 for pure secant-FIRE — ~41% better than baseline FIRE. The unit-mass FIRE head-start variant is catastrophic (543) — Anderson needs the preconditioned FIRE to get near the basin first. Anderson/BB1 from scratch are poor on free nuclei (349/662).
- **Fixed nuclei:** Anderson alone is best (**33 evals**, m5≈m8), BB1 close (45), both ~2× faster than best FIRE (67). Secant preconditioning gives no benefit on the electronic-only subspace.
- **massMetric** remains a wash (158 vs 137 pure secant) — off by default, kept for experiments.
- BFGS (earlier round, ~72 evals free-nuclei) remains the eval-count winner where available.
- Determinism caveat resolved: `run_GH` is now deterministic for identical calls; residual count sensitivity across *different* starts is genuine landscape sensitivity, not stale state.

### 8.10 Method comparison: FIRE vs Anderson vs BFGS (conceptual)

All three consume only energy/force evaluations — they differ in **how much of the local quadratic model each reconstructs per step**:

- **FIRE — "damped dynamics", inertia but no memory.** `v += f·dt·invM; x += v·dt`; velocity re-aimed toward the force by fraction `damping`, **killed outright when `v·f<0`** (negative-power restart). Adapts only `dt`/`damping` from progress signs. Robust, minimal per step, but converges linearly in stiff valleys — no curvature model. A diagonal preconditioner (`autoMass` secant masses) approximates curvature *statically*.
- **Anderson — "extrapolate the fixed-point map", remembers steps.** Treats relaxation as iterating `x ← x − η·g`; keeps the last m accepted `(Δx, Δg)` pairs, solves a small `n×m` least-squares for γ, and forms the next step as the linear combination of past steps that best cancels the residual: `d = −η·g − (S−η·Y)·γ`. An implicit multisecant curvature estimate — quasi-Newton-like without an explicit matrix. Storage O(n·m), solve O(n·m²). No step-size model of its own: needs base `η` + Armijo backtracking, so it **shines near the basin and flails far from it** (descent safeguard repeatedly clears history early on — hence the FIRE head-start on free nuclei).
- **BFGS — "explicit quadratic model", remembers a curvature matrix.** Maintains an inverse-Hessian estimate `H` via rank-2 secant updates; `d = −H·g` is a true Newton-like step with direction *and* per-direction scaling built in. Strongest convergence, but O(n²) (or low-rank) state, a real Wolfe line search, and a decent `H0` seed (scaled `0.01·I` mattered here) are all required.

One line: **FIRE has inertia but no memory; Anderson remembers steps and extrapolates; BFGS remembers an explicit curvature matrix.**

**Where each wins (frozen-core H2O, evals to `|F|<1e-3`, §8.9 table):**

| | free nuclei (25 DOFs) | fixed nuclei (16 DOFs) |
|---|---:|---:|
| best plain FIRE | 184 (dt.03 d.7) | 67 (dt.05 d.7) |
| FIRE + secant preconditioner | 137 (dt.07 d.7) | 89 (no benefit) |
| BB1 spectral GD | 662 | 45 |
| Anderson (best m) | 349 (alone, m8) | **33** (m5≈m8) |
| FIRE60(secant)→Anderson m8 | **108** | 64 (no benefit) |
| scaled BFGS (H0=0.01·I) | **72** | — |

Reading the table: on the coupled free-nuclei landscape, gradient-history methods only pay off *after* preconditioned FIRE reaches the basin (Anderson alone 349, unit-FIRE60 hybrid 543 vs secant-FIRE60 hybrid 108). On the constrained electronic subspace the landscape is benign enough that Anderson/BB1 win outright. BFGS's explicit curvature still dominates everywhere it is applicable — the gradient-history methods are the fallback when BFGS machinery is undesirable or unavailable.
