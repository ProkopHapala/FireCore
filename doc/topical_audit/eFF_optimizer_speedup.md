---
type: TopicalAudit
title: eFF relaxation optimizer — FIRE speedup / stiffness / mass preconditioning study
description: Why eFF relaxation under FIRE is slow or stalls; systematic study of dt, damping, masses, diagonal-curvature preconditioning; discussion of fixed-mass guesses, on-the-fly mass adaptation, and alternative algorithms (BFGS et al.).
tags: [eFF, FIRE, optimizer, preconditioning, masses, stiffness, DynamicOpt, hessian]
timestamp: 2026-11-08
---

# eFF optimizer speedup study

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

## 5. Open questions

- **Physics regression scope:** which published/expected results used coreMode 'f'/'e' with pre-fix (swapped-`fsi`/`fsj`) forces? Frozen-core goes through `evalAE` (was buggy); check whether 'e'/ECP goes through `evalAE_ECP` (was correct) or `evalAE` — needs a per-case look.
- **OpenCL regression:** `eFF.cl` received both the swap fix and the refit coefficients — rerun the CPU↔GPU parity test (`tEFF_ocl` or equivalent) before shipping.
- **Residual ~3e-4 esize f-vs-E mismatch** (`gauss_p8` fs vs `erfx_e6` E): worth removing if BFGS/L-BFGS becomes the production polish — derive `fs` from the same approximant. Also could now re-derive `fs = dy/ds` directly since the refit returns a consistent `dy`.
- Does preconditioning change the *physics* outcome on the production scan set (`distscan/angdistscan_H2O`, bCoreElectrons variants)? Compare basin-by-basin energies vs unit-mass at fixed protocol before adopting.
- How large can `k_ref` get before stiff-DOF descent is starved? (k_ref=30 beat 100 on H2O_spins; crossover somewhere between.)
- BFGS prototype: now converges standalone on both test systems (49/156 iters) — worth porting a bounded dense/L-BFGS + trust-region into `DynamicOpt` for mode-'a' polish; stall-detection still useful for `bPaired=False` collapse cases.
