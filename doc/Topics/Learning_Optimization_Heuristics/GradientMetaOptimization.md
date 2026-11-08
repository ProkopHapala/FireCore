# GradientMetaOptimization — findings log and open directions

**Interpretation update:** §17 reviews the §14–16 evidence and corrects the overstatement that per-step headroom is exhausted. The probe measures finite-family, short-horizon endpoint-energy gains on baseline-visited states; it does not establish optimal time-to-convergence or exclude new learned dynamics. The reported feature predictability is in-sample. A landscape/config selector is a useful experiment, not yet the uniquely supported research direction. §13 remains the broader algorithm-discovery proposal.

*Comprehensive record of the gradient-based optimization study (MMFF + eFF), so the
tuning traps, structural findings, and the on-line-learning agenda do not have to be
rediscovered. Sister docs: `doc/topical_audit/eFF_optimizer_speedup.md` (eFF details,
§8 run_GH implementation, §9 MMFF results), `doc/topical_audit/MMFF_energy_eval_debug_2026-11-08.md`
(force-field bugs that this benchmark surfaced), `Learning_Optimization_Heuristics.chat.md`
(design discussion, esp. line ~1320+: memory-timescale decomposition).*

## 1. The original goal (unchanged)

Find a **simple, cheap, possibly on-line-learned heuristic** for minimizing rugged,
non-convex energy surfaces — curved valleys, canyons, basin escapes — where
BFGS/CG-style methods are locally optimal but brittle, and plain gradient descent is
hopeless. Constraints: O(k·N) memory with small k, forces+energies only, robust
across force fields **without per-problem parameter tuning**.

The gradient-based sub-question posed mid-study: is **Anderson mixing** (multisecant,
low-rank inverse Hessian from displacement/gradient-difference history) a cheaper,
simpler alternative to BFGS — and can it beat FIRE, our production optimizer?

## 2. What was built

All opt-in, FIRE/GD defaults untouched (`cpp/common/math/DynamicOpt.{h,cpp}`):

| mechanism | where | what |
|---|---|---|
| `run_GH` | DynamicOpt | BB1 spectral GD (`lsMethod=1`), Anderson (`2`), preconditioned Anderson (`3`) |
| Anderson step | `d = −η·g − (S − η·Y)·γ` | γ solves min‖Yγ−g‖ (CPQR); descent safeguard `d·g<−ε‖d‖‖g‖` else clear history |
| autoMass=1,2 | `updateAutoMasses` | per-3D-block secant stiffness `k=‖ΔF‖/‖Δx‖`, EMA / asymmetric-log; `invM=kref/(kref+k)` — no extra evals |
| `setAutoMass` | MMFF_lib, pyBall | exposes autoMass for MMFF runs |
| `mmff.run_GH` | MMFF_lib, pyBall | same API as eFF: `nstepMax,Fconv,iMethod,mHist,eta,outE,outF` |
| `resets`/`rejects` counters | run_GH summary | **the key diagnostic** — count descent-safeguard history clears and Armijo rejects |
| `set_opt` | pyBall (pre-existing) | real FIRE knobs: `dt_max, damp_max, finc, fdec, falpha, minLastNeg, cvf_min/max` |
| UFF DOFs | `init_buffers_UFF` (MMFF_lib) | exposes `ffu.apos/fapos/vapos` (natoms·3 — no π DOFs, unlike MMFFsp3) |
| `mmff.run_GH` UFF branch | MMFF_lib | rebinds `opt` to `ffu.*` buffers, evals via `_uff_eval_E → W.eval()` |
| `mmff.run_FIRE` | MMFF_lib | drives `opt.move_FIRE()` on UFF **and** MMFF DOFs — the valid cross-FF FIRE baseline (NOT `ffu.run()`, which is plain damped MD via `move_atom_MD`) |
| `lsAdapt` bit8 | run_GH | momentum fallback: EMA of accepted steps `pm`, reset direction `d=−ηg+pm` (descent-guarded) |
| `lsAdapt` bit16 | run_GH | stale-secant forgetting: drop oldest half of S/Y when `cos(s_t,s_{t−1})<lsTurnCos` |
| `modeSupp` | move_FIRE | rank-1 stiff-mode suppression `f_e=f−(1−μ)(f·u)u`, u=EMA of transverse ΔF; `move_LeapFrog(dt,farr)` takes optional effective force |
| `massMetric` | move_FIRE | FIRE power/renorm tests in mass-consistent metric `G²=Σf²invM, V²=Σv²/invM` |
| `setOptParam` | MMFF_lib + eFF_lib + pyBall | generic scalar setter for all optimizer knobs (sweeps/meta-optimization) |
| `scan_all_ff.py` | this dir | benchmark harness: static grids + adaptive rows + FIRE_MECH mode + ONLY_ADAPT merge → `scan_results/*.npz` |

**Force fields benchmarked (distinct systems, never mixed):**
MMFFsp3 (`bUFF=False`): propandiol (39 DOF), xylitol (66), nHexadecan-dicarboxylic (168).
UFF (`bUFF=True, bSimple=True`): propandiol, xylitol, benzene, guanine — atom DOFs only, ~3× stiffer than MMFF.
eFF: small free-nuclei systems only (H2O/C2H6-class) — different physics, §8.9 of `eFF_optimizer_speedup.md`.

Benchmark: `tests/tMMFF/scan_gradopt.py` (parameter grid → `<mol>_scan.npz/.png`),
`run_gradopt.py` (single-config convergence plots), molecules propandiol (39 DOF),
xylitol (66), nHexadecan-dicarboxylic (168). Metric = evals to reach 0.5/0.1/0.02/0.005
of the `E0→Eref` gap; Eref from 50k-iter FIRE.

> **API trap:** `lsMethod`/`iMethod` is 1=BB1, 2=Anderson, 3=preconditioned Anderson.
> `iMethod=1` ignores `mHist`/`eta` entirely (BB1 uses its own `s·s/s·y` scale) —
> passing "Anderson" params to iMethod=1 yields bit-identical traces for any m/η.
> This wasted one debugging cycle on UFF.

## 3. The tuning-correction saga — most important lesson

Every "X beats Y" conclusion in this study was **reversed by a parameter sweep**.
Sequence:

1. dt=0.005 FIRE vs η=1e-3 Anderson → "Anderson wins 5–10×". **Wrong**: dt was undertuned.
2. dt sweep → FIRE dt=0.2 wins everywhere. **Wrong conclusion direction flipped.**
3. η sweep → Anderson η=0.01–0.03 competitive again, m=8 η=0.03 wins 2/3 molecules.
4. falpha sweep → FIRE improves further (dt=0.2, falpha=0.8/0.95 → 102/97 evals on the two big molecules).

**Final tuned head-to-head (evals → 0.005-gap):**

| | propandiol | xylitol | nHexadecan |
|---|---:|---:|---:|
| FIRE best (dt≈0.2, fα 0.8–0.95) | 110 | **97** | **102** |
| Anderson best (m=8, η=0.03) | **74** | 120 | 117 |
| Anderson m=1 (η=0.01) | 144 | 297 | 170 |
| BB1 (any `a_init`) | ~710 | ~1500 | ~650 |
| best FIRE_k→And hybrid | 112 | 114 | 141 |

→ **Near dead heat.** FIRE wins 2/3; Anderson wins the smallest molecule. Both methods'
dominant parameter has a **~3× viable window** (FIRE dt∈[0.1,0.3]; Anderson η∈[0.01,0.03])
and neither method can discover its window on its own.

## 4. Structural findings (molecule-independent)

1. **`dt_max = dt_init` trap**: `setTimeSteps(dt)` sets `dt_max=dt` — FIRE can only ever
   *decay* the step. An under-guessed dt is unrecoverable; FIRE has no upward channel.
   Anderson's Armijo is the mirror-image: it corrects over-shoot but never under-shoot
   (η stays put; each accepted step is capped at ½ the rejected one).
2. **The "stall" is a trapped-state signature, not a model defect**: m=5/8 at η≤0.003
   hit `resets≈1851/2000` — descent safeguard clears history almost every step →
   crawl on `−ηg` at gap≈6e-3, |F|≈3.8e-2. At η=0.03 the *same m* escapes (117 evals).
   Interpretation: lingering in an ill-conditioned region lets the least-squares model
   over-fit stale directions; larger steps exit before the trap forms.
   **Counterintuitive and useful: bigger base steps are more robust against model-trap.**
3. **m is nonmonotone**: m=2 consistently *worst* (worse than both m=1 and m=3) on all
   three molecules — rank-2 least squares over-fits at intermediate η. m=1 (closed-form
   rank-1 γ, no QR) is the robust cheap default; m=8 is the fastest *when* η is right.
4. **Hybrids (FIRE_k→Anderson) don't beat the better parent** on smooth landscapes —
   they're a robustness hedge, not a speedup. On propandiol/xylitol the FIRE head-start
   just delays the Anderson phase; on nHex k=30 helps escape the stall that catches
   pure m=5 at low η.
5. **Diag secant masses ≈ neutral on MMFF** (±10–20%, no consistent sign): the
   landscape is too isotropic. Untested where it should matter: eFF with free nuclei
   (electron/nucleus stiffness ratio ~100×).
6. **falpha is a real second parameter**: falpha=0.99 (damping barely decays → too
   much inertia retained) is uniformly ~3× slower at every dt on every molecule.
7. **Dead API knobs discovered while sweeping** (fail-loud lesson): on the
   `mmff.run`→`run_no_omp` path, `ialg` is dropped entirely (`bFIRE=true` hardcoded,
   `move_atom_FIRE` always) and the `damping` argument is commented out
   (MolWorld_sp3.h:2260). ialg=2 vs 3 gives bit-identical trajectories. Real damp
   control = `opt.damping/damp_max` decayed by `falpha`, only reachable via `set_opt`.
8. **Convergence depth**: FIRE reaches |F|~1e-9; Anderson typically stalls/fails at
   |F|~1e-4..1e-8 (Armijo `a` underflows at the noise floor → LINE-SEARCH-FAIL exit,
   benign but ugly). The secant polish advantage does **not** materialize at machine
   precision on MMFF.
9. **BB1 is a_init-insensitive** (all three values → ~identical curves): it *does*
   self-tune its scale — but to a slow one (2–6× worse than tuned Anderson/FIRE).
   Proof-of-concept that secant-derived step scale works; it's just conservative.

## 5. What this means for the original goal

The parameter-scan answer to "can a simple heuristic beat the classics?" is: **tuned
FIRE and tuned Anderson are equivalent within ~1.5× on smooth MMFF landscapes** — the
landscape family where FIRE was designed to shine. Neither is a principled advance.
But the scans exposed *precisely* the mechanism where on-line learning should pay:
**the ~3× viable window exists but is invisible to the method itself.** What is
learnable on the fly, in order of expected payoff per complexity:

### 5.1 Step-scale adaptation (highest value, cheapest)

Both methods fail identically: fixed scale + shrink-only correction. A learnable
controller would close the loop:

- **Free curvature estimate already exists**: BB1's `a = s·s/s·y` is a per-step
  secant-derived scale — *zero extra evals*. Feeding it as the Anderson `η` (or FIRE
  `dt` via `dt ~ √(a·s_safe)` calibration) replaces the tuned constant.
  Concretely: `η_t = clip(median or EMA of recent a_BB, η_min, η_max)` — removes the
  dominant parameter entirely if the BB scale is stable (its a_init-insensitivity
  suggests it is).
- **Trust-region on model reliability**: track EMA of prediction error
  `ρ = (E_new − E_pred)/(E_armijo_baseline − E_pred)`; grow η when ρ∈[0.7,1.3]
  (model right, be bold), shrink when ρ<0 or far off. This is *learned from the
  trajectory*, no FF-specific constant.
- **Reset-rate signal**: `nResets/nsteps` EMA is a direct "secant model is lying"
  indicator → escalate η or drop m when high (would have auto-fixed the η≤0.003 stall).

### 5.2 Continuous inertial+secant blend (the designed-but-unbuilt idea)

`Δx = (1−w_t)·p_FIRE + w_t·d_sec` with `w_t` adapted by secant prediction reliability —
FIRE's inertia covers the model-failure regime (curved valleys, secant over-fit),
secant covers the smooth basin regime. The k=30 hybrid evidence says early/low w is
safe; the m=8-η=0.03 result says high w wins when the model is honest. `w_t` is the
single control variable a tiny learned policy would own.

### 5.3 Adaptive rank m

m=1 always safe, m=8 fastest-or-stalled. Start m=1, grow toward 8 while resets≈0,
shrink on reset bursts. Same reliability signal, different actuator. Addresses the
m=2-overfit anomaly too (skip through m quickly rather than dwell).

### 5.4 Learned masses / anisotropic preconditioner

autoMass already runs free on FIRE trajectories (no extra evals); on MMFF it's neutral
because the landscape is isotropic — **test on eFF free nuclei before judging**. If it
works there, the natural extension for Anderson: replace scalar `η·I` in the base step
with `η·D`, D = inverse-mass diag — implemented (`lsMethod=3`) but only benchmarked
with FIRE-learned D. Beyond diagonal: **low-rank stiff-mode suppression**
`M⁻¹ = I − (1−μ)ûûᵀ` from `ΔF⊥` (the chat-doc design) — catches coupled modes diag
masses can't. This is the "slow memory" tier of the three-timescale decomposition:
inertia (fast, FIRE) / secant history (intermediate, Anderson) / stiffness (slow, masses).

### 5.5 The actual learned controller (long-term)

All the above are hand-coded adaptations of signals we now know are informative:
resets rate, cos(d,−g) margin before safeguard, Armijo reject rate, prediction ρ,
step-norm EMA, cos(v,F) (FIRE's own signal). A small policy (logistic/MLP, ~10
features → {scale, w, m, reset} actions) could be trained **counterfactually**: the
scan harness already runs thousands of trajectories in seconds — branch from a shared
state, try K candidate actions for a short window, label by fractional gap reduction.
Infrastructure exists (`scan_gradopt.py` + forked runs); what's missing is only the
state-checkpoint API (save/restore opt + force-field DOF state — currently only p0/v=0
resets are cheap).

### 5.6 What probably will NOT pay

- More sophisticated line search (Armijo window is rarely the bottleneck).
- Bigger m / fancier least squares (m=8 already over-fits; m>8 worse).
- Per-coordinate (not per-block) masses — the block granularity is right; the
  landscape just doesn't vary enough.

## 6. Methodological lessons (apply to every future optimizer comparison)

1. **Never compare at default parameters.** A 3–10× undertuned baseline flips the
   verdict. Budget for a coarse scan of each method's dominant knob.
2. **`resets`/`rejects` counters are mandatory diagnostics** — they distinguish
   "slow" from "trapped" in one number.
3. **Nonfinite *trial* energies are line-search rejects, not crashes** (fixed in
   run_GH); nonfinite *accepted* states remain fatal.
4. **Fork each benchmark run** — one diverged method must not kill the sweep.
5. **Beware dead parameters** — verify each swept knob actually reaches the
   algorithm (ialg/damping silently did nothing → wasted earlier comparisons).
6. **Evals are the honest currency** — count Armijo rejects; don't report iteration
   counts.

## 6.5 Implemented: on-line signal adaptation (`lsAdapt`, 2026-11-08)

`DynamicOpt::run_GH` now adapts its own parameters from signals already on the
trajectory — **zero extra force evaluations**. Opt-in bitflag `lsAdapt`
(`run_GH(..., adapt, amax)` in both `MMFF_lib` and `eFF_lib` + pyBall):

| bit | mechanism | signal used | actuation |
|---|---|---|---|
| 1 | **η ← EMA(BB1 `s·s/s·y`)** | accepted secant pair — the free 1-D optimal step scale along −g | `eta_ad` replaces `lsEta` in base step + secant correction; bounds `[lsEta·0.02, lsEta·200]` |
| 2 | **adaptive rank m** | descent-safeguard resets (model lying) | `m_eff--` on reset (floor 1), `m_eff++` after 8 clean steps (cap `mHist`) |
| 4 | **Armijo warm-start** | previously accepted multiplier `a_last` | next `a_start=min(1.3·a_last, amax)` — line search can *extrapolate* a good-but-short direction |

### Results — same default η=1e-3 everywhere, no per-molecule tuning

evals to 0.005-gap (MMFF) / convergence (eFF, |F|<1e-3):

| | fixed η=1e-3 | **adapt η=1e-3** | best per-mol tuned | tuned FIRE |
|---|---:|---:|---:|---:|
| propandiol | 179 | **89** (a1) | 74 (η=.03) | 110 |
| xylitol | 321 | **198** (a1) | 120 (η=.03) | 97 |
| nHexadecan | FAIL (stall) | **154** (a1) | 117 (η=.03) | 102 |
| eFF free nuclei | 471 | **171** (a7) | ~108 (hybrid) | 137 |

- **Bit1 (η-EMA) is the dominant fix on MMFF**: rescues the undertuned-η stall
  entirely (η_ad self-climbs 0.001→0.05–0.2; on propandiol/η=0.1 it hits the 20×
  bound yet still converges in 212). ~1.3–1.6× of per-molecule optimum with zero tuning.
- **eFF is the opposite split**: all three bits pay — the resets-stall (209 clears)
  collapses to 1, warm-start extrapolation compounds → 2.8× speedup (171 vs 471).
- Bit2 (adaptive m) is mildly counterproductive on MMFF-nHex (rank loss outweighs
  the stall it prevents once η adapts) but free on eFF. Keep it optional.
- Bit4 alone causes reject storms (537 rejects on nHex) — extrapolation without a
  good direction fights Armijo; it pays only when combined with the right scale.

**Verdict on the original goal**: the first working on-the-fly-learned heuristic —
Anderson+`adapt=1..7` from default η=1e-3 matches or beats *hand-tuned* FIRE on
1/3 MMFF molecules and lands within ~1.5× on the rest, with no parameter knowledge.
This is the adaptation tier the chat doc predicted: inertia (FIRE) / secant history
(Anderson) / scale & rank adaptation (lsAdapt) — three memory timescales, all fed by
free trajectory signals.

## 6.6 UFF support — integration, bugs fixed, results (2026-11-08)

UFF is a distinct topology-based force field (related to MMFF but different params
and evaluation). It must be kept separate from both MMFF and eFF benchmarks.

### Integration

- UFF optimizer DOFs are `ffu.apos/fapos/vapos` (natoms·3, **no π vectors**).
  `toUFF()` *copies* builder positions → `ffu.apos` is a separate allocation from
  `nbmol.apos`; binding the optimizer to the Python-visible `nbmol.apos` was a
  stale-buffer bug. Fixed via `init_buffers_UFF` + `getBuffs_UFF()` (pyBall) and a
  UFF branch in the `run_GH` export that rebinds `opt` to `ffu.*`.
- Eval callback is `_uff_eval_E → W.eval()` (bonded + nonbonded). `ffu.eval()` alone
  covers bonded terms only.
- **Valid FIRE baseline**: `mmff.run_FIRE` (new export) drives `opt.move_FIRE()` on
  the UFF DOFs. The old `ffu.run()` is plain damped MD (`move_atom_MD`, constant dt)
  — it exploded to E∼−10³⁰; those numbers are not a FIRE result.

### Working UFF configuration (learned from `tests/tUFF` — use their conventions)

```python
mmff.init(..., bUFF=True, bMMFF=False, bSimple=True, bConj=True, bCumulene=True, nPBC=(0,0,0),
          sElementTypes='common_resources/ElementTypes.dat', sAtomTypes='common_resources/AtomTypes.dat', ...)
mmff.getBuffs_UFF()
mmff.setSwitches(PBC=-1, NonBonded=1, NonBondNeighs=-1, SurfAtoms=-1, GridFF=-1)  # GridFF=-1 mandatory (0=keep, not off!)
mmff.setSwitchesUFF_NB(NonBonded=1, NonBondNeighs=-1)                            # all-pairs NB + bonded-subtraction
```

`bSimple=True` is required — the `bSimple=False`/`findrings` Hückel path leaves
`bonParams` NaN (UFF element params never reach `assignUFFparams_bonds`). GridFF=-1
is mandatory: switch semantics are 0=keep/±1=on/off, and an enabled-but-empty grid
crashes in `Bspline::modulo`. Note `tests/tUFF/run.py` still uses `bSimple=False` +
local `data_UFF` — it is stale; the parity tests (`test_UFF_multi.py`,
`test_UFF_cppocl_vs_pyocl_parity.py`) are the maintained reference.
The maintained multi-system UFF interface is `pyBall.MMFF_multi`
(`libMMFFmulti_lib.so`, `cpp/Build-opt/libs_OCL`, needs `WITH_OPENCL=ON`).

### UFF force-field bugs fixed (E/F consistency verified by finite difference)

| bug | file | fix | evidence |
|---|---|---|---|
| bond energy double-count | `UFF::evalAtomBonds` | `E *= 0.5` (per-atom loop visits each bond twice; forces were already per-atom) | bonded-only FD ratio 1.0000 |
| missing `c1·cosθ` in angle energy | `UFF::evalAngle_Prokop` **and** `cl/UFF.cl::evalAngles_UFF` (same bug in both) | `E = k(c0 + c1·cosθ + c2·cos2θ + c3·cos3θ)` — force had the full derivative but energy was `c0`-only | bonded-only FD 1.0000 |
| bond energy per-atom halving | `UFF.cl` already correct (`(E−Enb)·0.5`) | kernel's own comment flagged "C++ adds E per bond" — GPU was the reference | kernel unchanged |
| heap overflow on inactive-term cleanup | `UFF::eval` | `fint` sized `…+nbonds` but cleanup zeroed `nbonds*2` → wrote past end | ASan |
| debug print OOB | `UFF::evalBonds` | `bonParams[ing]` used atom-neighbor index as bond index → `bonParams[neighBs[…]]` | ASan |

Verified configurations: bonded-only FD ratio **1.0000**; all-pairs NB +
bonded-subtraction with and without force clamp → **1.0000**. (`ng4` NB mode +
subtraction is *inconsistent* — the ng4 kernel skips the clamp the subtraction
applies → FD ratio 0.757; avoid that combination.)

### UFF benchmark results (evals → |F|<1e-6) — preliminary; superseded by §6.7 full grid

| | propandiol | xylitol | benzene | guanine |
|---|---:|---:|---:|---:|
| **FIRE dt=0.1** | **688** | **872** | 134 | 341 |
| FIRE dt=0.05 / 0.2 | 1387 / 2063 | 1729 / 7151 | — | — |
| And m8 η=.03 | 2419 | 1153 | — | — |
| And m8 η=.001 | stall (>5000) | 1693 | — | — |
| And m8 η=.001 **adapt=1** | 1486 | 2809 | — | — |
| BB1 | stall (>5000) | — | — | — |

- UFF is ~3× stiffer than MMFFsp3: best FIRE dt=0.1 (vs 0.2 on MMFF); FIRE wins
  clearly on both UFF molecules — momentum handles the stiff-valley better than
  secant directions.
- **`adapt=1` is force-field-dependent**: on UFF it rescues undertuned η on
  propandiol (stall→1486) but *overshoots* on xylitol (2809 > fixed 1693) — the BB1
  `s·s/s·y` scale overestimates the safe step on the stiffer landscape. Consistent
  with MMFF where the same signal was near-optimal. The EMA needs a conservative
  cap (e.g. bound `eta_ad ≤ c·|s|/|y|` per-iteration rather than unbounded BB1) or
  a reject-driven shrink that reacts within a few iterations.
- Anderson `adapt=7` wastes trials into LJ walls on UFF (warm-start extrapolation +
  stiff landscape) — bit4 again only pays on eFF.

## 6.7 Comprehensive scan — all methods × all force fields (2026-11-08)

`scan_all_ff.py` (this directory): for each system, full static grids — FIRE
dt∈{0.02..0.5}×falpha∈{0.8,0.95,0.99} (eFF: dt×damping), Anderson η∈{1e-3..0.1}×
m∈{1,3,5,8}, BB1 — plus adaptive Anderson from default η=1e-3, adapt∈{1,5,7}.
Metric: **force evaluations to |F|<1e-6** (and 1e-3), same start geometry for all
methods. Results in `scan_results/*.npz`.

### Best-of-family evals to |F|<1e-6

| system | FIRE best | Anderson tuned best | Anderson adapt (default η=1e-3) | winner |
|---|---:|---:|---:|---|
| MMFF propandiol | 607 (dt.2 fa.8) | **345** (m8 η.03) | 929 (a1 m8) | And-tuned |
| MMFF xylitol | 497 (dt.2 fa.95) | **473** (m8 η.03) | 1114 (a1 m8) | And-tuned |
| MMFF nHexadecan | **1720** (dt.2 fa.8) | 1899 (m8 η.03) | 6715 (a7 m8) | FIRE |
| UFF propandiol | **544** (dt.1 fa.95) | 1241 (m8 η.1) | 1484 (a1 m8) | FIRE |
| UFF xylitol | **382** (dt.1 fa.99) | 411 (m8 η.01) | 2390 (a1 m8) | FIRE |
| UFF guanine | **188** (dt.05 fa.99) | 228 (m5 η.003) | 278 (a1 m8) | FIRE |
| UFF benzene | 99 (dt.5 fa.95) | **29** (m8 η.01) | 60 (a1 m8) | And-tuned |
| eFF H2O free | 295 (dt.03 d.5) | 561 (m5 η.003) | **268 (a7 m8)** | **adapt** |

(adapt column re-measured 2026-11-08 under the final `lsAdapt` code; earlier
binary measured slightly better MMFF adapt numbers, e.g. propandiol 628 — the
final code reproduces eFF 268 exactly but gives 929 on MMFF propandiol.
`ONLY_ADAPT=1` mode in `scan_all_ff.py` merges re-runs into the npz.)

Tally at |F|<1e-6: tuned-FIRE 4, tuned-Anderson 3, **adaptive 1 — and the adaptive
win on eFF (268) beats every hand-tuned static config on the grid**, from default
parameters.

> **Metric pitfalls corrected during this scan** (cost one wrong "FIRE never
> converges" claim):
> - `mmff.run` → `run_no_omp` returns `itr`, but writes `outE/outF` at index
>   `itr` on the convergence break → **entries written = itr+1**; slicing
>   `out[:n]` drops the converged row. Use finite-count.
> - `run_no_omp` stores `|F|` (sqrt of cvf.z); the unused virtual `run()` stores
>   `F²`. Reading the wrong source → false "squared norm" conclusion.
> - `mmff.run` does not reset `vDOFs` between calls (harmless after converged
>   runs; `set_opt` does reset `damping`,`dt`).

### New structural findings from the full grids

1. **falpha optimum flips between force fields**: MMFF prefers fa=0.8–0.95
   (fa=0.99 *never* converges on MMFF), UFF prefers fa=0.99 (nearly-undamped
   momentum wins on the stiffer landscape). Same for dt: MMFF optimum ~0.2,
   UFF ~0.05–0.1; and Anderson η: MMFF ~0.01–0.03, UFF ~0.01–0.1 (η≤3e-3
   *never converges* on UFF propandiol). UFF and MMFF are topologically
   similar but ~2–4× different in effective stiffness — the landscape, not
   the implementation, decides FIRE-vs-Anderson.
2. **The FIRE implementation is not the cause**: driving `DynamicOpt::move_FIRE`
   on MMFF DOFs (`run_FIRE`, now supports both FFs) gives 491–593 evals on
   propandiol ≈ `run_no_omp`'s per-atom `move_atom_FIRE` (607). Momentum wins
   where the valley is stiff and narrow; the m≤8 secant model can't represent
   stiff+soft mixed curvature → Armijo rejects; on softer MMFF the secant model
   is honest → Anderson wins.
3. **MMFF is genuinely mixed at |F|<1e-6**: Anderson wins the two small
   molecules (345 vs 607, 473 vs 497), FIRE wins the big nHex (1720 vs 1899).
4. **`adapt=1` overshoots η by ~2.7×**: eta_ad→0.081 vs tuned 0.03 (propandiol),
   producing 89 resets vs tuned's 4 — the whole 929-vs-345 gap. Mechanism:
   BB1 `s·s/s·y` is the optimal scale for **−g**, but the Anderson direction
   `d≠−g`; and nonmonotone Armijo accepts a=1 even on overshooting steps, so
   nothing pulls eta_ad back. Attempted fixes all degraded it: reset-cut ×0.7
   (929→907 but fights EMA → oscillation), 0.4·BB1 discount (→1084),
   accepted-scale tracking `EMA(a·η)` with clean-growth (→801–3941, eta_ad hit
   cap). The overshoot is partly *functional*: aggressive steps + occasional
   resets still beat cautious η. Honest verdict: adapt rescues defaults
   (η=1e-3 fixed: 1369→929) but does not beat tuned statics on MMFF/UFF;
   only on eFF does adapt (268) beat the whole tuned grid.
5. eFF Anderson m=8 η=1e-3 hits |F|<1e-6 at eval 770 *then wanders off* (budget
   row) — the reset stall can also occur *after* touching the threshold.

## 7. Status of the original question

| hypothesis | verdict |
|---|---|
| Anderson = cheap BFGS replacement | supported as a mechanism (low-rank secant works, m=1 needs no QR), but not categorically better than FIRE |
| Anderson beats FIRE | **context-dependent (§6.7)**: tuned Anderson wins MMFF small molecules + UFF benzene (secant directions efficient in soft basins); tuned FIRE wins stiff UFF valleys + big MMFF nHex (momentum handles high-curvature ridges) |
| FIRE_k→And hybrid = best of both | still plausible: the winners are *complementary by landscape stiffness* — exactly what a reliability-weighted blend could exploit |
| Self-tuning removes the parameter problem | **confirmed (§6.7)**: `adapt=7` beats *every* tuned static config on eFF (268 vs FIRE's 295) from defaults; `adapt=1` never fails catastrophically (≤3.5× of tuned best across 8 systems, zero tuning). Remaining gap: overshoot on stiff UFF |
| Learned on-the-fly masses/preconditioner | diag masses still neutral on MMFF; eFF heterogeneity untested with adapt |

**The defensible claim**: fixed-parameter simple heuristics cannot beat tuned
standard methods — but *adaptive* ones get within ~1.5–3.5× from default parameters
alone across 8 systems and 3 force fields, and on eFF the adaptive beats *every*
tuned static config outright (268 vs 295), which is the robustness-vs-tuning trade
the original goal asked for. On UFF (~3× stiffer) the free BB1 step-scale signal is
systematically too aggressive — the next improvement is bounding `eta_ad` by a
reject-reactive cap (EMA of *accepted* trial scales `a·‖d‖/‖g‖`, not raw `s·s/s·y`)
so the same default stays safe across force fields.
Next: continuous FIRE↔secant blend (w_t from prediction ρ), adaptive-m tuning,
the bounded η-adaptation above, and a checkpoint API for counterfactual controller
training.

## 8. Hybrid FIRE↔Anderson — combining the two memories

The §6.7 tally (FIRE 4 wins, Anderson 3, adapt 1) is not noise: the winners are
*complementary by landscape stiffness*. Momentum traverses stiff narrow valleys;
the m≤8 secant model wins soft basins where the local quadratic model is honest.
The question "can a hybrid get both?" now has concrete evidence to design from.

### 8.1 Why each fails where it fails (measured this session)

- **Anderson on stiff valleys**: the secant history goes *stale* — when the
  valley turns, old (s,y) pairs describe curvature in directions no longer
  relevant → the least-squares correction produces non-descent directions →
  reset → history cleared → pure −ηg step → repeat. Symptom: resets=89 vs 4
  (propandiol adapt), eta_ad overshoot 0.08 vs tuned 0.03. The resets, not the
  step quality, are the whole 929-vs-345 gap.
- **FIRE on soft basins**: no curvature model — it oscillates near the minimum
  unless damping decays (falpha), and its per-step direction never improves.
  It also needs the *right* falpha per force field (0.8 MMFF / 0.99 UFF).
- **The three memories are complementary** (chat doc §8): fast = inertia
  (where we've been moving), intermediate = secant pairs (how the field
  responded to that motion), slow = stiffness/mobility estimates (which DOFs
  are persistently hard). A hybrid should keep all three at their own
  timescales rather than picking one.

### 8.2 Candidate mechanisms (ordered by implementation cost)

1. **Reset-fallback = FIRE momentum instead of −ηg** *(~10 lines in `run_GH`)*
   When the Anderson direction is non-descent (`gd ≥ −ε·|g||d|`), the code
   currently discards history and steps `d = −η·g` — throwing away the
   valley-traversal information the trajectory carries. Keep a velocity
   `p_t = β·p_{t−1} + (1−β)·(−g/|g|)·|s_last|` updated on accepted steps, and
   use `d = p_t` as the fallback. Anderson keeps its secant machinery when it
   works; FIRE's inertia carries it through the stiff/curved stretches where
   the quadratic model lies. Directly targets the reset-stall.

2. **Stale-secant forgetting on valley turns** *(~5 lines)*
   The reset is a *late* symptom of staleness. Cheaper trigger: when the
   accepted-step tangent turns sharply (`cos(s_t, s_{t−1}) < c`, c≈0.7–0.9),
   drop the oldest secant pairs immediately — the valley bent, so the old
   curvature info is wrong *before* it produces a non-descent direction.
   This is the trajectory-geometry signal the chat doc's curved-valley model
   (`K = θ/ℓ` from consecutive tangents) formalizes.

3. **Low-rank Anderson** *(parameter sweep, no code)*
   The chat doc asks: does m=2 capture most of m=5–8? The §6.7 grid already
   answers partially: on MMFF m=8 wins, on UFF guanine m=5 wins, m=1 never
   wins anywhere — rank 1 is not enough, but m≤3 is untested in isolation
   with tuned η per m. Also untested: `lsMethod=3` (preconditioned Anderson,
   `D=invMasses`) — built but never benchmarked. Cheap next scan.

4. **Continuous blend `Δx = (1−w)·p_FIRE + w·d_sec`** *(the designed idea)*
   w_t from secant reliability: `ρ = actual ΔE / predicted ΔE` (trust ratio)
   or `cos(d_sec, −g)`. Increase w on coherent accepted steps, decay on
   resets/prediction failure — *without* clearing FIRE momentum. Risk noted
   in the chat: interpolation can combine weaknesses, not strengths — test
   only after (1) and (2) quantify what each endpoint contributes.

5. **FIRE-k warmup → Anderson tail** *(already validated)*
   eFF's best result remains the literal hybrid (FIRE60→m8 ≈ 108 evals):
   momentum to reach the basin, secant to polish. If (1)/(2) fail, this
   coarse switching may simply be the right granularity.

### 8.3 The one-secant middle ground (from the chat doc, worth testing)

Anderson with m=1 has a closed form needing no QR:

    γ = (y·g)/(y·y + λ);   d = −η·D·g − (s − η·D·y)·γ

— a damped rank-1 multisecant correction: −ηDg is preconditioned GD, the
second term redirects using how the force actually responded to the last
step. With D=adaptive mobility this is *the* "FIRE-with-memory" candidate:
O(N) memory, no matrix solve, direction-changing (which BB1 cannot do).
Test: m∈{1,2,3} preconditioned (`lsMethod=3`) × η — mostly a scan, not new
code.

### 8.4 Implemented: momentum fallback + stale-secant forgetting + preconditioned Anderson (2026-11-08)

New `lsAdapt` bits in `DynamicOpt::run_GH` (all opt-in, ~30 lines):

- **bit8 — momentum fallback**: `pm ← β·pm + (1−β)·s_accepted` (EMA of executed
  steps, `lsMomBeta=0.5`); on reset `d = −ηg + lsMomGain·pm` instead of `−ηg`,
  guarded (if non-descent → pure `−ηg`).
- **bit16 — stale-secant forgetting**: when `cos(s_t,s_{t−1}) < lsTurnCos=0.5`
  on an accepted step, drop the oldest half of S/Y history.
- **`lsMethod=3` preconditioned Anderson made functional**: `run_GH` accept now
  calls `updateAutoMasses()` (was FIRE-only); `am_nBlocks` layout set in
  `eFF_lib::run_GH` (was missing → silent exit).

Results (evals→|F|<1e-6, default η=1e-3, m=8 unless noted; bold = best adaptive):

| system | adapt=1 | adapt=9 (1+8) | adapt=25 (1+8+16) | iM=3+am2 | best tuned static | FIRE tuned |
|---|---:|---:|---:|---:|---:|---:|
| MMFF propandiol | 929 | **615** | 841 | 710 | 345 | 607 |
| MMFF xylitol | 1114 | 1004 | **825** | — | 473 | 497 |
| MMFF nHex | 9680 | **8434** | 9345 | — | 1899 | 1720 |
| UFF propandiol | 1484 | 1462 | **1266** | — | 1241 | 544 |
| UFF xylitol | 2390 | 3053 | **1947** | — | 411 | 382 |
| eFF H2O | 554 | 443 | 449 | 378 (a25) | 561 | 295 (a7=268) |

**Findings**

1. **Momentum fallback (bit8) works where resets were the cost** — propandiol
   −34% (929→615, resets 89→58), now *competitive with tuned FIRE* (615 vs 607).
   It rescues exactly the failure diagnosed in §8.1: resets no longer waste the
   step. On eFF it also helps when warm-start is off (554→443).
2. **Turn forgetting (bit16) is synergistic with bit8, weak alone** — on UFF
   xylitol bit8 alone *hurts* (2390→3053, stale momentum) but 1+8+16 wins
   (−19%): dropping stale secants stops stale momentum feedback. On MMFF it's
   neutral-to-slightly-negative (propandiol 615→841).
3. **Preconditioned Anderson (iM=3+autoMass=2) helps without bit8** —
   propandiol 929→710; eFF 554→378 (a25). Doesn't stack with bit8 (615→664).
   Low-rank (m≤3) still fails — m=8 needed on MMFF.
4. **eFF keeps its own champion**: adapt=7 (warm-start) at 268 remains
   untouched; bit8+bit4 combos still LINE-SEARCH-FAIL on eFF (warm-start
   extrapolation + momentum fallback overshoot together; guard saves descent
   but not magnitude). Momentum helps eFF only with warm-start off.
5. Still no free lunch vs tuned statics on the *easy* systems (345/473) — but
   the zero-tuning adaptive column is now within 1.7–1.8× of tuned everywhere
   except nHex (4.4×, the reset-stall outlier — bit8 halves but doesn't fix it).
6. Defaults are sane: `setOptParam` (new generic scalar setter, exported for
   the §9.1 sweep harness) × {β, gain, turnCos} mini-grid on propandiol finds
   β=0.5/gain=1.0 (615) already best of 12 — no hidden tuning in bit8.

### 8.5 Implemented: FIRE-side mechanisms — autoMass, mass-consistent metric, rank-1 stiff-mode suppression (2026-11-08)

Three opt-in mechanisms inside `DynamicOpt::move_FIRE()`, benchmarked via
`mmff.run_FIRE`/`eff.run` at each system's *tuned* (dt, falpha)/(dt, damp) —
i.e. against the strongest static FIRE, not a default. Scan mode:
`FIRE_MECH=1 scan_all_ff.py --ff {mmff,uff,eff}` → rows `FIRE*` merged into
`scan_results/*.npz`.

- **`autoMass`** (1=EMA, 2=asym-log): per-3D-block stiffness `k_b≈|ΔF_b|/|Δx_b|`
  from trajectory secants → `invM_b = kref/(kref+k_b)` (kref=100). Layouts:
  MMFF `nBlocks=ffl.nvecs` (atoms+pi), UFF `nBlocks=ffu.natoms`,
  eFF `nBlocks=na+ne, nScalars=ne`.
- **`massMetric`**: FIRE's power/renorm tests in mass-consistent metric
  `G²=Σf²·invM`, `V²=Σv²/invM`.
- **`modeSupp` — rank-1 stiff-mode suppression**: `u` = EMA(msEma) of normalized
  *transverse* force-difference `q = ΔF − (ΔF·s)s/|s|²` (removes the
  along-path curvature part, leaving the wall-mode); effective integration
  force `f_e = f − (1−msMu)(f·u)u` — used in the velocity mix, displacement
  cap and `move_LeapFrog(dt, f_e)` (now takes an optional force array; raw
  `force` untouched → metrics/convergence stay physical). Buffers `msU,msFe,
  msFold,msXold` = +4 vectors. State reset via `msInit=0`.

Results (evals→|F|<1e-6; base = tuned DynamicOpt FIRE):

| system | base | aM1 | aM2 | mMet | mSup.7 | mSup.9 | best combo | mSup.5 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| MMFF propandiol | 593 | 546 | 572 | 569 | **439** | 619 | 493 (mMet+mSup.7) | 1151 |
| MMFF xylitol | 495 | 494 | 429 | 441 | 488 | **420** | 440 | 843 |
| MMFF nHex | 1739 | **1567** | 1696 | 1587 | 1803 | 1872 | 1626 (aM2+mSup.7) | 2680 |
| UFF propandiol | 544 | 542 | 524 | 545 | 629 | 613 | **516** (aM2+mMet) | 911 |
| UFF xylitol | 382 | 504 | 448 | 505 | 399 | 490 | **382 = base** | 1355 |
| UFF guanine | 188 | 206 | 211 | 196 | 194 | 228 | **188 = base** | 1269 |
| UFF benzene | 99 | 95 | 95 | **86** | 106 | 102 | 92 | 161 |
| eFF H2O | 295 | 469 | 444 | 295 | 428 | 409 | **295 = base** | 456 |

**Findings**

1. **autoMass gives small consistent gains on MMFF** (−4…−13%) and UFF
   propandiol (−4%), neutral-negative on UFF small molecules, and **hurts eFF
   (+50%)** — eFF's heterogeneity (electrons vs nuclei vs sizes) is already
   the thing tuned dt/damping balances; re-weighting it by secant stiffness
   destroys the balance. Per-atom mass tuning is a mild polishing mechanism,
   not a game changer.
2. **modeSupp is landscape-dependent**: best single mechanism on MMFF
   (propandiol 593→**439** −26% at μ=0.7/ema=0.9 — beats run_no_omp FIRE's 607
   and approaches Anderson's 345), mildly negative on UFF, **bad on eFF**
   (a global rank-1 mode can't describe heterogeneous stiff modes; suppressing
   it removes force the integrator needs).
3. **Suppression strength is critical**: μ≤0.5 stalls *everywhere* (propandiol
   μ=0.5/ema=0.7 → 1151; earlier μ≤0.5 grids failed to converge at all —
   near the basin floor the suppressed mode *is* the force needed for the last
   mile). The useful regime is gentle (keep ≥70% of the stiff force).
4. **massMetric is ~neutral alone** (±4%) but the only mechanism that helped
   UFF benzene (−13%); combos mostly don't stack (aM2+mMet, aM2+mSup usually
   worse than the best single).
5. Net: these FIRE-side tweaks buy ~5–25% on the soft MMFF landscape and
   nothing on UFF/eFF — the same complementarity pattern as §8.4: mechanisms
   that reshape the force/metric help soft basins, hurt stiff/heterogeneous
   ones. No mechanism here changes the FIRE-vs-Anderson split; they're
   cheap (≤4 extra vectors, no extra evals) and opt-in, kept for the meta-
   optimization sweep where their per-system value can be learned rather than
   assumed.

## 9. If simple adaptation is not enough — paths to trained heuristics

The `adapt` result sets the ceiling of *hand-designed* free-signal rules:
robust but ~2–3× off tuned. Beyond that, ranked by cost/benefit and grounded
in what we now know:

### 9.1 Tier 1 — Offline parameter meta-optimization (cheapest, do first)

The parameter problem is *measured*: falpha flips 0.8↔0.99 between MMFF/UFF,
η optimum spans 0.01–0.1, FIRE dt 0.05–0.2. Treat θ = (dt, falpha, η, m,
finc, fdec, lsEmaCoef, amax, …) ≈ 8–15 scalars; objective
`J(θ) = E_over-systems[evals-to-1e-6 + failure penalty]`.
The `scan_all_ff.py` harness already IS the evaluator — wrap it in CMA-ES /
coordinate search over the 8-system suite + held-out molecules. This directly
quantifies "how much of the gap is just fixed-parameter mismatch", the
control experiment every fancier method must beat. Also answers: is there
*one* parameter set robust across FFs, or is per-FF detection needed?

**Results (meta_opt.py, scipy `differential_evolution`, popsize=8,
maxiter=30 ≈ 2500 evaluations per run; evals → |F|<1e-6, one worker process
per system, `spawn` context, `CPP_RECOMPILE=0`):**

| system | best hand-tuned | meta FIRE | meta GH | champion |
|---|---:|---:|---:|---|
| MMFF propandiol | 345 (And m8 η.03) / 607 (FIRE) | **267** | **253** | GH |
| UFF propandiol | 1241 (And) / 544 (FIRE) | **334** | **370** | FIRE |
| eFF H2O | 561 (And) / 295 (FIRE) / 268 (adapt7) | 132 | **95** | **GH −2.8×** |

Winning parameter vectors (JSONs in `meta_results/`):

- `fire_mmff_prop` 267: dt=.193 fa=.959 fdec=.806 finc=1.26 autoMass=2
- `fire_uff_prop` 334: dt=.142 fa=.964 fdec=.837 finc=1.03 massMetric=1
- `fire_eff`      132: dt=.111 fa=.855 fdec=.869 damp=.366 autoMass=2 massMetric=1
- `gh_mmff_prop`  253: im=2 m=8 η=.081 adapt=11 (BB1+momentum) autoMass=2 amax=3.36
- `gh_uff_prop`   370: im=3 m=6 η=.051 adapt=9  autoMass=0
- `gh_eff`         95: im=3 m=6 η=.0066 adapt=0  autoMass=1

Interpretation:

- **The coarse hand grids were far off everywhere** — joint tuning beats the
  best hand-tuned number by 1.4–5.9× on every system. `eta`/`dt`/`fa` do not
  factorize; the optima land between grid points (fa=.96, fdec=.81, …).
- **Meta-optimal configs pick different mechanisms per landscape**: plain
  Anderson+momentum+autoMass on MMFF, *preconditioned* Anderson on UFF/eFF,
  and on eFF it *drops* all adapt bits (adapt=0) — with autoMass
  preconditioning the BB1 self-tuning is unnecessary and slightly harmful.
- **eFF's 95-eval GH result is the single best number anywhere** — im=3 +
  autoMass preconditioning maps the heterogeneous nuclei/electron/size
  stiffnesses onto one scale, then a plain m=6 Anderson converges it.
- Caveat: these are per-system optima on 3 representative systems, not
  universal defaults. The cross-system question (one θ for all vs per-FF
  detection) is exactly what the control-law tier addresses.

Infra fixes that made the sweep reliable: `spawn` multiprocessing (forked
workers + mp.Queue feeder threads deadlocked on respawn), `CPP_RECOMPILE=0`
(else every worker import runs `make clean` and wipes sibling .so files —
the "disappearing libMMFF_lib.so" bug), eFF worker `chdir(tests/tEFF)`,
GH `m` bound capped at 8 (DynamicOpt mMax=8; m>8 killed workers).

### 9.2 Tier 2 — Symbolic/scalar control law (the sweet spot)

Per the chat doc's design: keep an explicit update rule, let a controller
adjust a few coefficients from dimensionless observables computed from
trajectory data at ~zero eval cost:

- inputs (already measurable): `cos(F,v)` downhill alignment; `cos(F_t,F_{t−1})`
  force persistence; `C⊥` transverse-force persistence (steering vs
  oscillation — the valley signal!); longitudinal curvature `−s·ΔF/|s|²`
  (free BB1); `|F|/F̄`; prediction ratio ρ; reject/reset rates
- outputs: bounded corrections to (momentum β, step η/dt, transverse gain c,
  damping) — *small corrections to baseline, not arbitrary outputs*
- train on **optimization rollouts** (not MD trajectories — the distributions
  differ; a good optimization step isn't inferable from static forces)

Two discovery modes: (a) fit the scalar mapping (tiny MLP ~90 params, or just
search over symbolic expressions — Lion-style program search over
{dot, norm, normalize, project, EMA, clip}); (b) the chat doc's preference —
find an *explicit formula* like `c_t = c0 + c1·max(0,C⊥) + c2·max(0,−C⊥)` so
the result stays interpretable. Our `resets`/`rejects`/`eta_ad` diagnostics
are already the observables it needs.

### 9.3 Tier 3 — Moving-valley geometric model (different axis of sophistication)

Not learned weights — a *better local model*: estimate the valley's tangent
t, turning normal n, geometric curvature `K ≈ θ_t/ℓ` from consecutive
accepted-step tangents (EMA-smoothed to reject zigzag noise), then propose

    Δx = ℓ·t + ½·K·ℓ²·n + α·v_random

— extrapolate along the *curved* centerline instead of the straight tangent.
The chat doc shows the transverse wall penalty grows as ℓ⁴ — predicting the
bend removes the leading error. This attacks exactly the failure mode in
§8.1 (stale secants on turns) and is complementary to any controller above:
it supplies the *proposal direction*, the controller supplies trust.

### 9.4 Tier 4 — Full learned optimizer (only if evidence demands)

Counterfactual training needs a checkpoint API (branch the trajectory,
replay alternatives) — the scan harness already generates thousands of
forked trajectories as raw material. Training hierarchy per the chat doc:
synthetic valleys (Rosenbrock/rotated anisotropic/curved canyons,
randomized stiffness) → real FF rollouts (eFF/UFF/MMFF) → held-out
*molecules*, not just held-out initial states, to test generalization.
Safeguard: fallback to FIRE when the learned rule misbehaves (the
"continually adapting learned optimizer" precedent).

### 9.5 Recommended order

1. ~~§8.2.1 fallback-momentum + §8.2.2 stale-secant forgetting~~ — **done
   (§8.4): bit8+bit16 implemented; momentum fallback is a real win**
   (−34% on MMFF propandiol, now ≈ tuned FIRE)
2. ~~§8.3 low-rank preconditioned Anderson~~ — **done (§8.4): iM=3+autoMass
   functional; helps without bit8, low-rank still fails (m=8 needed)**
3. ~~FIRE-side per-atom masses / mass metric / rank-1 suppression~~ — **done
   (§8.5): implemented and scanned on all 8 systems; mild gains on MMFF only,
   μ≤0.5 suppression stalls — kept as opt-in knobs for the sweep**
4. ~~§9.1 meta-optimization over `scan_all_ff.py`~~ — **done (§9.1):
   differential-evolution sweeps beat hand-tuned baselines 1.4–5.9× on all
   three representative systems; per-FF mechanism choice confirmed**
5. §9.2 scalar control law on the observables we already log — **next**:
   can a fixed low-dim rule recover most of these per-system gains?
6. §9.3/9.4 only if the gap justifies them

## 10. Artifact inventory — which files belong to this study

`git status` on this repo shows many unrelated untracked files (other studies in
`tests/pyFireball`, `tests/tMMFF`, `doc/py`). **The files that belong to THIS
study — the ones another agent needs — are only these:**

### This directory (`doc/Topics/Learning_Optimization_Heuristics/`, all untracked → commit it)

| file | role |
|---|---|
| `GradientMetaOptimization.md` | **this report** — single source of truth for methods + numbers |
| `Learning_Optimization_Heuristics.chat.md` | design discussion / idea backlog (esp. ~line 1320+: memory-timescale decomposition) |
| `scan_all_ff.py` | benchmark harness — static grids, `ONLY_ADAPT` merge mode, `FIRE_MECH` mechanism rows, `--table` renderer |
| `scan_results/*.npz` | per-system result rows (8 files: `mmff_*`, `uff_*`, `eff_H2O_free`) |
| `valleyopt/` | energy-only optimizer package: `valleyopt.py` (algorithms), `bench.py`, `train.py` (offline param search), `firecore_adapter.py` + `FIRECORE_INTEGRATION.md` (FireCore bridge), `test_valleyopt.py`, `README.md`, `IMPLEMENTATION_BRIEF.md` |

### Modified source files (optimizer machinery)

| file | what was added |
|---|---|
| `cpp/common/math/DynamicOpt.{h,cpp}` | `run_GH` (BB1/Anderson/precond), `updateAutoMasses`, `lsAdapt` bits 8+16 + knobs (`lsMomBeta/Gain`, `lsTurnCos`, `lsEmaCoef`, `lsEtaLo/Hi`, `lsAmax`), `modeSupp` + `msU/msFe/msFold/msXold` buffers, `massMetric`, `move_LeapFrog(dt,farr)` |
| `cpp/libs/Molecular/MMFF_lib.cpp` | `run_GH` + `run_FIRE` exports (MMFF+UFF eval bridges), `setAutoMass` (per-FF block layout), `setOptParam` |
| `cpp/libs/Molecular/eFF_lib.cpp` | `run_GH` autoMass layout fix, `setMassMetric`, `setOptParam` |
| `pyBall/MMFF.py` | `run_GH`, `run_FIRE`, `setAutoMass`, `setOptParam` wrappers |
| `pyBall/eFF.py` | `setOptParam` wrapper |
| `tests/tEFF/fire_study_lib.py` | eFF H2O fixture helpers (`set_atom_params`, `load_e`, `snapshot/restore`) used by the scan |

### Modified force-field files (bug fixes this study surfaced — keep them)

| file | fix |
|---|---|
| `cpp/common/molecular/UFF.h` | 0.5× bond-energy sharing (E/F consistency) |
| `cpp/common/molecular/UFF.cl` | same missing `c1·cosθ` term in `evalAngles_UFF` — GPU had the bug too |
| `cpp/common/molecular/{NBFF,MMFFsp3_loc,MolWorld_sp3,libMMFF}.h` | supporting fixes; see `doc/topical_audit/MMFF_energy_eval_debug_2026-11-08.md` for the full FD-consistency evidence |
| `cpp/common_resources/xyz/{propandiol,butandiol,nHexadecan_dicarboxylic}.xyz` | benchmark geometries |
| `tests/tMMFF/data/xyz/*.xyz` | MMFF benchmark geometries |

### Related docs

`doc/topical_audit/eFF_optimizer_speedup.md` (eFF FIRE study, §8 run_GH),
`doc/topical_audit/MMFF_energy_eval_debug_2026-11-08.md` (FF bug forensics),
`tests/tEFF/run_valleyopt.py`, `tests/tMMFF/run_valleyopt.py` +
`scan_gradopt.py` (energy-only valleyopt drivers).

### Build note

Mechanisms live in `libMMFF_lib.so` / `libeFF_lib.so` — rebuild via
`cd cpp/Build && make MMFF_lib eFF_lib`. Other build dirs (`Build-opt`,
`Build-dbg`, `Build-asan`) are stale w.r.t. these changes.

Everything else in `git status` untracked output is **not** part of this study.

## 11. Deep review: improve the information used, not the number of knobs

This review read this findings log, the linked design discussion, the current optimizer implementation, and the scan/meta-search drivers. It preserves the discussion's central constraint: the deployed rule sees **optimization dynamics only**, not molecule identity, topology, or a learned potential. The gradient-free bath discussion is useful context, but the immediate proposal below uses forces already available and needs no random probes or Hessian calls.

### 11.1 What the new results establish, and what they do not

The strongest new result is §9.1: joint tuning beats the coarse grids. It shows that the proposed mechanisms interact and that coarse-grid champions were not strong enough baselines. It does **not** yet show which mechanism caused the gains, whether one policy transfers, or whether more mechanisms are needed. The 95-eval eFF configuration, for example, uses adapt=0: its momentum, turn-forgetting, extrapolation and BB-EMA settings are inactive. Its active ingredients are preconditioned Anderson, memory 6, eta≈.00665 and autoMass=1. Do not interpret the inactive numbers in its JSON as discovered physical laws.

Several earlier explanations should remain hypotheses:

- Frequent non-descent resets do not prove that the valley curved or old secants became stale. Indefinite curvature, incomplete subspace coverage, ill-conditioned coefficients, constraints, coordinate projection, or an incorrect scale can produce the same symptom.
- An accepted-step turn is not proof of changing Hessian: conjugate directions turn on an exact fixed quadratic. Conversely, a nearly straight path can enter a region with a changed Hessian.
- Coordinate-space block masses may help without equalizing the spectrum. Report their measured benefit, not inferred spectral equalization unless independently diagnosed.
- Larger falpha means mixing alpha decays **more slowly**, so velocity is re-aimed toward force for longer. The §4 phrase "too much inertia retained" has the sign backwards. The observed FF-dependent optimum can be real while that explanation is wrong.
- The early claim that gradient-free trajectories beat stalled eFF does not by itself establish barrier-crossing superiority: the chat predates several eFF force/paired-energy fixes. Recheck those examples on the current Hamiltonian before using them as training evidence.

### 11.2 Repair the experimental oracle before spending on more search

**A. Accepted-state success versus a low-force rejected trial.** `run_GH` writes every trial into outE/outF before Armijo acceptance. `meta_opt.py::_cost` and `scan_all_ff.py::milestones` search these traces for the first force-threshold crossing. A rejected trial can therefore count as "converged" even when the returned accepted state is not. The report itself notes a threshold touch followed by wandering. Trial cost must be counted, but trial feasibility/acceptance/status must be separate fields.

Agent task: return explicit accepted-state events and terminal status, plus evaluation indices. Score a valid terminal accepted state (or stop immediately when an accepted state meets the declared criterion); distinguish threshold-touch, accepted convergence, budget exhaustion and failure. On failed search, coordinates are restored to xc but current forces still describe the rejected trial; restore/recompute the accepted energy/forces before exposing a terminal result, and charge any added eval. Require finite active forces before accepting a finite-energy trial. Test deliberately rejected low-force candidates.

**B. Candidate-order dependence in the MMFF/UFF persistent workers.** `_eval_mmff` resets autoMass mode and several knobs, but not invMasses. `set_opt` changes parameters, not masses; `setAutoMass(0)` disables updates without restoring old adaptive masses. `run_FIRE` clears velocity/history but not invMasses; `run_GH` also consumes them. Thus a mode-0 candidate can inherit the last mode-1/2 candidate's D, especially iMethod=3 with autoMass=0. This is a source-level path to contamination, not a claim that every saved winner is wrong. eFF's `initOpt(bMass=False)` explicitly resets its masses and avoids this particular problem.

Agent task: define a complete fresh optimizer-state initializer, preserving explicitly supplied static preconditioners only when requested. Test A, B, A in one worker against A in fresh workers; include autoMass on→off, FIRE↔GH, masks and bounds. Re-evaluate saved winners and baselines after repair. A counterfactual checkpoint must include all of this state, not just positions/zero velocity.

**C. Zero secant columns currently poison the QR helper.** A compiled direct call to `Lingebra::leastSquareFit_QR` with one all-zero column returned rank=1, a NaN coefficient and a nonfinite basis. `nv < rcond*nyj` is false when both norms are zero, so the code divides by zero. Reject zero-norm columns explicitly. Add zero/duplicate/nearly-dependent column tests and compare residuals/rank with a reference solver. The current implementation sorts columns by their initial norm, then performs one modified-Gram-Schmidt pass; it is not full dynamic residual-norm-pivoted CPQR despite its comment. Do not assume a change of rank/memory fixes numerical defects.

**D. Constraints must apply to histories as well as norms.** `set_g` projects fixmask entries to zero, but accepted `Y=-force-gc` and the BB-adaptation y use raw forces. A masked reaction force can enter the secant fit and create a nonzero masked correction. Construct y from projected new/old gradients and explicitly enforce zero direction on constrained components. Test nonzero reactions, not just unconstrained fixtures. This does not assert that the current free-nuclei winner used constraints.

These are small correctness/reproducibility tasks. They take precedence over another large differential-evolution run. The earlier fresh-autoMass initialization issue has already been fixed in `updateAutoMasses`; retain its unit-mass reset in the fresh-state contract.

### 11.3 Key algebra: eta controls the unexplained gradient, not the learned subspace

For the current Anderson proposal, define the least-squares residual:

```text
gamma = argmin norm(Y*gamma - g)
r     = g - Y*gamma
d     = -eta*D*g - (S - eta*D*Y)*gamma
      = -S*gamma - eta*D*r
```

This is an exact identity for the coefficients returned by the solve; a tiny independent algebra check verified it. It remains true with changing D at the current step. In particular, when r=0, changing eta cannot change the direction at all. When history is incomplete, eta sets the move **outside the gradient span the history explains**. It is therefore unsurprising that adapting eta from full-step `s dot s / (s dot y)` can fight the secant model: it measures a different combination of directions. BB1 is not generally the exact line-minimizing step along the current negative gradient; its reciprocal Rayleigh curvature comes from the previous displacement.

Highest-value new diagnostics (reuse the current QR and vector loops):

- `residual_fraction = norm(r)/norm(g)` — how much is left to the base rule;
- norm of `S*gamma`, norm of `eta*D*r`, and their angle;
- descent contribution of each term: `-g dot (-S*gamma)` versus `eta*g dot D*r` (either component can be misleading alone);
- effective numerical rank, coefficient/correction amplification, accepted scale, cap reason, resets and rejects.

If residual_fraction is small but the proposal fails, eta tuning is unlikely to fix the cause. If residual_fraction is large, increasing m blindly may not help unless the new history explores relevant directions. This identifies which actuator is appropriate before a learned controller is asked to choose it. Do not infer unseen residual-direction curvature from a full-step secant without checking excitation; coupling makes that inference underdetermined.

### 11.4 First new mechanism to test: bound the secant correction instead of erasing it

Current failure handling is discontinuous: one non-descent correction clears all history and jumps to the base gradient/momentum rule. A cheaper middle ground needs only scalar dot products. Write:

```text
d0 = -eta*D*g
c  = -(S - eta*D*Y)*gamma
d  = d0 + theta*c
q  = eta*(g dot D*g)                 # positive on the free subspace
```

Require `g dot d <= -eps*q`, with 0<eps<1. If `g dot c>0`, use:

```text
theta_max = min(1, (1-eps)*q/(g dot c))
```

Otherwise theta_max=1. Optionally reduce theta further with a correction-length bound in the chosen fixed coordinate metric, then retain the existing actual-step cap/line search. Start with the existing descent margin scale (eps≈1e-3) as a **test value**, not a universal optimum. For theta=0 the base direction descends if D is positive; theta=1 recovers the current proposal. Project constrained/bound-active components first and recheck the final direction. Log theta and the reason it was reduced.

This changes failure handling, not the LS solver or its memory budget. It preserves usable portions of history instead of converting every doubtful direction into complete amnesia. Persistent small theta must trigger diagnosis/forgetting, not an endless accumulation of bad secants. A correction that is already descending can still be excessively large, so the length safeguard is separate. This is an **untested candidate**, not a reported speedup.

Second ablation, only after that baseline: regularize the small LS fit. Damping coefficients merely by `lambda*norm(gamma)^2` depends on arbitrary column scaling; normalize history consistently or regularize the resulting correction/step in a defined metric. Because m≤8, QR with an explicit regularization block is practical and avoids normal-equation inversion. Do not simultaneously add learned rank, momentum and new forgetting rules to this test.

### 11.5 Replace turn-triggered amnesia with evidence of model error

Before appending a new accepted pair, test whether the old model predicts it. If the new s is well represented by the old displacement span, fit `S_old*a≈s`, and compare `Y_old*a` with the observed y. Use a normalized error and report displacement-span coverage. If s is mostly outside that span, a poor prediction is **missing information**, not proof that all old information is stale. This requires only a small auxiliary fit; sample it periodically if that cost matters.

An even cheaper scalar diagnostic is the endpoint trapezoidal work defect:

```text
work_defect = abs((Enew-Eold) - 0.5*(gold+gnew) dot s)
```

Normalize by energy-change/work magnitudes plus a declared numerical floor. It is exactly zero for a consistent quadratic model, irrespective of turning direction; it increases with nonlinearity or force/energy inconsistency. It is an after-step consistency signal, **not** a pre-step trust prediction, and cannot distinguish those two causes without additional tests. eFF's known force/energy approximant mismatch makes that distinction important near convergence.

The report's proposed trust-ratio formula also needs an explicit model: standard actual/predicted reduction is `(Eold-Enew)/(Eold-Epred)` when the predicted reduction is positive and numerically resolved. A nonmonotone Armijo window maximum is not the current model anchor. `-g dot s` is only a linear prediction, so a ratio below 1 can indicate ordinary positive quadratic curvature rather than a broken model. Do not train a policy on a quantity called "reliability" before defining what it predicts.

### 11.6 Fix the geometric meaning of mode suppression before extending it

For `y=H*s`, the transverse change `(I-t*t^T)*y` measures **coupling of the sampled displacement to other directions**. It is not automatically a stiff eigenvector, a wall normal or a soft-valley turning normal. It vanishes when s lies along an eigenvector, even the stiffest one. Current modeSupp is therefore a heuristic learned direction, not a measured stiff-mode projector.

Two concrete issues deserve isolated tests:

1. **Axis sign symmetry.** Suppression depends on `u*u^T`, so u and -u are the same model. Current normalized-vector EMA does not respect this symmetry. At msEma=.5, two opposite observations of the same axis cancel to zero; a tiny vector check reproduced it. Align a new normalized q with the stored axis (flip q when its dot with u is negative) before blending, or use a projector/Oja-style update. Test invariance to sign flips of equivalent observations. This is not evidence that all .9-EMA runs were affected equally.
2. **Positive factors do not guarantee a positive product.** Integration applies diagonal mobility D after suppression P, i.e. effective mobility D*P. Two SPD matrices need not have an SPD product in the Euclidean descent sense. A 2D check with `D=diag(1,.01)`, `P=I-.5*u*u^T`, `u=(1,1)/sqrt(2)` gives a symmetric-part eigenvalue −0.01338: some gradients yield `g dot D*P*g<0`. If combining them, use a consistent symmetric form such as `D^(1/2)*P*D^(1/2)` in the corresponding scaled coordinates. Power/mixing tests must use that same effective operator, not a diagonal-only metric while integration uses another one. Respect masks in the learned direction and normalization.

Recheck the combination on rotated 2D quadratics before interpreting its molecular stalls as a fundamental failure of low-rank preconditioning. Neither this correction nor adding more modes guarantees a gain. Estimating a few scalar bath statistics, as discussed in the chat, is a separate problem: trajectory directions and Armijo failures are selected, not isotropically random. Failed collinear backtracking trials do not identify an unexplored transverse spectrum. Keep synthetic observability checks before any Bayesian bath estimator.

### 11.7 What a minimal transferable controller should actually control

After the above corrections, use only **two initial actions**:

1. bounded increment of log base mobility/eta;
2. strength theta of the retained secant correction (or a switch between two already-tested rules).

Hold rank/memory fixed initially. Extra actions for every existing knob make counterfactual labels expensive and confounded. The controller can use dimensionless dynamic signals: residual_fraction, descent angle, accepted-step scale/cap status, correction amplification, held-out prediction error where covered, reset/reject rates, and work_defect above its numerical floor. No coordinates, topology, element identity, or molecule name.

**Scale invariance is a cheap falsification test.** For `x'=b*x`, `E'=a*E+c`, gradients transform as `g'=(a/b)*g`, and mobility/eta as `b²/a`. Fixed physical metrics/bounds/tolerances must transform too. A controller trained on raw force norm, fixed eta ranges, or unscaled covariance thresholds can memorize unit conventions instead of learning dynamics. For preconditioned Anderson, eta and an overall factor in D are redundant: only eta*D enters the direction. Fix a normalization convention for D and compensate eta so the effective operator is unchanged before learning two separate scales. With FIRE's changing dt, converting velocity to displacement momentum also introduces the factor `dt_new/dt_old` in carried inertia; dt²*invM alone does not characterize the full update when dt changes.

Use common random rotations, translations, coordinate/energy unit rescalings and initial perturbations during validation. MMFF pi-vector normalization makes the interface a constrained-coordinate problem; verify that saved s equals the **actual evaluated displacement**, not merely the requested a*d, when the callback normalizes/projects coordinates. Histories and counterfactual checkpoints must use that same coordinate contract.

### 11.8 Learning objective and benchmark split

The latest chat clarifies that the goal is efficient local descent through nonlinear/corrugated valleys, not a separately rewarded global search. But scoring only evals-to-1e-6 makes the parabolic tail dominate and can miss the desired improvement. Keep two explicit scores:

- accepted-state local convergence at several force tolerances, with failure rate and terminal energy/geometry checks;
- anytime best **feasible, validated** energy progress over a fixed evaluation-budget distribution, normalized by an independently fixed task scale.

A reference is an offline scoring aid, never an input to the deployed controller. If a candidate beats Eref, retain that information (or update/freeze references for all methods in a new benchmark version); clipping every better result to zero discards progress. Check that low force is not an invalid geometry, wrong Hamiltonian, or shallow/high-energy stationary state before treating it as comparable success. No special reward for barrier crossing is needed unless global search becomes an explicit new objective.

For short-horizon counterfactual training, branch from **identical complete states** and use common future randomness. Choose horizons from the parent state, not separately for each candidate. Label improvement per *force evaluation*, including rejection/restoration calls. Use two horizons or random budgets to discourage a greedy policy that wins the next step but destroys its useful memory. Do not penalize resets as an independent training objective initially: resets are diagnostics and can be a productive cost of progress.

Separate train/validation/test by molecule and initial geometry; include synthetic positive/indefinite quadratics and curved/corrugated valleys to make causal claims possible. Use matched tuning budgets for strong fixed baselines. Report median log evaluation ratio, a tail/failure statistic, and wall time rather than an average dominated by one molecule. Tune the shared rule on multiple systems simultaneously; per-system optima are training/oracle comparisons, not universal defaults. Search conditional active knobs, not dead dimensions: adapt=0 disables BB/momentum/turn/warm-start parameters, and modeSupp=0 disables its parameters.

### 11.9 Bounded next-agent plan

1. **Correctness first:** zero-column/rank tests, projected histories, accepted-state scoring/terminal state, and A→B→A reset tests. Revalidate existing JSON champions before claiming more speed. No new broad meta-run yet.
2. **Instrument the existing solver:** expose residual_fraction, correction scale, rank, theta/descent margin, cap reason, accepted/trial flag and status as structured output. Validate quadratic work_defect and gradient/energy scaling identities. Do not parse printf summaries as the long-term API.
3. **One mechanism experiment:** current Anderson versus descent-bounded correction retention, same eta/m/D/line search, on existing H2O plus a small MMFF/UFF set and one exact quadratic. If it helps, try held-out-error forgetting separately. Charge actual force calls and wall time.
4. **Only then fit the two-action scalar rule:** small counterfactual windows on the existing trajectories; use an explicit clipped/log-linear formula first. Compare against a jointly tuned shared constant rule and an offline per-system oracle. An MLP is justified only if those measured action labels are predictable and a simpler rule leaves useful gains.
5. **Defer new geometric models:** repair/test the sign/metric meaning of modeSupp separately. Multi-mode baths, bent-centerline extrapolation and trained exploration add unidentifiable state unless the simpler controller fails for a diagnosed information reason.

No production source was changed, no molecular sweep/training was run, and no speedup is claimed by this section. Checks performed: current-source zero-column QR probe; algebraic Anderson residual identity; 2D SPD-product counterexample; opposite-axis EMA cancellation. The document backup is `/tmp/GradientMetaOptimization_before_deep_review.md`.

## 12. Executing the §11 plan — fixes, instrumentation, and what re-validation actually showed

Scope: all work below follows the review's ordered plan — correctness first, then instrumentation, then the one mechanism experiment (θ-bounded correction), then re-validation. No new broad meta-search was run. All numbers are evals (energy+force calls, including rejected Armijo trials) to |F|<1e-6 from the standard start, measured on the post-fix build in fresh processes unless stated.

### 12.1 Correctness fixes — what was changed, where, how verified

| review item | fix (files) | verification result |
|---|---|---|
| C: zero-column QR NaN | `Lingebra.h::leastSquareFit_QR` — explicit rejection `nv < 1e-300` before normalization; zero/degenerate columns skipped, `kc` (effective numerical rank) reported to caller | `test_lingebra_qr.cpp` (in this dir): zero column, duplicate columns, exact rank-2-of-4 input → no NaNs, finite γ, correct residual. Before the fix, a zero secant column produced `nv=0` → `q=vwork/nv` → NaN in γ → NaN direction → Armijo rejects storm |
| D: projected histories | `run_GH`: S/Y writes and all inner-product sums (gd, gg, ss, dqs, df2…) skip masked DOFs; Y built from `set_g`-projected gradients so masked comps are exactly 0; `modeSupp` q/u sums made mask-aware | H2O has no masked DOFs → inert there by construction; correct-by-construction elsewhere (reaction forces no longer leak into secants) |
| A: accepted-state traces | rejected Armijo trials write `outF=DBL_MAX` sentinel — can never pass `<Fconv`; `outE` untouched on reject; Armijo accept now requires `isfinite(Ft)` in addition to `isfinite(Et)` | a rejected low-force trial can no longer masquerade as convergence in `milestones`/scoring; eFF's transient NaN-force trials are rejected instead of poisoning the accepted state |
| B: mass contamination | `DynamicOpt::saveMassesOnce()/restoreMasses()` (`am_mass0` buffer); `setAutoMass(mode>0)` saves base masses once, `mode=0` restores them; `run_GH`/`run_FIRE` entry calls `restoreMasses()` when `autoMass==0` | **A→B→A test PASSED** (`test_state_reset.py`): A with autoMass ran 456 evals; B(fixed) fresh = 359, B after A = 359 — bit-identical. Also the eFF block layout for `run_GH` (`am_nBlocks=na+ne`, `am_nScalars=ne`) added — previously preconditioned GH ran with uninitialized block map |
| modeSupp sign | `move_FIRE`: normalized q is sign-aligned to current u (`sgn = sign(q·u)`) before the EMA — u and −u span the same projector u·uᵀ so sign is gauge, but mixing signs in an EMA cancels the mode | opposite-direction transverse-ΔF observations now accumulate instead of cancelling; masked comps excluded from u |

Descent-test margin is **unchanged** (`gd >= -0.001*gg*dd`, line ~499 of `DynamicOpt.cpp`) — `lsDescEps` (default 1e-3) only controls θ inside the bit32 branch. The legacy path is trajectory-identical to the pre-fix build for adapt<32 configs.

### 12.2 Instrumentation — what run_GH now reports

The summary line (one per run) is:

```
run_GH(lsMethod=3 m=6 eta=0.006645 adapt=0) <STATUS> evals=N E=... |F|=...
    resets=R bounded=B rejects=J eta_ad=.. m_eff=.. kc=K rfrac=.. theta=..
```

- `resets` — non-descent events that cleared all history (legacy path)
- `bounded` — non-descent events handled by θ-shrink instead (bit32)
- `rejects` — Armijo trials evaluated but not accepted (each costs a force call; this is where evals go when the line search grinds)
- `eta_ad` — final adaptive base mobility (bit1)
- `m_eff` — final effective Anderson memory (bit2)
- `kc` — effective QR rank at last solve (0–m; how much usable history existed)
- `rfrac` — residual fraction |r|/|g| at last solve; **rfrac→1 means the secant model explains nothing** — the direction is pure −ηDg
- `theta` — mean correction-bound factor over bounded events (θ≈1 = correction retained at full strength, θ≈0 = correction nearly killed)

This is printf-level, not the structured API the review wants long-term — flagged as remaining work. `lsDescEps` exposed via `setOptParam` in both libs.

### 12.3 §11.4 experiment — θ-bounded correction vs full reset (bit 32)

Mechanism: when the proposed Anderson direction fails descent, instead of `nhist=0` + `d=−ηg`, keep the LS solution but split `d = −Sγ − ηDr` into correction and base parts, scale the correction by `θ = min(1, (1−ε)qq/gcc)` so that `g·d ≤ −lsDescEps·|g|²` holds, and **retain the history**. A second variant additionally drip-forgets (drops the oldest half of S/Y) on each such event.

Full results (all runs from η=1e-3 default, adapt=1 = BB1-EMA baseline, +8 momentum fallback, +32 θ-bound):

| config | MMFF propandiol | eFF H2O | UFF propandiol |
|---|---:|---:|---:|
| adapt=1 (legacy clear) | 716 | 436 | 1464 |
| adapt=9 (+momentum) | 615→1083† | 443 | 1462 |
| adapt=33 (+θ-bound) | 528 | **285** | 2979 |
| adapt=41 (+mom +θ-bound) | **509** | 389 | 3541 |
| adapt=33+drip-forget | 528 | 285 | 3541 |

† run-to-run variance on MMFF adapt=9 observed earlier in the study (615 vs 1083); the θ-bounded variants are the stable winners either way.

Diagnostics from representative runs:
- MMFF adapt=41: `resets=5, bounded=96` — 122 resets under legacy become 5 resets + 96 bounded corrections; the evaluator pays θ-shrinks once instead of paying a full-clear and several re-learning steps.
- UFF adapt=33: `θ≈0.9` mean — the bound almost never bites, so the bad secant correction is kept at nearly full strength; `rejects` explodes (5.8k in one run). Tightening `lsDescEps` did not rescue it — the problem is retained wrong history, not the margin.

**Verdict: landscape-dependent, exactly as the review predicted.** MMFF/eFF secant models are *occasionally* wrong → shrunken retention wins. UFF's model is *systematically* wrong → only full clears recover. The retain-vs-clear decision is itself the learnable signal (§11.7's "secant-correction strength" control).

### 12.4 Champion re-validation — the review's warning was correct

All six `meta_results/*.json` champions replayed from fresh state on the fixed build:

| champion | θ (key fields) | claimed | clean replay | verdict |
|---|---|---:|---:|---|
| fire_mmff_prop | dt=.193 fa=.959 fdec=.806 am=2 | 267 | **267** | reproducible |
| gh_mmff_prop | im=2 m=8 adapt=11 am=2 η=.081 | 253 | **305** | inflated by contamination |
| fire_uff_prop | dt=.142 fa=.964 fdec=.837 mm=1 am=0 | 334 | **445** | inflated — ran with leftover adaptive state |
| gh_uff_prop | im=3 m≈8 adapt=8+ am=0 | 370 | **540** | inflated — same |
| fire_eff | dt=.111 fa=.855 am=2 mm=1 | 132 | **132** | reproducible |
| gh_eff | im=3 m=6 adapt=0 am=1 η=.00664515 | 95 | **FAIL** (9198 evals, stalls \|F\|=8e-5) | **invalid — §12.5** |

Contamination mechanics: during the DE sweep a worker evaluated ~200 candidates sharing one `opt`. Candidates decoded with `am=0` skipped `setAutoMass` entirely → `opt.autoMass` kept the *previous* candidate's value and `invMasses` retained adapted entries → the champion was measured with mass preconditioning it did not ask for. The two UFF champions (both `am=0`) show the biggest inflation, consistent with that mechanism. The MMFF/eFF champions that used `am>0` explicitly and the `am=0` MMFF FIRE champion are clean.

Honest current standings (evals, fresh state):

| system | best FIRE | best GH-family |
|---|---:|---:|
| MMFF propandiol | **267** (meta, verified) | 305 (verified) / 509 adaptive-default |
| UFF propandiol | **445** (meta, verified) | 540 (verified) / 1266 adaptive |
| eFF H2O | **132** (meta, verified) | 251 adapt=7 (verified) / 285 θ-bound |

### 12.5 The knife-edge discovery — the most important finding of this pass

`gh_eff` is not merely contaminated — its recorded 95 is not reproducible *even with contamination*, because the objective itself is discontinuous. Fresh-process determinism confirmed (two identical replays → identical evals, E, |F|). Then, clean state, im=3 m=6 adapt=0 am=1, amax=2.377:

| η | result |
|---:|---|
| 0.0066 | **CONVERGED 204** |
| 0.00664515 | **CONVERGED 128** |
| 0.006645146757907484 (the saved θ) | **FAIL** — 9198 evals, stalls \|F\|=8.25e-5 |
| 0.0067 | fails (1.2e-6 @3000) |

am=2, im=3, m=6, adapt=0 — η scan:

| η | 0.004 | 0.005 | 0.006 | 0.0066 | 0.007 | 0.008 | 0.01 |
|---|---|---|---|---|---|---|---|
| result | fail 3.8e-5 | fail 4.9e-5 | fail 4.4e-5 | fail | fail 1.0e-5 | **conv 212** | fail 2.4e-5 |

`evals(θ)` is **discontinuous on a razor scale**: a 1e-9 relative change in η flips converge↔stall, and convergence islands (η≈0.00664515, η=0.008) sit inside a sea of stalls. A single Armijo accept/reject flip early in the trajectory redirects everything after it. DE climbed one razor peak.

**Stall signature** (diagnosed via new instrumentation): near \|F\|~1e-4 the secant model collapses — `kc→0–1`, `rfrac→1.0` — after which ~90% of all further evals are Armijo rejects bisecting `a` (e.g. rejects=8824 of 9198 evals; resets=35–166). Pure −ηDg against autoMass-rescaled DOFs cannot finish.

**Reframe**: point-optima of a discontinuous objective are noise. The *robust* champions are adaptive — eFF adapt=7 re-validated at **251** (previously 268), converging across a wide η range; MMFF adapt=41 at **509** from default η=1e-3. Self-tuning buys **robustness**, not just speed — arguably its real purpose, and the strongest argument yet for the learned two-action controller.

### 12.6 Implications and remaining work

1. **Meta-optimization scoring must be robustness-aware.** Concrete proposal: score θ at ±ε perturbations of the continuous params (η, dt, fa; ε≈1–3%) and take median or max evals; report "neighborhood convergence fraction". Otherwise any point-optimizer (DE, CMA-ES, Bayesian) will keep harvesting razor peaks.
2. **Per-system mechanism choice confirmed on the fixed build**: UFF → momentum-FIRE or full-clear Anderson; MMFF → θ-bounded Anderson; eFF → self-tuning (adapt bits). The heterogeneity conclusion stands and is now backed by clean numbers.
3. **The near-convergence reject-storm** is now the clearest remaining inefficiency on eFF: when `rfrac` stays ~1 for several consecutive accepts, the loop wastes thousands of evals. Candidate targeted fix: fall back to a local quadratic/trust-region step (or freeze autoMass) when the secant model demonstrably carries no information — likely worth more than controller tuning.
4. **Still open from §11**: structured diagnostics API (beyond printf), quadratic work-defect identity validation, prediction-error-based forgetting (turn-based bit16 was neutral in these tests), and the two-action controller fit — now with the requirement that training labels be robustness-weighted.

### 12.7 Files touched in this pass (for the inventory, §10)

- `cpp/common/math/Lingebra.h` — zero-column rejection, kc rank reporting, corrected QR ordering comment
- `cpp/common/math/DynamicOpt.h` — `lsDescEps`, `am_mass0`/`am_bMassSaved` fields, `saveMassesOnce/restoreMasses` decls
- `cpp/common/math/DynamicOpt.cpp` — run_GH rewrite (θ-bound, drip-forget, masked sums, sentinel traces, instrumentation), move_FIRE modeSupp sign fix + mask handling, mass save/restore impl
- `cpp/libs/Molecular/MMFF_lib.cpp`, `eFF_lib.cpp` — `saveMassesOnce/restoreMasses` wiring in `setAutoMass`, `lsDescEps` (+`lsEtaLo/Hi`, `lsAmax` for eFF) in `setOptParam`, eFF run_GH am-block layout
- `doc/Topics/Learning_Optimization_Heuristics/test_lingebra_qr.cpp` — QR rank/zero-col test (compile: `g++ -I cpp/common/math test_lingebra_qr.cpp`)
- `doc/Topics/Learning_Optimization_Heuristics/test_state_reset.py` — A→B→A mass-contamination test
- this file — §12

## 13. Research reset: use offline compute to discover a better algorithm

### 13.1 The intended advance is larger than parameter adaptation

The user's correction is right: §11 answered "how do we improve today's implementation safely?" rather than the more ambitious question in the discussion: **can massive offline experience discover cheap optimization dynamics that humans would not find by hand?** The narrow two-action controller proposed there is a control experiment, not a sufficient scope for the main research program.

§12 strengthens the case for changing direction. Clean replays give FIRE 267/445/132 versus the GH-family's 305/540/251 on the three representative tasks. The 95-eval GH point was not a robust champion. Those results are useful negative evidence: adding handcrafted mechanisms and tuning their constants has not yet produced a substantial, general algorithmic improvement. They do not show that learned nonlinear state evolution cannot do better.

The new objective is therefore **a learned update rule with richer, state-dependent memory**, followed if possible by a short interpretable formula. Thousands of actual relaxation rollouts and many more counterfactual continuations are justified offline. Deployment must still require only a few vectors, dot products and scalar calculations per force evaluation. No molecule coordinates/topology, potential fingerprint, task identity or reference solution is exposed to the policy.

Large offline compute is not evidence of a guaranteed gain. Define success before training: for example, seek a robust 1.5–2× gain over strong baselines on held-out difficult problems, or a major reduction of the slow/failure tail without sacrificing easy cases. These are proposed research targets, not predicted results. A small improvement over one undertuned baseline would not establish the intended advance.

### 13.2 My main research bet: learn the vector-memory recurrence

Do not restrict the learned model to outputs such as dt, falpha, or an Anderson mixing weight. Let it learn **how its small collection of directional memories evolves**. A controller that only interpolates existing proposals is largely confined to the information and failure modes those proposals already supply.

One illustrative family has r=4–8 vector registers, a few scalar memories, and a shared small controller. The registers can express fast/slow momentum, force-response memory, turning direction, or other learned filters; they need not have predetermined semantic roles after training. In displacement units, write:

```text
q_t       = B_t * F_t                         # force converted to a displacement-scale input
z_t       = invariant_features(F_t, q_t, R_t, recent scalar history)
(A,b,c,c0,scale_updates) = controller_theta(z_t)
R_(t+1)   = R_t*A + q_t*b^T
d_t       = R_(t+1)*c + c0*q_t
```

R is N×r. B is a scalar mobility initially, or a bounded block-diagonal mobility already learned from dynamics. All mixed vectors must have consistent units. A sparse recurrence in which each register combines at most two old registers plus q keeps vector work O(r*N); an unrestricted dense A costs O(r²*N) and must be budgeted accordingly. This is a design family for training, not a finalized new optimizer.

Useful inputs are normalized inner products/Gram entries, relative norms, work/progress ratios, recent accept/reject outcomes, memory ages, and running scale estimates. The controller can be a small two-layer network (e.g. 16–32 hidden units) evaluated **once per configuration**, not per atom. No arbitrary N-dimensional neural output is needed: scalar coefficients combine equivariant vectors, so rotating the problem rotates the proposed step. Scale conventions and masks/bounds must be explicit and tested.

Why this could yield more than parameter tuning:

- A learned recurrence can retain forward momentum while damping an oscillating component, rather than globally killing velocity.
- Several filters at learned timescales can approximate useful inverse-curvature responses without storing an N×N matrix or doing a QR solve every step.
- A force-change or turning-memory register can steer into a bend before the accumulated momentum becomes badly misaligned.
- The same recurrence can change its behavior by regime, rather than choosing one globally compromised set of coefficients.

These are hypotheses the architecture makes expressible, not claims that any one mechanism is the answer. The gradient-history span still limits what can be inferred: a learner cannot recover an unexcited direction's curvature merely by having more scalar parameters. Include controlled synthetic problems to test what information is actually sufficient.

The smallest policy need not be the first successful teacher. Train a somewhat richer but still bounded recurrence first; then prune registers, remove actions and simplify. Requiring the final explicit formula before collecting evidence may exclude the very mechanism we hope to discover.

### 13.3 An explicit-rule track in parallel: search programs, not just constants

The complementary route is a typed, cost-limited program space. Programs update a few vector/scalar registers using:

```text
axpy, dot, norm, projection, EMA, bounded scalar nonlinearities,
normalization, force/displacement differences, conditional memory refresh
```

Search the **recurrence structure, memory wiring and control expressions**, with constants optimized inside each candidate. That is different from running differential evolution over the current 10 FIRE/GH knobs. Mutations can replace an update equation, change which register it reads, introduce/remove a projection, or change a memory timescale rule. Respect vector/scalar types and physical units; exclude component-wise sign operations if rotation covariance is required.

Seed the population with actual FIRE and simple inertial/gradient-history rules, but allow rules outside their convex mixtures. Keep baselines as competitors, not as the only allowable proposals. Penalize memory traffic, register count and scalar/vector operations, and retain a Pareto frontier of performance versus deployment cost. Evaluate promoted programs on complete held-out rollouts, not only a short proxy where an aggressive rule can look good before failing.

Two valid discovery paths should be compared:

1. learn a richer controller/recurrence, then distill its behavior into a symbolic rule;
2. directly evolve short symbolic recurrences with offline constant fitting.

Distillation is a hypothesis, not a promised free simplification. If the compact neural recurrence is materially better and meets the runtime budget, it is itself a valid result. Conversely, a searched formula that matches its performance is preferable for inspectability. The relevant methodological precedent is [Symbolic Discovery of Optimization Algorithms](https://arxiv.org/abs/2302.06675): program search and simplification, not copying Lion's coordinate-wise sign update. [VeLO](https://arxiv.org/abs/2211.09760) demonstrates large-scale meta-training of learned update rules, but neither paper establishes transfer or a speedup for these force fields.

### 13.4 First decisive experiment: measure learnable headroom with an offline teacher

Before committing to a large training run, ask a more informative question than "does another fixed formula win?": **from the same dynamical state, are there substantially better continuations, and can our allowed observations predict them?**

1. Collect states along strong FIRE, GH, exploratory recurrent policies and deliberately varied starts. Include productive descent, sharp turns, oscillation, stagnation and final approach. Do not fill the dataset mostly with nearly converged copies of one molecule.
2. Checkpoint the complete optimizer and force-field state. Branch from the exact same parent; include memory registers, masses, histories, constraints, scales, acceptance windows and randomness. Verify uninterrupted versus checkpoint/resume equivalence.
3. Offline, try structured mutations of several controls **together**: memory retention times, longitudinal/transverse response, force-response correction, global scale, momentum quenching and memory refresh. Include a small number of new recurrence structures, not just FIRE-versus-Anderson selection. Optimize short control sequences or piecewise-constant decisions rather than only the next instantaneous move.
4. Run candidate continuations at short and longer horizons, with a common continuation policy and common future randomness where appropriate. Use roughly 16–32 proposals at a short horizon, then promote a few to longer horizons. This can use random search/CEM or another black-box planner; differentiating through the real force field is unnecessary.
5. Record the costs/advantages of **all** actions, not only a hard winner label. Ties and multimodal good actions should remain visible. A choice that wins one step and loses after its momentum/memory evolves should not be called optimal.

A feasible pilot scale is about 1,000 sampled parent states: 16 branches×16 evaluations, with four promoted branches extended to 64, is approximately 448,000 branch evaluations before overhead. This is an illustrative offline budget, not a job launched here. Later stages can use thousands of base trajectories and millions of branch evaluations if the pilot shows useful signal.

The teacher uses privileged branch outcomes **offline**; deployed inference does not run these branches. Its headroom is conditional on the action/proposal family and horizon, not a theoretical global optimum. Measure separately:

- teacher advantage over the strong baseline;
- fraction of that advantage predicted from the permitted dynamical features;
- fraction recovered by the deployed recurrent rule on complete unseen rollouts.

If teacher advantage is tiny, enlarge/change the proposal or memory class before training a bigger predictor. If advantage is large but not predictable, test missing dynamical memory or a longer observation window, not molecule identity. If predictions are accurate but rollouts fail, the problem is distribution shift or short-horizon credit assignment. These outcomes would tell us **why** a learned optimizer fails, instead of producing another inconclusive parameter table.

### 13.5 Training: use all the offline statistics to change the algorithm

**Stage A — supervised warm start.** Train the controller to predict action values/relative advantages, or regress a distribution over good control settings, from the invariant state features. Train vector-memory updates as part of the recurrent policy, rather than only fitting a static mapping from existing FIRE statistics. Handcrafted feature-controller training is a useful smaller baseline.

Store trajectory prefixes/sequences, not just isolated feature snapshots, when training the recurrence: its memory must be generated by its own preceding updates. Counterfactual actions within one fixed architecture share identical parent memory. Comparing different recurrence structures requires full rollouts or an explicit, validated state-initialization/replay protocol; copying a FIRE velocity into an unrelated learned register is not a fair common-state experiment.

**Stage B — on-policy data aggregation.** Run the learned optimizer, collect its new states, query the expensive teacher on selected uncertain/poor states, and retrain. Baseline trajectories alone do not cover the states a new optimizer creates. Keep failure states and recovery trajectories rather than removing them from the dataset.

**Stage C — full-rollout meta-optimization.** Fine-tune using complete and randomized-horizon optimization costs, including delays introduced by memory decisions. Black-box population/antithetic perturbation methods can optimize a compact controller on nondifferentiable force-field drivers. Differentiable synthetic tasks can also support unrolled training, but their advantage must transfer to real rollouts; force-field derivatives/Hessian information used offline must not silently become a deployment requirement.

**Stage D — simplify and search around successful rules.** Distill the trained actions/recurrence, prune registers, and mutate symbolic programs around the successful structure. Re-evaluate the simplified program on complete rollouts. Keep the richer teacher/policy only when its improvement justifies its measured deployment cost.

This sequence is intentionally more ambitious than §11's "fit two coefficients first." The latter remains a baseline for checking whether the richer learner adds value. Do not let it become an indefinite prerequisite that consumes the entire project in manual ablations.

### 13.6 Dataset and loss must match the intended difficult regime

Thousands of trajectories from one H2O start are not thousands of independent optimization tasks. Split **before generating variants**, by molecule/function family and geometry source. Train across independently varied stiffness, coupling, valley curvature, corrugation, dimension and energy/coordinate scales. Include MMFF/UFF/eFF tasks, with many distinct molecules/starts where feasible, and synthetic curved/nonconvex problems where the relevant mechanisms can be controlled. Hold out whole families, not just rotations of training instances.

Use multi-budget anytime progress plus accepted-state convergence/reliability. Keep the force-accurate local tail as a benchmark, but do not let reaching 1e-6 on smooth small molecules define the entire objective when the research target is nonlinear valley traversal. Offline references/scales are scoring aids only; the deployed rule must never see them. Include all actual force calls, rejected trials and safety/recovery work. Best feasible energy progress can reward useful traversal of small corrugations without turning the task into a separately optimized global search.

The safety/acceptance envelope should be common across competitors and explicitly logged. Keep hard physical/domain constraints and finite accepted states; do not impose an extra "every direction must be steepest-descent-like" restriction that prevents a learned inertial method from retaining momentum in the target regime. Bounded nonmonotone steps are already part of the existing design. A recovery branch is part of the evaluated algorithm and its cost, not hidden free assistance.

Evaluate neighborhoods of learned parameters, numerical precision changes, rotations, unit rescalings and perturbed starts. These tests target §12's knife-edge behavior. Do not train every discrete acceptance flip away or claim differentiability where none exists; optimize expected rollout performance and the failure tail over that distribution.

### 13.7 What would constitute a real result

Compare against:

- jointly tuned **shared** FIRE/GH settings using the same training-task access;
- robust per-task tuned baselines as an offline oracle comparison;
- a learned selector/controller over existing algorithms, to isolate new recurrence learning from mere algorithm selection;
- the learned recurrence with fixed coefficients, to isolate state-dependent control;
- memory/register ablations, to establish what information drives the gain.

Report force calls and wall time separately; a tiny network or r-vector update can still be expensive relative to an extremely cheap classical force kernel. Set a deployment budget in measured overhead/memory traffic, not just parameter count. Start with a practical target such as ≤10% optimizer overhead on representative systems, then tune the architecture/cadence or distill if it misses that target. This is a proposed engineering budget, not an established property of the architecture.

A substantial win must survive held-out tasks and matched tuning effort, not only beat an average default. If improvements occur only on the difficult tail, report that explicitly; it can still be valuable. If the learner merely reproduces FIRE or selects Anderson when appropriate, that is informative but not the intended new-algorithm result. If no headroom appears after broadening the proposal class, that is also a valid negative outcome.

### 13.8 Concrete project handoff

1. **Freeze trustworthy evaluation/checkpoint contracts** from §12. Treat this as finite infrastructure work, not the research program itself. Save baseline traces and accepted-state metrics so future comparisons remain comparable.
2. **Build a reusable trajectory/branch dataset generator**, sharing the existing force-field adapters and optimizer state rather than copying evaluation code into scripts. Pilot a bounded set of diverse states and short/long continuations; estimate teacher headroom and recoverability.
3. **Prototype one richer recurrent vector-memory family and one symbolic program grammar.** Document which novel dynamics they can express, their units/symmetries, register count and measured operation budget. Allow more than two learned controls where the coupled recurrence needs them; regularize complexity rather than forbidding it in advance.
4. **Train the pilot policy and run it end to end on held-out tasks.** Use on-policy data aggregation when it visits unfamiliar states. Only scale to thousands of trajectories/millions of branches after the pilot establishes an advantage that a cheap policy can recover.
5. **Distill and perform decisive ablations.** Separate gains due to a novel memory update, state-dependent control, better scale calibration and algorithm selection. Prefer the shortest deployable rule that retains the actual benefit.

No new training or benchmark was launched for this research reset, and no speedup is claimed. This section is a proposal for using substantial offline compute to make an algorithmic discovery, rather than another round of local manual tuning. Backup: `/tmp/GradientMetaOptimization_before_algorithm_discovery.md`.

## 14. §13 execution — checkpoint contract + branch-lab headroom pilot (v1 results)

### 14.1 Implemented infrastructure

**Fat checkpoint** (`DynamicOpt::stateSize/saveState/restoreState`): serializes every optimizer field (FIRE scalars, autoMass buffers incl. `am_mass0`, modeSupp buffers, full GH/Anderson workspace, pos/vel/**force**/invMasses contents) into a flat double blob. Restore copies contents into the *currently bound* arrays — ff-owned buffers are never re-bound. Fails loudly on magic/n mismatch. `force` must be in the blob: `move_FIRE` reads it on the first post-restore step; omitting it silently fed stale forces from the previous branch (caught by the equivalence check below, not by eyeballing).

**Step API** (both `MMFF_lib` and `eFF_lib`): `optBind()` (same binding `run_GH` does, incl. `getForce` bridge and UFF `vapos` aliasing), `optInitOpt(dt,damping)` (fresh-start semantics identical to `run()`), `optEval()` (evaluates at current pos through the lib's eval bridge), `optMoveFIRE()`, `optNDOFs()`, `optBuffer(which)` → numpy views of pos/vel/force/invMasses, `optStateSize/Save/Restore`, plus `dt/damping/lastNeg` added to `setOptParam`. Segmented Python stepping is trajectory-identical to the internal `run()` loop because all FIRE state lives in `opt` fields.

**`branch_lab.py`** (study dir): collects a stepwise FIRE trajectory, snapshots every accepted step, picks ~K parents, applies the action family to restored state, continues each branch H evals, records (features, per-action ΔE). Actions mutate **only optimizer memory/scalars** — never positions — so a branch is a counterfactual decision, not a restart.

Action family v1 (17): baseline; `v{0,0.5,1.5,2}` velocity scaling; `steer{0.25,0.5,1.0}` rotate v toward f̂; `killmode{0.5,1.0}` remove v-component along transverse-ΔF direction; `dt{0.5,2}`; `dampreset`; `fresh` (v=0+dampreset); `im{0.5,2}` global invMass scale; `kick` (v += 0.5|v|ŝ along last accepted step).

Features per parent (all dimensionless, computed from trajectory history only): `cos(F,v)`, `log10|F|`, `|v|/|F|`, BB1 `s·y/s·s`, `|F|/EMA(|F|)`, last-step progress per unit force.

**Checkpoint-equivalence is verified, not assumed**: 'base' branch from parent ix reproduces the recorded `E[ix+H]` to `diff=0.00e+00` on all three systems.

### 14.2 Pilot results (17 actions, H=24 evals, 20 parents per system)

| system | baseline ΔE | oracle ΔE | headroom mean | parents with any better action | best actions (win-rate) |
|---|---:|---:|---:|---:|---|
| MMFF propandiol | 0.003 | 0.004 | 0.001 (max 0.008) | 90% | v1.5/dt2.0 20% |
| UFF propandiol | 0.004 | 0.006 | 0.002 (max 0.012) | 95% | dt0.5 30%, v2/im0.5 15% |
| eFF H2O | 1.276 | 1.285 | 0.009 (max 0.128) | 45% | v2/steer1.0 15% |

Feature-predicted best action recovers 85–99% of the (small) oracle gap.

### 14.3 One-shot results (17 actions, H=24, 20 parents/system)

| system | baseline ΔE | oracle ΔE | headroom mean | improvable parents | best actions |
|---|---:|---:|---:|---:|---|
| MMFF propandiol | 0.003 | 0.004 | 0.001 (max 0.008) | 90% | v1.5/dt2.0 20% |
| UFF propandiol | 0.004 | 0.006 | 0.002 (max 0.012) | 95% | dt0.5 30%, v2/im0.5 15% |
| eFF H2O | 1.276 | 1.285 | 0.009 (max 0.128) | 45% | v2/steer1.0 15% |

Feature-predicted best action recovers 85–99% of the (small) oracle gap.

### 14.4 Persistent-policy branches (`--persist`: rule re-applied every step — the actual recurrence hypothesis)

Same actions, but each is a *policy* applied every step of H. This is where real signal appears:

| action | UFF win% | MMFF win% | eFF win% | mean ΔE (UFF/MMFF/eFF) |
|---|---:|---:|---:|---|
| base | 10 | 15 | 35 | +0.004 / +0.002 / +1.276 |
| **v1.5 persistent** | **70** | **40** | 25 | **+0.007** / +0.002 / +1.223 |
| v2.0 | 10 | 20 | 10 | +0.005 / +0.001 / +1.261 |
| v0.5 | 5 | 0 | 5 | +0.002 / +0.001 / +1.268 |
| kick | 5 | 15 | 5 | +0.005 / +0.001 / +1.248 |
| v0.0 (quench) | 0 | 0 | 0 | **−2.2 / −1.8 / +0.86** |
| steer1.0 | 0 | 0 | 0 | −3.3 / −1.3 / **−30.6** |
| dt2.0 | 0 | 0 | 0 | −3.2 / −1.9 / −9.3 |
| fresh | 0 | 0 | 0 | −0.26 / −0.75 / +0.90 |
| dt0.5 / im0.5 / steer0.5 / killmode / dampreset / steer0.25 / im2 | ~0 | ~0 | ~0 | ≈ baseline |

### 14.5 What the pilot established

1. **Headroom lives in persistent dynamics, not interventions** — *for an undertuned baseline*. With `falpha=0.8/fdec=0.5` defaults, `v←1.5v` each step (extra momentum injection) beat baseline on 70% of UFF, 40% of MMFF, 25% of eFF parents, and persistent quench/alignment/overshoot were catastrophic (−2 to −30 ΔE).
2. **Calibration control result — most of that was parameter tuning, not new dynamics.** Re-running UFF with the verified champion params (dt=0.142, fa=0.964, fdec=0.837): `v1.5` drops to **0% wins**, residual headroom concentrates in global `invMass` scaling (`im{0.5,2}` 45% each) — i.e., still effective step-size. The tuned FIRE dynamics leaves only ~2–8% relative ΔE headroom over 24 evals for the entire single-knob perturbation family.
3. **That is the pilot's real answer**: once the optimizer is well-calibrated, *no single scalar intervention* unlocks much — supporting the §13.2 thesis that better continuations require richer state (multiple timescales, mode-resolved memory), which this proposal family cannot express. A negative-but-informative result: the headroom question is now scoped to genuinely new dynamics.
4. **eFF is the hard case for predictability**: feature→best-action fit recovers −154% of gap (worse than always predicting base) — 6 dynamical features cannot identify *when* momentum boosting helps on eFF; UFF-tuned dropped to 19%. More memory is the identified missing ingredient.
5. **Infrastructure verified**: exact checkpoint equivalence (0.00e+00 on all systems), per-branch eval accounting, dimensionless features, action-value tables in `branch_results/*.npz` — ready as the recurrent-rule training set.

### 14.6 Next increments

1. ~~Calibration control~~ — done above: `v1.5` was an undertuned-baseline artifact; only invMass/step-scale slack remains.
2. **Richer proposal family**: multi-timescale/stateful branches — e.g., two-timescale momentum (fast+slow v), mode-resolved damping (per-block or per-u), steer-then-hold sequences. If scalar families show no headroom, the next test must add *memory dimension*, not more scalars.
3. **Harder state sampling**: parents at resets/stalls/turns (`cos(F,v)<0`, `rfrac→1`), not the converging tail.
4. **GH checkpointing**: `run_GH` accept-loop locals (`nhist`, `Ewin`, `a_last`, `eta_ad`) need a step-function refactor to branch Anderson fairly.
5. **Metric**: evals-to-threshold alongside ΔE@H (near-converged parents make ΔE degenerate).
6. Then the §13.2 recurrent family (r=4–8 registers, small controller) — the branch dataset is its training input.

## 15. Register-memory policies — PyFire replica, branch wins, and the rollout reality check

### 15.1 Python FIRE engine (`PyFire` in `branch_lab.py`)

Reason: the C++ step API executes `move_FIRE` exactly but cannot express new state (registers). A Python replica implements the identical math — `eval_cos_vf` (plain metric), damping update `cf=renorm·damping, cv=1−damping` / `cos<0 → cv=cv_kill, cf=0`, `lastNeg>minLastNeg → dt·=finc, damping·=falpha`, `cos<0 → dt·=fdec, damping=damp_max, lastNeg=0`, force-limit `dt` shrink, `dr_limit` rescale, then an open **policy hook** between the velocity mix and the leapfrog kick (`hook(d,st,ctx)` mutates `v` before `x+=v·dt`).

**Verification (`--verify`)**: `max|Δpos| ≤ 5e-14` over 30 steps on all three systems (UFF 1.8e-15, MMFF 4.9e-14, eFF 5.2e-18) — bit-equivalent to `move_FIRE`. The first verify attempt appeared to diverge; root cause was *test-side*: `reset()` restored positions but not FIRE scalars, so the Python run inherited the cpp run's residual `dt/lastNeg`. Fixed by snapshotting the post-init blob and `restore()`-ing before each engine — the checkpoint contract again earning its keep.

### 15.2 Register policy family (each = a few vector ops + private memory)

| policy | mechanism |
|---|---|
| `2tslow{g}` | slow-EMA velocity register `vs←0.9vs+0.1v`, boost `v+=g·vs` (two-timescale momentum) |
| `mkill{g}` | online transverse-mode unit vector `u` (EMA of sign-aligned `ΔF⊥s`), damp `v−=g(v·u)u` |
| `psteer{b}` | rotate `v` toward `f` by fraction `b`, preserving ‖v‖ |
| `vcurve{k}` | curved-valley extrapolation: extra `½k‖s‖²` displacement along turning normal `n=ŝ−ŝ0` |
| `gmkill` | `mkill` gated on `cos(F,v)<0.7` — only kill transverse modes while misaligned |
| `g2ts` | `2tslow` gated on `cos(F,v)>0.8` — only boost while aligned |

### 15.3 Branch results (H=24, 20 parents, tuned baselines)

| system | top actions (% parent wins) | headroom mean | feature predictability |
|---|---|---:|---:|
| UFF | mkill0.5 45%, 2tslow1.0 35% | ~1% of ΔE | 90% recovered |
| MMFF | 2tslow1.0 60% | ~0 | 49% |
| eFF | 2tslow0.5 45%, mkill0.5 25% | ~1% | 96% |

**Registers do beat baseline at branch granularity** — mode-resolved (`mkill`) and multi-timescale (`2tslow`) memories are the consistently preferred family, steering/curvature extrapolation is not.

### 15.4 Full rollouts — evals to |F|<1e-4 (`--rollout`, 5000 budget)

| policy | eFF | MMFF | UFF |
|---|---:|---:|---:|
| base FIRE | **520** | **1621** | 5000 stall* |
| 2tslow0.5 | 1636 | **1292 (−20%)** | 5000 |
| mkill0.5 | 619 | 2500 | 5000 |
| psteer0.1 | 884 | 2570 | 5000 |
| vcurve1.0 | 541 | 1628 | 5000 |
| gmkill | 565 | 1753 | — |
| g2ts | 563 | **1361 (−16%)** | — |

*UFF reaches only ~5e-4 in 5000 evals under these params — the tuned UFF champion used `massMetric=1` which PyFire's plain metric doesn't replicate; a handicapped but *common* baseline for all policies.

### 15.5 Reading — the §13.5 diagnostic resolves cleanly

- **Branch wins do not imply rollout wins.** `2tslow0.5` won 45% of eFF parents yet is 3× *worse* over a full trajectory (1636 vs 520). Short-horizon ΔE rewards momentum injection; over a trajectory the same boost must then be paid back in damping/resets.
- **Gating is the missing piece, and it works mechanically**: `g2ts` recovers 91% of the unconditional-2tslow damage on eFF (1636→563) and on MMFF converts a −20% mechanism into a −16% *deployable* policy. `cos(F,v)` — a free, dimensionless dynamical statistic — is already an adequate gate.
- **But the honest headline stands**: against calibrated FIRE, exactly one register policy wins (2tslow0.5 on MMFF, −20%) — and MMFF's baseline is itself handicapped (no `autoMass`; the real champion ran 267 to 1e-6). Every mechanism that "wins" so far recreates something tuned FIRE already encodes (slow-EMA momentum ≈ higher effective `falpha`; `mkill` ≈ what `massMetric`/`autoMass` preconditioning achieves).
- **Conclusion for §13**: the proposal space is now right (memory + conditional updates), the verification stack is airtight, and the limiting factor is *policy fitting against complete-rollout objectives* — i.e., learn `gate(features)` coefficients and register gains against evals-to-convergence on-policy, rather than branch ΔE. That is precisely §13.2 steps 2–3: the branch dataset (predictability ≥90% on UFF/eFF) says the actions are *learnable*; the rollouts say only the *timing* of applying them needs learning.

### 15.6 Gate×gain scan on rollouts (`--gatescan`) — the first cross-system win

2-parameter grid over the two gated registers, scored on **complete-rollout evals to |F|<1e-4** (5000 budget):

| system | base | best `g2ts` | Δ | `gmkill` best |
|---|---:|---:|---:|---:|
| eFF H2O (1e-4) | 520 | **378** c0.5 g0.25 | **−27%** | 495 c0.3 g1.0 |
| eFF H2O (1e-6) | 742 | **528** c0.5 g0.25 | **−29%** | 701 c0.3 g1.0 |
| MMFF propandiol | 1621 | **809** c0.9 g0.25 | **−50%** | ~1650 flat |
| MMFF **butandiol** (held-out) | 1566 | **761** c0.5 g0.25 | **−51%** | ~1550 flat |
| MMFF **formic_dimer** (held-out) | 4581 | **2178** c0.7 g0.25 | **−52%** | ~4550 flat |
| UFF propandiol | 5000 stall | **4416** CONV g1.0 | escapes stall | ~5000 flat |

**The rule**: `vs ← 0.9·vs + 0.1·v; if cos(F,v) > c_gate: v += g·vs` — boost velocity along its slow-EMA direction only while moving downhill. 4 vector ops + one EMA register + one scalar gate. No geometry, topology, or system identity.

**Response surface**: plateau at g=0.25 across all gates 0.5–0.9 on every system (e.g. MMFF: 809/911/926 vs base 1621 — not knife-edge). Cliff at g≥1.0 (diverges/stalls everywhere) — the safe region is unambiguous. `gmkill` (transverse-mode kill) is flat-neutral everywhere: once momentum is gated correctly, extra mode suppression adds nothing on these landscapes.

### 15.7 Interpretation and caveats

- **Why it works**: FIRE's momentum mixing `v←cv·v+cf·f̂` rescales velocity toward the *instantaneous* force direction every step; the gated slow-EMA register adds a *persistent* directional component that survives the resets. When `cos(F,v)` is high the optimizer is in a coherent descent phase — exactly when extra inertia along the established direction pays off (long valley). When cos drops (turning/wall) the gate closes and FIRE's own damping handles it. It is a *conditional* momentum timescale — the simplest possible version of §13's "learned memory" that actually helps.
- **vs the undertuned-baseline artifact** (§14.5): this is different in kind — g=0.25 is a *conservative* additive term gated by alignment, not a global scalar multiplication of v. It wins on held-out molecules never used for tuning, at gains far below the divergence cliff.
- **Honest caveats**: (a) all baselines here are handicapped — no `autoMass`, plain `massMetric=0`. The eFF meta-champion (132 evals) used `am=2`+preconditioning; whether g2ts *stacks* with autoMass or is subsumed by it is the decisive next test. (b) gains vs. *handicapped* FIRE, not vs champions — the fair test is g2ts-in-C++ vs the full champion configs. (c) UFF under `massMetric=0` stalls regardless — can't read UFF yet.
- **Status vs §13 targets**: 1.5–2× improvement on held-out problems — **met vs the common baseline on 4/5 systems**, pending the un-handicapped-baseline check.

### 15.8 g2ts ported to C++ `move_FIRE` (exact parity)

Fields `gate_cos`, `gate_g` (0=off), `vs_ema`, register buffer `vSlow` + `vSlowInit` in `DynamicOpt`; params exposed via `setOptParam` in both libs; `vSlow` included in the state checkpoint blob (flag bit 16, NSCAL 69). Gate uses **post-mix** cos(F,v) recomputed after the velocity update — required for parity with the PyFire hook (first version used pre-mix `cos_vf` → 4e-2 divergence; post-mix version matches PyFire to **1.1e-12 over 200 steps**).

### 15.9 The decisive test — champion replay ± gate (`test_gate_champion.py`)

Same eval path as meta-opt (`_eval_mmff`/`_eval_eff`, |F|<1e-6), saved champion θ, gate_cos=0.7:

| system | champion off | gate_g=0.25 | gate_g=0.5 |
|---|---:|---:|---:|
| MMFF propandiol | **267** | 3044 | 546 |
| UFF propandiol | **445** | 521 | 556 |
| eFF H₂O | **132** | 339 | 380 |

**Result: the gate is uniformly harmful under the champion configs** — the mechanism is *subsumed*, not stacked. `autoMass` preconditioning + `massMetric` + tuned FIRE parameters already deliver (better than) whatever the gated slow-momentum register was providing on the handicapped baseline.

### 15.10 The pattern, stated plainly

Three independent mechanisms now show the same lifecycle:

1. `v×1.5` persistent momentum: 70% UFF parent wins → **0%** after baseline tuning (§14.5).
2. Global invMass scaling: residual ~45% wins → effective step-size only.
3. `g2ts` gated register: −27 to −52% vs handicapped baselines **on 5/5 systems including 2 held-out molecules** → +20% to +1140% evals vs champions.

**Interpretation**: the information budget of a weakly-tuned FIRE is large — almost any reasonable memory mechanism improves it. The tuned champion's information advantage is real but *already exploited*: adaptive inverse-mass estimates are themselves a learned per-DOF memory register encoding directional mobility — precisely the role `vSlow` was playing. What remains genuinely unexplored is whether a *switching/policy* layer ON TOP of the full-featured optimizer (autoMass+massMetric+momentum registers) can adapt these mechanisms over the trajectory — i.e., the next controller must operate on the **strong** baseline, not replace it.

### 15.11 Updated next steps

1. ~~Port g2ts + champion replay~~ — done: subsumed. The C++ implementation remains (off by default; useful as a primitive for learned policies).
2. **Re-run headroom/predictability analysis against the STRONG baseline** (champion θ + autoMass): if branches from champion states show ~zero headroom for every register proposal, then the per-step policy question is closed for this action family and effort should go to §13's third leg — *program search over update rules* evaluated on complete rollouts of the strong optimizer.
3. If some headroom remains: it will be conditional (predictability was already ≥90% on 2/3 systems) — a small learned gate, not a fixed rule.
4. **Held-out systems still needed for the champions themselves** — all champion tuning was per-system; the robustness of e.g. am=2 across molecules is unknown and is itself a research question.

### 15.12 Held-out check on the champion itself

Replaying `fire_mmff_prop` champion θ on held-out MMFF molecules (same `_eval_mmff` path, |F|<1e-6):

| molecule | default FIRE | propandiol champion |
|---|---:|---:|
| butandiol | **448** | 471 |
| formic_dimer | **860** | 1317 (+53% slower) |

The champion's 267 is partly **overfit to propandiol's landscape** — on unseen molecules its dt/falpha/autoMass combination loses to plain defaults. This reframes the whole study's target once more: neither per-system constant tuning NOR a fixed mechanism generalizes. The only family that demonstrably adapted across the board so far is the *self-tuning* kind (adaptive GH bits, robust across η — §12.5), suggesting the real prize is an **online-adaptive** rule: a mechanism that senses trajectory statistics and adjusts — precisely the §13.2 recurrent-controller idea, but now with the correct strong baseline and a verified implementation stack to test it on.

## 16. Headroom probe on the STRONG baseline (`--strong` mode, C++ engine)

Setup: champion θ + `setAutoMass(ch.am)` + `massMetric` loaded into the step driver; parent checkpoints along the champion's own trajectory; 18 param-switch actions (gate on at (cos,g) grid, modeSupp/msKill, massMetric toggles, dt×0.67/1.5, falpha/fdec kicks, v×0.5/1.5, amOff) applied persistently for H=16 evals per branch. Checkpoint equivalence under full autoMass state: **exact (0.00e+00)** on all parents — the blob covers `am_k/am_xold/am_fold/am_mass0/invMasses`.

Infra fix during this: `optStateSize` was captured before lazy am-buffer allocation → heap overflow; solved by eager-allocating `am_k/am_xold/am_fold` in `setAutoMass` (both libs) — blob layout now stable from step 0.

### 16.1 Result: per-step headroom is exhausted on tuned trajectories

| system (champion config) | baseline ΔE/16ev | oracle headroom | detail |
|---|---:|---:|---|
| eFF (am=2, mm=1) | 2.117 | **~0.14%** | base wins 27%, all actions ≤13%; `amOff`→1.246 (−41%, autoMass IS the champion) |
| MMFF (am=2) | 0.002 | **~0** | noise-level; predictability 97% of ~nothing |
| UFF (am=0, mm=1) | 0.009 | ~20% rel., tiny abs | g2ts_c0.5g0.5 wins 47% of parents — but §15.9 full-run says it doesn't transfer (556 vs 445) |

**Qualification from §17:** this measures short-horizon endpoint-energy headroom for the tested interventions on sampled champion states. It does not establish per-step optimality for time-to-convergence. Current saved artifacts also contain phase-specific gains and in-sample prediction scores that require the narrower interpretation below.

### 16.2 autoMass is landscape-conditional too

`def` FIRE vs `def+am1/am2/am2+mm1` on held-out MMFF molecules (|F|<1e-6):

| molecule | def | +am1 | +am2 | +am2+mm1 |
|---|---:|---:|---:|---:|
| butandiol | **448** | 525 | 493 | 538 |
| formic_dimer | **860** | 842 | 877 | 853 |

The mechanism that defined the propandiol/eFF champions is **neutral-to-negative on unseen molecules** — autoMass is just another per-landscape tunable, not a free upgrade.

### 16.3 Synthesis — what the whole branch-lab campaign established

1. **Per-step headroom on a tuned baseline: ~0** (three systems, two action families, exact checkpoints). Per-step policy learning on these actions is closed.
2. **Every mechanism is landscape-conditional**: g2ts (weak FIRE, −50%) ⊂ autoMass (champion) ⊂ landscape-dependent too. massMetric helps UFF, modeSupp marginal, θ-bounded Anderson helps MMFF/eFF but not UFF (§12.3).
3. **Feature→action predictability: 90–98%** whenever any headroom existed — the *conditional* structure is real; what's missing is the thing to condition ON that has any headroom left per-step.
4. **The remaining exploitable axis is cross-landscape selection/adaptation**: champions don't transfer between molecules (+53% evals), mechanisms don't transfer (+10–20%), but features DO predict which mechanism wins. A deployable algorithm shaped by this evidence = **a landscape classifier on early trajectory statistics that selects the right tuned mechanism/config** — a coarse-grained version of the §13.2 controller (select among proven configs rather than synthesize per-step updates).
5. **Robustness of tuning itself**: champions beat defaults *in-place* (267 vs ~600 MMFF, 132 vs ~300 eFF) but lose out-of-place — meaning "the right config" is a property of (system, mechanism) joint structure that a selector could capture where any single default cannot.

### 16.4 Where effort should go next (updated)

- **(a) Mechanism-selector experiment**: meta-opt per molecule × mechanism-family on a molecule suite → fit trajectory-feature → config map → evaluate held-out. Uses existing infra (`meta_opt._eval_*`, features from `branch_lab`), answers "is the win in selection?".
- **(b) Program search (§13.3)** over update rules evaluated on *complete rollouts of the strong baseline* — only remaining place per-step novelty could exist.
- **(c) If selector wins**: distill it into an online rule = the adaptive algorithm the original brief asked for.

### 16.5 First selector data — config×molecule matrix and an |F|₀ split

evals to |F|<1e-6, MMFF, configs = {default, propandiol-champion θ, def+am2, def+mm1}:

| molecule | def | prop-champ | def+am2 | def+mm1 | 30-step warmup stats (def run) |
|---|---:|---:|---:|---:|---|
| propandiol | 593 | **267** | 572 | 593 | resets 4/30, mean\|cos\| 0.72, \|F\|₀=1.56 |
| benzene | 101 | **62** | 124 | 101 | resets 4/30, mean\|cos\| 0.90, \|F\|₀=2.59 |
| pyridine | 134 | **84** | 145 | 134 | resets 4/30, mean\|cos\| 0.82, \|F\|₀=4.31 |
| butandiol | **448** | 471 | 493 | 448 | resets 4/30, mean\|cos\| 0.50, \|F\|₀=12.27 |
| formic_dimer | **860** | 1317 | 877 | 860 | resets 3/30, mean\|cos\| 0.84, \|F\|₀=10.95 |

**Observed rule (5/5, small n)**: \|F\|₀>5 → defaults win; \|F\|₀<5 → champion wins. Plausible mechanism: high initial force = far-from-basin/rugged start → the champion's aggressive `dt`/`finc` and fast-adapting masses overshoot into walls and spend evals recovering; smooth-basin starts let the same settings glide. `def+am2`/`def+mm1` never win — the champion is a *joint* config, not decomposable.

Selector value if it holds on a bigger suite: oracle-per-molecule mean ≈ 344 vs always-default 427 vs always-champion 444 — a ~20% mean gain plus tail-risk removal (formic 1317→860). Caveats: n=5, one FF family, thresholds un-fitted — needs ≥15 molecules + held-out split before claiming.

**Historical hypothesis, qualified by §17:** a regime-selecting optimizer (probe → choose config → run) is a plausible and inexpensive baseline. The present data do not establish that per-step learning is exhausted, that unseen-landscape selection is predictable, or that selection is the only remaining avenue for algorithm discovery.

## 17. What local headroom means, why rollout performance differs, and what to learn next

### 17.1 The apparent contradiction is mostly a mismatch of definitions

The new results do **not** demonstrate "the steps are optimal, but the molecule is not." They demonstrate a much narrower property:

> From a small collection of complete states visited by one tuned FIRE trajectory, the tested interventions usually do not improve the **endpoint energy after 16 force evaluations** by much in the aggregate.

That is conditional on the action family, duration, scoring rule and sampled states. Moreover, a persistent action is a 16-evaluation control policy, not a single step. "Near-best tested short continuation" is the right description; "locally optimal step" suggests a stronger result that was not measured.

There are three distinct optimization questions:

| question | what is optimized | what the current experiment answers |
|---|---|---|
| short continuation | endpoint energy after H evaluations, starting at one fixed optimizer state | yes, for sampled states and 18 tested interventions |
| full relaxation | total evaluations to an accepted force-converged state, with acceptable energy/geometry | only the separate full-rollout comparisons |
| transferable algorithm | expected full-relaxation performance over unseen geometries/molecules | not established by one-start champions or the five-molecule split |

A champion jointly tuned on propandiol can be hard to improve by short interventions on **its own propandiol states**, yet transfer poorly to another initial geometry. Those are different states and different tasks. There is no logical conflict.

### 17.2 Why better short-term energy need not mean earlier convergence

Even energy and force norm do not rank states identically on a fixed convex quadratic. For `E(x,y)=0.5*(100*x²+y²)`:

| endpoint | energy | force norm | converged at force tolerance .5? |
|---|---:|---:|---|
| A=(.02,0) | .020 | 2.0 | no |
| B=(0,.3) | .045 | .3 | yes |

This algebraic example was checked directly. A has the lower energy, but B already meets the declared force tolerance. On a nonlinear valley, a branch also changes position, velocity, masses, secant information and future resets. A larger immediate drop can leave poorly aligned momentum or enter a difficult region; a smaller drop can build memory or arrive on a much faster continuation. Neither endpoint energy nor force alone captures all remaining cost.

Let z denote the **complete optimizer state**, not just coordinates. The current score is:

```text
G_H(z,a) = E(z) - E(z_H after intervention a)
```

A score aligned with a local time-to-convergence objective instead asks:

```text
Q_H^pi(z,a) = evaluations already spent in the branch
              + remaining evaluations under continuation policy pi from z_H
```

If convergence happens during the branch, stop and count its actual cost. Keeping the same continuation pi across actions makes the comparison meaningful, but it is still not the true optimal remaining cost. The energy-drop oracle is an oracle for G_H over the finite action list, **not** for Q or all possible algorithms.

If one had the exact optimal remaining-cost function, the full action space and every relevant state, a best next decision would account for its future and the local/global distinction would largely disappear. That is not what the current H=16 probe does. A finite-energy-horizon negative result therefore cannot rule out delayed-benefit actions, action sequences or a new recurrence that visits different states.

The UFF result already supplies an empirical example of the distinction: g2ts wins short branches at some parent states, but applying it through the full champion rollout loses (556 versus 445 in the reported comparison). This demonstrates a horizon/policy-transfer problem, not a violation of optimality.

### 17.3 Reanalysis of the currently saved strong-branch tables

Read-only NumPy analysis was run on `branch_results/mmff_strong.npz`, `branch_results/uff_strong.npz`, and `tests/tEFF/branch_results/eff_strong.npz`. Each contains **15 parents and 18 actions**. The files save action values, features, indices and parent energies/forces, but not horizon/command/config/build metadata; verify their provenance before equating them with a particular printed run. The discussion labels these strong probes H=16.

| current saved file | mean baseline energy drop | mean oracle extra drop | extra / mean baseline drop |
|---|---:|---:|---:|
| MMFF strong | .00217548 | .000191458 | **8.80%** |
| UFF strong | .00874675 | .00170066 | **19.44%** |
| eFF strong | 2.11746 | .00300768 | **0.142%** |

The MMFF aggregate is small in absolute energy, but "zero" is not its relative value. Absolute energy changes from different force fields are not a common performance scale. None of these percentages is a convergence-time speedup.

More importantly, aggregate energy-drop weighting hides trajectory phases. Specific current-file examples:

- eFF parent at evaluation index 50: baseline drop .0752129, g2ts adds .0200631 — **26.7% more short-horizon drop**, parent force norm .5073.
- eFF index 34: baseline drop .0138949, dt0.67 adds .0200563 — **144% more**, parent force norm .7907. Its actual intervention semantics need checking as described below.
- MMFF index 98: baseline drop .0000840338, g2ts adds .0000437979 — **52.1% more**, parent force norm .004733.

The eFF mean is dominated by large-drop portions of the trajectory. A small ratio of mean gains to mean drops cannot establish that all intermediate decisions have negligible headroom. Conversely, a large percentage near an almost-zero drop is not evidence of useful speedup. Report absolute gain, a numerically resolved relative gain, phase/remaining-cost and terminal convergence together.

### 17.4 "90–98% predictable" is currently an in-sample fit, not learned transfer

`branch_lab.py` fits ridge coefficients W using all parent rows and then chooses actions from `X@W` on **those same rows**. The printed quantity is recovered training-set oracle reward, not classification accuracy and not held-out predictability. Adjacent states from one trajectory are also strongly correlated.

The current saved files give:

| file | training-set recovered oracle gain | leave-one-parent-out recovered gain |
|---|---:|---:|
| MMFF | 96.7% | **−180.8%** |
| UFF | 97.6% | 88.8% |
| eFF | **−113.9%** | **−173.4%** |

These were computed with the script's same raw features and ridge coefficient 1e-4. Negative recovery means the predictor selects actions that are worse than always taking baseline on average. Recovery can be extremely unstable when the oracle gap is tiny. Leave-one-parent-out is merely a diagnostic and remains optimistic about transfer because other rows are from the same trajectory; it is not an unseen-molecule test. The eFF file does not support a blanket 90–98% claim even in sample.

Thus there are at least two different failure modes already visible: small short-horizon reward gaps, and poor prediction of useful interventions. They must not be merged into "the actions are learnable but only landscape-level selection remains." The UFF data do contain a promising within-trajectory signal; whether it predicts longer-term value is still open.

### 17.5 Action semantics limit the negative result further

The current strong test uses `make_param_actions`, not unrestricted new register programs or all combinations of controls. Several actions also differ from their informal names:

- With `--persist`, fn is called once before the loop and again before the first move. A multiplicative velocity/dt action is therefore applied twice before that first move.
- Reapplying dt×.67 or ×1.5 every evaluation compounds the intervention; it is not simply selecting a new fixed timestep ceiling. Both are legitimate experimental policies if named/scored accordingly, but answer different questions.
- Persistent modeSupp/msKill calls reset `msInit=0` every evaluation. `move_FIRE` consequently reinitializes its mode estimator instead of accumulating a direction. This does not fairly test a persistent learned-mode rule; initialization belongs at branch entry, separate from the per-step control.
- Newly activated memory registers may start cold at the parent. That tests cold-start activation, not a policy whose register has been accumulating throughout its preceding trajectory. Use explicit prefix replay or a declared cold-start protocol for each new recurrence.

Checkpoint equivalence is a substantial achievement, but it verifies baseline state restoration, not that each intervention implements the intended policy. Before interpreting a no-headroom result, trace actual dt/masses/mode state for a few steps of every action. Save H, flags, parent trajectory identity, action definitions, configuration and build/version with the result file.

### 17.6 What can legitimately be inferred from the molecule/config matrix?

It demonstrates **configuration transfer is imperfect** and supplies a useful selector benchmark. The different winners are real observations for these particular starts. It does not establish that the systems fall into stable classes called "smooth" and "rugged," or that all fixed mechanisms are universally useless outside their winning molecule.

Large initial force can mean initial strain, a stiff but entirely convex mode, a different energy scale or more DOFs. It does not measure nonconvexity, barrier structure or valley curvature. The threshold at 5 is a post-hoc separation of five observations; there is no untouched test set. The state can also change regime during a single relaxation, so "per-molecule versus per-step" is a false binary. A selector may actually be detecting an **initial-state regime**, not an intrinsic property of the molecule.

Expand with multiple initial distortions of the **same** molecule as well as held-out molecules. If the winner flips with distortion, a state-dependent controller is more appropriate than a fixed molecule-level label. Normalize/test features against system size and consistent unit rescaling; the raw-force threshold may still be useful within a fixed FF/unit setup, but is not a universal physical rule.

Arithmetic check on the five listed rows: always-default mean **427.2**, always-champion **440.2**, hindsight per-molecule best **344.2** evaluations. The ~19.4% saving versus always-default is a **post-hoc selection opportunity**, not a demonstrated deployable gain. It also cannot exceed this finite portfolio's oracle on the same task set; a novel update rule is a different search space.

`|F|_0` itself requires the initial force evaluation, not 30 warmup steps. If the eventual selector uses reset/alignment statistics from a warmup, the cost and state transition matter. Thirty probe evaluations followed by a restart and a 62-eval champion run costs 92, not 62. Continuing from the warmed state has a different trajectory and must be measured. A baseline that would finish in 75 evaluations makes an unnecessary 30-evaluation probe a substantial expense.

### 17.7 What to do next for Learning to Optimize

**First, test whether the observed mismatch is a horizon problem.** Use a small representative subset of parent states and correctly defined actions. Record rewards at H=16, H=64 and completion under the same continuation. Compare energy improvement, best feasible energy curve, raw-force progress and actual remaining evaluations. Compute whether action rankings agree across horizons. If short-horizon winners lose later, train longer-horizon value/advantage, not immediate energy-drop labels. If the same actions have no useful long-horizon advantage, that closes **those actions**, not all learned optimizers.

**Second, test policy sequences rather than only one held action.** A few two-segment decisions can reveal delayed benefits: conservative memory-building then acceleration; aggressive traversal then cleanup; activate a mode estimate once, allow it to learn, then use it. Always compare with the matched constant-policy baseline. This is a useful offline teacher experiment because it directly tests the credit-assignment hypothesis without yet training a large network.

**Third, collect states visited by the new policies.** Current parents come from the champion's state distribution. A better algorithm may take a different early route and never visit those states. Repeatedly query only the baseline path cannot establish absence of opportunity elsewhere. Train on complete rollouts/on-policy data aggregation, retaining the strong baseline as competitor. Hold out molecules and trajectory starts; do not fit and score the same branch rows.

**Fourth, keep selection as a baseline, not as a premature endpoint.** Test default versus champion portfolio selection on an untouched suite, including probe cost and multiple starts. In parallel retain §13's broader recurrence/program search, scored on complete trajectories. Current branch families have not tested that space. Only if a richer, properly scored teacher repeatedly finds no useful advantage should the algorithm-discovery objective be reconsidered.

For training labels, use remaining evaluation cost under a declared continuation or two-horizon normalized progress with a terminal accepted-state value. Short energy gains are still useful auxiliary signals, but not the sole target. Partial observability also matters: an instantaneous force norm may not distinguish momentum preparation, a turn, and an approaching instability. Give a recurrent policy the dynamical history needed to predict longer-term value, and measure its held-out regret relative to the finite-action oracle—not just action-label accuracy.

### 17.8 Revised conclusion

The evidence supports: **the tested champion is robust against many sampled short interventions on its own trajectory; some short interventions help locally but do not help when applied throughout the relaxation; configuration choice matters across starts; a selector is worth testing.**

It does not support: **per-step optimality, exhaustion of all learnable dynamics, proven unseen-landscape predictability, or a proof that program search can only rediscover FIRE.** The strongest next scientific question is whether the missed opportunity is longer-horizon control, missing memory/proposals, or genuinely little remaining headroom. The experiments above separate those explanations.

No molecular relaxation or training was run for this review. Read-only analyses used the three existing strong-branch NPZ files; the quadratic example was evaluated algebraically. Backup: `/tmp/GradientMetaOptimization_before_headroom_interpretation.md`.

## 18. §17 experiments executed — corrected baseline, Q-scan, sequences, on-policy states

### 18.0 Blocking bug found first: step-API vs `run_FIRE` parity (367 vs 267)

Before interpreting any sequence result, the step-API baseline had to reproduce `run_FIRE`. It did not: MMFF-propandiol champion converged at **367 evals** through the step API vs **267** via `run_FIRE`, with trajectories identical early and diverging late. Root cause: `lib.set_opt` (called by every `_eval_mmff` replay) sets `minLastNeg=5`, while `DynamicOpt`'s field default is **3** — the step API never set it, so it ran a *different* FIRE variant. (`eFF::run` never sets `minLastNeg` → eFF champions genuinely ran 3; only MMFF/UFF evals used 5.) Fix: `minLastNeg` added to `setOptParam` in both libs + `fp['minLastNeg']=5` in `branch_lab`'s MMFF/UFF driver. **Verified: step API now converges at exactly 267 = run_FIRE.**

Consequence: all previous `--strong`/`--qscan` MMFF/UFF numbers were measured on a minLastNeg=3 pseudo-champion and had to be re-run.

### 18.1 Q-scan re-run on the corrected baseline (answers §17.7-first)

Per-action evals-to-convergence from champion-path parents, plus ΔE@16/ΔE@64 for rank comparison:

| system | baseline remaining | best action | verdict vs previous run |
|---|---:|---|---|
| MMFF (minLN=5, am=2) | mean 140 | massMet1 **−63** | was **+18** — the apparent headroom was entirely the minLastNeg artifact |
| UFF (minLN=5, am=0) | mean 233 | v1.5 −3, g2ts_c0.9g0.25 −51 | was −3..−82 — same conclusion |
| eFF (minLN=3, am=2, mm=1) | mean 70 | all actions +47..+479 | unchanged (eFF baseline was already correct) |

Rank-consistency rows confirm §17.2 again: e.g. on MMFF `v0.5` has the best ΔE@64 (+0.0048) yet loses 123 evals; `dt1.5` has ΔE@16 below baseline but loses only 20. On the corrected baseline **no tested action beats the champion's own continuation** on any system — the §16 "near-locally-optimal" statement now survives the narrower Q-objective too, *within this 18-action family on champion-path states*.

### 18.2 Sequence tests (answers §17.7-second)

`massMetric 0→1` switch on MMFF-propandiol champion: baseline 267; switching at steps {25..225} gives 295–540 (all worse), at 250 gives 237 (already past the interesting regime). `g2ts` warmup→baseline at {10..200}: 295–635 (all worse). Baseline→`g2ts` late: 524–2000 (all worse). **Every two-segment schedule tested is worse than the constant champion policy** — the delayed-benefit hypothesis fails for these action pairs. Memory built under one rule does not transfer profitably to the other, likely because each policy's memory (vSlow EMA, metric history) is only optimal under its own continuation.

### 18.3 On-policy states (answers §17.7-third) — the one positive result

Parents collected on the **default-FIRE** trajectory (`set_opt` defaults, ~870-eval continuation), then branched: (a) stay default vs (b) **champify** — switch `finc/fdec/falpha/dt_max/massMetric` + enable `autoMass=2` mid-run:

```
ix= 10: default 860 -> champify 331   saved +529
ix=100: default 770 -> champify 407   saved +363
ix=250: default 620 -> champify 248   saved +372
ix=500: default 370 -> champify 206   saved +164
```

**Config-switching mid-run recovers most of the champion advantage from every visited state** (+103 to +529 evals) — autoMass re-seeds its secant estimate within a few steps and works even on states it didn't generate. This is the first verified instance of a *non-champion-path* state where a learned decision has real value — and it directly supports the §16.5 regime-selector direction: an online rule that detects "this landscape wants adaptive masses" captures most of the champion's benefit without knowing it a priori.

### 18.4 What the corrected evidence now supports

- On champion-path states, the 18-action family has **zero convergence-time headroom** on all three systems (was already true for eFF/UFF; MMFF's apparent +18 was the minLastNeg artifact).
- Two-segment action sequences of the tested pairs do **not** unlock delayed benefits.
- The exploitable signal is at **config-selection granularity**, and it is **recoverable online** (champify from any default-path state saves 100–500 evals).
- The selector hypothesis (|F|₀ split, §16.5) is now the only mechanism with verified headroom — still awaiting the ≥15-molecule + multiple-start suite.
- Open per §17.8: states from *alternative trained policies* (not just default-FIRE) remain untested; program search over update rules is untouched; and nothing here rules out a learned controller finding value that 18 hand-picked actions missed.

### 18.5 Methodological note (keep for all future branch work)

`minLastNeg` shows that **the same field default silently differs between `lib.set_opt` and `DynamicOpt`** — any step-API experiment that doesn't byte-compare its baseline against `run_FIRE`'s trajectory can measure artifacts for weeks without noticing. The blob-diff check (`save → restore → continue → compare`) caught the shift artifact, but scalar-field parity needed the eval-count comparison to surface. Both checks are now standard: trajectory-shift check *and* fresh-process eval-count parity before any branch experiment is trusted.
