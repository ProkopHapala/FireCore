# GradientMetaOptimization — findings log and open directions

**Current research direction:** §13 supersedes §11's conservative two-action recommendation as the main project: discover a new cheap vector-memory update rule through large-scale offline trajectory experiments, learned control and symbolic program search. §12's clean replays show why more small FIRE/Anderson patches are not a sufficient research objective. Correctness fixes remain prerequisites for trustworthy training labels, not the principal research contribution.

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
