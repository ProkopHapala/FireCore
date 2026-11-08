# Research/development brief for coding agent: topology-blind gradient-free optimization

## Scientific goal and non-negotiable constraints

Build an optimizer that **only receives a scalar energy sequence and the optimization trajectory**. It must NOT read molecular geometry/topology, atom types, chemical graph, force-field details, or gradients. FireCore eFF / UFF / RARFF are *black-box objectives* to which the exact same optimizer must apply. Training can be expensive and use arbitrary optimization rollouts, but inference must be O(kN) arithmetic and memory with k ~2–4, independent of trajectory length. A scalar log for benchmarking is fine; full-length-vector history is not. The desired output of offline learning is **a small, interpretable symbolic update or a fixed numerical parameter set**. A tiny MLP can later *control scalar decisions* (step, trust, exploration), but must never predict a full N-dimensional displacement or see the physical system's features.

Our research hypothesis: on correlated optimization trajectories in smooth curved valleys, a small stateful model can infer a useful moving descent subspace from very few energy probes, while an isotropic gradient-free method wastes most samples. This will not work for arbitrary high-dimensional black boxes with unrelated gradients: a new energy difference gives only one scalar gradient projection. Measure this limitation rather than conceal it.

**Target regime (explicit):** we do NOT aim to beat BFGS/FIRE near a parabolic minimum — gradient methods are already near-optimal there. The goal is to approach their performance on **complex non-convex long-curved-valley potentials**. A stochastic energy-only walker can jump over small corrugations and barriers, so it can be *faster* than gradient methods on non-smooth potentials, especially at small DOF counts. Confirmed on eFF `load_a` starts (tests/tEFF/run_valleyopt.py): FIRE stalls on a plateau for 20k iters while stochastic walks reach a much deeper basin in 3k evals.

## Existing executable starting points

- `valleyopt.py`: common `ask/tell`, five zero-order algorithms, scalar trace.
- `bench.py`: rotated curved valley and 2D Rosenbrock.
- `train.py`: simple offline parameter tuning, normalized area under improvement curve.
- `firecore_adapter.py`: a writable-flat-DOF energy wrapper.
- `test_valleyopt.py`: smoke and state-consistency tests.
- `FIRECORE_INTEGRATION.md`: verified paths, explicit adapter example.

None of the five classes faithfully implements a named paper. Start by extending/benchmarking these rather than rewriting the whole codebase.

## Five algorithms to compare and improve

### A. `MomentumWalk`: cheap control baseline

Maintain one unit tangent `t` (EMA of accepted steps). Propose `d = h * normalize(t + a * v)` where `v` is random perpendicular to `t`; if unknown, pure random with optional opposite probe. On success update `t` and grow h, otherwise shrink h and increase exploration. This is *not* enough to follow a curved valley well, but it is the bare-minimum ablation.

### B. `SecantSketch`: gradient from successful and rejected energy samples

Store `g_hat` (one N-vector). For any trial displacement `d`, predict `Delta E_hat = g_hat.d`, observe `Delta E`, then apply:

```
g_hat <- (1-forget) * g_hat + gain * (Delta E - g_hat.d) * d / (d.d + eps)
```

This is normalized LMS / Kaczmarz-style online identification of **gradient projections**. It is not a Hessian method. Keep a trust radius and optionally use opposite probes to reduce nonlinear bias. Candidate is approximately `-h * normalize(g_hat)`. Important warning: if the center moves or a step is large, the secant observation is not a true gradient at the *new* position. Add curvature- or residual-aware forgetting. Compare sketch accuracy to known synthetic gradient **only for analysis**, never as the optimizer's input.

### C. `CurvedWalk`: curvature of the trajectory

Keep recent unit tangent `t` and turn normal `n`. From consecutive accepted steps `t_old,t_new`, estimate `kappa ~ acos(t_old.t_new)/step_length`; propose:

```
d = h*t + 0.5*kappa*h*h*n
```

Treat it as a predictor. It measures curvature of the *accepted path*, which can differ from the centerline of the low-energy valley (zigzag transverse relaxation is the main confounder). Add low-pass filtering of the turn and confidence diagnostics. It currently does not estimate an energy Hessian.

### D. `WallProbeWalk`: exploit -> hit wall -> transverse exploration

If extrapolated point `y = x + d_forward` fails, hold the incumbent `x` fixed. Let `v` be the remembered normal or a random direction orthogonal to the tangent. Evaluate `E(y+qv)` and `E(y-qv)` and use already measured `E(y)`:

```
gv = (E(y+qv) - E(y-qv)) / (2q)
kv = (E(y+qv) + E(y-qv) - 2E(y)) / q**2
z  = clip(-gv/kv, -2q, 2q)     # when kv>0; otherwise choose lower endpoint
candidate = y + z*v
```

Accept a candidate only if it improves the incumbent (but always save best *evaluated* trial). If the correction fails, test another normal or shrink h. Enhance with a persistent exponentially averaged **bath stiffness** `k_bath <- EMA(kv)` from random orthogonal probes, and a better exploration policy. Test whether probing around *failed forward point* is more informative than probing at current x.

**Important:** ± probes are centered at `y` and are *not independent moves from current x*. They must never update the incumbent or shift the finite-difference anchor before the probe pair is finished. The current code tests this distinction.

### E. `RidgeTR`: local k-dimensional quadratic energy model

Maintain U (N×k, k=2 or 3), g (k), H (k×k), model/trust radius h. At the current anchor x, estimate g and H by energy-only finite differences in local coordinates `a=U.T d`:

```
E(x+Ua) ~= E(x) + g.a + 0.5*a.T H a
```

Use ± probes for diagonal entries and one combined probe per off-diagonal term; k=2 requires 5 probes if E(x) already known. Trust-region step can use regularized eigenvalue-clipped Newton direction within ||a||<=h. After accepted step `Ua`, keep the same frame and propagate `g <- g+H*a`, while evaluating actual-vs-predicted improvement `rho = (Enew-Eold)/(model_delta)` and rebuilding if rho is unreliable. A clearly novel *soft-subspace discovery* mechanism is still required when known directions are exhausted. The current code restarts with the previous accepted tangent and one newly sampled orthogonal vector; this is crude and the highest-priority weakness to address.

**Advanced option:** represent a moving valley itself with tangent t, normal n, longitudinal slope/curvature, transverse stiffness, and geometric curvature K:

```
E(s,z) = E0 + gs*s + 0.5*ks*s*s + 0.5*kperp*(z - 0.5*K*s*s)**2
```

**All fitted quantities stay SECOND-ORDER.** Along the soft direction: slope `gs` and quadratic `ks`. Along the normal(s): quadratic `kperp` — much stiffer, but still quadratic. The geometric path curvature `K` describes where the valley *goes*, not an energy term. The quartic `s^4` penalty is implicit in `(z - K*s^2/2)^2` and must NOT be fitted directly (ill-conditioned); the bend enters only through motion of the frame. In practice: estimate `K` from accepted-trajectory turns, then fit `gs, ks, kperp` from local observations. Do not fit `kperp` using only pure tangent samples: those cannot identify transverse stiffness reliably. Compare against `RidgeTR` using the **same number of energy calls**. This model has not yet been coded; it is the leading next research variant.

## Extracted from the papers: actionable parts, not claimed reproductions

**Lion — Chen et al., Symbolic Discovery of Optimization Algorithms, NeurIPS 2023.** Search over *short typed mathematical programs*, score on cheap proxy tasks, then simplify/select for transfer. Its actual discovered gradient-based rule uses mixed gradient/momentum, coordinate `sign` updates and momentum EMA. **We want the meta-method, not `sign` updates**, because sign is not equivariant under rotation of physical coordinates. Use typed primitives: scalar `{+, -, *, /safe, clip, min, max, EMA, sqrt, log1p}`; vector `{+, scalar_mul, dot, norm, normalize, orthogonal_projection, low_rank_project, sample_orth}`. Restrict expression size, operation count, vectors stored; penalize extra objective calls and numerical failures. https://arxiv.org/abs/2302.06675

**OMoRF — Gross & Parks, Optimization by Moving Ridge Functions, 2021.** Trust-region method with a local low-dimensional active subspace, quadratic model within it, and *separate interpolation sets*: one for determining the subspace (n+1 function samples in the described method), one for fitting a low-dimensional quadratic (for subspace dimension d, about (d+1)(d+2)/2 values). The subspace moves with the trust region. Our `RidgeTR` only implements a simplified local kD quadratic and does **not** implement the full n-dimensional subspace identification or interpolation-point geometry. This matters because n+1 energy evaluations are already unaffordable for n=1e6. https://arxiv.org/abs/2007.04893

**ASEBO — Choromanski et al., From Complexity to Simplicity: Adaptive ES-Active Subspaces for Blackbox Optimization, NeurIPS 2019.** Learn a low-dimensional active subspace from black-box gradient-estimator history, sample inside it and in its complement, and adapt the mixture based on how productive those samples are (their method uses more elaborate online bandit/smoothing-estimator machinery). Our natural adaptation is to measure improvement per evaluation / model residual inside and outside U, then adjust the probability of exploration. Do not label a simple mixture of probes as a faithful ASEBO implementation. https://arxiv.org/abs/1903.04268

**LM-CMA-ES — Loshchilov, GECCO 2014.** Approximate covariance-adapted sampling using a small set of historical direction vectors, rather than storing a dense N×N covariance. Its claimed sampling and memory scale as O(mN). Provides a serious limited-memory baseline and suggests an `adaptive covariance of exploration directions` rather than isotropic random exploration. Our prototypes do not implement CMA updates; compare to a maintained baseline implementation. https://arxiv.org/abs/1404.5520

**Powell direction-set and NEWUOA/BOBYQA.** Powell's original direction-set technique replaces a direction using cumulative movement after line minimizations, so a direction retains information from previous successes; NEWUOA/BOBYQA use underdetermined quadratic interpolation and trust-region steps, maintaining poised sample geometry. Use SciPy Powell and PRIMA COBYQA / Py-BOBYQA as actual baselines, not rough pseudocode imitations. https://doi.org/10.1093/comjnl/7.2.155 and https://www.libprima.net/

## Highest priority next changes for agent

1. **Correctness and deterministic tests:** fail-safe non-finite energies, anchor invariance of all probes, unit tests for diagonal/full quadratic, consistent count of initial/failed evaluations, best-evaluated vs incumbent states. Check energy adapter on each physical backend.
2. **Better subspace tracking:** when a new accepted step d is outside U, update U by modified Gram-Schmidt with a bounded change; rotate/transform or reset g,H consistently. Do not simply project into a stale subspace forever. If old samples are reused, store only a bounded cache of k-dimensional coefficients *in a fixed frame* or transform them correctly when recentering/rotating. Memory O(kN+k^2), not O(TN).
3. **Hard-bath statistics (moments, not eigensystem):** estimate low-order moments of the (n-k)-dimensional stiff complement — NOT exact per-direction stiffness, NOT the Hessian eigenbasis. A random unit probe v gives k_v = vᵀ H_bath v, a sample from a distribution with mean tr(H_bath)/(n-k). If the stiff space is dominated by one direction, observed k_i = Kmax·cos²(a_i), so Kmax ≈ max_i k_i/cos²(a_i) is a usable upper-bound statistic even with unknown a_i. Accumulate samples mostly from ORGANIC rejected evaluations (see item 4). Forget stale samples two ways: by iteration count (cheap EMA) AND by distance traveled — a sample's validity decays once the incumbent leaves its trust region. This is deliberately the cheap, realistic, still-useful goal: bath stiffness sets exploration amplitude and helps distinguish "hit a stiff wall" from "hit an undiscovered soft direction".
4. **Exploration rides on organic failures — no dedicated probes during descent:** while we move along a good descent direction, most trial samples are ACCEPTED; when the valley kinks in an unknown direction we get stuck and the necessary random trials automatically start failing around the incumbent — those failed samples ARE the transverse probes, obtained for free. Exploit this implicit dynamics: record rejected samples as subspace/curvature measurements anchored at the failed point, trigger extra probing ONLY after repeated failures, never burn evaluations on systematic probes while descending correctly. Learned decisions, if any: forward step vs mirror vs probe pair vs replace a subspace vector — scored by *expected reduction per energy call*, not acceptance rate.
5. **Curved valley surrogate:** explicit `K` and its turning normal; compare with local quadratic at same budget, vary curvature and k_perp independently.
6. **Meta-optimization:** first tune ~5–15 scalar parameters using cheap offline random search / CMA-ES. Then optional tiny MLP with 5–10 scalar features, <100 weights, one evaluation *per optimization step*. Inputs only dynamical observables: trust ratio, energy improvement, rejection EMA, tangent change, kappa*h, transverse finite-difference curvature, gradient-sketch prediction error. Outputs only strategy parameters and choice of next probe type. Finally try symbolic regression to replace MLP with an equation.
7. **Ablations:** disable momentum, valley curvature, wall probes, low-rank model, trust-ratio adaptation separately. Ablate rank k=1/2/3/4 and probes per stall.

## Benchmark design

Synthetic cases: positive/negative/rotated quadratics (with condition number 1..1e6), 2D Rosenbrock, nD *rotated* curved valleys with one soft tangent and hard perpendicular bath, paths with 90° turns, S-shaped valleys, branching bifurcations, locally flat saddle/plateau, periodic multimodal potentials. Vary n=2, 10, 100, 1000, ideally 1e6 for memory/runtime (do not construct N×N rotation for n>~100). Embed soft subspaces using a small number of random orthogonal vectors instead. Keep Hessian/gradient available only to benchmark evaluators, not optimizer.

Physical cases: from FireCore eFF, UFF, RARFF; use held-out molecules and perturbations; do not train and test only on the same molecule with different starts. Include atom-only vs electron-position vs electron-size degrees and mixed DOFs. Test with/without coordinate scaling and constraints. Be aware rigid-atom rotation needs a retraction adapter.

Log `nfev`, successful evaluations, best observed energy, best incumbent energy, wall time, final test-only max force, parameter trajectory, predicted vs actual delta E, trust ratio, rank, basis change angle, fraction of probes exploring outside U, and actual effective stiffness. Primary metrics: median evaluations to reference energy threshold, performance profiles across tasks, AUC of normalized best-E-vs-NFEV, and success fraction. Report percentiles, not just means. For global exploration report distinct best basins separately rather than disguising a local solver as a global solver.

### Training cost function (starting point)

On tasks with independent reference E_ref and E_scale = max(E_initial-E_ref, eps):

```
regret_t = max(0, E_best(t)-E_ref)/E_scale
Loss = mean_over_tasks[ mean_over_evals(log1p(30*regret_t)) + log1p(30*regret_T) ]
```

The constants/weight are provisional. Test generalization on new landscapes/seeds/potentials. Meta-training may use force trajectories as *sources of initial states*, but NEVER supply forces to the energy-only optimizer; it learns from scalar evaluations obtained during its own rollouts.

**Dense feedback — the sparse-signal problem.** Evals-to-convergence yields ONE scalar per rollout: too sparse to train even a handful of coefficients. Use dense per-step proxies instead:

- **Sliding-window relative drop:** `R_i = regret_i / regret_{i+w}` (w ~ 10–50 evals). Rewards a sustained descent rate and punishes stalls; cheap to log every step. This is the primary training signal.
- **Stuck penalty:** fraction of windows with no accepted-step progress / no best-energy improvement.
- **Instability penalty (secondary, less well-defined):** zigzag — successive accepted steps nearly anti-parallel; step-size blowup; divergence to non-finite energy.

Starting combined form: `Loss = mean_i w_i·log1p(30·regret_i) + c_stuck·(stuck fraction) + c_osc·(reversal fraction)`, with the task normalization above and `w_i` optionally up-weighting later progress.

## Critical known scientific limitations

- On general high-D problems, arbitrary unknown gradients cannot be recovered from k<<N scalar probes. Any high-dimensional efficiency gain must come from exploitable temporal/spatial structure or pretrained task-distribution priors.
- Approximate geometric curvature from accepted steps can reflect zigzag noise, not valley geometry.
- Near a strict local minimum, an energy-only method cannot certify tiny residual forces cheaply, and monotone acceptance cannot cross energy barriers. Optional final force diagnostic is fine, but not counted as gradient-free optimization progress.
- An explicit hard-bath *mean curvature* does not reveal the direction of an undiscovered rare soft mode.
- Some classical force-field implementations compute forces together with energy. Benchmark **real** CPU/GPU cost; zero-order evaluations may not save much wall time compared to analytic gradients.
