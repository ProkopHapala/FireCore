# ValleyOpt: minimal gradient-free optimizer discovery playground

**Status:** research skeleton / runnable prototypes, **not** vetted numerical solvers. The algorithms are original simplified experiments inspired by the references below; they are NOT exact reproductions of those publications. The current prototypes are deliberately energy-only and topology-blind.

```bash
python -m unittest -v test_valleyopt.py
python bench.py --n 12 --neval 1000
python train.py --method ridge --trials 100 --neval 300
```

Requires Python 3 and NumPy. Optional: use SciPy / PRIMA / Nevergrad externally for fair baselines; not imported by these files. See `IMPLEMENTATION_BRIEF.md` and `FIRECORE_INTEGRATION.md`.

## Interface

```python
import numpy as np
from valleyopt import RidgeTR, minimize

energy = lambda x: float(100*(x[1]-x[0]**2)**2+(1-x[0])**2)
s = minimize(energy, RidgeTR, np.array([-1.2,1.0]), max_evals=500, step=.1, seed=4)
print(s.best_e, s.best_x, s.trace.nfev)
```

Or use `ask()/tell()` directly for a FireCore C++/Python energy evaluator. No neural network and no problem-specific molecular descriptors are permitted in any optimizer. Coordinate transformation is the adapter's responsibility.

| Prototype | Concept | Stored N-dimensional vectors | Evaluations per informative action |
|---|---|---:|---|
| `MomentumWalk` | successful-direction EMA + random perpendicular exploration | 1 | 1 |
| `SecantSketch` | online normalized gradient sketch using all energy differences (NLMS/Kaczmarz style) | 1 (gradient estimate) | 1; mirror on selected failures |
| `CurvedWalk` | extrapolate turning of previous successful displacements | 2 (+ pending move) | 1 |
| `WallProbeWalk` | curved extrapolation; on wall hit probe ± transverse displacement, fit slope and curvature | ~5 (+ pending move) | usually 1, occasionally 3–4 |
| `RidgeTR` | quadratic model in rank-k moving subspace; trust ratio; active new direction | k + tangent | `2k+k(k-1)/2` probes per rebuilt model; otherwise 1 |

All algorithms retain scalar traces for diagnostics. Scalar trace length is O(T); **no optimizer stores full N-vectors for the last T iterations**. In production, traces should be optionally disabled.

### What is experimental/incomplete

- Many decisions are currently hard-coded; offline training tunes only fixed scalar hyperparameters. No MLP controller or symbolic search engine yet.
- `RidgeTR` is a simple **local quadratic model inside an adaptive subspace**, not OMoRF. It deliberately omits OMoRF's full-space subspace interpolation and geometry management.
- `WallProbeWalk` is a new heuristic based on the discussion, not a paper implementation. It estimates local transverse slope/curvature at a *failed forward point*. It may need many probes to discover an unrelated soft mode in large N.
- `SecantSketch` is based on local *linear* response; distant probes cause bias and state staleness. Add trust-aware forgetting/curvature before treating it as reliable.
- `CurvedWalk` estimates the curvature of the **accepted trajectory**, not necessarily the valley centerline. Don't conflate the two.
- `RidgeTR` does not model the hard bath explicitly yet, though the wall probes in `WallProbeWalk` provide the relevant measurements.
- All variants are monotonic **for their incumbent**. This is a local minimization benchmark, not barrier-crossing global optimization or true simulated annealing.
- There is no mathematically justified convergence guarantee. In some landscapes they can be much worse than Powell, COBYQA, FIRE, and L-BFGS.

## Smoke test observations (seed 3, 350 energy evaluations)

| Problem | Momentum | Secant | Curved | Wall | Ridge |
|---|---:|---:|---:|---:|---:|
| 2D Rosenbrock, final incumbent E | 0.065 | 0.130 | 0.166 | 5.5e-10 | 9.1e-25 |
| 12D rotated curved valley, final incumbent E | 1.189 | 1.051 | 1.169 | 1.176 | 1.104 |

For the rotated valley, E(x0) = 1.2 and E_min = 0. Landscape and optimizer random seeds are deliberately independent; a shared seed artificially aligns the first random trial with the hidden valley and corrupts the benchmark. Results do **not** establish performance; they demonstrate that successful 2D local fitting does not automatically imply efficient discovery of useful directions in high-dimensional space. Results can change significantly with scaling, seeding and hyperparameters. Details in `bench.py`.

## MMFFsp3_loc observations (tests/tMMFF/run_valleyopt.py, energy-only, equal eval budgets)

Run from `tests/tMMFF`: `LD_PRELOAD=$(gcc -print-file-name=libasan.so) ASAN_OPTIONS=detect_leaks=0 python3 run_valleyopt.py --xyz data/xyz/xylitol`. `ArrayEnergy` wraps `mmff.DOFs`/`mmff.eval_no_omp` (≡ `mmff.eval()` ≡ `run()` outE after fixes below); eval ~0.01–0.05 ms. Isolated-molecule runs use `nPBC=(0,0,0)` which now explicitly disables PBC even when the xyz carries `lvs`.

**Bugs found & fixed while validating the scalar energy** (finite-difference checks; this work's real product so far):

1. `MolWorld_sp3::eval()` called `evalLJQs_ng4_PBC_omp()` unconditionally — wrong NB path for non-PBC, and that routine never accumulated `E` → E/F inconsistent by construction. `eval()` now mirrors `eval_no_omp()` branching (`bPBC`/`bNonBondNeighs`/`bExclusion2`).
2. `setNonBondStrategy` was called after `ffl.eval()` inside `eval()` → the *first* `eval()` call computed the bonded NB-subtraction with stale flags → garbage E (e.g. -16983); moved before, matching `eval_no_omp` order.
3. Per-atom `evalLJQs_*_atom*` routines sum every pair's energy twice (once per endpoint) while forces are correct → `E_reported = E_b + 2·E_NB`. Fixed by returning the half-share `E*0.5` (standard partition).
4. **`ndims` bound to `W.ff` while DOF buffers point to `W.ffl`** — different DOF layouts (`natoms+npi` vs `natoms+nnode` pi-vectors): Python saw only 162 of ffl's 210 DOFs on nHexadecan → invisible `pipos[4..19]` mutated by `run()` → E depended on hidden state (restoring `DOFs` didn't restore E). Fixed `init_buffers` to report ffl dims; `nDOFs` is now 54/96/210 for the three molecules below.
5. **`eval_no_omp`/`run_no_omp` never normalized pi-vectors** (`ffl.eval()` does): pi-terms scale with `|pipos|` → unphysical unbounded DOF; scipy Powell legitimately drove `|pipos|≈5e4` → E=-54716. Fixed by `normalizePis()` in both paths. The optimizer now sees the projected unit-sphere manifold; tangent forces keep `|pi|≈1` in FIRE.

Finite-difference consistency after fixes: `numgrad(E) = -fDOFs` to rel ~6e-8 at relaxed states (non-PBC); at heavily clashed starts FD truncation on the r⁻¹² wall dominates (not a bug). Under PBC, `bClampNonBonded` caps close-contact forces at `FmaxNonBonded` while E keeps the true LJ → E/F deliberately disagree on clashed pairs (documented caveat, PBC path only).

**Metrics:** `E_inc` = best energy found within the eval budget. `dE%` = `100·(best−Eref)/(E0−Eref)` — fraction of the gap to the reference remaining; negative = deeper than reference. `FIRE ref` = minimum *visited* energy over a long FIRE run (iters shown; early exit = converged to Fconv). `FIRE*` = FIRE best-visited at the same eval count — force-based, not information-equivalent. Milestone columns = first eval index reaching `Eref + L·(E0−Eref)`. The run also writes `<mol>_convergence.npz`/`.png` in `tests/tMMFF/` (gitignored): best-visited relative gap vs eval index on log-y — the actual convergence-speed curves, regenerated with `python run_valleyopt.py --plot <name>`.

2000 evals, step=0.05, seed=3, non-PBC, corrected energy. `nDOFs` counts only **free** DOFs — the ffl layout allocates a pi-vector per node, but all-sp3 molecules have `Ksp=Kpp=0` (verified: E invariant under pi scramble) so dead pi DOFs are masked out via `ArrayEnergy(free_mask)`; nHexadecan keeps 6 active pi nodes (2 COOH groups):

| molecule | free nDOFs | E0 | FIRE ref (iters) | FIRE* 2k | powell | secant | ridge | momentum | curved | wall |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| propandiol | 39 (of 54) | 0.176 | 0.1075 (18669, conv) | 0.108 (0.9%) | **0.117** (13.4%) | 0.127 (27.8%) | 0.132 (35.2%) | 0.147 (57.9%) | 0.139 (45.7%) | 0.143 (51.7%) |
| xylitol | 66 (of 96) | 0.927 | 0.3338 (21669, conv) | 0.360 (4.4%) | **0.476** (24.0%) | 0.516 (30.7%) | 0.557 (37.6%) | 0.761 (72.0%) | 0.656 (54.4%) | 0.692 (60.4%) |
| nHexadecan-2COOH | 168 (of 210) | 11.499 | -1.227 (50000, ~best) | -1.138 (0.7%) | **0.299** (12.0%) | 1.982 (25.2%) | 4.280 (43.3%) | 9.903 (87.5%) | 7.444 (68.1%) | 8.778 (78.6%) |

Honest picture with a *correct* smooth energy: the landscape is far less "glassy" than the buggy energy suggested — FIRE now converges on propandiol/xylitol and is strong per-eval everywhere (0.7–4.4% gap at 2k). Among energy-only methods Powell is the benchmark to beat; `secant` remains the best prototype (best on xylitol, second on the others) — consistent with its amortized gradient-sketch design. The prototypes' remaining niche is the rugged/poor-start regime, not smooth descents; a valleyopt→FIRE hybrid is the realistic deployment.

## eFF observations (tests/tEFF/run_valleyopt.py, energy-only, equal eval budgets)

*(Historical/exploratory — eFF's physical soundness is questionable; MMFFsp3_loc above is the trusted validation target.)*

Run from `tests/tEFF` against the ASan build: `LD_PRELOAD=$(gcc -print-file-name=libasan.so) CPP_RECOMPILE=0 ASAN_OPTIONS=detect_leaks=0 python3 run_valleyopt.py` (Build-opt libs are stale — missing `setAutoMass`). `ArrayEnergy` wraps `eff.pDOFs`/`eff.eval`; eval ~15 µs (forces computed anyway, so energy-only saves no wall time here).

| start (atoms fixed, ~36 free DOFs) | E0 | FIRE, 20k iters | best energy-only, 3k evals | FIRE polish after |
|---|---:|---:|---:|---:|
| H2O, load_a (auto electrons) | -749.7 | -1259.0 stalled | ridge **-1616.3** | -1715.3 |
| CH4, load_a | -350.9 | -555.0 stalled | secant **-860.7** | -905.6 |
| CO, load_a | -822.4 | -2007.7 | secant -2033.3 | — |
| H2O_spins, load_e (good start) | -1473.9 | -1718.7 converged | ~-1670 (20–27% gap left) | — |

From badly-placed auto-generated electrons the stochastic walkers cross barriers/plateaus that stall FIRE for 20k iters, and land in basins FIRE-polishable to deep minima — the target regime of this project. From a good start on a smooth descent, SciPy Powell still dominates every prototype (0.64% of gap at 2000 evals on H2O_spins): the prototypes are not yet competitive in the smooth/quadratic regime.

## Dimensionality scaling on eFF (load_a, atoms fixed, 3000 evals, seed 3)

Percent of energy gap (E0→Eref) REMAINING after 3000 evals; negative = went below FIRE's 20k reference. `secant` is the best prototype on every real eFF surface tested — consistent with its design (amortized gradient sketch rather than per-step random probing).

| system | nfree | FIRE (forces) | powell | momentum | secant | curved | ridge | wall |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| H2O_spins (load_e, good start) | 36 | 0.0 | 0.64 | 24.8 | 21.1 | 27.4 | 20.0 | 21.3 |
| H2O (load_a) | 36 | ~100 stalled | -15.8 | 15.5 | -64.0 | -7.3 | -70.1 | -58.0 |
| CH4 (load_a) | 36 | ~stalled | — | -35.1 | -149.7 | -51.6 | -55.4 | -129.4 |
| CO (load_a) | 32 | 0.0 | — | 37.8 | -2.2 | 20.1 | 22.2 | 31.3 |
| adipic C6H10O4 (load_a) | 256 | 0.67 | -12.5 | 71.1 | 29.2 | 57.7 | 38.7 | 61.4 |
| xylitol C5H12O5 (load_a) | 288 | 5.5 | 37.2 | 79.4 | 42.4 | 69.4 | 53.2 | 78.0 |

Clear dimensionality dependence: at ~36 DOFs the stochastic methods escape FIRE stalls and win; at ~280 DOFs FIRE and Powell pull far ahead — random energy-only directions in 300-D space mostly fail, and the information cost of probing grows with the stiff subspace. `secant`'s amortized sketch degrades most gracefully. The prototypes' value on large systems is escaping bad basins, not beating FIRE on open descents.

## Scientifically meaningful comparisons

Prefer area under the best-energy-vs-NFEV curve, time/evaluations to fixed energy thresholds, and convergence success rate across many held-out landscapes. Evaluate at equal energy-evaluation budgets; also measure wall time (energy eval sometimes computes forces anyway in FireCore). Preserve a final high-accuracy force check as a diagnostic even for energy-only solvers. For model learning use distinct *molecules/potential parameters* in the held-out set, not only new starts from molecules used for training.

## References (concepts, not reproductions)

- Xiangning Chen et al., **Symbolic Discovery of Optimization Algorithms**, NeurIPS 2023 — program search, simplification, resulting Lion rule: https://arxiv.org/abs/2302.06675
- James Gross & Geoffrey Parks, **Optimization by moving ridge functions**, Engineering Optimization (2021) — local reduced subspaces, quadratic ridge trust-region models, two interpolation sets: https://arxiv.org/abs/2007.04893
- Krzysztof Choromanski et al., **From Complexity to Simplicity: Adaptive ES-Active Subspaces for Blackbox Optimization**, NeurIPS 2019 — online active subspaces and adaptive probing inside/outside: https://arxiv.org/abs/1903.04268
- Ilya Loshchilov, **A Computationally Efficient Limited-Memory CMA-ES for Large Scale Optimization**, GECCO 2014 — a small set of historical directions achieves O(mN) sampling/storage: https://arxiv.org/abs/1404.5520
- M. J. D. Powell, **An efficient method for finding the minimum ... without calculating derivatives**, Computer Journal (1964) — derivative-free direction updates and 1D minimization: https://doi.org/10.1093/comjnl/7.2.155
- Powell's BOBYQA, NEWUOA; modern PRIMA project — interpolating quadratic trust region baseline (not the same as the above low-rank prototype): https://www.libprima.net/
