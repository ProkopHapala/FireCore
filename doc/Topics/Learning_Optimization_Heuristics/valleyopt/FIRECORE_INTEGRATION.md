# FireCore `prokop` branch integration notes for coding agent

Repository: https://github.com/ProkopHapala/FireCore/tree/prokop

**Verified from branch source on 2026-11-08:**

- `pyBall/eFF.py`: `eval()` calls the C++ evaluator and returns `double`; `getBuffs()` binds `pDOFs` and `fDOFs` via writable NumPy views, with `nDOFs`. It also provides `initOpt`, `run`, and native FIRE baseline `ialg=2` in example code.
- `tests/tEFF/fire_study_lib.py`: instrumented FIRE for eFF, with `load_e`, `load_a`, `snapshot`, `restore`. The eFF DOF order is `[apos(3*na) | epos(3*ne) | esize(ne)]`; its `load_e` returns a fixed-DOF mask. This is the most straightforward first physical test.
- `pyBall/MMFF.py`: `eval()` returns a `double`, `getBuffs()` exposes `DOFs`, `fDOFs` for MMFF, and `getBuffs_UFF()` exposes `apos`, `fapos` for the UFF layout. `run()` supports the built-in force-based baseline. **Do not assume the MMFF and UFF DOF layouts are identical.**
- `pyBall/OCL/UFF.py`: GPU UFF supports `nSystems` batches and `upload_positions` for float32 XYZ positions. Use later for parallel training; it needs a proper energy-readback adapter, not the CPU buffer wrapper here.
- `cpp/common/math/DynamicOpt.h`: `DynamicOpt` is the CPU optimization class with FIRE parameters, velocities, force, timestep and optional fixing masks.
- `cpp/common/molecular/RARFF.h`: rigid-atom reactive FF includes translational and quaternion rotational state (`RigidAtom.pos`, `qrot`). `cpp/apps/RFF/RFFapp.cpp` is an application. I have **not identified a generic Python RARFF energy-evaluator binding**. Add a small adapter exposing a flat reduced DOF state and scalar energy, with rotation handled by a quaternion retraction (outside the optimizer).

## Minimal eFF usage (after native initialization)

```python
import sys, numpy as np
from pyBall import eFF as eff
from valleyopt import RidgeTR, minimize
from firecore_adapter import ArrayEnergy

# Example test setup, from repository:
# sys.path.insert(0, '/path/to/FireCore/tests/tEFF')
# import fire_study_lib as fs
# na, ne, fixed = fs.load_e('H2O_spins.xyz')
# The surrounding test harness should choose exact cwd and input path.
# eff.getBuffs() must already have been called after system initialization.
fixed = np.zeros(eff.nDOFs, dtype=bool)  # REPLACE using actual fixed indices
positive = np.arange(3*(eff.na+eff.ne), eff.nDOFs)  # electron sizes
objective = ArrayEnergy(eff.pDOFs, eff.eval, free_mask=~fixed,
                        scale=0.1, positive_indices=positive, min_positive=1.e-3)
z0 = np.zeros(objective.ndim)
s = minimize(objective, RidgeTR, z0, rank=2, step=.3, max_evals=1000)
objective.commit(s.best_x)   # restore best energy-evaluated state
```

**Critical:** the example `fixed` mask is placeholder, not a physically correct setting. Obtain it from `fire_study_lib.load_e()` / `load_a()` or the test's actual constraints. `positive_indices` refers to indices in the *full flattened* pDOFs array, not free-subspace indices. The example coordinate transform is linear; near positivity bounds use an alternative log parametrization for electron sizes.

`ArrayEnergy` assumes `eval()` reads `pDOFs` without remapping their ownership. This should be verified in the exact eFF executable mode before benchmarking; also reset any hidden optimizer velocities/constraints if evaluations depend on them. The optimized vector `z` lives in scaled free-DOF coordinates; the observer only sees `z` and energies.

## Proposed adapter API for C++ integration

```cpp
struct EnergyAdapter {
    int nfree();                       // optimization dimension after constraints
    void set_coordinates(const double* z); // map optimizer z into physical state
    double eval_energy();             // side-effect-safe scalar evaluation
    void commit_best(const double* z);     // restore best configuration
};
```

Build a Python adapter first for correctness. Then port only the O(kN) vector operations, scalar state and ask/tell logic to C++, or fuse into GPU kernels. The energy-only optimizer must not use `fDOFs` even if the engine computes forces internally.

## Physical constraints / fair cost measurement

- eFF has atom positions, electron positions **and electron sizes**. Sizes must remain positive; electron/core overlap, fixed coordinates and other native constraints may require validation or reparametrization. A naive vector step can produce invalid states.
- RARFF rotations are on SO(3), not unconstrained quaternion components; map a 3-vector of local rotations through an exponential-map retraction and compose with the current quaternion. Quaternion normalization alone is not an ideal unconstrained coordinate chart. The optimizer sees abstract coordinates only.
- UFF/MMFF need proper topology preparation, exclusion/fixed handling, energy mode, and masses (for the FIRE baseline); do not use the topology as input to the optimizer.
- Use one fixed, repeatable molecular structure per objective instance. Every `energy(z)` evaluation must set the complete state from an immutable reference, not add deltas to the last trial's mutated state.
- Compare **energy-call count and wall time**. FireCore `eval()` may compute forces at essentially no extra cost, so energy-only comparisons could be less impressive than a purely expensive gradient case.
- Use identical starting geometry across methods. Save original x before native FIRE runs because its native integrator mutates the state.
- Coordinate scaling matters: raw eFF coordinate subsets can have very different physical stiffnesses. Scaling is allowed in the objective adapter (or use a simple common length scale); do not let the learned controller inspect atom types or topology.

## Recommended first physical benchmark sequence

1. eFF: use `tests/tEFF/fire_study_lib.py`; minimize electrons with fixed nuclei first, then allow atoms. Compare to native FIRE (`ialg=2`) and independent SciPy Powell/PRIMA COBYQA baselines. Validate final forces using `eff.fDOFs` only for scoring.
2. MMFF/UFF: small molecules with nontrivial torsions, then adsorbed molecules, then flexible molecular assemblies. Compare fixed evaluation budgets and the same perturbation seeds.
3. RARFF: add a stable energy interface / quaternion retraction; begin with rigid-atom positions-only, add orientations later.
4. Batch candidate evaluation on GPU only after confirming single-evaluation state semantics. Batch *independent optimizer replicas*, or implement a proper frozen-center `ask_batch`; do not parallelize sequentially dependent ask/tell steps naively.

No repository files have been modified by this package.
