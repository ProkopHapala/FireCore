# CPU eFF follow-up review — 2026-11-08

Follow-up to [eFF_optimizer_speedup.md](eFF_optimizer_speedup.md). Scope: current CPU `eFF.h`, `DynamicOpt`, C API and `pyBall/eFF.py`. No production code changed. Findings below distinguish small direct probes from source inspection; no relaxation sweep was run. Existing benchmark numbers were not remeasured.

## 1. Correctness blockers before optimizer tuning

### P0 — paired orbitals omit their internal electron–electron interaction [probed]

`cpp/common/molecular/eFF.h:527` scales kinetic energy by `echarge`; `evalEE():547` visits only distinct orbitals (`j<i`). An `e2` orbital therefore has two electrons' kinetic and nuclear energies but no Coulomb repulsion between those two electrons. This is a physical interaction between opposite spins, not the spurious self-interaction of a single electron.

Direct He probe: one nucleus Z=2, point nucleus, electrons centered on it, s=0.5 Å:

| representation | E [eV] | force on shared size [eV/Å] |
|---|---:|---:|
| one `e2`, Q=2, spin=0 | −92.3885 | −1.89786 |
| two identical orbitals, Q=1 each, spins +1/−1 | −59.8920 | +63.0952 (sum) |

The missing **32.4965 eV** equals the explicit Coulomb pair energy. At coincidence `J(s,s)=2*k_e/(sqrt(pi)*s)` in the current width convention. The opposite-spin Pauli term tends to zero for identical orbitals (up to the current distance regularization). Missing J is size dependent: it changes equilibrium sizes, stiffness and molecular energies, not just an energy reference.

**Agent task:** define `e2` as two opposite-spin electrons constrained to the same center/size, then add its internal pair energy and both size derivatives exactly once, using the existing pair helper. Do not add Coulomb self-energy to Q=1 orbitals. Verify pair–single and pair–pair interactions against explicitly expanded spins too; do not blindly scale every Pauli term by Q_i*Q_j because spin multiplicities already enter `spinij`/`qq`.

**Tests:** isolated He and H2; compare `e2` to the same constrained explicit spin pair over several positive widths, energies and shared-size forces. For two spatially separated paired orbitals, compare with all four expanded electrons. Then rerun CH4/H2O. The earlier claim that `bPaired=True` establishes physical correctness is unsupported: it changes the variational space AND currently drops an interaction. Explicit core orbitals are always paired in `builder2EFFstatic`, even when `bPaired=False`; all-electron study cases containing `e2` are affected too.

### P0 — ECP Coulomb sign, occupancy and width argument [sign/occupancy probed]

`eFF.h:777` (`evalAE_ECP`) passes `(aPar.x-aPar.z)` as Coulomb charge product, with no `-qj`. The ordinary AE path uses `-qj*qCore`. ECP currently repels electrons.

Direct H probe, r=10 Å, s=0.5 Å, Coulomb only:

| path | E [eV] | electron Fx [eV/Å] |
|---|---:|---:|
| ordinary AE, Q=1 | −1.43996 | −0.143996 |
| ECP AE, Q=1 | +1.43996 | +0.143996 |
| ECP AE, Q=2 | +1.43996 | unchanged |

The call also supplies `aPar.z` (documented/core-table **charge**, Zcore_eff) as a Gaussian **length**; `aPar.y` contains R_eff. The standalone `evalECP():755` duplicates the suspect call.

**Agent task:** correct the charge product to `-qj*(Z-Zcore)`. Resolve intended smoothing radius from the parameter convention before replacing the width; do not merely repair the sign and assume ECP physics is validated. Audit both callers and record whether radii/Pauli fit coefficients are in Å or Bohr. The wrapper multiplies Pauli energy by Hartree_to_eV but does not visibly convert lengths; units need independent verification.

**Tests:** H reduction to the ordinary bare-nucleus path with Pauli disabled; C/O large-distance limit `E*r -> -k_e*qj*(Z-Zcore)`; occupancy doubling; Cartesian/size finite differences with Coulomb alone and Pauli alone.

### P0 — ECP Pauli forces have reversed signs [probed for s-type model 0]

`eFF_LAMMPS_funcs.h:485` returns `frc=-dE/dr`, `fre2=-dE/ds`. `evalAE_ECP` subtracts `fe_` from `fsj`, and adds `dR*fr_/r` into `f`, then subtracts `f` from the electron force. Both conversions reverse the required electron force. The standalone `evalECP` repeats the convention.

Direct s-type probe, model 0, Pauli only, r=0.4 Å, s=0.5 Å, A=22.7, B=0.728, C=1.103 (synthetic parameters for derivative testing):

| component | returned force | central FD `-dE/dx`, h=1e-5 |
|---|---:|---:|
| electron x | −243.956 | +243.956 |
| electron size | +36.0615 | −36.0615 |

**Agent task:** establish one radial-force convention and fix wrapper conversions; verify all three `iECPmodel` variants individually. At r=0 the wrapper divides by `rc`; derive the appropriate limit. For p-type formulas containing `(r-C*s)^2`, the radial derivative can be nonzero at zero: that can be a genuine Cartesian cusp, requiring an explicit model/domain decision, not an arbitrary epsilon direction.

### P0 — FGO/default initialization uses negative electron occupancy [probed]

`eFF.h:342` initializes `echarge` with `QE=-1`, but current compute paths treat it as **positive occupancy**. `from_xyz_line` uses `fabs(Q)`; `readSingleFGO():1336` never assigns it. Thus legacy FGO input and a fresh `init()` without explicit charge initialization can give negative kinetic energy and wrong AE signs.

Loading the supplied `cpp/sketches_SDL/Molecular/data/H_eFF.fgo` produced `echarge=-1`, `Ek=-45.7198 eV`, and NaN total energy. The NaN also involves the legacy core fields/zero width: correcting occupancy alone is not sufficient to certify this fixture.

**Agent task:** use positive occupancy consistently; assign it explicitly in every loader and initialization path. Define how FGO serializes paired occupancy: the writer currently writes a literal coefficient 1.0 and spin, without charge. Its velocity format also omits vsize despite passing it to fprintf; the closed-shell header has four sscanf arguments but three conversions.

**Tests:** fresh one-electron initialization; same H system loaded through XYZ/FGO; positive kinetic energy; round-trip singles and pairs, including widths, spins, occupancies and vsize. If the old format cannot encode a state, reject the export explicitly rather than silently losing occupancy.

## 2. Optimizer findings that affect the previous conclusions

### P1 — stopping criterion measures acceleration, not physical force [probed]

`DynamicOpt.h:137` calculates `ff=sum((F_i*invM_i)^2)`. `DynamicOpt.cpp:306` returns it from `move_FIRE`; `eFF_lib.cpp:228` compares it with `Fconv²` and exports it as outF. `fire`, `fire_groups`, and `fire_secant` in `tests/tEFF/fire_study_lib.py` use the same weighted stopping criterion.

One-DOF probe: F=1, invM=1e-6 reports |F|=1e-6 and satisfies Fconv=1e-3 while the true force remains 1. Increasing mass must not relax the convergence tolerance.

**Agent task:** keep the internal acceleration norm for existing mixing if desired, but stop/report using the raw physical force over free DOFs. Use the same criterion for GD/MD/FIRE and all Python prototypes. Record both L2 and maximum absolute component; never include constrained reaction forces. Reevaluate at the returned coordinates: current run evaluates, moves, then tests the old force and returns the new positions with old energies/forces. Check convergence before moving, and refresh the final state on budget exhaustion.

**Benchmark interpretation:** `summarize()` already measures crossings using raw `Fu`; those individual recorded crossings are meaningful. However weighted early termination can censor traces, `converged` is unreliable, and outF from C++ is weighted. Recheck each published number's actual source instead of either accepting or discarding the entire table.

### P1 — `dr_limit` is inactive; two inconsistent width floors [source inspection + width probe]

Active `move_FIRE` only reduces dt using an acceleration norm and `f_limit`. Its displacement-limit block is commented out in the older implementation; `dr_limit=0.2` is not enforced by active FIRE. Existing velocity can therefore cause a large displacement even when the current force is small.

`run` clamps sizes to 1e-3 after moving, but `evalKinetic` clamps to `min_esize=0.1` **and mutates positions inside energy evaluation**. Probe: setting s=0.05 then evaluating returns s=0.1. This invalidates unconstrained finite differences and optimizer secants across the boundary. A clipped width also retains its incoming velocity.

**Agent task:** choose/document one width domain. Apply step admissibility before evaluation, not inside the energy function; report an active boundary and use constrained stationarity if boundary minima are allowed. For relaxation, limit the actual proposed displacement (including velocity) and keep position/velocity consistent. Test a relative width-step cap as a cheap option. Do not transplant relaxation mass/clipping heuristics into physical MD.

**Tests:** eval does not change valid input; invalid widths fail explicitly; stiff one-dimensional oscillator with preexisting velocity; maximal displacement and positive widths at aggressive dt; unclipped interior trajectories remain unchanged.

### P1 — optimizer/constraint state can survive a load incorrectly [source inspection]

`processXYZ_e():792` reallocates ff on each invocation, but calls `initOpt` only if `!opt_initialized`. Otherwise it changes dt without rebinding opt's pos/force/vel/invMasses/fixmask or its n. After another load these may point to freed allocations or retain an old size. Updating the buffer registry does not repair optimizer bindings. `processXYZ` binds arrays but does not visibly refresh opt.fixmask there. `ff.realloc` also leaves `nfix`/fixed constraint arrays from a prior system in place.

`run` zeros velocity but preserves `damping`/`lastNeg` while resetting dt to its ceiling. `initOpt` likewise does not reset those history counters. This mixes restart and continuation semantics and can make scans depend on the preceding system's optimizer history.

**Agent task:** rebind all arrays and masks after every realloc; explicitly reset/rebuild constraints. Define fresh-start versus continuation behavior and reset FIRE history consistently. Preserve the already-documented build-only atom-constraint fix as a separate task.

**Tests:** ASan, sequential loads with different DOF counts/masks and constraints; compare a fresh process with the same system loaded after another system. Test both Python entry points and one/multiple run calls.

### P1 — nonfinite checks print errors and continue [probed + inspected]

`evalEE` and `eval` print/flush on nonfinite values but do not terminate; ECP can have finite energy with NaN Cartesian force at rc=0. Check all active coordinates, widths, forces and energies and terminate with context. Do not classify a persistent negative-power/dt_min stall as convergence unless the physical force criterion also passes; return a distinct stalled status.

## 3. Cheap convergence work, after correctness gates

1. **Retest existing secant masses at a common physical-force tolerance.** This remains the cheapest promising candidate; no new algorithm is needed. Keep unit masses for the baseline and physical masses for MD only. Test H2O in 'a', 'f', 'e' separately; repair the selected path's physics first. Compare within the same representation/basin, not raw energies between different Hamiltonians.
2. **Treat secant masses as a heuristic, not an estimated Hessian diagonal.** Along a coupled trajectory, `-ΔF_i/Δx_i ≈ sum_j H_ij*Δx_j/Δx_i`. Small Δx_i can amplify cross-coupling arbitrarily. EMA/clipping do not remove that bias. Record cap saturation and reject updates with insignificant displacement. Try one isotropic mass per electron-position triple to avoid arbitrary x/y/z and molecular-orientation dependence; retain a separate width mass. A block secant `abs(-Δx·ΔF)/(Δx·Δx)` is an inexpensive experiment, not guaranteed curvature recovery.
3. **Do not force dt_min to stay large when raising dt_max.** Current dt_min=0.1*dt_max couples the two. Expose a separate minimum or floor fraction and test a small grid. First add the actual displacement safeguard; reducing the floor alone can merely turn divergence into a very slow stall.
4. **Benchmark a physically meaningful FIRE power test as an isolated ablation.** Current power is v·(invM F), whereas energy rate is −v·F. Nonuniform masses can give opposite signs. A consistent fixed-mass coordinate transform uses power v·F, velocity norm sqrt(vᵀMv), and force norm sqrt(FᵀM⁻¹F), mixing velocity toward M⁻¹F. This costs only reductions. Keep frozen DOFs out and handle changing masses deliberately; do not combine this change with secant tuning in the first comparison. Existing acceleration-space FIRE may still work empirically, but its power is not an energy-descent diagnostic.
5. **Remove the known Coulomb force/energy mismatch before fine polishing.** For the same approximant `E=A*k*g(k*r)`, its width derivative follows from scaling: pair-helper coefficient `fs=(E+r*dE/dr)/(si²+sj²)`, then accumulate `fs*si`/`fs*sj`. This can reuse the existing radial derivative and eliminate the independent gauss_p8 call; validate against FD at r=0, moderate separation and the tail. Beware cancellation of E+r*dE/dr in the far tail. No new fit or optimizer is required.

The previous finite-difference refresh cost estimate needs correction: one diagonal costs `2*n_free` full force evaluations. Relative extra cost at refresh interval R is approximately `2*n_free/R`, not automatically 1–4%. For n_free=100 and R=200 this is **100% extra evaluations**, before Python overhead. Rank methods by total evaluations and wall time, not iteration count alone.

## 4. Bounded handoff and acceptance gates

- **Physics agent:** paired internal interactions + occupancy/FGO first; separate ECP sign/derivative/units patch. Use independent FD and expanded-spin equivalence, not CPU/GPU agreement as the physics oracle.
- **Optimizer/API agent:** raw-force convergence and final-state consistency; pointer/mask lifecycle; deterministic restart; width domain and real displacement cap. Avoid broad DynamicOpt changes until its other callers are inventoried.
- **Benchmark agent:** after the above, at most a small H2O matrix: unit versus existing secant versus isotropic block variant; a few dt ceilings; common free-DOF physical tolerance (e.g. 1e-3, then 1e-4). Report terminal raw norms, energy, basin identity, boundary hits, resets, iterations, force-evaluation count and elapsed time. Repeat one rotated geometry to expose coordinate dependence. No long molecular sweep or new quasi-Newton implementation yet.

The temporary probe was compiled against the current headers and `DynamicOpt.cpp`, without using a possibly stale shared library. Source/executable: `/tmp/eff_cpu_review/probe.cpp`, `/tmp/eff_cpu_review/probe` (ephemeral). Commands from the repository root:

```sh
g++ -std=c++17 -O0 -Icpp/common -Icpp/common/math -Icpp/common/molecular -Icpp/common/dataStructures -Icpp/common/algorithms /tmp/eff_cpu_review/probe.cpp cpp/common/math/DynamicOpt.cpp -o /tmp/eff_cpu_review/probe
/tmp/eff_cpu_review/probe
```

The first compile needed `<cmath>` before `eFF.h` because `datatypes.h` uses NAN without a preceding math include. Recompilation succeeded with existing header warnings. All numerical results above are diagnostic observations of broken behavior, not passing regression tests. Convert the probes into assertions in the sanctioned test location when implementing repairs; the checkout currently has `tests/tEFF` and no `/test` directory, so resolve that documented path discrepancy before adding permanent tests.
