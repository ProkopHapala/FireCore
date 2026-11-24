# Physics differences: prokop vs Tomas Zverina eFF

Three-way reference:
- **origZ** = Tomas's original (commit `201caa6a` / tag `zverina-orig`, sources in `cpp/libs/Molecular/origZ/`)
- **prokop** = `/home/prokop/git/FireCore` @ `85cf4f32`
- **ported** = current `FireCore-zverina` working tree = Tomas's physics + prokop's bugfixes & optimizer

Energy unit eV, length Å throughout. `aPar = {x:Z_nuc, y:R_eff(sQ), z:Zcore_eff(sP), w:cP(A)}`.

---

## 1. `cpp/common/molecular/InteractionsGauss.h` — pair-interaction kernels

### 1a. Different-spin Pauli term, `addPauliGauss_New` (spin<=0 branch)
| | formula | notes |
|--|--|--|
| origZ | `E = -rho2·DT·S22·si2sj2 / (S22-1)` | **attractive** term; extra `si2sj2` factor; **singular at S22=1**; uses file-global `double rho2 = -0.15` |
| prokop | `E = -rho·DT·S22 / (S22+1)` | finite everywhere; `rho=KRSrho.z` |
| ported | **origZ kept** (file untouched) |

### 1b. Same-spin Pauli (spin>=0 branch)
Same formula in both: `E = S22·DT·(-rho·S22 + rho-2)/(S22²-1)`.
- prokop guards `|S22²-1| < 1e-12` (clamps denominator); origZ unguarded → pole at S22=1.
- ported = origZ (**unguarded** — candidate source of the force spikes seen in H2O_spins).
- prokop also adds an `idebug`-gated pair debug print; cosmetic.

### 1c. `addPauliMultiGauss` — Tomas-only
```cpp
E = eConst * sum_i S(si, sj)   // overlap of electron sj with SET of core Gaussians paramsSi
```
Fields `paramsSi = {1.0,0.5,0.3}`, `eConst=1` live in `eFF.h` (origZ/ported). Not in prokop.

### 1d. `erfx_e6` (fastmath.h, not InteractionsGauss.h but same physics role)
| origZ | old fit — **wrong derivative at r=0** (|r|-cusp ⇒ finite-force kink at coincident charges) |
| prokop | 2026-11-08 constrained refit, `dy(0)=0` |
| ported | **prokop's refit** (CPU fastmath.h + OpenCL eFF.cl both) |

---

## 2. `cpp/common/molecular/eFF.h` — energy assembly

### 2a. Atom–electron Coulomb (`evalAE`, `evalAE_MultiGauss`)
| term | origZ / ported | prokop |
|--|--|--|
| nucleus↔e attraction | `addCoulombGauss(dR, 0, sj, ..., qCore*-qj)` — **point nucleus** | `addCoulombGauss(dR, aPar.y, sj, ...)` — smeared by core width |
| core↔e repulsion (bCoreCoul) | `addCoulombGauss(dR, aPar.y, sj, ..., qj*aPar.z)` | same |
| core Pauli | `addPauliGauss_New(dR, aPar.y, sj, ..., qj*aPar.z*0.5)` | same + `aPar.y>1e-8` zero-width guard (ported has it) |
| `qCore` | `bCoreCoul ? Z : Z-Zc` | same |

Tomas's comment in origZ: *"there should be theoretically 0 ... However, it works better with aPar.y — FIX IT"* — he knows point-nucleus is a hack; prokop's version uses `aPar.y`.

### 2b. Atom–atom repulsion (`evalAA`) — **structurally different**
| origZ / ported | explicit 4-term screened interaction: `Z_i·Z_j` point-point **plus** nucleus↔core-cloud terms `-Z_j·Zc_i`, `-Z_i·Zc_j` (point↔Gaussian, width `aPar.y`) **plus** core↔core `Zc_i·Zc_j` (Gaussian↔Gaussian) |
| prokop | single screened point charge `(Z_i-Zc_i)(Z_j-Zc_j)` — no core-cloud smearing |

⇒ origZ softens the short-range repulsion wall by the finite core density; prokop's is stiffer.

### 2c. Electron–electron (`evalEE`)
Identical structure (`qij=qi*qj` Coulomb, `addPauliGauss_New` with `spinij`, `sc=qq`).
**ported adds prokop's fix**: internal self-Coulomb of paired orbitals
`E += CoulombGauss(0, si, si, qi*(qi-1)/2)` — ~+27–32 eV per e2 orbital, **not in origZ or prokop-evalEE-j-loop** (it was simply missing). This alone shifts all paired-orbital energies.

### 2d. `evalAE_MultiGauss` — Tomas-only path
Dispatched when `bMultiGauss` (coreMode `'m'`): uses `addPauliMultiGauss` instead of `addPauliGauss_New`. Kept in ported, with prokop's fsj/fs_junk swap fix applied to its core-Coulomb term (origZ had the swap bug here too; Tomas had only fixed `evalAE`).

### 2e. ECP paths (`evalECP`, `evalAE_ECP`) — prokop's fixes now in ported
| bug in origZ | fix in ported |
|--|--|
| Coulomb width `aPar.z` (charge number as size!) | `aPar.y` (R_eff) |
| charge `+(Z-Zeff)` → **repulsive** | `-qj*(Z-Zeff)` attractive |
| `fsj -= fe_`, `f += fr_` (both signs wrong vs FD) | `fsj += fe_`, `f -= fr_` |
| `else` → `addPauliGauss_New` for ALL non-eCP atoms | `aPar.z>0` guard (bare H ⇒ 0/0 NaN before) |

⚠ Tomas's `paramsSi`/`eConst` and Si ECP fits were tuned **against the buggy signs** — ECP results will differ and need refit.

### 2f. `echarge` semantics
| origZ | `_realloc0(echarge, ne, QE)` QE=-1 default; `readSingleFGO` never assigns echarge |
| ported/prokop | occupancy magnitude, init `1.0`; `readSingleFGO` sets `fabs(c)` incl. paired-e block |

### 2g. Constants / params
| | origZ / ported | prokop |
|--|--|--|
| `default_AtomParams2` R_eff | **C 0.179, N 0.2, O 0.1** | C 0.621427, N 0.1, O 0.167813 |
| `min_esize` | **0.0** (floor inactive) | 0.1 |
| `KRSrho` | member Vec3d {1.125,0.9,-0.2}, `setKRSrho` setter | `constexpr static` same values |

### 2h. XYZ / FGO IO
- Header: origZ `na,ne,core %i %i %c` (+coreMode, required by `processXYZ_e` parser); prokop `na,ne %i %i`. Ported = origZ.
- Electron line: origZ writes integer charge (`-2`,`-1`); spin parsed from `e+/e-/e2` name suffix. Prokop writes `%4.1f` float where **sign = spin** and parses spin from charge sign. Ported = origZ (name-suffix convention).
- origZ has `bOnlyPairs` filter + `Setrho2`, `copyAtomForces`, `to_xyz` element-name table — kept.

### 2i. Safety additions (prokop → ported)
`!isfinite(Etot)` print in `eval()`, `!isfinite` checks in evalEE Pauli, `verbosity>3` per-pair debug prints, `eff_dbg_*` pair hooks. Non-physics (diagnostics only).

---

## 3. `cpp/libs/Molecular/eFF_lib.cpp` — exports & drivers

### Tomas-only (all kept in ported)
`setKRSrho`, `setGaussPauliSizes`, `setEnergyConstant`, `setParsECandPS`,
`EvalTwoElectrons`, `EvalHE`, `FindMinrho2` (rho2 scan driver), `processMol2`,
`setAtomParams`, `printSwitches`, `relaxed_scan`, 9-arg `setSwitches`,
`bFixedAtoms` global + fixed-atom logic in `processXYZ_e`, `na,ne,core` XYZ parser.

### prokop-only (added in ported)
Optimizer API: `setAutoMass/setMassMetric/setOptParam/run_GH`, `opt*` step API,
`initOpt(dt,damping,f_limit,bMass)`; run() fixes: pre-move convergence test,
stale-buffer rebind, NaN abort, esize floor via `ff.min_esize` (=0.0 on Zverina!).

### Physics-neutral but behavioral
- `processXYZ` EOF last-conf fix (ported).
- `init_buffers()` `buffers[k]=v` — stale Python pointer fix.

---

## 4. Parameter consistency warning (pre-existing, not port-introduced)

`ATOM_PARAMS` in `eFF_ocl.py` hardcodes R_eff O 0.167813 / C 0.621427 — while
Tomas's CPU `default_AtomParams2` uses O 0.1, C 0.179, N 0.2. This mismatch
already existed on Tomas's branch (his old `eFF_ocl.py` had the same values);
the GPU path therefore never matched his CPU defaults. Decide which R_eff set
is canonical for his fits.

## 5. What the ported build actually is
`ported CPU = Tomas's Hamiltonian` (rho2 term, singular Pauli denominators,
point-nucleus attraction, 4-term core-cloud evalAA, MultiGauss path, his params)
`+ prokop's correctness fixes` (e2 self-Coulomb, ECP sign/width, echarge, erfx refit, swap bugs, NaN guards)
`+ prokop's optimizer machinery` (autoMass, massMetric, run_GH, DynamicOpt v2).

Deliberately NOT taken from prokop: his different-spin Pauli formula, his S22-1
guard, his point→smeared evalAE/evalAA conventions, min_esize=0.1 — those would
invalidate Tomas's fitted parameters.
