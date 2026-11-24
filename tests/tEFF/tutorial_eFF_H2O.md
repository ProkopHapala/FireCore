# Relaxing Water with eFF — a Guided Tour of Optimization on a Noisy Energy Landscape

*A didactic tutorial + reproducible experiments. Audience: a clever student who knows what an
atom and an energy minimum are; the "Developer notes" sections go deeper.*

All experiments below are reproduced by one script (run from `tests/tEFF/`):

```bash
python3 tutorial_H2O_relax.py dimer dimer_bad forces water15 dtscan
```

It writes figures + `.xyz` structures to `artifacts/` and prints every number quoted here.

**Two libraries, side by side.** This branch contains *two* versions of the eFF library built
from the same Build directory: `eFF_lib` (the ported version — Tomas's physics + Prokop's
optimizer and bug fixes) and `eFF_lib_origZ` (**Tomas's untouched original**, its sources kept
in `cpp/libs/Molecular/origZ/`). Which one loads is chosen at runtime:

```bash
EFF_LIB=eFF_lib        python3 tutorial_H2O_relax.py dimer   # ported (default)
EFF_LIB=eFF_lib_origZ  python3 tutorial_H2O_relax.py dimer   # original backup
```

This lets you A/B the same script against both physics/optimizer stacks — used in §7.5.

---

## 1. What is eFF, in one paragraph

In ordinary molecular mechanics, atoms are balls connected by springs; electrons are hidden
inside the springs' parameters. **eFF (electron Force Field) does something bolder: it makes
the electrons themselves particles.** Each electron is a spherical floating Gaussian wavepacket
with three numbers to optimize: its center position `x, y, z` and its width `s` (size). Pauli
repulsion between electrons is computed from the quantum-mechanical overlap of the Gaussians —
same-spin electrons avoid each other, opposite-spin pairs can share a region and form what we
*see* as bond pairs and lone pairs.

Consequence: a water molecule in eFF is not 3 atoms, it is **3 nuclei + 8 electron wavepackets**
(4 e-pairs: 2 bond pairs + 2 lone pairs). Everything — covalent bonds, lone pairs, hydrogen
bonds — emerges from one energy function `E(positions, sizes)`. If the physics is right, the
hydrogen bond should form *by itself* when we minimize E. That is what we test here.

## 2. The thing we want to see: the hydrogen bond

A hydrogen bond is not a stick between atoms — it is a directional, mostly electrostatic
attraction: the (partially positive) **H of a donor water points into a lone pair of the
acceptor's oxygen**. In eFF language: the donor's O–H bond axis should point at one of the
acceptor's e-pair Gaussians sitting *opposite* to its own O–H bonds.

So our test is very concrete: start two waters **misoriented** and check whether relaxation
rotates the acceptor so that an e-pair ends up facing the donor's H.

## 3. Why is this hard for an optimizer?

Three reasons, in increasing subtlety:

**(a) The forces span ~4 orders of magnitude.** The forces that stretch an O–H bond or push an
electron into place are ~eV/Å. The torque that rotates a molecule into H-bond alignment is
~meV-scale — a thousand times weaker. An optimizer that treats all degrees of freedom equally
spends most of its time polishing strong modes and barely notices the weak ones.

**(b) The energy landscape has many minima.** The electron wavepackets can settle into several
different arrangements, and each arrangement creates a different valley for the nuclei. Which
valley you fall into depends on the path you took — relaxation is **path-dependent**. Starting
from a "wrong" geometry can trap you in a wrong-looking minimum forever (we show this below).

**(c) There is a force jitter floor.** Some e-pair directions are nearly flat — the electron
can wander without changing E much. A dynamical optimizer keeps rattling along those
directions: after the atoms are fully converged, `|F|` on the electron DOFs still jitters
around ~0.5. **Convergence measured on total |F| is therefore unreachable** — you must measure
convergence on *atom* forces `|Fa|` only. Measured directly:

```
iters      E [eV]     |Fa|     |Fepos|  |Fesize|      <- python3 tutorial_H2O_relax.py forces
   500    -891.449   0.17271   34.834    0.4307
  4500    -898.726   0.00280    7.093    0.9126
  8500    -899.113   0.00039    2.420    0.0449
 12500    -899.217   0.00023    0.978    0.0090
 19500    -899.383   0.00021    0.467    0.0120       <- atoms done at ~10k; epos jitters forever
```

The atoms are converged to |Fa|≈2e-4 while the electron-center force keeps bouncing around 0.5–1.
That residual is noise on flat directions, not a signal that atoms should still move.

## 4. The engine: FIRE, explained with a ball

Optimization = finding the bottom of a valley in `E(x)`. Steepest descent (roll straight
downhill, step ∝ force) is slow in curved valleys — it zig-zags. **FIRE** instead lets the
point behave like a **ball**: it accumulates velocity downhill (inertia carries it smoothly
along the valley), with two controls:

- **Steering (`damping` = mixing coefficient α):** each step, velocity is rotated a fraction α
  toward the force direction: `v = (1−α)·v + α·f̂·|v|`. Speed is kept, direction is corrected.
- **Emergency brake:** if `<v|f> < 0` — the ball started moving *uphill* — velocity is zeroed
  and it restarts (with dt shrunk, then grown back).

**Counterintuitive finding of this study:** on ordinary smooth energy surfaces, small α (~0.1,
nearly free-rolling ball + occasional brake) is fastest. On the eFF landscape α≈0.3 converges
~2× faster. Why? The electron jitter floor from §3(c) acts like random noise constantly
knocking the ball sideways; a larger α scrubs the mis-aimed part of the velocity every step
instead of letting it orbit. Pure inertia is optimal only when the valley floor is smooth.

## 5. The preconditioner: `autoMass`

FIRE uses a single time step for all DOFs, but the eFF DOFs have wildly different stiffnesses:
e-sizes and e-pair separations are stiff (curvature k ~ 10³), atom translations are soft.
A single dt must serve the stiffest mode → everything else crawls.

`setAutoMass(1, kref)` fixes this with a **diagonal preconditioner** — *not* physical masses.
It watches the trajectory, estimates each block's stiffness as `k_b = |ΔF_b| / |Δx_b|`
(a secant estimate — how much the force changed when this block moved), and sets

```
invMasses[i] = kref / (kref + k_b)
```

Stiff blocks get small invMass → they accelerate slowly → take small effective steps;
soft blocks keep invMass ≈ 1 → stride normally. One dt now serves all modes. Same energy
landscape, same minima — it only rescales the effective step size per block.
(This is distinct from `bMass=True` in `initOpt`, which assigns real atomic masses for MD.)

## 6. Putting it together: `relax_auto`

`eff.relax_auto()` (pyBall/eFF.py) wraps FIRE + autoMass into a self-tuning driver:

1. `setAutoMass(1, 100)` — normalize stiffness across all blocks.
2. Run FIRE in 2000-iter chunks; after each chunk measure `|Fa|` (atom forces).
3. If `|Fa|` explodes (> 10·√na) → the dt guess was too hot → `dt ×= 0.3`.
4. If `|Fa|` stops improving (less than ×0.7 per chunk) → we are at the jitter-limited
   plateau for this dt → `dt ×= 0.5` (automatic polish phase).
5. Stop when `|Fa| < FconvA` (default 1e-3), confirmed once more at `dt × 0.3`.
6. `damping=0.3` — the α from §4, ~2× faster on this landscape than the FIRE default.

Zero hand-tuned parameters: it finds the stable dt by probing, and it knows *what* converged
(atoms) vs *what* jitters (electrons).

## 7. Experiments

### 7.1 The good case: `dimer_90.mol2` — a correct hydrogen bond

The acceptor water starts rotated 90° away from the donor. Run
`python3 tutorial_H2O_relax.py dimer`:

```
=== dimer_90 ===
E: -687.66 -> -899.30 eV in 8000 iters (dt_final=0.03)
|Fa|=0.00060  |Fepos|=5.0967  |Fesize|=0.0062
O...O = 2.849 A,  H(donor)...O(acc) = 2.032 A
  e-pair 15 on acceptor: r(O-e)=0.44 A, cos(lp, H-bond dir)=+0.87
  e-pair 14 on acceptor: r(O-e)=0.44 A, cos(lp, H-bond dir)=+0.87
```

![dimer_90](artifacts/dimer_90_relaxed.png)

The acceptor **rotated**: O···O = 2.85 Å, donor H sits 2.03 Å from the acceptor O, and an
acceptor e-pair (indices 14/15 — the two spins of one lone pair) points at the incoming H with
**cos = +0.87**, nearly dead-on. This is the textbook H-bond geometry, produced by the model
itself, not imposed.

### 7.2 The trap: `dimer_bad.mol2` — path dependence is real

Same two waters, but the acceptor starts *fully reversed* — its own H's pointing at the donor.
`python3 tutorial_H2O_relax.py dimer_bad`:

```
=== dimer_bad ===
E: -699.32 -> -899.10 eV in 10000 iters
O...O = 2.999 A,  H(donor)...O(acc) = 2.086 A
  best acceptor e-pair cos = +0.74 (partial), closest pairs cos = -0.44 (facing AWAY)
```

![dimer_bad](artifacts/dimer_bad_relaxed.png)

Fully converged (|Fa|=4e-4) — but into a **wrong local minimum**: the molecule refined the bad
orientation rather than rotating. It is even ~0.2 eV shallower, yet the barrier to rotate out
is too big for a local optimizer. **Lesson: convergence ≠ correctness.** A converged relaxation
only proves you found *a* minimum, not *the* minimum. Tomas's original library converges to the
same trap — this is eFF landscape physics, not an optimizer bug.

### 7.3 dt is not just a step size — it picks the basin

`python3 tutorial_H2O_relax.py dtscan` — identical start, only dt changes:

```
    dt     E [eV]      |Fa|   O...O
  0.02   -899.594   0.00061   2.741
  0.05   -899.462   0.00014   2.837
  0.10   -899.040   0.00168   2.859
  0.30   -899.160   0.00067   2.843
```

Every run converges atoms — but to *different* minima ~0.5 eV apart. A larger dt doesn't just
walk faster; it takes different paths through the landscape and can hop over small barriers.
This is why "tune dt" is a legitimate (if unsatisfying) exploration tool — and why
`relax_auto` treats dt as a schedule, not a constant. For *global* optimization you want the
opposite: deliberately hot dynamics to explore, then a quench (Tomas's `dt=5` heating +
`dt=0.1` anneal protocol in `minimum.py` is exactly that).

### 7.4 The big one: 15 waters

`EMolecules/H2O/15.mol2` = 3 parallel chains of 5 waters (45 atoms, 120 electrons, 615 DOFs).
`python3 tutorial_H2O_relax.py water15`:

```
15 waters: E -5706.0 -> -6741.842 eV in 12000 iters
|Fa|=0.00062 |Fepos|=0.9310 |Fesize|=0.0607  (na=45, ne=120)
H-bonds (H...O < 2.6 A, not covalent): 12
```

![water15](artifacts/water15_relaxed.png)

12 directional H-bonds form *along* the chains (O···O ≈ 3.0 Å, H···O ≈ 2.1 Å). Note the chains
do **not** bond between each other: every water's free H points perpendicular to the chain
plane (a symmetric local minimum — no cross-chain bond is geometrically accessible without
breaking the symmetry). Same trap as the dimer, at larger scale: the input topology decides
which minima are reachable.

### 7.5 A/B: the same dimer on Tomas's original library

Same start (`dimer_90.mol2`), same driver (`relax_auto`), only the library changes —
`EFF_LIB=eFF_lib_origZ`. On the original library `setAutoMass` does not exist, so
`relax_auto` automatically skips the preconditioner (the `EFF_LIB` env var + `TAG` suffix
keep the outputs separate):

```
EFF_LIB=eFF_lib_origZ TAG=dimer_90_origZ python3 tutorial_H2O_relax.py dimer
```

| | ported `eFF_lib` | original `eFF_lib_origZ` |
|---|---|---|
| iters used | **8 000** | 60 000 (hit the cap, **not converged**) |
| final \|Fa\| | 6e-4 | ~0.2, still bouncing |
| dt trajectory | 0.1 → 0.03, done | decayed to the 2e-4 floor, still jittering |
| E [eV] | −899.30 | −899.91 (still sinking) |
| O···O / H···O [Å] | 2.85 / 2.03 | **2.49 / 1.68** |
| lone pair → donor H | cos +0.87 | **cos +1.00** |

Two opposite lessons in one table:

1. **The preconditioner is what makes relaxation practical.** Without autoMass, the original
   library burns 7× more iterations and *still* has not settled — dt decays to the floor
   while |Fa| hovers ~0.2. This is exactly the hand-tuning pain (hundreds of thousands of
   steps) that motivated the port.
2. **The original physics finds a *deeper, tighter* H-bond.** Once it gets there, origZ's
   minimum is −899.9 eV with O···O = 2.49 Å and a *perfectly* aligned lone pair (cos = 1.00),
   vs the ported −899.3 eV at 2.85 Å. The two Hamiltonians are not identical (see
   PHYSICS_DIFFS.md: the ported version adds the e2 self-Coulomb term and different Pauli/ECP
   parameters), so their minima differ — part of this gap is also that origZ's run was still
   descending at the cutoff. The qualitative physics — "H into the lone pair" — is confirmed
   on **both**.

## 8. What we learned, condensed

| Claim | Evidence |
|---|---|
| eFF forms proper directional H-bonds | §7.1: lp→H cos = +0.87, O···O = 2.85 Å |
| Relaxation is path-dependent; wrong starts trap | §7.2: anti-aligned minimum, fully converged |
| Total-\|F\| convergence is unreachable | §3(c): \|Fepos\| floor ≈ 0.5 forever |
| Converge on \|Fa\| (atoms) instead | \|Fa\| → 1e-4…6e-4 in every experiment |
| autoMass is a step-size preconditioner, not masses | §5: `invMasses = kref/(kref+k_b)` |
| Larger FIRE α helps on noisy landscapes | damping=0.3: ~2× fewer iters (5k vs 8k to 1e-3 on 15-mer) |
| dt picks basins; hot→quench explores | §7.3 + Tomas's anneal protocol |
| Typical cost with `relax_auto` | dimer ~8–10k, 15-mer ~12k, stiff H2O_spins ~14k iters |

## 9. Developer notes

- **Entry points** (pyBall/eFF.py): `processMol2`/`processXYZ_e` (load+build, use `nstepMax=0`
  to skip built-in relax), `initOpt` (dt, damping, f_limit), `setAutoMass`, `run`,
  `relax_auto`, `eval()` (E+forces at current `pDOFs`), `getBuffs` (numpy views on live
  `apos/epos/esize/fDOFs/pDOFs/invMasses` buffers — set `bVel=True` first for `invMasses`).
- **Line-search methods fail on this landscape.** Both `run_GH` (BB1/Anderson) and a Python
  L-BFGS over the live `pDOFs`/`fDOFs` views were tried: every descent direction needs step
  size → 0 because the e-size/Pauli walls make the directional derivative huge and rough.
  Damped dynamics (FIRE) is the right engine here; do not reach for quasi-Newton.
- **`_mol2_sizes` heap-corruption fix** (this branch): the C++ builder adds lone-pair
  electrons (`nepair=2` for O) that the old Python size guess missed → `epos` buffer overrun.
  `processMol2` now does a NULL-output dry-run and reads true `na`/`ne` from `ndims`.
- **Convergence criterion** is on `|Fa|` for a reason — see §3(c). If you ever see code
  converging on total |F|, it will run forever (or "converge" on a fluke).
- **Basin depth vs speed trade-off:** `relax_auto` quenches as soon as |Fa| stalls. A
  deliberately hot phase (dt≈0.3, ~10–20k iters) before the decay finds deeper minima
  (15-mer: −6745.5 eV vs −6741.8 quenched). For global search, wrap `relax_auto` in an
  anneal loop (perturb/heat → relax → keep best) — the primitive is ready.
- **Physics caveats in this branch's Hamiltonian** (see PHYSICS_DIFFS.md): Tomas's
  different-spin Pauli term `−rho2·DT·S22·si2sj2/(S22−1)` has an unguarded pole at S22=1 —
  a plausible source of force spikes. His ECP/`rho2`/`paramsSi` parameters were fitted
  against a different version of the energy function; refitting may be needed.
- **A/B testing vs Tomas's original**: `EFF_LIB=eFF_lib_origZ` selects the original-library
  build; its untouched sources live in `cpp/libs/Molecular/origZ/` and build as target
  `eFF_lib_origZ` in the same Build dir (`make eFF_lib_origZ`). A git tag `zverina-orig` also
  marks his last commit. `relax_auto` detects the missing `setAutoMass` symbol and degrades
  gracefully (dt-probe only). See §7.5 for the dimer A/B; the dimer trap and jitter floor
  reproduce identically on both libs — they are physics, not port bugs.

## 10. Reproduce

```bash
cd tests/tEFF
python3 tutorial_H2O_relax.py dimer        # ~10 s   -> figs + "=== dimer_90 ===" table
python3 tutorial_H2O_relax.py dimer_bad    # ~15 s   -> the trap minimum
python3 tutorial_H2O_relax.py forces       # ~20 s   -> |Fa| vs |Fepos| jitter table
python3 tutorial_H2O_relax.py dtscan       # ~20 s   -> dt -> basin scatter
python3 tutorial_H2O_relax.py water15      # ~60 s   -> 15-mer lattice + 12 H-bonds
EFF_LIB=eFF_lib_origZ TAG=dimer_90_origZ python3 tutorial_H2O_relax.py dimer
                                           # ~60 s   -> same dimer on the backup lib (no autoMass)
```

Outputs: `artifacts/*.png` (structures), `*.xyz` (relaxed geometries for `xyz_view`).
All numbers in this document come from these runs.
