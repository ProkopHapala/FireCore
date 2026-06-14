---
 description: Flat surface MMFF coupling notes (Hamaker LJ9-3 / Morse)
---

## Purpose
Implementation and debugging notes for the flat surface interaction in `MolWorld_sp3` used by MMFF relaxations (Hamaker LJ 9-3 and Morse surface potentials). Includes where to configure from Python/CLI and test scripts.

## Key code locations
- C++ core
  - `cpp/common/molecular/MolWorld_sp3.h`: flat-surface state (`surfFlat_pos0`, `surfFlat_normal`, `surfFlat_REQ`, `surfFlat_K`, `surfFlat_mode`, `bPlaneSurfForce`), setters `setSurfFlatPlane`, `setSurfFlatParams`, evaluation `evalSurfFlat`, trajectory saving in `run`/`run_omp`.
  - `cpp/common/math/Forces.h`: surface kernels `getHamakerLJ93` (LJ 9-3, 1D), `getMorseSurface` (1D Morse). Force is along the surface normal; minimum at `z = REQij.x` (see combine rule below).
  - `cpp/common/molecular/libMMFF.h`: C API exports for `setSurfFlatPlane`, `setSurfFlatParams`, `setTrjName`, `run`, `saveXYZ`.
  - `cpp/libs/Molecular/MMFF_lib.cpp`: sample surface kinds, C API wrapper dispatch to `MolWorld_sp3::run`/`run_omp`.
- Python bindings
  - `pyBall/MMFF.py`: ctypes bindings for plane/params (`setSurfFlatPlane`, `setSurfFlatParams`), trajectory (`setTrjName`), `run` wrapper, switch configuration.
- Tests/examples
  - `tests/tMMFF/run_relax_surf.py`: minimal relaxation of H2O over flat plane; sets plane at z=-5, Hamaker params, enables `SurfAtoms`, saves trajectory `h2o_flat_trj.xyz`.
  - `tests/tMMFF/run_sample_surf.py`: sampling harness including flat-surface kinds.
  - `tests/tMMFF/run.sh`: builds `libMMFF_lib.so` and runs the test scripts.

## Parameter / mixing behavior
- `setSurfFlatParams(mode, REQ, K)` sets:
  - `mode`: 1 = Hamaker LJ 9-3, 2 = Morse surface.
  - `REQ.x` (z0/R), `REQ.y` (epsilon), `K` (Morse stiffness).
- `evalSurfFlat` calls `combineREQ(surfFlat_REQ, atomREQ, REQij)`; current combine rule sums radii: `REQij.x = R_plane + R_atom`, epsilon multiplied, charges/Hb combined. The LJ9-3 minimum is at `z = REQij.x`; using a nonzero plane radius lifts the minimum for larger atoms.
- To make the minimum depend only on atom radius (smaller atoms closer), pass `REQ.x = 0` from Python (no C++ change needed), e.g. `mmff.setSurfFlatParams(mode=1, REQ=(0.0, eps, 0, 0))`.

## Force direction and logging
- Forces from `getHamakerLJ93` / `getMorseSurface` are strictly along the normalized surface normal; any x/y components in debug prints come from other force terms in `fi` before adding `f`.
- Debug print in `evalSurfFlat` shows `pi`, surface-only `f`, and combined REQij.

## Trajectory saving
- `MolWorld_sp3::run` and `run_omp` write frames when `trj_fname` is set and `savePerNsteps>0`; `setTrjName(path, savePerNsteps, bDel, nPBC)` from Python/C sets it. OMP path now mirrors serial path.

## Known pitfalls / fixes
- If plane at z=0 and atom z=0: LJ9-3 singularity; `getHamakerLJ93` early-return guard was commented out for testing—avoid z≈0 or raise plane.
- If you see zero atoms: ensure `bMMFF=True` in `mmff.init` and `SurfAtoms=1` in switches.
- If surface pushes larger atoms up: plane radius contributes to `REQij.x`; set plane `REQ.x=0` to avoid it.
- ASan double-free traces point to `MolWorld_sp3` cleanup; unresolved—ignore during current surface testing.

## Quick usage (Python)
```python
from pyBall import MMFF as mmff
mmff.setVerbosity(verbosity=1, idebug=1)
mmff.setSwitches(NonBonded=1, MMFF=1)
mmff.init(xyz_name="data/xyz/H2O", surf_name=None, bMMFF=True)
mmff.setSurfFlatPlane(pos0=(0,0,-5), normal=(0,0,1))
mmff.setSurfFlatParams(mode=1, REQ=(0.0, 1.0, 0, 0), K=1.6)  # plane radius 0 so H < O
mmff.setSwitches(SurfAtoms=1, GridFF=0, NonBonded=1, MMFF=1)
mmff.setTrjName("/home/prokop/git/FireCore/tests/tMMFF/h2o_flat_trj.xyz", savePerNsteps=1, bDel=True, nPBC=(1,1,1))
mmff.run(nstepMax=5000, dt=0.05, Fconv=1e-3, ialg=2, omp=True)
mmff.saveXYZ("h2o_flat_relaxed.xyz", "relaxed flat surf", 1)
```

## Outstanding items
- ASan double-free in `~MolWorld_sp3` still present; needs proper ownership cleanup.
- Consider optional clamp/softening for z→0 in LJ9-3 to avoid singularity.
- Add GUI plane visualization (transparent quad) in `MolGUI` (planned).

## Recent fixes (mol/mol2 loading & init)
- `load_mol` now tolerates V2000 counts line on line 3 or 4, parses bond lines as (atom1, atom2, order), assigns conformers to atoms, adds bonds to confs, and backfills `conf.nbond` when zero so `checkNumberOfBonds` passes for V2000 inputs (e.g., `xylitol.mol`).
- `MolWorld_sp3::init` guard now checks `builder.atoms.size()` (not `ffl.natoms`) so it doesn’t abort before `makeFFs()` populates FF structures; covers both mol and mol2 paths.
