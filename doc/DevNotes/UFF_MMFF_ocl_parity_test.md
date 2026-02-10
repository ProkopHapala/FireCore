# UFF/MMFF OpenCL Parity Test (Bonded-only, CPU vs GPU)

## Goal
Provide a single, repeatable procedure to check **force parity** between CPU (C++) and GPU (OpenCL) for **bonded-only** terms of UFF and MMFF. Energies and non-bonded terms are out of scope here.

## Scope and assumptions
- Parity target: forces (F) only, bonded terms only (bonds/angles/dihedrals/inversions or MMFF bonded set).
- Non-bonded and GridFF **disabled**.
- CPU reference: C++ implementations exposed via `pyBall.MMFF_multi`.
- GPU path: C++ OpenCL via `iParalel=2` (UFF/MMFF) or PyOpenCL (`pyBall/OCL/UFF.py`) for kernel isolation.
- Use small molecules (e.g., formic acid, water, methane) for fast runs and clear diffs.

## Quick commands (UFF, bonded-only)
Using the existing C++-backed driver:
```bash
# UFF: CPU vs GPU forces, bonded-only
python3 tests/tUFF/test_UFF_multi.py \
  --molecule cpp/common_resources/mol/formic_acid.mol2 \
  --tol 1e-3
```
Notes:
- The script already calls `setSwitches2(NonBonded=-1, GridFF=-1)` and sets UFF component flags; only bonded terms run.
- CPU path uses `iParalel=0`; GPU path uses `iParalel=2`. The script compares `F_gpu - F_cpu` and reports max |ΔF|.
- Energies are not compared in this script (forces only).

## Quick commands (PyOpenCL UFF kernel isolation)
For kernel-by-kernel checks without the C++ layer:
```bash
python3 tests/tUFF/test_UFF_ocl.py \
  --molecule cpp/common_resources/mol/formic_acid.mol2 \
  --components bonds,angles,dihedrals,inversions \
  --tolerance 1e-6
```
Notes:
- Uses `pyBall/OCL/UFF.py` to upload topology/params and run bonded kernels.
- Good for debugging kernel argument order/packing before C++ integration.

## MMFF (bonded-only) using the same C++ interface
There is no dedicated MMFF parity script yet, but you can reuse `pyBall.MMFF_multi` directly:

1) Initialize with MMFF, UFF off:
```python
from pyBall import MMFF_multi as ff
ff.init(bMMFF=True, bUFF=False)
```

2) Disable non-bonded / grid and enable bonded MMFF pieces:
```python
ff.setSwitches2(NonBonded=-1, GridFF=-1, MMFF=1, Angles=1, PiSigma=0, PiPiI=0)
```

3) Run CPU vs GPU and compare forces (pseudo-flow):
```python
ff.setTrjName("trj_mmff.xyz", savePerNsteps=1)
# CPU
ff.fapos[:,:] = 0.0
ff.run(nstepMax=1, dt=0.02, Fconv=1e-6, ialg=2, damping=0.1, iParalel=0)
F_cpu = ff.fapos.copy()
# GPU
ff.fapos[:,:] = 0.0
ff.run(nstepMax=1, dt=0.02, Fconv=1e-6, ialg=2, damping=0.1, iParalel=2)
F_gpu = ff.fapos.copy()
# Compare
import numpy as np
max_diff = np.max(np.abs(F_gpu - F_cpu))
print("max |ΔF| =", max_diff)
```

## Toggles and conventions
- **Non-bonded off**: `setSwitches2(NonBonded=-1, GridFF=-1)`; for UFF also keep `SubtractBondNonBond=-1`, `ClampNonBonded=-1` in `setSwitchesUFF`.
- **Bonded components**:
  - UFF: `DoBond/DoAngle/DoDihedral/DoInversion = 1` for terms you want; `DoAssemble=1` to sum `fint` -> `fapos`.
  - MMFF: `Angles=1` (and other Pi terms as needed) via `setSwitches2`.
- **Parallel flag**: `iParalel=0` (CPU serial), `iParalel=1` (CPU OMP where available), `iParalel=2` (GPU OpenCL).
- **Buffers**: UFF tests pad CPU angAtoms to stride 4 for GPU comparison; ensure you don’t repack when reusing helpers.

## Known gaps / cautions
- `tests/tUFF/test_UFF_multi.py` does **not** compare energies—forces only.
- UFF non-bonded/GridFF kernels are not bound/run in the current OpenCL setup; keep them disabled for parity.
- MMFF parity isn’t scripted; use the manual snippet above until a consolidated script is added.
- Angle/dihedral/inversion buffer packing must match kernel expectations (see `doc/DevNotes/UFF_muli_plan.md` for historical mismatches).

## Suggested consolidation (future work)
- Add `--ff {uff,mmff}` and `--tol` to `tests/tUFF/test_UFF_multi.py`, reuse the same CPU-vs-GPU force diff logic, and branch the init/switches accordingly.
- Optionally add energy comparison after ensuring GPU energies are written/downloaded.

## References
- `tests/tUFF/test_UFF_multi.py` – C++ UFF CPU vs GPU force scan
- `tests/tUFF/test_UFF_ocl.py` – PyOpenCL UFF vs CPU buffers/forces
- `doc/DevNotes/UFF-GPU-vs-CPU-Debug-Guide.md` – kernel arg/offset checklist
- `doc/DevNotes/UFF_muli_plan.md`, `UFF_ocl_plan_new.md`, `UFF_ocl_plan_updated.md` – wiring plans and risk notes
- `doc/DevNotes/MD_MMFF_OCL_notes.md` – MMFFsp3 OpenCL MD (rotational/RATTLE) learnings and toggles (useful for MMFF setup)
