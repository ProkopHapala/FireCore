# Practical workflow to validate GPU vs CPU UFF

- End‑to‑end steps

  - Build and run the test:
    - `python3 tests/tUFF/test_UFF_multi.py`
  - Ensure correct switches before run:
    - In Python: `uff.setSwitchesUFF(DoBond=1, DoAngle=1, DoDihedral=1, DoInversion=1, DoAssemble=1, SubtractBondNonBond=-1, ClampNonBonded=-1)`
    - In C++: `MMFFmulti_lib.cpp::setSwitchesUFF()` sets both host and OCL flags.
  - Fetch buffers:
    - `getBuffs_UFF()` prints sizes and exposes arrays, useful for pre/post inspection.
  - For GPU forces:
    - Call `uff.download(bForces=True)` then use `uff.fapos.copy()` for comparison.

- Where to add targeted prints for UFF

  - In `UFF.cl`:
    - Temporarily set `#define DBG_UFF 1` and set `IDBG_*` to specific interaction indices to get per‑DOF traces.
    - Watch headers `[GPU][BOND]`, `[GPU][ANGL]`, `[GPU][DIH ]`, `[GPU][INV ]` and per‑DOF lines for mismatch clues.
  - In `OCL_UFF.h::setup_kernels()`:
    - Leave `printf("OCL_UFF::setup_kernels().task_...")` banners.
    - Keep `OCL_checkError(err, "kernel.argN name")` after every binding to immediately pinpoint arg order/type errors (as with the angles fix.
  - In `MolWorld_sp3_multi.h::eval_UFF_ocl()`:
    - Prints already in place to confirm flow, kernel preparation, download, and final energy.

- Expected argument orders (audit checklist)

  - Bonds: 17 args total after scalars, ending with `fint`.
  - Angles: After `pbc_shifts` pass `neighs`, `neighCell`, `npbc`, then `fint`, then `Ea_contrib`.
  - Dihedrals: After `pbc_shifts` pass `neighs`, `neighCell`, `npbc`, then `fint`, then `Ed_contrib`.
  - Inversions: After `hneigh` pass `fint`, then `Ei_contrib`.
  - Any `CL_INVALID_ARG_SIZE` points to a position where a buffer/scalar got swapped – recheck against `UFF.cl`.

- Common pitfalls and resolutions (UFF)

  - Wrong arg order for angles/dihedrals/inversions:
    - Symptom: `CL_INVALID_ARG_SIZE` with arg index in `OCL_checkError`.
    - Fix: Reorder bindings in `OCL_UFF.h::setup_kernels()` to match `UFF.cl`.
  - Missing energy buffers:
    - Symptom: Build or runtime errors when kernel expects `__global float* Ea_contrib` etc.
    - Fix: Allocate `ibuff_Ea/Ed/Ei` in `OCL_UFF.h::realloc()` and bind them (no NULL).
  - Wrong `fint` offsets:
    - Symptom: Shape mismatch or forces scrambled vs CPU.
    - Fix: Recompute `i0bon/i0ang/i0dih/i0inv` = `{0, +2*B, +3*A, +4*D}`, bind consistently, and ensure `nf` total matches CPU planning.
  - No forces after GPU:
    - Symptom: CPU vs GPU differs with zeros on GPU.
    - Fix: Confirm `OCL_UFF.h::eval()` enqueues assembler if `bUFF_assemble` is true; confirm `OCL.h::upload()` and `OCL.h::download()` are called in run path.
  - Python retrieval mismatch:
    - Symptom: AttributeError for `gpu_aforces`.
    - Fix: Use `uff.fapos` after `uff.download(bForces=True)` as per `MMFFmulti_lib.cpp::init_buffers_UFF()` exposure.

- Minimal quick‑test recipe

  - In Python:
    - `uff.init(..., bMMFF=True, bUFF=True)`
    - `uff.setSwitches2(NonBonded=-1, SurfAtoms=-1, GridFF=-1)`
    - `uff.setSwitchesUFF(DoBond=1, DoAngle=1, DoDihedral=1, DoInversion=1, DoAssemble=1, SubtractBondNonBond=-1, ClampNonBonded=-1)`
    - `uff.getBuffs_UFF()`
    - `uff.run(nstepMax=1, dt=0.02, Fconv=1e-6, ialg=2, iParalel=2)`
    - `uff.download(bForces=True)`
    - `forces_gpu = uff.fapos.copy()`
  - Compare to CPU by running `iParalel=0` similarly, then `compare_results(...)`.

If you want, I can save these as two files, e.g.:

- `doc/Markdown/GPU-CPU-Debugging-and-OCL-Interface.md`
- `doc/Markdown/UFF-GPU-vs-CPU-Debug-Guide.md`

Let me know the preferred filenames/locations and I’ll add them.