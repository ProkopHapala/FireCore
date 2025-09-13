
# General GPU vs CPU debugging and OCL.h interface

- Scope: How to debug GPU (OpenCL) vs CPU (C++) implementations, and how to build robust OpenCL interfaces against `OCL.h`.
- Audience: Any subsystem doing CPU/GPU parity checks, not limited to UFF.

1) Implementing debug prints in GPU (OpenCL) and CPU (C++)

- OpenCL Kernel Debugging (printf)

  - Use compile‑time macros to gate prints. In `cpp/common_resources/cl/*.cl`, add at top:
    - `#ifndef DBG_<MODULE>; #define DBG_<MODULE> 0; #endif`
    - ID filters like `#define IDBG_<KIND> (-1)` to print only one work‑item.
  - In kernels, print only when a single work‑item runs to avoid flooding:
    - Header prints: `if ((DBG_UFF!=0) && get_global_id(0)==0) { printf(...); }`
    - Focused DOF prints keyed by index:
      - Example in `UFF.cl`:
        - Bonds: inside loop, `if (DBG_UFF && ib_dbg==IDBG_BOND) printf("[GPU][BOND-DOF] ...");`
        - Angles: `if (DBG_UFF && iang==IDBG_ANGLE) printf("[GPU][ANGL-DOF] ...");`
        - Dihedrals: `if (DBG_UFF && id==IDBG_DIH) printf("[GPU][DIH -DOF] ...");`
        - Inversions: `if (DBG_UFF && ii==IDBG_INV) printf("[GPU][INV -DOF] ...");`
  - Print what helps parity checks:
    - Input parameters for first 64 interactions (as already done in UFF).
    - Local variables: indices, computed vectors, energies (E) and Enb (if subtracted), forces fi/fj/fk/fl, intermediate coefficients.

- C++ CPU Debugging (printf, asserts)

  - Prefer loud failures and minimal masking:
    - Use `printf` in hot paths (as per Global Coding Rules) for quick traces.
    - Use `assert()` for invariants (e.g., valid index bounds, non‑null buffer pointers after allocation).
  - For OpenCL host binding errors, annotate each step:
    - After every `OCL.h::useArg()/useArgBuff()/useArg_()` call, use `OCL_checkError(err, "tag")`.
    - Include kernel name and positional index in the tag to quickly locate mismatches.
  - Example:
    - In `OCL_UFF.h::setup_kernels()`, bind each arg and immediately do:
      - `OCL_checkError(err, "evalAngles.arg14 neighs")`
    - This immediately showed the arg index causing `CL_INVALID_ARG_SIZE`.

- How to toggle debug prints

  - Add CMake or compile flags to define macros:
    - For GPU: build program with `-DDBG_UFF=1 -DIDBG_ANGLE=0` etc.
    - For CPU: guard prints with a module‑level verbosity flag or compile flag.

2) Building robust OpenCL interfaces with OCL.h

- Required steps and the exact functions/macros to use

  - Build program and create tasks:
    - `OCL.h::buildProgram(path, program)`, then `OCL.h::newTask("kernel_name", program, ...)`.
  - Cache task pointers:
    - `OCLtask* task = getTask("kernel_name");`
  - Allocate GPU buffers using the OCLsystem API:
    - `ibuff_X = newBuffer("name", count, sizeof(type), 0, CL_MEM_READ_ONLY/WRITE_ONLY/READ_WRITE);`
    - For optional outputs, allocate actual buffers (do not pass NULL) if the kernel signature expects `__global float*` – e.g., `Ea_contrib`, `Ed_contrib`, `Ei_contrib`.
  - Bind kernel args with error granularity and strict order:
    - Set sizes: `task->local.x = nloc; task->global.x = N + nloc - (N % nloc); task->global.y = nSystems;`
    - `OCL.h::useKernel(task->ikernel);`
    - For scalars: `err |= useArg(int_or_float)`.
    - For buffers by index: `err |= useArgBuff(ibuff_X)`.
    - For raw pointers when absolutely needed: `err |= useArg_(ptr, nbytes)` (rare; avoid for OpenCL buffers).
    - Always `OCL_checkError(err, "kernel.argN <label>")` after each.
  - Caveats and how to avoid them
    - Arg order must exactly match kernel signature.
      - Otherwise you get `CL_INVALID_ARG_SIZE` or silent incorrect binding.
    - Don’t pass NULL placeholder pointers for `__global T*` parameters.
      - Either remove the argument in kernel (if optional) or allocate a real device buffer and pass it.
    - Ensure nf/offsets are consistent between CPU buffer layout and GPU kernel logic.
      - Compute `i0bon`, `i0ang`, `i0dih`, `i0inv` once and reuse consistently.
    - Buffers with N=0: allocate at least size 1 to satisfy OpenCL (e.g., PBC shifts).
    - Set tasks only after buffers and dimensions are ready.
      - Prepare tasks in `OCL_UFF.h::makeKernels()`, bind args in `OCL_UFF.h::setup_kernels()` called after `OCL_UFF.h::realloc()`.
    - Do not mutate `task->args` directly if you use OCL_MM style:
      - Prefer `OCL.h::useKernel()/useArg*()` so you can check each binding step with `OCL_checkError`.
  - Catching and pinpointing errors
    - Always include the positional arg index in the error label: `"kernel.arg<N> name"`.
    - If you hit an error, compare the binding order with the exact order in the `.cl` file.
    - Use the kernel header prints to confirm the values received by the GPU.

- Typical failure signatures and fixes

  - `CL_INVALID_ARG_SIZE`: Wrong order or wrong arg type (e.g., scalar where buffer expected).
  - All zeros in output forces:
    - Kernel not enqueued; check flags and `task->enque()` calls.
    - Assembler kernel missing; check `bUFF_assemble` and that assemble task runs last.
  - Shape mismatches in Python:
    - Reconcile host buffer shapes and nf/index offsets; update `MMFF_multi.py::getBuffs_UFF()` shapes.

Document 2: UFF‑specific GPU vs CPU validation workflow

- Scope: How to validate and debug UFF GPU against CPU using your current files and scripts. Includes exact functions/paths and expected arg orders.

1) Relevant files and functions

- GPU Kernels and Signatures
  - `cpp/common_resources/cl/UFF.cl`
    - `evalBondsAndHNeigh_UFF(natoms, npbc, i0bon, bSubtractBondNonBond, Rdamp, FmaxNonBonded, apos, fapos, neighs, neighCell, pbc_shifts, neighBs, bonParams, REQs, bonAtoms, hneigh, fint)`
    - `evalAngles_UFF(nangles, i0ang, bSubtractAngleNonBond, Rdamp, FmaxNonBonded, angAtoms, angNgs, angParams1, angParams2_w, hneigh, REQs, apos, pbc_shifts, neighs, neighCell, npbc, fint, Ea_contrib)`
    - `evalDihedrals_UFF(ndihedrals, i0dih, SubNBTorsionFactor, Rdamp, FmaxNonBonded, dihAtoms, dihNgs, dihParams, hneigh, REQs, apos, pbc_shifts, neighs, neighCell, npbc, fint, Ed_contrib)`
    - `evalInversions_UFF(ninversions, i0inv, invAtoms, invNgs, invParams, hneigh, fint, Ei_contrib)`
    - Print macros and DOF prints already present.

- OpenCL Host Interface
  - `cpp/common/OpenCL/OCL_UFF.h`
    - `makeKernels(cl_src_dir)` creates tasks via `newTask()`.
    - `realloc(...)` allocates all buffers, including optional `ibuff_Ea/Ed/Ei`, safe `ibuff_pbcshifts`.
    - `setup_kernels()` binds args using OCL_MM style:
      - Uses `useKernel()`, `useArg()`, `useArgBuff()`, `OCL_checkError()` per arg.
      - Computes `i0bon=0; i0ang=i0bon+nBonds*2; i0dih=i0ang+nAngles*3; i0inv=i0dih+nDihedrals*4`.
      - Sets local/global sizes.
    - `eval()` enqueues in order:
      - Bonds, Angles, Dihedrals, Inversions, Assemble (guarded by `bUFF_*` flags).
    - Flags:
      - `bUFF_bonds`, `bUFF_angles`, `bUFF_dihedrals`, `bUFF_inversions`, `bUFF_assemble`, `bSubtractNB`, `bClampNonBonded`.
    - Lifecycle:
      - After `realloc()`, call `setup_kernels()` to bind args. Use `bKernelPrepared` to avoid re‑binding unless data changes.

- Integration with World/Library
  - `cpp/common/molecular/MolWorld_sp3_multi.h`
    - `eval_UFF_ocl(int niter)`:
      - If `!uff_ocl->bKernelPrepared`, calls `uff_ocl->setup_kernels()`.
      - Loops `uff_ocl->eval()`, then `uff_ocl->download(ibuff_energies)` and returns `E_tot` of system 0.
    - `pack_uff_system(...)`, `upload(...)`, `download(...)` used by library.
  - `cpp/libs_OCL/MMFFmulti_lib.cpp`
    - `init_buffers_UFF()` exposes buffers to Python (`apos`, `fapos`, `REQs`, `hneigh`, `fint`, param/topology buffers).
    - `setSwitchesUFF()`: sets both host flags (`W.ffu`) and OCL flags (`W.uff_ocl`) for bonds/angles/dihedrals/inversions, NB subtraction, clamp.
    - `run(...)`:
      - iParalel=2 path:
        - `W.pack_uff_system(0, ...), W.upload(true, false, false, true), W.eval_UFF_ocl(nstepMax), W.download(true,false)`.

- Python test harness
  - `tests/tUFF/test_UFF_multi.py`
    - `init()` – sets `bMMFF=True`, `bUFF=True`.
    - `setSwitches2()` – disables unrelated subsystems.
    - `setSwitchesUFF()` – enables UFF component flags, NB subtraction/clamp.
    - `getBuffs_UFF()` – reads CPU buffers (canonical shapes) via `getBuff()/getIBuff()`.
    - `run_uff(use_gpu)` – selects iParalel (GPU=2), runs library, then for GPU calls `uff.download(bForces=True)` and reads forces from `uff.fapos`.
    - `compare_results(...)` – compares CPU vs GPU forces.

