

## USER

look at @test_UFF_multi.py @test_MMFFsp3_pyOCL.py @test_UFF_ocl.py @MD_test_utils.py @UFF.py @UFFbuilder.py @MMFF.py @MolecularDynamics.py 

@relax_multi_mini.cl @UFF.h 

We have there two forcefields UFF and MMFF, we have them implemented/interfaced using pyOpenCL and we can also call it over C++ interface to @MMFF_multi.py 

I'm a bit lost what is the status, of testing, did we achieved parity? What are the relevant files and code pathways?

Analyse the relevant codes and tell me what is complete (probably) and what is incomplete. 

Do we have anywhere in this repo andy .md files which may contain andy notes on testing parity between different implementation of these forcefields?

I think I did it some time a go, but I do not remember much 

---

## GPT-5.1-Codex

Here’s the current state from the referenced files and docs:

**What’s implemented and apparently complete**
- **C++ UFF/MMFF multi wrapper**: `pyBall/MMFF_multi.py` exposes all C++ switches (MMFF vs UFF, CPU serial/OMP vs GPU OpenCL) and buffer accessors; run/scan paths are wired to the C++ `MMFFmulti_lib.cpp` dispatcher (CPU paths and GPU path `iParalel=2`). @pyBall/MMFF_multi.py#424-329 shows the exposed init/switch/run/scan hooks that route to C++.
- **C++ OpenCL UFF backend**: `cpp/libs_OCL/MMFFmulti_lib.cpp` sets up buffers, synchronizes host/GPU flags, and dispatches CPU vs GPU UFF via `iParalel` (0/-1/1 CPU, 2 GPU). GPU eval uses pack/upload → `eval_UFF_ocl()` → download/unpack, and scan batching is implemented. @cpp/libs_OCL/MMFFmulti_lib.cpp#143-329
- **GPU UFF kernel wrapper**: `cpp/common/OpenCL/OCL_UFF.h` has full buffer allocation, kernel setup/binding (bonds/angles/dihedrals/inversions/assembly, nonbond, grid), and eval. It mirrors CPU offsets (`i0dih/i0inv/i0ang/i0bon`) and allocates per-component energy buffers. @cpp/common/OpenCL/OCL_UFF.h#211-582
- **PyOpenCL UFF path (standalone)**: `pyBall/OCL/UFF.py` builds topology via `UFF_Builder`, allocates/upload buffers, maps atom→force indices, and runs kernels for bonded pieces (and optional nonbond). Buffer composition for angParams is handled. @pyBall/OCL/UFF.py#51-244

**Tests and parity coverage**
- **CPU vs GPU via C++ wrapper**: [tests/tUFF/test_UFF_multi.py](cci:7://file:///home/prokop/git/FireCore/tests/tUFF/test_UFF_multi.py:0:0-0:0) exercises the unified `pyBall.MMFF_multi` interface, sets UFF switches, runs `run()` with `iParalel=0` vs `2`, and has a `scan()` CPU vs GPU batch force comparison with tolerance/diagnostics. Buffer-print helpers pad CPU to GPU stride for fair diff. This is the main end-to-end parity script. @tests/tUFF/test_UFF_multi.py#23-245
- **CPU vs PyOpenCL UFF**: [tests/tUFF/test_UFF_ocl.py](cci:7://file:///home/prokop/git/FireCore/tests/tUFF/test_UFF_ocl.py:0:0-0:0) compares C++ UFF (CPU) to PyOpenCL UFF (`UFF_CL`), including buffer dumps and an energy/force scan utility. Good for cross-validating the Python OpenCL stack. @tests/tUFF/test_UFF_ocl.py#34-240
- **MMFFsp3 PyOpenCL test**: `tests/tUFF/test_MMFFsp3_pyOCL.py` exists (not re-read here), but is separate from UFF; parity status there is not established in the shown files.

**What appears incomplete or still in-progress**
- **Energy parity in [test_UFF_multi.py](cci:7://file:///home/prokop/git/FireCore/tests/tUFF/test_UFF_multi.py:0:0-0:0)**: The GPU/CPU run path sets `energy = 0` placeholder; forces are compared, but energies are not currently downloaded/compared in that script. @tests/tUFF/test_UFF_multi.py#183-190
- **PyOpenCL UFF nonbond / grid parity**: `pyBall/OCL/UFF.py` has a `bDoNonBonded` flag and kernel arg scaffolding, but the shown run path only enqueues bonded kernels; nonbond/grid kernels aren’t wired in that class.
- **C++ GPU nonbond/grid kernels**: In `OCL_UFF.h::setup_kernels` the nonbond/grid setup calls are commented; flags exist but kernel binding/enqueue is gated out. @cpp/common/OpenCL/OCL_UFF.h#513-525
- **Documentation of achieved parity**: No explicit “parity achieved” note found; energy parity not asserted in tests.

**Relevant doc notes**
- [doc/DevNotes/UFF-GPU-vs-CPU-Debug-Guide.md](cci:7://file:///home/prokop/git/FireCore/doc/DevNotes/UFF-GPU-vs-CPU-Debug-Guide.md:0:0-0:0) documents the debug workflow, switches, and common pitfalls for UFF GPU vs CPU validation (arg ordering, offsets, energy buffers). @doc/DevNotes/UFF-GPU-vs-CPU-Debug-Guide.md#1-63
- Other UFF planning notes exist: [doc/DevNotes/UFF_muli_plan.md](cci:7://file:///home/prokop/git/FireCore/doc/DevNotes/UFF_muli_plan.md:0:0-0:0), `UFF_ocl_plan*.md`, and summaries under `doc/Topics/ForceFields/UFF_multiSystem*.md` (not re-opened here).

**Code pathways (summary)**
- **CPU MMFF/UFF**: `pyBall/MMFF_multi.py` → C FFI → `cpp/libs_OCL/MMFFmulti_lib.cpp::run()` → `W.ffu.run` / `run_omp` (UFF) or `W.run_*` (MMFF).
- **GPU UFF (C++ OpenCL)**: same entry, `iParalel=2` → `pack_uff_system/upload` → `eval_UFF_ocl()` which uses `OCL_UFF` kernels; forces downloaded and unpacked back to CPU buffers.
- **GPU UFF (PyOpenCL standalone)**: `pyBall/OCL/UFF.py` (UFF_CL) builds topology via `UFFbuilder.py`, allocates GPU buffers, uploads, then enqueues bonded kernels; used by [tests/tUFF/test_UFF_ocl.py](cci:7://file:///home/prokop/git/FireCore/tests/tUFF/test_UFF_ocl.py:0:0-0:0).
- **Utilities**: `pyBall/MD_test_utils.py` provides scan/plot helpers used by [test_UFF_ocl.py](cci:7://file:///home/prokop/git/FireCore/tests/tUFF/test_UFF_ocl.py:0:0-0:0) for force/energy scans (not parity-critical by itself).

**Parity status**
- Force parity is partially tested (CPU vs GPU) via [test_UFF_multi.py](cci:7://file:///home/prokop/git/FireCore/tests/tUFF/test_UFF_multi.py:0:0-0:0) (C++ path) and CPU vs PyOpenCL via [test_UFF_ocl.py](cci:7://file:///home/prokop/git/FireCore/tests/tUFF/test_UFF_ocl.py:0:0-0:0); energies are not validated in [test_UFF_multi.py](cci:7://file:///home/prokop/git/FireCore/tests/tUFF/test_UFF_multi.py:0:0-0:0).
- Nonbond/grid kernels parity is not exercised or bound in the current OpenCL UFF setup (both C++ and PyOpenCL).
- MMFF (non-UFF) parity on GPU is not covered by the shown tests; only UFF paths are compared.

If you want next steps, suggest: (1) add energy download/compare in [test_UFF_multi.py](cci:7://file:///home/prokop/git/FireCore/tests/tUFF/test_UFF_multi.py:0:0-0:0), (2) wire and bind nonbond/grid kernels in `OCL_UFF.h` and enable in tests, (3) run `test_MMFFsp3_pyOCL.py` to assess MMFF parity separately.

---

## USER

make as much as possible to complete the parity chekcs. Make MD instance, check it, and eventually add parity check for MMFF following what we do for UFF, we want to minimize code duplication, reuse the same functionality. We want to do small changes to code when possible (always analysing exising code before adding new or overriging it), but when it is necessary you can add new functionality but reusing exising as much as possible. Also run thest when woringin on it. This is Test driven development. Goal is to make sure bonding interaction for both MMFF and UFF in pyOpenCL implementations work (reproduce same forces, but this means we should forst check if the bonding topology and the parameters arrays are exactly the same). Make detailed plan and notes what we need to do and which files/function to analyse to achive that. At sime point we will ned to add debug print in both OpenCL kenrel files *.cl and in C++ in exactly the same format, so we can print them and analyse (eventually diff-them) to check parity 

make detailed notes and plan

---

## GPT-5.1-Codex


Here’s the concrete plan to finish bonded parity (UFF + MMFF) with minimal changes and test-driven steps:

## What’s missing
- UFF PyOpenCL parity: [test_UFF_ocl.py](cci:7://file:///home/prokop/git/FireCore/tests/tUFF/test_UFF_ocl.py:0:0-0:0) still crashes because [scan_energy_force](cci:1://file:///home/prokop/git/FireCore/pyBall/MD_test_utils.py:686:0-730:5) needs an `md` (MolecularDynamics) instance.
- MMFF parity: no CPU↔GPU harness; only the MD diagnostics `test_MMFFsp3_pyOCL.py`.
- Topology/parameter parity verification (bonded) is not explicitly diffed across CPU/GPU for MMFF (and only implicitly for UFF).

## Minimal changes to implement
1) **Provide `md` in [test_UFF_ocl.py](cci:7://file:///home/prokop/git/FireCore/tests/tUFF/test_UFF_ocl.py:0:0-0:0)** for the scan:
   - Instantiate [MolecularDynamics](cci:2://file:///home/prokop/git/FireCore/pyBall/OCL/MolecularDynamics.py:36:0-640:21) from [pyBall/OCL/MolecularDynamics.py](cci:7://file:///home/prokop/git/FireCore/pyBall/OCL/MolecularDynamics.py:0:0-0:0).
   - Allocate buffers with the same UFF template (natoms/nvecs).
   - Pass `md` to [scan_energy_force(mm, md, …)](cci:1://file:///home/prokop/git/FireCore/pyBall/MD_test_utils.py:686:0-730:5) or disable the scan and do direct CPU-vs-GPU force diff (simpler).
   - Run: `cd tests/tUFF && python3 test_UFF_ocl.py --components bonds,angles,dihedrals,inversions --tolerance 1e-3`.

2) **Add MMFF bonded parity harness** reusing UFF logic to avoid duplication:
   - Option A: extend [tests/tUFF/test_UFF_multi.py](cci:7://file:///home/prokop/git/FireCore/tests/tUFF/test_UFF_multi.py:0:0-0:0) with `--ff {uff,mmff}` to reuse scan/diff logic:
     - `bUFF=False, bMMFF=True`, disable NonBonded/GridFF, run CPU (`iParalel=0`) vs GPU (`iParalel=2`), compare `fapos`.
   - Option B: a new lightweight script calling `pyBall.MMFF_multi` similarly.
   - Include topology/param dumps (bond/angle/dihedral/inversion arrays) for both CPU/GPU to verify identity before force diff.

3) **Debug/trace hooks for topology/params parity (bonded only):**
   - Host side: print/dump packed arrays just before upload:
     - UFF: [pack_uff_system](cci:1://file:///home/prokop/git/FireCore/cpp/common/molecular/MolWorld_sp3_multi.h:801:0-870:1) host slices (`REQs`, `bonAtoms/Params`, `angAtoms/Params`, `dihParams`, `invParams`, `a2f`).
     - MMFF: equivalent bonded arrays from `MMFF_multi` / [MolecularDynamics.pack_system](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/MolecularDynamics.py:169:4-207:122).
   - Kernel side (if needed): guarded prints in [UFF.cl](cci:7://file:///home/prokop/git/FireCore/cpp/common_resources/cl/UFF.cl:0:0-0:0) and `relax_multi_mini.cl` with per-system `IDBG_SYS` and a single work-item to dump offsets and a sample force accumulation path.

## Files/functions to touch
- [tests/tUFF/test_UFF_ocl.py](cci:7://file:///home/prokop/git/FireCore/tests/tUFF/test_UFF_ocl.py:0:0-0:0): add MD instantiation and pass `md` to [scan_energy_force](cci:1://file:///home/prokop/git/FireCore/pyBall/MD_test_utils.py:686:0-730:5), or skip scan and do force diff only.
- [tests/tUFF/test_UFF_multi.py](cci:7://file:///home/prokop/git/FireCore/tests/tUFF/test_UFF_multi.py:0:0-0:0) (or a new small script): add `--ff mmff` branch to run bonded-only MMFF CPU vs GPU using `pyBall.MMFF_multi` with `iParalel=0/2`.
- Optional debug prints:
  - Host: [MolWorld_sp3_multi::pack_uff_system](cci:1://file:///home/prokop/git/FireCore/cpp/common/molecular/MolWorld_sp3_multi.h:801:0-870:1) / [pack_system](cci:1://file:///home/prokop/git/FireCore/cpp/common/molecular/MolWorld_sp3_multi.h:588:0-676:1) for bonded params.
  - Kernels: [UFF.cl](cci:7://file:///home/prokop/git/FireCore/cpp/common_resources/cl/UFF.cl:0:0-0:0), `relax_multi_mini.cl` (bonded kernels) with `IDBG_SYS` gate.

## Test commands (bonded only, correct CWD)
- UFF C++ parity (already good): `cd tests/tUFF && python3 test_UFF_multi.py`
- UFF PyOCL (after md fix or scan skip): `cd tests/tUFF && python3 test_UFF_ocl.py --components bonds,angles,dihedrals,inversions --tolerance 1e-3`
- MMFF parity (after harness added): e.g., `cd tests/tUFF && python3 test_UFF_multi.py --ff mmff --tol 1e-5` (or dedicated script)

## Current status
- UFF C++ CPU vs GPU bonded: forces match (~5e-6), teardown double-free persists.
- UFF PyOCL: blocked by missing `md` in scan.
- MMFF: no parity test yet.

I can implement the `md` wiring for UFF PyOCL and add the MMFF parity branch with shared logic next—minimal edits, no duplication.


---

## USER

if scan_energy_force is working on MMFF rathet then UFF then make ne version from UFF still do int in @MD_test_utils.py so we consolidate the testing utilities. 

Notice that we are able to run both UFF and MMFF over C++ thise are the drivers
@OCL_MM.h @OCL_UFF.h and it should be possible to initialize both in @MMFF_multi.py @MMFFmulti_lib.cpp and also get relevant buffers see
init_buffers
init_buffers_UFF
MMFFmulti_lib::init() 

note common mistake it so thing that bMMFF should be false to run UFF that is not tru bMMFF must be always tru, wna we swich the two only by bUFF (MMFF is default when bUFF=false)

note that @OCL_MM.h is interfacing relax_multi.cl similarly as @MolecularDynamics.py 
and in the same way
@UFF.cl is interfaced by both @OCL_UFF.h and @UFF.py ,
that is important information for meating parity. If we prepare same inputs (topology, params) the two must give same results as thay use same code (same *.cl files) but crucial is to provide same impots (topology and params), for that is esential to make sure @UFFbuilder.py and @MMFF.py works properly and reporuce @MMFFBuilder.h 

we have also purely CPU / C++ based implementation not using OpenCL which is not paralelized over replicas of system, that is in 
@MMFFsp3_loc.h for MMFF and in @UFF.h for UFF


analyse the referenced files note importaint functions and variables, and make clar notes, this is very importaint insight you should alsways keep in mind wrilte debugging


---

## GPT-5.1-Codex

Below are consolidated, file-by-file notes focusing on how topology/params/positions are prepared and fed into C++ vs PyOpenCL for UFF/MMFF parity. I’ve highlighted key functions, required flags (bMMFF/bUFF), buffer layouts, and where to verify equality.

## cpp/common/OpenCL/OCL_UFF.h
- Wraps UFF.cl for C++ OpenCL path (used by MMFFmulti_lib when bUFF=True).
- Buffers: apos (float4), fapos (float4), REQs (float4), bonAtoms (int2), bonParams (float2), angAtoms (int4), angNgs (int2), angParams (float4 + K), dihAtoms (int4), dihNgs (int4), dihParams (float4), invAtoms (int4), invNgs (int3/4), invParams (float4), neighs (int4), neighBs (int4), a2f (offsets/counts/indices), fint (force pieces), hneigh.
- Key: **nf_per_system** = 2*nbonds + 3*nangles + 4*ndihedrals + 4*ninversions; fint sized by this, not natoms.
- Kernel args set in [setup_kernels()](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/MolecularDynamics.py:268:4-290:79): evalBondsAndHNeigh_UFF, evalAngles_UFF, evalDihedrals_UFF, evalInversions_UFF, assembleForces_UFF. Angles/dihedrals/inversions write into fint then assemble -> fapos.
- Uses 2D NDRange (global.y = nSystems) offsets: i0a, i0h, i0f, i0a2f, i0bon/i0ang/i0dih/i0inv.

## cpp/common/OpenCL/OCL_MM.h
- Interface to relax_multi_mini.cl for MMFF MD path (C++ OpenCL). Mirrors MolecularDynamics.py role.
- Buffers: apos/fapos, neighs/bs/cell, params, MDparams, etc.; updates via updateAtomsMMFFf4 kernels. Uses same kernels as PyOpenCL MD.

## cpp/common_resources/cl/UFF.cl
- UFF kernels shared by C++ OCL and PyOpenCL: clear_fapos_UFF, clear_fint_UFF, evalBondsAndHNeigh_UFF, evalAngles_UFF, evalDihedrals_UFF, evalInversions_UFF, assembleForces_UFF (+ NB kernels).
- Parameter order requires scalars: natoms, npbc, i0bon/i0ang/i0dih/i0inv, nf_per_system, Rdamp, FmaxNonBonded, SubNBTorsionFactor, bSubtractBond/Angle NB, etc.
- Expects hneigh sized natoms*4 float4 per system; fint sized nf_per_system per system.
- Uses a2f map to assemble angle/dih/inv forces.

## cpp/common/molecular/UFF.h (CPU)
- Reference UFF implementation (serial/OpenMP).
- Buffers: apos, fapos, REQs, bonAtoms/Params, angAtoms/Params, dihAtoms/Params, invAtoms/Params, neighs/bs, hneigh, fint, a2f.
- Flow: bonds fill fapos + hneigh; angles/dihedrals/inversions fill fint; assembleAtomForce uses a2f to accumulate into fapos.
- Flags: bUFF must be true; bMMFF must stay true (base class) even for UFF.

## cpp/common/molecular/MMFFsp3_loc.h (CPU MMFF)
- CPU MMFF reference; bonded + pi terms. Topology/params assembled by MMFFBuilder.

## cpp/common/molecular/MMFFBuilder.h (CPU builder)
- Builds bonded/nonbonded topology for MMFF (and UFF when bUFF on). Source of truth for how bon/ang/dih/inv arrays and params are formed.

## cpp/libs_OCL/MMFFmulti_lib.cpp
- Exposes C++ multi-system MMFF/UFF to Python (pyBall/MMFF_multi.py).
- init(): builds MolWorld_sp3_multi; bMMFF must be true; bUFF toggles UFF kernels.
- init_buffers_UFF(): exposes UFF-specific buffers to Python.
- run(): dispatch iParalel 0/1 CPU, 2 GPU OCL (UFF via OCL_UFF, MMFF via OCL_MM).
- scan(): UFF-only batched GPU vs CPU forces.

## pyBall/MMFF_multi.py
- Python ctypes wrapper to MMFFmulti_lib. Flags: bMMFF must be True always; bUFF toggles UFF.
- setSwitchesUFF controls bonded components and assemble flag; setSwitches2 controls MMFF/NonBonded/PBC/Grid.
- scan() for UFF GPU vs CPU.

## pyBall/OCL/UFFbuilder.py
- Builds UFF topology/params for PyOpenCL: bonAtoms/Params, angAtoms/Params, dih/inv, REQs, neighs/bs, a2f map.

## pyBall/OCL/UFF.py (PyOpenCL UFF)
- Loads UFF.cl. [toUFF(mol)](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/UFF.py:50:4-67:23) uses UFF_Builder, allocates buffers via realloc_buffers (natoms, nbonds, nangles, ndihedrals, ninversions, nf_per_system, hneigh sized natoms*4 float4, fint sized nf_per_system).
- upload_topology_params uploads all bonded arrays plus a2f_offsets/counts/indices, neighs/bs, REQs, atype.
- kernel_params include natoms… npbc, n, nf_per_system, i0bon/i0ang/i0dih/i0inv, Rdamp, Fmax, SubNBTorsionFactor, bSubtractBond/Angle NB.
- run_eval_step: clears fapos; clears fint if any bonded terms; runs bonds->angles->dihedrals->inversions; assembleForces_UFF; optional NonBonded.

## pyBall/OCL/MMFF.py (PyOpenCL MMFF)
- Builder for MMFF (pi orbitals, etc.), prepares topology/params for relax_multi_mini.cl; used with MolecularDynamics.py.

## pyBall/OCL/MolecularDynamics.py
- PyOpenCL MD driver for relax_multi_mini.cl (MMFF path). Manages apos/fapos, neighs/params, MDparams, runs getMMFFf4 / updateAtoms kernels. Not UFF.

## pyBall/MD_test_utils.py
- MMFF utilities: configure_md, run_step, evaluate_mmff_gpu, scan_energy_force(mm, md, …). This pipeline calls relax_multi_mini.cl and requires both mm and md (MMFF only). Not suitable for UFF kernels.

## tests/tUFF/test_UFF_ocl.py
- Intended for UFF CPU vs PyOpenCL UFF parity (topology/params/forces). Currently broken if calling scan_energy_force (MMFF-only). Should instead compare CPU UFF vs UFF_CL (buffers + fapos).

## Key parity rules / common pitfalls
- bMMFF must remain True; toggle UFF via bUFF to enter UFF paths.
- Match topology/params exactly: bonAtoms/Params, angAtoms/Params, dihParams, invParams, neighs/bs, REQs, a2f (offsets/counts/indices). Verify shapes/strides match C++ expectations.
- Buffers: fint sized by nf_per_system; hneigh sized natoms*4 float4 per system; a2f map per atom; use 2D NDRange offsets i0a/i0f/i0a2f.
- UFF PyOpenCL parity: do not use scan_energy_force (MMFF). Use UFF_CL.run_eval_step then compare buffers/forces vs CPU UFF.
- MMFF PyOpenCL parity: use scan_energy_force with mm+md (relax_multi_mini.cl), compare to C++ MMFF_multi CPU forces.

I have read and extracted the critical pieces from the referenced files; will proceed to (a) fix the UFF test to do buffer+force parity, and (b) add a UFF-specific scan helper in MD_test_utils, then (c) add MMFF parity harness.