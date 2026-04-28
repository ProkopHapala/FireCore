
# UFF Throughput Tutorial ([tests/tUFF/run_throughput_UFF.py](cci:7://file:///home/prokophapala/git/prokop_and_master/FireCore/tests/tUFF/run_throughput_UFF.py:0:0-0:0))

## Overview
This guide explains how to launch the UFF throughput molecular dynamics driver in [tests/tUFF/run_throughput_UFF.py](cci:7://file:///home/prokophapala/git/prokop_and_master/FireCore/tests/tUFF/run_throughput_UFF.py:0:0-0:0), highlights the Python ⇄ C++ ⇄ OpenCL stack it exercises, and outlines how to enable the GridFF accelerated non-bonded path. Use it as both user-facing instructions and a developer-facing reference to the relevant code.

## Prerequisites
- **Python** CPython ≥3.8 with NumPy available (only stdlib + NumPy are required by the script).
- **Compiled library** `libMMFFmulti_lib.so` produced under `cpp/Build/libs_OCL/` (the script loads it via `pyBall/cpp_utils_.py`).
- **GPU runtime** Valid OpenCL runtime, device drivers, and permissions for the target GPU.
- **UFF resources** Files in `tests/tUFF/common_resources/` (`ElementTypes.dat`, `AtomTypes.dat`, etc.) plus molecular `.xyz` inputs under `tests/tUFF/data_UFF/`.
- **Environment** Run from repository root or ensure working directory contains the data folders referenced above.

## Key Files and Interfaces
- **[tests/tUFF/run_throughput_UFF.py](cci:7://file:///home/prokophapala/git/prokop_and_master/FireCore/tests/tUFF/run_throughput_UFF.py:0:0-0:0)** main CLI driving UFF throughput loops.
- **[pyBall/MMFF_multi.py](cci:7://file:///home/prokophapala/git/prokop_and_master/FireCore/pyBall/MMFF_multi.py:0:0-0:0)** Python ctypes bindings; exposes [uff.init()](cci:1://file:///home/prokophapala/git/prokop_and_master/FireCore/cpp/libs_OCL/MMFFmulti_lib.cpp:145:0-187:1), [uff.MDloop()](cci:1://file:///home/prokophapala/git/prokop_and_master/FireCore/cpp/libs_OCL/MMFFmulti_lib.cpp:211:0-216:1), [uff.setSwitches2()](cci:1://file:///home/prokophapala/git/prokop_and_master/FireCore/pyBall/MMFF_multi.py:498:0-499:124), [uff.setSwitchesUFF()](cci:1://file:///home/prokophapala/git/prokop_and_master/FireCore/cpp/libs_OCL/MMFFmulti_lib.cpp:218:0-267:1), etc.
- **[cpp/libs_OCL/MMFFmulti_lib.cpp](cci:7://file:///home/prokophapala/git/prokop_and_master/FireCore/cpp/libs_OCL/MMFFmulti_lib.cpp:0:0-0:0)** C-entry points ([init](cci:1://file:///home/prokophapala/git/prokop_and_master/FireCore/cpp/libs_OCL/MMFFmulti_lib.cpp:145:0-187:1), [MDloop](cci:1://file:///home/prokophapala/git/prokop_and_master/FireCore/cpp/libs_OCL/MMFFmulti_lib.cpp:211:0-216:1), [setSwitches2](cci:1://file:///home/prokophapala/git/prokop_and_master/FireCore/pyBall/MMFF_multi.py:498:0-499:124), [setSwitchesUFF](cci:1://file:///home/prokophapala/git/prokop_and_master/FireCore/cpp/libs_OCL/MMFFmulti_lib.cpp:218:0-267:1)) gluing to the simulation world.
- **[cpp/common/molecular/MolWorld_sp3_multi.h](cci:7://file:///home/prokophapala/git/prokop_and_master/FireCore/cpp/common/molecular/MolWorld_sp3_multi.h:0:0-0:0)** host-side engine orchestrating replicas, OpenCL upload/download, and GridFF integration.
- **[cpp/common/OpenCL/OCL_UFF.h](cci:7://file:///home/prokophapala/git/prokop_and_master/FireCore/cpp/common/OpenCL/OCL_UFF.h:0:0-0:0)** GPU driver managing OpenCL buffers, task setup, and kernel dispatch.
- **`tests/tDFT/data/cl/UFF.cl`** OpenCL kernels for bonds/angles/dihedrals/inversions, GridFF evaluation, and non-bonded interactions.

## Workflow Overview
1. **CLI parsing:** [run_throughput_UFF.py](cci:7://file:///home/prokophapala/git/prokop_and_master/FireCore/tests/tUFF/run_throughput_UFF.py:0:0-0:0) collects flags such as `--nSys`, `--loops`, and GridFF toggles.
2. **Python init:** [uff.init()](cci:1://file:///home/prokophapala/git/prokop_and_master/FireCore/cpp/libs_OCL/MMFFmulti_lib.cpp:145:0-187:1) ([MMFF_multi.py](cci:7://file:///home/prokophapala/git/prokop_and_master/FireCore/pyBall/MMFF_multi.py:0:0-0:0)) forwards parameters (UFF flags, thermostat values, GridFF settings) to [MMFFmulti_lib::init()](cci:1://file:///home/prokophapala/git/prokop_and_master/FireCore/cpp/libs_OCL/MMFFmulti_lib.cpp:145:0-187:1).
3. **World setup:** [MolWorld_sp3_multi::init()](cci:1://file:///home/prokophapala/git/prokop_and_master/FireCore/cpp/libs_OCL/MMFFmulti_lib.cpp:145:0-187:1) loads structures, prepares UFF force fields, allocates multi-system replicas, and, when requested, schedules GridFF data.
4. **OpenCL preparation:** [OCL_UFF::realloc()](cci:1://file:///home/prokophapala/git/prokop_and_master/FireCore/cpp/common/OpenCL/OCL_UFF.h:232:4-310:5) and [OCL_UFF::setup_kernels()](cci:1://file:///home/prokophapala/git/prokop_and_master/FireCore/cpp/common/OpenCL/OCL_UFF.h:312:4-503:5) lay out device buffers and kernel arguments, referencing `tests/tDFT/data/cl/UFF.cl`.
5. **MD throughput loop:** The Python script repeatedly invokes [uff.MDloop()](cci:1://file:///home/prokophapala/git/prokop_and_master/FireCore/cpp/libs_OCL/MMFFmulti_lib.cpp:211:0-216:1) which calls [MMFFmulti_lib::MDloop()](cci:1://file:///home/prokophapala/git/prokop_and_master/FireCore/cpp/libs_OCL/MMFFmulti_lib.cpp:211:0-216:1); this drives [MolWorld_sp3_multi::MDloop()](cci:1://file:///home/prokophapala/git/prokop_and_master/FireCore/cpp/libs_OCL/MMFFmulti_lib.cpp:211:0-216:1) to run GPU kernels for force evaluation and integration.
6. **Trajectory output:** [uff.setTrjName()](cci:1://file:///home/prokophapala/git/prokop_and_master/FireCore/pyBall/MMFF_multi.py:534:0-538:53) keeps per-system `.xyz` trajectories (`traj_UFF_###.xyz`), refreshed per 1000 steps by default.

## Running the Script
```bash
python tests/tUFF/run_throughput_UFF.py  --nSys 512 --xyz_name data_UFF/xyz/HHO-h.xyz --loops 500 --perframe 1000 --perVF 100 --Fconv 1e-4 --bNonBonded 1 --bGridFF -1
```
- **`--nSys`** number of replicas packed for GPU evaluation.
- **`--xyz_name`** base name (without `.xyz`) in `tests/tUFF/`.
- **`--perframe`** MD steps between host syncs.
- **`--perVF`** inner kernel repetition count (`MolWorld_sp3_multi::nPerVFs`).
- **`--loops`** total [MDloop()](cci:1://file:///home/prokophapala/git/prokop_and_master/FireCore/cpp/libs_OCL/MMFFmulti_lib.cpp:211:0-216:1) invocations; progress logs every 50 iterations.
- **`--gridnPBC`** tuple string passed through to [OCL_UFF::realloc()](cci:1://file:///home/prokophapala/git/prokop_and_master/FireCore/cpp/common/OpenCL/OCL_UFF.h:232:4-310:5) for periodic images in GridFF mode.
- **`--bSaveToDatabase`** placeholder for database logging; `-1` keeps default off.

## Enabling GridFF
- **dovdW**: set `--dovdW 1` to enable non-bonded interactions (currently it is important only for labeling the simulation)
- **doSurfaceAtoms**: set `--doSurfAtoms 1` if surface should be computed
- **CLI flag**: set `--bGridFF 6` number 6 must be to choose bSpline mode, if negative GridFF is disabled
- **Periodic images**: adjust `--gridnPBC "(1,1,0)"` it is for adjusting periodic surfaces imporatnt for surfaces calculated with full atom
- **Surface data**: provide `--surf_name data_UFF/surf/YourSurface` if coupling molecules to surfaces.
- **Kernel path**: GridFF kernels are compiled from [OCL_UFF::makeKernels("common_resources/cl")](cci:1://file:///home/prokophapala/git/prokop_and_master/FireCore/cpp/common/OpenCL/OCL_UFF.h:207:4-230:5), so ensure matching grid assets under `tests/tUFF/common_resources/`.

Example GridFF run:
```bash
python tests/tUFF/run_throughput_UFF.py \
    --nSys 256 \
    --xyz_name data_UFF/xyz/benzene \
    --surf_name data_UFF/surf/graphene \
    --bGridFF 1 \
    --gridnPBC "(2,2,1)" \
    --T 1200 --gamma 0.05 --nExplore 500 --nRelax 80000
```

## Developer Reference

- **[uff.init()](cci:1://file:///home/prokophapala/git/prokop_and_master/FireCore/cpp/libs_OCL/MMFFmulti_lib.cpp:145:0-187:1)** (`pyBall/MMFF_multi.py:455`) marshals Python types to C signatures, sets global `nSys`, and keeps track of UFF toggles.
- **[uff.setSwitches2()](cci:1://file:///home/prokophapala/git/prokop_and_master/FireCore/pyBall/MMFF_multi.py:498:0-499:124)** (`pyBall/MMFF_multi.py:499`) toggles world-level flags such as `NonBonded`, `SurfAtoms`, `GridFF`; C++ handler [MMFFmulti_lib::setSwitches2()](cci:1://file:///home/prokophapala/git/prokop_and_master/FireCore/pyBall/MMFF_multi.py:498:0-499:124) updates [MolWorld_sp3_multi](cci:2://file:///home/prokophapala/git/prokop_and_master/FireCore/cpp/common/molecular/MolWorld_sp3_multi.h:106:0-3595:1) state.
- **[uff.setSwitchesUFF()](cci:1://file:///home/prokophapala/git/prokop_and_master/FireCore/cpp/libs_OCL/MMFFmulti_lib.cpp:218:0-267:1)** (`pyBall/MMFF_multi.py:505`) controls component-level masks and non-bond subtraction/clamping; [MMFFmulti_lib::setSwitchesUFF()](cci:1://file:///home/prokophapala/git/prokop_and_master/FireCore/cpp/libs_OCL/MMFFmulti_lib.cpp:218:0-267:1) synchronizes CPU and GPU flags and rebinds kernels.
- **[uff.MDloop()](cci:1://file:///home/prokophapala/git/prokop_and_master/FireCore/cpp/libs_OCL/MMFFmulti_lib.cpp:211:0-216:1)** (`pyBall/MMFF_multi.py:576`) wraps [lib.MDloop](cci:1://file:///home/prokophapala/git/prokop_and_master/FireCore/cpp/libs_OCL/MMFFmulti_lib.cpp:211:0-216:1); C++ implementation sets kernel iteration counts and calls [MolWorld_sp3_multi::MDloop()](cci:1://file:///home/prokophapala/git/prokop_and_master/FireCore/cpp/libs_OCL/MMFFmulti_lib.cpp:211:0-216:1).
- **Trajectory handling**: [uff.setTrjName()](cci:1://file:///home/prokophapala/git/prokop_and_master/FireCore/pyBall/MMFF_multi.py:534:0-538:53) (`pyBall/MMFF_multi.py:535`) and [MolWorld_sp3_multi](cci:2://file:///home/prokophapala/git/prokop_and_master/FireCore/cpp/common/molecular/MolWorld_sp3_multi.h:106:0-3595:1) cooperate to maintain `traj_UFF_###.xyz`.
- **Kernel tasks**: [OCL_UFF::setup_kernels()](cci:1://file:///home/prokophapala/git/prokop_and_master/FireCore/cpp/common/OpenCL/OCL_UFF.h:312:4-503:5) (`OCL_UFF.h:313`) defines workgroup sizes and binds arguments for `evalBondsAndHNeigh_UFF`, `evalAngles_UFF`, `getNonBond`, etc.
- **Multi-system layout**: [MolWorld_sp3_multi::realloc()](cci:1://file:///home/prokophapala/git/prokop_and_master/FireCore/cpp/common/OpenCL/OCL_UFF.h:232:4-310:5) (`MolWorld_sp3_multi.h:254`) sizes host and device buffers, including GridFF-specific arrays (`host_neighs_UFF`, `host_a2f_indices`).

## Troubleshooting & Tips
- **Library mismatch**: If Python raises `OSError` while loading `libMMFFmulti_lib.so`, rebuild via the project’s CMake/Bash flow producing `cpp/Build/libs_OCL/libMMFFmulti_lib.so`.
- **Kernel rebuild**: After editing `tests/tDFT/data/cl/UFF.cl`, clear cached binaries so [OCL_UFF::makeKernels()](cci:1://file:///home/prokophapala/git/prokop_and_master/FireCore/cpp/common/OpenCL/OCL_UFF.h:207:4-230:5) compiles fresh code.
- **Verbose diagnostics**: Temporarily call [uff.setVerbosity(verbosity=2, idebug=1)](cci:1://file:///home/prokophapala/git/prokop_and_master/FireCore/pyBall/MMFF_multi.py:420:0-421:48) before [uff.init()](cci:1://file:///home/prokophapala/git/prokop_and_master/FireCore/cpp/libs_OCL/MMFFmulti_lib.cpp:145:0-187:1) for detailed world/kernels logging.
- **Buffer inspection**: Use [uff.printBuffNames()](cci:1://file:///home/prokophapala/git/prokop_and_master/FireCore/pyBall/MMFF_multi.py:181:0-182:24) and [uff.getBuff()](cci:1://file:///home/prokophapala/git/prokop_and_master/FireCore/pyBall/MMFF_multi.py:196:0-200:48) helpers ([pyBall/MMFF_multi.py](cci:7://file:///home/prokophapala/git/prokop_and_master/FireCore/pyBall/MMFF_multi.py:0:0-0:0)) to introspect GPU-host buffer mappings during debugging.

# Task Status
- **Status** Completed: generated full Markdown tutorial covering usage, GridFF workflow, and developer references for the specified files.