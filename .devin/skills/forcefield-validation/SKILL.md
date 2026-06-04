---
name: forcefield-validation
description: Use when implementing or debugging interatomic force-fields (e.g., UFF, MMFF)
trigger:
  glob:
    - "**/MMFF*"
    - "**/UFF*"
    - "**/common_resources/**/*.dat"
    - "**/common_resources/**/*.xyz"
    - "**/common_resources/**/*.mol"
---

## Architecture

- Reference: `pyBall/MMFF_multi.py` (wrapper for C++ `MMFFmulti_lib.cpp`)
- Backends: `iParalel=0` (C++ CPU), `iParalel=2` (C++ OpenCL), `pyBall/OCL/UFF.py` (Python OpenCL)
- Data files: `AtomTypes.dat`, `BondTypes.dat`, `AngleTypes.dat` in `common_resources/`

## Switch Semantics

0=keep, >0=force-on, <0=force-off. Passing 0 does nothing.
- UFF: `setSwitchesUFF(DoBond, DoAngle, DoDihedral, DoInversion, DoAssemble, ...)`. `bSubtractBondNonBond` must align with `setSwitchesUFF_NB`.
- MMFF: `setSwitches2(..., MMFF, Angles, PiSigma, PiPiI)`. `bMMFF` is base class flag (True even for UFF). Check `bUFF` to distinguish.

## Buffer Parity

- Call `init_buffers(bUFF=...)` to populate C++ `buffers`/`ibuffers` maps. Read `ndims` first for shapes.
- UFF critical: Check `bonParams`, `angParams`, `a2f`. CPU `angParams` layout is `[K, c0, c1, c2, c3]`; kernels expect split `angParams1=[c0..c3]` and `angParams2_w=[K]`.
- MMFF critical: Check `apos` (Pi-orbitals at `natoms:natoms+nnode`), `bkNeighs` (Back Neighbors), `Ksp`/`Kpp`.

## Pitfalls

- Pi-orbitals: In MMFF, `apos` contains atoms AND pi-nodes. Loop bounds: `natoms` vs `natoms+nnode`.
- Node/cap layout: Current C++ builder sets `nnode=natoms`, allocates one pi slot per atom. Caps have `Ksp/Kpp=0` but occupy `nvecs`. `bkNeighs` sized `nSystems*nvecs`.

## Bonded-Only Parity Flow

`cleanF → getMMFFf4 → updateAtomsMMFFf4(dt=0)` to assemble recoil from `fneigh` via `bkNeighs` without moving atoms. Set `bSubtractVdW=0` when NonBonded off. Propagate switches to `ffls[isys]`. Copy `fapos` back to shared buffers.