## USER

I would like to make code map for my C++ Classical Forcefield system, it contains two forcefields UFF (universal forcefield) and MMFFsp3 (which is not Merk molecular forcefield but my own forcefield)

CPU C++ implementations of forcefields itself are in 
@MMFFsp3_loc.h for MMFFsp3
@UFF.h  for UFF
@NBFF.h for non covalent interactions

all of these calsses are derived from common parent @ForceField.h which is derived from @Atoms.h 
so they share common datamember and functions

these are integrated into molecular simulation world @MolWorld_sp3.h and can be imported into python library over
@MMFF_lib.cpp @MMFF.py 

To build molecular topology and assing types @MolWorld_sp3.h uses @MMFFBuilder.h @MMFFparams.h 

beside that there is also GPU OpenCL version of these forcefields, the kernels are implemented in
@relax_multi_mini.cl for MMFFsp3 and NBFF, and in @UFF.cl fr UFF, there are C++ interfaces
@OCL_MM.h @OCL_UFF.h the whole think is than integrated to @MolWorld_sp3_multi.h which is derived from 
@MolWorld_sp3.h , and the  main difference that @MolWorld_sp3_multi.h is optimized for simulation of many replicas of the system in parallel on GPU (this is because our systems are rather small 100-1000 atoms, and GPU has >10k threads, so we need to have up to thousands replicas to saturate GPU threads)

we also want to have interface to run GPU version directly from python using pyOpenCL withut the C++ intermediate layer. We have this here but maybe it is not completely finished
@MolecularDynamics.py @MMparams.py @MMFF.py @UFF.py @UFFbuilder.py 

molecular geometry is read from .xyz, .mol, .mol2 using the Builder class

Parameters itself are stored in this in @AtomTypes.dat @AngleTypes.dat @BondTypes.dat @DihedralTypes.dat @ElementTypes.dat 

now your job is read all these files, not most importaint functions, and note relations between the files.
Good would be make a table of modules and table of functions which are kinda equvalent between the different variants, I mean we have UFF, MMFFsp3, we have C++ CPU implementation, then OpenCL GPU implementation with C++ driver, and finally pyOpenCL GPU implemetation with python driver. It maybe be hard to navigate. 

What we often do that we are doing pariti testing if one version reporduce results of the other version, for that it is importaint to follow the execution path from reading geometry and forecfild parameters, preparing the topology (type assingment), preparing axuiliary buffers (neighbor lists, PBC mapping ...) then eventually preparing OpenCL buffers, and running OpenCL kernels in proper sequence downloading results every few cycles and chekcing convergence (we need to clean some buffers every cycle, adjust some parameters etc. ... )

Plese make comprehensive notes analysing and listing all this in systematic manner cmparing the differentc versions

---

## GPT-5.1-Codex

Here’s a concise code-map across CPU, GPU (C++ OpenCL), and pyOpenCL for the classical forcefields (MMFFsp3, UFF, NBFF), focusing on flow from geometry/params to execution. I drew on the retrieved code (notably MMFFsp3_loc, MMFFsp3, NBFF, and the pyOpenCL builder path) to outline equivalences and data paths.

## Module inventory (by layer)

**CPU C++ core**
- Atoms/base: `cpp/common/molecular/Atoms.h`, parent `ForceField.h` (shared members, PBC flags, REQs, etc.).
- Nonbonded: `NBFF.h` — Coulomb/LJ, PBC shifts, cell/AABB helpers @`NBFF.h#24-93`.
- Bonded (localized): `MMFFsp3_loc.h` — per-atom localized bonds/angles/pi, neighbor forces in `fneigh`, optional nonbond subtraction, torsion handled via pi terms @`MMFFsp3_loc.h#29-574`.
- Bonded (bond-list): `MMFFsp3.h` — bond-list-based variant with `bond2atom`, `bond_l0`, `bond_k`, `NeighParams`, etc.; uses `realloc` to size DOFs, bond arrays, neighbors @`MMFFsp3.h#15-149`.
- UFF: `UFF.h` — UFF bonded + NB terms (not opened here, parallel role to MMFFsp3_loc).
- Builders/params: `MMFFBuilder.h`, `MMFFparams.h` (topology/type assignment, reading *.dat param tables).
- World integration: [MolWorld_sp3.h](cci:7://file:///home/prokop/git/FireCore/cpp/common/molecular/MolWorld_sp3.h:0:0-0:0) (single system, CPU) ties Atoms, FFs, builder, neighbor lists.

**GPU OpenCL (C++ driver + kernels)**
- Kernels: `cpp/common_resources/cl/relax_multi_mini.cl` (MMFFsp3_loc + NBFF), `cpp/common_resources/cl/UFF.cl` (UFF bonded+NB).
- C++ drivers: `cpp/common/OpenCL/OCL_MM.h` (MMFFsp3/NB), `cpp/common/OpenCL/OCL_UFF.h` (UFF).
- Multi-replica world: `cpp/common/molecular/MolWorld_sp3_multi.h` (inherits CPU world, orchestrates buffer uploads per replica, kernel sequencing).

**pyOpenCL path (Python driver + kernels above)**
- Driver/topology: `pyBall/OCL/MMFF.py` (MMFFsp3_loc structure, host packing), `pyBall/OCL/UFF.py`, `pyBall/OCL/UFFbuilder.py`.
- Params helper: `pyBall/OCL/MMparams.py`.
- MD loop: `pyBall/OCL/MolecularDynamics.py` (kernel sequencing, buffer management).
- Uses same kernels (`relax_multi_mini.cl`, `UFF.cl`) via pyOpenCL.

**Parameter data**
- `cpp/common_resources/AtomTypes.dat`, `AngleTypes.dat`, `BondTypes.dat`, `DihedralTypes.dat`, `ElementTypes.dat` — read by builders/params on both CPU and pyOpenCL paths.

## Execution path (common logical phases)

1. **Geometry import**: Builder reads .xyz/.mol/.mol2, populates atoms (`apos`, `atypes`, `enames`), bonds list, neighbor lists.
2. **Parameter/type assignment**:
   - CPU: `MMFFBuilder` + `MMFFparams` map elements to Atom/Bond/Angle/Dihedral types (from *.dat), fill per-atom/bond params.
   - pyOpenCL: `UFFbuilder.py`/`MMFF.py.toMMFFsp3_loc()` ingest AtomicSystem and atom_types dict (from MMparams) to fill arrays.
3. **Topology/buffers**:
   - MMFFsp3_loc (CPU): `realloc(nnode,ncap,ntors)` allocates `apos/fapos`, `pipos/fpipos`, `neighs`, `bkneighs`, `bLs/bKs/Ksp/Kpp`, `apars`, `angles`, torsion data, constraints @`MMFFsp3_loc.h#133-171`.
   - MMFFsp3 (CPU bond-list): `realloc(nnode,nbonds,npi,ncap)` allocates bond list arrays `bond2atom`, `bond_l0/k/kPi`, neighbor arrays if enabled @`MMFFsp3.h#76-109`.
   - NBFF: allocates PBC shifts (`makePBCshifts_`), buckets/AABB for cutoffs.
   - GPU C++ driver: allocates mirrored OpenCL buffers for all above arrays (per replica slices).
   - pyOpenCL: allocates host numpy buffers matching kernel layout, uploads to device buffers; may reorder atoms (node-first) for kernel expectations.
4. **Aux data**:
   - Neighbor lists (`neighs` per atom up to 4) and back-neighbors (`bkneighs`) for assembly.
   - PBC: `invLvec/lvec`, optional `pbcShifts` per bond; pyOpenCL expects `neighCell` indices too.
   - REQ (LJ/Q) per atom for NBFF and bond/angle subtraction.
5. **Kernel/compute sequencing**:
   - CPU MMFFsp3_loc: `eval_atom(ia)` loops bonds -> angles -> pi interactions; accumulates into `fneigh/fneighpi`, then assembled to `fapos/fpipos`. Supports collision damping, subtracting bonded from NB @`MMFFsp3_loc.h#247-574`.
   - CPU NBFF: Coulomb/LJ with PBC; can subtract bonded pairs/angles when requested.
   - GPU (relax_multi_mini.cl): mirrored steps — bond stretch, angles, pi-sigma, pi-pi, NBFF; runs over replicas (2D NDRange in multi-replica).
   - GPU (UFF.cl): UFF bond/angle/dihedral/nonbond kernels similarly sequenced.
   - C++ driver (`OCL_MM`, `OCL_UFF`): uploads inputs, sets kernel args, enqueues kernel sequence per iteration; handles occasional downloads for convergence checks.
   - pyOpenCL (`MolecularDynamics.py`): same sequencing in Python; cleans buffers each cycle as needed, downloads energies/forces for parity tests.
6. **Parity testing**: Compare CPU vs C++-OpenCL vs pyOpenCL by running identical topology/params, inspecting energies/forces after fixed iterations; ensure cleanup/reset of temporary buffers (fneigh, accumulators) and consistent REQ/mixing rules.

## Cross-version equivalence (major structures/functions)

| Concept / Data             | CPU (C++)                                | GPU C++ driver / kernels                 | pyOpenCL (Python)                                  |
|----------------------------|------------------------------------------|------------------------------------------|----------------------------------------------------|
| Atomic arrays (pos/force)  | `apos/fapos` (Atoms)                     | CL buffers `apos`, `fapos`               | numpy `apos` host → CL buffer                      |
| Pi vectors                 | `pipos/fpipos` (MMFFsp3_loc)             | `pipos` in CL buffer (tail of apos in kernels) | stored then written to tail of `apos` @`pyBall/OCL/MMFF.py#556-560` |
| Neighbors (up to 4)        | `neighs` (`Quat4i` per atom)             | `fneigh`/`neighs` buffers                | `self.neighs` (numpy int32)                        |
| Back-neighbors             | `bkneighs` (inverse)                     | kernel expects back-neigh indices        | `make_back_neighs()` placeholder @`pyBall/OCL/MMFF.py#562-563` |
| Bond params                | `bLs`, `bKs` per atom (local) or `bond_l0/k` (bond-list) | CL buffers `bLs`, `bKs`                  | numpy `bLs`, `bKs` assigned in `toMMFFsp3_loc`     |
| Angle params               | `apars` (cos(ang0/2), Kss, Ksp, c0_e)    | CL buffer `apars`                        | numpy `apars` filled per atom                      |
| Pi interaction params      | `Ksp`, `Kpp`                             | CL buffers                               | numpy `Ksp`, `Kpp`                                 |
| REQ (NB params)            | `REQs` (ForceField/NBFF)                 | CL buffer `REQs`                         | numpy `REQs` (generated from atom types)           |
| PBC shifts/cells           | `neighCell`, `pbcShifts`, `invLvec/lvec` | CL buffers + kernel args                 | `neighCell`, lattice vectors passed to kernels     |
| Evaluation routine         | `MMFFsp3_loc::eval_atom`                 | `relax_multi_mini.cl` kernels            | kernel dispatch order in `MolecularDynamics.py`    |
| UFF bonded/NB              | `UFF.h`                                  | `UFF.cl` via `OCL_UFF`                   | `pyBall/OCL/UFF.py` driver                         |

## Notable call flows (CPU vs pyOpenCL) for MMFFsp3_loc

**CPU (conceptual)**
1. Builder loads geometry + params → fills `atypes`, `neighs`, `bLs/bKs`, `apars`, `Ksp/Kpp`, `REQs`.
2. `MMFFsp3_loc.realloc(...)` allocates DOFs/forces, neighbor param arrays.
3. Loop `eval_atom` for each atom; accumulates `fneigh/fneighpi` then assembled to `fapos/fpipos`; NBFF optionally subtracts bonded/angle pairs.

**pyOpenCL (in `MMFF.py.toMMFFsp3_loc`)**
1. Input `AtomicSystem` + `atom_types` (from MMparams / *.dat).
2. Optional reorder nodes-first; compute `npi`, `nep`, `nnode`, `ncap`, `nb`.
3. `realloc(...)` on Python side to size numpy buffers.
4. Fill `apos`, `atypes`, neighbor indices (`neighs`, `ngbs`), per-atom angle/ bond/ pi parameters (`apars`, `bLs/bKs`, `Ksp/Kpp`); build pi vectors and normalize @`pyBall/OCL/MMFF.py#247-560`.
5. Write pi vectors into tail of `apos` to match kernels @`#556-560`; create back-neighbor placeholder.
6. Upload buffers to device; `MolecularDynamics.py` runs kernel sequence.

## Parameter/data files usage

- **AtomTypes.dat**: element-level radii/charges, valence, Ksp/Kpp defaults.
- **BondTypes.dat**: bond equilibrium lengths/stiffness.
- **AngleTypes.dat**: angle equilibrium cos(ang0/2) and Kss.
- **DihedralTypes.dat**: torsion params (used in UFF).
- **ElementTypes.dat**: basic element properties (mass, cov radii, etc.).
- Both C++ builders and pyOpenCL MMparams/Builder read these to populate atom_types dicts or parameter tables.

## Suggested parity checklist

1. Ensure identical atom ordering (or consistent reorder-to-node-first) across CPU/GPU.
2. Confirm REQ mixing rules match between CPU (`_mixREQ`) and kernels.
3. Verify pi vector storage: pipos written at `apos[nAtoms : nAtoms+nnode]` for kernels.
4. Align neighbor/back-neighbor indexing (GPU expects inverse map).
5. PBC: same `neighCell` indices and lattice vectors; same cutoff/clamp flags.
6. Clean buffers each iteration (`fneigh`, accumulators) both CPU and GPU before kernels.
7. Kernel sequence parity: bond → angles → pi-sigma → pi-pi → NBFF; subtract bonded NB when enabled.

Let me know if you want deeper per-file notes (e.g., UFF.h specifics, MolWorld_sp3*_ orchestration, or the exact kernel argument order in OCL_MM/UFF), and I can extend the tables.


---

## GPT-5.2-Low

Got it. I’ll build you a **cross-variant code map** (CPU C++ vs C++/OpenCL vs pure pyOpenCL) focusing on the **execution path**: geometry+params → topology/type assignment → aux buffers (neighbors/PBC/replicas) → OpenCL buffer prep → kernel sequence → readback/convergence.  

Because this requires systematically reading a bunch of files, I’ll:
1. Locate the *entry points* and “main orchestration” functions in each layer (MolWorld, ForceField, OCL drivers, pyBall wrappers).
2. Extract the key datatypes/buffers and the call graphs.
3. Produce:
   - a **module table** (responsibility + key structs + dependencies)
   - an **equivalence table** (CPU vs GPU vs pyOCL “same concept” functions/kernels)
   - an **execution pipeline** diagram (stages + where implemented)
   - notes for **parity testing hooks** (what to compare, when to clear buffers, etc.)

I’ll start by searching within your worktree for these files (paths in the metadata point to main repo, but I must use the worktree copy).

# Code Map: Classical Forcefield System (CPU C++ / C++ OpenCL / pure pyOpenCL)

This is a **systematic cross-variant map** of your classical FF stack, focusing on:
- **Geometry + parameter loading**
- **Topology + type assignment**
- **Aux buffers** (neighbors, PBC shifts, exclusions, back-neighbors, atom→force assembly maps, replica packing)
- **Evaluation loop / kernel sequencing**
- **Parity-testing checkpoints** (what to compare between implementations)

I’m basing this strictly on the files you listed (from the worktree copy).

## 1) Big-picture architecture (layers and data flow)

## A. CPU C++ “single-system” path
- **Data model**
  - [Atoms](cci:2://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/Atoms.h:27:0-194:1) → `ForceField` → [NBFF](cci:2://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/NBFF.h:96:0-1312:1) → (bonded FFs)
  - Two bonded FF “families”:
    - **MMFFsp3_loc**: bonded (bonds/angles/pi terms) + uses NBFF for nonbonded
    - **UFF**: bonded (bonds/angles/dihedrals/inversions) + can subtract 1–2/1–3 nonbonded contributions + can also evaluate explicit nonbonded via [NBFF](cci:2://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/NBFF.h:96:0-1312:1) methods (in `eval_omp` path it does NB inside bond loop)
- **World orchestration**
  - [MolWorld_sp3](cci:2://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3.h:76:0-3258:1) owns:
    - builder/topology: `MM::Builder builder` (from [MMFFBuilder.h](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MMFFBuilder.h:0:0-0:0))
    - parameters: `MMFFparams params` (from `MMFFparams.h`)
    - forcefields: `MMFFsp3_loc ffl`, `UFF ffu`, `NBFF nbmol`, `NBFF surf`, `GridFF gridFF`, etc.
  - [MolWorld_sp3::init()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3_multi.h:371:0-464:1) is the central entry point: reads params, reads geometry/surface/SMILES, builds topology, calls [makeFFs()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3.h:1385:0-1485:5) (not in excerpt but referenced), etc.

## B. C++ OpenCL “many replicas” path
- `MolWorld_sp3_multi : MolWorld_sp3, MultiSolverInterface`
  - Calls **base** [MolWorld_sp3::init()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3_multi.h:371:0-464:1) to do *the same* geometry+topology+CPU FF build.
  - Then allocates/initializes OpenCL interface(s):
    - For **MMFF**: `OCL_MM ocl` using kernels from `relax_multi.cl` (or mini kernel file in repo)
    - For **UFF**: `OCL_UFF* uff_ocl` using kernels from [UFF.cl](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common_resources/cl/UFF.cl:0:0-0:0)
  - Packs **per-system replicas** into big SoA-ish arrays and uploads to GPU.
  - Runs kernel sequences for relaxation / evaluation.

## C. pure pyOpenCL path (no C++ driver)
- Two separate “tracks”:
  1) **MMFFsp3_loc-like GPU MD** using [pyBall/OCL/MolecularDynamics.py](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/MolecularDynamics.py:0:0-0:0) + [pyBall/OCL/MMFF.py](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/MMFF.py:0:0-0:0) and kernel [cpp/common_resources/cl/relax_multi_mini.cl](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common_resources/cl/relax_multi_mini.cl:0:0-0:0)
  2) **UFF GPU evaluation** using [pyBall/OCL/UFF.py](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/UFF.py:0:0-0:0) + [pyBall/OCL/UFFbuilder.py](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/UFFbuilder.py:0:0-0:0) and kernel [cpp/common_resources/cl/UFF.cl](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common_resources/cl/UFF.cl:0:0-0:0)
- This path replicates the C++ OpenCL buffer model (but some parts are still “in progress” / inconsistent—see notes below).

## 2) Module table (responsibilities + main dependencies)

## CPU / shared C++
| Module | Role | Key data it owns | Key dependencies / notes |
|---|---|---|---|
| [cpp/common/molecular/Atoms.h](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/Atoms.h:0:0-0:0) | Minimal atom container + XYZ IO | `natoms`, `atypes`, `apos`, optional `lvec*`, `charge*` | [atomsFromXYZ()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/Atoms.h:91:4-118:5) reads `atype` as integer directly; `lvec` can be read from comment line beginning with `lvs` |
| `cpp/common/molecular/ForceField.h` | MD/relaxation utilities (integrators, velocity, force buffers) | `fapos`, `vapos`, PBC flags, `shifts`, nonbond strategy flags | `setNonBondStrategy()` configures subtract/clamp logic; `move_*` integrators; `cleanForce()` |
| [cpp/common/molecular/NBFF.h](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/NBFF.h:0:0-0:0) | Nonbonded core (LJQH/MorseQH etc), PBC shifts, exclusions | `REQs`, `neighs`, `neighCell`, `excl`, `PLQs`, `shifts`, `nPBC`, etc. | [makePBCshifts()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/NBFF.h:1258:4-1272:5) and [makeSecondNeighs()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/NBFF.h:205:4-230:5) (EXCL table), `evalLJQs_*` variants including `ex2` |
| `cpp/common/molecular/MMFFsp3_loc.h` | **Bonded** “sp3 local” FF + uses NBFF | `apars,bLs,bKs,Ksp,Kpp`, `pipos`, `fneigh`, `bkneighs`, `neighs`, `neighCell`, constraints | [realloc()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3_multi.h:252:0-326:1) binds `apos/fapos` into `DOFs/fDOFs` layout; recoil stored in `fneigh` then assembled |
| `cpp/common/molecular/UFF.h` | UFF bonded terms + (optionally) nonbond subtract | `bonAtoms/bonParams`, `angAtoms/angParams/angNgs`, `dih*`, `inv*`, `hneigh`, `fint`, `a2f` mapping | [eval()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/MMFF.py:1139:0-1140:21) calls components, then `assembleAtomsForces()` via atom→force map; `evalAtomBonds()` fills `hneigh` as it goes |
| `cpp/common/molecular/MMFFparams.h` | Reads parameter tables; provides typed lookup (bond/angle/dihedral types, etc.) | atom/bond/angle/dihedral type tables | Used by `builder` and world init; contains the schema used by `AtomTypes.dat` etc. |
| [cpp/common/molecular/MMFFBuilder.h](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MMFFBuilder.h:0:0-0:0) | Builds topology + assigns types + builds FF instances | atoms/bonds/angles/dihedrals, confs, fragments, type assignment logic | Contains UFF bond formula helper [assignBondParamsUFF()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/MMFF.py:664:4-753:25) and many MMFF assignment utilities |
| [cpp/common/molecular/MolWorld_sp3.h](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3.h:0:0-0:0) | “simulation world” orchestrator on CPU | owns builder, params, FF objects, flags, IO names, run loops | [init()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3_multi.h:371:0-464:1) is key pipeline; selects between `bMMFF` and `bUFF` workflows |
| [cpp/libs/Molecular/MMFF_lib.cpp](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/libs/Molecular/MMFF_lib.cpp:0:0-0:0) | C ABI / Python bridge to `MolWorld_sp3 W` | global `MolWorld_sp3 W`, buffer registry | [init()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3_multi.h:371:0-464:1) configures W flags, calls [W.initParams()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/MMFF.py:1020:0-1027:126) then [W.init()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3_multi.h:371:0-464:1); [init_buffers_UFF()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/libs/Molecular/MMFF_lib.cpp:23:0-63:1) exposes UFF arrays |
| [pyBall/MMFF.py](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/MMFF.py:0:0-0:0) | Python ctypes wrapper to `MMFF_lib` | interfaces to [init](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3_multi.h:371:0-464:1), [setSwitches](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/MMFF.py:1073:0-1074:97), [getArrayPointer](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/libs/Molecular/MMFF_lib.cpp:180:0-191:1), etc. | Provides [makeGridFF](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/MMFF.py:306:0-316:15), Ewald utilities, and various sampling helpers |

## GPU OpenCL C++ driver
| Module | Role | Key data/buffers | Notes |
|---|---|---|---|
| [cpp/common/OpenCL/OCL_MM.h](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/OpenCL/OCL_MM.h:0:0-0:0) | OpenCL system for MMFFf4/NBFF/GridFF/MD update | `ibuff_atoms/aforces/avel/cvf`, `ibuff_neighs`, `ibuff_REQs`, `ibuff_*params`, `ibuff_pbcshifts`, etc. | Builds tasks from `relax_multi.cl` (or variants). Most kernels use **2D NDRange**: `global.y=nSystems` |
| `cpp/common/OpenCL/OCL_UFF.h` | OpenCL system for UFF | buffers for apos/fapos/REQs, topology (bond/angle/dih/inv), `hneigh`, `fint`, plus a2f mapping buffers | Builds tasks from [UFF.cl](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common_resources/cl/UFF.cl:0:0-0:0): `evalBondsAndHNeigh_UFF`, `evalAngles_UFF`, `evalDihedrals_UFF`, `evalInversions_UFF`, `assembleForces_UFF`, `getNonBond*`, `updateAtomsMMFFf4` |
| [cpp/common/molecular/MolWorld_sp3_multi.h](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3_multi.h:0:0-0:0) | Multi-replica orchestrator | host-packed arrays: `atoms,aforces,avel,REQs,...`, topology vectors, per-system lvec/pbc shifts | [init()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3_multi.h:371:0-464:1) calls base init then OpenCL init; [pack_system()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3_multi.h:587:0-675:1) / [pack_uff_system()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3_multi.h:800:0-869:1) fill host arrays; `upload_*` methods upload slices |

## pure pyOpenCL driver
| Module | Role | Notes |
|---|---|---|
| [pyBall/OCL/MolecularDynamics.py](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/MolecularDynamics.py:0:0-0:0) | OpenCL runner for MMFF kernels (MD/relax) | Mirrors [OCL_MM](cci:2://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/OpenCL/OCL_MM.h:54:0-1517:1) buffer set, loads [relax_multi_mini.cl](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common_resources/cl/relax_multi_mini.cl:0:0-0:0), builds kernel args by parsing headers |
| [pyBall/OCL/MMFF.py](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/MMFF.py:0:0-0:0) | Builds MMFFsp3_loc-like topology/params in Python | Has node-first reorder machinery, computes `neighs`, `bLs/bKs`, `apars`, `pipos`, `back_neighs` etc. |
| [pyBall/OCL/UFFbuilder.py](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/UFFbuilder.py:0:0-0:0) | Python implementation of UFF type assignment + parameter assignment | Attempts to mirror `assignUFFtypes*` logic; returns arrays similar to C++ `UFF` |
| [pyBall/OCL/UFF.py](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/UFF.py:0:0-0:0) | PyOpenCL runner for UFF kernels | Allocates buffers mirroring `OCL_UFF`/[UFF.cl](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common_resources/cl/UFF.cl:0:0-0:0); runs kernel sequence |

## 3) Equivalent concepts/functions table (CPU vs C++/OCL vs pyOpenCL)

This is the “Rosetta stone” you asked for: same conceptual step, different implementation points.

## A. Geometry / system state
| Concept | CPU C++ | C++ OpenCL multi | pyOpenCL |
|---|---|---|---|
| Atom container | [Atoms](cci:2://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/Atoms.h:27:0-194:1) / `ForceField` | packed `Quat4f* atoms` (per system) | `apos` float4 buffers ([MolecularDynamics](cci:2://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/MolecularDynamics.py:36:0-640:21), [UFF_CL](cci:2://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/UFF.py:12:0-322:51)) |
| Set positions | [MolWorld_sp3::setGeom()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3_multi.h:1428:0-1488:1) writes `ffl.apos` etc | [MolWorld_sp3_multi::setGeom(isys,...)](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3_multi.h:1428:0-1488:1) updates per replica and packs pbc shifts | [pack_system()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3_multi.h:587:0-675:1) uploads to GPU buffers |
| PBC shifts | [NBFF::makePBCshifts()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/NBFF.h:1258:4-1272:5) and/or [MolWorld_sp3::evalPBCshifts()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3.h:744:0-753:1) | stored in `pbcshifts[]` and uploaded to `ibuff_pbcshifts` | `pbc_shifts` buffer in [MolecularDynamics](cci:2://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/MolecularDynamics.py:36:0-640:21) and [UFF_CL](cci:2://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/UFF.py:12:0-322:51) |

## B. Nonbonded (NBFF)
| Concept | CPU C++ | C++ OpenCL | pyOpenCL |
|---|---|---|---|
| LJQH pair | `NBFF::getLJQH()` via `Forces.h` | `getNonBond*` kernels in `relax_multi*.cl` / [UFF.cl](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common_resources/cl/UFF.cl:0:0-0:0) | uses same kernel code (shared `.cl`) |
| Exclusions (1-2, 1-3) | `bSubtractBondNonBond`, `bSubtractAngleNonBond`, `excl` table for ex2 | `getNonBond_ex2` uses `ibuff_excl` | [MolecularDynamics](cci:2://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/MolecularDynamics.py:36:0-640:21) allocates `excl` but pack path not fully wired (MMFF side) |
| Second-neighbor exclusion table | [NBFF::makeSecondNeighs()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/NBFF.h:205:4-230:5) builds `excl[natoms*EXCL_MAX]` | [MolWorld_sp3_multi::pack_system()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3_multi.h:587:0-675:1) uploads `excl` when `bExclusion2` | not clearly finished for MMFF pyOCL path |

## C. Bonded: MMFFsp3_loc
| Concept | CPU C++ (`MMFFsp3_loc`) | C++ OpenCL ([OCL_MM](cci:2://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/OpenCL/OCL_MM.h:54:0-1517:1) + `relax_multi*.cl`) | pyOpenCL ([MolecularDynamics](cci:2://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/MolecularDynamics.py:36:0-640:21) + [pyBall/OCL/MMFF.py](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/MMFF.py:0:0-0:0)) |
|---|---|---|---|
| Per-node local param arrays | `apars,bLs,bKs,Ksp,Kpp` | `ibuff_MMpars, ibuff_BLs, ibuff_BKs, ibuff_Ksp, ibuff_Kpp` | same buffers allocated; uploaded in [pack_system()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3_multi.h:587:0-675:1) |
| Neighbor list | `neighs[natoms]` + `neighCell[natoms]` | `ibuff_neighs`, `ibuff_neighCell` | same |
| Back-neighbor mapping | `bkneighs` / `bkNeighs_new` on GPU | `ibuff_bkNeighs`, `ibuff_bkNeighs_new` | py side builds `back_neighs` and uploads as `bkNeighs` |
| Force accumulation strategy | recoil in `fneigh[]` then assemble to `fapos` | same via `fneigh` buffer; kernels compute recoil then update kernel assembles | `fneigh` allocated in [MolecularDynamics](cci:2://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/MolecularDynamics.py:36:0-640:21), but you should verify `runMD` uses it consistently |
| Core bonded kernel | `eval_atom()` / [eval()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/MMFF.py:1139:0-1140:21) style functions (not fully excerpted) | `__kernel getMMFFf4` (and `evalMMFFf4_local*`) | `prg.getMMFFf4` invoked in [run_MD_py](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/MolecularDynamics.py:392:4-407:68)/[run_step_basic](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/MolecularDynamics.py:409:4-417:27) |

## D. Bonded: UFF
| Concept | CPU C++ (`UFF`) | C++ OpenCL (`OCL_UFF` + [UFF.cl](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common_resources/cl/UFF.cl:0:0-0:0)) | pyOpenCL ([UFF_CL](cci:2://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/UFF.py:12:0-322:51) + [UFF.cl](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common_resources/cl/UFF.cl:0:0-0:0)) |
|---|---|---|---|
| Bonds | `evalBonds()` → `evalAtomBonds()` | `evalBondsAndHNeigh_UFF` | same kernel |
| Angles | `evalAngles()` → `evalAngle_Prokop()` | `evalAngles_UFF` | same kernel |
| Dihedrals | `evalDihedrals()` → `evalDihedral_Prokop()` | `evalDihedrals_UFF` | same kernel |
| Inversions | `evalInversions()` | `evalInversions_UFF` | same kernel |
| H-neigh cache | `hneigh[natoms*4]` filled during bond eval | `hneigh` buffer filled by bond kernel | same |
| Assembly model | CPU has `fint[]` and atom→force map `a2f` and `assembleAtomsForces()` | GPU has `fint` + `assembleForces_UFF` with a2f buffers | pyOpenCL has a2f buffers **allocated** but current [run_eval_step()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/UFF.py:157:4-168:29) does **not call assemble** (important gap) |

## 4) End-to-end execution pipeline (compare CPU vs GPU vs pyOCL)

## Stage 0: Read parameter tables
- **C++**
  - [MolWorld_sp3::init()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3_multi.h:371:0-464:1) calls [initParams(...)](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/MMFF.py:1020:0-1027:126) if `params.atypes.size()==0`
  - Default paths in C++ init: `common_resources/ElementTypes.dat`, `AtomTypes.dat`, `BondTypes.dat`, `AngleTypes.dat`, `DihedralTypes.dat`
- **Python ctypes (CPU C++ backend)**
  - [MMFF_lib.cpp::init()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3_multi.h:371:0-464:1) calls [W.initParams(sElementTypes, sAtomTypes, ...)](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/MMFF.py:1020:0-1027:126)
  - `pyBall/MMFF.py:init()` passes default `"data/*.dat"` paths (note: your tests may use different relative path roots)
- **pyOpenCL**
  - [pyBall/OCL/MMparams.py](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/MMparams.py:0:0-0:0) reads the same `ElementTypes.dat` and `AtomTypes.dat` and builds dicts.

**Parity checkpoint**
- Verify that **same dat files** are used (and relative path roots match). This is a common source of “mysterious differences”.

## Stage 1: Read geometry / construct initial molecule
- **C++ world**
  - [MolWorld_sp3::init()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3_multi.h:371:0-464:1):
    - If `surf_name`: [loadSurf()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3_multi.h:2925:0-2939:1) and sets `bGridFF=true` (also note `z0` handling differs: CPU uses `NAN` in code excerpt)
    - If `smile_name`: parse SMILES, `builder.addAllCapTopo()`, randomize, sets `bMMFF=true`
    - Else if `xyz_name` and (`bMMFF` or `bUFF`): [buildMolecule_xyz(xyz_name)](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3.h:1279:0-1314:5)
    - Else: [loadNBmol(xyz_name)](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3.h:1005:0-1018:5) as purely nonbonded-only system
- **pyOpenCL**
  - Depends on `AtomicSystem` object in pyBall (not in your file list), but [pyBall/OCL/MMFF.py](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/MMFF.py:0:0-0:0) expects:
    - `mol.apos`, `mol.enames`, and `mol.bonds` or ability to create them (`mol.findBonds()` is used in UFFbuilder)

**Parity checkpoint**
- Ensure atom ordering is consistent.
  - pyOpenCL MMFF builder often **reorders nodes first** (`reorder_nodes_first=True`). That means *indices no longer match C++* unless you apply the same permutation in C++ or map back when comparing.

## Stage 2: Topology building + type assignment
### MMFFsp3 path (your “MMFFsp3_loc”)
- **C++**: done inside builder + [makeFFs()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3.h:1385:0-1485:5) in [MolWorld_sp3](cci:2://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3.h:76:0-3258:1)
  - [MM::Builder](cci:2://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MMFFBuilder.h:92:0-3923:1) owns atoms/bonds/confs and assigns:
    - atom types (sp1/sp2/sp3 subTypes)
    - bond parameters via `params->getBondParams(...)`
    - angle parameters (Kss etc. from `AtomTypes.dat` columns)
- **pyOpenCL**: [pyBall/OCL/MMFF.py::toMMFFsp3_loc()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/MMFF.py:254:4-562:31)
  - Constructs:
    - `neighs` from `mol.ngs`
    - `apars` from `AtomType.Ass/Kss`
    - `bLs/bKs` from either simple rule or UFF-derived rule
    - `Ksp/Kpp` assignments and pi directions `pipos`
    - `back_neighs` via [make_back_neighs()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/MMFF.py:168:4-191:65)

### UFF path
- **C++**
  - UFF typing logic is in [MolWorld_sp3](cci:2://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3.h:76:0-3258:1) + `UFFbuilder.h` (included), and the actual runtime FF is `UFF` in `UFF.h`.
- **pyOpenCL**
  - [pyBall/OCL/UFFbuilder.py](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/UFFbuilder.py:0:0-0:0) implements:
    - `assign_uff_types_*` (trivial/nitro/treewalk/simplerule/findrings/assignrest/fixsaturation/cumulene/conjugation/amide/checks)
    - `assign_uff_params_*` to produce `bonParams/angParams/dihParams/invParams`
  - Returns arrays: `bonAtoms`, `bonParams`, `angAtoms`, `angParams`, `dihAtoms`, `dihParams`, `invAtoms`, `invParams`, plus `neighs`, `neighBs`.

**Parity checkpoint**
- Compare **typed atom-type indices** and bond orders between:
  - C++ UFFbuilder output and Python [UFFbuilder.py](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/UFFbuilder.py:0:0-0:0)
- For MMFFsp3_loc, compare:
  - `neighs` arrays
  - `bLs/bKs`
  - `apars` (especially `cos(Ass/2), sin(Ass/2), 4*Kss, sin(Ass)` conventions in pyOpenCL code)

## Stage 3: Prepare aux buffers (PBC, exclusions, back-neigh, assembly maps)
### PBC shifts
- C++:
  - [MolWorld_sp3::makePBCshifts()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/NBFF.h:1258:4-1272:5) allocates `pbc_shifts` and calls [evalPBCshifts()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3.h:744:0-753:1)
  - [NBFF::makePBCshifts()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/NBFF.h:1258:4-1272:5) exists too; both patterns appear in the codebase.
- GPU:
  - [OCL_MM](cci:2://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/OpenCL/OCL_MM.h:54:0-1517:1) expects `pbcshifts` as `float4` buffer sized `nSystems*npbc`
  - [MolWorld_sp3_multi::pack_system()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3_multi.h:587:0-675:1) packs `pbc_shifts` → `pbcshifts` and uploads
- pyOpenCL MD:
  - [MolecularDynamics](cci:2://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/MolecularDynamics.py:36:0-640:21) allocates `pbc_shifts` buffer but you must ensure it’s populated; current [pack_system()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3_multi.h:587:0-675:1) doesn’t show it being uploaded.

### Exclusions (`excl` = second neighbor list)
- CPU: [NBFF::makeSecondNeighs()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/NBFF.h:205:4-230:5) writes `excl[natoms*EXCL_MAX]`
- C++ OpenCL: [MolWorld_sp3_multi::pack_system()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3_multi.h:587:0-675:1) uploads `excl` if `bExclusion2`
- pyOpenCL: buffer exists but I don’t see consistent population for MMFF yet.

### UFF atom→force assembly map (`a2f`)
- CPU UFF: `assembleAtomForce()` uses `a2f.cellI0s/cellNs/cell2obj` to sum `fint` pieces into `fapos`
- C++ OpenCL UFF: `OCL_UFF` has `ibuff_a2f_offsets/counts/indices` and a dedicated `assembleForces_UFF` kernel
- pyOpenCL UFF:
  - [UFF_CL.mapAtomInteractions()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/UFF.py:243:4-322:51) builds arrays, and [set_a2f_map_size()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/UFF.py:227:4-241:44) allocates buffers
  - **But** [run_eval_step()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/UFF.py:157:4-168:29) currently does not call an assembly kernel, and `fint` buffer size is wrong-ish (it allocates `fint` as `nA*4` floats, but in C++ `fint` is sized to “force pieces”, not atoms). So expect parity problems unless this has been fixed elsewhere.

**Parity checkpoint**
- For UFF GPU parity, the most robust compare is:
  - Compare `hneigh` after bond kernel
  - Compare `fint` entries for angles/dihedrals/inversions **before assembly**
  - Then compare assembled `fapos`
- For MMFF GPU parity:
  - Compare recoil buffers `fneigh` and final `aforce`

## Stage 4: GPU buffer packing and offsets (multi replicas)

This is where many subtle bugs live.

## C++ OpenCL multi: canonical offset scheme
- **MMFF path** ([MolWorld_sp3_multi::pack_system()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3_multi.h:587:0-675:1)):
  - `i0a = isys * ocl.nAtoms` for atom-sized arrays (neighs, REQs, constraints)
  - `i0v = isys * ocl.nvecs` for vector-sized arrays (atoms + pi)
  - `i0n = isys * ocl.nnode` for node-sized arrays (apars/bLs/bKs/Ksp/Kpp)
  - `i0pbc = isys * ocl.npbc`
- **UFF path** ([pack_uff_system()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3_multi.h:800:0-869:1)):
  - `i0a = isys * natoms`
  - topology arrays are per-interaction counts: `i0b = isys*nBonds`, etc.
- **Kernel launch** is typically **2D NDRange** with `(global.x ~ nAtomsRounded, global.y = nSystems)`.

## pyOpenCL
- [MolecularDynamics.pack_system()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/MolecularDynamics.py:169:4-207:122) uses byte offsets consistent with `iSys * nvecs * sizeof(float4)` etc. That’s good.
- [UFF_CL.upload_topology_params()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/UFF.py:129:4-155:82) uses offsets in **elements**, then multiplies by element size—also good conceptually.

**Parity checkpoint**
- Always confirm that:
  - host array slice offsets == device buffer offsets
  - kernel uses the same per-system base `i0 = isys * stride`
This is consistent with your own retrieved memory rule about offsets; the C++ multi code follows that pattern.

## Stage 5: Kernel sequence / evaluation loop

## MMFF (GPU) typical per-step sequence
From [pyBall/OCL/MolecularDynamics.py](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/MolecularDynamics.py:0:0-0:0) (and consistent with [OCL_MM.h](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/OpenCL/OCL_MM.h:0:0-0:0) tasks):
1. `cleanForceMMFFf4` (or nonbond kernel overwriting, depending on mode)
2. `getNonBond` / `getNonBond_ex2` / `getNonBond_GridFF_*`
3. `getMMFFf4` (bonded)
4. `updateAtomsMMFFf4` (integrator / FIRE-like update)

In C++ multi, there are multiple run modes ([run_omp_ocl](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3_multi.h:2355:0-2490:1), [run_ocl_opt](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3_multi.h:2101:0-2263:1), etc.), but the kernel building in [OCL_MM::makeKrenels_MM()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/OpenCL/OCL_MM.h:167:4-200:5) suggests the same building blocks.

## UFF (GPU) typical evaluation sequence
From `OCL_UFF::makeKernels()` / [UFF.cl](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common_resources/cl/UFF.cl:0:0-0:0) kernels:
1. `clear_fapos_UFF` (and maybe `clear_fint_UFF`)
2. `evalBondsAndHNeigh_UFF` (produces `hneigh` + writes/accumulates into `fapos`)
3. `evalAngles_UFF` (writes angle force pieces to `fint`)
4. `evalDihedrals_UFF` (writes torsion pieces to `fint`)
5. `evalInversions_UFF` (writes inversion pieces to `fint`)
6. `assembleForces_UFF` (reduces `fint` via a2f mapping into `fapos`)
7. optional nonbond: `getNonBond*` or gridFF variants
8. optional `updateAtomsMMFFf4` (used as generic integrator kernel)

**Important mismatch to note**
- `pyBall/OCL/UFF.py:run_eval_step()` currently runs 2–5 but **does not call** step 6 (`assembleForces_UFF`). If so, `fapos` will miss angle/dihedral/inversion contributions unless the kernels directly accumulate into `fapos` (they don’t; they store into `fint`).

## 5) Where parameters live (dat files) and who consumes them

## File roles
- `ElementTypes.dat`
  - element-level radii, vdW, UFF electronegativity (`Eaff`), etc.
- `AtomTypes.dat`
  - atom-type resolved `Ruff/RvdW/EvdW/Qbase/Hb` and also MMFF-like parameters `Ass/Asp/Kss/Ksp/Kep/Kpp` (optional columns)
- `BondTypes.dat`, `AngleTypes.dat`, `DihedralTypes.dat`
  - tables used by C++ [MMFFparams](cci:2://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/MMparams.py:244:0-261:51) for topology-based typed terms

## Consumption patterns
- **MMFFsp3_loc** mostly uses:
  - `AtomTypes.dat` for angle + pi stiffness constants
  - `BondTypes.dat` (and maybe AngleTypes etc.) via `MMFFparams::getBondParams(...)`
- **UFF** uses:
  - `AtomTypes.dat` and `ElementTypes.dat` to derive bond lengths/force constants (see [Builder::assignBondParamsUFF()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/MMFF.py:664:4-753:25) formula)
  - `UFF` also stores `REQs` for nonbond mixing (LJQH)

## 6) Practical parity-testing checklist (what to compare, when)

This is the “systematic” part that tends to save time.

## A. Geometry/topology invariants (before any force eval)
- **Atom order**
  - If Python reorders nodes first, record `perm` and compare in mapped index space.
- **Neighbor tables**
  - Compare `neighs` and `neighCell` (and `neighBs` for UFF) exactly.
- **Parameter arrays**
  - MMFF:
    - `apars`, `bLs`, `bKs`, `Ksp`, `Kpp`
    - `pipos` (direction vectors) and whether they’re normalized
  - UFF:
    - `bonAtoms/bonParams`
    - `angAtoms/angParams` and `angNgs`
    - `dihAtoms/dihParams` and `dihNgs`
    - `invAtoms/invParams` and `invNgs`

## B. One-step evaluation parity (forces only; no integration)
### MMFF
- CPU: call bonded+nonbond evaluation once ([MolWorld_sp3::eval_no_omp()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3.h:2078:0-2113:1) likely; not shown but used in Hessian code)
- GPU: run `cleanForce` → `getNonBond*` → `getMMFFf4` (skip `updateAtoms`)
- Compare:
  - per-atom forces
  - optionally `fneigh` recoil and per-node contributions

### UFF
- CPU: [UFF::eval(bClean=true)](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/MMFF.py:1139:0-1140:21) (includes assembly)
- GPU: `clear_*` → `evalBondsAndHNeigh_UFF` → `evalAngles_UFF` → `evalDihedrals_UFF` → `evalInversions_UFF` → `assembleForces_UFF`
- Compare:
  - `hneigh` (after bond step)
  - `fint` slices (after each kernel)
  - final `fapos`

## C. Relaxation parity (multi-step)
- Run a few steps and compare:
  - force norms decreasing similarly
  - energies by component (Eb/Ea/Ed/Ei) if available on both sides
- Ensure buffers are cleared properly each iteration:
  - MMFF: `cleanForceMMFFf4` resets `aforce` and `fneigh`
  - UFF: need both `clear_fapos_UFF` and `clear_fint_UFF` (depending on whether kernels overwrite or accumulate)

## 7) Gaps / “hard to navigate” hotspots you should flag in the map

These are places where the variants diverge and parity testing often breaks.

## A. Node ordering differences (MMFF pyOCL)
- [pyBall/OCL/MMFF.py](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/MMFF.py:0:0-0:0) may reorder nodes first.
- C++ CPU builder likely has its own ordering ([builder.sortConfAtomsFirst()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MMFFBuilder.h:1859:4-1869:5) etc. appears in places).
- If you compare forces atom-by-atom without mapping indices, you’ll get false mismatches.

## B. UFF pyOCL assembly path looks incomplete
- [UFF_CL.run_eval_step()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/UFF.py:157:4-168:29) doesn’t call `assembleForces_UFF` (in the excerpt).
- [UFF_CL](cci:2://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/UFF.py:12:0-322:51) allocates `fint` as `nA*4*f32` which is not consistent with “number of force pieces” logic in `OCL_UFF` (`nf_per_system = 2*nBonds + 3*nAngles + 4*nDihedrals + 4*nInversions`).
- This likely needs finishing before UFF pyOCL can be parity-tested against CPU/GPU-C++ reliably.

## C. PBC shifts upload in pyOCL MMFF MD
- [MolecularDynamics.allocate_cl_buffers()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/MolecularDynamics.py:78:4-149:58) allocates `pbc_shifts`, but [pack_system()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3_multi.h:587:0-675:1) doesn’t upload it (in excerpt).
- If you rely on PBC in MMFF pyOCL kernels, make sure that buffer is populated.

# Completion status
- **Completed**
  - **CPU C++ map** (Atoms→ForceField→NBFF→MMFFsp3_loc/UFF, orchestrated by MolWorld_sp3, exposed via MMFF_lib.cpp + pyBall/MMFF.py)
  - **C++ OpenCL map** (OCL_MM/OCL_UFF + kernels + MolWorld_sp3_multi packing & upload)
  - **pyOpenCL map** (MolecularDynamics/MMFF and UFF_CL/UFFbuilder, including buffer logic and kernel sequence)
  - Produced **module table**, **equivalence table**, and **pipeline/parity checklist**
- **Open question (for you)**
  - Do you want this refined into a *single “call graph style” diagram* per variant (with exact function names like [buildMolecule_xyz()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3.h:1279:0-1314:5) → [makeFFs()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3.h:1385:0-1485:5) → `builder.toMMFF...()` etc.)? I can do that next, but I’d want you to confirm which entry point(s) you use most:
    - [pyBall/MMFF.py](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/MMFF.py:0:0-0:0) via `MMFF_lib` (ctypes),
    - `MolWorld_sp3_multi` GPU C++ path,
    - or pure `pyBall/OCL/*` path.


---

## USER

Yes, lets go deeper, esential taks is to make parity check betwen different versions, in particular between newest pyOpenCL path and the original C++/OpenCL path. We want to do it using python interface.

I forgot mention there is this
@MMFFmulti_lib.cpp and @MMFF_multi.py for interface to C++/OpenCL path from python. This is essential for comparing to pyOpenCL


There are three main strategy how to do such tests
1) we download and compare the internal buffers (list of parameters, neibhors - topology, positions, forces) and compare them after single iterations first, then multiple iterations
2) if there are discrepancies in ordering/indexing, we need to go inside the program and introduce some debug prints to understain the exactution flow and where they start to deviate. This debug prints mut be gated because we typically execute large number of item (many iterations, many atoms etc), and without carefill gating and filtering (selected atoms etc.) the output would be overhelmed and it may even freeze, definitely difficult to analyse. To easy comparison (eventually with diff) we should make sure debug prints have exactly he same format, Preferebly the formating should be with long lines (no line leng limit), e.g. all controbitions from one atom on single line, and with partfecly aligned number in colum (no number lengh variation due to +/- or different number of digits).

in all cases we run the same simple system (start with e.g. H2O, then HCOOH, then pyrrolo, cytosine, guanine...) and small number of iterations

we already have some test for UFF here:
@test_MMFFsp3_pyOCL.py @test_UFF_multi.py @test_UFF_ocl.py , not sure if they cover it perfectly and if we have something like that for MMFF, but perhaps we can make something analogous

this debugging can be easily automatically analysed by LLM from the terminal output and using diff / grep / sed

3) visual debugging - by user - if we need to understand some higher level emergent behaviour (e.g. conservation of momentum, drift, non-phsyical dynamics), the user need to see what the molecue is moving (animation of simulation like in game) and how it respond e.g. when atoms are pulled by mouse. We already have used such debuging tools in python here:
@XPTB_utils.py @test_rigid_XPBD_molecules.py @XPBD_2D_utils.py @test_XPBD_2D.py 

analyse it and figure out how we can reuse these debugging tools , especially LiveViz2D, LivePortViz with mouse-picking and update callbacks in simulation loop, for efficient visual debbuging. We do not want to reimplemnt it, but rhater just reuse it / import it, we can also modify it to fit different structure of the problems (instead of ports we have pi-orbitals, angles and dihedrasl...). But the general thing (positions of atoms, atomic types, bonds between atoms) are shared for all molecular forcefields

---

## GPT-5.1-Codex

Here’s a focused parity plan between **pyOpenCL** and **C++/OpenCL** using the Python interfaces ([pyBall/MMFF_multi.py](cci:7://file:///home/prokop/git/FireCore/pyBall/MMFF_multi.py:0:0-0:0) → `cpp/libs_OCL/MMFFmulti_lib.cpp`) plus debug/visual tooling. I’m keeping it structured so you can implement stepwise.

## A) Parity harness: buffer download + single/multi-iteration checks
Use [pyBall/MMFF_multi.py](cci:7://file:///home/prokop/git/FireCore/pyBall/MMFF_multi.py:0:0-0:0) as the bridge to the C++/OpenCL path:

1) **Initialize both stacks with the same inputs**
   - Same `.xyz` (H2O → HCOOH → pyrrole → cytosine → guanine).
   - Same param files (`ElementTypes.dat`, `AtomTypes.dat`, `BondTypes.dat`, `AngleTypes.dat`, `DihedralTypes.dat`).
   - Same PBC (`nPBC`), same flags (`SubtractBondNonBond`, `ClampNonBonded`, etc.).

2) **Expose buffers for C++/OpenCL**
   - After [init()](cci:1://file:///home/prokop/git/FireCore/pyBall/MMFF_multi.py:426:0-443:248) and [getBuffs()](cci:1://file:///home/prokop/git/FireCore/pyBall/MMFF_multi.py:221:0-268:66) (MMFF) or [getBuffs_UFF()](cci:1://file:///home/prokop/git/FireCore/pyBall/MMFF_multi.py:272:0-330:15) (UFF), pull:
     - Topology: `neighs`, `neighCell`, `bkNeighs`, `REQs`, `apars/bLs/bKs/Ksp/Kpp` (MMFF); `bonAtoms/bonParams`, `angAtoms/angParams`, `dihAtoms/dihParams`, `invAtoms/invParams`, `neighBs`, `hneigh`, `fint` (UFF).
     - Geometry: `apos`, `pipos` (tail), `lvec/ilvec` if PBC.
   - For GPU-side slices: use `gpu_*` buffers in [getBuffs()](cci:1://file:///home/prokop/git/FireCore/pyBall/MMFF_multi.py:221:0-268:66) to inspect packed replicas (stride: atoms/nvecs/nnode per system). Confirm offsets match `isys * stride`.

3) **Mirror pyOpenCL packing**
   - Run the pyOpenCL builder (`pyBall/OCL/MMFF.py.toMMFFsp3_loc` or `UFFbuilder.py`) on the same AtomicSystem.
   - Collect numpy buffers before upload: `neighs`, `neighCell`, `back_neighs`, `apars`, `bLs/bKs`, `Ksp/Kpp`, `pipos`, `REQs` (MMFF), and UFF topology/params.
   - If node-first reorder is used, record `perm`/`perm_inv`; remap for comparison.

4) **Single-iteration force parity**
   - C++/OCL: [upload()](cci:1://file:///home/prokop/git/FireCore/pyBall/MMFF_multi.py:569:0-570:52) → run one bonded+NB kernel pass (no integrator) → [download(bForces=True)](cci:1://file:///home/prokop/git/FireCore/pyBall/MMFF_multi.py:575:0-576:38).
   - pyOpenCL: enqueue `cleanForceMMFFf4` → `getNonBond`/`getNonBond_ex2` → `getMMFFf4` (MMFF) or UFF kernel sequence (ensure `assembleForces_UFF` is called) → read `fapos`.
   - Compare `fapos` (and optionally recoil buffers `fneigh` / `fint`).

5) **Multi-iteration drift check**
   - Run N steps (e.g., 10) with same `dt`, damping, clamp flags. Download forces/energies every step; check norms and energy components (Eb/Ea/Epp/Eps for MMFF; Eb/Ea/Ed/Ei for UFF).

## B) Debug-print gating (for ordering/index mismatches)
Add gated prints in both C++ and pyOpenCL with identical line formats (single-line per atom/interaction; fixed-width, no wrapping). Gate by:
- System index (`isys`), atom index list, and/or iteration modulus.
- Verbosity/idebug flags passed through [setVerbosity(verbosity, idebug)](cci:1://file:///home/prokop/git/FireCore/pyBall/MMFF_multi.py:420:0-421:48) or dedicated env.

Recommended formats:
- Atom state: `AT i=%4d sys=%3d pos=(% .6f,% .6f,% .6f) f=(% .6f,% .6f,% .6f) pi=(% .6f,% .6f,% .6f)`
- Bond contrib: `BOND sys=%3d i=%4d j=%4d k=% .6f l0=% .6f l=% .6f fb=(% .6f,% .6f,% .6f)`
- Angle contrib: `ANGL sys=%3d i=%4d j=%4d k=%4d K=% .6f c0h=(% .6f,% .6f) fij=(% .6f,% .6f,% .6f) fjk=(% .6f,% .6f,% .6f)`
- Pi terms: `PPI sys=%3d i=%4d j=%4d kpp=% .6f dot=% .6f fpi=(% .6f,% .6f,% .6f) fpj=(% .6f,% .6f,% .6f)`
- Nonbond pair: `NB sys=%3d i=%4d j=%4d REQ=(% .6f,% .6f,% .6f) r=% .6f f=(% .6f,% .6f,% .6f)`

Keep identical ordering of fields in pyOpenCL (host-side debug wrappers) and C++ (printf) so diff works.

## C) Specific gaps to address for reliable parity
1) **UFF pyOpenCL assembly**: Ensure `run_eval_step()` calls `assembleForces_UFF` and that `fint` is sized by number of force pieces (`2*nbonds + 3*nangles + 4*ndihedrals + 4*ninversions`), not `natoms*4`.
2) **PBC shifts upload in pyOpenCL MMFF**: Confirm `pbc_shifts` buffer is populated and passed to kernels.
3) **Node ordering**: If pyOpenCL reorders nodes-first, apply the same permutation to C++ buffers (or disable reorder) before comparing.
4) **REQ mixing rule**: Confirm `_mixREQ` in CPU matches kernel path.

## D) Test cases to add/extend
- Clone existing UFF tests (`tests/tUFF/test_UFF_ocl.py`, `test_UFF_multi.py`) to:
  - Compare CPU vs C++/OCL vs pyOpenCL on H2O, HCOOH, pyrrole with 1-step and 10-step runs.
  - Dump buffers (topology + forces) to JSON/npz for diff.
- MMFF parity test akin to `test_MMFFsp3_pyOCL.py`:
  - Build with [MMFF_multi.py](cci:7://file:///home/prokop/git/FireCore/pyBall/MMFF_multi.py:0:0-0:0) (C++/OCL) and pyOpenCL builder, compare topology arrays and 1-step forces.
  - Include a perm map check if reorder is on.

## E) Visual debugging reuse (LiveViz2D / LivePortViz)
- [pyBall/XPBD_2D/XPBD_2D_utils.py](cci:7://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:0:0-0:0) and `pyBall/XPDB_AVBD/XPTB_utils.py` provide:
  - LiveViz2D canvas with mouse pick and per-frame callback.
  - Port visualization (local neighbor vectors).
- Reuse pattern:
  - Adapt the “ports” to show bonds/neighbors or pi-vectors: supply `pos`, `types`, `neighs`, optional `pi_dirs` as “ports”.
  - Hook simulation step to your MD loop: on each iteration, update positions/forces and redraw.
  - Keep mouse-pick to drag atoms and see response; add overlay to show per-atom force magnitude and neighbor indices.
- Minimal integration sketch:
  - Build a lightweight adapter that converts current MMFF/UFF state (`apos`, `neighs`, maybe `pipos`) into the expected [LiveViz2D](cci:2://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:30:0-178:29) data structure.
  - Use the existing `sim.upload_state()` and `sim.upload_topology()` style API; replace “stiffness/ports” with bond lines and pi arrows.

## F) Next steps I recommend
1) Wire a parity script:
   - Load system → build C++/OCL via [MMFF_multi.py](cci:7://file:///home/prokop/git/FireCore/pyBall/MMFF_multi.py:0:0-0:0) and pyOpenCL via `MMFF.py`/`UFF.py`.
   - Dump topology/params/positions to npz from both; compare with permutation if needed.
   - Run 1-step force eval on both; diff forces.
2) Patch pyOpenCL UFF to assemble forces (if not already).
3) Add gated debug-print utilities (shared format) on both sides.
4) Add a small LiveViz adapter to visualize a running pyOpenCL or C++/OCL stepper for manual inspection.

Tell me which piece you want first (UFF assemble fix, parity script scaffold, or LiveViz adapter), and I’ll outline concrete code changes accordingly.

---

## Gemini-3-Pro High

# Classical Forcefield Code Map & Parity Testing Strategy

This document maps the implementation variants of your forcefield system and provides a concrete strategy for parity testing, focusing on the interface between **C++ OpenCL** (accessed via `MMFF_multi`) and **pure pyOpenCL** (`pyBall/OCL`).

## 1. Module & Interface Map

### A. The Three Variants
1.  **CPU C++**:
    *   **Core**: [MolWorld_sp3](cci:2://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/cpp/common/molecular/MolWorld_sp3.h:76:0-3258:1) (owns `MMFFsp3_loc`, `UFF`, [NBFF](cci:2://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/cpp/common/molecular/NBFF.h:96:0-1312:1)).
    *   **Entry**: [MMFF_lib.cpp](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/cpp/libs/Molecular/MMFF_lib.cpp:0:0-0:0) (exposed as [pyBall.MMFF](cci:2://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/pyBall/OCL/MMFF.py:61:0-937:34)).
    *   **Role**: Reference "ground truth" for physics, but serial/OpenMP only.
2.  **GPU C++ Driver**:
    *   **Core**: `MolWorld_sp3_multi` (owns [OCL_MM](cci:2://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/cpp/common/OpenCL/OCL_MM.h:54:0-1517:1), `OCL_UFF`).
    *   **Entry**: [MMFFmulti_lib.cpp](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/cpp/libs_OCL/MMFFmulti_lib.cpp:0:0-0:0) (exposed as `pyBall.MMFF_multi`).
    *   **Role**: High-performance production code. Uses [relax_multi_mini.cl](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/cpp/common_resources/cl/relax_multi_mini.cl:0:0-0:0) and [UFF.cl](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/cpp/common_resources/cl/UFF.cl:0:0-0:0).
3.  **GPU Python Driver (pyOpenCL)**:
    *   **Core**: [pyBall.OCL.MolecularDynamics](cci:2://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/pyBall/OCL/MolecularDynamics.py:36:0-640:21) (MMFF), `pyBall.OCL.UFF` (UFF).
    *   **Entry**: Direct Python usage.
    *   **Role**: Rapid prototyping, easy introspection. Re-implements buffer packing and kernel management in Python.

### B. Interface Equivalence Table
This table maps how to access the *same data* across the two GPU interfaces we want to compare.

| Data Concept | **C++/OpenCL Wrapper** ([MMFF_multi.py](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/pyBall/MMFF_multi.py:0:0-0:0)) | **Pure pyOpenCL** (`pyBall.OCL.MMFF/UFF`) |
| :--- | :--- | :--- |
| **Wrapper Class** | [MMFF_multi.py](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/pyBall/MMFF_multi.py:0:0-0:0) calling [lib.init_buffers](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/pyBall/MMFF.py:917:0-927:26) | [MolecularDynamics](cci:2://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/pyBall/OCL/MolecularDynamics.py:36:0-640:21) / [UFF_CL](cci:2://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/pyBall/OCL/UFF.py:12:0-322:51) class |
| **Atom Positions** | `buffers['apos']` (host pointer to `W.nbmol.apos`) | `self.apos` (numpy array) |
| **Atom Forces** | `buffers['fapos']` | `self.aforce` (numpy array) |
| **Neighbors** | `ibuffers['neighs']` (host `W.ffl.neighs`) | `self.neighs` |
| **Bond Params** | `fbuffers['gpu_BLs']` / `gpu_BKs` | `self.bLs`, `self.bKs` |
| **Pi-Orbitals** | `buffers['pipos']` | `self.apos[natoms:natoms+nnode]` |
| **UFF Topology** | `ibuffers['bonAtoms']`, `ibuffers['angAtoms']` | `self.bonAtoms`, `self.angAtoms` |
| **UFF Params** | `buffers['bonParams']`, `buffers['angParams']` | `self.bonParams`, `self.angParams` |

## 2. Parity Testing Strategy

The goal is to verify `pyOpenCL_implementation == C++_OpenCL_implementation`.

### Strategy 1: Buffer Comparison (Static & Single-Step)
Use [MMFF_multi.py](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/pyBall/MMFF_multi.py:0:0-0:0) to expose the internal C++ state and compare it against `pyBall.OCL` objects.

**Workflow:**
1.  **Setup**: Load the same `.xyz/.mol2` in both systems.
2.  **Upload**: Call [pack_system](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/cpp/common/molecular/MolWorld_sp3_multi.h:587:0-675:1) (C++) and [upload_positions](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/pyBall/OCL/UFF.py:117:4-127:119) (pyOCL).
3.  **Inspect (Static)**: Compare topology and parameters **before** any run.
    *   *Check*: Neighbor lists (indices must match, watch out for sorting/permutations).
    *   *Check*: Bond lengths/stiffness (float tolerance).
4.  **Run 1 Step**: Run **forces only** (no move).
    *   C++: [W.run(1, ...)](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/pyBall/MMFF_multi.py:537:0-538:166) or [scan_uff_relaxed](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/tests/tUFF/test_UFF_multi.py:246:0-295:27) with 0 steps.
    *   Py: [run_eval_step()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/pyBall/OCL/UFF.py:157:4-168:29) or [run_step_basic()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/pyBall/OCL/MolecularDynamics.py:409:4-417:27) with `damp=0` and `dt=0`.
5.  **Download & Diff**:
    *   Compare `fapos` (forces). This is the most sensitive metric.
    *   If forces match, MD trajectories will likely match (divergence is expected later due to chaotic nature/float precision).

**Common Pitfalls:**
*   **Atom Sorting**: [pyBall.OCL.MMFF](cci:2://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/pyBall/OCL/MMFF.py:61:0-937:34) often reorders nodes to the front. [MolWorld_sp3](cci:2://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/cpp/common/molecular/MolWorld_sp3.h:76:0-3258:1) might not. You **must** apply the permutation `perm` to indices when comparing.
*   **Buffer Padding**: C++ often aligns data to multiples of 32/64. Python buffers might be tight. Compare valid regions only.

### Strategy 2: Precision Logging (The "Diff" Approach)
When buffers mismatch, use "printf debugging" in kernels.

**Guideline for Kernel Printfs:**
To make `diff` work, the output format must be **identical** and **atomic**.
*   **Format**: `[ID] [Tag] Val1 Val2 Val3 ...`
*   **Alignment**: Use `%10.5f` or `%20.12e` to ensure columns align.
*   **Gating**: Only print for `if(global_id == 0)` to prevent terminal flooding.

**Example [relax_multi_mini.cl](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/cpp/common_resources/cl/relax_multi_mini.cl:0:0-0:0) / [UFF.cl](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/cpp/common_resources/cl/UFF.cl:0:0-0:0) Snippet:**
```c
// Enable via #define DEBUG_PRINT in kernel header
if(get_global_id(0) == 0){ // Watch only atom 0
    printf("KERN_BOND i=%4d n=%2d f=(%+12.8f,%+12.8f,%+12.8f) E=%+12.8f\n", 
           iatom, n, f.x, f.y, f.z, E);
}
```

**Implementation in Python:**
Run both scripts, redirect stdout to files, and use `diff` or `meld`.
```bash
python test_C_plus_plus.py > log_cpp.txt
python test_pyOpenCL.py    > log_py.txt
meld log_cpp.txt log_py.txt
```

### Strategy 3: Visual Debugging (Reusing XPDB Tools)
The user wants to reuse `XPDB_AVBD` / `XPDB_2D` tools ([LiveViz2D](cci:2://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:30:0-178:29), [LivePortViz](cci:2://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/pyBall/XPDB_AVBD/XPTB_utils.py:414:0-467:29)). These tools are excellent for seeing **instability** (explosions, drift).

**Adapter Concept:**
We need a small adapter to translate MMFF data structures into what [LivePortViz](cci:2://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/pyBall/XPDB_AVBD/XPTB_utils.py:414:0-467:29) expects.

| Feature | MMFF (pyBall) | XPDB Viz (LivePortViz) | Conversion Logic |
| :--- | :--- | :--- | :--- |
| **Atoms** | `apos[:,:3]` | `pos` | Direct copy. |
| **Types** | `atypes` (int) | `elems` (str) | Map int ID to string (e.g. 6->"C"). |
| **Bonds** | `neighs` (adj list) | `bonds` (pairs) | Iterate `neighs`, if `j!=-1` and `i<j`, add `(i,j)`. |
| **Ports** | `pipos` (vectors) | `port_local` (rel pos) | `port_local = pipos * scale_factor`. |
| **Forces** | `aforce` | `force` | Direct copy (for quiver plot). |

**Proposed Python Adapter Class (Sketch):**

```python
from pyBall.XPDB_AVBD.XPTB_utils import plot_state_with_ports

class MMFFVizAdapter:
    def __init__(self, mmff_sys):
        self.sys = mmff_sys
        self.elems = [str(t) for t in mmff_sys.atypes] # Or look up real names
        
    def render(self):
        # 1. Extract Positions
        pos = self.sys.apos[:self.sys.natoms, :3]
        
        # 2. Extract Pi-Orbitals (Ports)
        # MMFF stores pi-vecs as atoms [natoms : natoms+nnode]
        pipos = self.sys.apos[self.sys.natoms : self.sys.natoms + self.sys.nnode, :3]
        
        # 3. Create fake "port neighbors" for visualization lines
        # This is just for drawing; we make the "port" look like a stick 1.0 unit away
        port_n = np.full(self.sys.natoms, 1 if self.sys.nnode > 0 else 0, dtype=np.int32)
        # XPDB expects (natoms, max_ports, 3). We have 1 pi-orb per node.
        pneigh = np.zeros((self.sys.natoms, 1, 3)) 
        pneigh[:self.sys.nnode, 0, :] = pipos * 1.0 # Scale vector for visibility
        
        # 4. Update Viz
        plot_state_with_ports(
            self.elems, 
            pos, 
            pneigh, 
            port_n, 
            force=self.sys.aforce[:self.sys.natoms], 
            title=f"Step {self.sys.step}"
        )
```

**Integration:**
Insert `adapter.render()` inside your Python simulation loop (e.g., every 10 steps).

## 3. Recommended Actions
1.  **Standardize Buffers**: Update [test_MMFFsp3_pyOCL.py](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/tests/tUFF/test_MMFFsp3_pyOCL.py:0:0-0:0) to use a `BUF_SPECS` dictionary similar to [test_UFF_ocl.py](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/tests/tUFF/test_UFF_ocl.py:0:0-0:0) for consistent inspection.
2.  **Fix MMFF_multi**: Ensure [MMFFmulti_lib.cpp](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/cpp/libs_OCL/MMFFmulti_lib.cpp:0:0-0:0) [init_buffers](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/pyBall/MMFF.py:917:0-927:26) covers *all* comparison targets (it currently misses some detailed UFF params if they aren't explicitly uncommented).
3.  **Create "Parity Runner"**: A single script that:
    *   Initializes both `MMFF_multi` (C++) and [pyBall.OCL.MMFF](cci:2://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/pyBall/OCL/MMFF.py:61:0-937:34) (Python).
    *   Steps them in lockstep.
    *   Asserts `norm(force_cpp - force_py) < tol`.
    *   Launches [LivePortViz](cci:2://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/pyBall/XPDB_AVBD/XPTB_utils.py:414:0-467:29) if the error exceeds tolerance, showing the specific atom with the discrepancy.

---

## Claude Opus 4.5

Now let me read the existing UFF parity tests to understand current coverage:

# Unified Parity-Testing Strategy: C++/OpenCL vs pyOpenCL

Based on my analysis of all the relevant files, here is a comprehensive strategy for parity testing between the **C++/OpenCL path** (accessed via [MMFF_multi.py](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/MMFF_multi.py:0:0-0:0) → [MMFFmulti_lib.cpp](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/libs_OCL/MMFFmulti_lib.cpp:0:0-0:0)) and the **pure pyOpenCL path** (`pyBall/OCL/*`).

## 1) Python Interface Summary

### C++/OpenCL Path (reference implementation)
| Python module | C++ library | Key functions |
|---|---|---|
| `@/home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/MMFF_multi.py` | [MMFFmulti_lib.cpp](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/libs_OCL/MMFFmulti_lib.cpp:0:0-0:0) | [init()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/MMFF_multi.py:426:0-443:248), [getBuffs_UFF()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/MMFF_multi.py:272:0-330:15), [run()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/MMFF_multi.py:537:0-538:166), [scan()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/libs_OCL/MMFFmulti_lib.cpp:279:0-328:1), [scan_relaxed()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/libs_OCL/MMFFmulti_lib.cpp:330:0-383:1), [pack_system()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/MMFF_multi.py:545:0-546:63), [upload()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/libs_OCL/MMFFmulti_lib.cpp:511:0-513:1), [download()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/libs_OCL/MMFFmulti_lib.cpp:514:0-516:1) |
| Buffers exposed | via [getBuff()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/MMFF_multi.py:196:0-200:48), [getIBuff()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/MMFF_multi.py:187:0-191:48), [getfBuff()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/MMFF_multi.py:214:0-218:49) | `apos`, `fapos`, `REQs`, `hneigh`, `fint`, `bonAtoms`, `angAtoms`, `dihAtoms`, `invAtoms`, `bonParams`, `angParams`, `dihParams`, `invParams`, `neighs`, `neighBs` |

### pyOpenCL Path (to be validated)
| Python module | Kernel file | Key classes/functions |
|---|---|---|
| `@/home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/UFF.py` | [UFF.cl](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common_resources/cl/UFF.cl:0:0-0:0) | [UFF_CL.toUFF()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/UFF.py:50:4-67:23), [realloc_buffers()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/UFF.py:69:4-115:87), [upload_positions()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/UFF.py:117:4-127:119), [upload_topology_params()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/UFF.py:129:4-155:82), [run_eval_step()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/UFF.py:157:4-168:29), [get_forces()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/UFF.py:170:4-179:56) |
| `@/home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/UFFbuilder.py` | — | [UFF_Builder.build()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/UFFbuilder.py:113:4-116:32), [assign_uff_types()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/UFFbuilder.py:118:4-175:34), [assign_uff_params()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/UFFbuilder.py:177:4-185:49), [get_arrays()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/UFFbuilder.py:187:4-204:23) |
| `@/home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/MolecularDynamics.py` | [relax_multi_mini.cl](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common_resources/cl/relax_multi_mini.cl:0:0-0:0) | [MolecularDynamics.pack_system()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/MolecularDynamics.py:169:4-207:122), [run_step_basic()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/MolecularDynamics.py:409:4-417:27), [run_step_rot()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/MolecularDynamics.py:419:4-427:27) |
| `@/home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/MMFF.py` | — | [MMFF.toMMFFsp3_loc()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/MMFF.py:254:4-562:31), [realloc()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3_multi.h:252:0-326:1), [make_back_neighs()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/MMFF.py:168:4-191:65) |

## 2) Existing Test Coverage

### UFF Tests
| Test file | What it tests | Gap |
|---|---|---|
| `@/home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/tests/tUFF/test_UFF_multi.py` | C++ CPU vs C++ GPU (via `MMFF_multi`) | Does **not** compare to pyOpenCL |
| `@/home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/tests/tUFF/test_UFF_ocl.py` | C++ CPU vs pyOpenCL UFF | Compares topology/params buffers, but [run_eval_step()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/UFF.py:157:4-168:29) doesn't call assembly kernel |
| `@/home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/tests/tUFF/test_MMFFsp3_pyOCL.py` | pyOpenCL MMFF MD dynamics | No direct comparison to C++ path |

### Gap: No unified test that runs **both** C++/OpenCL and pyOpenCL on the same molecule and compares forces/positions step-by-step.

## 3) Three-Level Parity Testing Strategy

### Level 1: Buffer Comparison (automated, scriptable)

**Goal**: After initialization and single-step evaluation, compare all internal buffers between implementations.

**Buffers to compare** (UFF example):

| Buffer | C++ accessor | pyOpenCL accessor | Notes |
|---|---|---|---|
| `bonAtoms` | `uff.bonAtoms` | `uff_cl.download_buf('bonAtoms')` | shape `(nbonds, 2)` |
| `angAtoms` | `uff.angAtoms` | `uff_cl.download_buf('angAtoms')` | C++ has 3 cols, GPU has 4 (pad with -1) |
| `dihAtoms` | `uff.dihAtoms` | `uff_cl.download_buf('dihAtoms')` | shape `(ndihedrals, 4)` |
| `invAtoms` | `uff.invAtoms` | `uff_cl.download_buf('invAtoms')` | shape `(ninversions, 4)` |
| `bonParams` | `uff.bonParams` | `uff_cl.download_buf('bonParams')` | shape `(nbonds, 2)` |
| `angParams` | `uff.angParams` | `uff_cl.download_buf('angParams1/2')` | C++ has 5 cols, GPU splits into 4+1 |
| `dihParams` | `uff.dihParams` | `uff_cl.download_buf('dihParams')` | shape `(ndihedrals, 3)` |
| `invParams` | `uff.invParams` | `uff_cl.download_buf('invParams')` | shape `(ninversions, 4)` |
| `neighs` | `uff.neighs` | `uff_cl.download_buf('neighs')` | shape `(natoms, 4)` |
| `neighBs` | `uff.neighBs` | `uff_cl.download_buf('neighBs')` | shape `(natoms, 4)` |
| `apos` | `uff.apos` | `uff_cl.download_buf('apos')` | shape `(natoms, 3)` |
| `fapos` | `uff.fapos` | [uff_cl.get_forces()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/UFF.py:170:4-179:56) | shape `(natoms, 3)` |
| `hneigh` | `uff.hneigh` | `uff_cl.download_buf('hneigh')` | shape `(natoms*4, 4)` |
| `fint` | `uff.fint` | `uff_cl.download_buf('fint')` | **size mismatch likely** |

**Comparison helper** (already exists in [test_UFF_ocl.py](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/tests/tUFF/test_UFF_ocl.py:0:0-0:0)):
```python
def compare_bufs(cpu_bufs, gpu_bufs, tol=1e-6, buf_type='Buffers'):
    # ... prints detailed diffs for any mismatched buffers
```

**Recommended test flow**:
```python
# 1. Initialize both paths with same molecule
mol_path = "cpp/common_resources/mol/formic_acid.mol2"

# C++/OpenCL path
import pyBall.MMFF_multi as uff_cpp
uff_cpp.init(nSys_=1, xyz_name=mol_path, bMMFF=True, bUFF=True, ...)
uff_cpp.getBuffs_UFF()

# pyOpenCL path
from pyBall.AtomicSystem import AtomicSystem
from pyBall.OCL import UFF as uff_ocl
mol = AtomicSystem(fname=mol_path)
uff_cl = uff_ocl.UFF_CL(nloc=32)
uff_cl.toUFF(mol)
uff_cl.upload_positions(mol.apos)
uff_cl.upload_topology_params(...)

# 2. Compare topology buffers BEFORE any force eval
cpu_topo = get_cpu_bufs(uff_cpp, TOPOLOGY_SPECS)
gpu_topo = get_gpu_bufs(uff_cl, TOPOLOGY_SPECS)
assert compare_bufs(cpu_topo, gpu_topo)

# 3. Run single eval step on both
uff_cpp.run(nstepMax=1, iParalel=2, ...)  # GPU
uff_cl.run_eval_step(bClearForce=True)

# 4. Compare forces
F_cpp = uff_cpp.fapos.copy()
F_ocl = uff_cl.get_forces()
max_diff = np.max(np.abs(F_cpp - F_ocl))
print(f"max |ΔF| = {max_diff:.3e}")
```

### Level 2: Debug Prints (gated, formatted for diff)

**Goal**: When buffer comparison fails, insert debug prints in both C++ and Python to trace execution flow.

**Gating strategy**:
```cpp
// C++ (in UFF.h or kernels)
extern int idebug;  // global, set from Python via setVerbosity()
#define DBG_ATOM 0  // which atom to trace
#define DBG_ITER 0  // which iteration

if(idebug > 0 && ia == DBG_ATOM && iter == DBG_ITER) {
    printf("UFF_BOND ia=%3d jn=%3d  r=%12.6f l0=%12.6f k=%12.6f  f=(%12.6f,%12.6f,%12.6f)\n",
           ia, jn, r, l0, k, f.x, f.y, f.z);
}
```

```python
# Python (in UFF.py or UFFbuilder.py)
DBG_ATOM = 0
DBG_ITER = 0

if self.debug and ia == DBG_ATOM and iter == DBG_ITER:
    print(f"UFF_BOND ia={ia:3d} jn={jn:3d}  r={r:12.6f} l0={l0:12.6f} k={k:12.6f}  f=({f[0]:12.6f},{f[1]:12.6f},{f[2]:12.6f})")
```

**Formatting rules** (for easy diff):
- Fixed-width columns: `%12.6f` for floats, `%3d` for ints
- All contributions from one atom on **single line**
- No line-length limit
- Prefix with component name: `UFF_BOND`, `UFF_ANGLE`, `UFF_DIH`, `UFF_INV`, `UFF_NB`

**Comparison workflow**:
```bash
# Run C++ path with debug
python test_UFF_multi.py --molecule formic_acid.mol2 --verbose 1 > cpp_debug.log 2>&1

# Run pyOpenCL path with debug
python test_UFF_ocl.py --molecule formic_acid.mol2 --verbose 1 > pyocl_debug.log 2>&1

# Diff
diff -u cpp_debug.log pyocl_debug.log | head -100
# Or use grep to filter specific atoms
grep "UFF_BOND ia=  0" cpp_debug.log > cpp_atom0.log
grep "UFF_BOND ia=  0" pyocl_debug.log > pyocl_atom0.log
diff cpp_atom0.log pyocl_atom0.log
```
### Level 3: Visual Debugging (interactive, user-driven)

**Goal**: Observe emergent behavior (momentum conservation, drift, non-physical dynamics) by animating the simulation.

**Existing tools to reuse**:

| Tool | Location | Features |
|---|---|---|
| [LiveViz2D](cci:2://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:30:0-178:29) | `@/home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/XPBD_2D/XPBD_2D_utils.py` | 2D matplotlib animation, mouse picking, bond/port visualization |
| [LivePortViz](cci:2://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/XPDB_AVBD/XPTB_utils.py:414:0-467:29) | `@/home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/XPDB_AVBD/XPTB_utils.py` | 3D matplotlib animation, mouse picking, port/force quivers |
| [attach_picker_2d](cci:1://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:388:0-434:15) | `XPBD_2D_utils.py:389` | Mouse drag to move atoms interactively |
| [attach_picker_3d](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/XPDB_AVBD/XPTB_utils.py:470:0-537:15) | `XPTB_utils.py:471` | 3D mouse picking with projection handling |

**Adaptation for UFF/MMFF**:

The existing tools are designed for XPBD (ports = bond directions). For UFF/MMFF:
- **Atoms** → same (positions, types, bonds)
- **Ports** → can be repurposed for:
  - **Pi-orbitals** (MMFF): `pipos` vectors
  - **Neighbor directions** (UFF): `hneigh` vectors
  - **Force vectors**: show `fapos` as quivers

**Minimal integration example**:
```python
from pyBall.XPBD_2D.XPBD_2D_utils import LiveViz2D, attach_picker_2d

# Create viz
viz = LiveViz2D(elems=mol.enames, view_scale=10.0)

# Simulation loop
for step in range(nsteps):
    # Run one step
    md.run_step_basic()
    
    # Download positions
    apos = md.download_results()['apos'][:, :2]  # 2D projection
    
    # Update viz
    viz.update(apos, neighs=mm.neighs, nnode=mm.nnode, 
               title=f"Step {step}", info=f"E={energy:.3f}")
    
    # Handle mouse picking
    if pick['active']:
        ia = pick['idx']
        # Pin atom to mouse position
        md.set_atom_pos(ia, pick['mouse'])
```

**Key callbacks to implement**:
- `on_step(step, apos, fapos, energy)` → update viz, check invariants
- `on_pick(ia, mouse_pos)` → pin atom, apply external force
- [on_release(ia)](cci:1://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:418:4-422:26) → unpin atom

## 4) Recommended Test Molecules (progressive complexity)

| Molecule | natoms | Features | File |
|---|---|---|---|
| H₂O | 3 | 2 bonds, 1 angle, no dihedrals | `H2O.xyz` |
| HCOOH | 5 | 4 bonds, 4 angles, 1 dihedral | `formic_acid.mol2` |
| Pyrrole | 10 | Ring, aromaticity | `pyrrol.xyz` |
| Cytosine | 13 | Heterocycle, H-bonding | `cytosine.xyz` |
| Guanine | 16 | Fused rings | `guanine.xyz` |

## 5) Known Gaps to Fix Before Parity Testing

### UFF pyOpenCL
1. **[run_eval_step()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/UFF.py:157:4-168:29) doesn't call `assembleForces_UFF`** → angle/dihedral/inversion forces missing from `fapos`
2. **`fint` buffer size** → allocated as `nA*4` but should be `2*nBonds + 3*nAngles + 4*nDihedrals + 4*nInversions`
3. **`a2f` mapping** → [mapAtomInteractions()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/UFF.py:243:4-322:51) builds arrays but assembly kernel not invoked

### MMFF pyOpenCL
1. **`pbc_shifts` not uploaded** in [MolecularDynamics.pack_system()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/OCL/MolecularDynamics.py:169:4-207:122)
2. **`excl` (exclusion table)** not populated for 1-3 exclusions
3. **`bkNeighs`** uploaded as -1 placeholders → recoil accumulation disabled

## 6) Proposed New Test Script Structure

```
tests/tUFF/
├── test_UFF_parity.py          # NEW: unified C++/OpenCL vs pyOpenCL comparison
├── test_UFF_multi.py           # existing: C++ CPU vs GPU
├── test_UFF_ocl.py             # existing: C++ vs pyOpenCL (partial)
├── test_MMFFsp3_pyOCL.py       # existing: pyOpenCL MMFF dynamics
└── test_MMFFsp3_parity.py      # NEW: unified MMFF comparison
```

**`test_UFF_parity.py` outline**:
```python
def test_topology_parity(mol_path):
    """Compare topology buffers after builder."""
    
def test_single_eval_parity(mol_path, components=['bonds','angles','dihedrals','inversions']):
    """Compare forces after single evaluation."""
    
def test_relaxation_parity(mol_path, nsteps=100):
    """Compare trajectories after relaxation."""
    
def test_with_debug_prints(mol_path, atom_idx=0):
    """Enable debug prints and save logs for diff."""
```
## Summary

| Strategy | When to use | Tools |
|---|---|---|
| **Buffer comparison** | First pass; automated CI | [compare_bufs()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/tests/tUFF/test_UFF_ocl.py:101:0-153:19), [getBuff()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/MMFF_multi.py:196:0-200:48)/`download_buf()` |
| **Debug prints** | When buffers differ; trace execution | Gated `printf`/[print](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/cpp/common/molecular/MolWorld_sp3_multi.h:92:4-92:103), `diff`, `grep` |
| **Visual debugging** | Emergent behavior; user-driven | [LiveViz2D](cci:2://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:30:0-178:29), [LivePortViz](cci:2://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-ce550402/pyBall/XPDB_AVBD/XPTB_utils.py:414:0-467:29), mouse picking |

The existing test infrastructure in `tests/tUFF/` provides a good foundation. The main work needed is:
1. Fix the pyOpenCL UFF assembly path
2. Create a unified parity test that runs both implementations
3. Add gated debug prints with consistent formatting


