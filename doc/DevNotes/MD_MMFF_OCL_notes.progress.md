
# Notes: MD with MMFFsp3 in OpenCL


Goal: Run MD for molecules using MMFFsp3 force field via pyOpenCL (`pyBall/OCL/`). Implement and expose two new propagator kernels: `updateAtomsMMFFf4_rot` (rotational PI) and `updateAtomsMMFFf4_RATTLE` (unit constraint for PI)

**CRITICAL:** Kernel and resource paths are **relative to the working directory**. Do **not** change source paths in code; instead run scripts from the proper folder (e.g., `tests/tUFF/` for UFF/MMFF tests) so `common_resources/cl/*.cl` are found. Running from the wrong CWD causes “Could not open kernel source file”.

Related to script `/tests/tUFF/test_MMFFsp3_pyOCL.py`

###  Relevant files/functions
- `cpp/common_resources/cl/relax_multi_mini.cl`
  - `getMMFFf4()`: bonded forces + pi contributions
  - `updateAtomsMMFFf4()`: existing leap-frog
  - `updateAtomsMMFFf4_rot()`: rotational update for pi-orbitals (new)
  - `updateAtomsMMFFf4_RATTLE()`: RATTLE-like pi constraint (new)
  - `runMD()`: calls `getMMFFf4` + update (not yet switched to variants)
- `pyBall/OCL/MolecularDynamics.py`
  - `allocate_cl_buffers()`: OpenCL buffers
  - `pack_system()`: upload per-system buffers; now also uploads bkNeighs
  - `setup_kernels()`: added arg gens for new kernels
  - `run_updateAtomsMMFFf4_{rot,RATTLE}()`: wrappers
- `pyBall/OCL/MMFF.py`
  - `MMFF.toMMFFsp3_loc()`: builds MMFF arrays from `AtomicSystem`
  - Now writes pi-orbital orientations into `apos[natoms:natoms+nnode]`
  - Builds placeholder `back_neighs` (atom-to-atom) and we upload -1s to GPU for `bkNeighs`
- `pyBall/AtomicSystem.py`
  - Loads `.mol2` (`au.loadMol2`) and builds bonds, neighbors

###  Important offsets and 2D-NDRange
- Kernels compute per-system base indices: `i0a=iS*natoms`, `i0v=iS*(natoms+nnode)`
- Host uploads must match: byte offsets applied consistently for each buffer slice per system (REQs, neighs, apars, ...)
- Work sizes: `sz_na=(roundup(natoms,nloc), nSystems)`, `sz_nvec=(roundup(nvecs,nloc), nSystems)`

###  PI-orbital placement
- Kernels use `apos[iav + nAtoms]` as pi direction. We now fill `MMFF.apos[natoms : natoms+nnode, :3] = pipos`
- `MolecularDynamics.pack_system()` uploads full `nvecs` slice; therefore pi vectors are on device

###  Back-neighbor forces
- Kernels accumulate recoil forces via `bkNeighs` indices addressing `fneigh[]`
- Python now mirrors `MMFFsp3_loc::makeBackNeighs()`:
  - `pyBall/OCL/MMFF.py::make_back_neighs()` packs node bonds as `ia*4+ib` and assigns capping atom neighbors.
  - `pyBall/OCL/MolecularDynamics.py::pack_system()` uploads this packed table into the `bkNeighs` buffer.
  - Issue observed: without valid back-neighbors hydrogen caps drifted away. After populating them, trajectories remained stable.
- `tests/tUFF/test_MD_OCL_formic_acid.py --print-params 1` now dumps `neighs`, `back_neighs`, and bonded parameters to verify correct setup before GPU upload.

###  Plan / Next steps
- Hook runMD to choose propagator variant (rot or RATTLE) via a mode flag or separate run methods
- Implement proper `bkNeighs` mapping (CPU indices to `fneigh` layout) if needed
- Add thermostats/switching in Python control as needed
- Validate on small molecules and cross-check with CPU baseline

###  Invariant monitoring (`tests/tUFF/test_MMFFsp3_pyOCL.py`)
- Added COM/momentum diagnostics and energy tracking per MD step with CLI toggles `--monitor`, `--monitor-props`, `--monitor-plot`, `--save-monitor`, `--save-monitor-data`.
- `elements.index_mass` now supplies non-zero masses for kinetic/COM evaluation; zero or missing masses raise explicitly.
- Stored history supports optional PNG/NPZ export plus drift summary for `F`, `T`, `P`, `L`, and energies.
- Plot helper renders each requested invariant (vector components split per axis) to ease conservation checks without extra scripting.

###  VdW subtraction toggle (`bSubtractVdW`)
- `getMMFFf4()` subtracts Lennard-Jones/Coulomb overlap for bonded neighbors when `bSubtractVdW=1`.
- Rule: enable subtraction only when the non-bonded kernel runs (`--do-nb 1` / `subtract_vdw=1`). Disable it when skipping non-bonded forces, otherwise we remove interactions that were never added, leading to NaNs.
- `tests/tUFF/test_MD_OCL_formic_acid.py` exposes `--subtract-vdw` so CLI configs clearly reflect this coupling.
- During debugging, NaNs reproduced when `--do-nb 0` with subtraction still on; turning it off restored stability.

### Debugging helpers / learnings

- `--print-params 1` gives early visibility into neighbor mappings without needing GPU downloads.
- Label overlays (`--plot-labels {number,type}`) help identify which atoms deviate in trajectory plots.


### 2025-09-30 session notes  

- **[Pack-system diagnostics]** Added `--print-params` passthrough in `tests/tUFF/test_MMFFsp3_pyOCL.py` so that `pyBall/MD_test_utils.configure_md()` calls `MolecularDynamics.set_pack_system_debug(True)`, exposing `pack_system()` dumps before GPU upload for fast inspection of `neighs`, `bkNeighs`, `Ksp`, etc.
- **[Recoil force assembly]** Confirmed that valid `bkNeighs` are mandatory; when we forget to upload them the hydrogen caps drift because recoil forces are never accumulated back. The current workaround updates `pyBall/OCL/MMFF.make_back_neighs()` and `MolecularDynamics.pack_system()` to mirror `MMFFsp3_loc::makeBackNeighs()`, but we still owe a follow-up to assemble the OpenCL recoil forces on-device instead of relying on host verification logs.
- **[Methanol Ksp/Kpp check]** Discovered that `MMFF.toMMFFsp3_loc()` copies `AtomTypes.dat` values even for sp³ heteroatoms, so methanol’s `O_3` produced non-zero `Ksp`. For saturated systems we now override these entries to zero in `tests/tmp/data/AtomTypes.dat` (setting both `Ksp` and `Kpp` to 0.0) before building the MMFF tables.
- **[Validation]** After zeroing the `Ksp/Kpp` override the methanol test (`tests/tUFF/test_MMFFsp3_pyOCL.py --steps 200`) keeps angular momentum and forces at machine zero, confirming the parameter tweak resolves the spurious torques.

### 2025-10-02 session notes  

- **[Distorted starts for Pi testing]** `tests/tUFF/test_MMFFsp3_pyOCL.py` now exposes `--distort` and `--distort-seed` so we can inject random Cartesian and pi-vector perturbations before MD. This is essential to bend the pi plane and trigger pi-dependent forces without editing the input geometry by hand.
- **[Recoil buffer inspections]** Added `--dump-fneigh` which downloads `fneigh` after a run via `fetch_arrays()` (forcing `dtype=np.float32`). Dumps print the atom half and the pi half separately, matching the `nnode*4` offset used on the device.
- **[Pi back-neighbor mapping]** `pyBall/OCL/MMFF.make_back_neighs()` now fills `back_neighs` for both atoms and the corresponding pi orbitals. For each sigma bond we store the packed index `ia*4+ib` for the atom recoils and `ia*4+ib + nnode*4` for the pi slot (`natoms + ja`). `MolecularDynamics.pack_system()` uploads these tables unchanged, so `bkNeighs` finally points at the correct half of `fneigh`.
- **[Sigma-pi torque check]** After fixing the Rodrigues update in `updateAtomsMMFFf4_rot` and providing valid pi back-neighbors, sigma-pi alignment terms now conserve angular momentum in the rotational integrator runs.
- **[Pi-pi still problematic]** Pi-pi alignment continues to leak angular momentum. The torque imbalance tracks the pi section of `fneigh`, so the remaining defect is likely in how pi recoil torques are accumulated and read back (the indices are correct now, but magnitudes look off).
- **[Two-node baseline parity]** Matched `tests/tUFF/test_MMFFsp3_pyOCL.py` two-node builder to `doc/py/pi_dynamics/pipi_dynamics.py`: atoms start separated by `bL0`, pi vectors normalized, `Ksp/Kpp` CLI-tunable, and `apars[:,3]` forced to `0` so sigma–pi equilibrium cosines align. Upload logic now calls `md.pack_system(0, mm)` after overrides to push parameter edits to the GPU.
- **[Parameter verification]** With `--print-params 1`, first-step GPU prints show `Ks/Kp = 1` and correct node orientations. The standalone driver and OpenCL kernel produce identical pi-pi and pi-sigma forces for step 0 when fed the synchronized geometry.
- **[Recoil accumulation gap]** Despite matching instantaneous torques, GPU propagation doubles the torque on pi node B (`0.24010 → 0.48021`). Analysis points to recoil aggregation via `fneigh`+`bkNeighs`: pi orbitals receive both direct torque and neighbor recoil, while sigma atoms may consume pi recoil components. Need to audit `bkNeighs` indexing and ensure sigma entries never read pi slots.
- **[Diagnostics and tooling]** Added heavy debug prints in `relax_multi_mini.cl` (`getMMFFf4_rot` and `updateAtomsMMFFf4_rot`) plus in `doc/py/pi_dynamics/pipi_dynamics.cl` to log torques, bond frames, and recoil vectors each step. Both test harnesses now default to `--steps 5` for manageable trace lengths.
- **[Next actions]** Dump `bkNeighs` and `fneigh` immediately after `pack_system()` to confirm layout (`sigma` block first, `pi` block offset by `nnode*4`). Consider disabling recoil addition in propagation temporarily to isolate the defect. Once indices verified, adjust host packing or kernel readout so pi recoil mirrors the standalone driver without double counting.

### 2025-10-03 session notes  

- **[Pi-pi angular momentum fix]** Resolved the torque leak by dropping pi-orbital recoil accumulation (`fneigh` reads) and instead evaluating pi-pi alignment twice (explicit `i→j` and `j→i` contributions). This mirrors the standalone driver and keeps angular momentum conserved without relying on pi back-neighbor bookkeeping.
- **[Symmetric mixing audit]** Verified that `pyBall/OCL/MMFF.toMMFFsp3_loc()` derives bond parameters symmetrically: `bLs` use `ti.Ruff + tj.Ruff`, bond stiffness `bKs` take `sqrt(ti.Kss * tj.Kss)`, and pi-pi stiffness `Kpp` similarly applies `sqrt(atom_type.Kpp * jtyp.Kpp)`. Each node loop therefore writes identical values into both `i-j` and `j-i` slots, ensuring the explicit pair evaluation above remains balanced.

## Current parity test status (UFF/MMFF)

### Summary matrix (6-molecule sweep: methanol.mol2, HCONH2.xyz, uracil.xyz, xylitol.mol2, guanine.xyz, Si10_H.xyz)

| Forcefield | Driver pair | Scope | Result | Notes |
| --- | --- | --- | --- | --- |
| MMFF | C++ CPU ↔ C++ OpenCL | bonded-only | **PASS (all 6)** | `test_MMFF_multi_parity.py` with `--nonbond 0 --angles 1 --pisigma 1 --pipii 1`; tolerances 1e-4/0.0. |
| MMFF | C++ OpenCL ↔ PyOpenCL | bonded-only | **PASS (all 6)** | `test_MMFF_cpp_vs_pyocl_parity.py` with `--nonbond 0 --angles 1 --pisigma 1 --pipii 1`; tolerances 1e-4/0.0. |
| UFF | C++ CPU ↔ C++ OpenCL | bonds-only | **PASS (all 6)** | `test_parity_suite.py` stages `--components bonds`; max |ΔF| ~1e-6. |
| UFF | C++ CPU ↔ C++ OpenCL | bonds+angles | **PASS (all 6)** | Stage `--components bonds,angles`; float noise only. |
| UFF | C++ CPU ↔ C++ OpenCL | bonds+angles+dihedrals | **FAIL (methanol, HCONH2, uracil, xylitol, guanine); PASS (Si10_H)** | Mismatch isolated to dihedral term; max |ΔF| up to ~5e-3 for methanol, larger for others earlier; likely `dihNgs/hneigh` sign/assembly issue. |
| UFF | C++ CPU ↔ C++ OpenCL | full (incl. inversions) | **same as above** | Si10_H trivial/pass; others fail once dihedrals enabled. |

### Key scripts and roles (with dump/fast-exit options)
- `tests/tUFF/test_parity_suite.py`: Orchestrates multi-molecule runs, stages UFF components, reruns failures with `--dump`/`--dump-n`, uses `--fast-exit` to avoid C++ teardown aborts; logs in `OUT_parity_suite/`.
- `tests/tUFF/run_parity_uff.sh`: Rebuilds libMMFFmulti_lib then runs UFF PyOCL vs C++ OCL parity across 6 molecules and staged components (bonds → bonds+angles → +dihedrals → +inversions); prints SUMMARY lines and exits nonzero on failure.
- `tests/tUFF/test_MMFF_multi_parity.py`: C++ CPU vs C++ OpenCL MMFF bonded parity; flags `--dump/--dump-n/--save-npz/--fast-exit`.
- `tests/tUFF/test_MMFF_cpp_vs_pyocl_parity.py`: C++ OpenCL vs PyOpenCL MMFF parity; same dump flags.
- `tests/tUFF/test_MMFF_ocl_parity.py`: C++ CPU vs PyOpenCL MMFF parity; same dump flags (named ocl but compares CPU↔PyOCL).
- `tests/tUFF/run_parity_mmff.sh`: Rebuilds libMMFFmulti_lib then runs MMFF PyOCL vs C++ OCL parity over 6 molecules and 3 presets (bonds; bonds+angles+pi; bonds+angles+pi+nonbond); prints SUMMARY lines and exits nonzero on failure.
- `tests/tUFF/test_UFF_ocl.py`: C++ CPU vs C++ OpenCL UFF; supports `--components` staging plus `--dump/--dump-n/--fast-exit`; now uses `uff_cpp.eval()` to keep geometry fixed.
- `tests/tUFF/test_UFF_multi.py`: Legacy C++ UFF multi-system test (CPU vs GPU) — kept as reference.
- `tests/tUFF/test_MMFFsp3_pyOCL.py`: Rotational/pi-focused MMFF PyOCL harness with pack/dump options (not part of suite).

### Insights / open items
- MMFF parity achieved; all inputs and forces match across CPU/OCL/PyOCL for bonded terms on the 6-molecule set.
- UFF parity clean until dihedrals; staging proves bonds/angles OK. Next step: compare/dump `dihNgs`, `hneigh`, and `a2f_*` to find the sign/assembly mismatch in `evalDihedrals_UFF` vs CPU `bakeDihedralNeighs()` convention.
- Tolerances: MMFF 1e-4 (forces), UFF staged runs currently 1e-4; dihedral discrepancy exceeds float noise, so do not loosen tolerance until the sign issue is fixed.

### 2026-02-11 UFF PyOpenCL parity recap

- **Root causes fixed**
  - CPU was still running nonbonded/subtraction paths because `setSwitches2/setSwitchesUFF` treat `0` as “keep”. Added C++ API `setSwitchesUFF_NB` and Python binding to force-disable UFF internal `bNonBonded/bNonBondNeighs/bSubtractAngleNonBond` with `-1` for bonds-only parity.
  - GPU `angParams` upload was wrong: CPU layout `[K,c0,c1,c2,c3]` vs kernel expectation `angParams1=[c0..c3], angParams2_w=[K]`. Fixed upload/recompose; angles parity then matched.
  - Minor assemble debug OOB in `evalBondsAndHNeigh_UFF` fixed by indexing `bonParams` via `neighBs` instead of neighbor atom id.
  - `assembleForces_UFF` a2f base offset corrected for single-system map (a2f_indices packed by ref-count, not nf_per_system).

- **Results**
  - Bonds-only max |ΔF| ~4e-6.
  - Angles-only max |ΔF| ~5e-6 after angParams fix.
  - Dihedrals-only and inversions-only within 1e-5–4e-5 (float32 vs CPU float64 trig), default tol set to 5e-5.

- **Debugging strategy / insistence**
  - Always instrument kernels with deterministic `DBG_UFF` prints (default **on** in kernel sources because C++ OpenCL can’t pass build flags). Use global-id 0/IDBG_SYS gating to prove which kernels run.
  - Exact single-pass eval on both CPU/GPU; no scans/loops unless explicitly requested.
  - Disable nonbonded and subtraction paths when comparing bonded terms; match switches explicitly (use `-1` to force off).
  - Print and compare buffers (topology/params) before force comparison; when in doubt, dump per-atom/bond in kernels.
  - Accept float32 GPU vs float64 CPU differences by setting tolerance accordingly; don’t chase noise once physics and buffers align.

- **[Follow-up]** Keep Kpp/Ks overrides synchronized with AtomTypes overrides (e.g. methanol tweak) now that recoil is disabled; any future mixing changes must maintain symmetry so the double-counted pi-pi evaluation does not introduce bias.

### 2026-02-11 UFF PyOCL vs C++ OCL (6-molecule sweep, staged components)

| molecule      | bonds | bonds+angles | +dihedrals | +dihedrals+inversions |
|---------------|-------|--------------|------------|-----------------------|
| methanol.mol2 | PASS  | PASS         | PASS       | PASS                  |
| HCONH2.xyz    | PASS  | PASS         | PASS       | PASS                  |
| uracil.xyz    | PASS  | PASS         | PASS       | PASS                  |
| xylitol.mol2  | PASS  | PASS         | PASS       | PASS                  |
| guanine.xyz   | PASS  | PASS         | PASS       | PASS                  |
| Si10_H.xyz    | PASS  | PASS         | PASS       | PASS                  |

Notes
- Root cause fixed: PyOCL packed `dihParams` as float3 and `invNgs` as int3; kernel expects float4/int4. PyOCL now pads (dihParams w=0, invNgs w=-1), kernel uses `__global int4* invNgs`.
- C++ harness already allocated `invNgs` as int4, so only kernel signature + PyOCL upload needed change.
- DBG_UFF defaults to 2 in UFF.cl with IDBG gating; parity confirmed across all components on the 6-molecule sweep.

## MMFF GPU/CPU Parity Findings (bonded-only, Feb 2026)

- Flow (mirror run_ocl_opt for forces-only): cleanF -> getMMFFf4 -> updateAtomsMMFFf4 with dt=0 -> download/unpack. The dt=0 call assembles recoil forces from fneigh via bkNeighs but skips dynamics (see relax_multi_mini.cl: updateAtomsMMFFf4 early return at MDpars.x<1e-16).
- Switch propagation: setSwitches2 touches W.ffl; propagate flags to W.ffls[isys] before CPU eval so CPU/GPU match.
- VdW subtraction flag: set bSubtractVdW=0 when NonBonded is off and pass it to kernel arg 16 of getMMFFf4.
- Force readback: both CPU and GPU paths must copy ffls[isys].fapos into nbmol.fapos so Python sees the same buffer.
- bkNeighs upload sizing: ibuff_bkNeighs/bkNeighs_new are allocated as nSystems*nvecs (nvecs=natoms+nnode). Upload must use size=nvecs and offset=i0v. Uploading nbkng overruns the buffer (nbkng=nnode*4*2). In the current builder, nnode=natoms and each atom gets a pi slot; caps have Ksp/Kpp=0 so their pi entries are inert but still occupy nvecs slots.
- Parity result: HCOOH bonded-only (nonbond off) now matches CPU vs GPU with Max|dF|~1.9e-6, RMS~7e-7, energy match.

### 2026-02-11 MMFF C++/OCL vs PyOCL parity status
- Inputs: **PASS** — neighs, bkNeighs, REQs, MMpars/apars, BLs, BKs, Ksp, Kpp, and `gpu_atoms` (xyz + pi tail) match within tolerance.
- Forces: **PENDING** — bonded-only HCOOH still shows |ΔF|≈3.7 in the parity harness; C++ `gpu_aforces.w` is zero (kernels don’t accumulate energy), PyOCL sums nonzero `.w`.
- Kernels: both paths build from the same `cpp/common_resources/cl/relax_multi_mini.cl`; intended sequence is `cleanForceMMFFf4 -> getMMFFf4 -> updateAtomsMMFFf4(dt=0)`. Need to re-check flags (nPBC, bSubtractVdW, mask) and buffer clearing (aforce/fneigh) to ensure identical flow.
- Next checks: compare forces only (ignore `.w`), gate kernel DBG prints for one atom/system after `getMMFFf4` and after `updateAtomsMMFFf4`, verify work sizes and per-system offsets, and decide whether to add matching energy accumulation or skip energy comparison.