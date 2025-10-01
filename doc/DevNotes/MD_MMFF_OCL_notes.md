
# Notes: MD with MMFFsp3 in OpenCL


Goal: Run MD for molecules using MMFFsp3 force field via pyOpenCL (`pyBall/OCL/`). Implement and expose two new propagator kernels: `updateAtomsMMFFf4_rot` (rotational PI) and `updateAtomsMMFFf4_RATTLE` (unit constraint for PI)

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

### 2025-10-01 session notes
- **[Rotational kernels working]** Verified `getMMFFf4_rot()` and `updateAtomsMMFFf4_rot()` in `cpp/common_resources/cl/relax_multi_mini.cl` conserve angular momentum when driven from `tests/tUFF/test_MMFFsp3_pyOCL.py`; added targeted GPU debug prints (`GPU cleanForceMMFFf4`, etc.) during bring-up.
- **[Python driver cleanup]** Simplified the test harness to call `MolecularDynamics.run_step_rot()` / `run_step_basic()` directly, removed the `finalize_recording()` helper, and inlined trajectory export/plotting with optional XYZ dumps plus plot suptitles indicating molecule + dynamics mode.
- **[XYZ/plot exports]** Introduced `write_xyz_trajectory()` in `pyBall/MD_test_utils.py` and wired CLI flags so runs optionally emit `trj.npy`, `trj.xyz`, and PNG plots (`--save-plot`, `--save-xyz`).
- **[Pending energy fix]** Noted residual total-energy drift traceable to missing rotational kinetic term for pi orbitals; follow-up implements angular KE contribution inside `pyBall/MD_test_utils.compute_energies()`