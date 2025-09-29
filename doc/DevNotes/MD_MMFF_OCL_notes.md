---
# Notes: MD with MMFFsp3 in OpenCL
---

Goal: Run MD for molecules using MMFFsp3 force field via pyOpenCL (`pyBall/OCL/`). Implement and expose two new propagator kernels: `updateAtomsMMFFf4_rot` (rotational PI) and `updateAtomsMMFFf4_RATTLE` (unit constraint for PI)

Relevant files/functions
- cpp/common_resources/cl/relax_multi_mini.cl
  - getMMFFf4(): bonded forces + pi contributions
  - updateAtomsMMFFf4(): existing leap-frog
  - updateAtomsMMFFf4_rot(): rotational update for pi-orbitals (new)
  - updateAtomsMMFFf4_RATTLE(): RATTLE-like pi constraint (new)
  - runMD(): calls getMMFFf4 + update (not yet switched to variants)
- pyBall/OCL/MolecularDynamics.py
  - allocate_cl_buffers(): OpenCL buffers
  - pack_system(): upload per-system buffers; now also uploads bkNeighs
  - setup_kernels(): added arg gens for new kernels
  - run_updateAtomsMMFFf4_{rot,RATTLE}(): wrappers
- pyBall/OCL/MMFF.py
  - MMFF.toMMFFsp3_loc(): builds MMFF arrays from `AtomicSystem`
  - Now writes pi-orbital orientations into `apos[natoms:natoms+nnode]`
  - Builds placeholder `back_neighs` (atom-to-atom) and we upload -1s to GPU for `bkNeighs`
- pyBall/AtomicSystem.py
  - Loads `.mol2` (`au.loadMol2`) and builds bonds, neighbors

Important offsets and 2D-NDRange
- Kernels compute per-system base indices: `i0a=iS*natoms`, `i0v=iS*(natoms+nnode)`
- Host uploads must match: byte offsets applied consistently for each buffer slice per system (REQs, neighs, apars, ...)
- Work sizes: `sz_na=(roundup(natoms,nloc), nSystems)`, `sz_nvec=(roundup(nvecs,nloc), nSystems)`

PI-orbital placement
- Kernels use `apos[iav + nAtoms]` as pi direction. We now fill `MMFF.apos[natoms : natoms+nnode, :3] = pipos`
- `MolecularDynamics.pack_system()` uploads full `nvecs` slice; therefore pi vectors are on device

Back-neighbor forces
- Kernels accumulate recoil forces via `bkNeighs` indices addressing `fneigh[]`
- Python now mirrors `MMFFsp3_loc::makeBackNeighs()`:
  - `pyBall/OCL/MMFF.py::make_back_neighs()` packs node bonds as `ia*4+ib` and assigns capping atom neighbors.
  - `pyBall/OCL/MolecularDynamics.py::pack_system()` uploads this packed table into the `bkNeighs` buffer.
  - Issue observed: without valid back-neighbors hydrogen caps drifted away. After populating them, trajectories remained stable.
- `tests/tUFF/test_MD_OCL_formic_acid.py --print-params 1` now dumps `neighs`, `back_neighs`, and bonded parameters to verify correct setup before GPU upload.

Plan / Next steps
- Hook runMD to choose propagator variant (rot or RATTLE) via a mode flag or separate run methods
- Implement proper `bkNeighs` mapping (CPU indices to `fneigh` layout) if needed
- Add thermostats/switching in Python control as needed
- Validate on small molecules and cross-check with CPU baseline

VdW subtraction toggle (`bSubtractVdW`)
- `getMMFFf4()` subtracts Lennard-Jones/Coulomb overlap for bonded neighbors when `bSubtractVdW=1`.
- Rule: enable subtraction only when the non-bonded kernel runs (`--do-nb 1` / `subtract_vdw=1`). Disable it when skipping non-bonded forces, otherwise we remove interactions that were never added, leading to NaNs.
- `tests/tUFF/test_MD_OCL_formic_acid.py` exposes `--subtract-vdw` so CLI configs clearly reflect this coupling.
- During debugging, NaNs reproduced when `--do-nb 0` with subtraction still on; turning it off restored stability.

Debugging helpers / learnings
- `--print-params 1` gives early visibility into neighbor mappings without needing GPU downloads.
- Label overlays (`--plot-labels {number,type}`) help identify which atoms deviate in trajectory plots.
