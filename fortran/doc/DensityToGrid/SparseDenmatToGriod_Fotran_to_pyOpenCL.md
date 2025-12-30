# Sparse Density Matrix to Grid (Fortran → pyOpenCL)

## Principles
- Keep Fortran and Python arrays in the native Fortran-style shapes; **do not flatten** on host. Only the OpenCL kernel views the buffers as flat pointers and uses explicit stride/index arithmetic.
- Export neighbor and dimension metadata separately from matrix blocks to avoid overstuffed interfaces and to enable reuse across H/S/rho.
- Preserve sparsity: rho/H/S are stored per (iatom, ineigh) blocks sized `(numorb_max, numorb_max)`; neighbors are per-atom lists (`neighn`, `neigh_j`, `neigh_b`, `xl`).
- Load balance GPU: build a task list of active grid blocks; sort by estimated pair count; skip empty blocks.

### Python array ordering note
- We now allocate **C-contiguous** NumPy arrays on the Python side (default order). Fortran will fill them using its column-major expectations; this effectively reverses index order relative to Fortran declarations. When reading in Python, treat the **last index as fastest** and keep the index order consistent (e.g., Fortran `A(i,j,k,l)` corresponds to Python linear order `A[l,k,j,i]`). This avoids ctypes `CONTIGUOUS` flag errors while keeping the data layout predictable. Do **not** set `order='F'` in Python allocations.

## Fortran exports (`fortran/MAIN/libFireCore.f90`)
- `firecore_get_HS_neighs(...)` — outputs dims/meta + neighbor arrays: `neighn`, `neigh_j`, `neigh_b`, `xl`, orbital meta (`num_orb`, `degelec`, `iatyp`, `lssh`, `mu`, `nu`, `mvalue`, `nssh`, `nzx`).
- `firecore_get_HS_sparse(...)` — outputs H/S blocks only; shape `(numorb_max, numorb_max, neigh_max, natoms)`.
- `firecore_get_rho_sparse(...)` — outputs rho blocks in the same shape as H/S.
- Allocation reminders (existing in Fortran):
  - `rho(numorb_max, numorb_max, neigh_max, natoms)`
  - `h_mat(numorb_max, numorb_max, neigh_max, natoms)`
  - `s_mat(numorb_max, numorb_max, neigh_max, natoms)`
  - `neighn(natoms)`, `neigh_j(neigh_max, natoms)`, `neigh_b(neigh_max, natoms)`, `xl(3, neigh_max, natoms)`

## Python bindings (`pyBall/FireCore.py`)
- Add C bindings/wrappers for the three exports above.
- Return NumPy arrays with the same n-d layout/ordering as Fortran; **no flattening in Python**.
- Provide a convenience `FireballDims` struct (natoms, neigh_max, numorb_max, per-atom `num_orb`, `degelec`, etc.).

## Host plan (pyBall/FireballOCL/Grid.py)
- Inputs: atom positions/types, per-atom Rcut, neighbor arrays, rho blocks (native shape), grid spec (origin, dA/dB/dC, ngrid).
- Build macro-grid (~max Rcut); map atoms → macro-cells (each atom touches ≤4 in 2D / ≤8 in 3D).
- For each fine block inside occupied macro-cells:
  - Collect overlapping atoms (sphere-AABB).
  - Count overlapping rho pairs using precomputed adjacency (|rij| < Ri+Rj).
  - Append active task `{block origin, atom index list/slice, pair_count}`.
- Optionally sort tasks by `pair_count` to reduce tail effects; skip empty blocks.
- Upload atoms, neighbors, rho blocks, and task list to GPU; pass buffers as-is (contiguous in Fortran order) to OpenCL.

## Kernel plan (pyBall/FireballOCL/cl/Grid.cl)
- Buffers received as flat pointers but originally Fortran-ordered; kernel uses explicit strides:
  - Element index: `idx = ((iatom * neigh_max + ineigh) * numorb_max + imu) * numorb_max + inu` for rho/H/S.
- Per task (grid block):
  - Optionally build local pair list in shared memory from rho block sparsity.
  - One thread per voxel; iterate pair list, accumulate `rho_ij * phi_i * phi_j`.
  - Write once (no atomics). For extreme dense blocks, consider splitting pair lists (atomics only then).
  - Interleaved voxel indexing inside workgroup for better balance.

## Test script outline
1. Build/load pentacene geometry (existing test data). 
2. Run `firecore_SCF` with 1 step.
3. Fetch dims/neighs (`get_HS_neighs`) and rho (`get_rho_sparse`).
4. Define grid (ngrid, dCell, origin) and run Grid.py host to project density.
5. Compare integrated charge to `getGridDens`; save XSF or `imshow` slices for inspection.

## Notes
- Respect host shapes everywhere; only the OpenCL kernel treats buffers as flat memory.
- Keep rho/H/S interfaces minimal and orthogonal: neighbor meta vs matrix payloads.
- Task-based launch (active blocks only) is the primary load-balancing mechanism; sorting by pair count is optional but recommended.
