# Density Projection (Fireball → pyOpenCL) — TDD Plan

## 1. Scope & Goals
- Reimplement projection of Fireball electron density onto arbitrary points and 3D grids using pyOpenCL.
- Start with CPU (NumPy/Numba) prototypes, then port kernels to GPU with load-balanced execution.
- Inputs come from existing Fortran/Python APIs (orbital coefficients, sparse H/S, neighbor lists).
- Deliver test-driven milestones with reproducible scripts (acene generator → Fireball SCF → density projection comparisons).

## 2. Key References (current code)
- Fortran exports
  - `firecore_dens2points` (density at points) and `firecore_get_wfcoef`, `firecore_get_HS_sparse`, `firecore_get_HS_dims` in `fortran/MAIN/libFireCore.f90` @/home/prokophapala/git/FireCore/fortran/MAIN/libFireCore.f90#331-399.
- Python bindings
  - `dens2points`, `get_wfcoef`, `get_HS_dims`, `get_HS_sparse` in `pyBall/FireCore.py` @/home/prokophapala/git/FireCore/pyBall/FireCore.py#331-399.
- Example usage
  - Density along a bond using `dens2points` in `tests/pyFireball/density_along_line.py` @/home/prokophapala/git/FireCore/tests/pyFireball/density_along_line.py#67-77.
- Design notes and optimization ideas
  - `fortran/doc/DensityToGrid/DensToGrid_Opet.md` (B-splines over cubic coeffs, super-voxel broad-phase, atom-centric DM batching, Morton/Z-curve clustering ideas, load balancing discussions).

## 3. TDD Milestones
1) **Geometry/Test Harness**
   - Script to build linear acenes with arbitrary length (positions, atom types).
   - Deterministic Fireball SCF run via Python API; fixture saves positions, `wfcoef`, `H/S` sparse export, neighbor lists.
2) **CPU Density Matrix Prototype (Numba)**
   - Build sparse DM: ρ_{μν} = Σ_occ f_n C_{nμ} C_{nν} over neighbor pairs only.
   - Validate vs Fireball `dens2points` on small grids/lines.
3) **CPU Grid Projection Prototype (Numba)**
   - Project sparse DM onto arbitrary points and small 3D grids.
   - Compare voxel/point results to `dens2points` (and `firecore_getGridDens` if applicable).
4) **GPU Step 1: DM Construction (pyOpenCL)**
   - Atom-centric batching kernel with local caching of atom + neighbor coefficients (per DensToGrid notes).
   - Validate DM blocks vs CPU NumPy/Numba reference.
5) **GPU Step 2: Grid Projection (pyOpenCL)**
   - Super-voxel broad-phase + thread-level sub-voxel filtering; B-spline evaluation (nodes, not cubic segments).
   - Validate density vs CPU reference on identical grids.
6) **Load Balancing & Scaling Tests**
   - Sweep block sizes, neighbor batch sizes, Morton/Z-order (or alternative) atom sorting; measure occupancy/bandwidth.
7) **Integration & Regression**
   - CLI/pytest entrypoints to regenerate references and compare; hook into `tests/pyFireball` style.

## 4. Data & Interfaces Needed
- From Fireball (Python layer):
  - `get_HS_dims()` → dimensions (natoms, norbitals, neighbor caps, etc.).
  - `get_HS_sparse(...)` → neighbor lists (`neigh_j`, `neigh_b`), orbital block sizes, H/S blocks (shape: natoms × neigh_max × numorb_max × numorb_max).
  - `get_wfcoef(ikp=1)` → MO coefficients C[n, μ]; occupations from nelec and Fermi filling (or Fireball-provided).
  - `dens2points(points, f_den, f_den0)` → ground-truth density for tests.
- Derived structures for projection kernels:
  - CSR-like neighbor list: `neigh_offsets`, `neigh_indices`, `block_offsets` for DM blocks.
  - B-spline radial nodes per (species, l) for basis evaluation (prefer node values over 4×coeffs per segment).
  - Morton/Z-order (optional) atom ordering to improve cache reuse between adjacent work-groups.

## 5. Algorithms (CPU first)
- **Sparse DM build (Numba)**
  - Loop bands in batches; outer products on neighbor pairs only (use neighbor list to skip non-overlaps).
  - Store DM blocks as contiguous `[pair][mu][nu]` with strides matching GPU plan.
- **Projection onto points/grid (Numba)**
  - For each point, traverse atoms in cutoff using neighbor list; evaluate basis via B-spline nodes; accumulate ρ = Σ_{μν} ρ_{μν} φ_μ φ_ν.
  - Grid version: tile grid, reuse filtered atom sets across voxels inside tile.

## 6. GPU Design (pyOpenCL) — mapped to notes in DensToGrid_Opet.md
- **DM kernel (Atom-cluster caching)**
  - Work-group per atom; load atom + neighbors’ C coefficients for a band batch into local memory; compute all μ,ν pairs (4×4 sp3 or 9×9) and atomic-add to global DM buffer.
  - If neighbors don’t fit, process in neighbor batches; reuse atom coefficients across all neighbors.
- **Grid kernel (3-level filtering)**
  - Work-group size 32–64; each thread handles a 2×2×2 sub-block of voxels.
  - Level 1: Super-voxel broad-phase using all threads; compact active atoms into local list via local atomics (dense list, no -1 gaps).
  - Level 2: Thread-box filter in registers for its sub-block (Rcut + margin); build per-thread small active list.
  - Level 3: Per-voxel projection; evaluate B-spline; accumulate density.
  - If active-atom list exceeds local mem, iterate batches.
- **Basis evaluation**
  - Use cubic B-spline nodes (1 float/node); per-eval fetch 4 nodes (i-1…i+2) and apply fixed polynomial; prefer buffer reads over hardware interpolation (accuracy requirement).

## 7. Testing Strategy
- Unit tests (pytest):
  - Density at random points vs `dens2points` (scalar tolerance).
  - Small grid (e.g., 8³) vs CPU/Fireball reference.
  - DM block correctness vs NumPy reference for selected atom pairs.
- Property tests:
  - Translational invariance (shift molecule → shift grid results).
  - Cutoff behavior (points outside Rcut yield ~0 within tolerance).
- Performance checks (not gating):
  - Kernel occupancy/logging of active-atom counts per tile; timing of DM vs projection kernels.

## 8. Open Questions / Decisions
- Occupations source: use Fireball occupations/fermi level or reconstruct from eigenvalues + nelec.
- Pair clustering heuristic: start with Morton/Z-order on atoms; evaluate if pair-level clustering is needed; consider greedy neighbor-sharing if hotspots persist.
- Grid tile size and cutoff margin: start with tile size ~2–3 Å (e.g., 16³ voxels at 0.125 Å) and margin Rcut+0.5×dgrid; tune empirically.

## 9. Next Actions
1) Write acene generator + Fireball SCF harness to dump `wfcoef`, `H/S` sparse, neighbor lists, and reference `dens2points` data.
2) Implement CPU Numba DM builder + point/grid projector; add tests comparing to references.
3) Define buffer layouts (CSR neighbor, DM blocks, B-spline nodes) matching planned kernels.
4) Port DM kernel to pyOpenCL; validate against CPU.
5) Port grid kernel with 3-level filtering; validate vs CPU and Fireball `dens2points`.
