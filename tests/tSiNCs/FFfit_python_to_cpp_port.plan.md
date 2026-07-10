# C++ Porting Plan: FFfit Pipeline Hotspots

> **Status (updated 2025-07):** Phases 1–3 COMPLETE. All Tier 1 graph algorithms,
> Tier 3.3–3.5 dihedral/Wilson functions, and batch dihedral sensitivity are ported
> to C++ with parity tests (14 tests in `test_parity_graph_cpp.py`, 7 in
> `test_parity_py_cpp.py`). Optimizations: CSR bond-graph, bounded ring BFS,
> thread_local buffers, in-place FD with symmetry exploitation, batch typed
> accumulation. Phase 4 (hybrid assembly orchestration) and Tier 2 (typing reform)
> remain as future work.

> **Note (refactoring):** The monolithic `test_FFfit.py` has been split into:
> - `pyBall/FFfit_utils.py` — all reusable logic (type system, topology, dihedral physics, parameter mapping, sensitivity, fitting, frequency analysis, C++ bridge)
> - `pyBall/FFfit_plots.py` — visualization (spectra, equilibrium distributions, stiffness HTML maps)
> - `tests/tSiNCs/test_FFfit.py` — thin CLI wrapper importing from both modules
>
> Line numbers below refer to the **original** `test_FFfit.py` before refactoring.
> Functions are now in `pyBall/FFfit_utils.py`; search by function name.

## Analysis Method
Profiled Python loop patterns in `test_FFfit.py` and `pyBall/FFfit.py` against
existing C++ infrastructure (`FFfit.h`, `FFfit_lib.cpp`, `MMFFBuilderBase.h`,
`MMFFparams.h`). Categorized by loop type, data dependency, and existing C++ support.

---

## Tier 1: Graph/Topology Algorithms (Pure integer loops, 10-100x speedup) — ✅ COMPLETE

### 1.1 `shortest_path_distances(bond_pairs, natoms)` — BFS all-pairs — ✅ DONE
- **Location**: `test_FFfit.py:117-133`
- **Pattern**: N BFS traversals over adjacency list → O(N * (N+E)) Python loop
- **Called by**: `build_3rd_neighbor_bonds`, `local_hessian_mask`
- **C++ infrastructure**: `MMFFBuilderBase` already has adjacency via `AtomConf.neighs[]`
  and bond graph traversal. `autoBonds`/`autoBondsPBC` do similar neighbor finding.
- **Port plan**: Add `int* compute_bond_graph_distances(natoms, bonds, nbonds)` to FFfit.h
  — returns flattened (natoms*natoms) int array. BFS per source, same as Python.
  Expose via `fffit_compute_graph_distances()` in `FFfit_lib.cpp`.
- **Atom type ID opportunity**: Replace string-based `symbols[i]` with integer atom type
  IDs from `MMFFparams::atomTypeDict`. Already available as `Atom.type` (int).

### 1.2 `build_3rd_neighbor_bonds(symbols, positions, bond_pairs, max_dist)` — 1-4 pairs — ✅ DONE
- **Location**: `test_FFfit.py:136-151`
- **Pattern**: O(N²) scan over distance matrix + distance computation
- **Depends on**: `shortest_path_distances` (1.1)
- **Port plan**: Fold into 1.1 — after BFS, scan dist==3 pairs, compute |r_j-r_i|,
  filter by max_dist. Single C++ function.

### 1.3 `build_dihedrals(symbols, positions, bonds, d, n, dihedral)` — torsion enumeration — ✅ DONE
- **Location**: `test_FFfit.py:249-276`
- **Pattern**: Triple-nested loop over bonds × neighbors × neighbors
- **C++ infrastructure**: `MMFFBuilderBase` has `AtomConf.neighs[]` (bond-based adjacency)
  and `getBondToNeighbor()`. Dihedral enumeration is standard topology building.
- **Port plan**: Add `void enumerate_dihedrals(bonds, natoms, dihedrals_out)` to FFfit.h.
  Uses adjacency from bond list. Returns (i,j,k,l,d,n) tuples as flat int array.

### 1.4 `local_hessian_mask(natoms, bonds, max_graph_distance)` — BFS with cutoff — ✅ DONE
- **Location**: `pyBall/FFfit.py:555-594`
- **Pattern**: N BFS traversals with early termination at max_graph_distance
- **Called by**: `assemble_hybrid_hessian_system` (every hybrid fit)
- **Port plan**: Same BFS core as 1.1, just with cutoff. Return boolean mask as
  `uint8_t*` array of size (3N * 3N). Expose via `fffit_local_hessian_mask()`.

### 1.5 `build_topology(symbols, positions, bond_cutoff, ...)` — bond/angle discovery — ⬜ NOT DONE
- **Location**: `test_FFfit.py:307-370`
- **Pattern**: O(N²) distance scan for bonds, then O(bonds²) for angles
- **C++ infrastructure**: `MMFFBuilderBase::autoBonds()` already does distance-based
  bond finding with type-aware radii. Angle building is in `MMFFBuilderBase` (not shown
  but standard: for each atom pair sharing a bond, find common neighbors).
- **Port plan**: Use `autoBonds()` for bond finding. Add angle enumeration to FFfit.h
  or reuse from MMFFBuilderBase. Return bonds/angles as flat arrays.

---

## Tier 2: Type Assignment & Parameter Mapping (Dict loops, 5-20x speedup)

### 2.1 `assign_si_environment_types(symbols, bonds, enabled)` — Si subtyping
- **Location**: `test_FFfit.py:35-49`
- **Pattern**: Loop over bonds counting H neighbors, then classify Si atoms
- **Issue**: Uses string comparisons (`symbols[i] == 'Si'`)
- **Port plan**: Replace strings with integer atom type IDs. Use `MMFFparams::atypes`
  to get element Z from type. Count H neighbors by Z==1. Assign subtype as integer
  enum: `enum SiSubtype { Si_bulk=0, SiH=1, SiH2=2, SiH3=3 }`.
  Add to FFfit.h as `void assign_environment_types(atom_types_out, natoms, bonds, nbonds, Z_array)`.

### 2.2 `interaction_type_counts(systems)` — type census
- **Location**: `test_FFfit.py:52-69`
- **Pattern**: Loop over all systems × all interactions, dict counting
- **Port plan**: With integer type IDs, this becomes array histogramming —
  `int counts[n_interaction_types]` incremented per term. Trivially fast in C++.
  Add as `void count_interaction_types(systems, n_systems, counts_out)`.

### 2.3 `build_global_param_map(...)` — global type indexing
- **Location**: `test_FFfit.py:746-794`
- **Pattern**: 4 passes over all systems × all interactions, dict-based type discovery
- **C++ infrastructure**: `ParamMap::assign_bond_types()` and `assign_angle_types()`
  already do this for single-system. Need multi-system version.
- **Port plan**: Extend `ParamMap` with `assign_bond_types_multi()` and
  `assign_angle_types_multi()` that take arrays of systems. Use `std::map` with
  integer keys (packed atom type pairs) instead of string keys.
  Expose via `fffit_assign_global_param_map()` in FFfit_lib.cpp.

### 2.4 `compute_averaged_equilibrium(...)` — type-averaged l0, c0
- **Location**: `test_FFfit.py:796-870`
- **Pattern**: Loop over all systems × all bonds/angles, accumulate by type, then average
- **Port plan**: Single C++ pass accumulating sums and counts per type, then divide.
  Use integer type indices as array indices. Add to FFfit.h.

---

## Tier 3: Sensitivity & Hessian Assembly (Numerical loops, already partially ported)

### 3.1 `build_sensitivity_matrices(...)` — Python reference (DONE)
- **Location**: `test_FFfit.py:449-532`
- **Status**: ✅ Already ported to C++ (`FFfit::build_sensitivity_matrices()`).
  Parity verified to machine precision (test_parity_py_cpp.py).
- **Keep Python**: As reference for parity testing.

### 3.2 `compute_gradient_term_blocks(...)` — Python reference (DONE)
- **Location**: `test_FFfit.py:584-674`
- **Status**: ✅ Already in C++ (`FFfit::compute_gradient()`).
- **Keep Python**: As reference.

### 3.3 `dihedral_hessian(pos, d, n, h)` — Finite-difference Hessian — ✅ DONE
- **Location**: `test_FFfit.py:226-240`
- **Pattern**: 12×12 Hessian via 24 gradient evaluations (central FD), each gradient
  involves cross products, norms, trig — O(12 * 24 * ~50) flops in Python
- **C++ infrastructure**: `UFF.h` has `evalDihedral_Prokop()` with analytical gradient.
  Analytical Hessian can be derived from the gradient (chain rule on the 4 atom forces).
- **Port plan**: Port `dihedral_energy_gradient()` to C++ (mirrors `evalDihedral_Prokop`).
  Then either: (a) use analytical Hessian derivation, or (b) use FD in C++ (still 100x
  faster than Python FD). Option (b) is simpler and sufficient for fitting.
  Add `dihedral_dHdk()` to FFfit.h alongside bond_dHdk/angle_dHdk.

### 3.4 `compute_dihedral_sensitivity(...)` — Full (3N,3N) dihedral A_p — ✅ DONE
- **Location**: `test_FFfit.py:279-298`
- **Pattern**: Loop over dihedrals, compute 12×12 Hessian, scatter into (3N,3N)
- **Depends on**: 3.3
- **Port plan**: After 3.3, this is just scatter-add into flat array — same pattern
  as `bond_dHdk` / `angle_dHdk`. Add to `FFfit::build_sensitivity_matrices()`.

### 3.5 `build_wilson_matrix(positions, bonds, angles)` — B matrix — ✅ DONE
- **Location**: `pyBall/FFfit.py:424-480`
- **Pattern**: Loop over bonds (2 atoms each) + angles (3 atoms each), fill sparse rows
- **C++ infrastructure**: `bond_wilson()` and `angle_wilson_cos()` already exist in FFfit.h
- **Port plan**: Add `void build_wilson_matrix(apos, bonds, angles, B_out, labels_out)`
  to FFfit.h. Loop calling existing wilson functions, scatter into B flat array.
  Expose via `fffit_build_wilson_matrix()`.

### 3.6 `internal_hessian_projection(H, B, masses, ...)` — Wilson GF method
- **Location**: `pyBall/FFfit.py:511-552`
- **Pattern**: SVD of C = B*M^{-1/2}, then F = C^{+T} D C^{+}
- **Status**: The matmul-heavy parts were already optimized (einsum→matmul, 2500x speedup).
  The remaining cost is SVD (LAPACK) which is already fast.
- **Port plan**: LOW PRIORITY. The Python version uses `np.linalg.svd` (LAPACK) which
  is already C-backed. The only Python overhead is the wrapper. Port only if profiling
  shows significant time here for large systems.

### 3.7 `assemble_hybrid_hessian_system(...)` — Hybrid objective assembly
- **Location**: `pyBall/FFfit.py:609-712`
- **Pattern**: Calls build_wilson_matrix, internal_hessian_projection, mass_weight_hessian,
  reference_vibrational_modes, local_hessian_mask — all already vectorized or ported.
- **Status**: The matmul parts are already optimized. Remaining Python overhead is
  orchestration (calling order, normalization).
- **Port plan**: MEDIUM PRIORITY. Once 1.4 and 3.5 are ported, the assembly can be
  a single C++ function that calls all sub-steps without Python round-trips.
  Add `fffit_assemble_hybrid_system()` to FFfit_lib.cpp.

---

## Tier 4: I/O and Statistics (Low priority, minimal loops)

### 4.1 `equilibrium_distributions(systems, ...)` — collect relaxed values
- **Location**: `test_FFfit.py:881-1001`
- **Pattern**: Loop over systems × interactions, append to dict
- **Port plan**: LOW PRIORITY. This is I/O-bound (reading files) and the loop is
  over systems, not atoms. Python is fine here.

### 4.2 `dft_stiffness_distributions(systems)` — Wilson GF diagonal extraction
- **Location**: `test_FFfit.py:1004-1075` (user-added)
- **Pattern**: Loop over systems, call build_wilson_matrix + internal_hessian_projection
- **Port plan**: Benefits from 3.5 porting, but the per-system loop is small.
  LOW PRIORITY on its own.

### 4.3 `load_pyscf_case(case_dir, ...)` — file I/O
- **Location**: `test_FFfit.py:82-111`
- **Pattern**: File reading + string parsing
- **Port plan**: NOT NEEDED. I/O-bound, Python is appropriate.

---

## Cross-Cutting: Atom Type ID Reform

### Problem
Current Python code uses string comparisons everywhere:
```python
if symbols[i] == 'Si' and symbols[j] == 'H': nH[i] += 1
key = tuple(sorted([symbols[i], symbols[j]]))
```

### Solution: Integer atom type IDs
`MMFFparams` already has `AtomType` with `uint8_t iZ` (atomic number) and
`int type` (atom type index). Use these throughout:

```cpp
// Instead of string comparison:
if (Z[i] == 14 && Z[j] == 1) nH[i]++;  // Si=14, H=1

// Instead of string-sorted tuple keys:
uint64_t bond_key = pack64(min(Z[i],Z[j]), max(Z[i],Z[j]), 0, 0);
// Or for environment types:
uint16_t type_i = atom_types[i];  // already an int from MMFFparams
uint64_t key = BondType::getId(type_i, type_j, order);
```

### Migration plan
1. Add `int* atom_type_ids` and `int* Z_array` to FFfit alongside `symbols`
2. Python side: convert symbols→type IDs once at load time using `MMFFparams::getAtomType()`
3. All new C++ functions accept `int* atom_types` instead of `char** symbols`
4. Keep `symbols` for human-readable output only

### Enum for Si subtypes
```cpp
enum class SiEnv : uint8_t {
    Bulk = 0,  // Si with 0 H neighbors
    SiH  = 1,  // Si with 1 H neighbor
    SiH2 = 2,  // Si with 2 H neighbors
    SiH3 = 3,  // Si with 3+ H neighbors
};
```
Assign as `atom_types[i] = base_type + SiEnv::SiH2` or use separate `env_type` array.

---

## Implementation Priority & Effort/Impact

| Priority | Function | Effort | Impact | Notes |
|----------|----------|--------|--------|-------|
| **P0** ✅ | 1.1 shortest_path_distances | Small | High | DONE — `bond_graph_distances` in FFfit.h, CSR+BFS |
| **P0** ✅ | 1.4 local_hessian_mask | Small | High | DONE — `local_hessian_mask` + combined `local_mask_and_14pairs` |
| **P0** ✅ | 3.5 build_wilson_matrix | Small | High | DONE — `build_wilson_matrix` in FFfit.h, sparse fill |
| **P1** | 1.5 build_topology | Medium | High | Replaces O(N²) Python |
| **P1** | 2.3 build_global_param_map | Medium | Medium | Multi-system, use int keys |
| **P1** ✅ | 3.3 dihedral_hessian | Medium | Medium | DONE — `dihedral_hessian_fd` with in-place FD + symmetry |
| **P2** | 2.1 assign_si_env_types | Small | Low | Simple but string→int reform |
| **P2** | 2.2 interaction_type_counts | Small | Low | Trivial with int types |
| **P2** | 2.4 compute_averaged_equilibrium | Medium | Low | Accumulate+average |
| **P2** ✅ | 1.2 build_3rd_neighbor_bonds | Small | Low | DONE — `find_3rd_neighbor_bonds` in FFfit.h |
| **P2** ✅ | 1.3 build_dihedrals | Small | Low | DONE — `enumerate_dihedrals` in FFfit.h |
| **P2** ✅ | 3.4 compute_dihedral_sensitivity | Medium | Medium | DONE — `dihedral_dHdk_batch_typed` replaces Python loop |
| **P3** | 3.7 assemble_hybrid_system | Large | Medium | Orchestration, needs 1.4+3.5 |
| **P3** | 3.6 internal_hessian_proj | Large | Low | Already LAPACK-backed |
| **P3** | 4.1-4.3 I/O & stats | — | Low | Python-appropriate |

---

## Suggested Implementation Order

### Phase 1: Graph algorithms + Wilson matrix (P0) — ✅ COMPLETE
1. ✅ Port `shortest_path_distances` → `bond_graph_distances` (CSR bond-graph + BFS)
2. ✅ Port `local_hessian_mask` → `local_hessian_mask` + combined `local_mask_and_14pairs`
3. ✅ Port `build_wilson_matrix` → `build_wilson_matrix` (sparse fill, pre-zeroed buffer)
4. ✅ Add C wrappers in `FFfit_lib.cpp`
5. ✅ Add Python bindings in `FFfit.py`
6. ✅ Add parity tests in `test_parity_graph_cpp.py` (12 tests)

### Phase 2: Topology + typing reform (P1-P2) — PARTIALLY COMPLETE
1. ⬜ Add `atom_type_ids` / `Z_array` to FFfit class
2. ⬜ Port `build_topology` using `autoBonds`-style distance scan
3. ⬜ Port `build_global_param_map` with integer keys
4. ⬜ Port `assign_si_environment_types` with enum
5. ✅ Port `build_dihedrals` → `enumerate_dihedrals`
6. ✅ Port `build_3rd_neighbor_bonds` → `find_3rd_neighbor_bonds`
7. ⬜ Port `compute_averaged_equilibrium`
8. ✅ Parity tests for completed items

### Phase 3: Dihedral sensitivity (P1) — ✅ COMPLETE
1. ✅ Port `dihedral_energy_gradient` to C++ (mirrors `evalDihedral_Prokop`)
2. ✅ Add `dihedral_hessian_fd` (in-place FD, symmetry exploitation, precomputed inv_2h)
3. ✅ Add `dihedral_dHdk` (single) + `dihedral_dHdk_batch` + `dihedral_dHdk_batch_typed`
4. ✅ Parity tests (2 batch tests: vs single=0.00e+00, vs Python=1.90e-06)

### Phase 4: Hybrid assembly orchestration (P3, optional) — NOT STARTED
1. ⬜ Single `fffit_assemble_hybrid_system()` C++ function
2. ⬜ Calls all sub-steps without Python round-trips
3. ⬜ Only needed if Python orchestration overhead is significant for large systems

---

## What NOT to Port (Keep in Python)
- `load_pyscf_case` — file I/O
- `equilibrium_distributions` / `dft_stiffness_distributions` — I/O + plotting
- `solve_regularized_lsq` — uses `np.linalg.lstsq` (LAPACK), already C-backed
- `reference_vibrational_modes` — uses `np.linalg.eigh` (LAPACK)
- `rigid_and_vibrational_bases` — SVD-based, already C-backed
- `internal_hessian_projection` — SVD + matmul, already optimized
- All plotting functions
- Python reference implementations (keep for parity testing)
