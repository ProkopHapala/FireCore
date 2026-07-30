# MQCA / Hubbard / Ising Many-Body Solvers — Export Document

## Purpose & Goals

This project implements **ground-state** and **non-equilibrium steady-state** solvers for an Ising-like Hamiltonian on a lattice of interacting charge sites (quantum dots / molecular sites). The Hamiltonian is:

$$H = \sum_i n_i \left[ \epsilon_i + \sum_{j<i} W_{ij} n_j \right]$$

where $n_i \in \{0,1\}$ is the occupation of site $i$, $\epsilon_i$ is the on-site energy (including tip potential), and $W_{ij}$ is the inter-site Coulomb coupling (sparse, max 8 neighbors).

The project has **two solver families**, each with different scope:

1. **MQCA (Molecular Quantum Cellular Automata)** — Ground-state only, small systems (≤16 sites), brute-force Gray code enumeration. OpenCL only. Designed for logic gate analysis and degeneracy detection in small clusters.

2. **Hubbard/Ising** — Both ground state and non-equilibrium PME (Pauli Master Equation), large systems (up to 64+ sites). Uses Monte Carlo for ground state, kinetic subspace selection + dense PME for non-equilibrium currents. OpenCL only.

3. **PME reference solver** — Full Pauli Master Equation solver for small systems (4 sites, 16 states). Used as reference for parity validation of the Hubbard dense PME solver. **Identical copy** exists in ppafm repo (`cl/PME.cl`).

**No C++ versions** of the MQCA or Hubbard solvers exist in FireCore. The C++ PME solver (`pauli.hpp`) lives in the ppafm repo only. All solvers here are OpenCL + Python (PyOpenCL).

## Solver Comparison Matrix

| Feature | MQCA (Gray Code) | Hubbard (MC + PME) | PME.cl (Reference) |
|---------|-------------------|---------------------|---------------------|
| **System size** | ≤16 sites (2^16 states) | ≤64 sites (MC), ≤64 basis states (PME) | 4 sites (16 states, hardcoded) |
| **Ground state** | ✅ Exact (brute force) | ✅ Approximate (Monte Carlo) | ❌ (PME only) |
| **Non-equilibrium** | ❌ | ✅ (PME on reduced basis) | ✅ (full PME) |
| **Degeneracy analysis** | ✅ (top-8 variant) | ❌ | ❌ |
| **Current calculation** | ❌ | ✅ (tip current from PME probabilities) | ✅ |
| **Implementation** | OpenCL (MQCA.cl, MQCA_top8.cl) | OpenCL (hubbard.cl) | OpenCL (PME.cl) |
| **Python wrapper** | MQCASolver (MQCA.py) | HubbardSolver (HubbardSolver.py) | PauliSolverCL (pauli_ocl.py) |
| **Parallelization** | 1 workgroup per instance, Gray code sequential | 1 workgroup per tip position, threads = sites | 1 workgroup per tip position, threads cooperate on matrix |
| **W matrix format** | Sparse (W_val, W_idx, nNeigh), max 8 neigh | Sparse (W_val, W_idx, nNeigh), configurable | Dense (N_SITES=4) |
| **Logic gate analysis** | ✅ (truth tables, phase diagrams) | ❌ | ❌ |
| **Scans** | W1×W2 parameter scan | x-position, V-bias, 2D x-V maps | x-position, V-bias, 2D x-V maps |

## General Mechanics

### Ground-State Solvers (MQCA + Hubbard MC)

1. **Geometry**: N sites on a lattice (square grid, ring, or arbitrary). Each site has position (x,y) and on-site energy ε_i.
2. **Coupling**: Sparse W_ij matrix with max 8 neighbors per site. Fixed (non-active) neighbors' effect is pre-folded into ε_i as a constant bias.
3. **State encoding**: Occupation bitmask (uint16 for MQCA, uint4=128-bit for Hubbard).
4. **MQCA**: Iterates all 2^N states via Gray code (single bit flip per step), incrementally updating energy. Finds exact ground state.
5. **Hubbard MC**: Iterative local updates (flip one site at a time), with Monte Carlo acceptance. 2-phase solver with global exploration (random reset, neighbor copy). Finds approximate ground state for large N.
6. **Output**: Ground-state energy, occupation bitmask, (optionally) top-8 lowest states for degeneracy analysis.

### Non-Equilibrium PME Solvers (Hubbard + PME.cl)

1. **Tip interaction**: Site energies shifted by tip electrostatic potential (multipole + mirror charge). Tunneling amplitudes decay exponentially with tip-site distance.
2. **Many-body states**: 2^N states. State energy = Σ ε_i n_i + Σ W_ij n_i n_j.
3. **Transition rates**: Fermi's golden rule. Single-electron tunneling between site and leads (substrate μ=0, tip μ=eV). Rate depends on Fermi function and tunneling coupling.
4. **PME**: Steady-state dP/dt = W·P = 0. Rate matrix built from transition rates, normalized (ΣP=1), solved via Gauss-Jordan elimination.
5. **Current**: I = Σ_{a,b} (P_a · W_tip(a→b) - P_b · W_tip(b→a)).
6. **Hubbard approach**: MC ground state → kinetic basis scanner (selects ~48 single-flip + double-flip states) → dense PME on reduced basis (≤64 states). Scalable to large N.
7. **PME.cl approach**: Full 2^N rate matrix (N=4 only). Reference solver.

## File Listing

### OpenCL Kernels

- `/home/prokop/git/FireCore/pyBall/OCL/cl/MQCA.cl` — **MQCA ground-state solver**: Gray code brute-force, ≤16 sites, sparse W (max 8 neigh). Kernels: `mqca_groundstate` (shared geometry), `mqca_groundstate_batch_W` (per-instance W). 266 lines.
- `/home/prokop/git/FireCore/pyBall/OCL/cl/MQCA_top8.cl` — **MQCA top-8 variant**: Tracks 8 lowest energy states for degeneracy detection. Same Gray code traversal, bubble-sort maintenance of top-8 list. 149 lines.
- `/home/prokop/git/FireCore/pyBall/OCL/cl/hubbard.cl` — **Hubbard/Ising solver** (main, 1711 lines): Multi-kernel file with:
  - `eval_coupling_matrix` — compute W_ij from site positions
  - `eval_Oriented_Hopping` — tunneling amplitudes with orbital orientation
  - `solve_minBrute_fly` — brute-force ground state, on-the-fly W_ij (≤~20 sites)
  - `solve_minBrute_boltzmann` — Boltzmann-weighted thermal averaging
  - `solve_local_updates` — iterative local update / Monte Carlo (≤64 sites)
  - `solve_MC_neigh` — 2-phase MC with neighbor exchange
  - `solve_MC_2phase` — 2-phase MC with ping-pong buffers
  - `calculate_currents` — post-MC current from tunneling factors
  - `precalc_Esite_Thop` — precalculate site energies + tunneling from tip
  - `solve_pme_star_analytic` — analytic star PME (no lateral hopping)
  - `kinetic_basis_scanner` — subspace selection (single + double excitations)
  - `solve_pme_dense_batch` — dense PME on reduced basis (≤64 states), Gauss-Jordan with partial pivoting
- `/home/prokop/git/FireCore/pyBall/OCL/cl/PME.cl` — **Reference full PME solver** (447 lines): Hardcoded N_STATES=16, N_SITES=4. Kernels: `compute_tip_interaction` (tip field + tunneling), `solve_pme` (full rate matrix, parallel Gauss-Jordan). **Identical copy** in ppafm repo.
- `/home/prokop/git/FireCore/pyBall/OCL/cl/hubbard_all.cl` — Earlier version of hubbard.cl (431 lines, pre-PME integration).
- `/home/prokop/git/FireCore/pyBall/OCL/cl/hubbard_bak_bfore_PME.cl` — Backup of hubbard.cl before PME kernels were added (1681 lines).
- `/home/prokop/git/FireCore/pyBall/OCL/cl/hubbard_flash.cl` — Flash variant of hubbard.cl (869 lines, intermediate version).

### Python Wrappers

- `/home/prokop/git/FireCore/pyBall/OCL/MQCA.py` — **MQCASolver** class (810 lines): PyOpenCL wrapper for MQCA.cl + MQCA_top8.cl. Methods: `solve`, `solve_batch_W`, `solve_batch_W_top8`. Helper functions: `sq_lattice_sparse` (build W for square lattice), `apply_input_bias`, `eval_logic_table`, `identify_logic`, `scan_W1_W2`, `scan_W1_W2_top8`, `check_ground_state_uniqueness`. Visualization: `plot_ground_states`, `plot_logic_map`, `plot_logic_fraction_map`. Constants: `MAX_SITES=16`, `MAX_NEIGH=8`.
- `/home/prokop/git/FireCore/pyBall/OCL/HubbardSolver.py` — **HubbardSolver** class (1306 lines): PyOpenCL wrapper for hubbard.cl. Methods: `solve_local_updates`, `solve_mc`, `precalc_esite_thop`, `solve_pme_star`, `solve_pme_dense`, `solve_pme_dense_fullbasis`. Buffer management: `realloc_buffers`, `realloc_pme_buffers`, `realloc_local_update_buffers`, `realloc_precalc_buffers`, `realloc_mc_buffers`. Helper: `make_sparse_W`.
- `/home/prokop/git/FireCore/pyBall/OCL/pauli_ocl.py` — **PauliSolverCL** class (290 lines): PyOpenCL wrapper for PME.cl. Handles device setup, buffer management, kernel dispatch. Method: `scan_current_tip`.
- `/home/prokop/git/FireCore/pyBall/OCL/OpenCLBase.py` — **OpenCLBase** base class (741 lines): Common infrastructure for all PyOpenCL solvers. Device selection, kernel loading/caching, buffer allocation/reuse (`try_make_buffers`, `toGPU_`, `fromGPU_`).

### Application / Demo Scripts

- `/home/prokop/git/FireCore/pyBall/OCL/run_Hubbard.py` — **Demo script** for HubbardSolver (820+ lines): Functions for setting up site geometries (rings, grids), `demo_local_update` (2D scan with MC optimization + PME post-processing), plotting utilities for site occupancy, energy, tunneling patterns, and current maps.
- `/home/prokop/git/FireCore/pyBall/OCL/test_pme_parity_runall.py` — **Parity test suite** (512 lines): Compares PauliSolverCL (PME.cl full) vs HubbardSolver (hubbard.cl dense PME). Tests: 2-site and 4-site x-scans, V-scans, 2D x-V maps. Outputs CSV, TXT summaries, PNG plots.
- `/home/prokop/git/FireCore/pyBall/OCL/test_pme_parity_single.py` — Single-point parity test.
- `/home/prokop/git/FireCore/pyBall/OCL/test_pme_parity_sweep.py` — 1D sweep parity test.

### Test Scripts (tests/tMQCA/)

- `/home/prokop/git/FireCore/tests/tMQCA/test_mqca.py` — **MQCA test suite** (839 lines): Tests ground state of small clusters, logic map (W1×W2 phase diagram), geometry scan across multiple cluster shapes. Defines cluster geometries: T-shape, cross, fan, etc.
- `/home/prokop/git/FireCore/tests/tMQCA/test_top8.py` — Quick test of top-8 kernel for degeneracy detection.
- `/home/prokop/git/FireCore/tests/tMQCA/analyze_mqca_degeneracy.py` — Degeneracy analysis using top-8 kernel.
- `/home/prokop/git/FireCore/tests/tMQCA/analyze_degeneracy.py` — Earlier degeneracy analysis script.
- `/home/prokop/git/FireCore/tests/tMQCA/plot_all_degenerate.py` — Plot all degenerate states found by top-8 scan.
- `/home/prokop/git/FireCore/tests/tMQCA/plot_degenerate_states.py` — Plot individual degenerate state configurations.

### Documentation

- `/home/prokop/git/FireCore/doc/Topics/ManyBody/MQCA_GrayCodeSolver.chat.md` — **Design discussion** (1040 lines): Chat log documenting the design of the MQCA Gray code solver. Covers Hamiltonian formulation, Gray code traversal, sparse W matrix, workgroup parallelization strategy.
- `/home/prokop/git/FireCore/doc/Topics/ManyBody/HubbardSolver.md` — **Kernel reference** (314 lines): Documents each OpenCL kernel in hubbard.cl: arguments, parallelization strategy, local memory usage, buffer bindings.
- `/home/prokop/git/FireCore/doc/Topics/ManyBody/HubbardSolver_new.md` — **Updated documentation** (25.5 KB): Revised and expanded kernel documentation.
- `/home/prokop/git/FireCore/doc/Topics/ManyBody/HubbardSolver_Discussion.md` — **Discussion notes** (6.7 KB): Design decisions and open questions.
- `/home/prokop/git/FireCore/doc/Topics/ManyBody/MC_PME.progress.md` — **Progress report** (291 lines): Documents the MC + PME pipeline implementation, including physical formulation, two-solver approach, 4 critical bugs found and fixed, parity results table (max|dI| ≈ 1e-12), input parameters, and output directories.
- `/home/prokop/git/FireCore/doc/Topics/ManyBody/MC_PME_dicussion.chat.md` — **Discussion chat** (82.5 KB): Extended conversation about MC + PME integration, subspace selection, kinetic basis scanner, and debugging.

### Data Files

- `/home/prokop/git/FireCore/pyBall/OCL/Ruslan_long.txt` — 2-site geometry (long separation).
- `/home/prokop/git/FireCore/pyBall/OCL/Ruslan_kite.txt` — 4-site kite geometry.
- `/home/prokop/git/FireCore/pyBall/OCL/test_systems.xyz` — Test system XYZ file.

## Dependency Graph

```
OpenCL Kernels                Python Wrappers               Application Layer
──────────────                ───────────────               ──────────────────
MQCA.cl                       MQCASolver (MQCA.py)          tests/tMQCA/
  └─ Gray code ground state     ├─ sq_lattice_sparse         ├─ test_mqca.py
     (≤16 sites)                ├─ eval_logic_table          ├─ test_top8.py
                                ├─ scan_W1_W2               ├─ analyze_*degeneracy.py
MQCA_top8.cl                   ├─ plot_ground_states        └─ plot_degenerate*.py
  └─ Top-8 degeneracy           └─ plot_logic_map
     tracking
                              HubbardSolver                 run_Hubbard.py
hubbard.cl                      (HubbardSolver.py)           ├─ demo_local_update
  ├─ solve_minBrute_fly         ├─ solve_local_updates       ├─ plotting utils
  ├─ solve_minBrute_boltzmann   ├─ solve_mc                  └─ PME post-processing
  ├─ solve_local_updates        ├─ precalc_esite_thop
  ├─ solve_MC_neigh             ├─ solve_pme_star            test_pme_parity_*.py
  ├─ solve_MC_2phase            ├─ solve_pme_dense           ├─ parity vs PME.cl
  ├─ precalc_Esite_Thop         └─ solve_pme_dense_fullbasis  └─ 2/4-site scans
  ├─ kinetic_basis_scanner
  ├─ solve_pme_dense_batch     PauliSolverCL
  └─ calculate_currents         (pauli_ocl.py)
                                └─ scan_current_tip
PME.cl (reference)
  ├─ compute_tip_interaction   OpenCLBase (base class)
  └─ solve_pme                   ├─ device management
                                 ├─ kernel loading
                                 └─ buffer management
```

## Interconnections

### MQCA ↔ Hubbard

- **Same Hamiltonian**: Both solve the Ising-like H = Σ n_i[ε_i + Σ W_ij n_j].
- **Same sparse W format**: (W_val, W_idx, nNeigh) with max 8 neighbors.
- **MQCA is the small-system exact solver**: Designed for ≤16 sites where 2^16=65536 states are tractable by brute force. Used for logic gate design and degeneracy analysis.
- **Hubbard is the large-system scalable solver**: Monte Carlo for ground state (handles 64+ sites), plus PME for non-equilibrium currents. The `solve_local_updates` kernel is the direct analog of MQCA but with iterative updates instead of exhaustive enumeration.
- **MQCA was designed after Hubbard**: The `MQCA_GrayCodeSolver.chat.md` explicitly references hubbard.cl as the existing large-system solver and motivates MQCA as a simpler, exact solver for small systems.

### Hubbard ↔ PME.cl

- **PME.cl is the reference**: Full PME solver for 4 sites (16 states). Used to validate the Hubbard dense PME solver (`solve_pme_dense_batch` in hubbard.cl).
- **Parity achieved**: max|dI| ≈ 1e-12, max|dP| ≈ 3e-6 across all test configurations (see `MC_PME.progress.md`).
- **PME.cl is identical** in both FireCore and ppafm repos (same md5sum: `6247330ab98c3f2fcc6f073b9d6db2e4`).
- **Hubbard's PME approach**: MC ground state → kinetic basis scanner (selects relevant states) → dense PME on ≤64-state basis. This makes PME tractable for large N where full 2^N enumeration is impossible.
- **4 bugs found during parity testing** (documented in `MC_PME.progress.md`):
  1. `mpolCs` buffer 4x overallocation (float4 vs float)
  2. Active-sites basis generation including spectator sites
  3. Missing partial pivoting in Gauss-Jordan (hubbard.cl)
  4. Rate matrix dE sign convention bug (critical physics fix)

### PME.cl ↔ ppafm charge_rings

- **PME.cl is the same file** in both repos. The ppafm repo has the full charge rings project with C++ CPU solver (`pauli.hpp`), Python bindings, GUIs, and application scripts. FireCore has only the OpenCL kernel + Python wrapper (`pauli_ocl.py`).
- **ppafm is the primary home** for charge rings / PME simulation infrastructure. FireCore's copy is for parity validation of the Hubbard dense PME solver.

## Key Differences Summary

### Small vs Large Systems

| Aspect | Small (≤16 sites) | Large (up to 64+ sites) |
|--------|-------------------|-------------------------|
| **Ground state** | MQCA.cl (exact, Gray code) | hubbard.cl `solve_local_updates` / `solve_MC_neigh` (approximate, MC) |
| **PME** | PME.cl (full 2^N, N≤4) | hubbard.cl `solve_pme_dense_batch` (reduced basis ≤64 states) |
| **Degeneracy** | MQCA_top8.cl (top-8 tracking) | Not supported |
| **Logic gates** | MQCA.py helpers | Not supported |

### Ground State vs Non-Equilibrium

| Aspect | Ground State Only | Non-Equilibrium (PME) |
|--------|-------------------|-----------------------|
| **Solvers** | MQCA.cl, hubbard.cl (MC kernels) | PME.cl, hubbard.cl (PME kernels) |
| **Physics** | Find min-energy occupation | Steady-state probabilities + currents |
| **Tip model** | Static bias on ε_i | Full tip interaction (multipole, tunneling, bias) |
| **Output** | Energy + occupation bitmask | Current + probabilities + per-site current |
| **Temperature** | Not in MQCA; kT in MC acceptance | Fermi-Dirac in transition rates |

### C++ vs OpenCL

| Aspect | C++ | OpenCL |
|--------|-----|--------|
| **MQCA** | ❌ Not implemented | ✅ MQCA.cl, MQCA_top8.cl |
| **Hubbard MC** | ❌ Not implemented | ✅ hubbard.cl |
| **Hubbard PME** | ❌ Not implemented | ✅ hubbard.cl (dense PME) |
| **Full PME** | ✅ ppafm: `pauli.hpp` | ✅ PME.cl (both repos) |
| **Python bindings** | ✅ ppafm: `pauli.py` (ctypes) | ✅ FireCore: PyOpenCL wrappers |

## External Dependencies

- **NumPy** — array operations throughout
- **Matplotlib** — all plotting and visualization
- **PyOpenCL** — all GPU solvers
- **OpenCL** runtime (NVIDIA CUDA, AMD ROCm, or Intel)

## Build / Run Instructions

No compilation needed — OpenCL kernels are compiled at runtime when solver classes are instantiated.

### Running MQCA tests
```bash
cd /home/prokop/git/FireCore/tests/tMQCA
PYTHONPATH=../.. python -u test_mqca.py
PYTHONPATH=../.. python -u test_top8.py
```

### Running Hubbard demos
```bash
cd /home/prokop/git/FireCore
python -u -m pyBall.OCL.run_Hubbard
```

### Running PME parity tests
```bash
cd /home/prokop/git/FireCore
python -u -m pyBall.OCL.test_pme_parity_runall
```

## Related Topical Audit Entry

In `doc/topical_audit/topical_audit.md`, this work falls under the broader "AFM/STM Simulation" topic (Priority 4), specifically the charge rings / many-body PME subtopic. The topical audit does not yet have a dedicated entry for MQCA or Hubbard solvers — they are part of the "Missing Topics" that need auditing. The charge rings export document in ppafm (`docs/export/charge_rings.export.md`) covers the C++ and application layer; this document covers the OpenCL solver layer in FireCore.
