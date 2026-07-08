CODEMAP:
Hydrogen Transfer & Kekule Structures (DFTB+ Integration)
https://windsurf.com/codemaps/07d864ef-4614-40e5-a45e-d6be773de1c4-fe86ab10a43f3d18


---

# USER

Can you review this problem of Kekule structures and Nudged_Elstic Band and other such calculations with DFTB, I would like to extract them to different repo in some consolidated form. What exactly we need to extract? which python modules, and also where are the results so I can review it ? I think we have quite many script related to this Kekule ahdn hydrogen tranfer between nitrogens, both for molecules and for ribbons (extended systems with k-points)

do not code just make audit of existing code and results

---

# SWE-1.6

# Kekule/NEB/H-transfer Code Audit

## 1. Core Python Modules to Extract

### Kekule Structure Explorer
- **[pyBall/KekuleBackend.py](pyBall/KekuleBackend.py)**
  - Backend for Kekule Structure Explorer
  - Manages AtomicSystem state with hexagonal grid metadata
  - Hydrogen passivation logic (sp3 H orientation fix recently applied)
  - Functions: [honeycomb_ring_nodes()](pyBall/KekuleBackend.py:31-50), [_calc_h_directions_atom()](pyBall/KekuleBackend.py:822-900), [adjust_h()](pyBall/KekuleBackend.py:808-820)

- **[pyBall/KekuleExplorerGUI.py](pyBall/KekuleExplorerGUI.py)**
  - PyQt5 GUI for Kekule explorer
  - Vispy-based 3D rendering
  - Extension system (ExtensionManager)

- **[pyBall/Kekule.py](pyBall/Kekule.py)**
  - Python ctypes wrapper for C++ `Kekule_lib` (bond-order optimizer)
  - Functions: `init()`, `relax()`, `setDefaultBondOrders()`, `pinBondOrders()`

### Kekule Topology Analysis (doc/Topics/Kekule_Topology/)
- **[GrapheneRibbonBuilder.py](doc/Topics/Kekule_Topology/GrapheneRibbonBuilder.py)** - Build graphene ribbons with various passivations
- **[KekuleSolver.py](doc/Topics/Kekule_Topology/KekuleSolver.py)** - Kekule pattern solver
- **[KekuleCLI.py](doc/Topics/Kekule_Topology/KekuleCLI.py)** - Command-line interface
- **[KekuleGUI.py](doc/Topics/Kekule_Topology/KekuleGUI.py)** - GUI for Kekule topology
- **[Kekule_Optimizer.md](doc/Topics/Kekule_Topology/Kekule_Optimizer.md)** - Documentation
- **[Kekule_Topology.md](doc/Topics/Kekule_Topology/Kekule_Topology.md)** - Theoretical background (Jackiw-Rebbi model, grain boundaries)

### DFTB Integration
- **[pyBall/dftb_utils.py](pyBall/dftb_utils.py)**
  - Subprocess-based DFTB+ interface
  - HSD generation: [makeDFTBjob()](pyBall/dftb_utils.py:607-700), [makeDFTBjob_pbc()](pyBall/dftb_utils.py:702-801)
  - Waveplot integration: [run_waveplot()](pyBall/dftb_utils.py:939-1016), [read_cube()](pyBall/dftb_utils.py:1019-1038)
  - Hessian/vibrations: [get_hessian_ase()](pyBall/dftb_utils.py:830-864), [hessian_to_frequencies()](pyBall/dftb_utils.py:890-934)

- **[pyBall/FireCore.py](pyBall/FireCore.py)**
  - C bindings for Fireball/DFTB
  - PBC support: `firecore_set_lvs()`, `firecore_set_kpoints()`
  - k-point sampling for periodic systems

## 2. H-transfer NEB Scripts

### Molecular Systems (gamma-point)
- **[tests/pyFireball/neb_h_transfer_molecules.py](tests/pyFireball/neb_h_transfer_molecules.py)**
  - NEB for H-transfer between hydro/dehydro molecules (phenazine, pyrazine)
  - Modes: `rigid_scan`, `relax_scan`, `dry_run`
  - H-bond detection: [find_hydrogen_bonds()](tests/pyFireball/neb_h_transfer_molecules.py:230-284)
  - PBC-aware placement: [build_molecular_cell()](tests/pyFireball/neb_h_transfer_molecules.py:75-175) with NN-distance mode
  - Input files: [Phenazine/phenazine_hydro.xyz](tests/pyFireball/Phenazine/phenazine_hydro.xyz), [phenazine_dehydro.xyz](tests/pyFireball/Phenazine/phenazine_dehydro.xyz), `pyrazine_*.xyz`

### Ribbon Systems (periodic with k-points)
- **[tests/pyFireball/neb_h_transfer.py](tests/pyFireball/neb_h_transfer.py)**
  - NEB for H-transfer between N-passivated and NH-passivated ribbons
  - Periodic along x with k-point sampling (nk_x=16)
  - Uses GenFormat for DFTB+ input

- **[tests/pyFireball/build_two_ribbons.py](tests/pyFireball/build_two_ribbons.py)**
  - Build two-ribbon unit cell for H-transfer
  - Uses `GrapheneRibbonBuilder` from Kekule topology module
  - N and NH passivation

- **[tests/pyFireball/scan_ribbon_dftb.py](tests/pyFireball/scan_ribbon_dftb.py)**
  - Lattice scan for zigzag ribbons using DFTB+
  - PBC along x with k-points

- **[tests/pyFireball/scan_constrained.py](tests/pyFireball/scan_constrained.py)**
  - Constrained H-transfer scan (move H between ribbons)
  - Fix N atoms, relax others

- **[tests/pyFireball/scan_LHb.py](tests/pyFireball/scan_LHb.py)**
  - Scan energy vs H-bond length

- **[tests/pyFireball/relax_ribbon.py](tests/pyFireball/relax_ribbon.py)**
  - Lattice scan for N-doped graphene ribbons (pyridinic, pyrrolic) using FireCore relax

## 3. Test Scripts and Results

### Kekule Explorer Tests
- **[tests/tKekuleExplorer/test_kekule_topology.py](tests/tKekuleExplorer/test_kekule_topology.py)**
- **[tests/tKekuleExplorer/test_kekule_hbonds.py](tests/tKekuleExplorer/test_kekule_hbonds.py)**
- **[tests/tKekuleExplorer/test_kekule_relax.py](tests/tKekuleExplorer/test_kekule_relax.py)**
- **[tests/tKekuleExplorer/test_ribbon_parity.py](tests/tKekuleExplorer/test_ribbon_parity.py)**

### Kekule Bond-Order Optimizer Tests
- **[tests/tKekule/run.py](tests/tKekule/run.py)** - Test script for `pyBall.Kekule` bond-order relaxation on anthracene fragments
- **[tests/tKekule/generate_2.py](tests/tKekule/generate_2.py)** - Generate Kekule structure variants with donor/acceptor groups (-OH, -NH2, =O, =N-)

**Results directories:**
- [tests/tKekuleExplorer/out_topology/](tests/tKekuleExplorer/out_topology) - 18 geometry snapshots (benzene, naphthalene, anthracene, etc.)
- [tests/tKekuleExplorer/out_relax/](tests/tKekuleExplorer/out_relax) - DFTB relaxation results

### DFTB Test Scripts
- **`tests/dftb/example_hessian.py`** - Hessian calculation
- **`tests/dftb/example_orbitals.py`** - Orbital visualization
- **`tests/dftb/test_*.py`** - C-API tests

### PyFireball Results
**Energy profiles:**
- [neb_energy_profile_mol.png](tests/pyFireball/neb_energy_profile_mol.png) - NEB energy profile for molecules
- [relax_scan_energy.png](tests/pyFireball/relax_scan_energy.png) - Relax scan results
- [rigid_scan_energy.png](tests/pyFireball/rigid_scan_energy.png) - Rigid scan results
- `LHb_scan*.png` - H-bond length scans
- `lattice_scan_*.png` - Lattice optimization scans

**Geometry movies:**
- [neb_results_mol.npz](tests/pyFireball/neb_results_mol.npz) - NEB trajectory data
- [relax_scan_movie.xyz](tests/pyFireball/relax_scan_movie.xyz) - Relax scan trajectory
- [rigid_scan_movie.xyz](tests/pyFireball/rigid_scan_movie.xyz) - Rigid scan trajectory
- [dry_run_movie.xyz](tests/pyFireball/dry_run_movie.xyz) - Geometry-only check

**Input files:**
- [Phenazine/phenazine_hydro.xyz](tests/pyFireball/Phenazine/phenazine_hydro.xyz), [phenazine_dehydro.xyz](tests/pyFireball/Phenazine/phenazine_dehydro.xyz)
- [Phenazine/pyrazine_hydro.xyz](tests/pyFireball/Phenazine/pyrazine_hydro.xyz), [pyrazine_dehydro.xyz](tests/pyFireball/Phenazine/pyrazine_dehydro.xyz)
- [two_ribbons.gen](tests/pyFireball/two_ribbons.gen), [two_ribbons.xyz](tests/pyFireball/two_ribbons.xyz)

**Temporary calculation directories:**
- `temp_neb_mol_*/` - NEB image calculations (9 directories)
- `temp_relax_*/` - Relax scan calculations (15 directories)
- `temp_scan_*/` - Rigid scan calculations (11 directories)
- `scan_constrained_*/` - Constrained scan results (multiple subdirectories)

## 4. Dependencies

### Required for extraction:
- **PyQt5** - GUI framework
- **Vispy** - 3D rendering
- **NumPy, Matplotlib** - Numerical/plotting
- **ASE** - Atomic Simulation Environment (for cube reading, vibrations)
- **DFTB+ executable** - External QM package
- **Slater-Koster parameters** - 3ob-3-1 or similar

### Internal FireCore dependencies:
- `pyBall/AtomicSystem.py` - Atomic structure management
- `pyBall/AtomicGraph.py` - Graph representation
- `pyBall/atomicUtils.py` - Atomic utilities
- `pyBall/elements.py` - Element data
- `pyBall/plotUtils.py` - Plotting utilities
- `pyBall/VispyUtils.py` - Vispy helpers
- `pyBall/ExtensionManager.py` - Extension system
- `pyBall/config_utils.py` - Configuration

## 5. Consolidation Plan

### Recommended extraction structure:
```
kekule-neb-repo/
├── src/
│   ├── kekule/
│   │   ├── KekuleBackend.py
│   │   ├── KekuleExplorerGUI.py
│   │   ├── Kekule.py
│   │   └── __init__.py
│   ├── topology/
│   │   ├── GrapheneRibbonBuilder.py
│   │   ├── KekuleSolver.py
│   │   ├── KekuleCLI.py
│   │   └── __init__.py
│   ├── dftb/
│   │   ├── dftb_utils.py
│   │   ├── dftb_lib.py
│   │   └── __init__.py
│   └── neb/
│       ├── neb_h_transfer_molecules.py
│       ├── neb_h_transfer_ribbons.py
│       ├── build_two_ribbons.py
│       └── __init__.py
├── tests/
│   ├── test_kekule_topology.py
│   ├── test_kekule_hbonds.py
│   ├── test_neb_molecules.py
│   └── test_neb_ribbons.py
├── data/
│   ├── molecules/
│   │   ├── phenazine_hydro.xyz
│   │   ├── phenazine_dehydro.xyz
│   │   ├── pyrazine_hydro.xyz
│   │   └── pyrazine_dehydro.xyz
│   └── ribbons/
│       └── two_ribbons.gen
├── docs/
│   ├── Kekule_Topology.md
│   ├── Kekule_Optimizer.md
│   └── DFTB_consolidation_notes.md
└── requirements.txt
```

### Files to leave in FireCore:
- `pyBall/AtomicSystem.py` - General-purpose, used elsewhere
- `pyBall/AtomicGraph.py` - General-purpose
- `pyBall/atomicUtils.py` - General-purpose
- `pyBall/elements.py` - General-purpose
- `pyBall/plotUtils.py` - General-purpose
- `pyBall/VispyUtils.py` - General-purpose
- [pyBall/FireCore.py](pyBall/FireCore.py) - Fireball-specific (separate from DFTB+)

### Refactoring needed (per DFTB_consolidation_notes.md):
1. Move duplicated functions to [dftb_utils.py](pyBall/dftb_utils.py):
   - `load_molecule_xyz()` → [dftb_utils.load_molecule()](pyBall/dftb_utils.py:129-136)
   - [parse_forces()](pyBall/dftb_utils.py:263-284) → [dftb_utils.parse_forces()](pyBall/dftb_utils.py:263-284)
   - [save_xyz_movie()](tests/pyFireball/neb_h_transfer_molecules.py:26-26) → [dftb_utils.save_xyz_movie()](tests/pyFireball/neb_h_transfer_molecules.py:26-26)
   - Unified [run_pbc()](pyBall/dftb_utils.py:309-339) function

2. Extract ribbon builders to shared module:
   - [build_ribbon()](pyBall/KekuleBackend.py:1876-1909), [build_two_ribbon_cell()](tests/pyFireball/build_two_ribbons.py:28-106) → `ribbon_utils.py`

3. Extend `plotUtils.py`:
   - `plot_cube_slice()` for orbital visualization
   - `plotGeometryWithForces()` for geometry+force plots
   - `plot_scan_profile()` for energy profiles

## 6. Key Features to Preserve

- **Kekule Explorer:** Interactive GUI for building/editing Kekule structures
- **sp3 H orientation:** Recent fix for tetrahedral hydrogen placement
- **H-bond detection:** PBC-aware X-H...Y detection for molecules
- **NEB modes:** rigid_scan, relax_scan, dry_run for both molecules and ribbons
- **PBC support:** k-point sampling for periodic ribbon systems
- **DFTB+ integration:** Both subprocess and C-API interfaces
- **Orbital visualization:** Waveplot integration for molecular orbitals

---

## See Also

- [Topical Audit Index](topical_audit.md) — priority ranking, dependency graph, missing topics
- [Molecular Topology Editors](molecular_topology_editors.md) — KekuleExplorerGUI architecture, crystal building gap analysis
- [Molecular Topology Types](molecular_topology_types.md) — atom type assignment used in Kekule structures
- [AFM/STM Simulation](afm_stm_simulation.md) — DFTB+ integration details, STM orbital projection
- [GUI Feature Audit](gui_audit.md) — visualization & editor feature matrices, VisPy consolidation plan