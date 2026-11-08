# USER

Can you review this problem of Kekule structures and Nudged_Elstic Band and other such calculations with DFTB, I would like to extract them to different repo in some consolidated form. What exactly we need to extract? which python modules, and also where are the results so I can review it ? I think we have quite many script related to this Kekule ahdn hydrogen tranfer between nitrogens, both for molecules and for ribbons (extended systems with k-points)

do not code just make audit of existing code and results

---

# SWE-1.6

# Kekule/NEB/H-transfer Code Audit

## 1. Core Python Modules to Extract

### Kekule Structure Explorer
- **[pyBall/KekuleBackend.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/KekuleBackend.py:0:0-0:0)** (1975 lines)
  - Backend for Kekule Structure Explorer
  - Manages AtomicSystem state with hexagonal grid metadata
  - Hydrogen passivation logic (sp3 H orientation fix recently applied)
  - Functions: [honeycomb_ring_nodes()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/KekuleBackend.py:31:0-50:16), [_calc_h_directions_atom()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/KekuleBackend.py:822:4-900:17), [adjust_h()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/KekuleBackend.py:808:4-820:78)

- **[pyBall/KekuleExplorerGUI.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/KekuleExplorerGUI.py:0:0-0:0)** (1413 lines)
  - PyQt5 GUI for Kekule explorer
  - Vispy-based 3D rendering
  - Extension system (ExtensionManager)

- **[pyBall/Kekule.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/Kekule.py:0:0-0:0)**
  - Core Kekule module (needs inspection)

### Kekule Topology Analysis (doc/Topics/Kekule_Topology/)
- **[GrapheneRibbonBuilder.py](cci:7://file:///home/prokophapala/git/FireCore/doc/Topics/Kekule_Topology/GrapheneRibbonBuilder.py:0:0-0:0)** (19733 lines) - Build graphene ribbons with various passivations
- **[KekuleSolver.py](cci:7://file:///home/prokophapala/git/FireCore/doc/Topics/Kekule_Topology/KekuleSolver.py:0:0-0:0)** (5994 lines) - Kekule pattern solver
- **[KekuleCLI.py](cci:7://file:///home/prokophapala/git/FireCore/doc/Topics/Kekule_Topology/KekuleCLI.py:0:0-0:0)** (2635 lines) - Command-line interface
- **[KekuleGUI.py](cci:7://file:///home/prokophapala/git/FireCore/doc/Topics/Kekule_Topology/KekuleGUI.py:0:0-0:0)** (8701 lines) - GUI for Kekule topology
- **[Kekule_Optimizer.md](cci:7://file:///home/prokophapala/git/FireCore/doc/Topics/Kekule_Topology/Kekule_Optimizer.md:0:0-0:0)** (21756 lines) - Documentation
- **[Kekule_Topology.md](cci:7://file:///home/prokophapala/git/FireCore/doc/Topics/Kekule_Topology/Kekule_Topology.md:0:0-0:0)** (39721 lines) - Theoretical background (Jackiw-Rebbi model, grain boundaries)

### DFTB Integration
- **[pyBall/dftb_utils.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/dftb_utils.py:0:0-0:0)** (1294 lines)
  - Subprocess-based DFTB+ interface
  - HSD generation: [makeDFTBjob()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/dftb_utils.py:607:0-700:37), [makeDFTBjob_pbc()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/dftb_utils.py:702:0-801:24)
  - Waveplot integration: [run_waveplot()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/dftb_utils.py:939:0-1016:16), [read_cube()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/dftb_utils.py:1019:0-1038:22)
  - Hessian/vibrations: [get_hessian_ase()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/dftb_utils.py:830:0-864:23), [hessian_to_frequencies()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/dftb_utils.py:890:0-934:29)

- **[pyBall/FireCore.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/FireCore.py:0:0-0:0)** (1488 lines)
  - C bindings for Fireball/DFTB
  - PBC support: `firecore_set_lvs()`, `firecore_set_kpoints()`
  - k-point sampling for periodic systems

## 2. H-transfer NEB Scripts

### Molecular Systems (gamma-point)
- **[tests/pyFireball/neb_h_transfer_molecules.py](cci:7://file:///home/prokophapala/git/FireCore/tests/pyFireball/neb_h_transfer_molecules.py:0:0-0:0)** (967 lines)
  - NEB for H-transfer between hydro/dehydro molecules (phenazine, pyrazine)
  - Modes: `rigid_scan`, `relax_scan`, `dry_run`
  - H-bond detection: [find_hydrogen_bonds()](cci:1://file:///home/prokophapala/git/FireCore/tests/pyFireball/neb_h_transfer_molecules.py:230:0-284:17)
  - PBC-aware placement: [build_molecular_cell()](cci:1://file:///home/prokophapala/git/FireCore/tests/pyFireball/neb_h_transfer_molecules.py:75:0-175:35) with NN-distance mode
  - Input files: [Phenazine/phenazine_hydro.xyz](cci:7://file:///home/prokophapala/git/FireCore/tests/pyFireball/Phenazine/phenazine_hydro.xyz:0:0-0:0), [phenazine_dehydro.xyz](cci:7://file:///home/prokophapala/git/FireCore/tests/pyFireball/Phenazine/phenazine_dehydro.xyz:0:0-0:0), `pyrazine_*.xyz`

### Ribbon Systems (periodic with k-points)
- **[tests/pyFireball/neb_h_transfer.py](cci:7://file:///home/prokophapala/git/FireCore/tests/pyFireball/neb_h_transfer.py:0:0-0:0)** (329 lines)
  - NEB for H-transfer between N-passivated and NH-passivated ribbons
  - Periodic along x with k-point sampling (nk_x=16)
  - Uses GenFormat for DFTB+ input

- **[tests/pyFireball/build_two_ribbons.py](cci:7://file:///home/prokophapala/git/FireCore/tests/pyFireball/build_two_ribbons.py:0:0-0:0)** (173 lines)
  - Build two-ribbon unit cell for H-transfer
  - Uses `GrapheneRibbonBuilder` from Kekule topology module
  - N and NH passivation

- **[tests/pyFireball/scan_ribbon_dftb.py](cci:7://file:///home/prokophapala/git/FireCore/tests/pyFireball/scan_ribbon_dftb.py:0:0-0:0)** (245 lines)
  - Lattice scan for zigzag ribbons using DFTB+
  - PBC along x with k-points

- **[tests/pyFireball/scan_constrained.py](cci:7://file:///home/prokophapala/git/FireCore/tests/pyFireball/scan_constrained.py:0:0-0:0)** (453 lines)
  - Constrained H-transfer scan (move H between ribbons)
  - Fix N atoms, relax others

- **[tests/pyFireball/scan_LHb.py](cci:7://file:///home/prokophapala/git/FireCore/tests/pyFireball/scan_LHb.py:0:0-0:0)** (300 lines)
  - Scan energy vs H-bond length

## 3. Test Scripts and Results

### Kekule Explorer Tests
- **[tests/tKekuleExplorer/test_kekule_topology.py](cci:7://file:///home/prokophapala/git/FireCore/tests/tKekuleExplorer/test_kekule_topology.py:0:0-0:0)** (16177 lines)
- **[tests/tKekuleExplorer/test_kekule_hbonds.py](cci:7://file:///home/prokophapala/git/FireCore/tests/tKekuleExplorer/test_kekule_hbonds.py:0:0-0:0)** (4456 lines)
- **[tests/tKekuleExplorer/test_kekule_relax.py](cci:7://file:///home/prokophapala/git/FireCore/tests/tKekuleExplorer/test_kekule_relax.py:0:0-0:0)** (2520 lines)
- **[tests/tKekuleExplorer/test_ribbon_parity.py](cci:7://file:///home/prokophapala/git/FireCore/tests/tKekuleExplorer/test_ribbon_parity.py:0:0-0:0)** (22665 lines)

**Results directories:**
- [tests/tKekuleExplorer/out_topology/](cci:9://file:///home/prokophapala/git/FireCore/tests/tKekuleExplorer/out_topology:0:0-0:0) - 18 geometry snapshots (benzene, naphthalene, anthracene, etc.)
- [tests/tKekuleExplorer/out_relax/](cci:9://file:///home/prokophapala/git/FireCore/tests/tKekuleExplorer/out_relax:0:0-0:0) - DFTB relaxation results

### DFTB Test Scripts
- **`tests/dftb/example_hessian.py`** (283 lines) - Hessian calculation
- **`tests/dftb/example_orbitals.py`** (501 lines) - Orbital visualization
- **`tests/dftb/test_*.py`** - C-API tests

### PyFireball Results
**Energy profiles:**
- [neb_energy_profile_mol.png](cci:7://file:///home/prokophapala/git/FireCore/tests/pyFireball/neb_energy_profile_mol.png:0:0-0:0) - NEB energy profile for molecules
- [relax_scan_energy.png](cci:7://file:///home/prokophapala/git/FireCore/tests/pyFireball/relax_scan_energy.png:0:0-0:0) - Relax scan results
- [rigid_scan_energy.png](cci:7://file:///home/prokophapala/git/FireCore/tests/pyFireball/rigid_scan_energy.png:0:0-0:0) - Rigid scan results
- `LHb_scan*.png` - H-bond length scans
- `lattice_scan_*.png` - Lattice optimization scans

**Geometry movies:**
- [neb_results_mol.npz](cci:7://file:///home/prokophapala/git/FireCore/tests/pyFireball/neb_results_mol.npz:0:0-0:0) - NEB trajectory data
- [relax_scan_movie.xyz](cci:7://file:///home/prokophapala/git/FireCore/tests/pyFireball/relax_scan_movie.xyz:0:0-0:0) - Relax scan trajectory
- [rigid_scan_movie.xyz](cci:7://file:///home/prokophapala/git/FireCore/tests/pyFireball/rigid_scan_movie.xyz:0:0-0:0) - Rigid scan trajectory
- [dry_run_movie.xyz](cci:7://file:///home/prokophapala/git/FireCore/tests/pyFireball/dry_run_movie.xyz:0:0-0:0) - Geometry-only check

**Input files:**
- [Phenazine/phenazine_hydro.xyz](cci:7://file:///home/prokophapala/git/FireCore/tests/pyFireball/Phenazine/phenazine_hydro.xyz:0:0-0:0), [phenazine_dehydro.xyz](cci:7://file:///home/prokophapala/git/FireCore/tests/pyFireball/Phenazine/phenazine_dehydro.xyz:0:0-0:0)
- [Phenazine/pyrazine_hydro.xyz](cci:7://file:///home/prokophapala/git/FireCore/tests/pyFireball/Phenazine/pyrazine_hydro.xyz:0:0-0:0), [pyrazine_dehydro.xyz](cci:7://file:///home/prokophapala/git/FireCore/tests/pyFireball/Phenazine/pyrazine_dehydro.xyz:0:0-0:0)
- [two_ribbons.gen](cci:7://file:///home/prokophapala/git/FireCore/tests/pyFireball/two_ribbons.gen:0:0-0:0), [two_ribbons.xyz](cci:7://file:///home/prokophapala/git/FireCore/tests/pyFireball/two_ribbons.xyz:0:0-0:0)

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
- [pyBall/FireCore.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/FireCore.py:0:0-0:0) - Fireball-specific (separate from DFTB+)

### Refactoring needed (per DFTB_consolidation_notes.md):
1. Move duplicated functions to [dftb_utils.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/dftb_utils.py:0:0-0:0):
   - `load_molecule_xyz()` → [dftb_utils.load_molecule()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/dftb_utils.py:129:0-136:133)
   - [parse_forces()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/dftb_utils.py:263:0-284:17) → [dftb_utils.parse_forces()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/dftb_utils.py:263:0-284:17)
   - [save_xyz_movie()](cci:1://file:///home/prokophapala/git/FireCore/tests/pyFireball/neb_h_transfer_molecules.py:26:0-26:57) → [dftb_utils.save_xyz_movie()](cci:1://file:///home/prokophapala/git/FireCore/tests/pyFireball/neb_h_transfer_molecules.py:26:0-26:57)
   - Unified [run_pbc()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/dftb_utils.py:309:0-339:25) function

2. Extract ribbon builders to shared module:
   - [build_ribbon()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/KekuleBackend.py:1876:0-1909:31), [build_two_ribbon_cell()](cci:1://file:///home/prokophapala/git/FireCore/tests/pyFireball/build_two_ribbons.py:28:0-106:35) → `ribbon_utils.py`

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