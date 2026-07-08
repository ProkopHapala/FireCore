# Topical Audit Manifesto

## Purpose

FireCore is a >5 year computational science project that has evolved as a lab-book/journal rather than a polished product. Many topics have been started, abandoned, restarted, or implemented multiple times across different languages and frameworks. This folder serves as a consolidation effort to map and organize these scattered implementations.

## Problem Statement

The repository contains multiple implementations of the same scientific concepts across:
- **Languages**: Fortran, C++, Python, JavaScript, Julia
- **Acceleration frameworks**: OpenCL, CUDA, WebGPU, WebGL
- **Interfaces**: CLI, GUI (SDL), web-based
- **Abstraction levels**: Low-level kernels, high-level APIs, educational scripts

Without systematic documentation, it is difficult to:
- Find all implementations of a given algorithm
- Understand which implementation is current/active vs experimental
- Identify dependencies between scattered modules
- Avoid reimplementing existing functionality

## Goals

1. **Map scattered implementations** - Create topic-based indexes showing where each concept is implemented across the codebase
2. **Identify active vs experimental code** - Mark which implementations are production-ready, deprecated, or unfinished
3. **Document cross-language relationships** - Explain how C++ libraries, Python bindings, and Fortran modules relate to each other
4. **Consolidate duplicate efforts** - Identify when the same functionality exists in multiple places and recommend consolidation
5. **Preserve experimental work** - Keep track of promising but unfinished experiments for future reference

## Structure

Each markdown file in this folder should focus on a specific scientific topic or computational method, such as:

- `gridff.md` - Grid Force Field implementations (C++, OpenCL, Python, JavaScript)
- `afm_simulation.md` - AFM surface science simulation modules
- `qm_mm.md` - QM/MM coupling across Fortran/C++/Python
- `molecular_dynamics.md` - MD implementations in various frameworks
- `surface_adsorption.md` - Molecule-surface interaction calculations
- `kpoints_bands.md` - k-point sampling and band structure (recently added)

## Content Guidelines

For each topic, document:

- **Overview** - What is the scientific problem being solved?
- **Implementations** - List all locations where this topic appears (file paths, languages)
- **Status** - Active, experimental, deprecated, or unfinished
- **Relationships** - How different implementations relate (e.g., Python wrapper calling C++ library)
- **Notes** - Context about why multiple implementations exist, which to use, consolidation plans

## Usage

When working on a feature:
1. Check this folder first to see if related work already exists
2. Update the relevant topic file when adding new implementations
3. Mark old implementations as deprecated when superseded
4. Add new topic files for novel scientific areas

## Scope

This is **not** a replacement for:
- API documentation (use Doxygen/docstrings)
- Build instructions (use README/CMakeLists.txt)
- Tutorial examples (use `/tests/` directory)

This **is** a supplement to:
- Help navigate the complex, multi-language codebase
- Preserve context for long-term project continuity
- Enable systematic consolidation of duplicate efforts

---

# TOPIC HIERARCHY WITH PRIORITY RANKING

Topics are ranked by **consolidation priority** — how much duplicate/scattered code exists and how much would benefit from unification.

## Priority 1 (Critical) — Molecular Topology & GUI Consolidation

**Why critical:** 4+ GUIs with overlapping editor operations, 3 languages with duplicate topology logic, crystal building exists only in JS (no Python equivalent). This is the biggest consolidation opportunity.

### 1a. Molecular Topology — Base Graph & Connectivity
- Graph representations, bond finding, neighbor lists, ring/bridge detection across C++/Python/JavaScript
- **Gap identified:** Crystal building (CIF parsing, lattice replication, slab cutting) exists only in JavaScript `CrystalUtils.js` (1188 lines) — no Python equivalent
- **Audit Document:** [molecular_topology.md](molecular_topology.md)
- **Key files:** `AtomicGraph.py` (392 lines), `AtomicSystem.py` (1314 lines), `EditableMolecule.js` (1057 lines), `MMFFBuilderBase.h` (808 lines)

### 1b. Molecular Topology — Type Assignment
- Atom type assignment (sp1/sp2/sp3), VSEPR geometry, parameter loading
- **Audit Document:** [molecular_topology_types.md](molecular_topology_types.md)

### 1c. Molecular Topology — Editors & GUIs
- Interactive editors: KekuleExplorerGUI (Python/VisPy), molgui_webgpu (JS/WebGPU), MoleculeEditor2D (deprecated), molgui_web (legacy)
- **Consolidation target:** One unified Python VisPy GUI
- **Audit Document:** [molecular_topology_editors.md](molecular_topology_editors.md)
- **Feature audit:** [gui_audit.md](gui_audit.md) — detailed visualization & editor feature matrices
- **Key gap:** Crystal builder (`CrystalUtils.js`) needs Python port for unified GUI

### 1d. Crystal Building (Gap — No Audit Document Yet)
- **Current state:** `CrystalUtils.js` (1188 lines) — CIF parsing, symmetry ops, lattice vectors, frac↔cart, slab cutting (HKL), plane cutting, nanocrystal generation, dedup, bond building across cell boundaries
- **Python gap:** No equivalent. `KekuleBackend` does hex grid for graphene only
- **C++ gap:** `MMFFBuilder` has `insertAtoms`/`autoBonds` but no CIF/crystallography
- **Action needed:** Create `crystal_building.md` audit document; port `CrystalUtils.js` → Python
- **Related files:** `web/molgui_webgpu/CrystalUtils.js`, `web/molgui_webgpu/Nanocrystals.js` (631 lines), `web/molgui_webgpu/BuildersGUI.js` (1130 lines)

## Priority 2 (High) — Classical Force Fields

**Why high:** Multiple FF implementations (MMFF, UFF, XPBD, ProjectiveDynamics) with CPU/GPU parity concerns across 3 languages. Core simulation capability.

### 2a. Intramolecular Force Fields
- UFF, MMFFsp3, ProjectiveDynamics, XPBD/XPDB, RigidBody FF
- **Audit Documents:**
  - [forcefields_overview.md](forcefields_overview.md) — High-level taxonomy
  - [intramolecular_forcefields.md](intramolecular_forcefields.md) — Detailed per-FF analysis
  - [Forcefields_Audit.md](Forcefields_Audit.md) — Legacy consolidated audit
- **Key integrator:** `MolWorld_sp3.h` — orchestrates all FF components

### 2b. Non-Bonding Force Fields
- NBFF (LJ/Morse/Coulomb), FMM, PME, exclusion schemes
- **Audit Document:** [nonbonding_forcefields.md](nonbonding_forcefields.md)

### 2c. Web Force Fields
- XPBD/MMFFL in WebGPU compute shaders, WebGL legacy
- **Audit Document:** [forcefields_web_implementation.md](forcefields_web_implementation.md)

### 2d. Parameter Fitting
- FitREQ: H-bond parameter fitting, Monte Carlo optimization, GPU-accelerated
- **Related:** `FitREQ.md`, `FitREQ_CPU_Tutorial.guide.md` in `doc/DevNotes/`

## Priority 3 (High) — Surface Interactions & GridFF

**Why high:** Core FireCore capability (molecule-on-surface), multiple representations (GridFF, FAF, Ewald2D), well-documented but needs consolidation of GridFF variants.

### 3a. GridFF & Surface Potentials
- GridFF (B-spline interpolation), FoldedAtomicFunctions (compact basis), Surface.cl (unified kernel), Ewald2D
- **Audit Document:** [surface_interactions.md](surface_interactions.md)
- **Consolidation needed:** `GridFF.py` vs `GridFF_new.py` vs `GridFFRelaxedScan.py`

## Priority 4 (Medium) — AFM/STM Simulation

**Why medium:** Well-organized pipeline (ModularAFMPipeline), but depends on many other topics (GridFF, RigidBody, DFTB, STM). Consolidation mostly about interface cleanup.

### 4a. AFM Simulation
- AFMulator, RigidBodyAFM, ModularAFMPipeline (S1-S6 stages), AFMExtension GUI integration
- **Audit Document:** [afm_stm_simulation.md](afm_stm_simulation.md)

### 4b. STM Simulation
- Fireball SCF → H/S matrices → spectral function → DOS/STM current, NEGF Caroli formula, orbital projection
- **Audit Document:** [afm_stm_simulation.md](afm_stm_simulation.md) (sections 7-8)
- **Critical conventions:** Orbital ordering (Fortran vs OpenCL Hamiltonian vs OpenCL Grid) — see [memory ccfa7062]

## Priority 5 (Medium) — QM Integration (DFTB/Fireball)

**Why medium:** Specialized interfaces, not duplicated, but complex integration points.

- Fireball Hamiltonian Assembly: Fortran → PyOpenCL
- DFTB+ Integration: subprocess, C-API, parsers, OpenCL grid projection
- **Related audit:** [afm_stm_simulation.md](afm_stm_simulation.md) (DFTB sections)

## Priority 6 (Medium) — Crystallography & Nanocrystals

**Why medium:** Niche functionality but important for surface science workflows.

- Nanocrystal generation, force field setup, vibration spectroscopy
- **Audit Document:** [Nanocrystal_Vibrations.md](Nanocrystal_Vibrations.md)
- **Related:** `Nanocrystals.js` (631 lines) — JS nanocrystal builder using CrystalUtils

## Priority 7 (Low) — Kekule Topology & H-transfer NEB

**Why low:** Specialized research workflow, not duplicated.

- Kekule Structure Explorer, Graphene Ribbon Builder, H-transfer NEB
- **Audit Document:** [Htransfer_Kekule_DFTB.md](Htransfer_Kekule_DFTB.md)

## Priority 8 (Low) — Debugging & Utilities

- FFutils visual debugging consolidation
- XPBD_2D, XPDB_AVBD, C++ buffer systems

---

# TOPIC DEPENDENCY GRAPH

```
                    ┌─────────────────────────────────┐
                    │  Molecular Topology (1a, 1b)     │
                    │  AtomicGraph / EditableMolecule  │
                    └──────────┬──────────┬────────────┘
                               │          │
                    ┌──────────▼──┐  ┌────▼──────────────┐
                    │  Editors/GUI │  │  Crystal Building │
                    │  (1c, 1d)    │  │  (1d — GAP)       │
                    └──────┬───────┘  └────┬──────────────┘
                           │               │
                    ┌──────▼───────────────▼──────────────┐
                    │  Type Assignment / VSEPR (1b)        │
                    │  → Force Field Parameter Loading     │
                    └──────────────┬───────────────────────┘
                                   │
              ┌────────────────────┼────────────────────┐
              │                    │                    │
    ┌─────────▼────────┐ ┌────────▼─────────┐ ┌────────▼────────┐
    │ Intramolecular FF │ │ Non-Bonding FF   │ │ Surface Interac.│
    │ (2a: MMFF, UFF,  │ │ (2b: NBFF, FMM,  │ │ (3: GridFF, FAF │
    │  XPBD, RigidBody)│ │  PME)            │ │  Ewald2D)       │
    └────────┬─────────┘ └────────┬─────────┘ └────────┬────────┘
             │                    │                    │
             └────────────────────┼────────────────────┘
                                  │
                    ┌─────────────▼──────────────┐
                    │  AFM/STM Simulation (4)    │
                    │  AFMulator → RigidBody →   │
                    │  GridFF → Tip Relaxation → │
                    │  AFM Image / STM Current   │
                    └─────────────┬──────────────┘
                                  │
                    ┌─────────────▼──────────────┐
                    │  QM Integration (5)        │
                    │  Fireball SCF / DFTB+      │
                    │  → H/S matrices → orbitals │
                    └────────────────────────────┘
```

**Key dependencies:**
- Crystal Building (1d) → Topology (1a): crystal generation produces topology
- Topology (1a) → Type Assignment (1b): graph structure determines hybridization
- Type Assignment (1b) → Force Fields (2): types determine FF parameters
- Force Fields (2) + Surface Interactions (3) → AFM/STM (4): AFM uses both
- QM Integration (5) → AFM/STM (4): STM requires Fireball/DFTB orbital data
- Editors/GUI (1c) → All topics: GUI needs to visualize and edit everything

---

# MISSING TOPICS (Not Yet Audited)

| Topic | Why Important | Related Files | Priority |
|-------|--------------|---------------|----------|
| **Crystal Building** | CIF parsing, lattice replication, slab cutting — exists only in JS, no Python port | `CrystalUtils.js` (1188 lines), `Nanocrystals.js` (631 lines), `BuildersGUI.js` (1130 lines) | **Critical** |
| **Molecular Dynamics (MD)** | XPBD MD in PyOpenCL, MD in WebGPU compute shaders, integrators | `pyBall/OCL/MolecularDynamics.py`, `web/molgui_webgpu/` (compute shaders) | High |
| **File I/O & Format Parsing** | XYZ, MOL2, CIF, extended XYZ with Lattice — scattered across languages | `atomicUtils.py`, `CrystalUtils.js`, `IO_utils.h` | Medium |
| **Visualization & Rendering** | VisPy (Python), OpenGL (C++), WebGPU/WebGL (JS) — 3 separate rendering stacks | `KekuleExplorerGUI.py`, `MolGUI.h`, `MoleculeRenderer.js` | Medium |
| **Selection Systems** | 3 different selection implementations (C++ set, Python set, JS Selection class) | `Selection.js`, `AtomicSystem.py`, `MMFFBuilder.h` | Low |

---

# CONSOLIDATION PRIORITY SUMMARY

| # | What | Lines Affected | Effort | Impact |
|---|------|---------------|--------|--------|
| 1 | Port `CrystalUtils.js` → Python `CrystalBuilder.py` | ~1200 new | Medium | High — enables unified GUI |
| 2 | Port VSEPR capping from `EditableMolecule.js` → Python | ~300 new | Low | High — Python GUI parity |
| 3 | Unify Python graph: `AtomicGraph` vs `AtomicSystem` | ~1700 refactor | Medium | High — clear API |
| 4 | Consolidate `EditableMolecule.js` (WebGL vs WebGPU duplicates) | ~1050 dedup | Low | Medium — maintenance |
| 5 | Consolidate `MMParams.js` (WebGL vs WebGPU) | ~500 dedup | Low | Medium — maintenance |
| 6 | Archive `MoleculeEditor2D.py` (deprecated) | ~1600 removal | Trivial | Low — cleanup |
| 7 | Consolidate `GridFF.py` variants (3 files → 1) | ~1000 refactor | Medium | Medium — clarity |
| 8 | Build unified VisPy GUI from KekuleExplorer + crystal builder | ~2000 new | High | High — single GUI |

---

# UNIQUE CODEMAPS (Deduplicated)

* [All Windsurf Codemaps](https://windsurf.com/codemaps)

### AFM/STM Simulation
* [AFM Simulation: GPU Rigid Body Dynamics, CPU GridFF Relaxation, and Interactive GUI](https://windsurf.com/codemaps/594f7eaf-c3ab-4139-8f20-d1d2d7f8d401-fe86ab10a43f3d18)
* [AFM PyOpenCL System: Morse/LJ Path and FDBM Density-Based Path](https://windsurf.com/codemaps/9bb4c2a5-0c38-4943-abe9-254cfdcc75af-8796fe608a7d71c1)
* [AFM FDBM Pipeline: DFTB Backend & pySCF Integration Points](https://windsurf.com/codemaps/02d559c9-de47-4058-b07b-3318664b454e-fe86ab10a43f3d18)
* [Rigid Body Dynamics on Surfaces (pyOpenCL)](https://windsurf.com/codemaps/b5d9c2d2-50f0-4ba7-bc65-60db6e06e423-8796fe608a7d71c1)
* [Rigid Body Dynamics System for AFM Simulation](https://windsurf.com/codemaps/c9f13e1f-edfa-4702-814f-5036d03ea6c9-fe86ab10a43f3d18)

### STM Simulation
* [GPU Green's Function STM Implementation: Current Orbital Projection System & Planned GF Solver Integration](https://windsurf.com/codemaps/f398c2cf-5ff8-4d75-a398-c83e788e27b4-fe86ab10a43f3d18)
* [STM Simulation Pipeline: Orbital Projection & Quantum Transport](https://windsurf.com/codemaps/d0242216-c415-4f38-98f9-4c88b5dfeeb8-fe86ab10a43f3d18)
* [STM QMMM: Fireball DFTB Integration with GPU Density Projection](https://windsurf.com/codemaps/9fa40c64-e78c-42f2-9573-574936c8040d-fe86ab10a43f3d18)

### Surface Interactions
* [Interactive GridFF Scanning: PTCDA-on-CaF2 Constrained Relaxation System](https://windsurf.com/codemaps/99d506e2-223b-4ae7-bb60-8c2498fedfb9-8796fe608a7d71c1)
* [Surface Potential Evaluation: GridFF B-spline and XYZ Rigid Kernels](https://windsurf.com/codemaps/2a639fae-c9cb-407a-9d45-7b806c90c749-8796fe608a7d71c1)
* [FoldedAtomicFunctions: Surface Potential Basis Fitting System](https://windsurf.com/codemaps/c9fc44a7-57a2-47c5-906f-886fa301ccc7-8796fe608a7d71c1)
* [Molecule-Substrate Interaction Energy Scanning: Assembly, GUI Placement, Force Fields & Surface Evaluation](https://windsurf.com/codemaps/38bd3cb6-31c0-45b6-9e09-fda94257999c-8796fe608a7d71c1)
* [Molecule-on-Surface Systems: GridFF, XYZ Scanning, Surface Sampling, and Assembly](https://windsurf.com/codemaps/f8407e23-3a2e-41f1-abcf-9c15f3644c41-8796fe608a7d71c1)

### Force Fields
* [FireCore Classical Forcefields: MMFFsp3 & UFF (CPU/GPU/Python)](https://windsurf.com/codemaps/53f2fe2c-ac5c-4c0b-b905-af6653adde97-8796fe608a7d71c1)
* [MMFF/UFF CPU vs GPU Testing: C++ OpenCL and PyOpenCL Parity Infrastructure](https://windsurf.com/codemaps/8d1b056f-1502-4363-b52d-8257de4be453-8796fe608a7d71c1)
* [FitREQ_PN: Hydrogen-Bond Parameter Fitting System](https://windsurf.com/codemaps/d977d597-94b4-42c3-a92a-0cefe34a3e82-8796fe608a7d71c1)
* [FitREQ Interactive GUI: Monte Carlo Optimization & Energy Decomposition Integration](https://windsurf.com/codemaps/e25a0dfc-f9a8-42ab-b8bb-1d959037ca68-fe86ab10a43f3d18)
* [FitREQ Hydrogen Bond Fitting System - GPU-Accelerated Parameter Optimization](https://windsurf.com/codemaps/bf59a960-ac6c-4eea-b828-9bd18c3d44ac-fe86ab10a43f3d18)

### QM Integration
* [Fireball Hamiltonian Assembly: Fortran → PyOpenCL](https://windsurf.com/codemaps/92089c9f-b536-4b78-955a-915f4363f656-8796fe608a7d71c1)
* [DFTB+ Python Integration: Library Interfaces, Parsers, and OpenCL Grid Projection](https://windsurf.com/codemaps/1d6b4b7c-04de-49ef-b581-12cf5bfef54a-fe86ab10a43f3d18)
* [DFTB+ Eigenvector Export for OpenCL Orbital Projection](https://windsurf.com/codemaps/845d1373-d23e-4f7d-a109-c0d8eccebea9-fe86ab10a43f3d18)
* [DFTB+ Calculation Flow: Standalone Program, C API, and Python Wrapper](https://windsurf.com/codemaps/2c157118-9d28-4a7c-a234-a49a3d464424-fe86ab10a43f3d18)
* [DFTB Reference Calculation & FDBM AFM Forcefield Comparison System](https://windsurf.com/codemaps/1153fe89-ff29-4d4b-b4a6-e97d8f37047f-fe86ab10a43f3d18)

### Molecular Dynamics & Visualization
* [XPBD Molecular Dynamics pyOpenCL](https://windsurf.com/codemaps/2e558e51-fdbe-4bd4-8732-7818724d4ced-8796fe608a7d71c1)
* [WebGPU Molecular Visualization & Physics Simulation](https://windsurf.com/codemaps/65e0669d-bbfa-4f58-87f9-04050f2cdced-8796fe608a7d71c1)

### Debugging
* [FFutils Visual Debugging Consolidation: XPBD_2D, XPDB_AVBD, and C++ Buffer Systems](https://windsurf.com/codemaps/cc3069ab-a83a-4948-8daa-39dbd7d6464f-8796fe608a7d71c1)

---






















































