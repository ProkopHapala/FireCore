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

# TOPIC HIERARCHY

## AFM (Atomic Force Microscopy) Simulation
- GPU Rigid Body Dynamics, CPU GridFF Relaxation, and Interactive GUI
- PyOpenCL System: Morse/LJ Path and FDBM Density-Based Path
- Rigid Body Dynamics System
- FDBM Pipeline: DFTB Backend & pySCF Integration Points
- DFTB Reference Calculation & FDBM AFM Forcefield Comparison System

## STM (Scanning Tunneling Microscopy) Simulation
- GPU Green's Function STM Implementation: Current Orbital Projection System & Planned GF Solver Integration
- STM Simulation Pipeline: Orbital Projection & Quantum Transport
- STM QMMM: Fireball DFTB Integration with GPU Density Projection

## GridFF (Grid Force Field) & Surface Interactions
- Interactive GridFF Scanning: PTCDA-on-CaF2 Constrained Relaxation System
- Surface Potential Evaluation: GridFF B-spline and XYZ Rigid Kernels
- Molecule-on-Surface Systems: GridFF, XYZ Scanning, Surface Sampling, and Assembly
- Molecule-Substrate Interaction Energy Scanning: Assembly, GUI Placement, Force Fields & Surface Evaluation
- FoldedAtomicFunctions: Surface Potential Basis Fitting System

## Kekule Topology & H-transfer NEB
- Kekule Structure Explorer: PyQt5/Vispy GUI for building/editing Kekule patterns
- Graphene Ribbon Builder: Zigzag ribbons with N/NH passivation for H-bond studies
- H-transfer NEB: Molecular (gamma-point) and ribbon (k-point) hydrogen transfer calculations
- DFTB+ Integration: Subprocess and C-API interfaces with PBC support
- See [Htransfer_Kekule_DFTB.md](Htransfer_Kekule_DFTB.md) for detailed file inventory and consolidation plan

## Classical Force Fields (MMFF/UFF)
- FireCore Classical Forcefields: MMFFsp3 & UFF (CPU/GPU/Python)
- MMFF/UFF CPU vs GPU Testing: C++ OpenCL and PyOpenCL Parity Infrastructure
- FitREQ_PN: Hydrogen-Bond Parameter Fitting System
- FitREQ Interactive GUI: Monte Carlo Optimization & Energy Decomposition Integration
- FitREQ Hydrogen Bond Fitting System - GPU-Accelerated Parameter Optimization
- **Hierarchical Audit Documents:**
  - [forcefields_overview.md](forcefields_overview.md) — High-level taxonomy of all force field classes
  - [intramolecular_forcefields.md](intramolecular_forcefields.md) — UFF, MMFFsp3, ProjectiveDynamics, XPBD, RigidBody
  - [nonbonding_forcefields.md](nonbonding_forcefields.md) — NBFF, exclusion schemes, FMM
  - [surface_interactions.md](surface_interactions.md) — GridFF, FoldedAtomicFunctions
  - [forcefields_web_implementation.md](forcefields_web_implementation.md) — WebGL/WebGPU shader implementations
  - [Forcefields_Audit.md](Forcefields_Audit.md) — Consolidated comprehensive audit (legacy single-file version)

## Nanocrystal Vibrations
- Silicon/Diamond nanocrystal generation, force field setup, and vibration spectroscopy
- See [Nanocrystal_Vibrations.md](Nanocrystal_Vibrations.md) for detailed file inventory and workflow

## Rigid Body Dynamics
- Rigid Body Dynamics on Surfaces (pyOpenCL)
- Rigid Body Dynamics System for AFM Simulation

## DFTB / QM Integration
- Fireball Hamiltonian Assembly: Fortran → PyOpenCL
- DFTB+ Python Integration: Library Interfaces, Parsers, and OpenCL Grid Projection
- DFTB+ Eigenvector Export for OpenCL Orbital Projection
- DFTB+ Calculation Flow: Standalone Program, C API, and Python Wrapper
- DFTB Reference Calculation & FDBM AFM Forcefield Comparison System

## Molecular Dynamics
- XPBD Molecular Dynamics pyOpenCL

## Web & Visualization
- WebGPU Molecular Visualization & Physics Simulation

## Debugging & Utilities
- FFutils Visual Debugging Consolidation: XPBD_2D, XPDB_AVBD, and C++ Buffer Systems

---

# CODEMAPS

* [All Windsurf Codempas](https://windsurf.com/codemaps)

* [AFM Simulation: GPU Rigid Body Dynamics, CPU GridFF Relaxation, and Interactive GUI]
(https://windsurf.com/codemaps/594f7eaf-c3ab-4139-8f20-d1d2d7f8d401-fe86ab10a43f3d18)
* [GPU Green's Function STM Implementation: Current Orbital Projection System & Planned GF Solver Integration]
(https://windsurf.com/codemaps/f398c2cf-5ff8-4d75-a398-c83e788e27b4-fe86ab10a43f3d18)
* [STM Simulation Pipeline: Orbital Projection & Quantum Transport]
(https://windsurf.com/codemaps/d0242216-c415-4f38-98f9-4c88b5dfeeb8-fe86ab10a43f3d18)
* [STM QMMM: Fireball DFTB Integration with GPU Density Projection]
(https://windsurf.com/codemaps/9fa40c64-e78c-42f2-9573-574936c8040d-fe86ab10a43f3d18)
* [Interactive GridFF Scanning: PTCDA-on-CaF2 Constrained Relaxation System]
(https://windsurf.com/codemaps/99d506e2-223b-4ae7-bb60-8c2498fedfb9-8796fe608a7d71c1)
* [FireCore Classical Forcefields: MMFFsp3 & UFF (CPU/GPU/Python)]
(https://windsurf.com/codemaps/53f2fe2c-ac5c-4c0b-b905-af6653adde97-8796fe608a7d71c1)
* [MMFF/UFF CPU vs GPU Testing: C++ OpenCL and PyOpenCL Parity Infrastructure]
(https://windsurf.com/codemaps/8d1b056f-1502-4363-b52d-8257de4be453-8796fe608a7d71c1)
* [Surface Potential Evaluation: GridFF B-spline and XYZ Rigid Kernels]
(https://windsurf.com/codemaps/2a639fae-c9cb-407a-9d45-7b806c90c749-8796fe608a7d71c1)
* [Rigid Body Dynamics on Surfaces (pyOpenCL)]
(https://windsurf.com/codemaps/b5d9c2d2-50f0-4ba7-bc65-60db6e06e423-8796fe608a7d71c1)
* [FoldedAtomicFunctions: Surface Potential Basis Fitting System]
(https://windsurf.com/codemaps/c9fc44a7-57a2-47c5-906f-886fa301ccc7-8796fe608a7d71c1)
* [Molecule-Substrate Interaction Energy Scanning: Assembly, GUI Placement, Force Fields & Surface Evaluation]
(https://windsurf.com/codemaps/38bd3cb6-31c0-45b6-9e09-fda94257999c-8796fe608a7d71c1)
* [Molecule-on-Surface Systems: GridFF, XYZ Scanning, Surface Sampling, and Assembly]
(https://windsurf.com/codemaps/f8407e23-3a2e-41f1-abcf-9c15f3644c41-8796fe608a7d71c1)
* [AFM PyOpenCL System: Morse/LJ Path and FDBM Density-Based Path]
(https://windsurf.com/codemaps/9bb4c2a5-0c38-4943-abe9-254cfdcc75af-8796fe608a7d71c1)
* [AFM Simulation: GPU Rigid Body Dynamics, CPU GridFF Relaxation, and Interactive GUI]
(https://windsurf.com/codemaps/594f7eaf-c3ab-4139-8f20-d1d2d7f8d401-fe86ab10a43f3d18)
* [GPU Green's Function STM Implementation: Current Orbital Projection System & Planned GF Solver Integration]
(https://windsurf.com/codemaps/f398c2cf-5ff8-4d75-a398-c83e788e27b4-fe86ab10a43f3d18)
* [STM Simulation Pipeline: Orbital Projection & Quantum Transport]
(https://windsurf.com/codemaps/d0242216-c415-4f38-98f9-4c88b5dfeeb8-fe86ab10a43f3d18)
* [STM QMMM: Fireball DFTB Integration with GPU Density Projection]
(https://windsurf.com/codemaps/9fa40c64-e78c-42f2-9573-574936c8040d-fe86ab10a43f3d18)
* [Interactive GridFF Scanning: PTCDA-on-CaF2 Constrained Relaxation System]
(https://windsurf.com/codemaps/99d506e2-223b-4ae7-bb60-8c2498fedfb9-8796fe608a7d71c1)
* [FireCore Classical Forcefields: MMFFsp3 & UFF (CPU/GPU/Python)]
(https://windsurf.com/codemaps/53f2fe2c-ac5c-4c0b-b905-af6653adde97-8796fe608a7d71c1)
* [MMFF/UFF CPU vs GPU Testing: C++ OpenCL and PyOpenCL Parity Infrastructure]
(https://windsurf.com/codemaps/8d1b056f-1502-4363-b52d-8257de4be453-8796fe608a7d71c1)
* [Surface Potential Evaluation: GridFF B-spline and XYZ Rigid Kernels]
(https://windsurf.com/codemaps/2a639fae-c9cb-407a-9d45-7b806c90c749-8796fe608a7d71c1)
* [Rigid Body Dynamics on Surfaces (pyOpenCL)]
(https://windsurf.com/codemaps/b5d9c2d2-50f0-4ba7-bc65-60db6e06e423-8796fe608a7d71c1)
* [FoldedAtomicFunctions: Surface Potential Basis Fitting System]
(https://windsurf.com/codemaps/c9fc44a7-57a2-47c5-906f-886fa301ccc7-8796fe608a7d71c1)
* [Molecule-Substrate Interaction Energy Scanning: Assembly, GUI Placement, Force Fields & Surface Evaluation]
(https://windsurf.com/codemaps/38bd3cb6-31c0-45b6-9e09-fda94257999c-8796fe608a7d71c1)
* [Molecule-on-Surface Systems: GridFF, XYZ Scanning, Surface Sampling, and Assembly]
(https://windsurf.com/codemaps/f8407e23-3a2e-41f1-abcf-9c15f3644c41-8796fe608a7d71c1)
* [AFM PyOpenCL System: Morse/LJ Path and FDBM Density-Based Path]
(https://windsurf.com/codemaps/9bb4c2a5-0c38-4943-abe9-254cfdcc75af-8796fe608a7d71c1)
* [FitREQ_PN: Hydrogen-Bond Parameter Fitting System]
(https://windsurf.com/codemaps/d977d597-94b4-42c3-a92a-0cefe34a3e82-8796fe608a7d71c1)
* [Fireball Hamiltonian Assembly: Fortran → PyOpenCL]
(https://windsurf.com/codemaps/92089c9f-b536-4b78-955a-915f4363f656-8796fe608a7d71c1)
* [XPBD Molecular Dynamics pyOpenCL]
(https://windsurf.com/codemaps/2e558e51-fdbe-4bd4-8732-7818724d4ced-8796fe608a7d71c1)
* [AFM FDBM Pipeline: DFTB Backend & pySCF Integration Points]
(https://windsurf.com/codemaps/02d559c9-de47-4058-b07b-3318664b454e-fe86ab10a43f3d18)
* [FitREQ Interactive GUI: Monte Carlo Optimization & Energy Decomposition Integration]
(https://windsurf.com/codemaps/e25a0dfc-f9a8-42ab-b8bb-1d959037ca68-fe86ab10a43f3d18)
* [FitREQ Hydrogen Bond Fitting System - GPU-Accelerated Parameter Optimization]
(https://windsurf.com/codemaps/bf59a960-ac6c-4eea-b828-9bd18c3d44ac-fe86ab10a43f3d18)
* [DFTB Reference Calculation & FDBM AFM Forcefield Comparison System]
(https://windsurf.com/codemaps/1153fe89-ff29-4d4b-b4a6-e97d8f37047f-fe86ab10a43f3d18)
* [DFTB+ Python Integration: Library Interfaces, Parsers, and OpenCL Grid Projection]
(https://windsurf.com/codemaps/1d6b4b7c-04de-49ef-b581-12cf5bfef54a-fe86ab10a43f3d18)
* [DFTB+ Eigenvector Export for OpenCL Orbital Projection]
(https://windsurf.com/codemaps/845d1373-d23e-4f7d-a109-c0d8eccebea9-fe86ab10a43f3d18)
* [DFTB+ Calculation Flow: Standalone Program, C API, and Python Wrapper]
(https://windsurf.com/codemaps/2c157118-9d28-4a7c-a234-a49a3d464424-fe86ab10a43f3d18)
* [AFM PyOpenCL System: Morse/LJ Path and FDBM Density-Based Path]
(https://windsurf.com/codemaps/9bb4c2a5-0c38-4943-abe9-254cfdcc75af-fe86ab10a43f3d18)
* [Rigid Body Dynamics System for AFM Simulation]
(https://windsurf.com/codemaps/c9f13e1f-edfa-4702-814f-5036d03ea6c9-fe86ab10a43f3d18)
* [AFM Simulation: GPU Rigid Body Dynamics, CPU GridFF Relaxation, and Interactive GUI]
(https://windsurf.com/codemaps/594f7eaf-c3ab-4139-8f20-d1d2d7f8d401-fe86ab10a43f3d18)
* [GPU Green's Function STM Implementation: Current Orbital Projection System & Planned GF Solver Integration]
(https://windsurf.com/codemaps/f398c2cf-5ff8-4d75-a398-c83e788e27b4-fe86ab10a43f3d18)
* [STM Simulation Pipeline: Orbital Projection & Quantum Transport]
(https://windsurf.com/codemaps/d0242216-c415-4f38-98f9-4c88b5dfeeb8-fe86ab10a43f3d18)
* [STM QMMM: Fireball DFTB Integration with GPU Density Projection]
(https://windsurf.com/codemaps/9fa40c64-e78c-42f2-9573-574936c8040d-fe86ab10a43f3d18)
* [Interactive GridFF Scanning: PTCDA-on-CaF2 Constrained Relaxation System]
(https://windsurf.com/codemaps/99d506e2-223b-4ae7-bb60-8c2498fedfb9-8796fe608a7d71c1)
* [FireCore Classical Forcefields: MMFFsp3 & UFF (CPU/GPU/Python)]
(https://windsurf.com/codemaps/53f2fe2c-ac5c-4c0b-b905-af6653adde97-8796fe608a7d71c1)
* [MMFF/UFF CPU vs GPU Testing: C++ OpenCL and PyOpenCL Parity Infrastructure]
(https://windsurf.com/codemaps/8d1b056f-1502-4363-b52d-8257de4be453-8796fe608a7d71c1)
* [Surface Potential Evaluation: GridFF B-spline and XYZ Rigid Kernels]
(https://windsurf.com/codemaps/2a639fae-c9cb-407a-9d45-7b806c90c749-8796fe608a7d71c1)
* [Rigid Body Dynamics on Surfaces (pyOpenCL)]
(https://windsurf.com/codemaps/b5d9c2d2-50f0-4ba7-bc65-60db6e06e423-8796fe608a7d71c1)
* [FoldedAtomicFunctions: Surface Potential Basis Fitting System]
(https://windsurf.com/codemaps/c9fc44a7-57a2-47c5-906f-886fa301ccc7-8796fe608a7d71c1)
* [Molecule-Substrate Interaction Energy Scanning: Assembly, GUI Placement, Force Fields & Surface Evaluation]
(https://windsurf.com/codemaps/38bd3cb6-31c0-45b6-9e09-fda94257999c-8796fe608a7d71c1)
* [Molecule-on-Surface Systems: GridFF, XYZ Scanning, Surface Sampling, and Assembly]
(https://windsurf.com/codemaps/f8407e23-3a2e-41f1-abcf-9c15f3644c41-8796fe608a7d71c1)
* [AFM PyOpenCL System: Morse/LJ Path and FDBM Density-Based Path]
(https://windsurf.com/codemaps/9bb4c2a5-0c38-4943-abe9-254cfdcc75af-8796fe608a7d71c1)
* [FitREQ_PN: Hydrogen-Bond Parameter Fitting System]
(https://windsurf.com/codemaps/d977d597-94b4-42c3-a92a-0cefe34a3e82-8796fe608a7d71c1)
* [Fireball Hamiltonian Assembly: Fortran → PyOpenCL]
(https://windsurf.com/codemaps/92089c9f-b536-4b78-955a-915f4363f656-8796fe608a7d71c1)
* [XPBD Molecular Dynamics pyOpenCL]
(https://windsurf.com/codemaps/2e558e51-fdbe-4bd4-8732-7818724d4ced-8796fe608a7d71c1)
* [WebGPU Molecular Visualization & Physics Simulation]
(https://windsurf.com/codemaps/65e0669d-bbfa-4f58-87f9-04050f2cdced-8796fe608a7d71c1)
* [FFutils Visual Debugging Consolidation: XPBD_2D, XPDB_AVBD, and C++ Buffer Systems]
(https://windsurf.com/codemaps/cc3069ab-a83a-4948-8daa-39dbd7d6464f-8796fe608a7d71c1)
* [Rigid Body Dynamics System for AFM Simulation]
(https://windsurf.com/codemaps/c9f13e1f-edfa-4702-814f-5036d03ea6c9-fe86ab10a43f3d18)

---
























































