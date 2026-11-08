# FireCore Hermes Codemap

Auto-generated from `doc/`. Each file has a one-line summary plus top headings.

- Total markdown files: 535


## `doc/DevNotes/`

- `DFTB_consolidation_notes.md` (15 KB)
  - **Summary:** Systematize the DFTB+ Python scripts and modules without breaking any functionality. Move reusable logic into `pyBall/dftb_utils.py` (subpro 
  - **Headings:** DFTB+ Python Code Consolidation Notes,   - Goal,   - 1. Current Inventory,     - 1.1 Library Modules,     - 1.2 Test Scripts,   - 2. Duplication Map,   - 3. Proposed Refactoring Plan,     - 3.1 `pyBall/dftb_utils.py` — New/Modified Functions
    ... (+20 more)

- `FitREQ.md` (37 KB)
  - **Summary:** This document provides a comprehensive overview of the FitREQ system, a library designed for fitting non-covalent interaction parameters bet 
  - **Headings:** FitREQ Molecular Fitting Library Documentation,   - 1. Introduction to FitREQ,   - 2. Physical Principles and Fitting Objective,     - Hydrogen Bond Corrections,   - 3. System Architecture and Class Design,     - 3.1. C++ Core (`cpp/common/molecular/FitREQ.h`),     - 3.2. C-Interface Library (`cpp/libs/Molecular/FitREQ_lib.cpp`),     - 3.3. Python Bindings (`pyBall/FitREQ.py`)
    ... (+23 more)

- `FitREQ_CPU_Tutorial.md` (10 KB)
  - **Summary:** This tutorial explains how to use the CPU version of FitREQ for H-bond 2D scans using the script `tests/tFitREQ/opt_2D.py`. It focuses on wh 
  - **Headings:** FitREQ CPU Tutorial (student guide),   - 1) Physical model (REQ + H),   - 2) Inputs and where you control things,   - 3) What `opt_2D.py` does (step-by-step),   - 4) How to modify for your system,   - 5) Chemical notes for H-bonds,   - 6) Key functions you will touch (Python),   - 7) Minimal example (Morse-SR)
    ... (+2 more)

- `FitREQ_ocl.md` (32 KB)
  - **Summary:** This program is a high-performance tool designed to optimize the parameters of a classical force field for describing non-covalent interacti 
  - **Headings:** Documentation of OpenCL-accelerated Non-Covalent Interaction Fitting Program (`NonBondFitting.py`, `FitREQ.cl`),   - Part 1: Conceptual Overview,     - 1.1 Introduction,     - 1.2 Core Physical and Mathematical Principles,   - Part 2: User Guide,     - 2.1 Overview of Input Files,     - 2.2 Input File Formats, n0 3 Etot 3.4667 ...
    ... (+10 more)

- `FitREQ_ocl_vispy.md` (31 KB)
  - **Summary:** OK, what we want to do now is some interactive fitting system using VisPy, we have @VispyUtils.py and pyQt with @BaseGUI.py example where it 
  - **Headings:** USER, SWE-1.6, Interactive pyOpenCL Fitting System - Plan,   - Overview,   - Components,     - 1. Reference Implementations to Study,     - 2. GUI Layout,     - 3. Parameter Formats
    ... (+39 more)

- `FitREQ_Paolo.md` (2 KB)
  - **Summary:** What was exactly the problem with Host atom of free electron pair ? 
  - **Headings:**     - Charge conservation $Q_{ep} + Q_{host} = Q_{resp}$,     - Variation of H1 depends on Q ?

- `GPU-CPU-Debugging-and-OCL-Interface.md` (13 KB)
  - **Summary:** - Scope: How to debug GPU (OpenCL) vs CPU (C++) implementations, and how to build robust OpenCL interfaces against `OCL.h`. - Audience: Any subsystem doing CPU/GPU parity checks, not limited to UFF.
  - **Headings:** General GPU vs CPU debugging and OCL.h interface,   - Implementation Log - 2025-09-13: CPU/GPU UFF parity achieved,     - Summary,     - Symptoms Observed,     - Root Causes,     - Fixes Implemented,     - Validation Results (after fixes),     - Key File Touchpoints
    ... (+3 more)

- `Limit_Language_Server_Memory.md` (18 KB)
  - **Summary:** Do you have any idea how to solve this problem with winsurf language server eating all my memory? I want some genral hard solution which bas **Summary**
  - **Headings:** USER, DeepSeek-V4,   - Step 1 – Create a cgroup with a memory limit, Create a cgroup named "limit_lsp", Set memory limit to 5 GB (use 5G, 5000M, or whatever you prefer), (Optional) also limit swap usage to prevent swapping out – same value is fine,   - Step 2 – Automatically put the language server into that cgroup,   - Step 3 – Run the script as a systemd user service (auto‑start & resilient)
    ... (+54 more)

- `LLCAO_1D.md` (11 KB)
  - **Summary:** I'm investigating possibilities to create approximative quantum mechanical method which rely on localization of molecular orbitals. That is  
  - **Headings:** Localized Linar Combination of Atomic Orbitals (LLCAO),   - User Prompt,   - Design Document: 1D LCAO Quantum Solver,     - 1. Introduction,     - 2. Theoretical Framework,     - 3. Basis Set: 1D Gaussian Functions,     - 4. Matrix Elements Calculation,     - 5. Step-by-Step Implementation Plan

- `LocalyOrthogonalizeOrbitalForcefield.md` (11 KB)
  - **Summary:** We want to create something between reactive forcefield, electron forcefield, and approximative quantum mechanical model of chemical bonding 
  - **Headings:** Localy Orthogonalized Orbital Forcefield,   - Motivation,   - Nature of chemical bonding,   - The basic idea,   - On-site interaction - exact iterative orthogonalization,   - Inter site interactions,   - Electrostatics,   - Electron-electron interaction
    ... (+1 more)

- `MD_MMFF_OCL_notes.md` (18 KB)
  - **Summary:** Goal: Run MD for molecules using MMFFsp3 force field via pyOpenCL (`pyBall/OCL/`). Implement and expose two new propagator kernels: `updateA 
  - **Headings:** Notes: MD with MMFFsp3 in OpenCL,     - Relevant files/functions,     - Important offsets and 2D-NDRange,     - PI-orbital placement,     - Back-neighbor forces,     - Plan / Next steps,     - Invariant monitoring (`tests/tUFF/test_MMFFsp3_pyOCL.py`),     - VdW subtraction toggle (`bSubtractVdW`)
    ... (+12 more)

- `MD_MMFF_OCL_parity_discussion.md` (89 KB)
  - **Summary:** look at @test_UFF_multi.py @test_MMFFsp3_pyOCL.py @test_UFF_ocl.py @MD_test_utils.py @UFF.py @UFFbuilder.py @MMFF.py @MolecularDynamics.py 
  - **Headings:**   - USER,   - GPT-5.1-Codex,   - USER,   - GPT-5.1-Codex,   - What’s missing,   - Minimal changes to implement,   - Files/functions to touch,   - Test commands (bonded only, correct CWD)
    ... (+136 more)

- `MMFFsp3_pi_angular_momentum_problem.md` (44 KB)
  - **Summary:** The problem is that pi-orbital dynamics in MMFFsp3 does not conserve angular momentum causing some angular momentum drift. 
  - **Headings:** MMFFsp3 pi-DOF angular momentum problem, GPT-5 (ChatGPT), Gemini 2.5 pro (Google AI studio),     - Executive Summary of the Issue,     - Detailed Physical and Code Analysis,     - The Root Cause of the Problem,     - Detailed Derivation and The Bug,     - The Discrepancy
    ... (+3 more)

- `MolecularEditor_2D_py.md` (13 KB)
  - **Summary:** - **Purpose** Build a PyQt5 desktop tool for 2D planar molecular sketching/editing while maintaining reusable 3D data via [AtomicSystem](cci - **Primary Users** Surface-chemistry researchers prototyping adsorbed molecules.
  - **Headings:** (2D) Molecule Editor in python,   - Scope & Goals,   - Functional Requirements,   - Non-Functional Requirements,   - System Architecture,   - Interaction Design,   - Data Flow,   - AtomicSystem Integration
    ... (+5 more)

- `OpenCLKernelMetaprogramming.md` (3 KB)
  - **Summary:** This document describes a general approach for dynamic OpenCL kernel generation using Python-based textual substitution. The system enables: - Runtime code substitution in kernel templates
  - **Headings:** OpenCL Kernel Metaprogramming Guide,   - 1. Introduction,   - 2. Core Workflow,   - 3. Key Components (General),     - 3.1 Host Manager Class,     - 3.2 Kernel Template,     - 3.3 Test Harness,   - 4. Implementation Patterns
    ... (+15 more)

- `rotational_dynamics_implementation.md` (5 KB)
  - **Summary:** This document describes the implementation of rotational dynamics for pi-orbitals in the pyOpenCL MD framework, designed to properly conserv 
  - **Headings:** Rotational Dynamics Implementation for Pi-Orbitals,   - Overview,   - Changes Made,     - 1. MolecularDynamics.py - Added Rotational Force Kernel Driver,     - 2. MD_test_utils.py - Updated to Use Optimized Method,     - 3. Test Script Updates,     - 4. Comparison Script,   - Usage Examples
    ... (+12 more)

- `TestContect.md` (7 KB)
  - **Summary:** path: fortran/MAIN/ktransform.f90 range: (whole file)
  - **Headings:**   - ktransform.f90,   - libFireCore.f90

- `ToDo_FastCollision_2.md` (21 KB)
  - **Summary:** [Position based dynamics](https://matthias-research.github.io/pages/publications/posBasedDyn.pdf) and [projective dynamics](https://www.proj 
  - **Headings:** Particle collisions by linear solvers,   - Problems:,   - Split Quasi Linear Poential,   - The Linear poential cutoff,   - The smoothening potential

- `ToDo_FastCollision_3.md` (14 KB)
  - **Summary:** This document outlines the detailed implementation plan for integrating fast particle collisions using projective dynamics. The key innovati 1. Linear (harmonic) part - handled by iterative linear solver (e.g. Jacobi, Gauss-Seidel)
  - **Headings:** Fast Particle Collisions by Projective dynamics,   - Background,   - Problems:,   - Split Quasi Linear Poential,   - The Linear poential cutoff,   - The smoothening potential, Implementation plan,   - Key Components and Files
    ... (+3 more)

- `ToDo_Fit_Hbond_GPU.md` (3 KB)
  - **Summary:** * Polarizability by fitting force-constants of springs attached to the electron pairs to host-atom * Try non-linear fit with relaxation of the e-pairs and sigma-holes
  - **Headings:** **H-Bond fitting**,   - Ideas,   - Implementation on GPU,   - Practical problems

- `ToDo_GUI.md` (5 KB)
  - **Summary:** Generally, it makes sense to draw inspiration from the GUIs of other 3D editors and molecular editors. In particular: 
  - **Headings:** GUI Functionality to be Added, Topology Modification, Camera, Selection, Pose Manipulation, Topology-Based Pose Manipulation, Measurement, Visualization
    ... (+1 more)

- `UFF-GPU-vs-CPU-Debug-Guide.md` (3 KB)
  - **Summary:** - End‑to‑end steps 
  - **Headings:** Practical workflow to validate GPU vs CPU UFF

- `UFF_MMFF_ocl_parity_test.md` (4 KB)
  - **Summary:** Provide a single, repeatable procedure to check **force parity** between CPU (C++) and GPU (OpenCL) for **bonded-only** terms of UFF and MMF 
  - **Headings:** UFF/MMFF OpenCL Parity Test (Bonded-only, CPU vs GPU),   - Goal,   - Scope and assumptions,   - Quick commands (UFF, bonded-only), UFF: CPU vs GPU forces, bonded-only,   - Quick commands (PyOpenCL UFF kernel isolation),   - MMFF (bonded-only) using the same C++ interface, CPU
    ... (+6 more)

- `UFF_muli_plan.md` (18 KB)
  - **Summary:** I'll start by inspecting the referenced C++ and Python files to see how MMFFsp3 is implemented and what needs to be added for UFF. I'll open 
  - **Headings:** Findings, Recommended Actions, Validation Plan, Next Steps, Status, Findings, What I changed, Why this matches the architecture
    ... (+2 more)

- `UFF_ocl_plan.md` (29 KB)
  - **Summary:** I want to implement UFF on GPU using OpenCL and test it. Making simple test script which check forces on atom calculated using CPU version o 
  - **Headings:** GUI app integration, python / library integration,   - Phase 1: Analysis and Planning,   - Phase 2: Test Script Implementation, UFF OpenCL Implementation Plan,   - Current State Analysis,     - Existing UFF Components,     - Reference MMFF Implementation (for comparison)
    ... (+36 more)

- `UFF_ocl_plan_new copy.md` (14 KB)
  - **Summary:** I want to implement UFF on GPU using OpenCL and test it. Making simple test script which check forces on atom calculated using CPU version o 
  - **Headings:** GUI app integration, python / library integration,   - Phase 1: Analysis and Planning,   - Phase 2: Test Script Implementation,   - UFF OpenCL Implementation Plan,     - 1. High-Level Goal,     - 2. System Architecture,     - 3. Component Analysis
    ... (+1 more)

- `UFF_ocl_plan_new.md` (8 KB)
  - **Summary:** The primary objective is to create a GPU-accelerated version of the UFF force field using OpenCL, fully integrated into the `FireCore` multi 
  - **Headings:**   - UFF OpenCL Implementation Plan,     - 1. High-Level Goal,     - 2. System Architecture,     - 3. Component Analysis,     - 4. Detailed Implementation Plan

- `UFF_ocl_plan_updated copy.md` (6 KB)
  - **Summary:** **Priority Order:** 1. **`evalBondsAndHNeigh_UFF`** - Foundation kernel that:
  - **Headings:** UFF OpenCL Implementation - Updated Technical Plan,   - Critical Implementation Insights,     - 1. Kernel Development Strategy,     - 2. Buffer Management Strategy,     - 3. Integration with MMFFBuilder Workflow,   - Implementation Roadmap Refinement,     - Phase 1: Core Kernel Development (High Priority),     - Phase 2: C++ Integration Layer
    ... (+11 more)

- `UFF_ocl_plan_updated.md` (8 KB)
  - **Summary:** The project will proceed along two parallel but complementary paths. This strategy allows for rapid, isolated testing of the low-level OpenC 
  - **Headings:** UFF OpenCL Implementation - Updated Technical Plan,   - Refined Development Strategy: A Two-Path Approach,     - Path A: Direct PyOpenCL Kernel Validation,     - Path B: High-Performance C++ Integration,   - Critical Implementation Insights,     - 1. Kernel Development Strategy,     - 2. Buffer Management Strategy,     - 3. Integration with MMFFBuilder Workflow
    ... (+14 more)


## `doc/FMM/`

- `FMM.md` (9 KB)
  - **Summary:** https://aistudio.google.com/app/prompts?state=%7B%22ids%22:%5B%221jo9BKPgkz9Ro_gR_T5_naDegd5wm2n-e%22%5D,%22action%22:%22open%22,%22userId%2 
  - **Headings:**   - USER,   - Gemini,     - 1. The Strategy: Energy Blending,     - 2. The Multipole Math (Charge + Dipole),     - 3. Deriving the Forces (The "Messy" Part),     - 4. Explicit Formulas for Implementation,     - 5. Quadrupole: Yes or No?,     - 6. Summary for GPU Kernel Design


## `doc/Markdown/`

- `eFF_documentation.md` (19 KB)
  - **Summary:** 1.  The core eFF logic is implemented in C++ (`eFF_lib.cpp`, `eFF.h`, `InteractionsGauss.h`). 2.  The Python file (`eFF.py`) provides an interface to this C++ code, allowing it to be used from Python.
  - **Headings:** Electron Forcefield (EFF) Documentation,   - 1. Files relevant for the  Electron forcefield (EFF) project,     - Description of the Files:,   - Functions implemented in each file:,     - `tests/tEFF/run_tests.py`,     - `pyBall/eFF.py`,     - `cpp/libs/Molecular/eFF_lib.cpp`,     - `cpp/common/molecular/eFF.h`
    ... (+1 more)

- `graphen_ribbon_lattice_scan.md` (45 KB)
  - **Summary:** This document describes the refactoring of the lattice scan functionality for graphene ribbons in FireCore. The refactoring makes the lattic 
  - **Headings:** Graphene Ribbon Lattice Scan - Refactoring Documentation,   - Overview,   - Motivation,     - Original Problems,     - Goals,   - Changes Made,     - 1. `pyBall/FireCore.py`,     - 2. `pyBall/plotUtils.py`
    ... (+71 more)

- `GUI.md` (14 KB)
  - **Summary:** This document provides a comprehensive guide to the GUI system implemented in `GUI.h` and `GUI.cpp`, with usage examples from `MolGUI.h`. It 
  - **Headings:** GUI System Documentation,   - Design Philosophy,   - Key Classes,     - `GUI`,     - `GUIAbstractPanel`,     - `GUIPanel`,     - `MultiPanel`,     - `CheckBoxList`
    ... (+15 more)

- `kekule_research_notes copy 2.md` (33 KB)
  - **Summary:** We want to make nice visual app GUI for easy and fast exploration of Kekule structures of relaed aromatic system an that either finite or in 
  - **Headings:** USER,     - DFTB interface:,     - plotting and GUI,     - Ribbonds generation and optimization, Gemini-3-flash, Kekule Explorer: Code Research & Architecture Plan,   - 1. Component Research,     - 1.1 DFTB+ Interfacing (`pyBall/dftb_utils.py`)
    ... (+46 more)

- `kekule_research_notes copy.md` (15 KB)
  - **Summary:** We want to make nice visual app GUI for easy and fast exploration of Kekule structures of relaed aromatic system an that either finite or in 
  - **Headings:** USER,     - DFTB interface:,     - plotting and GUI,     - Ribbonds generation and optimization, Gemini-3-flash, Kekule Explorer: Code Research & Architecture Plan,   - 1. Component Research,     - 1.1 DFTB+ Interfacing (`pyBall/dftb_utils.py`)
    ... (+30 more)

- `kekule_research_notes.md` (129 KB)
  - **Summary:** We want to make nice visual app GUI for easy and fast exploration of Kekule structures of relaed aromatic system an that either finite or in 
  - **Headings:** USER,     - DFTB interface:,     - plotting and GUI,     - Ribbonds generation and optimization, Gemini-3-flash, Kekule Explorer: Code Research & Architecture Plan,   - 1. Component Research,     - 1.1 DFTB+ Interfacing (`pyBall/dftb_utils.py`)
    ... (+216 more)

- `KekuleGUI_Extension_System.md` (49 KB)
  - **Summary:** plese figure out how to make some extensiton manager in @KekuleExplorerGUI.py currently it use two external program/libraries DFTB and Fireball, since these progam mey be not installed we must load them lazzyliy, set s
  - **Headings:** USER, SWE-1.6,   - Current Extension Handling Analysis,     - KekuleExplorerGUI.py Current State,     - Extension File Analysis,   - Specific Requirements Identified,     - 1. **Lazy Loading Strategy**,     - 2. **Path Configuration**
    ... (+58 more)

- `MMFFBuilder.md` (21 KB)
  - **Summary:** This document provides a comprehensive overview of the `MMFFBuilder.h` header file, focusing on the `MM::Builder` class and related structur 
  - **Headings:** `MMFFBuilder.h` Documentation,   - 1. Introduction,   - Helper Classes,   - class `Atom`,   - class `AtomConf`,     - class `Bond`,     - class `Angle`,     - class `Dihedral`
    ... (+8 more)


## `doc/Markdown/cpp/common/molecular/`

- `Atoms.h.md` (4 KB)
  - **Summary:** The `Atoms.h` header file defines the `Atoms` class, which is used to store and manipulate atomic data in molecular dynamics simulations. It 
  - **Headings:** Atoms.h,   - Includes,   - Types (classes and structs),     - Class `MetaData`,     - Class `Points`,   - Class `Atoms`,   - Notes

- `EwaldGrid.h.md` (4 KB)
  - **Summary:** `EwaldGrid.h` implements the Ewald summation method for calculating electrostatic potentials in molecular dynamics simulations. It projects  
  - **Headings:** EwaldGrid.h,   - Includes,   - Free functions,     - Index Manipulation,     - FFTW Array Conversion,   - Types (classes and structs),     - class `EwaldGrid`

- `ForceField.h.md` (6 KB)
  - **Summary:** This header file defines the `ForceField` class, which is responsible for managing and calculating forces in a molecular dynamics simulation 
  - **Headings:** ForceField.h,   - Includes,   - Free Functions,   - Types (Classes and Structs),     - class `ForceField`,     - class `CollisionDamping`

- `GridFF.h.md` (7 KB)
  - **Summary:** `GridFF.h` defines the `GridFF` class, which provides functionality for evaluating force fields on a grid structure. This class is used to c 
  - **Headings:** GridFF.h,   - Includes,   - Free Functions,   - Types (Classes and Structs),     - Class `GridFF`

- `MMFFsp3_loc.h.md` (5 KB)
  - **Summary:** This file defines the `MMFFsp3_loc` class, which implements a localized version of the MMFF force field for classical molecular mechanics. T 
  - **Headings:** MMFFsp3_loc.h,   - Includes,   - Types (Classes and Structs),     - Class `MMFFsp3_loc`

- `MolWorld_sp3.h.md` (20 KB)
  - **Summary:** `MolWorld_sp3.h` is a header file that defines the `MolWorld_sp3` class, which serves as a comprehensive container for managing molecular si 
  - **Headings:** MolWorld_sp3.h,   - File Purpose,   - Includes,   - Types (Classes and Structs),     - Class `MolWorld_sp3`

- `NBFF.h.md` (10 KB)
  - **Summary:** The `NBFF.h` file defines a class that implements non-bonded force-field interactions between particles using potentials such as Lenard-Jone 
  - **Headings:** NBFF.h,   - Includes,   - Free Functions,   - Types (Classes and Structs),     - Class `NBFF`

- `UFF.h.md` (5 KB)
  - **Summary:** The `UFF.h` file implements a universal force field (UFF) for molecular systems based on the article "J. Am. Chem. Soc. 1992, 114, 25, 10024 
  - **Headings:** UFF.h,   - Includes,   - Free Functions,   - Types (Classes and Structs),     - Class `UFF`


## `doc/Markdown/cpp/common_SDL/SDL2OGL/`

- `EditorGizmo.h.md` (4 KB)
  - **Summary:** This header defines the `EditorGizmo` class, a Blender-like 3D transform gizmo for interactive translation, scaling, and rotation of selecte 
  - **Headings:** EditorGizmo.h,   - Includes,   - Free Functions,   - Types,   - Class `EditorGizmo`,     - Important members,     - Binding,     - Transform application
    ... (+7 more)

- `FragmentManipulation_Design.md` (15 KB)
  - **Summary:** Status: draft notes (2025-08-10) 
  - **Headings:** Molecular Editor: Fragment Selection & Manipulation (Design Scratchpad),   - Goals,   - Relevant APIs and Data,   - Reuse Plan,   - Missing Pieces / Proposed Additions (C++; not implemented yet),   - UI/UX Flow,   - Data Ownership / Synchronization,   - Edge Cases / Constraints
    ... (+7 more)

- `MolecularDraw.h.md` (2 KB)
  - **Summary:** Utility drawing helpers for molecular visualization and field-based rendering using OpenGL. Provides color maps, neighbor/trajectory visuali 
  - **Headings:** MolecularDraw.h,   - Includes,   - Color helpers,   - Debug/printing,   - Neighborhood and mapping,   - Grid-force interaction and trajectories,   - Planes and substrates,   - Notes
    ... (+1 more)

- `MolGUI.h.md` (18 KB)
  - **Summary:** This header defines the `MolGUI` class, a 3D molecular editor and visualization app built on SDL2/OpenGL. It manages molecular data binding, 
  - **Headings:** MolGUI.h,   - Includes,   - Free Functions (helpers),   - Class `MolGUI`,     - Construction,     - Core virtuals (override),     - Event modes and mouse,     - Binding / data setup
    ... (+28 more)


## `doc/Maxima/`

- `GaussianIntegrals.md` (23 KB)
  - **Summary:** Define Gaussian Functions shifted along x-axis and Laplacian of Gaussian 
  - **Headings:** Gaussian Integrals,   - Kinetic Energy,   - eFF Pauli potential in Gaussain Basis, Spherical Harmonics / Multiploles, Matching Lorenz like function,     - Bspline General Recursive formulas

- `MaximaTutorial.md` (13 KB)
  - **Summary:** This tutorial introduces Maxima, a powerful open-source computer algebra system, focusing on commands and techniques relevant to scientific  
  - **Headings:** Maxima Tutorial: Symbolic Computation for Scientists,   - 1. Getting Started with Maxima,     - 1.1. Basic Interaction,     - 1.2. Assignments,     - 1.3. Assumptions,   - 2. Defining Functions,   - 3. Symbolic Calculus,     - 3.1. Differentiation
    ... (+11 more)


## `doc/Topics/`

- `MoleculeBuilder.md` (8 KB)
  - **Summary:** I want to make simple graphical user interface for placing molecules on ionic substrate. My usin substrate fro starter is NaCl crystal, but  
  - **Headings:**   - USER,   - GPT-5.1-Codex,     - Differences between NaCl step scripts,     - Useful pieces in AtomicSystem,     - Useful pieces in MolecularPlacer (template for GUI features),     - Proposed goals (systematic),     - High-level plan


## `doc/Topics/AFM/`

- `AFM.md` (21 KB)
  - **Summary:** This document provides a comprehensive overview of the Atomic Force Microscopy (AFM) simulation package, detailing its structure, workflow,  
  - **Headings:** AFM Simulation Package Documentation,   - 1. Code Map: Key Files and Modules,   - 2. Workflow Review,     - 0. Mathematical Background for Density-Based Models,     - Step 1: DFT Calculation (Fireball),     - Step 2: Density Projection to Grid,     - Step 3: Potential & Energy Calculation,     - Step 4: Probe Particle Relaxation
    ... (+25 more)

- `AFM_FDBM_DFTB.md` (60 KB)
  - **Summary:** perfect! we did it. Now we should investigate how we can use this density for comouting afm images using Fill-densty based model (FDBM), 
  - **Headings:** USER, SWE-1.6, FDBM AFM with DFTB - Implementation Plan,   - Current Implementation Overview,   - What Needs to Change,     - 1. **Import Changes**,     - 2. **Step 1: Density Projection** (MAJOR CHANGE), Run Fireball SCF
    ... (+131 more)

- `AFM_FDBM_fitting.md` (56 KB)
  - **Summary:** CODEMAP: DFTB Reference Calculation & FDBM AFM Forcefield Comparison System https://windsurf.com/codemaps/1153fe89-ff29-4d4b-b4a6-e97d8f37047f-fe86ab10a43f3d18
  - **Headings:** USER, SWE-1.6, Comprehensive Plan for DFTB Reference Calculation and Grid Point Identification,   - Task 1: DFTB-based Z-Scan for CO Tip Approaching Pentacene,     - Objective,     - Current Capabilities Analysis,     - Recommended Approach: Use DFTBcore with Coordinate Updates,     - Alternative Approach: dftb_utils.constrained_scan()
    ... (+154 more)

- `AFM_FDBM_optimization.md` (43 KB)
  - **Summary:** I do not think that is opencl conplialtion which is slow. There must be something else you are doing in python which is slow, maybe you agai 
  - **Headings:** USER, SWE-1.6, USER, SWE-1.6, USER, SWE-1.6,   - Current Inefficiencies:,   - Optimization Opportunities:
    ... (+72 more)

- `AFM_FDBM_profiling_optimization.md` (150 KB)
  - **Summary:** OK, we have this project for AFM simulation @tests/tAFM/pyocl_fdbm/test_full_pipeline.py @pyBall/OCL/AFM_utils.py @pyBall/OCL/AFM.py @pyBall/DFTB/DFTBcore.py @pyBall/DFTB/Grid_dftb.py 
  - **Headings:** USER, SWE-1.6, AFM Pipeline Performance Analysis,   - 1. Python-Level Overheads Identified,     - **Major Overheads:**,     - **Minor Overheads:**,   - 2. Profiling Approach (No Code Modification),     - **Recommended Profiling Tools:**
    ... (+284 more)

- `AFM_FDBM_pySCF.md` (56 KB)
  - **Summary:** CODEMAP: AFM FDBM Pipeline: DFTB Backend & pySCF Integration Points https://windsurf.com/codemaps/02d559c9-de47-4058-b07b-3318664b454e-fe86ab10a43f3d18
  - **Headings:** USER, Gemini-2.5-flash,     - Why PySCF fits your requirements:,     - Python Code Implementation, --- Example Usage (Anthracene-like minimal test geometry) ---,     - Performance expectations for your workflow:, Kimi 2.6,     - 1. Density on a precisely specified grid in pySCF
    ... (+95 more)

- `AFM_migration_discusion.md` (83 KB)
  - **Summary:** read in @AFM.md and @AFM_migration_plan.md , sumarize what is it about and mention all relevant files in thich chat. 
  - **Headings:** USER, GPT-5.1-Codex, USER, GPT-5.1-Codex, USER, GPT-5.1-Codex, USER, GPT-5.1-Codex
    ... (+104 more)

- `AFM_migration_plan.md` (38 KB)
  - **Summary:** AFM Density-Based Simulation: Classical Force Fields and Density Projection Workflows https://windsurf.com/codemaps/7a2c010f-2075-4a19-893c-1b54ec8df995-fe86ab10a43f3d18
  - **Headings:** AFM PyOpenCL Migration Plan (2025 Update),   - 1. Scope and Goals,   - 2. Summary of the Current Workflow,   - 3. Target Python Package Layout,   - 4. Detailed Migration Tasks,     - 4.1 Stage 0 – Preparatory Work,     - 4.2 Stage 1 – Context & Resource Management (C++: `OCL_GridFF.cpp`, `OCL_DFT` core),     - 4.3 Stage 2 – FFT Infrastructure (C++: `initFFT`, `runFFT`, `convolution` scaffolding)
    ... (+20 more)

- `AFM_migration_progress.md` (12 KB)
  - **Summary:** Below is a concise report (pre-FDBM), based on [doc/Topics/AFM/AFM_migration_discusion.md](cci:7://file:///home/prokophapala/git/FireCore/do 
  - **Headings:**   - Goal,   - What we implemented (pre-FDBM),   - Key problems and fixes,   - Remaining work,   - Insights / takeaways for future debugging,   - Phase 3: Full Density-Based Model (FDBM) AFM (Fireball DFT),     - What we achieved,     - Problems we hit
    ... (+13 more)


## `doc/Topics/FTIR_Nanocrystals/`

- `gen_nanocrystals.md` (51 KB)
  - **Summary:** description: Silicon nanocrystal generator (scripts/gen_nanocrystals.mjs) 
  - **Headings:** Silicon Nanocrystal Generator,   - What it does,   - CLI (key options),   - Current surface-bridge logic (Si),   - Usage examples, Generate a set of nanocrystals,   - Linear Response Vibration Spectroscopy (FTIR),     - Key Components
    ... (+99 more)

- `Hessian_fitting.md` (54 KB)
  - **Summary:** Consider following physical problem - fitting parameters of classical forcefield to reproduce vibration Spectra resp Hessian. For simplicity 
  - **Headings:** USER, DeepSeek,     - ⚛️ From Hessians to Parameters,     - 🔧 Adding Regularization,     - 🚀 Off-the-Shelf Solutions,   - Which Python-friendly package exists?,   - Programming from scratch: a robust blueprint,     - 1. Energy and Hessian in internal coordinates
    ... (+88 more)


## `doc/Topics/ForceFields/`

- `GridFF_fit_CG.md` (3 KB)
  - **Summary:** - Implemented faster fitting for 3D B-spline grid potentials using Conjugate Gradient (CG) alongside existing damped MD/Jacobi. - Applies to Morse-only fitting (no Coulomb/Ewald) in the OpenCL GridFF path.
  - **Headings:** GridFF B-spline Fitting: CG vs Damped MD,   - Scope,   - Key Implementations,   - Design Decisions,   - BsplineConv3D coefficient conventions,   - Performance (NaCl_1x1_L3, Morse-only),   - How to Run, CG Morse-only fit with saved plot
    ... (+3 more)

- `GridFF_fit_jacobi.md` (44 KB)
  - **Summary:** I have cubic B-spline fitting in 3D on GPU in opencl / pyOpenCL which is currently using gradient descent (resp. dynamical relaxation with i 
  - **Headings:**   - User 1 :,   - LLM 1 :,     - Summary of the Jacobi Method for Your Problem,     - 1. New OpenCL Kernel for the Jacobi Update,     - 2. Updated Python Driver with Jacobi Solver, You would have your existing class structure here., class GridFF_cl:, def __init__(self, ...):
    ... (+29 more)

- `NeighborExclusion_Summary.md` (7 KB)
  - **Summary:** Computing non-covalent interactions on the GPU requires excluding 1-2 (bonded) and 1-3 (angle) pairs while preserving high throughput. Each  
  - **Headings:** Neighbor Exclusions for Non-Covalent OpenCL Kernels,   - Context and Goals,   - Assumptions and Preconditions,   - Evaluated Alternatives,   - Final Data Layout and Algorithm,   - Host-Side Preparation (C++),   - Kernel Fragment (OpenCL C),   - Key Takeaways

- `UFF_multiSystem.md` (61 KB)
  - **Summary:** These notes capture the target design for running UFF over `nSystems` replicas in parallel on the GPU, and the action items to fix the curre 
  - **Headings:** UFF multi‑system (GPU/OpenCL) design notes,   - Problem recap,   - Relevant Files,   - Target execution model,   - NDRange and OCL task configuration,   - Host/Device buffer layout and uploads,   - Kernel specifics,   - Debug instrumentation
    ... (+99 more)

- `UFF_multiSystem_summary.md` (8 KB)
  - **Summary:** These notes capture the target design for running UFF over `nSystems` replicas in parallel on the GPU, and the action items to fix the curre 
  - **Headings:** UFF multi‑system (GPU/OpenCL) design notes,   - Problem recap,   - Relevant Files, UFF OpenCL Multi-System Implementation: Final Summary and Future Reference,   - 1. Core Architecture: The Multi-System Execution Model,     - A. Kernel Launch Configuration (NDRange),     - B. Universal Kernel Indexing,     - C. Per-System Base Offsets (The Golden Rule)
    ... (+12 more)


## `doc/Topics/Kekule_Topology/`

- `Kekule_Optimizer.md` (21 KB)
  - **Summary:** Now lets combine chemistry and this pierls distortion effective hamiltionina picture. Notmally the kekule texture (configuration of double-s 
  - **Headings:** USER,   - Gemini-3-flash,   - 1. The \"Golden Rule\": Bond Dissociation Free Energy (BDFE),   - 2. Using the Thermodynamic Square Scheme,     - The Balancing Act:,   - 3. Design Principles for Energy Storage (Flow Cells),   - 4. Summary Table for Quick Estimation,     - Practical Rule:
    ... (+21 more)

- `Kekule_Topology.md` (38 KB)
  - **Summary:** can you read this article? https://pubs.acs.org/doi/10.1021/acs.nanolett.3c04703
  - **Headings:** USER, Gemini-3-flash,     - 1. The Effective Hamiltonian Method,     - 2. How to determine conductivity for your \"Hand-Drawn\" structure,     - 3. Practical Computation Steps, Deepseek-V4,     - 🔍 How the Article's Methodology Works,     - ⚛️ Applying the Method to Quinonic and Aromatic Domains
    ... (+92 more)


## `doc/Topics/LinearLagebra/`

- `Chebyshev_Acceleration.md` (5 KB)
  - **Summary:** This document provides a complete implementation guide for the Chebyshev-accelerated Jacobi solver as described in [Wang 2015] (Section 2 of 
  - **Headings:** Chebyshev-Accelerated Jacobi Solver for Projective Dynamics,   - 1. Introduction,   - 2. Core Algorithm,   - 2.1 Chebyshev Acceleration Formula,     - 2.2 Chebyshev Weights Calculation,     - 2.3 Complete Implementation,   - 3. Important Implementation Notes,     - 3.1 Spectral Radius (ρ)
    ... (+9 more)

- `ChebyshevAndInterialAccelerationOfLinearSolvers.md` (26 KB)
  - **Summary:** The core idea is to accelerate the convergence of an existing iterative solver (like Projective Dynamics or Position-Based Dynamics). Standa 
  - **Headings:**   - Design Document: Chebyshev Acceleration for Iterative Solvers,     - 1. High-Level Concept and Goal,     - 2. Theoretical Background (The Linear Case),     - 3. The Core Update Formula,     - 4. Adapting to Nonlinear Systems (The Paper's Key Insight),     - 5. The Critical Parameter: Estimating the Spectral Radius ($\rho$),     - 6. Practical Implementation (Pseudocode and Caveats),     - 7. Summary of Key Parameters and Ranges
    ... (+9 more)

- `NumricalStabilityLinearSolvers.md` (31 KB)
  - **Summary:** When simulating truss or spring-mass systems with implicit integration (e.g. implicit Euler), we must solve a large linear system at every t 
  - **Headings:** Iterative Solvers for Truss Dynamics: Jacobi and Gauss–Seidel on displacements,   - 1. Motivation,     - 1.1 Update loop in C++,     - 1.2 Iterative Linear solver,   - 2. Working with Displacements,   - 3. Jacobi vs. Gauss–Seidel,     - Jacobi,     - Gauss–Seidel
    ... (+24 more)

- `Parallel_GaussSeidel.md` (16 KB)
  - **Summary:** Nice — you’re on the right track and your high-level kernel loop is exactly what people use for GPU Gauss–Seidel style solvers. Below I: 
  - **Headings:**   - Short critique / checklist (things to think about),   - Assumptions for the kernel below,   - Kernel: OpenCL implementation (single workgroup handles the whole mesh / chunk),   - Notes, tuning knobs and follow-ups


## `doc/Topics/ManyBody/`

- `HubbardSolver.md` (28 KB)
  - **Summary:** This document provides an overview of the OpenCL kernels defined in `cl/hubbard.cl` and their interaction with the buffer management in `pyB 
  - **Headings:** Hubbard Solver OpenCL Kernels and Buffer Management,   - OpenCL Kernels Overview,     - 1. `eval_coupling_matrix`,     - 2. `eval_Oriented_Hopping`,     - 3. `solve_minBrute_fly`,     - 4. `solve_minBrute_boltzmann`,     - 5. `solve_local_updates`,     - 6. `solve_MC_neigh`
    ... (+7 more)

- `HubbardSolver_new.md` (24 KB)
  - **Summary:** This document provides a comprehensive overview of the Hubbard Solver project, focusing on its OpenCL kernels for GPU acceleration and their 
  - **Headings:** Hubbard Solver: Comprehensive Documentation,   - 1. Physical Principles and Computational Challenges,     - 1.1. The Hubbard Model in STM Simulations,     - 1.2. Limitations of Greedy Descent,     - 1.3. Configuration Representation,   - 2. The Hybrid Monte Carlo Framework (Iterated Local Search - ILS),     - 2.1. Motivation and Core Idea,     - 2.2. Global Mutation Strategies
    ... (+21 more)

- `MC_PME_dicussion.md` (80 KB)
  - **Summary:** I have pyOpenCL Monte-Carlo solver of Hubard/Ising model searchin minimum energy occupancy of array of quantum dot depening of position of S The problem is that this is just ground-state solver. But In experiments we see that shard ptansitions beweeen states are non-equalibirum an
  - **Headings:** USER, Gemini-3-pro,     - 1. Physics & Topology of the Sub-Space,     - 2. The Numerical Solver Strategy (GPU),     - 3. Proposed Pipeline,     - 4. Mathematical Refinement for the Basis,     - 5. Summary of the Plan, USER
    ... (+67 more)

- `MC_PME_progress.md` (10 KB)
  - **Summary:** The system consists of **N quantum dots** (sites) that can be occupied by electrons. The state of the system is defined by an occupancy vect ```
  - **Headings:** MC + Pauli Master Equation Solver: Implementation Progress Report,   - 1. General Formulation of the Problem,     - 1.1 Physical System,     - 1.2 Hamiltonian,     - 1.3 Tunneling & Rates (Fermi's Golden Rule),     - 1.4 Pauli Master Equation (PME),     - 1.5 The Challenge: Exponential Scaling,     - 1.6 Subspace Selection Strategy
    ... (+24 more)

- `MQCA_GrayCodeSolver.md` (39 KB)
  - **Summary:** What we want to do is to solve molecualr cellular automata MQCA This is basically lattice of sites with some onsite energy and neighbor capacitve coupling described by Ising-like hamiltonina
  - **Headings:** MQCA Ground-State Solver: Gray Code Brute-Force Search, USER, Plan (Kimi K2.5), MQCA Ground-State Solver: Gray Code Brute-Force Search,   - Overview,     - Key Design Goals,   - Physical Model,     - Hamiltonian
    ... (+63 more)


## `doc/Topics/MolecularGUI/`

- `MoleculeEditor2D.md` (8 KB)
  - **Summary:** title: MoleculeEditor2D User Guide 
  - **Headings:** Overview, Quick Start, Core Components, Usage Tips, Integration Points, Future Work

- `MolViewer_PyQt.md` (17 KB)
  - **Summary:** I want to make simple viewer of molecular geometry using python, pyQt5 and moderngl library 
  - **Headings:**     - 1. Project Structure,     - 2. OpenGL Utilities (`gl_utils.py`), gl_utils.py,     - 3. Shader Files (`shaders/`),     - 4. Detailed 3D Molecule Viewer (`gl_molecule_viewer.py`), gl_molecule_viewer.py,     - 5. Thumbnail Grid Widget (`thumbnail_grid.py`), thumbnail_grid.py
    ... (+11 more)

- `MolViewer_PyQt_2.md` (23 KB)
  - **Summary:** Excellent. You've provided a well-structured and powerful set of base classes. Reusing `GLGUI.py`, `AtomicSystem.py`, and the concepts from  
  - **Headings:**     - Project Goal Recap,     - New File Structure,     - Step 1: Shaders for Ball-and-Stick Rendering,     - Step 2: The Detailed 3D Viewer (`detailed_viewer.py`), detailed_viewer.py,     - Step 3: The Main Browser Window (`molecular_browser_app.py`), molecular_browser_app.py,     - Step 4: The Main Entry Point (`main_browser.py`)
    ... (+2 more)

- `MolViewer_PyQt_toGemini.md` (97 KB)
  - **Summary:** This is good start but I ralized I already have several classes which I would like to reuse (use a base for our new viewer) so we do not rei 
  - **Headings:** OpenGL Viewer Documentation,   - 1. Core OpenGL GUI Framework (`pyBall/GUI/GLGUI.py`),     - 1.1. `GLGUI` Class,     - 1.2. `Mesh` Class,     - 1.3. `InstancedData` Class,     - 1.4. Utility Functions,     - 1.5. `AppWindow` Class,   - 2. Molecular Viewer Application (`tests/tEFF/xyz_view_new.py`)
    ... (+36 more)

- `OpenGL_Viewer_Documentation.md` (14 KB)
  - **Summary:** This document provides an overview of the OpenGL viewer components within the `FireCore` project, focusing on `pyBall/GUI/GLGUI.py`, `tests/ 
  - **Headings:** OpenGL Viewer Documentation,   - 1. Core OpenGL GUI Framework (`pyBall/GUI/GLGUI.py`),     - 1.1. `GLGUI` Class,     - 1.2. `Mesh` Class,     - 1.3. `InstancedData` Class,     - 1.4. Utility Functions,     - 1.5. `AppWindow` Class,   - 2. Molecular Viewer Application (`tests/tEFF/xyz_view_new.py`)
    ... (+9 more)


## `doc/Topics/OnSurfaceAssembly/`

- `AFM_File_Analysis_and_Architecture.md` (22 KB)
  - **Summary:** Analyze existing FireCore files related to visualization, relaxation scanning, and rigid body dynamics to design a modular architecture for  
  - **Headings:** AFM Simulation: File Analysis and Architecture Design,   - **📋 OBJECTIVE**,   - **🔍 FILE CATEGORIZATION AND MAPPING**,     - **Category 1: GPU Rigid Body Dynamics (Core Backend)**,     - **Category 2: CPU GridFF Relaxation (Reference Backend)**,     - **Category 3: Visualization and GUI Frontends**,     - **Category 4: Testing and Methodology**,   - **🏗️ PROPOSED MODULAR ARCHITECTURE**
    ... (+34 more)

- `Comprehensive_Overview_Molecule_Surface_Systems.md` (25 KB)
  - **Summary:** **Purpose**: This document consolidates all implemented modules, scripts, and visualization programs for generating rigid-body configuration 
  - **Headings:** Comprehensive Overview: Molecule-on-Surface Systems in FireCore,   - 1. High-Level Architecture Overview,     - 1.1 Core Subsystems,     - 1.2 Key Directories,   - 2. Subsystem Deep Dives,     - 2.1 GridFF (Precomputed Substrate Potentials),     - 2.2 InteractionScanner (XYZ Rigid Scanning),     - 2.3 FoldedAtomicFunctions (Basis Fitting for Surface Potentials)
    ... (+41 more)

- `ElectrostaticContinuumEmbeding_report.md` (25 KB)
  - **Summary:** Description: Electrostatic continuum embedding report 
  - **Headings:** Electrostatic Continuum Embedding for Surface Interaction Scans,   - Problem Statement,   - Analytic Model,     - Rectangle Bounds (critical fix),     - Layer Moments,   - Python Reference Implementation,   - OpenCL Implementation,     - Host (packing / launch)
    ... (+32 more)

- `Ewald_2D.md` (84 KB)
  - **Summary:** I have surface of ionic crystal 2D, and want to compute electrostatic potential. I want to represent it by small basis of sin and cosine in  
  - **Headings:** USER, DeepSeek-V4,     - 1. Vacuum potential above the slab,     - 2. Potential inside the slab,     - 3. Relation to the 2D Ewald sum,     - 4. Practical algorithm to get the basis coefficients,     - 5. Efficiency and implementation tips, Gemini-3-flash
    ... (+106 more)

- `ExplorerVisPy_GridFF_Interactive_Plan.md` (28 KB)
  - **Summary:** This document now includes the **completed implementation** of GridFF relaxed scanning with proper physics, caching, and visualization. The  
  - **Headings:** Interactive GridFF Scanning Implementation Plan,   - **Implementation Summary - COMPLETED**,   - **✅ COMPLETED IMPLEMENTATION**,     - **Core GridFF Infrastructure**,   - **🐛 PROBLEMS ENCOUNTERED AND SOLUTIONS**,     - **Problem 1: GridFF Generation Crashes**,     - **Problem 2: GridFF Cache Not Saved**,     - **Problem 3: Incorrect PLQ Scaling**
    ... (+68 more)

- `FoldedSubstratePotential_OpenCL.md` (16 KB)
  - **Summary:** OK, so now we would like to implement The Folded atomif function into the forcefield, study the test which we have about it overview is here @RigidSurfPotential_GridFF.md @FoldedAtomicFunction.md @FoldedAtomicFunction_optimization.md but most importaint are kernels
  - **Headings:** USER, SWE-1.6 : Folded Atomic Function Implementation Plan,   - Overview,   - Existing Code Review,     - 1. **Existing Kernels** (cpp/common_resources/cl/Surface.cl),     - 2. **Python Framework** (doc/py/FoldedAtomicFunctions/FoldedAtomicFunction.md),     - 3. **OpenCL Ewald2D** (pyBall/OCL/SurfaceEwald.py, cpp/common_resources/cl/Surface.cl),   - Implementation Plan
    ... (+25 more)

- `GridFF_atoms_alignment.md` (30 KB)
  - **Summary:** I have persistin problem in my implementation of GridFF (grid projected forcefield) https://pubs.acs.org/doi/10.1021/acs.jctc.5c01223
  - **Headings:** USER, Gemini-3.5-flash,     - Part 1: Retroactive Alignment via Grid Analysis, Load your grid (assuming shape Nx, Ny, Nz), 1. Project 3D -> 1D by taking the minimum along axes, (Good for Lennard-Jones wells; use np.max for repulsive cores or electrostatic peaks), 2. Find indices of peaks/valleys in grid index-space, 'distance' ensures we don't pick adjacent pixels for the same atom
    ... (+62 more)

- `GridFF_FDBM_Fitting.md` (131 KB)
  - **Summary:** CODEMAP: FDBM Fitting Infrastructure: GridFF Sampling, Ewald Electrostatics, and Folded Basis Fitting https://windsurf.com/codemaps/a8b020b7-8af8-46db-bb73-e6a0627021f6-fe86ab10a43f3d18
  - **Headings:** USER, Gemini-3.5-flash,   - Phase 1: Atom Type Classification,   - Phase 2: Grid Sampling & Feature Matrix Generation,     - 1. Prerequisite Cubes from Substrate DFT,     - 2. Math Formulation for the Linear System,   - Phase 3: Step-by-Step Execution Pipeline,   - Phase 4: Refinement Strategies
    ... (+247 more)

- `GridFF_RelaxedScan_cpp_notes.md` (62 KB)
  - **Summary:** Document evidence-based analysis of C++ frameworks for implementing high-performance parallel relaxed scanning with spline-based constraint  
  - **Headings:** GridFF Relaxed Scan C++ Implementation Notes,   - **🎯 OBJECTIVE**,   - **🚀 MOLECULAR OPTIMIZATION SYSTEM**,     - **📋 OVERVIEW**,     - **🎯 KEY ACHIEVEMENTS**,   - **📁 CORE FILES AND FUNCTIONALITY**,     - **🔧 TipSplineOptimizer.py** - Main Optimization Engine,     - **🎮 run_tipSpline_scan.py** - Command Interface
    ... (+134 more)

- `RigidBodyAFM.md` (54 KB)
  - **Summary:** See following code maps: /@AFM Simulation: GPU Rigid Body Dynamics, CPU GridFF Relaxation, and Interactive GUI
  - **Headings:** USER, SWE-1.6,   - File Relationships,   - Existing Components Analysis,     - RigidBodyDynamics.py (Core Wrapper),     - Test Scripts Comparison,   - What You Need for AFM Simulation,   - Recommended Starting Point
    ... (+94 more)

- `SurfacePotentialVisualization_discussion.md` (31 KB)
  - **Summary:** We have this explorer of interaction of molecule with rigid substrate. @ExplorerVisPy.py
  - **Headings:** USER, GPT-5.1-Codex,   - Findings (from referenced files & docs),   - What’s needed (per user),   - Proposed implementation plan,   - Key design choices / reuse points,   - Next steps (if you agree), USER
    ... (+31 more)

- `SurfacePotentialVisualization_report.md` (4 KB)
  - **Summary:** - Render a seamless periodic isoheight surface z(x,y) for substrate–probe interaction, using the existing, validated backends (GridFF npy, I 
  - **Headings:** Post-mortem (surface potential visualization),   - Goal,   - What we did,   - Problems encountered (root causes),   - Fixes applied,   - Validation results (rigid, NaCl_1x1_L3),   - Key takeaways / guardrails to prevent repeat mistakes,   - Next actions (if needed)


## `doc/Topics/PathDiffusion/`

- `ManipulationPathOpt_Report.md` (9 KB)
  - **Summary:** Based on the `PathDiffusion_design.md` and related discussion notes, the goal was to simulate the manipulation of a molecule over a surface  
  - **Headings:** ManipulationPathOpt Implementation Report,   - 1. Overview and Physics Context,   - 2. The Algorithmic Strategy,     - Asymmetric Rest-Length Causal Tethering,     - Annealing and Bifurcation (Snap) Detection,   - 3. Implementation Details,     - `ManipulationPathOpt.py` Orchestrator,     - GPU Kernel Mapping (`relax_multi.cl`)
    ... (+10 more)

- `PathDiffusion_design.md` (15 KB)
  - **Summary:** Here is a detailed architectural design document and implementation plan for the **"Asymmetric Causal Band Softening & Dynamic n-section"**  
  - **Headings:** System Architecture & Implementation Plan,   - 1. Global Memory Map & 2D Grid Setup,   - 2. New OpenCL Kernels to Implement,     - A. `init_trajectory_spline_kernel`,     - B. `getTipMorse_kernel`,     - C. `reduce_trajectory_gaps_kernel`,   - 3. Modifications to Existing Kernels,     - A. `updateAtomsMMFFf4` (The Propagator & Causal Band)
    ... (+11 more)

- `PathDiffusion_discussion.md` (80 KB)
  - **Summary:** I'm thing about simulation molecular manipulation of molecule (such as PTCDA - perylene dicarboxilic anhydride) on ionic substrate by AFM (i 
  - **Headings:** USER,   - From the paper,     - 4. Configuration Sampling on GPU,   - Your job, Gemini-3.1-pro,     - Part 1: Physical Rigor & The "Teleportation vs. Snap" Paradox,     - Part 2: Numerical Efficiency vs. RL and Shooting Methods,     - Part 3: Implementation Plan / Design Document
    ... (+64 more)

- `PathDiffusion_discussion_old.md` (16 KB)
  - **Summary:** * Diffusion-policy for coupled tip–molecule trajectories in joint state $\mathbf{P} = (\mathbf{x}_{\text{tip}}, \mathbf{x}_{\text{mol}}, \ma * Force-balance guidance: enforce quasi-static equilibrium along the path
  - **Headings:** Summary,     - Trajectory Optimization / Machine-Learning, USER, Gemini-3-flash,   - 1. Does the logic hold up?,   - 2. Is it suitable for Diffusion?,     - Why it works:,   - 3. Challenges & Potential Problems
    ... (+23 more)


## `doc/Topics/RandomGPU/`

- `NormalDistributionGPU.md` (49 KB)
  - **Summary:** What is the fastest way how to generate random number with normal distribution on GPU? Box miller transformation use giniometric functions making it slow. Perhaps buffer with normally distributed numbers?
  - **Headings:** USER, Gemini-3-flash,     - 1. The Ziggurat Algorithm (Fastest Computation),     - 2. Inverse Transform Sampling with Fast Approximations,     - 3. Marsaglia Polar Method (Optimized Box-Muller),     - 4. The "Buffer" Approach (Memory-Bound),     - Summary Comparison Table,     - 1. Why Ziggurat Wins in OpenCL
    ... (+71 more)


## `doc/Topics/RigidBodyAssembly/`

- `RigidBodyAssemblyDiscussion.md` (86 KB)
  - **Summary:** What are possible and efficient ways to generate all possible rigid body rotations If I want to have homogeneous low discrepancy distributio 
  - **Headings:** USER, Grok,     - 1. Random Uniform Sampling (Monte Carlo),     - 2. Low-Discrepancy Deterministic Sampling (Quasi-Monte Carlo),     - 3. Systematic Grid Sampling (Deterministic Lattice/Grid), USER, Gemini-3.1-pro,     - 1. The Mathematical Pipeline (For each OpenCL Thread)
    ... (+83 more)

- `RigidBodyAssemblyDiscussion_tutorial.md` (8 KB)
  - **Summary:** - We sample a large set of rigid-body poses of a molecule inside a periodic lattice, enforcing a 6-fold hexagonal “flower” symmetry. The eva - **Rotation Sampling**: Can be full 3D (Super-Fibonacci sequence), perfectly flat 2D in-plane rotations, or 2D in-plane rotations with smal
  - **Headings:**   - 1) What the module does (scientific overview),   - 2) Key components,     - Host harness: `AssemblyOCL`,     - Transform generation utilities & Diagnostics,     - CLI driver: `tests/tMMFF/test_assembly.py`,   - 3) How to run (examples),     - Triptycene preset (GPU) with Tilt Mode and Diagnostics,     - Helicine with tighter shifts and stricter clash (Full 3D)
    ... (+4 more)

- `RigidGridFF_H2O_NaCl_report.md` (15 KB)
  - **Summary:** description: Rigid GridFF H2O-on-NaCl and PTCDA-on-NaCl multi-replica relaxation report and tutorial 
  - **Headings:** Rigid GridFF — H2O on NaCl and PTCDA on NaCl (multi-replica) — Report and Tutorial,   - 1) Objectives,   - 2) Implemented changes (code),     - Core rigid wrapper,     - H2O runner (single-body baseline),     - PTCDA multi-replica runner (production),   - 3) Workflows and commands,     - A) H2O on NaCl (single body, fast drop + settle)
    ... (+26 more)


## `doc/Topics/STM/`

- `STM_GF_new.md` (105 KB)
  - **Summary:** https://aistudio.google.com/app/prompts?state=%7B%22ids%22:%5B%221u6ceIskV52PNfeEJXMqyAWBGy-xFgKTW%22%5D,%22action%22:%22open%22,%22userId%2 
  - **Headings:** USER, Gemini-3-flash,     - 1. System Partitioning,     - 2. Modeling Source and Drain (Self-Energies),     - 3. The Retarded Green's Function,     - 4. The Path Analysis (Dyson Equation),     - 5. Transmission and Current,     - 6. Computational Efficiency for Scanning
    ... (+135 more)

- `STM_GPU_QMMM.md` (174 KB)
  - **Summary:** I'm doing simulation of STM junction composed of multiple molecules, some molecules are on surface, and one molecule is hanging on a tip (PT 
  - **Headings:** USER, Gemini Flash,     - 1. Formal Framework: The Hamiltonian and Transport,     - 2. Integration with FireCore and Fireball,     - 3. Fast Perturbative Approach (The "Propeller" Model),     - 4. Critical Analysis and Refinements,     - Recommended Next Steps for Implementation:, USER
    ... (+363 more)


## `doc/Topics/TrussSolver/`

- `High_Precision_Physics_on_singlepoint_GPU_hardware.md` (11 KB)
  - **Summary:** When simulating physics, the accuracy of geometric calculations is paramount. On a standard 32-bit `float` GPU, a critical problem arises wh 
  - **Headings:** High-Precision Physics on Single-Precision GPUs: A Design Document,   - Motivation: The Precision Problem in Large Worlds,   - The Core Idea: Split Representation,   - 5. Comparison and Recommendations,     - Recommendation, Method 1: Split-Float (Compensated Arithmetic),   - Integration into the MD Loop,     - Pre-computation Kernel
    ... (+6 more)

- `ProjectiveDynamics.md` (6 KB)
  - **Summary:** Projective Dynamics (PD) is a fast and robust numerical solver for physics-based animation, particularly useful for simulating deformable ob 
  - **Headings:** Projective Dynamics: A Fast and Stable Physics Solver,   - Introduction,   - Mathematical Foundation,     - System Description,     - Energy Formulation,   - Implementation Details,     - System Matrix Construction,     - Right-Hand Side Construction
    ... (+4 more)

- `ProjectiveDynamics_wegl_glsl.md` (86 KB)
  - **Summary:** https://aistudio.google.com/app/prompts?state=%7B%22ids%22:%5B%221SAaaB62dnOueLvFU0AJPwLI9JtRJkFtx%22%5D,%22action%22:%22open%22,%22userId%2 
  - **Headings:**   - USER,   - Gemini,     - 1. Key Architectural Decisions,     - 2. GLSL Implementation,     - 3. JavaScript / WebGL Setup Logic,     - 4. Key Implementation Nuances,   - USER,   - Gemini
    ... (+44 more)

- `VertexBlockDescent.md` (18 KB)
  - **Summary:** This document rewrites article [Vertex Block Descent, A.Chen,(2024)](https://dl.acm.org/doi/10.1145/3658179) focusing on the main VBD optimi 
  - **Headings:** Vertex Block Descent (VBD),   - Purpose,   - Design notes: modular implicit dynamics in `Truss`,   - Notation,   - 3.1 Global optimization,   - 3.2 Local system solver,   - 3.5 Collisions,   - Friction
    ... (+15 more)


## `doc/Topics/molgui_web/`

- `Buckets_GridNeighborDesign.md` (14 KB)
  - **Summary:** - [x] `web/common_js/Buckets.js`: JS bucket graph for storing membership + neighbor buckets (array-based, not typed arrays) - [x] `web/common_js/GridIndexer.js`: grid helpers for mapping positions to cells (if/when used)
  - **Headings:** Design notes for flat Buckets + Grid/AABB neighbor search in MolGUI/CrystalUtils, Status checklist (current JS implementation), Implemented notes and “don’t forget” items (important for future work),   - Index vs ID storage (topology robustness),   - Bucket bounds and visualization correctness,   - Empty bucket pruning,   - MolGUI debug UI / rendering hooks,   - Node test scripts in `web/molgui_web/scripts`
    ... (+37 more)

- `Crystal_Generator.md` (53 KB)
  - **Summary:** I wan to make simple crystal / surface genrator in javascript, something a bit like in ASE (atomistic simulation environment) It should be based on lattice vectors (or angles), few atoms in elementary unit cell, and number of copies. we have then number of replicas 
  - **Headings:**   - USER 1,   - Gemini,   - Prompt for the Coder LLM,     - Part 1: Unit Cell Definitions (The Database),     - Part 2: Surface Generation Logic (Miller Indices),     - Part 3: Vicinal Steps (The "Python Script" Logic),     - Summary of Inputs/Outputs for the Function,   - USER 2
    ... (+44 more)

- `Editor_for_Polymers_on_Surfaces_and_Edges.md` (25 KB)
  - **Summary:** I have 3D editor of molecules in javascript / webGL. What is specific about my use-case is that I focus on simulation of molecules on surfac 
  - **Headings:**   - USER 1, Design Document: Surface-Supported Polymer Assembly System,   - 1. Project Overview,   - 1.1 Implementation decomposition (agreed plan),     - Subtask 1: Surface / lattice replication (pure geometry, no topology),     - Subtask 2: Rigid-body placement/orientation (rotation + translation),     - Subtask 3: Topology operations (bond editing + atom replacement),   - 2. Functional Modules
    ... (+24 more)

- `Molecular_topology_editor_with_explicit_rings.md` (29 KB)
  - **Summary:** I'm thinking to make simple editor of molecules, resp moleculae drawing app. For the begigning consider we work in top view (moleculae are p 
  - **Headings:**   - USER 1,     - 1. Data Structure Definitions,     - 2. Solving the "Membership" Problem,     - 3. Key Algorithms,     - 4. Implementation Workflow Example,     - 5. Why this is efficient,   - Status update (Dec 2025): what is implemented in MolGUI Web now,     - 1) Faster bond search (spatial subdivision buckets)
    ... (+9 more)

- `Molecular_topology_js_plan.md` (14 KB)
  - **Summary:** - [x] Use **explicit hydrogens** as normal atoms (with neighbors etc.). - [x] Selection/reference policy: store **atom IDs** everywhere; only use indices as transient dense-array slots.
  - **Headings:** TODO, Status update (Dec 2025),   - What was implemented,   - Still missing / planned next, Best design (after reading the ring/chain editor doc + [MMFFBuilderBase.h](cci:7://file:///home/prokop/git/FireCore/cpp/common/molecular/MMFFBuilderBase.h:0:0-0:0)),   - Design decision: **Editable graph is AoS (Atom/Bond classes), export/pack is SoA (into existing MoleculeSystem)**, Core requirements distilled from your doc, Proposed JS architecture
    ... (+15 more)

- `molgui_web_dev.md` (15 KB)
  - **Summary:** - **`EditableMolecule` is now the authoritative molecular model** (AoS: `atoms[]` with `pos: Vec3`, `Z`, stable `id`). - Legacy `MoleculeSystem` usage was pushed into a **packed render buffer** only.
  - **Headings:** Session Summary (molgui_web refactor + performance + interaction fixes),   - Achieved / Implemented Changes,     - 1) **Authoritative model migration groundwork**,     - 2) **Legacy IO cleanup**,     - 3) **Math refactor: Vec3 + Mat3**,     - 4) **GUI boilerplate reduction**,     - 5) **Builder panels extracted from GUI.js**,     - 5b) **Polymer + attachment tooling (ported from tests/tAttach)**
    ... (+20 more)

- `Step_Edge_Editor_ideas.md` (0 KB)
  - **Summary:** submodules - drawing 2D top view
  - **Headings:** Step Edge Editor ideas


## `doc/dev_notes/`

- `EFF_otevrena_veda.md` (8 KB)
  - **Summary:** * Elektronove silove pole ktere nahrazuje kvantove interakce mezi elektrony jednoduchym modelovym potencialem mezi (1) elektronem-a-elektron * Silové pole má oproti kvantovým metodám obroskou výhodu v rychlosti a názornosti.
  - **Headings:** Projekt Otevrena Veda, Úvod / Pozadí, Cíl / myšlenka projektu, Problém / Výzva, Postup / Kroky, Kratkodoby vyhled, English Background

- `EFF_otevrena_veda_gpt.md` (12 KB)
  - **Summary:** Elektronové silové pole (Electron Force Field, EFF) je model, který nahrazuje pravé kvantové interakce mezi elektrony jednoduchým potenciále 
  - **Headings:** Projekt Otevřená Věda,   - Úvod / Pozadí,     - Výhody a nevýhody EFF,   - Cíl projektu,     - Jak nato,   - Problémy a výzvy,     - 1. Skryté parametry ($R_i, s_i$),     - 2. Regularizace
    ... (+8 more)

- `xyz_viever.md` (7 KB)
  - **Summary:** You're right to think about modernizing the OpenGL approach! Using shaders, VBOs, and instancing can lead to significant performance improve 


## `doc/dev_notes/EFF/`

- `ECP.md` (18 KB)
  - **Summary:** **Effective Core Pseudopotentials (ECP) in eFF** 
  - **Headings:**   - Particle Types and Interaction Matrix in eFF-ECP,     - Spin Values and Particle Types,   - Detailed Interaction Matrix,     - Abbreviations,     - Function Calls Key,   - Cases

- `eFF_epairs.md` (8 KB)
  - **Summary:** The standard eFF model treats each electron individually. This new mode aims to: 1.  Represent certain electrons as **spin-paired entities (electron pairs)** occupying a single Gaussian orbital.
  - **Headings:** Electron-Pair Mode in eFF: Design and Considerations,   - Motivation,   - Core Concept: `echarge` and `espin`,   - Key Implementation Considerations in Core `eval*` Functions,     - 1. `echarge` and `espin` Initialization,     - 2. Kinetic Energy (`evalKinetic`),     - 3. Coulomb Interactions,     - 4. Pauli Interactions
    ... (+4 more)


## `doc/dev_notes/MMFF/`

- `FFutils_consolidation_summary.md` (8 KB)
  - **Summary:** Successfully consolidated scattered debugging and visualization utilities from XPBD_2D, XPDB_AVBD, and MMFF/UFF into a single reusable modul 
  - **Headings:** FFutils Consolidation - Summary,   - Overview,   - Problem Solved,   - Created Files,     - Core Module,     - Documentation,     - Examples & Tests,   - Key Features
    ... (+27 more)

- `FFutils_migration_guide.md` (11 KB)
  - **Summary:** The `pyBall/FFutils.py` module consolidates scattered debugging and visualization utilities from: - `pyBall/XPBD_2D/XPBD_2D_utils.py`
  - **Headings:** FFutils Migration Guide,   - Overview,   - Key Features,     - 1. Memory-Safe C++ Buffer Wrapping, OLD (scattered in MMFF.py, MMFF_multi.py), NEW (consolidated, memory-safe), Subsequent calls return same object (cached),     - 2. Unified Visualization (2D/3D)
    ... (+58 more)

- `FFutils_visual_debugging_consolidation.md` (48 KB)
  - **Summary:** look we have this two forcefields /home/prokop/git/FireCore/pyBall/XPDB_AVBD
  - **Headings:** USER, GPT-5.1-Codex,     - Where things live now,     - Overlaps & common patterns ripe for refactor,     - Proposed FFutil.py structure (general, numpy-first, no hard dependency on specific sims),     - Refactor targets (what to move into FFutil),     - Interfaces to keep generic (no sim coupling),     - Suggested next steps
    ... (+29 more)

- `MMFF_flat_surface_potential.md` (4 KB)
  - **Summary:** description: Flat surface MMFF coupling notes (Hamaker LJ9-3 / Morse) 
  - **Headings:**   - Purpose,   - Key code locations,   - Parameter / mixing behavior,   - Force direction and logging,   - Trajectory saving,   - Known pitfalls / fixes,   - Quick usage (Python),   - Outstanding items
    ... (+1 more)

- `MMFF_flatsurf_pyOpenCL_progress.md` (9 KB)
  - **Summary:** The objective was to implement and test GPU flat surface potentials (Hamaker LJ9-3 and Morse) in MMFF simulations using PyOpenCL. This invol 
  - **Headings:** GPU Flat Surface Potential (Hamaker + Morse) — PyOpenCL Implementation,   - Background,   - Relevant Files and Functions,     - OpenCL Kernels,     - Python Bindings,     - Test Scripts,   - Problems Encountered,     - 1. Python IndentationError in MolecularDynamics.py
    ... (+23 more)

- `MMFF_variants.md` (87 KB)
  - **Summary:** I would like to make code map for my C++ Classical Forcefield system, it contains two forcefields UFF (universal forcefield) and MMFFsp3 (wh 
  - **Headings:**   - USER,   - GPT-5.1-Codex,   - Module inventory (by layer),   - Execution path (common logical phases),   - Cross-version equivalence (major structures/functions),   - Notable call flows (CPU vs pyOpenCL) for MMFFsp3_loc,   - Parameter/data files usage,   - Suggested parity checklist
    ... (+94 more)


## `doc/dev_notes/Windsurf_Ruls_Workflows/Backup/workflows/`

- `fireball-debugging-guidelines.md` (3 KB)
  - **Summary:** description: Evidence-based Fireball Fortran↔OpenCL parity and debugging 
  - **Headings:**   - Purpose,   - Rules (must follow),   - Workflow,   - Pitfalls to avoid,   - When to stop and fix prerequisites

- `forcefiled-parity-general.md` (4 KB)
  - **Summary:** **Goal:** Achieve numerical parity between a Reference Implementation (e.g., C++/CPU) and a Target Implementation (e.g., OpenCL/CUDA/WebGPU) 
  - **Headings:**   - Skill: scientific_debug_protocol

- `python-pyopencl-workflow-ws.md` (5 KB)
  - **Summary:** description: Python + C/C++ (ctypes) : Maintain and synchronize a C/C++ shared library and its Python `ctypes` wrapper with tests and plots  auto_execution_mode: 1
  - **Headings:** Workflow: Python + C/C++ (ctypes),   - Purpose,   - Workflow Steps,   - Repository pattern (generic) + examples,   - Mapped buffers and switches (general pattern),   - ABI / interop rules,   - Notes,   - Sync / change workflow

- `review.md` (1 KB)
  - **Summary:** auto_execution_mode: 0 description: Review code changes for bugs, security issues, and improvements

- `UFF-MMFF-parity.md` (7 KB)
  - **Summary:** 1. Build CPU reference and GPU targets. Rebuild C++ lib if switches/API changed (`cpp/Build-opt/libs/Molecular/make MMFF_lib`). 2. Disable nonbonded/subtraction: `setSwitches(NonBonded=-1, NonBondNeighs=-1, SurfAtoms=-1, GridFF=-1)` and `setSwitchesUFF_NB(NonBonded=-1
  - **Headings:**     - 5. Parity Run Checklist (UFF),     - 6. Parity Run Checklist (MMFF, MD path),   - Skill: firecore_parity_expert,     - 1. Controlling Physics (The Switches),     - 2. Accessing C++ Buffers (Input Parity),     - 3. OpenCL Debugging (Tracing),     - 4. Common FireCore Pitfalls,     - 5. **MMFF parity (CPU ↔ C++ OCL ↔ PyOpenCL)**


## `doc/dev_notes/Windsurf_Ruls_Workflows/skills/forcefield-validation/`

- `SKILL.md` (1 KB)
  - **Summary:** name: forcefield-validation description: Use when implementing or debugging interatomic force-fields (e.g., UFF, MMFF)
  - **Headings:**   - Architecture,   - Switch Semantics,   - Buffer Parity,   - Pitfalls,   - Bonded-Only Parity Flow


## `doc/dev_notes/Windsurf_Ruls_Workflows/skills/gpu-debugging/`

- `SKILL.md` (1 KB)
  - **Summary:** name: gpu-debugging description: Use when debugging GPU/OpenCL kernels or CUDA code
  - **Headings:**   - Gated Debug Macros,   - Synchronized Tracing,   - Barrier Deadlock,   - GPU-Specific Pitfalls


## `doc/dev_notes/Windsurf_Ruls_Workflows/skills/numerical-parity/`

- `SKILL.md` (2 KB)
  - **Summary:** name: numerical-parity description: Use when comparing two calculations or simulation methods
  - **Headings:**   - Protocol,   - Bisect Strategy,   - Structural vs. Numeric Divergence,   - Discrepancy Tracking


## `doc/dev_notes/Windsurf_Ruls_Workflows/skills/port-to-opencl/`

- `SKILL.md` (2 KB)
  - **Summary:** name: port-to-opencl description: Use when porting Fortran/C++ compute kernels to PyOpenCL
  - **Headings:**   - Why PyOpenCL First,   - Base Class: OpenCLBase,   - Porting Principles (Fortran→OpenCL),   - GPU Kernel Perforamnce Guidelines


## `doc/dev_notes/Windsurf_Ruls_Workflows/skills/python-native-bindings/`

- `SKILL.md` (2 KB)
  - **Summary:** name: python-native-bindings description: Use when maintaining Python bindings to C/C++/Fortran shared libraries via ctypes
  - **Headings:**   - Setup,   - Memory Layout: Fortran vs C/C++/Python,   - Buffer Wrapping (C++ → NumPy zero-copy), Usage: live view into C++ memory,   - Existing Bindings Examples, Disable nonbonded for parity testing, Then read ndims to know shapes, then getBuff(),   - Testing


## `doc/dev_notes/Windsurf_Ruls_Workflows/skills/python-performance/`

- `SKILL.md` (2 KB)
  - **Summary:** name: python-performance description: Use when writing Python code for scientific computing
  - **Headings:**   - Core Principle,   - Rules,     - 1. Batch Everything,     - 2. Vectorized Operations Only,     - 3. Minimal Allocation,     - 4. Advanced Slicing,   - Common Anti-Patterns,     - ❌ Grid Iteration
    ... (+13 more)


## `doc/dev_notes/Windsurf_Ruls_Workflows/skills/visual-debugging/`

- `SKILL.md` (2 KB)
  - **Summary:** name: visual-debugging description: Use when creating diagnostic plots or visualizations for debugging
  - **Headings:**   - Shared Utilities,   - Test Artifacts,   - Visual Review,   - Diagnostics,   - Consolidation Principle


## `doc/dev_notes/Windsurf_Ruls_Workflows/workflows/`

- `post-inventory.md` (1 KB)
  - **Summary:** description: Consolidation and refactoring review after coding to remove duplication and enforce separation of concerns 
  - **Headings:**   - Purpose,   - Workflow,   - Stop Conditions

- `pre-inventory.md` (1 KB)
  - **Summary:** description: Inventory existing code before writing new code to avoid duplication and find reusable components 
  - **Headings:**   - Purpose,   - Workflow


## `doc/dev_notes/pyocl_MMFF/`

- `pyocl_MMFF.md` (49 KB)
  - **Summary:** * I'm trying to move my C/C++ OpenCL interface to python pyOpenCL interface fro Molecular dynamics simulations on GPU. * Our approach is a bit special because it simulates multiple replicas of the same system, at the same time in order to exploit the GPU para
  - **Headings:**   - Task,   - class AtomicSystem in atomicUtils.py,   - MMFFsp3_loc,   - MMFFBuilder::toMMFFsp3_loc(),   - MolWorld_sp3_multi (multi system buffers),   - MolWorld_sp3_multi::run_ocl_opt(),   - C++ interface to allocation of OpenCL buffers,   - OpenCL kernel headers


## `doc/py/CoarseGrainedFF/`

- `GayBerne.md` (12 KB)
  - **Summary:** https://aistudio.google.com/app/prompts?state=%7B%22ids%22:%5B%221Mkp6QAPY5Gc34eGEVS7U0CS1WscFQdhb%22%5D,%22action%22:%22open%22,%22userId%2 
  - **Headings:**   - USER,   - Shadertory Adjugate  https://www.shadertoy.com/view/3s33zj,   - Gemini-3-pro,     - 1. The Standard Solution: The Gay-Berne Potential,     - 2. The "Rounded Disk" SDF Approach (Using your Shader idea),     - 3. The Best Practical Compromise: The "Three-Sphere" Model,     - 4. Implementation Idea: Analytical "Disk-Point" Potential,     - Summary Recommendation
    ... (+5 more)


## `doc/py/FoldedAtomicFunctions/`

- `FoldedAtomicFunction.md` (11 KB)
  - **Summary:** This document describes a novel computational framework for efficiently representing non-covalent interactions between molecules and crystal 
  - **Headings:** Surface Potential Fitting Framework: Design Document,   - Executive Summary,   - 1. Background and Motivation,     - 1.1 Problem Statement,     - 1.2 Current Limitations,   - 2. Proposed Method: Basis Set Compression,     - 2.1 Core Concept,     - 2.2 Mathematical Framework
    ... (+22 more)

- `FoldedAtomicFunction_optimization.md` (41 KB)
  - **Summary:** Look this is how we evaluate folded basis to represent surface ptoential 
  - **Headings:** USER, Gemini-3-Flash,     - 1. The Strategy: Fused Harmonic Generation,     - 2. Fast Polynomial Approximations (The "Base"),     - 3. Optimized Implementation (OpenCL),     - 4. Efficient Gradient Calculation,     - Recommended Workflow:, USER
    ... (+17 more)

- `FoldedAtomicFunction_rigorous_test_plan.md` (35 KB)
  - **Summary:** OK, what we want to do now: 
  - **Headings:** USER, GPT-5.1-Codex, USER, GPT-5.1-Codex,   - Code map (reuse),   - Proposed harness,   - Component separation (minimal changes),   - Outputs (rigorous/publishable)
    ... (+25 more)

- `OptimalAnalyticalBasiset.md` (6 KB)
  - **Summary:** We are given: 1.  A set of $M$ sample functions, $\{y_j(x)\}_{j=1}^M$. Each $y_j(x)$ is typically known at $N_x$ discrete points $x_i$, so we have data $(
  - **Headings:**   - Problem: Finding an Optimal Analytical Basis Set,   - Proposed Solution: Two-Step Approach using SVD on Coefficients

- `Optimize_basis.md` (10 KB)
  - **Summary:** This document outlines the principles, implementation details, and considerations for the basis set optimizer found in `Optimize_basis_new.p 
  - **Headings:** Basis Set Optimizer for Folded Atomic Functions,   - 1. Core Principle: Monte Carlo Optimization,   - 2. Problem Specifics: Basis Set for Morse Potentials,     - 2.1. Basis Definition,     - 2.2. Cost Function,   - 3. Mutation Strategies,   - 4. Important Considerations and Insights,     - 4.1. Two-Phase Optimization
    ... (+7 more)

- `optimize_z_basis.md` (16 KB)
  - **Summary:** This document describes a method and its implementation in `optimize_z_basis.py` for finding an optimal, compact set of analytical basis fun 
  - **Headings:** Optimal Analytical Basis along Z-Direction,   - 1. Motivation and Goals,   - 2. Theoretical Background,     - 2.1. Problem Statement,     - 2.2. Library Basis Representation,     - 2.3. Optimal Basis via Singular Value Decomposition (SVD),     - 2.4. Constructing the Optimal Basis,   - 3. Implementation Details and Workflow (`optimize_z_basis.py`)
    ... (+5 more)

- `tutorial_folded_basis_pareto.md` (6 KB)
  - **Summary:** "description": "Folded basis Pareto search tutorial" 
  - **Headings:** Folded Atomic Functions: Pareto Search Tutorial,   - Files and roles,   - Quickstart: run a batch and plot the Pareto front, Run a batch sweep (writes/updates results.jsonl by default), Plot Pareto front and per-basis fits from the DB, Outputs: pareto_scatter.png, pf_<idx>_* images in the current directory,   - Inputs,   - Outputs
    ... (+5 more)


## `doc/py/LLCAO1D/`

- `Localized_Solver_1D.md` (14 KB)
  - **Summary:** In standard Linear Combination of Atomic Orbitals (LCAO) methods, molecular orbitals (MOs) $\psi_i$ are typically found by solving the gener 
  - **Headings:** Finding Localized Molecular Orbitals via Energy Minimization and Constraints,   - Introduction,   - Orbital Representation,   - Energy Minimization Objective,   - Localization potential vs. hard cutoff,   - Constraints,     - 1. Normalization,     - 2. Orthogonality
    ... (+5 more)

- `OrbitalOrthogonalization.md` (20 KB)
  - **Summary:** This document outlines an efficient method for computing the orthogonalization update (or gradient) for a set of orbitals expressed in a loc 
  - **Headings:**   - Efficient Computation of Orthogonalization Gradients for Localized Orbitals,     - 1. Problem Definition,     - 2. The Algebraic View: Two Computational Paths,     - 3. Complexity Analysis: Why Path Matters,     - 4. The Practical Implementation: A Linear-Scaling Algorithm,     - 5. Illustrative Python Code, --- 1. Setup: Create sparse mock data ---, C: (N_orb, N_basis) coefficient matrix, sparse
    ... (+20 more)


## `doc/py/OrderN_QM/`

- `Approximate_Overlap.md` (31 KB)
  - **Summary:** https://aistudio.google.com/app/prompts?state=%7B%22ids%22:%5B%2213BfBgHtNbkazBS-RZIc0wpSbfLnUqEKW%22%5D,%22action%22:%22open%22,%22userId%2 
  - **Headings:**   - USER,     - 1. The Core Concept: Zero Differential Overlap (ZDO),     - 2. The Hierarchy: CNDO vs. INDO vs. NDDO,     - 3. Comparison of Modern NDDO-based Methods,     - 4. Beyond NDDO: The Linear Scaling Perspective,     - Summary Advice for Your Implementation,   - USER,     - 1. Clarification: Differential Overlap vs. Density Matrix
    ... (+15 more)

- `BlockCholesky.md` (7 KB)
  - **Summary:** https://aistudio.google.com/app/prompts?state=%7B%22ids%22:%5B%221xADxx5mysl_VQtkkYMRydiAlTU-2FRnS%22%5D,%22action%22:%22open%22,%22userId%2 
  - **Headings:**   - USER,   - Gemini,     - 1. Factorization: The "Tiled" Approach (Your "Cluster" Idea),     - 2. The Solve Step: Parallelizing Substitution,     - 3. Your "Smoothing Residuals" Idea (Domain Decomposition),     - Summary of recommendations for your implementation

- `CheFSI_Frontier_Orbitals.md` (66 KB)
  - **Summary:** https://aistudio.google.com/app/prompts?state=%7B%22ids%22:%5B%221BbJlsY2lhVkragBxr38FduksEM3FM7wc%22%5D,%22action%22:%22open%22,%22userId%2 
  - **Headings:**   - USER,   - Gemini,     - 1. The Core Problem: Shift-and-Invert,     - 2. The Solution: Chebyshev Filtering (CheFSI),     - 3. The "Folded Spectrum" Method (Gradient Descent approach),     - 4. LOBPCG (Locally Optimal Block Preconditioned Conjugate Gradient),     - Summary: Which one should you use?,     - Recommendation for your implementation
    ... (+48 more)

- `DensityMatrix_Idenpotency_NearestNeighbor.md` (44 KB)
  - **Summary:** https://gemini.google.com/share/245e57bba284 
  - **Headings:**   - USER,   - Gemini-Flash,   - 1. The Pauli Principle as Idempotency,   - 2. Linear Scaling and \"Nearsightedness\",     - How it is enforced locally:,   - 3. Local Constraints vs. Global Trace,     - Summary Table: MO vs. Density Matrix,   - USER
    ... (+57 more)

- `KekuleOrderN.md` (70 KB)
  - **Summary:** I'm thinking about connection between Haris functional and pierls distortion. Piers distortion in chain of carbon atoms make them dimerize w 
  - **Headings:** USER, Gemini-3-flash,   - 1. The Paradox: What Makes the Bond if $\rho_{ij} = 0$?,     - But what if we don't diagonalize at all?,   - 2. Applying this to the Hückel $\pi$-System,   - 3. How to Fix It Without Diagonalization: Fragment Selection,   - 4. Can We Estimate $\rho_{ij}(r_{ij})$ Directly (No Diagonalization)?,     - Method A: The Moment Expansion (Bond-Order Potentials)
    ... (+132 more)

- `OMM_Kernels.md` (76 KB)
  - **Summary:** https://aistudio.google.com/app/prompts?state=%7B%22ids%22:%5B%221J80UActrLop6Z4WDxUdlGt75mU9bFA-0%22%5D,%22action%22:%22open%22,%22userId%2 
  - **Headings:**   - USER,   - Gemmini,     - 1. Analysis of the Design,     - 2. PyOpenCL Kernel Implementation,     - 3. Optimization Notes for the Future (Python Wrapper Context),   - USER,   - Gemini,     - 1. Data Structure Design
    ... (+42 more)

- `Orbital_Minimization_Methond.md` (125 KB)
  - **Summary:** https://drive.google.com/file/d/10B4ojUt-DkFaRUYlPv91nnVFsZuO9N1h/view?usp=sharing, https://aistudio.google.com/app/prompts?state=%7B%22ids% 
  - **Headings:**   - USER 1,     - 1. Literature & Similar Methods,     - 2. The Algorithm: Projected Gradient Descent with Localization,     - 3. Recommendation for your Implementation,     - 4. Python Prototype (Constraint Projection Solver), --- Run ---, Plot, Visualize Matrix Sparsity
    ... (+94 more)

- `OrderN_Electronic_Structure_Solver_Gemini.md` (67 KB)
  - **Summary:** https://gemini.google.com/share/8e4c23a2107b 
  - **Headings:**   - USER 1,   - 1. Minimal Basis & Slater-Koster Framework,   - 2. Grid-Bin Neighbor Searching,   - 3. Sigma-Pi Orthogonality,   - 4. Sigma Space: Polynomial Expansion,   - 5. Pi-System: Green's Functions,     - Potential Challenges to Consider,   - USER 2
    ... (+76 more)

- `wrong_hopping.md` (6 KB)
  - **Summary:** The error is not in the implementation of the solver logic (which correctly handles the generalized eigenvalue problem), but rather in the * 
  - **Headings:**     - The Explanation,     - The Fix, ... (around line 81), ...,   - USER,     - 1. The Scaling Law (Harrison's Rule),     - 2. The Constants,     - 3. Calculating Reasonable Parameters
    ... (+8 more)


## `doc/py/PairLoadBalancing/`

- `PairLoadBalancing.md` (8 KB)
  - **Summary:** OK, I modified it like this. It works to some degree but the memory reuse it not so great, it seems it could be better. I think the Z-curve  try to brainstorm a bit. I guess this is rather common task and many strategies / improvements were proposed
  - **Headings:**   - USER,   - Gemini,     - Why the Midpoint Strategy is "Sub-optimal" for DM,     - Brainstorming Better Strategies for DM Pair-Chunking,     - Comparison Script: Midpoint vs. Lexicographical vs. 2D Block, =============================================================================, PAIR SORTING HEURISTICS, =============================================================================
    ... (+7 more)


## `doc/py/ParticleOnSpring/`

- `README.md` (8 KB)
  - **Summary:** This directory contains educational Python scripts that demonstrate the physics of particle manipulation using spring-based force microscopy 
  - **Headings:** Particle on Spring: Hysteresis and Manipulation Trajectory Demonstrations,   - Overview,   - Physical Background,     - Double-Well Potential,     - Spring-Particle System,     - Total Energy,     - Hysteresis Origin,   - Scripts Description
    ... (+24 more)


## `doc/py/SoftClampFit/`

- `SoftClampFit.md` (7 KB)
  - **Summary:** The primary objective when parameterizing an empirical potential (or "force field") is to create a model that accurately reproduces referenc 
  - **Headings:**   - Robust and Chemically-Relevant Fitting of Non-Covalent Interaction Potentials,     - 1. The Challenge of Fitting Empirical Potentials,     - 2. A Robust Error Function using Soft Clamping,     - 3. Rigorous Optimization: The Variational Derivative


## `doc/py/SparseMatrixToGrid/`

- `FittingLCAO.md` (59 KB)
  - **Summary:** https://aistudio.google.com/app/prompts?state=%7B%22ids%22:%5B%221NHWTpNMbbj3K5MUwVTvTFvxUP62UqZHi%22%5D,%22action%22:%22open%22,%22userId%2 
  - **Headings:**   - USER,   - Gemini,     - Part 1: The Operators (Projection and Back-Projection),     - Part 2: Optimization Algorithms (Solver Layer),     - Summary Recommendation,   - USER,   - Gemini,     - Didactic Introduction
    ... (+20 more)

- `Quadrature_Weights_optimization.md` (26 KB)
  - **Summary:** https://aistudio.google.com/app/prompts?state=%7B%22ids%22:%5B%221uS4jr2NwrA-QJIimtJND-_vEWi4jWSlr%22%5D,%22action%22:%22open%22,%22userId%2 
  - **Headings:**   - USER,   - Gemini,     - Mathematical Formulation,     - Python Implementation, Import the geometry generator from the previous script, ==========================================, 1. Basis Function Definitions, ==========================================
    ... (+28 more)


## `doc/py/SpatialPartitioning/`

- `Hierarchical_PCA_Particle_Clustering.md` (28 KB)
  - **Summary:** https://gemini.google.com/share/66aeeecb1c9d 
  - **Headings:**   - USER 1,     - The Algorithmic Flow,     - Refining the "Outlier" and "Equal Size" Problem,     - Comparison: PCA-Splitting vs. Standard K-Means,     - Implementation Tip: Memory Management, Conceptual logic for the split,   - USER 2,     - The Parallel Implementation Pipeline
    ... (+39 more)


## `doc/py/WrapedSamplinggGrid/`

- `WrapedSampling_tutorial.md` (5 KB)
  - **Summary:** We need training/diagnostic samples to **fit a force field against DFT references** for a molecule on a stepped NaCl surface. The workflow b 
  - **Headings:**   - Tutorial,     - Why we are doing this (purpose & motivation),     - Physical & mathematical background (quick recap),     - What this does,     - Main script for combined sampling: [sample_driver.py](cci:7://file:///home/prokophapala/git/FireCore/doc/py/WrapedSamplinggGrid/sample_driver.py:0:0-0:0),     - Angular-only visualization / orientation movie: [tilt_sampling_fan.py](cci:7://file:///home/prokophapala/git/FireCore/doc/py/WrapedSamplinggGrid/tilt_sampling_fan.py:0:0-0:0),     - Output naming conventions,     - Controlling u-levels with `u_list`
    ... (+2 more)

- `WrapedSamplingGrid_discussion.md` (17 KB)
  - **Summary:** https://aistudio.google.com/app/prompts?state=%7B%22ids%22:%5B%221g7hpcRmxV35Ge1Im4rC4NNQanI9SJ_ON%22%5D,%22action%22:%22open%22,%22userId%2 
  - **Headings:** USER, Gemini-3.1-pro,     - 1. The Riemannian Geometry Approach (Inspired by General Relativity),     - 2. The Optimal Transport Approach (The Ultimate Jacobian Mapping),     - 3. Gradient Flow / Continuous Normalizing Flows (Aerodynamics / Diffeomorphisms),     - Summary of Recommendation, USER, Gemini-3.1-pro
    ... (+29 more)


## `doc/py/wavelets/`

- `Wavelet_Gemini.md` (29 KB)
  - **Summary:** I'm thinking about using wavelets to represent molecular orbitals in my fast linear scalling code. I know that expansion of MOs in wavelets  
  - **Headings:**   - USER,     - 1. The Choice of Wavelet,     - 2. The Transform Algorithm (Grid $\leftrightarrow$ Wavelet),     - 3. Grid Strategy: Uniform vs. Multi-resolution,     - 4. Special Considerations for STM/AFM,     - 5. Implementation Roadmap (The "Fast" Path),   - USER 2,     - Why your "Mexican Hat" intuition leads to infinite spread
    ... (+30 more)


## `doc/skeletons/apps/EFF/`

- `EFFapp.cpp.skeleton.md` (1 KB)
  - **Summary:** * `E0` * `Efunc`
  - **Headings:** `EFFapp.cpp`, Classes,     - `TestAppRARFF`, Free Functions

- `test_CLCFGO_OGL3.cpp.skeleton.md` (0 KB)
  - **Summary:** * `die` * `draw`
  - **Headings:** `test_CLCFGO_OGL3.cpp`, Free Functions


## `doc/skeletons/apps/MolecularEditor/`

- `ConfSearch.cpp.skeleton.md` (0 KB)
  - **Summary:** * `atomSize` * `atomdist`
  - **Headings:** `ConfSearch.cpp`, Classes,     - `AppMolecularEditor2`, Free Functions

- `FireCoreMMFFmini.cpp.skeleton.md` (1 KB)
  - **Summary:** * `MOgrid` * `MOgrid`
  - **Headings:** `FireCoreMMFFmini.cpp`, Classes,     - `TestAppMMFFmini`, Free Functions

- `FireCoreVisual.cpp.skeleton.md` (2 KB)
  - **Summary:** * `MOgrid` * `MOgrid`
  - **Headings:** `FireCoreVisual.cpp`, Classes,     - `TestAppFireCoreVisual`, Free Functions

- `MolecularBrowser.cpp.skeleton.md` (1 KB)
  - **Summary:** * `atypeZ` * `atypes`
  - **Headings:** `MolecularBrowser.cpp`, Classes,     - `TestAppMolecularBrowser`, Free Functions

- `MolecularEditor_main.cpp.skeleton.md` (0 KB)
  - **Summary:** * `atomSize` * `builder`
  - **Headings:** `MolecularEditor_main.cpp`, Classes,     - `AppMolecularEditor2`, Free Functions

- `MolecularEditor_old_main.cpp.skeleton.md` (0 KB)
  - **Summary:** * `converged` * `fmaxConverg`
  - **Headings:** `MolecularEditor_old_main.cpp`, Classes,     - `MolecularEditorApp`, Free Functions

- `MolGUIapp.cpp.skeleton.md` (0 KB)
  - **Summary:** * `main`
  - **Headings:** `MolGUIapp.cpp`, Free Functions

- `MolGUIapp_argv.h.skeleton.md` (0 KB)
  - **Summary:** (no prose lines)
  - **Headings:** `MolGUIapp_argv.h`

- `MolGUIapp_Lua.h.skeleton.md` (0 KB)
  - **Summary:** * `initMyLua` * `l_addGUIpanel`
  - **Headings:** `MolGUIapp_Lua.h`, Free Functions

- `MolGUIapp_QMMM.cpp.skeleton.md` (0 KB)
  - **Summary:** * `main`
  - **Headings:** `MolGUIapp_QMMM.cpp`, Free Functions

- `test_MMFF.cpp.skeleton.md` (0 KB)
  - **Summary:** * `bakeMMFF` * `initRigidSubstrate`
  - **Headings:** `test_MMFF.cpp`, Free Functions


## `doc/skeletons/apps/RFF/`

- `RFFapp.cpp.skeleton.md` (1 KB)
  - **Summary:** * `bBlockAddAtom` * `bConsole`
  - **Headings:** `RFFapp.cpp`, Classes,     - `TestAppRARFF`, Free Functions


## `doc/skeletons/apps_OCL/MolecularEditorOCL/`

- `FireCoreVisualOCL.cpp.skeleton.md` (2 KB)
  - **Summary:** * `MOgrid` * `MOgrid`
  - **Headings:** `FireCoreVisualOCL.cpp`, Classes,     - `TestAppFireCoreVisual`, Free Functions

- `MolecularEditorOCL_main.cpp.skeleton.md` (1 KB)
  - **Summary:** * `atomSize` * `atom_count`
  - **Headings:** `MolecularEditorOCL_main.cpp`, Classes,     - `AppMolecularEditorOCL`, Free Functions

- `MolecularEditorOCL_scanner.cpp.skeleton.md` (1 KB)
  - **Summary:** * `atomSize` * `atom_count`
  - **Headings:** `MolecularEditorOCL_scanner.cpp`, Classes,     - `AppMolecularEditorOCL`, Free Functions

- `MolecularWorl2OCL.h.skeleton.md` (2 KB)
  - **Summary:** * `alpha` * `cl`
  - **Headings:** `MolecularWorl2OCL.h`, Classes,     - `GridFF_OCL`,     - `RigidMolecularWorldOCL`, Free Functions

- `MolGUIapp_multi.cpp.skeleton.md` (0 KB)
  - **Summary:** * `main`
  - **Headings:** `MolGUIapp_multi.cpp`, Free Functions

- `MolGUIapp_multi_Lua.h.skeleton.md` (0 KB)
  - **Summary:** * `initMyLua_multi` * `l_getReplicaCount`
  - **Headings:** `MolGUIapp_multi_Lua.h`, Free Functions

- `MolGUIapp_QMMM_multi.cpp.skeleton.md` (0 KB)
  - **Summary:** * `main`
  - **Headings:** `MolGUIapp_QMMM_multi.cpp`, Free Functions


## `doc/skeletons/common/`

- `argparse.h.skeleton.md` (0 KB)
  - **Summary:** * `func` * `n`
  - **Headings:** `argparse.h`, Classes,     - `ArgFunc`, Free Functions

- `arrayAlgs.h.skeleton.md` (0 KB)
  - **Summary:** * `i0` * `n`
  - **Headings:** `arrayAlgs.h`, Classes,     - `I0n`, Free Functions

- `browser.h.skeleton.md` (0 KB)
  - **Summary:** * `bMyDownScale` * `curCol`
  - **Headings:** `browser.h`, Classes,     - `Browser`

- `commandTree.h.skeleton.md` (0 KB)
  - **Summary:** * `func` * `info`
  - **Headings:** `commandTree.h`, Classes,     - `CommandNode`,     - `CommandTree`

- `constants.h.skeleton.md` (0 KB)
  - **Summary:** (no prose lines)
  - **Headings:** `constants.h`

- `datatypes.h.skeleton.md` (1 KB)
  - **Summary:** * `x` * `y`
  - **Headings:** `datatypes.h`, Classes,     - `double2`,     - `double4`,     - `double5`,     - `double8`,     - `float16`,     - `float2`
    ... (+7 more)

- `datatypes_utils.h.skeleton.md` (0 KB)
  - **Summary:** * `add_except` * `copy`
  - **Headings:** `datatypes_utils.h`, Free Functions

- `globals.h.skeleton.md` (0 KB)
  - **Summary:** (no prose lines)
  - **Headings:** `globals.h`

- `IO_utils.h.skeleton.md` (0 KB)
  - **Summary:** * `data` * `dtype`
  - **Headings:** `IO_utils.h`, Classes,     - `NumpyFile`, Free Functions

- `libUtils.h.skeleton.md` (0 KB)
  - **Summary:** * `getBBuff` * `getBuff`
  - **Headings:** `libUtils.h`, Free Functions

- `LuaClass.h.skeleton.md` (0 KB)
  - **Summary:** * `L` * `L`
  - **Headings:** `LuaClass.h`, Classes,     - `LuaClass`

- `LuaHelpers.h.skeleton.md` (0 KB)
  - **Summary:** * `Lua::checkError` * `Lua::clean`
  - **Headings:** `LuaHelpers.h`, Free Functions

- `macroUtils.h.skeleton.md` (0 KB)
  - **Summary:** * `data` * `n`
  - **Headings:** `macroUtils.h`, Classes,     - `Arr`,     - `Buf`, Free Functions

- `repl.h.skeleton.md` (0 KB)
  - **Summary:** * `functions` * `functions`
  - **Headings:** `repl.h`, Classes,     - `Interpreter`,     - `REPL::Interpreter`, Free Functions

- `script_parser.h.skeleton.md` (0 KB)
  - **Summary:** * `functions` * `nlevels`
  - **Headings:** `script_parser.h`, Classes,     - `ScripParser`,     - `code_segment`

- `ScriptLang.h.skeleton.md` (0 KB)
  - **Summary:** * `args` * `body`
  - **Headings:** `ScriptLang.h`, Classes,     - `Function`,     - `ScriptLang`,     - `Variable`

- `test_globals.h.skeleton.md` (0 KB)
  - **Summary:** * `data` * `data`
  - **Headings:** `test_globals.h`, Classes,     - `NDArray`

- `testUtils.h.skeleton.md` (0 KB)
  - **Summary:** * `T` * `t1`
  - **Headings:** `testUtils.h`, Classes,     - `StopWatch`, Free Functions


## `doc/skeletons/common/OpenCL/`

- `OCL.h.skeleton.md` (2 KB)
  - **Summary:** * `flags` * `imageFormat`
  - **Headings:** `OCL.h`, Classes,     - `OCLBuffer`,     - `OCLarg`,     - `OCLsystem`,     - `OCLtask`

- `OCL_device_picker.h.skeleton.md` (0 KB)
  - **Summary:** * `getDeviceList` * `getDeviceName`
  - **Headings:** `OCL_device_picker.h`, Free Functions

- `OCL_DFT.h.skeleton.md` (3 KB)
  - **Summary:** * `dim` * `global`
  - **Headings:** `OCL_DFT.h`, Classes,     - `KernelDims`,     - `OCL_DFT`, Free Functions

- `OCL_EFF.h.skeleton.md` (0 KB)
  - **Summary:** * `KRSrho` * `KRSrho`
  - **Headings:** `OCL_EFF.h`, Classes,     - `OCL_EFF`

- `OCL_MM.h.skeleton.md` (4 KB)
  - **Summary:** * `GFFparams` * `GFFparams`
  - **Headings:** `OCL_MM.h`, Classes,     - `OCL_MM`,     - `cl_Mat3`, Free Functions

- `OCL_PP.h.skeleton.md` (1 KB)
  - **Summary:** * `Rdamp` * `cl_invLvec`
  - **Headings:** `OCL_PP.h`, Classes,     - `OCL_PP`,     - `cl_Mat3`, Free Functions

- `OCLerrors.h.skeleton.md` (0 KB)
  - **Summary:** * `OCL_buildProgramFailure` * `OCL_check_error`
  - **Headings:** `OCLerrors.h`, Free Functions

- `OCLfft_errors.h.skeleton.md` (0 KB)
  - **Summary:** * `OCLfft_buildProgramFailure` * `OCLfft_check_error`
  - **Headings:** `OCLfft_errors.h`, Free Functions


## `doc/skeletons/common/dataStructures/`

- `arrayAlloc.h.skeleton.md` (0 KB)
  - **Summary:** * `arrays` * `nbtot`
  - **Headings:** `arrayAlloc.h`, Classes,     - `ArrayAllocator`

- `Buckets.h.skeleton.md` (0 KB)
  - **Summary:** * `cell2obj` * `cellI0s`
  - **Headings:** `Buckets.h`, Classes,     - `Buckets`

- `Buckets2D.h.skeleton.md` (0 KB)
  - **Summary:** * `Rcut` * `epsilon`
  - **Headings:** `Buckets2D.h`, Classes,     - `Buckets2D`

- `Buckets3D.h.skeleton.md` (0 KB)
  - **Summary:** * `bResizing` * `neighs_in`
  - **Headings:** `Buckets3D.h`, Classes,     - `Buckets3D`

- `BucketsND.h.skeleton.md` (0 KB)
  - **Summary:** * `n` * `pmax`
  - **Headings:** `BucketsND.h`, Classes,     - `BoxND`,     - `BucketsND`

- `CircularArray.h.skeleton.md` (0 KB)
  - **Summary:** * `bOwner` * `data`
  - **Headings:** `CircularArray.h`, Classes,     - `CircularArray`

- `CMesh.h.skeleton.md` (0 KB)
  - **Summary:** * `edges` * `faces`
  - **Headings:** `CMesh.h`, Classes,     - `CMesh`

- `containers.h.skeleton.md` (0 KB)
  - **Summary:** * `map` * `vec`
  - **Headings:** `containers.h`, Classes,     - `Dict`, Free Functions

- `CubicRuler.h.skeleton.md` (0 KB)
  - **Summary:** * `diRay` * `hRay`
  - **Headings:** `CubicRuler.h`, Classes,     - `CubicRuler`, Free Functions

- `DistanceHierarchy.h.skeleton.md` (0 KB)
  - **Summary:** * `R` * `n`
  - **Headings:** `DistanceHierarchy.h`, Classes,     - `DHLevel`,     - `DHNode`,     - `DHSearch`,     - `DistanceHierarchy`

- `Graph.h.skeleton.md` (0 KB)
  - **Summary:** * `V` * `adj`
  - **Headings:** `Graph.h`, Classes,     - `Graph`, Free Functions

- `Grid.h.skeleton.md` (1 KB)
  - **Summary:** * `cell` * `dCell`
  - **Headings:** `Grid.h`, Classes,     - `GridShape`, Free Functions

- `Grid3D.h.skeleton.md` (0 KB)
  - **Summary:** * `NX` * `NXY`
  - **Headings:** `Grid3D.h`, Classes,     - `Grid3D`,     - `TILE_vector`

- `grids3D.h.skeleton.md` (1 KB)
  - **Summary:** * `invStep` * `n`
  - **Headings:** `grids3D.h`, Classes,     - `CubeGridRuler`,     - `CubeGridRulerUnbound`,     - `GridRulerInterface`, Free Functions

- `GridUtils.h.skeleton.md` (0 KB)
  - **Summary:** * `coulombGrid_brute`
  - **Headings:** `GridUtils.h`, Free Functions

- `HashMap2D.h.skeleton.md` (0 KB)
  - **Summary:** * `invStep` * `step`
  - **Headings:** `HashMap2D.h`, Classes,     - `HashMap2D`

- `HashMap_temp.h.skeleton.md` (0 KB)
  - **Summary:** * `DEBUG_counter` * `DEBUG_counter`
  - **Headings:** `HashMap_temp.h`, Classes,     - `HashMap`,     - `HashMapField`, Free Functions

- `LimitedGraph.h.skeleton.md` (0 KB)
  - **Summary:** * `bPrint` * `disc`
  - **Headings:** `LimitedGraph.h`, Classes,     - `LimitedGraph`, Free Functions

- `Map2D.h.skeleton.md` (0 KB)
  - **Summary:** * `invStep` * `nx`
  - **Headings:** `Map2D.h`, Classes,     - `Map2D`

- `MethodDict.h.skeleton.md` (0 KB)
  - **Summary:** * `actionDispatch` * `actionDispatch`
  - **Headings:** `MethodDict.h`, Classes,     - `MethodDict`

- `parsing.h.skeleton.md` (1 KB)
  - **Summary:** * `cb` * `cl`
  - **Headings:** `parsing.h`, Classes,     - `Break`,     - `Node`,     - `ParserItem`,     - `Tuple`,     - `parsing::Break`,     - `parsing::Node`
    ... (+2 more)

- `PointCloudComparator.h.skeleton.md` (1 KB)
  - **Summary:** * `i` * `r2`
  - **Headings:** `PointCloudComparator.h`, Classes,     - `Neigh`,     - `PointCloudComparator`,     - `TypePointsComparator`, Free Functions

- `RBFInterpolator.h.skeleton.md` (0 KB)
  - **Summary:** * `Rcut` * `buckets2D`
  - **Headings:** `RBFInterpolator.h`, Classes,     - `RBFInterpolator`

- `SimplexRuler.h.skeleton.md` (0 KB)
  - **Summary:** * `MAP_OFFSET` * `ROUND_OFFSET`
  - **Headings:** `SimplexRuler.h`, Classes,     - `SimplexRuler`, Free Functions

- `SMILESparser.h.skeleton.md` (0 KB)
  - **Summary:** * `anchors` * `builder`
  - **Headings:** `SMILESparser.h`, Classes,     - `SMILESparser`

- `SortedStrings.h.skeleton.md` (0 KB)
  - **Summary:** * `table` * `table`
  - **Headings:** `SortedStrings.h`, Classes,     - `SortedStrings`

- `SphereTreeND.h.skeleton.md` (0 KB)
  - **Summary:** * `center` * `leafs`
  - **Headings:** `SphereTreeND.h`, Classes,     - `SphereNodeND`,     - `SphereTreeND`, Free Functions

- `StructParser.h.skeleton.md` (0 KB)
  - **Summary:** * `cCLOSE` * `cNEXT`
  - **Headings:** `StructParser.h`, Classes,     - `StructParser`,     - `parsing::StructParser`

- `Table.h.skeleton.md` (0 KB)
  - **Summary:** * `nsub` * `offset`
  - **Headings:** `Table.h`, Classes,     - `Atribute`,     - `Table`

- `TileBuffer3D.h.skeleton.md` (0 KB)
  - **Summary:** * `M` * `MNX`
  - **Headings:** `TileBuffer3D.h`, Classes,     - `TileBuffer3D`

- `Tree.h.skeleton.md` (0 KB)
  - **Summary:** * `content` * `parrent`
  - **Headings:** `Tree.h`, Classes,     - `PTree`,     - `Tree`,     - `TreeGen`


## `doc/skeletons/common/math/`

- `approximation.h.skeleton.md` (1 KB)
  - **Summary:** * `coefs` * `cost`
  - **Headings:** `approximation.h`, Classes,     - `Approx::ApproxBlock`,     - `Approx::ApproxVariant`,     - `Approx::AutoApprox`,     - `Approx::PolyFit`,     - `ApproxBlock`,     - `ApproxVariant`
    ... (+3 more)

- `Besier.h.skeleton.md` (0 KB)
  - **Summary:** * `ns` * `ps`
  - **Headings:** `Besier.h`, Classes,     - `BesierTriangle`, Free Functions

- `Body.cpp.skeleton.md` (0 KB)
  - **Summary:** (no prose lines)
  - **Headings:** `Body.cpp`

- `Body.h.skeleton.md` (1 KB)
  - **Summary:** * `force` * `mass`
  - **Headings:** `Body.h`, Classes,     - `BodyQ`,     - `KinematicBody`,     - `Particle3D`,     - `PointBody`,     - `PoseQ`,     - `RigidBody`
    ... (+2 more)

- `Body2D.cpp.skeleton.md` (0 KB)
  - **Summary:** (no prose lines)
  - **Headings:** `Body2D.cpp`

- `Body2D.h.skeleton.md` (1 KB)
  - **Summary:** * `phi` * `pos`
  - **Headings:** `Body2D.h`, Classes,     - `KinematicBody2D`,     - `PointBody2D`,     - `RigidBody2D`,     - `SpringConstrain2D`

- `Bspline.h.skeleton.md` (1 KB)
  - **Summary:** * `Bspline::basis` * `Bspline::basis5`
  - **Headings:** `Bspline.h`, Free Functions

- `Bspline_fit.h.skeleton.md` (0 KB)
  - **Summary:** * `Bspline::fit1D` * `Bspline::fit1D_EF`
  - **Headings:** `Bspline_fit.h`, Free Functions

- `Bspline_fit_2D.h.skeleton.md` (0 KB)
  - **Summary:** * `Bspline::assembleBound2D_pbc` * `Bspline::assemleBound2D`
  - **Headings:** `Bspline_fit_2D.h`, Free Functions

- `Bspline_fit_3D.h.skeleton.md` (0 KB)
  - **Summary:** * `Bspline::error_iz` * `Bspline::error_iz_pbc`
  - **Headings:** `Bspline_fit_3D.h`, Free Functions

- `Camera.h.skeleton.md` (0 KB)
  - **Summary:** * `aspect` * `persp`
  - **Headings:** `Camera.h`, Classes,     - `CameraT`

- `CG.h.skeleton.md` (0 KB)
  - **Summary:** * `A` * `Ap`
  - **Headings:** `CG.h`, Classes,     - `CG`

- `Convex2d.cpp.skeleton.md` (0 KB)
  - **Summary:** (no prose lines)
  - **Headings:** `Convex2d.cpp`

- `Convex2d.h.skeleton.md` (0 KB)
  - **Summary:** * `corners` * `id`
  - **Headings:** `Convex2d.h`, Classes,     - `Convex2d`

- `Convex3d.h.skeleton.md` (0 KB)
  - **Summary:** * `edges` * `planes`
  - **Headings:** `Convex3d.h`, Classes,     - `Convex2d`

- `CubicBSpline.h.skeleton.md` (0 KB)
  - **Summary:** * `CubicBSpline::basis` * `CubicBSpline::dbasis`
  - **Headings:** `CubicBSpline.h`, Free Functions

- `distribution_generator.h.skeleton.md` (0 KB)
  - **Summary:** * `Eps` * `distributionType`
  - **Headings:** `distribution_generator.h`, Classes,     - `DistributionGenerator`

- `DynamicOpt.cpp.skeleton.md` (0 KB)
  - **Summary:** (no prose lines)
  - **Headings:** `DynamicOpt.cpp`

- `DynamicOpt.h.skeleton.md` (1 KB)
  - **Summary:** * `avs` * `avs2`
  - **Headings:** `DynamicOpt.h`, Classes,     - `DynamicOpt`,     - `OptLog`

- `fastmath.h.skeleton.md` (0 KB)
  - **Summary:** * `clamp` * `clamp_abs`
  - **Headings:** `fastmath.h`, Free Functions

- `Forces.h.skeleton.md` (0 KB)
  - **Summary:** * `R4blob` * `R8down`
  - **Headings:** `Forces.h`, Free Functions

- `Fourier.h.skeleton.md` (0 KB)
  - **Summary:** * `FFT` * `FFT_2D`
  - **Headings:** `Fourier.h`, Free Functions

- `functions.h.skeleton.md` (0 KB)
  - **Summary:** * `cross_valey` * `gridValy`
  - **Headings:** `functions.h`, Free Functions

- `GaussianBasis.h.skeleton.md` (2 KB)
  - **Summary:** * `c` * `p`
  - **Headings:** `GaussianBasis.h`, Classes,     - `Blob`,     - `Gauss::Blob`,     - `Gauss::PairDeriv`,     - `Gauss::PairDerivs`,     - `Gauss::PairInt`,     - `PairDeriv`
    ... (+3 more)

- `general_splines.h.skeleton.md` (0 KB)
  - **Summary:** * `get_splineR2_0246` * `splineR2_0246`
  - **Headings:** `general_splines.h`, Free Functions

- `geom2D.h.skeleton.md` (2 KB)
  - **Summary:** * `ang0` * `dang`
  - **Headings:** `geom2D.h`, Classes,     - `Arc2d`,     - `Circle2d`,     - `Ellipse2D`,     - `Line2d`,     - `Ray2d`,     - `Rect2d`
    ... (+3 more)

- `geom3D.h.skeleton.md` (2 KB)
  - **Summary:** * `a` * `b`
  - **Headings:** `geom3D.h`, Classes,     - `Box`,     - `Capsula3D`,     - `Disk3D`,     - `Ellipsoide`,     - `Line3d`,     - `LineInterval3d`
    ... (+9 more)

- `gonioApprox.h.skeleton.md` (0 KB)
  - **Summary:** * `acos_nVidia` * `acos_table`
  - **Headings:** `gonioApprox.h`, Free Functions

- `grids2D.h.skeleton.md` (0 KB)
  - **Summary:** * `subLoopTile` * `subdivideLoopGrid`
  - **Headings:** `grids2D.h`, Free Functions

- `integerOps.h.skeleton.md` (0 KB)
  - **Summary:** * `clamp_index_fast` * `fhash_Wang`
  - **Headings:** `integerOps.h`, Free Functions

- `integration.h.skeleton.md` (0 KB)
  - **Summary:** * `AdaptiveIntegral::simpson` * `AdaptiveIntegral::simpsonRecur`
  - **Headings:** `integration.h`, Free Functions

- `IntersectionCurve.h.skeleton.md` (0 KB)
  - **Summary:** * `Cs` * `Econv`
  - **Headings:** `IntersectionCurve.h`, Classes,     - `IntersectionCurve`

- `lineSearch.h.skeleton.md` (0 KB)
  - **Summary:** * `bound` * `side`
  - **Headings:** `lineSearch.h`, Classes,     - `LineSearch`, Free Functions

- `Lingebra.cpp.skeleton.md` (1 KB)
  - **Summary:** * `GaussElimination` * `Lingebra::GaussElimination`
  - **Headings:** `Lingebra.cpp`, Free Functions

- `Lingebra.h.skeleton.md` (0 KB)
  - **Summary:** * `Ap` * `Ap`
  - **Headings:** `Lingebra.h`, Classes,     - `LinSolver`, Free Functions

- `Mat3.h.skeleton.md` (2 KB)
  - **Summary:** * `a` * `array`
  - **Headings:** `Mat3.h`, Classes,     - `Mat3T`, Free Functions

- `Mat4.h.skeleton.md` (1 KB)
  - **Summary:** * `a` * `arr2d`
  - **Headings:** `Mat4.h`, Classes,     - `Mat4T`

- `Mesh.cpp.skeleton.md` (0 KB)
  - **Summary:** (no prose lines)
  - **Headings:** `Mesh.cpp`

- `Mesh.h.skeleton.md` (1 KB)
  - **Summary:** * `Rbound` * `center`
  - **Headings:** `Mesh.h`, Classes,     - `Mesh`, Free Functions

- `MultipoleGrid.h.skeleton.md` (0 KB)
  - **Summary:** * `Q` * `pos`
  - **Headings:** `MultipoleGrid.h`, Classes,     - `Atom`,     - `MultipoleGrid`

- `Multipoles.h.skeleton.md` (0 KB)
  - **Summary:** * `Coulomb` * `EF_brute`
  - **Headings:** `Multipoles.h`, Free Functions

- `NN.h.skeleton.md` (0 KB)
  - **Summary:** * `layers` * `nLayers`
  - **Headings:** `NN.h`, Classes,     - `NN`,     - `NN_layer`

- `Noise.cpp.skeleton.md` (0 KB)
  - **Summary:** (no prose lines)
  - **Headings:** `Noise.cpp`

- `Noise.h.skeleton.md` (0 KB)
  - **Summary:** * `Noise::craterProfile` * `Noise::getCraterHeight`
  - **Headings:** `Noise.h`, Free Functions

- `NumT.h.skeleton.md` (2 KB)
  - **Summary:** * `from` * `to`
  - **Headings:** `NumT.h`, Classes,     - `NumT::Range`,     - `NumT::Scan`,     - `NumT::ScanN`,     - `NumT::Tensor`,     - `NumT::Tuple`,     - `NumT::VecN`
    ... (+11 more)

- `NURBS.h.skeleton.md` (0 KB)
  - **Summary:** * `NURBS::fe1d` * `NURBS::sample1D`
  - **Headings:** `NURBS.h`, Free Functions

- `Opt3d.h.skeleton.md` (0 KB)
  - **Summary:** * `damp` * `dampMax`
  - **Headings:** `Opt3d.h`, Classes,     - `Opt3d`

- `optimizer_random.h.skeleton.md` (0 KB)
  - **Summary:** * `bestFitness` * `dxMax`
  - **Headings:** `optimizer_random.h`, Classes,     - `OptimizerRandom`,     - `OptimizerRandom_2`,     - `OptimizerRandom_3`

- `OptRandomWalk.h.skeleton.md` (0 KB)
  - **Summary:** * `E` * `Ebest`
  - **Headings:** `OptRandomWalk.h`, Classes,     - `OptRandomWalk`

- `physics_constants.h.skeleton.md` (0 KB)
  - **Summary:** (no prose lines)
  - **Headings:** `physics_constants.h`

- `pixMapFractal.h.skeleton.md` (0 KB)
  - **Summary:** * `binaryPixMapFrac`
  - **Headings:** `pixMapFractal.h`, Free Functions

- `PolyLinear1d.h.skeleton.md` (0 KB)
  - **Summary:** * `n` * `xs`
  - **Headings:** `PolyLinear1d.h`, Classes,     - `PolyLinear1d`

- `PolyLinear2d.h.skeleton.md` (0 KB)
  - **Summary:** * `Lin2Dcross::val` * `val`
  - **Headings:** `PolyLinear2d.h`, Free Functions

- `PotentialFlow.h.skeleton.md` (0 KB)
  - **Summary:** * `CurrentLoop` * `ILineFinite`
  - **Headings:** `PotentialFlow.h`, Free Functions

- `QuadratureCoefsOptimizer.h.skeleton.md` (1 KB)
  - **Summary:** * `bRationalTrainingFunction` * `nparP`
  - **Headings:** `QuadratureCoefsOptimizer.h`, Classes,     - `QuadratureCoefsOptimizer`, Free Functions

- `quaternion.h.skeleton.md` (2 KB)
  - **Summary:** * `array` * `e`
  - **Headings:** `quaternion.h`, Classes,     - `Quat4T`,     - `Quat8T`, Free Functions

- `radial_splines.h.skeleton.md` (0 KB)
  - **Summary:** * `coefs` * `dm`
  - **Headings:** `radial_splines.h`, Classes,     - `GeneralSpline`

- `raytrace.h.skeleton.md` (0 KB)
  - **Summary:** * `capsulaIntersect` * `linePointDistance2`
  - **Headings:** `raytrace.h`, Free Functions

- `RBF2D.h.skeleton.md` (0 KB)
  - **Summary:** * `height` * `pos`
  - **Headings:** `RBF2D.h`, Classes,     - `RBF2D`

- `RotationMesh.h.skeleton.md` (0 KB)
  - **Summary:** * `mrots` * `n`
  - **Headings:** `RotationMesh.h`, Classes,     - `RotationMesh`

- `SchroedingerGreen1D.h.skeleton.md` (0 KB)
  - **Summary:** * `Apsi` * `E`
  - **Headings:** `SchroedingerGreen1D.h`, Classes,     - `SchroedingerGreen1D`

- `SchroedingerGreen2D.h.skeleton.md` (0 KB)
  - **Summary:** * `Apsi` * `E`
  - **Headings:** `SchroedingerGreen2D.h`, Classes,     - `SchroedingerGreen2D`

- `simd.h.skeleton.md` (0 KB)
  - **Summary:** * `hsum_double_avx` * `pack_Vec3d_simd`
  - **Headings:** `simd.h`, Free Functions

- `SMat3.h.skeleton.md` (0 KB)
  - **Summary:** * `array` * `diag`
  - **Headings:** `SMat3.h`, Classes,     - `SMat3`, Free Functions

- `Solids.h.skeleton.md` (0 KB)
  - **Summary:** (no prose lines)
  - **Headings:** `Solids.h`

- `spaceFillingCurves.h.skeleton.md` (0 KB)
  - **Summary:** * `HiblertCurve::d2xy` * `HiblertCurve::rot`
  - **Headings:** `spaceFillingCurves.h`, Free Functions

- `Spline2d.h.skeleton.md` (0 KB)
  - **Summary:** * `CPs` * `CPs`
  - **Headings:** `Spline2d.h`, Classes,     - `Spline2d`,     - `Spline2dSampler`,     - `SplinePoint2d`

- `spline_Circle.h.skeleton.md` (0 KB)
  - **Summary:** * `CPs` * `CPs`
  - **Headings:** `spline_Circle.h`, Classes,     - `CircleSpline`

- `spline_hermite.h.skeleton.md` (1 KB)
  - **Summary:** * `dx` * `dy0`
  - **Headings:** `spline_hermite.h`, Classes,     - `Sampler`,     - `Spline_Hermite::Sampler`, Free Functions

- `spline_triC1.h.skeleton.md` (0 KB)
  - **Summary:** * `Spline_triC1::pi` * `Spline_triC1::pitau`
  - **Headings:** `spline_triC1.h`, Free Functions

- `SplineManager.cpp.skeleton.md` (0 KB)
  - **Summary:** (no prose lines)
  - **Headings:** `SplineManager.cpp`

- `SplineManager.h.skeleton.md` (0 KB)
  - **Summary:** * `CPs` * `CPs`
  - **Headings:** `SplineManager.h`, Classes,     - `SplineManager`

- `TestFuncND.h.skeleton.md` (0 KB)
  - **Summary:** * `coefs` * `freq0`
  - **Headings:** `TestFuncND.h`, Classes,     - `TestFuncND`

- `tresholdFunctions.h.skeleton.md` (0 KB)
  - **Summary:** * `Treshold::p0` * `Treshold::p0i`
  - **Headings:** `tresholdFunctions.h`, Free Functions

- `UVfuncs.h.skeleton.md` (0 KB)
  - **Summary:** * `ConeUVfunc` * `HarmonicTubeUVfunc`
  - **Headings:** `UVfuncs.h`, Free Functions

- `Vec2.h.skeleton.md` (2 KB)
  - **Summary:** * `a` * `array`
  - **Headings:** `Vec2.h`, Classes,     - `Vec2T`, Free Functions

- `Vec3.h.skeleton.md` (3 KB)
  - **Summary:** * `a` * `array`
  - **Headings:** `Vec3.h`, Classes,     - `Vec3T`,     - `Vec6T`, Free Functions

- `Vec3Utils.h.skeleton.md` (0 KB)
  - **Summary:** * `affineTransform` * `average`
  - **Headings:** `Vec3Utils.h`, Free Functions

- `VecN.h.skeleton.md` (0 KB)
  - **Summary:** * `VecN::absmax` * `VecN::add`
  - **Headings:** `VecN.h`, Free Functions

- `VecNm.h.skeleton.md` (0 KB)
  - **Summary:** * `data` * `m`
  - **Headings:** `VecNm.h`, Classes,     - `VecNm`

- `voronoi.cpp.skeleton.md` (0 KB)
  - **Summary:** * `Vec2dCompare`
  - **Headings:** `voronoi.cpp`, Free Functions

- `Voronoi.cpp.skeleton.md` (0 KB)
  - **Summary:** (no prose lines)
  - **Headings:** `Voronoi.cpp`

- `voronoi.h.skeleton.md` (2 KB)
  - **Summary:** * `Sites` * `Vertices`
  - **Headings:** `voronoi.h`, Classes,     - `Edge`,     - `FreeNodeArrayList`,     - `Freelist`,     - `Freenode`,     - `GraphEdge`,     - `Halfedge`
    ... (+3 more)

- `Voronoi.h.skeleton.md` (1 KB)
  - **Summary:** * `direction` * `end`
  - **Headings:** `Voronoi.h`, Classes,     - `CompareEvent`,     - `VEdge`,     - `VEvent`,     - `VParabola`,     - `Voronoi`,     - `VoronoiNamespace::VEdge`
    ... (+4 more)

- `warpFunction2D.h.skeleton.md` (0 KB)
  - **Summary:** * `warp_sin` * `warp_x1period`
  - **Headings:** `warpFunction2D.h`, Free Functions


## `doc/skeletons/common/molecular/`

- `AOIntegrals.h.skeleton.md` (0 KB)
  - **Summary:** * `dot_rolled` * `dot_shifted`
  - **Headings:** `AOIntegrals.h`, Free Functions

- `AOrotations.h.skeleton.md` (0 KB)
  - **Summary:** * `Iss` * `Isz`
  - **Headings:** `AOrotations.h`, Classes,     - `AOPairType`

- `AtomicConfiguration.h.skeleton.md` (0 KB)
  - **Summary:** * `natoms` * `pos`
  - **Headings:** `AtomicConfiguration.h`, Classes,     - `AtomicConfiguration`,     - `FastAtomicMetric`

- `AtomicSystemDescriptors.h.skeleton.md` (0 KB)
  - **Summary:** * `Rmax` * `Rmin`
  - **Headings:** `AtomicSystemDescriptors.h`, Classes,     - `AtomicSystemDescriptors`, Free Functions

- `Atoms.h.skeleton.md` (0 KB)
  - **Summary:** * `Energy` * `apos`
  - **Headings:** `Atoms.h`, Classes,     - `Atoms`,     - `MetaData`,     - `Points`

- `AtomsInGrid.h.skeleton.md` (0 KB)
  - **Summary:** * `neighs` * `type`
  - **Headings:** `AtomsInGrid.h`, Classes,     - `AtomInBox`,     - `AtomsInGrid`,     - `GridPointDynamics`

- `AuxFF.h.skeleton.md` (0 KB)
  - **Summary:** * `REQs` * `bond`
  - **Headings:** `AuxFF.h`, Classes,     - `AuxFF`

- `BehlerSymFunc.h.skeleton.md` (0 KB)
  - **Summary:** * `Rcut` * `SFA`
  - **Headings:** `BehlerSymFunc.h`, Classes,     - `BehlerNNFF`, Free Functions

- `BOFF.h.skeleton.md` (0 KB)
  - **Summary:** * `dval` * `i`
  - **Headings:** `BOFF.h`, Classes,     - `BOFF`,     - `BOFF::BObond`,     - `BObond`

- `BondAdaptedMesh.h.skeleton.md` (0 KB)
  - **Summary:** * `mesh` * `tetrahedrons`
  - **Headings:** `BondAdaptedMesh.h`, Classes,     - `AtomGrid`

- `CarParrinello.h.skeleton.md` (0 KB)
  - **Summary:** (no prose lines)
  - **Headings:** `CarParrinello.h`, Classes,     - `CarParrinello`

- `CLCFGO.h.skeleton.md` (4 KB)
  - **Summary:** * `Eaa` * `Eae`
  - **Headings:** `CLCFGO.h`, Classes,     - `CLCFGO`

- `CLCFGO_Structs.h.skeleton.md` (3 KB)
  - **Summary:** * `Eaa` * `Eae`
  - **Headings:** `CLCFGO_Structs.h`, Classes,     - `CLCFGO`,     - `Nucleus`,     - `Orbital`

- `CLCFGO_tests.h.skeleton.md` (0 KB)
  - **Summary:** * `getNumDeriv` * `initTestElectrons`
  - **Headings:** `CLCFGO_tests.h`, Free Functions

- `CLCFSF.h.skeleton.md` (1 KB)
  - **Summary:** * `ICoulomb` * `IKinetic`
  - **Headings:** `CLCFSF.h`, Classes,     - `CLCFSF`,     - `PairInteractionTable`

- `Confs.h.skeleton.md` (0 KB)
  - **Summary:** * `cs` * `iatoms`
  - **Headings:** `Confs.h`, Classes,     - `CGroup`,     - `Confs`

- `constrains.h.skeleton.md` (1 KB)
  - **Summary:** * `acell` * `active`
  - **Headings:** `constrains.h`, Classes,     - `AngleConstr`,     - `Constrains`,     - `DistConstr`,     - `ForcedDriver`,     - `SplineConstr`,     - `TorsionConstr`

- `DipoleMap.h.skeleton.md` (0 KB)
  - **Summary:** * `FE` * `FE2`
  - **Headings:** `DipoleMap.h`, Classes,     - `DipoleMap`

- `DirectionStiffness.h.skeleton.md` (1 KB)
  - **Summary:** * `Ks` * `L`
  - **Headings:** `DirectionStiffness.h`, Classes,     - `Deformer`,     - `DirStiff`

- `e2FF.h.skeleton.md` (0 KB)
  - **Summary:** * `Eaa` * `Eae`
  - **Headings:** `e2FF.h`, Classes,     - `E2FF`,     - `EFF2AtomType`, Free Functions

- `eFF.h.skeleton.md` (2 KB)
  - **Summary:** * `Eaa` * `Eae`
  - **Headings:** `eFF.h`, Classes,     - `EFF`,     - `EFFAtomType`, Free Functions

- `eFF_bak.h.skeleton.md` (1 KB)
  - **Summary:** * `Eaa` * `Eae`
  - **Headings:** `eFF_bak.h`, Classes,     - `EFF`,     - `EFFAtomType`, Free Functions

- `eFF_old.h.skeleton.md` (1 KB)
  - **Summary:** * `Eaa` * `Eae`
  - **Headings:** `eFF_old.h`, Classes,     - `EFF`,     - `EFFAtomType`, Free Functions

- `eFF_old_aPar.h.skeleton.md` (1 KB)
  - **Summary:** * `Eaa` * `Eae`
  - **Headings:** `eFF_old_aPar.h`, Classes,     - `EFF`,     - `EFFAtomType`, Free Functions

- `eFF_plots.h.skeleton.md` (0 KB)
  - **Summary:** * `checkDerivs` * `checkFinite`
  - **Headings:** `eFF_plots.h`, Free Functions

- `elFF.h.skeleton.md` (0 KB)
  - **Summary:** * `a2w` * `acsums`
  - **Headings:** `elFF.h`, Classes,     - `ElFF`

- `EOFF.h.skeleton.md` (0 KB)
  - **Summary:** * `E0` * `E0`
  - **Headings:** `EOFF.h`, Classes,     - `EAFFelement`,     - `EOFF`,     - `EOFFAtom`, Free Functions

- `ESFF.h.skeleton.md` (0 KB)
  - **Summary:** * `aZs` * `aee`
  - **Headings:** `ESFF.h`, Classes,     - `ESFF`,     - `ESFFelement`

- `EwaldGrid.h.skeleton.md` (1 KB)
  - **Summary:** * `Qabs` * `Qabs_g`
  - **Headings:** `EwaldGrid.h`, Classes,     - `EwaldGrid`, Free Functions

- `FF2D.h.skeleton.md` (1 KB)
  - **Summary:** * `nb` * `ne`
  - **Headings:** `FF2D.h`, Classes,     - `Atom2D`,     - `FF2D`

- `FireCoreAPI.h.skeleton.md` (1 KB)
  - **Summary:** * `evalForce` * `getCharges`
  - **Headings:** `FireCoreAPI.h`, Classes,     - `FireCore::Lib`,     - `FireCore::QMMM`,     - `Lib`,     - `QMMM`

- `FitFF.h.skeleton.md` (0 KB)
  - **Summary:** * `ij` * `paramOffset`
  - **Headings:** `FitFF.h`, Classes,     - `AtomPair`,     - `FitStructure`, Free Functions

- `FitREQ.h.skeleton.md` (2 KB)
  - **Summary:** * `DOFs` * `Eilines`
  - **Headings:** `FitREQ.h`, Classes,     - `FitREQ`, Free Functions

- `FlexibleAtomReactiveFF.h.skeleton.md` (1 KB)
  - **Summary:** * `Eatoms` * `Epairs`
  - **Headings:** `FlexibleAtomReactiveFF.h`, Classes,     - `FARFF`,     - `FlexibleAtomType`,     - `FlexiblePairType`

- `FlexibleAtomReactiveFF_dyn.h.skeleton.md` (1 KB)
  - **Summary:** * `conf` * `energy`
  - **Headings:** `FlexibleAtomReactiveFF_dyn.h`, Classes,     - `FARFFAtom`,     - `FARFF_dyn`

- `ForceField.h.skeleton.md` (1 KB)
  - **Summary:** * `ang` * `bAng`
  - **Headings:** `ForceField.h`, Classes,     - `CollisionDamping`,     - `ForceField`, Free Functions

- `FspFF.h.skeleton.md` (1 KB)
  - **Summary:** * `Kea` * `Kee`
  - **Headings:** `FspFF.h`, Classes,     - `FspFF`,     - `Mat3Sd`, Free Functions

- `FspFFclean.h.skeleton.md` (1 KB)
  - **Summary:** * `Kae` * `Kee`
  - **Headings:** `FspFFclean.h`, Classes,     - `FspFF`

- `FTRFF.h.skeleton.md` (0 KB)
  - **Summary:** * `fbond` * `force`
  - **Headings:** `FTRFF.h`, Classes,     - `FTAtom`,     - `FTRFF`,     - `FTtype`, Free Functions

- `GlobalOptimizer.h.skeleton.md` (1 KB)
  - **Summary:** * `DoF` * `F_opt`
  - **Headings:** `GlobalOptimizer.h`, Classes,     - `GlobalOptimizer`

- `GOpt.h.skeleton.md` (0 KB)
  - **Summary:** * `T_current` * `T_target`
  - **Headings:** `GOpt.h`, Classes,     - `GOpt`

- `GridFF.h.skeleton.md` (3 KB)
  - **Summary:** * `Bspline_Coulomb` * `Bspline_London`
  - **Headings:** `GridFF.h`, Classes,     - `GridFF`, Free Functions

- `Groups.h.skeleton.md` (0 KB)
  - **Summary:** * `cog` * `fDrive`
  - **Headings:** `Groups.h`, Classes,     - `Group`,     - `Groups`

- `InteractionsGauss.h.skeleton.md` (0 KB)
  - **Summary:** * `CoulombGauss` * `CoulombGauss_FixSize`
  - **Headings:** `InteractionsGauss.h`, Free Functions

- `InteractionsGauss_old.h.skeleton.md` (0 KB)
  - **Summary:** * `CoulombGauss` * `CoulombGauss_FixSize`
  - **Headings:** `InteractionsGauss_old.h`, Free Functions

- `InterpolateTricubic.h.skeleton.md` (1 KB)
  - **Summary:** * `Spline_Hermite::basis` * `Spline_Hermite::basis_val`
  - **Headings:** `InterpolateTricubic.h`, Free Functions

- `InterpolateTrilinear.h.skeleton.md` (0 KB)
  - **Summary:** * `Trilinear::basis` * `Trilinear::fe3d_3`
  - **Headings:** `InterpolateTrilinear.h`, Free Functions

- `Kekule.h.skeleton.md` (1 KB)
  - **Summary:** * `Ea` * `Eb`
  - **Headings:** `Kekule.h`, Classes,     - `Kekule`

- `LatticeMatch2D.h.skeleton.md` (0 KB)
  - **Summary:** * `i` * `val`
  - **Headings:** `LatticeMatch2D.h`, Classes,     - `IVal`,     - `Latmiss`,     - `LatticeMatch2D`, Free Functions

- `LCAOFF.h.skeleton.md` (0 KB)
  - **Summary:** * `atoms` * `atoms`
  - **Headings:** `LCAOFF.h`, Classes,     - `LCAOFF`,     - `TBparams`, Free Functions

- `libMMFF.h.skeleton.md` (1 KB)
  - **Summary:** * `addAngConstrain` * `addDistConstrain`
  - **Headings:** `libMMFF.h`, Free Functions

- `MMFF.h.skeleton.md` (3 KB)
  - **Summary:** * `Fmax` * `Rconv`
  - **Headings:** `MMFF.h`, Classes,     - `AtomicManipulator`,     - `MMFF`,     - `MolType`, Free Functions

- `MMFFBuilder.h.skeleton.md` (11 KB)
  - **Summary:** * `C0` * `C1`
  - **Headings:** `MMFFBuilder.h`, Classes,     - `Angle`,     - `Atom`,     - `AtomConf`,     - `Bond`,     - `Builder`,     - `Dihedral`
    ... (+11 more)

- `MMFFf4.h.skeleton.md` (1 KB)
  - **Summary:** * `DOFs` * `Ea`
  - **Headings:** `MMFFf4.h`, Classes,     - `MMFFf4`, Free Functions

- `MMFFmini.h.skeleton.md` (1 KB)
  - **Summary:** * `Ea` * `Eb`
  - **Headings:** `MMFFmini.h`, Classes,     - `MMFFmini`, Free Functions

- `MMFFparams.h.skeleton.md` (3 KB)
  - **Summary:** * `angle0` * `atoms`
  - **Headings:** `MMFFparams.h`, Classes,     - `AngleType`,     - `AtomType`,     - `BondType`,     - `DihedralType`,     - `ElementType`,     - `MMFFparams`

- `MMFFsp3.h.skeleton.md` (2 KB)
  - **Summary:** * `DOFs` * `Ea`
  - **Headings:** `MMFFsp3.h`, Classes,     - `MMFFsp3`

- `MMFFsp3_loc.h.skeleton.md` (2 KB)
  - **Summary:** * `DOFs` * `Ea`
  - **Headings:** `MMFFsp3_loc.h`, Classes,     - `MMFFsp3_loc`

- `molecular.cpp.skeleton.md` (0 KB)
  - **Summary:** * `check_disturb` * `disturb`
  - **Headings:** `molecular.cpp`, Free Functions

- `molecular.h.skeleton.md` (0 KB)
  - **Summary:** * `EiiLJ` * `RiiLJ`
  - **Headings:** `molecular.h`, Classes,     - `AtomType`

- `molecular_utils.h.skeleton.md` (0 KB)
  - **Summary:** * `checkNumRange` * `ckeckNaN`
  - **Headings:** `molecular_utils.h`, Free Functions

- `MolecularConfiguration.h.skeleton.md` (0 KB)
  - **Summary:** * `poses` * `poses`
  - **Headings:** `MolecularConfiguration.h`, Classes,     - `MolecularConfiguration`,     - `MolecularConfigurationMemory`, Free Functions

- `MolecularDatabase.h.skeleton.md` (1 KB)
  - **Summary:** * `descs` * `descs`
  - **Headings:** `MolecularDatabase.h`, Classes,     - `Descriptor`,     - `MolecularDatabase`,     - `MolecularDatabaseDescriptor`

- `MolecularGraph.h.skeleton.md` (0 KB)
  - **Summary:** * `amask` * `atom2bond`
  - **Headings:** `MolecularGraph.h`, Classes,     - `MolecularGraph`

- `Molecule.h.skeleton.md` (1 KB)
  - **Summary:** * `REQs` * `ang0s`
  - **Headings:** `Molecule.h`, Classes,     - `Molecule`, Free Functions

- `MoleculeType.cpp.skeleton.md` (0 KB)
  - **Summary:** (no prose lines)
  - **Headings:** `MoleculeType.cpp`

- `MoleculeWorld2D.cpp.skeleton.md` (0 KB)
  - **Summary:** * `forceTransform` * `interMolForce`
  - **Headings:** `MoleculeWorld2D.cpp`, Free Functions

- `MoleculeWorld2D.h.skeleton.md` (0 KB)
  - **Summary:** * `Qs` * `Rs`
  - **Headings:** `MoleculeWorld2D.h`, Classes,     - `MoleculeType2D`,     - `MoleculeWorld`, Free Functions

- `MolWorld_sp3.h.skeleton.md` (7 KB)
  - **Summary:** * `ConstrZ_k` * `ConstrZ_l`
  - **Headings:** `MolWorld_sp3.h`, Classes,     - `MolWorld_sp3`

- `MolWorld_sp3_multi.h.skeleton.md` (5 KB)
  - **Summary:** * `cf` * `cos_vf`
  - **Headings:** `MolWorld_sp3_multi.h`, Classes,     - `FIRE`,     - `FIRE_setup`,     - `MolWorld_sp3_multi`

- `MolWorld_sp3_QMMM.h.skeleton.md` (0 KB)
  - **Summary:** * `bPrepared_mm` * `bPrepared_qm`
  - **Headings:** `MolWorld_sp3_QMMM.h`, Classes,     - `MolWorld_sp3_QMMM`

- `MolWorld_sp3_QMMM_multi.h.skeleton.md` (0 KB)
  - **Summary:** * `bPrepared_mm` * `bPrepared_qm`
  - **Headings:** `MolWorld_sp3_QMMM_multi.h`, Classes,     - `MolWorld_sp3_QMMM_multi`

- `MolWorld_sp3_simple.h.skeleton.md` (2 KB)
  - **Summary:** * `Etot` * `Kmorse`
  - **Headings:** `MolWorld_sp3_simple.h`, Classes,     - `MolWorld_sp3_simple`

- `MultiSolverInterface.h.skeleton.md` (0 KB)
  - **Summary:** (no prose lines)
  - **Headings:** `MultiSolverInterface.h`, Classes,     - `MultiSolverInterface`,     - `SolverInterface`

- `NBFF.h.skeleton.md` (2 KB)
  - **Summary:** * `BBs` * `BBs`
  - **Headings:** `NBFF.h`, Classes,     - `NBFF`, Free Functions

- `NBFF_old.h.skeleton.md` (0 KB)
  - **Summary:** * `neighList` * `nmask`
  - **Headings:** `NBFF_old.h`, Classes,     - `NBFF_old`, Free Functions

- `NBFF_SR.h.skeleton.md` (0 KB)
  - **Summary:** * `a2g` * `bInGroupNB`
  - **Headings:** `NBFF_SR.h`, Classes,     - `NBFF_SR`, Free Functions

- `PBCsystem.h.skeleton.md` (0 KB)
  - **Summary:** * `R2cut` * `R2min`
  - **Headings:** `PBCsystem.h`, Classes,     - `AnyShortRangeFF`,     - `PBCsystem`, Free Functions

- `PiPlanesFF.h.skeleton.md` (0 KB)
  - **Summary:** * `K` * `K`
  - **Headings:** `PiPlanesFF.h`, Classes,     - `PiPlane`,     - `PiPlanes`

- `QEq.h.skeleton.md` (0 KB)
  - **Summary:** * `J` * `J`
  - **Headings:** `QEq.h`, Classes,     - `QEq`, Free Functions

- `RARFF.h.skeleton.md` (0 KB)
  - **Summary:** * `atoms` * `invRotMass`
  - **Headings:** `RARFF.h`, Classes,     - `RARFF`,     - `RigidAtom`,     - `RigidAtomType`, Free Functions

- `RARFF2.h.skeleton.md` (1 KB)
  - **Summary:** * `atoms` * `fbonds`
  - **Headings:** `RARFF2.h`, Classes,     - `RARFF2`,     - `RigidAtom`,     - `RigidAtomType`, Free Functions

- `RARFF_SR.h.skeleton.md` (2 KB)
  - **Summary:** * `EvdW` * `RvdW`
  - **Headings:** `RARFF_SR.h`, Classes,     - `CapType`,     - `PairList`,     - `RARFF_SR`,     - `RigidAtomType`, Free Functions

- `RARFFarr.h.skeleton.md` (1 KB)
  - **Summary:** * `EvdW` * `RvdW`
  - **Headings:** `RARFFarr.h`, Classes,     - `CapType`,     - `RARFF2arr`,     - `RigidAtomType`, Free Functions

- `RBodyConfDyn.cpp.skeleton.md` (0 KB)
  - **Summary:** (no prose lines)
  - **Headings:** `RBodyConfDyn.cpp`

- `RBodyConfDyn.h.skeleton.md` (1 KB)
  - **Summary:** * `E_last` * `FposScale`
  - **Headings:** `RBodyConfDyn.h`, Classes,     - `RBodyConfDyn`,     - `RBodyPose`

- `RigidBodyFF.h.skeleton.md` (0 KB)
  - **Summary:** * `apos0` * `frots`
  - **Headings:** `RigidBodyFF.h`, Classes,     - `RigidBodyFF`

- `RRFF.h.skeleton.md` (0 KB)
  - **Summary:** * `expar` * `hbond`
  - **Headings:** `RRFF.h`, Classes,     - `NeighAtom`,     - `RRFF`,     - `RRFFelement`,     - `RRFFparams`,     - `RRFFtype`, Free Functions

- `RspFF.h.skeleton.md` (0 KB)
  - **Summary:** * `REQs` * `b2a`
  - **Headings:** `RspFF.h`, Classes,     - `RARFF2arr`, Free Functions

- `SimpleForceField.h.skeleton.md` (0 KB)
  - **Summary:** * `SimpleForceField::forceLJE` * `SimpleForceField::init`
  - **Headings:** `SimpleForceField.h`, Free Functions

- `StructuralHistogram.h.skeleton.md` (0 KB)
  - **Summary:** * `bondHistBounds` * `bondHistType`
  - **Headings:** `StructuralHistogram.h`, Classes,     - `StructuralHistogram_h`, Free Functions

- `SVG_render.h.skeleton.md` (0 KB)
  - **Summary:** * `bFill` * `bStroke`
  - **Headings:** `SVG_render.h`, Classes,     - `SVG_render`

- `UFF.h.skeleton.md` (2 KB)
  - **Summary:** * `Ea` * `Eb`
  - **Headings:** `UFF.h`, Classes,     - `UFF`, Free Functions


## `doc/skeletons/common_SDL/`

- `browser_sdl.h.skeleton.md` (1 KB)
  - **Summary:** * `bMyDownScale` * `bd`
  - **Headings:** `browser_sdl.h`, Classes,     - `BrowserSDL`

- `GLUtils.h.skeleton.md` (0 KB)
  - **Summary:** * `countVerts` * `hardFace`
  - **Headings:** `GLUtils.h`, Free Functions

- `SDL_utils.h.skeleton.md` (0 KB)
  - **Summary:** * `loadDataImageFloat` * `makeTexture`
  - **Headings:** `SDL_utils.h`, Free Functions


## `doc/skeletons/common_SDL/SDL2/`

- `SDLplot.h.skeleton.md` (0 KB)
  - **Summary:** * `DestR` * `IX0`
  - **Headings:** `SDLplot.h`, Classes,     - `SDLplot`


## `doc/skeletons/common_SDL/SDL2OGL/`

- `AppSDL2OGL.cpp.skeleton.md` (0 KB)
  - **Summary:** (no prose lines)
  - **Headings:** `AppSDL2OGL.cpp`

- `AppSDL2OGL.h.skeleton.md` (0 KB)
  - **Summary:** * `LMB` * `RMB`
  - **Headings:** `AppSDL2OGL.h`, Classes,     - `AppSDL2OGL`

- `AppSDL2OGL_3D.cpp.skeleton.md` (0 KB)
  - **Summary:** (no prose lines)
  - **Headings:** `AppSDL2OGL_3D.cpp`

- `AppSDL2OGL_3D.h.skeleton.md` (0 KB)
  - **Summary:** * `bDragging` * `cam`
  - **Headings:** `AppSDL2OGL_3D.h`, Classes,     - `AppSDL2OGL_3D`

- `cameraOGL.h.skeleton.md` (0 KB)
  - **Summary:** * `Cam::ortho` * `Cam::perspective`
  - **Headings:** `cameraOGL.h`, Free Functions

- `Console.h.skeleton.md` (0 KB)
  - **Summary:** * `bShift` * `callback`
  - **Headings:** `Console.h`, Classes,     - `Console`

- `Draw.cpp.skeleton.md` (0 KB)
  - **Summary:** (no prose lines)
  - **Headings:** `Draw.cpp`

- `Draw.h.skeleton.md` (0 KB)
  - **Summary:** * `Draw::color_of_hash` * `Draw::float2RGBA`
  - **Headings:** `Draw.h`, Free Functions

- `Draw2D.cpp.skeleton.md` (0 KB)
  - **Summary:** (no prose lines)
  - **Headings:** `Draw2D.cpp`

- `Draw2D.h.skeleton.md` (0 KB)
  - **Summary:** * `Draw2D::toGLMat` * `Draw2D::toGLMatCam`
  - **Headings:** `Draw2D.h`, Free Functions

- `Draw3D.cpp.skeleton.md` (2 KB)
  - **Summary:** * `Draw3D::drawArrow` * `Draw3D::drawAxis`
  - **Headings:** `Draw3D.cpp`, Free Functions

- `Draw3D.h.skeleton.md` (1 KB)
  - **Summary:** * `Draw3D::color` * `Draw3D::drawArrow`
  - **Headings:** `Draw3D.h`, Free Functions

- `Draw3D_Molecular.h.skeleton.md` (0 KB)
  - **Summary:** * `Draw3D::angle` * `Draw3D::atomLabels`
  - **Headings:** `Draw3D_Molecular.h`, Free Functions

- `Draw3D_Surf.h.skeleton.md` (0 KB)
  - **Summary:** * `Draw3D::drawExtrudedWireUVFunc` * `Draw3D::drawSmoothUVFunc`
  - **Headings:** `Draw3D_Surf.h`, Free Functions

- `DrawIso.h.skeleton.md` (0 KB)
  - **Summary:** * `Draw3D::drawIso` * `drawIso`
  - **Headings:** `DrawIso.h`, Free Functions

- `DrawSphereMap.h.skeleton.md` (0 KB)
  - **Summary:** * `SphereSampling::drawDiTri` * `SphereSampling::drawDiTriWire`
  - **Headings:** `DrawSphereMap.h`, Free Functions

- `EditorGizmo.h.skeleton.md` (1 KB)
  - **Summary:** * `Rgizmo` * `Rhandle`
  - **Headings:** `EditorGizmo.h`, Classes,     - `EditorGizmo`,     - `Pose3d`, Free Functions

- `GUI.cpp.skeleton.md` (0 KB)
  - **Summary:** (no prose lines)
  - **Headings:** `GUI.cpp`

- `GUI.h.skeleton.md` (3 KB)
  - **Summary:** * `label` * `master`
  - **Headings:** `GUI.h`, Classes,     - `CheckBox`,     - `CheckBoxList`,     - `Command`,     - `CommandList`,     - `Commander`,     - `DropDownList`
    ... (+11 more)

- `icosahedron.h.skeleton.md` (0 KB)
  - **Summary:** * `drawIcosaSphere` * `drawtriangle`
  - **Headings:** `icosahedron.h`, Free Functions

- `icosahedron_2.h.skeleton.md` (0 KB)
  - **Summary:** * `drawIcosa` * `drawtriangle`
  - **Headings:** `icosahedron_2.h`, Free Functions

- `loadOBJ_SDL2.h.skeleton.md` (0 KB)
  - **Summary:** * `loadObjToList`
  - **Headings:** `loadOBJ_SDL2.h`, Free Functions

- `MarchingCubes.h.skeleton.md` (0 KB)
  - **Summary:** * `Draw3D::Grid2Points` * `Draw3D::LerpQ`
  - **Headings:** `MarchingCubes.h`, Free Functions

- `MolecularDraw.h.skeleton.md` (0 KB)
  - **Summary:** * `colorBW` * `colorRB`
  - **Headings:** `MolecularDraw.h`, Free Functions

- `MolGUI.h.skeleton.md` (5 KB)
  - **Summary:** * `BondLengh_max` * `BondLengh_min`
  - **Headings:** `MolGUI.h`, Classes,     - `MolGUI`, Free Functions

- `MolGUI_tests.h.skeleton.md` (0 KB)
  - **Summary:** * `torsion_Paolo` * `torsion_Prokop`
  - **Headings:** `MolGUI_tests.h`, Free Functions

- `Plot2D.cpp.skeleton.md` (0 KB)
  - **Summary:** (no prose lines)
  - **Headings:** `Plot2D.cpp`

- `Plot2D.h.skeleton.md` (1 KB)
  - **Summary:** * `bSharedX` * `bView`
  - **Headings:** `Plot2D.h`, Classes,     - `DataLine2D`,     - `Plot2D`,     - `QuePlot2D`, Free Functions

- `PlotScreen2D.cpp.skeleton.md` (0 KB)
  - **Summary:** (no prose lines)
  - **Headings:** `PlotScreen2D.cpp`

- `PlotScreen2D.h.skeleton.md` (0 KB)
  - **Summary:** * `defaultExtend` * `defaultExtend`
  - **Headings:** `PlotScreen2D.h`, Classes,     - `PlotScreen2D`

- `ScreenSDL2OGL.cpp.skeleton.md` (0 KB)
  - **Summary:** * `setLightingNormal` * `setLightingRGB`
  - **Headings:** `ScreenSDL2OGL.cpp`, Free Functions

- `ScreenSDL2OGL.h.skeleton.md` (0 KB)
  - **Summary:** * `ASPECT_RATIO` * `GL_LOCK`
  - **Headings:** `ScreenSDL2OGL.h`, Classes,     - `ScreenSDL2OGL`

- `ScreenSDL2OGL_3D.cpp.skeleton.md` (0 KB)
  - **Summary:** (no prose lines)
  - **Headings:** `ScreenSDL2OGL_3D.cpp`

- `ScreenSDL2OGL_3D.h.skeleton.md` (0 KB)
  - **Summary:** * `cam` * `camDist`
  - **Headings:** `ScreenSDL2OGL_3D.h`, Classes,     - `ScreenSDL2OGL_3D`

- `TiledView.cpp.skeleton.md` (0 KB)
  - **Summary:** (no prose lines)
  - **Headings:** `TiledView.cpp`

- `TiledView.h.skeleton.md` (0 KB)
  - **Summary:** * `tiles` * `tiles`
  - **Headings:** `TiledView.h`, Classes,     - `TiledView`


## `doc/skeletons/common_SDL/SDL2OGL/temp/`

- `drawMath.cpp.skeleton.md` (0 KB)
  - **Summary:** * `drawAxis` * `drawBox`
  - **Headings:** `drawMath.cpp`, Free Functions

- `drawMath.h.skeleton.md` (0 KB)
  - **Summary:** * `setColorInt32` * `toGLMat`
  - **Headings:** `drawMath.h`, Free Functions

- `drawUtils.h.skeleton.md` (0 KB)
  - **Summary:** * `drawAxis` * `drawBox`
  - **Headings:** `drawUtils.h`, Free Functions

- `ScreenSDL2OGL_2D.cpp.skeleton.md` (0 KB)
  - **Summary:** (no prose lines)
  - **Headings:** `ScreenSDL2OGL_2D.cpp`

- `ScreenSDL2OGL_2D.h.skeleton.md` (0 KB)
  - **Summary:** * `ASPECT_RATIO` * `HEIGHT`
  - **Headings:** `ScreenSDL2OGL_2D.h`, Classes,     - `ScreenSDL2OGL_2D`


## `doc/skeletons/libs/`

- `quadrature_lib.cpp.skeleton.md` (0 KB)
  - **Summary:** * `distributPointsTetrahedron` * `evaluateQuadrature`
  - **Headings:** `quadrature_lib.cpp`, Free Functions


## `doc/skeletons/libs/Molecular/`

- `CLCFGO_lib.cpp.skeleton.md` (0 KB)
  - **Summary:** * `atomsPotAtPoints` * `coulombOrbPair`
  - **Headings:** `CLCFGO_lib.cpp`, Free Functions

- `eFF_lib.cpp.skeleton.md` (0 KB)
  - **Summary:** * `eval` * `evalFuncDerivs`
  - **Headings:** `eFF_lib.cpp`, Free Functions

- `FF2D_lib.cpp.skeleton.md` (0 KB)
  - **Summary:** * `addAtom` * `addBond`
  - **Headings:** `FF2D_lib.cpp`, Free Functions

- `FitREQ_lib.cpp.skeleton.md` (0 KB)
  - **Summary:** * `getEs` * `getType`
  - **Headings:** `FitREQ_lib.cpp`, Free Functions

- `Kekule_lib.cpp.skeleton.md` (0 KB)
  - **Summary:** * `eval` * `init`
  - **Headings:** `Kekule_lib.cpp`, Free Functions

- `Lattice2D_lib.cpp.skeleton.md` (0 KB)
  - **Summary:** * `getMatches` * `getVecMatch`
  - **Headings:** `Lattice2D_lib.cpp`, Free Functions

- `MMFF_lib.cpp.skeleton.md` (0 KB)
  - **Summary:** * `EwaldGridSolveLaplace` * `EwaldGridSolveLaplaceDebug`
  - **Headings:** `MMFF_lib.cpp`, Free Functions

- `MMFFsp3_lib.cpp.skeleton.md` (0 KB)
  - **Summary:** * `init` * `init_buffers`
  - **Headings:** `MMFFsp3_lib.cpp`, Free Functions

- `Molecular.cpp.skeleton.md` (0 KB)
  - **Summary:** * `addFitStruct` * `atoms2map`
  - **Headings:** `Molecular.cpp`, Free Functions

- `ReactiveFF.cpp.skeleton.md` (0 KB)
  - **Summary:** * `getBondCaps` * `getEbonds`
  - **Headings:** `ReactiveFF.cpp`, Free Functions

- `RigidMol.cpp.skeleton.md` (0 KB)
  - **Summary:** * `bakeMMFF` * `clear`
  - **Headings:** `RigidMol.cpp`, Free Functions

- `SchroedingerGreen1D_lib.cpp.skeleton.md` (0 KB)
  - **Summary:** * `init` * `init_buffers`
  - **Headings:** `SchroedingerGreen1D_lib.cpp`, Free Functions

- `SchroedingerGreen2D_lib.cpp.skeleton.md` (0 KB)
  - **Summary:** * `init` * `init_buffers`
  - **Headings:** `SchroedingerGreen2D_lib.cpp`, Free Functions

- `utils.cpp.skeleton.md` (0 KB)
  - **Summary:** * `find_Alphabet`
  - **Headings:** `utils.cpp`, Free Functions


## `doc/skeletons/libs_OCL/`

- `MMFFmulti_lib.cpp.skeleton.md` (0 KB)
  - **Summary:** * `change_lvec` * `download`
  - **Headings:** `MMFFmulti_lib.cpp`, Free Functions

- `OCL_GridFF.cpp.skeleton.md` (0 KB)
  - **Summary:** * `cleanup` * `convCoefs`
  - **Headings:** `OCL_GridFF.cpp`, Free Functions


## `doc/skeletons/libs_SDL/`

- `MolGUIlib.cpp.skeleton.md` (0 KB)
  - **Summary:** * `init` * `run`
  - **Headings:** `MolGUIlib.cpp`, Free Functions


## `doc/skeletons/libs_SDL/GLView/`

- `GLView.cpp.skeleton.md` (0 KB)
  - **Summary:** * `Cam::ortho` * `Cam::perspective`
  - **Headings:** `GLView.cpp`, Free Functions

- `GLView.h.skeleton.md` (0 KB)
  - **Summary:** (no prose lines)
  - **Headings:** `GLView.h`

- `GLView.hpp.skeleton.md` (1 KB)
  - **Summary:** * `ASPECT_RATIO` * `ASPECT_RATIO`
  - **Headings:** `GLView.hpp`, Classes,     - `GLView`


## `doc/skeletons/libs_SDL/Vis3D/`

- `Vis3D.cpp.skeleton.md` (0 KB)
  - **Summary:** * `displayLists` * `dragging`
  - **Headings:** `Vis3D.cpp`, Classes,     - `Vis3DApp`, Free Functions


## `doc/skeletons/sketches_SDL/Molecular/`

- `test_BondAdaptedMesh.cpp.skeleton.md` (0 KB)
  - **Summary:** * `agrid` * `bRun`
  - **Headings:** `test_BondAdaptedMesh.cpp`, Classes,     - `TestAppRARFF`, Free Functions

- `test_CLCFGO.cpp.skeleton.md` (0 KB)
  - **Summary:** * `bDrawAtoms` * `bDrawObjects`
  - **Headings:** `test_CLCFGO.cpp`, Classes,     - `TestAppCLCFSF`, Free Functions

- `test_CLCFSF.cpp.skeleton.md` (0 KB)
  - **Summary:** * `bRun` * `fontTex`
  - **Headings:** `test_CLCFSF.cpp`, Classes,     - `TestAppCLCFSF`, Free Functions

- `test_ConfDynamics.cpp.skeleton.md` (0 KB)
  - **Summary:** * `confWorld` * `vobConf`
  - **Headings:** `test_ConfDynamics.cpp`, Classes,     - `TestConfDynamics`, Free Functions

- `test_DirectionStiffness.cpp.skeleton.md` (1 KB)
  - **Summary:** * `atypes` * `bConverged`
  - **Headings:** `test_DirectionStiffness.cpp`, Classes,     - `TestAppDirectionStiffness`, Free Functions

- `test_eFF.cpp.skeleton.md` (0 KB)
  - **Summary:** * `E0` * `Efunc`
  - **Headings:** `test_eFF.cpp`, Classes,     - `TestAppRARFF`, Free Functions

- `test_eFF_MD.cpp.skeleton.md` (0 KB)
  - **Summary:** * `E0` * `Efunc`
  - **Headings:** `test_eFF_MD.cpp`, Classes,     - `TestAppRARFF`, Free Functions

- `test_eFF_old.cpp.skeleton.md` (0 KB)
  - **Summary:** * `E0` * `Efunc`
  - **Headings:** `test_eFF_old.cpp`, Classes,     - `TestAppRARFF`, Free Functions

- `test_eFFMC.cpp.skeleton.md` (0 KB)
  - **Summary:** * `bDrawAtoms` * `bDrawObjects`
  - **Headings:** `test_eFFMC.cpp`, Classes,     - `TestAppCLCFSF`, Free Functions

- `test_EOFF.cpp.skeleton.md` (0 KB)
  - **Summary:** * `grid` * `grid`
  - **Headings:** `test_EOFF.cpp`, Classes,     - `TestAppEOFF`, Free Functions

- `test_ESFF.cpp.skeleton.md` (0 KB)
  - **Summary:** * `grid` * `grid`
  - **Headings:** `test_ESFF.cpp`, Classes,     - `TestAppEOFF`, Free Functions

- `test_FARFF.cpp.skeleton.md` (0 KB)
  - **Summary:** * `Emax` * `Emin`
  - **Headings:** `test_FARFF.cpp`, Classes,     - `TestAppFARFF`, Free Functions

- `test_FF2D.cpp.skeleton.md` (0 KB)
  - **Summary:** * `bRun` * `ff`
  - **Headings:** `test_FF2D.cpp`, Classes,     - `TestAppFF2D`, Free Functions

- `test_FIRE.cpp.skeleton.md` (0 KB)
  - **Summary:** * `Fmax` * `Fmax`
  - **Headings:** `test_FIRE.cpp`, Classes,     - `TestAppFIRE`, Free Functions

- `test_FTRFF.cpp.skeleton.md` (0 KB)
  - **Summary:** * `grid` * `grid`
  - **Headings:** `test_FTRFF.cpp`, Classes,     - `TestAppFTRFF`, Free Functions

- `test_MatchLattice2D.cpp.skeleton.md` (0 KB)
  - **Summary:** * `LM` * `LM`
  - **Headings:** `test_MatchLattice2D.cpp`, Classes,     - `TestAppLatticeMatch2D`, Free Functions

- `test_MMFFmini.cpp.skeleton.md` (1 KB)
  - **Summary:** * `atypes` * `bConverged`
  - **Headings:** `test_MMFFmini.cpp`, Classes,     - `TestAppMMFFmini`, Free Functions

- `test_MMFFsp3.cpp.skeleton.md` (1 KB)
  - **Summary:** * `atypes` * `bAtomsSpheres`
  - **Headings:** `test_MMFFsp3.cpp`, Classes,     - `TestAppMMFFsp3`, Free Functions

- `test_Multipoles.cpp.skeleton.md` (0 KB)
  - **Summary:** * `grid` * `system1`
  - **Headings:** `test_Multipoles.cpp`, Classes,     - `TestAppMultipoles`, Free Functions

- `test_RARFF.cpp.skeleton.md` (0 KB)
  - **Summary:** * `Emax` * `Emin`
  - **Headings:** `test_RARFF.cpp`, Classes,     - `TestAppRARFF`, Free Functions

- `test_RARFF2.cpp.skeleton.md` (0 KB)
  - **Summary:** * `Emax` * `Emin`
  - **Headings:** `test_RARFF2.cpp`, Classes,     - `TestAppRARFF`, Free Functions

- `test_RARFF_SR.cpp.skeleton.md` (0 KB)
  - **Summary:** * `bBlockAddAtom` * `bRun`
  - **Headings:** `test_RARFF_SR.cpp`, Classes,     - `TestAppRARFF`, Free Functions

- `test_RARFFarr.cpp.skeleton.md` (0 KB)
  - **Summary:** * `bBlockAddAtom` * `bRun`
  - **Headings:** `test_RARFFarr.cpp`, Classes,     - `TestAppRARFF`, Free Functions

- `test_RRFF.cpp.skeleton.md` (0 KB)
  - **Summary:** * `grid` * `grid`
  - **Headings:** `test_RRFF.cpp`, Classes,     - `TestAppRRFF`, Free Functions

- `test_RspFF.cpp.skeleton.md` (0 KB)
  - **Summary:** * `bRun` * `ff`
  - **Headings:** `test_RspFF.cpp`, Classes,     - `TestAppRARFF`, Free Functions

- `test_SoftMolecularDynamics.cpp.skeleton.md` (0 KB)
  - **Summary:** * `builder` * `commands`
  - **Headings:** `test_SoftMolecularDynamics.cpp`, Classes,     - `TestAppSoftMolDyn`, Free Functions

- `test_sp3space.cpp.skeleton.md` (0 KB)
  - **Summary:** * `bRun` * `fontTex`
  - **Headings:** `test_sp3space.cpp`, Classes,     - `TestAppSp3Space`, Free Functions

- `test_spRotations.cpp.skeleton.md` (0 KB)
  - **Summary:** * `bRun` * `fontTex`
  - **Headings:** `test_spRotations.cpp`, Classes,     - `TestAppSp3Space`, Free Functions


## `doc/topical_audit/`

- `hermes_topical_overview.md` (1 KB)
  - **Summary:** High-level problem map for FireCore, prepared from existing docs and README files. This file captures what the project is trying to solve, not every implementation detail.
  - **Headings:** Hermes Topical Overview,   - Dominant scientific problems,   - Corresponding computational systems,   - Experimental linkage,   - Cross-cutting code concerns

- `RigidSurfPotential_GridFF.md` (37 KB)
  - **Summary:** This document reviews the GridFF (Grid Force Field) and rigid surface potential evaluation systems in FireCore. These systems compute molecu 
  - **Headings:** Rigid Surface Potential & GridFF: Surface Interaction Evaluation Systems,   - Overview,   - Codemaps Analyzed,   - GridFF Subsystem,     - Purpose,     - Implementation Locations,     - GridFF Generation Pipeline,     - GridFF GPU Evaluation
    ... (+87 more)

- `topical_audit.md` (14 KB)
  - **Summary:** FireCore is a >5 year computational science project that has evolved as a lab-book/journal rather than a polished product. Many topics have  
  - **Headings:** Topical Audit Manifesto,   - Purpose,   - Problem Statement,   - Goals,   - Structure,   - Content Guidelines,   - Usage,   - Scope
    ... (+11 more)


## `doc/root/`

- `DFTB_API_migration.md` (52 KB)
  - **Summary:** CODEMAP: DFTB+ Python Integration: Library Interfaces, Parsers, and OpenCL Grid Projection https://windsurf.com/codemaps/1d6b4b7c-04de-49ef-b581-12cf5bfef54a-fe86ab10a43f3d18
  - **Headings:** USER, SWE-1.6, DFTB Python Modules Migration Documentation,   - 1. Core Python Files to Copy,     - Essential Modules (self-contained package):,     - Optional/Utility Modules (not required for core functionality):,   - 2. External Python Dependencies,     - Required (no external compilation needed):
    ... (+95 more)

- `Global_Coding_Rules.md` (3 KB)
  - **Summary:** Applies to all languages (Python, C/C++, OpenCL) and all workflows. These are the non‑negotiable baseline rules for this repository. 
  - **Headings:** Global Coding Rules,   - Core principles,   - NEVER DO THIS,   - Debugging first,   - Performance (general),   - Style (cross‑language),   - Visualization

- `MolecularViever_3.md` (25 KB)
  - **Summary:** I was testing this approach, but it was not really working, I had problem wih OpenGL context distributed among the browser and viewer. But a 
  - **Headings:**     - Architecture Overview,     - Step 1: Utility for Bitmap Font, make_font.py, Configuration,     - Step 2: New Shaders,     - Step 3: The Main Application (`MolecularBrowser.py`), MolecularBrowser.py, Assuming your project structure allows these imports
    ... (+4 more)

- `OpenCLBase.md` (11 KB)
  - **Summary:** File: `pyBall/OCL/OpenCLBase.py` 
  - **Headings:** OpenCLBase: Minimal OpenCL helper for FireCore,   - Quick start (minimal),   - Architecture and conventions,   - Function reference,   - Buffer management patterns,   - Data transfer, low-level,   - Kernel headers → arguments without boilerplate
    ... (+6 more)

- `PATH_CONFIGURATION.md` (7 KB)
  - **Summary:** This document explains how to configure paths for FireCore on different systems. FireCore uses a flexible configuration system that supports 
  - **Headings:** FireCore Path Configuration Guide,   - Problem,   - Solution,   - Configuration Methods,     - Method 1: Using config.json (Recommended),     - Method 2: Using Environment Variables, DFTB+ Slater-Koster files directory, DFTB+ executable (optional, only needed for DFTB backend)
    ... (+19 more)

- `Project_Structure.md` (3 KB)
  - **Summary:** FireCore is a comprehensive computational chemistry and physics repository containing multiple subprojects, programs, and test scripts. It i 
  - **Headings:** FireCore Project Structure and Development Guide,   - Project Overview,   - Directory Structure and Components,     - 1. Core C++ Implementation (`cpp/`),     - 2. Fortran DFT Implementation (`fortran/`, `fortran2/`),     - 3. Python Interface (`pyBall/`),     - 4. OpenCL Implementation (`pyBall/OCL/`),     - 5. Tests and Examples (`tests/`)
    ... (+4 more)

- `QM_DFT_DFTB_packages_for_integration.md` (82 KB)
  - **Summary:** I'm looking for some quantum simulation code which can be easily integrated to my qmmm package. It should be based on LCAO and I consider bo 
  - **Headings:**   - USER :,     - Recommended Open-Source Quantum Simulation Codes for QM/MM Integration,     - Key Recommendations,   - USER 2,     - PySCF: Forces and Periodic Systems,     - Alternatives for Periodic Systems,     - Semi-Empirical Codes (e.g., Grimme's xTB/GFN2-xTB),   - USER 3
    ... (+28 more)

- `QM_DFT_DFTB_packages_for_integration_sumarized.md` (13 KB)
  - **Summary:** This document summarizes LLM-collected notes about **open-source (and a few reference closed-source) electronic-structure codes** that may b 
  - **Headings:** QM / DFT / DFTB / Semiempirical codes for Python-driven QM/MM integration,   - Selection criteria emphasized in the notes,   - 1) DFTB / Semiempirical (fast),   - 2) DFT (LDA/GGA/Hybrid/meta-hybrid),   - 3) HF / post-HF (reference / high-accuracy, GPU-hybrid focused),   - To investigate next (missing info / tick-the-box gaps)

- `README.md` (9 KB)
  - **Summary:** This directory contains comprehensive documentation, development notes, mathematical derivations, and prototype implementations for the Fire 
  - **Headings:** Documentation and Development Notes,   - Directory Structure,   - Technical Documentation (Markdown/),     - GUI System (GUI.md),     - Molecular Mechanics (MMFFBuilder.md),     - Electron Force Field (eFF_documentation.md),   - Development Notes (DevNotes/),     - GUI Development (ToDo_GUI.md)
    ... (+34 more)
