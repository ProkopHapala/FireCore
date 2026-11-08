# File Descriptions
One sentence per file, added as files are read.

## OnSurfaceAssembly
- `/home/prokophapala/git/FireCore/doc/Topics/OnSurfaceAssembly/GridFF_FDBM_Fitting.md` — Records a user request and assistant analysis for fitting a Full Density-Based Model (FDBM) to small-molecule/NaCl DFT interaction data using linear least squares on pre-sampled GridFF features (Pauli, London, Coulomb).
- `/home/prokophapala/git/FireCore/doc/Topics/OnSurfaceAssembly/GridFF_atoms_alignment.md` — Discusses grid-to-atom alignment debugging for GridFF by using 1D/2D/3D extrema detection on forcefield grids, plus an `.npz` bundling strategy to preserve origin metadata.
- `/home/prokophapala/git/FireCore/doc/Topics/OnSurfaceAssembly/FoldedSubstratePotential_OpenCL.md` — Plans a separable Fourier/exponential folded-substrate potential (sin/cos in XY, exp in Z) fit with periodic boundary conditions and evaluated via OpenCL kernels including `getSurfFolded_workgroup` and `getSurfFolded_harmonics`.
- `/home/prokophapala/git/FireCore/doc/Topics/OnSurfaceAssembly/Ewald_2D.md` — Surveys analytical 2D Ewald and reciprocal-space formulas to represent the surface electrostatic potential as a compact cosine/exponential basis, plus a Python/NumPy program for generating XY and XZ slice plots from ion coordinates.
- `/home/prokophapala/git/FireCore/doc/Topics/OnSurfaceAssembly/RigidBodyAFM.md` — Analyzes existing GPU rigid-body modules (`RigidBodyDynamics.py`, `Rigid.cl`) and test scripts to identify the missing GridFF+anchor path needed for tip-attached molecule AFM simulations.
- `/home/prokophapala/git/FireCore/doc/Topics/OnSurfaceAssembly/AFM_File_Analysis_and_Architecture.md` — Maps the GPU/CPU backends, CLI/GUI frontends, visualization tools, and test scripts relevant to AFM simulation; proposes `RigidBodyAFM.py` backend plus separate CLI and VisPy frontends.
- `/home/prokophapala/git/FireCore/doc/Topics/OnSurfaceAssembly/GridFF_RelaxedScan_cpp_notes.md` — Describes the C++/Python spline-anchored manipulation system (`TipSplineSAOptimizer`, `run_tipSpline_scan.py`) and documents generated plots, optimization results, and planned GPU migration.
- `/home/prokophapala/git/FireCore/doc/Topics/OnSurfaceAssembly/ExplorerVisPy_GridFF_Interactive_Plan.md` — Tracks implementation and debugging of interactive GridFF scanning in `ExplorerVisPy.py`, including Morse-derived PLQ coefficients, Coulomb reference subtraction, GridFF caching fixes, and parity fixes between fast GPU and reference paths.
- `/home/prokophapala/git/FireCore/doc/Topics/OnSurfaceAssembly/ElectrostaticContinuumEmbeding_report.md` — Documents a GPU-ready macroscopic electrostatic correction for PBC scans using monopole/dipole/quadrupole sheet formulas, validated to <1e-4 eV primitive-translation equivalence.
- `/home/prokophapala/git/FireCore/doc/Topics/OnSurfaceAssembly/SurfacePotentialVisualization_report.md` — Post-mortem of substrate surface visualization reporting validation results across GridFF, XYZ reference, and GPU-accelerated backends in VisPy/headless renderers.

## AFM
- `/home/prokophapala/git/FireCore/doc/Topics/AFM/AFM.md` — Describes the full AFM simulation pipeline (5-step workflow) and maps which code files implement each stage.
- `/home/prokophapala/git/FireCore/doc/Topics/AFM/AFM_migration_plan.md` — Outlines a 2025 plan to migrate the C++/OpenCL AFM host code to pyOpenCL while preserving API compatibility.
- `/home/prokophapala/git/FireCore/doc/Topics/AFM/AFM_migration_progress.md` — Tracks the pyOpenCL migration progress, completed milestones, and current open issues.
- `/home/prokophapala/git/FireCore/doc/Topics/AFM/AFM_migration_discusion.md` — Documents an LLM-mediated discussion consolidating AFM kernel/host equivalences and proposing a safe non-duplication strategy before implementation.
- `/home/prokophapala/git/FireCore/doc/Topics/AFM/AFM_FDBM_DFTB.md` — Documents integration of the Full Density Based Model (FDBM) with DFTB-generated substrate grids for AFM simulations.
- `/home/prokophapala/git/FireCore/doc/Topics/AFM/AFM_FDBM_fitting.md` — Defines workflows for fitting FDBM parameters to DFT reference data using linear least squares.
- `/home/prokophapala/git/FireCore/doc/Topics/AFM/AFM_FDBM_optimization.md` — Details performance optimizations applied to the FDBM implementations.
- `/home/prokophapala/git/FireCore/doc/Topics/AFM/AFM_FDBM_profiling_optimization.md` — Reports profiling results and further optimization steps for FDBM code paths.
- `/home/prokophapala/git/FireCore/doc/Topics/AFM/AFM_FDBM_pySCF.md` — Covers PySCF integration for FDBM electron density and electrostatic potentials, plus embedded historical discussion of normalization fixes and Green’s-function STM approaches.

## STM
- `/home/prokophapala/git/FireCore/doc/Topics/STM/STM_GF_new.md` — Describes the Green’s-function-based STM simulation workflow for calculating tunnel currents and density of states.
- `/home/prokophapala/git/FireCore/doc/Topics/STM/STM_GPU_QMMM.md` — Long technical chat log planning a GPU-accelerated QM/MM STM junction simulator combining classical MD with Green’s-function transport and perturbative density-of-states coupling.

## ForceFields
- `/home/prokophapala/git/FireCore/doc/Topics/ForceFields/GridFF_fit_CG.md` — Describes adding a conjugate-gradient fitting path alongside damped MD for 3D B-spline GridFF grids, including OpenCL kernels and test harness.
- `/home/prokophapala/git/FireCore/doc/Topics/ForceFields/GridFF_fit_jacobi.md` — Long derivation and chat-driven transition from the original damped-MD Bspline grid fitter to a formally derived Jacobi iterative solver, including code samples for the OpenCL `jacobi_step` kernel and a matching pyOpenCL driver.
- `/home/prokophapala/git/FireCore/doc/Topics/ForceFields/NeighborExclusion_Summary.md` — Specifies an amortized O(1) neighbor-exclusion scheme for OpenCL non-bonded kernels using sorted per-atom exclusion lists packed into `uint32_t`.
- `/home/prokophapala/git/FireCore/doc/Topics/ForceFields/UFF_multiSystem.md` — Design notes for running UFF over multiple replicas on the GPU, including host-side per-system packing and 2D NDRange kernel launches.
- `/home/prokophapala/git/FireCore/doc/Topics/ForceFields/UFF_multiSystem_summary.md` — Consolidated reference rules from a multi-session UFF GPU debug: offset indexing, context isolation, and CPU/GPU parity protocol across bonds, angles, dihedrals, inversions, and GridFF interaction modes.

## FTIR_Nanocrystals
- `/home/prokophapala/git/FireCore/doc/Topics/FTIR_Nanocrystals/Hessian_fitting.md` — Explains fitting procedures and practical heuristics for constructing disordered nanocrystal Hessians using 3D B-spline grids.
- `/home/prokophapala/git/FireCore/doc/Topics/FTIR_Nanocrystals/gen_nanocrystals.md` — Generates random nanocrystal geometries by stacking atoms onto substrate layers with element typing.

## Kekule_Topology
- `/home/prokophapala/git/FireCore/doc/Topics/Kekule_Topology/Kekule_Optimizer.md` — Describes an optimizer for Kekulé bond patterns and related resonance topology tasks.
- `/home/prokophapala/git/FireCore/doc/Topics/Kekule_Topology/Kekule_Topology.md` — Defines graph-topology operations for aromatic bond patterns.
- `/home/prokophapala/git/FireCore/doc/Topics/Kekule_Topology/GrapheneRibbonBuilder.py` — Script that builds armchair or zigzag graphene ribbon structures from a width parameter.
- `/home/prokophapala/git/FireCore/doc/Topics/Kekule_Topology/GrapheneRibbonGUI.py` — Small GUI wrapper around graphene ribbon geometry generation.
- `/home/prokophapala/git/FireCore/doc/Topics/Kekule_Topology/KekuleCLI.py` — Command-line interface wrapper for calling the Kekulé topology solver from a terminal.
- `/home/prokophapala/git/FireCore/doc/Topics/Kekule_Topology/KekuleGUI.py` — Graphical user interface that visualizes and edits Kekulé bond arrangements on aromatic lattices.
- `/home/prokophapala/git/FireCore/doc/Topics/Kekule_Topology/KekuleSolver.py` — Core algorithm for solving Kekulé patterns on honeycomb lattices.
- `/home/prokophapala/git/FireCore/doc/Topics/Kekule_Topology/Kekule_DeepSeek.py` — GPT/DeepSeek-assisted exploration of Kekulé topology rules and bond-counting heuristics.
- `/home/prokophapala/git/FireCore/doc/Topics/Kekule_Topology/Kekule_Gemini.py` — GPT/Gemini-assisted session on topological constraints for Kekulé resonance structures.
- `/home/prokophapala/git/FireCore/doc/Topics/Kekule_Topology/Kekule_Opt_Deepseek.py` — LLM-driven session on optimization formulations for Kekulé bond assignments.
- `/home/prokophapala/git/FireCore/doc/Topics/Kekule_Topology/Kekule_Opt_Gemini.py` — LLM-driven session comparing optimization strategies for Kekulé pattern selection.

## LinearLagebra
- `/home/prokophapala/git/FireCore/doc/Topics/LinearLagebra/ChebyshevAndInterialAccelerationOfLinearSolvers.md` — Discusses Chebyshev and inertial acceleration techniques for iterative linear solvers.
- `/home/prokophapala/git/FireCore/doc/Topics/LinearLagebra/Chebyshev_Acceleration.md` — Derives and benchmarks Chebyshev-accelerated iterative solvers for grids.
- `/home/prokophapala/git/FireCore/doc/Topics/LinearLagebra/NumricalStabilityLinearSolvers.md` — Analyzes numerical stability issues appearing in parallel linear system solvers.
- `/home/prokophapala/git/FireCore/doc/Topics/LinearLagebra/Parallel_GaussSeidel.md` — Describes parallel Gauss-Seidel implementations suitable for GPU/OpenCL execution.

## ManyBody
- `/home/prokophapala/git/FireCore/doc/Topics/ManyBody/HubbardSolver.md` — Documents a Hubbard-model solver for strongly correlated lattice physics.
- `/home/prokophapala/git/FireCore/doc/Topics/ManyBody/HubbardSolver_new.md` — Updated implementation notes for the lattice Hubbard solver with revised numerics.
- `/home/prokophapala/git/FireCore/doc/Topics/ManyBody/MC_PME_dicussion.md` — Discussion of Ewald-summation-based Monte Carlo sampling for many-body systems.
- `/home/prokophapala/git/FireCore/doc/Topics/ManyBody/MC_PME_progress.md` — Progress log for implementing fast Ewald PME screening in a Monte Carlo workflow.
- `/home/prokophapala/git/FireCore/doc/Topics/ManyBody/MQCA_GrayCodeSolver.md` — Describes a Gray-code-based solver for many-body quantum or combinatorial optimization problems.

## MolecularGUI
- `/home/prokophapala/git/FireCore/doc/Topics/MolecularGUI/MolViewer_PyQt.md` — Documents a PyQt-based molecular viewer for structural inspection.
- `/home/prokophapala/git/FireCore/doc/Topics/MolecularGUI/MolViewer_PyQt_2.md` — Updated notes for the PyQt viewer with refined interaction and rendering.
- `/home/prokophapala/git/FireCore/doc/Topics/MolecularGUI/MolViewer_PyQt_toGemini.md` — Chat-style conversation guiding the design of PyQt molecular viewer features using a Gemini assist.
- `/home/prokophapala/git/FireCore/doc/Topics/MolecularGUI/MoleculeEditor2D.md` — Plans and notes for a 2D molecular topology editor GUI.
- `/home/prokophapala/git/FireCore/doc/Topics/MolecularGUI/OpenGL_Viewer_Documentation.md` — Documents the OpenGL viewer workflow for molecular scenes and screenshots.

## MoleculeBuilder
- `/home/prokophapala/git/FireCore/doc/Topics/MoleculeBuilder.md` — Captures discussions and design ideas for building and manipulating molecular structures.

## PathDiffusion
- `/home/prokophapala/git/FireCore/doc/Topics/PathDiffusion/ManipulationPathOpt_Report.md` — Reports optimized manipulation-path results combined with cheminformatics and diffusion-style optimization tools.
- `/home/prokophapala/git/FireCore/doc/Topics/PathDiffusion/PathDiffusion_design.md` — Describes the architecture of a path-diffusion module for molecular-manipulation trajectories.
- `/home/prokophapala/git/FireCore/doc/Topics/PathDiffusion/PathDiffusion_discussion.md` — Discussion log on trajectory optimization, sampling, and reuse.
- `/home/prokophapala/git/FireCore/doc/Topics/PathDiffusion/PathDiffusion_discussion_old.md` — Older discussion thread on path diffusion designs superseded by the current design document.

## RigidBodyAssembly
- `/home/prokophapala/git/FireCore/doc/Topics/RigidBodyAssembly/RigidBodyAssemblyDiscussion.md` — Discussion on general rigid-body assembly workflows and constraints.
- `/home/prokophapala/git/FireCore/doc/Topics/RigidBodyAssembly/RigidBodyAssemblyDiscussion_tutorial.md` — Tutorial-oriented notes on rigid-body assembly and configuration sampling.
- `/home/prokophapala/git/FireCore/doc/Topics/RigidBodyAssembly/RigidGridFF_H2O_NaCl_report.md` — Detailed multi-replica rigid GridFF relaxation report and tutorial for H2O on NaCl and PTCDA on NaCl.

## TrussSolver
- `/home/prokophapala/git/FireCore/doc/Topics/TrussSolver/High_Precision_Physics_on_singlepoint_GPU_hardware.md` — Outlines strategies for high-precision physics computations on consumer single-precision GPU hardware.
- `/home/prokophapala/git/FireCore/doc/Topics/TrussSolver/ProjectiveDynamics.md` — Describes a projective-dynamics solver for truss-like systems.
- `/home/prokophapala/git/FireCore/doc/Topics/TrussSolver/ProjectiveDynamics_wegl_glsl.md` — Documents WebGL/GLSL implementations of projective dynamics for web-based viewing.
- `/home/prokophapala/git/FireCore/doc/Topics/TrussSolver/VertexBlockDescent.md` — Describes a vertex-block descent solver for large truss or mesh optimization problems.

## molgui_web
- `/home/prokophapala/git/FireCore/doc/Topics/molgui_web/Buckets_GridNeighborDesign.md` — Designs bucket/grid neighbor-lookup structures for the web molecular GUI.
- `/home/prokophapala/git/FireCore/doc/Topics/molgui_web/Crystal_Generator.md` — Documents crystal-structure generation helpers in the web GUI context.
- `/home/prokophapala/git/FireCore/doc/Topics/molgui_web/Editor_for_Polymers_on_Surfaces_and_Edges.md` — Describes editor interaction patterns for polymers on surfaces and edge sites.
- `/home/prokophapala/git/FireCore/doc/Topics/molgui_web/Molecular_topology_editor_with_explicit_rings.md` — Discusses explicit ring support in the molecular topology editor.
- `/home/prokophapala/git/FireCore/doc/Topics/molgui_web/Molecular_topology_js_plan.md` — JavaScript-centric implementation plan for molecular topology editing in the browser.
- `/home/prokophapala/git/FireCore/doc/Topics/molgui_web/Step_Edge_Editor_ideas.md` — Brainstorming notes for a step-edge molecular editor UI.
- `/home/prokophapala/git/FireCore/doc/Topics/molgui_web/molgui_web_dev.md` — General web GUI development notes and pointers for the molecular editor frontend.
