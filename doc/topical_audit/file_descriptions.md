# File Descriptions
One sentence per file, added as files are read.

## OnSurfaceAssembly
- `doc/Topics/OnSurfaceAssembly/GridFF_FDBM_Fitting.md` — Records a user request and assistant analysis for fitting a Full Density-Based Model (FDBM) to small-molecule/NaCl DFT interaction data using linear least squares on pre-sampled GridFF features (Pauli, London, Coulomb).
- `doc/Topics/OnSurfaceAssembly/GridFF_atoms_alignment.md` — Discusses grid-to-atom alignment debugging for GridFF by using 1D/2D/3D extrema detection on forcefield grids, plus an `.npz` bundling strategy to preserve origin metadata.
- `doc/Topics/OnSurfaceAssembly/FoldedSubstratePotential_OpenCL.md` — Plans a separable Fourier/exponential folded-substrate potential (sin/cos in XY, exp in Z) fit with periodic boundary conditions and evaluated via OpenCL kernels including `getSurfFolded_workgroup` and `getSurfFolded_harmonics`.
- `doc/Topics/OnSurfaceAssembly/Ewald_2D.md` — Surveys analytical 2D Ewald and reciprocal-space formulas to represent the surface electrostatic potential as a compact cosine/exponential basis, plus a Python/NumPy program for generating XY and XZ slice plots from ion coordinates.
- `doc/Topics/OnSurfaceAssembly/RigidBodyAFM.md` — Analyzes existing GPU rigid-body modules (`RigidBodyDynamics.py`, `Rigid.cl`) and test scripts to identify the missing GridFF+anchor path needed for tip-attached molecule AFM simulations.
- `doc/Topics/OnSurfaceAssembly/AFM_File_Analysis_and_Architecture.md` — Maps the GPU/CPU backends, CLI/GUI frontends, visualization tools, and test scripts relevant to AFM simulation; proposes `RigidBodyAFM.py` backend plus separate CLI and VisPy frontends.
- `doc/Topics/OnSurfaceAssembly/GridFF_RelaxedScan_cpp_notes.md` — Describes the C++/Python spline-anchored manipulation system (`TipSplineSAOptimizer`, `run_tipSpline_scan.py`) and documents generated plots, optimization results, and planned GPU migration.
- `doc/Topics/OnSurfaceAssembly/ExplorerVisPy_GridFF_Interactive_Plan.md` — Tracks implementation and debugging of interactive GridFF scanning in `ExplorerVisPy.py`, including Morse-derived PLQ coefficients, Coulomb reference subtraction, GridFF caching fixes, and parity fixes between fast GPU and reference paths.
- `doc/Topics/OnSurfaceAssembly/ElectrostaticContinuumEmbeding_report.md` — Documents a GPU-ready macroscopic electrostatic correction for PBC scans using monopole/dipole/quadrupole sheet formulas, validated to <1e-4 eV primitive-translation equivalence.
- `doc/Topics/OnSurfaceAssembly/SurfacePotentialVisualization_report.md` — Post-mortem of substrate surface visualization reporting validation results across GridFF, XYZ reference, and GPU-accelerated backends in VisPy/headless renderers.

## AFM
- `doc/Topics/AFM/AFM.md` — Describes the full AFM simulation pipeline (5-step workflow) and maps which code files implement each stage.
- `doc/Topics/AFM/AFM_migration_plan.md` — Outlines a 2025 plan to migrate the C++/OpenCL AFM host code to pyOpenCL while preserving API compatibility.
- `doc/Topics/AFM/AFM_migration_progress.md` — Tracks the pyOpenCL migration progress, completed milestones, and current open issues.
- `doc/Topics/AFM/AFM_migration_discusion.md` — Documents an LLM-mediated discussion consolidating AFM kernel/host equivalences and proposing a safe non-duplication strategy before implementation.
- `doc/Topics/AFM/AFM_FDBM_DFTB.md` — Documents integration of the Full Density Based Model (FDBM) with DFTB-generated substrate grids for AFM simulations.
- `doc/Topics/AFM/AFM_FDBM_fitting.md` — Defines workflows for fitting FDBM parameters to DFT reference data using linear least squares.
- `doc/Topics/AFM/AFM_FDBM_optimization.md` — Details performance optimizations applied to the FDBM implementations.
- `doc/Topics/AFM/AFM_FDBM_profiling_optimization.md` — Reports profiling results and further optimization steps for FDBM code paths.
- `doc/Topics/AFM/AFM_FDBM_pySCF.md` — Covers PySCF integration for FDBM electron density and electrostatic potentials, plus embedded historical discussion of normalization fixes and Green’s-function STM approaches.

## STM
- `doc/Topics/STM/STM_GF_new.md` — Describes the Green’s-function-based STM simulation workflow for calculating tunnel currents and density of states.
- `doc/Topics/STM/STM_GPU_QMMM.md` — Long technical chat log planning a GPU-accelerated QM/MM STM junction simulator combining classical MD with Green’s-function transport and perturbative density-of-states coupling.

## ForceFields
- `doc/Topics/ForceFields/GridFF_fit_CG.md` — Describes adding a conjugate-gradient fitting path alongside damped MD for 3D B-spline GridFF grids, including OpenCL kernels and test harness.
- `doc/Topics/ForceFields/GridFF_fit_jacobi.md` — Long derivation and chat-driven transition from the original damped-MD Bspline grid fitter to a formally derived Jacobi iterative solver, including code samples for the OpenCL `jacobi_step` kernel and a matching pyOpenCL driver.
- `doc/Topics/ForceFields/NeighborExclusion_Summary.md` — Specifies an amortized O(1) neighbor-exclusion scheme for OpenCL non-bonded kernels using sorted per-atom exclusion lists packed into `uint32_t`.
- `doc/Topics/ForceFields/UFF_multiSystem.md` — Design notes for running UFF over multiple replicas on the GPU, including host-side per-system packing and 2D NDRange kernel launches.
- `doc/Topics/ForceFields/UFF_multiSystem_summary.md` — Consolidated reference rules from a multi-session UFF GPU debug: offset indexing, context isolation, and CPU/GPU parity protocol across bonds, angles, dihedrals, inversions, and GridFF interaction modes.

## FTIR_Nanocrystals

Topic index: `doc/Topics/FTIR_Nanocrystals/README.md`. Audit: `doc/topical_audit/Nanocrystal_Vibrations.md`. Hub: `tests/tSiNCs/README.md`. Open items: `tests/tSiNCs/ToDo_Nanocrystal.md`.

- `doc/Topics/FTIR_Nanocrystals/Nanocrystal_NPZ_Pipeline.guide.md` — NPZ stages 01→05, relax/Hessian/spectrum consumers.
- `doc/Topics/FTIR_Nanocrystals/NPZ_Crystal_Schema.md` — Per-key NPZ contract (v1.2).
- `doc/Topics/FTIR_Nanocrystals/CPP_MolecularBrowser_NPZ.md` — C++ SDL NPZ molecular browser.
- `doc/Topics/FTIR_Nanocrystals/Python_Vispy_MolBrowser_Plugins.md` — Python Vispy browser + vibration plugin.
- `scripts/export_nanocrystal_bundle.mjs` — Export `01_crystal.npz` + `03_topology.npz` with defects and surface AABBs.
- `web/molgui_webgpu/NanocrystalExport.js` — Unified JS export bundle API.
- `pyBall/io/crystal_npz.py` — Shared NPZ loaders (Python viewers + pipeline).
- `cpp/common/io/` — NPY/NPZ decode for C++ MolecularBrowser.
- `pyBall/nanocrystal_pipeline.py` — NPZ relax, topology-linear Hessian, spectrum stages.
- `doc/Topics/FTIR_Nanocrystals/gen_nanocrystals.chat.md` — Generation CLI, Miller cuts, capping, bridge defects.
- `tests/tSiNCs/nanocrystals.mjs` — Unified nanocrystal CLI (`generate`, `ensemble`, topology, audit). Canonical location; `scripts/` copies deprecated.
- `tests/tSiNCs/gen_nanocrystals.mjs` — Deprecated single-crystal wrapper; use `nanocrystals.mjs generate`.
- `tests/tSiNCs/export_nanocrystal_bundle.mjs` — Crystal + topology NPZ bundle for viewer fixtures.
- `doc/Topics/FTIR_Nanocrystals/Phonon_testing.guide.md` — Phonon testing with MMFF; PBC/ASR pitfalls.
- `pyBall/FTIR.py` — Linear response vibration spectroscopy: Green's function probing, topology-linear Hessian, rigid mode projection.
- `pyBall/MMFF.py` — Python bindings for MMFF/UFF force fields with `getHessian3Nx3N`, `getPhononPhiBlocks`, and generalized `setBondParamsByType`/`setAngleParamsByType`.
- `cpp/libs/Molecular/MMFF_lib.cpp` — C++ implementation of Hessian and phonon block computation via central finite differences, plus Python buffer exposure.
- `cpp/common/molecular/NBFF.h` — Base non-bonded force field class with reusable `checkPBCNeighCells()` validation.
- `tests/tMMFF/test_diamond_phonon_bands.py` — Diamond phonon band structure using supercell + Bloch extraction; supports PBC, ASR, UFF/MMFF toggle, and parameter scaling.
- `tests/tMMFF/test_vibration_spectra.py` — End-to-end linear response FTIR spectrum: MMFF → Hessian → mass matrix → rigid mode projection → frequency probe.
- `tests/tMMFF/test_hessian_fitting.py` — Fit UFF bond/angle stiffness parameters to reference Hessian using linear least squares.
- `tests/tMMFF/test_diamond_gamma.py` — Γ-point phonon frequencies for diamond (simplified, no k-path).
- `tests/tMMFF/test_ethane_gamma.py` — Γ-point frequencies for ethane (molecular, no PBC).
- `tests/tMMFF/test_diatomic_hessian.py` — Hessian calculation validation for diatomic molecules.
- `tests/tMMFF/test_diamond_phonon_bands.py` — Diamond phonon bands; plotting inline (standalone `plot_phonon_bands.py` not yet extracted).
- `tests/tMMFF/plot_mmff_comparison.py` — Compare MMFF vs reference phonon bands.
- `tests/tMMFF/plot_mmff_scaling.py` — Overlay default vs scaled parameter phonon bands.
- `tests/tMMFF/run_hessian.py` — Batch Hessian runner for multiple systems.

## Molecular_Topology
- `doc/topical_audit/molecular_topology.md` — Cross-language audit of molecular graph representations, bond finding, neighbor lists, ring/bridge detection, and hybridization geometry (C++/Python/JS).
- `doc/topical_audit/molecular_topology_types.md` — Cross-language audit of atom type assignment: parameter files, type resolution, hybridization detection, pi-fragments, and bond parameters.
- `doc/topical_audit/molecular_topology_editors.md` — Cross-language audit of interactive molecular editors, GUI architectures, advanced editing operations, and consolidation roadmap.
- `tests/tUFF/data_UFF/ElementTypes.dat` — Canonical element database: symbols, atomic numbers, covalent/vdW radii, colors, UFF/QEq params. Loaded by C++/Python/JS parameter parsers.
- `tests/tUFF/data_UFF/AtomTypes.dat` — Canonical atom type database: hybridized types (C_sp2, O_hydroxyl, etc.), valence, pi orbitals, electron pairs, MMFF/UFF parameters. Loaded by all three language implementations.
- `tests/tUFF/data_UFF/BondTypes.dat` — Canonical bond parameter table: bond length `l0` and stiffness `k` indexed by atom type pair and bond order.
- `tests/tUFF/data_UFF/AngleTypes.dat` — Canonical angle parameter table: equilibrium angle `ang0` and stiffness `k` indexed by atom type triple.
- `tests/tUFF/data_UFF/DihedralTypes.dat` — Canonical dihedral/torsion parameter table: barrier height `k`, phase `ang0`, and periodicity indexed by atom type quadruple and bond order.

## Kekule_Topology
- `doc/Topics/Kekule_Topology/Kekule_Optimizer.md` — Describes an optimizer for Kekulé bond patterns and related resonance topology tasks.
- `doc/Topics/Kekule_Topology/Kekule_Topology.md` — Defines graph-topology operations for aromatic bond patterns.
- `doc/Topics/Kekule_Topology/GrapheneRibbonBuilder.py` — Script that builds armchair or zigzag graphene ribbon structures from a width parameter.
- `doc/Topics/Kekule_Topology/GrapheneRibbonGUI.py` — Small GUI wrapper around graphene ribbon geometry generation.
- `doc/Topics/Kekule_Topology/KekuleCLI.py` — Command-line interface wrapper for calling the Kekulé topology solver from a terminal.
- `doc/Topics/Kekule_Topology/KekuleGUI.py` — Graphical user interface that visualizes and edits Kekulé bond arrangements on aromatic lattices.
- `doc/Topics/Kekule_Topology/KekuleSolver.py` — Core algorithm for solving Kekulé patterns on honeycomb lattices.
- `doc/Topics/Kekule_Topology/Kekule_DeepSeek.py` — GPT/DeepSeek-assisted exploration of Kekulé topology rules and bond-counting heuristics.
- `doc/Topics/Kekule_Topology/Kekule_Gemini.py` — GPT/Gemini-assisted session on topological constraints for Kekulé resonance structures.
- `doc/Topics/Kekule_Topology/Kekule_Opt_Deepseek.py` — LLM-driven session on optimization formulations for Kekulé bond assignments.
- `doc/Topics/Kekule_Topology/Kekule_Opt_Gemini.py` — LLM-driven session comparing optimization strategies for Kekulé pattern selection.
- `pyBall/KekuleBackend.py` — Backend for Kekule Structure Explorer using `AtomicGraph` as authoritative state; hex grid editing, ring detection, passivation groups, H-capping.
- `pyBall/KekuleExplorerGUI.py` — PyQt5/Vispy GUI frontend for the Kekule Structure Explorer with interactive 3D rendering and extension system.
- `pyBall/Kekule.py` — Python ctypes wrapper for C++ Kekule_lib (bond-order optimizer with init, relax, setDefaultBondOrders, pinBondOrders).

## Kekule_Tests
- `tests/tKekuleExplorer/test_kekule_topology.py` — Test suite for Kekule Structure Explorer topology operations.
- `tests/tKekuleExplorer/test_kekule_hbonds.py` — Test suite for hydrogen bond detection in Kekule structures.
- `tests/tKekuleExplorer/test_kekule_relax.py` — Test suite for geometry relaxation of Kekule structures.
- `tests/tKekuleExplorer/test_ribbon_parity.py` — Parity tests for graphene ribbon Kekule patterns.

## Kekule_BondOrder
- `tests/tKekule/run.py` — Test script for `pyBall.Kekule` bond-order relaxation on anthracene fragments.
- `tests/tKekule/generate_2.py` — Generate Kekule structure variants with donor/acceptor groups (-OH, -NH2, =O, =N-).

## H_transfer_NEB
- `tests/pyFireball/neb_h_transfer_molecules.py` — NEB calculation for H-transfer between hydro/dehydro molecules (phenazine, pyrazine) with rigid_scan, relax_scan, and dry_run modes.
- `tests/pyFireball/neb_h_transfer.py` — NEB for H-transfer between N-passivated and NH-passivated ribbons with periodic k-point sampling.
- `tests/pyFireball/build_two_ribbons.py` — Build two-ribbon unit cell with N/NH passivation for H-transfer studies using GrapheneRibbonBuilder.
- `tests/pyFireball/scan_ribbon_dftb.py` — Lattice scan for zigzag graphene ribbons using DFTB+ with PBC along x.
- `tests/pyFireball/scan_constrained.py` — Constrained H-transfer scan between ribbons fixing N atoms and relaxing others via DFTB+.
- `tests/pyFireball/scan_LHb.py` — Scan energy versus H-bond length for two-ribbon systems using DFTB+.
- `tests/pyFireball/relax_ribbon.py` — Lattice scan for N-doped graphene ribbons (pyridinic, pyrrolic) using FireCore relax.

## DFTB_Utils
- `pyBall/dftb_utils.py` — Subprocess-based DFTB+ interface with HSD generation, waveplot integration, Hessian/vibration utilities, and PBC support.
- `pyBall/dftb_lib.py` — C-API wrapper for libdftbplus.so via ctypes with Hamiltonian/overlap/density matrix access.
- `pyBall/FireCore.py` — Python bindings for the Fireball Fortran library with PBC and k-point support.
- `tests/dftb/example_hessian.py` — Compute Hessian via DFTB+ SecondDerivatives driver and extract vibrational frequencies.
- `tests/dftb/example_orbitals.py` — Run DFTB+ with eigenvector output, run waveplot, and plot orbital 2D slices from cube files.


## LinearLagebra
- `doc/Topics/LinearLagebra/ChebyshevAndInterialAccelerationOfLinearSolvers.md` — Discusses Chebyshev and inertial acceleration techniques for iterative linear solvers.
- `doc/Topics/LinearLagebra/Chebyshev_Acceleration.md` — Derives and benchmarks Chebyshev-accelerated iterative solvers for grids.
- `doc/Topics/LinearLagebra/NumricalStabilityLinearSolvers.md` — Analyzes numerical stability issues appearing in parallel linear system solvers.
- `doc/Topics/LinearLagebra/Parallel_GaussSeidel.md` — Describes parallel Gauss-Seidel implementations suitable for GPU/OpenCL execution.

## ManyBody
- `doc/Topics/ManyBody/HubbardSolver.md` — Documents a Hubbard-model solver for strongly correlated lattice physics.
- `doc/Topics/ManyBody/HubbardSolver_new.md` — Updated implementation notes for the lattice Hubbard solver with revised numerics.
- `doc/Topics/ManyBody/MC_PME_dicussion.md` — Discussion of Ewald-summation-based Monte Carlo sampling for many-body systems.
- `doc/Topics/ManyBody/MC_PME_progress.md` — Progress log for implementing fast Ewald PME screening in a Monte Carlo workflow.
- `doc/Topics/ManyBody/MQCA_GrayCodeSolver.md` — Describes a Gray-code-based solver for many-body quantum or combinatorial optimization problems.

## MolecularGUI
- `doc/Topics/MolecularGUI/MolViewer_PyQt.md` — Documents a PyQt-based molecular viewer for structural inspection.
- `doc/Topics/MolecularGUI/MolViewer_PyQt_2.md` — Updated notes for the PyQt viewer with refined interaction and rendering.
- `doc/Topics/MolecularGUI/MolViewer_PyQt_toGemini.md` — Chat-style conversation guiding the design of PyQt molecular viewer features using a Gemini assist.
- `doc/Topics/MolecularGUI/MoleculeEditor2D.md` — Plans and notes for a 2D molecular topology editor GUI.
- `doc/Topics/MolecularGUI/OpenGL_Viewer_Documentation.md` — Documents the OpenGL viewer workflow for molecular scenes and screenshots.

## MoleculeBuilder
- `doc/Topics/MoleculeBuilder.md` — Captures discussions and design ideas for building and manipulating molecular structures.

## PathDiffusion
- `doc/Topics/PathDiffusion/ManipulationPathOpt_Report.md` — Reports optimized manipulation-path results combined with cheminformatics and diffusion-style optimization tools.
- `doc/Topics/PathDiffusion/PathDiffusion_design.md` — Describes the architecture of a path-diffusion module for molecular-manipulation trajectories.
- `doc/Topics/PathDiffusion/PathDiffusion_discussion.md` — Discussion log on trajectory optimization, sampling, and reuse.
- `doc/Topics/PathDiffusion/PathDiffusion_discussion_old.md` — Older discussion thread on path diffusion designs superseded by the current design document.

## RigidBodyAssembly
- `doc/Topics/RigidBodyAssembly/RigidBodyAssemblyDiscussion.md` — Discussion on general rigid-body assembly workflows and constraints.
- `doc/Topics/RigidBodyAssembly/RigidBodyAssemblyDiscussion_tutorial.md` — Tutorial-oriented notes on rigid-body assembly and configuration sampling.
- `doc/Topics/RigidBodyAssembly/RigidGridFF_H2O_NaCl_report.md` — Detailed multi-replica rigid GridFF relaxation report and tutorial for H2O on NaCl and PTCDA on NaCl.

## TrussSolver
- `doc/Topics/TrussSolver/High_Precision_Physics_on_singlepoint_GPU_hardware.md` — Outlines strategies for high-precision physics computations on consumer single-precision GPU hardware.
- `doc/Topics/TrussSolver/ProjectiveDynamics.md` — Describes a projective-dynamics solver for truss-like systems.
- `doc/Topics/TrussSolver/ProjectiveDynamics_wegl_glsl.md` — Documents WebGL/GLSL implementations of projective dynamics for web-based viewing.
- `doc/Topics/TrussSolver/VertexBlockDescent.md` — Describes a vertex-block descent solver for large truss or mesh optimization problems.

## molgui_web
- `doc/Topics/molgui_web/Buckets_GridNeighborDesign.md` — Designs bucket/grid neighbor-lookup structures for the web molecular GUI.
- `doc/Topics/molgui_web/Crystal_Generator.md` — Documents crystal-structure generation helpers in the web GUI context.
- `doc/Topics/molgui_web/Editor_for_Polymers_on_Surfaces_and_Edges.md` — Describes editor interaction patterns for polymers on surfaces and edge sites.
- `doc/Topics/molgui_web/Molecular_topology_editor_with_explicit_rings.md` — Discusses explicit ring support in the molecular topology editor.
- `doc/Topics/molgui_web/Molecular_topology_js_plan.md` — JavaScript-centric implementation plan for molecular topology editing in the browser.
- `doc/Topics/molgui_web/Step_Edge_Editor_ideas.md` — Brainstorming notes for a step-edge molecular editor UI.
- `doc/Topics/molgui_web/molgui_web_dev.md` — General web GUI development notes and pointers for the molecular editor frontend.

## ForceFields_Audit
- `doc/topical_audit/forcefields_audit_plan.md` — Phased plan for auditing classical forcefield implementations across CPU C++, OpenCL, PyOpenCL, and WebGPU/WebGL.
- `doc/topical_audit/Forcefields_Audit.md` — Consolidated single-file audit of intramolecular, non-bonding, and molecule-surface force fields (legacy scratchpad).
- `doc/topical_audit/forcefields_overview.md` — High-level taxonomy mapping force field classes to implementation files and status.
- `doc/topical_audit/intramolecular_forcefields.md` — Detailed audit of UFF, MMFFsp3, ProjectiveDynamics, XPBD, and RigidBody force fields.
- `doc/topical_audit/nonbonding_forcefields.md` — Detailed audit of NBFF, exclusion schemes, and Fast Multipole Method.
- `doc/topical_audit/surface_interactions.md` — Detailed audit of GridFF, FoldedAtomicFunctions, Surface.cl, and Ewald2D electrostatics.
- `doc/topical_audit/forcefields_web_implementation.md` — Audit of WebGL/WebGPU shader-based force field ports.

## ForceFields_CPP
- `cpp/common/molecular/UFF.h` — Main UFF class with bond, angle, dihedral, inversion data structures, force assembly, and neighbor lists.
- `cpp/common/molecular/MMFFsp3_loc.h` — MMFFsp3 localized data layout with π–π and π–σ stiffness, pipos orientation, and compact neighbor lists.
- `cpp/common/molecular/MolWorld_sp3.h` — Central orchestrator that instantiates and routes to UFF, MMFFsp3, RigidBodyFF, and GridFF.
- `cpp/common/molecular/RigidBodyFF.h` — Rigid-body dynamics with pose (pos+quat) updates, angular velocity integration, and GridFF coupling.
- `cpp/common/molecular/EwaldGrid.h` — Ewald grid construction and Laplace solving for long-range electrostatics on periodic substrates.
- `cpp/common/molecular/NBFF_SR.h` — Short-range non-bonded variant with repulsive R⁻⁴ potential and AABB collision acceleration.
- `cpp/common/molecular/NBFF_old.h` — Legacy NBFF implementation superseded by `NBFF.h`.
- `cpp/common/math/ProjectiveDynamics_d.h` — ProjectiveDynamics solver with Jacobi, Gauss-Seidel, Cholesky, CG, and momentum-accelerated iterations.
- `cpp/common/math/Multipoles.h` — C++ multipole math: charge projection, energy/force evaluation, and center-of-charge computation.
- `cpp/libs_OCL/OCL_GridFF.cpp` — OpenCL wrapper for grid construction and sampling.
- `cpp/common/molecular/MMFFBuilderBase.h` — Base topology builder with `Atom`, `Bond`, `Angle`, `AtomConf` structures; `autoBonds()`, `insertAtom()`, `insertBond()`.
- `cpp/common/molecular/MMFFBuilder.h` — Advanced topology editing: `assignSp3Type()`, `assignPiFragments()`, `assignBondParams()`, `splitGraphs()`, `splitByBond()`.
- `cpp/common/dataStructures/LimitedGraph.h` — Generic fixed-degree graph template with Tarjan DFS bridge-finding (`bridge()`, `bridgeUtil()`).
- `cpp/common/molecular/MMFFparams.h` — Force field parameter database: `ElementType`, `AtomType`, `BondType`, `AngleType`; loads `.dat` files; `assignSubTypes()` hierarchical typing.

## ForceFields_OpenCL
- `cpp/common_resources/cl/relax_multi.cl` — Unified OpenCL kernel with evalBond, evalAngCos, evalPiAling, getLJQH, getMorsePLQH, and multi-system force assembly.
- `cpp/common_resources/cl/GridFF.cl` — B-spline interpolation kernels: basis(), dbasis(), fe3d_pbc(), sample3D(), sample3D_grid().
- `cpp/common_resources/cl/Surface.cl` — Unified surface kernel: getSurfMorse, getSurfFolded, getSurfFolded_workgroup, Ewald2D, and macroscopic layer corrections.
- `cpp/common_resources/cl/FMM.cl` — Fast Multipole Method tile-based kernel with monopole-dipole-quadrupole force functions.
- `cpp/common_resources/cl/Rigid.cl` — Quaternion integration, B-spline field sampling, and PBC index helpers for rigid-body dynamics.
- `pyBall/pyTruss/truss.cl` — Sparse Jacobi and Gauss-Seidel iterative solvers for truss/projective-dynamics systems.

## ForceFields_Python
- `pyBall/OCL/UFF.py` — PyOpenCL wrapper for UFF GPU force evaluation.
- `pyBall/OCL/UFFbuilder.py` — Topology builder for UFF PyOpenCL simulations.
- `pyBall/OCL/MMFF.py` — PyOpenCL MMFF implementation with `toMMFFsp3_loc()` topology→force-field conversion, `initAtomProperties()` type assignment, and back-neighbor building.
- `pyBall/OCL/GridFF.py` — PyOpenCL grid sampler: GridFF_cl with sample3D, sample3D_comb, and buffer management.
- `pyBall/OCL/GridFF_new.py` — Refactored GridFF wrapper with improved buffer management.
- `pyBall/OCL/GridFFRelaxedScan.py` — Relaxed potential energy surface (PES) scanning routines for GridFF.
- `pyBall/OCL/SurfaceEwald.py` — SurfaceEwaldCL PyOpenCL wrapper for 2D Ewald electrostatics.
- `pyBall/OCL/Surface_utils.py` — GridFF I/O, alignment verification, FDBM fitting, and electrostatics comparison utilities.
- `pyBall/OCL/RigidBodyDynamics.py` — PyOpenCL rigid-body wrapper with REQ→PLQ conversion and GridFF initialization.
- `pyBall/OCL/Assembly.py` / `pyBall/OCL/cl/Assembly.cl` — Rigid-body packing evaluation and clash detection kernels.
- `pyBall/XPBD_2D/XPBD_2D.py` / `pyBall/XPBD_2D/XPBD_2D.cl` — 2D position-based dynamics using complex numbers for rotation.
- `pyBall/XPDB_AVBD/XPDB.py` / `pyBall/XPDB_AVBD/XPDB.cl` — 3D XPBD with angular-velocity-based quaternion integration.
- `pyBall/XPDB_AVBD/RRsp3.py` / `pyBall/XPDB_AVBD/RRsp3.cl` — Rigid-atom sp³ solver with port-based topology.
- `pyBall/MMFF_multi.py` — Multi-system UFF Python interface.
- `pyBall/Ewald2D.py` — Python reference implementation of 2D Ewald summation for charged slabs.
- `pyBall/pyTruss/truss_solver_ocl.py` — OpenCL truss/projective-dynamics solver tests.
- `pyBall/atomicUtils.py` — Topology utilities: bond/angle/dihedral detection, neighbor finding, graph preprocessing, rotation matrices, clustering.
- `pyBall/AtomicSystem.py` — Array-based atomic system with `.mol`/`.mol2`/`.xyz` I/O, `findBonds()`, `neighs()`, `grow_selection()`, `select_all_connected()`, PBC cloning.
- `pyBall/AtomicGraph.py` — Object-graph representation with stable-ID `Atom`, `Bond`, `Ring` objects; `detect_rings()` DFS cycle detection; `add_atom()`, `add_bond()`, `to_arrays()`.
- `pyBall/GUI/MoleculeEditor2D.py` — Matplotlib-based 2D molecule editor (deprecated; superseded by KekuleExplorerGUI).

## ProjectiveDynamics
- `doc/py/ProjectiveDynamics/projective_dynamics.py` — Reference Python implementation of ProjectiveDynamics using dense NumPy solvers.
- `doc/py/ProjectiveDynamics/projective_dynamics_iterative.py` — Iterative solver variant of ProjectiveDynamics.
- `doc/py/ProjectiveDynamics/example_pd.py` — Basic truss dynamics demo.
- `doc/py/ProjectiveDynamics/truss.py` — Truss structure builder and solver.
- `doc/py/ProjectiveDynamics/test_Jacobi_Chebyshev_convergence.py` — Convergence tests for Jacobi and Chebyshev-accelerated solvers.

## ForceFields_Tests
- `tests/tUFF/test_UFF_multi.py` — CPU vs GPU parity for multi-system UFF evaluation.
- `tests/tUFF/test_UFF_ocl.py` — OpenCL-specific UFF tests.
- `tests/tUFF/run_parity_uff.sh` — Automated UFF CPU/GPU parity suite.
- `tests/tMMFF/test_MMFF_ocl_parity.py` — CPU vs GPU parity for MMFF.
- `tests/tMMFF/run_test_GridFF.py` — Basic GridFF construction, sampling, and direct-sum comparison.
- `tests/tMMFF/run_test_GridFF_CaF2.py` — CaF₂(111) surface GridFF validation against known adsorption sites.
- `tests/tMMFF/run_test_GridFF_gauss_smear.py` — Gaussian charge smearing for Coulomb grid convergence.
- `tests/tMMFF/run_test_GridFF_ocl.py` / `tests/tMMFF/run_test_GridFF_ocl_new.py` — OpenCL vs CPU GridFF parity.
- `tests/tMMFF/test_electrostatics_comparison.py` — Compares GridFF, Ewald2D, and brute-force electrostatics on NaCl surfaces.
- `tests/tMMFF/test_gridff_alignment.py` — GridFF origin-convention detection and alignment verification.
- `tests/tMMFF/test_folded_fit_nacl1x1.py` — FoldedAtomicFunctions fit on NaCl(001) with Coulomb reference.
- `tests/tMMFF/gen_gridff_nacl_gpu.py` — GPU-accelerated GridFF generation for NaCl with DFT electrostatics.
- `tests/tMMFF/test_fdbm_fit_dft.py` — FDBM linear fitting against DFT reference adsorption energies.
- `tests/tMMFF/test_fdbm_fit_gridff_mock.py` — Mock FDBM fitting with synthetic data.
- `tests/tMMFF/gui_fdbm_fit.py` — Interactive PyQt5 GUI for live FDBM parameter tuning.
- `tests/tMMFF/GridFF_CaF2_doc_tutorial.md` — Step-by-step GridFF tutorial for CaF₂ surfaces.

## FoldedAtomicFunctions
- `doc/py/FoldedAtomicFunctions/FoldedAtomicFunction.md` — Design document for the plane-wave × exponential basis surface potential.
- `doc/py/FoldedAtomicFunctions/potentials.py` — Morse potential generation and 2D FAF potential profile extraction.
- `doc/py/FoldedAtomicFunctions/faf_func.py` — Minimal functional re-implementation: makePotentialXZ, scan_Qs, scan_basis_a0.
- `doc/py/FoldedAtomicFunctions/tutorial_folded_basis_pareto.md` — Pareto analysis tutorial for basis size vs. accuracy.
- `doc/py/FoldedAtomicFunctions/optimize_z_basis.md` — Optimization of z-basis exponential decay parameters.

## FMM
- `doc/FMM/FMM.md` — Detailed mathematical derivation of energy blending, force distribution, and switching functions for tile-based FMM.

## ForceFields_Web
- `web/molgui_web/js/ProjectiveDynamics.js` — JavaScript port of ProjectiveDynamics solver for browser-based demos.
- `web/molgui_webgpu/Draw3D_webgpu.js` — WebGPU rendering pipeline for atoms and bonds.
- `web/molgui_webgpu/CrystalUtils.js` — Lattice and PBC handling for periodic substrates in WebGPU.
- `web/common_js/Buckets.js` — Generic `Bucket`/`BucketGraph` spatial partitioning with AABB bounds, `aabbOverlap3D()`, `dist2ToAabb()`, `findOverlapNeighbors()`, `exportFlat()`.
- `web/common_js/BucketGrid3D.js` — Uniform 3D grid partitioning for crystal cells: `buildCrystalCellBucketsFromMol()`, `buildWireframeCellVerts()`, `buildWireframeAABBVerts()`.
- `web/common_js/Buckets_SoA.js` — Structure-of-Arrays spatial hash variant.
- `web/common_js/BucketAABBs.js` — AABB construction and overlap tests for rigid-body packing.
- `web/common_js/Selection.js` — Generic index selection with `add/remove/toggle`, tombstones, `SelectionBanks`, and predicate-based selection.
- `web/common_js/MeshBuilder.js` — Mesh construction with vertex/edge de-duplication, `SelectionBanks` integration, SDF-based vertex/edge selection.
- `web/molgui_webgpu/EditableMolecule.js` — Molecular topology editor with `Atom`, `Bond`, `Fragment`, `Bounds`; `id2atom` Map, `topoVersion` caching, VSEPR geometry helpers.
- `web/molgui_webgpu/MMFFLTopology.js` — Topology → XPDB packing: `buildMMFFLTopology()`, `computePiOrientations()`, `buildAngleBonds()`, `buildXPDBInputsFromMol()`.
- `web/molgui_webgpu/MMParams.js` — Parameter loader for `tests/tUFF/data_UFF/{ElementTypes,AtomTypes,BondTypes,AngleTypes,DihedralTypes}.dat`; `resolveTypeNameTable()`, `bondCutoff2()`, `bondLengthEstimate()`.
- `web/molgui_webgpu/MoleculeRenderer.js` — WebGPU atom/bond rendering pipeline.
- `web/molgui_webgpu/GUI.js` — WebGPU molecular editor UI components.
- `web/molgui_web/js/EditableMolecule.js` — **Legacy duplicate** of `molgui_webgpu/EditableMolecule.js` (WebGL era).
- `web/molgui_web/js/MMParams.js` — **Legacy** parameter parser (WebGL era; superseded by `molgui_webgpu/MMParams.js`).
