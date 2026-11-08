FireCore Consolidated Topic Map
Purpose: topic-first index of what this repo contains and where relevant implementation files live.
Sources: doc/Topics/*, doc/DevNotes/*, doc/Markdown/*, doc/*.md, top-level AGENTS.md and README.

Top-level references
- AGENTS.md — coding rules, GPU rules, test protocol, performance rules, repo conventions
- README.md — top-level directory map, modules, examples
- Project_Structure.md — directory map and roles
- PATH_CONFIGURATION.md — config layout and path conventions

AFM simulation
Core docs
- doc/Topics/AFM/AFM.md
- doc/Topics/AFM/AFM_migration_plan.md
- doc/Topics/AFM/AFM_migration_progress.md
- doc/Topics/AFM/AFM_FDBM_DFTB.md
- doc/Topics/AFM/AFM_FDBM_fitting.md
- doc/Topics/AFM/AFM_FDBM_optimization.md
- doc/Topics/AFM/AFM_FDBM_profiling_optimization.md
- doc/Topics/AFM/AFM_FDBM_pySCF.md
- doc/Topics/AFM/AFM_migration_discusion.md
Key code
- pyBall/DFT/jobs.py
- pyBall/DFT/high_level.py
- pyBall/DFT/oclfft.py
- cpp/libs_OCL/OCL_GridFF.cpp
- cpp/common/OpenCL/OCL_DFT.h
- cpp/common_resources/cl/myprog.cl
- cpp/common_resources/cl/GridFF.cl
- cpp/common_resources/cl/relax.cl
- pyBall/OCL/GridFF.py
- pyBall/FireballOCL/Grid.py
- pyBall/FireballOCL/cl/Grid.cl
- pyBall/AtomicSystem.py
- pyBall/OCL/MMparams.py
- cpp/common/molecular/QEq.h
- cpp/common_resources/ElementTypes.dat
- fortran/MAIN/libFireCore.f90
- tests/tDFT*
- tests/tDFT_CO*
- tests/tDFT_pentacene*
- tests/tAFM*

STM simulation
Core docs
- doc/Topics/STM/STM_GF_new.md
- doc/Topics/STM/STM_GPU_QMMM.md
Key code
- pyBall/OCL/GridFF.py
- pyBall/FireballOCL/Grid.py
- cpp/common_resources/cl/myprog.cl
- tests/pyFireball/test_grid_projection.py

Surface interactions and scanning
Core docs
- doc/Topics/OnSurfaceAssembly/Comprehensive_Overview_Molecule_Surface_Systems.md
- doc/Topics/OnSurfaceAssembly/GridFF_RelaxedScan_cpp_notes.md
- doc/Topics/OnSurfaceAssembly/GridFF_atoms_alignment.md
- doc/Topics/OnSurfaceAssembly/GridFF_FDBM_Fitting.md
- doc/Topics/OnSurfaceAssembly/SurfacePotentialVisualization_discussion.md
- doc/Topics/OnSurfaceAssembly/SurfacePotentialVisualization_report.md
- doc/Topics/OnSurfaceAssembly/ExplorerVisPy_GridFF_Interactive_Plan.md
- doc/Topics/OnSurfaceAssembly/FoldedSubstratePotential_OpenCL.md
- doc/Topics/OnSurfaceAssembly/RigidBodyAFM.md
- doc/Topics/OnSurfaceAssembly/ElectrostaticContinuumEmbeding_report.md
- doc/Topics/OnSurfaceAssembly/Ewald_2D.md
- doc/Topics/topical_audit/RigidSurfPotential_GridFF.md
Key code
- pyBall/OCL/InteractionEnergy.py
- pyBall/OCL/MolecularDynamics.py
- pyBall/OCL/Assembly.py
- pyBall/SurfaceSampling.py
- pyBall/SequencePlacerVisPy.py
- pyBall/MolecularPlacerVisPy.py
- pyBall/MolGUI_common.py
- pyBall/VispyUtils.py
- cpp/common_resources/cl/GridFF.cl
- cpp/common_resources/cl/relax.cl
- cpp/common_resources/cl/Assembly.cl
- tests/tMMFF/*

Force fields
Core docs
- doc/Topics/ForceFields/UFF_multiSystem.md
- doc/Topics/ForceFields/UFF_multiSystem_summary.md
- doc/Topics/ForceFields/GridFF_fit_CG.md
- doc/Topics/ForceFields/GridFF_fit_jacobi.md
- doc/Topics/ForceFields/NeighborExclusion_Summary.md
- doc/DevNotes/MMFF/* and related fit/FF docs
- doc/DevNotes/FitREQ*.md
Key code
- pyBall/OCL/UFF.py
- pyBall/OCL/MMFFL.py
- cpp/common/molecular/FitREQ.h
- cpp/libs/Molecular/FitREQ_lib.cpp
- pyBall/FitREQ.py
- cpp/common_resources/cl/relax_multi.cl
- cpp/common_resources/cl/relax_multi_mini.cl

Rigid body assembly
Core docs
- doc/Topics/RigidBodyAssembly/RigidBodyAssemblyDiscussion.md
- doc/Topics/RigidBodyAssembly/RigidBodyAssemblyDiscussion_tutorial.md
- doc/Topics/RigidBodyAssembly/RigidGridFF_H2O_NaCl_report.md
Key code
- pyBall/OCL/Assembly.py
- cpp/common_resources/cl/Assembly.cl
- tests/tMMFF/test_assembly.py

DFTB and QM integration
Core docs
- doc/DFTB_API_migration.md
- doc/DevNotes/DFTB_consolidation_notes.md
- doc/QM_DFT_DFTB_packages_for_integration*.md
Key code
- pyBall/dftb_utils.py
- pyBall/DFT/oclfft.py
- pyBall/DFT/jobs.py
- pyBall/FireballOCL/Grid.py
- pyBall/FireballOCL/cl/Grid.cl
- fortran/MAIN/libFireCore.f90
- tests/tDFT*

Molecular dynamics
Core docs
- doc/DevNotes/MD_MMFF_OCL_notes.md
- doc/DevNotes/MD_MMFF_OCL_parity_discussion.md
- doc/DevNotes/rotational_dynamics_implementation.md
- doc/Topics/RigidBodyAssembly/*
Key code
- pyBall/OCL/MolecularDynamics.py
- cpp/common_resources/cl/relax_multi.cl
- cpp/common_resources/cl/relax_multi_mini.cl

Path optimization and manipulation
Core docs
- doc/Topics/PathDiffusion/PathDiffusion_design.md
- doc/Topics/PathDiffusion/PathDiffusion_discussion.md
- doc/Topics/PathDiffusion/PathDiffusion_discussion_old.md
- doc/Topics/PathDiffusion/ManipulationPathOpt_Report.md
Key code
- pyBall/OCL/ManipulationPathOpt.py
- cpp/common_resources/cl/relax_multi.cl

Many-body and quantum transport
Core docs
- doc/Topics/ManyBody/MQCA_GrayCodeSolver.md
- doc/Topics/ManyBody/HubbardSolver.md
- doc/Topics/ManyBody/HubbardSolver_new.md
- doc/Topics/ManyBody/MC_PME_dicussion.md
- doc/Topics/ManyBody/MC_PME_progress.md
Key code
- pyBall/OCL/HubbardSolver.py

Linear algebra and solvers
Core docs
- doc/Topics/LinearLagebra/Parallel_GaussSeidel.md
- doc/Topics/LinearLagebra/Chebyshev_Acceleration.md
- doc/Topics/LinearLagebra/ChebyshevAndInterialAccelerationOfLinearSolvers.md
- doc/Topics/LinearLagebra/NumricalStabilityLinearSolvers.md

Truss and projective dynamics
Core docs
- doc/Topics/TrussSolver/ProjectiveDynamics.md
- doc/Topics/TrussSolver/ProjectiveDynamics_wegl_glsl.md
- doc/Topics/TrussSolver/VertexBlockDescent.md
- doc/Topics/TrussSolver/High_Precision_Physics_on_singlepoint_GPU_hardware.md

Kekule topology and aromatic systems
Core docs
- doc/Topics/Kekule_Topology/Kekule_Topology.md
- doc/Topics/Kekule_Topology/Kekule_Optimizer.md
- doc/Markdown/kekule_research_notes.md
- doc/Markdown/KekuleGUI_Extension_System.md

Molecular builder and substrates
Core docs
- doc/Topics/MoleculeBuilder.md
- doc/Topics/molgui_web/Step_Edge_Editor_ideas.md
- doc/Topics/molgui_web/Editor_for_Polymers_on_Surfaces_and_Edges.md
- doc/Topics/molgui_web/Crystal_Generator.md
- doc/Topics/molgui_web/Molecular_topology_editor_with_explicit_rings.md
- doc/Topics/molgui_web/Molecular_topology_js_plan.md
- doc/Topics/molgui_web/Buckets_GridNeighborDesign.md
- doc/Topics/molgui_web/molgui_web_dev.md
Key code
- web/molgui_web/js/EditableMolecule.js
- web/common_js/Buckets.js
- web/common_js/GridIndexer.js
- pyBall/GUI/*

Web visualization and GUI
Core docs
- doc/Topics/molgui_web/*
- doc/Topics/MolecularGUI/MolViewer_PyQt.md
- doc/Topics/MolecularGUI/MolViewer_PyQt_2.md
- doc/Topics/MolecularGUI/MoleculeEditor2D.md
- doc/Topics/MolecularGUI/OpenGL_Viewer_Documentation.md
- doc/Markdown/GUI.md
Key code
- pyBall/GUI/GLGUI.py
- pyBall/GUI/shaders/*
- pyBall/VispyUtils.py
- pyBall/ExplorerVisPy.py
- pyBall/SequencePlacerVisPy.py
- pyBall/MolecularPlacerVisPy.py
- tests/EFF/xyz_view_new.py

Scientific computing and supporting topics
Core docs
- doc/Topics/FTIR_Nanocrystals/gen_nanocrystals.md
- doc/Topics/FTIR_Nanocrystals/Hessian_fitting.md
- doc/Topics/RandomGPU/NormalDistributionGPU.md
- doc/Topics/FMM/FMM.md
- doc/Markdown/MMFFBuilder.md
- doc/Markdown/eFF_documentation.md
- doc/py/CoarseGrainedFF/GayBerne.md
- doc/py/FoldedAtomicFunctions/*.md
- doc/py/OrderN_QM/*.md
- doc/py/SparseMatrixToGrid/*.md
- doc/py/SpatialPartitioning/*.md
- doc/py/wavelets/*.md
- doc/Maxima/GaussianIntegrals.md
- doc/Maxima/MaximaTutorial.md

Dev notes and active work areas
Core docs
- doc/DevNotes/UFF_ocl_plan*.md
- doc/DevNotes/UFF_muli_plan.md
- doc/DevNotes/FitREQ_ocl.md
- doc/DevNotes/FitREQ_ocl_vispy.md
- doc/DevNotes/FitREQ_CPU_Tutorial.md
- doc/DevNotes/LocalyOrthogonalizeOrbitalForcefield.md
- doc/DevNotes/LLCAO_1D.md
- doc/DevNotes/MolecularEditor_2D_py.md
- doc/DevNotes/GPU-CPU-Debugging-and-OCL-Interface.md
- doc/DevNotes/OpenCLKernelMetaprogramming.md
- doc/DevNotes/ToDo_FastCollision_*.md
- doc/DevNotes/ToDo_GUI.md
- doc/DevNotes/ToDo_Fit_Hbond_GPU.md
- doc/DevNotes/Limit_Language_Server_Memory.md
- doc/DevNotes/MMFFsp3_pi_angular_momentum_problem.md
- doc/DevNotes/MMFF/FFutils_consolidation_summary.md
- doc/DevNotes/MMFF/FFutils_migration_guide.md
- doc/DevNotes/MMFF/FFutils_visual_debugging_consolidation.md
- doc/DevNotes/MMFF/MMFF_variants.md
- doc/DevNotes/MMFF/MMFF_flat_surface_potential.md
- doc/DevNotes/MMFF/MMFF_flatsurf_pyOpenCL_progress.md
- doc/DevNotes/MMFF/pyocl_MMFF/pyocl_MMFF.md
- doc/DevNotes/EFF/*.md

Tests and Python packages adjacent to docs
- tests/tAFM/*
- tests/tDFT*
- tests/tDFT_CO*
- tests/tDFT_pentacene*
- tests/tMMFF/*
- tests/tUFF/*
- tests/tFitREQ/*
- tests/tSiNCs/*
- tests/tEFF/*
- tests/tAssembly/*
- tests/pyFireball/*
- pyBall/OCL/*
- pyBall/DFT/*
- pyBall/GUI/*
- pyBall/FireballOCL/*
- pyBall/*.py
- cpp/common_resources/cl/*
- cpp/common_resources/Substrates/*
- cpp/common_resources/ElementTypes.dat
- fortran/MAIN/*
