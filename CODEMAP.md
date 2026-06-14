# CODEMAP.md

# FireCore Repository Structure

FireCore is a multi-component computational chemistry and physics repository containing DFT/DFTB electronic structure methods, classical molecular mechanics and MD, GPU/OpenCL/CUDA acceleration, QM/MM coupling, force-field fitting, visualization, and numerical method prototypes.

**The repository is NOT a single build target.** Different subsystems have different workflows and build procedures.

# Main Components

## `/cpp/` — High-performance C/C++ implementation

### `/cpp/common/` — Core numerical and molecular libraries
- **`math/`** — Linear algebra (`Mat3.h`, `Vec3.h`, `Vec2.h`, `Lingebra.h`), optimization (`DynamicOpt.h`, `CG.h`, `lineSearch.h`), splines (`Bspline.h`, `spline_hermite.h`), geometry (`geom2D.h`, `geom3D.h`, `quaternion.h`), integration (`integration.h`), PDE solvers (`SchroedingerGreen1D.h`, `SchroedingerGreen2D.h`), and numerical utilities (`fastmath.h`, `NumT.h`, `approximation.h`).
- **`molecular/`** — Molecular mechanics force fields (`MMFF.h`, `MMFFsp3.h`, `MMFFsp3_loc.h`, `UFF.h`, `eFF.h`, `NBFF.h`, `NBFF_SR.h`, `GridFF.h`, `RigidBodyFF.h`, `Ewald2D.h`), builders (`MMFFBuilder.h`, `MMFFparams.h`, `UFFbuilder.h`), molecular dynamics worlds (`MolWorld_sp3.h`, `MolWorld_sp3_multi.h`, `MolWorld_sp3_cuda.h`), and fitting (`FitREQ.h`, `FitFF.h`).
- **`dataStructures/`** — Spatial containers (`Buckets.h`, `Buckets2D.h`, `Buckets3D.h`, `Grid.h`, `Grid3D.h`, `HashMap2D.h`, `SimplexRuler.h`), graphs (`Graph.h`, `LimitedGraph.h`), and generic utilities (`Table.h`, `Tree.h`, `containers.h`).
- **`OpenCL/`** — OpenCL harness and wrappers for molecular mechanics (`OCL_MM.h`, `OCL_UFF.h`), DFT (`OCL_DFT.h`), and FFT (`OCLfft_errors.h`).

### `/cpp/common_resources/` — Shared runtime assets
- **`mol/`**, **`xyz/`**, **`crystals/`** — Example molecules and crystals in .mol2, .xyz, .cif, .poscar.
- **`.dat` parameter files** — Canonical force field parameters (`ElementTypes.dat`, `AtomTypes.dat`, `BondTypes.dat`, `AngleTypes.dat`, `DihedralTypes.dat`) loaded by C++, Python, and JS parsers.
- **`cl/`** — Canonical OpenCL compute kernels (`relax_multi.cl`, `UFF.cl`, `GridFF.cl`, `Surface.cl`, `eFF.cl`, `FMM.cl`, `Rigid.cl`) compiled by both C++ and PyOpenCL harnesses. Note: `pyBall/*/cl/` contain domain-specific kernels (DFTB/Grid, PME, hubbard, MQCA, Assembly) and are not duplicates.
- **`shaders/`** — GLSL shaders for WebGL rendering (atom, bond, label, MD predictor/corrector, ProjectiveDynamics).
- **Precalculated grids** — Large `.npy` B-spline PLQd grids for substrate surfaces (`NaCl_1x1_L3/`, `NaCl_8x8_L3_ClHole/`, `CaF2_6L_Ni3_rect_nx2_nz1_L2_top/`).

### `/cpp/apps/` — SDL-based standalone applications
- **`MolecularEditor/`** — Interactive 2D/3D molecular editor with force field relaxation.
- **`EFF/`** — Electron force field visualization.

### `/cpp/apps_OCL/` — OpenCL-accelerated applications
- **`MolecularEditorOCL/`** — GPU-accelerated molecular editor.

### `/cpp/libs/`, `/cpp/libs_OCL/`, `/cpp/libs_SDL/` — Compiled libraries
- `libs/Molecular/` — Compiled molecular mechanics libraries.

### `/cpp/common_SDL/` — SDL/OpenGL shared utilities
- `SDL2OGL/` — 2D/3D rendering helpers, shader management, camera controls.

### `/cpp/sketches_SDL/` — Experimental prototypes and sandbox code
- **`Molecular/`** — 148 experimental molecular visualization and simulation sketches.

## `/fortran/` — Original Fireball DFT/DFTB reference implementation
- **`MAIN/`** — SCF drivers and program entry points.
- **`MODULES/`** — Global modules, data structures, and constants.
- **`ASSEMBLERS/`** — Hamiltonian and overlap matrix assembly.
- **`DASSEMBLERS/`** — Density matrix and post-processing assembly.
- **`INTERACTIONS/`** — Two-center interaction integrals (Coulomb, exchange, etc.).
- **`GRID/`** — FFT/grid operations and electrostatics.
- **`NEIGHBORS/`** — Neighbor list construction.
- **`INITIALIZERS/`** — System setup and basis set initialization.
- **`READFILES/`** — File I/O for atomic structures and parameters.
- **`ROTATIONS/`** — Rotation matrices and angular momentum utilities.
- **`MATH/`** — Fortran math utilities.
- **`doc/`** — Fortran-specific documentation.

## `/fortran2/` — Refactored/reorganized Fortran work-in-progress.

## `/pyBall/` — Python bindings, orchestration, and utilities
- **`FireCore.py`** — Fireball/DFTB Python interface.
- **`MMFF.py`**, **`MMFFsp3.py`**, **`MMFF_multi.py`** — Molecular mechanics orchestration.
- **`AtomicSystem.py`**, **`atomicUtils.py`** — Structure manipulation, file I/O, lattice operations.
- **`Forces.py`**, **`Forces_cpp.py`** — Force field interfaces (Python and C++ bindings).
- **`plotUtils.py`**, **`VispyUtils.py`** — Shared plotting and visualization helpers.
- **`buildUtils.py`**, **`gen_makefile.py`**, **`config_utils.py`** — Build system helpers.
- **`eFF.py`**, **`eFF_terms.py`** — Electron force field Python implementation.
- **`FitREQ.py`**, **`FFFit.py`** — Force field parameter fitting.
- **`Kekule.py`**, **`KekuleBackend.py`**, **`KekuleExplorerGUI.py`** — Kekule structure enumeration and GUI.
- **`MolecularPlacer.py`**, **`SequencePlacer.py`** — On-surface molecular assembly.
- **`AFMExtension.py`** — AFM simulation extensions.
- **`FTIR.py`** — Vibrational spectra post-processing.
- **`Ewald2D.py`**, **`SurfaceSampling.py`**, **`SubstrateBuilder.py`** — Surface electrostatics and substrate generation.
- **`pyTruss/`** — Truss/structure solver in Python.
- **`XPBD_2D/`**, **`XPDB_AVBD/`** — XPBD physics experiments.

### `/pyBall/OCL/` — Standalone pyOpenCL implementations
- **`MolecularDynamics.py`** — GPU molecular dynamics with OpenCL.
- **`GridFF.py`**, **`Surface_utils.py`** — GridFF and surface interaction on GPU.
- **`MMFF.py`**, **`UFF.py`** — PyOpenCL MMFF and UFF force evaluation.
- **`cl/`** — Domain-specific OpenCL kernels (hubbard, MQCA, Assembly, PME, covolve).

### `/pyBall/FireballOCL/` — OpenCL-accelerated Fireball DFTB
- **`OCL_Hamiltonian.py`** — OpenCL Hamiltonian assembly.
- **`STM.py`** — STM current and orbital projection.
- **`Grid.py`** — Real-space density grid projection.
- **`STM_utils.py`** — Sparse/dense orbital data packing utilities.
- **`cl/`** — OpenCL kernels for DFTB grid projection and Hamiltonian.

### `/pyBall/DFTB/` — DFTB+ standalone wrapper
- Python bindings for DFTB+ C API.

### `/pyBall/GLCL/`**, **`GLCL2/`** — OpenGL/OpenCL interop experiments.

## `/web/` — WebGL/WebGPU molecular visualization and physics
- **`molgui_webgpu/`** — Active WebGPU molecular editor (compute shaders, XPBD, MMFF). `EditableMolecule.js`, `MMFFLTopology.js`, `MMParams.js`, `MoleculeRenderer.js`.
- **`molgui_web/`** — Legacy WebGL molecular editor (being phased out).
- **`common_js/`** — Shared JS utilities (`MeshBuilder.js`, `Vec3.js`, `SelectionBanks`).
- **`NBody2D_WebGPU/`** — WebGPU N-body physics demo.
- **`ppstm_web/`** — Web-based STM visualization.

## `/tests/` — Primary reference for implemented functionality and expected workflows
Each test directory is a self-contained workspace with its own `run.sh` and `AGENTS.md`.

### Force Field Tests
- **`tests/tUFF/`** — UFF force field CPU/GPU parity validation.
- **`tests/tMMFF/`** — MMFFsp3 force field, GridFF surface potentials, molecule-surface interactions, rigid-body dynamics, AFM/STM simulation pipelines.
- **`tests/tMMFFsp3/`** — MMFFsp3-specific unit tests.
- **`tests/tMMFFmulti/`** — Multi-system GPU MMFF evaluation.
- **`tests/tEFF/`** — eFF (electron force field) and RARFF reactive force field CPU/GPU parity.
- **`tests/tFitFF/`**, **`tests/tFitREQ/`** — Force field parameter fitting validation.
- **`tests/tAFM/`** — AFM simulation tests, FDBM pipeline, rigid-body tip dynamics.

### QM/DFTB Tests
- **`tests/pyFireball/`** — Fireball/DFTB reference calculations, STM orbital projection, Green's function transport, NEB hydrogen transfer, GPU density projection parity.
- **`tests/dftb/`** — DFTB+ standalone program, C API, and Python wrapper tests.
- **`tests/tDFT/`**, **`tests/tDFT_CO/`**, **`tests/tDFT_pentacene/`** — DFT test calculations.
- **`tests/pySCF/`** — PySCF integration tests.
- **`tests/pyocl_dft/`** — PyOpenCL DFT experiments.
- **`tests/Fireball/`** — Fortran Fireball test inputs and outputs.

### Molecular Editor & GUI Tests
- **`tests/tMolGUIapp/`** — SDL molecular editor tests.
- **`tests/tMolGUIapp_multi/`** — Multi-system SDL editor.
- **`tests/tMolGUIapp_QMMM/`**, **`tests/tMolGUIapp_QMMM_multi/`** — QM/MM coupled editor tests.
- **`tests/tKekule/`** — C++ Kekule bond-order optimizer tests.
- **`tests/tKekuleExplorer/`** — Kekule Structure Explorer GUI backend tests.

### Specialized Physics Tests
- **`tests/tEwald2D/`** — 2D Ewald summation tests.
- **`tests/tFF2D/`** — 2D force field tests.
- **`tests/tLattice2D/`** — 2D lattice dynamics.
- **`tests/tSchroedinger1D/`**, **`tests/tSchroedinger2D/`** — Quantum scattering tests.
- **`tests/tSiNCs/`** — Silicon nanocrystal vibration and phonon tests.
- **`tests/tQuadrature/`** — Numerical quadrature tests.
- **`tests/tMQCA/`** — MQCA (molecular quantum cellular automata) tests.
- **`tests/tCUDA/`** — CUDA-specific tests.
- **`tests/tLammpsTrj/`** — LAMMPS trajectory analysis tests.
- **`tests/tQMMM_diacetylene/`** — QM/MM diacetylene reaction tests.
- **`tests/tAttach/`** — Molecular attachment/assembly tests.
- **`tests/tPsi4resp/`** — Psi4 RESP charge fitting tests.
- **`tests/NonBondSampling/`** — Non-bonded interaction sampling tests.
- **`tests/blender/`** — Blender visualization integration tests.
- **`tests/web/`** — Web frontend tests.

## `/doc/` — Technical notes, derivations, experiments, and developer documentation

### `/doc/topical_audit/` — Global topic map
Cross-language, cross-framework audits of scattered implementations. These are the primary navigation aid for agents working on force fields, topology, surface interactions, etc.
- **`topical_audit.md`** — Index and guide to all topical audit documents.
- **`forcefields_overview.md`** — High-level taxonomy of all force fields.
- **`intramolecular_forcefields.md`** — Bonded force fields (UFF, MMFFsp3, XPBD, ProjectiveDynamics).
- **`nonbonding_forcefields.md`** — Non-bonded interactions (NBFF, FMM, Coulomb, LJ).
- **`surface_interactions.md`** — Molecule-surface interactions (GridFF, FAF, Surface.cl).
- **`RigidSurfPotential_GridFF.md`** — GridFF and rigid surface potential systems.
- **`forcefields_web_implementation.md`** — WebGL/WebGPU force field kernels.
- **`molecular_topology.md`** — Molecular graph representations, bond finding, ring detection.
- **`molecular_topology_types.md`** — Atom type assignment and hybridization detection.
- **`molecular_topology_editors.md`** — Interactive molecular editors across languages.
- **`Htransfer_Kekule_DFTB.md`** — H-transfer NEB and Kekule structure calculations.
- **`Nanocrystal_Vibrations.md`** — Nanocrystal structure generation and vibrational properties.
- **`file_descriptions.md`** — Detailed file-to-topic mapping across the entire repository.

### `/doc/DevNotes/` — Deep developer notes and design discussions
- Agentic documentation design, GPU/CPU parity debugging, UFF/MMFF OpenCL plans, MD implementation notes, FitREQ tutorials, rotational dynamics.

### `/doc/Topics/` — Topic-specific deep documentation
- **`AFM/`** — AFM simulation theory and methods.
- **`ForceFields/`** — Force field derivations and theory.
- **`Kekule_Topology/`** — Kekule structure enumeration and optimization.
- **`OnSurfaceAssembly/`** — On-surface molecular assembly.
- **`RigidBodyAssembly/`** — Rigid body assembly algorithms.
- **`STM/`** — Scanning tunneling microscopy theory.
- **`ManyBody/`** — Many-body physics notes.
- **`PathDiffusion/`** — Path integral and diffusion methods.
- **`TrussSolver/`** — Truss solver documentation.
- **`MoleculeBuilder.md`** — Molecular builder design.

### `/doc/Markdown/` — Tutorials and how-to guides
### `/doc/Julia/` — Julia language experiments and prototypes
### `/doc/Maxima/` — Computer algebra derivations
### `/doc/py/` — Python prototype scripts and experiments
### `/doc/FMM/` — Fast Multipole Method documentation

# Build & Execution

## Critical Rule

Do NOT compile manually using ad-hoc `make` commands.

Use provided:
- `run.sh`
- `make.sh`
- project test scripts

These scripts configure paths, rebuild dependencies, launch correct binaries, and manage runtime resources.

# Typical Workflows

## Run molecular mechanics tests
```bash
cd tests/tMMFF
./run.sh
```

## Run UFF parity tests
```bash
cd tests/tUFF
./run_parity_uff.sh
```

## Run DFTB reference calculations
```bash
cd tests/dftb
./run.sh
```