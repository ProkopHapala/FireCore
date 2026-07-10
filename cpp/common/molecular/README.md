# cpp/common/molecular — Header-only molecular mechanics and force-field algorithms

Header-only C++ library providing force-field evaluation, parameter fitting,
topology building, and molecular structure manipulation. Implementations are
in `.h` files (inline/template); C wrappers for Python are in `cpp/libs/Molecular/`.

## Force-Field Fitting

- **FFfit.h** — Hessian-based FF parameter fitting: Wilson B-matrix, sensitivity matrices (bond/angle/dihedral), linear least-squares + gradient descent with ParamMap symmetry sharing. CSR bond-graph algorithms (bounded BFS, 1-4 pairs, dihedral enumeration) and batch dihedral sensitivity with FD symmetry exploitation.
- **FitREQ.h** — Charge equilibrium (QEq) parameter fitting.
- **FitFF.h** — Generic force-field fitting utilities.

## Force Fields

- **UFF.h** — Universal Force Field: bond, angle, torsion, vdW, electrostatics. Prokop signed-dihedral convention.
- **MMFF.h** / **MMFFsp3.h** — MMFF94 sp3 force field with neighbor lists, Hessian FD.
- **MMFFsp3_loc.h** — Local MMFF variant with OpenCL-compatible flat arrays.
- **MMFFmini.h** — Minimal MMFF for testing.
- **NBFF.h** / **NBFF_SR.h** — Non-bonded force field (short-range variant).
- **eFF.h** — Electron force field (core-electron interactions).
- **RigidBodyFF.h** — Rigid-body force field with quaternion dynamics.
- **RRFF.h** / **RARFF.h** — Reactive and adaptive reactive force fields.
- **SimpleForceField.h** — Minimal FF for unit tests.

## Topology & Structure

- **MolecularGraph.h** — Bond graph representation and traversal.
- **MMFFBuilderBase.h** / **MMFFBuilder.h** — Topology → flat FF array conversion (Builder pattern).
- **UFFbuilder.h** / **UFFbuilder_new_Si.h** — UFF-specific topology builders.
- **Molecule.h** / **MoleculeType.cpp** — Molecule data structures.
- **MolecularConfiguration.h** — Configuration state management.
- **MolecularDatabase.h** — Molecular fragment database.
- **Atoms.h** — Atomic position/type arrays.
- **AtomicConfiguration.h** — Atomic-level configuration.
- **AtomicSystemDescriptors.h** — Descriptors for atomic systems.

## Molecular Dynamics

- **MolWorld_sp3.h** — sp3 molecular world (positions, forces, integrator).
- **MolWorld_sp3_multi.h** — Multi-system variant for parallel replicas.
- **MolWorld_sp3_QMMM.h** — QM/MM hybrid molecular world.
- **CarParrinello.h** — Car-Parrinello molecular dynamics.
- **GlobalOptimizer.h** — Global optimization (basin hopping, simulated annealing).
- **GOpt.h** — Geometry optimization (LBFGS, FIRE).
- **PBCsystem.h** — Periodic boundary condition management.

## Interactions

- **InteractionsGauss.h** — Gaussian interaction potentials.
- **InterpolateTricubic.h** / **InterpolateTrilinear.h** — Grid interpolation.
- **GridFF.h** — Grid-based force field (precomputed potential on grid).
- **Ewald2D.h** / **EwaldGrid.h** — Ewald summation for 2D and grid systems.
- **DipoleMap.h** — Dipole moment mapping.
- **QEq.h** — Charge equilibrium solver.
- **DirectionStiffness.h** — Direction-dependent stiffness tensors.

## Auxiliary

- **Vec3.h** / **Mat3.h** — 3D vector and matrix math (used throughout).
- **AOIntegrals.h** — Atomic orbital integral utilities.
- **AOrotations.h** — AO rotation matrices.
- **Confs.h** — Atom configuration (hybridization, bond order).
- **Groups.h** — Atom grouping and selection.
- **Manipulation.h** — Structure manipulation (rotate, translate, mirror).
- **LatticeMatch2D.h** — 2D lattice matching for surface interfaces.
- **StructuralHistogram.h** — Structural analysis histograms.
- **SVG_render.h** — SVG molecular rendering.
- **molecular.cpp** / **MoleculeType.cpp** — Non-header implementation files.
