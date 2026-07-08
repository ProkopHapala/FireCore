# Force Fields Overview

## Taxonomy of Force Field Implementations in FireCore

FireCore contains a broad spectrum of force field implementations spanning **intramolecular**, **non-bonding**, and **molecule-surface** interactions. These are realized across C++, Python, OpenCL, and JavaScript/WebGPU, reflecting both production-grade solvers and experimental research prototypes.

This document provides the high-level taxonomy. Detailed discussions of physics, implementation files, and performance considerations are delegated to:

**Windsurf Codemaps:**
- [FireCore Force Field Navigation: From Audit Docs to Implementation](https://windsurf.com/codemaps/d550a435-7c6f-47b1-aeb5-efc9b098564f-fe86ab10a43f3d18) — Interactive trace from audit documents to source code (UFF, GridFF, Ewald2D, MolWorld orchestration, REQ→PLQ, test parity).

- **`intramolecular_forcefields.md`** — Bonds, angles, dihedrals, torsions (UFF, MMFFsp3, ProjectiveDynamics, XPBD, RigidBody).
- **`nonbonding_forcefields.md`** — Lennard-Jones, Morse, Coulomb, exclusion schemes, and Fast Multipole Method (FMM).
- **`surface_interactions.md`** — GridFF (grid-based interpolation) and FoldedAtomicFunctions (plane-wave basis).
- **`forcefields_web_implementation.md`** — WebGL/WebGPU shader implementations for browser-based visualization and simulation.

---

**Related Codemaps (full list in topical_audit.md):**
- [FireCore Classical Forcefields: MMFFsp3 & UFF](https://windsurf.com/codemaps/53f2fe2c-ac5c-4c0b-b905-af6653adde97-8796fe608a7d71c1)
- [MMFF/UFF CPU vs GPU Testing](https://windsurf.com/codemaps/8d1b056f-1502-4363-b52d-8257de4be453-8796fe608a7d71c1)
- [Surface Potential Evaluation: GridFF B-spline](https://windsurf.com/codemaps/2a639fae-c9cb-407a-9d45-7b806c90c749-8796fe608a7d71c1)
- [FoldedAtomicFunctions](https://windsurf.com/codemaps/c9fc44a7-57a2-47c5-906f-886fa301ccc7-8796fe608a7d71c1)
- [AFM PyOpenCL System: Morse/LJ and FDBM](https://windsurf.com/codemaps/9bb4c2a5-0c38-4943-abe9-254cfdcc75af-8796fe608a7d71c1)
- [XPBD Molecular Dynamics pyOpenCL](https://windsurf.com/codemaps/2e558e51-fdbe-4bd4-8732-7818724d4ced-8796fe608a7d71c1)
- [Rigid Body Dynamics on Surfaces](https://windsurf.com/codemaps/b5d9c2d2-50f0-4ba7-bc65-60db6e06e423-8796fe608a7d71c1)

---

## 1. Intramolecular Force Fields

These describe covalent bonding geometry: bond lengths, angles, and (optionally) dihedrals / torsions.

| Class | Method | Key Files | Status |
|-------|--------|-----------|--------|
| **UFF** | Traditional topology-based force field with bonds, angles, dihedrals, inversions. | `cpp/common/molecular/UFF.h`, `relax_multi.cl` | Production |
| **MMFFsp3** | GPU-optimized variant replacing 4-body terms with π–π and π–σ alignment terms motivated by quantum chemistry of sp-hybridization. | `cpp/common/molecular/MMFFsp3_loc.h`, `pyBall/OCL/MMFF.py` | Production |
| **ProjectiveDynamics** | Position-based implicit solver for stiff harmonic bonds using Jacobi/Gauss-Seidel iterations. | `cpp/common/math/ProjectiveDynamics_d.h`, `pyBall/pyTruss/truss.cl` | Experimental |
| **XPBD (2D/3D)** | eXtended Position-Based Dynamics with per-atom quaternions and port-based constraints. | `pyBall/XPBD_2D/`, `pyBall/XPDB_AVBD/` | Experimental |
| **RigidBodyFF** | Whole-molecule rigid body dynamics (6/7 DOFs) for large time steps and robust assembly/AFM/STM. | `cpp/common/molecular/RigidBodyFF.h`, `pyBall/OCL/RigidBodyDynamics.py` | Production |
| **ReactiveFF** | RARFF, EFF — semi-quantum force fields using localized orbital ansätze. | `cpp/common/molecular/RARFF.h`, `eFF.h` | Highly experimental |

### Design Philosophy
- **GPU Parallelism**: MMFFsp3 and UFF use a *force-assembly* scheme to avoid atomic writes: forces are first emitted into auxiliary buffers, then assembled into the global force array.
- **Node Atoms vs. Capping Atoms**: MMFFsp3 distinguishes "node" atoms (C, N, O with multiple neighbors) from "capping" atoms (H, lone pairs). Only node atoms carry angular terms; capping atoms feel only recoil from their host.
- **Pi-Orbital Alignment**: MMFFsp3 replaces expensive 4-body dihedral terms with cheaper 2-body π–π (`k_pp`) and π–σ (`k_sp`) alignment terms that capture conjugation and orthogonality effects.

---

## 2. Non-Bonding Force Fields

These describe van der Waals, electrostatic, and short-range repulsive interactions between non-bonded atoms.

| Class | Method | Key Files | Status |
|-------|--------|-----------|--------|
| **NBFF** | Base class for LJ (12-6), Morse, Coulomb, and H-bond pseudo-charge interactions. | `cpp/common/molecular/NBFF.h`, `relax_multi.cl` | Production |
| **NBFF_SR** | Short-range variant with AABB collision acceleration for steric hindrance in rigid-body assembly. | `cpp/common/molecular/NBFF_SR.h` | Production |
| **FMM** | Fast Multipole Method (tile-based, single-layer) for long-range electrostatics on GPU. | `cpp/common_resources/cl/FMM.cl`, `cpp/common/math/Multipoles.h` | Experimental |

### Design Philosophy
- **Pairwise Potentials**: `getLJQH()` (`relax_multi.cl:164`) combines Lennard-Jones, Coulomb, and H-bond terms in a single inline function.
- **Exclusion Schemes**: Bonded neighbors (1-2, 1-3) are excluded via explicit exclusion lists (`excl[ natoms*EXCL_MAX ]`) or by subtracting non-bonded contributions from bonded pairs.
- **PLQ Factorization**: For grid evaluation, REQ parameters are converted to PLQ (Pauli, London, Charge) form to enable fast linear combination of precomputed grids.

---

## 3. Molecule-Surface Interactions

Specialized treatment of rigid substrates as an effective medium, avoiding explicit substrate dynamics.

| Class | Method | Key Files | Status |
|-------|--------|-----------|--------|
| **GridFF** | 3D precomputed grids (Pauli, London, Coulomb) sampled via B-spline / trilinear interpolation. | `cpp/common/molecular/GridFF.h`, `cpp/common_resources/cl/GridFF.cl` | Production |
| **FoldedAtomicFunctions** | Compact plane-wave × exponential basis for 2D periodic potentials. | `doc/py/FoldedAtomicFunctions/`, `FoldedAtomicFunction.md` | Experimental |

### Design Philosophy
- **Grid Interpolation**: GridFF evaluates the substrate potential at O(1) per atom via B-spline interpolation, versus O(N_substrate) for direct pairwise summation.
- **Long-Range Electrostatics**: Coulomb grids typically require Poisson solving (FFT or Ewald summation) to converge the long-range tail.
- **Folded Basis**: FAF uses `cos(kx)·exp(-a·z)` basis functions, reducing storage from millions of grid points to a few hundred coefficients—ideal for GPU local memory.

---

## Cross-Cutting Concerns

### Periodic Boundary Conditions (PBC)
Most force fields support PBC via precomputed shift vectors (`shifts[npbc]`). The number of periodic images is automatically determined from cell dimensions (`autoNPBC()` in `GridFF.h`).

### Multi-System GPU Evaluation
The `relax_multi.cl` kernel uses a 2D NDRange where `global.y = nSystems`, allowing simultaneous relaxation of many molecular copies (e.g., for sampling or machine learning).

### REQ → PLQ Conversion
For grid-based evaluation, the standard REQ (Radius, Energy, Charge, H-bond) parameters are converted to PLQ (Pauli strength, London strength, Charge) form:

```cpp
PLQ.x = sqrt(REQ.x * REQ.y);  // Pauli
PLQ.y = REQ.x^6 * REQ.y;      // London
PLQ.z = REQ.z;                // Charge
```

This factorization allows the same grid to be reused with different atom types via simple scaling.

---

## References

- [Topical Audit Index](topical_audit.md) — priority ranking, dependency graph, missing topics
- `Forcefields_Audit.md` — Consolidated audit with exhaustive file listings (legacy).
- `intramolecular_forcefields.md` — Detailed physics and implementation of bonded terms.
- `nonbonding_forcefields.md` — Detailed physics and implementation of non-bonded terms.
- `surface_interactions.md` — Detailed physics and implementation of substrate interactions.
- `forcefields_web_implementation.md` — WebGL/WebGPU shader details.
- [Molecular Topology](molecular_topology.md) — topology and graph representations underlying FF evaluation

---

*Last updated: 2026-06-13*
