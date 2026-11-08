# Hermes Topical Overview

High-level problem map for FireCore, prepared from existing docs and README files.
This file captures what the project is trying to solve, not every implementation detail.

## Dominant scientific problems

1. Non-covalent self-assembly of molecules on inorganic crystal surfaces.
2. Chemical reaction pathways between organic molecules on surfaces.
3. Computational prototyping of molecular designs for programmable self-assembling.
4. High-resolution AFM / STM simulation of organic molecules.

## Corresponding computational systems

- **Surface interaction models** — molecule-substrate energetics via GridFF and rigid scanning.
- **Assembly simulation** — pose generation, multi-molecule packing, and relaxed surface configurations.
- **Force field representations** — classical MMFF/UFF plus folded basis alternatives and FitREQ parameter fitting.
- **AFM simulation** — classic probe particle model plus full density-based model based on DFT density.
- **STM transport modeling** — density projection, orbital coupling, and transport via master-equation or Green-function approaches.
- **QM/MM coupling** — Fireball DFTB as the quantum region interacting with classical force fields.
- **Geometry optimization / dynamics** — GPU-accelerated relaxation, rigid-body dynamics, surface MD.

## Experimental linkage

- AFM: df images and force spectroscopy.
- STM: orbital-resolved tunneling and light-assisted STM variants.
- Manipulation: path optimization on surfaces, tethered relaxation, and tip-induced motion.

## Cross-cutting code concerns

- Multi-language stack: Fortran, C++, Python, OpenCL, WebGPU, JavaScript.
- Multiple acceleration paths: CPU, OpenMP, OpenCL, CUDA.
- GUI, CLI, web, and test entry points for the same physics kernels.
- Ongoing consolidation challenge across student branches.
