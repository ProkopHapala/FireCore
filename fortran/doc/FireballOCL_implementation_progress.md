# FireballOCL Implementation Progress

This document tracks the progress of reimplementing Fireball's Hamiltonian assembly using PyOpenCL. It serves as a technical bridge between the original Fortran implementation and the new GPU-accelerated Python package.

## 1. Project Overview

The goal is to port the core Hamiltonian and Overlap matrix assembly to PyOpenCL for significant performance gains, especially for large molecular systems. The implementation replicates the list-directed data parsing and real-space integral interpolation logic found in Fireball.

### Sub-package: `pyBall.FireballOCL`
- `FdataParser.py`: Parses 2-center (1D splines) and 3-center (2D grids) integrals.
- `OCL_Hamiltonian.py`: Orchestrates PyOpenCL buffers, kernels, and sparse matrix assembly.
- `cl/hamiltonian.cl`: OpenCL kernels for interpolation and orbital rotations.

---

## 2. Hamiltonian Decomposition & Mapping

Fireball decomposes the Hamiltonian $H$ as:
$$H = T + V_{na} + V_{nl} + V_{xc} + V_{eb}$$

| Term | Description | Fortran Subroutine | PyOpenCL Status |
| :--- | :--- | :--- | :--- |
| $S$ | Overlap Matrix | `assemble_2c` | **Verified** (Tol: $10^{-7}$) |
| $T$ | Kinetic Energy | `assemble_2c` | **Verified** (Tol: $10^{-7}$) |
| $V_{na}$ | Neutral Atom Potential | `assemble_2c`, `assemble_3c` | **Mismatch** ($10\times$ difference) |
| $\bar{\rho}$ | Average Density (3C) | `average_ca_rho` | **Implemented** (Needs 3C) |
| $V_{xc}$ | Exchange-Correlation | `assemble_olsxc_off` | Pending |

---

## 3. H2 Verification Results (Summary)

Verification was performed on an H2 molecule ($d=0.74$ Å). Results compared PyOpenCL `OCL_Hamiltonian` against Fortran `libFireCore` sparse export.

| Component | Max Abs Diff | Status | Notes |
| :--- | :--- | :--- | :--- |
| **Overlap S** | $4.86 \times 10^{-8}$ | **PASSED** | Correct rotation and interpolation. |
| **Kinetic T** | $2.39 \times 10^{-7}$ | **PASSED** | Correct rotation and interpolation. |
| **Vna (Total)**| $4.16 \times 10^{2}$ | **FAILED** | Major magnitude and diagonal mismatch. |

---

## 4. Vna Discrepancy Analysis

The discrepancy in $V_{na}$ is attributed to three main factors:

### A. Missing 3-Center Assembly
In the McWEDA/OLSXC scheme, $V_{na}$ is accumulated from both 2-center and 3-center terms.
- **Fortran**: Calls `assemble_3c` which adds massive contributions to the `vna` matrix (e.g., diagonal shifts of ~411 for H2).
- **PyOpenCL**: Currently only implements the 2-center spline interpolation (`vna_atom`, `vna_ontopl`, `vna_ontopr`).

### B. Energy Scaling Factor (`eq2`)
The Fortran implementation uses a global constant `eq2 = 14.39975` (defined in `constants_fireball.f90`) to scale electrostatic potential terms, likely converting to common energy units (eV). PyOpenCL lacks this scaling.

### C. Interaction Mapping
The mapping between `interaction` integer codes and data files is complex:
- `interaction=1`: Overlap (2C) or `bcna` (3C).
- `interaction=4`: $V_{na}$ Atom (2C).
- `interaction=13`: Kinetic (2C).
*Note: In `assemble_2c.f90`, the $V_{na}$ atom block actually re-uses the value of `interaction=13` set by the previous Kinetic block. This behavior must be replicated exactly to match the reference.*

---

## 5. Technical Implementation Details

### Rotation Logic (Slater-Koster)
Replicated logic from `ROTATIONS/epsilon.f90` and `ROTATIONS/twister.f90`.
- **Local Frame**: Bonds are aligned along the Z-axis in the local frame.
- **Orbital Order**: $(s, p_y, p_z, p_x)$ following Fireball convention.
- **Transformations**: Implemented in OpenCL kernels to handle arbitrary bond orientations.

### Interpolation Schemes
- **2-Center**: Cubic spline interpolation (replicating `buildspline_1d.f90`). Spline coefficients $(y, y'')$ are pre-calculated in Python and passed as `float4` to GPU.
- **3-Center**: 2D grid interpolation ($x, y$ - distance to NA and bond length) followed by Legendre polynomial expansion in $\cos \theta$. Replicated from `trescentros.f90` and `interpolate_2d.f90`.

---

## 4. Current Status & Verification

### Completed
- [x] Robust `Fdata` parser handling inhomogeneous multi-line data.
- [x] 1D Spline generator (Python) and GPU interpolation kernel.
- [x] 2-center orbital rotations for $s$ and $p$ orbitals.
- [x] 3-center 2D interpolation kernel.
- [x] Basic verification for H2 (Overlap/Kinetic) and H3 (Neutral Atom 3C).

### Next Steps
- [ ] Implement full **McWEDA / OLSXC** logic for XC potentials.
- [ ] Implement Non-local pseudopotential ($V_{nl}$) assembly.
- [ ] Full integration of `assemble_full` with the main SCF loop.
- [ ] Rigorous benchmarking and verification against `FireCore` sparse export.

---

## 5. Development Notes & References

- **Sparse Matrix**: Uses an atom-neighbor block format `[n_pairs, norb, norb]`.
- **Fireball Convention**: The "Ortega convention" for $p$-orbitals: $p_y, p_z, p_x$ mapped to indices $1, 2, 3$ (0-indexed: 1:y, 2:z, 3:x).
- **Theoretical Basis**:
    - PRB 40, 3979 (1989) - OLSXC / McWEDA.
    - Sankey and Niklewski, PRB 40, 3979 (1989) - Basis functions and integrals.

---
*Created by Antigravity AI on 2025-12-21*
