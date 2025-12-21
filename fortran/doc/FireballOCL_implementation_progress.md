# FireballOCL Implementation Progress
## Latest status (2025-12-21)

- **2-center (finished)**  
  - PyOpenCL now loads 2c data with the same naming as Fortran `read_2c`: `overlap`, `kinetic`, and all vna variants `vna_atom_<isorp>`, `vna_ontopl_<isorp>`, `vna_ontopr_<isorp>`, iterating `isorp` up to `nssh`.  
  - Generic `vna` requests are satisfied via the `vna_atom_*` tables (matching Fortran interaction 4 usage).  
  - Sign conventions for py/pz rotation in OpenCL kernels match Fortran.  
  - `verify_scan_2c.py` compares Fortran vs PyOpenCL and matches (no NaNs, diffs within threshold).

- **3-center (in progress)**  
  - **Data parity achieved:** Added `firecore_export_bcna_table` (Fortran) and `export_bcna_table` (ctypes) + `verify_Fdata_read_3c.py` plots; Fortran vs Python tables match bitwise for H bcna_01..05 (max|diff|=0).  
  - **Remaining issue:** OpenCL raw scan still mismatches Fortran on L0 (bcna_01) despite identical tables; likely a host/kernel mapping/stride bug (component selection or triplet index). L1/L2 etc. match to ~1e-9.  
  - `verify_scan_3c.py --mode raw` highlights this gap; interpolation kernel logic itself matches Fortran (off-by-one fixed), so focus is on the bcna component selection in the OpenCL path.
  - **Data-loading bugs fixed:**  
    - 3c header parsing in Python skipped wrong lines (charge lines with `<===`), causing wrong reshapes; fixed by matching Fortran `readheader_3c`.  
    - `num_nonzero` for 3c now uses the 3c-specific count (not 2c).  
    - Verification script `tests/pyFireball/verify_Fdata_read_3c.py` visualizes Fortran vs Python grids (5×3 subplot: Fortran/Python/Diff) and prints min/max/diff; confirms tables load identically.

- **Problem root cause fixed**  
  - PyOpenCL was previously missing shell-resolved vna files and inventing ad-hoc names; Fortran is ground truth. Loader now follows Fortran naming exactly; no extra names are introduced.

- **Next steps**  
  - Fix OpenCL bcna L0 mismatch: verify type_idx/theta/i_nz mapping in `scan_3c_raw_points` vs host packing, then re-run `verify_scan_3c.py --mode raw`.  
  - After raw parity, implement Fortran-equivalent 3c rotation/assembly in OpenCL.  
  - Keep verifying with `verify_scan_3c.py`; add batched paths once scalar parity is reached.


This document tracks the progress of reimplementing Fireball's Hamiltonian assembly using PyOpenCL. It serves as a technical bridge between the original Fortran implementation and the new GPU-accelerated Python package.

## 1. Project Overview

The goal is to port the core Hamiltonian and Overlap matrix assembly to PyOpenCL for significant performance gains, especially for large molecular systems. The implementation replicates the list-directed data parsing and real-space integral interpolation logic found in Fireball.

### Sub-package: `pyBall.FireballOCL`
- `FdataParser.py`: Parses 2-center (1D splines) and 3-center (2D grids) integrals.
- `OCL_Hamiltonian.py`: Orchestrates PyOpenCL buffers, kernels, and sparse matrix assembly.
- `cl/hamiltonian.cl`: OpenCL kernels for interpolation and orbital rotations.

## 2. Rules of Engagement

> [!IMPORTANT]
> **No Editing of Reference Fortran Code**
> 1. The original Fortran code in `fortran/` (except for `MAIN/libFireCore.f90`) is our **ground truth reference**. It is well-tested and must not be functionally modified.
> 2. Do not reorder logic, change variable meanings, or redefine interaction/orbital indices.
> 3. Conditional blocks (e.g., `if (ioff_Vna .eq. 1)`) can be added for component export, but core physics/assembly logic is strictly read-only.
> 4. Debug prints are allowed but should be minimal.

## 3. Variable & Index Mapping Reference

### A. General Variables
- `natoms`: Total number of atoms in the system.
- `iatom`, `jatom`: Current atom and its neighbor index (1-based system indices).
- `matom`: Index of the atom in the central cell (unit cell).
- `ineigh`: Index of the neighbor in the local neighbor list of `iatom`.
- `in1`, `in2`, `in3`: Species indices (pointing to the species defined in `info.dat`).
- `xl(3, mbeta)`: Lattice translation vector for image `mbeta`.
- `num_orb(in)`: Total number of orbitals for species `in`.
- `nssh(in)`: Number of shells for species `in`.
- `lssh(issh, in)`: Angular momentum $L$ for shell `issh`.
- `mu(idx, in1, in2)`, `nu(idx, in1, in2)`: Indices in the local $(num\_orb \times num\_orb)$ block for a specific interaction.
- `mvalue(idx, in1, in2)`: Absolute value of $m$ difference ($|m_1 - m_2|$) for 3-center integrals.
- `h_mat`, `s_mat`: Sparse matrix blocks `[n_pairs, norb, norb]`.

### B. Orbital Indexing (Ortega Convention)
Fireball uses a specific ordering for $L > 0$ orbitals within a shell:

| Shell | $L$ | Orbital Sequence (1-indexed) | Quantum Numbers ($m$) |
| :--- | :--- | :--- | :--- |
| **S** | 0 | 1: $s$ | $0$ |
| **P** | 1 | 2: $p_y$, 3: $p_z$, 4: $p_x$ | $-1, 0, 1$ |
| **D** | 2 | 5: $xy$, 6: $yz$, 7: $z^2$, 8: $xz$, 9: $x^2-y^2$ | $-2, -1, 0, 1, 2$ |

*Note: In PyOpenCL (0-indexed), $p_y, p_z, p_x$ map to indices $1, 2, 3$.*

### C. Interaction Codes (`interaction`)
These integers determine the type of integral being calculated/interpolated:

| Code | Interaction | Type |
| :--- | :--- | :--- |
| **1** | Overlap ($S$) / `bcna` (3C) | 2-Center / 3-Center |
| **2** | $V_{na}$ Ontop Left ($V(1)$) | 2-Center |
| **3** | $V_{na}$ Ontop Right ($V(2)$) | 2-Center |
| **4** | $V_{na}$ Atom ($V(3)$) | 2-Center |
| **5** | Non-local Pseudopotential ($V_{nl}$) | 2-Center |
| **13**| Kinetic Energy ($T$) | 2-Center |

### D. Physical Constants
- `eq2 = 14.39975`: Scaling factor used for electrostatic potential terms ($V_{na}$) to maintain consistency in eV/Angstrom units.

## 4. Hamiltonian Decomposition & Mapping

Fireball decomposes the Hamiltonian $H$ as:
$$H = T + V_{na} + V_{nl} + V_{xc} + V_{eb}$$

| Term | Description | Fortran Subroutine | PyOpenCL Status |
| :--- | :--- | :--- | :--- |
| $S$ | Overlap Matrix | `assemble_2c` | **Verified** (Tol: $10^{-7}$) |
| $T$ | Kinetic Energy | `assemble_2c` | **Verified** (Tol: $10^{-7}$) |
| $V_{na}$ | Neutral Atom Potential | `assemble_2c`, `assemble_3c` | **Mismatch** ($10\times$ difference) |
| $\bar{\rho}$ | Average Density (3C) | `average_ca_rho` | **Implemented** (Needs 3C) |
| $V_{xc}$ | Exchange-Correlation | `assemble_olsxc_off` | Pending |

## 5. H2 Verification Results (Summary)

Verification was performed on an H2 molecule ($d=0.74$ Å). Results compared PyOpenCL `OCL_Hamiltonian` against Fortran `libFireCore` sparse export.

| Component | Max Abs Diff | Status | Notes |
| :--- | :--- | :--- | :--- |
| **Overlap S** | $4.86 \times 10^{-8}$ | **PASSED** | Correct rotation and interpolation. |
| **Kinetic T** | $2.39 \times 10^{-7}$ | **PASSED** | Correct rotation and interpolation. |
| **Vna (Total)**| $4.16 \times 10^{2}$ | **FAILED** | Major magnitude and diagonal mismatch. |

## 6. Vna Discrepancy Analysis

The discrepancy in $V_{na}$ is attributed to two main factors:

### A. Missing 3-Center Assembly
In the McWEDA/OLSXC scheme, $V_{na}$ is accumulated from both 2-center and 3-center terms.
- **Fortran**: Calls `assemble_3c` which adds massive contributions to the `vna` matrix (e.g., diagonal shifts of ~411 for H2).
- **PyOpenCL**: Currently only implements 2-center spline interpolation terms.

### B. Energy Scaling Factor (`eq2`)
The Fortran implementation uses a global constant `eq2 = 14.39975` (defined in `constants_fireball.f90`) to scale electrostatic potential terms, converting to common energy units (eV/Angstrom). PyOpenCL currently lacks this scaling.

## 7. Technical Implementation Details

### Rotation Logic (Slater-Koster)
Replicated logic from `ROTATIONS/epsilon.f90` and `ROTATIONS/twister.f90`.
- **Local Frame**: Bonds are aligned along the Z-axis in the local frame.
- **Orbital Order**: $(s, p_y, p_z, p_x)$ following Fireball convention.
- **Transformations**: Implemented in OpenCL kernels to handle arbitrary bond orientations.

### Interpolation Schemes
- **2-Center**: Cubic spline interpolation (replicating `buildspline_1d.f90`). Spline coefficients $(y, y'')$ are pre-calculated in Python and passed as `float4` to GPU.
- **3-Center**: 2D grid interpolation ($x, y$ - distance to NA and bond length) followed by Legendre polynomial expansion in $\cos \theta$. Replicated from `trescentros.f90` and `interpolate_2d.f90`.

## 8. Current Status & Verification

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

## 9. Development Notes & References

- **Sparse Matrix**: Uses an atom-neighbor block format `[n_pairs, norb, norb]`.
- **Fireball Convention**: The "Ortega convention" for $p$-orbitals: $p_y, p_z, p_x$ mapped to indices $1, 2, 3$ (0-indexed: 1:y, 2:z, 3:x).
- **Theoretical Basis**:
    - PRB 40, 3979 (1989) - OLSXC / McWEDA.
    - Sankey and Niklewski, PRB 40, 3979 (1989) - Basis functions and integrals.

---
*Created by Antigravity AI on 2025-12-21*

## Recent Implementations (2025)

- **Scanning API for Granular Verification**
  - Added Fortran subroutines `firecore_scanHamPiece2c` and `firecore_scanHamPiece3c` in `fortran/MAIN/libFireCore.f90` to expose 2‑center and 3‑center integral scans.
  - Updated `pyBall/FireCore.py` with ctypes bindings `scanHamPiece2c` and `scanHamPiece3c`.
  - Implemented corresponding methods in `pyBall/FireballOCL/OCL_Hamiltonian.py` (`scanHamPiece2c`, `scanHamPiece3c`) that call the OpenCL kernels.
  - Created verification script `tests/pyFireball/verify_scan.py` performing radial, angular, and 3‑center scans and generating Matplotlib plots.
  - Fixed compilation issues by adding `firecore_options` module for assembly toggles and correcting module imports.
- **3‑Center Kernel Enhancements**
  - Added `assemble_3c` method and kernel support in `cl/hamiltonian.cl`.
  - Integrated 3‑center scans into verification workflow.
- **Rotation Logic Verification**
  - Added angular scan to verify Slater‑Koster rotation against Fortran reference.
- **Documentation Updates**
  - Updated this implementation progress document to reflect the new scanning capabilities and their locations.

These additions enable detailed, component‑wise testing of interpolation and rotation logic, facilitating future development and debugging.
