# FireballOCL Implementation Progress
## Latest status (2025-12-21)

- **2-center (finished)**  
  - PyOpenCL now loads 2c data with the same naming as Fortran `read_2c`: `overlap`, `kinetic`, and all vna variants `vna_atom_<isorp>`, `vna_ontopl_<isorp>`, `vna_ontopr_<isorp>`, iterating `isorp` up to `nssh`.  
  - Generic `vna` requests are satisfied via the `vna_atom_*` tables (matching Fortran interaction 4 usage).  
  - Sign conventions for py/pz rotation in OpenCL kernels match Fortran.  
  - `verify_scan_2c.py` compares Fortran vs PyOpenCL and matches (no NaNs, diffs within threshold).
  - **H2 parity achieved (S, T, Vna):** `verify_H2.py` now passes; the missing +14.759 eV diagonal in Vna was traced to skipped `vna_atom` on self neighbors (`i==j`) and missing accumulation into `neigh_self`. OpenCL now evaluates `vna_atom_00` for both `i!=j` and `i==j` and accumulates into the `(i,i)` block exactly as Fortran `assemble_2c` does.
  - **s-px sign fix:** Based on inverted s-px curve in `scan_2c_vna.png`, flipped `ezx` to `-ez.x` in `assemble_2c` kernel to match Fortran convention.
  - **Vna s-pz sign fix:** Identified that Vna (interaction=4) sp_sig spline data has opposite sign to S/T. Fixed by flipping sp_sig spline coeffs for Vna pairs in `prepare_splines` and adding conditional `ezz` rotation in kernel (`ezz = ez.z` for Vna, `-ez.z` for others). This aligns with Fortran's data loading and rotation conventions, ensuring robustness across geometries.

- **Vna (Neutral Atom Potential) - COMPLETED (2025-12-22)**  
  - **C2 parity achieved:** `verify_C2.py` now passes Vna (max diff ~2e-6). The remaining s-pz sign flip in C2 was due to Vna (interaction=4) using a different rotation frame than S/T: Fortran uses `epsilon(r2, sighat, eps)` with absolute neighbor position `r2`, while OpenCL initially used only relative dR. Implemented Fortran-accurate `epsilon(r_j, sighat, eps)` in `assemble_2c` for `is_vna_pair` flagged blocks, keeping S/T unchanged. This fixed the odd-parity s-p components without affecting other interactions.

- **Vnl (Non-local Pseudopotential) - PARTIALLY FIXED (2025-12-22)**  
  - **Scans (radial + angular) now match Fortran:** Added Vnl to `verify_scan_2c.py` with both radial and angular scans, plus extended angular scans for kinetic and Vna. The key fix for scans was routing `root='vnl'` to use `assemble_pp` kernel (PP projector overlaps) instead of `assemble_2c` (wrong kernel). Radial and angular Vnl scans match Fortran within ~1e-6.
  - **Full C2 verification now passes:** `verify_C2.py` now reports `SUCCESS: Vnl matches` after fixing the **contraction orientation/layout** in the CPU assembly path.
  - **Why scans were correct but verify_C2 initially failed:** scans validate the raw `sVNL` blocks (output of `assemble_pp`). Full Vnl assembly additionally performs the contraction (Fortran `assemble_2c_PP.f90`):
    `vnl(imu,inu) += sum_cc cl(cc) * sVNL(imu,cc) * sVNL(inu,cc)`.
    Our `assemble_pp` kernel outputs blocks in (inu,imu) convention (like other kernels), but the initial CPU contraction treated them as (imu,cc) without consistent transposes, producing a large apparent mismatch. Fixing the transpose handling and returning Vnl blocks in (inu,imu) convention resolved this.

- **3-center (in progress)**  
  - **Data parity achieved:** Added `firecore_export_bcna_table` (Fortran) and `export_bcna_table` (ctypes) + `verify_Fdata_read_3c.py` plots; Fortran vs Python tables match bitwise for H bcna_01..05 (max|diff|=0).  
  - **Remaining issue:** OpenCL raw scan still mismatches Fortran on L0 (bcna_01) despite identical tables; likely a host/kernel mapping/stride bug (component selection or triplet index). L1/L2 etc. match to ~1e-9.  
  - `verify_scan_3c.py --mode raw` highlights this gap; interpolation kernel logic itself matches Fortran (off-by-one fixed), so focus is on the bcna component selection in the OpenCL path.
  - **Kernel/host fixes done:** Added missing `type_idx` arg to 3c kernels, propagated `isorp` stride through packing and kernel indexing, and fixed off-by-one in 3c interpolation grid access. These bring all components except L0 into agreement (~1e-9).
  - **Data-loading bugs fixed:**  
    - 3c header parsing in Python skipped wrong lines (charge lines with `<===`), causing wrong reshapes; fixed by matching Fortran `readheader_3c`.  
    - `num_nonzero` for 3c now uses the 3c-specific count (not 2c).  
    - Verification script `tests/pyFireball/verify_Fdata_read_3c.py` visualizes Fortran vs Python grids (5×3 subplot: Fortran/Python/Diff) and prints min/max/diff; confirms tables load identically.

- **Next steps**  
  - Fix Vnl full assembly discrepancy: Debug the CPU contraction in `assemble_full` by comparing Python contraction terms against Fortran `assemble_2c_PP` loop. Add prints for individual `cl(cc) * sVNL(imu,ncc) * sVNL(inu,ncc)` terms to isolate sign/indexing errors. Once fixed, Vnl will pass in `verify_C2.py`.
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
| $V_{na}$ | Neutral Atom Potential | `assemble_2c`, `assemble_3c` | **Verified** for 2c (H2 parity); 3c assembly still pending |
| $V_{nl}$ | Non-local Pseudopotential | `assemble_2c_PP` | **Verified** (Tol: $10^{-6}$) |
| $\bar{\rho}$ | Average Density (3C) | `average_ca_rho` | **Implemented** (Needs 3C) |
| $V_{xc}$ | Exchange-Correlation | `assemble_olsxc_off` | Pending |

## 5. H2 Verification Results (Summary)

Verification was performed on an H2 molecule ($d=0.74$ Å). Results compared PyOpenCL `OCL_Hamiltonian` against Fortran `libFireCore` sparse export.

| Component | Max Abs Diff | Status | Notes |
| :--- | :--- | :--- | :--- |
| **Overlap S** | $4.86 \times 10^{-8}$ | **PASSED** | Correct rotation and interpolation. |
| **Kinetic T** | $2.39 \times 10^{-7}$ | **PASSED** | Correct rotation and interpolation. |
| **Vna (Total)**| $2 \times 10^{-6}$ | **PASSED** | 2c Vna (ontop + atom) matches Fortran; 3c Vna not yet assembled in PyOpenCL. |

## 6. Vna Status

- **2-center Vna parity achieved (H2):** `vna_atom` (interaction=4) is now evaluated for both `i!=j` and `i==j` and accumulated into the self block `(i,i)` (Fortran `neigh_self` semantics). Off-diagonal ontop terms (interactions 2 & 3) also match. Scaling with `eq2` is applied identically.
- **3-center Vna:** Not yet assembled in PyOpenCL. Fortran adds 3c Vna in `assemble_3c`; OpenCL still needs equivalent rotation/assembly.

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

### C2 Verification (2-center, 4 orbitals on Carbon)

Verification was performed on a C2 dimer (two Carbon atoms, 4 orbitals per atom: s, py, pz, px). Dense 8x8 matrices were reconstructed from the Fortran sparse blocked export and from the OpenCL blocked output.

| Component | Status | Notes |
| :--- | :--- | :--- |
| **Overlap S** | **PASSED** | Matches within ~1e-7 after correct 8x8 reconstruction. |
| **Kinetic T** | **PASSED** | Matches within ~1e-6. |
| **Vna (Total)** | **PASSED** | Matches within ~1e-6 after fixing vna_atom direction convention (see below). |
| **Vnl** | **PASSED** | Matches within ~1e-6 after fixing PP contraction and block orientation. |

### Next Steps
- [ ] Implement full **McWEDA / OLSXC** logic for XC potentials.
- [ ] Full integration of `assemble_full` with the main SCF loop.
- [ ] Rigorous benchmarking and verification against `FireCore` sparse export.

## 9. Development Notes & References

- **Sparse Matrix (Fortran allocation)**: Hamiltonian/overlap blocks are stored as
  `h_mat(numorb_max, numorb_max, neigh_max, natoms)` and `s_mat(numorb_max, numorb_max, neigh_max, natoms)` (same for `t_mat`).
  The Python interface exports this blocked structure; verification scripts reconstruct dense matrices by mapping neighbor blocks into global orbital offsets.
- **Fireball Convention**: The "Ortega convention" for $p$-orbitals: $p_y, p_z, p_x$ mapped to indices $1, 2, 3$ (0-indexed: 1:y, 2:z, 3:x).

### Vna Convention Bug (C2)

For H2, 2-center Vna parity was already achieved. C2 initially showed a large sign flip specifically in odd-parity s-p couplings (e.g. `H(s,p)` terms), while S and T were correct.

Root cause:
- In Fortran (`assemble_2c.f90`), `vna_atom` (interaction=4) is evaluated for geometry (i,j) but accumulated into the on-site self block `(i,i)` via `neigh_self(iatom)`:
  `vna(imu,inu,neigh_self(i),i) += <i|v(j)|i> * eq2`.
- In OpenCL, we evaluated the same table but used the same bond direction convention as the off-diagonal 2-center blocks.
- For multi-orbital species (Carbon), the `vna_atom` contribution contains odd-parity components whose sign depends on the bond direction convention. This produced the observed sign flips in C2.

Fix:
- When accumulating `vna_atom_00` into the self block `(i,i)`, OpenCL now evaluates the `vna_atom` block with reversed direction (swap i/j in the 2-center evaluation) but still stores it into `(i,i)`. This matches the Fortran convention for C2 and restores parity.
- **Theoretical Basis**:
    - PRB 40, 3979 (1989) - OLSXC / McWEDA.
    - Sankey and Niklewski, PRB 40, 3979 (1989) - Basis functions and integrals.

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

### Key Insights from C2 Vnl Debugging Session (2025-12-22)

This session resolved the final C2 verification issue, achieving **full 2-center parity** (S/T/Vna/Vnl pass with max diffs ~1e-6 to 2e-6). Key learnings:

#### Vnl Contraction Details
- **Fortran Vnl Assembly** (from `assemble_2c_PP.f90`):
  - Computes `sVNL` blocks via `doscentrosPP` (interpolate + recover_PP + rotatePP).
  - Contracts into full Vnl via: `vnl(imu,inu) += sum_cc cl(cc) * sVNL(imu,cc) * sVNL(inu,cc)`
  - `cl(cc)` are pseudopotential coefficients from `info.dat` (stored in `cl_pseudo`).

- **OpenCL Implementation**:
  - `assemble_pp` kernel: Interpolates, recovers, and rotates `sVNL` blocks (4x4, (inu,imu) convention).
  - CPU `contract_AB`: Implements the sum via matrix operations (transposes for indexing), returns (inu,imu) blocks.
  - Assembly: Accumulates contracted blocks into full Vnl matrix.

- **Bug and Fix**: Initial contraction produced blocks in wrong orientation (imu,inu), causing mismatch in dense reconstruction. Fixed by transposing the final result to match other kernels.

#### Takeaways and Precautions
- **Rotation Frames Vary by Interaction**: Don't assume uniform rotation (e.g., Vna uses `epsilon(r2, sighat)` with absolute `r2`; S/T use relative dR). Always check Fortran's `assemble_2c` or specific subroutines for each interaction.
- **Block Orientations Matter**: Kernels output in (inu,imu) for Python dense reconstruction (which transposes blocks). Contractions must maintain this to avoid silent mismatches in full assembly.
- **Scans ≠ Full Assembly**: Scans test raw kernels (e.g., sVNL blocks); full assembly adds contractions and accumulations. Verify both paths.
- **Indexing Conventions**: Fortran uses 1-based indices; OpenCL 0-based. Double-check muPP/nuPP mappings and orbital orders (Ortega: s, py, pz, px).
- **Debugging Strategy**: Use targeted prints (e.g., element-wise comparisons) and isolate issues (e.g., kernel vs contraction). Don't guess—trace Fortran paths.
- **Avoid Repeating Mistakes**: Always verify full `verify_C2.py` after scan fixes. Check Fortran data loading (e.g., sign flips in splines). Use `git status` to track changes.
