# Nanocrystal Vibrations: Generation, Force Fields, and Spectroscopy

## Overview

This topic covers the complete pipeline for generating silicon/diamond nanocrystal structures and computing their vibrational properties using classical force fields (MMFF/UFF). The workflow spans structure generation (JavaScript), force field evaluation (C++/Python), and post-processing for vibration spectra (Python/NumPy).

**Primary application**: Compute FTIR absorption spectra and phonon band structures for Si and C nanocrystals to compare with experimental spectroscopy data.

---

## Implementations

### 1. Nanocrystal Generation (JavaScript)

| File | Status | Description |
|------|--------|-------------|
| `scripts/gen_nanocrystals.mjs` | **Active** | High-throughput generation of Si nanocrystals with configurable Miller-index plane cuts, hydrogen capping, and surface bridge operations |
| `doc/Topics/FTIR_Nanocrystals/gen_nanocrystals.md` | **Active** | Documentation for the generation script; covers CLI arguments, plane templates, and output formats |

**Key functions** (in `gen_nanocrystals.mjs`):
- `parseArgs()` — CLI argument parsing (CIF path, plane templates, capping element, pruning threshold)
- `main()` — Load CIF → build supercell → apply Miller cuts → recalc bonds → prune undercoordinated → H-cap dangling bonds → output `.mol2`

**Dependencies**: Uses `web/molgui_webgpu/` modules (`Vec3`, `MMParams`, `EditableMolecule`, `CrystalUtils`, `MoleculeIO`, `MoleculeSelection`, `MoleculeUtils`)

---

### 2. Force Field Evaluation & Hessian Computation (C++ / Python)

#### 2a. Python API (`pyBall/MMFF.py`)

| Function | Status | Description |
|----------|--------|-------------|
| `getHessian3Nx3N(inds, dx)` | **Active** | Compute 3N×3N Hessian via central finite differences on C++ side; returns numpy array |
| `getPhononPhiBlocks(inds_total, inds_disp, dx)` | **Active** | PBC-aware force-constant blocks for phonon band structure; displace central-cell atoms, read forces on all supercell atoms |
| `setBondParamsByType(ti, tj, k, r0, forcefield)` | **Active** | Generalized selective bond modification by atom types (UFF or MMFF auto-detected) |
| `setAngleParamsByType(ti, tj, tk, k, c0, forcefield)` | **Active** | Generalized selective angle modification by atom types (UFF or MMFF auto-detected) |
| `getBuffs()` / `getBuffs_UFF()` | **Active** | Expose C++ parameter arrays as numpy views (`bKs`, `bLs`, `apars`, `Ksp`, `Kpp`, `atypes`, `neighs`, `neighCell`, etc.) |
| `getHessianContext()` | **Active** | Return bond/angle basis matrices and atom indices for Hessian fitting |
| `assembleHessianFromParams()` | **Active** | Assemble Hessian from fitted parameters (validation) |

#### 2b. C++ Implementation (`cpp/libs/Molecular/MMFF_lib.cpp`)

| Function | Status | Description |
|----------|--------|-------------|
| `getHessian3Nx3N()` | **Active** | Central FD loop: perturb each atom ±dx, evaluate forces, assemble symmetric Hessian |
| `getPhononPhiBlocks()` | **Active** | PBC validation → FD on central cell → return Phi[o,p] blocks; calls `NBFF::checkPBCNeighCells()` |
| `init_buffers()` | **Active** | Expose MMFF arrays to Python (`bKs`, `bLs`, `apars`, `Ksp`, `Kpp`, `atypes`, etc.) |
| `init_buffers_UFF()` | **Active** | Expose UFF arrays to Python (`bonParams`, `angParams`, `bonAtoms`, `angAtoms`, `atypes`) |

#### 2c. C++ Validation (`cpp/common/molecular/NBFF.h`)

| Function | Status | Description |
|----------|--------|-------------|
| `checkPBCNeighCells(func_name)` | **Active** | Reusable PBC validation: ensures `neighCell >= 0` for all valid neighbors when `bPBC=true`; eliminates code duplication across `MMFF_lib.cpp` and `MolWorld_sp3.h` |

#### 2d. C++ Force Field Classes

| File | Status | Description |
|------|--------|-------------|
| `cpp/common/molecular/MMFFsp3_loc.h` | **Active** | MMFFsp3_loc force field: node-based bond/angle evaluation with pi-orbital tracking |
| `cpp/common/molecular/NBFF.h` | **Active** | Base class for non-bonded force fields; owns `neighs`, `neighCell`, `shifts`, PBC infrastructure |
| `cpp/common/molecular/MolWorld_sp3.h` | **Active** | High-level world container; `makeMMFFs()` sets up UFF or MMFF force field from builder |

---

### 3. Vibration Spectrum & Phonon Post-Processing (Python)

#### 3a. Linear Response FTIR (`pyBall/FTIR.py`)

| Function | Status | Description |
|----------|--------|-------------|
| `mechanical_greens_probing(K, M, omegas, ...)` | **Active** | Dipole-driven Green's function probing: solves `A(ω)U = f` at each frequency without full diagonalization |
| `dynamic_stiffness(K, M, omega, eta)` | **Active** | Build `A = K - (ω + iη)²M` with stabilization |
| `project_rigid_modes(H, M, pos, shift)` | **Active** | Project out 6 rigid-body modes (translation + rotation) by eigenvalue shifting |
| `get_mass_matrix(mmff, n_atoms)` | **Active** | Build diagonal mass matrix from MMFF atom types |
| `fit_hessian_parameters(...)` | **Active** | Linear least-squares fitting of bond/angle stiffnesses to target Hessian (e.g., DFT reference) |
| `build_design_matrix_from_basis(...)` | **Active** | Assemble design matrix from bond/angle basis matrices for fitting |
| `save_xyz_vib(...)` | **Active** | Save vibration vectors to XYZ for Jmol visualization |

#### 3b. Phonon Band Structure (Python scripts in `tests/tMMFF/`)

| Script | Status | Description |
|--------|--------|-------------|
| `test_diamond_phonon_bands.py` | **Active** | Diamond phonon bands via supercell + Bloch extraction; supports PBC, ASR, UFF/MMFF toggle, bond/angle scaling |
| `test_diamond_gamma.py` | **Active** | Γ-point phonon frequencies for diamond (simplified, no k-path) |
| `test_ethane_gamma.py` | **Active** | Γ-point frequencies for ethane (molecular, no PBC) |
| `plot_phonon_bands.py` | **Active** | Plot phonon dispersion curves from `.npz` output |
| `plot_mmff_comparison.py` | **Active** | Compare MMFF vs reference phonon bands |
| `plot_uff_comparison.py` | **Active** | Compare UFF vs reference phonon bands |
| `plot_mmff_scaling.py` | **Active** | Overlay default vs scaled parameter phonon bands |
| `plot_angle_comparison.py` | **Active** | Compare angle parameter effect on phonon bands |

#### 3c. Hessian Tests & Validation

| Script | Status | Description |
|--------|--------|-------------|
| `test_vibration_spectra.py` | **Active** | End-to-end linear response spectrum: init MMFF → Hessian → mass matrix → project rigid modes → probe spectrum → plot |
| `test_hessian_fitting.py` | **Active** | Fit UFF bond/angle parameters to reference Hessian; validates `fit_hessian_parameters()` |
| `test_diatomic_hessian.py` | **Active** | Hessian calculation for diatomic molecules (simple validation) |
| `run_hessian.py` | **Active** | Batch Hessian runner for multiple systems |
| `test_asan_minimal.py` | **Active** | AddressSanitizer memory test for MMFF library |

---

### 4. Documentation

| File | Status | Description |
|------|--------|-------------|
| `doc/Topics/FTIR_Nanocrystals/Phonon_testing_guide.md` | **Active** | Practical guide for phonon testing with MMFF; covers common pitfalls (imaginary modes, ASR, PBC) |
| `doc/Topics/FTIR_Nanocrystals/Hessian_Kspace.md` | **Active** | Theory: Hessian in k-space for periodic crystals, Bloch extraction, and dynamical matrix construction |
| `doc/Topics/FTIR_Nanocrystals/Hessian_fitting.md` | **Active** | Theory and practice: fitting disordered nanocrystal Hessians, 3D B-spline grids, parameter fitting heuristics |
| `doc/Topics/FTIR_Nanocrystals/Debug_negative_phonon_freqs.md` | **Active** | Debugging guide for imaginary phonon modes; checklist of causes and fixes |
| `doc/Topics/FTIR_Nanocrystals/gen_nanocrystals.md` | **Active** | Nanocrystal generation documentation |

---

## Relationships

```
Nanocrystal Generation              Force Field Setup                Vibration Calculation
─────────────────────              ─────────────────                ─────────────────────
gen_nanocrystals.mjs        ──►     MMFF.py::init()          ──►     test_vibration_spectra.py
(CIF → supercell → cuts)           (load .mol2 → build topo)         (Hessian → FTIR spectrum)
                                          │
                                          ▼
                                   MMFF_lib.cpp::getHessian3Nx3N()
                                   (C++ central FD)
                                          │
                                          ▼
                                   test_diamond_phonon_bands.py
                                   (supercell → Phi blocks → Bloch → bands)
```

**Cross-language flow**:
1. **JavaScript** (`gen_nanocrystals.mjs`) generates structures → `.mol2`/`.xyz`
2. **Python** (`MMFF.py`) loads structures, configures force field, calls C++ via ctypes
3. **C++** (`MMFF_lib.cpp`) computes Hessian/phonon blocks via finite differences
4. **Python** (`FTIR.py`, `test_diamond_phonon_bands.py`) post-processes: phonon bands, FTIR spectra, parameter fitting

---

## Test Matrix

| Test Script | What it tests | Force Field | PBC | Key Flags |
|-------------|---------------|-------------|-----|-----------|
| `test_diamond_phonon_bands.py` | Phonon dispersion | MMFF/UFF | Optional | `--pbc`, `--uff`, `--scale-bond`, `--scale-angle`, `--asr` |
| `test_diamond_gamma.py` | Γ-point frequencies | MMFF | No | — |
| `test_ethane_gamma.py` | Molecular vibrations | MMFF | No | — |
| `test_vibration_spectra.py` | FTIR spectrum | MMFF | No | `--fmin`, `--fmax`, `--eta`, `--dx` |
| `test_hessian_fitting.py` | Parameter fitting | UFF | No | — |
| `test_diatomic_hessian.py` | Hessian correctness | MMFF | No | — |
| `test_asan_minimal.py` | Memory safety | MMFF | No | — |

---

## Status Summary

| Component | Status | Notes |
|-----------|--------|-------|
| Nanocrystal generation | **Active** | JavaScript pipeline; used for high-throughput Si structures |
| MMFF force field (C++) | **Active** | Production-ready; GPU path exists (OpenCL) but CPU path is primary for vibrations |
| UFF force field (C++) | **Active** | Used for comparison/validation; parameter fitting uses UFF |
| Hessian computation | **Active** | Central FD in C++; accurate but O(N²) force evaluations |
| Phonon bands (PBC) | **Active** | Supercell method with Bloch extraction; `checkPBCNeighCells()` validation ensures correctness |
| FTIR spectrum | **Active** | Green's function probing avoids full diagonalization; good for large nanocrystals |
| Parameter fitting | **Active** | Linear LS fit of bond/angle stiffnesses to DFT reference Hessians |
| GPU acceleration | **Experimental** | OpenCL/CUDA paths exist for MD but not yet for Hessian/phonon calculations |

---

## Known Limitations & TODOs

1. **PBC optical modes**: Cluster mode (no PBC) gives correct Γ optical modes but acoustic branches may go imaginary at finite k due to surface under-coordination. Use `--asr` (phonopy symmetrization) or PBC with caution.
2. **Hessian cost**: Central FD requires 6N force evaluations. For nanocrystals >1000 atoms, this becomes expensive. No GPU-accelerated Hessian yet.
3. **Parameter fitting**: Currently only bond/angle stiffnesses fitted; non-bonded parameters and pi-orbital terms not included.
4. **Nanocrystal surface**: Hydrogen capping in `gen_nanocrystals.mjs` is static; no relaxation of surface H positions before vibration calculation.

---

## Related Topics

- [Classical Force Fields (MMFF/UFF)](topical_audit.md#classical-force-fields-mmffuff) — Broader force field overview
- [Htransfer_Kekule_DFTB.md](Htransfer_Kekule_DFTB.md) — Related DFTB/force field work
- Codemap: [Silicon/Diamond Nanocrystal Vibration Spectroscopy Workflow](https://windsurf.com/codemaps) — Visual trace of the complete pipeline
