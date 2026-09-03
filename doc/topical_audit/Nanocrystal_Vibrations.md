---
type: TopicalAudit
title: Nanocrystal Vibrations
tags: [nanocrystal, ftir, phonon, mmff]
timestamp: 2026-09-02
---

# Nanocrystal Vibrations: Generation, Force Fields, and Spectroscopy

## Documentation map

This file is the **topical audit** (code inventory and status). It is maintained in sync with:

| Location | Role |
|----------|------|
| [`SiNCs.md`](SiNCs.md) | **Project briefing** (FTIR + XRD + FFfit) — paste into ChatGPT |
| [`XRD.md`](XRD.md) | Powder XRD / PDF inventory + Kusová strain-gradient checklist |
| [`Hessian_fitting.md`](Hessian_fitting.md) | Hessian / Wilson FFfit inventory + method checklist |
| **This file** | APIs, implementations, test matrix, limitations |
| [`doc/Topics/FTIR_Nanocrystals/PySCF_L1_neighborhood_PDOS.md`](../Topics/FTIR_Nanocrystals/PySCF_L1_neighborhood_PDOS.md) | PySCF L1 stacked PDOS — **plot template**; vs DFTB+ L1; Si Hessians not minima |
| [`doc/Topics/FTIR_Nanocrystals/MMFF_C_CH_vs_CH2_kfit.md`](../Topics/FTIR_Nanocrystals/MMFF_C_CH_vs_CH2_kfit.md) | C L1 own-min \(k\)-fit vs PBE; Morse + split angles; SA did not unify cube/octa |
| [`tests/tSiNCs/AGENTS.md`](../../tests/tSiNCs/AGENTS.md) | DOX contract — script ownership, `REPO`/`TEST_DIR` paths |
| [`doc/Topics/FTIR_Nanocrystals/README.md`](../Topics/FTIR_Nanocrystals/README.md) | Topic docs index (guides, chats, progress logs) |
| [`doc/Topics/FTIR_Nanocrystals/ChemAtlas_MMFF_relax.md`](../Topics/FTIR_Nanocrystals/ChemAtlas_MMFF_relax.md) | Atlas results location + surface-NB relax |
| [`tests/tSiNCs/README.md`](../../tests/tSiNCs/README.md) | Working hub — quick start, fixtures, viewers |
| [`tests/tSiNCs/ToDo_Nanocrystal.md`](../../tests/tSiNCs/ToDo_Nanocrystal.md) | Open items |

---

## Overview

Pipeline for generating silicon / diamond nanocrystal structures and computing vibrational properties with classical force fields (MMFF/UFF): structure generation (JavaScript + Python), force field evaluation (C++/Python), post-processing for FTIR and phonons (Python/NumPy).

**Primary application:** FTIR absorption spectra and phonon band structures for Si and C nanocrystals, compared with experiment and QM references (adamantane, sila-adamantane in `tests/tSiNCs/`).

**Canonical production spectrum method:** dense diagonalization (`np.linalg.eigh`) + mode histogram via `pyBall/nanocrystal_pipeline.py` and `FTIR.vibration_spectrum_from_modes`. See [`Nanocrystal_pipeline.progress.md`](../Topics/FTIR_Nanocrystals/Nanocrystal_pipeline.progress.md).

**Alternate method:** Green's-function frequency probing (`FTIR.mechanical_greens_probing`) — linear response without full diagonalization; used in `tests/tMMFF/test_vibration_spectra.py`. Sparse/GPU frequency solvers are **de-prioritized** (see [`Sparse_vibration_solver.progress.md`](../Topics/FTIR_Nanocrystals/Sparse_vibration_solver.progress.md)).

---

## Implementations

### 1. Nanocrystal generation

| File | Status | Description |
|------|--------|-------------|
| `web/molgui_webgpu/Nanocrystals.js` | **Active** | Core library: CIF → supercell → Miller/sphere cuts → prune → H-cap → bridges |
| `tests/tSiNCs/nanocrystals.mjs` | **Active** | Unified CLI: `generate`, `ensemble`, `topology`, `audit`, `nonbond`, `rings` |
| `tests/tSiNCs/gen_nanocrystals.mjs` | **Deprecated** | Thin wrapper; use `nanocrystals.mjs generate` |
| `tests/tSiNCs/gen_afm_tip.mjs` | **Active** | AFM-tip nanocrystal generator: truncated tetrahedron with [111] apex, Si + diamond |
| `tests/tSiNCs/gen_nanocrystals.py` | **Active** | Python CLI: spherical cuts native; Miller planes via Node |
| `pyBall/nanocrystal_gen.py` | **Active** | Python sphere-cut builder; parity target for JS |
| `doc/Topics/FTIR_Nanocrystals/gen_nanocrystals.chat.md` | **Active** | Generation CLI and `Nanocrystals.js` reference |

JS is feature-complete (plane cuts, defects, bridges). Python covers spherical cuts; Miller planes delegate to Node.

**Parity:** `tests/tSiNCs/crosscheck_nanocrystal_generators.py` — see [`tests/tSiNCs/README.md`](../../tests/tSiNCs/README.md).

---

### 2. Force field evaluation and Hessian (C++ / Python)

#### 2a. Python API (`pyBall/MMFF.py`)

| Function | Status | Description |
|----------|--------|-------------|
| `getHessian3Nx3N(inds, dx)` | **Active** | 3N×3N Hessian via central FD in C++ |
| `getHessianSparseBlocks(...)` | **Active** | Neighbor-shell sparse 3×3 blocks |
| `getPhononPhiBlocks(inds_total, inds_disp, dx)` | **Active** | PBC force-constant blocks for phonon bands |
| `setBondParamsByType` / `setAngleParamsByType` | **Active** | Selective bond/angle modification (UFF or MMFF) |
| `setMorseNonBond(i, K)` | **Active** | Exclusion2 cluster Morse α (`NBFF::bMorseNonBond`); GridFF Morse is a different path |
| `setEachAngle(i)` | **Active** | Per-slot `angles[nnode×6]` \(K_{\mathrm{ss}}\) (XXX / XXH / HXH). Default kernel is **one Kss per node** (`apars[:,2]`) |
| `getBuffs()` / `getBuffs_UFF()` | **Active** | Expose C++ parameter arrays as numpy views (`angles` included) |
| `getHessianContext()` / `assembleHessianFromParams()` | **Active** | Hessian fitting basis and validation |

#### 2b. C++ (`cpp/libs/Molecular/MMFF_lib.cpp`)

| Function | Status | Description |
|----------|--------|-------------|
| `getHessian3Nx3N()` | **Active** | Central FD Hessian. Cluster: keep Exclusion2 if FIRE had it (`eval_no_omp`). PBC: bonded-only |
| `setMorseNonBond` / `setEachAngle` | **Active** | Cluster Morse vs LJ; per-angle \(K_{\mathrm{ss}}\) |
| `getHessianSparseBlocks()` | **Active** | Sparse block extraction |
| `getPhononPhiBlocks()` | **Active** | PBC phonon Φ blocks; uses `NBFF::checkPBCNeighCells()` |
| `init_buffers()` / `init_buffers_UFF()` | **Active** | Python buffer exposure |

#### 2c. Related C++

| File | Role |
|------|------|
| `cpp/common/molecular/MMFFsp3_loc.h` | MMFFsp3 bond/angle forces |
| `cpp/common/molecular/NBFF.h` | PBC neighbors; `checkPBCNeighCells()`; `bMorseNonBond` on Exclusion2 |
| `cpp/common/molecular/MolWorld_sp3.h` | World container; `makeMMFFs()` |
| `cpp/common_resources/crystals/` | `Si_primitive`, `diamond_primitive`, CIF/XYZ inputs |

---

### 3. Vibration spectrum and phonon post-processing (Python)

#### 3a. `pyBall/FTIR.py`

| Function | Status | Description |
|----------|--------|-------------|
| `vibration_spectrum_from_modes(...)` | **Active** | **Production:** dense `eigh` + mode summation |
| `build_hessian_from_linear_topology(...)` | **Active** | Analytical 3N×3N K from exported `03_topology.npz` springs |
| `harmonic_stick_hessian_blocks`, `mass_matrix_from_Z` | **Active** | Per-stick Hessian blocks; masses from `Z` |
| `mechanical_greens_probing(...)` | **Active** | Frequency-domain Green's probing (linear response FTIR) |
| `mechanical_greens_probing_sparse(...)` | **Deferred** | Sparse variant; not on ensemble critical path |
| `dynamic_stiffness`, `project_rigid_modes`, `get_mass_matrix` | **Active** | Dynamical matrix, rigid modes, masses |
| `fit_hessian_parameters`, `build_design_matrix_from_basis` | **Active** | Bond/angle stiffness fitting |
| `save_xyz_vib` | **Active** | Jmol-compatible mode export |

#### 3b. Ensemble pipeline

| File | Status | Description |
|------|--------|-------------|
| `pyBall/nanocrystal_pipeline.py` | **Active** | NPZ stages: `relax --nonbond surface` (collision-group atom set) → topology-linear Hessian → `eigh` → spectrum; atlas L0/L1 uses mol2 init |
| `tests/tSiNCs/chem_atlas.json` | **Active** | Wulff Si ≡ C + L1 5/7-ring; `relax_tiers: L0_ref, L1_dft` |
| `pyBall/io/crystal_npz.py` | **Active** | `load_crystal_npz`, `load_topology_npz`, `validate_topology_crystal_parity` |
| `tests/tSiNCs/export_nanocrystal_bundle.mjs` | **Active** | JS export `01_crystal.npz` + `03_topology.npz` |
| `tests/tSiNCs/run_nanocrystal_ensemble.mjs` | **Deprecated** | Legacy ensemble wrapper; use `nanocrystals.mjs ensemble` |
| `tests/tSiNCs/ensemble.example.json` | **Active** | Example ensemble config |
| `tests/tSiNCs/atlas_shapes.json` | **Active** | Shape atlas presets for `--atlas` |
| `tests/tSiNCs/make_small_symmetric_nc.mjs` | **Active** | Symmetric Si/C sphere batch (&lt;100 atoms) |

#### 3c. GPU / sparse solvers (experimental, deferred)

| File | Status | Description |
|------|--------|-------------|
| `pyBall/OCL/VibSolver.py` | **Experimental** | pyOpenCL wrapper for `vib_jacobi.cl` |
| `cpp/common_resources/cl/vib_jacobi.cl` | **Experimental** | Block-Jacobi frequency-domain solver |

#### 3d. Phonon and vibration tests (`tests/tMMFF/`)

| Script | Status | Description |
|--------|--------|-------------|
| `test_diamond_phonon_bands.py` | **Active** | Diamond phonon bands; PBC, ASR, UFF/MMFF, scaling |
| `test_diamond_gamma.py` | **Active** | Γ-point diamond frequencies |
| `test_ethane_gamma.py` | **Active** | Γ-point molecular (no PBC) |
| `test_vibration_spectra.py` | **Active** | End-to-end Green's-function FTIR |
| `test_nanocrystal_sparse_hessian.py` | **Active** | Nanocrystal sparse → dense/sparse spectrum |
| `test_hessian_fitting.py` | **Active** | UFF parameter fitting to reference Hessian |
| `test_diatomic_hessian.py` | **Active** | Simple Hessian validation |
| `test_iterative_vibration_solvers.py` | **Active** | Iterative solver comparison |
| `test_vibration_solver_ladder.py` | **Active** | Staged solver ladder M-S0..M-S2 |
| `test_vibration_jacobi_ocl.py` | **Active** | GPU Jacobi one-iteration parity |
| `run_hessian.py` | **Active** | Ad-hoc Hessian exploration |
| `plot_phonon_bands.py` | **Missing** | Referenced in guides; plotting is **inline** in `test_diamond_phonon_bands.py` (backlog: extract) |

Phonon workflow: [`Phonon_testing.guide.md`](../Topics/FTIR_Nanocrystals/Phonon_testing.guide.md).

#### 3e. QM reference (`tests/tSiNCs/`)

| Script | Status | Description |
|--------|--------|-------------|
| `run_vib_spectra.py` | **Active** | Multi-method spectra (DFTB+, CP2K, GPAW, Psi4, PySCF) |
| `plot_hydride_motif_spectra.py` | **Active** | Neighborhood PDOS; `--pyscf` → `OUT_pyscf_jobs/`; `--dftb-l1` → `OUT_dftb_vs_pyscf_l1/` |
| `pyBall/FFfit_plots.py` `plot_stacked_method_pdos` | **Active** | Stacked chemistry PDOS + total DOS + rug; per-row xy inset |
| `pyBall/FFfit_utils.py` `apply_ring_tags` / `heavy_cycles` | **Active** | 5/7-ring exclusive PDOS groups |
| `pyBall/FFfit_utils.py` `load_dftb_vibrations_tag` | **Active** | DFTB+ `vibrations.tag` → signed cm⁻¹ + Cartesian modes |
| `vib_utils.py` | **Active** | ASE dispatch and disk cache |
| `plot_vib_spectra.py` | **Active** | Overlay cached QM spectra |
| `crosscheck_nanocrystal_generators.py` | **Active** | JS vs Python generator parity |

Details: [`tests/tSiNCs/README.md`](../../tests/tSiNCs/README.md).

#### 3g. NPZ pipeline viewers

| File | Status | Description |
|------|--------|-------------|
| `cpp/common/io/NpyIO.h`, `NpzIO.h`, `CrystalNpz.h`, `TopologyNpz.h` | **Active** | C++ NPZ/NPY decode; skips unicode metadata arrays |
| `cpp/apps/MolecularEditor/MolecularBrowser.cpp` | **Active** | SDL browser: `01`/`02`/`03` NPZ, bond color map, AABB overlay |
| `pyBall/io/crystal_npz.py` | **Active** | Python mmap loaders; `is_viewer_listable_basename`, `find_nanocrystal_pipeline_stages` |
| `pyBall/GUI/VispyMolBrowser.py` | **Active** | PyQt5 + Vispy grid browser; plugin host; vibration mode overlay on `AtomScene` |
| `pyBall/GUI/mol_browser_plugins/` | **Active** | Extensible side panels — `VibrationSpectrumPlugin` (FTIR plot + eigenmode arrows/animation) |
| `pyBall/nanocrystal_pipeline.py` | **Active** | `load_spectrum_npz`, `solve_normal_modes_from_hessian_npz` (viewer + pipeline shared) |
| `tests/tSiNCs/test_mol_browser_plugins.py` | **Active** | Plugin registry, NPZ filter, normal-mode load |
| `tests/tSiNCs/run_vispy_mol_browser.sh` | **Active** | Python GUI launcher |
| `tests/tSiNCs/test_cpp_npz_load.sh` | **Active** | Headless C++ verify |
| `tests/tSiNCs/run_cpp_mol_browser.sh` | **Active** | C++ GUI launcher |

Guides: [`CPP_MolecularBrowser_NPZ.md`](../Topics/FTIR_Nanocrystals/CPP_MolecularBrowser_NPZ.md) (C++), [`Python_Vispy_MolBrowser_Plugins.md`](../Topics/FTIR_Nanocrystals/Python_Vispy_MolBrowser_Plugins.md) (Python plugins). Schema: [`NPZ_Crystal_Schema.md`](../Topics/FTIR_Nanocrystals/NPZ_Crystal_Schema.md) v1.2. Pipeline: [`Nanocrystal_NPZ_Pipeline.guide.md`](../Topics/FTIR_Nanocrystals/Nanocrystal_NPZ_Pipeline.guide.md).

**Parity gaps (C++ vs Python):** C++ has no spectrum/Hessian analysis panel; Python has vibration visualization. Topology bonds from `neigh_idx` fully wired in Python; C++ still uses distance fallback on `03_topology.npz`. Automated cross-viewer atom/bbox test pending.

#### 3h. XRD (`pyBall/XRD/`, `tests/tXRD/`)

Powder Debye = sinc transform of the pair-distance histogram (experimental powder / PDF). **Checklist + Kusová strain-gradient science:** [`XRD.md`](XRD.md). Engine writeup: [`XRD_progress.md`](../Topics/FTIR_Nanocrystals/XRD_progress.md). FTIR+XRD briefing: [`SiNCs.md`](SiNCs.md).

| File | Status | Description |
|------|--------|-------------|
| `pyBall/XRD/debye_histogram.py` | **Active** | `XRDDebye` OpenCL; ensemble histogram accumulate; Hessian \(\sigma\) |
| `pyBall/XRD/form_factors.py` | **Active** | Cromer–Mann H, C, Si |
| `cpp/common_resources/cl/XRDDebye.cl` | **Active** | Histogram + Gaussian smear + Debye transform |
| `scripts/generate_xrd_webgl.py` | **Active** | Single-crystal 2D diffraction HTML (GLSL embedded; not powder) |
| `tests/tXRD/test_debye_histogram.py` | **Active** | Static / constant-σ / Hessian-σ (C R=6 sphere) |
| `tests/tXRD/test_ensemble_exact.py` | **Active** | Multi-crystal histogram sum |
| `tests/tXRD/test_large_crystal.py` | **Active** | Larger C diamond NCs |
| `tests/tXRD/README.md` | **Active** | Folder index |

---

### 4. Topic documentation (`doc/Topics/FTIR_Nanocrystals/`)

Full index: [`README.md`](../Topics/FTIR_Nanocrystals/README.md).

| File | Type | Description |
|------|------|-------------|
| `Phonon_testing.guide.md` | guide | Phonon testing workflow and pitfalls |
| `Nanocrystal_NPZ_Pipeline.guide.md` | guide | **NPZ stage contract** (01→05), consumers, batch tutorial |
| `NPZ_Crystal_Schema.md` | contract | Per-key dtype/shape tables for all NPZ kinds |
| `Sparse_Hessian_Vibration_Spectra.guide.md` | guide | Sparse Hessian and spectrum methods |
| `CPP_MolecularBrowser_NPZ.md` | guide | C++ SDL browser: NPZ load, VIEW keys, bond color map |
| `Python_Vispy_MolBrowser_Plugins.md` | guide | Python Vispy browser: plugin system, vibration spectrum panel |
| `HydrideMotif_spectra.md` | report | L2 DFTB+ vs MMFF group PDOS |
| `PySCF_L1_neighborhood_PDOS.md` | report | PySCF L1 stacked PDOS template; vs DFTB+ L1 SK sets; Si imaginaries |
| `gen_nanocrystals.chat.md` | chat | Generation CLI reference |
| `Hessian_Kspace.chat.md` | chat | k-space / Bloch theory |
| `Hessian_fitting.chat.md` | chat | Hessian parameter fitting |
| `Nanocrystal_pipeline.progress.md` | progress | Ensemble pipeline plan (**current direction**) |
| `Nanocrystal_generation.progress.md` | progress | Generation milestones |
| `Linearized_topology.progress.md` | progress | MMFFL linearized topology |
| `Nanocrystal_vibration_sparse.progress.md` | progress | Sparse vibration staging |
| `Sparse_vibration_solver.progress.md` | progress | GPU/iterative solvers (deferred) |
| `XRD_progress.md` | progress | GPU Debye engine + thermal σ. Science checklist: [`XRD.md`](../../topical_audit/XRD.md) |
| `Debug_ASan_double_free_and_eval_atom.progress.md` | progress | MMFF ASan debugging |

**Not present:** `Debug_negative_phonon_freqs.md` — imaginary-mode checklist is in [`Phonon_testing.guide.md`](../Topics/FTIR_Nanocrystals/Phonon_testing.guide.md) (Known Issues).

---

## Relationships

```
Generation                    Force field                    Vibration / spectra
──────────                    ───────────                    ───────────────────
Nanocrystals.js          ──►  MMFF.py::init()           ──►  nanocrystal_pipeline.py  (production: eigh + histogram)
tests/tSiNCs/nanocrystals.mjs   getHessian3Nx3N()              FTIR.vibration_spectrum_from_modes
nanocrystal_gen.py                 │                         test_vibration_spectra.py  (Green's probing)
       │                           ▼
       │                    MMFF_lib.cpp (C++ FD)
       │                           │
       └─ crosscheck ──────────────┼──► test_diamond_phonon_bands.py (phonon bands)
                                  │
tests/tSiNCs/run_vib_spectra.py ◄─┘  (QM reference, adamantane / sila-adamantane)
```

**Cross-language flow:**

1. **JavaScript** (`Nanocrystals.js`, `tests/tSiNCs/nanocrystals.mjs`) → `.mol2` / `.xyz` / NPZ cache
2. **Python** (`MMFF.py`) loads structures, calls C++ via ctypes
3. **C++** (`MMFF_lib.cpp`) computes Hessian / phonon blocks (central FD)
4. **Python** (`FTIR.py`, `nanocrystal_pipeline.py`, `tests/tMMFF/`) → spectra, phonons, fitting

---

## Test matrix

| Test script | What it tests | Force field | PBC | Key flags |
|-------------|---------------|-------------|-----|-----------|
| `test_diamond_phonon_bands.py` | Phonon dispersion | MMFF/UFF | Optional | `--pbc`, `--uff`, `--scale-bond`, `--scale-angle`, `--asr` |
| `test_diamond_gamma.py` | Γ-point frequencies | MMFF | No | — |
| `test_ethane_gamma.py` | Molecular vibrations | MMFF | No | — |
| `test_vibration_spectra.py` | FTIR (Green's probing) | MMFF | No | `--fmin`, `--fmax`, `--eta`, `--dx` |
| `test_nanocrystal_sparse_hessian.py` | NC sparse Hessian spectra | MMFF | No | `--spectrum`, `--plot` |
| `test_hessian_fitting.py` | Parameter fitting | UFF | No | — |
| `test_diatomic_hessian.py` | Hessian correctness | MMFF | No | — |
| `crosscheck_nanocrystal_generators.py` | JS vs Python gen | — | — | (in `tests/tSiNCs/`) |
| `run_vib_spectra.py` | QM reference spectra | QM | No | `--methods` |
| `fit_mmff_kss_pyscf.py` | Own-min stretch RMSE vs PySCF | MMFF | No | `--mol pyscf_nc --ref joint\|nbscan\|morse2d\|anneal` |
| `test_cpp_npz_load.sh` | C++ NPZ parse | — | — | `--gui` optional |

Run classical tests from `tests/tMMFF/` via `bash run.sh <script>`. QM and parity: see [`tests/tSiNCs/README.md`](../../tests/tSiNCs/README.md).

---

## Status summary

| Component | Status | Notes |
|-----------|--------|-------|
| Nanocrystal generation | **Active** | JS primary; Python sphere parity |
| MMFF (C++) | **Active** | CPU FD Hessian is production path |
| UFF (C++) | **Active** | Comparison / fitting; angle issues on diamond (see phonon guide) |
| Hessian computation | **Active** | O(6N) force evals; bottleneck for large NCs |
| Ensemble / NPZ pipeline | **Working** | `nanocrystal_pipeline.py`; `fixtures/si_1nm_passivation/` 01→05 per crystal |
| Hessian (topology-linear) | **Active** | `FTIR.build_hessian_from_linear_topology`; default in pipeline |
| Hessian (MMFF FD) | **Diagnostic** | `--source mmff` / `--compare-mmff` only |
| Production spectra | **Active** | Dense `eigh` + mode histogram |
| Linear-response FTIR | **Active** | Green's probing in `FTIR.py` |
| Phonon bands (PBC) | **Active** | Supercell + Bloch; `checkPBCNeighCells()` |
| QM references | **Active** | Cages: `run_vib_spectra.py`. L1 PySCF npy: `OUT_pyscf_jobs/` — **Si not stationary** |
| Sparse/GPU freq solvers | **Deferred** | Implemented but not on critical path |
| GPU Hessian | **Not started** | MD GPU paths exist; Hessian is CPU-only |
| `plot_phonon_bands.py` | **Missing** | Extract from test script (backlog) |
| `bootstrap_vibration_parallel_fixtures.py` | **Missing** | Needed for some `tMMFF` solver tests (backlog) |
| C++ topology bond sticks from `neigh_idx` | **Missing** | `03_topology.npz` uses weak distance fallback |
| Cross-viewer NPZ parity test | **Pending** | Atom + AABB count automation |

---

## Known limitations

1. **PBC optical modes:** cluster mode gives correct Γ optical modes; finite-k acoustic branches may go imaginary on small clusters. Use `--asr` or PBC with care — [`Phonon_testing.guide.md`](../Topics/FTIR_Nanocrystals/Phonon_testing.guide.md).
2. **Hessian cost:** 6N force evaluations; expensive above ~1000 atoms. No GPU Hessian.
3. **Parameter fitting:** Hessian FFfit is bond/angle stiffness. Own-min stretch RMSE vs PBE L1 (`fit_mmff_kss_pyscf.py`) can also scan surface Morse α and XXX/XXH/HXH; **one C pack did not unify CH vs CH₂** ([`MMFF_C_CH_vs_CH2_kfit.md`](../Topics/FTIR_Nanocrystals/MMFF_C_CH_vs_CH2_kfit.md)).
4. **Surface H:** chem atlas L0/L1 MMFF-relaxes H with surface-only NB ([`ChemAtlas_MMFF_relax.md`](../Topics/FTIR_Nanocrystals/ChemAtlas_MMFF_relax.md)); other generators may still leave static caps. Fitter historically bonded-only; cubes clash without surface NB.
5. **SiO₂ / silica:** out of scope.
6. **Relax vs Hessian FF:** relax uses C++ MMFF today; Hessian uses exported MMFFL linear sticks — LFF relax on `03_topology.npz` is backlog for full parity.
7. **Eigenvectors:** not stored in `05_spectrum.npz` v1.2 (modes used internally for probe weights).
8. **PySCF L1 (2026-09-03):** canonical `/home/prokop/SIMULATIONS/SiNCs/pyscf_vib_results`. **10/12 0-imag** (five Si fixed by tight FMAX). Still saddles: `cube_C`, `octahedron_7ring_Si`. Legacy `pySCF/jobs/results` Si superseded. DFTB+ L1 of the same 12 crystals is 0-imag (`OUT_dftb_vs_pyscf_l1/` needs rebuild vs tight Si). [`PySCF_L1_neighborhood_PDOS.md`](../Topics/FTIR_Nanocrystals/PySCF_L1_neighborhood_PDOS.md).

---

## Open items

Tracked in [`tests/tSiNCs/ToDo_Nanocrystal.md`](../../tests/tSiNCs/ToDo_Nanocrystal.md).

---

## Related topics

- [Topical Audit Index](topical_audit.md) — priority ranking, dependency graph, missing topics
- [Classical Force Fields (MMFF/UFF)](topical_audit.md#priority-2-high--classical-force-fields)
- [Molecular Topology Editors](molecular_topology_editors.md) — crystal building gap analysis (Nanocrystals.js uses CrystalUtils.js)
- [Intramolecular Forcefields](intramolecular_forcefields.md) — UFF/MMFF used for nanocrystal vibrations
- [Htransfer_Kekule_DFTB.md](Htransfer_Kekule_DFTB.md)
- [`doc/Topics/FTIR_Nanocrystals/README.md`](../Topics/FTIR_Nanocrystals/README.md) — topic doc index
- [`tests/tSiNCs/README.md`](../../tests/tSiNCs/README.md) — working hub
- [`CODEMAP.md`](../../CODEMAP.md)
- Codemap: [Silicon/Diamond Nanocrystal Vibration Spectroscopy Workflow](https://windsurf.com/codemaps) — visual trace of the complete pipeline
