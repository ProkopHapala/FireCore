---
type: TopicalAudit
title: Si / Diamond Nanocrystals (FTIR + XRD)
tags: [nanocrystal, ftir, xrd, phonon, mmff, fffit, ensemble]
timestamp: 2026-09-01
---

# Si / Diamond Nanocrystals — Project Briefing (FTIR + XRD)

Self-contained inventory for discussion (paste this file into ChatGPT). Repo: **FireCore**. Working hub: `tests/tSiNCs/`. Canonical scripts live there, not under `scripts/` (legacy copies deprecated).

Related audits (do not duplicate here):
- [`Nanocrystal_Vibrations.md`](Nanocrystal_Vibrations.md) — API-level vibration/phonon inventory
- [`doc/Topics/FTIR_Nanocrystals/README.md`](../Topics/FTIR_Nanocrystals/README.md) — topic doc index
- [`doc/Topics/FFfit/README.md`](../Topics/FFfit/README.md) — Hessian fitting theory
- [`tests/tSiNCs/AGENTS.md`](../../tests/tSiNCs/AGENTS.md) — local DOX contract

---

## 1. Scientific goal

Compute **FTIR / phonon spectra** and **powder XRD** for hydrogen-capped **Si** and **C diamond** nanocrystals, then compare to experiment by **attributing sub-peaks to local chemical neighborhoods** (SiH / SiH₂ / SiH₃; CH on {111}/{110}, CH₂ on {100}, CH₃) and by mixing **ensembles** of shapes/sizes.

Experimental target (K. Kusová plots, symlink `tests/tSiNCs/SiH_CH_FTIRs_Kusova/`):

| Window | Chemistry | Typical assignment |
|--------|-----------|--------------------|
| ~2000 cm⁻¹ | amorphous Si–H | disordered / undercoordinated |
| 2040–2100 | **SiH** | monohydride, often {111} |
| 2100–2140 | **SiH₂** | dihydride, often {100} |
| 2140–2170 | **SiH₃** | corners / edges |
| 2180–2280 | **OxSiHᵧ** | backoxidized — **not implemented** (no oxygen) |
| 2815–2975 | **C–H** on diamond NCs | CH {111}/{110}, CH₂ {100}, CH₃, ridges |

Samples: Linde/Messer/bubbler SiH standards (T-series change **weights**, not new bands); ND set RND-H, MSY, BND, CD, Cheng (`_NDsamples.txt`).

**Raman:** same modes, different intensities — not a priority.

---

## 2. Architecture (what is production vs deferred)

```
CIF (Si-sym / diamond_primitive)
        │
        ▼
JS Nanocrystals.js  ── generate, H-cap, Miller/sphere cuts, optional bridges
        │  mol2 / xyz / 01_crystal.npz + 03_topology.npz
        ▼
C++ MMFF (pyBall.MMFF)  ── relax geometry
        │  02_relaxed.npz
        ├──────────────────────────────┐
        ▼                              ▼
FD Hessian 3N×3N                 Debye powder XRD (GPU)
        │  04_hessian.npz              pair histogram + sinc(Qr)
        ▼                              optional σ_ij from Hessian
dense eigh + histogram                 ensemble = sum histograms
        │  05_spectrum.npz
        ▼
accumulate over crystals  →  FTIR-like envelope  +  powder diffractogram
```

**Production vibration path:** C++ MMFF finite-difference Hessian → `np.linalg.eigh` → mode histogram (`pyBall/nanocrystal_pipeline.py`). Linearized topology (`03_topology.npz`) is **exported but not used** for production spectra (Path B deferred).

**What “FTIR” currently is:** vibrational **DOS** (optionally probe-weighted with dummy `charges=1`). True IR needs dipole derivatives; peak **positions** of typed hydrides are the honest first experimental target.

---

## 3. Solvers and libraries

### 3.1 Structure generation

| Location | Status | Role |
|----------|--------|------|
| `web/molgui_webgpu/Nanocrystals.js` | **Active** | CIF → supercell → sphere/Miller cuts → prune → H-cap → CH₂/SiH₂ bridges |
| `web/molgui_webgpu/NanocrystalExport.js` | **Active** | NPZ bundle export |
| `pyBall/nanocrystal_gen.py` | **Active** | Python sphere cuts; Miller planes delegate to Node |
| `tests/tSiNCs/nanocrystals.mjs` | **Active** | Unified CLI: `generate`, `ensemble`, `topology`, `audit`, `nonbond`, `rings` |
| `tests/tSiNCs/make_small_symmetric_nc.mjs` | **Active** | Deterministic Si/C spheres &lt;100 atoms |
| `tests/tSiNCs/gen_afm_tip.mjs` | **Active** | Truncated tetrahedron {111} tip (AFM, not FTIR) |
| `tests/tSiNCs/gen_nanocrystals.mjs/.py` | **Deprecated** | Use `nanocrystals.mjs generate` |

Passivation knobs: `caps: H`, `capHHBonds`, `resolveClashes`, `insertProb` / `collapseProb` (surface bridges). **No oxygen.**

### 3.2 Force field / Hessian (C++ via Python)

| Location | Status | Role |
|----------|--------|------|
| `cpp/libs/Molecular/MMFF_lib.cpp` | **Active** | `getHessian3Nx3N` (central FD), sparse 3×3 blocks, phonon Φ |
| `cpp/common/molecular/MMFFsp3_loc.h` | **Active** | Bond/angle forces |
| `cpp/common/molecular/NBFF.h` | **Active** | Nonbond + PBC neighbor check |
| `pyBall/MMFF.py` | **Active** | ctypes: Hessian, `setBondParamsByType`, `setAngleParamsByType` |
| `pyBall/FTIR.py` | **Active** | Spectra, Green's probing, topology Hessian, rigid-mode projection |

**Nonbond policy (open, recommended):** no LJ/Coulomb inside the covalent Si/C network (1–2–3 exclusions cover bulk). Surface **H···H steric** needed — angle terms + optional `capHHBonds` / collision groups; do **not** enable full-crystal nonbond.

### 3.3 Vibration / FTIR solvers

| Method | Where | Status | Use |
|--------|-------|--------|-----|
| Dense `eigh` + histogram | `nanocrystal_pipeline.py`, `FTIR.vibration_spectrum_from_modes` | **Production** | Ensemble spectra, N ≲ 1500 atoms |
| Green's frequency probing | `FTIR.mechanical_greens_probing` | **Active** | Linear-response FTIR; local-bond probes without full eig |
| Topology-linear Hessian | `FTIR.build_hessian_from_linear_topology` | **Built, unused in prod** | Path B from `03_topology.npz` springs |
| Sparse Green's | `mechanical_greens_probing_sparse` | **Deferred** | Large N |
| GPU Jacobi | `pyBall/OCL/VibSolver.py`, `cl/vib_jacobi.cl` | **Deferred** | `eigh` is not the bottleneck |
| PBC phonon bands | `getPhononPhiBlocks` + Bloch | **Active** | Bulk diamond validation, not NC FTIR |
| Raman | — | **Not started** | Polarizability |

Hessian cost is **O(6N) force evals** (CPU). No GPU Hessian.

### 3.4 Force-field fitting (stiffness constants)

Fit classical \(k\) to **PySCF Hessians** so Si–H / C–H land in experimental windows. Hybrid objective: modes + local Hessian + Wilson row-space.

| Location | Role |
|----------|------|
| `cpp/common/molecular/FFfit.h`, `cpp/libs/Molecular/FFfit_lib.cpp` | Analytic sensitivity, LSQ, GD |
| `pyBall/FFfit.py`, `FFfit_utils.py`, `FFfit_plots.py` | Python wrap + SiH/SiH₂/SiH₃ typing |
| `tests/tSiNCs/test_FFfit.py` | Multi-system PySCF fit CLI |
| `assign_si_environment_types()` | Classify Si by bonded-H count |

**Current all-Si result** (`tests/tSiNCs/SiNCs_FFfit_summary.md`): hierarchical Si subtypes cut mean frequency RMSE ~41 → ~32 cm⁻¹. Stretch–bend couplings help the Hessian; torsions are **not** recommended for tetrahedral Si/C.

### 3.5 Powder XRD

Debye scattering on a **finite cluster** (no PBC):

\[
I(Q)=\sum_{i,j} f_i(Q)f_j(Q)\,\mathrm{sinc}(Q r_{ij})
\]

Implemented as pair-distance **histogram** then sinc transform (GPU). Thermal smear: Gaussian in \(r\) with \(\sigma_{ij}=\sqrt{k_B T/k_{ij}}\) from local Hessian blocks (frozen-stiffness lower bound). **Static distortion** is automatic: histogram uses **relaxed** coordinates vs bulk lattice.

| Location | Role |
|----------|------|
| `pyBall/XRD/debye_histogram.py` | `XRDDebye` OpenCL engine, ensemble accumulate, σ from Hessian |
| `pyBall/XRD/form_factors.py` | Cromer–Mann H, C, Si |
| `cpp/common_resources/cl/XRDDebye.cl` | `pair_histogram`, `pair_histogram_gaussian`, `debye_transform_q` |
| `scripts/generate_xrd_webgl.py` | Single-crystal 2D diffraction HTML viewer |

Ensemble XRD = **sum pair histograms** over crystals, one Debye transform (size/shape polydispersity).

---

## 4. NPZ pipeline contract (01→05)

Guide: [`Nanocrystal_NPZ_Pipeline.guide.md`](../Topics/FTIR_Nanocrystals/Nanocrystal_NPZ_Pipeline.guide.md). Schema: [`NPZ_Crystal_Schema.md`](../Topics/FTIR_Nanocrystals/NPZ_Crystal_Schema.md) v1.2.

| File | Content |
|------|---------|
| `01_crystal.npz` | Generated positions, Z, bonds |
| `02_relaxed.npz` | MMFF-relaxed positions |
| `03_topology.npz` | Linearized springs K₁₂/K₁₃/K₁₄ (viewer + Path B) |
| `04_hessian.npz` | 3N×3N Hessian |
| `05_spectrum.npz` | Mode frequencies + histogram (eigenvectors **not** stored) |

Orchestration: `python3 -m pyBall.nanocrystal_pipeline {relax,hessian,spectrum,accumulate}` spawned from `nanocrystals.mjs ensemble`.

---

## 5. Test scripts and CLIs

Canonical: run Node from **repo root**; Python tests in `tests/tSiNCs` or `tests/tMMFF` via `run.sh`.

### 5.1 `tests/tSiNCs/` — nanocrystal hub

| Script | Role |
|--------|------|
| `nanocrystals.mjs` | Unified CLI (generate / ensemble / topology / audit / nonbond / rings) |
| `make_small_symmetric_nc.mjs` + `small_symmetric_nc.json` | Si/C spheres &lt;100 atoms → `OUT_small_nc/` |
| `ensemble.example.json` | 8 C-diamond spheres (example ensemble) |
| `atlas_shapes.json` | Fixed shapes: spheres + {111}/{100}/{110} mixes → `OUT_nc_atlas/` |
| `export_nanocrystal_bundle.mjs` | `01` + `03` NPZ for fixtures |
| `crosscheck_nanocrystal_generators.py` | JS vs Python sphere parity |
| `test_nanocrystal_geometry.mjs` | Passivation audit G0–G5 |
| `test_nanocrystal_defect_export.mjs` | Defect export |
| `generate_si_passivation_examples.mjs` | Passivation gallery |
| `run_vib_spectra.py`, `vib_utils.py`, `plot_vib_spectra.py` | QM refs (DFTB+/CP2K/GPAW/Psi4/PySCF) on adamantane, sila-adamantane |
| `run_small_np_pyscf_vib.py`, `gen_small_nps.py` | Small-NP PySCF vibrations |
| `test_FFfit.py`, `test_fffit_hybrid.py` | Hessian FF fit |
| `test_parity_py_cpp.py`, `test_parity_graph_cpp.py` | FFfit Python/C++ parity |
| `bootstrap_vibration_parallel_fixtures.py` | Solver-ladder fixtures |
| `export_vibration_benchmarks.py`, `plot_vibration_benchmark_sparsity.py` | Benchmark NPZ |
| `run_cpp_mol_browser.sh`, `run_vispy_mol_browser.sh` | NPZ viewers |
| `test_cpp_npz_load.sh`, `test_crystal_npz_load.py`, `test_mol_browser_plugins.py` | Viewer/NPZ tests |

Deprecated wrappers (same folder): `run_nanocrystal_ensemble.mjs`, `gen_nanocrystals.mjs`, `build_linearized_topology.mjs`, `debug_nanocrystal_nonbond_groups.mjs`.

### 5.2 `tests/tMMFF/` — classical vibration / phonon

Run: `cd tests/tMMFF && bash run.sh <script>`.

| Script | Role |
|--------|------|
| `test_vibration_spectra.py` | Green's-function FTIR end-to-end |
| `test_diamond_phonon_bands.py` | Diamond dispersion; `--pbc --asr --uff` |
| `test_diamond_gamma.py`, `test_ethane_gamma.py` | Γ-point checks |
| `test_nanocrystal_sparse_hessian.py` | Sparse Hessian → spectrum |
| `test_hessian_fitting.py` | UFF k fit to Hessian |
| `test_diatomic_hessian.py` | FD Hessian sanity |
| `test_iterative_vibration_solvers.py`, `test_vibration_solver_ladder.py` | Solver ladder |
| `test_vibration_jacobi_ocl.py` | GPU Jacobi parity |

### 5.3 `tests/tXRD/` — powder diffraction

| Script | Role |
|--------|------|
| `test_debye_histogram.py` | Static / constant-σ / Hessian-σ comparison |
| `test_ensemble_exact.py` | Histogram accumulation over copies |
| `test_large_crystal.py` | R=8, R=10 C diamond generate–relax–Hessian–XRD |

### 5.4 Fixtures (committed vs generated)

| Path | Role |
|------|------|
| `fixtures/si_1nm_passivation/` | Nine Si ~1 nm crystals, full 01→05 |
| `fixtures/npz_viewer/` | Tiny viewer smoke NPZ |
| `fixtures/vibration_benchmarks/` | Size-scaling Hessian NPZ |
| `fixtures/vibration_parallel/` | Solver-ladder structures |
| `OUT_*` | **Generated** — do not commit (`OUT_small_nc`, `OUT_nc_ensemble_v2`, `OUT_nc_atlas`, `OUT_XRD`, `OUT_FFfit_plots`) |

---

## 6. Markdown documents

### 6.1 Topic folder `doc/Topics/FTIR_Nanocrystals/`

| File | Kind | Content |
|------|------|---------|
| `README.md` | index | Topic map |
| `Nanocrystal_NPZ_Pipeline.guide.md` | guide | 01→05 contract |
| `NPZ_Crystal_Schema.md` | contract | dtypes/shapes |
| `Phonon_testing.guide.md` | guide | PBC/ASR pitfalls |
| `Sparse_Hessian_Vibration_Spectra.guide.md` | guide | Sparse Hessian methods |
| `Vibration_benchmark_npz.guide.md` | guide | Timing vs N |
| `CPP_MolecularBrowser_NPZ.md` | guide | C++ SDL viewer |
| `Python_Vispy_MolBrowser_Plugins.md` | guide | Vispy + vibration plugin |
| `gen_nanocrystals.chat.md` | chat | Generation CLI |
| `Hessian_Kspace.chat.md` | chat | Bloch Hessian |
| `Hessian_fitting.chat.md` | chat | Older fitting discussion (prefer FFfit theory) |
| `XrayDiffraction.chat.md` | chat | Debye + thermal σ design |
| `Nanocrystal_pipeline.progress.md` | progress | Ensemble milestones |
| `Nanocrystal_generation.progress.md` | progress | G0–G5 generation |
| `Linearized_topology.progress.md` | progress | MMFFL springs |
| `Nanocrystal_vibration_sparse.progress.md` | progress | Sparse staging |
| `Sparse_vibration_solver.progress.md` | progress | GPU solvers (**de-prioritized**) |
| `XRD_progress.md` | progress | XRD implementation + user guide |
| `Debug_ASan*.md`, `Debug_negative_phonon_freqs.md` | debug | MMFF/phonon debugging |

### 6.2 Force-field fit `doc/Topics/FFfit/`

| File | Content |
|------|---------|
| `HessianFitting_Theory.md` | **Authoritative** Wilson / hybrid / subtype hierarchy |
| `Dihedral_Torsion.md` | Torsions implemented; not production for Si/C |
| `HessianFit.chat.md` | Chronological log (may be stale vs theory) |
| `tests/tSiNCs/SiNCs_FFfit_summary.md` | Latest all-Si numbers |

### 6.3 Hub docs in `tests/tSiNCs/`

`README.md`, `AGENTS.md`, `ToDo_Nanocrystal.md`, `small_symmetric_nc.md`, fixture READMEs under `fixtures/*/`.

---

## 7. Status vs experiment

| Need | Have | Gap |
|------|------|-----|
| Generate faceted H-capped Si/C NCs | Yes (Miller + sphere + atlas) | Si atlas JSON (example atlas is C diamond) |
| Relax + Hessian + spectrum | Yes (NPZ 01→05, ensemble accumulate) | Spectra are DOS, not IR intensities |
| Split SiH / SiH₂ / SiH₃ frequencies | Typing + FFfit **exists** | Not yet plotted as tagged sub-spectra vs Kusová |
| Facet-typed C–H (111/100/110) | Generator can cut facets | No CH/CH₂/CH₃ FFfit types wired like Si subtypes |
| Ensemble mixing \(I=\sum w_c I_c\) | `accumulate_spectra` sums histograms | No chemistry tags; no fit of \(w_c\) to experiment |
| Powder XRD + thermal + distortion | Debye GPU + Hessian σ + relaxed coords | Not run on Si FTIR atlas yet |
| Surface H steric | Angles + optional `capHHBonds` | Full-crystal nonbond should stay **off** |
| Backoxidation 2180–2280 cm⁻¹ | — | **Out of scope** until oxygen in gen+FF |
| Raman | — | Not started |
| sp² carbon dots (CD sample) | — | Generator is sp³ only |

**Known physics limits:** surface H historically generated then not always re-relaxed before vibration; PBC optical modes need `--asr`; Hessian FD expensive above ~1000 atoms; eigenvectors omitted from `05_spectrum.npz`.

---

## 8. Suggested next experiments (discussion fodder)

1. **Calibrate \(k\)** — `test_FFfit.py --si-subtypes` on SiH₄ + small NCs until SiH/SiH₂/SiH₃ stretch order matches ~2080 / 2120 / 2150 cm⁻¹.
2. **Chemical atlas** — three Si crystals ({111}, {100}, sphere, 50–80 atoms): relax, Hessian, histogram **per hydride type** (reuse `assign_si_environment_types` on eigenvectors). Overlay on Kusová Si–H window.
3. **Same crystals → Debye XRD** (static relaxed + Hessian σ).
4. **Mix weights \(w_c\)** to mimic Linde vs Messer vs T-series (intensity ratios, not new peaks).
5. Repeat 2–4 for diamond C–H once CH typing exists.
6. **Do not** start with N=100 random spheres until typed sub-peaks exist.

Green's probing is the right tool to drive only one hydride class without a huge ensemble.

---

## 9. How to run (cheat sheet)

```bash
# Small symmetric Si + C diamond
node tests/tSiNCs/make_small_symmetric_nc.mjs

# Shape atlas (C diamond today)
node tests/tSiNCs/nanocrystals.mjs ensemble --atlas tests/tSiNCs/atlas_shapes.json \
  --output-dir tests/tSiNCs/OUT_nc_atlas

# Random ensemble example
node tests/tSiNCs/nanocrystals.mjs ensemble --config tests/tSiNCs/ensemble.example.json

# Single-crystal NPZ stages
python3 -m pyBall.nanocrystal_pipeline relax --init-npz DIR/01_crystal.npz --out-npz DIR/02_relaxed.npz
python3 -m pyBall.nanocrystal_pipeline hessian --relaxed-npz DIR/02_relaxed.npz --out-npz DIR/04_hessian.npz
python3 -m pyBall.nanocrystal_pipeline spectrum --hessian-npz DIR/04_hessian.npz --out-npz DIR/05_spectrum.npz

# FFfit
CPP_BUILD_PATH=$PWD/cpp/Build-opt/libs python3 tests/tSiNCs/test_FFfit.py --cases all_Si --compare-typing

# Phonon / Green's FTIR
cd tests/tMMFF && bash run.sh test_vibration_spectra.py
cd tests/tMMFF && bash run.sh test_diamond_phonon_bands.py --unit THz --super-n 3

# Powder XRD
python3 tests/tXRD/test_debye_histogram.py
```

---

## 10. Open items

Tracked in [`tests/tSiNCs/ToDo_Nanocrystal.md`](../../tests/tSiNCs/ToDo_Nanocrystal.md). Experiment-facing gaps: hydride-tagged spectra, Si atlas, IR intensities, oxygen, Raman, CH facet typing, \(w_c\) fit to Kusová envelopes.
