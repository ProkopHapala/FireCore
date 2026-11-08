---
type: TopicalAudit
title: Si / Diamond Nanocrystals (FTIR + XRD)
tags: [nanocrystal, ftir, xrd, phonon, mmff, fffit, ensemble]
timestamp: 2026-09-03
---

# Si / Diamond Nanocrystals — Project Briefing (FTIR + XRD)

Self-contained inventory for discussion (paste this file into ChatGPT). Repo: **FireCore**. Working hub: `tests/tSiNCs/`. Canonical scripts live there, not under `scripts/` (legacy copies deprecated).

Related audits (do not duplicate here):
- [`XRD.md`](XRD.md) — **Powder XRD / PDF inventory + Kusová strain-gradient checklist**
- [`Nanocrystal_Vibrations.md`](Nanocrystal_Vibrations.md) — API-level vibration/phonon inventory
- [`doc/Topics/FTIR_Nanocrystals/README.md`](../Topics/FTIR_Nanocrystals/README.md) — topic doc index
- [`doc/Topics/FTIR_Nanocrystals/PySCF_L1_neighborhood_PDOS.md`](../Topics/FTIR_Nanocrystals/PySCF_L1_neighborhood_PDOS.md) — PySCF L1 stacked PDOS (plot template; vs DFTB+ L1; Si PySCF not a minimum)
- [`doc/Topics/FTIR_Nanocrystals/ChemAtlas_MMFF_relax.md`](../Topics/FTIR_Nanocrystals/ChemAtlas_MMFF_relax.md) — Wulff atlas output paths + surface-NB policy
- [`doc/Topics/FFfit/README.md`](../Topics/FFfit/README.md) — Hessian fitting theory
- [`doc/topical_audit/Hessian_fitting.md`](Hessian_fitting.md) — FFfit / Wilson / own-min \(k\) inventory + method checklist
- [`doc/Topics/FTIR_Nanocrystals/MMFF_C_CH_vs_CH2_kfit.md`](../Topics/FTIR_Nanocrystals/MMFF_C_CH_vs_CH2_kfit.md) — C PBE L1 own-min \(k\)-fit: CH vs CH₂ + surface Morse; **not** a transferable pack
- [`tests/tSiNCs/AGENTS.md`](../../tests/tSiNCs/AGENTS.md) — local DOX contract

---

## 1. Scientific goal

Compute **FTIR / phonon spectra** and **powder XRD** for hydrogen-capped **Si** and **C diamond** nanocrystals.

- **FTIR:** attribute sub-peaks to local chemistry (SiH / SiH₂ / SiH₃; CH on {111}/{110}, CH₂ on {100}) and mix ensemble weights.
- **XRD (Kusová 2026-09-03):** powder \(I(Q)\) is the Debye transform of the pair-distance histogram. They fit an *effective crystalline size*. That is **not** a discrete crystal core in amorphous stuff — it is a crystalline core whose bond lengths **gradually** distort toward the passivated/defective surface. Correlate that fit with maps of \(a(\mathbf{r})\), the PDF, and \(I(Q)\) on **large** (L2/L3) relaxed crystals. Full checklist: [`XRD.md`](XRD.md).

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
| `tests/tSiNCs/chem_atlas.json` | **Active** | Wulff Si ≡ C atlas + L1 5/7-ring defects; `relax_nonbond: surface` |
| `tests/tSiNCs/nanocrystals.mjs` | **Active** | Unified CLI: `generate`, `ensemble` (`--atlas`), `topology`, `audit`, `nonbond`, `rings` |
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
| `pyBall/MMFF.py` | **Active** | ctypes: Hessian, `setExclusion2`, `setBondParamsByType`, `setAngleParamsByType` |
| `pyBall/nanocrystal_pipeline.py` | **Active** | `relax --nonbond surface`: REQs=0 on 4-coordinated X + Exclusion2; atlas L0/L1 |
| `pyBall/FTIR.py` | **Active** | Spectra, Green's probing, topology Hessian, rigid-mode projection |

**Nonbond policy (atlas L0/L1, implemented):** no LJ inside the covalent Si/C network. Surface set = H + heavies with nHeavy&lt;4 (`CollisionWorkgroups.surfaceAtomIndices`); bulk 4-coordinated X has `REQs=0`. 1–2 and 1–3 skipped (`setExclusion2`) so geminal XH₂ is angle-only. Do **not** enable full-crystal nonbond. Do **not** write H–H as covalent bonds (`recalculateBonds` skip + `stripHHBonds`). Report: [`ChemAtlas_MMFF_relax.md`](../Topics/FTIR_Nanocrystals/ChemAtlas_MMFF_relax.md).

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
| `tests/tSiNCs/test_FFfit.py` | Multi-system PySCF Hessian-fit CLI |
| `tests/tSiNCs/fit_mmff_kss_pyscf.py` | Own-min stretch RMSE vs PBE L1 (`joint`/`nbscan`/`morse2d`/`anneal`) |
| `assign_si_environment_types()` | Classify Si by bonded-H count |
| `tests/tSiNCs/plot_hydride_motif_spectra.py` | Neighborhood PDOS orchestration (L1/L2 MMFF+DFTB; `--pyscf` L1 DFT) |
| `pyBall/FFfit_plots.plot_stacked_method_pdos` | **Plot template:** stacked chemistry PDOS + total DOS + rug + optional xy inset |
| `pyBall/FFfit_utils.apply_ring_tags` | 5/7-ring atoms (from `heavy_cycles`) steal hydride tags so \(\sum w=1\) |

**Current all-Si result** (`tests/tSiNCs/SiNCs_FFfit_summary.md`): hierarchical Si subtypes cut mean frequency RMSE ~41 → ~32 cm⁻¹. Stretch–bend couplings help the Hessian; torsions are **not** recommended for tetrahedral Si/C. Full method checklist (restricted Wilson next): [`Hessian_fitting.md`](Hessian_fitting.md).

**Own-min stretch RMSE vs PBE L1** (`tests/tSiNCs/fit_mmff_kss_pyscf.py`, not Hessian FFfit): one C pack \(k_{\mathrm{XH}}=1.775\) is ~83 cm⁻¹ mean (bonded-only). Surface H···H + Morse α + split XXX/XXH/HXH **exist** (`setMorseNonBond`, `setEachAngle`); SA on cube vs octahedron did **not** unify (~72 cm⁻¹ mean, see-saw). Report: [`MMFF_C_CH_vs_CH2_kfit.md`](../Topics/FTIR_Nanocrystals/MMFF_C_CH_vs_CH2_kfit.md). Do not keep fitting Si PBE L1 as FTIR.

### 3.5 Powder XRD

**SSOT (files, checkboxes, Kusová science):** [`XRD.md`](XRD.md).

Debye scattering on a **finite cluster** (no PBC) is the sine transform of the pair histogram — the same object as the experimental powder / PDF:

\[
I(Q)=\sum_{i,j} f_i(Q)f_j(Q)\,\mathrm{sinc}(Q r_{ij})
\]

GPU: `pyBall/XRD/` + `cpp/common_resources/cl/XRDDebye.cl`. Thermal Gaussian in \(r\) from Hessian blocks. **Static strain is automatic on relaxed coords.** Engine demos exist on C diamond spheres. Atlas L2/L3 xyz exist but are **unrelaxed** (no strain gradient yet). L1 is too small for powder. Strain-gradient maps, Scherrer-vs-structure, and experimental \(2\theta\) overlay are **[ ]**.

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
| `plot_hydride_motif_spectra.py` | Neighborhood PDOS: L2 DFTB+/MMFF; PySCF L1 `--pyscf` → `OUT_pyscf_jobs/` |
| `test_FFfit.py`, `test_fffit_hybrid.py` | Hessian FF fit |
| `fit_mmff_kss_pyscf.py` | Own-min \(k\) vs PBE L1 (CH vs CH₂; not a shared pack) |
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

Folder index: [`tests/tXRD/README.md`](../../tests/tXRD/README.md). Checklist: [`XRD.md`](XRD.md).

| Script | Role |
|--------|------|
| `test_debye_histogram.py` | Static / constant-σ / Hessian-σ comparison (C R=6 sphere) |
| `test_ensemble_exact.py` | Histogram accumulation over copies |
| `test_large_crystal.py` | R=8, R=10 C diamond generate–relax–Hessian–XRD |

### 5.4 Fixtures (committed vs generated)

| Path | Role |
|------|------|
| `fixtures/si_1nm_passivation/` | Nine Si ~1 nm crystals, full 01→05 |
| `fixtures/npz_viewer/` | Tiny viewer smoke NPZ |
| `fixtures/vibration_benchmarks/` | Size-scaling Hessian NPZ |
| `fixtures/vibration_parallel/` | Solver-ladder structures |
| `OUT_*` | **Generated** — do not commit (`OUT_small_nc`, `OUT_nc_ensemble_v2`, `OUT_nc_atlas`, `OUT_XRD`, `OUT_FFfit_plots`, `OUT_mmff_kss_fit`, `OUT_pyscf_jobs`) |

---

## 6. Markdown documents

### 6.1 Topic folder `doc/Topics/FTIR_Nanocrystals/`

| File | Kind | Content |
|------|------|---------|
| `README.md` | index | Topic map |
| `HydrideMotif_spectra.md` | report | L2 DFTB+ vs MMFF group PDOS |
| `PySCF_L1_neighborhood_PDOS.md` | report | PySCF PBE L1 stacked PDOS — **plot template**; vs DFTB+ L1 SK sets; all PySCF Si jobs imaginary |
| `MMFF_stiffness_scaling.md` | SSOT | `bKs` / `apars[:,2]`; do not FIRE with `apars[:,1]×0.30` |
| `MMFF_C_CH_vs_CH2_kfit.md` | report | C L1 cube vs octa: Morse + split angles; SA did not unify (~72 cm⁻¹) |
| `Hessian_at_own_minimum.md` | HARD | Spectrum only at that method’s own minimum |
| `ChemAtlas.plan.md` | plan | Wulff motif atlas vs CompChemUtils |
| `ChemAtlas_MMFF_relax.md` | report | Where `OUT_chem_atlas` is; collision-group MMFF |
| `Nanocrystal_NPZ_Pipeline.guide.md` | guide | 01→05 contract |
| `NPZ_Crystal_Schema.md` | contract | dtypes/shapes |
| `Phonon_testing.guide.md` | guide | PBC/ASR pitfalls |
| `Sparse_Hessian_Vibration_Spectra.guide.md` | guide | Sparse Hessian methods |
| `Vibration_benchmark_npz.guide.md` | guide | Timing vs N |
| `CPP_MolecularBrowser_NPZ.md` | guide | C++ SDL viewer |
| `Python_Vispy_MolBrowser_Plugins.md` | guide | Vispy + vibration plugin |
| `gen_nanocrystals.chat.md` | chat | Generation CLI |
| `Hessian_Kspace.chat.md` | chat | Bloch Hessian |
| `Hessian_fitting.chat.md` | chat | GPT 5.6 2026-09-03 restricted Wilson / locality; older bond-fit notes. Checklist: [`Hessian_fitting.md`](../../topical_audit/Hessian_fitting.md) |
| `XrayDiffraction.chat.md` | chat | Debye + thermal σ design. Checklist: [`XRD.md`](XRD.md) |
| `Nanocrystal_pipeline.progress.md` | progress | Ensemble milestones |
| `Nanocrystal_generation.progress.md` | progress | G0–G5 generation |
| `Linearized_topology.progress.md` | progress | MMFFL springs |
| `Nanocrystal_vibration_sparse.progress.md` | progress | Sparse staging |
| `Sparse_vibration_solver.progress.md` | progress | GPU solvers (**de-prioritized**) |
| `XRD_progress.md` | progress | GPU Debye engine (2026-06). Science: [`XRD.md`](XRD.md) |
| `Debug_ASan*.md`, `Debug_negative_phonon_freqs.md` | debug | MMFF/phonon debugging |

### 6.2 Force-field fit `doc/Topics/FFfit/`

| File | Content |
|------|---------|
| [`Hessian_fitting.md`](Hessian_fitting.md) | **Inventory + checklist** |
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
| Split SiH / SiH₂ / SiH₃ frequencies | Typing + FFfit **exists**; PySCF L1 PDOS **plotted** | **10/12 PBE L1 at min** (`pyscf_vib_results`). Still saddles: `cube_C`, `octahedron_7ring_Si`. Plots are DOS, not IR |
| Facet-typed C–H (111/100/110) | Wulff tags + stacked PDOS on L1 PySCF (0-imag cases except `cube_C`) | One L1 particle, no ensemble \(w_c\), no IR intensity; `cube_C` has 1 imag |
| Ensemble mixing \(I=\sum w_c I_c\) | `accumulate_spectra` sums histograms | No chemistry tags; no fit of \(w_c\) to experiment |
| Powder XRD + thermal + distortion | Debye GPU + Hessian σ; L1 atlas Debye preview; L2/L3 xyz | L2/L3 **unrelaxed**; no \(a(r)\) maps; no Scherrer-vs-structure; no exp. \(2\theta\) file. [`XRD.md`](XRD.md) |
| Surface H steric | Atlas `--nonbond surface`; fitter `enable_surface_nonbond` + optional Morse | Helps cubes, hurts octahedra; one \(k_{\mathrm{XH}}\) still fails ([`MMFF_C_CH_vs_CH2_kfit.md`](../Topics/FTIR_Nanocrystals/MMFF_C_CH_vs_CH2_kfit.md)) |
| Backoxidation 2180–2280 cm⁻¹ | — | **Out of scope** until oxygen in gen+FF |
| Raman | — | Not started |
| sp² carbon dots (CD sample) | — | Generator is sp³ only |

**Known physics limits:** atlas L0/L1 H **is** MMFF-relaxed (surface NB); other ensembles may still use static caps; PBC optical modes need `--asr`; Hessian FD expensive above ~1000 atoms; eigenvectors omitted from `05_spectrum.npz`.

---

## 8. Suggested next experiments (discussion fodder)

1. **Finish PySCF mode-follow** for `cube_C` and `octahedron_7ring_Si`; re-plot `OUT_pyscf_jobs/stacked_*.png` from `pyscf_vib_results`. Five Si L1 are already 0-imag — overlay those on Kusová only as DOS, not IR.
2. **Calibrate \(k\)** — `test_FFfit.py --si-subtypes` on SiH₄ + small NCs until SiH/SiH₂/SiH₃ stretch order matches ~2080 / 2120 / 2150 cm⁻¹. For C, next lever is **typed CH vs CH₂ \(k\)**, not more SA on five shared scalars ([`MMFF_C_CH_vs_CH2_kfit.md`](../Topics/FTIR_Nanocrystals/MMFF_C_CH_vs_CH2_kfit.md)).
3. **Chemical atlas** — three Si crystals ({111}, {100}, sphere, 50–80 atoms): relax, Hessian, histogram **per hydride type** (reuse `assign_si_environment_types` on eigenvectors). Overlay on Kusová Si–H window **only** from stationary Hessians. Use the stacked PDOS template ([`PySCF_L1_neighborhood_PDOS.md`](../Topics/FTIR_Nanocrystals/PySCF_L1_neighborhood_PDOS.md)).
4. **XRD strain-gradient product** — relax L2/L3 Si, map \(a(r)\), PDF, \(I(Q)\), Scherrer vs geometric size ([`XRD.md`](XRD.md) § Next compute). Do not use L1 for powder overlay.
5. **Mix weights \(w_c\)** to mimic Linde vs Messer vs T-series (intensity ratios, not new peaks).
6. Repeat 2–5 for diamond C–H once CH typing exists.
7. **Do not** start with N=100 random spheres until typed sub-peaks exist.

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
python3 tests/tXRD/test_debye_histogram.py   # engine demo; science checklist: doc/topical_audit/XRD.md
```

---

## 10. Open items

Tracked in [`tests/tSiNCs/ToDo_Nanocrystal.md`](../../tests/tSiNCs/ToDo_Nanocrystal.md). Experiment-facing gaps: **PySCF remaining saddles**, hydride-tagged **IR**, Kusová FTIR overlay, **XRD strain-gradient maps on relaxed L2/L3** ([`XRD.md`](XRD.md)), oxygen, Raman, \(w_c\) fit.
