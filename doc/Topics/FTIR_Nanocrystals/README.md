# FTIR Nanocrystals — Topic Documentation

Design notes, guides, and progress logs for Si / diamond nanocrystal generation, MMFF/UFF vibrations, FTIR, phonons, and related XRD.

## Documentation map

| Doc | Role |
|-----|------|
| [`tests/tSiNCs/README.md`](../../../tests/tSiNCs/README.md) | **Working hub** — commands, fixtures, viewers |
| [`tests/tSiNCs/AGENTS.md`](../../../tests/tSiNCs/AGENTS.md) | **DOX contract** — ownership, paths, verification |
| [`ChemAtlas.plan.md`](ChemAtlas.plan.md) | **Working plan** — motif atlas, FireCore vs CompChemUtils, tool runbook |
| [`ChemAtlas_status.md`](ChemAtlas_status.md) | **Status + checklist + browse bank** (what exists vs planned; images for Katerina) |
| [`HydrideMotif_spectra.md`](HydrideMotif_spectra.md) | **DFTB+ vs MMFF group PDOS** on 3 L2 C crystals; L1 {100}/{111}, XH/XH₂, 5/7-ring |
| [`PySCF_L1_neighborhood_PDOS.md`](PySCF_L1_neighborhood_PDOS.md) | **PySCF PBE L1 stacked PDOS** — plot **template**; vs DFTB+ L1; 10/12 0-imag (`pyscf_vib_results`) |
| [`browse_for_Katerina.html`](browse_for_Katerina.html) | Local HTML gallery of morphology / experiment / existing spectra |
| [`Hessian_at_own_minimum.md`](Hessian_at_own_minimum.md) | **HARD RULE** — spectrum only at that method’s own minimum |
| [`MMFF_stiffness_scaling.md`](MMFF_stiffness_scaling.md) | Detective: CH×1.995 is real \(k\); angle×0.30 hit \(\sin(\theta_0/2)\) not \(K\) |
| [`MMFF_C_CH_vs_CH2_kfit.md`](MMFF_C_CH_vs_CH2_kfit.md) | **C L1 k-fit report** — cube CH₂ vs octa CH; surface Morse/LJ + split angles; SA ~72 cm⁻¹, not a shared FF |
| [`../../topical_audit/Hessian_fitting.md`](../../topical_audit/Hessian_fitting.md) | **FFfit / Wilson audit** — method checklist; next: restricted \(B_H\) |
| [`DFTB_vs_MMFF_L2.html`](DFTB_vs_MMFF_L2.html) | Didactic DFTB+ 3ob-3-1 vs MMFF (3 L2 diamond NCs; DOS not IR intensity) |
| [`ChemAtlas_MMFF_relax.md`](ChemAtlas_MMFF_relax.md) | **Relax report** — where `OUT_chem_atlas` files are, surface-NB policy, last-run fmax |
| [`NPZ_Crystal_Schema.md`](NPZ_Crystal_Schema.md) | **Schema contract** — per-key dtype/shape tables |
| [`doc/topical_audit/SiNCs.md`](../../topical_audit/SiNCs.md) | **Project briefing** — FTIR + XRD + FFfit (ChatGPT-ready) |
| [`doc/topical_audit/XRD.md`](../../topical_audit/XRD.md) | **Powder XRD / PDF** — Debye inventory, Kusová strain-gradient checklist |
| [`doc/topical_audit/Nanocrystal_Vibrations.md`](../../topical_audit/Nanocrystal_Vibrations.md) | **Audit** — vibration API inventory, test matrix |
| [`tests/tSiNCs/ToDo_Nanocrystal.md`](../../../tests/tSiNCs/ToDo_Nanocrystal.md) | Open items and open questions |

When filenames disagree, prefer **schema + pipeline guide** and on-disk paths.

---

## Integrated system (summary)

```mermaid
flowchart LR
  EXP[export_nanocrystal_bundle.mjs] --> C01[01_crystal.npz]
  EXP --> C03[03_topology.npz]
  C01 --> VIEW[C++ / Vispy viewers]
  C03 --> VIEW
  C01 --> RELAX[nanocrystal_pipeline relax]
  RELAX --> C02[02_relaxed.npz]
  C02 --> HESS[hessian topology-linear]
  C03 --> HESS
  HESS --> C04[04_hessian.npz]
  C04 --> SPEC[spectrum eigh]
  SPEC --> C05[05_spectrum.npz]
```

| Layer | Key modules |
|-------|-------------|
| **Generate + export** | `Nanocrystals.js`, `NanocrystalExport.js`, `tests/tSiNCs/nanocrystals.mjs`, `tests/tSiNCs/export_nanocrystal_bundle.mjs` |
| **View** | `MolecularBrowser.cpp` + `cpp/common/io/`; `VispyMolBrowser.py` + plugins |
| **Relax / spectrum** | `nanocrystal_pipeline.py`, `pyBall/io/crystal_npz.py`, `FTIR.py` |
| **Fixtures** | `tests/tSiNCs/fixtures/si_1nm_passivation/` (full pipeline), `fixtures/npz_viewer/` (smoke) |

Viewer guides: [`CPP_MolecularBrowser_NPZ.md`](CPP_MolecularBrowser_NPZ.md) · [`Python_Vispy_MolBrowser_Plugins.md`](Python_Vispy_MolBrowser_Plugins.md)

---

## Suggested reading order

1. [`tests/tSiNCs/README.md`](../../../tests/tSiNCs/README.md)  
2. [`ChemAtlas.plan.md`](ChemAtlas.plan.md) / [`ChemAtlas_MMFF_relax.md`](ChemAtlas_MMFF_relax.md) — Wulff atlas + MMFF pre-relax  
3. [`PySCF_L1_neighborhood_PDOS.md`](PySCF_L1_neighborhood_PDOS.md) — DFT L1 neighborhood PDOS (plot template; vs DFTB+; Si PySCF not yet a minimum)  
4. [`MMFF_C_CH_vs_CH2_kfit.md`](MMFF_C_CH_vs_CH2_kfit.md) — MMFF vs PBE L1 C stretch (negative shared-pack result)  
5. [`Nanocrystal_NPZ_Pipeline.guide.md`](Nanocrystal_NPZ_Pipeline.guide.md)  
6. [`NPZ_Crystal_Schema.md`](NPZ_Crystal_Schema.md)  
7. [`gen_nanocrystals.chat.md`](gen_nanocrystals.chat.md) — generation and defects  
8. [`CPP_MolecularBrowser_NPZ.md`](CPP_MolecularBrowser_NPZ.md) or [`Python_Vispy_MolBrowser_Plugins.md`](Python_Vispy_MolBrowser_Plugins.md) — viewers  
9. [`Nanocrystal_pipeline.progress.md`](Nanocrystal_pipeline.progress.md) — ensemble history  
10. [`../../topical_audit/Nanocrystal_Vibrations.md`](../../topical_audit/Nanocrystal_Vibrations.md) — full inventory  

---

## Guides (`.guide.md`)

| File | Content |
|------|---------|
| [`Nanocrystal_NPZ_Pipeline.guide.md`](Nanocrystal_NPZ_Pipeline.guide.md) | NPZ stage contract, batch tutorial |
| [`NPZ_Crystal_Schema.md`](NPZ_Crystal_Schema.md) | Authoritative NPZ key tables (v1.2) |
| [`CPP_MolecularBrowser_NPZ.md`](CPP_MolecularBrowser_NPZ.md) | C++ SDL browser: NPZ, bond color, AABB |
| [`Python_Vispy_MolBrowser_Plugins.md`](Python_Vispy_MolBrowser_Plugins.md) | Vispy browser + vibration plugin |
| [`Phonon_testing.guide.md`](Phonon_testing.guide.md) | Phonon workflow; PBC/ASR pitfalls |
| [`Sparse_Hessian_Vibration_Spectra.guide.md`](Sparse_Hessian_Vibration_Spectra.guide.md) | Sparse Hessian + spectrum methods |

---

## Reference chats (`.chat.md`)

| File | Content |
|------|---------|
| [`gen_nanocrystals.chat.md`](gen_nanocrystals.chat.md) | Generation CLI, plane cuts, capping, bridges |
| [`Hessian_Kspace.chat.md`](Hessian_Kspace.chat.md) | Bloch / k-space Hessian theory |
| [`Hessian_fitting.chat.md`](Hessian_fitting.chat.md) | Bond/angle stiffness fitting; GPT 5.6 sol 2026-09-03 (restricted Wilson) — checklist in [`Hessian_fitting.md`](../../topical_audit/Hessian_fitting.md) |
| [`XrayDiffraction.chat.md`](XrayDiffraction.chat.md) | Debye vs FFT vs 2D detector; Hessian σ in \(r\). Checklist: [`XRD.md`](../../topical_audit/XRD.md) |

---

## Progress logs (`.progress.md`)

| File | Content |
|------|---------|
| [`Nanocrystal_pipeline.progress.md`](Nanocrystal_pipeline.progress.md) | Ensemble pipeline milestones |
| [`Nanocrystal_generation.progress.md`](Nanocrystal_generation.progress.md) | Generation G0–G5 |
| [`Linearized_topology.progress.md`](Linearized_topology.progress.md) | MMFFL K₁₂/K₁₃/K₁₄ |
| [`Nanocrystal_vibration_sparse.progress.md`](Nanocrystal_vibration_sparse.progress.md) | Sparse vibration staging |
| [`Sparse_vibration_solver.progress.md`](Sparse_vibration_solver.progress.md) | GPU/iterative solvers (de-prioritized) |
| [`XRD_progress.md`](XRD_progress.md) | GPU Debye engine (2026-06). Science: [`XRD.md`](../../topical_audit/XRD.md) |
| [`Debug_ASan_double_free_and_eval_atom.progress.md`](Debug_ASan_double_free_and_eval_atom.progress.md) | MMFF ASan debugging |

---

## Filename convention

| Suffix | Meaning |
|--------|---------|
| `.guide.md` | Curated how-to / tutorial |
| `.chat.md` | Design discussion, CLI reference |
| `.progress.md` | Milestones, session notes |

---

## Cross-links

- [`doc/topical_audit/topical_audit.md`](../../topical_audit/topical_audit.md)  
- [`doc/topical_audit/file_descriptions.md`](../../topical_audit/file_descriptions.md)  
- [`tests/tMMFF/AGENTS.md`](../../../tests/tMMFF/AGENTS.md)  
- [`CODEMAP.md`](../../../CODEMAP.md)  
