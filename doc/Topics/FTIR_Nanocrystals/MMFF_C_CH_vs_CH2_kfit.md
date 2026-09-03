---
type: report
title: MMFF C L1 — one k_XH cannot fit CH and CH₂
description: Own-min stretch RMSE vs PySCF PBE; surface Morse/LJ; split XXX/XXH/HXH; SA did not unify cube vs octahedron
tags: [nanocrystal, mmff, ftir, k-fit, morse, ch2]
timestamp: 2026-09-02
---

# MMFF C L1 — CH vs CH₂, steric NB, split angles (2026-09-02)

**Status:** negative but useful. A single transferable MMFF pack does **not** put cube {100} CH₂ and octahedron {111} CH into the same stretch window. Functionality to test that claim **exists** and should be reused, not re-invented.

**Do not cite `cube_C` or `octahedron_7ring_Si` PBE frequencies as FTIR** (still 1 imaginary each; mode-follow pending). Five other Si L1 jobs in `pyscf_vib_results` are now 0-imag; the 2026-09-02 Si own-min RMSE in `OUT_mmff_kss_fit/pyscf_pbe_ncs/*_Si` used the **old** off-min Hessians — discard those Si numbers. C is the fitted reference except `cube_C`.

SSOT for the June angle×0.30 bug and XH₄ \(k\): [`MMFF_stiffness_scaling.md`](MMFF_stiffness_scaling.md). Spectrum protocol: [`Hessian_at_own_minimum.md`](Hessian_at_own_minimum.md). Atlas surface-NB policy: [`ChemAtlas_MMFF_relax.md`](ChemAtlas_MMFF_relax.md).

---

## What was missing in the first general-C fit

`--mol pyscf_nc --ref joint` used **bonded+angles only** (`NonBonded=-1`). Cluster Hessian historically **forced NB off** (PBC 27-image phonon bug). Atlas `--nonbond surface` (H + undercoordinated heavies, Exclusion2 skip 1–2/1–3) was never on in the fitter.

On cubes that is unphysical: bonded-only FIRE on `cube_5ring_C` at \(k_{\mathrm{XH}}=1.775\) put vicinal H···H at **0.67 Å**. Octahedron has no such contact (PBE min vicinal **2.52 Å**).

Default kernel: **one `Kss` per heavy atom** (`apars[ia].z = AtomTypes.Kss×4 = 25.28`). H–C–H, C–C–H, and C–C–C on the same carbon share it. `AngleTypes.dat` is `* C * 109.5 5.1` — same \(k\) for all C-centred triples. `bEachAngle` + `angles[nnode×6]` can split slots; we copy from `apars` so scale=1 is **not** a silent jump to AngleTypes 5.1.

Cluster pairwise NB is **LJ** (`getLJQH`). Morse α lives in GridFF; `bMorseNonBond` now switches Exclusion2 to `getMorseQH`. UFF H well depth ~10⁻⁶ eV, so α must be large (~3) before it is a Pauli wall. Scaling H EvdW is the more direct steric knob.

---

## How to run (from `tests/tMMFF`)

```bash
python3 ../tSiNCs/fit_mmff_kss_pyscf.py --mol pyscf_nc --ref joint      # bonded-only minimax k_XH
python3 ../tSiNCs/fit_mmff_kss_pyscf.py --mol pyscf_nc --ref nbscan     # H EvdW scale at frozen k
python3 ../tSiNCs/fit_mmff_kss_pyscf.py --mol pyscf_nc --ref morse2d    # k_XH × Morse α
python3 ../tSiNCs/fit_mmff_kss_pyscf.py --mol pyscf_nc --ref anneal     # SA: k_XH, α, K_XXX, K_XXH, K_HXH
python3 ../tSiNCs/fit_mmff_kss_pyscf.py --mol pyscf_nc --ref anneal_plot
```

One MMFF molecule per process. Rebuild `MMFF_lib` after C++ Hessian/NB changes (`cmake --build cpp/Build --target MMFF_lib`). Session: `CompChemUtils/examples/tSiNCs/mmff_molecular_session.py` (`enable_surface_nonbond`, `set_morse_nonbond`, `enable_each_angle`).

Outputs (gitignored): `tests/tSiNCs/OUT_mmff_kss_fit/pyscf_pbe_ncs/`.

---

## Results (C, own-min stretch RMSE vs PBE/ccECP-cc-pVDZ)

PBE itself already splits the hydride window: octahedron_C ⟨stretch⟩ **2861** cm⁻¹ (24 CH + 6 CH₂); cube_5ring_C **3011** cm⁻¹ (12 CH + 11 CH₂, 17 vicinal H···H < 2.2 Å).

### Bonded-only general FF (`--ref joint`)

| Pack | \(k_{\mathrm{XH}}\) | \(K_{\mathrm{ss}}\) | worst RMSE | mean | Note |
|------|--------------------:|--------------------:|-----------:|-----:|------|
| C all 6 L1 | 1.775 | 0.488 | 98 (`cube_C`) | 83 | bonded-only start |
| Si all 6 | 1.097 | 0.400 | 71 | 55 | **do not use as FTIR** (off-minimum) |

`joint1d_C.png`: octa V-min ~1.68–1.72; cube V-min ~1.85–1.90. Shared minimax ~100 cm⁻¹. Same MMFF mean stretch (~2938) on both families when NB is off — the DFT split is chemistry/geometry, not two \(k_{\mathrm{XH}}\) already in MMFF.

### H EvdW scale at frozen \(k_{\mathrm{XH}}=1.775\) (`nbscan_C.png`)

Surface LJ, Exclusion2. Cube improves (87 → **46** at ×4); octahedron **worsens** (83 → 176). Families move **apart**. Cube bonded-only clash opens when NB is on.

### 2D \(k_{\mathrm{XH}}\) × Morse α (`morse2d_C.png`)

| Crystal | best \(k_{\mathrm{XH}}\) | best α | RMSE | MMFF ⟨st⟩ vs PBE |
|---------|-------------------------:|-------:|-----:|------------------|
| octahedron_C | 1.683 | 0.80 | **30** | 2862 vs 2861 |
| cube_5ring_C | 1.817 | 3.00 | **35** | 3011 vs 3011 |

Each family can be fitted to ~30 cm⁻¹. They want **opposite corners**. Shared minimax on that 4×4 grid: **114** cm⁻¹ (\(k=1.683\), α=3). Soft Morse on the cube still clashes (H···H 0.83 Å).

### SA (13 evals, two crystals)

Knobs: \(k_{\mathrm{XH}}\), \(K_{\mathrm{Morse}}=\alpha\), \(K_{\mathrm{XXX}}\), \(K_{\mathrm{XXH}}\), \(K_{\mathrm{HXH}}\). Loss = mean stretch RMSE. Start = 2D minimax + equal \(K_{\mathrm{ss}}=0.488\).

| trial | mean | max | octa | cube | \(k_{\mathrm{XH}}\) | α | XXX | XXH | HXH |
|------:|-----:|----:|-----:|-----:|--------------------:|--:|----:|----:|----:|
| 0 | 75.0 | 114 | 36 | 114 | 1.683 | 3.00 | 0.488 | 0.488 | 0.488 |
| **3** | **71.8** | **89** | 55 | 89 | 1.715 | 3.14 | 0.353 | 0.357 | 0.541 |
| 9 | 80.2 | 125 | 125 | **35** | 1.809 | 3.12 | 0.200 | 0.461 | 0.548 |

Incumbent **froze at trial 3**. Δmean only **−3.2** cm⁻¹. Running-best last four trials: all 71.85. Later grey trials still explore (\(k_{\mathrm{XH}}\) to 1.85, XXX to the 0.20 bound) — not a trapped local well, a **see-saw**. Trial 9 matches the cube’s private 2D fit and destroys the octahedron.

Plots: `anneal_rmse.png`, `anneal_params.png`, `<case>/anneal/spectra_overlay.png`. Octahedron stretch **blue-shifts** (RMSE 36→55). Cube **red-shift shrinks** (114→89) but remains ~80–100 cm⁻¹ below PBE’s 3000–3150 cluster.

---

## Code that had to exist for this test

| Piece | Role |
|-------|------|
| `getHessian3Nx3N` | Cluster: keep `bNonBonded`; `eval_no_omp` so Exclusion2 matches FIRE. PBC still bonded-only. |
| `NBFF::bMorseNonBond` / `setMorseNonBond` | Exclusion2 Morse α vs LJ |
| `setEachAngle` + `angles` buffer | Per-slot \(K_{\mathrm{ss}}\) (XXX / XXH / HXH) |
| `fit_mmff_kss_pyscf.py` | `--ref joint\|nbscan\|morse2d\|anneal\|anneal_plot` |
| `mmff_molecular_session.py` | surface mask, Morse, each-angle copy from `apars` |

`nanocrystal_pipeline.relax --nonbond surface` is the atlas path (LJ, not Morse). `_mmff_hessian_fd` there still forces `NonBonded=-1` — atlas production spectra are still bonded-only unless that is changed on purpose.

---

## Conclusion

{100} CH₂ crowding is real (geometry + PBE frequencies). Surface H···H is required on cubes. Splitting HXH from CCC is the right extra \(K\). **None of that yields one pack** that is ~30 cm⁻¹ on both families. Next physics is typed hydrides (CH vs CH₂ \(k\)), not a longer SA on the same five scalars.

Do **not** FIRE with June 2026 `apars[:,1]×0.30`.
