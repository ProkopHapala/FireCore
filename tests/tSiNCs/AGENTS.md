# Si / Diamond Nanocrystal Vibrations & FTIR

## Purpose

Own the **end-to-end nanocrystal vibration workflow** in FireCore: structure generation (JS/Python), MMFF relax/Hessian/spectrum (Python + C++), QM reference spectra, NPZ fixtures, and C++/Vispy viewers. Canonical home for all nanocrystal **orchestration scripts** and test outputs (`OUT_*`).

## Ownership

- Node CLIs: `nanocrystals.mjs`, ensemble/atlas configs, small symmetric batch generator
- Python test drivers: `run_vib_spectra.py`, `crosscheck_nanocrystal_generators.py`, NPZ bootstrap/export helpers
- FFfit scripts: `test_FFfit.py` (thin CLI wrapper), `test_fffit_hybrid.py` (hybrid fitting tests), `test_parity_py_cpp.py`, `test_parity_graph_cpp.py` — reusable logic lives in `pyBall/FFfit_utils.py` and `pyBall/FFfit_plots.py`
- Own-min \(k\)-fit vs PySCF: `fit_mmff_kss_pyscf.py` (run from `tests/tMMFF`). C cube vs octa did **not** share one pack: [`MMFF_C_CH_vs_CH2_kfit.md`](../../doc/Topics/FTIR_Nanocrystals/MMFF_C_CH_vs_CH2_kfit.md).
- Committed fixtures: `fixtures/si_1nm_passivation/`, `fixtures/npz_viewer/`, `fixtures/vibration_*`
- Viewer launchers: `run_cpp_mol_browser.sh`, `run_vispy_mol_browser.sh`, `test_cpp_npz_load.sh`
- Generated output trees under this folder (`OUT_chem_atlas/`, `OUT_small_nc/`, `OUT_nc_ensemble_v2/`, `OUT_nc_atlas/`, `OUT_mmff_kss_fit/`, `OUT_pyscf_jobs/`, `OUT_XRD/`) — not committed unless promoted to fixtures

## Local Contracts

- **Harmonic spectrum = Hessian at that method’s own minimum.** Relax with the **same** potential (switches, CH/angle scales, nonbond) until \(f_{\max}<f_{\mathrm{conv}}\), then Hessian. DFTB geometry + MMFF Hessian is **FFfit only**, never a spectrum. Smoking gun: MMFF modes piled from ~18 cm⁻¹ or hundreds of large imaginaries. SSOT: [`doc/Topics/FTIR_Nanocrystals/Hessian_at_own_minimum.md`](../../doc/Topics/FTIR_Nanocrystals/Hessian_at_own_minimum.md). Invalid 2026-09-02 L2 MMFF files: `OUT_dftb_vs_mmff/WRONG_at_DFTB_geometry/`. Canonical PySCF L1: `/home/prokop/SIMULATIONS/SiNCs/pyscf_vib_results` (10/12 0-imag). Still not FTIR: `cube_C`, `octahedron_7ring_Si`. Legacy `pySCF/jobs/results` Si are superseded. DFTB+ L1 (`/home/prokop/SIMULATIONS/SiNCs/DFTB/L1`, all 24 jobs) **is** at a minimum (`n_imag=0`); comparison: `OUT_dftb_vs_pyscf_l1/`.
- **Stacked neighborhood PDOS is the Si-NC plot template** (2026-09-02): chemistry colors (not one color per crystal), black total DOS, eigenstate rug, xy inset in the empty \(\omega\)-gap with bonds, 5/7-ring as exclusive PDOS groups. API: `pyBall.FFfit_plots.plot_stacked_method_pdos`. Do not invent a second stacked plotter. SSOT: [`doc/Topics/FTIR_Nanocrystals/PySCF_L1_neighborhood_PDOS.md`](../../doc/Topics/FTIR_Nanocrystals/PySCF_L1_neighborhood_PDOS.md).
- **Scripts resolve `REPO` as FireCore root** (`path.resolve(__dirname, '../..')`); **outputs default under `TEST_DIR`** (this folder).
- **Libraries stay outside:** `web/molgui_webgpu/Nanocrystals.js`, `web/common_js/npzIO.js`, `pyBall/nanocrystal_pipeline.py`, `pyBall/FTIR.py`, `pyBall/FFfit_utils.py`, `pyBall/FFfit_plots.py`.
- **Production spectrum:** dense `eigh` + mode histogram via `pyBall/nanocrystal_pipeline.py` (full MMFF Hessian FD at relax geometry). Linearized `03_topology.npz` is built in ensemble runs but **not** used for production spectra yet (Path B deferred).
- **Legacy `scripts/` copies** of nanocrystal CLIs are deprecated; edit and run **`tests/tSiNCs/`** versions only.
- **CIF and MMFF `.dat` files** come from `cpp/common_resources/` (read-only from REPO).

## Work Guidance

### Entrypoints (from repo root)

```bash
node tests/tSiNCs/nanocrystals.mjs ensemble --atlas tests/tSiNCs/chem_atlas.json \
  --output-dir tests/tSiNCs/OUT_chem_atlas
python3 tests/tSiNCs/chem_atlas_report.py
./tests/tSiNCs/run_cpp_mol_browser.sh tests/tSiNCs/OUT_chem_atlas/atlas/L1_dft
```

Results: `tests/tSiNCs/OUT_chem_atlas/atlas/` (gitignored). Policy + file map: [`doc/Topics/FTIR_Nanocrystals/ChemAtlas_MMFF_relax.md`](../../doc/Topics/FTIR_Nanocrystals/ChemAtlas_MMFF_relax.md).

L2 MMFF spectrum (relax with that FF, then Hessian; never at DFTB *q*). Default and 3ob-fitted \(k\) (`CH×1.812`, `Kss×0.488`, `CC×1`):

```bash
cd tests/tMMFF
python3 ../tSiNCs/run_mmff_ownmin_l2.py rhombic_dodecahedron_C
```

CH₄ / C₂H₆ / SiH₄ / Si₂H₆ / adamantane / Si₁₀H₁₆ stiffness re-fit (`apars[:,2]` = \(K_{\mathrm{ss}}\), frozen QM geom):

```bash
cd tests/tMMFF
python3 ../tSiNCs/fit_mmff_kss_pyscf.py --mol all_si
python3 ../tSiNCs/fit_mmff_kss_pyscf.py --mol adamantane
python3 ../tSiNCs/fit_mmff_kss_pyscf.py --mol CH4 --ref dftb_3ob-3-1
python3 ../tSiNCs/fit_mmff_kss_pyscf.py --mol pyscf_nc --ref joint   # one C pack + one Si pack; FIRE+Hessian; whole-spectrum overlays
python3 ../tSiNCs/fit_mmff_kss_pyscf.py --mol pyscf_nc --ref nbscan  # H EvdW at frozen k (cube vs octa)
python3 ../tSiNCs/fit_mmff_kss_pyscf.py --mol pyscf_nc --ref morse2d
python3 ../tSiNCs/fit_mmff_kss_pyscf.py --mol pyscf_nc --ref anneal
python3 ../tSiNCs/fit_mmff_kss_pyscf.py --mol pyscf_nc --ref anneal_plot
```

See [`MMFF_stiffness_scaling.md`](../../doc/Topics/FTIR_Nanocrystals/MMFF_stiffness_scaling.md) (k tables) and [`MMFF_C_CH_vs_CH2_kfit.md`](../../doc/Topics/FTIR_Nanocrystals/MMFF_C_CH_vs_CH2_kfit.md) (CH vs CH₂ + Morse; **not** a recommended general FF).

Hessian pin-ladder (frozen QM \(q\), local \(K_{ab}\), hydride freeze; **not** FTIR):

```bash
CPP_BUILD_PATH=$PWD/cpp/Build-opt/libs python3 tests/tSiNCs/test_FFfit.py --ladder si
```

PySCF vs DFTB+ L1 (same 12 crystals; DFTB from `vibrations.tag`):

```bash
python3 tests/tSiNCs/plot_hydride_motif_spectra.py \
  --pyscf /home/prokop/SIMULATIONS/SiNCs/pyscf_vib_results \
  --dftb-l1 /home/prokop/SIMULATIONS/SiNCs/DFTB/L1
```

Output: `tests/tSiNCs/OUT_dftb_vs_pyscf_l1/` (`compare_<crystal>.png` + DOS overlay).

PySCF L1 neighborhood PDOS (read npy from MetaCentrum; does not re-run DFT):

```bash
python3 tests/tSiNCs/plot_hydride_motif_spectra.py \
  --pyscf /home/prokop/SIMULATIONS/SiNCs/pyscf_vib_results \
  --pyscf-out tests/tSiNCs/OUT_pyscf_jobs
```

### Subcommands (`nanocrystals.mjs`)

| Subcommand | Role |
|------------|------|
| `generate` | Single crystal mol2/xyz |
| `ensemble` | Batch gen → relax → topology → Hessian → spectrum → accumulate |
| `topology` | MMFFL linearized `03_topology.npz` |
| `audit` | Geometry presets G0–G5 |
| `nonbond` | Collision / nonbond group debug |
| `rings` | Ring detection SVG |

### Documentation triad (keep in sync)

| Doc | Role |
|-----|------|
| [`README.md`](README.md) | Folder index + run commands |
| [`doc/Topics/FTIR_Nanocrystals/README.md`](../../doc/Topics/FTIR_Nanocrystals/README.md) | Topic guides and progress logs |
| [`doc/Topics/FTIR_Nanocrystals/ChemAtlas_MMFF_relax.md`](../../doc/Topics/FTIR_Nanocrystals/ChemAtlas_MMFF_relax.md) | Atlas output paths + surface-NB policy |
| [`doc/Topics/FTIR_Nanocrystals/PySCF_L1_neighborhood_PDOS.md`](../../doc/Topics/FTIR_Nanocrystals/PySCF_L1_neighborhood_PDOS.md) | PySCF L1 PDOS plot template; imaginaries vs experiment |
| [`doc/topical_audit/SiNCs.md`](../../doc/topical_audit/SiNCs.md) | Project briefing (FTIR + XRD + FFfit) |
| [`doc/topical_audit/XRD.md`](../../doc/topical_audit/XRD.md) | Powder XRD / pair-distribution checklist |
| [`doc/topical_audit/Hessian_fitting.md`](../../doc/topical_audit/Hessian_fitting.md) | Hessian / Wilson FFfit checklist |
| [`doc/topical_audit/Nanocrystal_Vibrations.md`](../../doc/topical_audit/Nanocrystal_Vibrations.md) | Vibration API inventory |

7. **Expand adaptively.** If experimental residuals contain a feature that no current motif can reproduce, generate structures specifically containing plausible missing motifs. Do not blindly increase ensemble size.

## Verification

- Generator parity: `python3 crosscheck_nanocrystal_generators.py` (from this folder)
- Small NC batch: `node make_small_symmetric_nc.mjs` → `OUT_small_nc/index.md`
- Chem atlas: `node nanocrystals.mjs ensemble --atlas chem_atlas.json` then `python3 chem_atlas_report.py` (L0/L1 need `relaxed.xyz`)
- NPZ gallery: `fixtures/si_1nm_passivation/` stages 01→05 per crystal
- Classical ladder: `tests/tMMFF/run.sh test_vibration_spectra.py`
- FFfit hybrid tests: `python3 -m pytest tests/tSiNCs/test_fffit_hybrid.py -v` (5 test cases)
- Graph + batch dihedral parity: `PYTHONPATH=.. python3 test_parity_graph_cpp.py` (14 tests)

## Child DOX Index

| Path | Role |
|------|------|
| [`fixtures/si_1nm_passivation/README.md`](fixtures/si_1nm_passivation/README.md) | Nine-crystal NPZ pipeline gallery |
| [`fixtures/npz_viewer/README.md`](fixtures/npz_viewer/README.md) | Minimal viewer smoke fixtures |
| [`fixtures/vibration_benchmarks/README.md`](fixtures/vibration_benchmarks/README.md) | Benchmark NPZ structures |
