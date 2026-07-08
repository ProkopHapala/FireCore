# Nanocrystal / NPZ Pipeline — Open Items

Tracked follow-ups for `tests/tSiNCs/` and the Si/diamond vibration workflow.  
**Done work** is documented in guides and the topical audit — not here.

**Hub:** [`README.md`](README.md) · **Pipeline:** [`Nanocrystal_NPZ_Pipeline.guide.md`](../../doc/Topics/FTIR_Nanocrystals/Nanocrystal_NPZ_Pipeline.guide.md) · **Schema:** [`NPZ_Crystal_Schema.md`](../../doc/Topics/FTIR_Nanocrystals/NPZ_Crystal_Schema.md) · **Audit:** [`Nanocrystal_Vibrations.md`](../../doc/topical_audit/Nanocrystal_Vibrations.md)

---

## Viewers & NPZ

| Item | Priority | Notes |
|------|----------|-------|
| Cross-viewer parity test | Medium | Automated: same NPZ → C++ `--verify-npz` vs Python `load_crystal_npz`; compare `natoms`, `n_groups` |
| C++ bonds from `neigh_idx` | Medium | Port `bonds_from_neigh_idx` from `pyBall/io/crystal_npz.py`; `03_topology.npz` today uses distance fallback |
| C++ `bond_k` stiffness coloring | Low | Python Vispy has it; C++ has length-only map |
| C++ lazy thumbnail NPZ parse | Low | All files parsed at grid load |
| C++ vibration / spectrum panel | Low | Python `VibrationSpectrumPlugin` leads |
| Browser NPZ file picker (web) | Low | `nanocrystalViewer.html` still JSON/manifest only |
| Commit / gitignore fixtures | High | `fixtures/` largely untracked; add `.gitignore` for `.vispy_mol_browser_cache*`, `out_defect_export/` |

---

## Pipeline & force fields

| Item | Priority | Notes |
|------|----------|-------|
| Hessian MMFF FD vs topology-linear parity | Medium | `04_hessian.status.json` has `parity_max_rel_diff: NaN`; wire `--compare-mmff` reporting |
| LFF relax on `03_topology.npz` | Medium | Relax is C++ MMFF; Hessian uses exported MMFFL sticks — document or implement XPDB/LFF relax |
| `tests/tSiNCs/nanocrystals.mjs` | Medium | **Done** — unified CLI in `tests/tSiNCs/`; legacy `scripts/` copies deprecated |
| Store eigenvectors in `05_spectrum.npz` | Low | v1.2 omits modes; plugin re-eighs from `04` |
| Ensemble at scale (N=100+) | Low | Timing report; confirm Hessian vs relax bottleneck |
| `relax` without temp mol2 | Low | `--init-npz` works; could skip mol2 round-trip entirely in MMFF init |

---

## Validation & references

| Item | Priority | Notes |
|------|----------|-------|
| `bootstrap_vibration_parallel_fixtures.py` | Medium | Referenced by `tests/tMMFF` solver ladder; script missing |
| QM ↔ MMFF adamantane parity | Medium | Compare `run_vib_spectra.py` vs `test_vibration_spectra.py` with tolerances |
| Generator crosscheck beyond sphere | Low | Miller planes / bridges — extend `crosscheck_nanocrystal_generators.py` or document JS-only |
| DFTB SK paths in `vib_utils.py` | Low | Hardcoded `/home/prokop/SIMULATIONS/...`; use `DFTB_SK_ROOT` env |
| `plot_phonon_bands.py` | Low | Extract from `test_diamond_phonon_bands.py` |
| Nanocrystal-scale QM references | Low | QM track is molecular cages only |

---

## Physics & scope (open questions)

| Question | Context |
|----------|---------|
| Relax surface H before vibration? | H caps are static at generation; may affect surface modes |
| When to use MMFF FD vs topology Hessian? | Production default is topology-linear; FD for validation only? |
| Defect typing after Si bridge insert | Z reassignment only; full MMFF atom types not refreshed |
| PBC phonons on nanocrystals vs clusters | Cluster preferred for Γ optical; document per use case |
| SiO₂ / silica | Out of scope unless new generator |

---

## Resolved (2026-07-07)

- NPZ crystal + topology export (`export_nanocrystal_bundle.mjs`, `NanocrystalExport.js`)
- C++ / Python NPZ molecular browsers
- `--init-npz` relax, topology-linear Hessian in pipeline
- Si passivation gallery `fixtures/si_1nm_passivation/` (01→05)
- Schema v1.2, pipeline guide, viewer user guides
