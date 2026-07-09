---
type: TopicalAudit
title: Vibration Data Formats
tags: [vibration, data-format, npz, npy, pyscf, nanocrystal, viewer]
timestamp: 2026-07-09
---

# Vibration Data Formats: Input/Output Contracts for Vibration Pipelines

## Summary

Two distinct data formats coexist in FireCore for storing vibrational analysis
results (frequencies, Hessian, normal modes): the **nanocrystal pipeline NPZ**
format (stages 04/05, produced by `pyBall/nanocrystal_pipeline.py`) and the
**PySCF loose .npy** format (produced by `tests/tSiNCs/run_small_np_pyscf_vib.py`
via `vib_utils.py`). Both are consumed by the `VibrationSpectrumPlugin` in the
VispyMolBrowser and by plotting scripts. This audit documents the exact array
shapes, dtypes, units, and file naming conventions for each format so that
producers and consumers stay in sync.

---

## Implementations

### Format 1: Nanocrystal Pipeline NPZ

| File | Status | Description |
|------|--------|-------------|
| `04_hessian.npz` | **Active** | Mass-weighted Hessian + metadata for re-diagonalization |
| `05_spectrum.npz` | **Active** | Pre-computed histogram + eigenvalues for plotting |

#### `04_hessian.npz` keys

| Key | Shape | Dtype | Units | Notes |
|-----|-------|-------|-------|-------|
| `K` or `K_projected` | (3N, 3N) | float64 | eV/Å² (internal) | `K_projected` has rigid modes shifted; `K` is raw |
| `M` | (3N,) or (N,3) | float64 | amu | Diagonal mass matrix |
| `pos` | (N, 3) | float64 | Å | Relaxed geometry |
| `Z` | (N,) | int32 | — | Atomic numbers |
| `natoms` | () | int | — | Optional; inferred from `pos` if absent |

#### `05_spectrum.npz` keys

| Key | Shape | Dtype | Units | Notes |
|-----|-------|-------|-------|-------|
| `omega_centers` | (nbins,) | float64 | cm⁻¹ | Histogram bin centers |
| `hist` | (nbins,) | float64 | arb. | Mode density histogram |
| `omegas_modes` | (3N,) | float64 | sqrt(eV/amu)/Å (internal) | All eigenvalues; converted to cm⁻¹ via `omega_internal_to_cm1()` |
| `omegas_modes_vib` | (n_vib,) | float64 | internal | Vibrational modes only (rigid filtered) |
| `probe_weight_x/y/z` | (n_vib,) | float64 | arb. | Optional IR/probe projection weights |
| `sigma_bins` | () | float64 | cm⁻¹ | Optional broadening width |
| `units` | () | str | — | Optional unit label |
| `grid_meta` | () | str | — | Optional grid metadata |

**Unit conversion:** internal units → cm⁻¹ via `omega_internal_to_cm1()` in
`pyBall/nanocrystal_pipeline.py` (factor: `omega * C_CM_S / (2*pi)` with
appropriate mass-weighting).

**Eigenvectors:** Not stored in `05_spectrum.npz` v1.2. Re-diagonalized from
`04_hessian.npz` on load by `solve_normal_modes_from_hessian_npz()`.

**Producer:** `pyBall/nanocrystal_pipeline.py` (`hessian` and `spectrum` subcommands).

**Consumer:** `VibrationSpectrumPlugin.set_bundle()` in
`pyBall/GUI/mol_browser_plugins/vibration_spectrum.py`.

---

### Format 2: PySCF Loose .npy

| File | Status | Description |
|------|--------|-------------|
| `frequencies_cm1.npy` | **Active** | Vibrational frequencies in cm⁻¹ |
| `modes.npy` | **Active** | Normal mode displacement vectors |
| `hessian.npy` | **Active** | Full Cartesian Hessian (optional for viewer) |
| `masses.npy` | **Active** | Atomic masses (optional for viewer) |
| `relaxed.xyz` | **Active** | Relaxed geometry (loaded by browser, not plugin) |
| `status.json` | **Optional** | Metadata: `E_hf`, `status` |

#### Array specifications

| File | Shape | Dtype | Units | Notes |
|------|-------|-------|-------|-------|
| `frequencies_cm1.npy` | (3N,) | float64 or complex128 | cm⁻¹ | Complex = imaginary (unstable) modes; `.imag` > 0 |
| `modes.npy` | (n_vib, N, 3) | float64 | Å (displacement) | n_vib = 3N-6; vibrational only, no trans/rot |
| `hessian.npy` | (N, N, 3, 3) | float64 | Hartree/Bohr² | Full Cartesian Hessian; plotted as (3N, 3N) |
| `masses.npy` | (N,) | int64 | amu | Atomic masses |

**Key differences from NPZ format:**
- Frequencies are already in cm⁻¹ (no `omega_internal_to_cm1()` conversion)
- Modes are pre-computed and stored (no re-diagonalization needed)
- Modes contain only vibrational modes (3N-6), not all 3N eigenmodes
- Hessian is in (N,N,3,3) form, not flattened to (3N,3N)
- No histogram pre-computed; built on load from frequency bins

**Hessian plotting layout:** Reshaped as `hess.reshape(N*3, N*3)` which yields
atomic-block ordering `[atom0_xyz, atom1_xyz, ...]` with 3×3 blocks on the
diagonal. This is distinct from the alternative `[all_x, all_y, all_z]`
ordering.

**Producer:** `tests/tSiNCs/run_small_np_pyscf_vib.py` → `vib_utils.py:run_pyscf_vib_full()`.

**Consumer:** `VibrationSpectrumPlugin.set_pyscf_bundle()` in
`pyBall/GUI/mol_browser_plugins/vibration_spectrum.py`; also
`tests/tSiNCs/plot_pyscf_vib_results.py`.

---

## Format Detection Logic

The `VibrationSpectrumPlugin` auto-detects the format in `on_molecule_selected()`:

1. Try `find_nanocrystal_pipeline_stages(mol_dir)` → if `spectrum` or `hessian` key found, use NPZ path (`set_bundle`)
2. Else try `_find_pyscf_vib_dir(mol_dir)` → if `frequencies_cm1.npy` exists, use PySCF path (`set_pyscf_bundle`)
3. Else try `_find_pyscf_vib_dir(ctx.directory)` (parent directory)
4. Else show "missing" message

`is_relevant()` returns True if either NPZ stages or PySCF `.npy` files are found.

---

## Parity Status

| Aspect | NPZ | PySCF .npy | Parity |
|--------|-----|------------|--------|
| Frequency units | internal (needs conversion) | cm⁻¹ (direct) | N/A — different pipelines |
| Mode storage | Re-diagonalized on load | Pre-computed in file | N/A |
| Mode count | 3N (all, rigid filtered by mask) | 3N-6 (vibrational only) | Different indexing |
| Hessian shape | (3N, 3N) flattened | (N, N, 3, 3) block form | Reshape needed |
| Histogram | Pre-computed in 05 | Built on load (10 cm⁻¹ bins) | Different resolution |
| Viewer integration | `set_bundle()` | `set_pyscf_bundle()` | Both produce same internal dicts |

---

## Open Issues

- PySCF `modes.npy` contains only 3N-6 modes; `vib_indices` = [0..n_vib-1] maps
  1:1 to modes array columns. If PySCF ever includes trans/rot modes in the
  file, the mask logic (`omegas_cm > 10.0`) and indexing must be revisited.
- `masses.npy` and `hessian.npy` are not used by the viewer plugin (only by
  plotting scripts). Could be used for re-diagonalization if modes are missing.
- No formal schema validation; missing files cause FileNotFoundError or
  silent "missing" message in viewer.

---

## Related Topics

- [Nanocrystal Vibrations](Nanocrystal_Vibrations.md) — full pipeline audit
- [Topical Audit Index](topical_audit.md)
- [`tests/tSiNCs/README.md`](../../tests/tSiNCs/README.md) — working hub
- [`pyBall/GUI/mol_browser_plugins/vibration_spectrum.py`](../../pyBall/GUI/mol_browser_plugins/vibration_spectrum.py) — viewer plugin
- [`pyBall/nanocrystal_pipeline.py`](../../pyBall/nanocrystal_pipeline.py) — NPZ producer
- [`tests/tSiNCs/plot_pyscf_vib_results.py`](../../tests/tSiNCs/plot_pyscf_vib_results.py) — PySCF plotting
- [`tests/tSiNCs/vib_utils.py`](../../tests/tSiNCs/vib_utils.py) — PySCF producer
