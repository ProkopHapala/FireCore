# Hydride neighborhoods — DFTB+ vs MMFF, spectra **and** where they sit

Open **[tests/tSiNCs/OUT_motif_spectra/index.html](../../../tests/tSiNCs/OUT_motif_spectra/index.html)** (figures in `tests/tSiNCs/OUT_motif_spectra/plots/`).

**Plot template for further Si NC figures** is the PySCF L1 stack (same `plot_stacked_method_pdos`, per-row xy inset, 5/7-ring as PDOS groups): [`PySCF_L1_neighborhood_PDOS.md`](PySCF_L1_neighborhood_PDOS.md) · `tests/tSiNCs/OUT_pyscf_jobs/stacked_{C,Si}.png`. L2 figures below keep a 3-view map on **top** because they compare **methods** on one crystal.

Lead figure per crystal: **`l2_compare_*.png`**. Top: 3D atom colors (the PDOS legend). Then **one panel per method** (DFTB+ 3ob-3-1, MMFF default, MMFF CH×1.995): matplotlib `stackplot` of group PDOS. Black line = total DOS. **Σ groups = total DOS** (each mode’s mass-weighted |u|² is partitioned over atoms).

```bash
python3 tests/tSiNCs/plot_hydride_motif_spectra.py --skip-hessian   # L0+L1+L2, reuse Hessians
python3 tests/tSiNCs/plot_hydride_motif_spectra.py --l2-only --skip-hessian
python3 tests/tSiNCs/plot_hydride_motif_spectra.py --l2-only --facet ridge_110  # ⟨110⟩ extrema = edge (default)
python3 tests/tSiNCs/plot_hydride_motif_spectra.py --l2-only --facet miller_111
python3 tests/tSiNCs/plot_hydride_motif_spectra.py --l2-only --facet xh_align
python3 tests/tSiNCs/plot_hydride_motif_spectra.py --l2-only --facet wulff    # revert face/edge to COM support
python3 -m pytest tests/tSiNCs/test_facet_xh_align.py -s
```

---

## L2 C — three crystals that have DFTB+ 3ob-3-1

**L2 MMFF (2026-09-02):** the first bake used the **DFTB-relaxed geometry** (no MMFF minimisation) — those frequencies are **invalid** as harmonic MMFF spectra. See [`Hessian_at_own_minimum.md`](Hessian_at_own_minimum.md). Quarantine: `tests/tSiNCs/OUT_dftb_vs_mmff/WRONG_at_DFTB_geometry/`. L0/L1 MMFF from `relaxed.mol2` was already the correct protocol. L2 MMFF was **recomputed** with FIRE then Hessian (`run_mmff_ownmin_l2.py`); use `OUT_dftb_vs_mmff/<crystal>/mmff_*` plus `mmff_ownmin_protocol.json`. Scaled MMFF at its own min has **0** large imaginaries (the 349 imaginaries were the off-minimum Hessian).

| Crystal | Faces | What the stacked plot shows |
|---------|-------|-----------------------------|
| `truncated_octahedron_C` | {100} CH₂ + {111} CH, ridge edges | DFTB+: **CH₂@{100}** sits above **CH@{111}** (~2986 vs ~2917 cm⁻¹). xh_align: 52 CH@{111}, 72 CH₂@{100} (unchanged), 72 CH edge (Wulff had 84 CH@{110} + 12 CH@{100}). MMFF default: whole packet ~2210, split collapsed. |
| `octahedron_C` | {111} CH faces + ridge edges | DFTB+: one envelope ~2920 (almost pure CH). MMFF default ~2208; scaled ~3100. Few CH₂ (vertices). **xh_align:** 64 CH@{111} vs 120 CH edge (Wulff was 28 vs 156 — isolated CH on the Wulff ribbon are now the face). |
| `rhombic_dodecahedron_C` | Wulff faces are {110}; miller_111 tags CH whose X–H points along a crystal ⟨111⟩ | **`miller_111` (default):** X–H vs 8 axis-aligned ⟨111⟩ (cos≥0.90). This is hydride orientation, not a claim that the rhombic cut has {111} facets. `xh_align` still marks zigzag {110} CH as edge (only ~4 isolated). Jmol: `OUT_motif_spectra/xyz/l2_halogen_*.xyz`. |

**Robustness:** the CH vs CH₂ / {100} vs {111} split is a **DFTB+ 3ob-3-1** feature on these L2 diamonds. Default MMFF at **its own minimum** is still too soft (~2240 vs ~2930). Scaling CH×1.995 (adamantane) then re-relaxing moves the *count* into the stretch window but **does not restore the group substructure**. Scaled fingerprint is softer because angle×0.30, not because of leftover forces.

Kusová CH 2815–2900 / CH₂ 2900–2975 is shaded on the right. DFTB sits slightly high; default MMFF is empty in that window.

---

## Checklist (L1 MMFF atlas, below the L2 block)

| | Question | Spectrum | Spatial |
|---|----------|----------|---------|
| **Q1** | {100} vs {111} | L2 trunc `l2_compare_truncated_octahedron_C.png`; L1 `nbhd_cube_Si.png` vs `nbhd_octa_Si.png` | `map_chem_L2_truncated_octahedron_C.png`; `map_chem_cube_Si.png` / `map_chem_octa_Si.png` |
| **Q2** | XH vs XH₂ | L2: CH₂@{100} vs CH@{111} on trunc; L1 cube: **XH@{110} edge** vs **XH2@{100}** | `loc_cube_Si_SiHwin.png` vs `loc_cube_Si_SiH2win.png` |
| **Q3** | 5-ring / 7-ring | `nbhd_cube5_Si.png`, `nbhd_cube7_Si.png` | purple/orange on `map_chem_cube5_Si.png` / `map_chem_cube7_Si.png` |

**Tagging.** Hydride class from \(n_H\). Four face/edge criteria (switch `--facet`); **do not delete** the older ones:

- **`ridge_110` (default, L2 plots):** C/Si within 0.5 Å of the max or min projection on any of 6 unsigned ⟨110⟩ axes (`heavies_near_110_extrema`; 12 signed ends). Those hydrides are `@edge`. Leftover hydrides take Wulff among primary faces. Unrotated crystals. On a flat {110} rhombic face the whole facet is at the extremum — this trial may paint the face as edge. `--ridge-below` changes the slab.
- **`miller_111`:** sitting {111} = which of the 8 crystal-frame ⟨111⟩ octants the heavy occupies (argmax \(r_{\mathrm{COM}}\cdot\hat n\)). Face if that same \(\hat u_{\mathrm{XH}}\cdot\hat n_{\mathrm{sit}}\ge 0.90\) (`is_xh_on_miller_111`; `--align-cos 0.95` is 5%). Assumes the nanocrystal is **not rotated** off the cubic axes. Matching X–H to *any* of the 8 ⟨111⟩ would tag every diamond hydride (all C–H are tetrahedral ⟨111⟩).
- **`xh_align`:** local hydride chemistry, not COM support. Diamond {111} terrace CH are **not bonded to other hydrides** (only bulk C). Ridge CH are bonded to other CH with tetrahedral X–H (\(|\hat n\cdot\hat n'|\sim 1/3\)). Diamond {100} CH₂: one CH rim allowed. On rhombic {110} the zigzag chains make almost everything “edge”.
- **`wulff` (kept):** winner-take-all COM support `facet_kind_from_vec` among the crystal’s Wulff families.

L1 spatial maps use `ridge_110` on the stored geometry for **all** atlas crystals. L1 stacked PDOS only for Hessians that pass the own-minimum test (today: `cube_Si`, `cube5_Si`). The other L0/L1 `04_hessian.npz` files are **not spectra** (imaginaries or missing rigid modes). On L1-sized particles a 0.5 Å ⟨110⟩ slab covers the whole surface, so almost every hydride is `@edge` — the useful face/edge split is the L2 set. Stretch \(q=\hat r\cdot(u_H-u_X)\). `--facet wulff` still exists.

---

## Method / provenance

- L2 C DFTB+: `/home/prokop/SIMULATIONS/SiNCs/DFTB/{octahedron,truncated_octahedron,rhombic_dodecahedron}_C/` (`modes.npy`, `frequencies_cm1.npy`, `geo_end.xyz`). 3ob-3-1.
- L2 C MMFF: `tests/tSiNCs/OUT_dftb_vs_mmff/<crystal>/mmff_{default,scaled}_hessian.npy` after **MMFF FIRE** (`mmff_ownmin_protocol.json`). Atom order matches atlas mol2 / DFTB. Do not use `WRONG_at_DFTB_geometry/`. Modes from mass-weighted eigh (`cartesian_modes_from_hessian`).
- L0/L1: MMFF FD from `relaxed.mol2`. Unscaled.
- Recompute L2 block: `python3 tests/tSiNCs/plot_hydride_motif_spectra.py --l2-only`

Helpers: `pyBall/FFfit_utils.py` (`neighborhood_xh_groups`, `heavy_cycles`, `apply_ring_tags`, `heavies_near_110_extrema`, `is_xh_on_miller_111`, `is_xh_align_terrace`, `facet_kind_from_vec`, `nbhd_legend_label`, `write_hydride_halogen_xyz`, `cartesian_modes_from_hessian`); stacked plot: `pyBall/FFfit_plots.plot_stacked_method_pdos`; views: `pyBall.plotUtils.plot_nc_views`. Face ≠ edge hues: {111} CH blue, {110} face purple, ridge orange/brown. 5-ring green, 7-ring orange.
