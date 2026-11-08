# NEVER compute a harmonic spectrum off that method’s own minimum

**Status:** hard protocol. Violating it produces **scientifically invalid frequencies**.  
**Caught:** 2026-09-02, L2 diamond MMFF vs DFTB+ (`tests/tSiNCs/OUT_dftb_vs_mmff/`).  
**Quarantine:** `tests/tSiNCs/OUT_dftb_vs_mmff/WRONG_at_DFTB_geometry/` (do not cite those MMFF curves as spectra).

---

## The rule

For **each** energy model (MMFF, scaled MMFF, DFTB+, DFT, …):

1. Relax with **that** potential (same switches, same scales, same nonbond policy).
2. Confirm a stationary point: \(f_{\max} < f_{\mathrm{conv}}\).
3. Only then evaluate the Hessian / frequencies **at that geometry**.

A valid harmonic spectrum has **~6 rigid-body modes** near zero, then a **gap**, then the first vibration. DFTB+ 3ob-3-1 on these L2 diamonds: 6 modes \(|\nu|\le 10\,\mathrm{cm}^{-1}\), then first vibration \(\sim 135\)–\(162\,\mathrm{cm}^{-1}\), **zero** modes in \(10\)–\(150\,\mathrm{cm}^{-1}\).

If MMFF “starts at 0” on a 0–3300 plot, or scaled MMFF has **hundreds of large imaginary modes**, you are **not** at an MMFF minimum. The plot is not a spectrum. Stop. Relax. Recompute.

**Same switches for relax and Hessian.** Relaxing with surface nonbond then taking a bonded-only Hessian is the same class of error (residual forces of a different energy). Isolated-cluster `getHessian3Nx3N` now **keeps** Exclusion2 NB when FIRE had it (`eval_no_omp`); PBC Hessians stay bonded-only (27-image NB is not a Hessian). Atlas `_mmff_hessian_fd` still forces `NonBonded=-1` — change both or neither. CH vs CH₂ steric tests: [`MMFF_C_CH_vs_CH2_kfit.md`](MMFF_C_CH_vs_CH2_kfit.md).

---

## Two different questions — do not mix them

| Question | Geometry | Valid product | Invalid product |
|----------|----------|---------------|-----------------|
| **Harmonic spectrum / DOS / FTIR peak positions** | Each method’s **own** minimum | Frequencies, DOS, PDOS | — |
| **FFfit / force-constant match** | **One shared** geometry (often QM/DFTB-relaxed) | Hessian residual \(\|H_{\mathrm{MMFF}}-H_{\mathrm{ref}}\|\) at that \(q\) | Calling those eigenvalues a “spectrum” |

[`Matching_Exp_Atlas.md`](Matching_Exp_Atlas.md) already asked for **both** comparisons. Shared-geometry Hessians isolate \(k\)-error. Own-minimum spectra are the predictive workflow. Plotting MMFF eigenvalues at the DFTB geometry as if they were MMFF vibrations **collapses the two** and is wrong.

L0/L1 chem-atlas MMFF Hessians from `relaxed.mol2` follow the rule. The 2026-09-02 L2 overlay did not.

---

## What went wrong (L2 C, 2026-09-02)

MMFF Hessians for `rhombic_dodecahedron_C`, `truncated_octahedron_C`, `octahedron_C` were evaluated at **DFTB-relaxed** `relaxed.xyz`, with **no MMFF FIRE**. Switches: bonds+angles, `NonBonded=-1`.

Smoking gun, rhombic, all 1503 eigenvalues, no rigid-body projection (plot cut \(\nu>10\)):

| Method | Geometry | Rigid \(\|\nu\|\le 10\) | Modes in 10–150 | First vibration | Large imaginaries |
|--------|----------|------------------------:|----------------:|----------------:|------------------:|
| DFTB+ 3ob-3-1 | DFTB min | 6 | **0** | 161.8 cm⁻¹ (×3) | ~5 numerical \(\lesssim 3\) |
| MMFF default | **DFTB geom** | **3** (translations only) | **14** | **18 cm⁻¹** | — |
| MMFF CH×1.995, angle×0.30 | **DFTB geom** | — | — | — | **349** modes to **−558 cm⁻¹** |

Residual forces/torques at a non-stationary point of MMFF. Scaling CH/angles on a DFTB geometry **without re-relaxing** makes the Hessian indefinite. That is not “missing trans/rot removal”.

HTML [`DFTB_vs_MMFF_L2.html`](DFTB_vs_MMFF_L2.html) and L2 stacked PDOS in [`HydrideMotif_spectra.md`](HydrideMotif_spectra.md) repeated the same MMFF files. **Do not use those MMFF frequencies, CH RMSE, or “scaled MMFF fingerprint” claims.** DFTB-only overlay (`OUT_dftb_vs_mmff/dftb_only/`) is still valid.

---

## Never do this again

- Do **not** pass DFTB/DFT `relaxed.xyz` to `MMFFMolecularSession.get_hessian` / `compute_frequencies` and call the result a spectrum.
- Do **not** apply `set_scales` / `set_scales_per_bond_type` and Hessian **without** relaxing **after** the scale.
- Do **not** drop the first 6 eigenvalues blindly, or `sqrt(max(λ,0))`, to hide a pile of 18 cm⁻¹ modes or imaginaries.
- **Fail loud** if \(n_{\mathrm{imag}}(|\nu|>10)\ne 0\), or if fewer than 6 modes sit in \(|\nu|<10\) (missing rotations). Do **not** `sqrt(max(λ,0))`, drop imaginaries, or `if n_imag > 80`. Shared helper: `pyBall.FFfit_utils.assert_harmonic_spectrum_at_minimum` (called from `cartesian_modes_from_hessian`, `solve_normal_modes_from_hessian_npz`, L2 plot/runner, `vibration_spectrum_from_modes`).
- `MMFF` cannot re-`init` a second molecule in one process (`Buffers already initialized`). One process per crystal is required; it is **not** an excuse to skip relax.

Correct MMFF sequence (one process):

```text
init from mol2 (explicit bonds) → optional set_scales → FIRE relax (same switches as Hessian)
→ check fmax → FD Hessian → signed frequencies (imaginary = negative)
```

Atlas L2 start: `tests/tSiNCs/OUT_chem_atlas/atlas/L2_dftb/<shape>_C.mol2` (same lattice+H as DFTB `start.xyz`). Each method relaxes from that start with **its** energy.

Runner: `tests/tSiNCs/run_mmff_ownmin_l2.py` (one crystal per process).

---

## “Scaled MMFF” is stiffness `k`, not a size scale — and we hit the **wrong angle column**

**You do not scale geometry.** `l0` (`bLs`) was never multiplied. Bond stiffness `bKs` ×1.995 is a real `k` scale (CH4 stretch fit).

**Angle ×0.30 was not `k`.** `MMFFMolecularSession._angle_col = 1` multiplies `apars[:,1]`. In the live kernel (`MMFFsp3_loc.h`):

```text
cs0_ss = (apar.x, apar.y)   # equilibrium (cos θ0, sin θ0)  — or half-angle cs
ssK    = apar.z             # angle STIFFNESS Kss
```

Column **1 is sin(θ₀)**, the equilibrium angle, not K. Column **2 is K**. Python comments (`apars = [c0, Kss, Ksp, c0_e]`) and the adamantane report (`apars[:,1]`) copied a **stale layout**. Builder writes `apars = (cos θ0, sin θ0, K, …)`.

Harmonic bond/angle energy: minimum is at `r=l0`, `θ=θ0`, **independent of `k`**. Scaling only `k` cannot collapse a crystal. Scaling `sin(θ0)` by 0.30 **moves the angle equilibrium** (and denormalizes the (c,s) pair). FIRE then follows that broken `θ0`. That is why Rg went 6.8→1.2 Å and C–C → 0.2 Å. Not “stiffness made diamond implode”.

Source of 1.995 / 0.30: `CompChemUtils/examples/tSiNCs/MMFF_VIBRATION_FITTING_REPORT.md` (CH4 / C2H6 / adamantane), Hessian at **frozen** QM geometry. Those jobs never FIRE’d, so the wrong angle column only faked softer bends in the Hessian. Transferring that “scale” into FIRE on an L2 NC is what collapsed the particle.

**Do not FIRE with `apars[:,1]×0.30`.** Do not treat the collapsed orange curves as a spectrum. Default-MMFF own-minimum (blue) is still a diamond. Full trail + correct column: [`MMFF_stiffness_scaling.md`](MMFF_stiffness_scaling.md). `0.30` is **not** a valid \(K_{\mathrm{ss}}\) until CH₄ is re-fit.

---

## Recalculation (2026-09-02, after this note)

FIRE with the Hessian switches (`NonBonded=-1`, angles on), `Fconv=1e-4`, then FD Hessian. Scaled: CH×1.995 and angle×0.30 **then** FIRE.

| Crystal | MMFF default ν₁ / n_rigid / n_imag | MMFF scaled ν₁ / n_imag | RMSD vs DFTB default / scaled (Å) |
|---------|--------------------------------------|-------------------------|-----------------------------------|
| rhombic | 121.9 / 6 / 0 | 37.7 / 0 | 0.064 / 3.56 |
| trunc octa | 117.4 / 6 / 0 | 30.0 / 0 | 0.145 / 3.87 |
| octa | 99.5 / 6 / 0 | 39.4 / 0 | 0.052 / 3.53 |

Default MMFF C–H remains ~2240 cm⁻¹ (FF too soft). Scaled stretch count = nH, mean ~3105 vs DFTB ~2930–2945. Scaled fingerprint is genuinely softer (angle×0.30), **not** 349 imaginaries. Do not cite the quarantine folder as spectra.

---

## Where this is wired

- This file (topic SSOT for the pitfall).
- [`tests/tSiNCs/AGENTS.md`](../../../tests/tSiNCs/AGENTS.md) local contract.
- `pyBall/FFfit_utils.py` — `signed_frequencies_cm1`, `assert_harmonic_spectrum_at_minimum` (any |ν|>10 imaginary or <6 rigid modes → `RuntimeError`).
- `pyBall/nanocrystal_pipeline.py` — `mmff_assert_stationary` before Hessian; `solve_normal_modes_from_hessian_npz` asserts on unprojected `K`.
- `CompChemUtils/examples/tSiNCs/mmff_molecular_session.py` (`relax`, `assert_stationary`, `mol2_path`).
- `tests/tSiNCs/run_mmff_ownmin_l2.py` — one crystal per process; `freq_report` uses the shared assert (no `n_imag > 80` slop).
- `tests/tSiNCs/test_spectrum_at_minimum.py` — synthetic + L2 own-min must pass, `WRONG_at_DFTB_geometry` must fail.
- `tests/tSiNCs/fit_mmff_kss_pyscf.py` — own-min stretch RMSE vs PBE L1 (`joint` / `nbscan` / `morse2d` / `anneal`); report [`MMFF_C_CH_vs_CH2_kfit.md`](MMFF_C_CH_vs_CH2_kfit.md).
- Agent skill: recurring mistake in scientific / forcefield-validation skills.

## PySCF L1 (same rule, 2026-09-03)

Canonical jobs: `/home/prokop/SIMULATIONS/SiNCs/pyscf_vib_results` (PBE/ccECP-cc-pVDZ). **10 of 12** are minima. Still **not** spectra: `cube_C` (1 imag, Td saddle, −82.5 cm⁻¹) and `octahedron_7ring_Si` (1 imag, tip saddle, −25.8 cm⁻¹) — mode-follow jobs submitted. The other five Si (tight FMAX=0.001) and five C (luna) **are** valid own-min references. Legacy `.../pySCF/jobs/results` Si (5–11 imag) are superseded. Neighborhood PDOS: [`PySCF_L1_neighborhood_PDOS.md`](PySCF_L1_neighborhood_PDOS.md). Do not `sqrt(max(λ,0))`.

DFTB+ L1 (`/home/prokop/SIMULATIONS/SiNCs/DFTB/L1`, 24 jobs, C: 3ob/mio, Si: pbc/matsci) is **at a minimum** (`n_imag=0` in every `status.json`). Comparison plots: `tests/tSiNCs/OUT_dftb_vs_pyscf_l1/`. Still own-minimum: do not mix DFTB modes with a PySCF geometry.
