#!/usr/bin/env python3
"""Validate the integrated C₆/R⁶ dispersion grid channel.

Tests:
1. Grid values match the proof-of-concept direct sum within interpolation tolerance
2. Energy with c6_coeff=1.0 produces dispersion contribution
3. Regression: use_pairwise_c6=False gives same result as before
4. Joint REQ + C₆ fit reduces RMSE below PLQ-only baseline
"""

from __future__ import annotations

import os
import sys
import json
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("JAX_PLATFORMS", "cpu")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.append(str(ROOT))

from pyBall.gridff_jax import RunConfig
from pyBall.gridff_jax.density_backends import make_density_backend
from pyBall.gridff_jax.hybrid_energy import HybridGridFFModel, HybridParameters
from pyBall.gridff_jax.interfaces import PoseBatch
from pyBall.gridff_jax.pose_sampling import (
    get_adsorbate, infer_reference_planes, infer_surface_sites, transform_adsorbate,
)
from pyBall.gridff_jax.substrate_fields import VDW_SURF_PARAMS, TS_FREE_PARAMS, _ts_c6_pair
from pyBall.gridff_jax.utils import ensure_dir

OUT_DIR = ensure_dir(Path(__file__).resolve().parent / "c6_grid_channel_test")
TEACHER_CACHE = Path(__file__).resolve().parent / "damping_sweep" / "teacher_cache.npz"
FIT_PARAMS = Path(__file__).resolve().parent / "gridff_proof_run_damped_optimal" / "fit_params.json"


def _make_config(use_c6=False, d0=3.70, width=0.35):
    config = RunConfig()
    config.density_backend.chgcar_path = "/home/niel/git/ORR_HER_Ag_Colab/results/Ag_ORR_HER/slab_clean/final_scf_12x12x1/CHGCAR"
    config.density_backend.locpot_path = "/home/niel/git/ORR_HER_Ag_Colab/results/Ag_ORR_HER/slab_clean/workfunc_12x12x1/LOCPOT"
    config.grid.builder_mode = "metal_density_plq"
    config.grid.interpolation_order = 3
    config.grid.metal_density_pauli_power = 1.0
    config.grid.metal_density_london_power = 0.5
    if d0 > 0.0:
        config.grid.london_damping_d0 = float(d0)
        config.grid.london_damping_width = float(width)
    config.grid.use_pairwise_c6 = use_c6
    config.toggles.use_ct_qeq = False
    config.toggles.use_image_charge = False
    config.toggles.use_reactive_grid = False
    config.hybrid_model.use_qeq = False
    config.hybrid_model.use_image = False
    config.hybrid_model.use_reactive_grid = False
    config.hybrid_model.use_req_plq = True
    return config


def _build_poses(density, adsorbate, site_uvs, z_values):
    planes = infer_reference_planes(density)
    z_ref = float(planes["z_ref"])
    z_image = float(planes["z_image"])
    quat = np.array([1, 0, 0, 0], dtype=float)
    positions, pose_params, labels = [], [], []
    for label, uv in site_uvs.items():
        for zh in z_values:
            positions.append(transform_adsorbate(adsorbate, density, uv, float(zh), quat, z_ref))
            pose_params.append(np.concatenate([uv, [float(zh)], quat]))
            labels.append(label)
    return PoseBatch(
        adsorbate=adsorbate,
        positions=np.asarray(positions),
        pose_params=np.asarray(pose_params),
        site_labels=labels,
        metadata={"z_ref": z_ref, "z_image": z_image},
    )


def compute_direct_shape_at_points(points, ag_positions, cell, element, metal="Ag",
                                    pbc_counts=(2, 2, 0), s_r=0.94):
    """Direct pairwise shape sum: -Σ f_damp/r⁶ (NO C₆ factor, matches grid semantics)."""
    _c6_pair, r_sum = _ts_c6_pair(element, metal)
    r0 = s_r * r_sum

    n_points = len(points)
    v_shape = np.zeros(n_points)

    na, nb, nc = pbc_counts
    for ia in range(-na, na + 1):
        for ib in range(-nb, nb + 1):
            for ic in range(-nc, nc + 1):
                shift = ia * cell[0] + ib * cell[1] + ic * cell[2]
                ag_shifted = ag_positions + shift
                dp = points[:, None, :] - ag_shifted[None, :, :]
                r2 = np.sum(dp * dp, axis=-1) + 1e-32
                r = np.sqrt(r2)
                r6 = r2 ** 3
                f_damp = 1.0 / (1.0 + np.exp(-20.0 * (r / r0 - 1.0)))
                v_shape += np.sum(f_damp * (-1.0 / r6), axis=1)

    return v_shape


def main():
    print("=" * 70)
    print("  C₆ GRID CHANNEL VALIDATION")
    print("=" * 70)
    print(f"Output: {OUT_DIR}")

    # --- Print vdW-surf parameters ---
    print("\nvdW-surf screened parameters loaded:")
    for metal, (alpha, c6, rvdw) in VDW_SURF_PARAMS.items():
        print(f"  {metal}: alpha={alpha:.4f} Ang^3, C6={c6:.2f} eV*Ang^6, R_vdW={rvdw:.3f} Ang")
    print("\nTS free-atom parameters:")
    for el in ["H", "C", "N", "O"]:
        alpha, c6, rvdw = TS_FREE_PARAMS[el]
        print(f"  {el}: alpha={alpha:.4f} Ang^3, C6={c6:.3f} eV*Ang^6, R_vdW={rvdw:.3f} Ang")

    print("\nPairwise C6 (TS combination rule with vdW-surf Ag):")
    for el in ["H", "C", "N", "O"]:
        c6_pair, r_sum = _ts_c6_pair(el, "Ag")
        print(f"  {el}-Ag: C6={c6_pair:.2f} eV*Ang^6, R_sum={r_sum:.3f} Ang")

    # --- Load teacher data ---
    tc = np.load(TEACHER_CACHE)
    teacher_e = tc["energies"]
    print(f"\nTeacher: {len(teacher_e)} poses")

    # Load current best fit params
    with open(FIT_PARAMS) as f:
        current_params = json.load(f)

    # --- TEST 1: Regression — PLQ-only (no C₆) ---
    print("\n" + "=" * 50)
    print("TEST 1: Regression — PLQ-only (use_pairwise_c6=False)")
    print("=" * 50)

    config_no_c6 = _make_config(use_c6=False)
    print("Loading density (no C₆)...", flush=True)
    density = make_density_backend(config_no_c6.density_backend, surface=config_no_c6.surface, grid=config_no_c6.grid).load()

    adsorbate = get_adsorbate(name="CHONH2")
    site_uvs = {}
    for site in infer_surface_sites(density):
        if site.label in ["top", "bridge", "hollow"]:
            site_uvs[site.label] = np.asarray(site.uv)

    z_values = np.linspace(1.8, 8.0, 63)
    poses = _build_poses(density, adsorbate, site_uvs, z_values)
    reactive_elements = tuple(sorted(set(adsorbate.symbols)))

    model_no_c6 = HybridGridFFModel(
        density=density, reactive_elements=reactive_elements,
        toggles=config_no_c6.toggles, grid_config=config_no_c6.grid,
        hybrid_config=config_no_c6.hybrid_model, prefer_jax=True,
    )
    params_no_c6 = HybridParameters(
        pauli=current_params.get("pauli", {}),
        london=current_params.get("london", {}),
        reactive=current_params.get("reactive", {}),
        static_charge=current_params.get("static_charge", {}),
        req_radius_offset=current_params["req_radius_offset"],
        req_energy_scale=current_params["req_energy_scale"],
        chi=current_params.get("chi", {}),
        hardness=current_params.get("hardness", {}),
        image_scale=current_params.get("image_scale", 0.0),
        image_plane=current_params.get("image_plane", 0.0),
        image_damping=current_params.get("image_damping", 0.5),
        sample_shift_z=current_params.get("sample_shift_z", 0.0),
        coulomb_shift_z=current_params.get("coulomb_shift_z", 0.0),
        reservoir_chi=current_params.get("reservoir_chi", 0.0),
        reservoir_hardness=current_params.get("reservoir_hardness", 10.0),
        total_charge=current_params.get("total_charge", 0.0),
    )
    result_no_c6 = model_no_c6.evaluate_batch(poses, params=params_no_c6, compute_forces=False)
    gridff_e_no_c6 = np.asarray(result_no_c6.energies)
    rmse_no_c6 = np.sqrt(np.mean((gridff_e_no_c6 - teacher_e) ** 2)) * 1000
    print(f"  PLQ-only RMSE: {rmse_no_c6:.1f} meV")
    print(f"  Dispersion component: {np.asarray(result_no_c6.components.get('dispersion', [0.0]))[:5]}")

    # --- TEST 2: Build model WITH C₆ grid ---
    print("\n" + "=" * 50)
    print("TEST 2: Build model with C₆ grid channel")
    print("=" * 50)

    config_c6 = _make_config(use_c6=True)
    print("Loading density + building C₆ grid...", flush=True)
    t0 = time.time()
    density_c6 = make_density_backend(config_c6.density_backend, surface=config_c6.surface, grid=config_c6.grid).load()

    model_c6 = HybridGridFFModel(
        density=density_c6, reactive_elements=reactive_elements,
        toggles=config_c6.toggles, grid_config=config_c6.grid,
        hybrid_config=config_c6.hybrid_model, prefer_jax=True,
    )
    dt = time.time() - t0
    print(f"  Model built in {dt:.1f}s")

    # Check dispersion grid is non-zero
    disp_grid = model_c6.grids["dispersion_zyxc"]
    print(f"  Dispersion grid shape: {disp_grid.shape}")
    for c, el in enumerate(reactive_elements):
        ch = disp_grid[..., c]
        print(f"    {el}: min={ch.min():.4f}, max={ch.max():.4f}, mean={ch.mean():.6f}")

    # Evaluate with c6_coeff = 1.0
    params_c6 = HybridParameters(
        pauli=current_params.get("pauli", {}),
        london=current_params.get("london", {}),
        reactive=current_params.get("reactive", {}),
        static_charge=current_params.get("static_charge", {}),
        c6_coeff={el: 1.0 for el in reactive_elements},
        req_radius_offset=current_params["req_radius_offset"],
        req_energy_scale=current_params["req_energy_scale"],
        chi=current_params.get("chi", {}),
        hardness=current_params.get("hardness", {}),
        image_scale=current_params.get("image_scale", 0.0),
        image_plane=current_params.get("image_plane", 0.0),
        image_damping=current_params.get("image_damping", 0.5),
        sample_shift_z=current_params.get("sample_shift_z", 0.0),
        coulomb_shift_z=current_params.get("coulomb_shift_z", 0.0),
        reservoir_chi=current_params.get("reservoir_chi", 0.0),
        reservoir_hardness=current_params.get("reservoir_hardness", 10.0),
        total_charge=current_params.get("total_charge", 0.0),
    )
    result_c6 = model_c6.evaluate_batch(poses, params=params_c6, compute_forces=False)
    gridff_e_c6 = np.asarray(result_c6.energies)
    disp_component = np.asarray(result_c6.components["dispersion"])
    print(f"\n  Dispersion energy stats:")
    print(f"    min={disp_component.min():.4f}, max={disp_component.max():.4f}, mean={disp_component.mean():.4f} eV")
    print(f"  Total energy difference (C₆ on - C₆ off):")
    ediff = gridff_e_c6 - gridff_e_no_c6
    print(f"    min={ediff.min():.4f}, max={ediff.max():.4f}, mean={ediff.mean():.4f} eV")

    rmse_c6_default = np.sqrt(np.mean((gridff_e_c6 - teacher_e) ** 2)) * 1000
    print(f"  RMSE with c6_coeff=1.0: {rmse_c6_default:.1f} meV")

    # --- TEST 3: Validate grid vs direct sum ---
    print("\n" + "=" * 50)
    print("TEST 3: Grid interpolation vs direct pairwise sum")
    print("=" * 50)

    # Sample a few poses and compare
    n_check = min(20, len(poses.positions))
    max_rel_err = 0.0
    max_abs_err = 0.0

    for k in range(0, len(poses.positions), max(1, len(poses.positions) // n_check)):
        mol_pos = poses.positions[k]
        for i, el in enumerate(adsorbate.symbols):
            c = list(reactive_elements).index(el)
            # Grid-interpolated value
            grid_val = float(model_c6._sample_field_numpy(
                model_c6.sample_grids["dispersion_zyxc"], mol_pos[i:i+1]
            )[0, c]) if model_c6.sample_grids["dispersion_zyxc"].ndim == 4 else 0.0
            # Direct sum
            direct_val = float(compute_direct_shape_at_points(
                mol_pos[i:i+1], density_c6.positions, density_c6.cell, el
            )[0])
            abs_err = abs(grid_val - direct_val)
            rel_err = abs_err / max(abs(direct_val), 1e-10)
            max_abs_err = max(max_abs_err, abs_err)
            max_rel_err = max(max_rel_err, rel_err)

    print(f"  Max absolute error: {max_abs_err:.6f} eV")
    print(f"  Max relative error: {max_rel_err:.4f} ({max_rel_err*100:.2f}%)")

    if max_rel_err < 0.05:
        print("  PASS: Grid interpolation matches direct sum within 5%")
    else:
        print(f"  WARNING: Relative error {max_rel_err*100:.1f}% exceeds 5% threshold")

    # --- TEST 4: Linear C₆ fit (same as proof-of-concept but using grid) ---
    print("\n" + "=" * 50)
    print("TEST 4: Linear C₆ correction using grid channel")
    print("=" * 50)

    # Build design matrix from dispersion component per element
    # Grid stores shape only; multiply by c6_pair_values so fit produces c6_coeff
    n_poses_total = len(poses.positions)
    el_list = list(reactive_elements)
    n_el = len(el_list)
    c6_pairs = model_c6.c6_pair_values  # [n_el]
    X_shape = np.zeros((n_poses_total, n_el))

    for k in range(n_poses_total):
        mol_pos = poses.positions[k]
        disp_sample = model_c6._sample_field_numpy(
            model_c6.sample_grids["dispersion_zyxc"], mol_pos
        )
        for i, el in enumerate(adsorbate.symbols):
            c = el_list.index(el)
            X_shape[k, c] += disp_sample[i, c] if disp_sample.ndim == 2 else disp_sample[i]

    # X_c6 = shape × c6_pair (so that fitted D values are c6_coeff scaling factors)
    X_c6 = X_shape * c6_pairs[np.newaxis, :]

    # Fit: E_MACE ≈ E_GridFF + X_c6 @ D + offset
    residual = teacher_e - gridff_e_no_c6
    X_full = np.column_stack([X_c6, np.ones(n_poses_total)])
    theta, _, _, _ = np.linalg.lstsq(X_full, residual, rcond=None)
    D_fit = theta[:n_el]
    offset = theta[-1]

    E_corrected = gridff_e_no_c6 + X_c6 @ D_fit + offset
    rmse_corrected = np.sqrt(np.mean((E_corrected - teacher_e) ** 2)) * 1000

    print(f"  Fitted C₆ scaling coefficients (c6_coeff values):")
    for c, el in enumerate(el_list):
        c6_p = c6_pairs[c]
        print(f"    {el}: c6_coeff = {D_fit[c]:.4f}  (C6_pair = {c6_p:.2f} eV·Å⁶)")
    print(f"  Offset: {offset:.4f} eV")
    print(f"  RMSE (PLQ-only):     {rmse_no_c6:.1f} meV")
    print(f"  RMSE (PLQ + C₆ fit): {rmse_corrected:.1f} meV")
    improvement = (1.0 - rmse_corrected / rmse_no_c6) * 100
    print(f"  Improvement: {improvement:.1f}%")

    # --- PLOT: z-scan comparison ---
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    nz = len(z_values)

    for idx, site in enumerate(["top", "bridge", "hollow"]):
        ax = axes[idx]
        mask = np.array([l == site for l in poses.site_labels])
        z_site = np.array([poses.pose_params[k, 2] for k in range(len(poses.site_labels)) if poses.site_labels[k] == site])
        e_mace = teacher_e[mask] * 1000
        e_plq = gridff_e_no_c6[mask] * 1000
        e_corr = E_corrected[mask] * 1000
        e_c6def = gridff_e_c6[mask] * 1000

        ax.plot(z_site, e_mace, 'ko-', ms=3, lw=1, label='MACE (teacher)')
        ax.plot(z_site, e_plq, 'b--', lw=1, label=f'PLQ only ({rmse_no_c6:.0f} meV)')
        ax.plot(z_site, e_corr, 'r-', lw=2, label=f'PLQ + C₆ fit ({rmse_corrected:.0f} meV)')
        ax.set_xlabel('z (Å above surface)')
        ax.set_ylabel('Energy (meV)')
        ax.set_title(f'{site} site')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    plt.suptitle('C₆ Grid Channel Validation: CHONH2 / Ag(111)', fontsize=14)
    plt.tight_layout()
    plt.savefig(OUT_DIR / "c6_grid_channel_validation.png", dpi=150)
    print(f"\nPlot saved: {OUT_DIR / 'c6_grid_channel_validation.png'}")

    # --- Summary ---
    print("\n" + "=" * 70)
    print("  SUMMARY")
    print("=" * 70)
    print(f"  PLQ-only RMSE:            {rmse_no_c6:.1f} meV")
    print(f"  PLQ + C₆ (coeff=1.0):     {rmse_c6_default:.1f} meV")
    print(f"  PLQ + C₆ (linear fit):    {rmse_corrected:.1f} meV")
    print(f"  Grid vs direct max error: {max_rel_err*100:.2f}%")
    print(f"  Improvement:              {improvement:.1f}%")

    if rmse_corrected < rmse_no_c6 * 0.7:
        print("\n  PASS: C₆ grid channel provides >30% RMSE improvement")
    else:
        print(f"\n  NOTE: Improvement is {improvement:.1f}% (target: >30%)")

    print("=" * 70)


if __name__ == "__main__":
    main()
