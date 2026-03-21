#!/usr/bin/env python3
"""Proof-of-concept: Can pairwise C₆/R⁶ dispersion fix the hollow-site error?

The PLQ fields are laterally smooth → same well depth at all sites.
Pairwise 1/r⁶ naturally counts more Ag neighbors at hollow → deeper well.

This script:
1. Computes a pairwise 1/r⁶ dispersion field over the Ag slab (with PBC)
2. Samples it at each molecular atom position for all poses
3. Fits: E_MACE ≈ E_GridFF_current + Σ_i D_el × V_disp(r_i) + offset
4. Reports: does this site-differentiated correction reduce RMSE?

No pyBall code is modified — this is a standalone diagnostic.
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
from pyBall.gridff_jax.utils import ensure_dir

OUT_DIR = ensure_dir(Path(__file__).resolve().parent / "c6_correction_test")
TEACHER_CACHE = Path(__file__).resolve().parent / "damping_sweep" / "teacher_cache.npz"
FIT_PARAMS = Path(__file__).resolve().parent / "gridff_proof_run_damped_optimal" / "fit_params.json"

ELEMENTS = ['C', 'H', 'N', 'O']
EL_MAP = {e: i for i, e in enumerate(ELEMENTS)}


def compute_pairwise_dispersion_at_points(points, ag_positions, cell, pbc_counts=(2, 2, 0)):
    """Compute V_disp(r) = -Σ_j 1/|r - r_j|⁶ for each point r.

    Sums over all Ag atoms and their periodic images.
    Returns array of shape (n_points,).

    The 1/r⁶ is the universal shape of London dispersion.
    Element-specific C₆ scaling is handled by the linear fit.
    """
    n_points = len(points)
    v_disp = np.zeros(n_points)

    # Generate periodic images
    na, nb, nc = pbc_counts
    shifts = []
    for ia in range(-na, na + 1):
        for ib in range(-nb, nb + 1):
            for ic in range(-nc, nc + 1):
                shifts.append(ia * cell[0] + ib * cell[1] + ic * cell[2])
    shifts = np.array(shifts)  # (n_shifts, 3)

    # Sum over all Ag atoms and their images
    for shift in shifts:
        ag_shifted = ag_positions + shift  # (n_ag, 3)
        # Vectorized over points and atoms
        dp = points[:, None, :] - ag_shifted[None, :, :]  # (n_points, n_ag, 3)
        r2 = np.sum(dp * dp, axis=-1)  # (n_points, n_ag)
        r2 = np.maximum(r2, 1.0)  # Avoid singularities (1 A² minimum)
        r6 = r2 ** 3
        v_disp -= np.sum(1.0 / r6, axis=1)

    return v_disp


def compute_pairwise_c6_grid(ag_positions, cell, grid_shape, voxel, origin, pbc_counts=(2, 2, 0)):
    """Build a 3D grid of the pairwise 1/r⁶ dispersion potential.

    Returns grid of shape (nz, ny, nx).
    """
    nz, ny, nx = grid_shape
    print(f"  Building C₆/R⁶ grid: {nz}×{ny}×{nx} = {nz*ny*nx:,} voxels", flush=True)

    # Generate all grid points
    iz, iy, ix = np.mgrid[0:nz, 0:ny, 0:nx]
    points = np.stack([
        ix.ravel() * voxel[0] + origin[0],
        iy.ravel() * voxel[1] + origin[1],
        iz.ravel() * voxel[2] + origin[2],
    ], axis=-1)  # (N, 3)

    # Process in chunks to avoid memory issues
    chunk_size = 50000
    n_total = len(points)
    v_disp = np.zeros(n_total)

    t0 = time.time()
    for start in range(0, n_total, chunk_size):
        end = min(start + chunk_size, n_total)
        v_disp[start:end] = compute_pairwise_dispersion_at_points(
            points[start:end], ag_positions, cell, pbc_counts
        )
        if (start // chunk_size) % 20 == 0:
            elapsed = time.time() - t0
            pct = end / n_total * 100
            print(f"    {pct:.0f}% done ({elapsed:.1f}s)", flush=True)

    dt = time.time() - t0
    print(f"  Grid built in {dt:.1f}s")

    return v_disp.reshape(nz, ny, nx)


def _make_config(d0=3.70, width=0.35):
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


def main():
    print("=" * 70)
    print("  PAIRWISE C₆/R⁶ CORRECTION: PROOF OF CONCEPT")
    print("=" * 70)
    print(f"Output: {OUT_DIR}")

    # Load teacher (MACE)
    tc = np.load(TEACHER_CACHE)
    teacher_e = tc["energies"]
    print(f"Teacher: {len(teacher_e)} poses")

    # Load current best fit params
    with open(FIT_PARAMS) as f:
        current_params = json.load(f)

    # Build density and model
    config = _make_config()
    print("Loading density...", flush=True)
    density = make_density_backend(config.density_backend, surface=config.surface, grid=config.grid).load()

    adsorbate = get_adsorbate(name="CHONH2")
    site_uvs = {}
    for site in infer_surface_sites(density):
        if site.label in ["top", "bridge", "hollow"]:
            site_uvs[site.label] = np.asarray(site.uv)

    z_values = np.linspace(1.8, 8.0, 63)
    nz = len(z_values)
    poses = _build_poses(density, adsorbate, site_uvs, z_values)

    # Get current GridFF predictions
    reactive_elements = tuple(sorted(set(adsorbate.symbols)))
    model = HybridGridFFModel(
        density=density, reactive_elements=reactive_elements,
        toggles=config.toggles, grid_config=config.grid,
        hybrid_config=config.hybrid_model, prefer_jax=True,
    )
    params = HybridParameters(
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
    student = model.evaluate_batch(poses, params=params, compute_forces=True)
    gridff_e = np.asarray(student.energies)

    # ============================================================
    # COMPUTE PAIRWISE 1/R⁶ AT MOLECULAR ATOM POSITIONS
    # ============================================================
    print("\nComputing pairwise 1/r⁶ dispersion at atom positions...", flush=True)

    # Get Ag positions and cell
    ag_positions = density.positions  # (n_ag, 3)
    cell = density.cell  # (3, 3) lattice vectors
    print(f"  Ag atoms: {len(ag_positions)}")
    print(f"  Cell: {cell[0]}")

    # For each pose, compute 1/r⁶ sum at each molecular atom position
    n_poses = len(poses.positions)
    n_mol = len(adsorbate.symbols)
    n_el = len(ELEMENTS)

    # Design matrix: per-element sum of 1/r⁶ at atom positions
    X_c6 = np.zeros((n_poses, n_el))  # 4 columns: D_C, D_H, D_N, D_O

    t0 = time.time()
    for k in range(n_poses):
        mol_pos = poses.positions[k]  # (n_mol, 3)
        for i in range(n_mol):
            el = adsorbate.symbols[i]
            j = EL_MAP[el]
            pos = mol_pos[i:i+1]  # (1, 3)
            v = compute_pairwise_dispersion_at_points(pos, ag_positions, cell, pbc_counts=(2, 2, 0))
            X_c6[k, j] += v[0]

    dt = time.time() - t0
    print(f"  Computed in {dt:.1f}s ({n_poses} poses × {n_mol} atoms)")

    # ============================================================
    # LINEAR FIT: E_MACE = E_GridFF + X_c6 @ D + offset
    # ============================================================
    print("\n" + "=" * 70)
    print("  FIT: GridFF + C₆ CORRECTION")
    print("=" * 70)

    # Filter to z >= 2.0
    site_names = list(site_uvs.keys())
    z_mask = np.zeros(n_poses, dtype=bool)
    for si in range(len(site_names)):
        lo, hi = si * nz, (si + 1) * nz
        z_mask[lo:hi] = z_values >= 2.0

    # Residual: what GridFF misses
    residual = teacher_e[z_mask] - gridff_e[z_mask]

    # Fit: residual ≈ X_c6 @ D + offset
    X_fit = np.column_stack([X_c6[z_mask], np.ones(z_mask.sum())])  # 5 params
    D_opt, _, _, _ = np.linalg.lstsq(X_fit, residual, rcond=None)

    corrected_e = gridff_e.copy()
    X_full = np.column_stack([X_c6, np.ones(n_poses)])
    correction = X_full @ D_opt
    corrected_e += correction

    # Compute metrics
    rmse_before = np.sqrt(np.mean((gridff_e[z_mask] - teacher_e[z_mask]) ** 2)) * 1000
    rmse_after = np.sqrt(np.mean((corrected_e[z_mask] - teacher_e[z_mask]) ** 2)) * 1000

    print(f"\n  Before correction (GridFF only):  RMSE = {rmse_before:.1f} meV")
    print(f"  After C₆ correction:              RMSE = {rmse_after:.1f} meV")
    print(f"  Improvement: {rmse_before - rmse_after:.1f} meV ({(1-rmse_after/rmse_before)*100:.1f}%)")

    print(f"\n  Fitted correction coefficients:")
    for j, el in enumerate(ELEMENTS):
        print(f"    D_{el}: {D_opt[j]:.6e}")
    print(f"    Offset: {D_opt[-1]*1000:.1f} meV")

    # Per-site analysis
    print("\n  Per-site comparison:")
    print(f"  {'Site':>8s}  {'Before':>10s}  {'After':>10s}  {'Improve':>10s}  {'Well_MACE':>10s}  {'Well_corr':>10s}")
    for si, site in enumerate(site_names):
        lo, hi = si * nz, (si + 1) * nz
        mask = z_values >= 2.0

        te_s = teacher_e[lo:hi]
        gf_s = gridff_e[lo:hi]
        co_s = corrected_e[lo:hi]

        rmse_b = np.sqrt(np.mean((gf_s[mask] - te_s[mask]) ** 2)) * 1000
        rmse_a = np.sqrt(np.mean((co_s[mask] - te_s[mask]) ** 2)) * 1000

        well_mace = np.min(te_s[mask]) * 1000
        well_corr = np.min(co_s[mask]) * 1000

        print(f"  {site:>8s}  {rmse_b:8.1f} meV  {rmse_a:8.1f} meV  {rmse_b-rmse_a:+8.1f} meV"
              f"  {well_mace:8.1f} meV  {well_corr:8.1f} meV")

    # Check: are C₆ field values site-dependent?
    print("\n  C₆ field site-dependence at well (z=2.7 A):")
    well_z_idx = np.argmin(np.abs(z_values - 2.7))
    for si, site in enumerate(site_names):
        k = si * nz + well_z_idx
        c6_total = np.sum(X_c6[k])
        print(f"    {site:>8s}: V_disp_total = {c6_total:.4e}  "
              f"(per el: C={X_c6[k,0]:.3e} H={X_c6[k,1]:.3e} N={X_c6[k,2]:.3e} O={X_c6[k,3]:.3e})")

    # Z-region analysis
    print("\n  Z-region analysis:")
    for zlo, zhi in [(2.0, 2.5), (2.5, 3.0), (3.0, 4.0), (4.0, 6.0), (6.0, 8.0)]:
        m_all = np.zeros(n_poses, dtype=bool)
        for si in range(len(site_names)):
            lo, hi = si * nz, (si + 1) * nz
            m_z = (z_values >= zlo) & (z_values < zhi)
            m_all[lo:hi] = m_z
        if m_all.sum() == 0:
            continue
        rmse_b = np.sqrt(np.mean((gridff_e[m_all] - teacher_e[m_all]) ** 2)) * 1000
        rmse_a = np.sqrt(np.mean((corrected_e[m_all] - teacher_e[m_all]) ** 2)) * 1000
        print(f"    z=[{zlo:.0f},{zhi:.0f}): Before={rmse_b:6.1f}  After={rmse_a:6.1f}  Δ={rmse_b-rmse_a:+6.1f} meV")

    # ============================================================
    # PLOT
    # ============================================================
    print("\nGenerating plots...", flush=True)

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    for si, site in enumerate(site_names):
        lo, hi = si * nz, (si + 1) * nz
        z = z_values
        mask = z >= 2.0

        ax = axes[0, si]
        ax.plot(z[mask], teacher_e[lo:hi][mask] * 1000, "k-", lw=2.5, label="MACE")
        ax.plot(z[mask], gridff_e[lo:hi][mask] * 1000, "b--", lw=1.5, label="GridFF (before)")
        ax.plot(z[mask], corrected_e[lo:hi][mask] * 1000, "r-", lw=1.5, label="GridFF + C₆")
        ax.set_xlabel("z [A]")
        ax.set_ylabel("E_int [meV]")
        ax.set_title(f"{site.capitalize()}")
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)
        ax.axhline(0, color="gray", lw=0.5)

        ax2 = axes[1, si]
        ax2.plot(z[mask], (gridff_e[lo:hi][mask] - teacher_e[lo:hi][mask]) * 1000,
                 "b--", lw=1.2, label="Before")
        ax2.plot(z[mask], (corrected_e[lo:hi][mask] - teacher_e[lo:hi][mask]) * 1000,
                 "r-", lw=1.2, label="After C₆ correction")
        ax2.axhline(0, color="gray", lw=0.5)
        ax2.set_xlabel("z [A]")
        ax2.set_ylabel("Residual [meV]")
        ax2.set_title(f"{site.capitalize()}: Residuals")
        ax2.legend(fontsize=7)
        ax2.grid(True, alpha=0.3)

    plt.suptitle(f"Pairwise C₆/R⁶ Correction: CHONH2/Ag(111)\n"
                 f"Before: {rmse_before:.1f} meV → After: {rmse_after:.1f} meV "
                 f"({(1-rmse_after/rmse_before)*100:+.1f}%)", fontsize=12)
    plt.tight_layout()
    plt.savefig(OUT_DIR / "c6_correction_comparison.png", dpi=200)
    plt.close()

    # Plot C₆ field values per site
    fig, ax = plt.subplots(figsize=(10, 6))
    for si, site in enumerate(site_names):
        lo, hi = si * nz, (si + 1) * nz
        c6_sum = np.sum(X_c6[lo:hi], axis=1)  # Sum over all atoms
        ax.plot(z_values, c6_sum, "-", lw=2, label=f"{site}")
    ax.set_xlabel("z above Ag(111) [A]")
    ax.set_ylabel("V_disp = -Σ 1/r⁶ [A⁻⁶]")
    ax.set_title("Pairwise 1/r⁶ Dispersion Field: Site Dependence")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xlim(2.0, 6.0)
    plt.tight_layout()
    plt.savefig(OUT_DIR / "c6_site_dependence.png", dpi=200)
    plt.close()

    # Save results
    results = {
        "rmse_before_meV": float(rmse_before),
        "rmse_after_meV": float(rmse_after),
        "improvement_pct": float((1 - rmse_after / rmse_before) * 100),
        "correction_coefficients": {el: float(D_opt[j]) for j, el in enumerate(ELEMENTS)},
        "correction_offset_meV": float(D_opt[-1] * 1000),
    }
    with open(OUT_DIR / "c6_results.json", "w") as f:
        json.dump(results, f, indent=2)

    print(f"\nPlots and results saved in: {OUT_DIR}")


if __name__ == "__main__":
    main()
