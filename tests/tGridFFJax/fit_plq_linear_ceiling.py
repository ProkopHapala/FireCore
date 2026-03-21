#!/usr/bin/env python3
"""PLQ Theoretical Ceiling: Can ANY linear PLQ coefficients reproduce MACE?

Uses existing MACE teacher data from damping sweep to solve the linear system:
    E_int(k) = Σ_i [ A_{t(i)}·V_P(r_i^k) + B_{t(i)}·V_L(r_i^k) + C_{t(i)}·V_Q(r_i^k) ]

4 elements {C, H, N, O} × 3 fields {P, L, Q} → 12 parameters.
Solves via numpy.linalg.lstsq (linear least squares — global optimum guaranteed).

If the optimal linear RMSE is close to current GridFF RMSE (~37 meV):
   → PLQ functional form IS the ceiling; need new physics (pairwise C₆, image, etc.)
If much better (e.g. < 15 meV):
   → PLQ CAN do better; the REQ Morse parameterization is suboptimal.
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

OUT_DIR = ensure_dir(Path(__file__).resolve().parent / "plq_ceiling_analysis")
TEACHER_CACHE = Path(__file__).resolve().parent / "damping_sweep" / "teacher_cache.npz"
FIT_PARAMS = Path(__file__).resolve().parent / "gridff_proof_run_damped_optimal" / "fit_params.json"

ELEMENTS = ['C', 'H', 'N', 'O']
EL_MAP = {e: i for i, e in enumerate(ELEMENTS)}
N_EL = len(ELEMENTS)


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
    print("  PLQ THEORETICAL CEILING: LINEAR REGRESSION ON MACE DATA")
    print("=" * 70)
    print(f"Output: {OUT_DIR}")

    # Load teacher
    if not TEACHER_CACHE.exists():
        raise FileNotFoundError(f"Run sweep_london_damping.py first: {TEACHER_CACHE}")
    tc = np.load(TEACHER_CACHE)
    teacher_energies = tc["energies"]
    print(f"Teacher cache: {len(teacher_energies)} poses, {TEACHER_CACHE}")

    # Load current best fit params for comparison
    with open(FIT_PARAMS) as f:
        current_params = json.load(f)

    # Build model
    config = _make_config(d0=3.70, width=0.35)
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

    reactive_elements = tuple(sorted(set(adsorbate.symbols)))
    model = HybridGridFFModel(
        density=density, reactive_elements=reactive_elements,
        toggles=config.toggles, grid_config=config.grid,
        hybrid_config=config.hybrid_model, prefer_jax=True,
    )

    # Also get current GridFF predictions for comparison
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
    gridff_energies = np.asarray(student.energies)

    # ============================================================
    # BUILD DESIGN MATRIX: sample P, L, Q at each atom position
    # ============================================================
    print("\nBuilding design matrix (sampling PLQ grids at atom positions)...", flush=True)

    n_poses = len(poses.positions)
    n_mol = len(adsorbate.symbols)
    X = np.zeros((n_poses, 3 * N_EL))  # 12 columns: [P_C, P_H, P_N, P_O, L_C, ..., Q_O]

    # Get the PLQ grids — build_substrate_grids returns a dict with raw and bspline grids
    from pyBall.gridff_jax.substrate_fields import build_substrate_grids
    from scipy.ndimage import map_coordinates

    grids = build_substrate_grids(density, reactive_elements, grid_config=config.grid,
                                   toggles=config.toggles)
    # Use pre-filtered B-spline coefficients directly (already computed by build_substrate_grids)
    plq_bspline = np.stack([grids["pauli_coeff_zyx"],
                             grids["london_coeff_zyx"],
                             grids["coulomb_coeff_zyx"]])  # (3, nz, ny, nx)
    print(f"  PLQ B-spline coeff shape: {plq_bspline.shape}")

    # Grid coordinate system
    voxel = density.voxel
    if voxel.ndim == 2:
        voxel_diag = np.array([voxel[0, 0], voxel[1, 1], voxel[2, 2]])
    else:
        voxel_diag = voxel
    origin = density.origin

    print(f"  Voxel: {voxel_diag}")
    print(f"  Origin: {origin}")
    print(f"  Grid size: {plq_bspline.shape[1:]}")

    # Sample at each atom position for each pose
    t0 = time.time()
    for k in range(n_poses):
        mol_pos = poses.positions[k]  # (n_mol, 3)
        for i in range(n_mol):
            el = adsorbate.symbols[i]
            j = EL_MAP[el]
            pos = mol_pos[i]
            # Fractional grid coordinates
            frac = (pos - origin) / voxel_diag
            for ch in range(3):  # P=0, L=1, Q=2
                val = map_coordinates(plq_bspline[ch], frac.reshape(3, 1),
                                      order=3, mode='wrap', prefilter=False)[0]
                X[k, ch * N_EL + j] += val

    dt = time.time() - t0
    print(f"  Design matrix built in {dt:.1f}s ({n_poses} poses × {n_mol} atoms)")

    # ============================================================
    # LINEAR LEAST SQUARES: find optimal coefficients
    # ============================================================
    y = teacher_energies  # MACE interaction energies (eV)

    # Filter to z >= 2.0 for physical relevance
    site_names = list(site_uvs.keys())
    z_mask = np.zeros(n_poses, dtype=bool)
    for si in range(len(site_names)):
        lo, hi = si * nz, (si + 1) * nz
        z_mask[lo:hi] = z_values >= 2.0

    X_fit = X[z_mask]
    y_fit = y[z_mask]

    print(f"\n  Fitting {X_fit.shape[0]} poses (z >= 2.0 A), {X_fit.shape[1]} parameters")

    # Solve: theta_opt minimizes ||X @ theta - y||²
    theta_opt, residuals, rank, sv = np.linalg.lstsq(X_fit, y_fit, rcond=None)

    y_pred_linear = X_fit @ theta_opt
    rmse_linear = np.sqrt(np.mean((y_pred_linear - y_fit) ** 2)) * 1000  # meV
    r2_linear = 1 - np.sum((y_pred_linear - y_fit) ** 2) / np.sum((y_fit - y_fit.mean()) ** 2)

    # Current GridFF comparison
    gridff_fit = gridff_energies[z_mask]
    rmse_gridff = np.sqrt(np.mean((gridff_fit - y_fit) ** 2)) * 1000
    r2_gridff = 1 - np.sum((gridff_fit - y_fit) ** 2) / np.sum((y_fit - y_fit.mean()) ** 2)

    # Also try with a constant offset term (13 params)
    X_fit_const = np.column_stack([X_fit, np.ones(len(X_fit))])
    theta_opt_c, _, _, _ = np.linalg.lstsq(X_fit_const, y_fit, rcond=None)
    y_pred_const = X_fit_const @ theta_opt_c
    rmse_const = np.sqrt(np.mean((y_pred_const - y_fit) ** 2)) * 1000

    print("\n" + "=" * 70)
    print("  RESULTS: PLQ THEORETICAL CEILING")
    print("=" * 70)
    print(f"\n  Current GridFF (REQ Morse):  RMSE = {rmse_gridff:.1f} meV  R² = {r2_gridff:.4f}")
    print(f"  Optimal linear PLQ (12p):   RMSE = {rmse_linear:.1f} meV  R² = {r2_linear:.4f}")
    print(f"  Optimal linear PLQ + offset: RMSE = {rmse_const:.1f} meV")
    print(f"\n  Improvement potential: {rmse_gridff:.1f} → {rmse_linear:.1f} meV "
          f"({(1 - rmse_linear/rmse_gridff)*100:+.1f}%)")

    # Print optimal coefficients
    print("\n  Optimal linear coefficients:")
    print(f"  {'Element':>8s}  {'A (Pauli)':>12s}  {'B (London)':>12s}  {'C (Coulomb)':>12s}")
    for j, el in enumerate(ELEMENTS):
        print(f"  {el:>8s}  {theta_opt[j]:12.6f}  {theta_opt[N_EL+j]:12.6f}  {theta_opt[2*N_EL+j]:12.6f}")

    if X_fit_const.shape[1] > X_fit.shape[1]:
        print(f"\n  Constant offset: {theta_opt_c[-1]*1000:.1f} meV")

    # ============================================================
    # PER-SITE ANALYSIS
    # ============================================================
    print("\n" + "=" * 70)
    print("  PER-SITE ANALYSIS")
    print("=" * 70)

    # Full predictions (including z < 2.0 for plotting)
    y_pred_full = X @ theta_opt

    for si, site in enumerate(site_names):
        lo, hi = si * nz, (si + 1) * nz
        mask = z_values >= 2.0

        te_s = y[lo:hi]
        gf_s = gridff_energies[lo:hi]
        lp_s = y_pred_full[lo:hi]

        rmse_gf = np.sqrt(np.mean((gf_s[mask] - te_s[mask]) ** 2)) * 1000
        rmse_lp = np.sqrt(np.mean((lp_s[mask] - te_s[mask]) ** 2)) * 1000

        well_i_mace = np.argmin(te_s[mask])
        well_i_gf = np.argmin(gf_s[mask])
        well_i_lp = np.argmin(lp_s[mask])

        z_phys = z_values[mask]

        print(f"\n  {site.upper()} site:")
        print(f"    MACE well:       {te_s[mask][well_i_mace]*1000:.1f} meV at z={z_phys[well_i_mace]:.2f} A")
        print(f"    GridFF well:     {gf_s[mask][well_i_gf]*1000:.1f} meV at z={z_phys[well_i_gf]:.2f} A")
        print(f"    Linear PLQ well: {lp_s[mask][well_i_lp]*1000:.1f} meV at z={z_phys[well_i_lp]:.2f} A")
        print(f"    GridFF RMSE:     {rmse_gf:.1f} meV")
        print(f"    Linear PLQ RMSE: {rmse_lp:.1f} meV")

    # ============================================================
    # ALSO TEST WITHOUT DAMPING
    # ============================================================
    print("\n" + "=" * 70)
    print("  COMPARISON: WITH vs WITHOUT LONDON DAMPING")
    print("=" * 70)

    config_nodamp = _make_config(d0=0.0)
    density_nodamp = make_density_backend(config_nodamp.density_backend,
                                           surface=config_nodamp.surface,
                                           grid=config_nodamp.grid).load()
    grids_nodamp = build_substrate_grids(density_nodamp, reactive_elements,
                                          grid_config=config_nodamp.grid,
                                          toggles=config_nodamp.toggles)
    plq_bspline_nd = np.stack([grids_nodamp["pauli_coeff_zyx"],
                                grids_nodamp["london_coeff_zyx"],
                                grids_nodamp["coulomb_coeff_zyx"]])

    X_nd = np.zeros((n_poses, 3 * N_EL))
    for k in range(n_poses):
        mol_pos = poses.positions[k]
        for i in range(n_mol):
            el = adsorbate.symbols[i]
            j = EL_MAP[el]
            pos = mol_pos[i]
            frac = (pos - origin) / voxel_diag
            for ch in range(3):
                val = map_coordinates(plq_bspline_nd[ch], frac.reshape(3, 1),
                                      order=3, mode='wrap', prefilter=False)[0]
                X_nd[k, ch * N_EL + j] += val

    X_nd_fit = X_nd[z_mask]
    theta_nd, _, _, _ = np.linalg.lstsq(X_nd_fit, y_fit, rcond=None)
    y_pred_nd = X_nd_fit @ theta_nd
    rmse_nd = np.sqrt(np.mean((y_pred_nd - y_fit) ** 2)) * 1000

    print(f"\n  Optimal linear PLQ WITHOUT damping: RMSE = {rmse_nd:.1f} meV")
    print(f"  Optimal linear PLQ WITH damping:    RMSE = {rmse_linear:.1f} meV")
    print(f"  Improvement from damping:           {rmse_nd - rmse_linear:+.1f} meV")

    # ============================================================
    # RIDGE REGRESSION (regularized) to check overfitting
    # ============================================================
    # Leave-one-site-out cross-validation
    print("\n  Cross-validation (leave-one-site-out):")
    for si_test, site_test in enumerate(site_names):
        lo_test, hi_test = si_test * nz, (si_test + 1) * nz
        mask_test = z_values >= 2.0
        test_idx = np.arange(lo_test, hi_test)[mask_test]

        all_idx = np.where(z_mask)[0]
        train_idx = np.array([i for i in all_idx if i not in test_idx])

        theta_cv, _, _, _ = np.linalg.lstsq(X[train_idx], y[train_idx], rcond=None)
        y_cv_pred = X[test_idx] @ theta_cv
        rmse_cv = np.sqrt(np.mean((y_cv_pred - y[test_idx]) ** 2)) * 1000
        print(f"    Test={site_test:8s}: RMSE = {rmse_cv:.1f} meV (trained on other sites)")

    # ============================================================
    # PLOTS
    # ============================================================
    print("\nGenerating plots...", flush=True)

    # Plot 1: Three curves comparison (MACE, GridFF, optimal linear PLQ)
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    for si, site in enumerate(site_names):
        lo, hi = si * nz, (si + 1) * nz
        z = z_values
        mask = z >= 2.0

        ax = axes[0, si]
        ax.plot(z[mask], y[lo:hi][mask] * 1000, "k-", lw=2.5, label="MACE teacher")
        ax.plot(z[mask], gridff_energies[lo:hi][mask] * 1000, "b--", lw=1.5, label=f"GridFF REQ ({rmse_gridff:.0f} meV)")
        ax.plot(z[mask], y_pred_full[lo:hi][mask] * 1000, "r-", lw=1.5, label=f"Optimal linear PLQ ({rmse_linear:.0f} meV)")
        ax.set_xlabel("z [A]")
        ax.set_ylabel("E_int [meV]")
        ax.set_title(f"{site.capitalize()}: Energy")
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)
        ax.axhline(0, color="gray", lw=0.5)

        # Residuals
        ax2 = axes[1, si]
        ax2.plot(z[mask], (gridff_energies[lo:hi][mask] - y[lo:hi][mask]) * 1000,
                 "b--", lw=1.2, label="GridFF - MACE")
        ax2.plot(z[mask], (y_pred_full[lo:hi][mask] - y[lo:hi][mask]) * 1000,
                 "r-", lw=1.2, label="Linear PLQ - MACE")
        ax2.axhline(0, color="gray", lw=0.5)
        ax2.set_xlabel("z [A]")
        ax2.set_ylabel("Residual [meV]")
        ax2.set_title(f"{site.capitalize()}: Residuals")
        ax2.legend(fontsize=7)
        ax2.grid(True, alpha=0.3)

    plt.suptitle(f"PLQ Theoretical Ceiling: CHONH2/Ag(111)\n"
                 f"GridFF REQ: {rmse_gridff:.1f} meV | Optimal linear PLQ: {rmse_linear:.1f} meV",
                 fontsize=12)
    plt.tight_layout()
    plt.savefig(OUT_DIR / "plot1_ceiling_comparison.png", dpi=200)
    plt.close()

    # Plot 2: Parity plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    colors = {"top": "tab:blue", "bridge": "tab:green", "hollow": "tab:red"}

    for si, site in enumerate(site_names):
        lo, hi = si * nz, (si + 1) * nz
        mask = z_values >= 2.0

        axes[0].scatter(y[lo:hi][mask] * 1000, gridff_energies[lo:hi][mask] * 1000,
                        c=colors[site], s=10, alpha=0.6, label=site)
        axes[1].scatter(y[lo:hi][mask] * 1000, y_pred_full[lo:hi][mask] * 1000,
                        c=colors[site], s=10, alpha=0.6, label=site)

    for ax, title, rmse, r2 in [(axes[0], "GridFF REQ Morse", rmse_gridff, r2_gridff),
                                  (axes[1], "Optimal Linear PLQ", rmse_linear, r2_linear)]:
        lims = [-500, 1500]
        ax.plot(lims, lims, "k--", lw=0.5, alpha=0.5)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.set_xlabel("MACE E_int [meV]")
        ax.set_ylabel(f"{title} E_int [meV]")
        ax.set_title(f"{title}\nRMSE={rmse:.1f} meV, R²={r2:.4f}")
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_aspect("equal")

    plt.tight_layout()
    plt.savefig(OUT_DIR / "plot2_parity.png", dpi=200)
    plt.close()

    # Plot 3: Coefficient comparison
    fig, ax = plt.subplots(figsize=(10, 5))
    x_pos = np.arange(N_EL)
    w = 0.25

    for ch_i, (ch_name, color) in enumerate(zip(["Pauli (A)", "London (B)", "Coulomb (C)"],
                                                  ["tab:red", "tab:green", "tab:purple"])):
        vals = theta_opt[ch_i * N_EL:(ch_i + 1) * N_EL]
        ax.bar(x_pos + ch_i * w, vals, w, label=ch_name, color=color, alpha=0.7)

    ax.set_xticks(x_pos + w)
    ax.set_xticklabels(ELEMENTS)
    ax.set_ylabel("Coefficient value")
    ax.set_title("Optimal Linear PLQ Coefficients per Element")
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    ax.axhline(0, color="gray", lw=0.5)

    plt.tight_layout()
    plt.savefig(OUT_DIR / "plot3_coefficients.png", dpi=200)
    plt.close()

    # ============================================================
    # SUMMARY
    # ============================================================
    print("\n" + "=" * 70)
    print("  CONCLUSION")
    print("=" * 70)

    if rmse_linear < rmse_gridff * 0.7:
        print(f"\n  PLQ CAN do significantly better: {rmse_gridff:.1f} → {rmse_linear:.1f} meV")
        print(f"  The REQ Morse parameterization is suboptimal.")
        print(f"  Direct linear coefficients achieve {(1-rmse_linear/rmse_gridff)*100:.0f}% lower RMSE.")
        print(f"  → Consider using direct PLQ coefficients instead of REQ coupling.")
    elif rmse_linear < rmse_gridff * 0.9:
        print(f"\n  Modest improvement possible: {rmse_gridff:.1f} → {rmse_linear:.1f} meV")
        print(f"  The REQ coupling captures most of what PLQ can do.")
        print(f"  Remaining error is partly in the PLQ functional form.")
    else:
        print(f"\n  PLQ IS AT ITS CEILING: {rmse_gridff:.1f} ≈ {rmse_linear:.1f} meV")
        print(f"  The REQ Morse parameterization is near-optimal for PLQ.")
        print(f"  Further improvement REQUIRES new physics beyond P, L, Q fields.")
        print(f"  Candidates: pairwise C₆/R⁶, image charges, site-specific corrections.")

    print(f"\n  All plots saved in: {OUT_DIR}")

    # Save results as JSON
    results = {
        "rmse_gridff_meV": float(rmse_gridff),
        "rmse_linear_plq_meV": float(rmse_linear),
        "rmse_linear_plq_offset_meV": float(rmse_const),
        "rmse_linear_no_damping_meV": float(rmse_nd),
        "r2_gridff": float(r2_gridff),
        "r2_linear_plq": float(r2_linear),
        "optimal_coefficients": {
            el: {"pauli": float(theta_opt[j]),
                 "london": float(theta_opt[N_EL + j]),
                 "coulomb": float(theta_opt[2*N_EL + j])}
            for j, el in enumerate(ELEMENTS)
        },
        "n_poses_fitted": int(X_fit.shape[0]),
        "london_damping": {"d0": 3.70, "width": 0.35},
    }
    with open(OUT_DIR / "ceiling_results.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"  Results saved: {OUT_DIR / 'ceiling_results.json'}")


if __name__ == "__main__":
    main()
