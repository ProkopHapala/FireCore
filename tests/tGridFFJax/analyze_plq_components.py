#!/usr/bin/env python3
"""Component-wise PLQ analysis: decompose GridFF error into P, L, Q contributions.

For metals, we can't compare against pairwise LAMMPS (interaction isn't pairwise).
Instead we:
1. Show P, L, Q channel values at each z-point and site
2. Compare total GridFF against MACE teacher
3. Compute the "missing interaction" = MACE - GridFF(P+L+Q)
4. Analyze shape of each channel to identify which needs fixing
5. Test interpolation accuracy: B-spline cubic vs trilinear vs raw grid sampling
"""

from __future__ import annotations

import os
import sys
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
from pyBall.gridff_jax.hybrid_energy import HybridGridFFModel
from pyBall.gridff_jax.interfaces import PoseBatch, TeacherResult
from pyBall.gridff_jax.pose_sampling import (
    get_adsorbate, infer_reference_planes, infer_surface_sites, transform_adsorbate,
)
from pyBall.gridff_jax.utils import ensure_dir

OUT_DIR = ensure_dir(Path(__file__).resolve().parent / "component_analysis")

# === Load cached teacher data from the damping sweep ===
TEACHER_CACHE = Path(__file__).resolve().parent / "damping_sweep" / "teacher_cache.npz"


def _make_config(d0=3.70, width=0.35, interp_order=3):
    import os as _os
    config = RunConfig()
    chgcar = _os.environ.get("GRIDFF_JAX_CHGCAR")
    locpot = _os.environ.get("GRIDFF_JAX_LOCPOT")
    if not chgcar or not locpot:
        raise RuntimeError(
            "analyze_plq_components.py needs GRIDFF_JAX_CHGCAR and GRIDFF_JAX_LOCPOT "
            "env vars (point at your VASP slab outputs). Old hardcoded "
            "workstation paths have been removed."
        )
    config.density_backend.chgcar_path = chgcar
    config.density_backend.locpot_path = locpot
    config.grid.builder_mode = "metal_density_plq"
    config.grid.interpolation_order = interp_order
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


def _load_fit_params(path="tests/tGridFFJax/gridff_proof_run_damped_optimal/fit_params.json"):
    import json
    with open(ROOT / path) as f:
        return json.load(f)


def main():
    print("=" * 70)
    print("  COMPONENT-WISE PLQ ANALYSIS: CHONH2 / Ag(111)")
    print("=" * 70)
    print(f"Output: {OUT_DIR}")

    # Load teacher
    if TEACHER_CACHE.exists():
        tc = np.load(TEACHER_CACHE)
        teacher = TeacherResult(energies=tc["energies"], forces=tc["forces"], metadata={})
        print(f"Teacher cache loaded: {TEACHER_CACHE}")
    else:
        raise FileNotFoundError(f"Run sweep_london_damping.py first to generate {TEACHER_CACHE}")

    # Build config and density
    config = _make_config(d0=3.70, width=0.35, interp_order=3)
    print("Loading density (d0=3.70, w=0.35, bspline cubic)...", flush=True)
    density = make_density_backend(config.density_backend, surface=config.surface, grid=config.grid).load()

    adsorbate = get_adsorbate(name="CHONH2")
    site_uvs = {}
    for site in infer_surface_sites(density):
        if site.label in ["top", "bridge", "hollow"]:
            site_uvs[site.label] = np.asarray(site.uv)

    z_values = np.linspace(1.8, 8.0, 63)
    nz = len(z_values)
    poses = _build_poses(density, adsorbate, site_uvs, z_values)

    # Load optimal fit params
    fp = _load_fit_params()

    # Build model and reconstruct HybridParameters
    reactive_elements = tuple(sorted(set(adsorbate.symbols)))
    model = HybridGridFFModel(
        density=density, reactive_elements=reactive_elements,
        toggles=config.toggles, grid_config=config.grid,
        hybrid_config=config.hybrid_model, prefer_jax=True,
    )

    # Reconstruct params from the fit_params.json
    from pyBall.gridff_jax.hybrid_energy import HybridParameters
    params = HybridParameters(
        pauli=fp.get("pauli", {}),
        london=fp.get("london", {}),
        reactive=fp.get("reactive", {}),
        static_charge=fp.get("static_charge", {}),
        req_radius_offset=fp["req_radius_offset"],
        req_energy_scale=fp["req_energy_scale"],
        chi=fp.get("chi", {}),
        hardness=fp.get("hardness", {}),
        image_scale=fp.get("image_scale", 0.0),
        image_plane=fp.get("image_plane", 0.0),
        image_damping=fp.get("image_damping", 0.5),
        sample_shift_z=fp.get("sample_shift_z", 0.0),
        coulomb_shift_z=fp.get("coulomb_shift_z", 0.0),
        reservoir_chi=fp.get("reservoir_chi", 0.0),
        reservoir_hardness=fp.get("reservoir_hardness", 10.0),
        total_charge=fp.get("total_charge", 0.0),
    )

    # Evaluate with components
    print("Evaluating student (B-spline cubic, optimal damping)...", flush=True)
    student = model.evaluate_batch(poses, params=params, compute_forces=True)

    te = np.asarray(teacher.energies)
    se = np.asarray(student.energies)
    pauli_e = np.asarray(student.components["pauli"])
    london_e = np.asarray(student.components["london"])
    coulomb_e = np.asarray(student.components["coulomb"])

    # === Now test trilinear interpolation for comparison ===
    print("Building model with trilinear interpolation...", flush=True)
    config_lin = _make_config(d0=3.70, width=0.35, interp_order=1)
    density_lin = make_density_backend(config_lin.density_backend, surface=config_lin.surface, grid=config_lin.grid).load()
    model_lin = HybridGridFFModel(
        density=density_lin, reactive_elements=reactive_elements,
        toggles=config_lin.toggles, grid_config=config_lin.grid,
        hybrid_config=config_lin.hybrid_model, prefer_jax=True,
    )
    student_lin = model_lin.evaluate_batch(poses, params=params, compute_forces=True)
    se_lin = np.asarray(student_lin.energies)
    pauli_lin = np.asarray(student_lin.components["pauli"])
    london_lin = np.asarray(student_lin.components["london"])
    coulomb_lin = np.asarray(student_lin.components["coulomb"])

    # === Also test WITHOUT damping to see London field difference ===
    print("Building model WITHOUT London damping...", flush=True)
    config_nodamp = _make_config(d0=0.0, width=0.5, interp_order=3)
    density_nodamp = make_density_backend(config_nodamp.density_backend, surface=config_nodamp.surface, grid=config_nodamp.grid).load()
    model_nodamp = HybridGridFFModel(
        density=density_nodamp, reactive_elements=reactive_elements,
        toggles=config_nodamp.toggles, grid_config=config_nodamp.grid,
        hybrid_config=config_nodamp.hybrid_model, prefer_jax=True,
    )
    # Use same params for fair comparison of FIELD shapes
    student_nodamp = model_nodamp.evaluate_batch(poses, params=params, compute_forces=True)
    london_nodamp = np.asarray(student_nodamp.components["london"])

    # ============================================================
    # ANALYSIS AND PLOTTING
    # ============================================================
    site_names = list(site_uvs.keys())

    print("\n" + "=" * 70)
    print("  PART 1: COMPONENT-WISE ERROR BREAKDOWN")
    print("=" * 70)

    for si, site in enumerate(site_names):
        lo, hi = si * nz, (si + 1) * nz
        z = z_values
        mask = z >= 2.0

        te_s = te[lo:hi]
        se_s = se[lo:hi]
        p_s = pauli_e[lo:hi]
        l_s = london_e[lo:hi]
        q_s = coulomb_e[lo:hi]

        # The "missing interaction" = what MACE sees that PLQ doesn't
        missing = te_s - se_s  # MACE - GridFF
        plq_total = p_s + l_s + q_s

        print(f"\n  {site.upper()} site:")
        print(f"    At well minimum (z ≈ 2.7 A):")
        well_i = np.argmin(te_s[mask]) + np.argmax(mask)
        print(f"      MACE total:    {te_s[well_i]*1000:8.1f} meV")
        print(f"      GridFF total:  {se_s[well_i]*1000:8.1f} meV")
        print(f"      Pauli (P):     {p_s[well_i]*1000:8.1f} meV")
        print(f"      London (L):    {l_s[well_i]*1000:8.1f} meV")
        print(f"      Coulomb (Q):   {q_s[well_i]*1000:8.1f} meV")
        print(f"      Missing:       {missing[well_i]*1000:8.1f} meV")
        print(f"      P/(P+|L|):     {p_s[well_i]/(p_s[well_i]+abs(l_s[well_i]))*100:.1f}%")

        print(f"\n    Z-region breakdown:")
        for zlo, zhi in [(2.0, 2.5), (2.5, 3.0), (3.0, 4.0), (4.0, 6.0), (6.0, 8.0)]:
            m = (z >= zlo) & (z < zhi)
            if m.sum() == 0:
                continue
            p_rmse = np.sqrt(np.mean(p_s[m]**2)) * 1000
            l_rmse = np.sqrt(np.mean(l_s[m]**2)) * 1000
            q_rmse = np.sqrt(np.mean(q_s[m]**2)) * 1000
            total_rmse = np.sqrt(np.mean((se_s[m] - te_s[m])**2)) * 1000
            missing_mean = np.mean(missing[m]) * 1000
            print(f"      z=[{zlo:.0f},{zhi:.0f}): E_RMSE={total_rmse:6.1f}  Missing_mean={missing_mean:+7.1f}  "
                  f"|P|={p_rmse:7.1f}  |L|={l_rmse:7.1f}  |Q|={q_rmse:6.1f} meV")

    # PART 2: Interpolation accuracy
    print("\n" + "=" * 70)
    print("  PART 2: INTERPOLATION ACCURACY (B-spline cubic vs trilinear)")
    print("=" * 70)

    for si, site in enumerate(site_names):
        lo, hi = si * nz, (si + 1) * nz
        z = z_values
        mask = z >= 2.0

        e_diff = (se[lo:hi] - se_lin[lo:hi]) * 1000
        p_diff = (pauli_e[lo:hi] - pauli_lin[lo:hi]) * 1000
        l_diff = (london_e[lo:hi] - london_lin[lo:hi]) * 1000
        q_diff = (coulomb_e[lo:hi] - coulomb_lin[lo:hi]) * 1000

        print(f"\n  {site.upper()} (B-spline - trilinear):")
        print(f"    Total E:  RMSE={np.sqrt(np.mean(e_diff[mask]**2)):.3f} meV  max={np.max(np.abs(e_diff[mask])):.3f} meV")
        print(f"    Pauli:    RMSE={np.sqrt(np.mean(p_diff[mask]**2)):.3f} meV  max={np.max(np.abs(p_diff[mask])):.3f} meV")
        print(f"    London:   RMSE={np.sqrt(np.mean(l_diff[mask]**2)):.3f} meV  max={np.max(np.abs(l_diff[mask])):.3f} meV")
        print(f"    Coulomb:  RMSE={np.sqrt(np.mean(q_diff[mask]**2)):.3f} meV  max={np.max(np.abs(q_diff[mask])):.3f} meV")

    # PART 3: London damping effect on the field
    print("\n" + "=" * 70)
    print("  PART 3: LONDON FIELD WITH vs WITHOUT DAMPING")
    print("=" * 70)

    for si, site in enumerate(site_names):
        lo, hi = si * nz, (si + 1) * nz
        z = z_values
        l_damp = london_e[lo:hi] * 1000
        l_nodamp = london_nodamp[lo:hi] * 1000

        print(f"\n  {site.upper()}:")
        for zi in range(len(z)):
            if z[zi] >= 2.5 and z[zi] <= 6.0 and zi % 3 == 0:
                ratio = l_damp[zi] / l_nodamp[zi] if abs(l_nodamp[zi]) > 0.01 else float('nan')
                print(f"    z={z[zi]:.2f}: L_damped={l_damp[zi]:8.1f}  L_undamped={l_nodamp[zi]:8.1f}  ratio={ratio:.4f}")

    # ============================================================
    # PLOTS
    # ============================================================
    print("\nGenerating plots...", flush=True)

    # Plot 1: Component decomposition per site
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    for si, site in enumerate(site_names):
        lo, hi = si * nz, (si + 1) * nz
        z = z_values
        ax_top = axes[0, si]
        ax_bot = axes[1, si]

        ax_top.plot(z, te[lo:hi] * 1000, "k-", lw=2.5, label="MACE total")
        ax_top.plot(z, se[lo:hi] * 1000, "b-", lw=1.8, label="GridFF total")
        ax_top.plot(z, pauli_e[lo:hi] * 1000, "r:", lw=1.2, alpha=0.7, label="Pauli (P)")
        ax_top.plot(z, london_e[lo:hi] * 1000, "g:", lw=1.2, alpha=0.7, label="London (L)")
        ax_top.plot(z, coulomb_e[lo:hi] * 1000, "m:", lw=1.2, alpha=0.7, label="Coulomb (Q)")
        ax_top.set_xlim(2.0, 6.0)
        ax_top.set_ylim(-500, 600)
        ax_top.set_xlabel("z [A]")
        ax_top.set_ylabel("Energy [meV]")
        ax_top.set_title(f"{site.capitalize()}: PLQ Components")
        ax_top.legend(fontsize=7)
        ax_top.grid(True, alpha=0.3)
        ax_top.axhline(0, color="gray", lw=0.5)

        # Missing interaction
        missing = (te[lo:hi] - se[lo:hi]) * 1000
        ax_bot.fill_between(z, missing, alpha=0.4, color="orange", label="Missing (MACE - GridFF)")
        ax_bot.axhline(0, color="gray", lw=0.5)
        ax_bot.set_xlim(2.0, 6.0)
        ax_bot.set_xlabel("z [A]")
        ax_bot.set_ylabel("Missing interaction [meV]")
        ax_bot.set_title(f"{site.capitalize()}: Missing Physics")
        ax_bot.legend(fontsize=8)
        ax_bot.grid(True, alpha=0.3)

    plt.suptitle("CHONH2/Ag(111): PLQ Component Decomposition + Missing Interaction", fontsize=12)
    plt.tight_layout()
    plt.savefig(OUT_DIR / "plot1_component_decomposition.png", dpi=200)
    plt.close()

    # Plot 2: Interpolation accuracy
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    for si, site in enumerate(site_names):
        lo, hi = si * nz, (si + 1) * nz
        z = z_values
        ax = axes[si]
        ax.plot(z, (pauli_e[lo:hi] - pauli_lin[lo:hi]) * 1000, "r-", lw=1, label="P (bspline-trilinear)")
        ax.plot(z, (london_e[lo:hi] - london_lin[lo:hi]) * 1000, "g-", lw=1, label="L (bspline-trilinear)")
        ax.plot(z, (coulomb_e[lo:hi] - coulomb_lin[lo:hi]) * 1000, "m-", lw=1, label="Q (bspline-trilinear)")
        ax.plot(z, (se[lo:hi] - se_lin[lo:hi]) * 1000, "k-", lw=1.5, label="Total (bspline-trilinear)")
        ax.set_xlabel("z [A]")
        ax.set_ylabel("Interpolation diff [meV]")
        ax.set_title(f"{site.capitalize()}: B-spline vs Trilinear")
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)
        ax.axhline(0, color="gray", lw=0.5)
    plt.suptitle("Interpolation Accuracy: B-spline Cubic vs Trilinear", fontsize=12)
    plt.tight_layout()
    plt.savefig(OUT_DIR / "plot2_interpolation_accuracy.png", dpi=200)
    plt.close()

    # Plot 3: London field with vs without damping
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    for si, site in enumerate(site_names):
        lo, hi = si * nz, (si + 1) * nz
        z = z_values
        ax = axes[si]
        ax.plot(z, london_e[lo:hi] * 1000, "b-", lw=2, label="L with damping (d0=3.70)")
        ax.plot(z, london_nodamp[lo:hi] * 1000, "r--", lw=1.5, label="L without damping")
        ax.set_xlabel("z [A]")
        ax.set_ylabel("London energy [meV]")
        ax.set_title(f"{site.capitalize()}: London Damping Effect")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(2.0, 8.0)
    plt.suptitle("London Field: Effect of Fermi Damping", fontsize=12)
    plt.tight_layout()
    plt.savefig(OUT_DIR / "plot3_london_damping_effect.png", dpi=200)
    plt.close()

    # Plot 4: Missing interaction shape analysis
    fig, ax = plt.subplots(figsize=(10, 6))
    colors = {"top": "tab:blue", "bridge": "tab:green", "hollow": "tab:red"}
    for si, site in enumerate(site_names):
        lo, hi = si * nz, (si + 1) * nz
        z = z_values
        missing = (te[lo:hi] - se[lo:hi]) * 1000
        ax.plot(z, missing, "-", color=colors[site], lw=2, label=f"{site}: MACE - GridFF")

    # Overlay what -q²/4z would look like (image charge model)
    # CHONH2 has charges: N~-0.8, O~-0.5, H~+0.3, C~+0.2 → sum_q² ≈ 1.3 e²
    sum_q2 = 0.8**2 + 0.5**2 + 3*0.3**2 + 2*0.2**2  # approximate
    hartree_to_ev = 27.211
    bohr_to_ang = 0.529
    # E_image = -sum_q² / (4z) in atomic units → convert
    z_model = np.linspace(2.0, 8.0, 100)
    e_image_model = -sum_q2 * hartree_to_ev * bohr_to_ang / (4 * z_model) * 1000  # meV
    ax.plot(z_model, e_image_model, "k--", lw=1, alpha=0.5, label=f"Image charge model (-q²/4z, Σq²≈{sum_q2:.2f})")

    ax.set_xlabel("z above Ag(111) [A]")
    ax.set_ylabel("Missing interaction [meV]")
    ax.set_title("Missing Physics: MACE - GridFF(P+L+Q)")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color="gray", lw=0.5)
    ax.set_xlim(2.0, 8.0)
    plt.tight_layout()
    plt.savefig(OUT_DIR / "plot4_missing_interaction_analysis.png", dpi=200)
    plt.close()

    print(f"\nAll plots saved in: {OUT_DIR}")
    print("  plot1_component_decomposition.png — P, L, Q channels + missing interaction")
    print("  plot2_interpolation_accuracy.png — B-spline cubic vs trilinear difference")
    print("  plot3_london_damping_effect.png — London field with/without Fermi damping")
    print("  plot4_missing_interaction_analysis.png — Missing physics shape analysis")


if __name__ == "__main__":
    main()
