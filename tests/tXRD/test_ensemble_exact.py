#!/usr/bin/env python3
"""Test ensemble accumulation + exact thermal broadening (splu selected inversion).

Run from repo root:
    python tests/tXRD/test_ensemble_exact.py

Expects fixtures at tests/tSiNCs/fixtures/vibration_parallel/.
Outputs to OUT_XRD/.
"""
import os
import sys
import time
import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, REPO)

from pyBall.XRD import (
    XRDDebye, load_xyz_atoms, get_form_factor_table,
    compute_sigma_from_sparse_blocks, compute_sigma_exact,
)


def compute_spectrum_with_self(engine, atoms, unique_Z, Q, hist, bin_edges, n_atoms_override=None):
    """Debye transform + self-scattering term."""
    ff_table = get_form_factor_table(unique_Z, Q)
    I = engine.compute_spectrum(hist, bin_edges, Q, ff_table, unique_Z)
    n_species = len(unique_Z)
    if n_atoms_override is not None:
        # For ensemble: approximate self-term as N_total * <f^2> per species
        # For single-species (C only): exact
        if n_species == 1:
            I_self = n_atoms_override * ff_table[0] ** 2
        else:
            # Need per-species counts; skip for now
            I_self = np.zeros(len(Q), dtype=np.float32)
    else:
        elem_idx = atoms[:, 3].astype(np.int32)
        counts = np.bincount(elem_idx, minlength=n_species).astype(np.float32)
        I_self = np.sum(counts[:, None] * (ff_table ** 2), axis=0)
    return I + I_self


def test_ensemble_and_exact():
    xyz_path = os.path.join(REPO, 'tests/tSiNCs/fixtures/vibration_parallel/structures/diamond_nc_R6_relaxed.xyz')
    hess_path = os.path.join(REPO, 'tests/tSiNCs/fixtures/vibration_parallel/hessian_mmff/diamond_nc_R6_sparse_blocks.npz')
    if not os.path.exists(xyz_path):
        raise FileNotFoundError(f"Fixture not found: {xyz_path}")

    print(f"Loading {xyz_path} ...")
    atoms, unique_Z = load_xyz_atoms(xyz_path)
    n_atoms = atoms.shape[0]
    print(f"  natoms={n_atoms}  species={unique_Z}")

    engine = XRDDebye(preferred_vendor='nvidia')

    Q = np.linspace(0.5, 15.0, 800, dtype=np.float32)
    r_max = 20.0
    dr = 0.02

    # ── 1. Single-crystal reference (static) ──
    print(f"\n[1] Single crystal static histogram ...")
    hist_single, bin_edges = engine.compute_histogram(atoms, r_max=r_max, dr=dr)
    I_single = compute_spectrum_with_self(engine, atoms, unique_Z, Q, hist_single, bin_edges)
    print(f"    total_pairs={hist_single.sum():.1f}")

    # ── 2. Ensemble: N rotated copies of the same crystal ──
    # Simulate ensemble by rotating the crystal N_ensemble times
    # (realistic ensemble would have different sizes/shapes; this tests accumulation)
    N_ensemble = 8
    print(f"\n[2] Ensemble of {N_ensemble} randomly rotated copies ...")
    rng = np.random.default_rng(42)
    atoms_list = []
    for k in range(N_ensemble):
        # Random rotation matrix
        axis = rng.standard_normal(3)
        axis /= np.linalg.norm(axis)
        angle = rng.uniform(0, 2 * np.pi)
        c, s = np.cos(angle), np.sin(angle)
        R = np.array([[c + axis[0]**2*(1-c), axis[0]*axis[1]*(1-c) - axis[2]*s, axis[0]*axis[2]*(1-c) + axis[1]*s],
                       [axis[1]*axis[0]*(1-c) + axis[2]*s, c + axis[1]**2*(1-c), axis[1]*axis[2]*(1-c) - axis[0]*s],
                       [axis[2]*axis[0]*(1-c) - axis[1]*s, axis[2]*axis[1]*(1-c) + axis[0]*s, c + axis[2]**2*(1-c)]])
        pos_rot = atoms[:, :3] @ R.T
        atoms_rot = atoms.copy()
        atoms_rot[:, :3] = pos_rot.astype(np.float32)
        atoms_list.append(atoms_rot)

    t0 = time.time()
    hist_ens, bin_edges_ens, total_atoms = engine.ensemble_histogram(atoms_list, r_max=r_max, dr=dr)
    t1 = time.time()
    print(f"    ensemble hist time: {t1-t0:.2f}s  total_pairs={hist_ens.sum():.1f}  total_atoms={total_atoms}")
    I_ens = compute_spectrum_with_self(engine, atoms, unique_Z, Q, hist_ens, bin_edges_ens, n_atoms_override=total_atoms)

    # Verify: ensemble of identical (rotated) crystals should give same I(Q) as single / N_ensemble
    I_ens_norm = I_ens / N_ensemble
    max_diff = np.max(np.abs(I_ens_norm - I_single))
    print(f"    ensemble/N vs single max|diff| = {max_diff:.6e}")

    # ── 3. Frozen-stiffness sigma (existing method) ──
    I_frozen = None
    sigma_frozen = None
    if os.path.exists(hess_path):
        print(f"\n[3] Frozen-stiffness thermal broadening ...")
        hess = np.load(hess_path)
        pos = atoms[:, :3].astype(np.float64)
        t0 = time.time()
        sigma_frozen = compute_sigma_from_sparse_blocks(
            pos, hess['neigh_idx'], hess['neigh_count'], hess['blocks'],
            kBT=0.02585, default_sigma=0.02
        )
        t1 = time.time()
        triu = np.triu_indices(n_atoms, k=1)
        sig_vals = sigma_frozen[triu]
        print(f"    time: {t1-t0:.2f}s  sigma range: {sig_vals.min():.5f} .. {sig_vals.max():.5f}  median={np.median(sig_vals):.5f}")
        hist_frozen, _ = engine.compute_histogram_gaussian(atoms, sigma_frozen, r_max=r_max, dr=dr)
        I_frozen = compute_spectrum_with_self(engine, atoms, unique_Z, Q, hist_frozen, bin_edges)
    else:
        print(f"\n[3] Hessian fixture not found: {hess_path}")

    # ── 4. Exact sigma via splu selected inversion ──
    I_exact = None
    sigma_exact = None
    if os.path.exists(hess_path):
        print(f"\n[4] Exact thermal broadening via splu ...")
        t0 = time.time()
        # Use C mass for all atoms (fixture is carbon diamond NC)
        masses = np.full(n_atoms, 12.011, dtype=np.float64)
        sigma_exact = compute_sigma_exact(
            pos, hess['neigh_idx'], hess['neigh_count'], hess['blocks'],
            masses=masses, kBT=0.02585, default_sigma=0.02, rigid_proj=True
        )
        t1 = time.time()
        sig_vals_ex = sigma_exact[triu]
        print(f"    time: {t1-t0:.2f}s  sigma range: {sig_vals_ex.min():.5f} .. {sig_vals_ex.max():.5f}  median={np.median(sig_vals_ex):.5f}")
        hist_exact, _ = engine.compute_histogram_gaussian(atoms, sigma_exact, r_max=r_max, dr=dr)
        I_exact = compute_spectrum_with_self(engine, atoms, unique_Z, Q, hist_exact, bin_edges)

        # Compare frozen vs exact
        ratio = sigma_exact[triu] / sigma_frozen[triu]
        print(f"    exact/frozen ratio: median={np.median(ratio):.3f}  mean={np.mean(ratio):.3f}  max={np.max(ratio):.3f}")
    else:
        print(f"\n[4] Hessian fixture not found, skipping exact broadening")

    # ── 5. Ensemble + exact thermal broadening ──
    I_ens_exact = None
    if sigma_exact is not None:
        print(f"\n[5] Ensemble ({N_ensemble}) + exact thermal broadening ...")
        sigma_list = [sigma_exact] * N_ensemble  # same crystal, same sigma
        t0 = time.time()
        hist_ens_ex, _, total_atoms_ex = engine.ensemble_histogram_gaussian(atoms_list, sigma_list, r_max=r_max, dr=dr)
        t1 = time.time()
        print(f"    time: {t1-t0:.2f}s  total_pairs={hist_ens_ex.sum():.1f}")
        I_ens_exact = compute_spectrum_with_self(engine, atoms, unique_Z, Q, hist_ens_ex, bin_edges, n_atoms_override=total_atoms_ex)

    # ── Save & Plot ──
    outdir = os.path.join(REPO, 'OUT_XRD')
    os.makedirs(outdir, exist_ok=True)

    save_dict = dict(Q=Q, I_single=I_single, I_ens=I_ens, I_ens_norm=I_ens_norm,
                     hist_single=hist_single, hist_ens=hist_ens, bin_edges=bin_edges,
                     unique_Z=unique_Z, N_ensemble=N_ensemble)
    if I_frozen is not None:
        save_dict['I_frozen'] = I_frozen
    if I_exact is not None:
        save_dict['I_exact'] = I_exact
    if I_ens_exact is not None:
        save_dict['I_ens_exact'] = I_ens_exact
    if sigma_frozen is not None:
        save_dict['sigma_frozen'] = sigma_frozen
    if sigma_exact is not None:
        save_dict['sigma_exact'] = sigma_exact
    npz_path = os.path.join(outdir, 'ensemble_exact_xrd.npz')
    np.savez(npz_path, **save_dict)
    print(f"\nSaved NPZ: {npz_path}")

    try:
        import matplotlib.pyplot as plt

        fig = plt.figure(figsize=(14, 12))
        gs = fig.add_gridspec(4, 3)

        # Row 0: sigma comparison
        if sigma_frozen is not None and sigma_exact is not None:
            ax = fig.add_subplot(gs[0, 0])
            ax.hist(sig_vals, bins=50, alpha=0.5, color='C0', label='frozen')
            ax.hist(sig_vals_ex, bins=50, alpha=0.5, color='C1', label='exact')
            ax.set_xlabel('σ (Å)')
            ax.set_ylabel('count')
            ax.set_title('Per-pair σ distribution')
            ax.legend()

            ax = fig.add_subplot(gs[0, 1])
            ax.scatter(sig_vals, sig_vals_ex, s=1, alpha=0.3)
            lim = [0, max(sig_vals.max(), sig_vals_ex.max())]
            ax.plot(lim, lim, 'r--', lw=0.5)
            ax.set_xlabel('σ_frozen (Å)')
            ax.set_ylabel('σ_exact (Å)')
            ax.set_title('Frozen vs exact σ')

            ax = fig.add_subplot(gs[0, 2])
            ax.hist(ratio, bins=50, alpha=0.7, color='C2')
            ax.set_xlabel('σ_exact / σ_frozen')
            ax.set_ylabel('count')
            ax.set_title(f'Ratio (median={np.median(ratio):.2f})')

        # Row 1: I(Q) — ensemble vs single
        ax = fig.add_subplot(gs[1, :])
        ax.plot(Q, I_single, '-', lw=1.0, color='C0', label='single (static)')
        ax.plot(Q, I_ens_norm, '-', lw=1.0, color='C3', label=f'ensemble/{N_ensemble} (static)')
        ax.set_xlabel('Q (Å⁻¹)')
        ax.set_ylabel('I(Q)')
        ax.set_title(f'Ensemble accumulation check — {N_ensemble} rotated copies')
        ax.legend()

        # Row 2: I(Q) — thermal broadening comparison
        ax = fig.add_subplot(gs[2, :])
        ax.plot(Q, I_single, '-', lw=0.8, color='C0', label='static')
        if I_frozen is not None:
            ax.plot(Q, I_frozen, '-', lw=0.8, color='C1', label='frozen σ')
        if I_exact is not None:
            ax.plot(Q, I_exact, '-', lw=0.8, color='C2', label='exact σ (splu)')
        if I_ens_exact is not None:
            ax.plot(Q, I_ens_exact / N_ensemble, '-', lw=0.8, color='C4', label=f'ensemble/{N_ensemble} + exact σ')
        ax.set_xlabel('Q (Å⁻¹)')
        ax.set_ylabel('I(Q)')
        ax.set_title('Thermal broadening: static vs frozen vs exact')
        ax.legend()

        # Row 3: difference and log
        ax = fig.add_subplot(gs[3, 0])
        if I_frozen is not None:
            ax.plot(Q, I_frozen - I_single, '-', lw=0.8, color='C1', label='frozen − static')
        if I_exact is not None:
            ax.plot(Q, I_exact - I_single, '-', lw=0.8, color='C2', label='exact − static')
        ax.set_xlabel('Q (Å⁻¹)')
        ax.set_ylabel('ΔI(Q)')
        ax.set_title('Broadening difference')
        ax.legend()

        ax = fig.add_subplot(gs[3, 1])
        if I_exact is not None and I_frozen is not None:
            ax.plot(Q, I_exact - I_frozen, '-', lw=0.8, color='C3')
            ax.set_title('exact − frozen')
        ax.set_xlabel('Q (Å⁻¹)')
        ax.set_ylabel('ΔI(Q)')

        ax = fig.add_subplot(gs[3, 2])
        ax.semilogy(Q, I_single, '-', lw=0.8, color='C0', label='static')
        if I_frozen is not None:
            ax.semilogy(Q, I_frozen, '-', lw=0.8, color='C1', label='frozen')
        if I_exact is not None:
            ax.semilogy(Q, I_exact, '-', lw=0.8, color='C2', label='exact')
        ax.set_xlabel('Q (Å⁻¹)')
        ax.set_ylabel('log I(Q)')
        ax.set_title('Log scale')
        ax.legend()

        fig.tight_layout()
        png_path = os.path.join(outdir, 'ensemble_exact_xrd.png')
        fig.savefig(png_path, dpi=200)
        print(f"Saved plot: {png_path}")
        plt.close(fig)
    except Exception as e:
        print(f"Plotting skipped: {e}")

    return I_single, Q


if __name__ == '__main__':
    test_ensemble_and_exact()
