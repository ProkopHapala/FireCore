#!/usr/bin/env python3
"""Generate larger C diamond nanocrystals, relax, compute Hessian, run XRD + exact thermal broadening."""
import os, sys, time
import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, REPO)

from pyBall.XRD import (
    XRDDebye, load_xyz_atoms, get_form_factor_table,
    compute_sigma_from_sparse_blocks, compute_sigma_exact,
)

OUTDIR = REPO + '/OUT_XRD'
HESS_CACHE = OUTDIR + '/hessian_cache'
os.makedirs(HESS_CACHE, exist_ok=True)


def generate_and_relax(R):
    """Generate spherical C diamond NC of radius R, relax with MMFF, return atoms + hessian data."""
    from pyBall.nanocrystal_gen import build_spherical_nanoparticle
    from pyBall import MMFF, FTIR

    print(f"  Generating C diamond NC R={R} Å ...")
    t0 = time.time()
    prim_xyz = REPO + '/tests/tSiNCs/fixtures/vibration_parallel/structures/diamond_primitive.xyz'
    elems, apos = build_spherical_nanoparticle(prim_xyz=prim_xyz, R=R, nrep=8, heavy_z=6, cap='H')
    natoms = len(elems)
    pos = apos
    # Convert element strings to Z
    sym_to_Z = {'H': 1, 'C': 6, 'Si': 14}
    Z = np.array([sym_to_Z.get(e, 6) for e in elems], dtype=np.int32)
    print(f"    {natoms} atoms ({(Z==6).sum()} C + {(Z==1).sum()} H)  gen time={time.time()-t0:.2f}s")

    # Write XYZ for MMFF
    work_xyz = f'{HESS_CACHE}/diamond_R{R}.xyz'
    sym_map = {1: 'H', 6: 'C'}
    with open(work_xyz, 'w') as f:
        f.write(f"{natoms}\nR={R} diamond NC\n")
        for i in range(natoms):
            s = sym_map.get(int(Z[i]), 'C')
            f.write(f"{s} {pos[i,0]:.6f} {pos[i,1]:.6f} {pos[i,2]:.6f}\n")

    # Relax + compute Hessian with MMFF
    tmmff = REPO + '/tests/tMMFF'
    data_link = tmmff + '/data'
    if not os.path.exists(data_link):
        os.symlink('../../cpp/common_resources', data_link)

    orig = os.getcwd()
    try:
        os.chdir(tmmff)
        MMFF.setVerbosity(0, 0)
        print(f"  MMFF init + relax ...")
        t0 = time.time()
        MMFF.init(xyz_name=str(work_xyz), bEpairs=False, bMMFF=True, nPBC=(0, 0, 0))
        MMFF.getBuffs()
        natoms = MMFF.natoms
        MMFF.setSwitches(NonBonded=-1, MMFF=+1, Angles=+1, SurfAtoms=-1, GridFF=-1, PiSigma=-1, PiPiI=-1)
        MMFF.set_opt(dt_max=0.05, dt_min=0.01, damp_max=0.1, cvf_min=-0.1, cvf_max=+0.1)
        n_steps = MMFF.run(nstepMax=2000, dt=0.05, Fconv=1e-4, damping=0.1, ialg=2, omp=False)
        pos_relaxed = MMFF.apos[:natoms].copy()
        fmax = float(np.max(np.abs(MMFF.fapos[:natoms])))
        print(f"    relaxed in {n_steps} steps, fmax={fmax:.3e}, time={time.time()-t0:.2f}s")

        # Build neighbor shells + compute sparse Hessian
        print(f"  Computing sparse Hessian ...")
        t0 = time.time()
        neighs = MMFF.neighs
        neigh_idx, neigh_count = MMFF.buildNeighShells(neighs, n_shells=2, max_neigh=64, include_self=True)
        blocks = MMFF.getHessianSparseBlocks(neigh_idx, neigh_count, dx=1e-4)
        dt_hess = time.time() - t0
        print(f"    Hessian time={dt_hess:.2f}s  max_neigh={neigh_count.max()}")

        # Get mass matrix
        M = FTIR.get_mass_matrix(MMFF, natoms)
        masses = np.diag(M).reshape(natoms, 3)[:, 0]

        # Get Z from MMFF
        from pyBall import elements
        Z_mmff = np.zeros(natoms, dtype=np.int32)
        for i in range(natoms):
            tn = MMFF.getTypeName(i)
            sym = tn.split('_')[0] or 'C'
            Z_mmff[i] = int(elements.ELEMENT_DICT[sym][0]) if sym in elements.ELEMENT_DICT else 6

    finally:
        os.chdir(orig)

    # Cache
    cache_path = f'{HESS_CACHE}/diamond_R{R}_hessian.npz'
    np.savez(cache_path, pos=pos_relaxed, Z=Z_mmff, neigh_idx=neigh_idx, neigh_count=neigh_count,
             blocks=blocks, masses=masses)
    print(f"  Cached: {cache_path}")

    return pos_relaxed, Z_mmff, neigh_idx, neigh_count, blocks, masses


def run_xrd_exact(name, pos, Z, neigh_idx, neigh_count, blocks, masses):
    natoms = pos.shape[0]
    print(f"\n  XRD + thermal broadening for {name} ({natoms} atoms) ...")

    # Build atoms array for XRD
    unique_Z = np.unique(Z)
    elem_map = {z: i for i, z in enumerate(unique_Z)}
    elem_idx = np.array([elem_map[z] for z in Z], dtype=np.int32)
    atoms = np.column_stack([pos.astype(np.float32), elem_idx.astype(np.float32)]).astype(np.float32)

    engine = XRDDebye(preferred_vendor='nvidia')
    Q = np.linspace(0.5, 15.0, 800, dtype=np.float32)
    r_max = 20.0
    dr = 0.02

    # Static histogram
    t0 = time.time()
    hist_static, bin_edges = engine.compute_histogram(atoms, r_max=r_max, dr=dr)
    print(f"    histogram: {time.time()-t0:.2f}s  pairs={hist_static.sum():.0f}")

    ff_table = get_form_factor_table(unique_Z, Q)
    counts = np.bincount(elem_idx, minlength=len(unique_Z)).astype(np.float32)
    I_static = engine.compute_spectrum(hist_static, bin_edges, Q, ff_table, unique_Z)
    I_static += np.sum(counts[:, None] * (ff_table ** 2), axis=0)

    # Frozen sigma
    t0 = time.time()
    sigma_frozen = compute_sigma_from_sparse_blocks(
        pos, neigh_idx, neigh_count, blocks, kBT=0.02585, default_sigma=0.02)
    dt_frozen = time.time() - t0
    triu = np.triu_indices(natoms, k=1)
    sig_f = sigma_frozen[triu]
    print(f"    frozen sigma: {dt_frozen:.2f}s  range={sig_f.min():.5f}..{sig_f.max():.5f}  med={np.median(sig_f):.5f}")

    hist_frozen, _ = engine.compute_histogram_gaussian(atoms, sigma_frozen, r_max=r_max, dr=dr)
    I_frozen = engine.compute_spectrum(hist_frozen, bin_edges, Q, ff_table, unique_Z)
    I_frozen += np.sum(counts[:, None] * (ff_table ** 2), axis=0)

    # Exact sigma via splu
    t0 = time.time()
    sigma_exact = compute_sigma_exact(
        pos, neigh_idx, neigh_count, blocks, masses=masses,
        kBT=0.02585, default_sigma=0.02, rigid_proj=True)
    dt_exact = time.time() - t0
    sig_e = sigma_exact[triu]
    ratio = sig_e / sig_f
    print(f"    exact sigma: {dt_exact:.2f}s  range={sig_e.min():.5f}..{sig_e.max():.5f}  med={np.median(sig_e):.5f}")
    print(f"    ratio exact/frozen: med={np.median(ratio):.3f}  mean={np.mean(ratio):.3f}  max={np.max(ratio):.3f}")

    hist_exact, _ = engine.compute_histogram_gaussian(atoms, sigma_exact, r_max=r_max, dr=dr)
    I_exact = engine.compute_spectrum(hist_exact, bin_edges, Q, ff_table, unique_Z)
    I_exact += np.sum(counts[:, None] * (ff_table ** 2), axis=0)

    # Save
    npz_path = f'{OUTDIR}/{name}_xrd_exact.npz'
    np.savez(npz_path, Q=Q, I_static=I_static, I_frozen=I_frozen, I_exact=I_exact,
             sigma_frozen=sigma_frozen, sigma_exact=sigma_exact,
             hist_static=hist_static, hist_frozen=hist_frozen, hist_exact=hist_exact,
             bin_edges=bin_edges, unique_Z=unique_Z,
             timing_frozen_s=dt_frozen, timing_exact_s=dt_exact, natoms=natoms)
    print(f"    Saved: {npz_path}")

    # Plot
    try:
        import matplotlib; matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig = plt.figure(figsize=(14, 10))
        gs = fig.add_gridspec(3, 3)

        ax = fig.add_subplot(gs[0, 0])
        ax.hist(sig_f, bins=50, alpha=0.5, color='C0', label='frozen')
        ax.hist(sig_e, bins=50, alpha=0.5, color='C1', label='exact')
        ax.set_xlabel('σ (Å)'); ax.set_ylabel('count')
        ax.set_title(f'σ distribution ({natoms} atoms)'); ax.legend()

        ax = fig.add_subplot(gs[0, 1])
        ax.scatter(sig_f, sig_e, s=1, alpha=0.3)
        lim = [0, max(sig_f.max(), sig_e.max())]
        ax.plot(lim, lim, 'r--', lw=0.5)
        ax.set_xlabel('σ_frozen'); ax.set_ylabel('σ_exact')
        ax.set_title('Frozen vs exact σ')

        ax = fig.add_subplot(gs[0, 2])
        ax.hist(ratio, bins=50, alpha=0.7, color='C2')
        ax.set_xlabel('σ_exact / σ_frozen'); ax.set_ylabel('count')
        ax.set_title(f'Ratio (med={np.median(ratio):.2f})')

        ax = fig.add_subplot(gs[1, :])
        ax.plot(Q, I_static, '-', lw=0.8, color='C0', label='static')
        ax.plot(Q, I_frozen, '-', lw=0.8, color='C1', label='frozen σ')
        ax.plot(Q, I_exact, '-', lw=0.8, color='C2', label='exact σ (splu)')
        ax.set_xlabel('Q (Å⁻¹)'); ax.set_ylabel('I(Q)')
        ax.set_title(f'Powder XRD — {name} ({natoms} atoms)'); ax.legend()

        ax = fig.add_subplot(gs[2, 0])
        ax.plot(Q, I_frozen - I_static, '-', lw=0.8, color='C1', label='frozen − static')
        ax.plot(Q, I_exact - I_static, '-', lw=0.8, color='C2', label='exact − static')
        ax.set_xlabel('Q (Å⁻¹)'); ax.set_ylabel('ΔI(Q)'); ax.legend()

        ax = fig.add_subplot(gs[2, 1])
        ax.plot(Q, I_exact - I_frozen, '-', lw=0.8, color='C3')
        ax.set_xlabel('Q (Å⁻¹)'); ax.set_ylabel('ΔI(Q)'); ax.set_title('exact − frozen')

        ax = fig.add_subplot(gs[2, 2])
        ax.semilogy(Q, I_static, '-', lw=0.8, color='C0', label='static')
        ax.semilogy(Q, I_frozen, '-', lw=0.8, color='C1', label='frozen')
        ax.semilogy(Q, I_exact, '-', lw=0.8, color='C2', label='exact')
        ax.set_xlabel('Q (Å⁻¹)'); ax.set_ylabel('log I(Q)'); ax.legend()

        fig.tight_layout()
        png_path = f'{OUTDIR}/{name}_xrd_exact.png'
        fig.savefig(png_path, dpi=200)
        print(f"    Plot: {png_path}")
        plt.close(fig)
    except Exception as e:
        print(f"    Plot failed: {e}")

    return {'natoms': natoms, 'frozen_s': dt_frozen, 'exact_s': dt_exact,
            'sig_f_med': float(np.median(sig_f)), 'sig_e_med': float(np.median(sig_e)),
            'ratio_med': float(np.median(ratio)), 'ratio_max': float(np.max(ratio))}


if __name__ == '__main__':
    RADII = [8.0, 10.0]  # R=6 (270 atoms) already tested; these give ~500-1000 atoms
    results = {}

    for R in RADII:
        name = f'diamond_R{R}'
        cache = f'{HESS_CACHE}/diamond_R{R}_hessian.npz'
        print(f"\n{'='*60}")
        print(f"R={R} Å")
        print(f"{'='*60}")

        if os.path.exists(cache):
            print(f"  Loading cached: {cache}")
            d = np.load(cache)
            pos, Z = d['pos'], d['Z']
            neigh_idx, neigh_count, blocks = d['neigh_idx'], d['neigh_count'], d['blocks']
            masses = d['masses']
        else:
            pos, Z, neigh_idx, neigh_count, blocks, masses = generate_and_relax(R)

        r = run_xrd_exact(name, pos, Z, neigh_idx, neigh_count, blocks, masses)
        results[name] = r

    # Also run the existing R=6 fixture for comparison
    print(f"\n{'='*60}")
    print(f"R=6 Å (existing fixture)")
    print(f"{'='*60}")
    xyz_path = REPO + '/tests/tSiNCs/fixtures/vibration_parallel/structures/diamond_nc_R6_relaxed.xyz'
    hess_path = REPO + '/tests/tSiNCs/fixtures/vibration_parallel/hessian_mmff/diamond_nc_R6_sparse_blocks.npz'
    atoms, unique_Z = load_xyz_atoms(xyz_path)
    hess = np.load(hess_path)
    pos = atoms[:, :3].astype(np.float64)
    masses = np.full(atoms.shape[0], 12.011, dtype=np.float64)
    r = run_xrd_exact('diamond_R6_fixture', pos, atoms[:, 3].astype(np.int32),
                      hess['neigh_idx'], hess['neigh_count'], hess['blocks'], masses)
    results['diamond_R6_fixture'] = r

    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"{'Crystal':<25} {'N':>5} {'frozen(s)':>10} {'exact(s)':>10} {'σ_f_med':>10} {'σ_e_med':>10} {'ratio_med':>10} {'ratio_max':>10}")
    for name, r in results.items():
        print(f"{name:<25} {r['natoms']:>5} {r['frozen_s']:>10.2f} {r['exact_s']:>10.2f} {r['sig_f_med']:>10.5f} {r['sig_e_med']:>10.5f} {r['ratio_med']:>10.3f} {r['ratio_max']:>10.3f}")
