#!/usr/bin/env python3

import os, sys, argparse
import numpy as np

sys.path.append("../../")
from pyBall import MMFF, FTIR
from pyBall.nanocrystal_gen import build_diamond_nanoparticle, save_xyz


def reconstruct_dense_from_blocks(neigh_idx, neigh_count, blocks):
    natoms, max_neigh = neigh_idx.shape
    dim = natoms*3
    H = np.zeros((dim,dim), dtype=np.float64)
    for p in range(natoms):
        i0 = p*max_neigh
        for j in range(int(neigh_count[p])):
            o = int(neigh_idx[p,j])
            if o < 0: continue
            blk = blocks[p,j]
            H[o*3:(o+1)*3, p*3:(p+1)*3] = blk
    H = 0.5*(H + H.T)
    return H


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--R', type=float, default=10.0)
    ap.add_argument('--nrep', type=int, default=8)
    ap.add_argument('--dx', type=float, default=1e-4)
    ap.add_argument('--shells', type=int, default=2)
    ap.add_argument('--max-neigh', type=int, default=64)
    ap.add_argument('--outdir', type=str, default='OUT_nanoparticle')
    ap.add_argument('--plot', action='store_true')
    ap.add_argument('--diag', action='store_true')
    ap.add_argument('--spectrum', action='store_true', help='Compute vibration spectrum')
    ap.add_argument('--spectrum-method', type=str, default='dense', choices=['dense','sparse','both'],
                      help='Spectrum computation method: dense=eigenmode-based O(N^3)+O(nfreq*N), '
                           'sparse=sparse spsolve per freq O(nfreq*nnz), both=run both and compare')
    ap.add_argument('--fmin', type=float, default=0.01)
    ap.add_argument('--fmax', type=float, default=50.0)
    ap.add_argument('--nfreq', type=int, default=2000)
    ap.add_argument('--eta', type=float, default=0.05)
    ap.add_argument('--plot-spectrum', action='store_true')
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    elems, apos = build_diamond_nanoparticle(R=args.R, nrep=args.nrep)
    xyz0 = os.path.join(args.outdir, f"diamond_nc_R{args.R:.1f}_init.xyz")
    save_xyz(xyz0, elems, apos, comment=f"diamond nanoparticle R={args.R}")
    print(f"Saved init structure: {xyz0}  natoms={len(elems)}")

    MMFF.setVerbosity(0,0)
    MMFF.init(xyz_name=xyz0, bEpairs=False, bMMFF=True, nPBC=(0,0,0))
    MMFF.getBuffs()
    if MMFF.natoms != len(elems):
        raise ValueError(f"MMFF.init() natoms mismatch MMFF.natoms={MMFF.natoms} file={len(elems)}")
    MMFF.setSwitches(NonBonded=-1, MMFF=+1, Angles=+1, SurfAtoms=-1, GridFF=-1, PiSigma=-1, PiPiI=-1)

    # Relax
    MMFF.set_opt(dt_max=0.05, dt_min=0.01, damp_max=0.1, cvf_min=-0.1, cvf_max=+0.1)
    MMFF.run(nstepMax=2000, dt=0.05, Fconv=1e-4, damping=0.1, ialg=2, omp=False)

    xyz1 = os.path.join(args.outdir, f"diamond_nc_R{args.R:.1f}_relaxed.xyz")
    MMFF.saveXYZ(xyz1, "")
    print(f"Saved relaxed structure: {xyz1}")

    neigh1 = MMFF.neighs.copy()
    natoms = MMFF.natoms
    neigh_idx, neigh_count = MMFF.buildNeighShells(neigh1[:natoms], n_shells=args.shells, max_neigh=args.max_neigh, include_self=True)
    print(f"Neighbor shells built: shells={args.shells} max_neigh={neigh_idx.shape[1]} max_count={neigh_count.max()}")

    blocks = MMFF.getHessianSparseBlocks(neigh_idx, neigh_count, dx=args.dx)
    if np.isnan(blocks).any() or np.isinf(blocks).any():
        raise ValueError('NaN/Inf in sparse Hessian blocks')

    npz_path = os.path.join(args.outdir, f"diamond_nc_R{args.R:.1f}_sparse_hessian.npz")
    np.savez(npz_path, neigh_idx=neigh_idx, neigh_count=neigh_count, blocks=blocks)
    print(f"Saved sparse Hessian blocks: {npz_path}")

    H = reconstruct_dense_from_blocks(neigh_idx, neigh_count, blocks)
    if np.isnan(H).any() or np.isinf(H).any():
        raise ValueError('NaN/Inf in reconstructed dense Hessian')

    if args.diag:
        ndof = natoms * 3
        print(f"WARNING: diagonalizing dense {ndof}x{ndof} Hessian...")
        M = FTIR.get_mass_matrix(MMFF, natoms)
        H2 = FTIR.project_rigid_modes(H, M, MMFF.apos[:natoms].copy(), shift=1e6)
        w = np.linalg.eigvalsh(H2)
        w_sorted = np.sort(w)
        print(f"eigval min/max: {w_sorted[0]} {w_sorted[-1]}")

    if args.spectrum:
        ndof = natoms * 3
        M = FTIR.get_mass_matrix(MMFF, natoms)
        omegas = np.linspace(args.fmin, args.fmax, args.nfreq)
        charges = np.ones(natoms, dtype=np.float64)

        methods_to_run = []
        if args.spectrum_method in ('dense', 'both'):
            methods_to_run.append('dense')
        if args.spectrum_method in ('sparse', 'both'):
            methods_to_run.append('sparse')

        results = {}

        for method in methods_to_run:
            if method == 'dense':
                print(f"WARNING: computing vibration spectrum ({args.nfreq} freq points, {ndof} DOF) using DENSE eigenmode method...")
                H2 = FTIR.project_rigid_modes(H, M, MMFF.apos[:natoms].copy(), shift=1e6)
                res = FTIR.vibration_spectrum_from_modes(H2, M, omegas, eta=args.eta,
                                                          direction_vec=np.array([1.0,0.0,0.0]), charges=charges)
                results['dense'] = res
                spec_npz = os.path.join(args.outdir, f"diamond_nc_R{args.R:.1f}_spectrum_dense.npz")
                np.savez(spec_npz, omega=res['omega'], energy=res['energy'], dipole=res['dipole'])
                print(f"Saved dense spectrum: {spec_npz}")

            elif method == 'sparse':
                print(f"WARNING: computing vibration spectrum ({args.nfreq} freq points, {ndof} DOF) using SPARSE spsolve method...")
                # Build sparse Hessian directly from blocks (avoid dense reconstruction)
                H_sparse = FTIR.build_sparse_hessian_from_blocks(neigh_idx, neigh_count, blocks, symmetrize=True)
                res = FTIR.mechanical_greens_probing_sparse(H_sparse, M, omegas, eta=args.eta,
                                                             direction_vec=np.array([1.0,0.0,0.0]), charges=charges)
                results['sparse'] = res
                spec_npz = os.path.join(args.outdir, f"diamond_nc_R{args.R:.1f}_spectrum_sparse.npz")
                np.savez(spec_npz, omega=res['omega'], energy=res['energy'], dipole=res['dipole'])
                print(f"Saved sparse spectrum: {spec_npz}")

        if args.spectrum_method == 'both' and 'dense' in results and 'sparse' in results:
            # Compare results
            e_dense = results['dense']['energy']
            e_sparse = results['sparse']['energy']
            max_diff = np.max(np.abs(e_dense - e_sparse))
            rel_diff = max_diff / (np.max(np.abs(e_dense)) + 1e-12)
            print(f"SPECTRUM COMPARISON: max abs diff = {max_diff:.6e}, rel diff = {rel_diff:.6e}")

        # Plot first available result (or dense if both)
        if args.plot_spectrum:
            import matplotlib.pyplot as plt
            fig = plt.figure(figsize=(8,4))
            ax = fig.add_subplot(1,1,1)
            plot_method = 'dense' if 'dense' in results else list(results.keys())[0]
            res = results[plot_method]
            ax.plot(res['omega'], res['energy'], '-', lw=1.0, label=plot_method)
            if len(results) > 1:
                for method, res in results.items():
                    if method != plot_method:
                        ax.plot(res['omega'], res['energy'], '--', lw=1.0, label=method, alpha=0.7)
            ax.legend()
            ax.set_xlabel('omega')
            ax.set_ylabel('response energy')
            ax.set_title('Linear-response vibration spectrum')
            fig.tight_layout()
            png = os.path.join(args.outdir, f"diamond_nc_R{args.R:.1f}_spectrum.png")
            fig.savefig(png, dpi=200)
            print(f"Saved spectrum plot: {png}")

    if args.plot:
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(10,4))
        ax1 = fig.add_subplot(1,2,1)
        ax1.set_title('Dense Hessian (log|H|)')
        ax1.imshow(np.log10(np.abs(H) + 1e-12), cmap='viridis', origin='lower', interpolation='nearest')

        ax2 = fig.add_subplot(1,2,2)
        ax2.set_title('Sparse blocks rectangle')
        nat, max_neigh = neigh_idx.shape
        rect = np.zeros((nat*3, max_neigh*3), dtype=np.float64)
        for p in range(nat):
            for j in range(int(neigh_count[p])):
                rect[p*3:(p+1)*3, j*3:(j+1)*3] = blocks[p,j]
        ax2.imshow(np.log10(np.abs(rect) + 1e-12), cmap='viridis', origin='lower', interpolation='nearest', aspect='auto')

        fig.tight_layout()
        png_path = os.path.join(args.outdir, f"diamond_nc_R{args.R:.1f}_hessian_plots.png")
        fig.savefig(png_path, dpi=200)
        print(f"Saved plot: {png_path}")


if __name__ == "__main__":
    main()
