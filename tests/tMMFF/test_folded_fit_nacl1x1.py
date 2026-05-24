#!/usr/bin/env python3
import os, sys, json, argparse, time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from pyBall.OCL.MolecularDynamics import MolecularDynamics
from pyBall.OCL.InteractionEnergy import load_xyz_with_REQs

# CRITICAL: load_xyz_with_REQs returns REQ with E = sqrt(EvdW) for GridFF compatibility.
# This matches the GridFF generation convention in ocl_GridFF_new.py::make_atoms_arrays(bSqrtEvdw=True).
# When manually specifying probe REQs, ensure E is sqrt(EvdW), NOT raw EvdW.


def rmse(a, b):
    d = np.asarray(a, dtype=np.float64) - np.asarray(b, dtype=np.float64)
    return float(np.sqrt(np.mean(d*d)))


def rmse_masked(a, b, mask):
    a = np.asarray(a, dtype=np.float64); b = np.asarray(b, dtype=np.float64)
    m = np.asarray(mask, dtype=bool).reshape(-1)
    if len(a) != len(m):
        raise ValueError(f"rmse_masked(): len(a)={len(a)} len(mask)={len(m)}")
    if not np.any(m):
        return float('nan')
    d = a[m] - b[m]
    return float(np.sqrt(np.mean(d*d)))


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)


def pack_transforms(xyz):
    xyz = np.asarray(xyz, dtype=np.float32).reshape(-1, 3)
    T = np.zeros((len(xyz), 3, 4), dtype=np.float32)
    T[:, 0, 0] = 1.0; T[:, 1, 1] = 1.0; T[:, 2, 2] = 1.0
    T[:, :, 3] = xyz
    return T.reshape(-1, 12)


def eval_probe_component_reference(md_ref, xyz, component):
    T = pack_transforms(xyz)
    out = md_ref.eval_rigid_getSurfMorse_components(T, chunk_size=md_ref.nSystems, components=(component,))
    return np.asarray(out[component], dtype=np.float64)


def eval_probe_component_folded(md_fit, xyz, component):
    T = pack_transforms(xyz)
    out = md_fit.eval_rigid_getSurfFolded(T, chunk_size=md_fit.nSystems, component=component)
    return np.asarray(out['total'], dtype=np.float64)


def plot_1d_lines(md, xyz_1d, zs, tag, outdir, comps, line_label, line_x, line_y, ewald_info=None):
    """Plot 1D line scan: potential vs z for ref and folded at given (x,y)."""
    fig, axes = plt.subplots(3, 2, figsize=(14, 10))
    fig.suptitle(f'{tag} 1D Scan at ({line_x:.2f}, {line_y:.2f}) — {line_label}', fontsize=13)
    colors = {'pauli': '#e74c3c', 'london': '#3498db', 'coulomb': '#2ecc71', 'total': '#9b59b6'}
    for i, ck in enumerate(comps):
        ax = axes[i, 0]
        # Reference
        if ck == 'coulomb' and ewald_info is not None:
            ew, z_top, q_probe = ewald_info
            phi = ew.eval_full(xyz_1d[:, 0].reshape(-1,1), xyz_1d[:, 1].reshape(-1,1), (xyz_1d[:,2]-z_top).reshape(-1,1)).reshape(-1)
            y_ref = np.asarray(phi, dtype=np.float64) * float(q_probe)
        else:
            y_ref = eval_probe_component_reference(md, xyz_1d, ck)
        y_fit = eval_probe_component_folded(md, xyz_1d, ck)
        ax.plot(zs, y_ref, '-', color=colors[ck], lw=2, label=f'{ck} ref')
        ax.plot(zs, y_fit, '--', color=colors[ck], lw=1.5, label=f'{ck} folded')
        ax.set_ylabel('Energy [eV]')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        # Error panel
        ax2 = axes[i, 1]
        err = y_fit - y_ref
        ax2.plot(zs, err, '-', color='#e67e22', lw=1.5)
        ax2.axhline(0, color='k', ls='--', lw=0.5)
        ax2.set_ylabel('Error [eV]')
        ax2.set_xlabel('z - z_top [Å]')
        ax2.grid(True, alpha=0.3)
        rmse_val = float(np.sqrt(np.mean(err**2)))
        ax2.set_title(f'{ck} error (RMSE={rmse_val:.3e})', fontsize=9)
    # Total energy
    ax = axes[2, 0]
    y_tot_ref = np.zeros(len(xyz_1d))
    y_tot_fit = np.zeros(len(xyz_1d))
    for ck in comps:
        if ck == 'coulomb' and ewald_info is not None:
            ew, z_top, q_probe = ewald_info
            phi = ew.eval_full(xyz_1d[:,0].reshape(-1,1), xyz_1d[:,1].reshape(-1,1), (xyz_1d[:,2]-z_top).reshape(-1,1)).reshape(-1)
            y_tot_ref += np.asarray(phi, dtype=np.float64) * float(q_probe)
        else:
            y_tot_ref += eval_probe_component_reference(md, xyz_1d, ck)
        y_tot_fit += eval_probe_component_folded(md, xyz_1d, ck)
    ax.plot(zs, y_tot_ref, '-', color='#2c3e50', lw=2.5, label='Total ref')
    ax.plot(zs, y_tot_fit, '--', color='#e67e22', lw=2, label='Total folded')
    ax.set_xlabel('z - z_top [Å]')
    ax.set_ylabel('Energy [eV]')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax2 = axes[2, 1]
    err_tot = y_tot_fit - y_tot_ref
    ax2.plot(zs, err_tot, '-', color='#e67e22', lw=1.5)
    ax2.axhline(0, color='k', ls='--', lw=0.5)
    ax2.set_xlabel('z - z_top [Å]')
    ax2.set_ylabel('Error [eV]')
    ax2.grid(True, alpha=0.3)
    rmse_tot = float(np.sqrt(np.mean(err_tot**2)))
    ax2.set_title(f'Total error (RMSE={rmse_tot:.3e})', fontsize=9)
    fig.tight_layout()
    path = os.path.join(outdir, f'folded_fit_{tag}_{line_label.replace(" ","_")}_1d.png')
    fig.savefig(path, dpi=160, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved 1D line scan {path}')
    return rmse_tot


def plot_rmse_vs_z(md, xyz_3d, zs, nxy, nz, tag, outdir, comps, ewald_info=None):
    """Plot RMSE vs z height for each component."""
    fig, ax = plt.subplots(figsize=(10, 6))
    colors = {'pauli': '#e74c3c', 'london': '#3498db', 'coulomb': '#2ecc71'}
    # Extract z_top from ewald_info (used for x-axis labeling)
    _, z_top, _ = ewald_info if ewald_info is not None else (None, -1, None)
    for ck in comps:
        if ck == 'coulomb' and ewald_info is not None:
            ew, _, q_probe = ewald_info
            phi = ew.eval_full(xyz_3d[:,0].reshape(-1,1), xyz_3d[:,1].reshape(-1,1), (xyz_3d[:,2]-z_top).reshape(-1,1)).reshape(-1)
            y_ref = np.asarray(phi, dtype=np.float64) * float(q_probe)
        else:
            y_ref = eval_probe_component_reference(md, xyz_3d, ck)
        y_fit = eval_probe_component_folded(md, xyz_3d, ck)
        Eref = y_ref.reshape(nz, nxy, nxy)
        Efit = y_fit.reshape(nz, nxy, nxy)
        rmse_z = np.array([float(np.sqrt(np.mean((Efit[iz]-Eref[iz])**2))) for iz in range(nz)])
        ax.semilogy(zs - z_top, rmse_z, '-o', color=colors[ck], lw=1.5, ms=3, label=ck)
    ax.set_xlabel('z - z_top [Å]')
    ax.set_ylabel('RMSE [eV]')
    ax.set_title(f'{tag} RMSE vs height (each z-slice)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    path = os.path.join(outdir, f'folded_fit_{tag}_rmse_vs_z.png')
    fig.savefig(path, dpi=160)
    plt.close(fig)
    print(f'  Saved RMSE vs z plot {path}')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--surf', default=os.path.join('..', '..', 'cpp', 'common_resources', 'xyz', 'NaCl_1x1_L3.xyz'))
    parser.add_argument('--outdir', default=os.path.join('results_electrostatics', 'folded_fit_nacl1x1'))
    parser.add_argument('--nu', type=int, default=4)
    parser.add_argument('--nv', type=int, default=4)
    parser.add_argument('--nz', type=int, default=4)
    parser.add_argument('--nxy', type=int, default=32)
    parser.add_argument('--nz_samp', type=int, default=40)
    parser.add_argument('--z0', type=float, default=0.5)
    parser.add_argument('--z1', type=float, default=6.0)
    parser.add_argument('--nPBC', type=int, nargs=3, default=(1, 1, 0))
    parser.add_argument('--ewald_n_harm', type=int, default=6)
    parser.add_argument('--kernel', default='workgroup', choices=('orig', 'harmonics', 'workgroup'))
    parser.add_argument('--plot_nz', type=int, default=3)
    parser.add_argument('--z_fit_min', type=float, default=1.0)
    parser.add_argument('--repel_cut', type=float, default=1.0)
    parser.add_argument('--weight_power', type=float, default=6.0)
    args = parser.parse_args()

    ensure_dir(args.outdir)

    surf = os.path.normpath(os.path.join(os.path.dirname(__file__), args.surf))
    apos_s, REQs_s, enames_s, Zs_s, lvec_s = load_xyz_with_REQs(surf)
    z_top = float(np.max(apos_s[:, 2]))
    a = np.array(lvec_s[0, :3], dtype=np.float64)
    b = np.array(lvec_s[1, :3], dtype=np.float64)

    # Two probe types: +1 and -1 charges with dummy R,E
    probe_REQs = np.array([
        [1.5, 0.01, +1.0, 0.0],
        [1.5, 0.01, -1.0, 0.0],
    ], dtype=np.float32)
    probe_names = ['q+1', 'q-1']

    # Prepare sample grid in one cell
    us = np.linspace(0.0, 1.0, int(args.nxy), endpoint=False, dtype=np.float64)
    vs = np.linspace(0.0, 1.0, int(args.nxy), endpoint=False, dtype=np.float64)
    zs = np.linspace(z_top + float(args.z0), z_top + float(args.z1), int(args.nz_samp), endpoint=True, dtype=np.float64)

    xyz = np.empty((len(zs) * len(vs) * len(us), 3), dtype=np.float32)
    iu = 0
    for z in zs:
        for v in vs:
            for u in us:
                p = a*u + b*v
                xyz[iu, 0] = p[0]
                xyz[iu, 1] = p[1]
                xyz[iu, 2] = z
                iu += 1

    report = {
        'surf': surf,
        'z_top': z_top,
        'basis': {'nu': args.nu, 'nv': args.nv, 'nz': args.nz},
        'sample_grid': {'nxy': args.nxy, 'nz_samp': args.nz_samp, 'z_range_above_top': [args.z0, args.z1]},
        'kernel': args.kernel,
        'coulomb_solver': {'mode': 'ewald2d', 'ewald_n_harm': args.ewald_n_harm},
        'results': {},
    }

    for ip, (req, tag) in enumerate(zip(probe_REQs, probe_names)):
        print('='*70)
        print(f'Probe {tag} REQ={req.tolist()}')
        print('='*70)

        nSystems = min(len(xyz), 8192)
        md = MolecularDynamics(nloc=32, debug_build_options='-DDBG_UFF=0')
        md.init_rigid_molecule_batch(np.zeros((1, 3), dtype=np.float32), req[None, :], nSystems=nSystems)
        md.set_surface(surf, nPBC=tuple(args.nPBC), pos0=(0.0, 0.0, 0.0), alpha_morse=1.8, r_damp=0.0, bMacro=False)
        md.set_folded_kernel_kind(args.kernel)

        # Build fit mask based on physically relevant region (avoid huge repulsion)
        # Reference: Pauli+London from getSurfMorse, Coulomb from Ewald2D (same as in fitting)
        T_all = pack_transforms(xyz)
        ref_morse = md.eval_rigid_getSurfMorse_components(T_all, chunk_size=md.nSystems, components=('pauli', 'london'))
        ref_pauli = np.asarray(ref_morse['pauli'], dtype=np.float64)
        ref_london = np.asarray(ref_morse['london'], dtype=np.float64)
        # Coulomb ref via Ewald2D (quickly instantiate for mask calculation)
        from pyBall.OCL.SurfaceEwald import SurfaceEwaldCL
        ion_data = np.zeros((len(md.surface_atoms), 4), dtype=np.float32)
        ion_data[:, :3] = md.surface_atoms[:, :3]
        ion_data[:, 2] -= z_top
        ion_data[:, 3] = np.asarray(md.surface_REQs[:, 2], dtype=np.float32)
        ew = SurfaceEwaldCL(platform='nvidia')
        ew.prepare_system(ion_data, np.array(md.surface_lvec[0, :2], dtype=np.float32), np.array(md.surface_lvec[1, :2], dtype=np.float32), n_harm=int(args.ewald_n_harm))
        phi = ew.eval_full(xyz[:, 0].reshape(-1, 1), xyz[:, 1].reshape(-1, 1), (xyz[:, 2] - z_top).reshape(-1, 1)).reshape(-1)
        ref_coul = np.asarray(phi, dtype=np.float64) * float(req[2])
        ref_total = ref_pauli + ref_london + ref_coul
        z_rel = xyz[:, 2].astype(np.float64) - z_top
        fit_mask = (z_rel >= float(args.z_fit_min)) & (ref_total <= float(args.repel_cut))
        if not np.any(fit_mask):
            raise ValueError(f"Fit mask is empty: z_fit_min={args.z_fit_min} repel_cut={args.repel_cut} z_rel=[{z_rel.min():.3f},{z_rel.max():.3f}] ref_total=[{ref_total.min():.3f},{ref_total.max():.3f}]")
        print(f"[fit_mask] keep {int(np.count_nonzero(fit_mask))}/{len(fit_mask)} points | z_fit_min={args.z_fit_min} repel_cut={args.repel_cut}")

        t0 = time.perf_counter()
        md.fit_folded_surface_basis(surf, nPBC=tuple(args.nPBC), z_range=(z_top + args.z0, z_top + args.z1), nu=args.nu, nv=args.nv, nz=args.nz, nxy=args.nxy, nz_samp=args.nz_samp, r_damp=0.0, alpha_morse=1.8, bMacro=False, components=('pauli', 'london', 'coulomb'), fit_mask=fit_mask, weight_power=float(args.weight_power), coulomb_solver='ewald2d', ewald_n_harm=args.ewald_n_harm)
        t1 = time.perf_counter()
        print(f'[fit] time {t1-t0:.3f} s')

        # Evaluate reference and folded on same sample grid
        comps = ('pauli', 'london', 'coulomb')
        res = {}
        for ck in comps:
            # Reference
            if ck == 'coulomb':
                # Reference is already Ewald2D in fit, but here we compare folded-vs-reference stored in fit_info
                # We'll regenerate reference via fit_folded_surface_basis internals by calling it again is too heavy.
                # Instead: use current folded_fit_info refs.
                y_ref = np.asarray(md.folded_fit_info['refs'][ck][0], dtype=np.float64)
            else:
                y_ref = eval_probe_component_reference(md, xyz, ck)
            # Folded
            y_fit = eval_probe_component_folded(md, xyz, ck)
            err = y_fit - y_ref
            res[ck] = {
                'rmse_full': rmse(y_ref, y_fit),
                'rmse_fitmask': rmse_masked(y_ref, y_fit, fit_mask),
                'max_err_full': float(np.max(np.abs(err))),
                'max_err_fitmask': float(np.max(np.abs(err[fit_mask]))),
                'ref_min': float(np.min(y_ref)),
                'ref_max': float(np.max(y_ref)),
                'fit_min': float(np.min(y_fit)),
                'fit_max': float(np.max(y_fit)),
            }
            print(f'[{ck}] rmse(full)={res[ck]["rmse_full"]:.6e} max(full)={res[ck]["max_err_full"]:.6e}  rmse(mask)={res[ck]["rmse_fitmask"]:.6e} max(mask)={res[ck]["max_err_fitmask"]:.6e}')

        report['results'][tag] = res

        # Plot a few z-slices of errors for each component
        nxy = int(args.nxy)
        nz = int(args.nz_samp)
        nsel = min(int(args.plot_nz), nz)
        izs = np.linspace(0, nz-1, nsel, dtype=int)
        X = xyz[:, 0].reshape(nz, nxy, nxy)
        Y = xyz[:, 1].reshape(nz, nxy, nxy)

        for ck in comps:
            if ck == 'coulomb':
                y_ref = np.asarray(md.folded_fit_info['refs'][ck][0], dtype=np.float64)
            else:
                y_ref = eval_probe_component_reference(md, xyz, ck)
            y_fit = eval_probe_component_folded(md, xyz, ck)
            Eref = y_ref.reshape(nz, nxy, nxy)
            Efit = y_fit.reshape(nz, nxy, nxy)
            Err = (Efit - Eref)

            fig, axs = plt.subplots(nsel, 3, figsize=(12, 4*nsel))
            if nsel == 1:
                axs = axs.reshape(1, 3)
            for i, iz in enumerate(izs):
                vmin = float(min(Eref[iz].min(), Efit[iz].min()))
                vmax = float(max(Eref[iz].max(), Efit[iz].max()))
                im0 = axs[i, 0].pcolormesh(X[iz], Y[iz], Eref[iz], shading='auto', cmap='RdBu_r', vmin=vmin, vmax=vmax)
                axs[i, 0].set_title(f'{ck} ref z={float(zs[iz]-z_top):.2f}Å')
                plt.colorbar(im0, ax=axs[i, 0])
                im1 = axs[i, 1].pcolormesh(X[iz], Y[iz], Efit[iz], shading='auto', cmap='RdBu_r', vmin=vmin, vmax=vmax)
                axs[i, 1].set_title(f'{ck} folded')
                plt.colorbar(im1, ax=axs[i, 1])
                em = float(np.max(np.abs(Err[iz])))
                im2 = axs[i, 2].pcolormesh(X[iz], Y[iz], Err[iz], shading='auto', cmap='RdBu_r', vmin=-em, vmax=em)
                axs[i, 2].set_title(f'error (max {em:.2e})')
                plt.colorbar(im2, ax=axs[i, 2])
                for j in range(3):
                    axs[i, j].set_xlabel('x [Å]'); axs[i, j].set_ylabel('y [Å]')
            fig.tight_layout()
            path = os.path.join(args.outdir, f'folded_fit_{tag}_{ck}.png')
            fig.savefig(path, dpi=160)
            plt.close(fig)

        # === 1D Line Scans ===
        ewald_info = (ew, z_top, float(req[2]))
        nz_line = int(args.nz_samp)
        zs_line = np.linspace(z_top + float(args.z0), z_top + float(args.z1), nz_line, endpoint=True)
        # Line positions: above Na (0,0), above Cl (2,2), midpoint (1,1)
        line_positions = [
            (0.05, 0.05, 'on_Na'),
            (2.05, 2.05, 'on_Cl'),
            (1.0, 1.0, 'midpoint'),
        ]
        for lx, ly, llabel in line_positions:
            xyz_line = np.empty((nz_line, 3), dtype=np.float32)
            p = np.array([a[0]*lx + b[0]*ly, a[1]*lx + b[1]*ly, 0.0], dtype=np.float32)
            xyz_line[:, 0] = p[0]
            xyz_line[:, 1] = p[1]
            xyz_line[:, 2] = zs_line
            rmse_tot = plot_1d_lines(md, xyz_line, zs_line - z_top, tag, args.outdir, comps, llabel, lx, ly, ewald_info=ewald_info)
            res[f'1d_{llabel}_rmse_total'] = float(rmse_tot)

        # === RMSE vs z profile ===
        plot_rmse_vs_z(md, xyz, zs, nxy, nz, tag, args.outdir, comps, ewald_info=ewald_info)

    with open(os.path.join(args.outdir, 'report.json'), 'w') as f:
        json.dump(report, f, indent=2)

    print('\nSaved report to', os.path.join(args.outdir, 'report.json'))


if __name__ == '__main__':
    main()
