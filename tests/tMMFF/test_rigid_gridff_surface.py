#!/usr/bin/env python3
import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, BASE)

from pyBall.OCL.RigidBodyDynamics import RigidBodyDynamics


def write_xyz_frame(fout, enames, pos3, comment=""):
    nat = pos3.shape[0]
    fout.write(f"{nat}\n")
    fout.write(f"{comment}\n")
    for i in range(nat):
        fout.write(f"{enames[i]:3s} {pos3[i,0]:12.6f} {pos3[i,1]:12.6f} {pos3[i,2]:12.6f}\n")


def quat_from_euler(rx_deg, ry_deg, rz_deg):
    rx, ry, rz = np.radians([rx_deg, ry_deg, rz_deg]).astype(np.float32)
    cx, sx = np.cos(rx * 0.5), np.sin(rx * 0.5)
    cy, sy = np.cos(ry * 0.5), np.sin(ry * 0.5)
    cz, sz = np.cos(rz * 0.5), np.sin(rz * 0.5)
    qx = np.array([sx, 0.0, 0.0, cx], dtype=np.float32)
    qy = np.array([0.0, sy, 0.0, cy], dtype=np.float32)
    qz = np.array([0.0, 0.0, sz, cz], dtype=np.float32)
    def qmul(a, b):
        return np.array([
            a[3]*b[0] + a[0]*b[3] + a[1]*b[2] - a[2]*b[1],
            a[3]*b[1] - a[0]*b[2] + a[1]*b[3] + a[2]*b[0],
            a[3]*b[2] + a[0]*b[1] - a[1]*b[0] + a[2]*b[3],
            a[3]*b[3] - a[0]*b[0] - a[1]*b[1] - a[2]*b[2],
        ], dtype=np.float32)
    q = qmul(qmul(qz, qy), qx)
    q /= np.linalg.norm(q)
    return q


def quat_relative_angle_deg(q0, q1):
    q0 = np.asarray(q0, dtype=np.float64)
    q1 = np.asarray(q1, dtype=np.float64)
    q0 = q0 / np.linalg.norm(q0)
    q1 = q1 / np.linalg.norm(q1)
    dot = abs(float(np.dot(q0, q1)))
    dot = min(1.0, max(-1.0, dot))
    return 2.0 * np.degrees(np.arccos(dot))


def make_line_positions(n, p0, dp):
    ps = np.zeros((n, 3), dtype=np.float32)
    for i in range(n):
        ps[i] = np.asarray(p0, dtype=np.float32) + i * np.asarray(dp, dtype=np.float32)
    return ps


def resolve_default_paths(root, mol, sub):
    mol_path = mol if os.path.isabs(mol) else os.path.join(root, mol)
    sub_xyz = sub if os.path.isabs(sub) else os.path.join(root, sub)
    sub_name = os.path.splitext(os.path.basename(sub_xyz))[0]
    sub_grid = os.path.join(root, sub_name, 'Bspline_PLQd.npy')
    return mol_path, sub_xyz, sub_grid


def print_param_table(rbd):
    dbg = rbd.get_debug_dict()
    req = dbg['REQ']
    plq = dbg['PLQ']
    masses = dbg['masses']
    print("Assigned atom parameters:")
    print("i  type  mass      RvdW      sqrtE        Q         Hb        cP         cL         Q         cH")
    for i, en in enumerate(dbg['enames']):
        print(f"{i:1d}  {en:3s}  {masses[i]:7.4f}  {req[i,0]:8.5f}  {req[i,1]:8.5f}  {req[i,2]:9.5f}  {req[i,3]:8.5f}  {plq[i,0]:10.5f}  {plq[i,1]:10.5f}  {plq[i,2]:9.5f}  {plq[i,3]:10.5f}")
    print(f"Grid shape   : {dbg['grid_shape']}")
    print(f"Grid p0      : {dbg['grid_p0']}")
    print(f"Grid step    : {dbg['grid_step']}")


def save_zscan_plot(zs, Es, Fzs, labels, fout_png):
    fig, axs = plt.subplots(2, 1, figsize=(7, 8), sharex=True)
    for E, label in zip(Es, labels):
        axs[0].plot(zs, E, label=label)
    axs[0].set_ylabel('E [eV]')
    axs[0].grid(True)
    axs[0].legend()
    for Fz, label in zip(Fzs, labels):
        axs[1].plot(zs, Fz, label=label)
    axs[1].set_xlabel('COM z [A]')
    axs[1].set_ylabel('Fz [eV/A]')
    axs[1].grid(True)
    axs[1].legend()
    fig.tight_layout()
    fig.savefig(fout_png, dpi=150)
    plt.close(fig)


def save_traj_plot(hist, fout_png, enames):
    steps = np.asarray(hist['step'], dtype=np.int32)
    pos = np.asarray(hist['atom_pos'], dtype=np.float32)
    com = np.asarray(hist['com'], dtype=np.float32)
    fig, axs = plt.subplots(3, 1, figsize=(8, 10), sharex=True)
    for ia, en in enumerate(enames):
        axs[0].plot(steps, pos[:, ia, 2], label=f'{en}{ia}')
        axs[1].plot(steps, pos[:, ia, 0], label=f'{en}{ia}')
        axs[2].plot(steps, pos[:, ia, 1], label=f'{en}{ia}')
    axs[0].plot(steps, com[:, 2], 'k--', lw=1.0, label='COM')
    axs[1].plot(steps, com[:, 0], 'k--', lw=1.0, label='COM')
    axs[2].plot(steps, com[:, 1], 'k--', lw=1.0, label='COM')
    axs[0].set_ylabel('z [A]')
    axs[1].set_ylabel('x [A]')
    axs[2].set_ylabel('y [A]')
    axs[2].set_xlabel('step')
    for ax in axs:
        ax.grid(True)
        ax.legend(ncol=4, fontsize=8)
    fig.tight_layout()
    fig.savefig(fout_png, dpi=150)
    plt.close(fig)


def save_conv_plot(hist, fout_png):
    steps = np.asarray(hist['step'], dtype=np.int32)
    f = np.asarray(hist['f_norm'], dtype=np.float64)
    t = np.asarray(hist['t_norm'], dtype=np.float64)
    fig, axs = plt.subplots(2, 1, figsize=(7, 8), sharex=True)
    axs[0].semilogy(steps, np.maximum(f, 1e-16))
    axs[0].set_ylabel('|F_body| [eV/A]')
    axs[0].grid(True)
    axs[1].semilogy(steps, np.maximum(t, 1e-16))
    axs[1].set_ylabel('|T_body|')
    axs[1].set_xlabel('step')
    axs[1].grid(True)
    fig.tight_layout()
    fig.savefig(fout_png, dpi=150)
    plt.close(fig)


def do_zscan(args, rbd, out_dir):
    q0 = quat_from_euler(args.rx, args.ry, args.rz)
    zs = np.linspace(args.zscan_zmin, args.zscan_zmax, args.zscan_n, dtype=np.float32)
    sites = [
        ('Na_site', np.array([0.0, 0.0, 0.0], dtype=np.float32)),
        ('Cl_site', np.array([2.0, 2.0, 0.0], dtype=np.float32)),
    ]
    all_E = []
    all_Fz = []
    labels = []
    for label, site in sites:
        body_positions = np.zeros((1, 3), dtype=np.float32)
        body_positions[:, 0] = site[0]
        body_positions[:, 1] = site[1]
        quats = q0[None, :].copy()
        rbd_scan = RigidBodyDynamics.from_xyz_and_grid(args._mol_path, args._sub_grid, args._sub_xyz, n_bodies=1, body_positions=body_positions, quats=quats, alpha_morse=args.alpha, debug=args.debug_gpu)
        Es = np.zeros_like(zs)
        Fz = np.zeros_like(zs)
        for i, z in enumerate(zs):
            pos = np.zeros((1, 4), dtype=np.float32)
            pos[0, :3] = (site[0], site[1], z)
            pos[0, 3] = rbd_scan.download_outputs()['pos'][0, 3]
            outs = rbd_scan.download_outputs()
            cur_pos = outs['pos']
            cur_pos[0, :3] = (site[0], site[1], z)
            rbd_scan.toGPU('poss', cur_pos)
            rbd_scan.run_gridff(1, 0.0, lin_damp=1.0, ang_damp=1.0, force_scale=1.0, torque_scale=1.0)
            outs = rbd_scan.download_outputs()
            Es[i] = float(np.sum(outs['atom_force'][0, :, 3]))
            Fz[i] = float(outs['body_force'][0, 2])
        np.savetxt(os.path.join(out_dir, f'zscan_{label}.dat'), np.c_[zs, Es, Fz], header='z E Fz')
        all_E.append(Es)
        all_Fz.append(Fz)
        labels.append(label)
    save_zscan_plot(zs, all_E, all_Fz, labels, os.path.join(out_dir, 'zscan.png'))
    print(f'saved zscan    : {os.path.join(out_dir, "zscan.png")}')


def run_relax(args):
    data_root = os.path.join(BASE, 'cpp', 'common_resources')
    mol_path, sub_xyz, sub_grid = resolve_default_paths(data_root, args.mol, args.sub)
    if not os.path.exists(mol_path):
        raise FileNotFoundError(mol_path)
    if not os.path.exists(sub_xyz):
        raise FileNotFoundError(sub_xyz)
    if not os.path.exists(sub_grid):
        raise FileNotFoundError(sub_grid)
    if args.line_n < 1:
        raise ValueError(f'line_n must be >= 1, got {args.line_n}')
    args._mol_path = mol_path
    args._sub_xyz = sub_xyz
    args._sub_grid = sub_grid

    body_positions = make_line_positions(args.line_n, (args.x0, args.y0, args.z0), (args.dx, args.dy, args.dz))
    q0 = quat_from_euler(args.rx, args.ry, args.rz)
    quats = np.repeat(q0[None, :], args.line_n, axis=0)

    type_map = None
    if args.type_map:
        type_map = {}
        for item in args.type_map.split(','):
            item = item.strip()
            if not item:
                continue
            if ':' not in item:
                raise ValueError(f"Invalid type-map item '{item}', expected ELEM:TYPE")
            k, v = item.split(':', 1)
            type_map[k.strip()] = v.strip()
    rbd = RigidBodyDynamics.from_xyz_and_grid(mol_path, sub_grid, sub_xyz, n_bodies=args.line_n, body_positions=body_positions, quats=quats, alpha_morse=args.alpha, debug=args.debug_gpu, type_map=type_map)
    if args.init_from_npz:
        init_npz = np.load(args.init_from_npz)
        pos_init = np.array(init_npz['pos'], dtype=np.float32)
        quat_init = np.array(init_npz['quats'], dtype=np.float32)
        if pos_init.shape != (args.line_n, 4):
            raise ValueError(f"init_from_npz pos shape {pos_init.shape} != {(args.line_n,4)}")
        if quat_init.shape != (args.line_n, 4):
            raise ValueError(f"init_from_npz quats shape {quat_init.shape} != {(args.line_n,4)}")
        rbd.reset_pose(pos_init, quat_init)
        body_positions = pos_init[:, :3].copy()

    out_dir = os.path.abspath(args.outdir)
    os.makedirs(out_dir, exist_ok=True)
    traj_xyz = os.path.join(out_dir, args.traj)
    scan_xyz = os.path.join(out_dir, args.final_scan_xyz)
    final_npz = os.path.join(out_dir, args.final_npz)

    print(f'molecule     : {mol_path}')
    print(f'substrate    : {sub_xyz}')
    print(f'grid         : {sub_grid}')
    print(f'n_bodies     : {args.line_n}')
    print(f'initial pos0 : {(args.x0, args.y0, args.z0)}')
    print(f'line step    : {(args.dx, args.dy, args.dz)}')
    print(f'euler deg    : {(args.rx, args.ry, args.rz)}')
    print(f'dt           : {args.dt}')
    print(f'nsteps       : {args.nsteps}')
    print(f'dump_every   : {args.dump_every}')
    print_param_table(rbd)
    if args.do_zscan:
        do_zscan(args, rbd, out_dir)

    outputs = rbd.download_outputs()
    enames = rbd.enames
    nat = outputs['atom_positions'].shape[1]
    pos0 = outputs['atom_positions'][:, :, :3].copy()
    quat0 = outputs['quats'].copy()
    hist = {'step': [0], 'atom_pos': [outputs['atom_positions'][0, :, :3].copy()], 'com': [outputs['pos'][0, :3].copy()]}
    conv = {'step': [0], 'f_norm': [float(np.linalg.norm(outputs['body_force'][0, :3]))], 't_norm': [float(np.linalg.norm(outputs['body_torque'][0, :3]))]}
    reached_force_step = -1
    reached_torque_step = -1

    with open(traj_xyz, 'w') as fout:
        write_xyz_frame(fout, enames, outputs['atom_positions'][0, :, :3], comment='# step 0 body 0')
        for istep in range(1, args.nsteps + 1):
            rbd.run_gridff(1, args.dt, lin_damp=args.lin_damp, ang_damp=args.ang_damp, force_scale=args.force_scale, torque_scale=args.torque_scale)
            outputs = rbd.download_outputs()
            pos = outputs['atom_positions'][:, :, :3]
            com = outputs['pos'][:, :3]
            vel = outputs['lin_mom'][:, :3]
            omg = outputs['ang_mom'][:, :3]
            af = outputs['atom_force'][:, :, :3]
            aE = outputs['atom_force'][:, :, 3]
            bf = outputs['body_force'][:, :3]
            bt = outputs['body_torque'][:, :3]
            f_norm = float(np.sqrt(np.max(np.sum(bf * bf, axis=1))))
            t_norm = float(np.sqrt(np.max(np.sum(bt * bt, axis=1))))
            if not np.isfinite(pos).all() or not np.isfinite(com).all() or not np.isfinite(vel).all() or not np.isfinite(omg).all():
                raise RuntimeError(f'Non-finite state at step {istep}')
            vmax = float(np.sqrt(np.max(np.sum(vel * vel, axis=1)))) if len(vel) else 0.0
            wmax = float(np.sqrt(np.max(np.sum(omg * omg, axis=1)))) if len(omg) else 0.0
            zmin = float(pos[:, :, 2].min())
            zmax = float(pos[:, :, 2].max())
            disp = np.sqrt(np.sum((pos - pos0) ** 2, axis=2))
            dmin = float(disp.min())
            dmax = float(disp.max())
            if vmax > args.vmax_abort:
                raise RuntimeError(f'Linear blow-up at step {istep}: vmax={vmax}')
            if wmax > args.wmax_abort:
                raise RuntimeError(f'Angular blow-up at step {istep}: wmax={wmax}')
            if zmin < args.zmin_abort:
                raise RuntimeError(f'Molecule fell below zmin at step {istep}: zmin={zmin}')
            if zmax > args.zmax_abort:
                raise RuntimeError(f'Molecule escaped above zmax at step {istep}: zmax={zmax}')
            if reached_force_step < 0 and f_norm < args.fconv:
                reached_force_step = istep
            if reached_torque_step < 0 and t_norm < args.tconv:
                reached_torque_step = istep
            if (istep % args.dump_every) == 0 or istep == 1 or istep == args.nsteps:
                write_xyz_frame(fout, enames, pos[0], comment=f'# step {istep} body 0 com {com[0,0]:.6f} {com[0,1]:.6f} {com[0,2]:.6f} vmax {vmax:.6e} wmax {wmax:.6e}')
            if (istep % args.report_every) == 0 or istep == 1 or istep == args.nsteps:
                print(f'step {istep:6d}  com0=({com[0,0]:8.4f},{com[0,1]:8.4f},{com[0,2]:8.4f})  vmax={vmax:10.4e}  wmax={wmax:10.4e}  |F|={f_norm:10.4e}  |T|={t_norm:10.4e}  zspan=({zmin:8.4f},{zmax:8.4f})  disp[min,max]=({dmin:10.4e},{dmax:10.4e})')
                print(f'  body_force0 = ({bf[0,0]:11.4e},{bf[0,1]:11.4e},{bf[0,2]:11.4e})  body_torque0 = ({bt[0,0]:11.4e},{bt[0,1]:11.4e},{bt[0,2]:11.4e})')
                for ia, en in enumerate(enames):
                    print(f'  atom {ia:1d} {en:2s} pos=({pos[0,ia,0]:9.4f},{pos[0,ia,1]:9.4f},{pos[0,ia,2]:9.4f})  f=({af[0,ia,0]:11.4e},{af[0,ia,1]:11.4e},{af[0,ia,2]:11.4e})  E={aE[0,ia]:11.4e}')
            if (istep % args.plot_every) == 0 or istep == args.nsteps:
                hist['step'].append(istep)
                hist['atom_pos'].append(pos[0].copy())
                hist['com'].append(com[0].copy())
            if (istep % args.conv_every) == 0 or istep == 1 or istep == args.nsteps:
                conv['step'].append(istep)
                conv['f_norm'].append(f_norm)
                conv['t_norm'].append(t_norm)
            if args.stop_on_conv and (f_norm < args.fconv) and (t_norm < args.tconv):
                print(f'converged at step {istep} with |F|={f_norm:.6e} |T|={t_norm:.6e}')
                break

    with open(scan_xyz, 'w') as fout:
        for ib in range(args.line_n):
            write_xyz_frame(fout, enames, outputs['atom_positions'][ib, :, :3], comment=f'# final body {ib} com {outputs["pos"][ib,0]:.6f} {outputs["pos"][ib,1]:.6f} {outputs["pos"][ib,2]:.6f}')

    np.savez(
        final_npz,
        pos=outputs['pos'],
        quats=outputs['quats'],
        lin_mom=outputs['lin_mom'],
        ang_mom=outputs['ang_mom'],
        atom_positions=outputs['atom_positions'],
        mol_path=mol_path,
        sub_xyz=sub_xyz,
        sub_grid=sub_grid,
    )
    disp = np.sqrt(np.sum((outputs['atom_positions'][:, :, :3] - pos0) ** 2, axis=2))
    com_shift = outputs['pos'][:, :3] - body_positions
    rot_deg = np.array([quat_relative_angle_deg(quat0[i], outputs['quats'][i]) for i in range(args.line_n)], dtype=np.float64)
    save_traj_plot(hist, os.path.join(out_dir, 'trajectory.png'), enames)
    save_conv_plot(conv, os.path.join(out_dir, 'convergence.png'))
    print(f'saved trajectory : {traj_xyz}')
    print(f'saved final scan : {scan_xyz}')
    print(f'saved final npz  : {final_npz}')
    print(f'saved traj plot  : {os.path.join(out_dir, "trajectory.png")}')
    print(f'saved conv plot  : {os.path.join(out_dir, "convergence.png")}')
    print(f'final disp[min,max] : ({disp.min():.6f}, {disp.max():.6f}) A')
    print(f'final COM shift max : {np.sqrt(np.max(np.sum(com_shift * com_shift, axis=1))):.6f} A')
    print(f'final rotation max  : {rot_deg.max():.6f} deg')
    print(f'force threshold step: {reached_force_step}')
    print(f'torque threshold step: {reached_torque_step}')
    print(f'final |F| max       : {max(conv["f_norm"][-1], 0.0):.6e}')
    print(f'final |T| max       : {max(conv["t_norm"][-1], 0.0):.6e}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Rigid-body GridFF relaxation on substrate')
    parser.add_argument('--mol', default='xyz/H2O.xyz')
    parser.add_argument('--sub', default='xyz/NaCl_1x1_L3.xyz')
    parser.add_argument('--outdir', default='output_rigid_gridff')
    parser.add_argument('--traj', default='rigid_h2o_relax.xyz')
    parser.add_argument('--final-scan-xyz', dest='final_scan_xyz', default='rigid_scan_final.xyz')
    parser.add_argument('--final-npz', dest='final_npz', default='rigid_scan_final.npz')
    parser.add_argument('--nsteps', type=int, default=2400)
    parser.add_argument('--dt', type=float, default=0.06)
    parser.add_argument('--lin-damp', dest='lin_damp', type=float, default=0.994)
    parser.add_argument('--ang-damp', dest='ang_damp', type=float, default=0.988)
    parser.add_argument('--force-scale', dest='force_scale', type=float, default=1.0)
    parser.add_argument('--torque-scale', dest='torque_scale', type=float, default=1.0)
    parser.add_argument('--alpha', type=float, default=1.5)
    parser.add_argument('--x0', type=float, default=0.6)
    parser.add_argument('--y0', type=float, default=0.2)
    parser.add_argument('--z0', type=float, default=5.2)
    parser.add_argument('--rx', type=float, default=35.0)
    parser.add_argument('--ry', type=float, default=20.0)
    parser.add_argument('--rz', type=float, default=25.0)
    parser.add_argument('--line-n', dest='line_n', type=int, default=1)
    parser.add_argument('--dx', type=float, default=0.25)
    parser.add_argument('--dy', type=float, default=0.0)
    parser.add_argument('--dz', type=float, default=0.0)
    parser.add_argument('--dump-every', dest='dump_every', type=int, default=10)
    parser.add_argument('--plot-every', dest='plot_every', type=int, default=10)
    parser.add_argument('--report-every', dest='report_every', type=int, default=50)
    parser.add_argument('--conv-every', dest='conv_every', type=int, default=10)
    parser.add_argument('--vmax-abort', dest='vmax_abort', type=float, default=20.0)
    parser.add_argument('--wmax-abort', dest='wmax_abort', type=float, default=20.0)
    parser.add_argument('--zmin-abort', dest='zmin_abort', type=float, default=-1.0)
    parser.add_argument('--zmax-abort', dest='zmax_abort', type=float, default=15.0)
    parser.add_argument('--fconv', type=float, default=1e-3)
    parser.add_argument('--tconv', type=float, default=1e-3)
    parser.add_argument('--stop-on-conv', dest='stop_on_conv', action='store_true')
    parser.add_argument('--do-zscan', dest='do_zscan', action='store_true')
    parser.add_argument('--zscan-zmin', dest='zscan_zmin', type=float, default=0.5)
    parser.add_argument('--zscan-zmax', dest='zscan_zmax', type=float, default=6.0)
    parser.add_argument('--zscan-n', dest='zscan_n', type=int, default=80)
    parser.add_argument('--debug-gpu', dest='debug_gpu', action='store_true')
    parser.add_argument('--type-map', dest='type_map', default='O:O_3,H:H_OH')
    parser.add_argument('--init-from-npz', dest='init_from_npz', default='')
    args = parser.parse_args()
    run_relax(args)
