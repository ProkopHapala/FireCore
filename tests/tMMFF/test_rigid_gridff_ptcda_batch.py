#!/usr/bin/env python3
import os, sys, time, argparse
import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(BASE)

from pyBall.OCL.RigidBodyDynamics import RigidBodyDynamics
from pyBall.OCL.InteractionEnergy import load_xyz_with_REQs

np.set_printoptions(linewidth=200, suppress=True)


def write_xyz_frame(fout, enames, pos3, comment=""):
    fout.write(f"{len(enames)}\n")
    fout.write(f"{comment}\n")
    for e, p in zip(enames, pos3):
        fout.write(f"{e:3s} {p[0]:12.6f} {p[1]:12.6f} {p[2]:12.6f}\n")


def save_xyz_stack(fname, enames, atom_positions, comments):
    with open(fname, 'w') as fout:
        for i in range(len(atom_positions)):
            write_xyz_frame(fout, enames, atom_positions[i, :, :3], comment=comments[i])


def quat_mul(q1, q2):
    q1 = np.asarray(q1, dtype=np.float32)
    q2 = np.asarray(q2, dtype=np.float32)
    x1, y1, z1, w1 = q1[:, 0], q1[:, 1], q1[:, 2], q1[:, 3]
    x2, y2, z2, w2 = q2[:, 0], q2[:, 1], q2[:, 2], q2[:, 3]
    out = np.empty_like(q1)
    out[:, 0] = w1*x2 + x1*w2 + y1*z2 - z1*y2
    out[:, 1] = w1*y2 - x1*z2 + y1*w2 + z1*x2
    out[:, 2] = w1*z2 + x1*y2 - y1*x2 + z1*w2
    out[:, 3] = w1*w2 - x1*x2 - y1*y2 - z1*z2
    return out


def quat_from_axis_angle(axis, ang):
    ang = np.asarray(ang, dtype=np.float32).reshape(-1)
    axis = np.asarray(axis, dtype=np.float32)
    axis = axis / np.linalg.norm(axis)
    h = 0.5 * ang
    out = np.zeros((len(ang), 4), dtype=np.float32)
    out[:, :3] = axis[None, :] * np.sin(h)[:, None]
    out[:, 3] = np.cos(h)
    return out


def make_rotation_set(nrot, tilt_max, ntilt):
    phis = np.linspace(0.0, 2.0*np.pi, nrot, endpoint=False, dtype=np.float32)
    qz = quat_from_axis_angle((0.0, 0.0, 1.0), phis)
    if ntilt <= 1 or tilt_max <= 0.0:
        return qz
    tilts = np.linspace(-tilt_max, tilt_max, ntilt, dtype=np.float32)
    qx = quat_from_axis_angle((1.0, 0.0, 0.0), tilts)
    qy = quat_from_axis_angle((0.0, 1.0, 0.0), tilts)
    qs = [qz]
    for qa in qx:
        qs.append(quat_mul(qz, np.repeat(qa[None, :], len(qz), axis=0)))
    for qa in qy:
        qs.append(quat_mul(qz, np.repeat(qa[None, :], len(qz), axis=0)))
    q = np.concatenate(qs, axis=0)
    q /= np.linalg.norm(q, axis=1)[:, None]
    return q.astype(np.float32)


def make_translation_grid(lvec, nx, ny, z0, margin):
    fx = (np.arange(nx, dtype=np.float32) + 0.5) / float(nx)
    fy = (np.arange(ny, dtype=np.float32) + 0.5) / float(ny)
    gx, gy = np.meshgrid(fx, fy, indexing='xy')
    frac = np.stack([gx.reshape(-1), gy.reshape(-1)], axis=1)
    if margin > 0.0:
        frac = frac * (1.0 - 2.0*margin) + margin
    xy = frac[:, 0, None] * lvec[0][None, :3] + frac[:, 1, None] * lvec[1][None, :3]
    pos = np.zeros((len(frac), 3), dtype=np.float32)
    pos[:, :2] = xy[:, :2]
    pos[:, 2] = z0
    return pos, frac


def tile_transforms(pos_grid, frac_grid, quats):
    nt = len(pos_grid)
    nr = len(quats)
    pos = np.repeat(pos_grid, nr, axis=0).astype(np.float32)
    frac = np.repeat(frac_grid, nr, axis=0).astype(np.float32)
    q = np.tile(quats, (nt, 1)).astype(np.float32)
    rot_id = np.tile(np.arange(nr, dtype=np.int32), nt)
    return pos, frac, q, rot_id


def wrap_frac_delta(df):
    return df - np.round(df)


def quat_dist(q1, q2):
    d1 = np.linalg.norm(q1 - q2)
    d2 = np.linalg.norm(q1 + q2)
    return min(float(d1), float(d2))


def frac_to_cart_xy(df, lvec):
    return df[0] * lvec[0, :2] + df[1] * lvec[1, :2]


def cart_to_frac_xy(pos_xy, lvec):
    return np.linalg.solve(lvec[:2, :2].T, np.asarray(pos_xy, dtype=np.float32).T).T.astype(np.float32)


def cluster_minima(pos, frac, quats, f_norm, t_norm, e_tot, pos_eps, quat_eps, lvec):
    order = np.argsort(e_tot)
    reps = []
    labels = -np.ones(len(order), dtype=np.int32)
    for ii, idx in enumerate(order):
        assigned = -1
        for ic, rep in enumerate(reps):
            dfrac = wrap_frac_delta(frac[idx] - frac[rep])
            dxy = float(np.linalg.norm(frac_to_cart_xy(dfrac, lvec)))
            dz = float(abs(pos[idx, 2] - pos[rep, 2]))
            dpos = float(np.sqrt(dxy*dxy + dz*dz))
            dquat = quat_dist(quats[idx], quats[rep])
            if dpos < pos_eps and dquat < quat_eps:
                assigned = ic
                break
        if assigned < 0:
            assigned = len(reps)
            reps.append(idx)
        labels[ii] = assigned
    counts = np.zeros(len(reps), dtype=np.int32)
    for lab in labels:
        counts[lab] += 1
    clusters = []
    for ic, rep in enumerate(reps):
        clusters.append({
            'rep': rep,
            'count': int(counts[ic]),
            'pos': pos[rep].copy(),
            'quat': quats[rep].copy(),
            'f_norm': float(f_norm[rep]),
            't_norm': float(t_norm[rep]),
            'energy': float(e_tot[rep]),
        })
    clusters.sort(key=lambda c: c['energy'])
    return clusters


def save_cluster_xyz(fname, enames, atom_positions, clusters, nmax):
    with open(fname, 'w') as fout:
        for i, c in enumerate(clusters[:nmax]):
            ib = c['rep']
            write_xyz_frame(fout, enames, atom_positions[ib, :, :3], comment=f"# cluster {i} count {c['count']} E {c['energy']:.6e} F {c['f_norm']:.6e} T {c['t_norm']:.6e}")


def sampled_stability_report(atom_positions, prev_atom_positions=None, sample_stride=1, zmin_limit=None, jump_limit=None):
    idx = np.arange(0, len(atom_positions), max(1, int(sample_stride)), dtype=np.int32)
    apos = atom_positions[idx, :, :3]
    zmin = apos[:, :, 2].min(axis=1)
    out = {
        'sample_indices': idx,
        'zmin': zmin,
        'bad_z': np.where((zmin_limit is not None) & (zmin < zmin_limit))[0] if zmin_limit is not None else np.zeros(0, dtype=np.int32),
    }
    if prev_atom_positions is not None:
        d = apos - prev_atom_positions[idx, :, :3]
        jump = np.sqrt(np.sum(d * d, axis=2)).max(axis=1)
        out['jump'] = jump
        out['bad_jump'] = np.where((jump_limit is not None) & (jump > jump_limit))[0] if jump_limit is not None else np.zeros(0, dtype=np.int32)
    return out


def run_batch(args):
    data_root = os.path.join(BASE, 'cpp', 'common_resources')
    mol_path = os.path.join(data_root, args.mol)
    sub_xyz = os.path.join(data_root, args.sub)
    sub_grid = os.path.join(data_root, args.grid)
    if not os.path.exists(mol_path): raise FileNotFoundError(mol_path)
    if not os.path.exists(sub_xyz): raise FileNotFoundError(sub_xyz)
    if not os.path.exists(sub_grid): raise FileNotFoundError(sub_grid)

    _, _, _, _, lvec = load_xyz_with_REQs(sub_xyz)
    if lvec is None:
        raise RuntimeError(f"Substrate {sub_xyz} has no lattice vectors in XYZ comment")

    if args.init_from_npz:
        init = np.load(args.init_from_npz)
        pos_init = np.array(init['pos'], dtype=np.float32)
        quats = np.array(init['quats'], dtype=np.float32)
        n_bodies = len(pos_init)
        body_positions = pos_init[:, :3].copy()
        if 'frac0' in init:
            body_frac = np.array(init['frac0'], dtype=np.float32)
        else:
            frac_xy = np.linalg.solve(lvec[:2, :2].T, body_positions[:, :2].T).T.astype(np.float32)
            body_frac = frac_xy
        if 'rot_id' in init:
            rot_ids = np.array(init['rot_id'], dtype=np.int32)
        else:
            rot_ids = np.full(n_bodies, -1, dtype=np.int32)
    else:
        pos_grid, frac_grid = make_translation_grid(lvec, args.nx, args.ny, args.z0, args.margin)
        qrots = make_rotation_set(args.nrot, np.deg2rad(args.tilt_deg), args.ntilt)
        body_positions, body_frac, quats, rot_ids = tile_transforms(pos_grid, frac_grid, qrots)
        n_bodies = len(body_positions)
        if args.max_bodies > 0 and n_bodies > args.max_bodies:
            body_positions = body_positions[:args.max_bodies]
            body_frac = body_frac[:args.max_bodies]
            quats = quats[:args.max_bodies]
            rot_ids = rot_ids[:args.max_bodies]
            n_bodies = len(body_positions)

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

    out_dir = os.path.abspath(args.outdir)
    os.makedirs(out_dir, exist_ok=True)

    rbd = RigidBodyDynamics.from_xyz_and_grid(
        mol_path, sub_grid, sub_xyz,
        n_bodies=n_bodies,
        body_positions=body_positions,
        quats=quats,
        alpha_morse=args.alpha,
        debug=args.debug_gpu,
        type_map=type_map,
        mass_trans=args.mass_trans,
        mass_rot=args.mass_rot,
    )
    init_outputs = rbd.download_selected(('atom_positions',))

    t0 = time.time()
    print(f"molecule   : {mol_path}")
    print(f"substrate  : {sub_xyz}")
    print(f"grid       : {sub_grid}")
    print(f"n_grid     : {args.nx} x {args.ny}")
    print(f"n_rot      : {args.nrot if not args.init_from_npz else 'restart'}")
    print(f"n_bodies   : {n_bodies}")
    print(f"dt         : {args.dt}")
    print(f"nsteps     : {args.nsteps}")
    print(f"chunk      : {args.chunk_steps}")

    reached_force = np.full(n_bodies, -1, dtype=np.int32)
    reached_torque = np.full(n_bodies, -1, dtype=np.int32)
    unstable = np.zeros(n_bodies, dtype=np.bool_)
    n_done_prev = 0
    outputs = None
    prev_sample_atoms = init_outputs['atom_positions'].copy() if args.sample_check_stride > 0 else None

    for istep0 in range(0, args.nsteps, args.chunk_steps):
        nrun = min(args.chunk_steps, args.nsteps - istep0)
        rbd.run_gridff(nrun, args.dt, lin_damp=args.lin_damp, ang_damp=args.ang_damp, force_scale=args.force_scale, torque_scale=args.torque_scale)
        fields = ['pos', 'quats', 'body_force', 'body_torque']
        do_sample = (args.sample_check_stride > 0) and (((istep0 + nrun) % args.sample_every_steps) == 0 or ((istep0 + nrun) == args.nsteps))
        if do_sample:
            fields.append('atom_positions')
        outputs = rbd.download_selected(tuple(fields))
        bf = outputs['body_force'][:, :3] / args.force_scale
        bt = outputs['body_torque'][:, :3] / args.torque_scale
        f_norm = np.sqrt(np.sum(bf * bf, axis=1))
        t_norm = np.sqrt(np.sum(bt * bt, axis=1))
        mask_f = (reached_force < 0) & (f_norm < args.fconv)
        mask_t = (reached_torque < 0) & (t_norm < args.tconv)
        reached_force[mask_f] = istep0 + nrun
        reached_torque[mask_t] = istep0 + nrun
        converged = (f_norm < args.fconv) & (t_norm < args.tconv)
        n_done = int(np.count_nonzero(converged))
        if do_sample:
            rep = sampled_stability_report(outputs['atom_positions'], prev_atom_positions=prev_sample_atoms, sample_stride=args.sample_check_stride, zmin_limit=args.zmin_limit, jump_limit=args.jump_limit)
            if 'bad_jump' in rep and len(rep['bad_jump']) > 0:
                unstable[rep['sample_indices'][rep['bad_jump']]] = True
            if len(rep['bad_z']) > 0:
                unstable[rep['sample_indices'][rep['bad_z']]] = True
            prev_sample_atoms = outputs['atom_positions'].copy()
            print(f"sample-check step {istep0+nrun:6d}: zmin[min,max]=({rep['zmin'].min():.4f},{rep['zmin'].max():.4f}) bad_z={len(rep['bad_z'])}" + (f" bad_jump={len(rep['bad_jump'])} jump[max]={rep['jump'].max():.4f}" if 'jump' in rep else ""))
        if ((istep0 + nrun) % args.report_every_steps == 0) or ((istep0 + nrun) == args.nsteps):
            print(f"step {istep0+nrun:6d}/{args.nsteps} converged {n_done:6d}/{n_bodies} |F|max {f_norm.max():10.4e} |T|max {t_norm.max():10.4e} wall {time.time()-t0:8.2f}s")
        if args.stop_when_all and n_done == n_bodies and n_done > 0 and n_done != n_done_prev:
            print(f"all replicas converged by step {istep0+nrun}")
            break
        n_done_prev = n_done

    outputs = rbd.download_outputs()
    bf = outputs['body_force'][:, :3] / args.force_scale
    bt = outputs['body_torque'][:, :3] / args.torque_scale
    f_norm = np.sqrt(np.sum(bf * bf, axis=1))
    t_norm = np.sqrt(np.sum(bt * bt, axis=1))
    e_tot = np.sum(outputs['atom_force'][:, :, 3], axis=1)
    converged = (f_norm < args.fconv) & (t_norm < args.tconv)
    iconv = np.where(converged)[0]

    np.savez(
        os.path.join(out_dir, args.final_npz),
        pos=outputs['pos'],
        quats=outputs['quats'],
        atom_positions=outputs['atom_positions'],
        body_force=outputs['body_force'],
        body_torque=outputs['body_torque'],
        body_force_phys=bf,
        body_torque_phys=bt,
        f_norm=f_norm,
        t_norm=t_norm,
        energy=e_tot,
        converged=converged,
        unstable=unstable,
        reached_force=reached_force,
        reached_torque=reached_torque,
        frac0=body_frac,
        rot_id=rot_ids,
        mol_path=mol_path,
        sub_xyz=sub_xyz,
        sub_grid=sub_grid,
    )

    with open(os.path.join(out_dir, 'converged_indices.txt'), 'w') as fout:
        for i in iconv:
            fout.write(f"{i} {f_norm[i]:.8e} {t_norm[i]:.8e} {e_tot[i]:.8e} {body_frac[i,0]:.8f} {body_frac[i,1]:.8f} {rot_ids[i]}\n")
    with open(os.path.join(out_dir, 'unstable_indices.txt'), 'w') as fout:
        for i in np.where(unstable)[0]:
            fout.write(f"{i} {f_norm[i]:.8e} {t_norm[i]:.8e} {e_tot[i]:.8e} {body_frac[i,0]:.8f} {body_frac[i,1]:.8f} {rot_ids[i]}\n")

    frac_final = cart_to_frac_xy(outputs['pos'][:, :2], lvec)
    clusters = cluster_minima(outputs['pos'][iconv], frac_final[iconv], outputs['quats'][iconv], f_norm[iconv], t_norm[iconv], e_tot[iconv], args.cluster_pos_eps, args.cluster_quat_eps, lvec) if len(iconv) else []
    with open(os.path.join(out_dir, 'cluster_summary.txt'), 'w') as fout:
        fout.write(f"# n_bodies {n_bodies}\n")
        fout.write(f"# n_converged {len(iconv)}\n")
        fout.write(f"# n_clusters {len(clusters)}\n")
        fout.write(f"# fconv {args.fconv:.8e} tconv {args.tconv:.8e} pos_eps {args.cluster_pos_eps:.8e} quat_eps {args.cluster_quat_eps:.8e}\n")
        for i, c in enumerate(clusters):
            p = c['pos']
            q = c['quat']
            fout.write(f"{i:6d} rep {iconv[c['rep']]:6d} count {c['count']:6d} E {c['energy']: .8e} F {c['f_norm']: .8e} T {c['t_norm']: .8e} pos {p[0]: .6f} {p[1]: .6f} {p[2]: .6f} quat {q[0]: .6f} {q[1]: .6f} {q[2]: .6f} {q[3]: .6f}\n")

    if len(clusters):
        save_cluster_xyz(os.path.join(out_dir, 'cluster_representatives.xyz'), rbd.enames, outputs['atom_positions'][iconv], clusters, args.export_clusters)
    init_comments = [f"# replica {i} init rot_id {int(rot_ids[i])} frac {body_frac[i,0]:.8f} {body_frac[i,1]:.8f}" for i in range(n_bodies)]
    final_comments = [f"# replica {i} final E {e_tot[i]:.8e} F {f_norm[i]:.8e} T {t_norm[i]:.8e} converged {int(converged[i])} unstable {int(unstable[i])} rot_id {int(rot_ids[i])} frac0 {body_frac[i,0]:.8f} {body_frac[i,1]:.8f} fracf {frac_final[i,0]:.8f} {frac_final[i,1]:.8f}" for i in range(n_bodies)]
    save_xyz_stack(os.path.join(out_dir, 'initial_all.xyz'), rbd.enames, init_outputs['atom_positions'], init_comments)
    save_xyz_stack(os.path.join(out_dir, 'final_all.xyz'), rbd.enames, outputs['atom_positions'], final_comments)
    if len(iconv):
        conv_comments = [f"# converged replica {int(iconv[j])} E {e_tot[iconv[j]]:.8e} F {f_norm[iconv[j]]:.8e} T {t_norm[iconv[j]]:.8e} rot_id {int(rot_ids[iconv[j]])} frac0 {body_frac[iconv[j],0]:.8f} {body_frac[iconv[j],1]:.8f} fracf {frac_final[iconv[j],0]:.8f} {frac_final[iconv[j],1]:.8f}" for j in range(len(iconv))]
        save_xyz_stack(os.path.join(out_dir, 'converged_all.xyz'), rbd.enames, outputs['atom_positions'][iconv], conv_comments)

    print(f"saved npz      : {os.path.join(out_dir, args.final_npz)}")
    print(f"saved conv idx : {os.path.join(out_dir, 'converged_indices.txt')}")
    print(f"saved unstab   : {os.path.join(out_dir, 'unstable_indices.txt')}")
    print(f"saved clusters : {os.path.join(out_dir, 'cluster_summary.txt')}")
    print(f"saved xyz init : {os.path.join(out_dir, 'initial_all.xyz')}")
    print(f"saved xyz final: {os.path.join(out_dir, 'final_all.xyz')}")
    if len(iconv):
        print(f"saved xyz conv : {os.path.join(out_dir, 'converged_all.xyz')}")
    if len(clusters):
        print(f"saved reps xyz : {os.path.join(out_dir, 'cluster_representatives.xyz')}")
    print(f"converged      : {len(iconv)} / {n_bodies}")
    print(f"distinct minima: {len(clusters)}")


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Production multi-replica rigid GridFF relaxation for PTCDA on NaCl')
    ap.add_argument('--mol', default='xyz/PTCDA.xyz')
    ap.add_argument('--sub', default='xyz/NaCl_1x1_L3.xyz')
    ap.add_argument('--grid', default='NaCl_1x1_L3/Bspline_PLQd.npy')
    ap.add_argument('--outdir', default='output_rigid_gridff_ptcda_batch')
    ap.add_argument('--final-npz', dest='final_npz', default='rigid_ptcda_batch_final.npz')
    ap.add_argument('--nx', type=int, default=20)
    ap.add_argument('--ny', type=int, default=20)
    ap.add_argument('--nrot', type=int, default=12)
    ap.add_argument('--ntilt', type=int, default=1)
    ap.add_argument('--tilt-deg', dest='tilt_deg', type=float, default=0.0)
    ap.add_argument('--z0', type=float, default=3.2)
    ap.add_argument('--margin', type=float, default=0.0)
    ap.add_argument('--max-bodies', dest='max_bodies', type=int, default=5000)
    ap.add_argument('--nsteps', type=int, default=2000)
    ap.add_argument('--chunk-steps', dest='chunk_steps', type=int, default=100)
    ap.add_argument('--report-every-steps', dest='report_every_steps', type=int, default=500)
    ap.add_argument('--dt', type=float, default=0.05)
    ap.add_argument('--lin-damp', dest='lin_damp', type=float, default=0.99)
    ap.add_argument('--ang-damp', dest='ang_damp', type=float, default=0.99)
    ap.add_argument('--mass-trans', dest='mass_trans', type=float, default=1.0)
    ap.add_argument('--mass-rot', dest='mass_rot', type=float, default=4.0)
    ap.add_argument('--force-scale', dest='force_scale', type=float, default=1.0)
    ap.add_argument('--torque-scale', dest='torque_scale', type=float, default=1.0)
    ap.add_argument('--alpha', type=float, default=1.5)
    ap.add_argument('--fconv', type=float, default=1e-3)
    ap.add_argument('--tconv', type=float, default=1e-3)
    ap.add_argument('--cluster-pos-eps', dest='cluster_pos_eps', type=float, default=0.1)
    ap.add_argument('--cluster-quat-eps', dest='cluster_quat_eps', type=float, default=0.02)
    ap.add_argument('--export-clusters', dest='export_clusters', type=int, default=64)
    ap.add_argument('--type-map', dest='type_map', default='')
    ap.add_argument('--init-from-npz', dest='init_from_npz', default='')
    ap.add_argument('--sample-check-stride', dest='sample_check_stride', type=int, default=64)
    ap.add_argument('--sample-every-steps', dest='sample_every_steps', type=int, default=500)
    ap.add_argument('--zmin-limit', dest='zmin_limit', type=float, default=2.0)
    ap.add_argument('--jump-limit', dest='jump_limit', type=float, default=1.0)
    ap.add_argument('--stop-when-all', dest='stop_when_all', action='store_true')
    ap.add_argument('--debug-gpu', dest='debug_gpu', action='store_true')
    args = ap.parse_args()
    run_batch(args)
