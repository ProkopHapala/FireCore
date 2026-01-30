import numpy as np
import argparse
import os
import sys

import matplotlib.pyplot as plt
from numpy.random import default_rng

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from XPDB_new import XPDB_new


def load_xyz(fname):
    with open(fname, 'r') as f:
        lines_raw = [l.rstrip('\n') for l in f.readlines()]

    # XYZ: line0=nAtoms, line1=comment (can be empty), then n atom lines.
    # Some files here use "-lvs ..." as comment line.
    if len(lines_raw) < 2:
        raise ValueError(f"load_xyz: file too short '{fname}'")

    n = int(lines_raw[0].strip().split()[0])
    if n <= 0:
        raise ValueError(f"load_xyz: invalid nAtoms={n} in '{fname}'")

    comment = lines_raw[1].strip()
    i0 = 2
    if comment.startswith('-lvs'):
        i0 = 2

    # Skip any extra empty lines before atoms (some xyz have blank comment line + extra blank)
    while i0 < len(lines_raw) and (lines_raw[i0].strip() == ''):
        i0 += 1
    if i0 + n > len(lines_raw):
        raise ValueError(f"load_xyz: not enough atom lines in '{fname}' (need {n}, have {len(lines_raw)-i0})")
    elems = []
    xyz = np.zeros((n, 3), dtype=np.float32)
    q = np.zeros((n,), dtype=np.float32)
    for i in range(n):
        ws = lines_raw[i0 + i].split()
        elems.append(ws[0])
        xyz[i, 0] = float(ws[1]); xyz[i, 1] = float(ws[2]); xyz[i, 2] = float(ws[3])
        if len(ws) > 4:
            q[i] = float(ws[4])
    return elems, xyz, q


def masses_from_elems(elems):
    m = np.zeros((len(elems),), dtype=np.float32)
    for i, e in enumerate(elems):
        if e == 'H':
            m[i] = 1.0
        elif e == 'C':
            m[i] = 12.0
        elif e == 'N':
            m[i] = 14.0
        elif e == 'O':
            m[i] = 16.0
        else:
            raise ValueError(f"Unknown element '{e}'")
    return m


def build_neighs_bk_from_bonds(n, bonds, max_deg=4):
    neighs = np.full((n, max_deg), -1, dtype=np.int32)
    bks = np.full((n, max_deg), -1, dtype=np.int32)

    deg = np.zeros((n,), dtype=np.int32)

    for (i, j) in bonds:
        if deg[i] >= max_deg or deg[j] >= max_deg:
            raise RuntimeError(f"Degree>={max_deg} for bond {i}-{j}")
        si = int(deg[i]); sj = int(deg[j])
        neighs[i, si] = j
        neighs[j, sj] = i
        bks[i, si] = sj
        bks[j, sj] = si
        deg[i] += 1
        deg[j] += 1

    return neighs, bks


def bonds_for_molecule(name, elems):
    # Debug test only: hardcoded bonds for the specific xyz files
    if name.lower().startswith('h2o'):
        # O(0)-H(1), O(0)-H(2)
        return [(0, 1), (0, 2)], 0  # nnode=0 (no rigid ports for water)
    if name.lower().startswith('ch2nh'):
        # C(0)=N(1), C(0)-H(2), C(0)-H(3), N(1)-H(4)
        return [(0, 1), (0, 2), (0, 3), (1, 4)], 2  # nnode=2 (C,N)
    raise ValueError(f"No bond template for '{name}'")


def make_bLs_bKs_from_neighs(xyz, neighs, *, k_bond=200.0):
    n = xyz.shape[0]
    bLs = np.zeros((n, 4), dtype=np.float32)
    bKs = np.zeros((n, 4), dtype=np.float32)
    for i in range(n):
        for k in range(4):
            j = int(neighs[i, k])
            if j < 0:
                bLs[i, k] = 0.0
                bKs[i, k] = 0.0
                continue
            bLs[i, k] = float(np.linalg.norm(xyz[j] - xyz[i]))
            bKs[i, k] = float(k_bond)
    return bLs, bKs


def com(pos, m):
    M = float(np.sum(m))
    return np.sum(pos * m[:, None], axis=0) / M


def linear_momentum(vel, m):
    return np.sum(vel * m[:, None], axis=0)


def angular_momentum(pos, vel, m, omega, Iiso):
    c = com(pos, m)
    r = pos - c
    L_orb = np.sum(np.cross(r, vel * m[:, None]), axis=0)
    L_spin = np.sum(omega * Iiso[:, None], axis=0)
    return L_orb + L_spin


def quat_mul(a, b):
    out = np.empty_like(a)
    ax, ay, az, aw = a[:, 0], a[:, 1], a[:, 2], a[:, 3]
    bx, by, bz, bw = b[:, 0], b[:, 1], b[:, 2], b[:, 3]
    out[:, 0] = aw*bx + ax*bw + ay*bz - az*by
    out[:, 1] = aw*by - ax*bz + ay*bw + az*bx
    out[:, 2] = aw*bz + ax*by - ay*bx + az*bw
    out[:, 3] = aw*bw - ax*bx - ay*by - az*bz
    return out


def normalize_quat(q):
    norms = np.linalg.norm(q, axis=1)
    norms[norms == 0.0] = 1.0
    return q / norms[:, None]


def perturb_state(pos, quat, pos_scale, rot_scale, rng):
    pos_out = pos.copy()
    quat_out = quat.copy()
    if pos_scale > 0:
        pos_out += rng.normal(scale=pos_scale, size=pos_out.shape)
    if rot_scale > 0:
        axes = rng.normal(size=(len(quat_out), 3))
        axis_norm = np.linalg.norm(axes, axis=1) + 1e-12
        axes /= axis_norm[:, None]
        angles = rng.normal(scale=rot_scale, size=(len(quat_out),))
        half = 0.5 * angles
        sin_half = np.sin(half)
        dq = np.zeros_like(quat_out)
        dq[:, :3] = axes * sin_half[:, None]
        dq[:, 3] = np.cos(half)
        quat_out = normalize_quat(quat_mul(dq, quat_out))
    return pos_out, quat_out


def write_xyz_with_ports(fname, elems, pos, pneigh, port_n):
    n = pos.shape[0]
    nports = int(np.sum(port_n))
    with open(fname, 'a') as f:
        f.write(f"{n + nports}\n")
        f.write("\n")
        for i in range(n):
            f.write(f"{elems[i]} {pos[i,0]:.6f} {pos[i,1]:.6f} {pos[i,2]:.6f}\n")
        ip = 0
        for i in range(n):
            pi = pos[i]
            for k in range(int(port_n[i])):
                tip = pi + pneigh[i, k, :3]
                f.write(f"X {tip[0]:.6f} {tip[1]:.6f} {tip[2]:.6f}\n")
                ip += 1


def plot_state_with_ports(elems, pos, pneigh, port_n, force=None, *, title=""):
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title(title)
    ax.scatter(pos[:, 0], pos[:, 1], pos[:, 2], s=60, c='k')
    for i, e in enumerate(elems):
        ax.text(pos[i, 0], pos[i, 1], pos[i, 2], e)
    for i in range(pos.shape[0]):
        pi = pos[i]
        for k in range(int(port_n[i])):
            tip = pi + pneigh[i, k, :3]
            ax.plot([pi[0], tip[0]], [pi[1], tip[1]], [pi[2], tip[2]], '-', c='C0', lw=1)
            ax.scatter([tip[0]], [tip[1]], [tip[2]], s=20, c='C0')
    if force is not None:
        f = force[:, :3]
        ax.quiver(pos[:, 0], pos[:, 1], pos[:, 2], f[:, 0], f[:, 1], f[:, 2], length=0.1, normalize=True, color='r')
    ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
    plt.tight_layout()
    if not getattr(run_case, '_noshow', False):
        plt.show()


def run_case(xyz_path, *, dt=0.1, dt_force=0.01, iters=50, iters_force=2000, k_bond=200.0, k_rot=50.0, perturb_pos=0.1, perturb_rot=0.1, seed=0, damp_force=0.98,
             dump_xyz=None, dump_every=10, viz_force=False, viz_every=100):
    name = os.path.basename(xyz_path)
    elems, xyz0, _q = load_xyz(xyz_path)
    m = masses_from_elems(elems)

    bonds, nnode = bonds_for_molecule(name, elems)
    neighs, bks = build_neighs_bk_from_bonds(len(elems), bonds, max_deg=4)
    bLs, bKs = make_bLs_bKs_from_neighs(xyz0, neighs, k_bond=k_bond)

    quat0 = np.zeros((len(elems), 4), dtype=np.float32)
    quat0[:, 3] = 1.0

    vel0 = np.zeros((len(elems), 3), dtype=np.float32)
    omega0 = np.zeros((len(elems), 3), dtype=np.float32)

    rng = default_rng(seed)
    pos_init, quat_init = perturb_state(xyz0, quat0, perturb_pos, perturb_rot, rng)

    sim = XPDB_new(len(elems))
    sim.upload_rigid_state(pos_init, m, quat=quat_init, vel=vel0, omega=omega0)
    sim.upload_rigid_topology(neighs, bks, bLs, bKs)

    # atom port types: 0 none, 1 sp1, 2 sp2, 3 sp3
    atom_types = np.zeros((len(elems),), dtype=np.uint8)
    for i, e in enumerate(elems):
        if e in ('C', 'N'):
            atom_types[i] = 2  # sp2
        else:
            atom_types[i] = 0
    sim.upload_rigid_atom_types(atom_types)

    # Explicit-force local ports: per-slot vectors in body frame, length-scaled.
    # For this simple test we use the initial geometry as body frame (quat=identity).
    port_local = np.zeros((len(elems), 4, 4), dtype=np.float32)
    port_n = np.zeros((len(elems),), dtype=np.uint8)
    for i in range(len(elems)):
        nn = 0
        for k in range(4):
            j = int(neighs[i, k])
            if j < 0:
                continue
            v = xyz0[j] - xyz0[i]
            port_local[i, k, :3] = v
            nn += 1
        port_n[i] = nn
    sim.upload_rigid_ports_local(port_local, port_n)

    Iiso = 0.4 * m * 1.0 * 1.0

    pos4, q4, v4, om4 = sim.download_rigid_state()
    p0 = pos4[:, :3].copy()

    # 1) projective (position-based) single iteration
    sim.rigid_projective_step(nnode=nnode, dt=dt, iterations=1)
    pos4, q4, v4, om4 = sim.download_rigid_state()
    p1 = pos4[:, :3].copy()
    v1 = (p1 - p0) / float(dt)

    P1 = linear_momentum(v1, m)
    L1 = angular_momentum(p1, v1, m, om4[:, :3], Iiso)

    print(f"\n=== {name} ===")
    print(f"nAtoms={len(elems)} nnode={nnode} dt={dt} iters={iters}")
    print(f"P(after 1 iter) = {P1}")
    print(f"L(after 1 iter) = {L1}")

    # Reset to initial state
    sim.upload_rigid_state(pos_init, m, quat=quat_init, vel=vel0, omega=omega0)
    sim.upload_rigid_atom_types(atom_types)

    # A) Projective iterations
    pos_prev = pos_init.copy()
    for it in range(int(iters)):
        sim.rigid_projective_step(nnode=nnode, dt=dt, iterations=1)
        pos4, q4, v4, om4 = sim.download_rigid_state()
        p = pos4[:, :3]
        v = (p - pos_prev) / float(dt)
        pos_prev = p.copy()

        P = linear_momentum(v, m)
        L = angular_momentum(p, v, m, om4[:, :3], Iiso)

        if (it == 0) or (it == int(iters) - 1):
            print(f"proj iter={it:4d} |P|={np.linalg.norm(P):.6e} |L|={np.linalg.norm(L):.6e}")

    # convergence check: bond length residuals (host side)
    max_err = 0.0
    for i in range(len(elems)):
        for k in range(4):
            j = int(neighs[i, k])
            if j < 0:
                continue
            d = float(np.linalg.norm(pos_prev[j] - pos_prev[i]))
            err = abs(d - float(bLs[i, k]))
            if err > max_err:
                max_err = err
    print(f"max bond |d-L0| after proj iters = {max_err:.6e}")

    # B) Explicit force-based dynamics
    # Use nnode_force=natoms so every atom can have ports (e.g. O in H2O).
    nnode_force = len(elems)
    sim.upload_rigid_state(pos_init, m, quat=quat_init, vel=vel0, omega=omega0)
    sim.upload_rigid_atom_types(atom_types)
    sim.upload_rigid_ports_local(port_local, port_n)
    pos_prev = pos_init.copy()
    f_hist = []
    it_hist = []
    print_every = max(1, int(iters_force // 20))
    dump_xyz_i = None
    if dump_xyz is not None:
        base, ext = os.path.splitext(dump_xyz)
        if ext == '':
            ext = '.xyz'
        dump_xyz_i = f"{base}_{name}{ext}"
        open(dump_xyz_i, 'w').close()
    for it in range(int(iters_force)):
        sim.rigid_force_explicit_step(nnode=nnode_force, dt=dt_force, nsteps=1, damp=damp_force)
        if (it % print_every) == 0 or (it == int(iters_force) - 1):
            f4 = sim.download_buffer(sim.cl_rforce, (sim.num_atoms, 4))
            fnorm = float(np.linalg.norm(f4[:, :3]))
            f_hist.append(fnorm)
            it_hist.append(it)
            print(f"force it={it:6d} |F|={fnorm:.6e} log10|F|={np.log10(fnorm+1e-30): .3f}")

        if dump_xyz_i is not None and ((it % int(dump_every)) == 0 or (it == int(iters_force) - 1)):
            pos4_it, q4_it, v4_it, om4_it = sim.download_rigid_state()
            pneigh_it = sim.download_buffer(sim.cl_rpneigh, (sim.num_atoms * 4, 4)).reshape(sim.num_atoms, 4, 4)
            write_xyz_with_ports(dump_xyz_i, elems, pos4_it[:, :3], pneigh_it, port_n)

        if viz_force and ((it % int(viz_every)) == 0 or (it == int(iters_force) - 1)):
            pos4_it, q4_it, v4_it, om4_it = sim.download_rigid_state()
            pneigh_it = sim.download_buffer(sim.cl_rpneigh, (sim.num_atoms * 4, 4)).reshape(sim.num_atoms, 4, 4)
            f4 = sim.download_buffer(sim.cl_rforce, (sim.num_atoms, 4))
            plot_state_with_ports(elems, pos4_it[:, :3], pneigh_it, port_n, force=f4, title=f"{name} it={it}")

    pos4, q4, v4, om4 = sim.download_rigid_state()
    p = pos4[:, :3]
    v = v4[:, :3]
    P = linear_momentum(v, m)
    L = angular_momentum(p, v, m, om4[:, :3], Iiso)
    print(f"force final |P|={np.linalg.norm(P):.6e} |L|={np.linalg.norm(L):.6e}")
    pos_prev = p.copy()
    max_err = 0.0
    for i in range(len(elems)):
        for k in range(4):
            j = int(neighs[i, k])
            if j < 0:
                continue
            d = float(np.linalg.norm(pos_prev[j] - pos_prev[i]))
            err = abs(d - float(bLs[i, k]))
            if err > max_err:
                max_err = err
    print(f"max bond |d-L0| after force iters = {max_err:.6e}")

    if len(f_hist) > 1:
        plt.figure(figsize=(6, 4))
        plt.plot(it_hist, np.log10(np.array(f_hist) + 1e-30))
        plt.xlabel('iter')
        plt.ylabel('log10 |F|')
        plt.title(f'force convergence: {name}')
        plt.tight_layout()
        if not getattr(run_case, '_noshow', False):
            plt.show()


'''
python3 pyBall/XPDB_AVBD/test_rigid_XPBD_molecules.py  --dt_force 0.001 --iters_force 1000 --perturb_pos 0.1 --perturb_rot 0.1 --damp_force 0.98 --viz_force --viz_every 50

python3 pyBall/XPDB_AVBD/test_rigid_XPBD_molecules.py --dt_force 0.001 --iters_force 2000 --perturb_pos 0.1 --perturb_rot 0.1 --damp_force 0.98 --dump_xyz traj.xyz --dump_every 10 --noshow

'''


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dt', type=float, default=0.1)
    ap.add_argument('--dt_force', type=float, default=0.01)
    ap.add_argument('--iters', type=int, default=50)
    ap.add_argument('--iters_force', type=int, default=2000)
    ap.add_argument('--k_bond', type=float, default=200.0)
    ap.add_argument('--k_rot', type=float, default=50.0)
    ap.add_argument('--perturb_pos', type=float, default=0.1)
    ap.add_argument('--perturb_rot', type=float, default=0.1)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--damp_force', type=float, default=0.98)
    ap.add_argument('--dump_xyz', type=str, default=None)
    ap.add_argument('--dump_every', type=int, default=10)
    ap.add_argument('--viz_force', action='store_true')
    ap.add_argument('--viz_every', type=int, default=100)
    ap.add_argument('--noshow', action='store_true')
    args = ap.parse_args()

    base = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'cpp', 'common_resources', 'xyz'))
    xyz_h2o = os.path.join(base, 'H2O.xyz')
    xyz_ch2nh = os.path.join(base, 'CH2NH.xyz')
    run_case._noshow = bool(args.noshow)
    run_case(xyz_h2o, dt=args.dt, dt_force=args.dt_force, iters=args.iters, iters_force=args.iters_force,
             k_bond=args.k_bond, k_rot=args.k_rot,
             perturb_pos=args.perturb_pos, perturb_rot=args.perturb_rot, seed=args.seed, damp_force=args.damp_force,
             dump_xyz=args.dump_xyz, dump_every=args.dump_every, viz_force=args.viz_force, viz_every=args.viz_every)
    run_case(xyz_ch2nh, dt=args.dt, dt_force=args.dt_force, iters=args.iters, iters_force=args.iters_force,
             k_bond=args.k_bond, k_rot=args.k_rot,
             perturb_pos=args.perturb_pos, perturb_rot=args.perturb_rot, seed=args.seed + 1, damp_force=args.damp_force,
             dump_xyz=args.dump_xyz, dump_every=args.dump_every, viz_force=args.viz_force, viz_every=args.viz_every)


if __name__ == '__main__':
    main()
