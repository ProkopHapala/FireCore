import numpy as np
import argparse
import os
import sys

from numpy.random import default_rng

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

_PYBALL = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if _PYBALL not in sys.path:
    sys.path.insert(0, _PYBALL)

from XPDB_new import XPDB_new, build_neighs_bk_from_bonds, make_bLs_bKs_from_neighs
from XPTB_utils import load_xyz, masses_from_elems, perturb_state


def quat_to_mat3(q):
    q = np.asarray(q, dtype=np.float32)
    x, y, z, w = float(q[0]), float(q[1]), float(q[2]), float(q[3])
    n = x*x + y*y + z*z + w*w
    if n == 0.0:
        raise ValueError("quat_to_mat3: zero quaternion")
    s = 2.0 / n
    xx, yy, zz = x*x*s, y*y*s, z*z*s
    xy, xz, yz = x*y*s, x*z*s, y*z*s
    wx, wy, wz = w*x*s, w*y*s, w*z*s
    return np.array([
        [1.0 - (yy + zz), xy - wz, xz + wy],
        [xy + wz, 1.0 - (xx + zz), yz - wx],
        [xz - wy, yz + wx, 1.0 - (xx + yy)],
    ], dtype=np.float32)


def mat3_det(R):
    return float(np.linalg.det(np.asarray(R, dtype=np.float32)))


def procrustes_cpu(r_local, p_world, w):
    r_local = np.asarray(r_local, dtype=np.float32)
    p_world = np.asarray(p_world, dtype=np.float32)
    w = np.asarray(w, dtype=np.float32)
    if r_local.shape != p_world.shape:
        raise ValueError(f"procrustes_cpu: shape mismatch r_local={r_local.shape} p_world={p_world.shape}")
    if r_local.ndim != 2 or r_local.shape[1] != 3:
        raise ValueError(f"procrustes_cpu: expected (m,3), got {r_local.shape}")
    if w.shape != (r_local.shape[0],):
        raise ValueError(f"procrustes_cpu: w.shape={w.shape} expected ({r_local.shape[0]},)")

    sw = float(np.sum(w))
    if not np.isfinite(sw) or sw <= 0.0:
        raise ValueError(f"procrustes_cpu: invalid sum_w={sw}")

    cr = (w[:, None] * r_local).sum(axis=0) / sw
    cp = (w[:, None] * p_world).sum(axis=0) / sw

    r0 = r_local - cr[None, :]
    p0 = p_world - cp[None, :]

    H = (w[:, None, None] * p0[:, :, None] * r0[:, None, :]).sum(axis=0)

    U, S, Vt = np.linalg.svd(H, full_matrices=False)
    R = U @ Vt
    if mat3_det(R) < 0.0:
        U[:, 2] *= -1.0
        R = U @ Vt

    t = cp - R @ cr
    return t.astype(np.float32), R.astype(np.float32)


def energy(t, R, r_local, p_world, w):
    r_local = np.asarray(r_local, dtype=np.float32)
    p_world = np.asarray(p_world, dtype=np.float32)
    w = np.asarray(w, dtype=np.float32)
    Rt = (R @ r_local.T).T
    d = (t[None, :] + Rt) - p_world
    return float(np.sum(w * np.sum(d*d, axis=1)))


'''
python test_procrustes_ch4.py --perturb_pos 0.1 --perturb_rot 0.1 --seed 0 --method both
'''

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--perturb_pos', type=float, default=0.1)
    ap.add_argument('--perturb_rot', type=float, default=0.1)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--method', type=str, default='both', choices=['quat', 'mat', 'both'])
    args = ap.parse_args()

    base = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'cpp', 'common_resources', 'xyz'))
    xyz_path = os.path.join(base, 'CH4.xyz')

    elems, xyz0, _q = load_xyz(xyz_path)
    m = masses_from_elems(elems)

    if elems[0] != 'C' or any(e != 'H' for e in elems[1:]):
        raise RuntimeError(f"Unexpected atom order in CH4.xyz elems={elems}")

    bonds = [(0, 1), (0, 2), (0, 3), (0, 4)]
    neighs, bks = build_neighs_bk_from_bonds(len(elems), bonds, max_deg=4)

    # For this test, stiffness just provides weights; bond lengths are irrelevant.
    bLs, bKs = make_bLs_bKs_from_neighs(xyz0, neighs, k_bond=1.0)

    nnode = 1
    port_local = np.zeros((nnode, 4, 4), dtype=np.float32)
    port_n = np.zeros((nnode,), dtype=np.uint8)
    w = np.zeros((4,), dtype=np.float32)

    nn = 0
    for k in range(4):
        j = int(neighs[0, k])
        if j < 0:
            continue
        v = xyz0[j] - xyz0[0]
        port_local[0, k, :3] = v
        w[nn] = float(bKs[0, k])
        nn += 1
    if nn != 4:
        raise RuntimeError(f"CH4 expected 4 neighbors for node 0, got nn={nn}, neighs[0]={neighs[0]}")
    port_n[0] = 4

    # Ground truth pose: keep identity rotation and translation at origin (rest), then perturb everything.
    quat0 = np.zeros((len(elems), 4), dtype=np.float32)
    quat0[:, 3] = 1.0

    rng = default_rng(int(args.seed))
    pos_init, quat_init = perturb_state(xyz0, quat0, float(args.perturb_pos), float(args.perturb_rot), rng)

    sim = XPDB_new(len(elems))
    sim.upload_rigid_state(pos_init, m, quat=quat_init)
    sim.upload_rigid_topology(neighs, bks, bLs, bKs)
    sim.upload_rigid_ports_local(port_local, port_n, nnode=nnode)
    sim.upload_rigid_node_stiffness_flat(bKs, nnode=nnode)

    # CPU reference uses the *current* neighbor positions (after perturbation)
    r_local = port_local[0, :, :3].copy()
    p_world = pos_init[neighs[0, :], :3].copy()
    w_ref = np.ones((4,), dtype=np.float32)
    t_ref, R_ref = procrustes_cpu(r_local, p_world, w_ref)

    # Initial guess from current state
    pos4, q4, *_ = sim.download_rigid_state()
    t0 = pos4[0, :3].copy()
    R0 = quat_to_mat3(q4[0])

    E0 = energy(t0, R0, r_local, p_world, w_ref)
    Eref = energy(t_ref, R_ref, r_local, p_world, w_ref)

    print("=== CH4 single-node Procrustes test ===")
    print(f"perturb_pos={args.perturb_pos} perturb_rot={args.perturb_rot} seed={args.seed}")
    print(f"E(initial)={E0:.6e} E(cpu_ref)={Eref:.6e} ratio={E0/(Eref+1e-30):.6e}")

    if not np.isfinite(E0) or not np.isfinite(Eref):
        raise RuntimeError("Non-finite energy detected")

    if args.method in ('quat', 'both'):
        sim.rigid_procrustes_quat(nnode=nnode)
        tpos = sim.download_rigid_procrustes_targets(nnode=nnode, with_rot=False)
        pos4, q4, *_ = sim.download_rigid_state()  # quaternion updated in-place
        t_est = tpos[0, :3].copy()
        R_est = quat_to_mat3(q4[0])

        Eq = energy(t_est, R_est, r_local, p_world, w_ref)
        EqT = energy(t_est, R_est.T, r_local, p_world, w_ref)
        dt = float(np.linalg.norm(t_est - t_ref))
        dR = float(np.linalg.norm(R_est - R_ref))
        dRT = float(np.linalg.norm(R_est.T - R_ref))

        print("-- quat kernel --")
        print(f"E(quat)={Eq:.6e}  E/Eref={Eq/(Eref+1e-30):.6e}  |t-tref|={dt:.6e}  |R-Rref|_F={dR:.6e} det(R)={mat3_det(R_est):.6f}")
        print(f"E(quat,R^T)={EqT:.6e}  E/Eref={EqT/(Eref+1e-30):.6e}  |R^T-Rref|_F={dRT:.6e}")
        if not np.isfinite(Eq):
            raise RuntimeError("quat: non-finite energy")

    if args.method in ('mat', 'both'):
        sim.rigid_procrustes_mat(nnode=nnode)
        out = sim.download_rigid_procrustes_targets(nnode=nnode, with_rot=True)
        tpos, r0, r1, r2 = out
        t_est = tpos[0, :3].copy()
        R_est = np.stack([r0[0, :3], r1[0, :3], r2[0, :3]], axis=0).astype(np.float32)

        Em = energy(t_est, R_est, r_local, p_world, w_ref)
        dt = float(np.linalg.norm(t_est - t_ref))
        dR = float(np.linalg.norm(R_est - R_ref))

        print("-- mat kernel --")
        print(f"E(mat)={Em:.6e}  E/Eref={Em/(Eref+1e-30):.6e}  |t-tref|={dt:.6e}  |R-Rref|_F={dR:.6e} det(R)={mat3_det(R_est):.6f}")
        if not np.isfinite(Em):
            raise RuntimeError("mat: non-finite energy")

    print("OK")


if __name__ == '__main__':
    main()
