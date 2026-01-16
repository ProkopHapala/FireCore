#!/usr/bin/env python3

import os
import sys
import argparse
import numpy as np

import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(__file__), "../.."))

from pyBall.AtomicSystem import AtomicSystem
from pyBall import atomicUtils as au
from pyBall.OCL.eFF_ocl import EFF_OCL


def load_density_xyz(fname):
    with open(fname, 'r') as f:
        lines = f.readlines()
    n = int(lines[0].strip())
    pts = np.zeros((n, 4), dtype=np.float32)
    labels = []
    for i in range(n):
        w = lines[i + 2].split()
        labels.append(w[0])
        pts[i, 0] = float(w[1]); pts[i, 1] = float(w[2]); pts[i, 2] = float(w[3]); pts[i, 3] = float(w[4])
    return pts, labels


def infer_mol_from_density_samples(grid, labels):
    apos = []
    enames = []
    for p, lab in zip(grid, labels):
        if not lab.startswith('SMP.ATOM.'):
            continue
        # label e.g. SMP.ATOM.O0 or SMP.ATOM.H12
        tail = lab[len('SMP.ATOM.') :]
        if len(tail) == 0:
            continue
        e = tail[0]
        if e.isalpha() and len(tail) > 1 and tail[1].isalpha():
            e = tail[:2]
        apos.append(p[:3].astype(np.float64))
        enames.append(e)
    if len(apos) == 0:
        raise RuntimeError('No SMP.ATOM.* records found in density xyz; cannot infer molecule geometry')
    apos = np.array(apos, dtype=np.float64)
    atypes = np.zeros((len(apos),), dtype=np.int32)
    mol = AtomicSystem(apos=apos, atypes=atypes, enames=enames, bPreinit=True)
    return mol


def make_initial_basis_from_mol(mol, s0=1.0, q_atom=1.0, q_bond=1.0, q_epair=1.0, epair_dist=0.5, bond_fracs=(0.5,)):
    if mol.bonds is None:
        mol.findBonds()
    if mol.ngs is None:
        mol.neighs()
    n0 = len(mol.apos)

    centers = []
    amps = []

    # atom centers
    for i in range(n0):
        centers.append(np.array([mol.apos[i, 0], mol.apos[i, 1], mol.apos[i, 2], s0], dtype=np.float32))
        amps.append(q_atom)

    # bond midpoints
    if mol.bonds is not None and len(mol.bonds) > 0:
        for (i, j) in mol.bonds:
            ri = mol.apos[int(i)]; rj = mol.apos[int(j)]
            d = (rj - ri)
            for t in bond_fracs:
                p = ri + float(t) * d
                centers.append(np.array([p[0], p[1], p[2], s0], dtype=np.float32))
                amps.append(q_bond)

    # epairs (reuse AtomicSystem logic)
    molE = mol.clonePBC()
    molE.bonds = np.array(mol.bonds, dtype=np.int32) if mol.bonds is not None else None
    molE.neighs()
    molE.add_electron_pairs()
    n1 = len(molE.apos)
    if n1 > n0:
        for ia in range(n0, n1):
            p = molE.apos[ia]
            centers.append(np.array([p[0], p[1], p[2], s0], dtype=np.float32))
            amps.append(q_epair)

    centers = np.array(centers, dtype=np.float32)
    amps = np.array(amps, dtype=np.float32)
    return centers, amps


def eval_model_density(points_xyz, centers, amps):
    # matches kernel: gauss_density = exp(-2 r^2 / s^2)
    pts = points_xyz[:, :3].astype(np.float64)
    rho = np.zeros((pts.shape[0],), dtype=np.float64)
    for j in range(centers.shape[0]):
        c = centers[j, :3].astype(np.float64)
        s = float(centers[j, 3])
        q = float(amps[j])
        d = pts - c[None, :]
        r2 = np.sum(d * d, axis=1)
        rho += q * np.exp(-2.0 * r2 / (s * s + 1e-16))
    return rho


def eval_loss(points_xyz, rho_model):
    diff = rho_model - points_xyz[:, 3]
    E = float(np.dot(diff, diff))
    rms = float(np.sqrt(E / max(1, diff.size)))
    mae = float(np.mean(np.abs(diff)))
    return E, rms, mae


def eval_force_norms(f4, fq):
    # f4: [ne,4] for xyzs ; fq: [ne]
    if f4 is None or fq is None:
        return None, None
    f2 = float(np.sum(f4.astype(np.float64) ** 2) + np.sum(fq.astype(np.float64) ** 2))
    fnorm = float(np.sqrt(f2))
    finf = float(max(np.max(np.abs(f4)), np.max(np.abs(fq))))
    return fnorm, finf


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dens', default=os.path.join(os.path.dirname(__file__), 'dens_outputs', 'H2O_dens.xyz'))
    ap.add_argument('--mol',  default=os.path.join(os.path.dirname(__file__), 'dens_outputs', 'H2O.xyz'))
    ap.add_argument('--nsteps',  type=int,   default=1000)
    ap.add_argument('--m',       type=int,   default=1, help='Print diagnostics every m optimizer steps (implemented by chunking)')
    ap.add_argument('--dt',      type=float, default=0.005)
    ap.add_argument('--opt', default='fire', choices=['fire','md','gd'])
    ap.add_argument('--md_damp', type=float, default=0.9)
    ap.add_argument('--fire_inc', type=float, default=1.1)
    ap.add_argument('--fire_dec', type=float, default=0.5)
    ap.add_argument('--fire_alpha', type=float, default=0.01)
    ap.add_argument('--fire_dt_max', type=float, default=0.01)

    ap.add_argument('--k_pos',   type=float, default=0.0)
    ap.add_argument('--k_size',  type=float, default=0.0)
    ap.add_argument('--k_q',     type=float, default=0.0)
    ap.add_argument('--s0',      type=float, default=1.0)
    ap.add_argument('--q_atom',  type=float, default=1.0)
    ap.add_argument('--q_bond',  type=float, default=1.0)
    ap.add_argument('--q_epair', type=float, default=1.0)
    ap.add_argument('--local_size', type=int, default=32)
    ap.add_argument('--dump_xyz', default='')
    ap.add_argument('--plot', type=int, default=1)
    args = ap.parse_args()

    grid, labels = load_density_xyz(args.dens)

    if args.mol and os.path.exists(args.mol):
        mol = AtomicSystem(args.mol)
    else:
        mol = infer_mol_from_density_samples(grid, labels)

    centers0, amps0 = make_initial_basis_from_mol(mol, s0=args.s0, q_atom=args.q_atom, q_bond=args.q_bond, q_epair=args.q_epair)

    rho0 = eval_model_density(grid, centers0, amps0)
    E0, rms0, mae0 = eval_loss(grid, rho0)
    print(f"initial: n_basis={len(amps0)} n_grid={len(grid)} E={E0:.6e} rms={rms0:.6e} mae={mae0:.6e}")

    ocl = EFF_OCL(nloc=args.local_size)
    centers = centers0.copy()
    amps = amps0.copy()

    vel = np.zeros_like(centers, dtype=np.float32)
    vq  = np.zeros((amps.shape[0],), dtype=np.float32)
    fire_state = np.zeros((1, 4), dtype=np.float32)
    fire_np    = np.zeros((1,), dtype=np.int32)

    m = int(args.m)
    if m <= 0:
        raise ValueError(f"m must be > 0, got {m}")
    nsteps = int(args.nsteps)
    nhist = []
    Eh    = []

    istep = 0
    opt_mode = {'fire':0,'md':1,'gd':2}[args.opt]
    while istep < nsteps:
        chunk = min(m, nsteps - istep)
        centers, amps, vel, vq, fire_state, fire_np, f4, fq = ocl.fit_density_fire(
            centers,
            grid,
            amps=amps,
            vel=vel,
            vq=vq,
            fire_state=fire_state,
            fire_np=fire_np,
            nsteps=chunk,
            dt=args.dt,
            k_pos=args.k_pos,
            k_size=args.k_size,
            k_q=args.k_q,
            fire_inc=args.fire_inc,
            fire_dec=args.fire_dec,
            fire_alpha=args.fire_alpha,
            fire_dt_max=args.fire_dt_max,
            opt_mode=opt_mode,
            md_damp=args.md_damp,
            local_size=args.local_size,
            bForces=True,
        )

        rho = eval_model_density(grid, centers, amps)
        E, rms, mae = eval_loss(grid, rho)
        fnorm, finf = eval_force_norms(f4, fq)
        istep += chunk
        Eh.append(E)
        nhist.append((istep, fnorm))
        print(f"step {istep:6d}  E {E:12.6e}  |F| {fnorm:12.6e}  F_inf {finf:12.6e}  dt {float(fire_state[0,0]):9.3e}  alpha {float(fire_state[0,1]):9.3e}  n_pos {int(fire_np[0])}")

    rho1 = eval_model_density(grid, centers, amps)
    E1, rms1, mae1 = eval_loss(grid, rho1)
    print(f"final:   n_basis={len(amps)} n_grid={len(grid)} E={E1:.6e} rms={rms1:.6e} mae={mae1:.6e}")

    if args.dump_xyz:
        enames = ["G"] * centers.shape[0]
        apos = centers[:, :3].astype(np.float64)
        qs = amps.astype(np.float64)
        Rs = centers[:, 3].astype(np.float64)
        au.saveXYZ(enames, apos, args.dump_xyz, qs=qs, Rs=Rs, comment="gaussian_basis_fit")
        print(f"wrote: {args.dump_xyz}")

    if args.plot and len(nhist) > 0:
        steps = np.array([s for (s, _) in nhist], dtype=np.int32)
        fn = np.array([f for (_, f) in nhist], dtype=np.float64)
        plt.figure(figsize=(6, 4))
        plt.plot(steps, fn)
        plt.yscale('log')
        plt.xlabel('step')
        plt.ylabel('|F|')
        plt.title('FIRE convergence')
        plt.grid(True)
        plt.tight_layout()

        # 2D slice through nuclei plane (approx): span H-O-H with normal via cross
        apos = mol.apos.astype(np.float64)
        c0 = apos.mean(axis=0)
        a = apos[1] - apos[0]
        b = apos[2] - apos[0]
        n = np.cross(a, b)
        nn = np.linalg.norm(n)
        if nn < 1e-8:
            n = np.array([0.0, 0.0, 1.0])
        else:
            n /= nn
        ex = a / (np.linalg.norm(a) + 1e-16)
        ey = np.cross(n, ex); ey /= (np.linalg.norm(ey) + 1e-16)
        L = 3.0
        N = 200
        xs = np.linspace(-L, L, N)
        ys = np.linspace(-L, L, N)
        XY = np.stack(np.meshgrid(xs, ys, indexing='xy'), axis=-1).reshape(-1, 2)
        pts = c0[None, :] + XY[:, 0:1] * ex[None, :] + XY[:, 1:2] * ey[None, :]
        grid2d = np.zeros((pts.shape[0], 4), dtype=np.float32)
        grid2d[:, 0:3] = pts.astype(np.float32)
        rho2 = ocl.eval_density_grid(centers, grid2d, amps=amps, local_size=128).reshape(N, N)
        plt.figure(figsize=(6, 5))
        plt.imshow(rho2, origin='lower', extent=[xs[0], xs[-1], ys[0], ys[-1]], aspect='equal')
        plt.colorbar(label='rho')
        plt.title('Fitted density slice')
        plt.xlabel('x in plane')
        plt.ylabel('y in plane')
        plt.show()


if __name__ == '__main__':
    main()
