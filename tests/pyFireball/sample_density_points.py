#!/usr/bin/env python3

import os
import sys
import argparse
import numpy as np

sys.path.append(os.path.join(os.path.dirname(__file__), "../.."))
from pyBall.AtomicSystem import AtomicSystem
from pyBall import atomicUtils as au


def _parse_floats_csv(s):
    if s is None or s == "":
        return []
    return [float(x) for x in s.split(',') if x.strip() != ""]


def _normalize(v):
    return au.normalize(v)


def _dedup_points(points, labels, hosts, bonds, tol):
    if tol is None or tol <= 0:
        return points, labels, hosts, bonds
    inv = 1.0 / tol
    hmap = {}
    keep = []
    old2new = {}
    for i, p in enumerate(points):
        key = (int(np.floor(p[0] * inv + 0.5)), int(np.floor(p[1] * inv + 0.5)), int(np.floor(p[2] * inv + 0.5)))
        j = hmap.get(key)
        if j is None:
            hmap[key] = i
            old2new[i] = len(keep)
            keep.append(i)
        else:
            # keep first, drop duplicates
            old2new[i] = old2new[j]
    pts2 = points[keep]
    lab2 = [labels[i] for i in keep]
    host2 = [hosts[i] for i in keep]

    bonds2 = []
    for a, b in bonds:
        a2 = old2new.get(a)
        b2 = old2new.get(b)
        if a2 is None or b2 is None:
            continue
        if a2 == b2:
            continue
        bonds2.append((a2, b2))

    return pts2, lab2, host2, bonds2


def _enumerate_angles(ngs):
    # ngs is list of dicts: ngs[j].keys() are neighbors
    angs = []
    for j, gj in enumerate(ngs):
        neis = [int(i) for i in gj.keys() if int(i) >= 0]
        neis.sort()
        n = len(neis)
        if n < 2:
            continue
        for a in range(n):
            i = neis[a]
            for b in range(a + 1, n):
                k = neis[b]
                angs.append((i, j, k))
    return angs


def build_samples(mol, opts):
    """Return (points[n,3], labels[n], hosts[n]) where hosts are tuples of atom indices used to define the sample."""

    # Ensure topology
    if mol.bonds is None:
        mol.findBonds(RvdwCut=opts.bond_rvdwcut)
    if mol.ngs is None:
        mol.neighs()
    if mol.Rs is None:
        mol.preinitialize_atomic_properties()

    n0 = len(mol.apos)

    pts = []
    labels = []
    hosts = []

    # ---------------- Atom centers
    if opts.atoms:
        for i in range(n0):
            pts.append(mol.apos[i].copy())
            labels.append(f"SMP.ATOM.{mol.enames[i]}{i}")
            hosts.append((i,))

    # ---------------- Bond points
    if opts.bonds:
        fracs = opts.bond_fracs
        if len(fracs) == 0:
            fracs = [0.5]
        for (i, j) in mol.bonds:
            i = int(i)
            j = int(j)
            ri = mol.apos[i]
            rj = mol.apos[j]
            d = (rj - ri)
            for t in fracs:
                p = ri + t * d
                pts.append(p)
                labels.append(f"SMP.BOND.{mol.enames[i]}{i}-{mol.enames[j]}{j}@{t:.2f}")
                hosts.append((i, j))
                if opts.invert_bonds:
                    p2 = ri - t * d
                    pts.append(p2)
                    labels.append(f"SMP.BONDINV.{mol.enames[i]}{i}-{mol.enames[j]}{j}@{t:.2f}")
                    hosts.append((i, j))

    # ---------------- Epairs via AtomicSystem (reuse)
    epair_pts = []
    epair_hosts = []
    if opts.epairs or opts.far_epairs or opts.invert_epairs or opts.angles:
        # Clone to avoid mutating mol
        molE = mol.clonePBC()
        molE.bonds = np.array(mol.bonds, dtype=np.int32) if mol.bonds is not None else None
        molE.neighs()
        molE.add_electron_pairs()
        n1 = len(molE.apos)
        if n1 > n0:
            # infer epair->host from new bonds
            for a, b in molE.bonds:
                a = int(a)
                b = int(b)
                if a >= n0 and b < n0:
                    epair_pts.append(molE.apos[a].copy())
                    epair_hosts.append(b)
                elif b >= n0 and a < n0:
                    epair_pts.append(molE.apos[b].copy())
                    epair_hosts.append(a)

    if opts.epairs:
        for p, ih in zip(epair_pts, epair_hosts):
            ih = int(ih)
            v = p - mol.apos[ih]
            base_dist = float(np.linalg.norm(v))
            d = _normalize(v) if base_dist > 1e-12 else np.array([1.0, 0.0, 0.0])
            dist_list = opts.epair_dists if (opts.epair_dists and len(opts.epair_dists) > 0) else [base_dist]
            for dd in dist_list:
                p_new = mol.apos[ih] + d * dd
                pts.append(p_new)
                labels.append(f"SMP.EPAIR.{mol.enames[ih]}{ih}@{dd:.3f}")
                hosts.append((ih,))
                if opts.invert_epairs:
                    p_inv = mol.apos[ih] - d * dd
                    pts.append(p_inv)
                    labels.append(f"SMP.EPAIRINV.{mol.enames[ih]}{ih}@{dd:.3f}")
                    hosts.append((ih,))

    # ---------------- Far epair prolongations
    if opts.far_epairs:
        for p, ih in zip(epair_pts, epair_hosts):
            ih = int(ih)
            d = _normalize(p - mol.apos[ih])
            pf = mol.apos[ih] + d * opts.far_dist
            pts.append(pf)
            labels.append(f"SMP.FAREP.{mol.enames[ih]}{ih}@{opts.far_dist:.2f}")
            hosts.append((ih,))
            if opts.invert_far:
                pf2 = mol.apos[ih] - d * opts.far_dist
                pts.append(pf2)
                labels.append(f"SMP.FAREPINV.{mol.enames[ih]}{ih}@{opts.far_dist:.2f}")
                hosts.append((ih,))

    # ---------------- Apex points (sigma-hole like)
    if opts.apex:
        for i in range(n0):
            gi = mol.ngs[i]
            neis = [int(j) for j in gi.keys() if int(j) >= 0]
            if len(neis) == 0:
                continue
            ri = mol.apos[i]
            if opts.apex_mode == 'avg':
                v = np.zeros(3)
                for j in neis:
                    v += _normalize(mol.apos[j] - ri)
                d = _normalize(v)
                if np.linalg.norm(d) < 1e-8:
                    continue
                p = ri - d * opts.apex_dist
                pts.append(p)
                labels.append(f"SMP.APEXAVG.{mol.enames[i]}{i}@{opts.apex_dist:.2f}")
                hosts.append((i,))
            elif opts.apex_mode == 'per-bond':
                for j in neis:
                    d = _normalize(mol.apos[j] - ri)
                    p = ri - d * opts.apex_dist
                    pts.append(p)
                    labels.append(f"SMP.APEXBOND.{mol.enames[i]}{i}-{mol.enames[j]}{j}@{opts.apex_dist:.2f}")
                    hosts.append((i, j))
            else:
                raise ValueError(f"Unknown apex_mode={opts.apex_mode}")

    # ---------------- Angle bisectors (includes epair legs only if epairs were actually added into topology; here we use original mol.ngs)
    if opts.angles:
        angs = _enumerate_angles(mol.ngs)
        for (i, j, k) in angs:
            rj = mol.apos[j]
            v1 = _normalize(mol.apos[i] - rj)
            v2 = _normalize(mol.apos[k] - rj)
            b = _normalize(v1 + v2)
            if np.linalg.norm(b) < 1e-8:
                continue
            p = rj + b * opts.angle_dist
            pts.append(p)
            labels.append(f"SMP.ANG.{mol.enames[i]}{i}-{mol.enames[j]}{j}-{mol.enames[k]}{k}@{opts.angle_dist:.2f}")
            hosts.append((i, j, k))

    pts = np.array(pts, dtype=np.float64)

    # Visualization bonds among samples: none; we only connect each sample to its host atom when saving
    samp_bonds = []

    return pts, labels, hosts, samp_bonds


def save_augmented(mol, pts, labels, hosts, out_xyz, out_mol2):
    n0 = len(mol.apos)
    nS = len(pts)

    apos = np.vstack([mol.apos, pts]) if nS > 0 else mol.apos.copy()
    enames = list(mol.enames)
    atypes = np.array(mol.atypes, copy=True) if mol.atypes is not None else None
    qs = np.array(mol.qs, copy=True) if mol.qs is not None else None

    # Add sample “atoms” with label encoded in ename.
    for lab in labels:
        enames.append(lab)

    if atypes is not None:
        atypes = np.concatenate([atypes, np.zeros(nS, dtype=atypes.dtype)])
    if qs is not None:
        qs = np.concatenate([qs, np.zeros(nS, dtype=qs.dtype)])

    # Combine bonds: original + sample->host
    bonds = []
    if mol.bonds is not None:
        for (i, j) in mol.bonds:
            bonds.append((int(i), int(j)))

    for isamp, h in enumerate(hosts):
        ia = n0 + isamp
        if len(h) >= 1:
            bonds.append((int(h[0]), ia))

    out = AtomicSystem(apos=apos, atypes=atypes, enames=enames, bonds=np.array(bonds, dtype=np.int32), qs=qs)

    if out_xyz is not None:
        out.saveXYZ(out_xyz)

    if out_mol2 is not None:
        out.save_mol2(out_mol2)

    return out


def plot_3d(mol, pts, labels, no_plot=False):
    if no_plot:
        return
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(111, projection='3d')

    # atoms
    ax.scatter(mol.apos[:, 0], mol.apos[:, 1], mol.apos[:, 2], c='k', s=60, label='atoms')

    if len(pts) > 0:
        # group by prefix SMP.X
        cats = {}
        for p, lab in zip(pts, labels):
            key = lab.split('.')[1] if '.' in lab else lab
            cats.setdefault(key, []).append(p)
        colors = {
            'ATOM': 'gray',
            'BOND': 'tab:blue',
            'BONDINV': 'tab:cyan',
            'EPAIR': 'tab:orange',
            'EPAIRINV': 'tab:red',
            'FAREP': 'tab:green',
            'FAREPINV': 'tab:olive',
            'APEXAVG': 'tab:purple',
            'APEXBOND': 'tab:pink',
            'ANG': 'tab:brown',
        }
        for k, ps in cats.items():
            ps = np.array(ps)
            ax.scatter(ps[:, 0], ps[:, 1], ps[:, 2], s=20, c=colors.get(k, 'tab:gray'), label=k)

    # enforce equal aspect by setting symmetric limits
    all_pts = [mol.apos]
    if len(pts) > 0:
        all_pts.append(np.array(pts))
    coords = np.vstack(all_pts)
    mins = coords.min(axis=0)
    maxs = coords.max(axis=0)
    center = 0.5 * (mins + maxs)
    span = (maxs - mins).max()
    half = span * 0.55  # a bit of padding
    for a, c0 in zip(['x', 'y', 'z'], center):
        getattr(ax, f"set_{a}lim")(c0 - half, c0 + half)

    ax.set_xlabel('x [Å]')
    ax.set_ylabel('y [Å]')
    ax.set_zlabel('z [Å]')
    ax.legend(loc='best')
    ax.set_box_aspect([1, 1, 1])
    plt.tight_layout()
    plt.show()


def main():
    ap = argparse.ArgumentParser(description="Generate debug sampling points around a molecule for later Fireball density sampling")
    ap.add_argument('--mol',    default='H2O', help='Molecule name (without .xyz) in tests/pyFireball/molecules/ (default H2O)')
    #ap.add_argument('--mol',    default='CH4', help='Molecule name (without .xyz) in tests/pyFireball/molecules/ (default H2O)')
    #ap.add_argument('--mol',    default='NH3', help='Molecule name (without .xyz) in tests/pyFireball/molecules/ (default H2O)')


    ap.add_argument('--xyz',    default=None,  help='Explicit xyz path (overrides --mol)')
    ap.add_argument('--atoms',         type=int, default=1, help='Add atom-center sample points (duplicates positions)')
    ap.add_argument('--bonds',         type=int, default=1, help='Add bond fraction points')
    ap.add_argument('--epairs',        type=int, default=1, help='Add electron-pair points (generated by AtomicSystem.add_electron_pairs)')
    ap.add_argument('--epair-dists',   default='0.5,1.0', help='Comma-separated distances for epairs (Å); default = use inherent epair distance')
    ap.add_argument('--invert-bonds',  type=int, default=1, help='Add inverted bond points (from host atom opposite direction)')
    ap.add_argument('--angles',        type=int, default=1, help='Add angle-bisector points from bonded topology')
    ap.add_argument('--apex',          type=int, default=1, help='Add apex (sigma-hole-like) points opposite neighbor directions')
    ap.add_argument('--invert-far',    type=int, default=1, help='Also add inverted far points')
    ap.add_argument('--far-epairs',    type=int, default=1, help='Add far points by prolonging epair directions')
    ap.add_argument('--invert-epairs', type=int, default=1, help='Add inverted electron-pair points')
    ap.add_argument('--apex-dist',     type=float, default=0.8, help='Distance for apex points (Å)')
    ap.add_argument('--bond-fracs',    default='0.5', help='Comma-separated fractions along bond from atom i to j (default 0.5)')
    ap.add_argument('--bond-rvdwcut',  type=float, default=1.5, help='RvdwCut for findBonds (default 1.5)')
    ap.add_argument('--angle-dist',    type=float, default=0.7, help='Distance along bisector for angle points (Å)')
    ap.add_argument('--apex-mode',     choices=['avg', 'per-bond'], default='avg')
    ap.add_argument('--far-dist',      type=float, default=6.0, help='Distance for far points (Å)')
    ap.add_argument('--dedup',         type=float, default=1e-6, help='Deduplicate points by rounding with tolerance (0 disables)')
    ap.add_argument('--drop-if-nonhost-nearest', type=int, default=0, help='If 1, drop samples whose nearest atom is not among the listed hosts')
    ap.add_argument('--out-xyz',       default='out_samples.xyz', help='Output xyz path (default out_samples.xyz)')
    ap.add_argument('--out-mol2',      default='out_samples.mol2', help='Output mol2 path (default out_samples.mol2)')
    ap.add_argument('--noPlot',        action='store_true', help='Disable matplotlib visualization')

    opts = ap.parse_args()

    opts.bond_fracs = _parse_floats_csv(opts.bond_fracs)
    opts.epair_dists = _parse_floats_csv(opts.epair_dists)

    if opts.xyz is None:
        xyz = os.path.join(os.path.dirname(__file__), 'molecules', f"{opts.mol}.xyz")
    else:
        xyz = opts.xyz

    if not os.path.exists(xyz):
        raise FileNotFoundError(f"XYZ not found: {xyz}")

    mol = AtomicSystem(fname=xyz)
    mol.findBonds(RvdwCut=opts.bond_rvdwcut)
    mol.neighs()

    pts, labels, hosts, samp_bonds = build_samples(mol, opts)

    # dedup
    pts, labels, hosts, samp_bonds = _dedup_points(pts, labels, hosts, samp_bonds, opts.dedup)

    # nearest-atom distances
    atom_pos = mol.apos
    keep_mask = []
    dists = []
    for p, h in zip(pts, hosts):
        dv = atom_pos - p[None, :]
        r2 = np.einsum('ij,ij->i', dv, dv)
        imin = int(np.argmin(r2))
        dmin = float(np.sqrt(r2[imin]))
        dists.append(dmin)
        if opts.drop_if_nonhost_nearest:
            keep_mask.append(imin in h)
        else:
            keep_mask.append(True)

    if opts.drop_if_nonhost_nearest:
        pts = pts[keep_mask]
        dists = [d for d, k in zip(dists, keep_mask) if k]
        labels = [lab for lab, k in zip(labels, keep_mask) if k]
        hosts = [h for h, k in zip(hosts, keep_mask) if k]

    # print table
    print(f"Loaded: {xyz}")
    print(f"Atoms: {len(mol.apos)}  Samples: {len(pts)}")
    for i, (p, lab, h, d) in enumerate(zip(pts, labels, hosts, dists)):
        print(f"[{i:4d}] {lab:32s} host={h!s:12s}  d_near={d:7.4f}  p=({p[0]:9.4f},{p[1]:9.4f},{p[2]:9.4f})")

    out = save_augmented(mol, pts, labels, hosts, out_xyz=opts.out_xyz, out_mol2=opts.out_mol2)
    print(f"Wrote: {opts.out_xyz}")
    print(f"Wrote: {opts.out_mol2}")

    plot_3d(mol, pts, labels, no_plot=opts.noPlot)


if __name__ == '__main__':
    main()
