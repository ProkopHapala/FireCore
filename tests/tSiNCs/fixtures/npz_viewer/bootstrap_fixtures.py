#!/usr/bin/env python3
"""Bootstrap NPZ viewer fixtures (Agent B, before Agent C delivers full set)."""
import os
import sys
import numpy as np

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..'))
sys.path.insert(0, _ROOT)

from pyBall.AtomicSystem import AtomicSystem

OUT = os.path.dirname(os.path.abspath(__file__))

def _write_crystal_from_xyz(xyz_path, out_npz):
    sys = AtomicSystem(fname=xyz_path, bPreinit=True)
    pos = np.asarray(sys.apos, dtype=np.float64)
    Z = np.asarray(sys.atypes, dtype=np.int32)
    bonds = np.asarray(sys.bonds, dtype=np.int32) if sys.bonds is not None else np.zeros((0, 2), dtype=np.int32)
    natoms = np.int32(pos.shape[0])
    np.savez_compressed(out_npz, pos=pos, Z=Z, natoms=natoms, bonds_ij=bonds)
    return pos, Z, bonds

def _write_topology(pos, Z, bonds, out_npz, n_groups=3):
    N = pos.shape[0]
    M = 48
    neigh_idx = np.full((N, M), -1, dtype=np.int32)
    bond_l0 = np.zeros((N, M), dtype=np.float32)
    bond_k = np.zeros((N, M), dtype=np.float32)
    stick_class = np.zeros((N, M), dtype=np.int32)
    for bi, (i, j) in enumerate(bonds):
        if bi >= M:
            break
        neigh_idx[i, bi] = j; neigh_idx[j, bi] = i
        d = float(np.linalg.norm(pos[i] - pos[j]))
        bond_l0[i, bi] = d; bond_l0[j, bi] = d
        bond_k[i, bi] = 200.0 + 50.0 * bi; bond_k[j, bi] = bond_k[i, bi]
        stick_class[i, bi] = 1; stick_class[j, bi] = 1
    icolGroup = np.full(N, -1, dtype=np.int32)
    chunk = max(N // n_groups, 1)
    for g in range(n_groups):
        sl = slice(g * chunk, min((g + 1) * chunk, N))
        icolGroup[sl] = g
    mins = []; maxs = []
    for g in range(n_groups):
        m = icolGroup == g
        if not np.any(m):
            mins.append([0, 0, 0]); maxs.append([1, 1, 1]); continue
        p = pos[m]
        mins.append(p.min(axis=0)); maxs.append(p.max(axis=0))
    group_bbox_min = np.asarray(mins, dtype=np.float64)
    group_bbox_max = np.asarray(maxs, dtype=np.float64)
    radius = np.full(N, 1.5, dtype=np.float64)
    np.savez_compressed(
        out_npz,
        pos=pos, Z=Z,
        neigh_idx=neigh_idx, bond_l0=bond_l0, bond_k=bond_k, stick_class=stick_class,
        max_neighbors=np.int32(M), n_bond=np.int32(bonds.shape[0]),
        icolGroup=icolGroup, icol=icolGroup.copy(),
        group_bbox_min=group_bbox_min, group_bbox_max=group_bbox_max,
        n_groups=np.int32(n_groups), schema_version=np.int32(1),
        radius=radius,
    )

def _write_n500(out_npz):
    rng = np.random.default_rng(42)
    N = 500
    pos = rng.normal(size=(N, 3)).astype(np.float64) * 2.0
    Z = np.full(N, 14, dtype=np.int32)
    Z[::10] = 1
    bonds = []
    for i in range(N - 1):
        bonds.append([i, i + 1])
    bonds = np.asarray(bonds, dtype=np.int32)
    natoms = np.int32(N)
    np.savez_compressed(out_npz, pos=pos, Z=Z, natoms=natoms, bonds_ij=bonds)

def main():
    os.makedirs(OUT, exist_ok=True)
    xyz_src = os.path.join(_ROOT, 'cpp', 'common_resources', 'xyz', 'Si10_H.xyz')
    if not os.path.isfile(xyz_src):
        raise FileNotFoundError(xyz_src)
    pos, Z, bonds = _write_crystal_from_xyz(xyz_src, os.path.join(OUT, '01_init.npz'))
    _write_topology(pos, Z, bonds, os.path.join(OUT, '03_topology.npz'))
    _write_n500(os.path.join(OUT, 'bench_500.npz'))
    npy_dir = os.path.join(OUT, 'npy_fast')
    os.makedirs(npy_dir, exist_ok=True)
    np.save(os.path.join(npy_dir, 'pos.npy'), pos)
    np.save(os.path.join(npy_dir, 'Z.npy'), Z)
    np.save(os.path.join(npy_dir, 'bonds_ij.npy'), bonds)
    import shutil
    for name in ('CH4.xyz',):
        src = os.path.join(_ROOT, 'cpp', 'common_resources', 'xyz', name)
        if os.path.isfile(src):
            shutil.copy2(src, os.path.join(OUT, name))
    print(f'bootstrap_fixtures: wrote fixtures to {OUT}')

if __name__ == '__main__':
    main()
