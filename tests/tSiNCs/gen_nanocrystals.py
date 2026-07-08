#!/usr/bin/env python3
"""Nanocrystal generator — Python CLI. Sphere cut: native; Miller planes + bridges: node tests/tSiNCs/gen_nanocrystals.mjs."""
import argparse
import os
import subprocess
import sys

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, REPO)

from pyBall.nanocrystal_gen import build_spherical_nanoparticle, save_xyz, find_cap_hh_pairs, save_hh_bonds_json


def run_node_planes(args):
    cmd = ['node', os.path.join(TEST_DIR, 'gen_nanocrystals.mjs'),
           '--samples', '1', '--maxFiles', '1',
           '--cutMode', 'planes',
           '--cif', args.cif,
           '--caps', args.caps,
           '--outDir', args.out_dir,
           '--prefix', args.prefix,
           '--bondFactor', str(args.bond_factor),
           '--defaultRcut', str(args.default_rcut),
           '--minHeavyDegree', str(args.min_heavy_degree),
           '--pruneMaxIter', str(args.prune_max_iter),
           '--outwardBias', str(args.outward_bias),
           '--insertProb', str(args.insert_prob),
           '--collapseProb', str(args.collapse_prob),
           '--nx-range', args.nx_range,
           '--ny-range', args.ny_range,
           '--nz-range', args.nz_range,
           '--centered', '1' if args.centered else '0',
           '--planeTemplates', args.plane_templates,
           '--planeSymC', str(args.plane_sym_c),
           '--planeCScale', str(args.plane_c_scale),
           '--planeCJitter', str(args.plane_c_jitter),
           ]
    if args.cap_hh_bonds:
        cmd.extend(['--capHHBonds', '1', '--capHHBondDist', str(args.cap_hh_bond_dist)])
    if args.seed:
        cmd.extend(['--seed', str(args.seed)])
    print('[gen_nanocrystals.py] node', ' '.join(cmd[2:]))
    subprocess.run(cmd, cwd=REPO, check=True)


def run_sphere_native(args):
    heavy_z = 14 if args.element.lower() == 'si' else 6
    prim = args.primitive_xyz
    if not prim:
        prim = os.path.join(REPO, 'cpp', 'common_resources', 'crystals',
                            'Si_primitive.xyz' if heavy_z == 14 else 'diamond_primitive.xyz')
    elems, apos = build_spherical_nanoparticle(
        prim_xyz=prim, R=args.sphere_r, nrep=args.sphere_nrep,
        min_degree=args.min_heavy_degree, heavy_z=heavy_z,
        rcut_heavy=args.rcut_heavy, outward_bias=args.outward_bias,
        resolve_clashes=(args.caps == 'H'),
    )
    os.makedirs(args.out_dir, exist_ok=True)
    tag = f'{args.prefix}_sphere_R{args.sphere_r:.1f}_nrep{args.sphere_nrep}_nat{len(elems)}'
    xyz_path = os.path.join(args.out_dir, tag + '.xyz')
    save_xyz(xyz_path, elems, apos, comment=f'sphere R={args.sphere_r} nrep={args.sphere_nrep} element={args.element}')
    print(f'[gen_nanocrystals.py] wrote {xyz_path} atoms={len(elems)}')
    if args.cap_hh_bonds:
        pairs = find_cap_hh_pairs(apos, elems, args.cap_hh_bond_dist)
        bonds_path = xyz_path + '.hh_bonds.json'
        save_hh_bonds_json(bonds_path, pairs)
        print(f'[gen_nanocrystals.py] wrote {bonds_path} hhBonds={len(pairs)}')
    if args.insert_prob > 0 or args.collapse_prob > 0:
        print('[gen_nanocrystals.py] WARNING: bridge insert/collapse requires --cutMode planes (node engine)')
    return xyz_path


def main():
    ap = argparse.ArgumentParser(description='Nanocrystal generator (Python)')
    ap.add_argument('--cutMode', choices=['sphere', 'planes'], default='sphere')
    ap.add_argument('--element', choices=['C', 'Si', 'c', 'si'], default='C')
    ap.add_argument('--cif', default=os.path.join(REPO, 'cpp/common_resources/crystals/Si-sym.cif'))
    ap.add_argument('--primitive-xyz', default=None)
    ap.add_argument('--sphere-r', type=float, default=6.0, dest='sphere_r')
    ap.add_argument('--sphere-nrep', type=int, default=5, dest='sphere_nrep')
    ap.add_argument('--rcut-heavy', type=float, default=1.75, dest='rcut_heavy')
    ap.add_argument('--caps', default='H')
    ap.add_argument('--outDir', dest='out_dir', default=os.path.join(TEST_DIR, 'OUT_nanocrystals_py'))
    ap.add_argument('--prefix', default='nc_py')
    ap.add_argument('--bondFactor', type=float, default=1.30, dest='bond_factor')
    ap.add_argument('--defaultRcut', type=float, default=1.60, dest='default_rcut')
    ap.add_argument('--minHeavyDegree', type=int, default=2, dest='min_heavy_degree')
    ap.add_argument('--pruneMaxIter', type=int, default=10, dest='prune_max_iter')
    ap.add_argument('--outwardBias', type=float, default=0.35, dest='outward_bias')
    ap.add_argument('--insertProb', type=float, default=0.0, dest='insert_prob')
    ap.add_argument('--collapseProb', type=float, default=0.0, dest='collapse_prob')
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--nx-range', default='2,2', dest='nx_range')
    ap.add_argument('--ny-range', default='2,2', dest='ny_range')
    ap.add_argument('--nz-range', default='2,2', dest='nz_range')
    ap.add_argument('--centered', action='store_true', default=True)
    ap.add_argument('--planeTemplates', default='a111', dest='plane_templates')
    ap.add_argument('--planeSymC', type=float, default=6.0, dest='plane_sym_c')
    ap.add_argument('--planeCScale', type=float, default=0.40, dest='plane_c_scale')
    ap.add_argument('--planeCJitter', type=float, default=0.0, dest='plane_c_jitter')
    ap.add_argument('--capHHBonds', action='store_true', dest='cap_hh_bonds')
    ap.add_argument('--capHHBondDist', type=float, default=1.8, dest='cap_hh_bond_dist')
    args = ap.parse_args()

    if args.cutMode == 'planes':
        run_node_planes(args)
    else:
        run_sphere_native(args)


if __name__ == '__main__':
    main()
