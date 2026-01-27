#!/usr/bin/env python3
import argparse
import numpy as np
import os
from pyBall import AtomicSystem
from pyBall.OCL.MMFFL import MMFFL
from pyBall.XPDB_AVBD.test_TiledJacobi_molecules import _bonds_to_adj, _pack_fixed_bonds, dump_xpdb_inputs_text


def build_topology_only(fname, args):
    mol = AtomicSystem.AtomicSystem(fname=fname)
    if mol.bonds is None or len(mol.bonds) == 0:
        mol.findBonds()
    mol.neighs()

    ff = MMFFL(
        L_pi=args.L_pi,
        two_pi_dummies=bool(args.two_pi),
        Kang=(args.k_angle if args.k_angle is not None else args.bond_k),
        Kpi_host=args.k_pi,
        Kpi_orth=args.k_pi_orth,
        Kpi_align=args.k_pi_align,
        L_epair=args.L_epair,
        Kep_host=args.k_ep,
        Kep_orth=args.k_ep_orth,
        Kep_pair=args.k_ep_pair,
        align_pi_vectors=bool(args.align_pi_vectors),
    )

    topo = ff.build_topology(
        mol,
        type_source=args.type_source,
        add_angle=bool(args.enable_angles),
        add_pi=bool(args.add_pi),
        two_pi_dummies=bool(args.two_pi),
        add_pi_align=bool(args.add_pi_align),
        add_epair=bool(args.add_epair),
        add_epair_pairs=bool(args.add_epair_pairs),
    )

    # gather arrays
    apos_all = np.asarray(ff.apos[:ff.natoms], dtype=np.float32)
    bonds_primary = topo.get("bonds_primary", [])
    bonds_derived = topo.get("bonds_angle", []) + topo.get("bonds_pi", []) + topo.get("bonds_pi_align", []) + topo.get("bonds_epair", []) + topo.get("bonds_epair_pair", [])
    bonds_all = sorted(set(bonds_primary + bonds_derived))
    n_all = apos_all.shape[0]
    bonds_adj = _bonds_to_adj(
        n_all,
        bonds_all,
        default_L=None,
        default_K=args.bond_k,
        apos=apos_all,
    )
    bond_indices, bond_lengths, bond_stiffness, bond_type_mask = _pack_fixed_bonds(
        n_all,
        bonds_adj,
        n_max_bonded=args.max_bonded,
    )

    radius = np.full(n_all, args.atom_rad, dtype=np.float32)
    mass = np.full(n_all, args.atom_mass, dtype=np.float32)

    dump_xpdb_inputs_text(
        args.dump_inputs,
        xyz_path=fname,
        pos0=apos_all[:, :3],
        radius=radius,
        mass=mass,
        bond_indices=bond_indices,
        bond_lengths=bond_lengths,
        bond_stiffness=bond_stiffness,
        fixed=args.dump_fixed,
    )
    print(f"[DUMP] wrote {args.dump_inputs}")


def main():
    p = argparse.ArgumentParser(description="Headless XPDB input dump (no matplotlib)")
    p.add_argument("--molecule", type=str, required=True, help="Path to molecule (.xyz/.mol/.mol2)")
    p.add_argument("--type_source", type=str, default="table", help="MMFFL type source")
    p.add_argument("--enable_angles", type=int, default=1)
    p.add_argument("--add_pi", type=int, default=1)
    p.add_argument("--two_pi", type=int, default=1)
    p.add_argument("--add_pi_align", type=int, default=0)
    p.add_argument("--add_epair", type=int, default=1)
    p.add_argument("--add_epair_pairs", type=int, default=0)
    p.add_argument("--align_pi_vectors", type=int, default=0)
    p.add_argument("--L_pi", type=float, default=1.0)
    p.add_argument("--L_epair", type=float, default=0.5)
    p.add_argument("--k_angle", type=float, default=None)
    p.add_argument("--k_pi", type=float, default=50.0)
    p.add_argument("--k_pi_orth", type=float, default=30.0)
    p.add_argument("--k_pi_align", type=float, default=15.0)
    p.add_argument("--k_ep", type=float, default=40.0)
    p.add_argument("--k_ep_orth", type=float, default=25.0)
    p.add_argument("--k_ep_pair", type=float, default=10.0)
    p.add_argument("--bond_k", type=float, default=200.0)
    p.add_argument("--bond_len", type=float, default=1.3)
    p.add_argument("--atom_rad", type=float, default=0.2)
    p.add_argument("--atom_mass", type=float, default=1.0)
    p.add_argument("--max_bonded", type=int, default=16)
    p.add_argument("--dump_fixed", type=int, default=6, help="Pad adjacency to this size; matches JS nMaxBonded")
    p.add_argument("--dump_inputs", type=str, required=True, help="Output path for line-by-line XPDB inputs")
    args = p.parse_args()

    build_topology_only(args.molecule, args)


if __name__ == "__main__":
    main()
