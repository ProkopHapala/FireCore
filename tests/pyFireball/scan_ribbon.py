#!/usr/bin/env python3
"""Universal lattice scan script for zigzag graphene ribbons.

Usage examples:
  python scan_ribbon.py --passivation NH --width 6
  python scan_ribbon.py --passivation CH --width 4 --do_relax --reverse
  python scan_ribbon.py --passivation C=O --Lx_min 2.0 --Lx_max 3.0 --nLx 16 --do_relax

Passivation types: N, NH, CH, C=O, C-OH (see GrapheneRibbonBuilder)
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import argparse

sys.path.append("../../")
from pyBall import FireCore as fc

ELEM_MAP = {'H': 1, 'C': 6, 'N': 7, 'O': 8}


def make_ribbon_builder(passivation, width_chains, length_cells=1, a_CC=1.42):
    """Return a geometry_builder callable(Lx) -> (pos2d, atypes) for a zigzag ribbon."""
    from doc.Topics.Kekule_Topology.GrapheneRibbonBuilder import GrapheneRibbonBuilder
    xa_nom = a_CC * np.cos(np.pi / 6)

    def builder(Lx):
        b = GrapheneRibbonBuilder(a_CC=a_CC)
        scale_x = Lx / (2.0 * xa_nom)
        pos2d, elems, bonds = b.build_zigzag_ribbon(
            width_chains=width_chains, length_cells=length_cells,
            passivation=passivation, scale_x=scale_x)
        atypes = np.array([ELEM_MAP[e] for e in elems], dtype=np.int32)
        return np.array(pos2d), atypes

    return builder


def main():
    parser = argparse.ArgumentParser(description='Universal lattice scan for zigzag graphene ribbons')
    parser.add_argument('--passivation', type=str, default='NH',   help='Passivation type: N, NH, CH, C=O, C-OH')
    parser.add_argument('--width',       type=int, default=6,      help='width_chains (ribbon thickness in y)')
    parser.add_argument('--length',      type=int, default=1,      help='length_cells (unit cells along x)')
    parser.add_argument('--a_CC',        type=float, default=1.42, help='C-C bond length (Ang)')
    parser.add_argument('--nk',          type=int, default=16,     help='Number of k-points')
    parser.add_argument('--nmax_scf',    type=int, default=200,    help='Max SCF iterations')
    parser.add_argument('--do_relax',    action='store_true',      help='Use FIRE relaxation')
    parser.add_argument('--nstep_relax', type=int, default=100,    help='Max FIRE ionic steps')
    parser.add_argument('--Lx_min',      type=float, default=2.0,  help='Min Lx (Ang)')
    parser.add_argument('--Lx_max',      type=float, default=3.0,  help='Max Lx (Ang)')
    parser.add_argument('--nLx',         type=int, default=16,     help='Number of Lx points')
    parser.add_argument('--reverse',     action='store_true',      help='Scan from Lx_max down to Lx_min')
    parser.add_argument('--Ly',          type=float, default=20.0, help='Cell size y (vacuum)')
    parser.add_argument('--Lz',          type=float, default=20.0, help='Cell size z (vacuum)')
    args = parser.parse_args()

    passivation = args.passivation
    label = passivation  # used for file naming
    label_safe = label.replace('=', '').replace('-', '')
    geom_label = f"w{args.width}_l{args.length}"  # e.g., "w6_l1"
    mode = 'relax' if args.do_relax else 'fixed'

    if args.reverse:
        Lx_vals = np.linspace(args.Lx_max, args.Lx_min, args.nLx)
    else:
        Lx_vals = np.linspace(args.Lx_min, args.Lx_max, args.nLx)

    print(f"Ribbon lattice scan: passivation={passivation}, width_chains={args.width}, length_cells={args.length}")
    print(f"  nk={args.nk}, nmax_scf={args.nmax_scf}, do_relax={args.do_relax}, nstep_relax={args.nstep_relax}")
    print(f"  Lx range: {Lx_vals[0]:.3f} to {Lx_vals[-1]:.3f} ({args.nLx} points)")

    # Build geometry_builder callable
    geometry_builder = make_ribbon_builder(passivation, args.width, args.length, args.a_CC)

    # Run scan
    res, _ = fc.run_lattice_scan(
        label, geometry_builder, Lx_vals,
        nk=args.nk, nmax_scf=args.nmax_scf,
        do_relax=args.do_relax, nstep_relax=args.nstep_relax,
        use_continuous_path=args.do_relax,
        Ly=args.Ly, Lz=args.Lz, geom_label=geom_label)

    # Sort results by Lx for plotting (in case of reverse scan)
    sort_idx = np.argsort(res[:, 0])
    res = res[sort_idx]

    # Save npz
    out_npz = f"lattice_scan_{label_safe}_{geom_label}_{mode}.npz"
    np.savez(out_npz, res=res, Lx_vals=Lx_vals, nk=args.nk, nmax_scf=args.nmax_scf,
             passivation=passivation, width=args.width, length=args.length)

    # Plot
    plt.figure(figsize=(8, 6))
    valid = ~np.isnan(res[:, 1])
    plt.plot(res[valid, 0], res[valid, 1], 'rs-', lw=2, markersize=6, label=passivation)
    plt.xlabel("Lattice Constant Lx (Å)")
    plt.ylabel("Total Energy (eV)")
    plt.grid(True, alpha=0.3)
    if np.any(valid):
        idx = np.argmin(res[valid, 1])
        lx_opt = res[valid, 0][idx];  e_opt = res[valid, 1][idx]
        plt.axvline(lx_opt, color='r', ls='--', alpha=0.5)
        plt.title(f"Zigzag Ribbon [{passivation}] {geom_label}\nLx_opt={lx_opt:.3f} Å  E_min={e_opt:.4f} eV")
        print(f"Result: Lx_opt={lx_opt:.4f} Å, E_min={e_opt:.6f} eV")
    plt.tight_layout()
    out_png = f"lattice_scan_{label_safe}_{geom_label}_{mode}.png"
    plt.savefig(out_png, dpi=150)

    # Write log
    out_log = f"lattice_scan_{label_safe}_{geom_label}_{mode}.log"
    with open(out_log, 'w') as f:
        f.write(f"# Lattice scan: passivation={passivation} {geom_label}\n")
        f.write(f"# nk={args.nk}, nmax_scf={args.nmax_scf}, do_relax={args.do_relax}\n")
        f.write(f"# Lx range: {args.Lx_min} to {args.Lx_max} ({args.nLx} points)\n")
        f.write(f"# {'Lx':>8}  {'E_tot':>14}  {'E_bs':>14}  {'Fmax':>12}\n")
        for row in res:
            f.write(f"  {row[0]:8.4f}  {row[1]:14.6f}  {row[2]:14.6f}  {row[3]:12.6f}\n")

    print(f"\nOutput files:")
    print(f"  {out_png}")
    print(f"  {out_log}")
    print(f"  lattice_scan_{label_safe}_{geom_label}_{mode}.xyz")
    print(f"  {out_npz}")


if __name__ == "__main__":
    main()
