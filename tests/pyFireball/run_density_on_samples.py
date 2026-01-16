#!/usr/bin/env python3
"""
Generate sampling points using sample_density_points.py, then evaluate Fireball density at those points.
Outputs augmented XYZ with density as 5th column for each input molecule.
"""
import os
import sys
import argparse
import numpy as np

sys.path.append(os.path.join(os.path.dirname(__file__), "../.."))
from pyBall.AtomicSystem import AtomicSystem
from pyBall import FireCore as fc

# We will import the sampling builder from sample_density_points by execfile-like approach
# to avoid duplicating logic; alternatively we could refactor to a module.
from sample_density_points import build_samples, _parse_floats_csv


def save_xyz_with_density(enames, apos, density, fname):
    """Save XYZ with density in 5th column."""
    n = len(enames)
    with open(fname, "w") as f:
        f.write(f"{n}\n")
        f.write("samples with density\n")
        for e, p, d in zip(enames, apos, density):
            f.write(f"{e:<24s}{p[0]:14.6f}{p[1]:14.6f}{p[2]:14.6f}{d:16.6e}\n")


def process_one(xyz_path, opts):
    # Remember current working directory
    cwd0 = os.getcwd()

    # Load molecule
    mol = AtomicSystem(fname=xyz_path)
    mol.findBonds(RvdwCut=opts.bond_rvdwcut)
    mol.neighs()

    # Build sample points using existing builder
    pts, labels, hosts, _ = build_samples(mol, opts)

    # Initialize FireCore per molecule; ensure cwd contains Fdata (usually parent of molecules/)
    workdir = os.path.abspath(os.path.join(os.path.dirname(xyz_path), '..'))
    try:
        os.chdir(workdir)
        fc.setVerbosity(0)
        fc.initialize(atomType=mol.atypes, atomPos=mol.apos, verbosity=0)
        fc.evalForce(mol.apos, nmax_scf=opts.nmax_scf)
    finally:
        os.chdir(cwd0)

    # Compute density
    dens = fc.dens2points(pts, f_den=opts.f_den, f_den0=opts.f_den0)

    # Compose output enames/apos (only samples, not original atoms)
    out_enames = []
    out_apos = []
    for lab, p in zip(labels, pts):
        out_enames.append(lab)
        out_apos.append(p)
    out_apos = np.array(out_apos)

    out_basename = os.path.splitext(os.path.basename(xyz_path))[0]
    out_xyz = os.path.join(opts.out_dir, f"{out_basename}_dens.xyz")
    os.makedirs(opts.out_dir, exist_ok=True)
    save_xyz_with_density(out_enames, out_apos, dens, out_xyz)
    print(f"Wrote {out_xyz} with {len(out_enames)} samples")

default_molecules = [ 'H2O.xyz', 'NH3.xyz', 'CH4.xyz' ]

def main():
    ap = argparse.ArgumentParser(description="Generate sampling points and evaluate Fireball density at those points")
    ap.add_argument('--xyz-list', nargs='+', default=None, help='List of input xyz files')

    # sampling options forwarded to build_samples
    ap.add_argument('--atoms',         type=int, default=1)
    ap.add_argument('--bonds',         type=int, default=1)
    ap.add_argument('--epairs',        type=int, default=1)
    ap.add_argument('--epair-dists',   default='0.5,1.0')
    ap.add_argument('--invert-bonds',  type=int, default=1)
    ap.add_argument('--angles',        type=int, default=1)
    ap.add_argument('--apex',          type=int, default=1)
    ap.add_argument('--invert-far',    type=int, default=1)
    ap.add_argument('--far-epairs',    type=int, default=1)
    ap.add_argument('--invert-epairs', type=int, default=1)
    ap.add_argument('--apex-dist',     type=float, default=0.8)
    ap.add_argument('--bond-fracs',    default='0.5')
    ap.add_argument('--bond-rvdwcut',  type=float, default=1.5)
    ap.add_argument('--angle-dist',    type=float, default=0.7)
    ap.add_argument('--apex-mode',     choices=['avg', 'per-bond'], default='avg')
    ap.add_argument('--far-dist',      type=float, default=6.0)
    ap.add_argument('--dedup',         type=float, default=1e-6)
    ap.add_argument('--drop-if-nonhost-nearest', type=int, default=0)

    # Fireball eval options
    ap.add_argument('--nmax-scf', type=int, default=100)
    ap.add_argument('--f-den',    type=float, default=1.0)
    ap.add_argument('--f-den0',   type=float, default=0.0)

    ap.add_argument('--out-dir',  default='dens_outputs', help='Directory to write *_dens.xyz files')
    ap.add_argument('--spawn-per-mol', type=int, default=1, help='If 1, run each xyz in a fresh python process to avoid Fireball reinit issues')

    opts = ap.parse_args()

    # parse comma lists
    opts.bond_fracs = _parse_floats_csv(opts.bond_fracs)
    opts.epair_dists = _parse_floats_csv(opts.epair_dists)

    if opts.xyz_list is None:
        # Default list used when CLI not provided
        base = os.path.join(os.path.dirname(__file__), 'molecules')
        opts.xyz_list = [ os.path.join(base, m) for m in default_molecules ]

    # If requested, spawn new interpreter per molecule to avoid Fortran reallocation errors
    if opts.spawn_per_mol:
        import subprocess
        script_path = os.path.abspath(__file__)
        for xyz in opts.xyz_list:
            if not os.path.exists(xyz):
                print(f"[WARN] missing xyz: {xyz}")
                continue
            cmd = [sys.executable, script_path, '--xyz-list', xyz,
                   '--atoms', str(opts.atoms), '--bonds', str(opts.bonds), '--epairs', str(opts.epairs),
                   '--epair-dists', ','.join(str(x) for x in opts.epair_dists),
                   '--invert-bonds', str(opts.invert_bonds), '--angles', str(opts.angles), '--apex', str(opts.apex),
                   '--invert-far', str(opts.invert_far), '--far-epairs', str(opts.far_epairs), '--invert-epairs', str(opts.invert_epairs),
                   '--apex-dist', str(opts.apex_dist), '--bond-fracs', ','.join(str(x) for x in opts.bond_fracs),
                   '--bond-rvdwcut', str(opts.bond_rvdwcut), '--angle-dist', str(opts.angle_dist), '--apex-mode', opts.apex_mode,
                   '--far-dist', str(opts.far_dist), '--dedup', str(opts.dedup), '--drop-if-nonhost-nearest', str(opts.drop_if_nonhost_nearest),
                   '--nmax-scf', str(opts.nmax_scf), '--f-den', str(opts.f_den), '--f-den0', str(opts.f_den0),
                   '--out-dir', opts.out_dir, '--spawn-per-mol', '0']
            print(f"[spawn] {' '.join(cmd)}")
            subprocess.run(cmd, check=True)
        return

    for xyz in opts.xyz_list:
        if not os.path.exists(xyz):
            print(f"[WARN] missing xyz: {xyz}")
            continue
        print(f"Processing {xyz}")
        process_one(xyz, opts)


if __name__ == '__main__':
    main()
