#!/usr/bin/env python3
"""Scan energy vs hydrogen bond length L_Hb for two-ribbon system.

Computes energy for both rigid (single-point) and relaxed (geometry optimization)
calculations as a function of L_Hb.
"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import argparse

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))
from pyBall import dftb_utils as dftbu
from pyBall import plotUtils
from doc.Topics.Kekule_Topology.GrapheneRibbonBuilder import build_two_ribbon_cell

def run_dftb_calculation(apos, enames, lvs, basis_path, do_relax=False, nk_x=16, workdir=None,
                         Temperature=300, MixingParameter=0.1, MaxScc=500, SCCTolerance=1e-4, allow_unconverged_energy=False):
    """Run DFTB+ calculation (rigid or relaxed) with convergence parameters."""
    try:
        E, apos_final, forces = dftbu.run_pbc(apos, enames, lvs, basis_path=basis_path, do_relax=do_relax, nk=(nk_x, 1, 1), dftb_exe=dftbu.DFTB_EXE, workdir=workdir, Temperature=Temperature, MixingParameter=MixingParameter, MaxScc=MaxScc, SCCTolerance=SCCTolerance, allow_unconverged_energy=allow_unconverged_energy)
        return E, apos_final, enames
    except Exception as e:
        print(f"  ERROR: {e}")
        return None, None, None

def main():
    parser = argparse.ArgumentParser(description='Scan energy vs L_Hb for two-ribbon system')
    parser.add_argument('--width', type=int, default=4, help='Ribbon width in chains (default: 4)')
    parser.add_argument('--Lx', type=float, default=2.4, help='Lattice constant Lx (default: 2.4)')
    parser.add_argument('--shift_x', type=float, default=0.0, help='x-shift for second ribbon in factors of Lx (default: 0.0)')
    parser.add_argument('--LHb_min', type=float, default=1.5, help='Min L_Hb (default: 1.5)')
    parser.add_argument('--LHb_max', type=float, default=2.5, help='Max L_Hb (default: 2.5)')
    parser.add_argument('--n_steps', type=int, default=21, help='Number of L_Hb points (default: 21 => step=0.05 in 1.5-2.5)')
    parser.add_argument('--nk', type=int, default=16, help='k-points along x (default: 16)')
    parser.add_argument('--Temperature', type=float, default=300.0, help='Fermi filling temperature [K] (default: 300)')
    parser.add_argument('--MixingParameter', type=float, default=0.1, help='Broyden mixing parameter (default: 0.1)')
    parser.add_argument('--MaxScc', type=int, default=500, help='Max SCC iterations (default: 500)')
    parser.add_argument('--SCCTolerance', type=float, default=1e-4, help='SCC tolerance (default: 1e-4)')
    parser.add_argument('--allow_unconverged', action='store_true', help='If SCC does not converge, still store last SCC electronic energy (prints loud warning)')
    parser.add_argument('--skip_rigid', action='store_true', help='Skip rigid scan')
    parser.add_argument('--skip_relax', action='store_true', help='Skip relaxed scan')
    args = parser.parse_args()
    
    # Parameters
    width_chains = args.width
    length_cells = 1
    Lx = args.Lx
    shift_x = args.shift_x
    a_CC = 1.42
    basis_path = dftbu.DEFAULT_SK_PATH
    nk_x = args.nk
    Temperature = float(args.Temperature)
    MixingParameter = float(args.MixingParameter)
    MaxScc = int(args.MaxScc)
    SCCTolerance = float(args.SCCTolerance)
    allow_unconverged_energy = bool(args.allow_unconverged)
    
    L_Hb_min = args.LHb_min
    L_Hb_max = args.LHb_max
    assert args.n_steps >= 2, "n_steps must be >= 2"
    L_Hb_values = np.linspace(L_Hb_min, L_Hb_max, args.n_steps)
    
    # Results storage
    results_rigid = []
    results_relax = []
    
    print("=" * 60)
    print("L_Hb scan for two-ribbon system")
    print("=" * 60)
    print(f"L_Hb range: {L_Hb_min:.2f} to {L_Hb_max:.2f} Å ({len(L_Hb_values)} points)")
    print(f"Width: {width_chains}, Length: {length_cells}, Lx: {Lx} Å")
    print("=" * 60)
    
    if not args.skip_rigid:
        print("\n--- Rigid scan (single-point) ---")
        for i, L_Hb in enumerate(L_Hb_values):
            print(f"L_Hb={L_Hb:.3f} Å ({i+1}/{len(L_Hb_values)})...", end=' ')
            apos, atypes, elems, lvs = build_two_ribbon_cell(width_chains=width_chains, length_cells=length_cells, Lx=Lx, a_CC=a_CC, L_Hb=L_Hb, shift_x=shift_x)
            wdir = f"temp_LHb_{os.getpid()}_rigid_{i:03d}"
            E, apos_final, enames = run_dftb_calculation(apos, elems, lvs, basis_path, do_relax=False, nk_x=nk_x, workdir=wdir, Temperature=Temperature, MixingParameter=MixingParameter, MaxScc=MaxScc, SCCTolerance=SCCTolerance, allow_unconverged_energy=allow_unconverged_energy)
            if E is not None:
                print(f"E={E:.4f} eV")
                results_rigid.append({'L_Hb': float(L_Hb), 'E': float(E), 'apos': apos_final, 'enames': enames})
            else:
                print("FAILED")
    
    if not args.skip_relax:
        print("\n--- Relaxed scan (geometry optimization) ---")
        for i, L_Hb in enumerate(L_Hb_values):
            print(f"L_Hb={L_Hb:.3f} Å ({i+1}/{len(L_Hb_values)})...", end=' ')
            apos, atypes, elems, lvs = build_two_ribbon_cell(width_chains=width_chains, length_cells=length_cells, Lx=Lx, a_CC=a_CC, L_Hb=L_Hb, shift_x=shift_x)
            wdir = f"temp_LHb_{os.getpid()}_relax_{i:03d}"
            E, apos_final, enames = run_dftb_calculation(apos, elems, lvs, basis_path, do_relax=True, nk_x=nk_x, workdir=wdir, Temperature=Temperature, MixingParameter=MixingParameter, MaxScc=MaxScc, SCCTolerance=SCCTolerance, allow_unconverged_energy=allow_unconverged_energy)
            if E is not None:
                print(f"E={E:.4f} eV")
                results_relax.append({'L_Hb': float(L_Hb), 'E': float(E), 'apos': apos_final, 'enames': enames})
            else:
                print("FAILED")
    
    # Plot results
    if results_rigid:
        L_Hb_rigid = [r['L_Hb'] for r in results_rigid]
        E_rigid = [r['E'] for r in results_rigid]
    
    if results_relax:
        L_Hb_relax = [r['L_Hb'] for r in results_relax]
        E_relax = [r['E'] for r in results_relax]
    
    fig, ax = plt.subplots(figsize=(10, 6))
    if results_rigid:
        plotUtils.plot_scan_profile(L_Hb_rigid, E_rigid, xlabel='H-bond length L_Hb [Å]', ylabel='Energy [eV]', title=f'Energy vs L_Hb (w{width_chains}, Lx={Lx} Å)', label='Rigid (single-point)', color='blue', marker='o', ax=ax)
    if results_relax:
        plotUtils.plot_scan_profile(L_Hb_relax, E_relax, xlabel='H-bond length L_Hb [Å]', ylabel='Energy [eV]', title=f'Energy vs L_Hb (w{width_chains}, Lx={Lx} Å)', label='Relaxed (optimized)', color='red', marker='s', ax=ax)
    fname_plot = f'LHb_scan_w{width_chains}_Lx{Lx}.png'
    fig.savefig(fname_plot, dpi=150, bbox_inches='tight')
    print(f"\nSaved plot: {fname_plot}")
    plt.close(fig)
    
    # Save results to file
    fname_log = f'LHb_scan_w{width_chains}_Lx{Lx}.log'
    with open(fname_log, 'w') as f:
        f.write(f"# L_Hb [Å]  E_rigid [eV]  E_relax [eV]  (width={width_chains}, Lx={Lx})\n")
        for i, L_Hb in enumerate(L_Hb_values):
            E_rigid = results_rigid[i]['E'] if i < len(results_rigid) else 'nan'
            E_relax = results_relax[i]['E'] if i < len(results_relax) else 'nan'
            f.write(f"{L_Hb:.4f}  {E_rigid}  {E_relax}\n")
    print(f"Saved log: {fname_log}")
    
    # Save XYZ movies
    fname_xyz_rigid = f'LHb_scan_w{width_chains}_Lx{Lx}'
    fname_xyz_relax = f'LHb_scan_w{width_chains}_Lx{Lx}'
    if shift_x != 0.0:
        fname_xyz_rigid += f'_sx{shift_x:.2f}'
        fname_xyz_relax += f'_sx{shift_x:.2f}'
    fname_xyz_rigid += '_rigid.xyz'
    fname_xyz_relax += '_relax.xyz'
    if results_rigid:
        dftbu.save_xyz_movie(results_rigid, fname_xyz_rigid, lvs=lvs, label='Rigid', key_order=['L_Hb', 'E'])
    if results_relax:
        dftbu.save_xyz_movie(results_relax, fname_xyz_relax, lvs=lvs, label='Relaxed', key_order=['L_Hb', 'E'])
    
    print("\nDone!")

if __name__ == '__main__':
    main()
