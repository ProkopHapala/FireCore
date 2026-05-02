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

sys.path.append("../../")
from pyBall import dftb_utils as dftbu

# Import from build_two_ribbons
ELEM_MAP = {'H': 1, 'C': 6, 'N': 7, 'O': 8}
ELEM_MAP_INV = {v: k for k, v in ELEM_MAP.items()}

def build_ribbon(passivation, width_chains, length_cells, Lx, a_CC=1.42):
    """Build a single ribbon geometry."""
    from doc.Topics.Kekule_Topology.GrapheneRibbonBuilder import GrapheneRibbonBuilder
    xa_nom = a_CC * np.cos(np.pi / 6)
    b = GrapheneRibbonBuilder(a_CC=a_CC)
    scale_x = Lx / (2.0 * xa_nom)
    pos2d, elems, bonds = b.build_zigzag_ribbon(
        width_chains=width_chains, length_cells=length_cells,
        passivation=passivation, scale_x=scale_x)
    atypes = np.array([ELEM_MAP[e] for e in elems], dtype=np.int32)
    return np.array(pos2d), atypes, elems

def build_two_ribbon_cell(width_chains=4, length_cells=1, Lx=2.4, a_CC=1.42, L_Hb=2.0, shift_x=0.0):
    """Build two-ribbon unit cell with N and NH passivation for hydrogen bonding.
    
    Args:
        width_chains: ribbon width in chains
        length_cells: length in unit cells
        Lx: lattice constant along x
        a_CC: C-C bond length
        L_Hb: hydrogen bond length between ribbons (Å)
        shift_x: x-shift for second ribbon in factors of Lx
    """
    # Build N-passivated ribbon (bottom)
    pos2d_N, atypes_N, elems_N = build_ribbon('N', width_chains, length_cells, Lx, a_CC)
    
    # Build NH-passivated ribbon (top)
    pos2d_NH, atypes_NH, elems_NH = build_ribbon('NH', width_chains, length_cells, Lx, a_CC)
    
    # Convert to 3D
    apos_N = np.zeros((len(atypes_N), 3))
    apos_N[:, 0:2] = pos2d_N
    
    apos_NH = np.zeros((len(atypes_NH), 3))
    apos_NH[:, 0:2] = pos2d_NH
    
    # Center each ribbon in y
    apos_N[:, 1] -= np.mean(apos_N[:, 1])
    apos_NH[:, 1] -= np.mean(apos_NH[:, 1])
    
    # Apply x-shift to NH ribbon
    apos_NH[:, 0] += shift_x * Lx
    
    # Calculate ribbon spans in y
    y_min_N = np.min(apos_N[:, 1])
    y_max_N = np.max(apos_N[:, 1])
    y_span_N = y_max_N - y_min_N
    
    y_min_NH = np.min(apos_NH[:, 1])
    y_max_NH = np.max(apos_NH[:, 1])
    y_span_NH = y_max_NH - y_min_NH
    
    # Shift NH ribbon so: y_min_NH = y_max_N + L_Hb
    shift_y = (y_max_N + L_Hb) - y_min_NH
    apos_NH[:, 1] += shift_y
    
    # Recalculate extents after shift
    y_min_NH = np.min(apos_NH[:, 1])
    y_max_NH = np.max(apos_NH[:, 1])
    
    # Calculate cell size: Ly = y_span_N + y_span_NH + 2 * L_Hb
    Ly = y_span_N + y_span_NH + 2 * L_Hb
    
    # Combine ribbons
    apos = np.vstack([apos_N, apos_NH])
    atypes = np.concatenate([atypes_N, atypes_NH])
    elems = list(elems_N) + list(elems_NH)
    
    # Add z-separation (vacuum)
    z_center = 10.0
    apos[:, 2] = z_center
    
    # Lz is vacuum in z
    Lz = 20.0
    
    # Center in cell (y-direction centered on the combined structure)
    apos[:, 1] -= np.mean(apos[:, 1])
    apos[:, 2] -= np.mean(apos[:, 2])
    apos[:, 2] += 0.5 * Lz
    
    # Lattice vectors: periodic along x, vacuum along y, z
    lvs = np.array([[Lx, 0.0, 0.0], [0.0, Ly, 0.0], [0.0, 0.0, Lz]])
    
    return apos, atypes, elems, lvs

def run_dftb_calculation(apos, atypes, lvs, basis_path, do_relax=False, nk_x=16, 
                         Temperature=600, MixingParameter=0.1, MaxScc=500, SCCTolerance=1e-4):
    """Run DFTB+ calculation (rigid or relaxed) with convergence parameters."""
    # Convert atomic types to element names
    enames = [ELEM_MAP_INV[int(at)] for at in atypes]
    
    # Create temporary directory
    temp_dir = f"temp_LHb_{os.getpid()}"
    os.makedirs(temp_dir, exist_ok=True)
    cwd = os.getcwd()
    os.chdir(temp_dir)
    
    try:
        # Write periodic DFTB+ input with improved convergence parameters
        dftbu.makeDFTBjob_pbc(
            enames=enames, apos=apos, lvs=lvs, fname='dftb_in.hsd',
            basis_path=basis_path,
            nk=(nk_x, 1, 1), k_shift=(0.5, 0.0, 0.0),
            opt=do_relax, params=dftbu.default_params,
            Temperature=Temperature, MixingParameter=MixingParameter,
            MaxScc=MaxScc, SCCTolerance=SCCTolerance
        )
        
        # Run DFTB+
        dftb_path = '/home/prokophapala/miniconda3/bin/dftb+'
        os.system(f'{dftb_path} > OUT')
        
        # Parse energy
        Estr = os.popen('grep "Total Energy" OUT | tail -1 | cut -b 52-70').read().strip()
        if not Estr:
            raise ValueError("Could not parse energy from DFTB+ output")
        E = float(Estr)
        
        # Read final geometry if relaxed
        if do_relax and os.path.exists('geom.out.xyz'):
            with open('geom.out.xyz', 'r') as fxyz:
                natoms = int(fxyz.readline().strip())
                fxyz.readline()  # comment line
                apos_final = np.zeros((natoms, 3), dtype=np.float64)
                for j in range(natoms):
                    parts = fxyz.readline().split()
                    apos_final[j, 0] = float(parts[1])
                    apos_final[j, 1] = float(parts[2])
                    apos_final[j, 2] = float(parts[3])
            return E, apos_final, enames
        else:
            return E, apos, enames
            
    except Exception as e:
        print(f"  ERROR: {e}")
        return None, None, None
    finally:
        os.chdir(cwd)
        # Clean up
        # os.system(f'rm -rf {temp_dir}')

def save_xyz_movie(results, fname, label, lvs):
    """Save XYZ movie with properties in comments."""
    with open(fname, 'w') as f:
        for r in results:
            apos = r['apos']
            enames = r['enames']
            L_Hb = r['L_Hb']
            E = r['E']
            
            f.write(f"{len(enames)}\n")
            f.write(f"{label}: L_Hb={L_Hb:.4f} E={E:.6f} eV Lx={lvs[0,0]:.2f} Ly={lvs[1,1]:.2f} Lz={lvs[2,2]:.2f}\n")
            for ename, pos in zip(enames, apos):
                f.write(f"{ename} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}\n")
    print(f"Saved XYZ movie: {fname}")

def main():
    parser = argparse.ArgumentParser(description='Scan energy vs L_Hb for two-ribbon system')
    parser.add_argument('--width', type=int, default=4, help='Ribbon width in chains (default: 4)')
    parser.add_argument('--Lx', type=float, default=2.4, help='Lattice constant Lx (default: 2.4)')
    parser.add_argument('--shift_x', type=float, default=0.0, help='x-shift for second ribbon in factors of Lx (default: 0.0)')
    args = parser.parse_args()
    
    # Parameters
    width_chains = args.width
    length_cells = 1
    Lx = args.Lx
    shift_x = args.shift_x
    a_CC = 1.42
    basis_path = '/home/prokophapala/SIMULATIONS/dftbplus/slakos/3ob-3-1/'
    nk_x = 16
    
    # L_Hb scan range (fine sampling in 1.5-2.5 Å range)
    L_Hb_min = 1.5
    L_Hb_max = 2.5
    step = 0.05
    L_Hb_values = np.arange(L_Hb_min, L_Hb_max + step, step)
    
    # Results storage
    results_rigid = []
    results_relax = []
    
    print("=" * 60)
    print("L_Hb scan for two-ribbon system")
    print("=" * 60)
    print(f"L_Hb range: {L_Hb_min:.2f} to {L_Hb_max:.2f} Å ({len(L_Hb_values)} points)")
    print(f"Width: {width_chains}, Length: {length_cells}, Lx: {Lx} Å")
    print("=" * 60)
    
    # Rigid scan
    print("\n--- Rigid scan (single-point) ---")
    for i, L_Hb in enumerate(L_Hb_values):
        print(f"L_Hb={L_Hb:.3f} Å ({i+1}/{len(L_Hb_values)})...", end=' ')
        
        # Build geometry
        apos, atypes, elems, lvs = build_two_ribbon_cell(
            width_chains=width_chains, length_cells=length_cells,
            Lx=Lx, a_CC=a_CC, L_Hb=L_Hb, shift_x=shift_x)
        
        # Run rigid calculation
        E, apos_final, enames = run_dftb_calculation(apos, atypes, lvs, basis_path, do_relax=False, nk_x=nk_x)
        
        if E is not None:
            print(f"E={E:.4f} eV")
            results_rigid.append({'L_Hb': L_Hb, 'E': E, 'apos': apos_final, 'enames': enames})
        else:
            print("FAILED")
    
    # Relaxed scan
    print("\n--- Relaxed scan (geometry optimization) ---")
    for i, L_Hb in enumerate(L_Hb_values):
        print(f"L_Hb={L_Hb:.3f} Å ({i+1}/{len(L_Hb_values)})...", end=' ')
        
        # Build geometry
        apos, atypes, elems, lvs = build_two_ribbon_cell(
            width_chains=width_chains, length_cells=length_cells,
            Lx=Lx, a_CC=a_CC, L_Hb=L_Hb, shift_x=shift_x)
        
        # Run relaxed calculation
        E, apos_final, enames = run_dftb_calculation(apos, atypes, lvs, basis_path, do_relax=True, nk_x=nk_x)
        
        if E is not None:
            print(f"E={E:.4f} eV")
            results_relax.append({'L_Hb': L_Hb, 'E': E, 'apos': apos_final, 'enames': enames})
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
        ax.plot(L_Hb_rigid, E_rigid, 'o-', label='Rigid (single-point)', color='blue')
    
    if results_relax:
        ax.plot(L_Hb_relax, E_relax, 's-', label='Relaxed (optimized)', color='red')
    
    ax.set_xlabel('H-bond length L_Hb [Å]')
    ax.set_ylabel('Energy [eV]')
    ax.set_title(f'Energy vs H-bond length for two-ribbon system (w{width_chains}, Lx={Lx} Å)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    fname_plot = f'LHb_scan_w{width_chains}_Lx{Lx}.png'
    plt.savefig(fname_plot, dpi=150)
    print(f"\nSaved plot: {fname_plot}")
    plt.close()
    
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
    save_xyz_movie(results_rigid, fname_xyz_rigid, 'Rigid', lvs)
    save_xyz_movie(results_relax, fname_xyz_relax, 'Relaxed', lvs)
    
    print("\nDone!")

if __name__ == '__main__':
    main()
