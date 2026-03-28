#!/usr/bin/env python3

import os
import sys
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
if ROOT not in sys.path:
    sys.path.append(ROOT)

from pyBall.OCL.GridFFRelaxedScan import GridFFRelaxedScan, ensure_dir

def get_test_particle_plq(particle_type, charge_override):
    """Get PLQ coefficients for test particle based on atom type and charge"""
    # Actual REQ parameters from AtomTypes.dat
    atom_reqs = {
        'H': {'R': 1.4430, 'E': np.sqrt(0.00190802059), 'Q': 0.0},
        'O': {'R': 1.7500, 'E': np.sqrt(0.00260184625), 'Q': -0.4},
        'C': {'R': 1.9255, 'E': np.sqrt(0.00455323095), 'Q': 0.0}
    }
    
    req = atom_reqs.get(particle_type, atom_reqs['H'])
    if charge_override != 0.0:
        req['Q'] = charge_override
    
    # Proper Morse-derived PLQ coefficients
    # GridFF stores: V_P = Σ_i exp(-2α*(r-Ri)), V_L = Σ_i exp(-α*(r-Ri))  
    # Test particle provides: P = exp(2α*Rj), L = 2E * exp(α*Rj), Q = qj
    alpha = 1.5  # Default alphaMorse from GridFF generation
    P = np.exp(2 * alpha * req['R'])      # Pauli prefactor
    L = 2 * req['E'] * np.exp(alpha * req['R'])  # London prefactor  
    Q = req['Q']                           # Coulomb charge
    
    return (P, L, Q), req

def main():
    outdir = os.path.abspath('out_gridff_comparison')
    ensure_dir(outdir)
    
    # Use the same CaF2 GridFF cache
    gridff_path = os.path.join(ROOT, 'tests', 'tMMFF', 'data', 'CaF2_6L_Ni3_rect_nx2_nz1_L2_top', 'Bspline_PLQd_sigma_0p000.npy')
    
    if not os.path.exists(gridff_path):
        print(f'ERROR: GridFF cache not found at {gridff_path}')
        return 1
    
    # Initialize scanner (no need to prepare for plotting only)
    scanner = GridFFRelaxedScan(
        mol_path=os.path.join(ROOT, 'cpp', 'common_resources', 'xyz', 'PTCDA.xyz'),
        sub_xyz_path=os.path.join(ROOT, 'cpp', 'common_resources', 'Substrates', 'CaF2_6L_Ni3.xyz'),
        gridff_path=gridff_path,
        out_dir=outdir,
        mol_type_map={'C': 'C_R', 'O': 'O_2', 'H': 'H'},
        grid_alpha=1.5,
        grid_sigma=0.0,
        grid_step=[23.175/225.0, 20.070/200.0, 48.472/382.0],
    )
    
    scanner.prepare(anchor_z_above=4.0, lateral_shift=(2.0, 2.0), generate_grid=False)
    
    # Test particles and charges
    test_configs = [
        ('H', 0.0),    # H neutral
        ('H', +0.2),   # H positive
        ('O', 0.0),    # O neutral  
        ('O', -0.4),   # O negative (default)
    ]
    
    generated_files = []
    
    for particle_type, charge in test_configs:
        plq_plot, test_req = get_test_particle_plq(particle_type, charge)
        
        # Generate descriptive filename
        charge_str = f"{'p' if charge > 0 else 'm'}{abs(charge):.1f}" if charge != 0 else "0"
        basename = f"gridff_{particle_type}_{charge_str}"
        
        # XZ plot
        xz_png = os.path.join(outdir, f"{basename}_xz.png")
        scanner.plot_gridff_xz(xz_png, iy=-1, plq=plq_plot, vmax=0.1,
                              test_req=test_req, test_particle=particle_type)
        
        # XY plot
        xy_png = os.path.join(outdir, f"{basename}_xy.png")
        scanner.plot_gridff_xy(xy_png, z=7.0, plq=plq_plot, vmax=0.1,
                              test_req=test_req, test_particle=particle_type)
        
        generated_files.extend([xz_png, xy_png])
        
        # Print PLQ coefficients for reference
        print(f"{particle_type} (Q={charge:+.1f}): P={plq_plot[0]:.2e}, L={plq_plot[1]:.2e}, Q={plq_plot[2]:.2f}")
    
    print('\nGenerated plots:')
    for f in generated_files:
        print(f'  {f}')
    
    return 0

if __name__ == '__main__':
    sys.exit(main())
