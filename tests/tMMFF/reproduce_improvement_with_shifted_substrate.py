#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import sys, os
import argparse

sys.path.append("../../")
from pyBall import MMFF as mmff
from pyBall import plotUtils
from pyBall import atomicUtils
from TipSplineOptimizer import _plot_single_improvement

def reproduce_improvement_with_shifted_substrate(improvement_dir, improvement_num=59, substrate_shift=[-10.0, -12.0, 2.5]):
    """
    Reproduce improvement plot using existing _plot_single_improvement function
    """
    
    # Load trajectory data
    aps_file = os.path.join(improvement_dir, f'best_{improvement_num:05d}_APs.npy')
    fcon_file = os.path.join(improvement_dir, f'best_{improvement_num:05d}_Fcon.npy')
    
    if not os.path.exists(aps_file):
        print(f"Error: {aps_file} not found")
        return
    
    print(f"Loading improvement {improvement_num} data...")
    APs = np.load(aps_file)  # (n_steps, n_atoms, 3)
    Fcon = np.load(fcon_file)  # (n_steps, 3)
    
    print(f"APs shape: {APs.shape}")
    print(f"Fcon shape: {Fcon.shape}")
    
    # Get optimization parameters
    meta_file = os.path.join(improvement_dir, 'meta.json')
    if os.path.exists(meta_file):
        import json
        with open(meta_file, 'r') as f:
            meta = json.load(f)
        ia_anchor = meta['ia_anchor']
        ia_opposite = meta['ia_opposite']
        target_pos = meta['target_pos']
        print(f"Anchor atom: {ia_anchor}, Opposite atom: {ia_opposite}")
        print(f"Target position: {target_pos}")
    else:
        ia_anchor = 28
        ia_opposite = 26
        target_pos = [5.0, 5.0, 10.0]
        print("Using default parameters")
    
    # Load spline control points for this improvement
    spline_file = os.path.join(improvement_dir, f'best_{improvement_num:05d}.dat')
    current_control_pts = None
    if os.path.exists(spline_file):
        print(f"Loading spline control points from {spline_file}")
        with open(spline_file, 'r') as f:
            lines = f.readlines()
        # Parse control points (skip first comment line)
        control_pts = []
        for line in lines[1:]:
            if line.strip() and not line.startswith('#'):
                parts = line.split()
                if len(parts) >= 3:
                    control_pts.append([float(parts[0]), float(parts[1]), float(parts[2])])
        current_control_pts = control_pts
        print(f"Loaded {len(current_control_pts)} control points")
    
    # Load molecule structure
    xyz_file = "../../cpp/common_resources/xyz/PTCDA.xyz"
    if os.path.exists(xyz_file):
        mol_data = atomicUtils.load_xyz(xyz_file)
        mol_pos = mol_data[0]  # positions (n_atoms, 3)
        mol_names = mol_data[2]  # element names
        print(f"Loaded molecule with {len(mol_pos)} atoms")
    else:
        print("Warning: Molecule file not found")
        print(f"Tried: {xyz_file}")
        mol_pos = None
        mol_names = None
    
    # Create dummy MMFF instance (needed for function but not used for plotting)
    mmff_instance = None
    
    # Create molecule bonds using findAllBonds with proper element types
    if mol_pos is not None and mol_names is not None:
        # Convert element symbols to atomic numbers
        element_types = []
        for element_symbol in mol_names:
            element_num = 6  # Default to carbon
            if element_symbol == 'H':
                element_num = 1
            elif element_symbol == 'C':
                element_num = 6
            elif element_symbol == 'N':
                element_num = 7
            elif element_symbol == 'O':
                element_num = 8
            elif element_symbol == 'S':
                element_num = 16
            element_types.append(element_num)
        
        # Convert positions to atom format expected by findAllBonds
        # findAllBonds expects atoms with format: [element_type, x, y, z]
        natoms = len(mol_pos)
        atoms_for_bonds = np.zeros((natoms, 4))
        atoms_for_bonds[:, 0] = element_types  # Use proper element types
        atoms_for_bonds[:, 1:4] = mol_pos  # Use positions for bond detection
        
        bonds, bonds_vecs = atomicUtils.findAllBonds(atoms_for_bonds, Rcut=3.0, RvdwCut=0.6)
        print(f"Created {len(bonds)} bonds using findAllBonds with proper element types")
    else:
        bonds = None
    
    # Create dummy energy data (function expects it but we can use simple values)
    Es = np.linspace(0, 1, len(APs))  # Simple linear energy profile
    
    # Use the existing plotting function
    print("Calling _plot_single_improvement with shifted substrate...")
    _plot_single_improvement(
        APs=APs,
        Fcon=Fcon,
        target_pos=target_pos,
        ia_anchor=ia_anchor,
        ia_opposite=ia_opposite,
        attempt=improvement_num,
        out_dir=".",  # Current directory
        mmff_instance=mmff_instance,
        substrate_shift=substrate_shift,
        surface_name="../../cpp/common_resources/Substrates/generated_rect/CaF2_6L_Ni3_rect_nx2_nz1_L2_top",
        plot_substrate=True,
        plot_molecule_atoms=False,  # As in original
        plot_molecule_bonds=True,   # As in original
        opt_target=target_pos,
        current_control_pts=current_control_pts,
        suffix="_shifted_substrate",
        title_extra=" (Shifted Substrate)",
        Es=Es
    )
    
    print(f"Saved shifted substrate plot: improvement_{improvement_num:05d}_shifted_substrate.png")
    print(f"Used substrate shift: {substrate_shift}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Reproduce improvement plot with shifted substrate')
    parser.add_argument('--improvement-dir', type=str, default='opt_3d_target',  help='Directory containing improvement data')
    parser.add_argument('--improvement-num', type=int, default=59,  help='Improvement number to reproduce')
    parser.add_argument('--substrate-shift', type=float, nargs=3,  default=[-11.0, -11.8, 2.5], help='Substrate shift (x y z) in Angstroms')
    
    args = parser.parse_args()
    
    reproduce_improvement_with_shifted_substrate(
        args.improvement_dir, 
        args.improvement_num, 
        args.substrate_shift
    )
