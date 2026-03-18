#!/usr/bin/env python3
# coding: utf-8
import re
import sys
import os
import numpy as np
from scipy.spatial.transform import Rotation as R
import argparse
import matplotlib.pyplot as plt
import pubchempy as pcp
from ase.io import read, write
from ase import Atoms
from ase.build import add_adsorbate
from ase.optimize import BFGS
from ase.constraints import FixAtoms
import shutil
import torch


# --- PPAFM imports
from ppafm.ocl.AFMulator import AFMulator
from ppafm.ocl.field import HartreePotential
from ppafm.io import loadGeometryIN
from PIL import Image

def parse_from_CID(cid):
    """
    Helper function for accessing the PubChem API to obtain the 3D Conformer JSON
    and parse it into an ASE Atoms object.
    """
    conformer_json = pcp.get_json(cid, record_type='3d')
    elements = np.asarray(conformer_json['PC_Compounds'][0]['atoms']['element'])
    positions = np.zeros([len(elements), 3])
    positions[:, 0] = np.asarray(conformer_json['PC_Compounds'][0]['coords'][0]['conformers'][0]['x'])
    positions[:, 1] = np.asarray(conformer_json['PC_Compounds'][0]['coords'][0]['conformers'][0]['y'])
    positions[:, 2] = np.asarray(conformer_json['PC_Compounds'][0]['coords'][0]['conformers'][0]['z'])
    
    return Atoms(elements, positions=positions)

def load_calculator(use_mace, device):
    if use_mace:
        from mace.calculators import MACECalculator
        #model_path = 'models/MACE_lambda_10_scaleshift/MACE_model_stagetwo_compiled.model'
        model_path = "models/MACE_FM_finetuned/MACE_finetuned_compiled.model"
        try:
            return MACECalculator(model_paths=model_path, device=str(device), default_dtype="float32", cu_equivariance=True)
        except TypeError:
            print("cu_equivariance not supported. Loading MACE without it.")
            return MACECalculator(model_paths=model_path, device=str(device), default_dtype="float32")
    else:
        from nequip.ase import NequIPCalculator
        model_path = "models/nequip_full_dataset/deployed_model.nequip.pth"
        return NequIPCalculator.from_deployed_model(model_path, device=device)


def optimize_structure(atoms, calculator):
    """
    Runs structure optimization with BFGS until forces < 0.01 eV/Å.
    """
    atoms.set_calculator(calculator)
    dyn = BFGS(atoms)
    dyn.run(fmax=0.01)
    return atoms

def align_molecule_to_z(mol: Atoms):
    """
    Rotate the molecule so that the best-fit plane through its atoms is perpendicular to the z-axis.
    
    Parameters:
        mol (Atoms): ASE Atoms object (modified in place and returned)

    Returns:
        Atoms: Rotated ASE Atoms object
        np.ndarray: rotation axis
        float: rotation angle in degrees
    """
    positions = mol.get_positions()
    positions -= positions.mean(axis=0)  # Center molecule at origin

    # PCA: best-fit plane is defined by the first two principal components
    cov = np.cov(positions.T)
    eigvals, eigvecs = np.linalg.eigh(cov)
    
    # Normal to the best-fit plane is the eigenvector with the smallest eigenvalue
    normal_vec = eigvecs[:, np.argmin(eigvals)]
    
    # Calculate rotation from normal_vec to z-axis
    z_axis = np.array([0.0, 0.0, 1.0])
    axis = np.cross(normal_vec, z_axis)
    angle = np.arccos(np.clip(np.dot(normal_vec, z_axis), -1.0, 1.0))  # in radians

    if np.linalg.norm(axis) < 1e-6:
        return mol, axis, 0.0  # Already aligned

    axis /= np.linalg.norm(axis)
    rot = R.from_rotvec(axis * angle)
    rotated_positions = rot.apply(positions)

    mol.set_positions(rotated_positions + mol.get_positions().mean(axis=0))  # Re-center to original position
    return mol, axis, np.degrees(angle)


def get_folder_names():
    """
    Utility to return the grandparent and current folder name,
    for naming output files as in your original script.
    """
    current_directory = os.getcwd()
    current_folder = os.path.basename(current_directory)
    parent_directory = os.path.dirname(current_directory)
    grandparent_directory = os.path.dirname(parent_directory)
    grandparent_folder = os.path.basename(grandparent_directory)
    return grandparent_folder, current_folder

def run_ppm_simulation(geometry_path, mol_len):
    os.chdir(os.path.dirname(geometry_path))

    sample = read(geometry_path)
    
    # Get the scanning area of the molecule
    coords_min = sample[-mol_len:].positions.min(axis=0)
    coords_max = sample[-mol_len:].positions.max(axis=0)
    mean_molec_height = sample[-mol_len:].positions[:, 2].mean()

    
    grandparent_folder, current_folder = get_folder_names()
    filename = f"afm_{grandparent_folder}_{current_folder}.npz"

    print(f"Data will be saved to {filename}")


    # Initialize the AFMulator
    afmulator = AFMulator(
        scan_dim=(200, 200, 60),
        scan_window=((coords_min[0] - 3, coords_min[1] - 3, mean_molec_height + 6), 
                     (coords_max[0] + 3, coords_max[1] + 3, mean_molec_height + 9)),
        iZPP=8,
        # Amplitude is (scan_window[1][2] - scan_window[0][2]) / scan_dim * df_steps
        # Here (scan_window[1][2] - scan_window[0][2]) / scan_dim = 0.05 Angstrom
        df_steps=40, # 40 -> Amp = 2 A   
        tipStiffness=(0.37, 0.37, 0.0, 20.0),
        rho={"dz2": -0.10},
        tipR0=[0.0, 0.0, 4.0],
        npbc=(2, 1, 0),
    )

    print("Loading final geometry for AFM simulation...")
    xyzs, Zs, lvec = loadGeometryIN(geometry_path)

    # Charges set to zero
    qs = np.zeros(len(Zs), dtype=np.float64)

    print("Running AFM simulation...")
    X = afmulator(xyzs, Zs, qs, plot_to_dir="./simulated_images")

    # Create directory for images without axes
    output_directory = "simulated_without_axes"
    os.makedirs(output_directory, exist_ok=True)

    # Save grayscale images without axes
    for i in range(X.shape[-1]):
        pixels = X[:, :, i]
        pixels = ((pixels - pixels.min()) / np.ptp(pixels) * 255).astype(np.uint8)
        img = Image.fromarray(pixels[::-1].T, mode="L")
        img = img.rotate(180)
        img.save(os.path.join(output_directory, f"afm_{i}.png"))

    # Save the 3D AFM data volume
    np.savez(filename, afm=X)
    print("AFM simulation complete.")

def main():
    parser = argparse.ArgumentParser(description="Modular PPM MLIP pipeline")
    source_group = parser.add_mutually_exclusive_group(required=True)
    source_group.add_argument("--fromPubchem", type=int, help="PubChem CID of the molecule")
    source_group.add_argument("--fromGeom", type=str, help="Path to existing geometry file (e.g. .xyz, .in)")
    parser.add_argument("--surface", type=str, required=True, help="Surface folder under surfaces/")
    parser.add_argument("--no_align", action="store_true", help="Do not rotate molecule to XY plane")
    parser.add_argument("--adsorption", action="store_true", help="Run geometry optimization: gas + adsorbed")
    parser.add_argument("--afm", action="store_true", help="Run AFM simulations using relaxed geometries")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--nequip", action="store_true", help="Use NequIP model (default)")
    group.add_argument("--MACE", action="store_true", help="Use MACE model")



    args = parser.parse_args()

    if not args.adsorption and not args.afm:
        print("Nothing selected. Use --adsorption, --afm, or both.")
        sys.exit(0)
    
    if args.fromPubchem:
        CID = args.fromPubchem
    else:
        CID = re.sub(r'\W+', '_', os.path.splitext(os.path.basename(args.fromGeom))[0])

    surface_name = args.surface
    current_path = os.getcwd()
    out_folder = os.path.join(current_path, f"molecule_{CID}")
    pubchem_folder = os.path.join(out_folder, "from_pubchem")
    opt_gas_folder = os.path.join(out_folder, "optimized_gas_phase")
    opt_adsorbed = os.path.join(out_folder, "optimized_adsorbed")
    initial_mol_path = os.path.join(pubchem_folder, "molecule.in")
    final_mol_path = os.path.join(opt_gas_folder, "molecule.in")
    final_sys_path = os.path.join(opt_adsorbed, "molecule_on_surface.in")


    if args.adsorption:
        # 1) Load either MACE or Nequip calculator
        use_mace = args.MACE
        
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        print(f"Using device: {device}")
    
        print(f"Model selected: {'MACE' if use_mace else 'NequIP'}")
        calculator = load_calculator(use_mace, device)
    
        # 2) Read and prepare the surface
        surface_path = f"surfaces/{surface_name}/geometry.in"
        if not os.path.exists(surface_path):
            raise FileNotFoundError(f"Surface file not found: {surface_path}")
    
        surface = read(surface_path)
        # Translate so that the lowest z is at 0
        surface.translate([0, 0, -surface.positions[:, 2].min()])
    
        # Fix the bottom 2/3 of atoms in the surface
        n_atoms_fix = int(len(surface) / 3 * 2)
        constraint = FixAtoms(indices=[i for i in range(n_atoms_fix)])
    
        cell = surface.cell
    
        # 3) Fetch molecule from PubChem and set up folder structure
        os.makedirs(pubchem_folder, exist_ok=True)
        os.makedirs(opt_gas_folder, exist_ok=True)
        os.makedirs(opt_adsorbed, exist_ok=True)

        if args.fromPubchem:
            molecule = parse_from_CID(CID)
        else:
            molecule = read(args.fromGeom)
        
        molecule.cell = cell
        if args.no_align:
            print("--no_align flag was chosen. Molecules won't be fixed into the XY plane")
        else:
            molecule, axis, angle = align_molecule_to_z(molecule)
            print(f"Molecule was rotated by {angle:.2f} degrees around axis {axis} to align with the XY plane")

    
        # Center molecule in the surface cell
        cell_center = (
            (cell[0, 0] + cell[1, 0]) / 2.0,
            cell[1, 1] / 2.0
        )
        molecule.positions[:, 0] += -1 * molecule.positions[:, 0].mean() + cell_center[0]
        molecule.positions[:, 1] += -1 * molecule.positions[:, 1].mean() + cell_center[1]
    
        # Write initial molecule
        write(initial_mol_path, molecule)
    
        # 4) Relax the molecule alone
        print(f'Relaxing gas phase molecule with CID: {CID}')
        molecule_rlx = optimize_structure(molecule.copy(), calculator)
        write(final_mol_path, molecule_rlx)
    
        
        # 5) Adsorb the relaxed molecule onto the surface
        system = surface.copy()
        # Recompute cell center in case the surface cell changed
        cell_center = (
            (system.cell[0, 0] + system.cell[1, 0]) / 2.0,
            system.cell[1, 1] / 2.0
        )
        add_adsorbate(system, molecule_rlx, height=3.0, position=cell_center)
        system.set_constraint(constraint)
    
    
        # 6) Relax the combined system
        print('Relaxing adsorbed system')
        system_rlx = optimize_structure(system.copy(), calculator)
    
        # Write final optimized system
        write(final_sys_path, system_rlx)
    
    
        # 7) Print predicted adsorption height
        mean_molec_height = system_rlx[-len(molecule_rlx):].positions[:, 2].mean()
        max_surface_z = surface.positions[:, 2].max()
        adsorption_height = mean_molec_height - max_surface_z
        print(f"Predicted adsorption height: {adsorption_height:.3f} Å")
    
    # 8) AFM Simulation
    if args.afm:
        
        if not os.path.isfile(final_mol_path) or not os.path.isfile(final_sys_path):
            raise FileNotFoundError("AFM simulation requested, but optimized structures are missing.")

        molecule_rlx = read(final_mol_path)

        print(f'Computing PPM simulation of relaxed gas phase molecule in: {final_mol_path}')
        run_ppm_simulation(final_mol_path, len(molecule_rlx))
    
        print(f'Computing PPM simulation of relaxed adsorbed in: {final_sys_path}')
        run_ppm_simulation(final_sys_path, len(molecule_rlx))
       
    

if __name__ == "__main__":
    main()
