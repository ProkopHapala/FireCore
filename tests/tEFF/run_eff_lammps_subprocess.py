import sys
import numpy as np
import os
import shutil
import tempfile
from pathlib import Path

sys.path.append("../../")
#import lammps_utils as lu
from pyBall import lammps_utils as lu

# EFF LAMMPS script template with placeholders for formatting
EFF_SCRIPT_TEMPLATE = """
units           electron
newton          on
boundary        f f f
atom_style      electron

read_data       {data_file}

pair_style      eff/cut 100.0 ecp 1 O
pair_coeff      * *

compute         energies all pair eff/cut
variable        eke equal c_energies[1]
variable        epauli equal c_energies[2]
variable        ecoul equal c_energies[3]
variable        erres equal c_energies[4]

thermo          1
thermo_style    custom step etotal v_eke v_epauli v_ecoul v_erres
thermo_modify   format float %23.15g

group           atoms type 1 2
fix             freeze atoms setforce 0.0 0.0 0.0

min_style       cg
compute         1 all property/atom spin eradius erforce
dump            2 all custom 1 mini.lammpstrj id type q c_1[1] c_1[2] x y z fx fy fz c_1[3]
minimize        0 1e-6 2000 4000

write_data      {output_file}
"""

# EFF-specific field mapping for energy extraction
EFF_FIELD_MAP = {
    'total': {'pattern': 'TotEng', 'index': 1},
    'eke': {'pattern': 'v_eke', 'index': 2},
    'epauli': {'pattern': 'v_epauli', 'index': 3},
    'ecoul': {'pattern': 'v_ecoul', 'index': 4},
    'erres': {'pattern': 'v_erres', 'index': 5}
}

# Callback function to create eFF input file content
def create_eff_input(elements, positions, charges, radii, **kwargs):
    """Callback function to create eFF input from XYZ data"""
    return lu.create_lammps_data_eff(elements, positions, charges, radii, **kwargs)

# Callback function to process eFF output
def process_eff_output(output_file, stdout):
    """Callback function to process eFF output"""
    # Read final geometry
    final_elements, final_positions = lu.read_lammps_data(output_file)
    
    # Return processed data
    return {
        'elements': final_elements,
        'positions': final_positions
    }

# Function to write XYZ file with energy in comment line
def write_xyz_with_energy(xyz_file, elements, positions, energy=None):
    """Write XYZ file with energy in comment line"""
    with open(xyz_file, 'w') as f:
        f.write(f"{len(elements)}\n")
        if energy is not None:
            f.write(f"Etot({energy['total']:.6f}) Eke({energy['eke']:.6f}) Epauli({energy['epauli']:.6f}) Ecoul({energy['ecoul']:.6f}) Erres({energy['erres']:.6f})\n")
        else:
            f.write("\n")
        
        for elem, pos in zip(elements, positions):
            f.write(f"{elem} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}\n")

# Example usage
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Process XYZ file with eFF LAMMPS')
    parser.add_argument('xyz_file', help='XYZ file to process')
    parser.add_argument('--output-dir', default='.', help='Directory for final output files')
    parser.add_argument('--scratch-dir', default=None, help='Temporary directory for intermediate files')
    parser.add_argument('--box-size', type=float, default=500, help='Box size for simulation')
    parser.add_argument('--max-frames', type=int, default=None, help='Maximum number of frames to process')
    parser.add_argument('--lammps-exec', default='lmp_serial', help='LAMMPS executable name/path')
    parser.add_argument('--electron-type', type=int, default=None, help='Atom type to use for electrons')
    parser.add_argument('--keep-scratch', action='store_true', help='Keep scratch directory after completion')
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Set up scratch directory
    if args.scratch_dir:
        scratch_dir = args.scratch_dir
        os.makedirs(scratch_dir, exist_ok=True)
    else:
        scratch_dir = tempfile.mkdtemp(prefix='lammps_scratch_')
    
    print(f"Using scratch directory: {scratch_dir}")
    
    try:
        # Process the frames with LAMMPS using the generic function directly
        results = lu.process_xyz_with_lammps(
            xyz_file=args.xyz_file,
            create_input_fn=create_eff_input,
            process_output_fn=process_eff_output,
            script_template=EFF_SCRIPT_TEMPLATE,
            field_map=EFF_FIELD_MAP,
            work_dir=scratch_dir,
            box_size=args.box_size,
            max_frames=args.max_frames,
            lammps_exec=args.lammps_exec,
            electron_type=args.electron_type
        )
        
        # Output a summary
        print(f"\nSuccessfully processed {len(results)} frames")
        
        # Write summary of energies to file
        summary_file = os.path.join(args.output_dir, "energy_summary.txt")
        with open(summary_file, 'w') as f:
            f.write("Frame Total_Energy E_ke E_pauli E_coul E_rres\n")
            for result in results:
                energies = result['energies']
                if energies:
                    f.write(f"{result['frame']} {energies.get('total', 0):.6f} {energies.get('eke', 0):.6f} {energies.get('epauli', 0):.6f} {energies.get('ecoul', 0):.6f} {energies.get('erres', 0):.6f}\n")
        
        print(f"Energy summary written to {summary_file}")
        
        # Save the geometries to XYZ with energy in comment line
        xyz_file = os.path.join(args.output_dir, "final_geometries.xyz")
        with open(xyz_file, 'w') as f:
            for result in results:
                processed = result.get('processed_data', {})
                elements = processed.get('elements', [])
                positions = processed.get('positions', [])
                energy = result['energies'] if 'energies' in result else None
                
                if elements and positions and len(elements) == len(positions):
                    f.write(f"{len(elements)}\n")
                    if energy is not None:
                        f.write(f"Etot({energy['total']:.6f}) Eke({energy['eke']:.6f}) Epauli({energy['epauli']:.6f}) Ecoul({energy['ecoul']:.6f}) Erres({energy['erres']:.6f})\n")
                    else:
                        f.write("\n")
                    
                    for elem, pos in zip(elements, positions):
                        f.write(f"{elem} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}\n")
        
        print(f"Final geometries written to {xyz_file}")
        
    finally:
        # Clean up scratch directory unless requested to keep it
        if not args.keep_scratch and (args.scratch_dir is None):
            shutil.rmtree(scratch_dir)
            print(f"Removed scratch directory: {scratch_dir}")
        elif args.keep_scratch:
            print(f"Kept scratch directory: {scratch_dir}")
