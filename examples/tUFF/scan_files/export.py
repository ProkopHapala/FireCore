
import os
import datetime
import numpy as np

def export_scan_data(
    output_dir,
    bUFF,
    scan_atom_info,
    z_scan_params,
    z_scan_data,
    xy_scan_params,
    xy_scan_data,
    surface_name
):
    """
    Exports Z-scan and XY-scan data to separate files.

    Args:
        output_dir (str): The directory to save the export files in.
        bUFF (bool): Flag to determine if UFF or MMFF is used.
        scan_atom_info (dict): Information about the scanned atom.
        z_scan_params (dict): Parameters used for the Z-scan.
        z_scan_data (tuple): A tuple containing Z coordinates and Fz forces.
        xy_scan_params (dict): Parameters used for the XY-scan.
        xy_scan_data (tuple): A tuple containing X, Y meshgrids and Fz grid.
        surface_name (str): The name of the surface file.
    """
    # Create the output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Generate timestamp and force field string
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    ff_string = "UFF" if bUFF else "MMFF"

    # --- Export Z-scan data ---
    z_filename = os.path.join(output_dir, f"{timestamp}_{ff_string}_z-scan.dat")
    Z, Fz = z_scan_data

    with open(z_filename, 'w') as f:
        f.write(f"# Z-Scan Data Export\n")
        f.write(f"# Timestamp: {timestamp}\n")
        f.write(f"# Force Field: {ff_string}\n")
        f.write(f"# Surface: {surface_name}\n")
        f.write(f"#\n")
        f.write(f"# Scanned Atom:\n")
        f.write(f"#   - Index: {scan_atom_info['index']}\n")
        f.write(f"#   - Type: {scan_atom_info['type']}\n")
        f.write(f"#   - Charge: {scan_atom_info['charge']:.4f} e\n")
        f.write(f"#   - Position: ({scan_atom_info['pos'][0]:.3f}, {scan_atom_info['pos'][1]:.3f}, {scan_atom_info['pos'][2]:.3f}) Å\n")
        f.write(f"#\n")
        f.write(f"# Z-Scan Parameters:\n")
        f.write(f"#   - Scan Range: {z_scan_params['start']:.3f} to {z_scan_params['end']:.3f} Å\n")
        f.write(f"#   - Step: {z_scan_params['step']:.4f} Å\n")
        f.write(f"#   - Fixed Position: X={z_scan_params['scan_x']:.3f}, Y={z_scan_params['scan_y']:.3f} Å\n")
        f.write(f"#\n")
        f.write(f"# Data Columns:\n")
        f.write(f"# 1. Z [Å]\n")
        f.write(f"# 2. Fz [eV/Å]\n")
        f.write(f"#-------------------------------------\n")
        
        if Z is not None and Fz is not None and len(Z) == len(Fz):
            np.savetxt(f, np.c_[Z, Fz], fmt='%.6f')

    print(f"Z-scan data exported to: {z_filename}")

    # --- Export XY-scan data ---
    xy_filename = os.path.join(output_dir, f"{timestamp}_{ff_string}_xy-scan.dat")
    X, Y, Fz_grid = xy_scan_data

    with open(xy_filename, 'w') as f:
        f.write(f"# XY-Scan Data Export\n")
        f.write(f"# Timestamp: {timestamp}\n")
        f.write(f"# Force Field: {ff_string}\n")
        f.write(f"# Surface: {surface_name}\n")
        f.write(f"#\n")
        f.write(f"# Scanned Atom:\n")
        f.write(f"#   - Index: {scan_atom_info['index']}\n")
        f.write(f"#   - Type: {scan_atom_info['type']}\n")
        f.write(f"#   - Charge: {scan_atom_info['charge']:.4f} e\n")
        f.write(f"#   - Position: ({scan_atom_info['pos'][0]:.3f}, {scan_atom_info['pos'][1]:.3f}, {scan_atom_info['pos'][2]:.3f}) Å\n")
        f.write(f"#\n")
        f.write(f"# XY-Scan Parameters:\n")
        f.write(f"#   - Scan Range: {xy_scan_params['start']:.3f} to {xy_scan_params['end']:.3f} Å\n")
        f.write(f"#   - Step: {xy_scan_params['step']:.4f} Å\n")
        f.write(f"#   - Fixed Position: Z={xy_scan_params['scan_z']:.3f} Å\n")
        f.write(f"#\n")
        f.write(f"# Data Columns:\n")
        f.write(f"# 1. X [Å]\n")
        f.write(f"# 2. Y [Å]\n")
        f.write(f"# 3. Fz [eV/Å]\n")
        f.write(f"#-------------------------------------\n")

        if X is not None and Y is not None and Fz_grid is not None:
            # Flatten the meshgrid data for column-wise saving
            x_flat = X.flatten()
            y_flat = Y.flatten()
            fz_flat = Fz_grid.flatten()
            if len(x_flat) == len(y_flat) == len(fz_flat):
                np.savetxt(f, np.c_[x_flat, y_flat, fz_flat], fmt='%.6f')

    print(f"XY-scan data exported to: {xy_filename}")
