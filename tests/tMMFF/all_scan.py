#!/usr/bin/env python3

# Comprehensive scan module with all scan types
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF as mmff


def scanPlot(nscan=1000, span=(0.0,8.0), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="E", saveFig=None, saveData=None):
    """
    Standard scan using default MMFF settings
    """
    ts = np.linspace(span[0], span[1], nscan, endpoint=False)
  
    poss = np.zeros((nscan, 3))
    poss[:, 0] = p0[0] + ts * dir[0]
    poss[:, 1] = p0[1] + ts * dir[1]
    poss[:, 2] = p0[2] + ts * dir[2]

   
    Es, Fs, Ps = mmff.scan(poss, bF=True, bP=True)
    
    Es_end = Es[-1]
    # Es_end = 0

    if saveData is not None:
        np.savetxt(saveData, np.column_stack((ts, Es-Es_end)), header="ts\tEnergy", comments="# ")

    if saveFig is not None:
        plt.figure(figsize=(10, 6))
        plt.title(label)
        plt.plot(ts, Es-Es_end, '-', lw=1.5, label=label)
        plt.xlabel(f"Distance along ({dir[0]}_{dir[1]}_{dir[2]}) direction (\u00c5)")
        plt.ylabel(f"Energy (eV)")
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.legend()
        plt.tight_layout()
        plt.savefig(saveFig, dpi=300)
        plt.close()


def scanPlot_uff(nscan=1000, span=(0.0,8.0), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="E", saveFig=None, saveData=None):
    """
    Scan using UFF parameters (rigid scan)
    """
    ts = np.linspace(span[0], span[1], nscan, endpoint=False)
  
    poss = np.zeros((nscan, 3))
    poss[:, 0] = p0[0] + ts * dir[0]
    poss[:, 1] = p0[1] + ts * dir[1]
    poss[:, 2] = p0[2] + ts * dir[2]

   
    Es, Fs, Ps = mmff.scan_rigid_uff(poss, bF=True, bP=True)
    
    Es_end = Es[-1]
    # Es_end = 0

    if saveData is not None:
        np.savetxt(saveData, np.column_stack((ts, Es-Es_end)), header="ts\tEnergy", comments="# ")

    if saveFig is not None:
        plt.figure(figsize=(10, 6))
        plt.title(label)
        plt.plot(ts, Es-Es_end, '-', lw=1.5, label=label)
        plt.xlabel(f"Distance along ({dir[0]}_{dir[1]}_{dir[2]}) direction (\u00c5)")
        plt.ylabel(f"Energy (eV)")
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.legend()
        plt.tight_layout()
        plt.savefig(saveFig, dpi=300)
        plt.close()

def scanPlot2D_uff(nscan1=1000, nscan2=1000, span1=(0.0,4.0), span2=(0.0,4.0),p0=(0.0,0.0,0.0), dir1=(1.0,0.0,0.0), dir2=(0.0,1.0,0.0),label="E_2D_UFF", saveFig=None, saveData=None):
    """
    Perform a 2D scan using rigid UFF across two directions.
    
    Args:
        nscan1 (int): Number of scan points in first direction
        nscan2 (int): Number of scan points in second direction
        span1 (tuple): Distance range to scan in first direction (min, max)
        span2 (tuple): Distance range to scan in second direction (min, max)
        dir1 (tuple): First direction vector for the scan
        dir2 (tuple): Second direction vector for the scan
        p0 (tuple): Starting position
        label (str): Label for the plot
        saveFig (str): Path to save figure
        saveData (str): Path to save data
    
    Returns:
        tuple: (energies, x_grid, y_grid) - 2D energy array and coordinate grids
    """
    # Create coordinate arrays for the scan
    # Create linspace arrays for both scan directions
    t1 = np.linspace(span1[0], span1[1], nscan1, endpoint=False)
    t2 = np.linspace(span2[0], span2[1], nscan2, endpoint=False)
    
    # Generate a meshgrid for these parameters
    T1, T2 = np.meshgrid(t1, t2, indexing="ij")
    
    # Prepare positions array
    poss = np.zeros((nscan1*nscan2, 3))
    
    # Each scanned position is the starting point plus contributions along two directions
    poss[:, 0] = p0[0] + T1.ravel()*dir1[0] + T2.ravel()*dir2[0]
    poss[:, 1] = p0[1] + T1.ravel()*dir1[1] + T2.ravel()*dir2[1]
    poss[:, 2] = p0[2] + T1.ravel()*dir1[2] + T2.ravel()*dir2[2]
    
    # Call the scan function using the computed positions.
 
    Es, Fs, Ps = mmff.scan_rigid_uff(poss, bF=True, bP=True)
    
    # Reshape the energies into a 2D grid matching the T1, T2 shape
    Egrid = Es.reshape(nscan1, nscan2)
    
    # Create a contour plot for the 2D energy scan.
    plt.figure()
    cp = plt.contourf(T1, T2, Egrid, levels=20, cmap="viridis")
    plt.colorbar(cp)
    plt.title(label)
    plt.xlabel(f"Scan parameter along ({dir1[0]}_{dir1[1]}_{dir1[2]}) direction")
    plt.ylabel(f"Scan parameter along ({dir2[0]}_{dir2[1]}_{dir2[2]}) direction")

    if saveFig is not None:
        plt.savefig(saveFig)
    if saveData is not None:
        # Save in 3-column format: x, y, E
        data_out = np.column_stack((T1.ravel(), T2.ravel(), Es))
        header = f"# Scan coordinates along:\n# dir1: {dir1}\n# dir2: {dir2}\n# x\ty\tEnergy(eV)"
        np.savetxt(saveData, data_out, header=header, comments='')
    # plt.show()
    
    return Egrid, T1, T2



def relax_scanPlot1D(nscan=1000, span=(0.0,4.0),
               p0=(0.0,0.0,0.0), dir=(1.0,0.0,0.0),
               label="E_1D", saveFig=None, saveData=None,
               bRelax=False, niter_max=10000, dt=0.05, Fconv=1e-5):
    """
    Perform a 1D scan along a specified direction with optional relaxation.
    """
    # Create linspace array for scan direction
    t = np.linspace(span[0], span[1], nscan, endpoint=False)
    
    # Prepare positions array
    poss = np.zeros((nscan, 3))
    
    # Each scanned position is the starting point plus contribution along direction
    poss[:, 0] = p0[0] + t*dir[0]
    poss[:, 1] = p0[1] + t*dir[1]
    poss[:, 2] = p0[2] + t*dir[2]
    
    # Call the scan function using the computed positions.
    if bRelax:
        Es, Fs, Ps = mmff.scan(poss, bF=True, bP=True, bRelax=True, niter_max=niter_max, dt=dt, Fconv=Fconv)
    else:
        Es, Fs, Ps = mmff.scan(poss, bF=True, bP=True)
    
    print(f"Shape of returned arrays: Es={Es.shape}, Ps={Ps.shape}")
    
    # Create a line plot for the 1D energy scan
    plt.figure(figsize=(12, 8))
    plt.subplot(2, 1, 1)
    plt.plot(t, Es)
    plt.title(f"{label} - Energy Profile")
    plt.xlabel(f"Scan parameter along ({dir[0]:.3f}_{dir[1]:.3f}_{dir[2]:.3f}) direction")
    plt.ylabel("Energy (eV)")
    plt.grid(True)
    
    # Plot position information
    plt.subplot(2, 1, 2)
    
    # Handle multi-atom position data
    if len(Ps.shape) == 3:  # (n_scan, n_atoms, 3)
        n_atoms = Ps.shape[1]
        
        # Plot the center of mass movement for simplicity
        com = np.mean(Ps, axis=1)  # Average across atoms
        
        plt.plot(t, com[:, 0], 'r-', label="Center of mass - x")
        plt.plot(t, com[:, 1], 'g-', label="Center of mass - y")
        plt.plot(t, com[:, 2], 'b-', label="Center of mass - z")
        
        # Calculate and plot distance from initial position
        init_com = com[0]
        distances = np.sqrt(np.sum((com - init_com)**2, axis=1))
        plt.plot(t, distances, 'k--', label="Distance from start")
        
        title = "Center of mass trajectory"
    else:
        # Fallback for other shapes (unlikely with your data)
        for i in range(min(Ps.shape[1], 3)):
            plt.plot(t, Ps[:, i], label=f"Dimension {i}")
        title = "Position during scan"
    
    plt.title(title)
    plt.xlabel(f"Scan parameter")
    plt.ylabel("Position (Å)")
    plt.legend()
    plt.grid(True)
    
    plt.tight_layout()
    
    if saveFig is not None:
        plt.savefig(saveFig)
    if saveData is not None:
        # Remove file extension to use as base name
        base_name = saveData.rsplit('.', 1)[0]

        # Get atom types from MMFF
        atom_types = mmff.getAtomTypes()

        # Save energy data separately
        energy_file = f"{base_name}_energy.dat"
        energy_data = np.column_stack((t, Es))
        energy_header = f"# Scan along dir: {dir}\n# t\tEnergy(eV)"
        np.savetxt(energy_file, energy_data, header=energy_header, comments='')

        # Save combined trajectory (initial + relaxed structures)
        if len(Ps.shape) == 3:
            xyz_file = f"{base_name}_trajectory.xyz"
            n_atoms = Ps.shape[1]
            n_frames = len(t)

            with open(xyz_file, 'w') as f:
                # First frame: Write initial structure (before relaxation)
                f.write(f"{n_atoms}\n")
                f.write(f"Initial structure, E = {Es[0]:.6f} eV\n")
                # Write initial positions (first frame of Ps before relaxation)
                for j in range(n_atoms):
                    f.write(f"{atom_types[j]} {Ps[0,j,0]:12.6f} {Ps[0,j,1]:12.6f} {Ps[0,j,2]:12.6f}\n")

                # Then write all relaxed structures
                for i in range(n_frames):
                    f.write(f"{n_atoms}\n")
                    f.write(f"t = {t[i]:.3f}, E = {Es[i]:.6f} eV\n")
                    for j in range(n_atoms):
                        f.write(f"{atom_types[j]} {Ps[i,j,0]:12.6f} {Ps[i,j,1]:12.6f} {Ps[i,j,2]:12.6f}\n")
    
    plt.show()
    
    return t, Es, Ps


def scanPlot2D(nscan1=1000, nscan2=1000, span1=(0.0,4.0), span2=(0.0,4.0),
               p0=(0.0,0.0,0.0), dir1=(1.0,0.0,0.0), dir2=(0.0,1.0,0.0),
               label="E_2D", saveFig=None, saveData=None):
    """
    Perform a 2D scan along two specified directions.
    """
    # Create linspace arrays for each scan direction
    t1 = np.linspace(span1[0], span1[1], nscan1, endpoint=False)
    t2 = np.linspace(span2[0], span2[1], nscan2, endpoint=False)
    
    # Create a meshgrid for the two scan parameters
    T1, T2 = np.meshgrid(t1, t2, indexing='ij')
    
    # Prepare positions array
    poss = np.zeros((nscan1 * nscan2, 3))
    
    # Flatten the meshgrid for easier indexing
    T1_flat = T1.flatten()
    T2_flat = T2.flatten()
    
    # Each scanned position is the starting point plus contributions along both directions
    poss[:, 0] = p0[0] + T1_flat*dir1[0] + T2_flat*dir2[0]
    poss[:, 1] = p0[1] + T1_flat*dir1[1] + T2_flat*dir2[1]
    poss[:, 2] = p0[2] + T1_flat*dir1[2] + T2_flat*dir2[2]
    
    # Call the scan function using the computed positions.
    Es, Fs, Ps = mmff.scan(poss, bF=True, bP=True)
    
    # Reshape the energy array to match the 2D grid
    Es_2D = Es.reshape(nscan1, nscan2)
    
    # Create a 2D plot for the energy scan
    plt.figure(figsize=(10, 8))
    plt.contourf(T1, T2, Es_2D, 100, cmap='viridis')
    plt.colorbar(label='Energy (eV)')
    plt.title(f"{label} - Energy Surface")
    plt.xlabel(f"Scan parameter along ({dir1[0]:.3f}_{dir1[1]:.3f}_{dir1[2]:.3f}) direction")
    plt.ylabel(f"Scan parameter along ({dir2[0]:.3f}_{dir2[1]:.3f}_{dir2[2]:.3f}) direction")
    plt.grid(True)
    
    if saveFig is not None:
        plt.savefig(saveFig)
    if saveData is not None:
        # Save the 2D energy data
        np.savez(saveData, t1=t1, t2=t2, energy=Es_2D)
    
    plt.show()
    
    return t1, t2, Es_2D


def relax_scanPlot2D(nscan1=1000, nscan2=1000, span1=(0.0,4.0), span2=(0.0,4.0),
               p0=(0.0,0.0,0.0), dir1=(1.0,0.0,0.0), dir2=(0.0,1.0,0.0),
               label="E_2D", saveFig=None, saveData=None,
               bRelax=False, niter_max=10000, dt=0.05, Fconv=1e-5):
    """
    Perform a 2D scan along two specified directions with optional relaxation.
    """
    # Create linspace arrays for each scan direction
    t1 = np.linspace(span1[0], span1[1], nscan1, endpoint=False)
    t2 = np.linspace(span2[0], span2[1], nscan2, endpoint=False)
    
    # Create a meshgrid for the two scan parameters
    T1, T2 = np.meshgrid(t1, t2, indexing='ij')
    
    # Prepare positions array
    poss = np.zeros((nscan1 * nscan2, 3))
    
    # Flatten the meshgrid for easier indexing
    T1_flat = T1.flatten()
    T2_flat = T2.flatten()
    
    # Each scanned position is the starting point plus contributions along both directions
    poss[:, 0] = p0[0] + T1_flat*dir1[0] + T2_flat*dir2[0]
    poss[:, 1] = p0[1] + T1_flat*dir1[1] + T2_flat*dir2[1]
    poss[:, 2] = p0[2] + T1_flat*dir1[2] + T2_flat*dir2[2]
    
    # Call the scan function using the computed positions.
    if bRelax:
        Es, Fs, Ps = mmff.scan(poss, bF=True, bP=True, bRelax=True, niter_max=niter_max, dt=dt, Fconv=Fconv)
    else:
        Es, Fs, Ps = mmff.scan(poss, bF=True, bP=True)
    
    # Reshape the energy array to match the 2D grid
    Es_2D = Es.reshape(nscan1, nscan2)
    
    # Create a 2D plot for the energy scan
    plt.figure(figsize=(10, 8))
    plt.contourf(T1, T2, Es_2D, 100, cmap='viridis')
    plt.colorbar(label='Energy (eV)')
    plt.title(f"{label} - Energy Surface")
    plt.xlabel(f"Scan parameter along ({dir1[0]:.3f}_{dir1[1]:.3f}_{dir1[2]:.3f}) direction")
    plt.ylabel(f"Scan parameter along ({dir2[0]:.3f}_{dir2[1]:.3f}_{dir2[2]:.3f}) direction")
    plt.grid(True)
    
    if saveFig is not None:
        plt.savefig(saveFig)
    if saveData is not None:
        # Save the 2D energy data
        np.savez(saveData, t1=t1, t2=t2, energy=Es_2D)
    
    plt.show()
    
    return t1, t2, Es_2D


# Scan type functions - these set up the MMFF for different scan types
def setup_total_scan():
    """
    Set up MMFF for total potential scan (all components active)
    """
    mmff.setSwitches2(NonBonded=1, MMFF=1, SurfAtoms=1, GridFF=1)
    print(f"Using all force field components for total potential scan")
    print(f"PLQs shape: {mmff.PLQs.shape}, non-zero values: P={np.count_nonzero(mmff.PLQs[:,0])}, L={np.count_nonzero(mmff.PLQs[:,1])}, Q={np.count_nonzero(mmff.PLQs[:,2])}")
    return "total"


def setup_morse_scan():
    """
    Set up MMFF for London-Pauli only scan (no Coulomb)
    """
    mmff.setSwitches2(NonBonded=1, MMFF=1, SurfAtoms=1, GridFF=1)
    mmff.PLQs[:,2] = 0.0  # Zero out Coulomb component
    print(f"Disabling Coulomb interactions for London-Pauli only scan")
    print(f"Coulomb component zeroed: Q={np.all(mmff.PLQs[:,2] == 0.0)}")
    return "morse"


def setup_coulomb_scan():
    """
    Set up MMFF for Coulomb only scan (no London-Pauli)
    """
    mmff.setSwitches2(NonBonded=1, MMFF=1, SurfAtoms=1, GridFF=1)
    mmff.PLQs[:,0] = 0.0  # Zero out Pauli component
    mmff.PLQs[:,1] = 0.0  # Zero out London component
    print(f"Disabling London-Pauli components for Coulomb only scan")
    print(f"London-Pauli components zeroed: Pauli={np.all(mmff.PLQs[:,0] == 0.0)}, London={np.all(mmff.PLQs[:,1] == 0.0)}")
    return "coul"


def setup_pauli_scan():
    """
    Set up MMFF for Pauli-only scan (no London or Coulomb)
    """
    mmff.setSwitches2(NonBonded=1, MMFF=1, SurfAtoms=1, GridFF=1)
    mmff.PLQs[:,1] = 0.0  # Zero out London component
    mmff.PLQs[:,2] = 0.0  # Zero out Coulomb component
    print(f"Disabling London and Coulomb interactions for Pauli-only scan")
    print(f"London and Coulomb components zeroed: London={np.all(mmff.PLQs[:,1] == 0.0)}, Coulomb={np.all(mmff.PLQs[:,2] == 0.0)}")
    return "pauli"


def setup_london_scan():
    """
    Set up MMFF for London-only scan (no Pauli or Coulomb)
    """
    mmff.setSwitches2(NonBonded=1, MMFF=1, SurfAtoms=1, GridFF=1)
    mmff.PLQs[:,0] = 0.0  # Zero out Pauli component
    mmff.PLQs[:,2] = 0.0  # Zero out Coulomb component
    print(f"Disabling Pauli and Coulomb interactions for London-only scan")
    print(f"Pauli and Coulomb components zeroed: Pauli={np.all(mmff.PLQs[:,0] == 0.0)}, Coulomb={np.all(mmff.PLQs[:,2] == 0.0)}")
    return "london"


# Dictionary mapping scan type names to setup functions
# Note: We're keeping all scan types in the dictionary for reference,
# but we'll only allow total, morse, and coulomb to be used
SCAN_TYPES = {
    'total': setup_total_scan,
    'morse': setup_morse_scan,
    'coulomb': setup_coulomb_scan,
    'pauli': setup_pauli_scan,
    'london': setup_london_scan
}

# List of allowed scan types that can be run
ALLOWED_SCAN_TYPES = ['total', 'morse', 'coulomb']


def run_scan(molecule, substrate, output_dir, scan_type='total', scan_params=None, skip_init=False):
    """
    Run a scan with the specified scan type and parameters
    
    Args:
        molecule (str): Path to the molecule xyz file
        substrate (str): Path to the substrate xyz file
        output_dir (str): Directory to save output files
        scan_type (str): Type of scan to perform (only 'total', 'morse', 'coulomb' are supported)
        scan_params (dict): Parameters for the scan
        skip_init (bool): If True, skip MMFF initialization (assumes it's already initialized)
    
    Returns:
        bool: True if scan completed successfully
    
    Raises:
        ValueError: If scan_type is not a string or not one of the allowed types
        Exception: For other errors during scan execution
    """
    # Validate scan_type is a string
    if not isinstance(scan_type, str):
        raise ValueError(f"scan_type must be a string, got {type(scan_type).__name__}")
    
    # Validate scan_type is one of the allowed types
    if scan_type.lower() not in ALLOWED_SCAN_TYPES:
        raise ValueError(f"Invalid scan type: {scan_type}. Must be one of {ALLOWED_SCAN_TYPES}")
    
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Extract molecule name for file naming
    mol_name = os.path.basename(molecule)
    
    # Default scan parameters if not provided
    if scan_params is None:
        scan_params = {
            "nscan": 125,
            "span": (2.6, 15.1),
            "dir": (0.0, 0.0, 1.0),
            "p0": (0, 0, 0),
        }
    
    # Initialize MMFF with molecule and substrate only if not already initialized
    if not skip_init:
        print(f"Initializing MMFF for {scan_type.upper()} scan with {molecule} molecule and {substrate} substrate...")
        mmff.init(xyz_name=molecule, surf_name=substrate, bUFF=True, bSimple=True)
        mmff.getBuffs()

    # Set up scan based on scan_type
    scan_file_prefix = SCAN_TYPES[scan_type.lower()]()
    scan_label = f"{mol_name} {scan_type} scan"

    # Add label to scan parameters
    scan_params_with_label = scan_params.copy()
    scan_params_with_label["label"] = scan_label

    # Determine if this is a 1D or 2D scan based on parameters
    is_2d_scan = 'nscan1' in scan_params and 'nscan2' in scan_params

    # Generate data for the scan
    print(f"Generating {scan_type} potential {'2D' if is_2d_scan else '1D'} scan...")
    
    if is_2d_scan:
        # 2D scan
        scanPlot2D_uff(
            **scan_params_with_label,
            saveFig=f"{output_dir}/{mol_name}_{scan_file_prefix}_2d.png",
            saveData=f"{output_dir}/{mol_name}_{scan_file_prefix}_2d.dat"
        )
        print(f"{scan_type.capitalize()} 2D scan completed. Data saved in {output_dir}/{mol_name}_{scan_file_prefix}_2d.dat")
    else:
        # 1D scan
        scanPlot_uff(
            **scan_params_with_label,
            saveFig=f"{output_dir}/{mol_name}_{scan_file_prefix}.png",
            saveData=f"{output_dir}/{mol_name}_{scan_file_prefix}.dat"
        )
        print(f"{scan_type.capitalize()} 1D scan completed. Data saved in {output_dir}/{mol_name}_{scan_file_prefix}.dat")

    return True
