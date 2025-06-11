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
               bRelax=False, niter_max=10000, dt=0.05, Fconv=1e-5,cons_atom=26):
    """
    Perform a 1D scan along a specified direction with optional relaxation.
    """
    
    # Normalize the direction vector
    dir_array = np.array(dir, dtype=float)
    dir_norm = np.linalg.norm(dir_array)
    if dir_norm > 0:
        dir_normalized = dir_array / dir_norm
        print(f"Direction vector normalized: {dir} → {dir_normalized}")
    else:
        dir_normalized = np.array([1.0, 0.0, 0.0])  # Default to x-direction if zero vector
        print("Warning: Zero direction vector provided, defaulting to x-direction")
    
    # Create linspace array for scan direction
    t = np.linspace(span[0], span[1], nscan, endpoint=False)
   
    # Prepare positions array
    poss = np.zeros((nscan, 3))
    
    # Each scanned position is the starting point plus contribution along normalized direction
    poss[:, 0] = p0[0] + t*dir_normalized[0]
    poss[:, 1] = p0[1] + t*dir_normalized[1]
    poss[:, 2] = p0[2] + t*dir_normalized[2]
    
    # print(f"poss: {poss.shape}", poss )
    # Call the scan function using the computed positions
    # if bRelax:
    #     Es, Fs, Ps = mmff.scan(poss, bF=True, bP=True, bRelax=True, niter_max=niter_max, dt=dt, Fconv=Fconv)
    # else:
    #     Es, Fs, Ps = mmff.scan(poss, bF=True, bP=True)

    # scan_constr(nconf, ncontr, icontrs, contrs, Es=None, aforces=None, aposs=None, bHardConstr=False, omp=False, niter_max=10000, dt=0.05, Fconv=1e-5, Flim=100.0 ):
    iconstr = np.array( [cons_atom], np.int32)
    # iconstr = np.arange(1, 38, 1, dtype=np.int32)
    # iconstr = np.array( [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37], np.int32)
    # nconf = 10
    contrs = np.zeros( (nscan,len(iconstr),4), np.float64)
    contrs[:,0,0] = poss[:,0] # x
    contrs[:,0,1] = poss[:,1] # y
    contrs[:,0,2] = poss[:,2] # z
    contrs[:,:,3] = 1.0 # stiffness
    # print(f"contrs: {contrs.shape}", contrs )
    # print(f"iconstr: {iconstr.shape}", iconstr )
    with open( "scan_constr.xyz","w") as f: f.write("") # delete the file content to  old data

 
    Es, Fs, Ps = mmff.scan_constr( iconstr, contrs, bHardConstr=True, bF=True, bP=True, niter_max=niter_max, dt=dt, Fconv=Fconv)
    
    # print(f"Shape of returned arrays: Es={Es.shape}, Ps={Ps.shape}")
    
    # Deeply check Ps data
    # print(f"Ps data check:")
    # print(f"  Ps type: {type(Ps)}")
    # print(f"  Ps dtype: {Ps.dtype}")
    # print(f"  Ps min: {np.min(Ps)}, max: {np.max(Ps)}")
    # print(f"  Ps contains NaN: {np.isnan(Ps).any()}")
    # print(f"  Ps contains only zeros: {np.all(Ps == 0)}")
    # print(f"  Ps memory layout: {Ps.flags}")
    
    # Create the figure first
    fig = plt.figure(figsize=(18, 10))
    
    # Energy plot
    ax1 = fig.add_subplot(311)
    ax1.plot(t, Es)
    ax1.set_title(f"{label} - Energy Profile")
    ax1.set_xlabel(f"Scan parameter along ({dir[0]:.3f}_{dir[1]:.3f}_{dir[2]:.3f}) direction")
    ax1.set_ylabel("Energy (eV)")
    ax1.grid(True)
    
    # Position plot for both atoms
    ax2 = fig.add_subplot(312)
    
    # Plot the first atom coordinates
    if Ps.shape[1] >= 1:
        ax2.plot(t, Ps[:, 26, 0], 'r-', label="Atom 26 x")
        ax2.plot(t, Ps[:, 26, 1], 'g-', label="Atom 26 y")
        ax2.plot(t, Ps[:, 26, 2], 'b-', label="Atom 26 z")
    
    # Plot the second atom coordinates
    if Ps.shape[1] >= 2:
        ax2.plot(t, Ps[:, 29, 0], 'r--', label="Atom 29 x")
        ax2.plot(t, Ps[:, 29, 1], 'g--', label="Atom 29 y")
        ax2.plot(t, Ps[:, 29, 2], 'b--', label="Atom 29 z")
    
    ax2.set_title("Atom Positions During Scan")
    ax2.set_xlabel("Scan parameter")
    ax2.set_ylabel("Position (Å)")
    ax2.legend()
    ax2.grid(True)


    ax3 = fig.add_subplot(313)
    
    # Plot the first atom coordinates
    for i in range(Ps.shape[1]):
        ax3.plot(t, Ps[:, i, 0], 'r-', label=f"Atom {i} x")
        ax3.plot(t, Ps[:, i, 1], 'g-', label=f"Atom {i} y")
        ax3.plot(t, Ps[:, i, 2], 'b-', label=f"Atom {i} z")

      
    ax3.set_title("Atom Positions During Scan")
    ax3.set_xlabel("Scan parameter")
    ax3.set_ylabel("Position (Å)")
    ax3.legend()
    ax3.grid(True)

    plt.tight_layout()
    
    if saveFig is not None:
        plt.savefig(saveFig)
    
    if saveData is not None:
        # Remove file extension to use as base name
        base_name = saveData.rsplit('.', 1)[0]
        
        # Save energy data separately
        energy_file = f"{base_name}.dat"
        
        Es_end = Es[-1]
        # Es_end = 0
        energy_data = np.column_stack((t, Es-Es_end))
        # energy_data = np.column_stack((t, Es_end-Es))
        energy_header = f"# Scan along dir: {dir}\n# t\tEnergy(eV)"
        np.savetxt(energy_file, energy_data, header=energy_header, comments='')
        
         # Get atom type information for molecule
        atom_types = []

        for j in range(Ps.shape[1]):
            if j < 24:
                atom_types.append("C")  # Default to carbon if no specific type info
            elif j >= 24 and j < 30:
                atom_types.append("O")  # Default to hydrogen if no specific type info
            else:
                atom_types.append("H")  # Default to oxygen if no specific type info
                
            # atom_types.append(mmff.geteElements())

        # atom_types = mmff.getAtomTypes()
        # # Check what was returned and convert to element symbols
        # atom_symbols = []
        
        # if isinstance(atom_types, tuple) and len(atom_types) == 2:
        #     # If it returned a tuple (types, count), use the first element
        #     atom_types = atom_types[0]
        
        # for i in range(Ps.shape[1]):

        #     atom_symbols.append(atom_types[i])

        

        


        # Get substrate atom positions and types
        substrate_pos, _ = mmff.get_atom_positions()
        n_sub = substrate_pos.shape[0]
        
        # Create substrate atom types (default to Na and Cl alternating)
        substrate_types = []
        for j in range(n_sub):
            if j % 2 == 0:
                substrate_types.append("Na")
            else:
                substrate_types.append("Cl")
        
        # Save trajectory in XYZ format with both substrate and molecule
        if len(Ps.shape) == 3:
            xyz_file = f"{base_name}_trajectory.xyz"
            n_mol_atoms = Ps.shape[1]
            n_frames = len(t)
            total_atoms = n_mol_atoms + n_sub
            
            with open(xyz_file, 'w') as f:
                # Write all frames including substrate atoms
                for i in range(n_frames):
                    f.write(f"{total_atoms}\n")
                    f.write(f"t = {t[i]:.3f}, E = {Es[i]:.6f} eV\n")
                    
                    # First write molecule atoms (which move during relaxation)
                    for j in range(n_mol_atoms):
                        f.write(f"{atom_types[j]} {Ps[i,j,0]:12.6f} {Ps[i,j,1]:12.6f} {Ps[i,j,2]:12.6f}\n")

                    # Then write substrate atoms (fixed positions)
                    for j in range(n_sub):
                        f.write(f"{substrate_types[j]} {substrate_pos[j,0]:12.6f} {substrate_pos[j,1]:12.6f} {substrate_pos[j,2]:12.6f}\n")
                    
                    
    
   # plt.show()
    
    # Create a separate 3D trajectory plot
    plt.figure(figsize=(10, 8))
    ax3d = plt.axes(projection='3d')
    
    # Plot trajectory of atom 26
    if Ps.shape[1] >= 26:  # Check if atom 26 exists (0-based indexing)
        ax3d.plot3D(Ps[:, 26, 0], Ps[:, 26, 1], Ps[:, 26, 2], 'red', label="Atom 26", marker='o', linestyle='-', markersize=4)
    
    # Plot trajectory of atom 29
    if Ps.shape[1] >= 29:  # Check if atom 29 exists (0-based indexing)
        ax3d.plot3D(Ps[:, 29, 0], Ps[:, 29, 1], Ps[:, 29, 2], 'blue', label="Atom 29", marker='o', linestyle='-', markersize=4)
    
    ax3d.set_xlabel('X')
    ax3d.set_ylabel('Y')
    ax3d.set_zlabel('Z')
    ax3d.legend()
    plt.title("3D Atom Trajectories")
   # plt.show()
    
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


def run_scan(molecule, substrate, output_dir, scan_type, scan_params, skip_init=False):
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
    """
    # Initialize MMFF if not already initialized
    if not skip_init:
        print(f"Initializing MMFF for {scan_type.upper()} scan with {molecule} molecule and {substrate} substrate...")
        mmff.init(xyz_name=molecule, surf_name=substrate, bUFF=True, bSimple=True)
        mmff.getBuffs()
        
    # Get the molecule name for file naming
    mol_name = os.path.splitext(os.path.basename(molecule))[0]
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
        if scan_params.get('relaxed', False):
            # Run relaxed 1D scan
            relax_scanPlot1D(
                nscan=scan_params['nscan'],
                span=scan_params['span'],
                dir=scan_params['dir'],
                p0=scan_params['p0'],
                label=scan_label,
                saveFig=f"{output_dir}/{mol_name}_{scan_file_prefix}.png",
                saveData=f"{output_dir}/{mol_name}_{scan_file_prefix}.dat",
                bRelax=True,
                niter_max=scan_params['niter_max'],
                dt=scan_params['dt'],
                Fconv=scan_params['Fconv'],
                cons_atom=scan_params['cons_atom']
            )
        else:
            # Regular 1D scan
            scanPlot_uff(
                **scan_params_with_label,
            saveFig=f"{output_dir}/{mol_name}_{scan_file_prefix}.png",
            saveData=f"{output_dir}/{mol_name}_{scan_file_prefix}.dat"
            )
        print(f"{scan_type.capitalize()} 1D scan completed. Data saved in {output_dir}/{mol_name}_{scan_file_prefix}.dat")

    return True
