import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF as mmff


########## 1D Scan #####################################################
def scanPlot(nscan=1000, span=(0.0,8.0), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="E", saveFig=None, saveData=None):
    
    ts = np.linspace(span[0], span[1], nscan, endpoint=False)
  
    poss = np.zeros((nscan, 3))
    poss[:, 0] = p0[0] + ts * dir[0]
    poss[:, 1] = p0[1] + ts * dir[1]
    poss[:, 2] = p0[2] + ts * dir[2]

   
    Es, Fs, Ps = mmff.scan(poss, bF=True, bP=True)
    
    Es_end = Es[-1]

    if saveData is not None:
        np.savetxt(saveData, np.column_stack((ts, Es-Es_end)), header="ts\tEnergy", comments="# ")

    plt.title(label)
    plt.plot(ts, Es, '-', lw=0.5, label=label)
    plt.xlabel(f"Scaned along ({dir[0]}_{dir[1]}_{dir[2]}) direction ")
    plt.ylabel(f"Scaned Energy (eV)")
    
    # Optionally, save the figure to a file.
    if saveFig is not None:
        plt.savefig(saveFig)
    plt.show()

#/home/indranil/git/FireCore/tests/tMMFF

def relax_scanPlot1D_0(nscan=1000, span=(0.0,4.0),
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

def relax_scanPlot1D(nscan=1000, span=(0.0,4.0),
               p0=(0.0,0.0,0.0), dir=(1.0,0.0,0.0),
               label="E_1D", saveFig=None, saveData=None,
               bRelax=False, niter_max=10000, dt=0.05, Fconv=1e-5):
    """
    Perform a 1D scan along a specified direction with optional relaxation.
    """
    
    """
    Debug version with extensive logging
    """
    # print(f"===== DEBUGGING INFO =====")
    # print(f"Input parameters:")
    # print(f"  nscan: {nscan}")
    # print(f"  span: {span}")
    # print(f"  p0: {p0}")
    # print(f"  dir: {dir}")
    # print(f"  bRelax: {bRelax}")
    
    # Create linspace array for scan direction
    t = np.linspace(span[0], span[1], nscan, endpoint=False)

    #nscan = 10
    #p0=[0.0,0.0,8.6]
    #t = np.linspace( 0.0, 10.0, nscan, endpoint=False)

    
    # Prepare positions array
    poss = np.zeros((nscan, 3))
    
    # Each scanned position is the starting point plus contribution along direction
    poss[:, 0] = p0[0] + t*dir[0]
    poss[:, 1] = p0[1] + t*dir[1]
    poss[:, 2] = p0[2] + t*dir[2]
    
    # print(f"poss: {poss.shape}", poss )
    # Call the scan function using the computed positions
    # if bRelax:
    #     Es, Fs, Ps = mmff.scan(poss, bF=True, bP=True, bRelax=True, niter_max=niter_max, dt=dt, Fconv=Fconv)
    # else:
    #     Es, Fs, Ps = mmff.scan(poss, bF=True, bP=True)

    # scan_constr(nconf, ncontr, icontrs, contrs, Es=None, aforces=None, aposs=None, bHardConstr=False, omp=False, niter_max=10000, dt=0.05, Fconv=1e-5, Flim=100.0 ):
    iconstr = np.array( [29], np.int32)
    # nconf = 10
    contrs = np.zeros( (nscan,len(iconstr),4), np.float64)
    contrs[:,0,0] = poss[:,0] # x
    contrs[:,0,1] = poss[:,1] # y
    contrs[:,0,2] = poss[:,2] # z
    contrs[:,:,3] = 1.0 # stiffness
    print(f"contrs: {contrs.shape}", contrs )
    print(f"iconstr: {iconstr.shape}", iconstr )
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
        energy_file = f"{base_name}_energy.dat"
        energy_data = np.column_stack((t, Es))
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
                    
                    
    
    plt.show()
    
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
    plt.show()
    
    return t, Es, Ps

def relax_scanPlot1D_debug(nscan=1000, span=(0.0,4.0),
               p0=(0.0,0.0,0.0), dir=(1.0,0.0,0.0),
               label="E_1D", saveFig=None, saveData=None,
               bRelax=False, niter_max=10000, dt=0.05, Fconv=1e-5):
    """
    Debug version with extensive logging
    """
    # print(f"===== DEBUGGING INFO =====")
    # print(f"Input parameters:")
    # print(f"  nscan: {nscan}")
    # print(f"  span: {span}")
    # print(f"  p0: {p0}")
    # print(f"  dir: {dir}")
    # print(f"  bRelax: {bRelax}")
    
    # Create linspace array for scan direction
    t = np.linspace(span[0], span[1], nscan, endpoint=False)
    
    # Prepare positions array
    poss = np.zeros((nscan, 3))
    
    # Each scanned position is the starting point plus contribution along direction
    poss[:, 0] = p0[0] + t*dir[0]
    poss[:, 1] = p0[1] + t*dir[1]
    poss[:, 2] = p0[2] + t*dir[2]
    
    print(f"Input scan positions (first 5 points):")
    for i in range(min(5, nscan)):
        print(f"  Point {i}: ({poss[i, 0]:.3f}, {poss[i, 1]:.3f}, {poss[i, 2]:.3f})")
    
    # Call the scan function using the computed positions
    if bRelax:
        print(f"Calling mmff.scan with bRelax=True, niter_max={niter_max}")
        Es, Fs, Ps = mmff.scan(poss, bF=True, bP=True, bRelax=True, niter_max=niter_max, dt=dt, Fconv=Fconv)
    else:
        print("Calling mmff.scan with bRelax=False")
        Es, Fs, Ps = mmff.scan(poss, bF=True, bP=True)
    
    print(f"Shape of returned arrays: Es={Es.shape}, Ps={Ps.shape}")
    
    # Deeply check Ps data
    print(f"Ps data check:")
    print(f"  Ps type: {type(Ps)}")
    print(f"  Ps dtype: {Ps.dtype}")
    print(f"  Ps min: {np.min(Ps)}, max: {np.max(Ps)}")
    print(f"  Ps contains NaN: {np.isnan(Ps).any()}")
    print(f"  Ps contains only zeros: {np.all(Ps == 0)}")
    print(f"  Ps memory layout: {Ps.flags}")
    
    # Print some sample values from the array
    if len(Ps.shape) == 3:
        print(f"Sample atom positions (first 3 points):")
        for i in range(min(3, nscan)):
            for atom_idx in range(Ps.shape[1]):
                print(f"  Scan point {i}, Atom {atom_idx}: ({Ps[i, atom_idx, 0]:.6f}, {Ps[i, atom_idx, 1]:.6f}, {Ps[i, atom_idx, 2]:.6f})")
    
    # Create a simple 1D plot of a few key values to verify data
    plt.figure(figsize=(10, 6))
    plt.plot(t, Es, 'b-', label="Energy")
    
    # Try to extract and plot direct coordinate values to debug
    if len(Ps.shape) == 3 and Ps.shape[1] > 0:
        atom0_z = np.array([p[0][2] for p in Ps])  # Atom 0's z-coordinate
        plt.plot(t, atom0_z, 'r-', label="Atom 0 z-coord")
        
        if Ps.shape[1] > 1:  # If we have a second atom
            atom1_z = np.array([p[1][2] for p in Ps])  # Atom 1's z-coordinate
            plt.plot(t, atom1_z, 'g-', label="Atom 1 z-coord")
    
    plt.legend()
    plt.title("Raw data check")
    plt.grid(True)
    plt.show()
    
    # Proceed with original visualization if we have valid data
    if not np.all(Ps == 0) and not np.isnan(Ps).any():
        # Create the figure
        fig = plt.figure(figsize=(12, 10))
        
        # Energy plot
        ax1 = fig.add_subplot(211)
        ax1.plot(t, Es)
        ax1.set_title(f"{label} - Energy Profile")
        ax1.set_xlabel(f"Scan parameter")
        ax1.set_ylabel("Energy (eV)")
        ax1.grid(True)
        
        # Position plot
        ax2 = fig.add_subplot(212)
        
        if len(Ps.shape) == 3 and Ps.shape[1] >= 1:
            ax2.plot(t, Ps[:, 0, 0], 'r-', label="Atom 0 x")
            ax2.plot(t, Ps[:, 0, 1], 'g-', label="Atom 0 y")
            ax2.plot(t, Ps[:, 0, 2], 'b-', label="Atom 0 z")
            
            if Ps.shape[1] >= 2:
                ax2.plot(t, Ps[:, 1, 0], 'r--', label="Atom 1 x")
                ax2.plot(t, Ps[:, 1, 1], 'g--', label="Atom 1 y")
                ax2.plot(t, Ps[:, 1, 2], 'b--', label="Atom 1 z")
        
        ax2.set_title("Atom Positions During Scan")
        ax2.set_xlabel("Scan parameter")
        ax2.set_ylabel("Position (Å)")
        ax2.legend()
        ax2.grid(True)
        
        plt.tight_layout()
        plt.show()
    else:
        print("WARNING: Ps array contains invalid data - skipping main visualization")
    
    return t, Es, Ps

def visualize_molecular_trajectory(t, Ps, sample_indices=None):
    """
    Visualize the molecular trajectory during a scan.
    
    Parameters:
    - t: scan parameter values
    - Ps: position array with shape (n_scan_points, n_atoms, 3)
    - sample_indices: indices to plot (default: 5 evenly spaced points)
    """
    if sample_indices is None:
        # Choose 5 evenly spaced indices by default
        sample_indices = np.linspace(0, len(t)-1, 5, dtype=int)
    
    # Determine the shape
    n_scan_points, n_atoms, n_coords = Ps.shape
    
    # Create 3D plot
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot molecular structure at each sampled point
    cmap = plt.cm.jet
    colors = [cmap(i) for i in np.linspace(0, 1, len(sample_indices))]
    
    for idx, sample_idx in enumerate(sample_indices):
        # Plot atoms for this scan point with proper color
        c = colors[idx]
        ax.scatter(Ps[sample_idx, :, 0], 
                  Ps[sample_idx, :, 1], 
                  Ps[sample_idx, :, 2],
                  color=c, s=50, label=f"t={t[sample_idx]:.2f}")
        
        # Connect atoms with lines to show molecular structure
        for i in range(n_atoms):
            for j in range(i+1, n_atoms):
                ax.plot([Ps[sample_idx, i, 0], Ps[sample_idx, j, 0]],
                       [Ps[sample_idx, i, 1], Ps[sample_idx, j, 1]],
                       [Ps[sample_idx, i, 2], Ps[sample_idx, j, 2]],
                       color=c, alpha=0.5)
    
    # Also plot the trajectory of each atom across all scan points
    for atom_idx in range(n_atoms):
        ax.plot(Ps[:, atom_idx, 0], 
               Ps[:, atom_idx, 1], 
               Ps[:, atom_idx, 2], 
               'k-', alpha=0.2)
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.legend()
    plt.title("Molecular structure at different scan points")
    
    plt.tight_layout()
    plt.show()


def visualize_relaxed_structures(t, Ps, n_structures=None, bond_length_threshold=2.0):
    """
    Visualize the final relaxed molecular structures at different scan points.
    
    Parameters:
    - t: scan parameter values
    - Ps: position array with shape (n_scan_points, n_atoms, 3)
    - n_structures: number of structures to visualize (default: all structures)
    - bond_length_threshold: maximum distance to consider atoms as bonded (in Angstroms)
    """
    n_scan_points = Ps.shape[0]
    
    if n_structures is None:
        # If not specified, show all structures
        n_structures = n_scan_points
        indices = np.arange(n_scan_points)
    else:
        # Create evenly spaced indices based on desired number of structures
        indices = np.linspace(0, n_scan_points-1, n_structures, dtype=int)
    
    # Create 3D plot
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot molecular structure at each sampled point
    cmap = plt.cm.jet
    colors = [cmap(i) for i in np.linspace(0, 1, len(indices))]
    
    for idx, scan_idx in enumerate(indices):
        # Plot atoms for this scan point with proper color
        c = colors[idx]
        positions = Ps[scan_idx]
        
        # Plot atoms
        ax.scatter(positions[:, 0], 
                  positions[:, 1], 
                  positions[:, 2],
                  color=c, s=1, label=f"z = {t[scan_idx]:.2f} Å")
        
        # Connect atoms to all their first nearest neighbors
        for i in range(positions.shape[0]):
            # Calculate distances to all other atoms
            distances = np.sqrt(np.sum((positions - positions[i])**2, axis=1))
            distances[i] = np.inf  # Exclude self
            
            # Find all neighbors within the threshold distance
            neighbors = np.where(distances < bond_length_threshold)[0]
            
            # Draw bonds to all nearest neighbors
            for neighbor in neighbors:
                if i < neighbor:  # Avoid double drawing
                    ax.plot([positions[i, 0], positions[neighbor, 0]],
                           [positions[i, 1], positions[neighbor, 1]],
                           [positions[i, 2], positions[neighbor, 2]],
                           color=c, alpha=0.5)
    
    ax.set_xlabel('X (Å)')
    ax.set_ylabel('Y (Å)')
    ax.set_zlabel('Z (Å)')
    ax.legend()
    plt.title("Final relaxed structures at different scan points")
    
    # Set equal aspect ratio
    ax.set_box_aspect([1,1,1])
    
    plt.tight_layout()
    plt.show()

########## 2D Scan #####################################################
def relax_scanPlot2D(nscan1=1000, nscan2=1000, span1=(0.0,4.0), span2=(0.0,4.0),
               p0=(0.0,0.0,0.0), dir1=(1.0,0.0,0.0), dir2=(0.0,1.0,0.0),
               label="E_2D", saveFig=None, saveData=None,
               bRelax=False, niter_max=10000, dt=0.05, Fconv=1e-5):
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
    if bRelax:
        Es, Fs, Ps = mmff.scan(poss, bF=True, bP=True, bRelax=True, niter_max=niter_max, dt=dt, Fconv=Fconv)
    else:
        Es, Fs, Ps = mmff.scan(poss, bF=True, bP=True)
    
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
    plt.show()

def scanPlot2D(nscan1=1000, nscan2=1000, span1=(0.0,4.0), span2=(0.0,4.0),
               p0=(0.0,0.0,0.0), dir1=(1.0,0.0,0.0), dir2=(0.0,1.0,0.0),
               label="E_2D", saveFig=None, saveData=None):
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
 
    Es, Fs, Ps = mmff.scan(poss, bF=True, bP=True)
    
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
    plt.show()
#======== Body ###########

mmff.setVerbosity( verbosity=2, idebug=1 )


################################### Molecule and Substrate Selection##############################################

#mmff.init( xyz_name="data/xyz/pyridine", surf_name="data/NaCl_1x1_L2" )    
#mmff.init( xyz_name="data/xyz/nHexadecan_dicarboxylic", bMMFF=True  )     
# mmff.init( xyz_name="data/xyz/O", surf_name="data/xyz/NaCl_1x1_L3" )  
#mmff.init( xyz_name="data/xyz/H2O", surf_name="data/xyz/NaCl_1x1_L3" )    
# mmff.init( xyz_name="data/xyz/PTCDA", surf_name="data/xyz/NaCl_8x8_L3" )    
# mmff.init( xyz_name="data/xyz/PTCDA", surf_name="data/xyz/NaCl_8x8_L3_copy" )

# mmff.init( xyz_name="data/xyz/C.iz0", surf_name="data/xyz/C.iz0" )
# mmff.init( xyz_name="data/xyz/C.iz0", surf_name="data/xyz/Na.iz0" )
# mmff.init( xyz_name="data/xyz/C.iz0", surf_name="data/xyz/Cl.iz0" )
# mmff.init( xyz_name="data/xyz/C.iz0", surf_name="data/xyz/NaCl.iz0" )
# mmff.init( xyz_name="data/xyz/C_Cl_iz0", surf_name="data/xyz/NaCl.iz0" )

# mmff.init( xyz_name="data/xyz/Na+Cl-", surf_name="data/xyz/NaCl_8x8_L3_Coulumb" )

# mmff.init( xyz_name="data/xyz/PTCDA", surf_name="data/xyz/NaCl_coulomb.iz0" )

# mmff.init( xyz_name="data/xyz/PTCDA_paolo", surf_name="data/xyz/NaCl.ptcda" )

# mmff.init( xyz_name="data/xyz/PTCDA_charge_paolo", surf_name="data/xyz/NaCl_paolo" )
# mmff.init( xyz_name="data/xyz/ptcda_charge", surf_name="data/xyz/NaCl_paolo" )


# mmff.init( xyz_name="data/xyz/molNaCl_Na.iz0", surf_name="data/xyz/NaCl_coulomb.iz0" )
# mmff.init( xyz_name="data/xyz/molNaCl_Cl.iz0", surf_name="data/xyz/NaCl_coulomb.iz0" )
# mmff.init( xyz_name="data/xyz/molNaCl.ix0.iy0", surf_name="data/xyz/NaCl_coulomb.iz0" )

# mmff.init( xyz_name="data/xyz/2_atom_NaCl_mol", surf_name="data/xyz/2_atom_NaCl" )


###***** NaCl Molecule #################
# mmff.init( xyz_name="data/xyz/molNa_0.9_Cl_-0.9_Na", surf_name="data/xyz/Na_0.9_Cl_-0.9" )
# mmff.init( xyz_name="data/xyz/molNa_0.9_Cl_-0.9_Cl", surf_name="data/xyz/Na_0.9_Cl_-0.9" )
# mmff.init( xyz_name="data/xyz/molNa_0.9_Cl_-0.9", surf_name="data/xyz/Na_0.9_Cl_-0.9" )
# mmff.init( xyz_name="data/xyz/molNa_0.9_Cl_-0.9_xy", surf_name="data/xyz/Na_0.9_Cl_-0.9" )
# mmff.init( xyz_name="data/xyz/H2O", surf_name="data/xyz/Na_0.9_Cl_-0.9" )


############### PTCDA On NaCl +0.9 -0.9 ##############
# mmff.init( xyz_name="data/xyz/new_PTCDA_charge_on_Na", surf_name="data/xyz/Na_0.9_Cl_-0.9" , bUFF=True, bSimple=True )   ### For uff relaxed scan the bUFF has to be true


#mmff.init( xyz_name="data/xyz/PTCDA_charge_on_Na", surf_name="data/xyz/Na_0.9_Cl_-0.9", bUFF=True, bSimple=True )   
mmff.init( xyz_name="data/xyz/PTCDA_charge_on_Na", surf_name="data/xyz/Na_0.9_Cl_-0.9" )   
# mmff.init( xyz_name="data/xyz/PTCDA_charge_on_Cl", surf_name="data/xyz/Na_0.9_Cl_-0.9" )
# mmff.init( xyz_name="data/xyz/PTCDA_charge_on_hollow", surf_name="data/xyz/Na_0.9_Cl_-0.9" )
# mmff.init( xyz_name="data/xyz/PTCDA_charge_xy", surf_name="data/xyz/Na_0.9_Cl_-0.9" )
# print("After init: ", mmff.ndims if hasattr(mmff, 'ndims') else "ndims not set yet")


####################################################################################################################################

mmff.getBuffs()
# print("After getBuffs: ndims=", mmff.ndims)
# print("natoms=", mmff.natoms, "nnode=", mmff.nnode, "ncap=", mmff.ncap)
# mmff.ipicked = 30

#print( "ffflags ", mmff.ffflags )

# mmff.setSwitches( NonBonded=-1, MMFF=1, SurfAtoms=0, GridFF=1,PBC=-1 )   ### For Relaxed Scan MMFF has to be 1 
# mmff.setSwitches( NonBonded=-1, MMFF=1, SurfAtoms=0, GridFF=1 )   ### For Relaxed Scan MMFF has to be 1 
mmff.setSwitches( NonBonded=-1, MMFF=1, SurfAtoms=1, GridFF=1 )   #### For Rigid Scan to make ay of the flag noneffective eed to set -1 0 will not work 



################# Mode Decision Morse Coulomb #######################################################################################################################
# mmff.PLQs[:,0] = 0.0  # delete Pauli
# mmff.PLQs[:,1] = 0.0  # delete London

# mmff.PLQs[:,2 ] = 0.0 # delete Coulomb (charges)

###################################################################################################################################################

######################################################################################## 1D Scan Plotting #########################################
# scanPlot( nscan=10, span=(0.0,4.0), dir=(1.0,0.0,0.0), p0=(0.0,0.0,0.0),  label="E_x", saveFig="E_x_scan.png", saveData="E_x_scan.dat")
# scanPlot( nscan=1000, span=(0.0,4.0), dir=(0.0,1.0,0.0), p0=(0.0,0.0,0.0),  label="E_y", saveFig="E_y_scan.png", saveData="E_y_scan.dat" )

# scanPlot( nscan=100, span=(3.3,8.3), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="Atop_C", saveFig="E_z_scan_Atop_C.png", saveData="E_z_scan_Atop_C.dat" )
# scanPlot( nscan=100, span=(2.9,8.3), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="Atop_C on Na", saveFig="E_z_scan_Atop_Na.png", saveData="E_z_scan_Atop_Na.dat" )
# scanPlot( nscan=100, span=(3.3,8.3), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="Atop_C on Cl ", saveFig="E_z_scan_Atop_Cl.png", saveData="E_z_scan_Atop_Cl.dat" )
# scanPlot( nscan=100, span=(0.8,6.3), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="Atop_Na", saveFig="E_z_scan_Atop_Na_SYM.png", saveData="E_z_scan_Atop_Na_SYM.dat" )
# scanPlot( nscan=100, span=(2.3,8.3), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="Atop_Na", saveFig="E_z_scan_Atop_Na_Ewald.png", saveData="E_z_scan_Atop_Na_Ewald.dat" )
# scanPlot( nscan=100, span=(7.9,18.3), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="Atop_Na", saveFig="E_z_scan_Atop_Na_Ewald.png", saveData="E_z_scan_Atop_Na_Ewald.dat" )
# scanPlot( nscan=100, span=(3.1,8), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="Atop_Cl", saveFig="E_z_scan_Atop_Cl_No_SYM.png", saveData="E_z_scan_Atop_Cl_No_SYM.dat" )

# scanPlot( nscan=100, span=(3,10), dir=(0.0,0.0,1.0), p0=(0.0,0.0,5.656854), label="PTCDA z0=0", saveFig="E_z_scan_PTCDA_Ewald.png", saveData="E_z_scan_PTCDA_Ewald.dat" )

# scanPlot( nscan=100, span=(3,10), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0), label="PTCDA", saveFig="E_z_scan_PTCDA_Ewald_new_Direct.png", saveData="E_z_scan_PTCDA_Ewald_new_Direct.dat" )
# scanPlot( nscan=100, span=(3,10), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0), label="PTCDA", saveFig="E_z_scan_PTCDA_Ewald_new_trial.png", saveData="E_z_scan_PTCDA_Ewald_new_trial.dat" )
# scanPlot( nscan=100, span=(3,10), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0), label="PTCDA", saveFig="E_z_scan_PTCDA_Morse_new.png", saveData="E_z_scan_PTCDA_Morse_new.dat" )

# scanPlot( nscan=101, span=(3,8.3), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="Atop_Na", saveFig="E_z_scan_Atop_Na_Ewald.png", saveData="E_z_scan_Atop_Na_Ewald.dat" )
# scanPlot( nscan=101, span=(3,8.3), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="Atop_Cl", saveFig="E_z_scan_Atop_Cl_Ewald.png", saveData="E_z_scan_Atop_Cl_Ewald.dat" )

# scanPlot( nscan=100, span=(3.1,8.1), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="Atop_Na", saveFig="E_z_scan_Atop_Na_Morse.png", saveData="E_z_scan_Atop_Na_Morse.dat" )
# scanPlot( nscan=100, span=(2.5,8.5), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="Atop_Cl", saveFig="E_z_scan_Atop_Cl_Morse.png", saveData="E_z_scan_Atop_Cl_Morse.dat" )

# scanPlot( nscan=100, span=(0,20), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="Atop_Na", saveFig="E_z_scan_2atom_NaCl_Ewald.png", saveData="E_z_scan_2atom_NaCl_Ewald.dat" )

# scanPlot( nscan=100, span=(2.6,10), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0), label="PTCDA on Na", saveFig="E_z_scan_on_Na_PTCDA_NaCl.png", saveData="E_z_scan_on_Na_PTCDA_NaCl.dat" )
# scanPlot( nscan=100, span=(2.6,10), dir=(0.0,0.0,1.0), p0=(2.0,2.0,0), label="PTCDA on Cl", saveFig="E_z_scan_on_Cl_PTCDA_NaCl.png", saveData="E_z_scan_on_Cl_PTCDA_NaCl.dat" )

################******************* NaCl
# scanPlot( nscan=150, span=(1.5,16.5), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0), label="NaCl_mol on Na", saveFig=None, saveData="trial_E_z_scan_on_Na_mol_NaCl_Morse.dat" )
# scanPlot( nscan=150, span=(1.5,16.5), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0), label="NaCl_mol on Na", saveFig="E_z_scan_on_Na_mol_NaCl_Morse.png", saveData="E_z_scan_on_Na_mol_NaCl_Morse.dat" )
# scanPlot( nscan=150, span=(1.5,16.5), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0), label="NaCl_mol on Na", saveFig="E_z_scan_on_Na_mol_NaCl_Coulomb.png", saveData="E_z_scan_on_Na_mol_NaCl_Coulomb.dat" )
#scanPlot( nscan=150, span=(1.5,16.5), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0), label="NaCl_mol on Na", saveFig="E_z_scan_on_Na_mol_NaCl_Morse_Coulomb.png", saveData="E_z_scan_on_Na_mol_NaCl_Morse_Coulomb.dat" )

# scanPlot( nscan=150, span=(1.5,16.5), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0), label="NaCl_mol on Cl", saveFig="E_z_scan_on_Cl_mol_NaCl_Morse.png", saveData="E_z_scan_on_Cl_mol_NaCl_Morse.dat" )
# scanPlot( nscan=150, span=(1.5,16.5), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0), label="NaCl_mol on Cl", saveFig="E_z_scan_on_Cl_mol_NaCl_Coulomb.png", saveData="E_z_scan_on_Cl_mol_NaCl_Coulomb.dat" )
# scanPlot( nscan=150, span=(1.5,16.5), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0), label="NaCl_mol on Cl", saveFig="E_z_scan_on_Cl_mol_NaCl_Morse_Coulomb.png", saveData="E_z_scan_on_Cl_mol_NaCl_Morse_Coulomb.dat" )

# scanPlot( nscan=150, span=(1.5,16.5), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0), label="NaCl_mol on hollow", saveFig="E_z_scan_on_hollow_mol_NaCl_Morse.png", saveData="E_z_scan_on_hollow_mol_NaCl_Morse.dat" )
# scanPlot( nscan=150, span=(1.5,16.5), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0), label="NaCl_mol on hollow", saveFig="E_z_scan_on_hollow_mol_NaCl_Coulomb.png", saveData="E_z_scan_on_hollow_mol_NaCl_Coulomb.dat" )
# scanPlot( nscan=150, span=(1.5,16.5), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0), label="NaCl_mol on hollow", saveFig="E_z_scan_on_hollow_mol_NaCl_Morse_Coulomb.png", saveData="E_z_scan_on_hollow_mol_NaCl_Morse_Coulomb.dat" )


################*******************PTCDA
scanPlot( nscan=125, span=(2.6,15.1), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0), label="PTCDA on Na", saveFig="E_z_scan_on_Na_PTCDA_Morse.png", saveData="E_z_scan_on_Na_PTCDA_Morse.dat" )
# scanPlot( nscan=125, span=(2.6,15.1), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0), label="PTCDA on Na", saveFig="E_z_scan_on_Na_PTCDA_Coulomb.png", saveData="E_z_scan_on_Na_PTCDA_Coulomb.dat" )
# scanPlot( nscan=125, span=(2.6,15.1), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0), label="PTCDA on Na", saveFig="E_z_scan_on_Na_PTCDA_Morse_Coulomb.png", saveData="E_z_scan_on_Na_PTCDA_Morse_Coulomb.dat" )

# scanPlot( nscan=125, span=(2.6,15.1), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0), label="PTCDA on Cl", saveFig="E_z_scan_on_Cl_PTCDA_Morse.png", saveData="E_z_scan_on_Cl_PTCDA_Morse.dat" )
# scanPlot( nscan=125, span=(2.6,15.1), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0), label="PTCDA on Cl", saveFig="E_z_scan_on_Cl_PTCDA_Coulomb.png", saveData="E_z_scan_on_Cl_PTCDA_Coulomb.dat" )
# scanPlot( nscan=125, span=(2.6,15.1), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0), label="PTCDA on Cl", saveFig="E_z_scan_on_Cl_PTCDA_Morse_Coulomb.png", saveData="E_z_scan_on_Cl_PTCDA_Morse_Coulomb.dat" )

# scanPlot( nscan=125, span=(2.6,15.1), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0), label="PTCDA on hollow", saveFig="E_z_scan_on_hollow_PTCDA_Morse.png", saveData="E_z_scan_on_hollow_PTCDA_Morse.dat" )
# scanPlot( nscan=125, span=(2.6,15.1), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0), label="PTCDA on hollow", saveFig="E_z_scan_on_hollow_PTCDA_Coulomb.png", saveData="E_z_scan_on_hollow_PTCDA_Coulomb.dat" )
# scanPlot( nscan=125, span=(2.6,15.1), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0), label="PTCDA on hollow", saveFig="E_z_scan_on_hollow_PTCDA_Morse_Coulomb.png", saveData="E_z_scan_on_hollow_PTCDA_Morse_Coulomb.dat" )



import gc
gc.disable()

###########********************* Relax 1D
# t,Es,Ps=relax_scanPlot1D(bRelax=True, nscan=125, span=(2.6,15.1), dir=(0.0,0.0,1.0), p0=(0.0,0.0,(0+2.6)), label="PTCDA on Na", saveFig=None, saveData=None,niter_max=100 )   
# t,Es,Ps=relax_scanPlot1D(bRelax=True, nscan=125, span=(2.6,15.1),  dir=(0.0,0.0,1.0), p0=(0.0,0.0,(0+0)), label="PTCDA on Na", saveFig=None, saveData="new_trial_relax_scan_ptcda_test",
#                         niter_max=100000,Fconv=1e-3,dt=0.02 )  ### z scan dt 0.05 is giving energy in the order of less than 100 but for more smaller step it is giving absolute energy in the order of 1e7 and greater value of like 0.1 is giving random values 0.1 is to match with LAMMPs 0.001femto

# t,Es,Ps=relax_scanPlot1D(bRelax=True, nscan=120, span=(0,12), dir=(0.866,0.5,0.0), p0=(0.0,0.0,(0+3.1)), label="PTCDA on Na", saveFig=None, saveData="trial_relax_scan_ptcda_line_test",
#                         niter_max=50000,Fconv=1e-6,dt=0.1 )  # x y scan  and diagonal #dir=(0.866,0.5,0.0) for 30 degree  nscan=351, span=(0,35.1)

# # Visualize trajectory
# visualize_molecular_trajectory(t, Ps)
# visualize_relaxed_structures(t, Ps,n_structures=25)


######Add these lines at the end of your script
plt.close('all')  # Close all matplotlib figures
gc.enable()
gc.collect()      # Force garbage collection

"""
The time i FireCore is            1.0180506e-14 s   

In lammps using 0.001 pico second == 1e-15 second 
In lammps they are using fire algorithm which can gradually decrease the time step if it detects instability in the system dynamics.


sed -i 's/\xEF\xBB\xBF//g' /home/indranil/git/FireCore/cpp/common/molecular/MolWorld_sp3.h

real	77m33.584s
user	76m30.038s
sys	0m5.760s

9.467951 -5.65

"""



# LLAMPS TIME 1D relaxed scan o top of Na aotm 
# real	183m28.087s
# user	182m25.173s
# sys	0m6.427s

########################################### 2D Scan Plotting #########################################
# scanPlot2D(nscan1=41, nscan2=41,
#            span1=(0.0, 4.1), span2=(0.0, 4.1),
#            dir1=(1.0,0.0,0.0), dir2=(0.0,1.0,0.0),
#            p0=(0.0,0.0,2.0), label="E_xy for z=2.0", #) #,
#            saveFig="E_xy_scan_mol_NaCl_Morse.png", saveData="E_xy_scan_mol_NaCl_Morse.dat")

# scanPlot2D(nscan1=41, nscan2=41,
#            span1=(0.0, 4.1), span2=(0.0, 4.1),
#            dir1=(1.0,0.0,0.0), dir2=(0.0,1.0,0.0),
#            p0=(0.0,0.0,3.3), label="E_xy for z=3.3", #) #,
#            saveFig="E_xy_scan_PTCDA_Morse_Coulomb.png", saveData="E_xy_scan_PTCDA_Morse_Coulomb.dat")

# scanPlot2D(nscan1=41, nscan2=61,
#            span1=(0.0, 4.1), span2=(2.3, 8.3),
#            dir1=(1.0,0.0,0.0), dir2=(0.0,0.0,1.0),
#            p0=(0.0,0.0,0), label="E_xz for Y=3.4",
#            saveFig="E_xz_scan_C.png", saveData="E_xz_scan_C.dat")

# scanPlot2D(nscan1=41, nscan2=61,
#            span1=(0.0, 4.1), span2=(2.3, 8.3),
#            dir1=(0.0,1.0,0.0), dir2=(0.0,0.0,1.0),
#            p0=(0.0,0.0,0), label="E_xz for Y=3.4",
#            saveFig="E_yz_scan_C.png", saveData="E_yz_scan_C.dat") 




############################ Relaxed 2D Scan ***********************************************************
# relax_scanPlot2D(nscan1=64, nscan2=64,
#            span1=(0.0, 32), span2=(0.0, 32),
#            dir1=(1.0,0.0,0.0), dir2=(0.0,1.0,0.0),
#            p0=(0.0,0.0,3.1), label="E_xy for z=3.1", #) #,
#            bRelax=True,
#            niter_max=1000,
#            dt=0.05,
#            Fconv=1e-5, #)
#            saveFig="E_xy_Relax_scan_PTCDA_NaCl.png", saveData="E_xy_Relax_scan_PTCDA_NaCl.dat")

        

          

# plt.legend()
# plt.grid()
# plt.show()
# print(mmff.PLQs.shape)



# Initialize the system. Make sure you call mmff.init() and then mmff.getBuffs()
# mmff.setVerbosity(verbosity=1, idebug=1)
# mmff.init(xyz_name="data/xyz/PTCDA", surf_name="data/xyz/NaCl_8x8_L3_copy")
# mmff.getBuffs()

# (Optionally) set up any switches, etc.
# mmff.setSwitches(NonBonded=-1, MMFF=-1, SurfAtoms=0, GridFF=1)

# Now, call our new saveXSF_geometry function.
# We assume that: 
#   - The grid data is stored in an array registered in the C++ global dictionary with key "BsplinePaul_pbc"
#   - mmff.natoms holds the total number of atoms of the substrate plus PTCDA.
#   - mmff.AtomType is a pointer or array of atom types (or use mmff.params.atypes if that’s your interface)
#   - mmff.fapos holds the geometry (atomic positions)

# It might be that the name "BsplinePaul_pbc" is not registered; adjust the key if needed.
# grid_array = mmff.getArrayPointer("BsplinePaul_pbc")
# if grid_array is None:
#     raise RuntimeError('Grid array "BsplinePaul_pbc" not found. Ensure the grid is built.')

# # Save to XSF file.
# ret = mmff.saveXSF_geometry("rigid_scan.xsf", grid_array, pitch=1, offset=0, 
#                               natoms=mmff.natoms,
#                               atypes=mmff.params.atypes if hasattr(mmff, "params") else None,
#                               fapos=mmff.fapos, bPrimCoord=True)
# if ret != 0:
#     print("XSF file saved successfully!")
# else:
#     print("Error saving XSF file.")

# # Optionally, run the dynamics or further analysis
# mmff.run()



# Save the grid and geometry to an XSF file:
# mmff.saveXSF("rigid_scan.xsf", mmff.getArrayPointer("BsplinePaul_pbc"), 1, 0, mmff.natoms, mmff.AtomType, mmff.fapos)
# mmff.saveXSF("rigid_scan.xsf", mmff.getArrayPointer("BsplinePaul_pbc"))

# E = 0.0
# mmff.run(omp=True, nstepMax=20)
# print( "E=", E )
#for i in range(200):
#mmff.addSnapshot()
# with open("gopt_trajectory.xyz", "w") as file:pass
# mmff.printDatabase()

#mmff.init( xyz_name="dataxyz//pyridine", surf_name="dataxyz/NaCl_1x1_L2", bMMFF=False, gridStep=-1 )  # without gridFF
#mmff.getBuffs()
#mmff.eval()
#mmff.relax(1000)
#print( "FORCES:\n mmff.fapos:\n ", mmff.fapos )
#mmff.plot(bForce=True, Fscale=10.0 )
#plt.show()
# exit(0)


# '''
# mmff.init()
# #mmff.init_params( "data/AtomTypes.dat", "data/BondTypes.dat", "data/AngleTypes.dat" )
# #mmff.insertSMILES("CC");
# #mmff.insertSMILES("C=C");
# #mmff.insertSMILES("C#C");
# #mmff.insertSMILES("C#CCN=C", True );
# #mmff.insertSMILES("C1#CCN=C1", True );
# #mmff.insertSMILES("C=C1NC#CC1CO", True, True );

# #mmff.initWithSMILES( "C=C1NC#CC1CO" )
# mmff.initWithSMILES( "C=C" )
# mmff.getBuffs()
# mmff.relax(1000)
# mmff.plot()
# plt.show()

# exit()
# '''

# '''
# # ======== Oritent Molecule
# xyzs,Zs,enames,qs = au.loadAtomsNP( "data/xyz/Benzene_deriv.xyz" )
# au.orient( 2, (5,2), (1,3), xyzs, bFlipXZ=True )
# au.saveXYZ( enames, xyzs, "data/xyz/Benzene_deriv_.xyz", qs=qs, Rs=None )
# plt.plot( xyzs[:,0],xyzs[:,1], "o" )
# plt.axis('equal')
# plt.show()
# exit()
# '''

# '''
# # ============== C2H4,xyz
# #mmff.initWithMolFile( "C2H4.xyz", bNonBonded=False, bOptimizer=True)
# #mmff.printBuffNames()
# #mmff.getBuffs() #;print( mmff.ndims )
# #mmff.eval()
# #mmff.relax(1000, bWriteTrj=True )
# #Es=mmff.scanRotation( [1,4,5], 0, 0,1, np.pi*2, 100, bWriteTrj=True)   ;print("Es=", Es)
# #plt.plot(Es)
# #print( "Es(Etot,Eb,Ea,Eps,EppT,EppI):", mmff.Es )
# #nsel = mmff.splitAtBond(6-1)  ;print( "split to:\n", mmff.selection[:nsel],"\n", mmff.selection[nsel:] )
# '''



# # ============== Benzene_deriv.xyz
# #mmff.initWithMolFile( "data/xyz/Benzene_deriv.xyz", bNonBonded=False, bOptimizer=True)


# # mmff.initWithMolFile( "data/xyz/PTCDA.xyz", bNonBonded=False, bOptimizer=True)
# # mmff.getBuffs() #;print( mmff.ndims )

# #nsel = mmff.splitAtBond(5)  ;print( "split to:\n", mmff.selection[:nsel],"\n", mmff.selection[nsel:] )
# #nsel = mmff.splitAtBond(6)  ;print( "split to:\n", mmff.selection[:nsel],"\n", mmff.selection[nsel:] )
# #nsel = mmff.splitAtBond(10)  ;print( "split to:\n", mmff.selection[:nsel],"\n", mmff.selection[nsel:] )
# #nsel = mmff.splitAtBond(2)  ;print( "split to:\n", mmff.selection[:nsel],"\n", mmff.selection[nsel:] )
# #nsel = mmff.splitAtBond(4)  ;print( "split to:\n", mmff.selection[:nsel],"\n", mmff.selection[nsel:] )
# #print( "nsel ", nsel, len(mmff.selection)-nsel )
# #Es=mmff.scanRotation( 1, 1,11, np.pi*2, 100, bWriteTrj=True, _0=1) ;plt.plot(Es) ;print("Es=", Es)

# #Es = mmff.scanBondRotation( 6, np.pi*2, 100, bWriteTrj=True );  plt.plot(Es)
# #Es = mmff.scanBondRotation( 2, np.pi*2, 100, bWriteTrj=True );  plt.plot(Es); plt.grid()


# #mmff.eval()
# #mmff.relax(1000, Ftol=1e-4, bWriteTrj=True )
# #Es=mmff.scanRotation( 1, 1,11, np.pi*2, 100, sel=[11,13,14,20]+[29,30,31,32], bWriteTrj=True, _0=1) ;plt.plot(Es) ;print("Es=", Es)

# # plt.figure()
# # mmff.plot()
# # #mmff.plot_selection( mmff.selection[:nsel] )
# # #mmff.plot_selection( [1,2,3] )



# file_path = "/home/indranil/git/FireCore/cpp/common_resources/NaCl_8x8_L3_copy/Bspline_PLQd.npy"
# grid = np.load(file_path)

# print("Shape of grid:", grid.shape)  # expected (Nx, Ny, Nz, 3)

# channel = 0  # select the potential channel to view (e.g., Morse potential)

# # For this example we assume the grids share the same shape.
# Nx, Ny, Nz, _ = grid.shape

# # Use fixed indices (you can change these to whichever slice you want)
# mid_z = 10 + Nz//2  # for horizontal (xy) slice at constant z
# mid_x = 10 + Nx//2  # for vertical_x (yz) slice at constant x
# mid_y = 10 + Ny//2  # for vertical_y (xz) slice at constant y

# # Extract slices from the grid
# slice_xy = grid[:, :, mid_z, channel]
# slice_yz = grid[mid_x, :, :, channel]
# slice_xz = grid[:, mid_y, :, channel]


# # ---------------------------
# # Retrieve grid parameters and atom positions from the C++ extension
# gff_shift0, gff_pos0, gff_cell, gff_dCell, gff_natoms, gff_natoms_ = mmff.get_gridFF_info()
# # print("gff_pos0:", gff_pos0)
# print("gff_shift0:", gff_shift0)
# print("gff_cell:\n", gff_cell)
# print("gff_dCell:\n", gff_dCell)

# # For scaling purposes (we assume dCell is diagonal and defines the grid spacing)
# dx = gff_dCell[0, 0]
# dy = gff_dCell[1, 1]
# dz = gff_dCell[2, 2]

# # Compute physical extents for each slice.
# # Here our coordinate system will be shifted so that the grid origin becomes zero.
# extent_xy = [0, Nx * dx, 0, Ny * dy]     # x spans along axis0 and y axis1
# extent_yz = [0, Ny * dy, 0, Nz * dz]       # for a slice taken at constant x; horizontal axis: y, vertical: z
# extent_xz = [0, Nx * dx, 0, Nz * dz]       # for a slice taken at constant y; horizontal: x, vertical: z
# # extent_xy = [ gff_pos0[0], gff_pos0[0] + Nx * dx, gff_pos0[1], gff_pos0[1] + Ny * dy ]
# # extent_yz = [ gff_pos0[1], gff_pos0[1] + Ny * dy, gff_pos0[2], gff_pos0[2] + Nz * dz ]
# # extent_xz = [ gff_pos0[0], gff_pos0[0] + Nx * dx, gff_pos0[2], gff_pos0[2] + Nz * dz ]

# # ---------------------------
# # Get atom positions (substate and molecule)
# substrate_apos, molecule_apos = mmff.get_atom_positions()

# # --- Scaling and Shifting Atom Positions ---
# # The grid (and substrate) has a global origin gff_pos0.
# # For plotting in a coordinate system where the lower left is 0, we subtract gff_pos0.
# # Additionally, the molecule atoms are offset by gff_shift0.
# substrate_apos_scaled = substrate_apos #+ gff_pos0
# molecule_apos_scaled  = molecule_apos  - gff_shift0  #+ gff_pos0

# # ---------------------------
# # Plotting: Create a figure with three subplots for grid slices.
# fig, axs = plt.subplots(1, 3, figsize=(18, 6))

# # Horizontal slice (xy plane) at constant z
# im0 = axs[0].imshow(slice_xy.T, extent=extent_xy, origin="lower", cmap="viridis", aspect="auto")
# axs[0].set_title(f"Horizontal (xy) slice at z = {mid_z*dz}")
# axs[0].set_xlabel("x (physical dimesion)")
# axs[0].set_ylabel("y (physical dimension)")
# fig.colorbar(im0, ax=axs[0], fraction=0.046, pad=0.04)
# # Overlay atom positions (x-y)
# axs[0].scatter(substrate_apos_scaled[:, 0], substrate_apos_scaled[:, 1],
#                s=10, color='red', label='Substrate')
# axs[0].scatter(molecule_apos_scaled[:, 0], molecule_apos_scaled[:, 1],
#                s=10, color='blue', label='Molecule')
# axs[0].legend()

# # Vertical slice (yz plane) taken at constant x (mid_x)
# im1 = axs[1].imshow(slice_yz.T, extent=extent_yz, origin="lower", cmap="viridis", aspect="auto")
# axs[1].set_title(f"Vertical (yz) slice at x = {mid_x*dx}")
# axs[1].set_xlabel("y (physical dimension)")
# axs[1].set_ylabel("z (physical dimension)")
# fig.colorbar(im1, ax=axs[1], fraction=0.046, pad=0.04)
# # Overlay atom positions in y-z plane
# axs[1].scatter(substrate_apos_scaled[:, 1], substrate_apos_scaled[:, 2],
#                s=10, color='red')
# axs[1].scatter(molecule_apos_scaled[:, 1], molecule_apos_scaled[:, 2],
#                s=10, color='blue')

# # Vertical slice (xz plane) taken at constant y (mid_y)
# im2 = axs[2].imshow(slice_xz.T, extent=extent_xz, origin="lower", cmap="viridis", aspect="auto")
# axs[2].set_title(f"Vertical (xz) slice at y = {mid_y*dy}")
# axs[2].set_xlabel("x (physical dimesion)")
# axs[2].set_ylabel("z (physical dimension)")
# fig.colorbar(im2, ax=axs[2], fraction=0.046, pad=0.04)
# # Overlay atom positions in x-z plane
# axs[2].scatter(substrate_apos_scaled[:, 0], substrate_apos_scaled[:, 2],
#                s=10, color='red')
# axs[2].scatter(molecule_apos_scaled[:, 0], molecule_apos_scaled[:, 2],
#                s=10, color='blue')

# plt.tight_layout()
# plt.show()

exit()
