import numpy as np
import matplotlib.pyplot as plt
from pyscf import gto, scf, dft
import matplotlib.cm as cm

def calculate_molecule(mol_string, basis='sto-3g'):
    """
    Builds molecule and calculates Electron Density and Orbitals
    using Standard DFT (PySCF) as the ground truth.
    """
    mol = gto.M(atom=mol_string, basis=basis, symmetry=True)
    
    # Run simple RHF (Restricted Hartree Fock) to get orbitals
    # We use HF/DFT because EEM claims to reproduce these densities.
    mf = scf.RHF(mol)
    mf.verbose = 0
    mf.kernel()
    
    return mol, mf

def get_grid_data(mol, mf, orbital_index, slice_axis='z', slice_pos=0.0, range_val=4.0, resolution=100):
    """
    Evaluates the orbital and its gradients on a 2D slice to visualize 
    the 'Hydrodynamic' features (Vortices/Nodes).
    """
    # Create a 2D grid
    x = np.linspace(-range_val, range_val, resolution)
    y = np.linspace(-range_val, range_val, resolution)
    X, Y = np.meshgrid(x, y)
    
    # Construct 3D coordinates for the slice
    coords = np.zeros((resolution**2, 3))
    
    if slice_axis == 'z': # Cross section XY
        coords[:, 0] = X.ravel()
        coords[:, 1] = Y.ravel()
        coords[:, 2] = slice_pos
    elif slice_axis == 'x': # Cross section YZ
        coords[:, 0] = slice_pos
        coords[:, 1] = X.ravel()
        coords[:, 2] = Y.ravel()
        
    # Evaluate the specific Orbital (Wavefunction Psi)
    # NOTE: we need both values and first derivatives; "*_deriv1" returns
    # (4, ngrid, nao) -> [value, dx, dy, dz].
    ao = mol.eval_gto("GTOval_sph_deriv1", coords)
    coeff = mf.mo_coeff[:, orbital_index]
    
    # Psi (Scalar field)
    psi = ao[0].dot(coeff).reshape(resolution, resolution)
    
    # Gradient of Psi (Vector field -> The "Stress" or "Twist" direction)
    # This corresponds to the Hofer bivector orientation change
    grad_x = ao[1].dot(coeff).reshape(resolution, resolution)  # DEBUG: grad components from deriv1
    grad_y = ao[2].dot(coeff).reshape(resolution, resolution)
    grad_z = ao[3].dot(coeff).reshape(resolution, resolution)
    
    # --- The "Hofer/Hydrodynamic" Quantities ---
    
    # 1. Mass Density
    rho = psi**2
    
    # 2. The "Quantum Potential" / Kinetic Energy Density
    # In Hofer's model, KE is the energy of the twist.
    # KE_density ~ (grad_psi)^2
    ke_density = 0.5 * (grad_x**2 + grad_y**2 + grad_z**2)
    
    # 3. The "Centrifugal Force" (Bohm Potential Gradient)
    # This force pushes density away from the node.
    # We visualize the gradient vector field (grad_x, grad_y)
    
    return X, Y, psi, rho, ke_density, grad_x, grad_y

def plot_ethylene_pi_bond():
    """
    Visualizes C2H4 (Ethylene) Pi-Bond as a Vortex Structure
    """
    # Define Ethylene (C2H4) aligned on axes
    # C=C bond along Z axis for this visualization setup
    mol_str = """
    C 0.0 0.0  0.66
    C 0.0 0.0 -0.66
    H 0.0 0.92  1.23
    H 0.0 -0.92 1.23
    H 0.0 0.92 -1.23
    H 0.0 -0.92 -1.23
    """
    
    print("Calculating Ethylene (C2H4) electronic structure...")
    mol, mf = calculate_molecule(mol_str)
    
    # Identify HOMO (The Pi Bond)
    # In C2H4 with STO-3G, HOMO is usually index 7 (starting 0)
    homo_idx = mf.mo_occ.argmax() + 1 # A guess, finding the Pi orbital
    # Actually, let's just find the orbital with a node along the bond axis
    # For this specific geometry, the Pi bond (px) node is the YZ plane (x=0)
    # Let's visualize the XY plane slice to see the two lobes and the node between them.
    
    # NOTE: Depending on basis set, orbital ordering changes. 
    # For STO-3G C2H4:
    # 0-4: Core/Sigma
    # The Pi orbital is usually the highest occupied.
    orbital_to_plot = 7 # HOMO for C2H4 in STO-3G usually
    
    print(f"Visualizing Orbital #{orbital_to_plot} (Likely HOMO Pi-Bond)")

    # --- Slice 1: The 'Side View' showing the Node (XZ Plane) ---
    # We slice perpendicular to the molecule to see the "Dumbbell" cross section
    # The Pi bond lobes are above and below the plane.
    X, Y, psi, rho, ke, gx, gy = get_grid_data(mol, mf, orbital_to_plot, 
                                              slice_axis='z', slice_pos=0.0, range_val=3.0)

    fig, ax = plt.subplots(1, 3, figsize=(18, 5))
    
    # Plot 1: The Mass Density (Standard Chemistry View)
    ax[0].set_title("1. Mass Density (rho)\nStandard View")
    cf = ax[0].contourf(X, Y, rho, 20, cmap='inferno')
    plt.colorbar(cf, ax=ax[0])
    ax[0].text(0,0, "Node (Zero Density)", color='white', ha='center')
    
    # Plot 2: The "Rotational Stress" (Kinetic Energy Density)
    # In Hofer's model, this is the energy cost of the twist.
    # Notice it is HIGHEST exactly where the density is LOWEST (The Node).
    ax[1].set_title("2. 'Twist' Energy (Kinetic Density)\nHofer's Rotational Energy")
    cf2 = ax[1].contourf(X, Y, ke, 20, cmap='viridis')
    plt.colorbar(cf2, ax=ax[1])
    ax[1].scatter([0], [0], color='red', marker='x', label='Vortex Core')
    ax[1].legend()
    
    # Plot 3: The "Centrifugal Field" (Gradient Vectors)
    # This shows the "Mechanical Polarization". 
    # The vectors point violently away/across the node, representing the 
    # rapid rotation required to flip the bivector orientation.
    ax[2].set_title("3. Mechanical Polarization (Field Stress)\nThe 'Centrifugal' Pauli Force")
    
    # Normalize vectors for cleaner streamplot
    speed = np.sqrt(gx**2 + gy**2)
    lw = 2 * speed / speed.max()
    
    # Plot streamlines of the gradient
    strm = ax[2].streamplot(X, Y, gx, gy, color=speed, cmap='autumn', density=1.5)
    ax[2].contour(X, Y, psi, levels=[0], colors='black', linewidths=2, linestyles='--')
    ax[2].text(0,0.2, "Node/Vortex Wall", ha='center')
    
    plt.tight_layout()
    plt.show()

    # --- BONUS: Visualizing the "Spin Vortex" ---
    # To see the "Rotation" explicitly, we need a complex wavefunction.
    # In EEM, a standing wave (real orbital) is 2 counter-rotating waves.
    # Let's synthesize the "Forward Moving" component: Psi_complex = Psi + i * Grad_Psi
    
    print("Generating Vortex Filament Visualization...")
    
    fig2, ax2 = plt.subplots(figsize=(8, 8))
    ax2.set_title("The 'Vortex Filament' (Phase Circulation)\nHow Hofer's model sees the Node")
    
    # Synthesize a "Traveling Wave" version of the Pi-orbital 
    # to reveal the hidden topological winding number.
    # We add a phase gradient perpendicular to the node.
    complex_psi = psi + 1j * (gx * 0.5) 
    
    # Get phase
    phase = np.angle(complex_psi)
    
    # Plot Phase
    cp = ax2.contourf(X, Y, phase, levels=20, cmap='hsv') # HSV is cyclic (good for phase)
    plt.colorbar(cp, ax=ax2, label='Rotor Phase Angle')
    
    # Plot Density contour on top
    ax2.contour(X, Y, rho, levels=5, colors='black', alpha=0.5)
    
    ax2.set_xlabel("Space X")
    ax2.set_ylabel("Space Y")
    
    plt.show()

if __name__ == "__main__":
    plot_ethylene_pi_bond()