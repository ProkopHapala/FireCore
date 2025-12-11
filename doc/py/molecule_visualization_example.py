import numpy as np
import matplotlib.pyplot as plt
from plot_utils import MoleculeTrajectoryVisualizer

def main():
    """
    Example script demonstrating how to use the MoleculeTrajectoryVisualizer class
    to visualize molecular trajectories and structures.
    """
    # Create a visualizer instance
    visualizer = MoleculeTrajectoryVisualizer()
    
    # Path to your XYZ trajectory file
    # Replace with the actual path to your trajectory file
    trajectory_file = "/home/indranil/Documents/Project_1/FireCore/UFF/new_trial_relax_scan_ptcda_test_trajectory.xyz"
    
    # Read the trajectory file
    visualizer.read_xyz_trajectory(trajectory_file)
    
    # Example 1: Plot trajectories of specific atoms (atoms 26 and 29 as mentioned in the example)
    atom_indices = [26, 29]  # 0-based indexing
    fig1, axes1 = visualizer.plot_atom_trajectories(
        atom_indices=atom_indices,
        projections=['xy', 'xz', 'yz'],  # Show all three 2D projections in one row
        figsize=(15, 5),  # Wider figure for 1 row, 3 columns layout
        colors=['red', 'blue'],
        title="Selected Atom Trajectories"  # Global title
    )
    
    plt.savefig("/home/indranil/Documents/Project_1/FireCore/UFF/selected_atom_trajectories.png", dpi=300, bbox_inches='tight')
    plt.show()
    plt.close(fig1)
    
    # Example 2: Visualize molecular structures at evenly spaced intervals
    # Select 5 evenly spaced frames from the trajectory
    n_structures = 3  # Number of structures to visualize
    frame_indices = np.linspace(0, visualizer.num_frames-1, n_structures, dtype=int)
    
    fig2, axes2 = visualizer.visualize_structures(
        frame_indices=frame_indices,
        projections=['xy', 'xz', 'yz'],
        figsize=(15, 5),  # Wider figure for 1 row, 3 columns layout
        bond_length_threshold=2.0,  # Maximum distance to consider atoms as bonded (in Angstroms)
        atom_size=10,  # Smaller atom size to reduce clutter
        title="Molecular Structures"  # Global title
    )
    
    plt.savefig("/home/indranil/Documents/Project_1/FireCore/UFF/molecular_structures.png", dpi=300, bbox_inches='tight')
    plt.show()
    plt.close(fig2)
    
    # Example 3: Combined visualization - both atom trajectories and molecular structures
    fig3, axes3 = visualizer.visualize_combined(
        atom_indices=atom_indices,
        frame_indices=frame_indices,
        projections=['xy', 'xz', 'yz'],
        figsize=(15, 5),  # Wider figure for 1 row, 3 columns layout
        bond_length_threshold=2.0,
        atom_size=10,  # Smaller atom size to reduce clutter
        trajectory_colors=['red', 'blue'],
        title="Combined Atom Trajectories and Molecular Structures"  # Global title
    )
    
    plt.savefig("/home/indranil/Documents/Project_1/FireCore/UFF/combined_visualization.png", dpi=300, bbox_inches='tight')
    plt.show()
    plt.close(fig3)
    
    print("Visualizations completed and saved as PNG files.")

if __name__ == "__main__":
    main()