import numpy as np
import matplotlib.pyplot as plt
# Import from the new file
from quantum_solver_1D import (
    QuantumSolver1D,
    gaussian_overlap_wrapper,
    gaussian_kinetic_wrapper,
    gaussian_coulomb_wrapper,
    gaussian_basis_eval_wrapper
)

# trun of line-brakes when printing numpy arrays
np.set_printoptions(linewidth=1000)

def plot_results(solver, energies, coefficients, num_mos_to_plot=4, is_localized=False):
    """
    Plots the molecular orbitals and their corresponding energy levels.
    """
    fig, ax = plt.subplots(figsize=(10, 7)) # Create a new figure for each plot

    # Define a spatial grid for plotting the wavefunctions
    x_min = solver.basis_centers.min() - 5 # Use basis_centers or nuclear_positions
    x_max = solver.basis_centers.max() + 5
    x_grid = np.linspace(x_min, x_max, 1000) # Increased grid points for smoother plot
    num_mos_to_plot = min(num_mos_to_plot, solver.num_basis_functions)

    # Plot each molecular orbital
    for i in range(num_mos_to_plot):
        # Construct the i-th MO wavefunction on the grid
        if is_localized:
            psi = solver.get_localized_molecular_orbital(x_grid, i)
        else:
            psi = solver.get_molecular_orbital(x_grid, i)

        # Plot the wavefunction, shifted vertically by its energy for visualization
        # Scale psi for better visualization if needed (e.g., psi * 0.5)
        ax.plot(x_grid, psi * 0.5 + energies[i], label=f'MO {i}, E = {energies[i]:.3f} Ha')

        # Draw a dashed line indicating the energy level of the orbital
        ax.axhline(y=energies[i], color='gray', linestyle='--', linewidth=0.7)

    # Plot the positions of the nuclei as markers
    y_marker = energies.min() - 0.2 if energies.size > 0 else -0.2
    ax.plot(solver.nuclear_positions, [y_marker] * len(solver.nuclear_positions),
            'kv', markersize=8, label='Nuclei Positions', linestyle='none')
    ax.set_xlabel('Position (a.u.)') # Use atomic units (Bohr)
    title_suffix = " (Localized)" if is_localized else " (Standard)"
    ax.set_ylabel('Energy (Hartree) / Wavefunction Amplitude')
    ax.set_title(f'Molecular Orbitals for a {solver.num_basis_functions}-Atom System{title_suffix}')
    ax.legend()
    plt.show()


def main():
    num_atoms = 20 # Increased atoms for a slightly larger system
    spacing   = 1.5 # Reduced spacing to increase overlap

    nuclear_positions = np.linspace(-(num_atoms - 1) * spacing / 2, (num_atoms - 1) * spacing / 2, num_atoms)
    nuclear_charges = np.array([1.0] * num_atoms)

    # Define widths 'w' for Gaussian basis functions.
    # Example: if you want alpha = 0.5, then w = sqrt(1/0.5) = sqrt(2) ~ 1.414
    # Let's use w = 1.0 for all atoms for this example
    basis_widths_array = np.array([1.0] * num_atoms) # Smaller width for more localized basis

    atoms_config = {
        'positions': nuclear_positions, # For nuclear charges and basis centers
        'charges': nuclear_charges,
        'widths': basis_widths_array    # For basis functions
    }

    solver = QuantumSolver1D(
        atoms_config,
        overlap_fn=gaussian_overlap_wrapper,
        kinetic_fn=gaussian_kinetic_wrapper,
        coulomb_fn=gaussian_coulomb_wrapper,
        basis_eval_fn=gaussian_basis_eval_wrapper
    )

    # --- Localized-solver parameter setup ---
    # Adjust these knobs from the script without touching the library
    solver.loc_step_size   = 0.02       # gradient step
    solver.ortho_damp      = 0.5        # orthogonalisation damping
    solver.ortho_iter      = 5
    solver.loc_max_iter    = 1000
    solver.loc_coeff_tol   = 1e-7
    solver.loc_overlap_tol = 1e-5

    # localisation controls
    solver.alpha_loc         = 0.05      # quadratic potential strength (0 → none)
    solver.apply_hard_cutoff = False    # disable hard cutoff for delocalised test

    # --- Running Localized Solver ---
    print("\n--- Running Localized Solver ---")
    
    # Localize each MO around its corresponding atom/basis center
    localization_centers = solver.basis_centers.copy()

    # Choose a cutoff slightly larger than half the spacing, but smaller than the spacing
    # to encourage localization but allow some overlap.
    # to encourage localization but allow some overlap.
    localization_cutoff = spacing * 8.0 # Example cutoff distance
    print(f"Using localization cutoff R_loc = {localization_cutoff:.2f} a.u.")

    localized_energies, localized_coefficients = solver.solve_localized(
        localization_centers=localization_centers,
        localization_cutoff=localization_cutoff
    )

    if localized_energies is not None and localized_coefficients is not None:
        print("\n--- Localized Solver Results ---")
        print("\nCalculated Localized Orbital Energies (Hartree):")
        for i, E in enumerate(localized_energies):
            print(f"  E{i}: {E:.6f}")
        print("\nLocalized MO Coefficients (columns of the matrix):")
        # print(localized_coefficients)  # Keep commented out unless needed for debugging

        plot_results(solver, localized_energies, localized_coefficients, num_mos_to_plot=num_atoms, is_localized=True)
    else:
        print("Localized solver failed to converge or produced invalid results.")

if __name__ == '__main__':
    main()