# In your run_simulation.py (or a new script)

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

# ... (plot_results function can remain similar, but ensure it uses solver.basis_centers or solver.nuclear_positions correctly)

def plot_results(solver, energies, coefficients, num_mos_to_plot=4):
    fig, ax = plt.subplots(figsize=(10, 7))
    x_min = solver.basis_centers.min() - 5 # Use basis_centers or nuclear_positions
    x_max = solver.basis_centers.max() + 5
    x_grid = np.linspace(x_min, x_max, 500)
    num_mos_to_plot = min(num_mos_to_plot, solver.num_basis_functions)

    for i in range(num_mos_to_plot):
        psi = solver.get_molecular_orbital(x_grid, i)
        ax.plot(x_grid, psi * 0.5 + energies[i], label=f'MO {i}, E = {energies[i]:.3f} Ha')
        ax.axhline(y=energies[i], color='gray', linestyle='--', linewidth=0.7)

    y_marker = energies.min() - 0.2 if energies.size > 0 else -0.2
    ax.plot(solver.nuclear_positions, [y_marker] * len(solver.nuclear_positions),
            'kv', markersize=10, label='Nuclei Positions', linestyle='none')
    ax.set_xlabel('Position (Bohr)')
    ax.set_ylabel('Energy (Hartree) / Wavefunction Amplitude')
    ax.set_title(f'Molecular Orbitals for a {solver.num_basis_functions}-Atom System')
    ax.legend()
    plt.show()


def main():
    num_atoms = 4
    spacing = 2.0
    
    nuclear_positions = np.linspace(-(num_atoms - 1) * spacing / 2, (num_atoms - 1) * spacing / 2, num_atoms)
    nuclear_charges = np.array([1.0] * num_atoms)
    
    # Define widths 'w' for Gaussian basis functions.
    # Example: if you want alpha = 0.5, then w = sqrt(1/0.5) = sqrt(2) ~ 1.414
    # Let's use w = 1.5 for all atoms for this example
    basis_widths_array = np.array([1.5] * num_atoms) 

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
    energies, coefficients = solver.solve()

    print("--- 1D Quantum Solver Results (New Callback Version) ---")
    print(f"System: {num_atoms} atoms with spacing {spacing} Bohr")
    print(f"Nuclear Positions: {solver.nuclear_positions}")
    print(f"Basis Function Widths: {solver.basis_widths}")
    print("\nCalculated Orbital Energies (Hartree):")
    for i, E in enumerate(energies):
        print(f"  E{i}: {E:.6f}")

    print("\nMO Coefficients (Eigenvectors, columns of the matrix):")
    np.set_printoptions(precision=4, suppress=True)
    print(coefficients)

    plot_results(solver, energies, coefficients, num_mos_to_plot=min(4, num_atoms))

if __name__ == '__main__':
    main()
