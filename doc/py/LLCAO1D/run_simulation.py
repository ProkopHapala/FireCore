# run_simulation.py

import numpy as np
import matplotlib.pyplot as plt
from quantum_solver_1D import QuantumSolver1D

def plot_results(solver, energies, coefficients, num_mos_to_plot=4):
    """
    Plots the molecular orbitals and their corresponding energy levels.
    """
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(10, 7))

    # Define a spatial grid for plotting the wavefunctions
    x_min = solver.positions.min() - 5
    x_max = solver.positions.max() + 5
    x_grid = np.linspace(x_min, x_max, 500)

    num_mos_to_plot = min(num_mos_to_plot, solver.num_basis_functions)

    # Plot each molecular orbital
    for i in range(num_mos_to_plot):
        # Construct the i-th MO wavefunction on the grid
        psi = solver.get_molecular_orbital(x_grid, i)
        
        # Plot the wavefunction, shifted vertically by its energy for visualization
        ax.plot(x_grid, psi*0.5 + energies[i], label=f'MO {i}, E = {energies[i]:.3f} Ha')
        
        # Draw a dashed line indicating the energy level of the orbital
        ax.axhline(y=energies[i], color='gray', linestyle='--', linewidth=0.7)

    # Plot the positions of the nuclei as markers
    y_marker = energies.min() - 0.2
    ax.plot(solver.positions, [y_marker] * len(solver.positions), 
            'kv', markersize=10, label='Nuclei Positions', linestyle='none')

    ax.set_xlabel('Position (Bohr)')
    ax.set_ylabel('Energy (Hartree) / Wavefunction Amplitude')
    ax.set_title(f'Molecular Orbitals for a {solver.num_basis_functions}-Atom System')
    ax.legend()
    plt.show()

def main():
    """
    Main function to define the system, run the solver, and display results.
    """
    # --- 1. Define the System: A 1D chain of 4 "Hydrogen-like" atoms ---
    num_atoms = 4
    spacing = 2.0  # Separation between atoms in Bohr
    
    # Center the atoms around the origin
    positions = np.linspace(-(num_atoms - 1) * spacing / 2, (num_atoms - 1) * spacing / 2, num_atoms)
    
    atoms_config = {
        'positions': positions,
        'charges': np.array([1.0] * num_atoms),       # All nuclear charges are +1
        'exponents': np.array([0.5] * num_atoms)      # Gaussian exponent for all basis functions
    }

    # --- 2. Initialize and Run the Solver ---
    solver = QuantumSolver1D(atoms_config)
    energies, coefficients = solver.solve()

    # --- 3. Print and Plot Results ---
    print("--- 1D Quantum Solver Results ---")
    print(f"System: {num_atoms} atoms with spacing {spacing} Bohr")
    print(f"Atomic Positions: {atoms_config['positions']}")
    print("\nCalculated Orbital Energies (Hartree):")
    for i, E in enumerate(energies):
        print(f"  E{i}: {E:.6f}")

    print("\nMO Coefficients (Eigenvectors, columns of the matrix):")
    np.set_printoptions(precision=4, suppress=True)
    print(coefficients)

    # Plot the resulting molecular orbitals
    plot_results(solver, energies, coefficients, num_mos_to_plot=4)

if __name__ == '__main__':
    main()