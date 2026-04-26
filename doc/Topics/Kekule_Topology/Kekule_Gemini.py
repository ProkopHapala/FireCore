import numpy as np
import matplotlib.pyplot as plt

def solve_kekule_boundary():
    """
    DIDACTIC DERIVATION:
    1. Geometry: We use a honeycomb lattice (graphene).
    2. Physics: Hopping (t) represents the resonance integral (beta in Huckel).
    3. Perturbation: We vary t slightly (+delta or -delta). 
       - If the pattern shifts across a boundary, we create a 'Phase Slip'.
    4. Topology: This slip forces the Hamiltonian to close the gap at the boundary,
       creating a Jackiw-Rebbi soliton (a zero-energy conductive state).
    """
    
    # --- 1. Hyperparameters ---
    width = 10   # atoms across
    length = 40  # atoms long
    t0 = 2.7     # Base hopping (eV)
    delta = 0.4  # Strength of Kekule distortion (the 'Mass')
    
    # Generate coordinates for a standard Armchair Nanoribbon
    coords = []
    for x in range(length):
        for y in range(width):
            # Hexagonal offsets
            off = 0.3 if x % 2 == 0 else 0
            coords.append([x * 0.866, y + (0.5 if x % 2 else 0)])
    
    coords = np.array(coords)
    num_atoms = len(coords)
    H = np.zeros((num_atoms, num_atoms))
    
    # --- 2. Building the Hamiltonian with Domain Boundary ---
    # We define a boundary at length // 2
    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):
            dist = np.linalg.norm(coords[i] - coords[j])
            
            # Standard nearest-neighbor distance is ~1.0 in these units
            if 0.9 < dist < 1.1:
                # Determine which domain we are in
                mid_x = (coords[i][0] + coords[j][0]) / 2
                
                # Domain 1 (Left): Pattern A
                if mid_x < (length * 0.866) / 2:
                    # Logic: Arbitrarily strengthen bonds based on index parity
                    # This creates a specific 'Kekule' phase
                    t = t0 + delta if (i + j) % 4 == 0 else t0 - delta/2
                
                # Domain 2 (Right): Pattern B (The 'Phase Slip')
                else:
                    # We flip the logic! What was weak is now strong.
                    # This represents a grain boundary in the bond-order.
                    t = t0 + delta if (i + j) % 4 != 0 else t0 - delta/2
                
                H[i, j] = -t
                H[j, i] = -t

    # --- 3. Solving the System ---
    energies, wavefunctions = np.linalg.eigh(H)
    
    # The 'Soliton' or conductive state is the one closest to E=0 (Fermi Level)
    # In a perfect mid-gap system, this is the num_atoms // 2 state
    mid_idx = num_atoms // 2
    soliton_state = wavefunctions[:, mid_idx]
    probability_density = np.abs(soliton_state)**2
    
    # --- 4. Visualization ---
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))
    
    # Plot A: The Spectrum (DOS)
    ax1.scatter(range(len(energies)), energies, s=5, color='black')
    ax1.axhline(0, color='red', linestyle='--', alpha=0.5, label='Fermi Level')
    ax1.set_title("Energy Spectrum: Note the states at E=0 (The Solitons)")
    ax1.set_ylabel("Energy (eV)")
    ax1.legend()

    # Plot B: The LDOS Map (The 'Conductive Channel')
    scatter = ax2.scatter(coords[:, 0], coords[:, 1], c=probability_density, 
                          cmap='inferno', s=100, edgecolors='none')
    ax2.set_title("LDOS at E=0: Visualization of the Conductive Grain Boundary")
    ax2.set_aspect('equal')
    plt.colorbar(scatter, ax=ax2, label='Probability Density $|\psi|^2$')
    
    print(f"Gap size in bulk: ~{2*delta} eV")
    print(f"Found {np.sum(np.isclose(energies, 0, atol=1e-2))} state(s) at the Fermi level.")
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    solve_kekule_boundary()