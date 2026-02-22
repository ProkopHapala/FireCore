import numpy as np
import matplotlib.pyplot as plt

def plot_potential_decomposition(r, E_coul, E_vdw, E_hcorr, E_epair):
    """
    r: array of distances
    E_*: arrays of energy components
    """
    
    # 1. Define the cumulative sums
    base_energy = E_coul + E_vdw
    base_plus_H = base_energy + E_hcorr
    total_energy = base_plus_H + E_epair
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # 2. Plot the Total and Base as lines
    ax.plot(r, base_energy, 'k--', linewidth=1.5, label='Baseline (Morse+Coul)', alpha=0.6)
    ax.plot(r, total_energy, 'k-', linewidth=2.5, label='Total Energy')
    
    # 3. Fill the contributions
    # Contribution of H-bond corrections (Difference between Base and Base+H)
    ax.fill_between(r, base_energy, base_plus_H, color='cyan', alpha=0.4, label='H-bond Corr (H1/H2)')
    
    # Contribution of Epairs/SigmaHoles (Difference between Base+H and Total)
    ax.fill_between(r, base_plus_H, total_energy,  color='orange', alpha=0.4, label='E-pair/Sigma (SR terms)')

    # Optional: Plot the components individually as thin lines to see shape
    # ax.plot(r, E_hcorr, 'c:', linewidth=1)
    # ax.plot(r, E_epair, 'orange', linestyle=':', linewidth=1)

    # Styling
    ax.set_ylim(-10.0, 10.0) # Set reasonable limits to ignore extreme repulsion
    ax.axhline(0, color='gray', linewidth=0.5)
    ax.set_xlabel('Distance [Å]')
    ax.set_ylabel('Energy [kcal/mol]')
    ax.set_title('Decomposition of Interaction Energy')
    ax.legend()
    
    return fig, ax

# --- Dummy Data Generation for testing ---
r = np.linspace(1.5, 5.0, 100)
# Fake Morse
E_vdw = 0.1 * ((1 - np.exp(-1.5*(r-3.0)))**2 - 1) 
# Fake Coulomb (slightly attractive)
E_coul = -3.0 * 0.2 / r 
# Fake H-bond (short range attraction)
E_hcorr = -2.0 * np.exp(-2.0 * (r - 2.8)**2)
# Fake Epair (angular/short range)
E_epair = -1.5 * np.exp(-4.0 * (r - 2.5)**2)

plot_potential_decomposition(r, E_coul, E_vdw, E_hcorr, E_epair)
plt.show()