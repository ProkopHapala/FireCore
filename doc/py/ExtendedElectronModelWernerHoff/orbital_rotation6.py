import numpy as np
import matplotlib.pyplot as plt
from pyscf import gto, scf

def calculate_pauli_energy_structure():
    # 1. Define Ethylene Geometry (D2h Symmetry)
    # Using Angstroms
    cc_half = 1.34 / 2
    ch_len = 1.09
    
    # H positions (120 deg geometry)
    hx = cc_half + ch_len * 0.5   # cos(60)
    hy = ch_len * 0.866           # sin(60)
    
    mol_str = f"""
    C  {-cc_half}  0.0  0.0
    C  {cc_half}   0.0  0.0
    H  {-hx}   {hy}  0.0
    H  {-hx}  {-hy}  0.0
    H  {hx}    {hy}  0.0
    H  {hx}   {-hy}  0.0
    """
    
    # 2. Run DFT (Symmetry Enforced for cleanliness)
    print("Computing Orbitals (Symmetry Enabled)...")
    mol = gto.M(atom=mol_str, basis='6-31g', verbose=0, symmetry=True)
    mf = scf.RHF(mol)
    mf.kernel()
    
    # 3. Setup Grid (XY Plane)
    res = 200
    lim = 3.5
    x = np.linspace(-lim, lim, res)
    y = np.linspace(-lim, lim, res)
    X, Y = np.meshgrid(x, y)
    
    # Z-slice: We need to be slightly off-plane (z=0.1) 
    # because at z=0 exactly, the pi-orbitals are zero, causing division by zero issues.
    coords = np.zeros((res**2, 3))
    coords[:, 0] = X.ravel()
    coords[:, 1] = Y.ravel()
    coords[:, 2] = 0.15 
    
    # 4. Evaluate Orbitals and Gradients
    print("Calculating Kinetic Energy Densities...")
    ao_vals = mol.eval_gto("GTOval_sph_deriv1", coords)
    mo_coeffs = mf.mo_coeff
    n_occ = mol.nelectron // 2
    
    # Accumulators
    rho = np.zeros(res**2)
    grad_rho = np.zeros((3, res**2))
    tau_total = np.zeros(res**2)
    
    # Sum over occupied orbitals
    for i in range(n_occ):
        c = mo_coeffs[:, i]
        
        # Psi and Gradient
        psi = np.dot(ao_vals[0], c)
        dx  = np.dot(ao_vals[1], c)
        dy  = np.dot(ao_vals[2], c)
        dz  = np.dot(ao_vals[3], c)
        
        # Density (Factor of 2 for RHF spin pairing)
        rho_i = 2.0 * psi**2
        rho += rho_i
        
        # Gradient of Density: Sum(2 * 2 * psi * grad_psi)
        grad_rho[0] += 4.0 * psi * dx
        grad_rho[1] += 4.0 * psi * dy
        grad_rho[2] += 4.0 * psi * dz
        
        # Total Kinetic Energy Density: Sum(2 * 0.5 * |grad_psi|^2)
        # Factor 2 for spin, 0.5 for KE operator
        tau_i = (dx**2 + dy**2 + dz**2)
        tau_total += tau_i

    # 5. Calculate Pauli Kinetic Energy
    # T_Pauli = T_Total - T_VonWeizsacker
    # T_VW = (1/8) * |grad_rho|^2 / rho
    
    grad_rho_sq = grad_rho[0]**2 + grad_rho[1]**2 + grad_rho[2]**2
    
    # Avoid division by zero
    mask = rho > 1e-6
    tau_vw = np.zeros_like(rho)
    tau_vw[mask] = (1.0/8.0) * grad_rho_sq[mask] / rho[mask]
    
    tau_pauli = tau_total - tau_vw
    # Numerical noise can make it slightly negative, clamp to 0
    tau_pauli = np.maximum(tau_pauli, 0)

    # Reshape
    Rho_Grid = rho.reshape((res, res))
    Tau_P_Grid = tau_pauli.reshape((res, res))
    
    return X, Y, Rho_Grid, Tau_P_Grid, mol

# --- PLOTTING ---
X, Y, Rho, Tau_P, mol = calculate_pauli_energy_structure()

fig, ax = plt.subplots(1, 2, figsize=(16, 7), facecolor='#111111')

# Plot 1: Total Density (The "Blob")
ax[0].set_facecolor('black')
im0 = ax[0].imshow(Rho, extent=[-3.5,3.5,-3.5,3.5], origin='lower', cmap='bone')
ax[0].set_title("Standard Electron Density ($\\rho$)", color='white', fontsize=14)
ax[0].contour(X, Y, Rho, levels=10, colors='cyan', linewidths=0.5, alpha=0.3)

# Plot 2: Pauli Energy (The "Hidden Structure")
# We use Log scale because the core potentials are huge compared to valence bonds
from matplotlib.colors import LogNorm

ax[1].set_facecolor('black')
# Adding small epsilon for log scale
im1 = ax[1].imshow(Tau_P + 1e-3, extent=[-3.5,3.5,-3.5,3.5], origin='lower', 
                   cmap='inferno', norm=LogNorm(vmin=0.01, vmax=Tau_P.max()))

ax[1].set_title("Pauli Kinetic Energy Density ($\\tau_P$)\n(Hofer's Rotational Stress)", color='white', fontsize=14)

# Overlay Atoms
atom_coords = mol.atom_coords()
for ax_i in ax:
    ax_i.scatter(atom_coords[:,0], atom_coords[:,1], color='cyan', s=100, edgecolors='white', zorder=10)
    # Draw Bonds
    ax_i.plot([atom_coords[0,0], atom_coords[1,0]], [0,0], 'w-', alpha=0.3) 

# Highlight the Filaments (Ridges of Pauli Energy)
# We plot contours of high Pauli energy
ax[1].contour(X, Y, Tau_P, levels=[0.1, 0.5, 1.0, 2.0], colors='white', linewidths=0.5, alpha=0.5)

plt.tight_layout()
plt.show()