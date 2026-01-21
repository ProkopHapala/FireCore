import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from pyscf import gto, scf
import matplotlib.colors as mcolors

# --- 1. System Setup (Acetylene) ---
def setup_acetylene():
    # C-C bond along Z-axis
    mol = gto.M(
        atom="""
        C 0.0 0.0 0.6
        C 0.0 0.0 -0.6
        H 0.0 0.0 1.66
        H 0.0 0.0 -1.66
        """,
        basis='sto-3g',
        symmetry=False
    )
    mf = scf.RHF(mol)
    mf.verbose = 0
    mf.kernel()
    return mol, mf

# --- 2. Grid & Wavefunction Evaluator ---
def get_slice_grid(L=3.0, N=50, axis='xy', z_val=0.3):
    """
    Creates a 2D grid for plotting.
    z_val: The height at which we slice (crucial to see cross-sections)
    """
    pts = np.linspace(-L, L, N)
    u, v = np.meshgrid(pts, pts)
    
    coords = np.zeros((N*N, 3))
    if axis == 'xy': # Cross section perpendicular to bond
        coords[:, 0] = u.ravel()
        coords[:, 1] = v.ravel()
        coords[:, 2] = z_val
    elif axis == 'xz': # Cross section along the bond
        coords[:, 0] = u.ravel()
        coords[:, 1] = 0
        coords[:, 2] = v.ravel()
        
    return u, v, coords

def eval_wavefunction_dynamics(mol, mf, indices, coeffs, coords, t):
    """
    Computes Psi(r, t) = Sum c_n * phi_n(r) * exp(-i * E_n * t)
    Returns Density and Flux (Current Density)
    """
    n_points = coords.shape[0]
    
    # Get Atomic Orbitals and Gradients on grid
    # ao_val shape: (4, N_points, N_ao) -> [val, dx, dy, dz]
    ao_all = mol.eval_gto("GTOval_sph_deriv1", coords)
    
    # MO Coefficients and Energies
    mo_coeff = mf.mo_coeff
    mo_energies = mf.mo_energy
    
    # Initialize Psi (Complex) and Grad Psi (Complex Vector)
    psi_tot = np.zeros(n_points, dtype=complex)
    grad_psi_tot = np.zeros((3, n_points), dtype=complex)
    
    for idx, c in zip(indices, coeffs):
        # 1. Get Spatial Orbital phi_n
        # Expand MO in terms of AO
        coeff_vec = mo_coeff[:, idx]
        
        # Scalar value phi_n
        phi_n = np.dot(ao_all[0], coeff_vec)
        
        # Gradient vector (d/dx, d/dy, d/dz) phi_n
        grad_phi_n_x = np.dot(ao_all[1], coeff_vec)
        grad_phi_n_y = np.dot(ao_all[2], coeff_vec)
        grad_phi_n_z = np.dot(ao_all[3], coeff_vec)
        
        # 2. Time Evolution Factor
        # Energy E_n
        E = mo_energies[idx]
        phase = np.exp(-1j * E * t)
        
        # 3. Add to total Wavefunction
        term = c * phi_n * phase
        psi_tot += term
        
        # Add to total Gradient
        grad_psi_tot[0] += c * grad_phi_n_x * phase
        grad_psi_tot[1] += c * grad_phi_n_y * phase
        grad_psi_tot[2] += c * grad_phi_n_z * phase

    # --- Compute Observables ---
    
    # Probability Density: |Psi|^2
    density = np.abs(psi_tot)**2
    
    # Probability Flux (Current): j = Im( Psi* . Grad(Psi) )
    # (Using atomic units, h_bar/m = 1)
    # j_x = Im( conj(Psi) * dPsi/dx )
    jx = np.imag(np.conj(psi_tot) * grad_psi_tot[0])
    jy = np.imag(np.conj(psi_tot) * grad_psi_tot[1])
    jz = np.imag(np.conj(psi_tot) * grad_psi_tot[2])
    
    return density, jx, jy, jz, psi_tot

# --- 3. Animation Routine ---

def visualize_orbitals(mol, mf):
    """
    Helper to see which orbitals are which (energies)
    """
    print("\n--- Molecular Orbitals ---")
    for i, e in enumerate(mf.mo_energy):
        occ = "Occ" if i < mol.nelectron//2 else "Vir"
        print(f"MO {i}: E = {e:.4f} Ha ({occ})")

mol, mf = setup_acetylene()
visualize_orbitals(mol, mf)

# --- CONFIGURATION: Choose your Mixing Here ---

# Scenario A: The "Vortex" (Degenerate Pi mix)
# Acetylene (STO-3G): MO 5 and 6 are the degenerate Pi_u (HOMO)
# We mix them as (1, i) to create rotation.
# Because E5 ~= E6, the time dependence exp(-i(E5-E6)t) is constant (static shape).
indices_A = [5, 6]
coeffs_A = [1.0, 1.0j] # 1 * px + i * py

# Scenario B: "Quantum Sloshing" (Sigma - Pi mixing)
# Mix HOMO (Pi, idx 5) with a Sigma orbital (idx 4)
# Energies are different. The density will beat at freq (E_pi - E_sigma).
indices_B = [4, 5] 
coeffs_B = [1.0, 1.0]

# SELECT SCENARIO
SCENARIO = 'A' # Change to 'B' to see sloshing

if SCENARIO == 'A':
    print("\nVisualizing Scenario A: Degenerate Pi-Mixing (Static Vortex)")
    print("Mixing MO 5 (px) + i * MO 6 (py)")
    chosen_indices = indices_A
    chosen_coeffs = coeffs_A
    # For Acetylene vortex, look at XY cross section (cutting the bond)
    slice_ax = 'xy'
    slice_z = 0.5 # Slightly above carbon to slice the lobe
    time_scale = 1.0
elif SCENARIO == 'B':
    print("\nVisualizing Scenario B: Non-Degenerate Mixing (Charge Sloshing)")
    print("Mixing MO 4 (Sigma) + MO 5 (Pi)")
    chosen_indices = indices_B
    chosen_coeffs = coeffs_B
    # For sloshing, look side-on (XZ)
    slice_ax = 'xz'
    slice_z = 0.0
    # Beat frequency is dE. Scale time to see oscillation.
    dE = abs(mf.mo_energy[chosen_indices[0]] - mf.mo_energy[chosen_indices[1]])
    time_scale = 2 * np.pi / dE / 50.0 # Resolution
    print(f"Energy Gap: {dE:.4f} Ha. Time scaling applied.")

# --- Setup Animation ---
L_box = 3.0
u, v, coords = get_slice_grid(L=L_box, N=40, axis=slice_ax, z_val=slice_z)

fig, ax = plt.subplots(figsize=(7, 6))

# Initial Plot
den, jx, jy, jz, psi = eval_wavefunction_dynamics(mol, mf, chosen_indices, chosen_coeffs, coords, 0.0)

# Reshape for contour
shape = u.shape
den_grid = den.reshape(shape)

if slice_ax == 'xy':
    J_u, J_v = jx.reshape(shape), jy.reshape(shape)
else:
    J_u, J_v = jx.reshape(shape), jz.reshape(shape)

# Plot elements
contour = ax.contourf(u, v, den_grid, levels=30, cmap='inferno')
# Quiver (Current Flux)
# We mask small currents to avoid visual clutter
mask = den_grid > (np.max(den_grid) * 0.05)
quiver = ax.quiver(u[mask], v[mask], J_u[mask], J_v[mask], color='cyan', scale=20, width=0.005)
title = ax.set_title(f"Time t=0.00")

def update(frame):
    t = frame * time_scale
    
    # Recalculate
    den, jx, jy, jz, psi = eval_wavefunction_dynamics(mol, mf, chosen_indices, chosen_coeffs, coords, t)
    den_grid = den.reshape(shape)
    
    if slice_ax == 'xy':
        J_u, J_v = jx.reshape(shape), jy.reshape(shape)
    else:
        J_u, J_v = jx.reshape(shape), jz.reshape(shape)
    
    # Update Density Plot
    ax.clear()
    ax.contourf(u, v, den_grid, levels=30, cmap='inferno')
    
    # Update Flux Vectors
    mask = den_grid > (np.max(den_grid) * 0.05)
    if np.any(mask):
        ax.quiver(u[mask], v[mask], J_u[mask], J_v[mask], color='cyan', scale=None, width=0.005)
    
    # Highlight Zeros (Vortex Cores)
    # Plot contour near zero
    ax.contour(u, v, den_grid, levels=[0.001], colors='white', linewidths=1.5, linestyles='--')
    
    ax.set_title(f"Density & Current Flux | Time t={t:.2f} au")
    ax.set_xlim(-L_box, L_box)
    ax.set_ylim(-L_box, L_box)

ani = FuncAnimation(fig, update, frames=np.arange(0, 100), interval=50)
plt.show()