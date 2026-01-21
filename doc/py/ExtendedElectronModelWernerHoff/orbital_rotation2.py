import numpy as np
import matplotlib.pyplot as plt
from pyscf import gto, scf

def get_pauli_fields(mol_str, basis='sto-3g', slice_z=0.0, grid_res=100, box_L=4.0):
    """
    Calculates Density, Pauli Energy, and 'Vorticity Proxy' from DFT orbitals.
    """
    # 1. Build Molecule & Solve DFT
    mol = gto.M(atom=mol_str, basis=basis, verbose=0)
    mf = scf.RHF(mol)
    mf.kernel()
    
    # 2. Setup Grid
    x = np.linspace(-box_L, box_L, grid_res)
    y = np.linspace(-box_L, box_L, grid_res)
    X, Y = np.meshgrid(x, y)
    
    # Grid points in 3D (Z is fixed at slice_z)
    coords = np.zeros((grid_res**2, 3))
    coords[:, 0] = X.ravel()
    coords[:, 1] = Y.ravel()
    coords[:, 2] = slice_z
    
    # 3. Evaluate Orbitals and Gradients
    # ao_vals shape: [4, N_points, N_ao]. Index 0=Value, 1=dx, 2=dy, 3=dz
    ao_vals = mol.eval_gto("GTOval_sph_deriv1", coords)
    mo_coeffs = mf.mo_coeff
    occ = mf.mo_occ
    
    # Identify occupied orbitals
    occ_idxs = np.where(occ > 0)[0]
    
    # Initialize Accumulators
    rho = np.zeros(grid_res**2)
    grad_rho = np.zeros((3, grid_res**2))
    tau_total = np.zeros(grid_res**2) # Total Kinetic Energy Density
    
    # Store individual orbital gradients for Vorticity calc
    orb_grads = [] 
    
    for idx in occ_idxs:
        c = mo_coeffs[:, idx]
        
        # Orbital Psi and Gradient
        psi = np.dot(ao_vals[0], c)
        dpsi_x = np.dot(ao_vals[1], c)
        dpsi_y = np.dot(ao_vals[2], c)
        dpsi_z = np.dot(ao_vals[3], c)
        
        # Density contribution (2 electrons per orbital for RHF)
        dens_i = 2.0 * psi**2
        rho += dens_i
        
        # Gradient Density contribution: Grad(psi^2) = 2*psi*Grad(psi)
        grad_rho[0] += 2.0 * 2.0 * psi * dpsi_x
        grad_rho[1] += 2.0 * 2.0 * psi * dpsi_y
        grad_rho[2] += 2.0 * 2.0 * psi * dpsi_z
        
        # Kinetic Energy contribution: |Grad(psi)|^2
        # T = 1/2 * Sum |Grad Psi|^2
        grad_sq = dpsi_x**2 + dpsi_y**2 + dpsi_z**2
        tau_total += 2.0 * 0.5 * grad_sq
        
        # Store for vorticity (weighted by occupancy/density)
        orb_grads.append({
            'psi': psi,
            'grad': np.array([dpsi_x, dpsi_y, dpsi_z])
        })

    # 4. Calculate Von Weizsacker (Bosonic) Kinetic Energy
    # T_vw = 1/8 * |Grad Rho|^2 / Rho
    grad_rho_sq = np.sum(grad_rho**2, axis=0)
    tau_vw = (1.0/8.0) * grad_rho_sq / (rho + 1e-12)
    
    # 5. Calculate Pauli Kinetic Energy (The "Hidden" Energy)
    tau_pauli = tau_total - tau_vw
    tau_pauli = np.maximum(tau_pauli, 0) # Numerical cleanup
    
    # 6. Calculate "Hofer Vorticity Proxy"
    # Vector Field V = Sum_{i<j} (Grad_i x Grad_j)
    # This measures the non-colinearity of the orbital gradients.
    vorticity = np.zeros((3, grid_res**2))
    
    # We only take the Cross Product of the HOMO and HOMO-1 
    # (The most active valence electrons determining chemistry)
    # Using all pairs is computationally heavy and noisy.
    
    # Get last two occupied orbitals
    orb_A = orb_grads[-1] # HOMO (Pi bond)
    orb_B = orb_grads[-2] # HOMO-1 (Sigma skeleton)
    
    # Cross Product of Gradients
    # V = Grad(Psi_A) x Grad(Psi_B)
    gxA, gyA, gzA = orb_A['grad']
    gxB, gyB, gzB = orb_B['grad']
    
    vorticity[0] = gyA*gzB - gzA*gyB # X component
    vorticity[1] = gzA*gxB - gxA*gzB # Y component
    vorticity[2] = gxA*gyB - gyA*gxB # Z component
    
    # Reshape for plotting
    shape = (grid_res, grid_res)
    return X, Y, rho.reshape(shape), tau_pauli.reshape(shape), vorticity, shape

# --- MAIN EXECUTION ---

# Ethylene Geometry (C=C along Z axis)
# H atoms in XZ plane.
mol_str = """
C 0.0 0.0  0.66
C 0.0 0.0 -0.66
H 1.7 0.0  1.2
H -1.7 0.0  1.2
H 1.7 0.0 -1.2
H -1.7 0.0 -1.2
"""

# Slice at Z=0 (Cutting through the middle of the C=C bond)
# This is where the Pi-bond node is (the XY plane).
X, Y, rho, tau_p, vort_flat, shape = get_pauli_fields(mol_str, slice_z=0.0)

# Extract Vector Components for XY plane
Vx = vort_flat[0].reshape(shape)
Vy = vort_flat[1].reshape(shape)
Vz = vort_flat[2].reshape(shape) # This is the "Swirl" in the plane

# --- VISUALIZATION ---

fig, ax = plt.subplots(1, 3, figsize=(18, 6), facecolor='white')

# 1. Total Density (Standard)
ax[0].set_title("1. Total Electron Density $\\rho$")
im0 = ax[0].imshow(rho, extent=[-4,4,-4,4], origin='lower', cmap='Blues')
ax[0].contour(X, Y, rho, colors='k', linewidths=0.5)
plt.colorbar(im0, ax=ax[0])

# 2. Pauli Kinetic Energy (The Scalar "Pressure")
# This shows WHERE the "Hidden Rotation" energy is stored.
# Notice it is ZERO at the center (where density is high) 
# and HIGH at the nodal surfaces (where orbitals fight).
ax[1].set_title("2. Pauli Kinetic Energy $T_{Pauli}$")
im1 = ax[1].imshow(tau_p, extent=[-4,4,-4,4], origin='lower', cmap='inferno')
plt.colorbar(im1, ax=ax[1])

# 3. The Hofer Vorticity Field (Vector)
# This shows the AXIS of the rotation.
# We plot the Z-component (Spinning in the plane) as color, 
# and XY components as arrows.
ax[2].set_title("3. Hofer Vorticity Proxy\n$\\nabla \\psi_{HOMO} \\times \\nabla \\psi_{HOMO-1}$")
# Plot magnitude/direction of rotation
mag = np.sqrt(Vx**2 + Vy**2)
strm = ax[2].streamplot(X, Y, Vx, Vy, color=Vz, cmap='coolwarm', linewidth=1, density=1.5)
plt.colorbar(strm.lines, ax=ax[2], label='Z-Vorticity (Spin out of page)')

# Mark the center (Bond Axis)
ax[2].scatter([0], [0], color='black', marker='x')

plt.tight_layout()
plt.show()