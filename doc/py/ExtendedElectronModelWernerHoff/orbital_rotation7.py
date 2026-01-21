import numpy as np
import matplotlib.pyplot as plt
from pyscf import gto, scf

'''

https://drive.google.com/file/d/1CEoLgk0kOBtpHT93aIAJxkUDa5wW-4eB/view?usp=sharing, 
https://drive.google.com/file/d/1CRuiU2tAXgQJEYgyhmKVGw-PW1RE4lsX/view?usp=sharing, 
https://drive.google.com/file/d/1FL6k14iDRVhL3Fz4DgilARsBzQbL6_yv/view?usp=sharing, 
https://drive.google.com/file/d/1TJxL5GzhRrMHoiQvyop72fMdssevNz-3/view?usp=sharing, 
https://aistudio.google.com/app/prompts?state=%7B%22ids%22:%5B%221TnNS3SOZe2bek_mcFFNOJh03lBf5owmO%22%5D,%22action%22:%22open%22,%22userId%22:%22100958146796876347936%22,%22resourceKeys%22:%7B%7D%7D&usp=sharing, 
https://drive.google.com/file/d/1tbCYUF7Gat0Te75nA4Zjb1lWvVIrCL7z/view?usp=sharing
'''


def calculate_topological_flow():
    # 1. Define Ethylene (Symmetric D2h)
    cc_half = 0.67
    ch_len = 1.09
    hx = cc_half + ch_len * 0.5
    hy = ch_len * 0.866
    
    mol = gto.M(
        atom=f"""
        C  {-cc_half}  0.0  0.0
        C  {cc_half}   0.0  0.0
        H  {-hx}   {hy}  0.0
        H  {-hx}  {-hy}  0.0
        H  {hx}    {hy}  0.0
        H  {hx}   {-hy}  0.0
        """,
        basis='6-31g', verbose=0, symmetry=True) # Symmetry is crucial!
    
    mf = scf.RHF(mol)
    mf.kernel()
    
    # 2. Setup High-Res Grid
    res = 200
    lim = 3.5
    x = np.linspace(-lim, lim, res)
    y = np.linspace(-lim, lim, res)
    X, Y = np.meshgrid(x, y)
    
    # Slice slightly above Z=0 to define orientation
    coords = np.zeros((res**2, 3))
    coords[:, 0] = X.ravel()
    coords[:, 1] = Y.ravel()
    coords[:, 2] = 0.15
    
    # 3. Compute Pauli Kinetic Energy Density (Scalar)
    # This is our robust, invariant "Height Map"
    print("Calculating Pauli Energy Landscape...")
    ao_vals = mol.eval_gto("GTOval_sph_deriv1", coords)
    mo_coeffs = mf.mo_coeff
    n_occ = mol.nelectron // 2
    
    rho = np.zeros(res**2)
    tau_total = np.zeros(res**2)
    grad_rho = np.zeros((3, res**2))
    
    for i in range(n_occ):
        c = mo_coeffs[:, i]
        psi = np.dot(ao_vals[0], c)
        dx  = np.dot(ao_vals[1], c)
        dy  = np.dot(ao_vals[2], c)
        dz  = np.dot(ao_vals[3], c)
        
        rho_i = 2.0 * psi**2
        rho += rho_i
        
        # Gradient of Density
        grad_rho[0] += 4.0 * psi * dx
        grad_rho[1] += 4.0 * psi * dy
        grad_rho[2] += 4.0 * psi * dz
        
        # Kinetic Energy Density
        tau_total += (dx**2 + dy**2 + dz**2)

    # Von Weizsacker Term
    grad_rho_sq = np.sum(grad_rho**2, axis=0)
    tau_vw = (1.0/8.0) * grad_rho_sq / (rho + 1e-12)
    
    # Pauli Energy (Scalar)
    tau_p = np.maximum(tau_total - tau_vw, 0).reshape((res, res))
    
    # 4. Construct the Vector Field (The Topological Trick)
    # The flow runs along the isolines of the Pauli Energy.
    # Flow = Curl(Tau_P) => Rotate Gradient by 90 degrees.
    
    # Calculate Gradient of Tau_P numerically
    # (grad_y, grad_x) because numpy arrays are [row, col]
    gy, gx = np.gradient(tau_p)
    
    # Rotate 90 degrees: (x, y) -> (-y, x)
    # This creates a solenoid field circulating around the peaks
    Vx = -gy
    Vy = gx
    
    # 5. Extract the "Sigma-Pi" Twist (Specific Orbital Interaction)
    # To confirm the direction, we calculate the HOMO x HOMO-1 cross product
    # HOMO is B1u (Pi), HOMO-1 is Ag (Sigma)
    # (Indices depend on basis set size, usually last two occupied)
    c_pi = mo_coeffs[:, n_occ-1]
    c_sig = mo_coeffs[:, n_occ-2]
    
    psi_pi = np.dot(ao_vals[0], c_pi)
    grad_pi = np.array([np.dot(ao_vals[1], c_pi), np.dot(ao_vals[2], c_pi), np.dot(ao_vals[3], c_pi)])
    
    psi_sig = np.dot(ao_vals[0], c_sig)
    grad_sig = np.array([np.dot(ao_vals[1], c_sig), np.dot(ao_vals[2], c_sig), np.dot(ao_vals[3], c_sig)])
    
    # Cross Product (Sigma x Pi)
    # In XY plane, Pi grad is mostly Z, Sigma grad is XY.
    # Result is in XY plane.
    v_twist_x = grad_sig[1]*grad_pi[2] - grad_sig[2]*grad_pi[1]
    v_twist_y = grad_sig[2]*grad_pi[0] - grad_sig[0]*grad_pi[2]
    
    return X, Y, tau_p, Vx, Vy, v_twist_x.reshape((res, res)), v_twist_y.reshape((res, res)), mol

# --- PLOTTING ---
X, Y, Tau_P, Vx, Vy, Tx, Ty, mol = calculate_topological_flow()

fig, ax = plt.subplots(1, 2, figsize=(16, 8), facecolor='#111111')

# Plot 1: The "Topological Flow" (Rotated Gradient of Pauli Energy)
# This shows the global circulation of the electron fluid.
ax[0].set_facecolor('black')
ax[0].set_title("1. Global Pauli Flow Field\n(Calculated from $\\nabla \\tau_P$)", color='white', fontsize=14)

# Scalar Background (Filaments)
im0 = ax[0].imshow(np.log1p(Tau_P), extent=[-3.5,3.5,-3.5,3.5], origin='lower', cmap='magma')

# Streamlines (brighter/thicker)
speed = np.sqrt(Vx**2 + Vy**2) + 1e-9
lw = 0.8 + 2.5 * speed / speed.max()
ax[0].streamplot(X, Y, Vx, Vy, color='cyan', density=2.0, linewidth=lw, arrowsize=1.6)

# Quiver overlay for direction cues
step = max(1, Vx.shape[0]//22)
ax[0].quiver(X[::step, ::step], Y[::step, ::step],
             Vx[::step, ::step]/speed[::step, ::step],
             Vy[::step, ::step]/speed[::step, ::step],
             color='white', scale=18, width=0.0045, alpha=0.8)


# Plot 2: The "Sigma-Pi Twist" (Specific Interaction)
# This shows the specific mechanical rotation of the double bond.
ax[1].set_facecolor('black')
ax[1].set_title("2. The C=C Bond Twist\n(HOMO $\\times$ HOMO-1 Interaction)", color='white', fontsize=14)

im1 = ax[1].imshow(np.log1p(Tau_P), extent=[-3.5,3.5,-3.5,3.5], origin='lower', cmap='gray', alpha=0.5)

# Twist Vectors
tspeed = np.sqrt(Tx**2 + Ty**2) + 1e-9
tlw = 1.0 + 2.5 * tspeed / tspeed.max()
# Use bright color for flow lines
ax[1].streamplot(X, Y, Tx, Ty, color='yellow', density=1.8, linewidth=tlw, arrowsize=1.6)
# Quiver overlay
step2 = max(1, Tx.shape[0]//22)
ax[1].quiver(X[::step2, ::step2], Y[::step2, ::step2],
             Tx[::step2, ::step2]/tspeed[::step2, ::step2],
             Ty[::step2, ::step2]/tspeed[::step2, ::step2],
             color='white', scale=18, width=0.0045, alpha=0.8)

# Atoms Overlay
atom_coords = mol.atom_coords()
for a in ax:
    a.scatter(atom_coords[:,0], atom_coords[:,1], color='white', s=100, zorder=10)
    a.plot([atom_coords[0,0], atom_coords[1,0]], [0,0], 'w-', alpha=0.5)

plt.tight_layout()
plt.show()