import numpy as np
import matplotlib.pyplot as plt
from pyscf import gto, scf
import itertools

def calculate_ethylene_vortex_structure():
    # 1. Define Ethylene Geometry (In XY Plane)
    # C-C bond length approx 1.34 A
    # C-H bond length approx 1.09 A
    # Angles approx 120 degrees
    
    cc_half = 1.34 / 2
    ch_len = 1.09
    
    # Trigonometry for H positions (120 deg from X axis)
    # cos(60) = 0.5, sin(60) = 0.866
    h_x_offset = ch_len * 0.5
    h_y_offset = ch_len * 0.866
    
    mol_str = f"""
    C  {-cc_half}  0.0  0.0
    C  {cc_half}   0.0  0.0
    H  {-cc_half - h_x_offset}   {h_y_offset}  0.0
    H  {-cc_half - h_x_offset}  {-h_y_offset}  0.0
    H  {cc_half + h_x_offset}    {h_y_offset}  0.0
    H  {cc_half + h_x_offset}   {-h_y_offset}  0.0
    """
    
    # 2. Run DFT (RHF)
    print("Computing Molecular Orbitals...")
    mol = gto.M(atom=mol_str, basis='6-31g', verbose=0)
    mf = scf.RHF(mol)
    mf.kernel()
    
    # 3. Setup Grid (XY Plane, Z=0.001)
    # We use Z=0.001 slightly off-plane to catch the Pi-orbital gradients 
    # which are zero exactly at Z=0 but have strong vertical derivatives.
    res = 120
    lim = 3.5
    x = np.linspace(-lim, lim, res)
    y = np.linspace(-lim, lim, res)
    X, Y = np.meshgrid(x, y)
    
    coords = np.zeros((res**2, 3))
    coords[:, 0] = X.ravel()
    coords[:, 1] = Y.ravel()
    coords[:, 2] = 0.05 # Slightly above plane to capture Pi-interaction
    
    # 4. Evaluate Gradients of ALL Occupied Orbitals
    # ao_vals: [4, N_points, N_ao] -> (Val, dx, dy, dz)
    print("Evaluating Gradients on Grid...")
    ao_vals = mol.eval_gto("GTOval_sph_deriv1", coords)
    mo_coeffs = mf.mo_coeff
    n_occ = mol.nelectron // 2
    
    # Store gradients for all occupied orbitals
    # grads[i] = (3, N_points)
    grads = []
    
    for i in range(n_occ):
        c = mo_coeffs[:, i]
        # Calculate gradient field for orbital i
        dx = np.dot(ao_vals[1], c)
        dy = np.dot(ao_vals[2], c)
        dz = np.dot(ao_vals[3], c)
        grads.append(np.array([dx, dy, dz]))

    # 5. Compute "Total Pauli Vorticity"
    # Sum of Cross Products of all pairs
    V_net = np.zeros((3, res**2))
    Scalar_Stress = np.zeros(res**2)
    
    print(f"Summing interactions of {n_occ*(n_occ-1)//2} orbital pairs...")
    
    # Iterate all unique pairs
    for i, j in itertools.combinations(range(n_occ), 2):
        # Grad_i x Grad_j
        gi = grads[i]
        gj = grads[j]
        
        # Cross product components
        # cx = gi_y * gj_z - gi_z * gj_y
        cx = gi[1]*gj[2] - gi[2]*gj[1]
        cy = gi[2]*gj[0] - gi[0]*gj[2]
        cz = gi[0]*gj[1] - gi[1]*gj[0]
        
        # Accumulate Vector Field
        # Note: The sign is arbitrary (i,j vs j,i), but consistent summing 
        # reveals the topological structure.
        V_net[0] += cx
        V_net[1] += cy
        V_net[2] += cz
        
        # Accumulate Scalar "Stress" (Magnitude)
        # This shows "Where" the filaments are, regardless of direction
        mag = np.sqrt(cx**2 + cy**2 + cz**2)
        Scalar_Stress += mag

    # Reshape
    Vx = V_net[0].reshape((res, res))
    Vy = V_net[1].reshape((res, res))
    Vz = V_net[2].reshape((res, res))
    Stress = Scalar_Stress.reshape((res, res))
    
    return X, Y, Vx, Vy, Vz, Stress, mol

# --- VISUALIZATION ---
X, Y, Vx, Vy, Vz, Stress, mol = calculate_ethylene_vortex_structure()

fig, ax = plt.subplots(1, 2, figsize=(16, 7), facecolor='black')

# Plot 1: The "Filament Walls" (Scalar Stress)
# This shows where the gradients are fighting maximally.
ax[0].set_facecolor('black')
im1 = ax[0].imshow(Stress, extent=[-3.5,3.5,-3.5,3.5], origin='lower', cmap='magma', vmin=0)
ax[0].contour(X, Y, Stress, levels=15, colors='white', linewidths=0.5, alpha=0.6)
ax[0].set_title("1. Topological Stress Magnitude\nSum |$\\nabla \\psi_i \\times \\nabla \\psi_j$|", color='white')

# Overlay Atoms
atom_coords = mol.atom_coords()
# Convert Bohr to Angstrom for consistency if needed, but PySCF usually defaults to Bohr internally or Input.
# The input coordinates were Angstroms, PySCF converts to Bohr.
# Let's just plot based on grid limits.
ax[0].scatter(atom_coords[:,0], atom_coords[:,1], color='cyan', s=100, edgecolors='white')


# Plot 2: The "Current/Flow" (Vector Field)
# We plot the In-Plane components (Vx, Vy) as streamlines.
# We color them by the Vertical component (Vz) - Blue = CW, Red = CCW
ax[1].set_facecolor('black')
ax[1].set_title("2. Pauli Vorticity Flow\n(Streamlines = In-Plane Twist, Color = Vertical Twist)", color='white')

# Streamplot of the "Sigma-Pi" interaction (which generates in-plane vorticity)
speed = np.sqrt(Vx**2 + Vy**2) + 1e-9
lw = 2 * speed / speed.max()

strm = ax[1].streamplot(X, Y, Vx, Vy, color=Vz, cmap='coolwarm', density=1.5, linewidth=lw, arrowsize=1.5)

# Add sparse quiver for clearer direction cues
step = max(1, Vx.shape[0]//25)
mask = (speed > 0.2 * speed.max())
ax[1].quiver(X[::step, ::step], Y[::step, ::step],
             Vx[::step, ::step]/speed[::step, ::step],
             Vy[::step, ::step]/speed[::step, ::step],
             Vz[::step, ::step], cmap='coolwarm', scale=20, width=0.004, alpha=0.7)

# Overlay Atoms
ax[1].scatter(atom_coords[:,0], atom_coords[:,1], color='cyan', s=100, zorder=10)

# Add bond lines for reference
ax[1].plot([atom_coords[0,0], atom_coords[1,0]], [0,0], 'w-', alpha=0.3) # C-C

plt.tight_layout()
plt.show()