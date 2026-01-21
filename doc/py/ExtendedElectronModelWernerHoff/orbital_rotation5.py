import numpy as np
import matplotlib.pyplot as plt
from pyscf import gto, scf
import itertools

def calculate_symmetric_vortex_structure():
    # 1. Define Ethylene Geometry (In XY Plane)
    cc_half = 0.67  # 1.34 / 2
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
    
    # 2. Run DFT with SYMMETRY ENFORCED
    print("Computing Symmetric Molecular Orbitals...")
    # symmetry=True is the key fix!
    mol = gto.M(atom=mol_str, basis='6-31g', verbose=0, symmetry=True)
    mf = scf.RHF(mol)
    mf.kernel()
    
    # 3. Setup Grid (XY Plane)
    res = 140
    lim = 3.5
    x = np.linspace(-lim, lim, res)
    y = np.linspace(-lim, lim, res)
    X, Y = np.meshgrid(x, y)
    
    # Slice slightly above Z=0 to capture the Pi-orbital flow
    coords = np.zeros((res**2, 3))
    coords[:, 0] = X.ravel()
    coords[:, 1] = Y.ravel()
    coords[:, 2] = 0.1 
    
    # 4. Evaluate Gradients
    # GTOval_sph_deriv1 returns [Value, dX, dY, dZ]
    print("Evaluating Gradients...")
    ao_vals = mol.eval_gto("GTOval_sph_deriv1", coords)
    mo_coeffs = mf.mo_coeff
    
    # Only use Occupied Orbitals
    n_occ = mol.nelectron // 2
    grads = []
    psi_vals = []
    energies = mf.mo_energy[:n_occ]
    
    for i in range(n_occ):
        c = mo_coeffs[:, i]
        val = np.dot(ao_vals[0], c)
        dx  = np.dot(ao_vals[1], c)
        dy  = np.dot(ao_vals[2], c)
        dz  = np.dot(ao_vals[3], c)
        
        psi_vals.append(val)
        grads.append(np.array([dx, dy, dz]))

    # 5. Compute Pauli Vorticity
    V_net = np.zeros((3, res**2))
    Scalar_Stress = np.zeros(res**2)
    
    print("Summing interactions...")
    
    for i, j in itertools.combinations(range(n_occ), 2):
        # Weighting: Interactions are strongest between orbitals that are 
        # spatially overlapping. We weight by Density_i * Density_j.
        # This removes "ghost" vortices far from the molecule.
        weight = np.abs(psi_vals[i] * psi_vals[j])
        
        gi = grads[i]
        gj = grads[j]
        
        # Cross Product: v = gi x gj
        cx = gi[1]*gj[2] - gi[2]*gj[1]
        cy = gi[2]*gj[0] - gi[0]*gj[2]
        cz = gi[0]*gj[1] - gi[1]*gj[0]
        
        # Accumulate Vector Field
        # Note: Vector addition can still be destructive (canceling), 
        # but Scalar Stress is additive.
        V_net[0] += cx 
        V_net[1] += cy 
        V_net[2] += cz 
        
        # Accumulate Magnitude (Stress)
        mag = np.sqrt(cx**2 + cy**2 + cz**2)
        Scalar_Stress += mag # * weight # Uncomment weight to clean up edges

    # Reshape
    Vx = V_net[0].reshape((res, res))
    Vy = V_net[1].reshape((res, res))
    Vz = V_net[2].reshape((res, res))
    Stress = Scalar_Stress.reshape((res, res))
    
    return X, Y, Vx, Vy, Vz, Stress, mol

# --- PLOTTING ---
X, Y, Vx, Vy, Vz, Stress, mol = calculate_symmetric_vortex_structure()

fig, ax = plt.subplots(1, 2, figsize=(16, 8), facecolor='#111111')

# 1. Scalar Stress (The Filament Map)
ax[0].set_facecolor('black')
# Use a Log scale to see both strong core filaments and weak outer ones
im1 = ax[0].imshow(np.log1p(Stress), extent=[-3.5,3.5,-3.5,3.5], origin='lower', cmap='inferno')
ax[0].contour(X, Y, Stress, levels=12, colors='white', linewidths=0.5, alpha=0.4)
ax[0].set_title("Topological Stress (Filaments)", color='white', fontsize=14)

# 2. Vector Flow (The Hofer Spin)
ax[1].set_facecolor('black')
ax[1].set_title("Pauli Vorticity Flow (XY Plane)", color='white', fontsize=14)

# Streamplot colored by Z-vorticity (Red = CCW, Blue = CW)
speed = np.sqrt(Vx**2 + Vy**2) + 1e-9
lw = 2.0 * speed / speed.max()
ax[1].imshow(Stress, extent=[-3.5,3.5,-3.5,3.5], origin='lower', cmap='Greys', alpha=0.12)
strm = ax[1].streamplot(X, Y, Vx, Vy, color=Vz, cmap='coolwarm', 
                        density=1.6, linewidth=lw, arrowsize=1.5)
plt.colorbar(strm.lines, ax=ax[1], label='Vz (twist out of plane)')

# Add sparse quiver for clearer direction cues
step = max(1, Vx.shape[0]//24)
ax[1].quiver(X[::step, ::step], Y[::step, ::step],
             Vx[::step, ::step]/speed[::step, ::step],
             Vy[::step, ::step]/speed[::step, ::step],
             Vz[::step, ::step], cmap='coolwarm', scale=19, width=0.0042, alpha=0.75)

# Overlay Atoms
atom_coords = mol.atom_coords() # Returns Bohr (approx)
for i in range(2):
    ax[0].scatter(atom_coords[:,0], atom_coords[:,1], color='cyan', s=120, edgecolors='white', zorder=10)
    ax[1].scatter(atom_coords[:,0], atom_coords[:,1], color='cyan', s=120, edgecolors='white', zorder=10)

plt.tight_layout()
plt.show()