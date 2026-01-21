import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.colors import hsv_to_rgb
from pyscf import gto, scf

# --- 1. SYSTEM SETUP (Ethylene) ---
def build_molecule():
    # C2H4 Geometry (Au)
    # C=C aligned along Z, Molecule in XZ plane
    mol = gto.M(
        atom="""
        C  0.0  0.0   1.26
        C  0.0  0.0  -1.26
        H  1.75 0.0   2.3
        H -1.75 0.0   2.3
        H  1.75 0.0  -2.3
        H -1.75 0.0  -2.3
        """,
        basis='sto-3g',
        verbose=0
    )
    return mol

# --- 2. CALCULATE ORBITALS ---
def run_scf(mol):
    print("\n[.] Running Hartree-Fock calculation (STO-3G)...")
    mf = scf.RHF(mol)
    mf.kernel()
    return mf

def print_occupied_table(mf, mol):
    n_occ = mol.nelectron // 2
    print("\nOccupied Molecular Orbitals:")
    print(f"{'Idx':<4} | {'Energy (Ha)':<12}")
    print("-"*22)
    for i in range(n_occ):
        print(f"{i:<4} | {mf.mo_energy[i]:<12.4f}")

def parse_mix_string(mix_str, homo_index):
    """Parse CLI mix string like '6:1.0,7:1.0j' into [{'idx':6,'c':1+0j}, ...]"""
    components = []
    parts = [p for p in mix_str.split(',') if p.strip()]
    for p in parts:
        idx_str, coeff_str = p.split(':')
        idx = int(idx_str.strip())
        coeff = complex(coeff_str.strip())
        components.append({'idx': idx, 'c': coeff})
    if not components:
        raise ValueError("Empty mix string")
    return components

# --- 4. GRID GENERATION ---
def make_grid(plane='xz', res=100, L=6.0):
    # 1D arrays
    u = np.linspace(-L, L, res)
    v = np.linspace(-L, L, res)
    U, V = np.meshgrid(u, v)
    
    coords = np.zeros((res*res, 3))
    
    if plane == 'xz':
        # Plane Y=0 (Molecular plane)
        coords[:, 0] = U.ravel() # X
        coords[:, 1] = 0.0       # Y
        coords[:, 2] = V.ravel() # Z
        xlabel, ylabel = 'X', 'Z'
    else:
        # Plane Z=0 (Cutting through C-C bond center)
        # But wait, atoms are at Z +/- 1.26. 
        # Let's cut at Z=0 to see the cross section of the bond.
        coords[:, 0] = U.ravel() # X
        coords[:, 1] = V.ravel() # Y
        coords[:, 2] = 0.0       # Z (Center of bond)
        xlabel, ylabel = 'X', 'Y'
        
    return U, V, coords, xlabel, ylabel

# --- 5. ANIMATION ENGINE ---
def visualize(mol, mf, components, mode, plane, frames=200, res=120):
    
    # 1. Setup Grid
    U, V, coords, xl, yl = make_grid(plane, res=res)
    
    # 2. Precompute Spatial Orbitals on Grid
    print(f"\n[.] Evaluating Orbitals on {plane.upper()} grid...")
    # shape: (4, N_points, N_orbitals). Index 0 is value.
    ao_values = mol.eval_gto("GTOval", coords) 
    mo_coeffs = mf.mo_coeff
    
    # We only store the specific spatial MOs requested to save memory
    spatial_mos = {}
    energies = {}
    
    for item in components:
        idx = item['idx']
        # Transform AO to MO: sum( c_i * ao_i )
        mo_spatial = np.dot(ao_values, mo_coeffs[:, idx])
        spatial_mos[idx] = mo_spatial
        energies[idx] = mf.mo_energy[idx]
        
    # 3. Time Scaling
    # Find beat frequency to set animation speed
    e_vals = list(energies.values())
    if len(e_vals) > 1:
        dE = max(e_vals) - min(e_vals)
        if dE < 1e-6: dE = 1.0 # Degenerate
    else:
        dE = 1.0
    
    # If dE is 0.5 Hartree, period is 2pi/0.5 ~ 12 time units.
    # We want ~100 frames to cover a period.
    dt = (2 * np.pi / dE) / 50.0
    print(f"[.] Time step: {dt:.4f} au (based on Energy gap {dE:.4f})")

    # 4. Plot Setup
    fig, ax = plt.subplots(figsize=(7, 6), facecolor='black')
    
    # Initial Calculation for scale
    psi_0 = np.zeros(res*res, dtype=complex)
    for item in components:
        psi_0 += item['c'] * spatial_mos[item['idx']]
        
    rho_max = np.max(np.abs(psi_0)**2)
    
    # Image container
    img = ax.imshow(np.zeros((res, res)), origin='lower', extent=[-6,6,-6,6], animated=True)
    
    # Decoration
    ax.set_title(f"Ethylene | Mode: {mode.upper()}", color='white')
    ax.set_xlabel(xl, color='white')
    ax.set_ylabel(yl, color='white')
    ax.tick_params(colors='white')
    
    # Draw Atoms markers (Approximate projection)
    if plane == 'xz':
        ax.scatter([0,0], [1.26, -1.26], c='cyan', marker='o', label='C')
        ax.scatter([1.75, -1.75], [2.3, 2.3], c='white', s=10, label='H')
        ax.scatter([1.75, -1.75], [-2.3, -2.3], c='white', s=10)
    elif plane == 'xy':
        # Looking down Z-axis. C atoms overlap at origin (0,0) in projection?
        # No, we are slicing Z=0. Atoms are above/below.
        # So we just mark the center.
        ax.scatter([0], [0], c='cyan', marker='x', label='Bond Center')

    txt = ax.text(0.05, 0.95, "", transform=ax.transAxes, color='white')

    def update(frame):
        t = frame * dt
        
        # Superposition
        psi = np.zeros(res*res, dtype=complex)
        for item in components:
            idx = item['idx']
            amp = item['c']
            E = energies[idx]
            # Time Evolution
            psi += amp * spatial_mos[idx] * np.exp(-1j * E * t)
            
        psi_grid = psi.reshape((res, res))
        
        if mode == 'density':
            # Simple Density Plot
            data = np.abs(psi_grid)**2
            img.set_data(data)
            img.set_cmap('inferno')
            img.set_clim(0, rho_max)
            
        elif mode == 'phase':
            # Phase-Color Mapping (The Rainbow)
            # 1. Amplitude (Brightness)
            amp = np.abs(psi_grid)
            # Normalize amp for display
            val = np.clip(amp / np.sqrt(rho_max), 0, 1)
            
            # 2. Phase (Color)
            phase = np.angle(psi_grid) # -pi to pi
            hue = (phase + np.pi) / (2 * np.pi) # 0 to 1
            
            # 3. Saturation (Always 1)
            sat = np.ones_like(hue)
            
            # Stack HSV
            hsv = np.dstack((hue, sat, val))
            rgb = hsv_to_rgb(hsv)
            img.set_data(rgb)
            
        txt.set_text(f"t={t:.2f} au")
        return img, txt

    ani = FuncAnimation(fig, update, frames=frames, interval=50, blit=True)
    plt.show()

def parse_args(mf, mol):
    n_occ = mol.nelectron // 2
    homo_idx = n_occ - 1
    default_mix = f"{homo_idx}:1.0,{homo_idx-1}:1.0j"  # HOMO + HOMO-1 gives time variation
    parser = argparse.ArgumentParser(description="Visualize ethylene orbital rotations (no prompts)")
    parser.add_argument("--mix", type=str, default=default_mix,
                        help="Comma list of idx:coeff (e.g. '6:1.0,7:1.0j')")
    parser.add_argument("--mode", choices=["density", "phase"], default="density", help="Render density or phase hue")
    parser.add_argument("--plane", choices=["xz", "xy"], default="xz", help="Slice plane: xz (molecular plane) or xy (cross-section)")
    parser.add_argument("--frames", type=int, default=200, help="Animation frames")
    parser.add_argument("--res", type=int, default=120, help="Grid resolution")
    args = parser.parse_args()
    try:
        components = parse_mix_string(args.mix, homo_idx)
    except Exception as e:
        print(f"[WARN] Failed to parse mix '{args.mix}': {e}. Using HOMO.")
        components = [{'idx': homo_idx, 'c': 1.0}]
    return args, components


# --- MAIN EXECUTION ---
if __name__ == "__main__":
    mol = build_molecule()
    mf = run_scf(mol)
    print_occupied_table(mf, mol)
    args, comps = parse_args(mf, mol)
    visualize(mol, mf, comps, args.mode, args.plane, frames=args.frames, res=args.res)