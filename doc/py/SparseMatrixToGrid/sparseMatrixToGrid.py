import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def load_xyz(fname, rcut_default):
    with open(fname, "r") as f:
        lines = f.read().strip().splitlines()
    n = int(lines[0])
    atoms = []
    for line in lines[2:2+n]:
        parts = line.split()
        if len(parts) < 4: continue
        sym, x, y, z = parts[0], float(parts[1]), float(parts[2]), float(parts[3])
        atoms.append({'pos': np.array([x, y]), 'Rcut': rcut_default, 'type': sym})
    return atoms

def get_benzene_coords():
    # Standard Benzene: C-C dist ~1.4A, C-H dist ~1.1A
    cc = 1.40
    ch = 1.08
    angles = np.linspace(0, 2 * np.pi, 6, endpoint=False)
    
    # Carbons
    cx = cc * np.cos(angles)
    cy = cc * np.sin(angles)
    # Hydrogens
    hx = (cc + ch) * np.cos(angles)
    hy = (cc + ch) * np.sin(angles)
    
    atoms = []
    for i in range(6):
        atoms.append({'pos': np.array([cx[i], cy[i]]), 'Rcut': 3.5, 'type': 'C'})
    for i in range(6):
        atoms.append({'pos': np.array([hx[i], hy[i]]), 'Rcut': 2.5, 'type': 'H'})
    return atoms

def check_overlap_sphere_aabb(center, radius, box_min, box_max):
    """ Fast AABB-Sphere collision: Find closest point in box to sphere center """
    closest_p = np.clip(center, box_min, box_max)
    distance_sq = np.sum((center - closest_p)**2)
    return distance_sq < (radius**2)

def simulate_projection(grid_size=None, block_res=0.8, macro_res=None, margin=0.1, seed=None, verbose=True, xyz=None, rcut_default=3.0, bbox_margin=1.0):
    if seed is not None:
        np.random.seed(seed)
    atoms = load_xyz(xyz, rcut_default) if xyz is not None else get_benzene_coords()
    n_atoms = len(atoms)
    max_rcut = max(a["Rcut"] for a in atoms)
    if block_res is None:
        block_res = 0.8  # fine voxel default
        if verbose: print(f"[DEBUG] block_res auto-set to default fine voxel {block_res}")
    auto_macro_res = 2.0 * max_rcut + margin  # big boxes sized to sphere diameter + margin
    if macro_res is None:
        macro_res = auto_macro_res
        if verbose: print(f"[DEBUG] macro_res auto-set to {macro_res} (= 2*max_rcut + margin)")
    elif macro_res < auto_macro_res:
        if verbose: print(f"[DEBUG] macro_res bumped from {macro_res} to {auto_macro_res} to fit sphere diameter+margin")
        macro_res = auto_macro_res

    # bounding box based on atoms (planar, xy only)
    apos = np.array([a['pos'] for a in atoms])
    radii = np.array([a['Rcut'] for a in atoms])
    bb_min = (apos - radii[:, None]).min(axis=0) - bbox_margin
    bb_max = (apos + radii[:, None]).max(axis=0) + bbox_margin
    size_vec = bb_max - bb_min
    if grid_size is None:
        grid_size = float(size_vec.max())
        if verbose: print(f"[DEBUG] grid_size auto-set to {grid_size} from bbox")
    if verbose:
        print(f"[DEBUG] grid_size={grid_size} block_res={block_res} macro_res={macro_res} max_rcut={max_rcut} margin={margin} bbox_min={bb_min} bbox_max={bb_max}")
    
    # Setup Grid of Super-Voxels anchored at bbox_min
    ticks_x = np.arange(bb_min[0], bb_max[0], block_res)
    ticks_y = np.arange(bb_min[1], bb_max[1], block_res)
    grid_nx = len(ticks_x)
    grid_ny = len(ticks_y)
    
    # Data to store
    workload_map = np.zeros((grid_ny, grid_nx))   # pairs per fine block
    atom_count_map = np.zeros((grid_ny, grid_nx)) # atoms per fine block

    # Pre-calculate atom-atom adjacency (Density Matrix Sparsity)
    # In DFTB, Pi,j != 0 only if atoms are close
    adj_matrix = np.zeros((n_atoms, n_atoms), dtype=bool)
    for i in range(n_atoms):
        for j in range(i, n_atoms):
            dist = np.linalg.norm(atoms[i]['pos'] - atoms[j]['pos'])
            # Only consider pairs within combined cutoff
            if dist < (atoms[i]['Rcut'] + atoms[j]['Rcut']):
                adj_matrix[i, j] = adj_matrix[j, i] = True

    # Coarse macro mapping: atom -> macro cells it touches (big boxes)
    macro_grid = {}
    for idx, a in enumerate(atoms):
        x_range = [int((a['pos'][0]-a['Rcut'] - bb_min[0])/macro_res),
                   int((a['pos'][0]+a['Rcut'] - bb_min[0])/macro_res)]
        y_range = [int((a['pos'][1]-a['Rcut'] - bb_min[1])/macro_res),
                   int((a['pos'][1]+a['Rcut'] - bb_min[1])/macro_res)]
        for ix in range(x_range[0], x_range[1] + 1):
            for iy in range(y_range[0], y_range[1] + 1):
                macro_grid.setdefault((ix, iy), []).append(idx)

    # Process Blocks
    for ix, x in enumerate(ticks_x):
        for iy, y in enumerate(ticks_y):
            b_min = np.array([x, y])
            b_max = b_min + block_res
            macro_ix = int((x - bb_min[0]) / macro_res)
            macro_iy = int((y - bb_min[1]) / macro_res)
            macro_atoms = macro_grid.get((macro_ix, macro_iy), [])
            if verbose:
                print(f"[BLOCK] macro=({macro_ix},{macro_iy}) fine=({ix},{iy}) pos=({x:.2f},{y:.2f}) candidates={len(macro_atoms)}")
            if not macro_atoms:
                continue
            
            # 1. Find which atoms overlap this block
            overlapping_atoms = []
            for i in macro_atoms:
                if check_overlap_sphere_aabb(atoms[i]['pos'], atoms[i]['Rcut'], b_min, b_max):
                    overlapping_atoms.append(i)
            
            na = len(overlapping_atoms)
            atom_count_map[iy, ix] = na
            
            # 2. Count active density elements P_ij for this block
            ne = 0
            for idx_a in range(na):
                for idx_b in range(idx_a, na):
                    atom_i = overlapping_atoms[idx_a]
                    atom_j = overlapping_atoms[idx_b]
                    if adj_matrix[atom_i, atom_j]:
                        ne += 1
            
            workload_map[iy, ix] = ne

    return (ticks_x, ticks_y), atoms, atom_count_map, workload_map, block_res, macro_res, (bb_min, bb_max)
def parse_args():
    p = argparse.ArgumentParser(description="Sparse matrix to grid projection with coarse macro culling.")
    p.add_argument("--grid-size", type=float, default=None, dest="grid_size", help="Simulation box size; defaults to bbox span")
    p.add_argument("--block-res", type=float, default=0.8, dest="block_res", help="Fine voxel size")
    p.add_argument("--macro-res", type=float, default=None, dest="macro_res", help="Macro (super-box) size; default auto = 2 * maxRcut + margin")
    p.add_argument("--bbox-margin", type=float, default=1.0, dest="bbox_margin", help="Padding around molecular bbox")
    p.add_argument("--xyz", type=str, default="/home/prokophapala/git/FireCore/tests/tUFF/data/xyz/pentacene.xyz", help="Path to XYZ file; if absent, benzene is used")
    p.add_argument("--rcut-default", type=float, default=3.0, dest="rcut_default", help="Default Rcut assigned to atoms when loading XYZ")
    p.add_argument("--margin", type=float, default=0.1, help="Extra padding added to sphere diameter for block sizing")
    p.add_argument("--seed", type=int, default=None, help="Random seed for reproducibility")
    p.add_argument("--no-plot", action="store_true", help="Disable plotting")
    p.add_argument("--quiet", action="store_true", help="Reduce debug prints")
    return p.parse_args()

def main():
    args = parse_args()
    (ticks_x, ticks_y), atoms, atom_counts, workload, b_res, m_res, (bb_min, bb_max) = simulate_projection(
        grid_size=args.grid_size,
        block_res=args.block_res,
        macro_res=args.macro_res,
        bbox_margin=args.bbox_margin,
        xyz=args.xyz,
        rcut_default=args.rcut_default,
        margin=args.margin,
        seed=args.seed,
        verbose=not args.quiet,
    )
    if args.no_plot:
        return

    fig, ax = plt.subplots(1, 2, figsize=(14, 6))

    # Plot 1: Geometry and Atom Counts per fine block
    ax[0].set_title("Atoms & Fine Grid (counts in red)")
    for a in atoms:
        circle = patches.Circle(a['pos'], a['Rcut'], color='blue', alpha=0.05)
        ax[0].add_patch(circle)
        color = 'black' if a['type'] == 'C' else 'gray'
        ax[0].scatter(a['pos'][0], a['pos'][1], c=color, s=50, zorder=3)

    # Draw fine grid lines
    for t in ticks_x:
        ax[0].axvline(t, color='black', lw=0.5, alpha=0.2)
    for t in ticks_y:
        ax[0].axhline(t, color='black', lw=0.5, alpha=0.2)
    # Draw macro grid lines
    for t in np.arange(bb_min[0], bb_max[0] + m_res, m_res):
        ax[0].axvline(t, color='green', lw=1.0, alpha=0.3)
    for t in np.arange(bb_min[1], bb_max[1] + m_res, m_res):
        ax[0].axhline(t, color='green', lw=1.0, alpha=0.3)

    ax[0].set_xlim(bb_min[0], bb_max[0])
    ax[0].set_ylim(bb_min[1], bb_max[1])
    ax[0].set_aspect('equal')
    # annotate atom counts
    for i, x in enumerate(ticks_x):
        for j, y in enumerate(ticks_y):
            val = atom_counts[j, i]
            if val > 0:
                ax[0].text(x + b_res/2, y + b_res/2, int(val), color='red', ha='center', va='center', fontsize=7)

    # Plot 2: Heatmap of pair workload
    im = ax[1].imshow(workload, origin='lower', extent=[bb_min[0], bb_max[0], bb_min[1], bb_max[1]], cmap='Grays')
    plt.colorbar(im, ax=ax[1], label='Number of Density Matrix Pairs ($N_e$)')
    ax[1].set_title("Compute Cost Heatmap ($O(N_e)$)")

    # Annotate values on the heatmap
    for i, x in enumerate(ticks_x):
        for j, y in enumerate(ticks_y):
            val = workload[j, i]
            if val > 0:
                ax[1].text(x + b_res/2, y + b_res/2, int(val), color='r', ha='center', va='center', fontsize=7)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()