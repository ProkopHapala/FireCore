import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def get_system(n_atoms=20, size=20, rmin=2, rmax=2):
    # Random atoms for a larger "diffuse" system
    atoms = []
    for _ in range(n_atoms):
        pos = np.random.uniform(-size/2, size/2, 2)
        rcut = np.random.uniform(rmin, rmax)
        atoms.append({'pos': pos, 'Rcut': rcut})
    return atoms

def check_overlap(center, radius, b_min, b_max):
    closest_p = np.clip(center, b_min, b_max)
    return np.sum((center - closest_p)**2) < (radius**2)

def run_simulation(natom=10, size=24.0, block_res=None, macro_res=None, rmin=2.0, rmax=2.0, margin=0.1, seed=None, verbose=True):
    if seed is not None:
        np.random.seed(seed)
    atoms = get_system(natom, size, rmin=rmin, rmax=rmax)
    max_rcut = max(a['Rcut'] for a in atoms)
    auto_block_res = 2.0 * max_rcut + margin
    if block_res is None:
        block_res = auto_block_res
        if verbose: print(f"[DEBUG] block_res auto-set to {block_res} (= 2*max_rcut + margin)")
    elif block_res < auto_block_res:
        if verbose: print(f"[DEBUG] block_res bumped from {block_res} to {auto_block_res} to fit sphere diameter+margin")
        block_res = auto_block_res

    if macro_res is None:
        macro_res = 2.0 * block_res
        if verbose: print(f"[DEBUG] macro_res auto-set to {macro_res} (= 2 * block_res)")
    elif macro_res < block_res:
        if verbose: print(f"[DEBUG] macro_res bumped from {macro_res} to {block_res} to not be smaller than block_res")
        macro_res = block_res

    ratio = macro_res / block_res
    if not np.isclose(ratio, round(ratio)):
        macro_res = block_res * np.ceil(ratio)
        if verbose: print(f"[DEBUG] macro_res adjusted to {macro_res} to be integer multiple of block_res")

    if verbose:
        print(f"[DEBUG] natom={natom} size={size} block_res={block_res} macro_res={macro_res} rmin={rmin} rmax={rmax} max_rcut={max_rcut} margin={margin}")
    assert block_res >= auto_block_res, f"block_res={block_res} must be >= sphere diameter+margin={auto_block_res} to guarantee <=2x2 blocks per sphere"

    # 1. NAIVE SEARCH COUNT
    n_blocks_side = int(size / block_res)
    naive_checks = len(atoms) * (n_blocks_side**2)
    
    # 2. HIERARCHICAL SEARCH
    # Step A: Atom -> Macro Mapping
    macro_grid = {} # key: (ix, iy), val: [atom_indices]
    for i, a in enumerate(atoms):
        # Find macro-cells covered by the atom's AABB
        x_range = [int((a['pos'][0]-a['Rcut'] + size/2)/macro_res), 
                   int((a['pos'][0]+a['Rcut'] + size/2)/macro_res)]
        y_range = [int((a['pos'][1]-a['Rcut'] + size/2)/macro_res), 
                   int((a['pos'][1]+a['Rcut'] + size/2)/macro_res)]
        
        for ix in range(x_range[0], x_range[1] + 1):
            for iy in range(y_range[0], y_range[1] + 1):
                macro_grid.setdefault((ix, iy), []).append(i)
    
    # Step B: Refine Fine Blocks
    hierarchical_checks = 0
    active_blocks = []  # entries: (bx, by, count)
    
    for (m_ix, m_iy), atom_indices in macro_grid.items():
        # Only check fine blocks inside this occupied macro cell
        assert macro_res % block_res == 0, f"macro_res must be divisible by block_res (got macro_res={macro_res}, block_res={block_res})"
        sub_blocks_per_macro = int(macro_res / block_res)
        for fx in range(sub_blocks_per_macro):
            for fy in range(sub_blocks_per_macro):
                hierarchical_checks += len(atom_indices)
                
                # Calculate fine block bounds
                bx = m_ix * macro_res + fx * block_res - size/2
                by = m_iy * macro_res + fy * block_res - size/2
                b_min, b_max = np.array([bx, by]), np.array([bx+block_res, by+block_res])
                
                # Check actual intersection and count overlaps
                overlaps = [i for i in atom_indices if check_overlap(atoms[i]['pos'], atoms[i]['Rcut'], b_min, b_max)]
                if verbose:
                    print(f"[BLOCK] macro=({m_ix},{m_iy}) fine=({fx},{fy}) pos=({bx:.2f},{by:.2f}) count={len(overlaps)}")
                if overlaps:
                    active_blocks.append((bx, by, len(overlaps)))

    print(f"Naive Intersection Checks: {naive_checks}")
    print(f"Hierarchical Intersection Checks: {hierarchical_checks}")
    print(f"Efficiency Gain: {naive_checks / hierarchical_checks:.1f}x")
    
    return atoms, active_blocks, macro_res, block_res, size

# atoms, active_blocks, m_res, b_res, size = run_simulation()

def parse_args():
    p = argparse.ArgumentParser(description="Hierarchical grid culling demo (sphere vs AABB).")
    p.add_argument("--natom", type=int, default=10, help="Number of random atoms")
    p.add_argument("--size", type=float, default=24.0, help="Simulation box size")
    p.add_argument("--block-res", type=float, default=None, dest="block_res", help="Fine block size; default auto = 2*maxRcut + margin")
    p.add_argument("--macro-res", type=float, default=None, dest="macro_res", help="Macro cell size; default auto = 2 * block-res")
    p.add_argument("--margin", type=float, default=0.1, help="Extra padding added to sphere diameter for block sizing")
    p.add_argument("--rmin", type=float, default=2.0, help="Minimum atom cutoff (Rcut)")
    p.add_argument("--rmax", type=float, default=2.0, help="Maximum atom cutoff (Rcut)")
    p.add_argument("--seed", type=int, default=None, help="Random seed for reproducibility")
    p.add_argument("--no-plot", action="store_true", help="Disable plotting (CLI only metrics)")
    p.add_argument("--quiet", action="store_true", help="Reduce debug prints")
    return p.parse_args()

def main():
    args = parse_args()
    atoms, active_blocks, m_res, b_res, size = run_simulation(
        natom=args.natom,
        size=args.size,
        block_res=args.block_res,
        macro_res=args.macro_res,
        rmin=args.rmin,
        rmax=args.rmax,
        seed=args.seed,
        verbose=not args.quiet,
    )
    if args.no_plot:
        return

    fig, ax = plt.subplots(figsize=(8,8))
    ax.set_aspect('equal')
    ax.set_xlim(-size/2, size/2); ax.set_ylim(-size/2, size/2)
    
    # Draw Macro Grid
    for x in np.arange(-size/2, size/2 + m_res, m_res): ax.axvline(x, color='blue', alpha=0.3, lw=2)
    for y in np.arange(-size/2, size/2 + m_res, m_res): ax.axhline(y, color='blue', alpha=0.3, lw=2)
    
    # Draw Atoms
    for a in atoms:
        ax.add_patch(patches.Circle(a['pos'], a['Rcut'], color='green', alpha=0.1, zorder=1))
        ax.scatter(*a['pos'], color='black', s=10)
    
    # Draw Active Fine Blocks
    for bx, by, cnt in active_blocks:
        ax.add_patch(patches.Rectangle((bx, by), b_res, b_res, color='red', alpha=0.4, zorder=2))
        ax.text(bx + b_res/2, by + b_res/2, str(cnt), color='black', ha='center', va='center', fontsize=8, zorder=3)
    
    plt.title("Hierarchical Culling: Blue=Macro Grid, Red=Active Workgroup Blocks")
    plt.show()

if __name__ == "__main__":
    main()