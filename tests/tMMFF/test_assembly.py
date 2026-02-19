import sys
import os
import argparse
import time
import numpy as np
import matplotlib.pyplot as plt

# Add FireCore to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from pyBall.AtomicSystem import AtomicSystem
from pyBall import elements
from pyBall import plotUtils as pu
from pyBall.OCL.Assembly import AssemblyOCL, generate_transform_buffer, generate_transform_buffer_simple

def pack_atoms_with_radii(mol, radius_override=None):
    natoms = mol.natoms
    atoms = np.zeros((natoms, 4), dtype=np.float32)
    atoms[:, :3] = mol.apos
    if radius_override is not None:
        atoms[:, 3] = radius_override
    else:
        for i, ename in enumerate(mol.enames):
            elem = ename.split('_')[0] if '_' in ename else ename
            if elem in elements.ELEMENT_DICT:
                try:
                    atoms[i, 3] = elements.ELEMENT_DICT[elem][7]
                except:
                    atoms[i, 3] = 1.0
            else:
                atoms[i, 3] = 1.0
    return atoms

def parse_lvs(lvs_str):
    parts = lvs_str.replace('lvs', '').split()
    vals = [float(x) for x in parts]
    return np.array(vals).reshape(3, 3)

def main():
    parser = argparse.ArgumentParser(description="Test rigid body assembly using OpenCL")
    parser.add_argument('--preset',  type=str, choices=['small', 'helicene', 'triptycene', 'none'], default='none', help='Use preset configuration')
    parser.add_argument('--mol',     type=str,                 help='Path to molecule file (.xyz, .mol2)')
    parser.add_argument('--cell',    type=str,                 help='Lattice vectors string "lvs ax ay az bx by bz cx cy cz"')
    parser.add_argument('--nrot',    type=int,   default=500,  help='Number of rotations (Super-Fibonacci)')
    parser.add_argument('--nshift',  type=int,   default=10,    help='Number of shifts per 2D cell axis')
    parser.add_argument('--shift_range', type=float, default=0.5, help='Range of fractional shifts [-val, val]')
    parser.add_argument('--penalty', type=float, default=50.0, help='Maximum clash penalty for early exit')
    parser.add_argument('--zpenalty',type=float, default=2.0,  help='Penalty weight for Z-span to favor flat molecules')
    parser.add_argument('--zspan_max',type=float, default=10.0, help='Max allowed z-span for prefiltering and export')
    parser.add_argument('--clash_max',type=float, default=10.0, help='Max allowed clash penalty for export')
    parser.add_argument('--radius',  type=float, default=1.0,  help='Override collision radius for all atoms (default 1.0 A)')
    parser.add_argument('--wg',      type=int,   default=128,  help='Workgroup size')
    parser.add_argument('--device',  type=int,   default=0,    help='OpenCL device index')
    parser.add_argument('--dump',    action='store_true',      help='Dump best configuration to XYZ and render plots')
    parser.add_argument('--simple',  action='store_true',      help='Use simple 1-replica transform generation for smoke testing')
    parser.add_argument('--top_k',   type=int,   default=10,   help='Number of top configurations to save in movie')
    parser.add_argument('--outdir',  type=str,   default=None, help='Output directory (default tests/tMMFF/assembly_out_<preset>)')
    
    args = parser.parse_args()
    
    mol_path = args.mol
    cell_lvs = None
    
    if args.preset == 'small':
        mol_path = os.path.join(os.path.dirname(__file__), '../tUFF/data/xyz/HCCOCCH.xyz')
        cell_lvs = np.array([[15.0, 0, 0], [0, 15.0, 0], [0, 0, 15.0]])
        args.nrot = 10
        args.nshift = 4
    elif args.preset == 'helicene':
        mol_path = os.path.join(os.path.dirname(__file__), 'DiTetraceno_helicene_1a.xyz')
        cell_lvs = parse_lvs("32.7 0 0   16.35 28.319 0   0 0 40")
    elif args.preset == 'triptycene':
        mol_path = os.path.join(os.path.dirname(__file__), 'DiTriptyceno_helicene_3a.xyz')
        cell_lvs = parse_lvs("32.7 0 0   16.35 28.319 0   0 0 40")
        
    if mol_path is None:
        print("Please provide --mol or use a --preset")
        sys.exit(1)
        
    if not os.path.exists(mol_path):
        print(f"Error: Molecule file not found: {mol_path}")
        sys.exit(1)
        
    print(f"Loading molecule: {mol_path}")
    mol = AtomicSystem(fname=mol_path)
    
    if args.cell:
        cell_lvs = parse_lvs(args.cell)
    elif mol.lvec is not None:
        cell_lvs = mol.lvec
        
    if cell_lvs is None:
        print("Warning: No lattice vectors provided. Using default 20x20x20 box.")
        cell_lvs = np.array([[20.0, 0, 0], [0, 20.0, 0], [0, 0, 20.0]])
        
    print(f"Lattice Vectors:\n{cell_lvs}")
    
    base_atoms = pack_atoms_with_radii(mol, radius_override=args.radius)
    
    print(f"Generating transforms... (nrot={args.nrot}, nshift={args.nshift}, shift_range=[- {args.shift_range}, {args.shift_range}])")
    t0 = time.time()
    
    if args.simple:
        transforms, n_confs, R_conf, T_conf = generate_transform_buffer_simple(cell_lvs[0], cell_lvs[1], args.nrot, args.nshift)
        nmols = 1
    else:
        # Pass shift limits to generator if possible, else modify generator to accept it. We'll modify it.
        from pyBall.OCL.Assembly import R
        def custom_generate_transform_buffer(lattice_a, lattice_b, n_rot, n_shift, shift_range):
            angles_60 = np.linspace(0, 2 * np.pi, 6, endpoint=False)
            S_sym = R.from_euler('z', angles_60).as_matrix()
            u = np.array([0, -1, 1])
            uu, vv = np.meshgrid(u, u)
            L_lat = np.zeros((9, 3))
            L_lat[:, 0] = uu.flatten() * lattice_a[0] + vv.flatten() * lattice_b[0]
            L_lat[:, 1] = uu.flatten() * lattice_a[1] + vv.flatten() * lattice_b[1]
            L_lat[:, 2] = uu.flatten() * lattice_a[2] + vv.flatten() * lattice_b[2]
            
            from pyBall.OCL.Assembly import super_fibonacci_rotations, pack_transforms
            R_base = super_fibonacci_rotations(n_rot)
            
            fa = np.linspace(-shift_range, shift_range, n_shift)
            fb = np.linspace(-shift_range, shift_range, n_shift)
            FA, FB = np.meshgrid(fa, fb, indexing='ij')
            
            T_base = np.zeros((n_shift**2, 3))
            T_base[:, 0] = FA.flatten() * lattice_a[0] + FB.flatten() * lattice_b[0]
            T_base[:, 1] = FA.flatten() * lattice_a[1] + FB.flatten() * lattice_b[1]
            T_base[:, 2] = FA.flatten() * lattice_a[2] + FB.flatten() * lattice_b[2]
            
            N_rot = len(R_base)
            N_shift = len(T_base)
            N_confs = N_rot * N_shift
            
            R_conf = np.repeat(R_base, N_shift, axis=0)
            T_conf = np.tile(T_base, (N_rot, 1))
            
            R_sym = np.einsum('kij,cjl->ckil', S_sym, R_conf)
            T_sym = np.einsum('kij,cj->cki', S_sym, T_conf)
            
            R_all = np.tile(R_sym, (1, 9, 1, 1, 1))
            T_all = np.zeros((N_confs, 9, 6, 3))
            for l in range(9):
                T_all[:, l, :, :] = T_sym + L_lat[l].reshape(1, 1, 3)
                
            R_all = R_all.reshape(N_confs, 54, 3, 3)
            T_all = T_all.reshape(N_confs, 54, 3)
            
            cl_transforms = pack_transforms(R_all, T_all)
            return cl_transforms, N_confs, R_conf, T_conf

        transforms, n_confs, R_conf, T_conf = custom_generate_transform_buffer(cell_lvs[0], cell_lvs[1], args.nrot, args.nshift, args.shift_range)
        nmols = 54
        
    t1 = time.time()
    print(f"Generated {n_confs} configurations ({n_confs * nmols} total replicas) in {t1-t0:.3f} s")
    
    z_coords = np.einsum('cij,aj->cai', R_conf, mol.apos)
    z_spans = z_coords[:, :, 2].max(axis=1) - z_coords[:, :, 2].min(axis=1)
    z_ok = z_spans < args.zspan_max
    n_pre = int(np.sum(z_ok))
    print(f"Z-span prefilter: keep {n_pre}/{len(z_spans)} configs with zspan<{args.zspan_max}")

    transforms = transforms[z_ok]
    R_conf     = R_conf[z_ok]
    T_conf     = T_conf[z_ok]
    z_spans    = z_spans[z_ok]
    n_confs    = transforms.shape[0]

    print("Initializing OpenCL AssemblyOCL...")
    ocl = AssemblyOCL(nloc=args.wg, device_index=args.device)
    ocl.upload_base_atoms(base_atoms)

    print(f"Running evaluate_packing_3d kernel on {n_confs} configs...")
    t0 = time.time()
    scores = ocl.evaluate_packing(transforms, nmols=nmols, max_clash_penalty=args.penalty)
    t1 = time.time()

    print(f"Kernel execution finished in {t1-t0:.3f} s")
    
    # Total composite score
    total_scores = scores + args.zpenalty * z_spans
    
    valid_mask = scores < args.penalty
    valid_scores = scores[valid_mask]
    
    if len(valid_scores) > 0:
        print(f"Found {len(valid_scores)} valid configurations (clash < {args.penalty})")
        print(f"Min clash score: {np.min(valid_scores):.4f}")
        print(f"Median clash score: {np.median(valid_scores):.4f}")
    else:
        print(f"No valid configurations found. All configurations exceeded max penalty {args.penalty}.")
        
    export_mask = (z_spans < args.zspan_max) & (scores < args.clash_max)
    export_indices = np.where(export_mask)[0]
    if export_indices.size > 0:
        export_sorted = export_indices[np.argsort(total_scores[export_indices])]
        best_indices = export_sorted[:args.top_k]
    elif len(valid_scores) == 0:
        best_indices = [int(np.argmin(total_scores))]
    else:
        valid_indices = np.where(valid_mask)[0]
        sorted_valid = valid_indices[np.argsort(total_scores[valid_indices])]
        best_indices = sorted_valid[:args.top_k]
        
    best_idx = best_indices[0]
    print(f"Best configuration index: {best_idx} with Clash {scores[best_idx]:.4f}, Z-span {z_spans[best_idx]:.4f}, Total Score {total_scores[best_idx]:.4f}")
    
    if args.outdir is None:
        args.outdir = os.path.join(os.path.dirname(__file__), f"assembly_out_{args.preset}")
    os.makedirs(args.outdir, exist_ok=True)

    if args.dump:
        # Pareto Plot
        plt.figure(figsize=(8, 6))
        plot_mask = export_mask
        sc = plt.scatter(scores[plot_mask], z_spans[plot_mask], c=total_scores[plot_mask], cmap='viridis', alpha=0.6, s=10)
        plt.colorbar(sc, label=f'Total Score (clash + {args.zpenalty}*zspan)')
        plt.xlabel('Clash Penalty')
        plt.ylabel('Z-Span (A)')
        plt.title('Pareto Frontier: Clash vs Z-Span')
        plt.scatter(scores[best_idx], z_spans[best_idx], color='red', marker='x', s=150, linewidths=3, label='Best Config')
        plt.xlim(-1, min(args.penalty, float(np.max(scores[plot_mask]) * 1.1) if np.any(plot_mask) else args.penalty))
        plt.legend()
        pareto_name = os.path.join(args.outdir, f'pareto_{args.preset}.png')
        plt.savefig(pareto_name)
        print(f"Saved Pareto plot to {pareto_name}")
        plt.close()
        
        # Dump configuration movie
        out_name = os.path.join(args.outdir, f"assembly_best_{args.preset}_movie.xyz")
        export_sorted = export_indices[np.argsort(total_scores[export_indices])] if export_indices.size > 0 else np.array(best_indices, dtype=int)
        print(f"Dumping {len(export_sorted)} configurations (zspan<{args.zspan_max}, clash<{args.clash_max}) to {out_name}...")
        with open(out_name, 'w') as f:
            for i_rank, idx in enumerate(export_sorted):
                best_transforms = transforms[idx]
                out_atoms = ocl.emit_configuration(best_transforms, nmols)
                
                r_flat = R_conf[idx].flatten()
                t_vec = T_conf[idx]
                transform_str = " ".join([f"{x:.4f}" for x in r_flat]) + " " + " ".join([f"{x:.4f}" for x in t_vec])
                
                f.write(f"{len(out_atoms)}\n")
                f.write(f"# idx={idx} rank={i_rank} clash={scores[idx]:.4f} zspan={z_spans[idx]:.4f} transform: {transform_str}\n")
                for i in range(len(out_atoms)):
                    orig_idx = i % mol.natoms
                    elem = mol.enames[orig_idx]
                    f.write(f"{elem} {out_atoms[i,0]:.6f} {out_atoms[i,1]:.6f} {out_atoms[i,2]:.6f}\n")
        print(f"Saved movie to {out_name}")
        
        # Render Atoms Plot for the very best config
        plt.figure(figsize=(10, 10))
        best_transforms = transforms[best_idx]
        out_atoms = ocl.emit_configuration(best_transforms, nmols)
        apos = out_atoms[:, :3]
        elems = [mol.enames[i % mol.natoms] for i in range(len(apos))]
        
        # Calculate bonds for plotting
        if mol.bonds is None:
            mol.findBonds(Rcut=3.0, RvdwCut=0.5)
            
        try:
            # We want small dots for atoms and thin lines for bonds
            # The plotUtils functions might not support lw out of the box easily without modification,
            # so we use matplotlib directly here for custom styling if needed, or pass args
            pu.plotAtoms(apos, es=np.array(elems), sizes=4.0, marker='.') # small dots
            
            if mol.bonds is not None and len(mol.bonds) > 0:
                # Need to map the bonds to all replicas
                all_links = []
                for k in range(nmols):
                    offset = k * mol.natoms
                    for b in mol.bonds:
                        all_links.append([b[0]+offset, b[1]+offset])
                pu.plotBonds(links=all_links, ps=apos, axes=(0,1), colors='k', lws=0.5)
            
            # Plot 3x3 lattice frame
            pts = np.array([
                [-cell_lvs[0,0]-cell_lvs[1,0], -cell_lvs[0,1]-cell_lvs[1,1]],
                [ 2*cell_lvs[0,0]-cell_lvs[1,0],  2*cell_lvs[0,1]-cell_lvs[1,1]],
                [ 2*cell_lvs[0,0]+2*cell_lvs[1,0], 2*cell_lvs[0,1]+2*cell_lvs[1,1]],
                [-cell_lvs[0,0]+2*cell_lvs[1,0], -cell_lvs[0,1]+2*cell_lvs[1,1]],
                [-cell_lvs[0,0]-cell_lvs[1,0], -cell_lvs[0,1]-cell_lvs[1,1]]
            ])
            plt.plot(pts[:,0], pts[:,1], 'k--', lw=2, label='3x3 Supercell Bounds')
            
            # Plot origin cell
            cpts = np.array([[0,0], cell_lvs[0,:2], cell_lvs[0,:2]+cell_lvs[1,:2], cell_lvs[1,:2], [0,0]])
            plt.plot(cpts[:,0], cpts[:,1], 'b-', lw=2, label='Unit Cell (0,0)')
            
            plt.axis('equal')
            plt.title(f'Best Config {best_idx} (Clash: {scores[best_idx]:.2f}, Z-span: {z_spans[best_idx]:.2f})')
            plt.legend()
            img_name = os.path.join(args.outdir, f'assembly_best_{args.preset}.png')
            plt.savefig(img_name, dpi=300)
            print(f"Saved configuration 2D plot to {img_name}")
            plt.close()
        except Exception as e:
            print(f"Warning: Failed to render 2D plot: {e}")

if __name__ == "__main__":
    main()
