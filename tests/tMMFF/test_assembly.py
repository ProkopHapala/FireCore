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
from pyBall.OCL.Assembly import AssemblyOCL, generate_transform_buffer, generate_transform_buffer_simple, pack_transforms
from scipy.spatial.transform import Rotation as R

def prune_duplicate_rotations_sym(R_arr, sym_mats, tol_rad, progress=False):
    """Remove rotational duplicates considering a set of symmetry operations.
    Warning: O(N^2) and slow for large grids; disable with --dedup=false.
    """
    keep = []
    report_every = max(500, len(R_arr)//10) if progress else len(R_arr)+1
    for i, R_new in enumerate(R_arr):
        duplicate = False
        for R_ref in keep:
            for S in sym_mats:
                # Check equivalence: R_ref ~ S @ R_new
                M = R_ref.T @ (S @ R_new)
                cosang = 0.5 * (np.trace(M) - 1.0)
                cosang = np.clip(cosang, -1.0, 1.0)
                ang = np.arccos(cosang)
                if ang < tol_rad:
                    duplicate = True
                    break
            if duplicate:
                break
        if not duplicate:
            keep.append(R_new)
        if progress and (i+1) % report_every == 0:
            print(f"[rot] dedup progress: {i+1}/{len(R_arr)} processed, kept {len(keep)}")
    return np.array(keep)


def generate_rotations(mode, nrot, tilt_range=0.1, n_tilt=3, rot_tol=1e-2, progress=False, do_dedup=True):
    if mode == 'full3d':
        from pyBall.OCL.Assembly import super_fibonacci_rotations
        R_raw = super_fibonacci_rotations(nrot)
    elif mode == 'inplane':
        angles = np.linspace(0, 2*np.pi, nrot, endpoint=False)
        R_raw = R.from_euler('z', angles).as_matrix()
    elif mode == 'tilt':
        angles_z = np.linspace(0, 2*np.pi, nrot, endpoint=False)
        tilt_x   = np.linspace(-tilt_range, tilt_range, n_tilt)
        tilt_y   = np.linspace(-tilt_range, tilt_range, n_tilt)
        # yaw (Z) should be fastest index -> X, Y, Z with indexing='ij'
        X, Y, Z = np.meshgrid(tilt_x, tilt_y, angles_z, indexing='ij')
        euler = np.column_stack((X.flatten(), Y.flatten(), Z.flatten()))
        R_raw = R.from_euler('xyz', euler).as_matrix()
    else:
        raise ValueError(f"Unknown rotation mode: {mode}")
    if not do_dedup:
        if progress:
            print(f"[rot] dedup disabled (--no-dedup); using {len(R_raw)} raw rotations")
        return R_raw

    angles_60 = np.linspace(0, 2*np.pi, 6, endpoint=False)
    S_sym = R.from_euler('z', angles_60).as_matrix()
    R_unique = prune_duplicate_rotations_sym(R_raw, S_sym, rot_tol, progress=progress)
    if len(R_unique) < len(R_raw):
        print(f"Rotation dedup: kept {len(R_unique)} / {len(R_raw)} (tol={rot_tol} rad) considering 6-fold Z-symmetry")
    return R_unique

def plot_translations(T_conf, cell_lvs, outpath):
    plt.figure(figsize=(6,6))
    T_unique = np.unique(T_conf, axis=0)
    plt.scatter(T_unique[:,0], T_unique[:,1], s=10, c='b')
    cpts = np.array([[0,0], cell_lvs[0,:2], cell_lvs[0,:2]+cell_lvs[1,:2], cell_lvs[1,:2], [0,0]])
    plt.plot(cpts[:,0], cpts[:,1], 'k--', lw=2, label='Unit Cell')
    plt.axis('equal')
    plt.title(f'Translation Sampling ({len(T_unique)} unique shifts)')
    plt.legend()
    plt.savefig(outpath)
    plt.close()

def plot_rotations(R_conf, outpath):
    fig = plt.figure(figsize=(15, 5))
    R_u = np.unique(R_conf, axis=0)
    axes_names = ['a (x-axis)', 'b (y-axis)', 'c (z-axis)']
    colors = ['r', 'g', 'b']
    
    for i in range(3):
        ax = fig.add_subplot(1, 3, i+1, projection='3d')
        vecs = R_u[:, :, i]
        ax.scatter(vecs[:,0], vecs[:,1], vecs[:,2], c=colors[i], s=5)
        ax.set_xlim([-1.1, 1.1]); ax.set_ylim([-1.1, 1.1]); ax.set_zlim([-1.1, 1.1])
        ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
        ax.set_title(f'Basis {axes_names[i]}')
        
    plt.suptitle(f'Rotation Sampling ({len(R_u)} unique orientations)')
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()

def make_toy_molecule_movie(R_base, outpath, scale=2.0):
    with open(outpath, 'w') as f:
        for i, R_mat in enumerate(R_base):
            f.write("4\n")
            f.write(f"Toy molecule orientation {i}\n")
            f.write("C 0.000 0.000 0.000\n")
            ax = R_mat[:,0] * scale; f.write(f"O {ax[0]:.3f} {ax[1]:.3f} {ax[2]:.3f}\n")
            ay = R_mat[:,1] * scale; f.write(f"N {ay[0]:.3f} {ay[1]:.3f} {ay[2]:.3f}\n")
            az = R_mat[:,2] * scale; f.write(f"F {az[0]:.3f} {az[1]:.3f} {az[2]:.3f}\n")

def make_sym_toy_movie(R_base, cell_lvs, outpath, scale=2.0, replicas=6):
    angles_60 = np.linspace(0, 2 * np.pi, 6, endpoint=False)
    S_sym = R.from_euler('z', angles_60).as_matrix()
    
    T0 = 0.5 * cell_lvs[0] + 0.5 * cell_lvs[1]
    
    u = np.array([0, -1, 1])
    uu, vv = np.meshgrid(u, u)
    L_lat = np.zeros((9, 3))
    L_lat[:, 0] = uu.flatten() * cell_lvs[0, 0] + vv.flatten() * cell_lvs[1, 0]
    L_lat[:, 1] = uu.flatten() * cell_lvs[0, 1] + vv.flatten() * cell_lvs[1, 1]
    L_lat[:, 2] = uu.flatten() * cell_lvs[0, 2] + vv.flatten() * cell_lvs[1, 2]
    
    with open(outpath, 'w') as f:
        for i, R_mat in enumerate(R_base):
            atoms_to_write = []
            
            # 6-fold symmetry
            for s in range(6):
                S = S_sym[s]
                R_s = S @ R_mat
                T_s = S @ T0
                
                if replicas == 54:
                    for l in range(9):
                        T_final = T_s + L_lat[l]
                        atoms_to_write.append((T_final, R_s))
                else:
                    atoms_to_write.append((T_s, R_s))
            
            f.write(f"{len(atoms_to_write) * 4}\n")
            f.write(f"Symmetric Toy molecule orientation {i}\n")
            for T_final, R_s in atoms_to_write:
                f.write(f"C {T_final[0]:.3f} {T_final[1]:.3f} {T_final[2]:.3f}\n")
                ax = T_final + R_s[:,0] * scale; f.write(f"O {ax[0]:.3f} {ax[1]:.3f} {ax[2]:.3f}\n")
                ay = T_final + R_s[:,1] * scale; f.write(f"N {ay[0]:.3f} {ay[1]:.3f} {ay[2]:.3f}\n")
                az = T_final + R_s[:,2] * scale; f.write(f"F {az[0]:.3f} {az[1]:.3f} {az[2]:.3f}\n")

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

'''

# 1) Triptycene, in-plane rotations (flat), with diagnostics (plots + toy XYZ), lenient z-span 
python3 test_assembly.py --preset triptycene   --rot_mode inplane --nrot 12 --nshift 4 --zspan_max 40 --plot_trans --plot_rot --plot_toy

# 2) Triptycene, tilt mode (slight out-of-plane), with diagnostics
python3 test_assembly.py --preset triptycene --rot_mode tilt --nrot 12 --n_tilt 3 --tilt_range 0.1 --nshift 4 --zspan_max 40 --plot_trans --plot_rot --plot_toy

# 3) Same as (2) but also dump Pareto, XYZ movie, and best-config PNG
python3 test_assembly.py --preset triptycene  --rot_mode tilt --nrot 12 --n_tilt 3 --tilt_range 0.1 --nshift 4 --zspan_max 40 --plot_trans --plot_rot --plot_toy --dump

# 4) Helicene example, full 3D sampling with stricter spans/clashes (optional)
python3 test_assembly.py --preset helicene --rot_mode full3d --nrot 800 --nshift 12 --shift_range 0.4 --zspan_max 10 --clash_max 5 --dump --plot_trans --plot_rot --plot_toy


python3 test_assembly.py --preset triptycene --rot_mode tilt --nrot 32 --n_tilt 7 --tilt_range 0.25 --nshift 20 --shift_region triangle --shift_sum_max 1.0 --zspan_max 40 --clash_max 2 --plot_trans --plot_rot --plot_best_k 10 --z_highlight 0.4 --dump

python3 test_assembly.py --preset triptycene --rot_mode tilt --nrot 16 --n_tilt 5 --tilt_range 0.4 --nshift 10 --clash_max 5.0 --plot_trans --plot_rot --plot_best_k 0 --z_highlight 0.4 --dump --nPBC_test 2 --nPBC_xyz 1

python3 test_assembly.py --preset triptycene --rot_mode tilt --nrot 32 --n_tilt 9 --tilt_range 0.35 --nshift 10 --clash_max 5.0 --plot_trans --plot_rot --plot_best_k 0 --z_highlight 0.4 --dump --nPBC_test 2 --nPBC_xyz 1

'''


def main():
    parser = argparse.ArgumentParser(description="Test rigid body assembly using OpenCL")
    parser.add_argument('--preset',  type=str, choices=['small', 'helicene', 'triptycene', 'none'], default='none', help='Use preset configuration')
    parser.add_argument('--mol',     type=str,                 help='Path to molecule file (.xyz, .mol2)')
    parser.add_argument('--cell',    type=str,                 help='Lattice vectors string "lvs ax ay az bx by bz cx cy cz"')
    parser.add_argument('--nrot',    type=int,   default=500,  help='Number of rotations')
    parser.add_argument('--rot_mode',type=str,   choices=['full3d', 'inplane', 'tilt'], default='full3d', help='Rotation sampling mode')
    parser.add_argument('--tilt_range',type=float, default=0.1, help='Range of tilt angles (radians) for tilt mode')
    parser.add_argument('--n_tilt',  type=int,   default=3,    help='Number of tilt steps per axis for tilt mode')
    parser.add_argument('--nshift',  type=int,   default=10,    help='Number of shifts per 2D cell axis')
    parser.add_argument('--shift_range', type=float, default=1.0, help='Range of fractional shifts [-val, val]')
    parser.add_argument('--shift_region', type=str, choices=['square', 'triangle'], default='triangle', help='Shape of the translation sampling region')
    parser.add_argument('--shift_sum_max', type=float, default=0.8001, help='Max sum of fractional coordinates (ca+cb) for triangle region')
    parser.add_argument('--penalty', type=float, default=50.0, help='Maximum clash penalty for early exit')
    parser.add_argument('--zpenalty',type=float, default=2.0,  help='Penalty weight for Z-span to favor flat molecules')
    parser.add_argument('--zspan_max',type=float, default=20.0, help='Max allowed z-span for prefiltering and export')
    parser.add_argument('--clash_max',type=float, default=10.0, help='Max allowed clash penalty for export')
    parser.add_argument('--export_max',type=int, default=500,  help='Maximum number of configurations to write to the movie to prevent massive dumps')
    parser.add_argument('--nPBC_test', type=int, default=2,    help='Number of PBC lattice replicas used in OpenCL evaluation (0=central cell, 1=3x3, 2=5x5)')
    parser.add_argument('--nPBC_xyz', type=int, default=1,     help='Number of PBC replicas to write to XYZ and PNGs')
    parser.add_argument('--dist_min', type=float, default=1.0, help='Minimum allowed closest-approach distance between molecules')
    parser.add_argument('--progress', action='store_true',     help='Print progress during pose generation to show it is running')
    parser.add_argument('--dedup',    type=int, default=0, help='Deduplicate rotations (slow, O(N^2))')
    parser.add_argument('--xyz_simple',action='store_true',    help='Strip elemental suffixes (like _ar) in XYZ output', default=True)
    parser.add_argument('--radius',  type=float, default=1.0,  help='Override collision radius for all atoms (default 1.0 A)')
    parser.add_argument('--wg',      type=int,   default=128,  help='Workgroup size')
    parser.add_argument('--device',  type=int,   default=0,    help='OpenCL device index')
    parser.add_argument('--dump',    action='store_true',      help='Dump best configuration to XYZ and render plots')
    parser.add_argument('--plot_trans', action='store_true',   help='Plot translation sampling')
    parser.add_argument('--plot_rot', action='store_true',     help='Plot rotation sampling on sphere')
    parser.add_argument('--plot_toy', action='store_true',     help='Generate toy 4-atom orientation movie')
    parser.add_argument('--plot_toy_sym', action='store_true', help='Generate symmetry-multiplied toy molecule movie for first config')
    parser.add_argument('--toy_sym_replicas', type=int, choices=[6, 54], default=6, help='Number of replicas for symmetry toy movie')
    parser.add_argument('--toy_scale', type=float, default=1.5,help='Scale of toy molecule basis vectors')
    parser.add_argument('--plot_best_k', type=int, default=1,  help='Number of top configurations to plot as 2D PNGs')
    parser.add_argument('--z_highlight', type=float, default=0.4, help='Z threshold from top to highlight atoms in red')
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
        mol_path = os.path.join(os.path.dirname(__file__), '../tMolGUIapp/DiTriptyceno_helicene_3a_xy.mol2')
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
        def custom_generate_transform_buffer(lattice_a, lattice_b, n_rot, n_shift, shift_range, rot_mode, tilt_range, n_tilt, shift_region, shift_sum_max, n_pbc_test):
            angles_60 = np.linspace(0, 2 * np.pi, 6, endpoint=False)
            S_sym = R.from_euler('z', angles_60).as_matrix()
            
            if n_pbc_test == 0:
                u = np.array([0])
            elif n_pbc_test == 1:
                u = np.array([0, -1, 1])
            elif n_pbc_test == 2:
                u = np.array([0, -1, 1, -2, 2])
            elif n_pbc_test == 3:
                u = np.array([0, -1, 1, -2, 2, -3, 3])
            else:
                raise ValueError(f"nPBC_test={n_pbc_test} not supported (use 0, 1, 2, 3)")
                
            uu, vv = np.meshgrid(u, u)
            n_lat = len(uu.flatten())
            L_lat = np.zeros((n_lat, 3))
            L_lat[:, 0] = uu.flatten() * lattice_a[0] + vv.flatten() * lattice_b[0]
            L_lat[:, 1] = uu.flatten() * lattice_a[1] + vv.flatten() * lattice_b[1]
            L_lat[:, 2] = uu.flatten() * lattice_a[2] + vv.flatten() * lattice_b[2]
            if args.progress:
                print(f"[gen] lattice grid: n_pbc_test={n_pbc_test}, cells={n_lat} n_rot={n_rot} n_tilt={n_tilt}  ({(2*n_pbc_test+1)}x{(2*n_pbc_test+1)})")
            

            R_base = generate_rotations(rot_mode, n_rot, tilt_range, n_tilt, progress=args.progress, do_dedup=args.dedup)
            if args.progress:
                print(f"[gen] rotations: {len(R_base)}")
            
            if shift_region == 'triangle':
                # Strict single-cell triangle: fractional coords in [0, shift_range], mask by fa+fb<=shift_sum_max
                fa = np.linspace(0.0, shift_range, n_shift)
                fb = np.linspace(0.0, shift_range, n_shift)
            else:
                fa = np.linspace(-shift_range, shift_range, n_shift)
                fb = np.linspace(-shift_range, shift_range, n_shift)
            FA, FB = np.meshgrid(fa, fb, indexing='ij')
            fa_flat = FA.flatten()
            fb_flat = FB.flatten()
            
            if shift_region == 'triangle':
                mask = (fa_flat >= 0) & (fb_flat >= 0) & ((fa_flat + fb_flat) <= shift_sum_max)
                fa_flat = fa_flat[mask]
                fb_flat = fb_flat[mask]
                if len(fa_flat) == 0:
                    raise ValueError(f"Shift region mask removed all translations! Try increasing --shift_range or --shift_sum_max")
            if args.progress:
                print(f"[gen] translations: {len(fa_flat)} kept")
            
            n_actual_shifts = len(fa_flat)
            T_base = np.zeros((n_actual_shifts, 3))
            T_base[:, 0] = fa_flat * lattice_a[0] + fb_flat * lattice_b[0]
            T_base[:, 1] = fa_flat * lattice_a[1] + fb_flat * lattice_b[1]
            T_base[:, 2] = fa_flat * lattice_a[2] + fb_flat * lattice_b[2]
            
            N_rot = len(R_base)
            N_shift = len(T_base)
            N_confs = N_rot * N_shift
            if args.progress:
                print(f"[gen] pose grid: N_rot={N_rot}, N_shift={N_shift}, N_confs={N_confs}")
            
            R_conf = np.repeat(R_base, N_shift, axis=0)
            T_conf = np.tile(T_base, (N_rot, 1))
            
            R_sym = np.einsum('kij,cjl->ckil', S_sym, R_conf)
            T_sym = np.einsum('kij,cj->cki', S_sym, T_conf)
            
            R_all = np.tile(R_sym[:, None, :, :, :], (1, n_lat, 1, 1, 1))
            T_all = np.zeros((N_confs, n_lat, 6, 3))
            for l in range(n_lat):
                T_all[:, l, :, :] = T_sym + L_lat[l].reshape(1, 1, 3)
                
            R_all = R_all.reshape(N_confs, 6*n_lat, 3, 3)
            T_all = T_all.reshape(N_confs, 6*n_lat, 3)
            if args.progress:
                print(f"[gen] sym/lat tiled: replicas_per_conf={6*n_lat}")
            
            cl_transforms = pack_transforms(R_all, T_all)
            return cl_transforms, N_confs, R_conf, T_conf, 6*n_lat

        transforms, n_confs, R_conf, T_conf, nmols_eval = custom_generate_transform_buffer(
            cell_lvs[0], cell_lvs[1], args.nrot, args.nshift, args.shift_range, 
            args.rot_mode, args.tilt_range, args.n_tilt, args.shift_region, args.shift_sum_max, args.nPBC_test
        )
        nmols = nmols_eval
        
        # Calculate requested output replicas
        if args.nPBC_xyz == 0:
            nmols_out_req = 6
        else:
            nmols_out_req = 6 * ((2 * args.nPBC_xyz + 1)**2)
            
        if nmols_out_req > nmols_eval:
            raise RuntimeError(f"--nPBC_xyz {args.nPBC_xyz} requires {nmols_out_req} replicas, "
                               f"but --nPBC_test {args.nPBC_test} only generated {nmols_eval} replicas. "
                               "Export replicas cannot exceed evaluation replicas. Increase nPBC_test or decrease nPBC_xyz.")
        
        # Build lattice grid consistent with nPBC_test
        if args.nPBC_test == 0:
            u = np.array([0])
        elif args.nPBC_test == 1:
            u = np.array([0, -1, 1])
        elif args.nPBC_test == 2:
            u = np.array([0, -1, 1, -2, 2])
        elif args.nPBC_test == 3:
            u = np.array([0, -1, 1, -2, 2, -3, 3])
        else:
            raise ValueError(f"nPBC_test={args.nPBC_test} not supported (use 0, 1, 2, 3)")
            
        uu, vv = np.meshgrid(u, u)
        uu_flat, vv_flat = uu.flatten(), vv.flatten()
        
        # Select centered block for export based on nPBC_xyz
        sel_lat = [i for i, (a, b) in enumerate(zip(uu_flat, vv_flat)) if max(abs(a), abs(b)) <= args.nPBC_xyz]
        expected_cells = (2 * args.nPBC_xyz + 1)**2
        if len(sel_lat) != expected_cells:
            raise RuntimeError(f"Expected {expected_cells} lattice cells for nPBC_xyz={args.nPBC_xyz}, got {len(sel_lat)}")
            
        sel_indices = []
        for l in sel_lat:
            sel_indices.extend([l*6 + s for s in range(6)])  # keep all 6 symmetry copies per selected cell
        nmols_out = len(sel_indices)
        
    t1 = time.time()
    print(f"Generated {n_confs} configurations ({n_confs * nmols} total replicas) in {t1-t0:.3f} s")
    
    if args.outdir is None:
        args.outdir = os.path.join(os.path.dirname(__file__), f"assembly_out_{args.preset}")
    os.makedirs(args.outdir, exist_ok=True)
    
    if args.plot_trans:
        trans_plot_path = os.path.join(args.outdir, f"trans_sampling_{args.preset}.png")
        plot_translations(T_conf, cell_lvs, trans_plot_path)
        print(f"Saved translation plot to {trans_plot_path}")
        
    if args.plot_rot:
        rot_plot_path = os.path.join(args.outdir, f"rot_sampling_{args.preset}.png")
        plot_rotations(R_conf, rot_plot_path)
        print(f"Saved rotation plot to {rot_plot_path}")
        
    if args.plot_toy:
        toy_plot_path = os.path.join(args.outdir, f"toy_rotations_{args.preset}.xyz")
        make_toy_molecule_movie(R_conf, toy_plot_path, scale=args.toy_scale)
        print(f"Saved toy molecule orientations to {toy_plot_path}")
    
    if args.plot_toy_sym:
        toy_sym_plot_path = os.path.join(args.outdir, f"toy_rotations_sym_{args.preset}.xyz")
        n_toy = min(args.top_k, len(R_conf))
        make_sym_toy_movie(R_conf[:n_toy], cell_lvs, toy_sym_plot_path, scale=args.toy_scale, replicas=args.toy_sym_replicas)
        print(f"Saved symmetric toy molecule orientations to {toy_sym_plot_path}")
    
    z_coords = np.einsum('cij,aj->cai', R_conf, mol.apos)
    z_spans = z_coords[:, :, 2].max(axis=1) - z_coords[:, :, 2].min(axis=1)
    z_ok = z_spans < args.zspan_max
    n_pre = int(np.sum(z_ok))
    print(f"Z-span prefilter: keep {n_pre}/{len(z_spans)} configs with zspan<{args.zspan_max}")

    if n_pre == 0:
        print(f"Error: All configurations were filtered out by --zspan_max {args.zspan_max}")
        print(f"Actual z-spans range from {z_spans.min():.2f} to {z_spans.max():.2f}")
        sys.exit(1)

    transforms = transforms[z_ok]
    R_conf     = R_conf[z_ok]
    T_conf     = T_conf[z_ok]
    z_spans    = z_spans[z_ok]
    n_confs    = transforms.shape[0]

    print("Initializing OpenCL AssemblyOCL...")
    ocl = AssemblyOCL(nloc=args.wg, device_index=args.device)
    ocl.upload_base_atoms(base_atoms)

    print(f"Running evaluate_packing_3d kernel on {n_confs} configs (nmols={nmols})...")
    t0 = time.time()
    scores, min_dists = ocl.evaluate_packing(transforms, nmols=nmols, max_clash_penalty=args.penalty)
    t1 = time.time()

    print(f"Kernel execution finished in {t1-t0:.3f} s")
    
    # Total composite score (using clash penalty as before)
    total_scores = scores + args.zpenalty * z_spans
    
    valid_mask = scores < args.penalty
    valid_scores = scores[valid_mask]
    
    if len(valid_scores) > 0:
        print(f"Found {len(valid_scores)} valid configurations (clash < {args.penalty})")
        print(f"Min clash score: {np.min(valid_scores):.4f}")
        print(f"Median clash score: {np.median(valid_scores):.4f}")
    else:
        print(f"No valid configurations found. All configurations exceeded max penalty {args.penalty}.")
        
    export_mask = (z_spans < args.zspan_max) & (scores < args.clash_max) & (min_dists >= args.dist_min)
    export_indices = np.where(export_mask)[0]
    
    print(f"Configurations passing export criteria (zspan < {args.zspan_max}, clash < {args.clash_max}, dist_min >= {args.dist_min}): {len(export_indices)}")
    
    if export_indices.size == 0:
        print(f"WARNING: No configurations passed the strict export criteria (zspan < {args.zspan_max}, clash < {args.clash_max}, dist_min >= {args.dist_min}).")
        print("         The Pareto plot will still be generated, but no XYZ movie or 2D plots will be exported.")
        export_sorted = np.array([], dtype=int)
        best_indices = np.array([], dtype=int)
        best_idx = int(np.argmin(total_scores)) # best overall, even if invalid, for plot reference
    else:
        export_sorted = export_indices[np.argsort(total_scores[export_indices])]
        export_sorted = export_sorted[:args.export_max] # Cap export size
        best_indices = export_sorted[:args.top_k]
        best_idx = best_indices[0]
        
    print(f"Best configuration overall: index {best_idx} with Clash {scores[best_idx]:.4f}, Min-Dist {min_dists[best_idx]:.4f}, Z-span {z_spans[best_idx]:.4f}, Total Score {total_scores[best_idx]:.4f}")
    
    if args.dump:
        
        # DOUBLE Pareto Plot (Clash vs Z-Span AND MinDist vs Z-Span)
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
        
        # --- Plot 1: Clash vs Z-span ---
        ax1.scatter(scores, z_spans, c='black', alpha=0.3, s=5, label='All Configs')
        if len(export_sorted) > 0:
            ax1.scatter(scores[export_sorted], z_spans[export_sorted], c='red', alpha=0.8, s=20, label='Exported')
            ax1.scatter(scores[best_idx], z_spans[best_idx], color='lime', marker='*', s=150, linewidths=2, label='Best Exported')
        else:
            ax1.scatter(scores[best_idx], z_spans[best_idx], color='orange', marker='x', s=150, linewidths=2, label='Best Overall (Failing)')
            
        pts1 = np.vstack((scores, z_spans)).T
        is_pareto1 = np.ones(pts1.shape[0], dtype=bool)
        for i, p in enumerate(pts1):
            if is_pareto1[i]:
                dominating = np.all(pts1 <= p, axis=1) & np.any(pts1 < p, axis=1)
                is_pareto1[i] = not np.any(dominating)
                
        pareto_idx1 = np.where(is_pareto1)[0]
        pareto_sorted1 = pareto_idx1[np.argsort(scores[pareto_idx1])]
        ax1.plot(scores[pareto_sorted1], z_spans[pareto_sorted1], 'b-', alpha=0.5, lw=2, label='Pareto Front')
        
        ax1.axvline(x=args.clash_max, color='r', linestyle='--', alpha=0.3, label='Clash Max')
        ax1.axhline(y=args.zspan_max, color='b', linestyle='--', alpha=0.3, label='Z-Span Max')
        
        ax1.set_xlabel('Clash Penalty (sum)')
        ax1.set_ylabel('Z-Span (A)')
        ax1.set_title('Clash vs Z-Span')
        
        max_clash_plot = max(args.clash_max * 2, np.percentile(scores, 10))
        max_zspan_plot = max(args.zspan_max * 1.2, np.percentile(z_spans, 50))
        ax1.set_xlim(-0.1, min(args.penalty, max_clash_plot))
        ax1.set_ylim(0, max_zspan_plot)
        ax1.legend()
        
        # --- Plot 2: Min-Dist vs Z-span ---
        # Note: For Min-Dist, we want to MAXIMIZE it, and MINIMIZE Z-span. 
        # So Pareto domination is reversed for Min-Dist (higher is better).
        ax2.scatter(min_dists, z_spans, c='black', alpha=0.3, s=5, label='All Configs')
        if len(export_sorted) > 0:
            ax2.scatter(min_dists[export_sorted], z_spans[export_sorted], c='red', alpha=0.8, s=20, label='Exported')
            ax2.scatter(min_dists[best_idx], z_spans[best_idx], color='lime', marker='*', s=150, linewidths=2, label='Best Exported')
        else:
            ax2.scatter(min_dists[best_idx], z_spans[best_idx], color='orange', marker='x', s=150, linewidths=2, label='Best Overall (Failing)')
            
        pts2 = np.vstack((-min_dists, z_spans)).T # Negate min_dist so we can minimize both for pareto logic
        is_pareto2 = np.ones(pts2.shape[0], dtype=bool)
        for i, p in enumerate(pts2):
            if is_pareto2[i]:
                dominating = np.all(pts2 <= p, axis=1) & np.any(pts2 < p, axis=1)
                is_pareto2[i] = not np.any(dominating)
                
        pareto_idx2 = np.where(is_pareto2)[0]
        # Sort by min_dist (descending) for line plotting
        pareto_sorted2 = pareto_idx2[np.argsort(-min_dists[pareto_idx2])]
        ax2.plot(min_dists[pareto_sorted2], z_spans[pareto_sorted2], 'b-', alpha=0.5, lw=2, label='Pareto Front')
        
        ax2.axvline(x=args.dist_min, color='r', linestyle='--', alpha=0.3, label='Min-Dist Min')
        ax2.axhline(y=args.zspan_max, color='b', linestyle='--', alpha=0.3, label='Z-Span Max')
        
        ax2.set_xlabel('Minimum Inter-Molecular Distance (A)')
        ax2.set_ylabel('Z-Span (A)')
        ax2.set_title('Min-Dist vs Z-Span')
        
        # Avoid huge ranges
        # Use a generous upper bound so high min-distance configs are visible
        max_dist_plot = max(np.percentile(min_dists, 99), args.dist_min * 2.0)
        ax2.set_xlim(0, max_dist_plot * 1.05)
        ax2.set_ylim(0, max_zspan_plot)
        ax2.legend()
        
        plt.tight_layout()
        pareto_name = os.path.join(args.outdir, f'pareto_{args.preset}.png')
        plt.savefig(pareto_name, dpi=150)
        print(f"Saved double Pareto plot to {pareto_name}")
        plt.close()
        
        # Dump configuration movie
        if len(export_sorted) > 0:
            out_name = os.path.join(args.outdir, f"assembly_best_{args.preset}_movie.xyz")
            print(f"Dumping {len(export_sorted)} configurations (capped at {args.export_max}) to {out_name} with {nmols_out} replicas...")
            with open(out_name, 'w') as f:
                for i_rank, idx in enumerate(export_sorted):
                    best_transforms = transforms[idx][sel_indices]
                    out_atoms = ocl.emit_configuration(best_transforms, nmols_out)
                    
                    r_flat = R_conf[idx].flatten()
                    t_vec = T_conf[idx]
                    transform_str = " ".join([f"{x:.4f}" for x in r_flat]) + " " + " ".join([f"{x:.4f}" for x in t_vec])
                    
                    f.write(f"{len(out_atoms)}\n")
                    f.write(f"# idx={idx} rank={i_rank} clash={scores[idx]:.4f} mindist={min_dists[idx]:.4f} zspan={z_spans[idx]:.4f} transform: {transform_str}\n")
                    for i in range(len(out_atoms)):
                        orig_idx = i % mol.natoms
                        elem = mol.enames[orig_idx]
                        if args.xyz_simple:
                            elem = elem.split('_')[0] if '_' in elem else elem
                        f.write(f"{elem} {out_atoms[i,0]:.6f} {out_atoms[i,1]:.6f} {out_atoms[i,2]:.6f}\n")
            print(f"Saved movie to {out_name}")
        
        # Render Atoms Plot for the top K configs
        plot_k = min(args.plot_best_k, len(best_indices))
        for rank in range(plot_k):
            idx = best_indices[rank]
            plt.figure(figsize=(10, 10))
            best_transforms = transforms[idx][sel_indices]
            out_atoms = ocl.emit_configuration(best_transforms, nmols_out)
            apos = out_atoms[:, :3]
            
            # Height-aware coloring and highlights
            z_vals = apos[:, 2]
            z_max = z_vals.max()
            z_min = z_vals.min()
            z_range = z_max - z_min if (z_max - z_min) > 0 else 1.0
            
            # Normalize z for alpha (top = 1.0, bottom = 0.1)
            alphas = 0.1 + 0.9 * (z_vals - z_min) / z_range
            
            # Identify highlight atoms
            highlight_mask = z_vals > (z_max - args.z_highlight)
            
            if mol.bonds is None:
                mol.findBonds(Rcut=3.0, RvdwCut=0.5)
                
            try:
                # Plot normal atoms
                for j in range(len(apos)):
                    if highlight_mask[j]:
                        plt.scatter(apos[j,0], apos[j,1], s=20, c='red', alpha=1.0, zorder=10, edgecolors='none')
                    else:
                        plt.scatter(apos[j,0], apos[j,1], s=10, c='gray', alpha=alphas[j], zorder=5, edgecolors='none')
                
                # Calculate bonds for plotting
                if mol.bonds is not None and len(mol.bonds) > 0:
                    all_links = []
                    for k in range(nmols_out):
                        offset = k * mol.natoms
                        for b in mol.bonds:
                            all_links.append([b[0]+offset, b[1]+offset])
                    # draw all bonds with low alpha
                    for link in all_links:
                        p1 = apos[link[0]]
                        p2 = apos[link[1]]
                        plt.plot([p1[0], p2[0]], [p1[1], p2[1]], color='black', alpha=0.3, lw=0.5, zorder=1)
                
                # Plot 3x3 lattice frame if full PBC, else just central cell bounds
                if args.nPBC_xyz == 1:
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
                plt.title(f'Rank {rank+1}: Config {idx} (Clash: {scores[idx]:.2f}, Z-span: {z_spans[idx]:.2f})')
                plt.legend()
                img_name = os.path.join(args.outdir, f'assembly_best_{args.preset}_{rank}.png')
                plt.savefig(img_name, dpi=300)
                print(f"Saved configuration 2D plot to {img_name}")
                plt.close()
            except Exception as e:
                print(f"Warning: Failed to render 2D plot: {e}")

if __name__ == "__main__":
    main()
