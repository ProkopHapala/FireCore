import sys
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
if ROOT not in sys.path:
    sys.path.append(ROOT)

from pyBall.AtomicSystem import AtomicSystem
from pyBall.OCL.ManipulationPathOpt import ManipulationPathOpt

def create_h2o():
    mol = AtomicSystem()
    mol.enames = ["O", "H", "H"]
    mol.apos = np.array([
        [0.0, 0.0, 0.0],
        [0.757, 0.586, 0.0],
        [-0.757, 0.586, 0.0]
    ], dtype=np.float32)
    mol.qs = np.array([-0.82, +0.41, +0.41], dtype=np.float32)
    mol.natoms = 3
    mol.atypes = [8, 1, 1]
    mol.bonds = [(0, 1), (0, 2)]
    mol.neighs() # builds ngs, nngs
    
    mol.atom_types_mmff = ["O_3", "H_", "H_"] # standard mmff names
    return mol

def parse_args():
    parser = argparse.ArgumentParser(description="Test ManipulationPathOpt with Asymmetric Causal Tethers")
    parser.add_argument("--mol",        type=str,    default="", help="Path to molecule .xyz or .mol2. If empty, uses built-in H2O.")
    parser.add_argument("--grid",       type=str,    default=os.path.join(ROOT, "tests/tUFF/data/NaCl_1x1_L3/Bspline_PLQd.npy"), help="Path to substrate GridFF Bspline_PLQd .npy. Empty or missing → vacuum.")
    parser.add_argument("--grid_p0",    type=float, nargs=3, default=[-2.0, -2.0, 0.0], help="Grid origin p0 (Angstrom), used for Bspline_PLQd")
    parser.add_argument("--grid_step",  type=float, nargs=3, default=[0.1, 0.1, 0.1], help="Grid step (Angstrom), used for Bspline_PLQd")
    parser.add_argument("--grid_tex",   type=int,   default=0,   help="Use 3D texture for GridFF sampling")
    parser.add_argument("--grid_alpha", type=float, default=2.2, help="alpha_morse used to factorize Morse in GridFF (see relax_multi.cl)")
    parser.add_argument("--z_start",    type=float, default=4.0, help="Initial molecule z shift above surface (Angstrom). Debugging default: 4.0")
    parser.add_argument("--anchor",     type=int,   default=1,   help="Index of tip anchor atom (default 1, which is H in H2O)")
    parser.add_argument("--plot_atom",  type=int,   default=2,   help="Index of atom to plot trajectory for (default 0, which is O in H2O)")
    parser.add_argument("--De",         type=float, default=1.0, help="Morse De for tip constraint")
    parser.add_argument("--r0",         type=float, default=2.0, help="Morse r0 for tip constraint")
    parser.add_argument("--alpha",      type=float, default=1.5, help="Morse alpha for tip constraint")
    parser.add_argument("--tip_dz",             type=float, default=0.0, help="Tip trajectory z offset relative to anchor atom initial z (Angstrom). Use 0.0 for 'hold atom at its initial height'.")
    parser.add_argument("--shift_xy",           type=float, nargs=2, default=[0.0, 0.0], help="Rigid shift of initial molecule in x,y (Angstrom).")
    parser.add_argument("--auto_xy",            type=int,   default=0, help="If 1, scan a small dx,dy grid and pick shift giving most attractive Fz on atom2 (H2) from GridFF before relaxation.")
    parser.add_argument("--auto_xy_span",       type=float, default=1.5, help="Half-span (Angstrom) of auto_xy scan in dx,dy.")
    parser.add_argument("--auto_xy_n",          type=int,   default=7, help="Number of samples per axis in auto_xy scan.")
    parser.add_argument("--scan_z",             type=int,   default=0, help="If 1, scan GridFF Fz vs z_start for current xy to locate attractive region (Fz<0).")
    parser.add_argument("--scan_z0",            type=float, default=0.5, help="scan_z: starting z_start (Angstrom)")
    parser.add_argument("--scan_z1",            type=float, default=6.0, help="scan_z: ending z_start (Angstrom)")
    parser.add_argument("--scan_nz",            type=int,   default=23, help="scan_z: number of samples")
    parser.add_argument("--use_tip_constraint", type=int,   default=1, help="If 1, hard-constrain anchor atom to tip trajectory using constr/constrK (recommended for debugging rotation).")
    parser.add_argument("--K_constraint",       type=float, default=2000.0, help="Constraint stiffness for hard tip constraint (eV/A^2-ish units in this integrator).")
    parser.add_argument("--n_steps_path",       type=int,   default=41, help="Number of replicas along manipulation path. Debugging default: 41 (matches reference runs). Use 1 for single-point relaxation debug.")
    parser.add_argument("--dx",                 type=float, default=0.1, help="Tip step in x per replica (A). Debugging default: 0.1 (matches sequential relaxed scan).")
    parser.add_argument("--n_steps_md",         type=int,   default=20000, help="Number of MD integration steps per relaxation run. Debugging default: 20000.")
    parser.add_argument("--dt",                 type=float, default=0.02, help="Integrator timestep")
    parser.add_argument("--damp",               type=float, default=0.01, help="Damping fraction per step (same semantics as test_relaxed_scan_tip.py): vfac=1-damp")
    parser.add_argument("--K_band_sched",       type=str,   default="0.5,0.2,0.05,0", help="Comma-separated K_band schedule for causal tethers. Use '0' to disable. Keep K<1 to avoid dissociating bonds under large inter-replica gaps.")
    parser.add_argument("--L_allowed",          type=float, default=0.2, help="Allowed per-atom gap before tethers pull (A). Set 0 with K=0 to fully disable.")
    parser.add_argument("--zero_vel",           type=int,   default=1, help="If 1, set initial velocities to zero before relaxation (recommended).")
    parser.add_argument("--shift_init",         type=int,   default=1, help="If 1, rigidly shift initial geometry for each replica to match tip position (avoids initial bond strain).")
    parser.add_argument("--trace_interval",     type=int,   default=5, help="Store trace frame every N MD steps")
    parser.add_argument("--out_dir",            type=str,   default="out_pathopt", help="Directory for output files")
    parser.add_argument("--init_seq_npz",       type=str,   default="", help="If set, initialize replicas from sequential relaxed scan npz (expects A_scan with shape [nscan+1,natoms,3]). This makes band-off PathOpt directly comparable to sequential scan (warm-start).")
    return parser.parse_args()

def _print_z_stats(z, label):
    z = np.asarray(z, dtype=float)
    print(f"{label}: z0={z[0]:.6f} z1={z[-1]:.6f} zmin={z.min():.6f} zmax={z.max():.6f} dz={z[-1]-z[0]:.6f}")

def _plot_xz_traj(A0, enames, out_path, title="x-z trajectories"):
    # A0: (nFrames, natoms, 3)
    nF, nA, _ = A0.shape
    plt.figure(figsize=(6, 6))
    for ia in range(nA):
        x = A0[:, ia, 0]
        z = A0[:, ia, 2]
        name = enames[ia] if (enames is not None and ia < len(enames)) else f"a{ia}"
        plt.plot(x, z, '-', label=f"{ia}:{name}")
        plt.plot([x[0]], [z[0]], 'o', ms=4)
        plt.plot([x[-1]], [z[-1]], 's', ms=4)
    plt.xlabel('x (A)')
    plt.ylabel('z (A)')
    plt.title(title)
    plt.grid(True)
    plt.axis('equal')
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()

if __name__ == "__main__":
    args = parse_args()
    np.set_printoptions(precision=4, suppress=True)
    
    os.makedirs(args.out_dir, exist_ok=True)
    
    if args.mol:
        mol = AtomicSystem(args.mol)
        # Assuming ngs, atypes etc are properly set or need setting
        # For a robust test, user might need to ensure mol is ready for MMFF
    else:
        mol = create_h2o()

    if args.z_start != 0.0:
        mol.apos[:, 2] += np.float32(args.z_start)

    if (args.shift_xy is not None) and (float(args.shift_xy[0]) != 0.0 or float(args.shift_xy[1]) != 0.0):
        mol.apos[:, 0] += np.float32(args.shift_xy[0])
        mol.apos[:, 1] += np.float32(args.shift_xy[1])

    tip_handle_idx = args.anchor 
    plot_atom_idx = args.plot_atom
    
    n_pop = 1
    n_steps = int(args.n_steps_path)
    
    # Tip trajectory: move across x
    tip_traj = np.zeros((n_pop, n_steps, 3), dtype=np.float32)
    # Match sequential relaxed scan sampling: tip.x = x0 + dx*i
    for i in range(n_steps):
        tip_traj[0, i, 0] = mol.apos[tip_handle_idx, 0] + np.float32(float(args.dx) * i)
        tip_traj[0, i, 1] = mol.apos[tip_handle_idx, 1]
        tip_traj[0, i, 2] = mol.apos[tip_handle_idx, 2] + float(args.tip_dz)
        
    print(f"Initializing ManipulationPathOpt with {n_pop} trajectories, {n_steps} steps...")
    grid_ok = False
    grid_path_in = args.grid
    if grid_path_in:
        grid_path = grid_path_in
        if not os.path.isabs(grid_path):
            grid_path = os.path.abspath(os.path.join(os.getcwd(), grid_path))
        if os.path.exists(grid_path):
            grid_ok = True
        else:
            raise FileNotFoundError(f"Grid file not found: {grid_path}")

    opt = ManipulationPathOpt(mol, tip_handle_idx=tip_handle_idx, n_pop=n_pop, n_steps=n_steps, enable_nonbond=bool(grid_ok))

    if grid_ok:
        grid_path = grid_path
        print(f"Loading GridFF Bspline from: {grid_path}")
        V = np.load(grid_path)
        if V.ndim != 4:
            raise ValueError(f"Expected Bspline_PLQd.npy to be 4D [nx,ny,nz,nch]; got {V.shape}")
        if V.shape[3] == 3:
            V4 = np.zeros((V.shape[0], V.shape[1], V.shape[2], 4), dtype=np.float32)
            V4[:, :, :, :3] = V.astype(np.float32)
        elif V.shape[3] == 4:
            V4 = V.astype(np.float32)
        else:
            raise ValueError(f"Expected Bspline_PLQd.npy channels 3 or 4; got {V.shape}")

        opt.set_substrate_gridff(
            V4,
            grid_p0=tuple(args.grid_p0),
            grid_step=tuple(args.grid_step),
            use_texture=bool(args.grid_tex),
            alpha_morse=float(args.grid_alpha),
        )

    def eval_gridff_fz(apos):
        opt.init_replica_states(apos, tip_traj)
        opt.md.run_cleanForceMMFFf4()
        if bool(args.grid_tex):
            opt.md.run_getNonBond_GridFF_Bspline_tex()
        else:
            # Use the same exclusion-aware kernel as the sequential relaxed scan when available
            if getattr(opt.md, 'kernel_args_getNonBond_GridFF_Bspline_ex2', None) is not None:
                opt.md.run_getNonBond_GridFF_Bspline_ex2()
            else:
                opt.md.run_getNonBond_GridFF_Bspline()
        f0 = opt.download_forces()[0]
        return f0[:, 2].copy()
    
    print(f"Setting tip Morse params: De={args.De}, r0={args.r0}, alpha={args.alpha}")
    opt.set_morse_params(args.De, args.r0, args.alpha)
    
    print("Setting up tip constraints...")
    opt.set_tip_trajectory(tip_traj)
    if int(args.use_tip_constraint) != 0:
        print(f"Enabling hard tip positional constraint: K={args.K_constraint}")
        opt.set_tip_pos_constraints(K=float(args.K_constraint))
    
    print("Initializing replica states (GPU interpolation)...")
    opt.init_replica_states(mol.apos, tip_traj, bShift=(int(args.shift_init) != 0))

    bWarmStart = (str(args.init_seq_npz).strip() != "")
    if str(args.init_seq_npz).strip() != "":
        D = np.load(args.init_seq_npz)
        if 'A_scan' not in D:
            raise RuntimeError(f"--init_seq_npz={args.init_seq_npz} does not contain key 'A_scan'")
        A_scan = np.asarray(D['A_scan'], dtype=np.float32)
        if A_scan.ndim != 3 or A_scan.shape[1] != mol.natoms or A_scan.shape[2] != 3:
            raise RuntimeError(f"A_scan has unexpected shape {A_scan.shape}; expected (nFrames,{mol.natoms},3)")
        if A_scan.shape[0] < n_steps:
            raise RuntimeError(f"A_scan has only {A_scan.shape[0]} frames, but n_steps_path={n_steps} replicas requested")

        # Fill per-replica positions from sequential scan frames (warm-start)
        apos_all = np.zeros((opt.n_rep, opt.mm.nvecs, 4), dtype=np.float32)
        for i in range(n_steps):
            apos_all[i, :mol.natoms, :3] = A_scan[i, :, :]
        opt.md.toGPU('apos', apos_all)

        # If hard tip constraints are enabled, refresh them (constr depends on tip positions but not on apos; still safe)
        if int(args.use_tip_constraint) != 0:
            opt.set_tip_pos_constraints(K=float(args.K_constraint))

    if int(args.zero_vel) != 0:
        v0 = np.zeros((opt.n_rep, opt.mm.nvecs, 4), dtype=np.float32)
        opt.md.toGPU('avel', v0)

    if args.grid and (not bWarmStart):
        fz = eval_gridff_fz(mol.apos)
        print(f"Initial GridFF Fz (replica0) O,H1,H2 = {fz}")

        if int(args.scan_z) != 0:
            z0 = float(args.scan_z0)
            z1 = float(args.scan_z1)
            nz = int(args.scan_nz)
            zs = np.linspace(z0, z1, nz)
            print("scan_z: z_start, Fz(O),Fz(H1),Fz(H2)")
            apos0 = mol.apos.copy()
            zref = float(apos0[0, 2] - float(args.z_start))
            for z in zs:
                apos = apos0.copy()
                apos[:, 2] = np.float32(zref + z)
                fz = eval_gridff_fz(apos)
                print(f"{z:8.3f}  {fz[0]:+12.6f} {fz[1]:+12.6f} {fz[2]:+12.6f}")

    if bWarmStart and (int(args.scan_z) != 0 or int(args.auto_xy) != 0):
        raise RuntimeError("--init_seq_npz is incompatible with --scan_z/--auto_xy because those diagnostics reinitialize replicas and overwrite the warm-start geometry")

    if int(args.auto_xy) != 0:
        # Pick best dx,dy for attraction of H2 (atom index 2) to the substrate.
        apos0 = mol.apos.copy()
        span = float(args.auto_xy_span)
        nxy = int(args.auto_xy_n)
        dxs = np.linspace(-span, span, nxy, dtype=np.float32)
        dys = np.linspace(-span, span, nxy, dtype=np.float32)
        best = None
        best_xy = (0.0, 0.0)
        for dx in dxs:
            for dy in dys:
                apos = apos0.copy()
                apos[:, 0] += dx
                apos[:, 1] += dy
                opt.init_replica_states(apos, tip_traj)
                # Only GridFF force (no MMFF), evaluate once
                opt.md.run_cleanForceMMFFf4()
                if bool(args.grid_tex):
                    opt.md.run_getNonBond_GridFF_Bspline_tex()
                else:
                    opt.md.run_getNonBond_GridFF_Bspline()
                f = opt.download_forces()[0]  # (natoms,4)
                fz = float(f[2, 2])
                if (best is None) or (fz < best):
                    best = fz
                    best_xy = (float(dx), float(dy))
        print(f"auto_xy: best shift dx,dy={best_xy} gives Fz(H2)={best:+.6f}")
        mol.apos = apos0
        mol.apos[:, 0] += np.float32(best_xy[0])
        mol.apos[:, 1] += np.float32(best_xy[1])
        opt.init_replica_states(mol.apos, tip_traj)

        if args.grid:
            opt.md.run_cleanForceMMFFf4()
            if bool(args.grid_tex):
                opt.md.run_getNonBond_GridFF_Bspline_tex()
            else:
                opt.md.run_getNonBond_GridFF_Bspline()
            f0 = opt.download_forces()[0]
            fz = f0[:, 2]
            print(f"After auto_xy GridFF Fz (replica0) O,H1,H2 = {fz}")
    
    # Save initial unrelaxed interpolation
    apos_initial = opt.download_states()
    opt.save_xyz_movie(apos_initial, os.path.join(args.out_dir, "traj_initial.xyz"))
    
    print("Running relaxation with annealing...")
    K_band_sched = [float(x) for x in str(args.K_band_sched).split(',') if len(x.strip()) > 0]
    L_allowed = float(args.L_allowed)
    
    history_apos, history_gaps, trace_apos = opt.run_relax(
        n_steps_md=int(args.n_steps_md), 
        K_band_sched=K_band_sched, 
        L_allowed=L_allowed,
        dt=float(args.dt),
        damp=float(args.damp),
        out_dir=args.out_dir,
        bTrace=True,
        trace_interval=int(args.trace_interval)
    )

    # ---- Diagnostics: z(t) for both hydrogens in the single-point relaxation trace
    if len(trace_apos) > 0:
        A = np.array(trace_apos, dtype=np.float32)           # (nFrames, n_rep, natoms, 3)
        A0 = A[:, 0, :, :]                                  # (nFrames, natoms, 3)
        zH1 = A0[:, 1, 2]
        zH2 = A0[:, 2, 2]
        _print_z_stats(zH1, "H1")
        _print_z_stats(zH2, "H2")
        t = np.arange(len(zH1)) * float(args.trace_interval) * float(args.dt)
        plt.figure(figsize=(8, 5))
        plt.plot(t, zH1, '-', label='H1 z(t)')
        plt.plot(t, zH2, '-', label='H2 z(t)')
        plt.xlabel('time')
        plt.ylabel('z (A)')
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(args.out_dir, 'hydrogens_z_trace.png'))
        plt.close()

        _plot_xz_traj(A0, getattr(mol, 'enames', None), os.path.join(args.out_dir, 'atoms_xz_trace.png'), title='Relaxation trace x-z (replica 0)')
    
    # Save a movie of the relaxation of the first replica (pop 0, step 0) 
    # to see how it "hangs" and relaxes on the surface.
    trace_xyz_path = os.path.join(args.out_dir, "relax_trace_p0_s0.xyz")
    print(f"Saving relaxation trace movie to {trace_xyz_path} ...")
    with open(trace_xyz_path, "w") as f:
        for i, apos_t in enumerate(trace_apos):
            f.write(f"{mol.natoms}\n")
            f.write(f"Trace frame {i} pop 0 step 0\n")
            # apos_t is (n_rep, natoms, 3)
            for ia in range(mol.natoms):
                name = mol.enames[ia] if hasattr(mol, 'enames') else "X"
                r = apos_t[0, ia, :] # just taking replica 0
                f.write(f"{name} {r[0]:.5f} {r[1]:.5f} {r[2]:.5f}\n")

    print("Plotting trajectories and penalties...")
    opt.plot_trajectories(history_apos, history_gaps, K_band_sched, plot_atom_idx=plot_atom_idx, out_dir=args.out_dir)
    
    print(f"Done. Outputs saved to {args.out_dir}/")

