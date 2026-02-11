#!/usr/bin/env python3

import os, sys, argparse, time
import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(BASE)

from pyBall.AtomicSystem import AtomicSystem
from pyBall.OCL.MMFF import MMFF as MMFF_pyocl
from pyBall.OCL.MolecularDynamics import MolecularDynamics
from pyBall.MD_test_utils import clean_element_symbols_list, apply_random_displacement, max_norm3, write_xyz_frame, iter_should_write, conv_plot_path_from_xyz, plot_fmax_history_png, format_relax_summary, effective_verbosity

np.set_printoptions(linewidth=np.inf, suppress=False)


def main():
    ap = argparse.ArgumentParser(description='PyOpenCL MMFF relaxation (forces+integration on GPU) with safety checks and XYZ trajectory output')
    ap.add_argument('--molecule', required=True, type=str)
    ap.add_argument('--dt', type=float, default=1e-3)
    ap.add_argument('--damp', type=float, default=0.05, help='damping in [0,1]; velocity factor = 1-damp')
    ap.add_argument('--dmax', type=float, default=0.2, help='random initial displacement in [-dmax,dmax] for each coordinate; set 0 to disable')
    ap.add_argument('--seed', type=int, default=0, help='RNG seed for initial displacement; 0 means deterministic default generator seed')
    ap.add_argument('--steps', type=int, default=500)
    ap.add_argument('--stride', type=int, default=1, help='(deprecated) write XYZ every N steps; use --write-stride')
    ap.add_argument('--write-stride', type=int, default=1, help='write XYZ every N steps (also writes step 0 and final)')
    ap.add_argument('--fconv', type=float, default=1e-3)
    ap.add_argument('--fmax-abort', type=float, default=1e+3)
    ap.add_argument('--rmax-abort', type=float, default=1e+3)
    ap.add_argument('--out', type=str, default='')
    ap.add_argument('--no-traj', type=int, default=0)
    ap.add_argument('--nloc', type=int, default=32)
    ap.add_argument('--quiet', type=int, default=0)
    ap.add_argument('--verbosity', type=int, default=1, help='0=silent, 1=progress, 2=verbose progress')
    ap.add_argument('--nonbond', type=int, default=0, help='keep off by default; 1 enables nonbond kernel path')
    ap.add_argument('--angles', type=int, default=1)
    ap.add_argument('--pisigma', type=int, default=1)
    ap.add_argument('--pipii', type=int, default=1)
    ap.add_argument('--force-node-all', type=int, default=1)
    ap.add_argument('--fast-exit', type=int, default=0)
    args = ap.parse_args()

    assert args.steps > 0
    assert args.dt > 0
    assert 0.0 <= args.damp < 1.0

    os.environ.setdefault('PYBALL_DEBUG_ATOMICSYSTEM', '0')

    mol = AtomicSystem(fname=args.molecule)
    if getattr(mol, 'ngs', None) is None:
        mol.neighs()

    enames = clean_element_symbols_list(mol.enames)
    apply_random_displacement(mol.apos, args.dmax, seed=(None if args.seed == 0 else int(args.seed)))

    out_path = args.out if args.out else (os.path.splitext(os.path.basename(args.molecule))[0] + '_MMFF_pyocl_relax.xyz')

    write_stride = int(args.write_stride if args.write_stride is not None else args.stride)
    if write_stride <= 0:
        write_stride = int(args.stride)
    if write_stride <= 0:
        write_stride = 1
    verbosity = effective_verbosity(args.quiet, args.verbosity)

    mm = MMFF_pyocl(bTorsion=False, verbosity=0)
    if args.force_node_all:
        mm.capping_atoms = set()
        mm.reorder_nodes_first = False

    mm.toMMFFsp3_loc(mol=mol, atom_types=mm.atom_types, bRealloc=True, bEPairs=False, bUFF=False)
    mm.make_back_neighs(b_cap_neighs=not args.force_node_all)

    # Build relax_multi.cl with debug prints disabled to avoid huge slowdowns
    md = MolecularDynamics(nloc=args.nloc, debug_build_options='-DDBG_UFF=0', enable_nonbond=bool(args.nonbond))
    md.realloc(mm, nSystems=1)
    md.setup_kernels()
    md.set_pack_system_debug(False)
    md.pack_system(0, mm)

    # Force-only / relaxation params (dt, Flim, vel_damp_factor, 0)
    dt = np.float32(args.dt)
    vfac = np.float32(1.0 - args.damp)
    Flim = np.float32(1e+6)
    md.toGPU('MDparams', np.array([dt, Flim, vfac, 0.0], dtype=np.float32), byte_offset=0)

    # Nonbond off by default
    md.kernel_params['bSubtractVdW'] = np.int32(1 if args.nonbond else 0)

    t0 = time.time()
    n_written = 0

    fout = None
    if not args.no_traj:
        fout = open(out_path, 'w')

    import pyopencl as cl
    apos_buf = np.empty((mm.nvecs, 4), dtype=np.float32)
    af_buf   = np.empty((mm.nvecs, 4), dtype=np.float32)
    def download_aforce_apos():
        cl.enqueue_copy(md.queue, apos_buf, md.buffer_dict['apos'])
        cl.enqueue_copy(md.queue, af_buf,   md.buffer_dict['aforce'])
        md.queue.finish()
        return apos_buf[:mm.natoms, :3].copy(), af_buf[:mm.natoms, :3].copy()

    steps_hist = []
    fmax_hist  = []
    converged = False
    abort_reason = None
    last_fmax = 0.0
    last_step = 0

    try:
        pos0, F0 = download_aforce_apos()
        if fout is not None:
            write_xyz_frame(fout, enames, pos0, comment="# step 0")
            n_written += 1

        for istep in range(1, args.steps + 1):
            # bonded-only path: clean -> getMMFF -> updateAtoms
            md.run_cleanForceMMFFf4()
            md.run_getMMFFf4()
            md.run_updateAtomsMMFFf4()

            pos, F = download_aforce_apos()
            if not np.isfinite(pos).all():
                raise RuntimeError(f"NaN/Inf in positions at step {istep}")
            if not np.isfinite(F).all():
                raise RuntimeError(f"NaN/Inf in forces at step {istep}")

            fmax = max_norm3(F)
            rmax = max_norm3(pos)
            steps_hist.append(int(istep))
            fmax_hist.append(float(fmax))
            last_fmax = float(fmax)
            last_step = int(istep)

            if fmax > args.fmax_abort:
                raise RuntimeError(f"Blow-up: max|F|={fmax:.6e} > fmax-abort={args.fmax_abort:.6e} at step {istep}")
            if rmax > args.rmax_abort:
                raise RuntimeError(f"Blow-up: max|r|={rmax:.6e} > rmax-abort={args.rmax_abort:.6e} at step {istep}")

            if fout is not None and iter_should_write(istep, write_stride, args.steps):
                write_xyz_frame(fout, enames, pos, comment=f"# step {istep}  maxF {fmax:.6e}  maxR {rmax:.6e}")
                n_written += 1

            if verbosity > 0:
                k = 1 if verbosity > 1 else 10
                if (istep % max(1, (k * write_stride)) == 0) or (istep == 1):
                    dt_wall = time.time() - t0
                    print(f"MMFF step {istep:6d}/{args.steps}  max|F|={fmax:.3e}  max|r|={rmax:.3e}  wall={dt_wall:.2f}s")

            if fmax < args.fconv:
                converged = True
                if verbosity > 0:
                    print(f"MMFF converged at step {istep}: max|F|={fmax:.6e} < fconv={args.fconv:.6e}")
                break

    finally:
        if fout is not None:
            fout.close()

    out_png = conv_plot_path_from_xyz(out_path)
    plot_fmax_history_png(out_png, steps_hist, fmax_hist, fconv=args.fconv, title=f"MMFF {os.path.basename(args.molecule)}")

    dt_wall = time.time() - t0
    if verbosity > 0:
        print(f"MMFF done: wrote {n_written} frames to {out_path}  wall={dt_wall:.2f}s")

    print(format_relax_summary('MMFF', args.molecule, converged, last_step, last_fmax, out_path, out_png, dt_wall, abort_reason=abort_reason))

    if args.fast_exit:
        sys.stdout.flush(); sys.stderr.flush(); os._exit(0)

    return 0


if __name__ == '__main__':
    raise SystemExit(main())
