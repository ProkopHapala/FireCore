#!/usr/bin/env python3

import os, sys, argparse, math, time
import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(BASE)

from pyBall.AtomicSystem import AtomicSystem
from pyBall.OCL.UFF import UFF_CL
from pyBall.MD_test_utils import clean_element_symbols_list, apply_random_displacement, max_norm3, write_xyz_frame, iter_should_write, conv_plot_path_from_xyz, plot_fmax_history_png, format_relax_summary, effective_verbosity

np.set_printoptions(linewidth=np.inf, suppress=False)


def main():
    ap = argparse.ArgumentParser(description='PyOpenCL UFF relaxation (forces on GPU, integration on CPU) with safety checks and XYZ trajectory output')
    ap.add_argument('--molecule', required=True, type=str)
    ap.add_argument('--dt', type=float, default=1e-3)
    ap.add_argument('--damp', type=float, default=0.05, help='damping in [0,1]; velocity factor = 1-damp')
    ap.add_argument('--dmax', type=float, default=0.2, help='random initial displacement in [-dmax,dmax] for each coordinate; set 0 to disable')
    ap.add_argument('--seed', type=int, default=0, help='RNG seed for initial displacement; 0 means deterministic default generator seed')
    ap.add_argument('--steps', type=int, default=500)
    ap.add_argument('--stride', type=int, default=1, help='(deprecated) write XYZ every N steps; use --write-stride')
    ap.add_argument('--write-stride', type=int, default=1, help='write XYZ every N steps (also writes step 0 and final)')
    ap.add_argument('--fconv', type=float, default=1e-3, help='stop if max|F| (per-atom norm) < fconv')
    ap.add_argument('--fmax-abort', type=float, default=1e+3, help='abort if max|F| exceeds this')
    ap.add_argument('--rmax-abort', type=float, default=1e+3, help='abort if max|r| exceeds this (from origin)')
    ap.add_argument('--out', type=str, default='', help='output XYZ trajectory path')
    ap.add_argument('--no-traj', type=int, default=0)
    ap.add_argument('--nloc', type=int, default=32)
    ap.add_argument('--quiet', type=int, default=0)
    ap.add_argument('--verbosity', type=int, default=1, help='0=silent, 1=progress, 2=verbose progress')
    ap.add_argument('--fast-exit', type=int, default=0)
    args = ap.parse_args()

    assert args.steps > 0
    assert args.dt > 0
    assert 0.0 <= args.damp < 1.0

    # Disable noisy AtomicSystem debug unless user explicitly enables it
    os.environ.setdefault('PYBALL_DEBUG_ATOMICSYSTEM', '0')

    mol = AtomicSystem(fname=args.molecule)
    if mol.enames is None or mol.apos is None:
        raise RuntimeError('Failed to load molecule')

    enames = clean_element_symbols_list(mol.enames)

    nat = len(mol.apos)
    pos = np.array(mol.apos, dtype=np.float32, copy=True)
    vel = np.zeros_like(pos, dtype=np.float32)

    apply_random_displacement(pos, args.dmax, seed=(None if args.seed == 0 else int(args.seed)))

    out_path = args.out if args.out else (os.path.splitext(os.path.basename(args.molecule))[0] + '_UFF_pyocl_relax.xyz')

    # Compile UFF kernels with debug prints off (PyOpenCL build options)
    uff = UFF_CL(nloc=args.nloc, bPrint=False, debug_build_options='-DDBG_UFF=0')
    uff_data = uff.toUFF(mol, bRealloc=True)
    uff.upload_positions(pos)
    uff.upload_topology_params(uff_data)

    # Ensure bonded-only for now
    uff.bDoNonBonded = False

    vfac = np.float32(1.0 - args.damp)
    dt = np.float32(args.dt)

    write_stride = int(args.write_stride if args.write_stride is not None else args.stride)
    if write_stride <= 0:
        write_stride = int(args.stride)
    if write_stride <= 0:
        write_stride = 1
    verbosity = effective_verbosity(args.quiet, args.verbosity)

    t0 = time.time()
    n_written = 0

    fout = None
    if not args.no_traj:
        fout = open(out_path, 'w')
        write_xyz_frame(fout, enames, pos, comment=f"# step 0")
        n_written += 1

    steps_hist = []
    fmax_hist  = []
    converged = False
    abort_reason = None
    last_fmax = 0.0
    last_step = 0

    try:
        for istep in range(1, args.steps + 1):
            uff.run_eval_step(bClearForce=True)
            F = np.array(uff.get_forces(iSys=0), dtype=np.float32, copy=True)

            if not np.isfinite(F).all():
                raise RuntimeError(f"NaN/Inf in forces at step {istep}")

            fmax = max_norm3(F)
            steps_hist.append(int(istep))
            fmax_hist.append(float(fmax))
            last_fmax = float(fmax)
            last_step = int(istep)
            if fmax > args.fmax_abort:
                raise RuntimeError(f"Blow-up: max|F|={fmax:.6e} > fmax-abort={args.fmax_abort:.6e} at step {istep}")

            # Integrate (simple damped explicit Euler / velocity form)
            vel[:] = (vel + dt * F) * vfac
            pos[:] = pos + dt * vel

            if not np.isfinite(pos).all():
                raise RuntimeError(f"NaN/Inf in positions at step {istep}")

            rmax = max_norm3(pos)
            if rmax > args.rmax_abort:
                raise RuntimeError(f"Blow-up: max|r|={rmax:.6e} > rmax-abort={args.rmax_abort:.6e} at step {istep}")

            # Upload updated positions for next force evaluation
            uff.upload_positions(pos)
            uff.queue.finish()

            if (fout is not None) and iter_should_write(istep, write_stride, args.steps):
                write_xyz_frame(fout, enames, pos, comment=f"# step {istep}  maxF {fmax:.6e}  maxR {rmax:.6e}")
                n_written += 1

            if verbosity > 0:
                k = 1 if verbosity > 1 else 10
                if (istep % max(1, k * write_stride) == 0) or (istep == 1):
                    dt_wall = time.time() - t0
                    print(f"UFF step {istep:6d}/{args.steps}  max|F|={fmax:.3e}  max|r|={rmax:.3e}  wall={dt_wall:.2f}s")

            dt_wall = time.time() - t0
            if fmax < args.fconv:
                converged = True
                if verbosity > 0:
                    print(f"UFF converged at step {istep}: max|F|={fmax:.6e} < fconv={args.fconv:.6e}")
                break

    finally:
        if fout is not None:
            fout.close()

    out_png = conv_plot_path_from_xyz(out_path)
    plot_fmax_history_png(out_png, steps_hist, fmax_hist, fconv=args.fconv, title=f"UFF {os.path.basename(args.molecule)}")

    dt_wall = time.time() - t0
    if verbosity > 0:
        print(f"UFF done: wrote {n_written} frames to {out_path}  wall={dt_wall:.2f}s")

    print(format_relax_summary('UFF', args.molecule, converged, last_step, last_fmax, out_path, out_png, dt_wall, abort_reason=abort_reason))

    if args.fast_exit:
        sys.stdout.flush(); sys.stderr.flush(); os._exit(0)

    return 0


if __name__ == '__main__':
    raise SystemExit(main())
