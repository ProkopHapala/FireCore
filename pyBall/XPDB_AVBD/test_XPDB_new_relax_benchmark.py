import argparse
import os
import sys

import numpy as np
import pyopencl as cl

# Silence PyOpenCL compiler output unless explicitly enabled
os.environ.setdefault("PYOPENCL_COMPILER_OUTPUT", "0")

# Ensure local imports
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from XPDB_new import XPDB_new
import XPTB_utils as xu


def backup_device_positions(sim):
    import pyopencl as cl
    mf = cl.mem_flags
    cl_pos_backup = cl.Buffer(sim.ctx, mf.READ_WRITE, sim.num_atoms * 16)
    cl_prev_backup = cl.Buffer(sim.ctx, mf.READ_WRITE, sim.num_atoms * 16)
    cl.enqueue_copy(sim.queue, cl_pos_backup, sim.cl_pos).wait()
    cl.enqueue_copy(sim.queue, cl_prev_backup, sim.cl_prev_pos).wait()
    return cl_pos_backup, cl_prev_backup


def restore_device_positions(sim, cl_pos_backup, cl_prev_backup):
    import pyopencl as cl
    cl.enqueue_copy(sim.queue, sim.cl_pos, cl_pos_backup).wait()
    cl.enqueue_copy(sim.queue, sim.cl_prev_pos, cl_prev_backup).wait()


def prepare_h2o_sim(args):
    pos_ref, bonds_adj = xu.make_h2o_geometry(add_angle=not args.no_angle)

    d = np.array([float(x) for x in args.direction.split(",")], dtype=np.float32)
    if d.size != 3:
        raise ValueError("--direction must have 3 comma-separated components")

    if args.deform == "shift":
        pos_def = xu.deform_shift_atom(pos_ref, atom_idx=args.atom, direction=d, shift=args.shift)
    elif args.deform == "scale":
        pos_def = xu.deform_scale_along_direction(pos_ref, direction=d, scale=args.scale)
    else:
        raise ValueError(f"unknown deform {args.deform}")

    n = len(pos_def)
    vel0 = np.zeros_like(pos_def, dtype=np.float32)
    radius = np.full(n, args.radius, dtype=np.float32)
    mass = np.full(n, args.mass, dtype=np.float32)

    sim = XPDB_new(n, group_size=args.group_size, prefer_gpu=bool(args.prefer_gpu), device_idx=args.device_idx)
    sim.upload_data(pos_def, vel0, radius, mass)
    sim.upload_bonds_fixed(bonds_adj)

    # Build topology once for molecule-only case (no collisions) and rebuild each iter if collisions enabled
    Rmax = max(float(args.radius), float(xu.bonds_to_max_L0(bonds_adj, default=1.0)))
    sim.build_local_topology(Rmax=Rmax, coll_scale=args.coll_scale, bbox_scale=args.bbox_scale)
    sim.reset_prev_pos()

    return sim, bonds_adj


def run_single(sim, args):
    xu.ensure_outdir(args.outdir)

    xu.print_run_header(args.tag, {
        "mode": "single",
        "deform": args.deform,
        "dt": args.dt,
        "omega": args.omega,
        "beta": args.momentum_beta,
        "k_coll": args.k_coll,
        "iters": args.iters,
        "tol": args.stop_tol,
    })

    iters = int(args.iters)
    traj = np.empty((iters + 1, sim.num_atoms, 3), dtype=np.float32) if args.save_traj else None
    res = []

    if traj is not None:
        traj[0] = sim.get_positions()

    for i in range(iters):
        if not args.use_inertia:
            cl.enqueue_copy(sim.queue, sim.cl_pred_pos, sim.cl_pos).wait()
        if args.k_coll != 0.0:
            # Collisions depend on positions; rebuild local topology each outer step.
            # (Bond reindexing is cheap for small molecules.)
            # Keep Rmax the same.
            sim.build_local_topology(Rmax=args.Rmax, coll_scale=args.coll_scale, bbox_scale=args.bbox_scale)

        if args.solver == "global":
            sim.solve_cluster_jacobi(dt=args.dt, iterations=1, k_coll=args.k_coll, omega=args.omega, momentum_beta=args.momentum_beta)
        elif args.solver == "local":
            sim.solve_cluster_jacobi_local(dt=args.dt, inner_iters=1, k_coll=args.k_coll, omega=args.omega, momentum_beta=args.momentum_beta)
        elif args.solver == "local_nocoll":
            sim.solve_cluster_jacobi_local_nocoll(dt=args.dt, inner_iters=1, omega=args.omega, momentum_beta=args.momentum_beta)
        else:
            raise ValueError(f"unknown solver {args.solver}")
        r = sim.bond_residual_norms()
        res.append((i + 1, r["linf"], r["l2"]))

        if traj is not None:
            traj[i + 1] = sim.get_positions(out=None)

        if args.stop_tol > 0.0 and r["linf"] < args.stop_tol:
            print(f"[STOP] reached tol at iter={i+1} linf={r['linf']:.3e} l2={r['l2']:.3e}")
            if traj is not None and i + 1 < iters:
                traj[i + 2 :] = traj[i + 1]
            break

        if (i % args.print_every) == 0 or i == iters - 1:
            print(f"[ITER {i+1:4d}] linf={r['linf']:.6e} l2={r['l2']:.6e}")

    res = np.array(res, dtype=np.float32)

    np.savetxt(os.path.join(args.outdir, args.tag + "_residuals.txt"), res, header="iter linf l2")
    if traj is not None:
        np.save(os.path.join(args.outdir, args.tag + "_traj.npy"), traj)

    if args.plot and (not args.noshow):
        xu.plot_residual_series(res, title=args.tag, noshow=args.noshow, outpath=None)

    return res


def run_sweep(sim, args):
    xu.ensure_outdir(args.outdir)

    dt_vals = np.array([float(x) for x in args.dt_list.split(",") if x.strip() != ""], dtype=np.float32)
    omega_vals = np.array([float(x) for x in args.omega_list.split(",") if x.strip() != ""], dtype=np.float32)
    beta_vals = np.array([float(x) for x in args.beta_list.split(",") if x.strip() != ""], dtype=np.float32)
    if dt_vals.size == 0 or omega_vals.size == 0 or beta_vals.size == 0:
        raise ValueError("dt_list, omega_list, beta_list must be non-empty")

    it_req = np.full((len(dt_vals), len(omega_vals), len(beta_vals)), int(args.iters) + 1, dtype=np.int32)
    linf_fin = np.full_like(it_req, np.nan, dtype=np.float32)

    pos_bak, prev_bak = backup_device_positions(sim)

    for idt, dt in enumerate(dt_vals):
        for io, omega in enumerate(omega_vals):
            for ib, beta in enumerate(beta_vals):
                restore_device_positions(sim, pos_bak, prev_bak)

                xu.print_run_header(args.tag, {
                    "mode": "sweep",
                    "dt": float(dt),
                    "omega": float(omega),
                    "beta": float(beta),
                    "k_coll": float(args.k_coll),
                })

                last_linf = None
                for i in range(int(args.iters)):
                    if not args.use_inertia:
                        cl.enqueue_copy(sim.queue, sim.cl_pred_pos, sim.cl_pos).wait()
                    if args.k_coll != 0.0:
                        sim.build_local_topology(Rmax=args.Rmax, coll_scale=args.coll_scale, bbox_scale=args.bbox_scale)

                    if args.solver == "global":
                        sim.solve_cluster_jacobi(dt=float(dt), iterations=1, k_coll=float(args.k_coll), omega=float(omega), momentum_beta=float(beta))
                    elif args.solver == "local":
                        sim.solve_cluster_jacobi_local(dt=float(dt), inner_iters=1, k_coll=float(args.k_coll), omega=float(omega), momentum_beta=float(beta))
                    elif args.solver == "local_nocoll":
                        sim.solve_cluster_jacobi_local_nocoll(dt=float(dt), inner_iters=1, omega=float(omega), momentum_beta=float(beta))
                    else:
                        raise ValueError(f"unknown solver {args.solver}")
                    r = sim.bond_residual_norms()
                    last_linf = float(r["linf"])
                    if args.stop_tol > 0.0 and r["linf"] < args.stop_tol:
                        it_req[idt, io, ib] = i + 1
                        linf_fin[idt, io, ib] = r["linf"]
                        break
                linf_fin[idt, io, ib] = float(last_linf) if last_linf is not None else np.nan

    np.save(os.path.join(args.outdir, args.tag + "_itreq.npy"), it_req)
    np.save(os.path.join(args.outdir, args.tag + "_linf_final.npy"), linf_fin)

    # Print best
    best = int(np.min(it_req))
    idx = np.unravel_index(np.argmin(it_req), it_req.shape)
    print(f"[BEST] iters={best} idx(dt,omega,beta)={idx} dt={dt_vals[idx[0]]} omega={omega_vals[idx[1]]} beta={beta_vals[idx[2]]} linf={linf_fin[idx]:.3e}")


def main():
    parser = argparse.ArgumentParser(description="XPDB_new tiled Jacobi relaxation benchmark (fixed-slot bonds)")
    parser.add_argument("--mode", type=str, default="single", choices=["single", "sweep"])

    parser.add_argument("--solver", type=str, default="global", choices=["global", "local", "local_nocoll"], help="Solver variant: global=kernel-per-iter (default), local=in-kernel loop (isolated systems only), local_nocoll=in-kernel loop without collisions")

    parser.add_argument("--tag", type=str, default="h2o_new")
    parser.add_argument("--outdir", type=str, default="/tmp/xpdb_new_benchmark")

    parser.add_argument("--deform", type=str, default="shift", choices=["shift", "scale"])
    parser.add_argument("--direction", type=str, default="1,0,0")
    parser.add_argument("--shift", type=float, default=2.0)
    parser.add_argument("--scale", type=float, default=1.2)
    parser.add_argument("--atom", type=int, default=1)

    parser.add_argument("--no_angle", type=int, default=1, help="1 disables angle-derived H-H constraint (for H2O)")

    parser.add_argument("--dt", type=float, default=0.05)
    parser.add_argument("--omega", type=float, default=0.8)
    parser.add_argument("--momentum_beta", type=float, default=0.0)
    parser.add_argument("--k_coll", type=float, default=0.0)
    parser.add_argument("--use_inertia", type=int, default=0, help="1 keeps pred_pos fixed (inertial anchor); 0 resets pred_pos=pos each iteration")

    parser.add_argument("--iters", type=int, default=30)
    parser.add_argument("--stop_tol", type=float, default=1e-4)
    parser.add_argument("--print_every", type=int, default=1)

    parser.add_argument("--dt_list", type=str, default="0.02,0.05,0.1")
    parser.add_argument("--omega_list", type=str, default="0.5,0.8,1.0")
    parser.add_argument("--beta_list", type=str, default="0.0,0.3,0.5,0.7")

    parser.add_argument("--plot", type=int, default=0)
    parser.add_argument("--noshow", type=int, default=1)
    parser.add_argument("--save_traj", type=int, default=1)

    parser.add_argument("--group_size", type=int, default=64)
    parser.add_argument("--prefer_gpu", type=int, default=1)
    parser.add_argument("--device_idx", type=int, default=0)

    parser.add_argument("--radius", type=float, default=0.3)
    parser.add_argument("--mass", type=float, default=1.0)
    parser.add_argument("--coll_scale", type=float, default=2.0)
    parser.add_argument("--bbox_scale", type=float, default=2.0)
    parser.add_argument("--Rmax", type=float, default=1.0, help="override Rmax used for build_local_topology when collisions enabled")

    args = parser.parse_args()

    sim, bonds_adj = prepare_h2o_sim(args)

    # If collisions enabled and user did not override Rmax explicitly, choose based on bonds and radius
    if args.Rmax <= 0.0:
        args.Rmax = max(float(args.radius), float(xu.bonds_to_max_L0(bonds_adj, default=1.0)))

    if args.mode == "single":
        run_single(sim, args)
    else:
        run_sweep(sim, args)


if __name__ == "__main__":
    main()
