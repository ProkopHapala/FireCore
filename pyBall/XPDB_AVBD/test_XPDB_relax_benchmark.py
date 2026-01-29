import argparse
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

# Silence PyOpenCL compiler output unless explicitly enabled
os.environ.setdefault("PYOPENCL_COMPILER_OUTPUT", "0")

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

# Ensure we can import XPDB from this directory
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from XPDB import XPDB
from pyBall import AtomicSystem


def _as_unit(v):
    v = np.array(v, dtype=np.float32)
    n = float(np.linalg.norm(v))
    if not np.isfinite(n) or n <= 0.0:
        raise ValueError(f"direction must be non-zero finite vector, got {v}")
    return v / n


def deform_shift_atom(pos, atom_idx=0, direction=(1.0, 0.0, 0.0), shift=2.0):
    """Shift a single atom by fixed displacement in a fixed direction."""
    d = _as_unit(direction) * float(shift)
    out = pos.copy()
    if atom_idx < 0 or atom_idx >= len(out):
        raise IndexError(f"atom_idx out of range: {atom_idx} vs n={len(out)}")
    out[atom_idx] += d
    return out


def deform_scale_along_direction(pos, direction=(1.0, 0.0, 0.0), scale=1.2, origin=None, atom_indices=None):
    """Scale coordinates along a direction by factor around an origin.

    p' = origin + (p-origin) + (scale-1) * dot(p-origin, d) * d
    """
    d = _as_unit(direction)
    out = pos.copy()
    if origin is None:
        origin = out.mean(axis=0)
    origin = np.array(origin, dtype=np.float32)
    if atom_indices is None:
        sel = np.arange(len(out), dtype=np.int32)
    else:
        sel = np.array(atom_indices, dtype=np.int32)
        if np.any(sel < 0) or np.any(sel >= len(out)):
            raise IndexError(f"atom_indices out of range: min={sel.min()} max={sel.max()} n={len(out)}")

    x = out[sel] - origin[None, :]
    s = np.dot(x, d)
    out[sel] = origin[None, :] + x + (float(scale) - 1.0) * s[:, None] * d[None, :]
    return out


def _parse_float_list(spec, name):
    vals = [float(x) for x in spec.split(",") if x.strip() != ""]
    if not vals:
        raise ValueError(f"{name} list must be non-empty")
    return np.array(vals, dtype=np.float32)


def make_h2o_geometry():
    """Return (pos, bonds_adj) for a simple H2O molecule.

    Atom order: 0=O, 1=H1, 2=H2.
    Geometry: O at origin, H atoms in XY plane.
    """
    # Approx equilibrium: r(OH)=0.9572 Å, angle(HOH)=104.52°
    r = 0.9572
    ang = np.deg2rad(104.52)
    h1 = np.array([r, 0.0, 0.0], dtype=np.float32)
    h2 = np.array([r * np.cos(ang), r * np.sin(ang), 0.0], dtype=np.float32)
    o = np.array([0.0, 0.0, 0.0], dtype=np.float32)
    pos = np.stack([o, h1, h2], axis=0).astype(np.float32)

    # Bond stiffness is arbitrary for convergence study; keep it high to see solver behavior
    k_oh = 500.0

    def dist(i, j):
        return float(np.linalg.norm(pos[i] - pos[j]))

    # Real bonds: O-H1 and O-H2
    bonds = [[] for _ in range(3)]
    L01 = dist(0, 1)
    L02 = dist(0, 2)
    bonds[0].append((1, L01, k_oh)); bonds[1].append((0, L01, k_oh))
    bonds[0].append((2, L02, k_oh)); bonds[2].append((0, L02, k_oh))

    # Angle-derived bond: H1-H2 distance implied by geometry
    k_hh = 200.0
    L12 = dist(1, 2)
    bonds[1].append((2, L12, k_hh)); bonds[2].append((1, L12, k_hh))

    return pos, bonds


def prepare_sim_for_molecule(pos0, bonds_adj, prefer_gpu=True, device_idx=0):
    n = len(pos0)
    sim = XPDB(n, group_size=64, prefer_gpu=prefer_gpu, device_idx=device_idx)

    vel0 = np.zeros_like(pos0, dtype=np.float32)
    radius = np.full(n, 0.3, dtype=np.float32)
    mass = np.full(n, 1.0, dtype=np.float32)

    sim.upload_data(pos0, vel0, radius, mass)
    sim.upload_bonds(bonds_adj)
    sim.ensure_neighbors(margin=1.0)

    return sim, vel0, radius, mass


def backup_device_positions(sim):
    import pyopencl as cl
    mf = cl.mem_flags
    cl_pos_backup = cl.Buffer(sim.ctx, mf.READ_WRITE, sim.num_atoms * 16)
    cl.enqueue_copy(sim.queue, cl_pos_backup, sim.cl_pos).wait()
    return cl_pos_backup


def restore_device_positions(sim, cl_pos_backup):
    import pyopencl as cl
    cl.enqueue_copy(sim.queue, sim.cl_pos, cl_pos_backup).wait()
    sim.reset_iteration_buffers(src_buffer=sim.cl_pos, reset_pred=True)


def run_single(sim, cl_pos_backup, args, title=""):
    iters = int(args.iters)
    traj = np.empty((iters + 1, sim.num_atoms, 3), dtype=np.float32) if args.save_traj else None

    # seed iteration buffers from backup
    restore_device_positions(sim, cl_pos_backup)

    # record initial
    if traj is not None:
        traj[0] = sim.get_positions()

    res = []
    # Use A as predictor => pure relaxation of constraints around current state
    pred_buf = sim.cl_iter_pos_A

    for i in range(iters):
        sim.jacobi_iteration(
            dt=args.dt,
            k_coll=args.k_coll,
            omega=args.omega,
            momentum_beta=args.momentum_beta,
            cheby_state=None,
            pred_buffer=pred_buf,
            collect_diag=False,
            out_pos_host=None,
        )
        r = sim.bond_residual_norms_from_pos(sim.cl_iter_pos_A)
        res.append((i + 1, r["linf"], r["l2"]))
        if traj is not None:
            traj[i + 1] = sim.get_positions(out=None)

        if args.stop_tol > 0.0 and r["linf"] < args.stop_tol:
            if traj is not None and i + 1 < iters:
                traj[i + 2 :] = traj[i + 1]
            break

    res = np.array(res, dtype=np.float32)

    if args.plot:
        plt.figure(figsize=(6, 4))
        plt.semilogy(res[:, 0], np.maximum(res[:, 1], 1e-12), label="L_inf")
        plt.semilogy(res[:, 0], np.maximum(res[:, 2], 1e-12), label="L2")
        plt.xlabel("Jacobi iteration")
        plt.ylabel("bond residual")
        plt.title(title)
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()

    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    np.savetxt(os.path.join(outdir, args.tag + "_residuals.txt"), res, header="iter linf l2")
    if traj is not None:
        np.save(os.path.join(outdir, args.tag + "_traj.npy"), traj)

    return res


def run_sweep(sim, cl_pos_backup, args):
    dt_vals = _parse_float_list(args.dt_list, "dt_list")
    beta_vals = _parse_float_list(args.beta_list, "beta_list")

    it_req = np.full((len(dt_vals), len(beta_vals)), int(args.iters) + 1, dtype=np.int32)
    linf_fin = np.full((len(dt_vals), len(beta_vals)), np.nan, dtype=np.float32)

    for idt, dt in enumerate(dt_vals):
        for ib, beta in enumerate(beta_vals):
            restore_device_positions(sim, cl_pos_backup)

            pred_buf = sim.cl_iter_pos_A
            last_linf = None
            for it in range(int(args.iters)):
                sim.jacobi_iteration(
                    dt=float(dt),
                    k_coll=args.k_coll,
                    omega=args.omega,
                    momentum_beta=float(beta),
                    cheby_state=None,
                    pred_buffer=pred_buf,
                    collect_diag=False,
                    out_pos_host=None,
                )
                r = sim.bond_residual_norms_from_pos(sim.cl_iter_pos_A)
                last_linf = r["linf"]
                if r["linf"] < args.stop_tol:
                    it_req[idt, ib] = it + 1
                    break
            linf_fin[idt, ib] = float(last_linf) if last_linf is not None else np.nan

    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)
    np.save(os.path.join(outdir, args.tag + "_itreq.npy"), it_req)
    np.save(os.path.join(outdir, args.tag + "_linf_final.npy"), linf_fin)

    if args.plot:
        plt.figure(figsize=(7, 4))
        plt.imshow(it_req, origin="lower", aspect="auto")
        plt.colorbar(label="iters to tol")
        plt.xticks(np.arange(len(beta_vals)), [f"{b:g}" for b in beta_vals], rotation=45)
        plt.yticks(np.arange(len(dt_vals)), [f"{d:g}" for d in dt_vals])
        plt.xlabel("momentum_beta")
        plt.ylabel("dt")
        plt.title(f"iters to reach L_inf < {args.stop_tol:g}")
        plt.tight_layout()

    return it_req, linf_fin


def run_beta_sweep(sim, cl_pos_backup, args):
    beta_vals = _parse_float_list(args.beta_list, "beta_list")
    tol_vals = _parse_float_list(args.accuracy_levels, "accuracy_levels")

    tol_vals = np.array(tol_vals, dtype=np.float32)
    sort_idx = np.argsort(tol_vals)
    tol_vals = tol_vals[sort_idx]

    it_limit = int(args.iters)
    it_req = np.full((len(beta_vals), len(tol_vals)), it_limit + 1, dtype=np.int32)
    linf_fin = np.full(len(beta_vals), np.nan, dtype=np.float32)

    for ib, beta in enumerate(beta_vals):
        restore_device_positions(sim, cl_pos_backup)
        pred_buf = sim.cl_iter_pos_A
        tol_hit = np.zeros(len(tol_vals), dtype=bool)

        for it in range(it_limit):
            sim.jacobi_iteration(
                dt=args.dt,
                k_coll=args.k_coll,
                omega=args.omega,
                momentum_beta=float(beta),
                cheby_state=None,
                pred_buffer=pred_buf,
                collect_diag=False,
                out_pos_host=None,
            )
            r = sim.bond_residual_norms_from_pos(sim.cl_iter_pos_A)
            linf = r["linf"]
            linf_fin[ib] = linf

            for jt, tol in enumerate(tol_vals):
                if not tol_hit[jt] and linf < tol:
                    it_req[ib, jt] = it + 1
                    tol_hit[jt] = True

            if tol_hit.all():
                break

    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)
    np.savez(
        os.path.join(outdir, args.tag + "_beta_sweep.npz"),
        beta_vals=beta_vals,
        tol_vals=tol_vals,
        iters=it_req,
        linf_final=linf_fin,
    )

    if args.plot:
        plt.figure(figsize=(7, 4))
        for jt, tol in enumerate(tol_vals):
            y = np.array(it_req[:, jt], dtype=np.float32)
            y[y > it_limit] = np.nan
            plt.plot(beta_vals, y, marker="o", label=f"tol={tol:g}")
        plt.xlabel("momentum_beta")
        plt.ylabel("iterations to reach tolerance")
        plt.title("Heavy-ball sweep")
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()

    return it_req, linf_fin


def main():
    parser = argparse.ArgumentParser(description="XPDB Jacobi relaxation benchmark (bond residual based)")
    parser.add_argument("--mode",      type=str,   default="beta_sweep", choices=["single", "sweep", "beta_sweep"], help="run mode")
    parser.add_argument("--deform",    type=str,   default="shift",  choices=["shift", "scale"], help="deformation type")
    parser.add_argument("--direction", type=str,   default="1,0,0",  help="direction vector as x,y,z")
    parser.add_argument("--shift",     type=float, default=2.0,      help="shift magnitude for deform=shift")
    parser.add_argument("--scale",     type=float, default=2.0,      help="scale factor for deform=scale")
    parser.add_argument("--atom",      type=int,   default=1,        help="atom index for deform=shift")
    parser.add_argument("--dt",        type=float, default=1.0)
    parser.add_argument("--omega",     type=float, default=1.0)
    parser.add_argument("--momentum_beta", type=float, default=0.5)
    parser.add_argument("--k_coll",     type=float, default=0.0)
    parser.add_argument("--iters",      type=int,   default=30)
    parser.add_argument("--stop_tol",   type=float, default=1e-4, help="stop when L_inf < tol (0 disables)")
    parser.add_argument("--dt_list",    type=str,   default="0.1,0.5,1.0,5.0")
    parser.add_argument("--beta_list",  type=str,   default="0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9")
    parser.add_argument("--accuracy_levels", type=str, default="0.1,0.01,0.001,0.0001")
    parser.add_argument("--plot",       type=int,   default=1)
    parser.add_argument("--save_traj",  type=int,   default=1)
    parser.add_argument("--outdir",     type=str,   default="/tmp/xpdb_benchmark")
    parser.add_argument("--tag",        type=str,   default="h2o")
    parser.add_argument("--molecule",   type=str,   default="../../cpp/common_resources/xyz/pentacene.xyz", help="Optional molecule file (.xyz/.mol/.mol2) to load instead of built-in H2O")
    parser.add_argument("--prefer_gpu", type=int,  default=1)
    parser.add_argument("--device_idx", type=int,  default=0)

    args = parser.parse_args()

    d = np.array([float(x) for x in args.direction.split(",")], dtype=np.float32)
    if d.size != 3:
        raise ValueError("--direction must have 3 comma-separated components")

    if args.molecule:
        mol_path = args.molecule
        if not os.path.isabs(mol_path):
            mol_path = os.path.join(os.path.dirname(__file__), mol_path)
        mol = AtomicSystem.AtomicSystem(fname=mol_path)
        if mol.bonds is None or len(mol.bonds) == 0:
            mol.findBonds()
        pos_def = np.array(mol.apos, dtype=np.float32)
        bonds = [[] for _ in range(len(pos_def))]
        for (ia, ib) in mol.bonds:
            ia = int(ia); ib = int(ib)
            if ia == ib:
                continue
            L = float(np.linalg.norm(pos_def[ia] - pos_def[ib]))
            bonds[ia].append((ib, L, 500.0))
            bonds[ib].append((ia, L, 500.0))
        if len(pos_def) > 0:
            atom_idx = max(0, min(int(args.atom), len(pos_def) - 1))
            if args.deform == "shift":
                pos_def = deform_shift_atom(pos_def, atom_idx=atom_idx, direction=d, shift=args.shift)
            elif args.deform == "scale":
                pos_def = deform_scale_along_direction(pos_def, direction=d, scale=args.scale)
        if args.tag == "h2o":
            args.tag = os.path.splitext(os.path.basename(mol_path))[0]
    else:
        pos_ref, bonds = make_h2o_geometry()
        if args.deform == "shift":
            pos_def = deform_shift_atom(pos_ref, atom_idx=args.atom, direction=d, shift=args.shift)
        elif args.deform == "scale":
            pos_def = deform_scale_along_direction(pos_ref, direction=d, scale=args.scale)
        else:
            raise ValueError(f"unknown deform {args.deform}")

    sim, vel0, radius, mass = prepare_sim_for_molecule(pos_def, bonds, prefer_gpu=bool(args.prefer_gpu), device_idx=args.device_idx)
    cl_pos_backup = backup_device_positions(sim)

    title = f"{args.tag} {args.deform} dt={args.dt:g} omega={args.omega:g} beta={args.momentum_beta:g}"

    if args.mode == "single":
        run_single(sim, cl_pos_backup, args, title=title)
    elif args.mode == "sweep":
        run_sweep(sim, cl_pos_backup, args)
    else:
        run_beta_sweep(sim, cl_pos_backup, args)

    if args.plot:
        plt.show()


if __name__ == "__main__":
    main()
