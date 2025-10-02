"""Minimal PyOpenCL driver for the `simulate_pi_pi_rot` kernel.

This script mirrors `pi_dynamics.py` but targets the pi–pi alignment system where
each node carries a rigid-body pi orbital. It initializes two nodes connected by
an elastic bond and integrates their coupled translational/rotational dynamics.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pyopencl as cl


def create_context() -> cl.Context:
    platforms = cl.get_platforms()
    for platform in platforms:
        gpus = platform.get_devices(device_type=cl.device_type.GPU)
        if gpus:
            return cl.Context(devices=[gpus[0]])
    for platform in platforms:
        devices = platform.get_devices()
        if devices:
            return cl.Context(devices=[devices[0]])
    raise RuntimeError("No OpenCL devices found.")


def build_program(ctx: cl.Context, kernel_path: Path) -> cl.Program:
    src = kernel_path.read_text()
    return cl.Program(ctx, src).build()

def urot( ang, i1=0, i2=1 ):
    c= np.cos(ang) 
    s=-np.sin(ang)
    arr = np.array([0.,0.,0.], dtype=np.float32)
    arr[i1] = c
    arr[i2] = s
    return arr

def normalized( arr ):
    return arr / np.linalg.norm(arr)

def run_kernel(
    prg: cl.Program,
    queue: cl.CommandQueue,
    steps: int,
    dt: float,
    damping: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n_records = steps * 4
    pos_host   = np.zeros((n_records, 4), dtype=np.float32)
    vel_host   = np.zeros((n_records, 4), dtype=np.float32)
    force_host = np.zeros((n_records, 4), dtype=np.float32)

    bL0 = 1.5

    pa = np.array([0.0, 0.0, 0.0, 0.0], dtype=np.float32)
    pb = np.array([bL0, 0.0, 0.0, 0.0], dtype=np.float32)
    # ha3 = urot(np.pi/2+0.2, 1, 2)
    # hb3 = urot(np.pi/2-0.3, 1, 2)

    ha3 = normalized(np.array([0.1,1.0, 0.15], dtype=np.float32))
    hb3 = normalized(np.array([0.2,1.0,-0.10], dtype=np.float32))
    ha = np.zeros(4, dtype=np.float32)
    hb = np.zeros(4, dtype=np.float32)
    ha[:3] = ha3
    hb[:3] = hb3

    pos_host[0] = pa
    pos_host[1] = pb
    pos_host[2] = ha
    pos_host[3] = hb

    vel_host[0] = np.zeros(4, dtype=np.float32)
    vel_host[1] = np.zeros(4, dtype=np.float32)
    vel_host[2] = np.zeros(4, dtype=np.float32)
    vel_host[3] = np.zeros(4, dtype=np.float32)

    mf = cl.mem_flags
    pos_buf   = cl.Buffer(queue.context, mf.READ_WRITE | mf.COPY_HOST_PTR, hostbuf=pos_host)
    vel_buf   = cl.Buffer(queue.context, mf.READ_WRITE | mf.COPY_HOST_PTR, hostbuf=vel_host)
    force_buf = cl.Buffer(queue.context, mf.READ_WRITE | mf.COPY_HOST_PTR, hostbuf=force_host)

    evt = prg.simulate_pi_pi_rot(queue, (1,), None,  pos_buf, vel_buf, force_buf, np.float32(dt), np.float32(damping), np.int32(steps))
    evt.wait()

    cl.enqueue_copy(queue, pos_host, pos_buf).wait()
    cl.enqueue_copy(queue, vel_host, vel_buf).wait()
    cl.enqueue_copy(queue, force_host, force_buf).wait()

    return (
        pos_host.reshape(steps, 4, 4),
        vel_host.reshape(steps, 4, 4),
        force_host.reshape(steps, 4, 4),
    )


def compute_invariants(pos: np.ndarray, vel: np.ndarray) -> dict[str, np.ndarray]:
    """Return linear/angular momentum and energies per step using simple numpy ops."""

    r_nodes = pos[:, :2, :3]      # positions of node atoms
    v_nodes = vel[:, :2, :3]
    omega   = vel[:, 2:, :3]

    masses   = np.array([1.0, 1.0], dtype=np.float32)
    inertias = np.array([1.0, 1.0], dtype=np.float32)

    p_nodes = v_nodes * masses.reshape(1, 2, 1)
    lin_mom = p_nodes.sum(axis=1)

    L_orb = np.cross(r_nodes, p_nodes).sum(axis=1)
    L_rot = (inertias.reshape(1, 2, 1) * omega).sum(axis=1)
    ang_mom = L_orb + L_rot

    ke_trans = 0.5 * (masses.reshape(1, 2) * (v_nodes * v_nodes).sum(axis=2)).sum(axis=1)
    ke_rot   = 0.5 * (inertias.reshape(1, 2) * (omega * omega).sum(axis=2)).sum(axis=1)

    pot = 0.5 * pos[:, :, 3].sum(axis=1)
    energy = ke_trans + ke_rot + pot

    return {
        "p": lin_mom,
        "L": ang_mom,
        "L_orb": L_orb,
        "L_rot": L_rot,
        "E_kin": ke_trans + ke_rot,
        "E_pot": pot,
        "E_tot": energy,
    }


def summarize_invariants(invars: dict[str, np.ndarray]) -> None:
    def drift(series: np.ndarray) -> float:
        return float(np.max(np.linalg.norm(series - series[0], axis=-1 if series.ndim == 2 else 0)))

    print("=== Invariant diagnostics ===")
    p = invars["p"]
    L = invars["L"]
    print(f"|Δp| max  {drift(p):.3e}")
    print(f"|ΔL| max  {drift(L):.3e}")
    for key in ("E_kin", "E_pot", "E_tot"):
        series = invars[key]
        print(f"{key} drift {float(np.max(np.abs(series - series[0]))):.3e}")


def plot_results(pos: np.ndarray, vel: np.ndarray, force: np.ndarray, dt: float, invariants: dict[str, np.ndarray]) -> None:
    steps = np.arange(pos.shape[0]) * dt

    energy = pos[:, 0, 3]
    torqueA = force[:, 2, :3]
    torqueB = force[:, 3, :3]

    fig_dofs, axes = plt.subplots(4, 1, figsize=(8, 12), sharex=True)

    axes[0].plot(steps, pos[:, 0, 0], label="nodeA_x")
    axes[0].plot(steps, pos[:, 0, 1], label="nodeA_y")
    axes[0].plot(steps, pos[:, 0, 2], label="nodeA_z")
    axes[0].plot(steps, pos[:, 1, 0], "--", label="nodeB_x")
    axes[0].plot(steps, pos[:, 1, 1], "--", label="nodeB_y")
    axes[0].plot(steps, pos[:, 1, 2], "--", label="nodeB_z")
    axes[0].set_ylabel("Position [Å]")
    axes[0].legend(loc="upper right")

    axes[1].plot(steps, pos[:, 2, 0], label="piA_x")
    axes[1].plot(steps, pos[:, 2, 1], label="piA_y")
    axes[1].plot(steps, pos[:, 2, 2], label="piA_z")
    axes[1].plot(steps, pos[:, 3, 0], "--", label="piB_x")
    axes[1].plot(steps, pos[:, 3, 1], "--", label="piB_y")
    axes[1].plot(steps, pos[:, 3, 2], "--", label="piB_z")
    axes[1].set_ylabel("Orientation")
    axes[1].legend(loc="upper right")

    axes[2].plot(steps, vel[:, 0, 0], label="vA_x")
    axes[2].plot(steps, vel[:, 0, 1], label="vA_y")
    axes[2].plot(steps, vel[:, 0, 2], label="vA_z")
    axes[2].plot(steps, vel[:, 1, 0], "--", label="vB_x")
    axes[2].plot(steps, vel[:, 1, 1], "--", label="vB_y")
    axes[2].plot(steps, vel[:, 1, 2], "--", label="vB_z")
    axes[2].set_ylabel("Velocity")
    axes[2].legend(loc="upper right")

    axes[3].plot(steps, np.linalg.norm(torqueA, axis=1), label="|tqA|")
    axes[3].plot(steps, np.linalg.norm(torqueB, axis=1), label="|tqB|")
    axes[3].plot(steps, energy, label="Bond energy")
    axes[3].set_ylabel("Energy / Torque")
    axes[3].legend(loc="upper right")

    fig_dofs.tight_layout()

    fig_inv, inv_axes = plt.subplots(3, 1, figsize=(8, 10), sharex=True)

    inv_axes[0].plot(steps, invariants["E_tot"], label="E_tot")
    inv_axes[0].plot(steps, invariants["E_kin"], label="E_kin")
    inv_axes[0].plot(steps, invariants["E_pot"], label="E_pot")
    inv_axes[0].set_ylabel("Energy")
    inv_axes[0].legend(loc="upper right")

    inv_axes[1].plot(steps, invariants["p"][:, 0], label="p_x")
    inv_axes[1].plot(steps, invariants["p"][:, 1], label="p_y")
    inv_axes[1].plot(steps, invariants["p"][:, 2], label="p_z")
    inv_axes[1].set_ylabel("Linear mom.")
    inv_axes[1].legend(loc="upper right")

    inv_axes[2].plot(steps, invariants["L"][:, 0], label="L_x")
    inv_axes[2].plot(steps, invariants["L"][:, 1], label="L_y")
    inv_axes[2].plot(steps, invariants["L"][:, 2], label="L_z")
    inv_axes[2].set_xlabel("Time [ps]")
    inv_axes[2].set_ylabel("Angular mom.")
    inv_axes[2].legend(loc="upper right")

    fig_inv.tight_layout()


def report_extrema(pos: np.ndarray, vel: np.ndarray, force: np.ndarray) -> None:
    labels = ("nodeA", "nodeB", "piA", "piB")

    for idx, label in enumerate(labels):
        p_min = pos[:, idx, :3].min(axis=0)
        p_max = pos[:, idx, :3].max(axis=0)
        v_min = vel[:, idx, :3].min(axis=0)
        v_max = vel[:, idx, :3].max(axis=0)
        f_slice = force[:, idx, :3]
        f_min = f_slice.min(axis=0)
        f_max = f_slice.max(axis=0)
        print(f"{label} position min {p_min} max {p_max}")
        print(f"{label} velocity min {v_min} max {v_max}")
        print(f"{label} force    min {f_min} max {f_max}")

def write_xyz(path: Path, pos: np.ndarray) -> None:
    with path.open("w") as fh:
        for frame in pos:
            pa = frame[0, :3]
            pb = frame[1, :3]
            tip_a = pa + frame[2, :3]
            tip_b = pb + frame[3, :3]

            fh.write("4\n")
            fh.write("pi-pi test frame\n")
            fh.write(f"C {pa[0]:.6f} {pa[1]:.6f} {pa[2]:.6f}\n")
            fh.write(f"C {pb[0]:.6f} {pb[1]:.6f} {pb[2]:.6f}\n")
            fh.write(f"E {tip_a[0]:.6f} {tip_a[1]:.6f} {tip_a[2]:.6f}\n")
            fh.write(f"E {tip_b[0]:.6f} {tip_b[1]:.6f} {tip_b[2]:.6f}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run pi–pi alignment OpenCL test and plot results.")
    parser.add_argument("--steps", type=int,   default=5, help="Number of integration steps")
    parser.add_argument("--dt",    type=float, default=0.01, help="Integration time step")
    parser.add_argument("--damp",  type=float, default=1.0, help="Damping factor applied each step")
    parser.add_argument("--xyz",   type=Path,  default=Path("pipi_trj.xyz"), help="XYZ trajectory output path (set to empty to skip)")
    args = parser.parse_args()

    kernel_path = Path(__file__).with_name("pipi_dynamics.cl")

    ctx = create_context()
    queue = cl.CommandQueue(ctx)
    program = build_program(ctx, kernel_path)

    pos_hist, vel_hist, force_hist = run_kernel(program, queue, args.steps, args.dt, args.damp)
    invariants = compute_invariants(pos_hist, vel_hist)

    report_extrema(pos_hist, vel_hist, force_hist)
    summarize_invariants(invariants)

    if args.xyz:
        write_xyz(args.xyz, pos_hist)
        print(f"Trajectory written to {args.xyz}")

    plot_results(pos_hist, vel_hist, force_hist, args.dt, invariants)
    plt.savefig("pipi_dynamics.png")
    plt.show()