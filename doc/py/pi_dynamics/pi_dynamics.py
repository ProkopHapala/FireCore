"""Minimal PyOpenCL driver for the `simulate_sigma_pi_single` kernel.

This script builds the OpenCL kernel found in `pi_dynamics.cl`, runs the
single-system sigma–pi alignment test, and plots the resulting trajectories.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pyopencl as cl


def create_context() -> cl.Context:
    """Create an OpenCL context preferring the first available GPU."""

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


def run_kernel(
    prg: cl.Program,
    queue: cl.CommandQueue,
    steps: int,
    dt: float,
    damping: float,
    bRot: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Execute the sigma–pi simulation kernel and return history arrays."""

    n_records  = steps * 3
    pos_host   = np.zeros((n_records, 4), dtype=np.float32)
    vel_host   = np.zeros((n_records, 4), dtype=np.float32)
    force_host = np.zeros((n_records, 4), dtype=np.float32)

    pix = np.sqrt(2) / 2
    bL0 = 1.5;

    # --- Initial conditions (first three entries are the live state)
    node_pos = np.array([0.0, 0.0, 0.0, 0.0], dtype=np.float32)
    cap_pos  = np.array([bL0, 0.0, 0.0, 0.0], dtype=np.float32)
    pi_vec   = np.array([pix, pix, 0.0, 0.0], dtype=np.float32)

    pos_host[0] = node_pos
    pos_host[1] = cap_pos
    pos_host[2] = pi_vec

    vel_host[0] = np.array([0.0, 0.0, 0.0, 0.0], dtype=np.float32)
    vel_host[1] = np.array([0.0, 0.0, 0.0, 0.0], dtype=np.float32)  # kick the cap slightly
    vel_host[2] = np.array([0.0, 0.0, 0.0, 0.0], dtype=np.float32)

    mf = cl.mem_flags
    pos_buf   = cl.Buffer(queue.context, mf.READ_WRITE | mf.COPY_HOST_PTR, hostbuf=pos_host)
    vel_buf   = cl.Buffer(queue.context, mf.READ_WRITE | mf.COPY_HOST_PTR, hostbuf=vel_host)
    force_buf = cl.Buffer(queue.context, mf.READ_WRITE | mf.COPY_HOST_PTR, hostbuf=force_host)

    if bRot:
        evt = prg.simulate_sigma_pi_rot( queue, (1,), None, pos_buf, vel_buf,force_buf, np.float32(dt), np.float32(damping), np.int32(steps),)
    else:
        evt = prg.simulate_sigma_pi    ( queue, (1,), None, pos_buf, vel_buf,force_buf, np.float32(dt), np.float32(damping), np.int32(steps),)
    evt.wait()
    cl.enqueue_copy(queue, pos_host, pos_buf).wait()
    cl.enqueue_copy(queue, vel_host, vel_buf).wait()
    cl.enqueue_copy(queue, force_host, force_buf).wait()
    return pos_host.reshape(steps, 3, 4), vel_host.reshape(steps, 3, 4), force_host.reshape(steps, 3, 4)


def plot_results(pos: np.ndarray, vel: np.ndarray, force: np.ndarray, dt: float, label: str) -> None:
    steps = np.arange(pos.shape[0]) * dt

    energy   = pos  [:, 0, 3]
    torque   = force[:, 2, :3]

    fig, axes = plt.subplots(3, 1, figsize=(8, 10), sharex=True)

    axes[0].plot(steps, pos[:,0,0], "-", label="node x")
    axes[0].plot(steps, pos[:,0,1], "-", label="node y")
    axes[0].plot(steps, pos[:,0,2], "-", label="node z")

    axes[0].plot(steps, pos[:,1,0], "--", label="cap x")
    axes[0].plot(steps, pos[:,1,1], "--", label="cap y")
    axes[0].plot(steps, pos[:,1,2], "--", label="cap z")

    axes[0].plot(steps, pos[:,2,0], ":", label="pi x")
    axes[0].plot(steps, pos[:,2,1], ":", label="pi y")
    axes[0].plot(steps, pos[:,2,2], ":", label="pi z")


    axes[0].set_ylabel("Position [Å]")
    axes[0].set_title("Node and cap positions")
    axes[0].legend(loc="upper right")

    axes[1].plot(steps, vel[:, 0, 0], "-", label="node_x")
    axes[1].plot(steps, vel[:, 0, 1], "-", label="node_y")
    axes[1].plot(steps, vel[:, 0, 2], "-", label="node_z")

    axes[1].plot(steps, vel[:, 1, 0], "--", label="cap_x")
    axes[1].plot(steps, vel[:, 1, 1], "--", label="cap_y")
    axes[1].plot(steps, vel[:, 1, 2], "--", label="cap_z")

    axes[1].plot(steps, vel[:, 2, 0], ":", label="pi_x")
    axes[1].plot(steps, vel[:, 2, 1], ":", label="pi_y")
    axes[1].plot(steps, vel[:, 2, 2], ":", label="pi_z")

    axes[1].set_ylabel("Velocity [Å/ps]")
    axes[1].set_title("Node and cap velocities")
    axes[1].legend(loc="upper right")

    axes[2].plot(steps, energy, label="Energy")
    axes[2].plot(steps, np.linalg.norm(torque, axis=1), label="|Torque|")
    axes[2].set_xlabel("Time [ps]")
    axes[2].set_ylabel("Energy / Torque")
    axes[2].set_title("Energy and torque magnitude")
    axes[2].legend(loc="upper right")

    fig.suptitle(label)
    fig.tight_layout()


def report_extrema(pos: np.ndarray, vel: np.ndarray, force: np.ndarray) -> None:
    def fmt(vec: np.ndarray) -> str:
        return "[" + ", ".join(f"{v:+.6f}" for v in vec) + "]"

    labels = ("node", "cap", "pi")

    for idx, label in enumerate(labels):
        p_min = pos[:, idx, :3].min(axis=0)
        p_max = pos[:, idx, :3].max(axis=0)
        v_min = vel[:, idx, :3].min(axis=0)
        v_max = vel[:, idx, :3].max(axis=0)
        f_min = force[:, idx if idx < 2 else 2, :3].min(axis=0)
        f_max = force[:, idx if idx < 2 else 2, :3].max(axis=0)

        print(f"{label} position min {fmt(p_min)} max {fmt(p_max)}")
        print(f"{label} velocity min {fmt(v_min)} max {fmt(v_max)}")
        if idx < 2:
            print(f"{label} force    min {fmt(f_min)} max {fmt(f_max)}")
        else:
            print(f"{label} torque   min {fmt(f_min)} max {fmt(f_max)}")

    energy = pos[:, 0, 3]
    torque_mag = np.linalg.norm(force[:, 2, :3], axis=1)
    print(f"Energy range      [{energy.min():+.6f}, {energy.max():+.6f}]")
    print(f"|Torque| range    [{torque_mag.min():+.6f}, {torque_mag.max():+.6f}]")


def write_xyz(path: Path, pos: np.ndarray) -> None:
    with path.open("w") as fh:
        for frame in pos:
            node = frame[0, :3]
            cap = frame[1, :3]
            pi_tip = node + frame[2, :3]

            fh.write("3\n")
            fh.write("sigma-pi test frame\n")
            fh.write(f"O {node[0]:.6f} {node[1]:.6f} {node[2]:.6f}\n")
            fh.write(f"H {cap[0]:.6f} {cap[1]:.6f} {cap[2]:.6f}\n")
            fh.write(f"E {pi_tip[0]:.6f} {pi_tip[1]:.6f} {pi_tip[2]:.6f}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run sigma–pi alignment OpenCL test and plot results.")
    parser.add_argument("--steps", type=int,   default=1000,      help="Number of integration steps")
    parser.add_argument("--dt",    type=float, default=0.01,      help="Integration time step")
    parser.add_argument("--damp",  type=float, default=1.0,       help="Damping factor applied each step")
    parser.add_argument("--xyz",   type=Path,  default="trj.xyz", help="Optional XYZ trajectory output path")
    parser.add_argument("--rot",   type=bool,  default=True,      help="Enable rotational dynamics")
    args = parser.parse_args()

    kernel_path = Path(__file__).with_suffix(".cl")

    ctx = create_context()
    queue = cl.CommandQueue(ctx)
    program = build_program(ctx, kernel_path)

    pos_hist, vel_hist, force_hist = run_kernel(program, queue, args.steps, args.dt, args.damp, args.rot)

    report_extrema(pos_hist, vel_hist, force_hist)

    print(f"Final node position: {pos_hist[-1, 0, :3]}")
    print(f"Final cap position:  {pos_hist[-1, 1, :3]}")
    print(f"Final pi vector:     {pos_hist[-1, 2, :3]}")
    print(f"Final energy:        {pos_hist[-1, 0, 3]:.6f}")

    if args.xyz is not None:
        write_xyz(args.xyz, pos_hist)
        print(f"Trajectory written to {args.xyz}")

    plot_results( pos_hist, vel_hist, force_hist, args.dt, label=( "rot" if args.rot else "no_rot" )  )

    plt.savefig("pi_dynamics.png")

    plt.show()
