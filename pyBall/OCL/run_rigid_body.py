import argparse
import os
import sys
import numpy as np

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.abspath(os.path.join(THIS_DIR, "..", ".."))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

from pyBall.AtomicSystem import AtomicSystem
from pyBall.OCL.RigidBodyDynamics import RigidBodyDynamics


def load_system(xyz_path):
    if not os.path.exists(xyz_path):
        raise FileNotFoundError(f"XYZ file not found: {xyz_path}")
    system = AtomicSystem(fname=xyz_path)
    if system.apos is None:
        raise ValueError("Atomic positions not loaded")
    return system


def _prepare_positions(system):
    apos = np.asarray(system.apos, dtype=np.float64)
    apos = np.asarray(system.apos, dtype=np.float64)
    if apos.ndim != 2 or apos.shape[1] != 3:
        raise ValueError(f"Unexpected apos shape {apos.shape}")
    return apos


def compute_mass_properties(rel_pos, mass_scale):
    n_atoms = rel_pos.shape[0]
    masses = np.full(n_atoms, mass_scale, dtype=np.float64)
    total_mass = masses.sum()
    if total_mass <= 0:
        raise ValueError("Total mass must be positive")
    # inertia tensor for point masses
    inertia = np.zeros((3, 3), dtype=np.float64)
    for m, r in zip(masses, rel_pos):
        inertia += m * ((r @ r) * np.eye(3) - np.outer(r, r))
    # regularize to avoid singular matrices for planar molecules
    eps = 1e-6
    inertia[np.diag_indices(3)] += eps
    inertia_inv = np.linalg.inv(inertia)
    return float(total_mass), inertia_inv.astype(np.float32)


def make_initial_state(com_pos, rel_pos, total_mass, inertia_inv):
    pos = np.array([[com_pos[0], com_pos[1], com_pos[2], 0.0]], dtype=np.float32)
    quat = np.array([[0.0, 0.0, 0.0, 1.0]], dtype=np.float32)
    lin_mom = np.zeros((1, 4), dtype=np.float32)
    ang_mom = np.zeros((1, 4), dtype=np.float32)
    mass = np.array([total_mass], dtype=np.float32)
    inv_mass = np.array([1.0 / total_mass], dtype=np.float32)
    inertia = inertia_inv.reshape(1, 3, 3)
    atoms_body = rel_pos.astype(np.float32).reshape(1, rel_pos.shape[0], 3)
    return pos, quat, lin_mom, ang_mom, mass, inv_mass, inertia, atoms_body


def write_xyz_frame(file_handle, enames, coords, comment=""):
    natoms = len(enames)
    file_handle.write(f"{natoms}\n")
    file_handle.write(f"{comment}\n")
    for name, (x, y, z) in zip(enames, coords):
        file_handle.write(f"{name:2s} {x:16.8f} {y:16.8f} {z:16.8f}\n")


def run_simulation(xyz_path, dt, steps, iterations, mass_scale, traj_path, verbose):
    system = load_system(xyz_path)
    apos = _prepare_positions(system)
    com = apos.mean(axis=0)
    rel = apos - com
    total_mass, inertia_inv = compute_mass_properties(rel, mass_scale)
    state = make_initial_state(com, rel, total_mass, inertia_inv)

    rbd = RigidBodyDynamics()
    rbd.realloc(1, rel.shape[0])
    rbd.upload_state(*state)

    print(f"Loaded {rel.shape[0]} atoms from {xyz_path}")
    print(f"Total mass: {total_mass:.3f}  dt: {dt}  steps/launch: {steps}  iterations: {iterations}")

    with open(traj_path, "w") as traj_file:
        for it in range(iterations):
            rbd.run(steps, dt)
            outputs = rbd.download_outputs()
            coords_now = outputs['atom_positions'][0, :rel.shape[0], :3]
            write_xyz_frame(traj_file, system.enames, coords_now, comment=f"iteration {it}")
            com = outputs['pos'][0, :3]
            quat = outputs['quats'][0]
            if verbose:
                print(f"iter {it:03d}  pos=({com[0]:8.4f} {com[1]:8.4f} {com[2]:8.4f})  quat=({quat[0]:.4f} {quat[1]:.4f} {quat[2]:.4f} {quat[3]:.4f}) |q|={np.linalg.norm(quat):.4f}")
            if it + 1 < iterations:
                rbd.sync_outputs_to_inputs()

    np.savetxt("rigid_body_final_positions.txt", outputs['atom_positions'][0, :, :3], fmt="%12.6f")
    print("Final atom positions (first 5):")
    print(outputs['atom_positions'][0, :5, :3])
    print("Saved full positions to rigid_body_final_positions.txt")
    print(f"Saved trajectory to {traj_path}")



"""
python3 -u -m pyBall.OCL.run_rigid_body --xyz tests/tDFT/data/xyz/PTCDA.xyz --dt 0.01 --steps 1 --iterations 3 --verbose --traj rigid_body_traj.xyz
"""
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run rigid-body dynamics on a single molecule using OpenCL")
    parser.add_argument("--xyz",        type=str,   default="tests/tDFT/data/xyz/PTCDA.xyz",           help="Path to input XYZ file")
    parser.add_argument("--traj",       type=str,   default="rigid_body_traj.xyz", help="Output XYZ trajectory path")
    parser.add_argument("--dt",         type=float, default=0.1,                   help="Integration timestep")
    parser.add_argument("--steps",      type=int,   default=50,                    help="Number of integration steps per kernel launch")
    parser.add_argument("--iterations", type=int,   default=10,                    help="Number of kernel launches to perform")
    parser.add_argument("--mass-scale", type=float, default=1.0,                   help="Uniform atomic mass scale (amu)")
    parser.add_argument("--verbose",    type=int,   default=1,                     help="Print per-iteration summaries")
    args = parser.parse_args()
    run_simulation(
        xyz_path=args.xyz,
        dt=float(args.dt),
        steps=int(args.steps),
        iterations=int(args.iterations),
        mass_scale=float(args.mass_scale),
        traj_path=args.traj,
        verbose=bool(args.verbose),
    )