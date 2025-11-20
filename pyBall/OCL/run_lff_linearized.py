import argparse
import sys
from pathlib import Path

import numpy as np

from pyBall.AtomicSystem import AtomicSystem
from pyBall.OCL.MMFFL import MMFFL
from pyBall.OCL.LFFSolver import LFFSolver, MAX_NEIGHBORS, DEFAULT_WORKGROUP_SIZE
from pyBall.OCL.mmffl_export import _build_linearized

def load_molecule(path: Path) -> AtomicSystem:
    if not path.exists(): raise FileNotFoundError(f"Input molecule '{path}' not found")
    return AtomicSystem(fname=str(path))


def add_bond(bond_map, i, j, l0, k, tag):
    if k is None or not np.isfinite(k) or k <= 0.0:  k = 1.0
    if l0 is None or not np.isfinite(l0) or l0 <= 0.0: l0 = 1.0
    key = tuple(sorted((int(i), int(j))))
    entry = bond_map.get(key)
    if entry is None:
        bond_map[key] = {"l0": float(l0), "k": float(k), "tags": [str(tag)]}
        return
    prev_k = entry["k"]
    new_k = prev_k + float(k)
    if new_k <= 0.0:
        new_k = 1.0
    entry["l0"] = (entry["l0"] * prev_k + float(l0) * float(k)) / new_k
    entry["k"] = new_k
    entry["tags"].append(str(tag))


def assemble_neighbor_data(mmffl: MMFFL, *, max_neighbors: int, include_linear: bool):
    natoms = int(mmffl.natoms)
    bond_map: dict[tuple[int, int], dict[str, object]] = {}

    # Base MMFF sigma bonds from neighbor table
    for ia in range(natoms):
        neigh_row = mmffl.neighs[ia]
        for slot, ja in enumerate(neigh_row):
            ja = int(ja)
            if ja < 0 or ja >= natoms:
                continue
            l0 = mmffl._equilibrium_bond_length(int(ia), ja)
            k = 0.0
            if ia < mmffl.bKs.shape[0] and slot < mmffl.bKs.shape[1]:
                k = float(mmffl.bKs[ia, slot])
            add_bond(bond_map, ia, ja, l0, k, tag="sigma")

    # Linearized supplemental bonds (angles, pi dummies, etc.)
    if include_linear:
        for (i, j, l0, k, tag) in mmffl.linear_bonds:
            add_bond(bond_map, i, j, l0, k, tag=tag)

    neighbors = [[] for _ in range(natoms)]
    for (i, j), info in sorted(bond_map.items()):
        l0 = float(info["l0"])
        k = float(info["k"])
        tags = tuple(info["tags"])
        neighbors[i].append((j, k, l0, tags))
        neighbors[j].append((i, k, l0, tags))

    for ia, entries in enumerate(neighbors):
        if len(entries) > max_neighbors:
            raise ValueError(f"Atom {ia} has {len(entries)} neighbors (> {max_neighbors}). Increase --max-neigh or prune bonds.")

    neigh_arr = np.full((natoms, max_neighbors), -1, dtype=np.int32)
    KLs = np.zeros((natoms, max_neighbors, 2), dtype=np.float32)
    for ia, entries in enumerate(neighbors):
        for slot, (ja, k, l0, _tags) in enumerate(entries):
            neigh_arr[ia, slot] = int(ja)
            KLs[ia, slot, 0] = np.float32(k)
            KLs[ia, slot, 1] = np.float32(l0)
    return neigh_arr, KLs, neighbors


def print_particle_table(names, apos, neighbors):
    print("\n== Linearized atoms and neighbors ==")
    for ia, (name, pos, entries) in enumerate(zip(names, apos, neighbors)):
        neighbor_strs = []
        for (ja, k, l0, tags) in entries:
            tag_txt = "+".join(tags)
            neighbor_strs.append(f"{ja}(k={k:.3f}, l0={l0:.3f}, tag={tag_txt})")
        joined = ", ".join(neighbor_strs) if neighbor_strs else "<none>"
        print(f"[{ia:3d}] {name:>3s} pos=({pos[0]:8.4f} {pos[1]:8.4f} {pos[2]:8.4f}) -> {joined}")


def write_xyz_frame(path: Path, names, coords, *, step_index: int, append: bool):
    coords = np.asarray(coords, dtype=np.float32)
    if coords.ndim != 2 or coords.shape[1] < 3:
        raise ValueError(f"coords must have shape (natoms, >=3), got {coords.shape}")
    if coords.shape[1] > 3:
        coords = coords[:, :3]
    if len(names) != coords.shape[0]:
        raise ValueError(f"names length {len(names)} does not match coords rows {coords.shape[0]}")
    mode = "a" if append else "w"
    with Path(path).open(mode) as fh:
        fh.write(f"{len(names)}\n")
        fh.write(f"step {step_index}\n")
        for name, vec in zip(names, coords):
            fh.write(f"{str(name):>2s} {vec[0]:12.6f} {vec[1]:12.6f} {vec[2]:12.6f}\n")


def build_state_arrays(mmffl: MMFFL, names, *, mass_default: float):
    natoms = int(mmffl.natoms)
    apos = np.array(mmffl.apos[:natoms, :4], dtype=np.float32, copy=True)
    pos = apos

    vel = np.zeros((natoms, 4), dtype=np.float32)
    vel[:, 3] = np.float32(mass_default)
    return pos, vel


def apply_random_displacement(pos, amplitude, *, seed=None):
    amp = float(amplitude)
    if amp <= 0.0:
        return pos
    rng = np.random.default_rng(seed)
    disp = rng.uniform(-amp, amp, size=pos[:, :3].shape).astype(np.float32)
    pos[:, :3] += disp
    return pos


def initialize_solver(mmffl: MMFFL, names, neigh_arr, KLs, args):
    pos, vel = build_state_arrays(mmffl, names, mass_default=float(args.mass))
    pos = apply_random_displacement(pos, args.randamp, seed=None if args.seed is None else int(args.seed))

    if args.verbosity >= 0 and float(args.randamp) > 0.0:
        print(f"Applied random displacement with amplitude ±{float(args.randamp):.4f}")
    if args.verbosity >= 1:
        print("Initial positions (first 10):")
        for ia in range(min(10, len(names))):
            print(f"[{ia:3d}] {names[ia]:>3s} pos=({pos[ia, 0]:8.4f} {pos[ia, 1]:8.4f} {pos[ia, 2]:8.4f})")

    solver = LFFSolver(workgroup_size=int(args.workgroup), max_neighbors=int(args.max_neigh))
    solver.realloc([int(mmffl.natoms)])
    solver.upload_state(pos=pos, vel=vel, neighs=neigh_arr, KLs=KLs)
    return solver, pos, vel


def integrate_chunks(solver: LFFSolver, *, args, names, traj_path: Path, state):
    total_steps = int(args.nStep)
    kernel_steps = int(args.kernelSteps)
    if total_steps < 0: raise ValueError("nStep must be non-negative")
    if kernel_steps <= 0: raise ValueError("kernelSteps must be positive")

    write_xyz_frame(traj_path, names, state['pos'], step_index=0, append=False)

    remaining = total_steps
    step_index = 0

    while remaining > 0:
        chunk = kernel_steps if kernel_steps < remaining else remaining
        solver.set_params(dt=float(args.dt), n_outer=int(chunk), n_inner=int(args.nInner), efield=tuple(args.efield))
        if args.verbosity >= 0:  print(f"\nRunning LFF solver chunk: steps {chunk}, remaining before run {remaining}")
        solver.run()
        state = solver.download_state()
        step_index += chunk
        remaining -= chunk
        write_xyz_frame(traj_path, names, state['pos'], step_index=step_index, append=True)
        if args.verbosity >= 0:  print(f"Completed step {step_index}/{total_steps}")
    return state


def log_neighbor_parameters(names, pos, neigh_arr, KLs):
    print("\nPython neighbor topology (j, L0ij, Kij | current)")
    pos = np.asarray(pos, dtype=np.float32)
    for ia, name in enumerate(names):
        tuples = []
        for slot in range(neigh_arr.shape[1]):
            ja = int(neigh_arr[ia, slot])
            if ja < 0:
                break
            current = float(np.linalg.norm(pos[ia, :3] - pos[ja, :3]))
            lij = float(KLs[ia, slot, 1])
            kij = float(KLs[ia, slot, 0])
            tuples.append(f"({ja},{lij:.4f},{kij:.4f}|{current:.4f})")
        joined = " ".join(tuples) if tuples else "<none>"
        print(f"[{ia:3d}] {name:>3s} {joined}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run linearized force field simulation")
    parser.add_argument("--input",     type=Path,  default=Path("cpp/common_resources/xyz/CH2F2.xyz"), help="Input molecule file (xyz/mol/mol2)")
    parser.add_argument("--dt",        type=float, default=1.0, help="Time step passed to the LFF kernel")
    parser.add_argument("--nStep",     type=int,   default=5,    help="Total outer iterations to execute")
    parser.add_argument("--kernelSteps", type=int, default=1,    help="Outer iterations executed per kernel launch")
    parser.add_argument("--nInner",    type=int,   default=4,    help="Inner Jacobi iterations inside the kernel")
    parser.add_argument("--efield",    type=float, nargs=3, default=(0.0, 0.0, 1.0), help="External electric field vector")
    parser.add_argument("--workgroup", type=int,   default=DEFAULT_WORKGROUP_SIZE,   help="Local workgroup size for the kernel launch")
    parser.add_argument("--max-neigh", type=int,   default=MAX_NEIGHBORS, help="Maximum neighbors per atom passed to the kernel")
    parser.add_argument("--Lpi",       type=float, default=1.0,  help="Distance of pi dummy atoms from the host")
    parser.add_argument("--twopi",     type=int,   default=0,    help="Place pi dummies on both sides of the host")
    parser.add_argument("--Kang",      type=float, default=0.0,  help="Override stiffness for angle replacement bonds")
    parser.add_argument("--KpiA",      type=float, default=0.0,  help="Stiffness for host-pi bonds")
    parser.add_argument("--KpiB",      type=float, default=0.0,  help="Stiffness for pi-orthogonality bonds")
    parser.add_argument("--no-linear", type=int,   default=0,    help="Disable extra linearized bonds (angles, pi dummies)")
    parser.add_argument("--no-reorder", type=int,  default=1,    help="Disable node-first reordering in MMFF preprocessing")
    parser.add_argument("--mass",      type=float, default=1.0,  help="Default mass assigned to every site (vel.w)")
    parser.add_argument("--randamp",   type=float, default=0.2,  help="Uniform random displacement amplitude applied to initial positions")
    parser.add_argument("--seed",      type=int,   default=3254, help="Random seed for initial displacement (None = nondeterministic)")
    parser.add_argument("--trajectory", type=Path, default=Path("lff_trj.xyz"), help="Output XYZ trajectory file")
    parser.add_argument("--verbosity", type=int,   default=1,   help="Verbosity level (0: quiet, >=1 print topology, >=2 dump KL tensor)")
    args = parser.parse_args()

    try:
        mol = load_molecule(args.input)
    except Exception as exc:
        print(f"Failed to load molecule: {exc}", file=sys.stderr)
        raise

    if not hasattr(mol, "natoms") or mol.natoms is None:  mol.natoms = int(mol.apos.shape[0])
    if getattr(mol, "enames", None) is not None and not isinstance(mol.enames, list):  mol.enames = list(mol.enames)

    mmffl = MMFFL(
        L_pi=float(args.Lpi),
        two_pi_dummies=bool(args.twopi),
        Kang=float(args.Kang),
        Kpi_host=float(args.KpiA),
        Kpi_orth=float(args.KpiB),
        verbosity=int(args.verbosity),
        lone_pairs_pi=True,
        align_pi_vectors=True,
        reorder_nodes_first=not bool(args.no_reorder),
    )
    include_linear = not bool(args.no_linear)
    mmffl.build_linearized(mol, bUFF=False, include_linear=include_linear)

    charges = np.zeros(int(mmffl.natoms), dtype=np.float32)
    if getattr(mol, "qs", None) is not None:
        src = np.asarray(mol.qs, dtype=np.float32).flatten()
        ncopy = min(src.shape[0], charges.shape[0])
        if ncopy > 0:
            charges[:ncopy] = src[:ncopy]
    mmffl.apos[:mmffl.natoms, 3] = charges

    apos_linear, names, bonds_export, n_dummies = _build_linearized(mmffl, mol)
    if args.verbosity >= 0:
        print(f"Input atoms: {mol.apos.shape[0]} | Dummy atoms added: {n_dummies} | Total linearized atoms: {len(names)}")
        print(f"Total bonds exported (sigma + linearized): {len(bonds_export)}")

    neigh_arr, KLs, neighbor_dbg = assemble_neighbor_data(mmffl, max_neighbors=int(args.max_neigh), include_linear=include_linear)
    if args.verbosity >= 1:
        print_particle_table(names, mmffl.apos[:mmffl.natoms, :3], neighbor_dbg)
    if args.verbosity >= 2:
        print("\nKLs tensor (stiffness, l0) per neighbor slot:")
        np.set_printoptions(precision=4, suppress=True)
        print(KLs)

    solver, pos, vel = initialize_solver(mmffl, names, neigh_arr, KLs, args)
    log_neighbor_parameters(names, pos, neigh_arr, KLs)
    traj_path = Path(args.trajectory)
    state = integrate_chunks(solver, args=args, names=names, traj_path=traj_path, state={'pos': pos, 'vel': vel})

    print("\nUpdated positions (first 10):")
    for ia in range(min(10, state['pos'].shape[0])):
        vec = state['pos'][ia]
        print(f"[{ia:3d}] ({vec[0]:8.4f} {vec[1]:8.4f} {vec[2]:8.4f}) q={vec[3]:8.4f}")
    print(f"Trajectory written to {traj_path}")

