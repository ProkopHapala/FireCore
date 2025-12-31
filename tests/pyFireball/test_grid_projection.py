import argparse
import numpy as np
import os
import sys
import matplotlib.pyplot as plt

# Adjust path to FireCore/pyBall directory
sys.path.append("../../")
from pyBall import FireCore as fc
from pyBall.FireballOCL import Grid as ocl_grid

def plot_density_slices(data, title="Density Slices", cmap='magma'):
    nx, ny, nz = data.shape
    fig, axes = plt.subplots(1, 3, figsize=(20, 5))
    fig.suptitle(title)
    print( "Density range:", data.min(), data.max(), data.shape )
    im0 = axes[0].imshow(data[nx//2, :, :].T, origin='lower', cmap=cmap); axes[0].set_title("X slice"); fig.colorbar(im0, ax=axes[0])
    im1 = axes[1].imshow(data[:, ny//2, :].T, origin='lower', cmap=cmap); axes[1].set_title("Y slice"); fig.colorbar(im1, ax=axes[1])
    im2 = axes[2].imshow(data[:, :, nz//2].T, origin='lower', cmap=cmap); axes[2].set_title("Z slice"); fig.colorbar(im2, ax=axes[2])
    
def plot_density_maxproj(data, title="Density Max Projections", cmap='magma'):
    nx, ny, nz = data.shape
    fig, axes = plt.subplots(1, 3, figsize=(20, 5))
    fig.suptitle(title)
    max_x = data.max(axis=0).T        # shape ny x nz -> shown as nz x ny after .T
    max_y = data.max(axis=1).T        # shape nx x nz -> shown as nz x nx after .T
    max_z = data.max(axis=2).T        # shape nx x ny -> shown as ny x nx after .T
    im0 = axes[0].imshow(max_x, origin='lower', cmap=cmap); axes[0].set_title("Max over X (ny x nz)"); fig.colorbar(im0, ax=axes[0])
    im1 = axes[1].imshow(max_y, origin='lower', cmap=cmap); axes[1].set_title("Max over Y (nx x nz)"); fig.colorbar(im1, ax=axes[1])
    im2 = axes[2].imshow(max_z, origin='lower', cmap=cmap); axes[2].set_title("Max over Z (nx x ny)"); fig.colorbar(im2, ax=axes[2])
    print("[DEBUG] maxproj shapes max_x", max_x.shape, "max_y", max_y.shape, "max_z", max_z.shape)

def test_pentacene_projection(args):
    # 1. Load Pentacene
    xyz_path = args.xyz
    # xyz_path = "../../tests/tUFF/data/xyz/benzene.xyz"
    if not os.path.exists(xyz_path):
        print(f"Error: XYZ file not found at {xyz_path}")
        return

    with open(xyz_path, 'r') as f:
        lines = f.readlines()
    natoms = int(lines[0])
    atomTypes = []
    atomPos = []
    for line in lines[2:2+natoms]:
        parts = line.split()
        sym = parts[0]
        z = 6 if sym == 'C' else 1
        atomTypes.append(z)
        atomPos.append([float(parts[1]), float(parts[2]), float(parts[3])])
    
    atomTypes = np.array(atomTypes, dtype=np.int32)
    atomPos = np.array(atomPos, dtype=np.float64)

    # 2. Initialize FireCore and run 1 SCF step
    print("Initializing FireCore...")
    fc.initialize(atomType=atomTypes, atomPos=atomPos)
    fc.setVerbosity(1)
    
    print("Running 1 SCF iteration...")
    fc.SCF(atomPos, nmax_scf=1)
    
    # 3. Get Sparse Data
    print("Fetching sparse data...")
    dims = fc.get_HS_dims()
    neighs = fc.get_HS_neighs(dims)
    rho = fc.get_rho_sparse(dims, data=neighs).rho
    # DEBUG: inspect density matrix blocks to ensure non-zero values
    print("[DEBUG] rho shape", rho.shape)
    print("[DEBUG] iatyp (Z per atom):", neighs.iatyp)
    print("[DEBUG] num_orb table len:", len(neighs.num_orb), "head:", neighs.num_orb[:min(16, len(neighs.num_orb))])

    # Optional hack: zero rho and set a single 4x4 diagonal block (i,j) with provided values
    if hack_block is not None:
        i_atom, j_atom, vals = hack_block
        rho = np.zeros_like(rho)
        ineigh = np.where(neighs.neigh_j[i_atom] == (j_atom + 1))[0]
        if len(ineigh) == 0:
            print(f"[WARN] No neighbor entry for i={i_atom} -> j={j_atom}; hack skipped")
        else:
            slot = int(ineigh[0])
            block = np.zeros((4, 4), dtype=rho.dtype)
            for k, v in enumerate(vals[:4]):
                block[k, k] = v
            rho[i_atom, slot, :4, :4] = block
            # fill reverse if present
            ineigh_rev = np.where(neighs.neigh_j[j_atom] == (i_atom + 1))[0]
            if len(ineigh_rev) > 0:
                slot_rev = int(ineigh_rev[0])
                rho[j_atom, slot_rev, :4, :4] = block
            print(f"[DEBUG] Hack block set for i={i_atom}, j={j_atom}, vals={vals}")
    # for i in range(neighs.iatyp.shape[0]):
    #     iatyp_z = int(neighs.iatyp[i])
    #     norb_i = int(neighs.num_orb[iatyp_z - 1]) if (iatyp_z - 1) < len(neighs.num_orb) else 0
    #     for ineigh in range(neighs.neigh_max):
    #         j_raw = int(neighs.neigh_j[i, ineigh])
    #         if j_raw <= 0:
    #             continue
    #         j = j_raw - 1  # Fortran -> Python index
    #         jatyp_z = int(neighs.iatyp[j])
    #         norb_j = int(neighs.num_orb[jatyp_z - 1]) if (jatyp_z - 1) < len(neighs.num_orb) else 0
    #         block = rho[i, ineigh, :norb_i, :norb_j]
    #         abs_sum = np.sum(np.abs(block))
    #         print(f"[DEBUG] rho block i={i} (Z={iatyp_z}, norb={norb_i}) "
    #               f"j={j} (Z={jatyp_z}, norb={norb_j}) ineigh={ineigh} "
    #               f"abs_sum={abs_sum:.6e}")
    #         print(block)
    
    #exit()
    
    # 4. Setup Grid Projector
    # Use local test Fdata (contains basis/*.wf1 for H, C)
    fdata_dir = os.path.join(os.path.dirname(__file__), "Fdata", "basis")
    projector = ocl_grid.GridProjector(fdata_dir)
    
    print("Loading basis functions...")
    species_nz = sorted(list(set(atomTypes)))
    projector.load_basis(species_nz)
    
    # 5. Define Grid (centered, fixed step, rounded to block size)
    margin = args.margin
    step = args.step
    block = args.block
    pos_min_raw = np.min(atomPos, axis=0) - margin
    pos_max_raw = np.max(atomPos, axis=0) + margin
    span = pos_max_raw - pos_min_raw
    ngrid = np.ceil(span / step).astype(int)
    ngrid = ((ngrid + block - 1) // block) * block  # round up to multiple of block
    dCell = np.array([step, step, step], dtype=np.float64)
    total_span = ngrid * dCell
    center = 0.5 * (pos_min_raw + pos_max_raw)
    origin = center - 0.5 * total_span
    print(f"Grid spec: origin={origin}, step={dCell}, ngrid={ngrid}, total_span={total_span}")
    grid_spec = {
        'origin': origin,
        'dA': [dCell[0], 0, 0],
        'dB': [0, dCell[1], 0],
        'dC': [0, 0, dCell[2]],
        'ngrid': ngrid
    }
    
    # 6. Project Density
    print("Projecting density to grid on GPU...")
    # Map atom data for projector
    atoms_dict = {
        'pos': atomPos,
        'Rcut': np.array([4.5 if z == 6 else 3.5 for z in atomTypes]), # Dummy Rcuts
        'type': atomTypes
    }
    
    dens = projector.project(rho, neighs, atoms_dict, grid_spec, nMaxAtom=args.nmaxatom, use_gpu_tasks=args.gpu_tasks)

    # Report block stats and optional histogram
    if hasattr(projector, "last_block_atom_counts"):
        counts = projector.last_block_atom_counts
        max_c = counts.max() if counts.size else 0
        empty = np.sum(counts == 0)
        ones = np.sum(counts == 1)
        multi = np.sum(counts > 1)
        print(f"[DEBUG] block atom stats (post-project): max={max_c}, empty={empty}, one={ones}, multi={multi}")
        if args.plot_block_hist:
            plt.figure()
            plt.hist(counts, bins=np.arange(0, counts.max()+2)-0.5, edgecolor='k')
            plt.xlabel("Atoms per block")
            plt.ylabel("Count of blocks")
            plt.title("Histogram of atoms per block")
    
    if dens is not None:
        print(f"Density range: {dens.min()} to {dens.max()}")
        print(f"Total charge (integrated): {np.sum(dens) * np.prod(dCell)}")
        plot_density_slices(dens, title="Pentacene Density (GPU) : Slices")
        plot_density_maxproj(dens, title="Pentacene Density (GPU) : Max Projections")
    else:
        print("Projection failed or returned None.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    #parser.add_argument("--xyz",      type=str,   default="../../tests/tUFF/data/xyz/pentacene.xyz", help="Path to XYZ geometry")
    #parser.add_argument("--xyz",      type=str,   default="../../tests/tUFF/data/xyz/citosine.xyz", help="Path to XYZ geometry")
    parser.add_argument("--xyz",      type=str,   default="../../tests/tUFF/data/xyz/guanine.xyz", help="Path to XYZ geometry")

    parser.add_argument("--margin",   type=float, default=4.0,    help="Grid margin (Angstrom)")
    parser.add_argument("--step",     type=float, default=0.1,    help="Grid spacing (Angstrom)")
    parser.add_argument("--block",    type=int,   default=8,      help="Block size for tasks (voxel edge count)")
    parser.add_argument("--nmaxatom", type=int,   default=64,     help="Max atoms per block/task")
    parser.add_argument("--gpu-tasks", action="store_true",      help="Use GPU-based task builder")
    parser.add_argument("--hack-block", nargs="+", type=float, default=None,   help="i j val_s [val_px val_py val_pz]; zeros rho and sets one 4x4 diagonal block")
    parser.add_argument("--plot-block-hist", type=int, default=0, help="Plot histogram of atoms per block (from build_tasks)")
    args = parser.parse_args()

    hack_block = None
    if args.hack_block is not None:
        if len(args.hack_block) < 2:
            raise SystemExit("hack-block requires at least i j")
        i_atom = int(args.hack_block[0])
        j_atom = int(args.hack_block[1])
        vals = [float(v) for v in args.hack_block[2:]] if len(args.hack_block) > 2 else [1.0]
        hack_block = (i_atom, j_atom, vals)
    args.hack_block = hack_block

    test_pentacene_projection(args)
    plt.show()
