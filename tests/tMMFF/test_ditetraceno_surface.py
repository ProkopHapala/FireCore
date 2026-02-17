#!/usr/bin/env python3
"""
DiTetracenoHelice relaxation on a flat surface using GPU (PyOpenCL MMFF).
Reads lattice vectors from mol2 #lvs comment, uses nPBC=(1,1,0).
Saves trajectory every 100 steps (molecule is large: 624 atoms).
"""

import os, sys, time, argparse
import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(BASE)

from pyBall.AtomicSystem import AtomicSystem
from pyBall.OCL.MMFF import MMFF as MMFF_pyocl
from pyBall.OCL.MolecularDynamics import MolecularDynamics
import pyopencl as cl

np.set_printoptions(linewidth=200, suppress=True)

def write_xyz_frame(fout, enames, pos3, comment=""):
    nat = pos3.shape[0]
    fout.write(f"{nat}\n")
    fout.write(f"{comment}\n")
    for i in range(nat):
        fout.write(f"{enames[i]:3s} {pos3[i,0]:12.6f} {pos3[i,1]:12.6f} {pos3[i,2]:12.6f}\n")

def download(md, natoms, nvecs):
    apos_buf = np.empty((nvecs, 4), dtype=np.float32)
    af_buf   = np.empty((nvecs, 4), dtype=np.float32)
    cl.enqueue_copy(md.queue, apos_buf, md.buffer_dict['apos'])
    cl.enqueue_copy(md.queue, af_buf,   md.buffer_dict['aforce'])
    md.queue.finish()
    return apos_buf[:natoms, :3].copy(), af_buf[:natoms, :3].copy(), af_buf[:natoms, 3].copy()

def run_ditetraceno(conf_name="conf1", nsteps=5000, dt=0.005, damp=0.1, surf_z=0.0, surf_REQ=(1.5, 1.0, 0.0, 0.0), K=1.6, mode=1, mode_name="Hamaker", write_stride=100, nPBC=(1,1,0)):
    mol2_path = os.path.join(BASE, f"tests/tMolGUIapp/DiTetracenoHelice_{conf_name}.mol2")
    print(f"\n{'='*70}")
    print(f"DiTetracenoHelice {conf_name} — {mode_name} surface, nPBC=(1,1,0)")
    print(f"  mol2: {mol2_path}")
    print(f"  surf_z={surf_z}, surf_REQ={surf_REQ}, K={K}, mode={mode}")
    print(f"  dt={dt}, damp={damp}, nsteps={nsteps}, write_stride={write_stride}")
    print(f"{'='*70}")

    # 1. Load molecule
    mol = AtomicSystem(fname=mol2_path)
    if getattr(mol, 'ngs', None) is None:
        mol.neighs()
    print(f"Loaded: natoms={mol.natoms}, enames[:5]={mol.enames[:5]}")

    # Check lattice vectors from mol2 #lvs
    lvec = getattr(mol, 'lvec', None)
    if lvec is not None:
        print(f"  lvec from mol2:\n{lvec}")
    else:
        print("  WARNING: No lattice vectors found in mol2, using default 50x50x50 box")
        lvec = np.diag([50.0, 50.0, 50.0])
        mol.lvec = lvec

    # Place molecule: shift z so center of mass is at z_start above surface
    z_min = mol.apos[:, 2].min()
    z_shift = surf_z + 4.0 - z_min  # place lowest atom ~4 Angstrom above surface
    mol.apos[:, 2] += z_shift
    print(f"  z_shift={z_shift:.3f}, z_range after shift: [{mol.apos[:,2].min():.2f}, {mol.apos[:,2].max():.2f}]")

    # 2. Build MMFF topology — force_node_all mode for stability with large aromatic systems
    mm = MMFF_pyocl(bTorsion=False, verbosity=0)
    mm.capping_atoms = set()          # treat all atoms as nodes
    mm.reorder_nodes_first = False
    mm.nPBC = nPBC
    mm.lvec = lvec.astype(np.float32)
    mm.toMMFFsp3_loc(mol=mol, atom_types=mm.atom_types, bRealloc=True, bEPairs=False, bUFF=False)
    mm.make_back_neighs(b_cap_neighs=False)

    print(f"MMFF built: natoms={mm.natoms}, nnode={mm.nnode}, nvecs={mm.nvecs}, ncap={mm.ncap}")

    # 3. Setup GPU
    md = MolecularDynamics(nloc=32, debug_build_options='-DDBG_UFF=0')
    md.realloc(mm, nSystems=1)
    md.setup_kernels()
    md.pack_system(0, mm)

    # Set MD params
    vfac = np.float32(1.0 - damp)
    md.toGPU('MDparams', np.array([dt, 1e+6, vfac, 0.0], dtype=np.float32), byte_offset=0)

    # 4. Set surface parameters
    md.kernel_params['surf_pos0']   = np.array([0.0, 0.0, surf_z, 0.0], dtype=np.float32)
    md.kernel_params['surf_normal'] = np.array([0.0, 0.0, 1.0, 0.0],    dtype=np.float32)
    md.kernel_params['surf_REQ']    = np.array(surf_REQ,                  dtype=np.float32)
    md.kernel_params['surf_param']  = np.array([K, float(mode), 0.0, 0.0], dtype=np.float32)

    enames = mol.enames
    out_name = f"trj_ditetraceno_{conf_name}_{mode_name.lower()}.xyz"
    print(f"Writing trajectory (every {write_stride} steps) to {out_name}")

    t0 = time.time()
    fout = open(out_name, 'w')

    # Initial frame
    pos, F, E = download(md, mm.natoms, mm.nvecs)
    write_xyz_frame(fout, enames, pos, comment=f"# step 0")

    for istep in range(1, nsteps + 1):
        # Force evaluation: clean → surface → bonded → integrate
        if md.kernel_args_cleanForceMMFFf4 is not None:
            md.run_cleanForceMMFFf4()
        else:
            md.toGPU('aforce', np.zeros((mm.nvecs, 4), dtype=np.float32))

        md.run_getSurfFlat()

        if md.kernel_args_getMMFFf4 is not None:
            md.run_getMMFFf4()

        if md.kernel_args_updateAtomsMMFFf4 is not None:
            md.run_updateAtomsMMFFf4()

        # Download only when needed (every write_stride steps or last step)
        if (istep % write_stride == 0) or (istep == nsteps) or (istep <= 5):
            pos, F, E = download(md, mm.natoms, mm.nvecs)
            fmax = float(np.sqrt(np.max(np.sum(F*F, axis=1)))) if F.size else 0.0
            Esurf_total = float(np.sum(E))

            if not np.isfinite(pos).all():
                print(f"ERROR: NaN/Inf in positions at step {istep}")
                break
            if fmax > 1e+5:
                print(f"ERROR: Blow-up at step {istep}, max|F|={fmax:.3e}")
                break

            write_xyz_frame(fout, enames, pos, comment=f"# step {istep} maxF {fmax:.6e} Esurf {Esurf_total:.6e}")

            z_min = pos[:, 2].min()
            z_mean = pos[:, 2].mean()
            dt_wall = time.time() - t0
            print(f"  step {istep:6d}/{nsteps}  max|F|={fmax:.3e}  Esurf={Esurf_total:.4e}  z_min={z_min:.3f}  z_mean={z_mean:.3f}  wall={dt_wall:.1f}s")

    fout.close()
    dt_wall = time.time() - t0
    print(f"Done: {nsteps} steps in {dt_wall:.1f}s, trajectory in {out_name}")
    return out_name

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('--conf', default='conf1', choices=['conf1','conf2'])
    ap.add_argument('--mode', type=int, default=1, help='1=Hamaker, 2=Morse')
    ap.add_argument('--nsteps', type=int, default=5000)
    ap.add_argument('--dt', type=float, default=0.001)
    ap.add_argument('--damp', type=float, default=0.1)
    ap.add_argument('--stride', type=int, default=100)
    args = ap.parse_args()

    mode_name = "Hamaker" if args.mode == 1 else "Morse"
    trj = run_ditetraceno(
        conf_name=args.conf, nsteps=args.nsteps, dt=args.dt, damp=args.damp,
        surf_z=0.0, surf_REQ=(1.5, 1.0, 0.0, 0.0), K=1.6,
        mode=args.mode, mode_name=mode_name, write_stride=args.stride
    )
    print(f"\nTrajectory saved: {trj}")
