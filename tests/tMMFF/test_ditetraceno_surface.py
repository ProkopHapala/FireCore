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

def _maybe_simplify(enames, simple_names=False):
    if not simple_names:
        return enames
    return [e.split('_')[0] for e in enames]

def write_xyz_frame(fout, enames, pos3, comment="", simple_names=False):
    enames = _maybe_simplify(enames, simple_names)
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

def dump_pbc_images_xyz(fname, enames, pos, lvec, nPBC=(1,1,0), simple_names=False):
    enames = _maybe_simplify(enames, simple_names)
    nx, ny, nz = nPBC
    shifts = []
    for iz in range(-nz, nz+1):
        for iy in range(-ny, ny+1):
            for ix in range(-nx, nx+1):
                shifts.append(ix*lvec[0] + iy*lvec[1] + iz*lvec[2])
    shifts = np.array(shifts, dtype=np.float32)
    nat = pos.shape[0]
    with open(fname, 'w') as f:
        f.write(f"{nat*len(shifts)}\n")
        f.write(f"# PBC images nPBC={nPBC} natoms={nat}\n")
        for ish, sh in enumerate(shifts):
            for i in range(nat):
                p = pos[i] + sh
                f.write(f"{enames[i]:3s} {p[0]:12.6f} {p[1]:12.6f} {p[2]:12.6f}\n")

def report_close_pairs_pbc(fname, pos, neighs, lvec, nPBC=(1,1,0), rcut=2.5):
    nat = pos.shape[0]
    nx, ny, nz = nPBC
    shifts = []
    for iz in range(-nz, nz+1):
        for iy in range(-ny, ny+1):
            for ix in range(-nx, nx+1):
                shifts.append((ix,iy,iz, ix*lvec[0] + iy*lvec[1] + iz*lvec[2]))
    neighs = np.asarray(neighs, dtype=np.int32)
    bonded = [set([int(j) for j in neighs[i] if j >= 0]) for i in range(nat)]
    r2cut = float(rcut*rcut)
    found = []
    for i in range(nat):
        pi = pos[i]
        for j in range(nat):
            if j == i: continue
            if j in bonded[i]: continue
            pj0 = pos[j]
            for (ix,iy,iz, sh) in shifts:
                dp = (pj0 + sh) - pi
                r2 = float(dp[0]*dp[0] + dp[1]*dp[1] + dp[2]*dp[2])
                if r2 < r2cut:
                    found.append((r2, i, j, ix, iy, iz, dp[0], dp[1], dp[2]))
    found.sort(key=lambda x: x[0])
    with open(fname, 'w') as f:
        f.write(f"# close nonbonded pairs under rcut={rcut} A, nPBC={nPBC}, nat={nat}\n")
        f.write("# columns: r  i  j  ix iy iz  dx dy dz\n")
        for (r2,i,j,ix,iy,iz,dx,dy,dz) in found[:2000]:
            f.write(f"{np.sqrt(r2):12.6f} {i:6d} {j:6d} {ix:3d} {iy:3d} {iz:3d} {dx:12.6f} {dy:12.6f} {dz:12.6f}\n")
    return len(found)

def run_ditetraceno(mol, conf_name="conf1", nsteps=5000, dt=0.01, damp=0.95, surf_z=0.0, surf_REQ=(1.5, 1.0, 0.0, 0.0), K=50.0, mode=1, mode_name="Hamaker", write_stride=100, nPBC=(1,1,0), simple_names=False):
    print(f"\n{'='*70}")
    print(f"DiTetracenoHelice {conf_name} — {mode_name} surface, nPBC=(1,1,0)")
    print(f"  surf_z={surf_z}, surf_REQ={surf_REQ}, K={K}, mode={mode}")
    print(f"  dt={dt}, damp={damp}, nsteps={nsteps}, write_stride={write_stride}")
    print(f"{'='*70}")

    # Ensure neighbors are available
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

    # Place molecule: shift z so lowest atom is very close to surface
    z_min = mol.apos[:, 2].min()
    z_shift = surf_z + 1.0 - z_min  # place lowest atom ~1 Angstrom above surface
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

    # Build packed exclusion list for non-bonded kernel getNonBond_ex2
    mm.excl = mm._make_excl_1_2_3(mm.neighs, neighCell=mm.neighCell, npbc=mm.npbc, EXCL_MAX=16)

    print(f"MMFF built: natoms={mm.natoms}, nnode={mm.nnode}, nvecs={mm.nvecs}, ncap={mm.ncap}")

    # 3. Setup GPU
    md = MolecularDynamics(nloc=32, debug_build_options='-DDBG_UFF=0', enable_nonbond=True)
    md.realloc(mm, nSystems=1)
    md.setup_kernels()
    md.pack_system(0, mm)

    # Set MD params
    vfac = np.float32(1.0 - damp)
    md.toGPU('MDparams', np.array([dt, 1e+6, vfac, 0.0], dtype=np.float32), byte_offset=0)

    # Initialize velocities to small random values to start dynamics
    np.random.seed(42)
    avel = np.random.normal(0, 0.1, (mm.nvecs, 4)).astype(np.float32)  # increased from 0.01 to 0.1
    md.toGPU('avel', avel)

    # 4. Set surface parameters
    md.kernel_params['surf_pos0']   = np.array([0.0, 0.0, surf_z, 0.0], dtype=np.float32)
    md.kernel_params['surf_normal'] = np.array([0.0, 0.0, 1.0, 0.0],    dtype=np.float32)
    md.kernel_params['surf_REQ']    = np.array(surf_REQ,                  dtype=np.float32)
    md.kernel_params['surf_param']  = np.array([K, float(mode), 0.0, 0.0], dtype=np.float32)

    enames = mol.enames
    out_name = f"trj_ditetraceno_{conf_name}_{mode_name.lower()}.xyz"
    dbg_img_name  = f"dbg_pbc_images_{conf_name}.xyz"
    dbg_pair_name = f"dbg_close_pairs_{conf_name}.txt"
    print(f"Writing trajectory (every {write_stride} steps) to {out_name}")

    t0 = time.time()
    fout = open(out_name, 'w')

    # Initial frame
    pos, F, E = download(md, mm.natoms, mm.nvecs)
    apos0 = pos.copy()  # store initial positions for RMSD calculation
    write_xyz_frame(fout, enames, pos, comment=f"# step 0", simple_names=simple_names)
    dump_pbc_images_xyz(dbg_img_name, enames, pos, lvec, nPBC=nPBC, simple_names=simple_names)
    nclose0 = report_close_pairs_pbc(dbg_pair_name, pos, mm.neighs, lvec, nPBC=nPBC, rcut=2.5)
    print(f"Initial close-pair report: nclose={nclose0} -> {dbg_pair_name}, images -> {dbg_img_name}")

    for istep in range(1, nsteps + 1):
        # Force evaluation: clean → nonbonded → surface → bonded → integrate
        want_out = (istep % write_stride == 0) or (istep == nsteps) or (istep <= 5)
        if md.kernel_args_cleanForceMMFFf4 is not None:
            md.run_cleanForceMMFFf4()
        else:
            md.toGPU('aforce', np.zeros((mm.nvecs, 4), dtype=np.float32))

        # Non-bonded with exclusions (Pauli/LJ + Coulomb)
        if md.kernel_args_getNonBond_ex2 is not None:
            md.run_getNonBond_ex2()
        elif md.kernel_args_getNonBond is not None:
            md.run_getNonBond()

        Enb_total = None
        Esurf_total = None
        Ebond_total = None
        Etot_total = None
        if want_out:
            _, _, E_nb = download(md, mm.natoms, mm.nvecs)
            Enb_total = float(np.sum(E_nb))

        md.run_getSurfFlat()

        if want_out:
            _, _, E_nbs = download(md, mm.natoms, mm.nvecs)
            Esurf_total = float(np.sum(E_nbs) - Enb_total)

        if md.kernel_args_getMMFFf4 is not None:
            md.run_getMMFFf4()

        if want_out:
            _, _, E_tot = download(md, mm.natoms, mm.nvecs)
            Etot_total = float(np.sum(E_tot))
            # NOTE: getMMFFf4() does not currently accumulate bonded energy into aforce.w
            # so Etot_total reflects (Enb + Esurf) only. Do not report Ebond by subtraction.
            Ebond_total = None

        if md.kernel_args_updateAtomsMMFFf4 is not None:
            md.run_updateAtomsMMFFf4()

        # Download only when needed (every write_stride steps or last step)
        if want_out:
            pos, F, _ = download(md, mm.natoms, mm.nvecs)
            fmax = float(np.sqrt(np.max(np.sum(F*F, axis=1)))) if F.size else 0.0

            if not np.isfinite(pos).all():
                print(f"ERROR: NaN/Inf in positions at step {istep}")
                dump_pbc_images_xyz(dbg_img_name, enames, pos, lvec, nPBC=nPBC)
                report_close_pairs_pbc(dbg_pair_name, pos, mm.neighs, lvec, nPBC=nPBC, rcut=2.5)
                break

            z_min = pos[:, 2].min()
            z_mean = pos[:, 2].mean()
            rmsd = float(np.sqrt(np.mean((pos - apos0)**2)))
            dt_wall = time.time() - t0
            if Ebond_total is None:
                write_xyz_frame(fout, enames, pos, comment=f"# step {istep} maxF {fmax:.6e} Enb {Enb_total:.6e} Esurf {Esurf_total:.6e} Etot {Etot_total:.6e} zmean {z_mean:.4f} zmin {z_min:.4f} rmsd {rmsd:.4f}", simple_names=simple_names)
                print(f"  step {istep:6d}/{nsteps}  max|F|={fmax:.3e}  Enb={Enb_total:.3e}  Esurf={Esurf_total:.3e}  Etot={Etot_total:.3e}  z_min={z_min:.3f}  z_mean={z_mean:.3f}  rmsd={rmsd:.3f}  wall={dt_wall:.1f}s")
            else:
                write_xyz_frame(fout, enames, pos, comment=f"# step {istep} maxF {fmax:.6e} Enb {Enb_total:.6e} Esurf {Esurf_total:.6e} Ebond {Ebond_total:.6e} Etot {Etot_total:.6e} zmean {z_mean:.4f} zmin {z_min:.4f} rmsd {rmsd:.4f}", simple_names=simple_names)
                print(f"  step {istep:6d}/{nsteps}  max|F|={fmax:.3e}  Enb={Enb_total:.3e}  Esurf={Esurf_total:.3e}  Ebond={Ebond_total:.3e}  Etot={Etot_total:.3e}  z_min={z_min:.3f}  z_mean={z_mean:.3f}  rmsd={rmsd:.3f}  wall={dt_wall:.1f}s")

    fout.close()
    dt_wall = time.time() - t0
    # Save final configuration as mol2 (keep full typing unless simple_names requested)
    final_mol2 = f"trj_ditetraceno_{conf_name}_{mode_name.lower()}_final.mol2"
    # update mol positions to final pos
    mol.apos = pos.copy()
    comment = ""
    if lvec is not None:
        lv = lvec
        comment = ("#lvs %g %g %g   %g %g %g   %g %g %g" % (lv[0,0],lv[0,1],lv[0,2], lv[1,0],lv[1,1],lv[1,2], lv[2,0],lv[2,1],lv[2,2]))
    mol.save_mol2(final_mol2, comment=comment, simple_names=simple_names)
    print(f"Done: {nsteps} steps in {dt_wall:.1f}s, trajectory in {out_name}, final mol2 in {final_mol2}")
    return out_name

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('--conf',               default='conf1', choices=['conf1','conf2'])
    ap.add_argument('--mode',   type=int,   default=1, help='1=Hamaker, 2=Morse')
    ap.add_argument('--nsteps', type=int,   default=10000)
    ap.add_argument('--dt',     type=float, default=0.02)
    ap.add_argument('--damp',   type=float, default=0.01)
    ap.add_argument('--K',      type=float, default=5.0, help='Hamaker/Morse strength')
    ap.add_argument('--stride', type=int, default=100)
    ap.add_argument('--simple-names', action='store_true', help='strip type suffixes (C_R -> C) in trajectory/debug XYZ output')
    ap.add_argument('--flip-z',       action='store_true', help='flip molecule along z-axis before placement')
    ap.add_argument('--shift',  type=float, nargs=3, default=(0.0,0.0,0.0), help='rigid shift (dx dy dz) applied before placement')
    args = ap.parse_args()

    mol2_path = os.path.join(BASE, f"tests/tMolGUIapp/DiTetracenoHelice_{args.conf}_out_preserve.mol2")
    mol = AtomicSystem(fname=mol2_path)
    if getattr(mol, 'ngs', None) is None:
        mol.neighs()
    if args.flip_z:
        mol.apos[:,2] *= -1.0
    if args.shift is not None:
        dx, dy, dz = args.shift
        if (dx!=0.0) or (dy!=0.0) or (dz!=0.0):
            mol.apos[:,0] += dx
            mol.apos[:,1] += dy
            mol.apos[:,2] += dz

    mode_name = "Hamaker" if args.mode == 1 else "Morse"
    trj = run_ditetraceno(
        mol=mol, conf_name=args.conf, nsteps=args.nsteps, dt=args.dt, damp=args.damp,
        surf_z=0.0, surf_REQ=(1.5, 1.0, 0.0, 0.0), K=args.K,
        mode=args.mode, mode_name=mode_name, write_stride=args.stride, simple_names=args.simple_names
    )
    print(f"\nTrajectory saved: {trj}")
