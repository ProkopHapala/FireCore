#!/usr/bin/env python3
"""
H2O relaxation on a flat surface using GPU (PyOpenCL MMFF).
Tests both Hamaker LJ9-3 and Morse surface potentials.
Saves full trajectory to XYZ for review.
"""

import os, sys, time
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

def run_h2o_surface_relax(mode=1, mode_name="Hamaker", nsteps=1000, dt=0.01, damp=0.1, surf_z=0.0, z_start=5.0, surf_REQ=(1.5, 0.5, 0.0, 0.0), K=1.6):
    print(f"\n{'='*60}")
    print(f"H2O relaxation on flat surface — {mode_name} (mode={mode})")
    print(f"  surf_z={surf_z}, z_start={z_start}, surf_REQ={surf_REQ}, K={K}, dt={dt}, damp={damp}, nsteps={nsteps}")
    print(f"{'='*60}")

    # 1. Load H2O
    mol2_path = os.path.join(BASE, "cpp/common_resources/mol/H2O.mol2")
    mol = AtomicSystem(fname=mol2_path)
    if getattr(mol, 'ngs', None) is None:
        mol.neighs()
    print(f"Loaded H2O: natoms={mol.natoms}, enames={mol.enames}")
    print(f"  apos:\n{mol.apos}")

    # Place molecule at z_start above the surface plane
    mol.apos[:, 2] = mol.apos[:, 2] - mol.apos[:, 2].mean() + z_start
    print(f"  apos after shift to z_start={z_start}:\n{mol.apos}")

    # 2. Build MMFF topology
    mm = MMFF_pyocl(bTorsion=False, verbosity=0)
    mm.capping_atoms = set()          # treat H as node for this small molecule
    mm.reorder_nodes_first = False
    mm.toMMFFsp3_loc(mol=mol, atom_types=mm.atom_types, bRealloc=True, bEPairs=False, bUFF=False)
    mm.make_back_neighs(b_cap_neighs=False)

    print(f"MMFF built: natoms={mm.natoms}, nnode={mm.nnode}, nvecs={mm.nvecs}, ncap={mm.ncap}")
    print(f"  neighs:\n{mm.neighs}")
    print(f"  REQs:\n{mm.REQs}")
    print(f"  bLs:\n{mm.bLs}")
    print(f"  bKs:\n{mm.bKs}")
    print(f"  apars:\n{mm.apars}")

    # 3. Setup GPU
    md = MolecularDynamics(nloc=32, debug_build_options='-DDBG_UFF=0')
    md.realloc(mm, nSystems=1)
    md.setup_kernels()
    md.pack_system(0, mm)

    # Set MD params: (dt, Flimit, vel_damp_factor, 0)
    vfac = np.float32(1.0 - damp)
    md.toGPU('MDparams', np.array([dt, 1e+6, vfac, 0.0], dtype=np.float32), byte_offset=0)

    # 4. Set surface parameters
    md.kernel_params['surf_pos0']   = np.array([0.0, 0.0, surf_z, 0.0], dtype=np.float32)
    md.kernel_params['surf_normal'] = np.array([0.0, 0.0, 1.0, 0.0], dtype=np.float32)
    md.kernel_params['surf_REQ']    = np.array(surf_REQ, dtype=np.float32)
    md.kernel_params['surf_param']  = np.array([K, float(mode), 0.0, 0.0], dtype=np.float32)

    enames = mol.enames
    out_name = f"trj_h2o_{mode_name.lower()}_surf.xyz"
    print(f"Writing trajectory to {out_name}")

    t0 = time.time()
    fout = open(out_name, 'w')

    # Initial frame
    pos, F, E = download(md, mm.natoms, mm.nvecs)
    write_xyz_frame(fout, enames, pos, comment=f"# step 0 maxF 0.0 Esurf 0.0")

    for istep in range(1, nsteps + 1):
        # Force evaluation: clean → surface → bonded → integrate
        if md.kernel_args_cleanForceMMFFf4 is not None:
            md.run_cleanForceMMFFf4()
        else:
            # Zero forces manually if cleanForce kernel unavailable
            md.toGPU('aforce', np.zeros((mm.nvecs, 4), dtype=np.float32))

        md.run_getSurfFlat()

        if md.kernel_args_getMMFFf4 is not None:
            md.run_getMMFFf4()

        if md.kernel_args_updateAtomsMMFFf4 is not None:
            md.run_updateAtomsMMFFf4()

        pos, F, E = download(md, mm.natoms, mm.nvecs)
        fmax = float(np.sqrt(np.max(np.sum(F*F, axis=1)))) if F.size else 0.0
        Esurf_total = float(np.sum(E))

        if not np.isfinite(pos).all():
            print(f"ERROR: NaN/Inf in positions at step {istep}")
            break
        if fmax > 1e+4:
            print(f"ERROR: Blow-up at step {istep}, max|F|={fmax:.3e}")
            break

        # Write every frame for H2O (small molecule)
        write_xyz_frame(fout, enames, pos, comment=f"# step {istep} maxF {fmax:.6e} Esurf {Esurf_total:.6e}")

        if istep % 50 == 0 or istep == 1:
            print(f"  step {istep:5d}/{nsteps}  max|F|={fmax:.3e}  Esurf={Esurf_total:.6e}  z_O={pos[0,2]:.4f}")

    fout.close()
    dt_wall = time.time() - t0
    print(f"Done: {nsteps} steps in {dt_wall:.2f}s, trajectory in {out_name}")
    print(f"  Final positions:\n{pos}")
    print(f"  Final max|F|={fmax:.3e}")
    return out_name

if __name__ == "__main__":
    # R_eff = R_surf + R_O ≈ 1.5 + 1.75 = 3.25 for O; 1.5+1.443=2.943 for H
    # Hamaker minimum is at z=R_eff, Morse minimum at z=R_eff
    # Start at z=4.5, well above equilibrium → molecule should drop to ~3.2

    # Test 1: Hamaker LJ 9-3
    trj1 = run_h2o_surface_relax(mode=1, mode_name="Hamaker", nsteps=2000, dt=0.05, damp=0.15, surf_z=0.0, z_start=4.5, surf_REQ=(1.5, 2.0, 0.0, 0.0), K=1.6)

    # Test 2: Morse
    trj2 = run_h2o_surface_relax(mode=2, mode_name="Morse", nsteps=2000, dt=0.05, damp=0.15, surf_z=0.0, z_start=4.5, surf_REQ=(1.5, 2.0, 0.0, 0.0), K=1.6)

    print(f"\nTrajectories saved: {trj1}, {trj2}")
    print("Review them in a viewer to check molecule settles onto surface.")
