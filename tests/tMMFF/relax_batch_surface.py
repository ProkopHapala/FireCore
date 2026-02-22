#!/usr/bin/env python3
"""
Batch relaxation of multiple MOL2 configurations on a flat surface.
Uses nSystems>1 to relax all configurations in parallel on the GPU.
Outputs final relaxed geometries to a specified directory.
"""

import os
import sys
import glob
import time
import argparse
import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(BASE)

from pyBall.AtomicSystem import AtomicSystem
from pyBall.OCL.MMFF import MMFF as MMFF_pyocl, MOL2_ATOMTYPE_ALIASES
from pyBall.OCL.MolecularDynamics import MolecularDynamics
import pyopencl as cl

np.set_printoptions(linewidth=200, suppress=True)

def download_system_pos(md, iSys, nvecs, natoms):
    """Download positions for a specific system from the packed buffer."""
    buf = np.empty((md.nSystems * nvecs, 4), dtype=np.float32)
    cl.enqueue_copy(md.queue, buf, md.buffer_dict['apos'])
    md.queue.finish()
    offset = iSys * nvecs
    return buf[offset:offset+natoms, :3].copy()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mol2-dir', type=str, required=True, help='Directory containing input MOL2 files')
    ap.add_argument('--out-dir', type=str, default='out_batch', help='Directory to save relaxed MOL2 files')
    ap.add_argument('--max-systems', type=int, default=1000, help='Maximum number of systems to process')
    ap.add_argument('--nsteps', type=int, default=10000, help='Number of MD steps')
    ap.add_argument('--dt', type=float, default=0.02, help='Time step')
    ap.add_argument('--damp', type=float, default=0.01, help='Damping factor')
    ap.add_argument('--K', type=float, default=5.0, help='Hamaker/Morse strength')
    ap.add_argument('--mode', type=int, default=1, help='1=Hamaker, 2=Morse')
    ap.add_argument('--nPBC', type=int, nargs=3, default=(1, 1, 0), help='PBC replications (x,y,z)')
    ap.add_argument('--stride', type=int, default=1000, help='Console output stride')
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    
    # Find and sort mol2 files
    mol2_files = sorted(glob.glob(os.path.join(args.mol2_dir, '*.mol2')))
    if len(mol2_files) == 0:
        print(f"No MOL2 files found in {args.mol2_dir}")
        return
    
    mol2_files = mol2_files[:args.max_systems]
    nSystems = len(mol2_files)
    print(f"Found {nSystems} MOL2 files to relax.")

    # Load all molecules
    mols = []
    for f in mol2_files:
        m = AtomicSystem(fname=f)
        # Shift z slightly up so they don't start exactly at 0 if needed (test_ditetraceno uses 1.0)
        z_shift = 1.0 - m.apos[:, 2].min()
        m.apos[:, 2] += z_shift
        mols.append(m)
        
    print(f"All {nSystems} systems loaded. First system z-range: [{mols[0].apos[:,2].min():.2f}, {mols[0].apos[:,2].max():.2f}]")

    # Build MMFF topologies (we build a list of MMFF objects to pack later)
    # We must ensure all systems share the same topological size limits for the GPU buffer
    # But currently `md.realloc` takes a single MMFF template and duplicates its size.
    # Therefore, we assume all configs in this batch have identical number of atoms/nodes/vecs
    # which is true for rigid body assembly sampling!
    
    mm_list = []
    nPBC = tuple(args.nPBC)
    
    for i, mol in enumerate(mols):
        mm = MMFF_pyocl(bTorsion=False, verbosity=0)
        mm.capping_atoms = set()
        mm.reorder_nodes_first = False
        mm.nPBC = nPBC
        mm.lvec = mol.lvec.astype(np.float32) if mol.lvec is not None else None
        mm.toMMFFsp3_loc(mol=mol, atom_types=mm.atom_types, bRealloc=True, bEPairs=False, bUFF=False)
        mm.make_back_neighs(b_cap_neighs=False)
        mm.excl = mm._make_excl_1_2_3(mm.neighs, neighCell=mm.neighCell, npbc=mm.npbc, EXCL_MAX=16)
        mm_list.append(mm)
        if i == 0:
            print(f"Topology built for template. natoms={mm.natoms}, nnode={mm.nnode}, nvecs={mm.nvecs}")

    # Initialize GPU MD with nSystems
    template_mm = mm_list[0]
    md = MolecularDynamics(nloc=32, debug_build_options='-DDBG_UFF=0', enable_nonbond=True)
    md.realloc(template_mm, nSystems=nSystems)
    md.setup_kernels()
    
    # Pack each system into the MD buffer
    for iSys in range(nSystems):
        md.pack_system(iSys, mm_list[iSys])
        
    print(f"All {nSystems} systems packed to GPU.")

    # Set MD params
    vfac = np.float32(1.0 - args.damp)
    md.toGPU('MDparams', np.array([args.dt, 1e+6, vfac, 0.0], dtype=np.float32), byte_offset=0)

    # Init velocities
    np.random.seed(42)
    avel = np.random.normal(0, 0.1, (nSystems * template_mm.nvecs, 4)).astype(np.float32)
    md.toGPU('avel', avel)

    # Set surface
    surf_REQ = (1.5, 1.0, 0.0, 0.0)
    md.kernel_params['surf_pos0']   = np.array([0.0, 0.0, 0.0, 0.0], dtype=np.float32)
    md.kernel_params['surf_normal'] = np.array([0.0, 0.0, 1.0, 0.0], dtype=np.float32)
    md.kernel_params['surf_REQ']    = np.array(surf_REQ, dtype=np.float32)
    md.kernel_params['surf_param']  = np.array([args.K, float(args.mode), 0.0, 0.0], dtype=np.float32)

    print(f"Starting parallel relaxation for {args.nsteps} steps...")
    t0 = time.time()
    
    for istep in range(1, args.nsteps + 1):
        md.run_cleanForceMMFFf4()
        
        if md.kernel_args_getNonBond_ex2 is not None:
            md.run_getNonBond_ex2()
        elif md.kernel_args_getNonBond is not None:
            md.run_getNonBond()
            
        md.run_getSurfFlat()
        
        if md.kernel_args_getMMFFf4 is not None:
            md.run_getMMFFf4()
            
        if md.kernel_args_updateAtomsMMFFf4 is not None:
            md.run_updateAtomsMMFFf4()

        if istep % args.stride == 0 or istep == args.nsteps:
            # Quick status - download just forces/energies to check stability
            # But downloading everything is slow, so we just print step count
            dt_wall = time.time() - t0
            print(f"  step {istep:6d}/{args.nsteps}  wall={dt_wall:.1f}s")
            
    dt_total = time.time() - t0
    print(f"Relaxation complete in {dt_total:.1f}s. Saving outputs...")
    
    # Download and save
    natoms = template_mm.natoms
    nvecs = template_mm.nvecs
    
    for iSys in range(nSystems):
        pos = download_system_pos(md, iSys, nvecs, natoms)
        mol = mols[iSys]
        mol.apos = pos.copy()
        
        fname_in = os.path.basename(mol2_files[iSys])
        fname_out = os.path.join(args.out_dir, fname_in.replace('.mol2', '_relaxed.mol2'))
        
        comment = f"# relaxed {args.nsteps} steps"
        if mol.lvec is not None:
            lv = mol.lvec
            comment = ("#lvs %g %g %g   %g %g %g   %g %g %g" % (lv[0,0],lv[0,1],lv[0,2], lv[1,0],lv[1,1],lv[1,2], lv[2,0],lv[2,1],lv[2,2]))
            
        # simple_names=False ensures we keep full Tripos types in column 6
        mol.save_mol2(fname_out, comment=comment, simple_names=False)
        
        if iSys < 5:
            z_min = pos[:, 2].min()
            z_mean = pos[:, 2].mean()
            print(f"  Saved {fname_out} (z_min={z_min:.3f}, z_mean={z_mean:.3f})")

    print("All done!")

if __name__ == "__main__":
    main()
