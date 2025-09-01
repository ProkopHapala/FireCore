#!/usr/bin/env python3
"""Simple UFF validation test comparing CPU and GPU implementations"""

import sys
import os
import argparse

# Add FireCore to path
base_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(base_path)

data_dir = os.path.join(base_path, "cpp/common_resources")

import numpy as np
from pyBall.AtomicSystem import AtomicSystem 
from pyBall import     MMFF as uff_cpp
from pyBall.OCL import UFF  as uff_ocl



def printCPU_bufs( uff_cpp ):
    # bonParams = getBuff ( "bonParams", (nbonds,2)      )
    # angParams = getBuff ( "angParams", (nangles,5)     )
    # dihParams = getBuff ( "dihParams", (ndihedrals,3)  )
    # invParams = getBuff ( "invParams", (ninversions,4) )

    # bonAtoms  = getIBuff( "bonAtoms",  (nbonds,2)      )
    # angAtoms  = getIBuff( "angAtoms",  (nangles,3)     )
    # dihAtoms  = getIBuff( "dihAtoms",  (ndihedrals,4)  )
    # invAtoms  = getIBuff( "invAtoms",  (ninversions,4) )

    # angNgs    = getIBuff( "angNgs",    (nangles,2)     )
    # dihNgs    = getIBuff( "dihNgs",    (ndihedrals,3)  )
    # invNgs    = getIBuff( "invNgs",    (ninversions,3) )

    # neighs    = getIBuff( "neighs",    (natoms,4)      )
    # neighBs   = getIBuff( "neighBs",   (natoms,4)      )

    print( "uff_cpp.bonParams\n", uff_cpp.bonParams )
    print( "uff_cpp.angParams\n", uff_cpp.angParams )
    print( "uff_cpp.dihParams\n", uff_cpp.dihParams )
    print( "uff_cpp.invParams\n", uff_cpp.invParams )

    print( "uff_cpp.bonAtoms\n", uff_cpp.bonAtoms )
    print( "uff_cpp.angAtoms\n", uff_cpp.angAtoms )
    print( "uff_cpp.dihAtoms\n", uff_cpp.dihAtoms )
    print( "uff_cpp.invAtoms\n", uff_cpp.invAtoms )

    print( "uff_cpp.angNgs\n", uff_cpp.angNgs )
    print( "uff_cpp.dihNgs\n", uff_cpp.dihNgs )
    print( "uff_cpp.invNgs\n", uff_cpp.invNgs )

    print( "uff_cpp.neighs\n", uff_cpp.neighs )
    print( "uff_cpp.neighBs\n", uff_cpp.neighBs )



def printCPU_bufs_2( uff_cl, bufs=None ):
    """Print GPU buffer contents for debugging"""
    if bufs is None:
        bufs = get_cpu_bufs(uff_cl)
    for name, buf in bufs.items():
        print(f"{name}\n", buf)


def printGPU_bufs( uff_cl, bufs=None ):
    """Print GPU buffer contents for debugging"""
    if bufs is None:
        bufs = get_gpu_bufs(uff_cl)
    for name, buf in bufs.items():
        print(f"{name}\n", buf)



def get_cpu_bufs(uff_cpp):
    """Automatically collect all CPU buffers"""
    buf_specs = {
        'bonAtoms':  ('bonAtoms', (-1,2)),
        'bonParams': ('bonParams', (-1,2)),
        'angAtoms':  ('angAtoms', (-1,4)),
        'angParams': ('angParams', (-1,5)),  # Will be processed specially
        'dihAtoms':  ('dihAtoms', (-1,4)),
        'dihParams': ('dihParams', (-1,3)),
        'invAtoms':  ('invAtoms', (-1,4)),
        'invParams': ('invParams', (-1,4)),
        'neighs':    ('neighs', (-1,4)),
        'neighBs':   ('neighBs', (-1,4))
    }
    
    bufs = {}
    for name, (attr, shape) in buf_specs.items():
        if hasattr(uff_cpp, attr):
            buf = getattr(uff_cpp, attr)
            if name == 'angParams':
                bufs[name] = np.hstack((buf[:,:4], buf[:,4:5]))
            else:
                bufs[name] = buf
    return bufs

def get_gpu_bufs(uff_cl):
    """Automatically collect all GPU buffers"""
    buf_specs = {
        'bonAtoms':  ('bonAtoms', (-1,2)),
        'bonParams': ('bonParams', (-1,2)),
        'angAtoms':  ('angAtoms', (-1,4)),
        'angParams': (('angParams1', 'angParams2_w'), (-1,5)),
        'dihAtoms':  ('dihAtoms', (-1,4)),
        'dihParams': ('dihParams', (-1,3)),
        'invAtoms':  ('invAtoms', (-1,4)),
        'invParams': ('invParams', (-1,4)),
        'neighs':    ('neighs', (-1,4)),
        'neighBs':   ('neighBs', (-1,4))
    }
    
    bufs = {}
    for name, (attrs, shape) in buf_specs.items():
        if isinstance(attrs, tuple):  # Special case for angParams
            buf1 = uff_cl.download_buf(attrs[0]).reshape(-1,4)
            buf2 = uff_cl.download_buf(attrs[1]).reshape(-1,1)
            bufs[name] = np.hstack((buf1, buf2))
        elif hasattr(uff_cl, 'download_buf'):
            bufs[name] = uff_cl.download_buf(attrs).reshape(shape)
    return bufs
    
def compare_bufs(uff_cpp, uff_cl, tol=1e-6):
    """Compare CPU and GPU buffers for consistency"""
    cpu_bufs = get_cpu_bufs(uff_cpp)
    gpu_bufs = get_gpu_bufs(uff_cl)
    
    errors = []
    for name in cpu_bufs:
        cpu = cpu_bufs[name]
        gpu = gpu_bufs[name]
        if cpu.shape != gpu.shape:
            errors.append(f"Shape mismatch in {name}: CPU {cpu.shape} vs GPU {gpu.shape}")
            continue
        if np.issubdtype(cpu.dtype, np.floating):
            diff = np.abs(cpu - gpu)
            max_diff = np.max(diff)
            if max_diff > tol:
                errors.append(f"Numerical diff in {name}: max_diff={max_diff:.2e}")
        else:
            if not np.array_equal(cpu, gpu):
                errors.append(f"Value mismatch in {name}")
    if errors:
        print("\nBuffer comparison errors:")
        for err in errors:
            print(f"- {err}")
        return False
    else:
        print("\nAll buffers match between CPU and GPU")
        return True



def run_uff_cpp( args ):        
    uff_cpp.init(
        xyz_name=args.molecule,
        sElementTypes=os.path.join(data_dir, "ElementTypes.dat"),
        sAtomTypes=os.path.join(data_dir, "AtomTypes.dat"),
        sBondTypes=os.path.join(data_dir, "BondTypes.dat"),
        sAngleTypes=os.path.join(data_dir, "AngleTypes.dat"),
        sDihedralTypes=os.path.join(data_dir, "DihedralTypes.dat"),
        bMMFF=True,
        bUFF=True,
        b141=True,
    )
    uff_cpp.setSwitches(NonBonded=-1, SurfAtoms=-1, GridFF=-1)
    #uff_cpp.setSwitchesUFF( DoAssemble=+1, DoBond=1, DoAngle=-1, DoDihedral=-1, DoInversion=-1,  SubtractBondNonBond=-1, ClampNonBonded=-1 )
    uff_cpp.setSwitchesUFF( DoAssemble=+1, DoBond=1, DoAngle=1, DoDihedral=1, DoInversion=1,  SubtractBondNonBond=1, ClampNonBonded=1 )
    print("-----------\n uff_cpp.getBuffs_UFF()  ")
    uff_cpp.getBuffs_UFF()
    print("-----------\n uff_cpp.setTrjName()  ")

    #printCPU_bufs( uff_cpp )
    printCPU_bufs_2( uff_cpp )

    uff_cpp.setTrjName("trj.xyz", savePerNsteps=1)
    print("-----------\n uff_cpp.run()  ")
    #uff_cpp.run( nstepMax=1000, dt=0.02, Fconv=1e-6, ialg=2, damping=0.1 )
    energy = uff_cpp.Es[0]
    forces = uff_cpp.fapos.copy()
    return energy, forces

def run_uff_ocl(args):
    mol      = AtomicSystem(fname=args.molecule)
    uff_cl   = uff_ocl.UFF_CL(nloc=32)
    uff_data = uff_cl.toUFF(mol)
    
    # Parse components string into dictionary
    components = {
        'bonds': 'bonds' in args.components,
        'angles': 'angles' in args.components,
        'dihedrals': 'dihedrals' in args.components,
        'inversions': 'inversions' in args.components
    }
    
    uff_cl.bDoBonds      = components['bonds']
    uff_cl.bDoAngles     = components['angles']
    uff_cl.bDoDihedrals  = components['dihedrals']
    uff_cl.bDoInversions = components['inversions']
    
    uff_cl.upload_topology_params(uff_data)
    uff_cl.upload_positions(mol.apos)
    printGPU_bufs(uff_cl)
    
    uff_cl.run_eval_step(bClearForce=True)
    
    energy = uff_cl.get_total_energy()
    forces = uff_cl.get_forces()
    return energy, forces, uff_cl

def compare_results(cpu_energy, cpu_forces, gpu_energy, gpu_forces):
    print(f"CPU Energy: {cpu_energy}")
    print(f"GPU Energy: {gpu_energy}")
    print(f"Energy diff: {abs(cpu_energy - gpu_energy)}")
    
    force_diff = np.linalg.norm(cpu_forces - gpu_forces, axis=1)
    print(f"Max force diff: {np.max(force_diff)}")
    print(f"Mean force diff: {np.mean(force_diff)}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='UFF Validation Test Script')
    # Build absolute default path to avoid CWD dependence
    default_mol = os.path.join(data_dir, 'mol', 'formic_acid.mol2')
    parser.add_argument('-m','--molecule',      type=str,     default=default_mol,             help='Molecule file (.mol, .xyz)')
    parser.add_argument('-g','--gpu',           type=int,     default=1,                     help='Use GPU/OpenCL implementation' )
    parser.add_argument('-c','--components',    default='bonds,angles,dihedrals,inversions', help='Comma-separated list of UFF components to test')
    parser.add_argument('-v','--verbose',       action='count', default=1,                   help='Increase verbosity level'    )
    parser.add_argument('-t','--tolerance',     type=float,   default=1e-6,                  help='Energy comparison tolerance' )
    args = parser.parse_args()

    cpu_energy, cpu_forces = run_uff_cpp(args)
    
    if args.gpu:
        gpu_energy, gpu_forces, uff_cl = run_uff_ocl(args)
        compare_results(cpu_energy, cpu_forces, gpu_energy, gpu_forces)
        printGPU_bufs(uff_cl)
        compare_bufs(uff_cpp, uff_cl)
    else:
        print(f"CPU Energy: {cpu_energy}")
        print(f"CPU Forces:\n{cpu_forces}")

    print("=========\nALL DONE")
