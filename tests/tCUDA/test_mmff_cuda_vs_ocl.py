import os
import sys
import numpy as np

sys.path.append("../../")
from pyBall.AtomicSystem import AtomicSystem
from pyBall.OCL import MMFF
#from pyBall.OCL.MMFF import 

from pyBall.OCL.MMparams import read_element_types, read_atom_types, generate_REQs_from_atom_types

os.environ['PYOPENCL_COMPILER_OUTPUT'] = '1'
os.environ['PYOPENCL_CTX'] = '0:1'  # Enables double precision on device 0

from pyBall.OCL import cuMMFF            as cuMD
from pyBall.OCL import MolecularDynamics as clMD

#mol = AtomicSystem( "common_resources/xyz/CH2NH.xyz" )

# Load forcefield parameters from data files
path = '../../cpp/common_resources/'
element_types_path = path + 'ElementTypes.dat'
atom_types_path    = path + 'AtomTypes.dat'

# Read element and atom types from data files
element_types = read_element_types(element_types_path)
atom_types = read_atom_types(atom_types_path, element_types)

# Create AtomTypeDict needed for MMFF.toMMFFsp3_loc function
#AtomTypeDict = create_atom_type_dict_for_mmff(atom_types, element_types)

mol = AtomicSystem(
    apos=np.array([[0.0, 0.0, 0.0],
                   [1.5, 0.0, 0.0],
                   [-.5,-1.0, 0.0],
                   [-.5, 1.0, 0.0],
                   [2.0, 1.0, 1.0]]),
    atypes=np.array([0,1, 2, 2, 2], dtype=np.int32),  # Example type indices
    enames=["C_2", "N_2", "H", "H", "H"],
    lvec=np.identity(3, dtype=np.float32),
    qs=np.array([-0.2, -0.3, +0.1, +0.1, +0.3], dtype=np.float32),
    bonds =[ (0,1), (0,2), (0,3), (1,4) ]
)

mol.neighs()
print( "mol.ngs ", mol.ngs)

# Set pi orbitals and electron pairs attributes after creation
mol.npi_list = np.array([1, 1, 0, 0,0], dtype=np.int32)
mol.nep_list = np.array([0, 1, 0, 0,0], dtype=np.int32)
mol.isNode   = np.array([1, 1, 0, 0,0], dtype=np.int32)
# Generate REQs from atom types
mol.REQs = generate_REQs_from_atom_types(mol, atom_types)


mmff = MMFF.MMFF(bTorsion=False, verbosity=1)
mmff.toMMFFsp3_loc( mol=mol, atom_types=atom_types)
for ia in range(mmff.natoms):
    mmff.printAtomConf(ia, mol)  # Replace 0 with desired atom index

# =========== Initialization of MMFF system

#exit()
#cuMD.init( mol.natoms, mol.natoms, mol.natoms, 0, 0 )

# Print MMFF dimensions to verify they are set correctly
print(f"\nMMFF Dimensions before MD:  natoms: {mmff.natoms}  nvecs: {mmff.nvecs}  nnode: {mmff.nnode}  ncap: {mmff.ncap}  ntors: {mmff.ntors}")

print("apos",   mmff.apos )
print("REQs",   mmff.REQs )
print("neighs", mmff.neighs )
print("bLs",    mmff.bLs )
print("bKs",    mmff.bKs )
print("apars",  mmff.apars )
print("Ksp",    mmff.Ksp )
print("Kpp",    mmff.Kpp )




# ===== RUN OpenCL Molecular Dynamics
print("\n\n\n################# RUN OpenCL MMFF #################")
mdcl = clMD.MolecularDynamics(nloc=32)
mdcl.realloc( mmff=mmff, nSystems=1,)   # Allocate memory for 1 system (nSystems=1) using the MMFF template
mdcl.pack_system(iSys=0, mmff=mmff)  # Pack the MMFF data into GPU buffers for system index 0
mdcl.upload_all_systems()            # Upload all system data to the GPU
mdcl.setup_kernels()                 # Set up kernels with their arguments

mdcl.run_getNonBond()
mdcl.run_getMMFFf4()
mdcl.run_updateAtomsMMFFf4()
mdcl.queue.finish()
mdcl.run_getNonBond()
mdcl.run_getMMFFf4()
mdcl.queue.finish()

#mdcl.make_MD_queue_batch(perBatch=10)
#mdcl.run_MD_batched( nsteps=1 )



#iter_done = mdcl.run_ocl_opt(niter=1, Fconv=1e-6, nPerVFs=1)
#final_pos, final_forces = mdcl.download_results()
print("################# END OpenCL MMFF #################")

#exit()

# # ===== RUN CUDA Molecular Dynamics
# print("\n\n\n################# RUN CUDA MMFF #################")
# cuMD.init( nAtoms=mmff.natoms, nnode=mmff.nnode, npbc=1, nMaxSysNeighs=4, nSystems=1 )
# cuMD.upload("apos",   mmff.apos)
# cuMD.upload("REQs",   mmff.REQs)
# cuMD.upload("neighs", mmff.neighs)
# cuMD.upload("BLs",    mmff.bLs)
# cuMD.upload("BKs",    mmff.bKs)
# cuMD.upload("MMpars", mmff.apars)
# cuMD.upload("Ksp",    mmff.Ksp)
# cuMD.upload("Kpp",    mmff.Kpp)
# cuMD.synchronize()
# #cuMD.upload("atypes", mmff.atypes)
# #cuMD.run_cleanForceMMFFf4()
# cuMD.run_getNonBond()
# #cuMD.run_getMMFFf4()
# #cuMD.run_updateAtomsMMFFf4()
# print("################# END CUDA MMFF #################")



