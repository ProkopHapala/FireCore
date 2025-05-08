import os
import sys
import numpy as np

sys.path.append("../../")
from pyBall.AtomicSystem import AtomicSystem
from pyBall.OCL import MMFF
#from pyBall.OCL.MMFF import 
from pyBall.OCL.MMparams import read_AtomAndElementTypes #read_element_types, read_atom_types, generate_REQs_from_atom_types
os.environ['PYOPENCL_COMPILER_OUTPUT'] = '1'
os.environ['PYOPENCL_CTX'] = '0:1'  # Enables double precision on device 0
from pyBall.OCL import cuMMFF            as cuMD
from pyBall.OCL import MolecularDynamics as clMD

# ========== Body

element_types, atom_types = read_AtomAndElementTypes(path='../../cpp/common_resources/')

#mol = AtomicSystem( "common_resources/xyz/CH2NH.xyz" )
mol = AtomicSystem    ( "./common_resources/xyz/nHexadecan_fold.xyz" )
#mol = AtomicSystem( "./common_resources/xyz/hydropentacene_cross.xyz" )

#; print( "mol.ngs ", mol.ngs)

# =========== Initialization of MMFF system
mol.neighs() 
mmff = MMFF.MMFF(bTorsion=False, verbosity=1)
mmff.toMMFFsp3_loc( mol=mol, atom_types=atom_types)
#for ia in range(mmff.natoms):   mmff.printAtomConf(ia, mol)  # Replace 0 with desired atom index
#mmff.printArrays()

# ===== RUN OpenCL Molecular Dynamics
print("\n\n\n################# RUN OpenCL MMFF #################")
mdcl = clMD.MolecularDynamics(nloc=32, perBatch=2)
#mdcl.realloc( mmff=mmff, nSystems=5,)   # Allocate memory for 1 system (nSystems=1) using the MMFF template
mdcl.realloc( mmff=mmff, nSystems=1,)   # Allocate memory for 1 system (nSystems=1) using the MMFF template

mdcl.setup_kernels()
mdcl.pack_system(iSys=0, mmff=mmff)  # Pack the MMFF data into GPU buffers for system index 0
mdcl.upload_all_systems()   
mdcl.init_kernel_params()         # Upload all system data to the GPU

#mdcl.run_getNonBond()
#mdcl.run_getMMFFf4()
#mdcl.run_updateAtomsMMFFf4()
mdcl.run_runMD()
mdcl.queue.finish()

# mdcl.run_getNonBond()
# mdcl.run_getMMFFf4()
# mdcl.queue.finish()

#mdcl.make_MD_queue_batch(perBatch=5)
#mdcl.run_MD_batched( nsteps=20 )



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



