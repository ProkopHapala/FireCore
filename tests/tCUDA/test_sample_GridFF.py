import os
import sys
import numpy as np
import time
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall.AtomicSystem import AtomicSystem
from pyBall.OCL import MMFF
#from pyBall.OCL.MMFF import 
from pyBall.OCL.MMparams import read_AtomAndElementTypes #read_element_types, read_atom_types, generate_REQs_from_atom_types
os.environ['PYOPENCL_COMPILER_OUTPUT'] = '1'
os.environ['PYOPENCL_CTX'] = '0:1'  # Enables double precision on device 0
from pyBall.OCL import MolecularDynamics as clMD

# ========== Body

element_types, atom_types = read_AtomAndElementTypes(path='../../cpp/common_resources/')

#mol = AtomicSystem( "common_resources/xyz/CH2NH.xyz" )
mol = AtomicSystem( "common_resources/xyz/O2.xyz" )
#mol = AtomicSystem    ( "./common_resources/xyz/nHexadecan_fold.xyz" )
#mol = AtomicSystem( "./common_resources/xyz/hydropentacene_cross.xyz" )

#; print( "mol.ngs ", mol.ngs)

# =========== Initialization of MMFF system
mol.neighs() 
mmff = MMFF.MMFF(bTorsion=False, verbosity=1)
mmff.toMMFFsp3_loc( mol=mol, atom_types=atom_types)
#for ia in range(mmff.natoms):   mmff.printAtomConf(ia, mol)  # Replace 0 with desired atom index
#mmff.printArrays()

# ===== RUN OpenCL Molecular Dynamics
# print("\n\n\n################# RUN OpenCL MMFF #################")
# nPerBatch = 10
# mdcl = clMD.MolecularDynamics(nloc=32, perBatch=nPerBatch)
# mdcl.realloc( mmff=mmff, nSystems=200,)   # Allocate memory for 1 system (nSystems=1) using the MMFF template
#mdcl.realloc( mmff=mmff, nSystems=1,)   # Allocate memory for 1 system (nSystems=1) using the MMFF template
# mdcl.realloc( mmff=mmff, nSystems=5,) 
# mdcl.setup_kernels()
# for iSys in range(mdcl.nSystems):
#     mdcl.pack_system(iSys=iSys, mmff=mmff)  # Pack the MMFF data into GPU buffers for system index 0
# mdcl.upload_all_systems()   
# mdcl.init_kernel_params()         # Upload all system data to the GPU

# # accurate performance time measurement
# t0 = time.perf_counter()
# mdcl.run_MD_py(nsteps=nPerBatch)
# #mdcl.run_runMD()
# #mdcl.queue.finish()

# t1 = time.perf_counter(); print("OpenCL MD time: ", t1-t0)

# mdcl.run_getNonBond()
# mdcl.run_getMMFFf4()
# mdcl.queue.finish()

#mdcl.make_MD_queue_batch(perBatch=5)
#mdcl.run_MD_batched( nsteps=20 )

#iter_done = mdcl.run_ocl_opt(niter=1, Fconv=1e-6, nPerVFs=1)
#final_pos, final_forces = mdcl.download_results()
print("################# END OpenCL MMFF #################")

#exit()

# ========== GridFF B-spline Test
# Load grid data (RGBA: Pauli, London, Coulomb, H-bond)

grid_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../cpp/common_resources/NaCl_1x1_L2/Bspline_PLQd.npy'))
#grid_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../cpp/common_resources/NaCl_1x1_L3/Bspline_PLQd.npy'))
print(f"Loading grid from {grid_path}")
arr = np.load(grid_path)
print(f"Loaded grid shape: {arr.shape}")


mid = arr.shape[2] // 2
plt.figure()
plt.imshow(arr[:,:,mid,0], origin='lower')
plt.colorbar()
plt.title('GridFF Pauli slice (Z='+str(mid)+')')
plt.savefig('gridff_pauli_slice.png')

# Initialize MD CL and GridFF
# Prepare atomic system (example: CH2NH molecule)
mol = AtomicSystem(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'cpp', 'common_resources', 'xyz', 'CH2NH.xyz')))
mol.neighs()
mmff = MMFF.MMFF(bTorsion=False, verbosity=1)
mmff.toMMFFsp3_loc(mol=mol, atom_types=atom_types)

mdcl = clMD.MolecularDynamics(nloc=32)
mdcl.realloc(mmff=mmff, nSystems=1)

# init GridFF (no texture)
grid_shape = arr.shape[:3]
grid_p0 = (0.0, 0.0, 0.0)
grid_step = (1.0, 1.0, 1.0)
mdcl.initGridFF(grid_shape, arr.astype(np.float32), grid_p0, grid_step, use_texture=False, r_damp=0.5, alpha_morse=2.0)

mdcl.setup_kernels()
# pack and upload
for i in range(mdcl.nSystems): mdcl.pack_system(i, mmff)
mdcl.upload_all_systems()

# run kernel
mdcl.run_getNonBond_GridFF_Bspline()
# fetch and display forces for atom 0
pos, frc = mdcl.download_results()
print('Atom 0 force:', frc[0])
