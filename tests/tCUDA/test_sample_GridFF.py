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
from pyBall.tests.utils import create_linear_texture

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
grid_data = np.load(grid_path)
print(f"Loaded grid shape: {grid_data.shape}")



# Create a simple system with one atom
print("\nCreating single atom system...")
# Initialize with apos and enames directly in constructor
apos = np.array([[0.0, 0.0, 0.0]], dtype=np.float32)
enames = ['C']
mol = AtomicSystem(apos=apos, enames=enames)
mol.neighs()
#print(f"Number of atoms: {mol.natoms}")

# Initialize MMFF with minimal settings
print("Initializing MMFF...")
mmff = MMFF.MMFF(bTorsion=False, verbosity=1)
mmff.toMMFFsp3_loc(mol=mol, atom_types=atom_types)

# Initialize OpenCL MD
print("Initializing OpenCL MD...")
mdcl = clMD.MolecularDynamics(nloc=32)

# Initialize with our system
print("Allocating memory...")
mdcl.realloc(mmff=mmff, nSystems=1)

# Setup GridFF
grid_shape = grid_data.shape[:3]
# We need to pad grid_data to 4 channels (RGBA. float32)
grid_data_ = np.zeros(grid_shape[:3] + (4,), dtype=np.float32)
grid_data_[:,:,:,:3] = grid_data
grid_data = grid_data_
#grid_data = grid_data.astype(np.float32)

# We need to pad grid_data to 4 channels (RGBA. float32)
#grid_data_ = np.zeros(grid_shape[:3] + (4,), dtype=np.float32)
#grid_data_[:,:,:,:3] = grid_data

grid_data = create_linear_texture(grid_shape, sizes=None, dtype=np.float32)


# Visualize a slice of the grid
iz  = grid_shape[2] // 2
plt.figure(figsize=(9, 3))
plt.subplot(1, 3, 1); plt.imshow(grid_data[:,:,iz,0], origin='lower'); plt.title('GridFF Pauli'); plt.colorbar()
plt.subplot(1, 3, 2); plt.imshow(grid_data[:,:,iz,1], origin='lower'); plt.title('GridFF London'); plt.colorbar()
plt.subplot(1, 3, 3); plt.imshow(grid_data[:,:,iz,2], origin='lower'); plt.title('GridFF Coulomb'); plt.colorbar()
plt.savefig('gridff_pauli_slice.png')


grid_p0   = (0.0, 0.0, 0.0)  # Grid origin
grid_step = (0.1, 0.1, 1.0)  # Grid spacing
print(f"Initializing GridFF with shape {grid_shape}, p0={grid_p0}, step={grid_step}")

use_texture = True
#use_texture = False

mdcl.initGridFF(grid_shape, grid_data_, grid_p0, grid_step,  use_texture=use_texture, r_damp=0.5, alpha_morse=2.0)

# Setup kernels
print("Setting up kernels...")
mdcl.setup_kernels()


# ----- Sample 1D

# Initialize position and force grid_dataays
nsteps = 100
pos    = np.zeros((nsteps, 4), dtype=np.float32)
forces = np.zeros((nsteps,4), dtype=np.float32)
x = np.linspace(-5.0, 5.0, nsteps)  # Scan from -5 to 5 Angstroms
print("Running 1D force scan...")
for i in range(nsteps):
    # Update atom position along x-axis
    #mmff.apos[0, 0] = x[i]  # Move along x-axis
    #mmff.apos[0, 1] = x[i]  # Move along y-axis
    mmff.apos[0, 2] = x[i]  # Move along z-axis
    #mdcl.pack_system(0, mmff)  # Update positions on device
    mdcl.upload_all_systems()
    # Run the kernel
    if use_texture:
        mdcl.run_getNonBond_GridFF_Bspline_tex()
    else:
        mdcl.run_getNonBond_GridFF_Bspline()
    # Download results
    pos_i, force_i = mdcl.download_results()
    forces[i] = force_i[0]  # Get force on first (and only) atom
    pos[i] = pos_i[0]  # Store position for reference
    #if i % 10 == 0:  print(f"Step {i+1}/{nsteps}: x = {x[i]:.2f} Å, F = ({forces[i,0]:.3f}, {forces[i,1]:.3f}, {forces[i,2]:.3f}) kJ/mol/Å")

# Plot the force components
plt.figure(figsize=(10, 6))
plt.plot(x, forces[:, 0], 'b-', label='Fx')
plt.plot(x, forces[:, 1], 'g-', label='Fy')
plt.plot(x, forces[:, 2], 'r-', label='Fz')
plt.plot(x, forces[:, 3], 'k-', label='E')

plt.axhline(0, color='k', linestyle='--', alpha=0.3)
plt.axvline(0, color='k', linestyle='--', alpha=0.3)
plt.xlabel('X Position (Å)')
plt.ylabel('Force / Energy ')
plt.title('GridFF Force Scan')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
#plt.savefig('gridff_force_scan.png')
#print("Saved force scan plot to gridff_force_scan.png")


'''
# Sample 2D
# Initialize position and force grid_dataays
nsteps = 50
pos    = np.zeros((nsteps, 4), dtype=np.float32)
forces = np.zeros((nsteps,nsteps,4), dtype=np.float32)
x = np.linspace(0.0, 5.0, nsteps)  # Scan from -5 to 5 Angstroms
y = np.linspace(0.0, 5.0, nsteps)  # Scan from -5 to 5 Angstroms
print("Running 2D force scan...")
for ix in range(nsteps):
    for iy in range(nsteps):
        #mmff.apos[0, 0] = x[ix]  # Move along x-axis
        mmff.apos[0, 1] = y[iy]  # Move along y-axis
        mmff.apos[0, 2] = x[ix]  # Move along x-axis
        mdcl.upload_all_systems()
        if use_texture:
            mdcl.run_getNonBond_GridFF_Bspline_tex()
        else:
            mdcl.run_getNonBond_GridFF_Bspline()
        pos_i, force_i = mdcl.download_results()
        forces[ix, iy,:] = force_i[0,:]  # Get force on first (and only) atom
        #pos   [ix, iy,:] = pos_i[0,:3]  # Store position for reference

plt.figure(figsize=(12,4))
plt.subplot(1, 3, 1); plt.imshow(forces[:, :, 0], origin='lower'); plt.title('GridFF Force (X component)'); plt.colorbar()
plt.subplot(1, 3, 2); plt.imshow(forces[:, :, 1], origin='lower'); plt.title('GridFF Force (Y component)'); plt.colorbar()
plt.subplot(1, 3, 3); plt.imshow(forces[:, :, 2], origin='lower'); plt.title('GridFF Force (Z component)'); plt.colorbar()
plt.show()
'''

# Show the plots
plt.show()
