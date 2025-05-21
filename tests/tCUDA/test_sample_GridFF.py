from ast import Import
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
import pyBall.tests.utils as ut

# ========== Fucntions

def scan_1D( nsteps=100, d=[0.1,0.0,0.0], p0=[0.0,0.0,0.0], use_texture=False ):
    # Initialize position and force grid_dataays
    pos    = np.zeros((nsteps,4), dtype=np.float32)
    forces = np.zeros((nsteps,4), dtype=np.float32)
    print("Running 1D force scan...")
    d  = np.array(d, dtype=np.float32)
    p0 = np.array(p0,  dtype=np.float32)
    for i in range(nsteps):
        # Update atom position along x-axis
        mmff.apos[0,:3] = p0 + d*i
        mdcl.upload_all_systems()
        if use_texture:
            mdcl.run_getNonBond_GridFF_Bspline_tex()
        else:
            mdcl.run_getNonBond_GridFF_Bspline()
        pos_i, force_i = mdcl.download_results()
        forces[i] = force_i[0]  # Get force on first (and only) atom
        pos[i] = pos_i[0]  # Store position for reference
        #if i % 10 == 0:  print(f"Step {i+1}/{nsteps}: x = {x[i]:.2f} Å, F = ({forces[i,0]:.3f}, {forces[i,1]:.3f}, {forces[i,2]:.3f}) kJ/mol/Å")
    return pos, forces

def scan_2D( ns=(50,50), du=[0.1,0.0,0.0], dv=[0.0,0.1,0.0],  p0=[0.0,0.0,0.0], use_texture=False ):
    pos    = np.zeros((ns[0],ns[1],4), dtype=np.float32)
    forces = np.zeros((ns[0],ns[1],4), dtype=np.float32)
    du = np.array(du, dtype=np.float32)
    dv = np.array(dv, dtype=np.float32)
    p0 = np.array(p0,  dtype=np.float32)
    print("Running 2D force scan...")
    for ix in range(ns[0]):
        for iy in range(ns[1]):
            p = p0 + du*ix + dv*iy
            mmff.apos[0,:3] = p
            mdcl.upload_all_systems()
            if use_texture:
                mdcl.run_getNonBond_GridFF_Bspline_tex()
            else:
                mdcl.run_getNonBond_GridFF_Bspline()
            pos_i, force_i = mdcl.download_results()
            forces[ix, iy,:] = force_i[0,: ]  # Get force on first (and only) atom
            pos   [ix, iy,:] = pos_i  [0,: ]  # Store position for reference
    return pos, forces

def plot_1d_fe(x, fe, mask=(1,1,1,1), ax=None, title=None):
    # Plot the force components
    if ax is None: fig,ax = plt.subplots(figsize=(5,5))
    print( x.shape, fe.shape)
    if mask[0]: ax.plot(x, fe[:,0], 'b-', label='Fx')
    if mask[1]: ax.plot(x, fe[:,1], 'g-', label='Fy')
    if mask[2]: ax.plot(x, fe[:,2], 'r-', label='Fz')
    if mask[3]: ax.plot(x, fe[:,3], 'k-', label='E')
    ax.set_xlabel('position')
    ax.set_ylabel('force / energy')
    if title is not None: ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.3)
    #ax.tight_layout()

    #plt.savefig('gridff_force_scan.png')
    #print("Saved force scan plot to gridff_force_scan.png")


def plot_2d_fe(fe):
    plt.figure(figsize=(12,4))
    plt.subplot(1, 3, 1); plt.imshow(fe[:, :, 0], origin='lower'); plt.title('Fx'); plt.colorbar()
    plt.subplot(1, 3, 2); plt.imshow(fe[:, :, 1], origin='lower'); plt.title('Fy'); plt.colorbar()
    plt.subplot(1, 3, 3); plt.imshow(fe[:, :, 2], origin='lower'); plt.title('Fz'); plt.colorbar()
    #plt.show()

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

# ========== Load GridFF data
# Load grid data (RGBA: Pauli, London, Coulomb, H-bond)

grid_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../cpp/common_resources/NaCl_1x1_L2/Bspline_PLQd.npy'))
#grid_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../cpp/common_resources/NaCl_1x1_L3/Bspline_PLQd.npy'))
print(f"Loading grid from {grid_path}")
grid_data = np.load(grid_path)
print(f"Loaded grid shape: {grid_data.shape}")

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

#grid_data = create_linear_texture(grid_shape, sizes=None, dtype=np.float32)

def func(X,Y,Z,fe):
    fe[:, :, :, 0] = np.sin(X*0.17)
    fe[:, :, :, 1] = np.cos(Y*0.2)
    fe[:, :, :, 2] = np.sin(Z*0.25)
    return fe
    
grid_data = ut.create_linear_func( func, grid_shape, sizes=None, dtype=np.float32)

grid_p0   = (0.0, 0.0, 0.0)  # Grid origin
grid_step = (0.1, 0.1, 1.0)  # Grid spacing
print(f"Initializing GridFF with shape {grid_shape}, p0={grid_p0}, step={grid_step}")

#plot_2d_fe(grid_data[:, :, 10])

# ======= Create mdcl

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


use_texture = True
#use_texture = False

mdcl.initGridFF(grid_shape, grid_data, grid_p0, grid_step,  use_texture=use_texture, r_damp=0.5, alpha_morse=2.0)
mdcl.setup_kernels()


# ------ 1D scan
pos_x, fe_x = scan_1D(nsteps=40, d=[0.1,0.0,0.0], p0=[0.0,0.0,0.0], use_texture=use_texture)
pos_y, fe_y = scan_1D(nsteps=40, d=[0.0,0.1,0.0], p0=[0.0,0.0,0.0], use_texture=use_texture)
pos_z, fe_z = scan_1D(nsteps=40, d=[0.0,0.0,0.1], p0=[0.0,0.0,0.0], use_texture=use_texture)
fig, (ax1,ax2,ax3) = plt.subplots(1,3, figsize=(15,5))
plot_1d_fe(pos_x[:,0], fe_x, ax=ax1, title="x-direction")
plot_1d_fe(pos_y[:,1], fe_y, ax=ax2, title="y-direction")
plot_1d_fe(pos_z[:,2], fe_z, ax=ax3, title="z-direction")

# pos, fe = scan_2D( ns=(50,50), du=[0.1,0.0,0.0], dv=[0.0,0.1,0.0],  p0=[0.0,0.0,0.0], use_texture=use_texture)
# plot_2d_fe(fe)

# Show the plots
plt.show()
