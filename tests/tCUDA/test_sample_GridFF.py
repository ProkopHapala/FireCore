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

def func(X, Y, Z, fe):
    sx = np.sin(X*0.25)
    sy = np.sin(Y*0.2)
    sz = np.sin(Z*0.15)
    fe[:, :, :, 0] =  sx*0  # Pauli
    fe[:, :, :, 1] =  sy*0  # London
    fe[:, :, :, 2] =  sx**2 + sy**2 + sz**2   # Coulomb
    fe[:, :, :, 3] = (sx**2 + sy**2 + sz**2)*0 
    return fe
    
grid_data = ut.create_linear_func( func, grid_shape, sizes=None, dtype=np.float32)

grid_p0   = (0.0, 0.0, 0.0)  # Grid origin
grid_step = (0.1, 0.1, 0.1)  # Grid spacing
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
mmff.REQs[:,2] = 1.0 # Q = 1.0

# Initialize OpenCL MD
print("Initializing OpenCL MD...")
mdcl = clMD.MolecularDynamics(nloc=32)

# Initialize with our system
print("Allocating memory...")
mdcl.realloc(mmff=mmff, nSystems=1)


#use_texture = True
use_texture = False

mdcl.initGridFF(grid_shape, grid_data, grid_p0, grid_step,  use_texture=use_texture, r_damp=0.5, alpha_morse=2.0)
mdcl.setup_kernels()

# ------ 1D scan
#pos_x, fe_x = mdcl.scan_1D(nsteps=40, d=[0.1,0.0,0.0], p0=[0.0,0.0,0.0], use_texture=use_texture)
#pos_y, fe_y = mdcl.scan_1D(nsteps=40, d=[0.0,0.1,0.0], p0=[0.0,0.0,0.0], use_texture=use_texture)
#pos_z, fe_z = mdcl.scan_1D(nsteps=40, d=[0.0,0.0,0.1], p0=[0.0,0.0,0.0], use_texture=use_texture)

d = 0.01
dd=1e-4
p0=[1.0+dd,1.0+dd,1.0+dd]
pos_x, fe_x = mdcl.scan_1D(nsteps=160, d=[d,0.0,0.0], p0=p0, use_texture=use_texture)
pos_y, fe_y = mdcl.scan_1D(nsteps=160, d=[0.0,d,0.0], p0=p0, use_texture=use_texture)
pos_z, fe_z = mdcl.scan_1D(nsteps=160, d=[0.0,0.0,d], p0=p0, use_texture=use_texture)

fig, (ax1,ax2,ax3) = plt.subplots(1,3, figsize=(15,5))
ut.plot_1d_fe(pos_x[:,0], fe_x, ax=ax1, title="x-direction")
ut.plot_1d_fe(pos_y[:,1], fe_y, ax=ax2, title="y-direction")
ut.plot_1d_fe(pos_z[:,2], fe_z, ax=ax3, title="z-direction")

# pos, fe = scan_2D( ns=(50,50), du=[0.1,0.0,0.0], dv=[0.0,0.1,0.0],  p0=[0.0,0.0,0.0], use_texture=use_texture)
# plot_2d_fe(fe)

# Show the plots
plt.show()
