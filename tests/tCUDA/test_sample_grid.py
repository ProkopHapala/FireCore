import os
import sys
import numpy as np
import time
import matplotlib.pyplot as plt
import pyopencl as cl


sys.path.append("../../")
from pyBall.AtomicSystem import AtomicSystem
from pyBall.OCL import MMFF
#from pyBall.OCL.MMFF import 
from pyBall.OCL.MMparams import read_AtomAndElementTypes #read_element_types, read_atom_types, generate_REQs_from_atom_types
os.environ['PYOPENCL_COMPILER_OUTPUT'] = '1'
os.environ['PYOPENCL_CTX'] = '0:1'  # Enables double precision on device 0
from pyBall.OCL import MolecularDynamics as clMD
import pyBall.tests.utils as ut


# Path to the OpenCL kernel
KERNEL_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '../tmp/data/cl/relax_multi_mini.cl'))

# Default REQ parameters if not provided


# ========== Helper Functions
import pyBall.tests.utils as ut

# Path to the OpenCL kernel
KERNEL_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '../tmp/data/cl/relax_multi_mini.cl'))

# Default REQ parameters if not provided
REQ_DEFAULT = np.array([1.7, 0.1, 0.0, 0.0], dtype=np.float32)  # R, E, Q, padding



# ========== Load GridFF data
# Load grid data (RGBA: Pauli, London, Coulomb, H-bond)
grid_path = os.path.abspath(os.path.join(os.path.dirname(__file__),   '../../cpp/common_resources/NaCl_1x1_L2/Bspline_PLQd.npy'))
print(f"Loading grid from {grid_path}")
grid_data = np.load(grid_path)
print(f"Loaded grid shape: {grid_data.shape}")

# Setup GridFF
grid_shape = grid_data.shape[:3]
# Pad grid_data to 4 channels (RGBA, float32)
grid_data_ = np.zeros(grid_shape[:3] + (4,), dtype=np.float32)
grid_data_[:,:,:,:3] = grid_data
grid_data = grid_data_

# Create a test function for the grid
def func(X, Y, Z, fe):
    fe[:, :, :, 0] = np.sin(X*0.17)  # Pauli
    fe[:, :, :, 1] = np.sin(Y*0.2)   # London
    fe[:, :, :, 2] = np.sin(Z*0.25)  # Coulomb
    fe[:, :, :, 3] = fe[:, :, :, 0]**2 + fe[:, :, :, 1]**2 + fe[:, :, :, 2]**2 
    return fe
    
grid_data = ut.create_linear_func(func, grid_shape, sizes=None, dtype=np.float32)

grid_p0   = (0.0, 0.0, 0.0)  # Grid origin
grid_step = (0.1, 0.1, 1.0)  # Grid spacing
print(f"Initializing GridFF with shape {grid_shape}, p0={grid_p0}, step={grid_step}")

# ======= Create atoms and REQ parameters

def point_along_line(  p0=(0.0,0.0,0.0), d=(1.0,0.0,0.0), n=10):
    p0 = np.array(p0, dtype=np.float32)
    d = np.array (d, dtype=np.float32)
    ps = np.zeros((n, 4), dtype=np.float32)
    for i in range(n):
        ps[i,:3] = p0 + d*i
    return ps

atoms = point_along_line()

# ======= Initialize MolecularDynamics with atoms
print("\nInitializing MolecularDynamics with atoms directly...")

mdcl = clMD.MolecularDynamics(nloc=32)
mdcl.init_with_atoms(atoms)

# Initialize the grid force field
print("Initializing GridFF...")
use_texture = True
mdcl.initGridFF(grid_shape, grid_data, grid_p0, grid_step, use_texture=use_texture, r_damp=0.5, alpha_morse=2.0, bKernels=False)

mdcl.get_work_sizes()
mdcl.init_kernel_params()
mdcl.kernel_args_sampleGrid_tex = mdcl.generate_kernel_args("sampleGrid_tex", bPrint=False) 
ps_x = point_along_line(p0=(0.0,0.0,0.0), d=(1.0,0.0,0.0), n=10); fe_x = mdcl.run_sampleGrid_tex( apos=ps_x); #print("fe_x",fe_x)
ps_y = point_along_line(p0=(0.0,0.0,0.0), d=(0.0,1.0,0.0), n=10); fe_y = mdcl.run_sampleGrid_tex( apos=ps_y); #print("fe_y",fe_y)
ps_z = point_along_line(p0=(0.0,0.0,0.0), d=(0.0,0.0,1.0), n=10); fe_z = mdcl.run_sampleGrid_tex( apos=ps_z); #print("fe_z",fe_z)
fig, (ax1,ax2,ax3) = plt.subplots(1,3, figsize=(15,5))
ut.plot_1d_fe(ps_x[:,0], fe_x, ax=ax1, title="x-direction")
ut.plot_1d_fe(ps_y[:,1], fe_y, ax=ax2, title="y-direction")
ut.plot_1d_fe(ps_z[:,2], fe_z, ax=ax3, title="z-direction")

plt.show()
