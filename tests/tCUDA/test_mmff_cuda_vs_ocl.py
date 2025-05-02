
import sys
import numpy as np

sys.path.append("../../")
from pyBall.AtomicSystem import AtomicSystem
from pyBall.OCL import MMFF
from pyBall.OCL.MMFF import AtomType, Bond, Dihedral

#from pyBall.OCL import cuMMFF            as cuMD
from pyBall.OCL import MolecularDynamics as clMD

#mol = AtomicSystem( "common_resources/xyz/CH2NH.xyz" )

# Define AtomType instances with npi and ne
AtomTypeDict = {
    "C": AtomType(name="C", Ruff= 0.757, Quff=1.912, Eaff=-4.528, Kss=300.0, Ass=109.5, Ksp=100.0, Kpp=150.0, npi=0, ne=0),
    "H": AtomType(name="H", Ruff=0.354 , Quff=0.712, Eaff=-4.528, Kss=200.0, Ass=109.5, Ksp=50.0,  Kpp=75.0,  npi=0, ne=0),
    "N": AtomType(name="N", Ruff=0.700 , Quff=2.544, Eaff=-6.899, Kss=350.0, Ass=109.5, Ksp=120.0, Kpp=180.0, npi=0, ne=0),
    "O": AtomType(name="O", Ruff=0.658 , Quff=2.300, Eaff=-8.741, Kss=350.0, Ass=109.5, Ksp=120.0, Kpp=180.0, npi=0, ne=0),
    "E": AtomType(name="E", Ruff=0.50  , Quff=0.0,   Eaff=0.0,    Kss=0.0,   Ass=0.0,   Ksp=0.0,   Kpp=0.0,   npi=0, ne=1)  # Electron pairs
}

mol = AtomicSystem(
    apos=np.array([[0.0, 0.0, 0.0],
                   [1.5, 0.0, 0.0],
                   [-.5,-1.0, 0.0],
                   [-.5, 1.0, 0.0],
                   [2.0, 1.0, 1.0]
                   ], dtype=np.float32),
    atypes=np.array([0,1, 2, 2, 2], dtype=np.int32),  # Example type indices
    enames=["C", "N", "H", "H", "H"],
    lvec=np.identity(3, dtype=np.float32),
    qs=np.array([0.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float32),
    bonds =[
        (0,1),
        (0,2),
        (0,3),
        (1,4)
    ]
    #dihedrals=[],          # No dihedrals in this simple example
)

mol.neighs()
print( "mol.ngs ", mol.ngs)


# Set pi orbitals and electron pairs attributes after creation
mol.npi_list = np.array([1, 1, 0, 0,0], dtype=np.int32)
mol.nep_list = np.array([0, 1, 0, 0,0], dtype=np.int32)
mol.isNode   = np.array([1, 1, 0, 0,0], dtype=np.int32)
mol.REQs=np.array(
    [[1.5, 0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0, 0.0]], dtype=np.float32),

# Initialize MMFF instance
mmff = MMFF.MMFF(bTorsion=False, verbosity=1)

# Convert AtomicSystem to MMFFsp3_loc representation
mmff.toMMFFsp3_loc(
    atomic_system=mol,
    AtomTypeDict=AtomTypeDict,
    bRealloc=True,
    bEPairs=True,
    bUFF=False
)

# Optionally, print atom configurations for debugging
for ia in range(mmff.natoms):
    mmff.printAtomConf(ia, mol)  # Replace 0 with desired atom index


# =========== Initialization of MMFF system


#exit()
#cuMD.init( mol.natoms, mol.natoms, mol.natoms, 0, 0 )

# Print MMFF dimensions to verify they are set correctly
print(f"\nMMFF Dimensions before MD:\n  natoms: {mmff.natoms}\n  nvecs: {mmff.nvecs}\n  nnode: {mmff.nnode}\n  ncap: {mmff.ncap}\n  ntors: {mmff.ntors}")

# Ensure MMFF dimensions are properly set if they are zero
if mmff.natoms == 0 or mmff.nvecs == 0 or mmff.nnode == 0:
    print("Fixing MMFF dimensions...")
    # Set dimensions based on the molecule structure
    mmff.natoms = mol.natoms  # Total atoms in the molecule
    mmff.nvecs = mol.natoms   # Vector elements (typically same as natoms)
    mmff.nnode = sum(mol.isNode)  # Number of nodes (atoms with configurations)
    mmff.ncap = mol.natoms - mmff.nnode  # Capping atoms (non-nodes)
    mmff.ntors = 0  # No torsions in this example
    print(f"Fixed MMFF dimensions:\n  natoms: {mmff.natoms}\n  nvecs: {mmff.nvecs}\n  nnode: {mmff.nnode}\n  ncap: {mmff.ncap}\n  ntors: {mmff.ntors}")

# Initialize MolecularDynamics with default nloc=32
md = clMD.MolecularDynamics(nloc=32)

# Allocate memory for 1 system (nSystems=1) using the MMFF template
md.realloc(nSystems=1, mmff=mmff)

# Pack the MMFF data into GPU buffers for system index 0
md.pack_system(iSys=0, mmff=mmff)

# Upload all system data to the GPU
md.upload_all_systems()

# Set up kernels with their arguments
md.setup_kernels()
md.setup_run_ocl_opt()

# Run optimization for 100 iterations with force convergence of 1e-6
print("\nRunning OpenCL optimization...")
iter_done = md.run_ocl_opt(niter=100, Fconv=1e-6, nPerVFs=10)
print(f"OpenCL optimization completed in {iter_done} iterations")

# Download results from GPU
final_pos, final_forces = md.download_results()
print("\nFinal positions:")
print(final_pos[0, :5, :3])  # Print positions of first 5 atoms
print("\nFinal forces:")
print(final_forces[0, :5, :3])  # Print forces of first 5 atoms

# Skip CUDA implementation for now
print("\nSkipping CUDA implementation comparison.")