
import sys
import numpy as np

sys.path.append("../../")
from pyBall.AtomicSystem import AtomicSystem
from pyBall.OCL import MMFF
from pyBall.OCL.MMFF import AtomType, Bond, Dihedral

from pyBall.OCL import cuMMFF            as cuMD
from pyBall.OCL import MolecularDynamics as clMD

mol = AtomicSystem( "common_resources/xyz/CH2NH.xyz" )

# Define AtomType instances with npi and ne
AtomTypeDict = {
    "C": AtomType(name="C", Kss=300.0, Asp=109.5, Ksp=100.0, Kpp=150.0, npi=0, ne=0),
    "H": AtomType(name="H", Kss=200.0, Asp=109.5, Ksp=50.0, Kpp=75.0, npi=0, ne=0),
    "O": AtomType(name="O", Kss=350.0, Asp=109.5, Ksp=120.0, Kpp=180.0, npi=0, ne=0),
    "N": AtomType(name="N", Kss=350.0, Asp=109.5, Ksp=120.0, Kpp=180.0, npi=0, ne=0),
    "E": AtomType(name="E", Kss=0.0, Asp=0.0, Ksp=0.0, Kpp=0.0, npi=0, ne=1)  # Electron pairs
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
mol.npi_list = np.array([0, 0, 0, 0,0], dtype=np.int32)
mol.nep_list = np.array([0, 0, 0, 0,0], dtype=np.int32)
mol.isNode   = np.array([1, 0, 0, 0,0], dtype=np.int32)
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


exit()


cuMD.init( mol.natoms, mol.natoms, mol.natoms, 0, 0 )