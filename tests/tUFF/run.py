
import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF as mmff

MYbSimple = False
MYbConj   = False

#======== Body
mmff.setVerbosity( verbosity=1, idebug=1 )
#mmff.setVerbosity( verbosity=2, idebug=1 )

# MolWorld_sp3::init
mmff.init( xyz_name=str(sys.argv[1]),
           surf_name =None, 
           smile_name=None, 
           sElementTypes  = "data_UFF/ElementTypes.dat",
           sAtomTypes     = "data_UFF/AtomTypes.dat", 
           sBondTypes     = "data_UFF/BondTypes.dat", 
           sAngleTypes    = "data_UFF/AngleTypes.dat",
           sDihedralTypes = "data_UFF/DihedralTypes.dat",
           bMMFF=True,
           bEpairs=False,
           nPBC=(1,1,1),
           gridStep=0.1,
           bUFF=True,
           b141=True,
           bSimple=MYbSimple,
           bConj=MYbConj,
           bCumulene=True)

# MolWorld_sp3::eval
mmff.eval()
exit(0)

#mmff.relax(1000)
print( "FORCES:\n mmff.fapos:\n ", mmff.fapos )
mmff.plot(bForce=True, Fscale=10.0 )
plt.show()
exit(0)
