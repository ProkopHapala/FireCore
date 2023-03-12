import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import plotUtils   as plu
from pyBall import MMFFsp3     as mmff


#======== Body

#mmff.setVerbosity( verbosity=1, idebug=0 )
mmff.setVerbosity( verbosity=0, idebug=0 )
#mmff.init( xyz_name="data/pyridine"  ) 

#------ Short Initialization
#mmff.init( xyz_name="data/HCOOH"  ) 

#------ Long Initialization
mmff.initParams()
print("=========== HCOOH ")
mmff.buildMolecule_xyz( xyz_name="data/HCOOH"  )
mmff.makeFFs()

mmff.getBuffs() 
labels = [ mmff.getType(i) for i in range(mmff.natoms) ]   ;print( labels ); print( mmff.apos )
plu.plotAtoms( mmff.apos, labels=labels, axes=(0,2) )
plt.show()

print("=========== pyridine ")
mmff.clear()
#mmff.buildMolecule_xyz( xyz_name="data/HCOOH"   )
mmff.buildMolecule_xyz( xyz_name="data/pyridine"  )
mmff.makeFFs()
print("=========== HCOOH ")
mmff.clear()
#mmff.buildMolecule_xyz( xyz_name="data/HCOOH"   )
mmff.buildMolecule_xyz( xyz_name="data/pyridine"  )
mmff.makeFFs()

'''
mmff.getBuffs()   # this has only sense with   makeFFs()                         
mmff.eval()                                   
mmff.setTrjName("relax.xyz",1)
mmff.run(1000)                               
'''


#print( "FORCES:\n mmff.fapos:\n ", mmff.fapos )
#mmff.plot(bForce=True, Fscale=10.0 )
#plt.show()
#exit(0)

