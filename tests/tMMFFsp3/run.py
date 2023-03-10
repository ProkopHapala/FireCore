import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFFsp3     as mmff


#======== Body

#mmff.setVerbosity( verbosity=1, idebug=0 )
mmff.setVerbosity( verbosity=0, idebug=0 )
#mmff.init( xyz_name="data/pyridine"  ) 
mmff.init( xyz_name="data/HCOOH"  ) 

'''
mmff.getBuffs()                            
mmff.eval()                                   
mmff.setTrjName("relax.xyz",1)
mmff.run(1000)                               
'''


#print( "FORCES:\n mmff.fapos:\n ", mmff.fapos )
#mmff.plot(bForce=True, Fscale=10.0 )
#plt.show()
#exit(0)

