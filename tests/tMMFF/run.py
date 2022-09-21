import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import MMFF as mmff

mmff.initWithMolFile( "C2H4.xyz", bNonBonded=False, bOptimizer=True)
#mmff.printBuffNames()
mmff.getBuffs() #;print( mmff.ndims )
#mmff.eval()
#mmff.relax(1000, bWriteTrj=True )
Es=mmff.scanRotation( [1,4,5], 0, 0,1, np.pi*2, 100, bWriteTrj=True)   ;print("Es=", Es)

plt.plot(Es)
print( "Es(Etot,Eb,Ea,Eps,EppT,EppI):", mmff.Es )

plt.show()