import sys
import numpy as np
import os

sys.path.append("../../")
from pyBall import MMFF as mmff

mmff.initWithMolFile( "C2H4.xyz", bNonBonded=False, bOptimizer=True)
#mmff.printBuffNames()
mmff.getBuffs() #;print( mmff.ndims )
mmff.eval()

print( "Es(Etot,Eb,Ea,Eps,EppT,EppI):", mmff.Es )