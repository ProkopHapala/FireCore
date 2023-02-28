import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

sys.path.append("../../")
from pyBall import FitREQ as fit



# just fit charge
typeMask = np.array([ [0,0,1],])
typREQs  = np.array([ 1.0,0.0, +0.2 ])   # charge = 0.2
  
fit.init_types( typeMask, typREQs )
#fit.setRigidSamples( int n, Es, poses, bool bCopy=True, bAlloc=True )
