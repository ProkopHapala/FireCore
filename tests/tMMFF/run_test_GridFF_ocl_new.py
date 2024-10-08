import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall.tests import ocl_GridFF_new as gff

# =============  Functions

R0 = 3.5
E0 = 0.1
a  = 1.8


Q = 0.4
p0 = [-2.0,-2.0,0.0]

#mol_name="PTCDA.xyz"
gff.test_gridFF_vocl( fname="data/xyz/NaCl_1x1_L2.xyz" )




plt.show()