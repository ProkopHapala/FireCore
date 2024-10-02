import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import MMFF as mmff
#from pyBall.tests import GridFF as gff
from pyBall.tests import ocl_GridFF as gff

# =============  Functions

R0 = 3.5
E0 = 1.0
a  = 1.8

#Q = 0.0
Q = 0.4
#p0 = [1.0,-5.05,2.0]
#p0 = [0.0,0.0,2.0]
p0 = [-2.0,-2.0,0.0]

mmff.initParams()

Emax=0.00001 
Fmax=0.00001

gff.test_gridFF_vs_ocl( mode=6, title="Bspline_o3 \n(z-cut)" , p0=[0.0,0.0,2.0],  Q=0.4, E0=0, bRefine=False, nPBC=[400,400,0], Emax=Emax, Fmax=Fmax )
gff.test_gridFF_vs_ocl( mode=6, title="Bspline_o3 \n(z-cut)" , p0=[2.0,2.0,2.0],  Q=0.4, E0=0, bRefine=False, nPBC=[400,400,0], Emax=Emax, Fmax=Fmax )
plt.show()