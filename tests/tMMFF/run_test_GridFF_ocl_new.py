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
#gff.test_gridFF_ocl( fname="data/xyz/NaCl_1x1_L2.xyz" )
gff.test_gridFF_ocl( fname="/home/prokop/git/FireCore/tests/pyutils/NaCl_8x8_L3.xyz", save_name="double3"  )
#gff.test_gridFF_ocl( fname="/home/prokop/git/FireCore/tests/pyutils/NaCl_8x8_L3.xyz" )

# ======== Ewald

d=0.6
apos=np.array([
    [-d,.0,0.],
    [+d,.0,0.],
    [0.,-d,0.],
    [0.,+d,0.],
])
qs = [ +1.,+1.,-1.,-1. ]

# d=0.6
# apos=np.array([
#     [0.,.0,-d],
#     [0.,.0,+d],
# ])
# qs = [ +1.,-1. ]

gff.test_Ewald( apos, qs,  ns=(100,100,100), dg=(0.10,0.10,0.10), order=3, bSlab=True,  nPBC=(100,100,0) )
# gff.test_Ewald( apos, qs,  ns=[100,100,150], dg=[0.10,0.10,0.10], order=3, bSlab=True,  nPBC=[100,100,0] )
# gff.test_Ewald( apos, qs,  ns=[100,100,200], dg=[0.10,0.10,0.10], order=3, bSlab=True,  nPBC=[100,100,0] )
# gff.test_Ewald( apos, qs,  ns=[100,100,300], dg=[0.10,0.10,0.10], order=3, bSlab=True,  nPBC=[100,100,0] )
# gff.test_Ewald( apos, qs,  ns=[100,100,400], dg=[0.10,0.10,0.10], order=3, bSlab=True,  nPBC=[100,100,0] )

plt.show()