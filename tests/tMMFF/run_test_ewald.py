import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches

sys.path.append("../../")
from pyBall import MMFF as mmff
from pyBall.tests import Ewald as ew

# d=0.6
# apos=np.array([
#     [-d,.0,0.],
#     [+d,.0,0.],
#     [0.,-d,0.],
#     [0.,+d,0.],
# ])
# qs = [ +1.,+1.,-1.,-1. ]

d=0.6
apos=np.array([
    [0.,.0,-d],
    [0.,.0,+d],
])
qs = [ +1.,-1. ]

# d=0.6
# apos = []
# qs   = [] 
# for ix in range(0,10):
#     for iy in range(0,10):
#         apos.append( [(ix-5)+0.5, (iy-5)+0.5, +d] ); qs.append(+0.01)
#         apos.append( [(ix-5)+0.5, (iy-5)+0.5, -d] ); qs.append(-0.01)


# ----- testing slab dipole correction
ew.test_vs_direct( apos, qs,  ns=[100,100,1000], dg=[0.10,0.10,0.10], order=3, bSlab=True,  nPBC=[200,200,0] )
# ew.test_vs_direct( apos, qs,  ns=[100,100,150], dg=[0.10,0.10,0.10], order=3, bSlab=True,  nPBC=[100,100,0] )
# ew.test_vs_direct( apos, qs,  ns=[100,100,200], dg=[0.10,0.10,0.10], order=3, bSlab=True,  nPBC=[100,100,0] )
# ew.test_vs_direct( apos, qs,  ns=[100,100,300], dg=[0.10,0.10,0.10], order=3, bSlab=True,  nPBC=[100,100,0] )
# ew.test_vs_direct( apos, qs,  ns=[100,100,400], dg=[0.10,0.10,0.10], order=3, bSlab=True,  nPBC=[100,100,0] )
# ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.10,0.10], order=3, bSlab=False, nPBC=[200,200,0] )

#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.10,0.10], order=3, bSlab=True,  nPBC=[100,100,0] )
#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.10,0.10], order=3, bSlab=False, nPBC=[100,100,0] )

#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.10,0.10], order=3, bPython=True )  
#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.10,0.10], order=3, bPython=True, pos0=[0,0,-5.0] )  
#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.10,0.10], order=2, bPython=True, pos0=[0,0,-5.0] )  

#ew.test_vs_direct( apos, qs,  ns=[16,16,16], dg=[0.10,0.10,0.10], order=2, bPython=True, pos0=[-0.8,-0.8,-0.8], bPlot1D=False )
#ew.test_vs_direct( apos, qs,  ns=[16,16,16], dg=[0.10,0.10,0.10], order=2, bPython=True, pos0=[ 0.0, 0.0,-0.8], bPlot1D=False )  

#ew.test_project_dens( apos, qs, ns=[16,16,16], pos0=[0.0,0.0,0.0],     order=2 )
#ew.test_project_dens( apos, qs, ns=[16,16,16], pos0=[-0.8,-0.8,-0.8], order=2 )
#ew.test_project_dens( apos, qs, ns=[16,16,16], pos0=[ 0.0, 0.0,-0.8], order=2 )

#ew.test_project2D( apos[:,:2].copy(), g0=[ 0.0, 0.0], order=3, ws=qs )
#ew.test_project2D( apos[:,:2].copy(), g0=[-0.8,-0.8], order=3, ws=qs )


# --- change voxel size  homogeneously in all directions
#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.15,0.15,0.15] )   # GOOD, This is perfect
#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.07,0.07,0.07] )   # GOOD, This is perfect

# --- change number of grid points homogeneously in all directions
#ew.test_vs_direct( apos, qs,  ns=[150,150,150], dg=[0.10,0.10,0.10] )    # GOOD, This is perfect
#ew.test_vs_direct( apos, qs,  ns=[80,80,80],    dg=[0.10,0.10,0.10] )    # GOOD, This is perfect

# ========= Changing step size in two directions

# --- change voxel size  homogeneously in xy-directions
#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.15,0.15,0.10] )    # GOOD, This is perfect
#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.07,0.07,0.10] )    # GOOD, This is perfect

# --- change voxel size  homogeneously in xz-directions
#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.07,0.10,0.07] )    # Quite good
#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.15,0.10,0.15] )    # Quite good

# --- change voxel size  homogeneously in yz-directions
#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.07,0.07] )    # Quite good
#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.15,0.15] )    # Quite good


# ========= Changing step size in one direction


# --- change voxel size  only in x-directions
#ew.test_vs_direct( apos, qs, ns=[100,100,100], dg=[0.07,0.10,0.10] )    # Quite good
#ew.test_vs_direct( apos, qs, ns=[100,100,100], dg=[0.15,0.10,0.10] )    # Quite good
#ew.test_vs_direct( apos, qs, ns=[100,100,100], dg=[0.20,0.10,0.10] )    # NOT SO PERFECT, iax0 goes to ~1.5 almost

#ew.test_vs_direct( apos, qs, ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[30,30,30] )    # NOT SO PERFECT, iax1 goes to ~1.5 almost
#ew.test_vs_direct( apos, qs, ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30] )    # NOT SO PERFECT, iax1 goes to ~1.5 almost

#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[30,30,30], order=3, bPython=True )    # NOT SO PERFECT, iax1 goes to ~1.5 almost

#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[100,100,100], order=3, bPython=True ) 

#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[30,30,30], order=3, bPython=True )    # Almost perfect
#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=3 )    # Almost perfect
#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[30,60,30], order=3 )    # Almost perfect
#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[30,30,60], order=3 )    # Almost perfect
#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.20,0.10,0.10], nPBC=[30,60,30], order=3 )    # Almost perfect

#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[120,60,60], order=3, bPython=True )    # Almost perfect

#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[30,30,30], order=3 )    # NOT SO PERFECT, iax1 goes to ~1.5 almost
#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=3 )    # Almost perfect


#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=3, nBlur=4, cV=0.85  )    # Almost perfect

#ew.test_vs_direct( apos, qs, ns=[100,100,100], dg=[0.10,0.20,0.10], order=3 )    # NOT SO PERFECT, iax1 goes to ~1.5 almost
#ew.test_vs_direct( apos, qs, ns=[100,100,100], dg=[0.10,0.10,0.20] )    # GOOD, This is perfect

# --- change voxel size  only in y-directions
#ew.test_vs_direct( apos, qs, ns=[100,100,100], dg=[0.10,0.07,0.10] )    # Quite good
#ew.test_vs_direct( apos, qs, ns=[100,100,100], dg=[0.10,0.15,0.10] )    # Quite good

# --- change voxel size  only in z-directions
#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.10,0.07] )    # GOOD, This is perfect
#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.10,0.15] )    # GOOD, This is perfect




# ========= Changing number of grid points in different directions

#ew.test_vs_direct( apos, qs,  ns=[150,100,100], dg=[0.10,0.10,0.10] )    # 
#ew.test_vs_direct( apos, qs,  ns=[100,150,100], dg=[0.10,0.10,0.10] )    # 
#ew.test_vs_direct( apos, qs,  ns=[100,100,150], dg=[0.10,0.10,0.10] )    # GOOD, This is perfect


#ew.test_vs_direct( apos, qs,  ns=[200,100,100], dg=[0.10,0.10,0.10] )    # NOT SO PERFECT - iax=0 goes up to 1.2
#ew.test_vs_direct( apos, qs,  ns=[100,200,100], dg=[0.10,0.10,0.10] )    # NOT SO PERFECT -  iax=1 goes up to 1.2


#ew.test_vs_direct( apos, qs,  ns=[50,100,100], dg=[0.10,0.10,0.10] )    # NOT SO PERFECT - iax=1 goes up to 1.2
#ew.test_vs_direct( apos, qs,  ns=[100,50,100], dg=[0.10,0.10,0.10] )    # NOT SO PERFECT -  iax=0 goes up to 1.2


#ew.test_vs_direct( apos, qs,  ns=[100,100,200], dg=[0.10,0.10,0.10] )   # GOOD, This is perfect
#ew.test_vs_direct( apos, qs,  ns=[100,100,50], dg=[0.10,0.10,0.10] )    # GOOD, This is perfect



# ========= Real space smoothening from Bspline(order=3)
#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=3, nBlur=2, cV=0.85, yrange=[0.98,1.02] )    # Almost perfect
#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=3, nBlur=2, cV=0.90, yrange=[0.98,1.02] )    # Almost perfect
#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=3, nBlur=2, cV=0.93, yrange=[0.98,1.02] )    # Almost perfect  # BEST
#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=3, nBlur=2, cV=0.95, yrange=[0.98,1.02] )    # Almost perfect
#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=3, nBlur=2, cV=0.97, yrange=[0.98,1.02] )    # Almost perfect
#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=3, nBlur=3, cV=0.85, yrange=[0.98,1.02] )    # Almost perfect
#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=3, nBlur=3, cV=0.90, yrange=[0.98,1.02] )    # Almost perfect
#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=3, nBlur=3, cV=0.94, yrange=[0.98,1.02] )    # Almost perfect
#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=3, nBlur=3, cV=0.95, yrange=[0.98,1.02] )    # Almost perfect  # BEST
#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=3, nBlur=3, cV=0.96, yrange=[0.98,1.02] )    # Almost perfect
#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=3, nBlur=4, cV=0.85, yrange=[0.98,1.02] )    # Almost perfect
#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=3, nBlur=4, cV=0.90, yrange=[0.98,1.02] )    # Almost perfect
#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=3, nBlur=4, cV=0.95, yrange=[0.98,1.02] )    # Almost perfect
#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=3, nBlur=4, cV=0.85, yrange=[0.99,1.01] )    # Almost perfect
#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=3, nBlur=4, cV=0.90, yrange=[0.99,1.01] )    # Almost perfect
#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=3, nBlur=4, cV=0.95, yrange=[0.99,1.01] )    # Almost perfect  # BEST

# ========= Real space smoothening from Bspline(order=2)

#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=2, nBlur=4, cV=0.850, yrange=[0.99,1.01] )   # Almost perfect
#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=2, nBlur=4, cV=0.900, yrange=[0.99,1.01] )   # Almost perfect
#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=2, nBlur=4, cV=0.930, yrange=[0.99,1.01] )   # Almost perfect
#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=2, nBlur=4, cV=0.935, yrange=[0.99,1.01] )   # Almost perfect  # BEST
#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=2, nBlur=4, cV=0.940, yrange=[0.99,1.01] )   # Almost perfect
#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=2, nBlur=4, cV=0.950, yrange=[0.99,1.01] )   # Almost perfect
#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=2, nBlur=4, cV=0.960, yrange=[0.99,1.01] )   # Almost perfect
#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.15,0.10], nPBC=[45,30,30], order=2, nBlur=4, cV=0.930, yrange=[0.99,1.01] )   # Almost perfect  # BEST
#ew.test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.15,0.10], nPBC=[45,30,30], order=2, nBlur=4, cV=0.935, yrange=[0.99,1.01] )   # Almost perfect  # BEST

plt.show()