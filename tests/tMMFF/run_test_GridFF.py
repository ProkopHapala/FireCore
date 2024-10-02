import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import MMFF as mmff
from pyBall.tests import GridFF as gff

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

#gff.test_gridFF    ( mode=6, title="Bspline_o3 \n(z-cut)" ,    Q=0.0, )
#gff.test_gridFF_lat( mode=6, title="Bspline_o3 \n(lat iax=0)", Q=0.0, )
#gff.test_gridFF_lat( mode=6, title="Bspline_o3 \n(lat iax=0)", Q=0.0, p0=p0, iax=0 )

#gff.test_gridFF    ( mode=1, title="tri-linar force \n(z-cut)"          )
#gff.test_gridFF_lat( mode=1, title="tri-Linear Force", Q=0.0, p0=p0, iax=0 )

#Emax=0.01 
#Fmax=0.01

Emax=0.00001 
Fmax=0.00001

#gff.test_gridFF    ( mode=6, title="Bspline_o3 \n(z-cut)" , p0=[0.0,0.0,2.0],  Q=0.4, E0=0, bRefine=False, nPBC=[50 ,50 ,0], Emax=Emax, Fmax=Fmax )
#gff.test_gridFF    ( mode=6, title="Bspline_o3 \n(z-cut)" , p0=[0.0,0.0,2.0],  Q=0.4, E0=0, bRefine=False, nPBC=[100,100,0], Emax=Emax, Fmax=Fmax )
#gff.test_gridFF    ( mode=6, title="Bspline_o3 \n(z-cut)" , p0=[0.0,0.0,2.0],  Q=0.4, E0=0, bRefine=False, nPBC=[150,150,0], Emax=Emax, Fmax=Fmax )

#gff.test_gridFF    ( mode=6, title="Bspline_o3 \n(z-cut)" , p0=[0.0,0.0,2.0],  Q=0.4, E0=0, bRefine=False, nPBC=[100,100,0], Emax=Emax, Fmax=Fmax )
#gff.test_gridFF    ( mode=6, title="Bspline_o3 \n(z-cut)" , p0=[0.0,0.0,2.0],  Q=0.4, E0=0, bRefine=False, nPBC=[300,300,0], Emax=Emax, Fmax=Fmax )
gff.test_gridFF    ( mode=6, title="Bspline_o3 \n(z-cut)" , p0=[0.0,0.0,2.0],  Q=0.4, E0=0, bRefine=False, nPBC=[400,400,0], Emax=Emax, Fmax=Fmax )
gff.test_gridFF    ( mode=6, title="Bspline_o3 \n(z-cut)" , p0=[2.0,2.0,2.0],  Q=0.4, E0=0, bRefine=False, nPBC=[400,400,0], Emax=Emax, Fmax=Fmax )

Emax=0.01 
Fmax=0.01

gff.test_gridFF_lat( mode=6, title="Bspline_o3 \n(lat iax=0)", p0=[0.0,0.0,2.0], nPBC=[400,400,0],  Q=0.4, E0=0.0, Emax=Emax, Fmax=Fmax  )
gff.test_gridFF_lat( mode=6, title="Bspline_o3 \n(lat iax=0)", p0=[2.0,2.0,2.0], nPBC=[400,400,0],  Q=0.4, E0=0.0, Emax=Emax, Fmax=Fmax  )

Emax=0.1 
Fmax=0.1

#gff.test_gridFF    ( mode=6, title="Bspline_o3 \n(z-cut)" , p0=[0.0,0.0,2.0],  Q=0.0, E0=0.1, bRefine=False, nPBC=[5,5,0], Emax=Emax, Fmax=Fmax )
#gff.test_gridFF    ( mode=6, title="Bspline_o3 \n(z-cut)" , p0=[2.0,2.0,2.0],  Q=0.0, E0=0.1, bRefine=False, nPBC=[5,5,0], Emax=Emax, Fmax=Fmax )
#gff.test_gridFF_lat( mode=6, title="Bspline_o3 \n(lat iax=0)", p0=[0.0,0.0,2.0], Q=0.0, E0=0.1 )
#gff.test_gridFF_lat( mode=6, title="Bspline_o3 \n(lat iax=0)", p0=[2.0,2.0,2.0], Q=0.0, E0=0.1 )

#gff.test_gridFF_npy( ps_xy=[(0.0,0.0),(0.0,0.5),(0.5,0.0),(0.5,0.5)], mode=6, title="" )

#gff.test_gridFF_npy( ps_xy=[(0.0,0.0)], mode=6, title="" )
#gff.test_gridFF_npy_lat( ps_zy=[(0.0,0.0)], mode=6, title="" )
#gff.test_gridFF_npy_lat( ps_zy=[(0.1,0.1)], mode=6, title="" )

plt.show()