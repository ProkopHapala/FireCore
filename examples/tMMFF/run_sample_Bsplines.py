import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
#from pyBall import atomicUtils as au
#from pyBall import MMFF as mmff
#from pyBall import FunctionSampling as fu

from pyBall import MMFF as mmff
from pyBall.tests import Bsplines as bsp

#mmff.setVerbosity( 2 )
mmff.setVerbosity( 3 )

#test_PBCindexes( order=3 )
#test_PBCindexes( order=5 )
#test_project1D( [ 0.050, 1.050] , g0=0.0, dg=0.1, ng=20, order=3 )
#test_project1D( [ 0.025, 1.025] , g0=0.0, dg=0.1, ng=20, order=3 )
# test_project1D( [ 0.000, 1.000] , g0=0.0, dg=0.1, ng=20, order=3 )
# test_project1D( [-0.025, 0.975] , g0=0.0, dg=0.1, ng=20, order=3 )
# test_project1D( [-0.050, 0.950] , g0=0.0, dg=0.1, ng=20, order=3 )
#test_project1D( [0.05, 1.0] , g0=0.0, dg=0.1, ng=20, order=5 )
#plt.legend()

d=0.6
apos=np.array([
    [-d,.0],
    [+d,.0],
    [0.,-d],
    [0.,+d],
])
qs = [ +1.,+1.,-1.,-1. ]

# test_project2D( apos, g0=[ 0.0, 0.0], order=3, ws=qs )
# test_project2D( apos, g0=[-0.8,-0.8], order=3, ws=qs )


#plt.figure(figsize=(5,10))
#test_eval_1D(order=3)
#test_eval_1D(order=5)

#test_NURBS( g0=0.0, dg=0.5, dsamp=0.05 )

#test_fit_1D( bUseForce=True )
#test_fit_1D( bUseForce=False, bHalf=False ,title="No-Half")
#test_fit_1D( bUseForce=False, bHalf=True  ,title="Half")
#test_fit_1D( g0=0.0, gmax=2.0, dg=0.1, bUseForce=False )

#test_fit_1D( g0=0.0, ng=10, dg=0.25 )
bsp.test_fit_1D( g0=0.0, ng=8, dg=0.25 )




#test_fit_2D(  )
#test_fit_2D( g0=(-1.0,-1.0), gmax=(1.0,1.0) )

#test_fit_2D_debug( title="2D fit run debug" )
#test_fit_3D_debug( title="3D fit run debug" )

#test_comb3_2D()

# Gs = test_comb3_3D(   iax=0, title="x-cut" );
# test_comb3_3D( Gs=Gs, iax=1, title="y-cut" );
# test_comb3_3D( Gs=Gs, iax=2, title="z-cut" );
#test_comb3_3D( iax=1 );
#test_comb3_3D( iax=2 );

#test_fit_3D(  )


plt.show()