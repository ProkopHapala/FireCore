import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall.OCL.splines import OCLSplines
from pyBall.tests import ocl_Bsplines as bsp

os.environ['PYOPENCL_CTX'] = '0'
ocl_splines = OCLSplines()

def test_fit_Bspline( fname="debug_VCoul.npy", path="./common_resources/NaCl_1x1_L2/" ):
    E_ref = np.load( path+fname )
    Gs,conv = ocl_splines.fit3D( E_ref=E_ref, nmaxiter=3000, dt=0.3, Ftol=1e-8, cdamp=0.95, bAlloc=True, nPerStep=10, bConvTrj=True )
    print( "Gs.min(),Gs.max() ", Gs.min(),Gs.max() )

    if conv is not None:
        plt.plot( conv[:,0], conv[:,1], label="|F|_max " )
        plt.plot( conv[:,0], conv[:,2], label="|E|_max " )
        plt.yscale('log')
        plt.grid()
        plt.legend()


#test_eval_1D( )

# bsp.test_eval_3D( p0=[0.0,0.0,2.0], Q=0.4, E0=0 )
# bsp.test_eval_3D( p0=[2.0,2.0,2.0], Q=0.4, E0=0 )
# bsp.test_eval_3D( p0=[4.0,4.0,2.0], Q=0.4, E0=0 )


test_fit_Bspline(  )





plt.show()