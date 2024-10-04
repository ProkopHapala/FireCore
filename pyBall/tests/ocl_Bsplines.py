import numpy as np
import os
import matplotlib.pyplot as plt

from . import utils as ut
from ..OCL.splines import OCLSplines
#from .. import atomicUtils as au
#from .. import MMFF as mmff
#from .. import FunctionSampling as fu

# set the environment variable PYOPENCL_CTX='0' to avoid being asked again.
#sys.env['PYOPENCL_CTX'] = '0'
os.environ['PYOPENCL_CTX'] = '0'
ocl_splines = OCLSplines()

def test_eval_1D( g0=0.0, ng=10, dg=0.1, dsamp=0.02, bUseForce=True, scErr=100.0, bHalf=False, title=None, order=3 ):
    gmax = g0 + dg*ng
    xs = np.linspace(0.0, dg*ng, ng, endpoint=False ) + g0   ;print("xs = ", xs)
    Gs = np.sin(2*np.pi*xs).astype(np.float32)               ;print("Gs = ", Gs)
    ps = np.linspace( -gmax, gmax*2, 100, endpoint=False ).astype(np.float32)

    result_1d = ocl_splines.sample1D_pbc(g0, dg, ng, Gs, ps)

    plt.figure(figsize=(10, 6))
    plt.plot(ps, result_1d[:, 0],'-k' ,lw=1.0, label='Interpolated')
    plt.plot(xs, Gs             ,'.-r',lw=0.5, label='Original')
    plt.title('1D Spline Interpolation with PBC')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.grid()
    plt.legend()
    plt.show()

def test_eval_3D( path="./data/NaCl_1x1_L2/", dsamp=0.02, p0=[0.0,0.0,2.0], R0=3.5, E0=0.1, a=1.6, Q=0.4, H=0.0, scErr=100.0, Emax=None, Fmax=None, maxSc=5.0, title=None, bSaveFig=True ):
    VPLQ = ut.load_potential_comb( path )
    #PLQH = (1.0,1.0,1.0,0.0)
    PLQH = ut.getPLQH( R0, E0, a, Q, H )

    g0   = (0.0,0.0,0.0)
    dg   = (0.2,0.2,0.2)

    ps, ts = ut.make_sample_points_f4(p0, dsamp=dsamp)

    ocl_splines.prepare_sample3D( g0, dg, VPLQ.shape[:3], VPLQ )
    fe = ocl_splines.sample3D_comb(  ps, PLQH )

    #plt.figure(figsize=(10, 6))
    plt.plot( ts, fe[:, 3],'-k' ,lw=1.0, label='Interpolated')
    #plt.plot(xs, Gs     ,'.-r',lw=0.5, label='Original')
    plt.title('1D Spline Interpolation with PBC')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.legend()
    plt.grid()
   
    

if __name__ == "__main__":

    #test_eval_1D( )

    test_eval_3D( p0=[0.0,0.0,2.0], Q=0.4, E0=0 )
    #test_eval_3D( p0=[2.0,2.0,2.0], Q=0.4, E0=0 )
    test_eval_3D( p0=[4.0,4.0,2.0], Q=0.4, E0=0 )


    plt.show()