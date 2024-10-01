import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall.OCL.splines import OCLSplines
#from pyBall import atomicUtils as au
#from pyBall import MMFF as mmff
#from pyBall import FunctionSampling as fu

# set the environment variable PYOPENCL_CTX='0' to avoid being asked again.
#sys.env['PYOPENCL_CTX'] = '0'
os.environ['PYOPENCL_CTX'] = '0'
ocl_splines = OCLSplines()

def getPLQH( R0, E0, a, Q, H ):
    e  = np.exp(a*R0);
    cL = e*E0;
    cP = e*cL;
    cH = e*e*H;
    return np.array([ cP, cL, Q, cH ])

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
    VPaul = np.load(path+"debug_BsplinePaul_pbc.npy")
    VLond = np.load(path+"debug_BsplineLond_pbc.npy")
    VCoul = np.load(path+"debug_BsplineCoul_pbc.npy")
    print("VCoul.shape = ", VCoul.shape)
    print("VLond.shape = ", VLond.shape)
    print("VPaul.shape = ", VPaul.shape)
    VPLQ = np.zeros(VCoul.shape+(4,), dtype=np.float32)
    VPLQ[:,:,:,0] = VPaul
    VPLQ[:,:,:,1] = VLond
    VPLQ[:,:,:,2] = VCoul
    VPLQ[:,:,:,3] = 0.0

    print("VPLQ.shape = ", VPLQ.shape)
    VPLQ = VPLQ.transpose( (2,1,0,3) ).copy()
    sh = VPLQ.shape; print("sh = ", sh)

    # xs = np.array( range( sh[0] ), dtype=np.float32 ); print("xs ", xs)
    # ys = np.array( range( sh[1] ), dtype=np.float32 ); print("ys ", ys)
    # zs = np.array( range( sh[2] ), dtype=np.float32 ); print("zs ", zs)
    # VPLQ[:,:,:,0] = xs[:,None,None]
    # VPLQ[:,:,:,1] = ys[None,:,None]
    # VPLQ[:,:,:,2] = zs[None,None,:]
    # VPLQ[:,:,:,3] = 0.0

    #PLQH = (1.0,1.0,1.0,0.0)
    PLQH = getPLQH( R0, E0, a, Q, H )

    g0   = (0.0,0.0,0.0)
    dg   = (0.2,0.2,0.2)

    zs = np.arange(0.0, 10.0, dsamp)
    ps = np.zeros( (len(zs), 4) )
    #ps[:,0] = 1.05
    #ps[:,1] = 1.05
    ps[:,0] = p0[0]
    ps[:,1] = p0[1]
    ps[:,2] = zs
    ps[:,3] = 0.0

    ocl_splines.prepare_sample3D( g0, dg, sh[:3], VPLQ )
    fe = ocl_splines.sample3D_comb(  ps, PLQH )

    #plt.figure(figsize=(10, 6))
    plt.plot(zs, fe[:, 3],'-k' ,lw=1.0, label='Interpolated')
    #plt.plot(xs, Gs     ,'.-r',lw=0.5, label='Original')
    plt.title('1D Spline Interpolation with PBC')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.legend()
    plt.grid()
   
    



# def test_eval_1D( g0=2.0, gmax=4.0, dg=0.2, dsamp=0.02, bUseForce=True, scErr=100.0, bHalf=False, title=None, order=3 ):
#     xs  = np.arange(g0, gmax+1e-8, dg)     ; ng=len(xs)
#     xs_ = np.arange(g0, gmax+1e-8, dsamp)  ; nsamp=len(xs_)
#     print("ng ", ng," nsamp ", nsamp)
    
#     Gs = np.zeros(ng)
#     Gs[ ng//2 ]=1.0

#     #FEout = mmff.sample_Bspline( xs_, Gs, x0=g0, dx=dg, order=3 )
#     FEout = clsp.sample_Bspline( xs_, Gs, x0=g0, dx=dg, order=order )
#     Es = FEout[:,0]
#     Fs = FEout[:,1]

#     Fnum = fu.numDeriv( Es, xs_ )

#     #plt.figure(figsize=(5,10))
#     plt.subplot(2,1,1); 
#     plt.plot( xs,  Gs, "o",           label="Gs poins" );
#     plt.plot( xs_, FEout[:,0], "-",  lw=0.5,  label="E_spline" );
#     plt.subplot(2,1,2); 
#     plt.plot( xs_, FEout[:,1], "-",  lw=0.5,  label="F spline" );
#     plt.plot( xs_[1:-1], Fnum, ":b",       lw=2.0,  label="F_num" );




#test_eval_1D( )

test_eval_3D( p0=[0.0,0.0,2.0], Q=0.4, E0=0 )
#test_eval_3D( p0=[2.0,2.0,2.0], Q=0.4, E0=0 )
test_eval_3D( p0=[4.0,4.0,2.0], Q=0.4, E0=0 )


plt.show()