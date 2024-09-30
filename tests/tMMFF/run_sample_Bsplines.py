import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF as mmff
from pyBall import FunctionSampling as fu

# =============

apos = [
    [-2.0,0.0,0.0],
    [ 2.0,0.0,0.0],
]
REs=[
    [3.5,1.0],
    [3.5,1.0],
]

# =============

### ==== ToDo: Extract these functions to GridUtils.py (?)



def test_NURBS( g0=0.0, dg=0.5, dsamp=0.05 ):
    Gs  = [ 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, -0.5, 0.0, 0.0, 0.0  ]
    Ws  = [ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0  ]
    Gs  = np.array(Gs); Ws=np.array(Ws)
    ng=len(Gs); gmax=ng*dg
    xgs = np.arange(g0, gmax, dg)     ;
    xs  = np.arange(g0, gmax+1e-8,  dsamp)  ; nsamp=len(xs)
    #Ws[3] = 10.1
    #Ws[4] = 10.1
    #Ws[5] = 10.1
    Ws[6] = 100.0
    #Ws[6] = .01

    #plt.figure(figsize=(5,5))
    #plt.plot( xgs, Gs, ":ok",  label="Gs" )
    #for iw,w in enumerate([ 0.2, 0.5, 1.0, 1.5, 2.0 ]):
    #    Ws[3] = w
    #    FEout = mmff.sample_NURBS( xs, Gs, Ws, x0=g0, dx=dg )
    #    #plt.subplot(2,1,1)    
    #    plt.plot( xs, FEout[:,0], "-",  lw=0.5,  label=("E_fit(w=%2.2f)" %w) )

    FEout = mmff.sample_NURBS( xs, Gs, Ws, x0=g0, dx=dg )

    plt.figure(figsize=(5,10))
    plt.subplot(2,1,1) 
    plt.plot( xgs, Gs, ":ok",  label="Gs" )
    plt.plot( xs, FEout[:,0], "-",  lw=0.5,  label="E_fit" )
    #plt.ylim(Emin*1.2,-Emin*1.2)
    plt.grid()
    plt.legend()
    #plt.title("Energy")
    
    Fnum = fu.numDeriv( FEout[:,0], xs ) #;print("Fnum ", Fnum)
    plt.subplot(2,1,2) 
    plt.plot( xs, -FEout[:,1], "-b", lw=0.5, label="F_fit" )
    plt.plot( xs[1:-1], -Fnum, ":b", lw=2.0, label="F_num" )
    #plt.ylim(Fmin*1.2,-Fmin*1.2)
    plt.legend()
    #plt.title("Force")
    plt.grid()

def test_eval_1D( g0=2.0, gmax=4.0, dg=0.2, dsamp=0.02, bUseForce=True, scErr=100.0, bHalf=False, title=None, order=3 ):
    xs  = np.arange(g0, gmax+1e-8, dg)     ; ng=len(xs)
    xs_ = np.arange(g0, gmax+1e-8, dsamp)  ; nsamp=len(xs_)
    print("ng ", ng," nsamp ", nsamp)
    
    Gs = np.zeros(ng)

    Gs[ ng//2 ]=1.0

    #FEout = mmff.sample_Bspline( xs_, Gs, x0=g0, dx=dg, order=3 )
    FEout = mmff.sample_Bspline( xs_, Gs, x0=g0, dx=dg, order=order )
    Es = FEout[:,0]
    Fs = FEout[:,1]

    Fnum = fu.numDeriv( Es, xs_ )

    #plt.figure(figsize=(5,10))
    plt.subplot(2,1,1); 
    plt.plot( xs,  Gs, "o",           label="Gs poins" );
    plt.plot( xs_, FEout[:,0], "-",  lw=0.5,  label="E_spline" );
    plt.subplot(2,1,2); 
    plt.plot( xs_, FEout[:,1], "-",  lw=0.5,  label="F spline" );
    plt.plot( xs_[1:-1], Fnum, ":b",       lw=2.0,  label="F_num" );



def test_fit_1D( g0=-2.0, ng=10, dg=0.2, dsamp=0.02, bUseForce=False, scErr=20.0, bHalf=False, bPBC=True, Kreg=0.01, title=None ):
    #x0 = 2.0
    #dx = 0.1
    gmax=g0+ng*dg 
    Dg = gmax-g0
    if bHalf: dg = dg*0.5
    xs  = np.arange(g0,    gmax,      dg    ); ng=len(xs)
    xs_ = np.arange(g0,    gmax+1e-8, dsamp ); nsamp=len(xs_)
    #xs_ = np.arange(g0-Dg, gmax+Dg, dsamp ); nsamp=len(xs_)
    xsg=xs; 
    if bHalf: xsg=xs[::2]
    
    # E,F         = fu.getLJ( xs,  3.5, 1.0 )
    # E_ref,F_ref = fu.getLJ( xs_, 3.5, 1.0 )
    
    E,F         = fu.getCos( xs,  np.pi )
    E_ref,F_ref = fu.getCos( xs_, np.pi )

    #print( "E_ref ", E_ref )
    #print( "F_ref ", F_ref )

    Emin =  E.min()
    Fmin = -F.max()
    #print( "Emin ", Emin," Fmin ", Fmin )
    Ecut = -2.0
    #Ws = 1/np.sqrt( (E/Ecut)**2 + 1 ); # EWs = E*Ws
    #Ws = 1/( E - Ecut );
    #E*=Ws
    #FEg[:,1]*=-1

    #Ecut2 =  1.0
    #mask    = E>Ecut2
    #E[mask] = Ecut2

    dgs = dg; 
    if bHalf: dgs=dg*2

    Gs = E.copy()
    
    #FEout0 = mmff.sample_Bspline( xs_, Gs, x0=g0, dx=dgs )

    if bUseForce:
        FEg = np.zeros( (len(xs),2) )
        FEg[:,0] = E[:]
        FEg[:,1] = F[:]
        Gs = E*0.9
        Ws     = np.zeros((ng,2));  Ws[:,0]=1.0; Ws[:,1]=0.0;
        Gs, Ws = mmff.fitEF_Bspline( dg, FEg, Gs=Gs, Ws=Ws, Ftol=1e-6, nmaxiter=1000, dt=0.1 )
    else:
        Gs, Ws = mmff.fit_Bspline( E, Ws=None, dt=0.4, nmaxiter=200, Ftol=1e-12, Kreg=Kreg, bHalf=bHalf, bPBC=bPBC )
        #Gs, Ws_ = mmff.fit_Bspline( E, Ws=Ws,   dt=1.0, nmaxiter=1000, Ftol=1e-9, bHalf=bHalf )
        #Gs, Ws = mmff.fit_Bspline( FEg[:,0].copy(), Ws=Ws,   dt=0.4, nmaxiter=1000, Ftol=1e-7 )

        dgb_vG = mmff.getArrayPointer("dgb_vG"); print(dgb_vG.shape)
        dgb_dF = mmff.getArrayPointer("dgb_dF"); print(dgb_dF.shape)


    FEout = mmff.sample_Bspline( xs_, Gs, x0=g0, dx=dgs )
    #FEout = mmff.sample_NURBS( xs_, Gs, Ws=np.ones(len(Gs)), x0=g0, dx=dgs )

    plt.figure(figsize=(10,15))
    plt.subplot(3,1,1)    
    plt.plot( xs, Ws*10.0, "-c", lw=1.0, label="Ws" )
    #plt.plot( xs, EWs, ".-k" )
    plt.plot( xsg, Gs,         ".-m", lw=0.25, label="Gs" )
    plt.plot( xs_, E_ref,      ":k",  lw=1.5,  label="E_ref" )
    plt.plot( xs_, FEout[:,0], "-g",  lw=0.5,  label="E_fit" )
    #plt.plot( xs_, FEout0[:,0],"--k",  lw=0.5,  label="E_fit0" )
    plt.plot( xs_, (FEout[:,0]-E_ref)*scErr, "-r", lw=0.5, label=("error*%g" % scErr) )
    #print( "Gs: ", Gs )
    plt.ylim(Emin*1.2,-Emin*1.2)
    plt.xlim(g0,gmax*1.5)
    plt.grid()
    plt.legend()
    plt.title("Energy")
    plt.subplot(3,1,2)
    #plt.plot( xs,  -F ,       ".-k", lw=0.5, label="F_ref" )    
    plt.plot( xs_, -F_ref ,     ":k", lw=1.5, label="F_ref" )    
    plt.plot( xs_, -FEout[:,1], "-g", lw=0.5, label="F_fit" )
    plt.plot( xs_, (FEout[:,1]-F_ref)*scErr, "-r", lw=0.5, label=("error*%g" % scErr) )
    plt.ylim(Fmin*1.2,-Fmin*1.2)
    plt.xlim(g0,gmax*1.5)
    plt.legend()
    plt.title("Force")
    plt.grid()

    plt.subplot(3,1,3)
    plt.plot( xs, dgb_vG[:,0],      ".-r",  lw=1.0,  label="dgb_vG[0]" )
    plt.plot( xs, dgb_vG[:,1],      ".-g",  lw=1.0,  label="dgb_vG[1]" )
    plt.plot( xs, dgb_vG[:,2],      ".-c",  lw=1.0,  label="dgb_vG[2]" )
    plt.plot( xs, dgb_vG[:,3],      ".-b",  lw=1.0,  label="dgb_vG[3]" )
    # plt.plot( xs, dgb_dF[:,0],      ":r",  lw=1.0,  label="dgb_dF[0]" )
    # plt.plot( xs, dgb_dF[:,1],      ":g",  lw=1.0,  label="dgb_dF[1]" )
    # plt.plot( xs, dgb_dF[:,2],      ":c",  lw=1.0,  label="dgb_dF[2]" )
    # plt.plot( xs, dgb_dF[:,3],      ":b",  lw=1.0,  label="dgb_dF[3]" )
    plt.xlim(g0,gmax*1.5)
    plt.legend()
    plt.title("Force")
    plt.title("Regularization")
    plt.grid()
    if title is not None: plt.suptitle(title)
    #plt.show()


def test_fit_2D( g0=(-5.0,2.0), gmax=(5.0,10.0), dg=(0.1,0.1), dsamp=(0.05,0.05), title=None, bPBC=True ):
    #cmap="RdBu_r"
    cmap="bwr"
    #x0 = 2.0
    #dx = 0.1
    Xs,Ys   = fu.make2Dsampling(  g0=g0, gmax=gmax, dg=dg )
    Xs_,Ys_ = fu.make2Dsampling(  g0=g0, gmax=gmax, dg=dsamp )

    Xs_*=0.999999; Ys_*=0.999999;

    sh_samp = Xs_.shape
    ps      = fu.pack_ps2D( Xs_, Ys_)
    #ps, sh_samp = make2Dsampling_ps(  g0=g0, gmax=gmax, dg=dsamp )

    print( "Xs.shape ", Xs.shape )
    
    #E, Fx,Fy,Fz = getLJ_atoms( apos, REs, Xs,Ys,Xs*0.0 )

    E,  Fx,Fy      =  fu.getCos2D( Xs, Ys   )
    E_r, Fx_r,Fy_r =  fu.getCos2D( Xs_, Ys_ )

    #xs_ = np.arange(g0, gmax, dsamp)  ; nsamp=len(xs_)
    Emin =  E.min()
    Fmin = -Fy.max()
    print( "Emin ", Emin," Fmin ", Fmin )
    #FEg = np.zeros( (len(xs),2) )
    #FEg[:,0] = E[:]
    #FEg[:,1] = F[:]
    Ecut = 100.0

    Gs, Ws = mmff.fit2D_Bspline( E, Ws=None, dt=0.4, nmaxiter=1000, Ftol=1e-7, bPBC=bPBC )
    E_f = mmff.sample_Bspline2D( ps, Gs, g0, dg, fes=None  ).reshape(sh_samp+(3,))

    Gmin = -np.abs(Gs).max()
    dG = Gs-E                    ; dGmin = -np.abs(dG).max()
    dE = E_f[:,:,2] - E_r[:,:]   ; dEmin = -np.abs(dE[4:-4,4:-4]).max()
    #dE = E_f[1:,1:,2] - E_r[:-1,:-1]   ; dEmin = -np.abs(dE).max()
    #dE = E_f[:-1,:-1,2] - E_r[1:,1:]   ; dEmin = -np.abs(dE).max()
    #dE = E_f[:-2,:-2,2] - E_r[2:,2:]   ; dEmin = -np.abs(dE).max()

    extent=(g0[0],gmax[0],g0[1],gmax[1])
    plt.figure(figsize=(15,10))
    plt.subplot(2,3,1); plt.imshow( E,         origin="lower", extent=extent, vmin=Emin,  vmax=-Emin,  cmap=cmap ) ;plt.colorbar(); plt.title("Eg")
    plt.subplot(2,3,2); plt.imshow( Gs,        origin="lower", extent=extent, vmin=Gmin,  vmax=-Gmin,  cmap=cmap ) ;plt.colorbar(); plt.title("Gs fit")
    plt.subplot(2,3,3); plt.imshow( dG,        origin="lower", extent=extent, vmin=dGmin, vmax=-dGmin, cmap=cmap ) ;plt.colorbar(); plt.title("Gs-Eg")
    
    plt.subplot(2,3,4); plt.imshow( E_r,        origin="lower", extent=extent, vmin=Emin, vmax=-Emin,   cmap=cmap ) ;plt.colorbar(); plt.title("E  ref")
    plt.subplot(2,3,5); plt.imshow( E_f[:,:,2], origin="lower", extent=extent, vmin=Emin, vmax=-Emin,   cmap=cmap ) ;plt.colorbar(); plt.title("E  fit")
    plt.subplot(2,3,6); plt.imshow( dE,         origin="lower", extent=extent, vmin=dEmin, vmax=-dEmin, cmap=cmap ) ;plt.colorbar(); plt.title("E(fit-ref)")
    plt.axis('equal')


def test_comb3_2D( g0=(-2.0,-2.0), gmax=(0.0,0.0), dg=(0.1,0.1), dsamp=(0.05,0.05), title=None, scErr=1000.0, bPBC=True, Ccomb=[1.0,0.0,0.0], bFit=True ):
    print( "test_comb3_2D START" )
    cmap="bwr"
    Xs,Ys   = fu.make2Dsampling(  g0=g0, gmax=gmax, dg=dg )
    Xs_,Ys_ = fu.make2Dsampling(  g0=g0, gmax=gmax, dg=dsamp )
    Xs_*=0.999999; Ys_*=0.999999;
    sh_samp = Xs_.shape
    ps      = fu.pack_ps2D( Xs_, Ys_)

    print( "Xs.shape ", Xs.shape )
    
    #E, Fx,Fy,Fz = getLJ_atoms( apos, REs, Xs,Ys,Xs*0.0 )

    E1,  Fx,Fy      =  fu.getCos2D( Xs,    Ys   )
    E2,  Fx,Fy      =  fu.getCos2D( Xs*2,  Ys   )
    E3,  Fx,Fy      =  fu.getCos2D( Xs,  Ys*2   )

    E1_,  Fx,Fy     =  fu.getCos2D( Xs_,    Ys_   )
    E2_,  Fx,Fy     =  fu.getCos2D( Xs_*2,  Ys_   )
    E3_,  Fx,Fy     =  fu.getCos2D( Xs_,    Ys_*2 )

    #E_r, Fx_r,Fy_r =  fu.getCos2D( Xs_, Ys_ )

    Gs = np.zeros( E1.shape+(3,) )

    if bFit:
        G1, Ws = mmff.fit2D_Bspline( E1, dt=0.4, nmaxiter=1000, Ftol=1e-10, bPBC=bPBC )
        G2, Ws = mmff.fit2D_Bspline( E2, dt=0.4, nmaxiter=1000, Ftol=1e-10, bPBC=bPBC )
        G3, Ws = mmff.fit2D_Bspline( E3, dt=0.4, nmaxiter=1000, Ftol=1e-10, bPBC=bPBC )
        Gs[:,:,0] = G1
        Gs[:,:,1] = G2
        Gs[:,:,2] = G3
    else:
        Gs[:,:,0] = E1
        Gs[:,:,1] = E2
        Gs[:,:,2] = E3

    #Gs[:] = np.roll(Gs, shift=-1, axis=0).copy()

    #ps[:,0] += -0.00003
    #ps[:,1] += -0.005375

    E_f = mmff.sample_Bspline2D_comb3( ps, Gs, g0, dg,  Cs=Ccomb  ).reshape(sh_samp+(3,))

    Eref = E1_*Ccomb[0] + E2_*Ccomb[1] + E3_*Ccomb[2]

    Err = (E_f[:,:,2]-Eref)

    Wmax=Ws.max()
    Emin=Eref.min()
    #Emin=E_f.min()
    dEmax = np.abs(Err[:,4:-8]).max()
    extent=(g0[0],gmax[0],g0[1],gmax[1])
    plt.figure(figsize=(20,10))
    plt.subplot(2,4,4); plt.imshow( Ws [:,:],   origin="lower", extent=extent, vmin=-Wmax,  vmax=Wmax,  cmap=cmap ) ;plt.colorbar(); plt.title("Ws"    )
    plt.subplot(2,4,1); plt.imshow( E_f[:,:,2], origin="lower", extent=extent, vmin=Emin,  vmax=-Emin,  cmap=cmap ) ;plt.colorbar(); plt.title("E_fit" )
    plt.subplot(2,4,2); plt.imshow( Eref,       origin="lower", extent=extent, vmin=Emin,  vmax=-Emin,  cmap=cmap ) ;plt.colorbar(); plt.title("E_ref" )
    plt.subplot(2,4,3); plt.imshow( Err,        origin="lower", extent=extent, vmin=-dEmax,vmax=dEmax,  cmap=cmap ) ;plt.colorbar(); plt.title("Error" ) #plt.title("Error *%i" %scErr)
    plt.axis('equal')

    ix0=Eref.shape[0]//2
    iy0=Eref.shape[1]//2
    plt.subplot(2,2,3); 
    plt.plot( Eref[:,iy0  ],                    'g-',lw=0.5, label="E_ref(x)" ); 
    plt.plot( E_f [:,iy0,2],                    'k-',lw=0.5, label="E_fit(x)" ); 
    plt.plot((E_f[:,iy0,2]-Eref[:,iy0])*scErr,  'r-',lw=1.0, label=("Err(x)*%g" %scErr) ); 
    plt.ylim(Emin,-Emin)
    plt.grid()
    plt.legend()

    plt.subplot(2,2,4); 
    plt.plot( Eref[ix0,:  ],                   'g-' ,lw=0.5, label="E_ref(y)" ); 
    plt.plot( E_f[ix0,:,2],                    'k-' ,lw=0.5, label="E_fit(y)" );
    plt.plot((E_f[ix0,:,2]-Eref[ix0,:])*scErr, 'r-' ,lw=1.0, label=("Err(y)*%g" %scErr) ); 
    plt.ylim(Emin,-Emin)
    plt.grid()
    plt.legend()
    print( "test_comb3_2D DONE" )

def test_fit_2D_debug( g0=(-2.0,-2.0), gmax=(2.0,2.0), dg=(0.2,0.2), title=None, bPBC=True ):
    cmap="bwr"
    Xs,Ys   = fu.make2Dsampling(  g0=g0, gmax=gmax, dg=dg )

    E,  Fx,Fy      =  fu.getCos2D( Xs, Ys   )
    #E,  Fx,Fy = fu.getGauss2D( Xs,  Ys-4.0, w=0.5  )
    Ws = np.zeros( E.shape )

    nplt = 5
    nitr = 10

    plt.figure(figsize=(5*nplt,5*3))
    Gs = None
    for i in range(nplt):
        itr = i*nitr
        Gs, Ws = mmff.fit2D_Bspline( E, Ws=Ws, Gs=Gs, dt=0.4, nmaxiter=nitr, Ftol=1e-12, bPBC=bPBC )
        plt.subplot(3,nplt,i       +1); plt.imshow( Ws[:,:],  origin="lower"             ) ;plt.colorbar(); plt.title("Ws[iter=%i]" %itr )
        plt.subplot(3,nplt,nplt  +i+1); plt.imshow( Gs[:,:],  origin="lower"             ) ;plt.colorbar(); plt.title("Gs[iter=%i]" %itr )
        plt.subplot(3,nplt,nplt*2+i+1); plt.imshow( E [:,:],  origin="lower", cmap='bwr' ) ;plt.colorbar(); plt.title("Gs[iter=%i]" %itr )
    if title is not None: plt.suptitle(title)

def test_fit_3D_debug( g0=(-2.0,-2.0,-2.0), gmax=(2.0,2.0,2.0), dg=(0.2,0.2,0.2), title=None, bPBC=True  ):
    cmap="bwr"
    Xs,Ys,Zs      = fu.make3Dsampling(  g0=g0, gmax=gmax, dg=dg )

    E,  Fx,Fy,Fz  =  fu.getCos3D( Xs, Ys, Zs   )
    #E,  Fx,Fy,Fz        = fu.getGauss3D( Xs,  Ys , Zs-4.0, w=0.5  )
    Ws = np.zeros( E.shape )
    
    iz0=10
    nplt = 5
    nitr = 100

    plt.figure(figsize=(5*nplt,10))
    Gs = None
    for i in range(nplt):
        itr = i*nitr
        Gs = mmff.fit3D_Bspline( E, Ws=Ws, Gs=Gs, dt=0.4, nmaxiter=nitr, Ftol=1e-12, bPBC=bPBC )
        plt.subplot(2,nplt,i+1     ); plt.imshow( Ws[iz0,:,:],  origin="lower" ) ;plt.colorbar(); plt.title("Ws[iter=%i]" %itr )
        plt.subplot(2,nplt,nplt+i+1); plt.imshow( Gs[iz0,:,:],  origin="lower" ) ;plt.colorbar(); plt.title("Gs[iter=%i]" %itr )
    if title is not None: plt.suptitle(title)


def test_fit_3D( g0=(-2.0,-2.0,-2.0), gmax=(2.0,2.0,2.0), dg=(0.2,0.2,0.2), dsamp=(0.05,0.05,0.05), bPBC=True  ):
    cmap="bwr"
    Xs,Ys,Zs    = fu.make3Dsampling(  g0=g0, gmax=gmax, dg=dg )
    Xs_,Ys_,Zs_ = fu.make3Dsampling(  g0=g0, gmax=gmax, dg=dsamp )

    #Xs_*=0.999999; Ys_*=0.999999; Zs_*=0.999999;

    sh_samp = Xs_.shape
    ps      = fu.pack_ps3D(Xs_,Ys_,Ys_)
    #ps, sh_samp = make2Dsampling_ps(  g0=g0, gmax=gmax, dg=dsamp )
    print( "Xs.shape ", Xs.shape )
    #E, Fx,Fy,Fz = getLJ_atoms( apos, REs, Xs,Ys,Xs*0.0 )
    E,  Fx,Fy,Fz        = fu.getCos3D( Xs,  Ys , Zs  )
    #E_r, Fx_r,Fy_r,Fz_r = fu.getCos3D( Xs_, Ys_, Zs_ )

    #E,  Fx,Fy,Fz        = fu.getGauss3D( Xs,  Ys , Zs-4.0, w=0.5  )
    # E_r, Fx_r,Fy_r,Fz_r = fu.getCos3D( Xs_, Ys_, Zs_ )

    #Emin =  E_r.min()
    #Fmin = -Fy_r.max()

    ix0=10;iy0=10;iz0=50;
    #E[:,:,:]=0.0
    #E[iz0,iy0,ix0]=1.0

    #Gs = mmff.fit3D_Bspline( E, dt=0.1, nmaxiter=100, Ftol=1e-6,  bOMP=False )
    Gs = mmff.fit3D_Bspline( E, dt=0.1, nmaxiter=1000, Ftol=1e-6, bOMP=True, bPBC=bPBC )

    #Bspline_Pauli = mmff.getArrayPointer( "Bspline_Pauli" ); print( "Bspline_Pauli ", Bspline_Pauli.shape, Bspline_Pauli[0,0,:] )

    plt.figure(figsize=(15,10))

    plt.subplot(2,3,1); plt.imshow( E[iz0  ,:,:],  origin="lower" ) ;plt.colorbar(); plt.title("Es[iz= 0]")
    plt.subplot(2,3,2); plt.imshow( E[iz0-1,:,:],  origin="lower" ) ;plt.colorbar(); plt.title("Es[iz=-1]")
    plt.subplot(2,3,3); plt.imshow( E[iz0+1,:,:],  origin="lower" ) ;plt.colorbar(); plt.title("Es[iz=+1]")

    plt.subplot(2,3,4); plt.imshow( Gs[iz0  ,:,:],  origin="lower" ) ;plt.colorbar(); plt.title("Gs[iz= 0]")
    plt.subplot(2,3,5); plt.imshow( Gs[iz0-1,:,:],  origin="lower" ) ;plt.colorbar(); plt.title("Gs[iz=-1]")
    plt.subplot(2,3,6); plt.imshow( Gs[iz0+1,:,:],  origin="lower" ) ;plt.colorbar(); plt.title("Gs[iz=+1]")

    #plt.subplot(2,3,1); plt.imshow( E[iz0  ,:,:],  origin="lower" ) ;plt.colorbar(); plt.title("Es[iz= 0]")
    #plt.subplot(2,3,2); plt.imshow( E[iz0-1,:,:],  origin="lower" ) ;plt.colorbar(); plt.title("Es[iz=-1]")
    #plt.subplot(2,3,3); plt.imshow( E[iz0+1,:,:],  origin="lower" ) ;plt.colorbar(); plt.title("Es[iz=+1]")

    #plt.subplot(2,3,4); plt.imshow( Gs[iz0  ,:,:],  origin="lower" ) ;plt.colorbar(); plt.title("Gs[iz= 0]")
    #plt.subplot(2,3,5); plt.imshow( Gs[iz0-1,:,:],  origin="lower" ) ;plt.colorbar(); plt.title("Gs[iz=-1]")
    #plt.subplot(2,3,6); plt.imshow( Gs[iz0+1,:,:],  origin="lower" ) ;plt.colorbar(); plt.title("Gs[iz=+1]")

    #E_f    = mmff.sample_Bspline3D( ps, Gs, g0, dg, fes=None  ).reshape(sh_samp+(3,))

    #Gmin = -np.abs(Gs).max()
    #dG = Gs-E                       ; dGmin = -np.abs(dG).max()
    #dE = E_f[:,:,:,2] - E_r[:,:,:]  ; dEmin = -np.abs(dE[4:-4,4:-4,4:-4]).max()

    # extent=(g0[0],gmax[0],g0[1],gmax[1])
    # plt.figure(figsize=(15,10))
    # plt.subplot(2,3,1); plt.imshow( E,         origin="lower", extent=extent, vmin=Emin,  vmax=-Emin,  cmap=cmap ) ;plt.colorbar(); plt.title("Eg")
    # plt.subplot(2,3,2); plt.imshow( Gs,        origin="lower", extent=extent, vmin=Gmin,  vmax=-Gmin,  cmap=cmap ) ;plt.colorbar(); plt.title("Gs fit")
    # plt.subplot(2,3,3); plt.imshow( dG,        origin="lower", extent=extent, vmin=dGmin, vmax=-dGmin, cmap=cmap ) ;plt.colorbar(); plt.title("Gs-Eg")
    
    # plt.subplot(2,3,4); plt.imshow( E_r,        origin="lower", extent=extent, vmin=Emin, vmax=-Emin,   cmap=cmap ) ;plt.colorbar(); plt.title("E  ref")
    # plt.subplot(2,3,5); plt.imshow( E_f[:,:,2], origin="lower", extent=extent, vmin=Emin, vmax=-Emin,   cmap=cmap ) ;plt.colorbar(); plt.title("E  fit")
    # plt.subplot(2,3,6); plt.imshow( dE,         origin="lower", extent=extent, vmin=dEmin, vmax=-dEmin, cmap=cmap ) ;plt.colorbar(); plt.title("E(fit-ref)")
    # plt.axis('equal')




def test_comb3_3D(g0=(-2.0,-2.0,-2.0), gmax=(2.0,2.0,2.5), dg=(0.1,0.1,0.1), dsamp=(0.05,0.05,0.05), title=None, scErr=1000.0, bPBC=True, Ccomb=[1.0,0.0,0.0], bFit=True, iax=2, Gs=None ):
    print("test_comb3_3D START")
    cmap = "bwr"
    
    x0=0.0; y0=0.0; z0=0.0;
    ts = np.arange( g0[iax]-1.0, gmax[iax]*2.0, dsamp[iax])
    ps = np.zeros( (len(ts), 3,) )
    
    ps[:,0] = x0
    ps[:,1] = y0
    ps[:,2] = z0
    ps[:,iax] = ts

    E1_, Fx1, Fy1, Fz1 = fu.getCos3D( ps[:,0],   ps[:,1],   ps[:,2] )
    E2_, Fx2, Fy2, Fz2 = fu.getCos3D( ps[:,0]*2, ps[:,1],   ps[:,2] )
    E3_, Fx3, Fy3, Fz3 = fu.getCos3D( ps[:,0],   ps[:,1]*2, ps[:,2] )

    E_r =  E1_ * Ccomb[0] + E2_ * Ccomb[1] + E3_ * Ccomb[2]   

    if  iax == 0:
        F_r =   Fx1 * Ccomb[0] + Fx2 * Ccomb[1] + Fx3 * Ccomb[2]
    elif iax == 1:
        F_r =  Fy1 * Ccomb[0] + Fy2 * Ccomb[1] + Fy3 * Ccomb[2]
    elif iax == 2:
        F_r =  Fz1 * Ccomb[0] + Fz2 * Ccomb[1] + Fz3 * Ccomb[2]


    if Gs is None:
        # Create 3D sampling grids
        Xs, Ys, Zs = fu.make3Dsampling(g0=g0, gmax=gmax, dg=dg)
        tgs = np.arange( g0[iax], gmax[iax], dg[iax])
        # Generate 3D cosine functions
        E1, Fx, Fy, Fz  = fu.getCos3D(Xs, Ys, Zs)
        E2, Fx, Fy, Fz  = fu.getCos3D(Xs * 2, Ys, Zs)
        E3, Fx, Fy, Fz  = fu.getCos3D(Xs, Ys * 2, Zs)

        # Initialize Gs tensor for 3D case
        Gs = np.zeros(E1.shape + (3,))

        if bFit:
            G1 = mmff.fit3D_Bspline(E1.transpose( (2,1,0) ).copy(), dt=0.4, nmaxiter=1000, Ftol=1e-10, bPBC=bPBC)
            G2 = mmff.fit3D_Bspline(E2.transpose( (2,1,0) ).copy(), dt=0.4, nmaxiter=1000, Ftol=1e-10, bPBC=bPBC)
            G3 = mmff.fit3D_Bspline(E3.transpose( (2,1,0) ).copy(), dt=0.4, nmaxiter=1000, Ftol=1e-10, bPBC=bPBC)
            Gs[:,:,:,0] = G1.transpose( (2,1,0) )
            Gs[:,:,:,1] = G2.transpose( (2,1,0) )
            Gs[:,:,:,2] = G3.transpose( (2,1,0) )
        else:
            Gs[:,:,:,0] = E1
            Gs[:,:,:,1] = E2
            Gs[:,:,:,2] = E3

        #print( "Gs.shape BEFOR ", Gs.shape )
        #Gs = Gs.transpose( (2,1,0,3) ).copy()
        #print( "Gs.shape AFTER ", Gs.shape )

        plt.figure()
        plt.imshow(Gs[:,:,Gs.shape[2]//2,0], origin="lower", extent=(g0[0],gmax[0],g0[1],gmax[1]))

        #dpx = 0.05

    ix0 = Gs.shape[0]//2
    iy0 = Gs.shape[1]//2
    iz0 = Gs.shape[2]//2
    print( "ix0,iy0,iax Gs.shape", ix0,iy0,iz0,iax, Gs.shape )
    
    plt.figure(figsize=(5,10))
    # plt.subplot(2,1,1)
    # plt.plot( Gs[:,iy0,iz0,0], '.-', label="G(x)")
    # plt.plot( Gs[ix0,:,iy0,0], '.-', label="G(y)")
    # plt.plot( Gs[ix0,iy0,:,0], '.-', label="G(z)")

    E_f = mmff.sample_Bspline3D_comb3(ps, Gs, g0, dg, Cs=Ccomb)
    Err = E_r - E_f[:,3]
    Frr = E_f[:,iax] - F_r

    plt.subplot(2,1,1)    
    plt.plot( ts,  E_f[:,3], '-' , lw=0.5, label=("E_fit(z)" ))
    plt.plot( ts,  E_r     , 'k:', lw=1.0, label=("E_ref(z)" )    )
    plt.plot( ts, Err*scErr, 'r-', lw=1.0, label=("Err(z)*%g" % scErr))
    Emin = E_r.min()*1.2
    plt.ylim(Emin,-Emin)
    plt.grid()
    plt.legend()

    plt.subplot(2,1,2)
    plt.plot( ts,  E_f[:,iax],'-' , lw=0.5, label=("F_fit(z)" ))
    plt.plot( ts,  F_r       ,'k:', lw=1.0, label=("F_ref(z)" )    )
    plt.plot( ts,  Frr*scErr, 'r-', lw=1.0, label=("F_err(z)*%g" % scErr))
    Fmin = F_r.min()*1.2
    plt.ylim(Fmin,-Fmin)
    plt.grid()
    plt.legend()



    if title is not None: plt.title(title)
    plt.savefig( "test_comb3_3D_iax_%i.png" %iax, bbox_inches='tight' )
    print("test_comb3_3D DONE")

    return Gs

def test_PBCindexes( n=60, ng=20, order=3 ):
    
    inds = np.array( range(n), dtype=np.int32 )-(n//2)
    iout = mmff.samplePBCindexes( inds, ng=ng, order=order )
    plt.figure()
    for i in range(order+1):
        plt.plot( inds, iout[:,i], '.-', label=("%i" %i ) )
    plt.grid()
    plt.suptitle("order="+str(order) )

def test_project1D( xs,  g0=0.0, dg=0.1, ng=20, order=3 ):
    L = ng*dg
    ys = mmff.projectBspline1D( xs, g0, dg, ng, order=order )
    Qtot = ys.sum(); print( "Qtot ",Qtot," order=",order," xs=", xs );

    xg = np.arange( ng )*dg + g0    #;print("xg ", xg )

    xg = np.concatenate([xg, xg+L ])  # Shift and concatenate x-array
    ys = np.concatenate([ys, ys   ]) 

    plt.figure()
    #plt.plot( xg, ys, '.-', label=("order=%i" %order) )
    plt.plot( xg, ys, '.-', label=("x[0]=%7.4f Qtot=%g" %(xs[0],Qtot) ) )
    plt.axvline(L*0,ls='--', c='r')
    plt.axvline(L*1,ls='--', c='r')
    plt.axvline(L*2,ls='--', c='r')
    plt.grid()
    plt.legend()
    #plt.suptitle("order="+str(order) )

def test_project2D( xs,  g0=[0.0,0.0], dg=[0.1,0.1], ng=[16,16], order=3, ws=None ):
    Ls = [ng[0]*dg[0], ng[1]*dg[1]  ]
    ys = mmff.projectBspline2D( xs, g0, dg, ng, order=order, ws=ws )
    Qtot = ys.sum(); QtotAbs=np.abs(ys).sum();
    #print( "Qtot ",Qtot," order=",order," xs=", xs );
    print( "Qtot ",QtotAbs," order=",order, " g0 ", g0 );

    #xg = np.arange( ng )*dg + g0      #;print("xg ", xg )
    #xg = np.concatenate([xg, xg+L ])  # Shift and concatenate x-array
    #ys = np.concatenate([ys, ys   ]) 

    plt.figure()
    #plt.plot( xg, ys, '.-', label=("order=%i" %order) )
    plt.imshow( ys, origin="lower", extent=(g0[0],g0[0]+Ls[0],g0[1],g0[1]+Ls[1]), vmin=-1.0,vmax=1.0,  cmap='bwr' )

    #plt.axvline(L*0,ls='--', c='r')
    #plt.axvline(L*1,ls='--', c='r')
    #plt.axvline(L*2,ls='--', c='r')
    #plt.grid()
    #plt.legend()
    plt.colorbar()
    plt.suptitle("g0="+str(g0)+"order="+str(order) )


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
test_fit_1D( g0=0.0, ng=8, dg=0.25 )




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