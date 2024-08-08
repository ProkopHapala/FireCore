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


def test_fit_1D( g0=2.0, gmax=10.0, dg=0.2, dsamp=0.02, bUseForce=True, scErr=100.0, bHalf=False, title=None ):
    #x0 = 2.0
    #dx = 0.1
    if bHalf: dg = dg*0.5
    xs  = np.arange(g0, gmax+1e-8, dg)     ; ng=len(xs)
    xs_ = np.arange(g0, gmax+1e-8, dsamp)  ; nsamp=len(xs_)
    xsg=xs; 
    if bHalf: xsg=xs[::2]
    
    E,F         = fu.getLJ( xs,  3.5, 1.0 )
    E_ref,F_ref = fu.getLJ( xs_, 3.5, 1.0 )
    
    #E,F         = getCos( xs,  np.pi )
    #E_ref,F_ref = getCos( xs_, np.pi )

    #print( "E_ref ", E_ref )
    #print( "F_ref ", F_ref )

    Emin =  E.min()
    Fmin = -F.max()
    #print( "Emin ", Emin," Fmin ", Fmin )
    Ecut = -2.0
    #Ws = 1/np.sqrt( (E/Ecut)**2 + 1 ); # EWs = E*Ws
    Ws = 1/( E - Ecut );
    #E*=Ws
    #FEg[:,1]*=-1

    #Ecut2 =  1.0
    #mask    = E>Ecut2
    #E[mask] = Ecut2

    if bUseForce:
        FEg = np.zeros( (len(xs),2) )
        FEg[:,0] = E[:]
        FEg[:,1] = F[:]
        Gs = E*0.9
        Ws     = np.zeros((ng,2));  Ws[:,0]=1.0; Ws[:,1]=0.0;
        Gs, Ws = mmff.fitEF_Bspline( dg, FEg, Gs=Gs, Ws=Ws, Ftol=1e-6, nmaxiter=1000, dt=0.1 )
    else:

        Gs, Ws_ = mmff.fit_Bspline( E, Ws=None, dt=0.4, nmaxiter=1000, Ftol=1e-7, bHalf=bHalf )

        #Gs, Ws_ = mmff.fit_Bspline( E, Ws=Ws,   dt=1.0, nmaxiter=1000, Ftol=1e-9, bHalf=bHalf )

        #Gs, Ws = mmff.fit_Bspline( FEg[:,0].copy(), Ws=Ws,   dt=0.4, nmaxiter=1000, Ftol=1e-7 )
    dgs = dg
    if bHalf: dgs=dg*2
    FEout = mmff.sample_Bspline( xs_, Gs, x0=g0, dx=dgs )
    plt.figure(figsize=(5,10))
    plt.subplot(2,1,1)    
    plt.plot( xs, Ws, ":m", lw=1.5, label="Ws" )
    #plt.plot( xs, EWs, ".-k" )
    #plt.plot( xs, E, ".-k", label="E_ref" )
    #plt.plot( xsg, Gs,          ".-m", lw=0.25, label="Gs" )
    plt.plot( xs_, E_ref,      "-k",  lw=0.5,  label="E_ref" )
    plt.plot( xs_, FEout[:,0], "-b",  lw=0.5,  label="E_fit" )
    plt.plot( xs_, (FEout[:,0]-E_ref)*scErr, "-r", lw=0.5, label=("error*%g" % scErr) )
    #print( "Gs: ", Gs )
    plt.ylim(Emin*1.2,-Emin*1.2)
    plt.grid()
    plt.legend()
    plt.title("Energy")
    plt.subplot(2,1,2)
    #plt.plot( xs,  -F ,       ".-k", lw=0.5, label="F_ref" )    
    plt.plot( xs_, -F_ref ,     "-k", lw=0.5, label="F_ref" )    
    plt.plot( xs_, -FEout[:,1], "-b", lw=0.5, label="F_fit" )
    plt.plot( xs_, (FEout[:,1]-F_ref)*scErr, "-r", lw=0.5, label=("error*%g" % scErr) )
    plt.ylim(Fmin*1.2,-Fmin*1.2)
    plt.legend()
    plt.title("Force")
    plt.grid()
    if title is not None: plt.suptitle(title)
    #plt.show()


def test_fit_2D( g0=(-5.0,2.0), gmax=(5.0,10.0), dg=(0.1,0.1), dsamp=(0.05,0.05) ):
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

    Gs, Ws = mmff.fit2D_Bspline( E, Ws=None, dt=0.4, nmaxiter=1000, Ftol=1e-7 )
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

def test_fit_3D( g0=(-5.0,-5.0,2.0), gmax=(5.0,-5.0,10.0), dg=(0.1,0.1,0.1), dsamp=(0.05,0.05,0.05) ):
    cmap="bwr"
    Xs,Ys,Ys    = fu.make3Dsampling(  g0=g0, gmax=gmax, dg=dg )
    Xs_,Ys_,Ys_ = fu.make3Dsampling(  g0=g0, gmax=gmax, dg=dsamp )

    Xs_*=0.999999; Ys_*=0.999999; Zs_*=0.999999;

    sh_samp = Xs_.shape
    ps      = fu.pack_ps3D(Xs_,Ys_,Ys_)
    #ps, sh_samp = make2Dsampling_ps(  g0=g0, gmax=gmax, dg=dsamp )
    print( "Xs.shape ", Xs.shape )
    #E, Fx,Fy,Fz = getLJ_atoms( apos, REs, Xs,Ys,Xs*0.0 )
    E,  Fx,Fy      = fu.getCos3D( Xs, Ys   )
    E_r, Fx_r,Fy_r = fu.getCos3D( Xs_, Ys_ )

    Emin =  E.min()
    Fmin = -Fy.max()

    Gs, Ws = mmff.fit3D_Bspline( E, Ws=None, dt=0.4, nmaxiter=1000, Ftol=1e-7 )
    E_f    = mmff.sample_Bspline3D( ps, Gs, g0, dg, fes=None  ).reshape(sh_samp+(3,))

    Gmin = -np.abs(Gs).max()
    dG = Gs-E                       ; dGmin = -np.abs(dG).max()
    dE = E_f[:,:,:,2] - E_r[:,:,:]  ; dEmin = -np.abs(dE[4:-4,4:-4,4:-4]).max()

    # extent=(g0[0],gmax[0],g0[1],gmax[1])
    # plt.figure(figsize=(15,10))
    # plt.subplot(2,3,1); plt.imshow( E,         origin="lower", extent=extent, vmin=Emin,  vmax=-Emin,  cmap=cmap ) ;plt.colorbar(); plt.title("Eg")
    # plt.subplot(2,3,2); plt.imshow( Gs,        origin="lower", extent=extent, vmin=Gmin,  vmax=-Gmin,  cmap=cmap ) ;plt.colorbar(); plt.title("Gs fit")
    # plt.subplot(2,3,3); plt.imshow( dG,        origin="lower", extent=extent, vmin=dGmin, vmax=-dGmin, cmap=cmap ) ;plt.colorbar(); plt.title("Gs-Eg")
    
    # plt.subplot(2,3,4); plt.imshow( E_r,        origin="lower", extent=extent, vmin=Emin, vmax=-Emin,   cmap=cmap ) ;plt.colorbar(); plt.title("E  ref")
    # plt.subplot(2,3,5); plt.imshow( E_f[:,:,2], origin="lower", extent=extent, vmin=Emin, vmax=-Emin,   cmap=cmap ) ;plt.colorbar(); plt.title("E  fit")
    # plt.subplot(2,3,6); plt.imshow( dE,         origin="lower", extent=extent, vmin=dEmin, vmax=-dEmin, cmap=cmap ) ;plt.colorbar(); plt.title("E(fit-ref)")
    # plt.axis('equal')

#mmff.setVerbosity( 2 )
mmff.setVerbosity( 3 )

#test_fit_1D( bUseForce=True )
test_fit_1D( bUseForce=False, bHalf=False ,title="No-Half")
test_fit_1D( bUseForce=False, bHalf=True  ,title="Half")
#test_fit_1D( g0=0.0, gmax=2.0, dg=0.1, bUseForce=False )

#test_fit_2D(  )
#test_fit_2D( g0=(-1.0,-1.0), gmax=(1.0,1.0) )

plt.show()