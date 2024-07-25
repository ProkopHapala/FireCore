import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF as mmff


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

def getCos( xs, freq ):
    E =       np.cos( freq*xs )
    F = -freq*np.sin( freq*xs )
    return E,F

def getLJ( r, R0, E0 ):
    #r = np.sqrt(x**2 + y**2 + z**2)
    E  = E0*    ( (R0/r)**12 - 2*(R0/r)**6 )
    fr = E0*-12*( (R0/r)**12 -   (R0/r)**6 )/(r*r)
    return E,fr*r

def getLJ_atoms( apos, REs, Xs,Ys,Zs ):
    ng  = len(Xs)
    #Es = np.zeros( ng )
    #Fx = np.zeros( ng )
    #Fy = np.zeros( ng )
    #Fz = np.zeros( ng )
    Es = Xs*0.0
    Fx = Xs*0.0
    Fy = Xs*0.0
    Fz = Xs*0.0
    for i,p in enumerate(apos):
        print( p )
        dx = Xs-p[0]
        dy = Ys-p[1]
        dz = Zs-p[2]
        r = np.sqrt( dx**2 + dy**2 + dz**2 )
        R0,E0 = REs[i]
        E, fr = getLJ( r, R0, E0 )
        Es += E
        Fx += fr*dx/r
        Fy += fr*dy/r
        Fz += fr*dz/r
    return Es, Fx,Fy,Fz

def test_fit_1D( g0=2.0, gmax=10.0, dg=0.1, dsamp=0.02, bUseForce=True ):
    #x0 = 2.0
    #dx = 0.1
    xs  = np.arange(g0, gmax+1e-8, dg)     ; ng=len(xs)
    xs_ = np.arange(g0, gmax+1e-8, dsamp)  ; nsamp=len(xs_)
    #E,F = getLJ( xs, 3.5, 1.0 )
    E,F         = getCos( xs,  np.pi )
    E_ref,F_ref = getCos( xs_, np.pi )

    print( "E_ref ", E_ref )
    #print( "F_ref ", F_ref )

    Emin =  E.min()
    Fmin = -F.max()
    #print( "Emin ", Emin," Fmin ", Fmin )
    Ecut = 100.0
    Ws = Ecut/np.sqrt( E**2 + Ecut**2 )  ; EWs = E*Ws
    #E*=Ws
    #FEg[:,1]*=-1
    if bUseForce:
        FEg = np.zeros( (len(xs),2) )
        FEg[:,0] = E[:]
        FEg[:,1] = F[:]
        Gs = E*0.9
        Ws     = np.zeros((ng,2));  Ws[:,0]=1.0; Ws[:,1]=0.0;
        Gs, Ws = mmff.fitEF_Bspline( dg, FEg, Gs=Gs, Ws=Ws, Ftol=1e-6, nmaxiter=1000, dt=0.1 )
    else:
        Gs, Ws = mmff.fit_Bspline( E, Ws=None, dt=0.4, nmaxiter=1000, Ftol=1e-7 )
        #Gs, Ws = mmff.fit_Bspline( FEg[:,0].copy(), Ws=Ws,   dt=0.4, nmaxiter=1000, Ftol=1e-7 )
    FEout = mmff.sample_Bspline( xs_, Gs, x0=g0, dx=dg )
    plt.figure(figsize=(5,10))
    plt.subplot(2,1,1)    
    #plt.plot( xs, EWs, ".-k" )
    #plt.plot( xs, E, ".-k", label="E_ref" )
    plt.plot( xs_, E_ref, "-k", label="E_ref" )
    plt.plot( xs, Gs, ".-m", label="Gs" )
    plt.plot( xs_, FEout[:,0], "-b", label="E_fit" )
    #print( "Gs: ", Gs )
    plt.ylim(Emin*1.2,-Emin*1.2)
    plt.grid()
    plt.legend()
    plt.title("Energy")
    plt.subplot(2,1,2)
    #plt.plot( xs,  -F ,     ".-k", label="F_ref" )    
    plt.plot( xs_, -F_ref , "-k", label="F_ref" )    
    plt.plot( xs_, -FEout[:,1], "-b", label="F_fit" )
    #plt.ylim(Fmin*1.2,-Fmin*1.2)
    plt.legend()
    plt.title("Force")
    plt.grid()
    plt.show()

def make2Dsampling(  g0=(-5.0,2.0), gmax=(5.0,10.0), dg=(0.1,0.1) ):
    xs  = np.arange(g0[0], gmax[0], dg[0])
    ys  = np.arange(g0[1], gmax[1], dg[0])
    Xs,Ys = np.meshgrid(xs,ys)
    return Xs,Ys

def make2Dsampling_ps(  g0=(-5.0,2.0), gmax=(5.0,10.0), dg=(0.1,0.1) ):
    Xs,Ys = make2Dsampling(  g0=g0, gmax=gmax, dg=dg )
    sh = Xs.shape
    ps = np.zeros( ( len(Xs.flat), 2) )
    ps[:,0] = Xs.flat
    ps[:,1] = Ys.flat
    return ps, sh

def test_fit_2D( g0=(-5.0,2.0), gmax=(5.0,10.0), dg=(0.1,0.1), dsamp=(0.05,0.05) ):
    #cmap="RdBu_r"
    cmap="bwr"
    #x0 = 2.0
    #dx = 0.1
    Xs,Ys = make2Dsampling(  g0=g0, gmax=gmax, dg=dg )
    ps, sh_samp = make2Dsampling_ps(  g0=g0, gmax=gmax, dg=dsamp )

    print( "Xs.shape ", Xs.shape )
    
    E, Fx,Fy,Fz = getLJ_atoms( apos, REs, Xs,Ys,Xs*0.0 )

    #xs_ = np.arange(g0, gmax, dsamp)  ; nsamp=len(xs_)
    Emin =  E.min()
    Fmin = -Fz.max()
    print( "Emin ", Emin," Fmin ", Fmin )
    #FEg = np.zeros( (len(xs),2) )
    #FEg[:,0] = E[:]
    #FEg[:,1] = F[:]
    Ecut = 100.0

    Gs, Ws = mmff.fit2D_Bspline( E, Ws=None, dt=0.4, nmaxiter=1, Ftol=1e-7 )

    E_ = mmff.sample_Bspline2D( ps, Gs, g0, dg, fes=None  )

    extent=(g0[0],gmax[0],g0[1],gmax[1])
    plt.figure(figsize=(15,5))
    plt.subplot(1,3,1); plt.imshow( E,                   origin="lower", extent=extent, vmin=Emin, vmax=-Emin, cmap=cmap ) ;plt.title("E  ref")
    plt.subplot(1,3,2); plt.imshow( Gs,                  origin="lower", extent=extent, vmin=Emin, vmax=-Emin, cmap=cmap ) ;plt.title("Gs fit")
    plt.subplot(1,3,3); plt.imshow( E_.reshape(sh_samp), origin="lower", extent=extent, vmin=Emin, vmax=-Emin, cmap=cmap ) ;plt.title("E  fit")
    plt.axis('equal')
    '''
    #Ws = Ecut/np.sqrt( E**2 + Ecut**2 )  ; EWs = E*Ws
    Gs, Ws = mmff.fit_Bspline2D( E, Ws=None, dt=0.4, nmaxiter=1000, Ftol=1e-7 )
    FEout = mmff.sample_Bspline( xs_, Gs, x0=g0, dx=dg )
    plt.figure(figsize=(5,10))
    plt.subplot(2,1,1)    
    #plt.plot( xs, EWs, ".-k" )
    plt.plot( xs, E, ".-k", label="E_ref" )
    plt.plot( xs, Gs, ".-m", label="Gs" )
    plt.plot( xs_, FEout[:,0], "-b", label="E_fit" )
    #print( "Gs: ", Gs )
    plt.ylim(Emin,-Emin)
    plt.grid()
    plt.legend()
    plt.title("Energy")
    plt.subplot(2,1,2)
    plt.plot( xs,  -F , ".-k", label="F_ref" )    
    plt.plot( xs_, -FEout[:,1], "-b", label="F_fit" )
    plt.ylim(Fmin,-Fmin)
    plt.legend()
    plt.title("Force")
    plt.grid()
    plt.show()
    '''


mmff.setVerbosity( 2 )
#mmff.setVerbosity( 3 )

#test_fit_1D( bUseForce=True )
#test_fit_1D( bUseForce=False )
test_fit_1D( g0=0.0, gmax=2.0, dg=0.1, bUseForce=False )

#test_fit_2D(  )

plt.show()