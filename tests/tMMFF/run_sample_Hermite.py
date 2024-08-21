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


def test_1D( func, funcS, g0=2.0, gmax=10.0, dg=0.2, dsamp=0.02, scErr=100.0, title=None, mode=1 ):
    xs  = np.arange(g0, gmax+1e-8, dg)     ; ng=len(xs)
    xs_ = np.arange(g0, gmax+1e-8 -3*dg, dsamp)  ; nsamp=len(xs_)
    
    #E0 = 1.0
    #R0 = 3.5
    #a  = 1.7
    E_,F_ = func( xs_ )
    if   ( mode==1 ):   # simple Hermite using just values
        E,F   = func( xs)
        Gs = E
        print( "xs ",xs )
        FEout = mmff.sample_SplineHermite( xs_, Gs, g0=g0, dg=dg )
    elif ( mode==2 ):   # simple Hermite using derivatives
        E,F   = func( xs )
        Gs    = np.zeros( (ng,2) )
        Gs[:,0] = E
        Gs[:,1] = F
        FEout = mmff.sample_SplineHermite1D_deriv( xs_, Gs, g0=g0, dg=dg )
    elif ( mode==3 ):   # hermite using linear combination of splited function
        Ep,El,Fp,Fl = funcS( xs )
        Gs = np.zeros( (ng,4) )
        Gs[:,0] = Ep
        Gs[:,1] = Fp*dg
        Gs[:,2] = El
        Gs[:,3] = Fl*dg
        FEout = mmff.sample_SplineHermite_comb( xs_, Gs, [1.0,1.0], ncomb=2, g0=g0, dg=dg )

    # ------ plotting

    Emin =  E_.min()
    Fmin = -F_.max()

    plt.figure(figsize=(5,10))
    plt.subplot(2,1,1)    
    plt.plot( xs_, E_,                    "-k", lw=0.5,  label="E_ref" )
    plt.plot( xs_,  FEout[:,0]          , "-b", lw=0.5, label="E_fit" )
    plt.plot( xs_, (FEout[:,0]-E_)*scErr, "-r", lw=0.5, label=("error*%g" % scErr) )
    #print( "Gs: ", Gs )
    plt.ylim(Emin*1.2,-Emin*1.2)
    plt.grid()
    plt.legend()
    plt.title("Energy")
    plt.subplot(2,1,2)   
    plt.plot( xs_, -F_ ,                  "-k", lw=0.5, label="F_ref" )    
    plt.plot( xs_, -FEout[:,1],           "-b", lw=0.5, label="F_fit" )
    plt.plot( xs_, (FEout[:,1]-F_)*scErr, "-r", lw=0.5, label=("error*%g" % scErr) )
    plt.ylim(Fmin*1.2,-Fmin*1.2)
    plt.legend()
    plt.title("Force")
    plt.grid()
    if title is not None: plt.suptitle(title)
    #plt.show()


def test_fit_2D( g0=(-3.0,2.0), gmax=(3.0,7.0), dg=(0.2,0.2), dsamp=(0.05,0.05), scErr=100.0, title=None, mode=1 ):
    #cmap="RdBu_r"
    cmap="bwr"
    #x0 = 2.0
    #dx = 0.1
    Xs,Ys   = fu.make2Dsampling(  g0=g0, gmax=gmax, dg=dg   , bSwapXY=(mode==3) )
    Xs_,Ys_ = fu.make2Dsampling(  g0=g0, gmax=gmax, dg=dsamp )

    Xs_*=0.999999; Ys_*=0.999999;

    sh_samp = Xs_.shape
    ps      = fu.pack_ps2D( Xs_, Ys_)
    #ps, sh_samp = make2Dsampling_ps(  g0=g0, gmax=gmax, dg=dsamp )
    
    E,   Fx,Fy,_     = fu.getLJ_atoms( apos, REs, Xs, Ys, Xs*0.0 )
    E_r, Fx_r,Fy_r,_ = fu.getLJ_atoms( apos, REs, Xs_,Ys_,Xs_*0.0 )

    # ---- Cos2D
    #E,   Fx,  Fy   =  fu.getCos2D( Xs,  Ys  )
    #E_r, Fx_r,Fy_r =  fu.getCos2D( Xs_, Ys_ )

    ng=E.shape


    #xs_ = np.arange(g0, gmax, dsamp)  ; nsamp=len(xs_)
    Emin =  E.min();   Emax=-Emin
    Fmin = -Fy.max()
    print( "Emin ", Emin," Fmin ", Fmin )

    if ( mode==1 ):   # simple Hermite using just values
        Gs = E.copy()
        print( "Gs.shape ", Gs.shape )
        FEout = mmff.sample_SplineHermite2D( ps, Gs, g0, dg ).reshape(sh_samp+(3,))
    elif ( mode==2 ):   # simple Hermite using derivatives
        Gs = np.zeros( ng+(3,) )
        Gs[:,:,0] = E
        Gs[:,:,1] = Fx
        Gs[:,:,2] = Fy
        FEout = mmff.sample_SplineHermite2D_deriv( ps, Gs, g0, dg ).reshape(sh_samp+(3,))
    elif ( mode==3 ):   # hermite using linear combination of splited function
        Gs = np.zeros( ng+(4,) )
        Gs[:,:,0] = E
        Gs[:,:,1] = Fy*dg[1]
        Gs[:,:,2] = 0
        Gs[:,:,3] = 0
        E = E.transpose()
        FEout = mmff.sample_SplineHermite2D_comb( ps, Gs, g0, dg ).reshape(sh_samp+(3,))

    dE = FEout[:,:,2] - E_r[:,:]; #dEmin = -np.abs(dE[4:-4,4:-4]).max(); dEmin = 0.1

    #Emin = None; Emax = None

    extent=(g0[0],gmax[0],g0[1],gmax[1])
    plt.figure(figsize=(10,10))
    plt.subplot(2,2,1); plt.imshow( E,            origin="lower", interpolation='nearest', extent=extent, vmin=Emin, vmax=Emax,  cmap=cmap ) ;plt.colorbar(); plt.title("Eg")    
    #plt.subplot(1,2,1); plt.imshow( Gs,            origin="lower", interpolation='nearest', extent=extent, vmin=Emin, vmax=Emax,  cmap=cmap ) ;plt.colorbar(); plt.title("Eg")    
    plt.subplot(2,2,2); plt.imshow( E_r,          origin="lower", interpolation='nearest', extent=extent, vmin=Emin, vmax=Emax,  cmap=cmap ) ;plt.colorbar(); plt.title("E  ref")
    plt.subplot(2,2,3); plt.imshow( FEout[:,:,2], origin="lower", interpolation='nearest', extent=extent, vmin=Emin, vmax=Emax,  cmap=cmap ) ;plt.colorbar(); plt.title("E  fit")
    plt.subplot(2,2,4); plt.imshow( dE,           origin="lower", interpolation='nearest', extent=extent, vmin=Emin/scErr, vmax=-Emin/scErr, cmap=cmap ) ;plt.colorbar(); plt.title("E(fit-ref)")
    plt.axis('equal')


def plotGridFF_1D( ff, ix=20,iy=20 ):
    # ------- Plot GridFF
    plt.figure(figsize=(5,10))
    plt.subplot(2,1,1);
    plt.plot( ff[ix,iy,:, 0], "-b", lw=0.5, label="E_Pauli" )
    plt.plot( ff[ix,iy,:, 2], "-g", lw=0.5, label="E_London" )
    plt.plot( ff[ix,iy,:, 4], "-r", lw=0.5, label="E_Coulomb" )
    Fp_num = (ff[ix,iy,2:, 0] - ff[20,20,:-2, 0])/(-0.2)
    Fl_num = (ff[ix,iy,2:, 2] - ff[20,20,:-2, 2])/(-0.2)
    Fc_num = (ff[ix,iy,2:, 4] - ff[20,20,:-2, 4])/(-0.2)
    zs = np.arange(0,ff.shape[2])
    plt.subplot(2,1,2);
    plt.plot( ff[ix,iy,:, 1], "-b", lw=0.5, label="F_Pauli" )
    plt.plot( ff[ix,iy,:, 3], "-g", lw=0.5, label="F_London" )
    plt.plot( ff[ix,iy,:, 5], "-r", lw=0.5, label="F_Coulomb" )
    plt.plot( zs[1:-1], Fp_num, ":b", lw=2.0, label="F_Pauli"   )
    plt.plot( zs[1:-1], Fl_num, ":g", lw=2.0, label="F_London"  )
    plt.plot( zs[1:-1], Fc_num, ":r", lw=2.0, label="F_Coulomb" )

def getPLQH( R0, E0, a, Q, H ):
    e  = np.exp(a*R0);
    cL = e*E0;
    cP = e*cL;
    cH = e*e*H;
    return np.array([ cP, cL, Q, cH ])

def test_gridFF( name="data/NaCl_1x1_L2", mode=4, dsamp=0.02,  R0=3.5, E0=0.1, a=1.6, Q=0.4, H=0.0, scErr=100.0, title=None, ):
    print( "test_gridFF() START" )
    #mode = 4
    #mode = 1
    mmff.makeGridFF( name=name, mode=mode )
    #plotGridFF_1D( EFg, ix=20,iy=20 )
    PLQH = getPLQH( R0, E0, a, Q, H )

    zs = np.arange(0.0, 10.0, dsamp)
    ps = np.zeros( (len(zs), 3) )
    ps[:,0] = 1.05
    ps[:,1] = 1.05
    ps[:,2] = zs
    
    FF_ref = mmff.evalGridFFAtPoints( ps, PLQH=PLQH, bSplit=False )
    FF_ref = mmff.evalGridFFAtPoints( ps, PLQH=PLQH, bSplit=True  )
    
    #FFout  = mmff.sample_SplineHermite3D_comb3( ps, EFg, g0=[0.0,0.0,0.0], dg=[0.1,0.1,0.1], fes=None, Cs=PLQH )
    # Es,Fs = sampleSurf( name, zs, Es=None, fs=None, kind=1, atyp=0, Q=0.0, K=-1.0, Rdamp=1.0, pos0=(0.,0.,0.), bSave=False )

    ps_ = ps.copy();
    if ( mode==1 ):  # Trilinear interpolation of force (Float, simple prec.)
        ps_[:,2]+=-2.0;  
    if ( mode==2 ):  # Trilinear interpolation of force (Double prec.)
        ps_[:,2]+=1.2;  
    elif ( mode==4 ):   # Hibrid-Hermite interpolation of potential
        ps_[:,2]+=-2.0;
    elif ( mode==6 ):   # Hibrid-Hermite interpolation of potential
        ps_[:,2]+=-2.0;
        #ps_[:,0]+= 0.1
        #ps_[:,1]+= 0.1
    
    FFout = mmff.sampleSurf_new( ps_, PLQH, mode=mode, Rdamp=1.0 )
    
    #Emin = FFout[:,3].min();
    #Fmin = FFout[:,2].min();

    Emin = FF_ref[:,3].min(); 
    Fmin = FF_ref[:,2].min()

    plt.figure(figsize=(5,10))
    plt.subplot(2,1,1);
    plt.plot( zs, FFout[:,3],  "-g", lw=0.5, label="Etot_fit" )
    plt.plot( zs, FF_ref[:,3], ":k", lw=2.0, label="Etot_ref" )
    plt.plot( zs, (FFout[:,3]-FF_ref[:,3])*scErr, "-r", lw=0.5, label=("Etot_err*%.2f" %scErr) )
    #plt.plot( zs, FFout[:,0], "-r", lw=0.5, label="Ftot_x" )
    #plt.plot( zs, FFout[:,1], "-g", lw=0.5, label="Ftot_y" )
    #plt.plot( zs, FFout[:,2], "-b", lw=0.5, label="Etot_z" )
    plt.axhline(0.0, c="k", ls='--', lw=0.5)
    plt.ylim( Emin, -Emin )
    plt.legend()
    plt.subplot(2,1,2);
    plt.plot( zs, FFout[:,2],  "-g", lw=0.5, label="Ftot_fit" )
    plt.plot( zs, FF_ref[:,2], ":k", lw=2.0, label="Ftot_ref" )
    plt.plot( zs, (FFout[:,2]-FF_ref[:,2])*scErr, "-r", lw=0.5, label=("Ftot_err*%.2f" %scErr) )
    plt.axhline(0.0, c="k", ls='--', lw=0.5)
    plt.ylim( Fmin, -Fmin )
    plt.legend()

    if ( title is not None ): plt.suptitle( title )
    
    #print( "ff.shape ", EFg.shape )
    print( "test_gridFF() DONE" )
    #return EFg

def test_gridFF_lat( name="data/NaCl_1x1_L2", iax=0, tmin=0.0,tmax=10.0, p0=[1.05,1.05,2.0], mode=4, dsamp=0.02,  R0=3.5, E0=0.1, a=1.6, Q=0.4, H=0.0, scErr=100.0, title=None, ):
    print( "test_gridFF() START" )
    #mode = 4
    #mode = 1
    mmff.makeGridFF( name=name, mode=mode )
    #plotGridFF_1D( EFg, ix=20,iy=20 )
    PLQH = getPLQH( R0, E0, a, Q, H )

    ts = np.arange( tmin, tmax, dsamp)
    ps = np.zeros( (len(ts), 3) )
    ps[:,0] = p0[0]
    ps[:,1] = p0[1]
    ps[:,2] = p0[2]
    #ps[:,iax] = ts
    
    ps[:,0] = ts
    ps[:,1] = ts

    FF_ref = mmff.evalGridFFAtPoints( ps, PLQH=PLQH )
    
    #FFout  = mmff.sample_SplineHermite3D_comb3( ps, EFg, g0=[0.0,0.0,0.0], dg=[0.1,0.1,0.1], fes=None, Cs=PLQH )
    # Es,Fs = sampleSurf( name, zs, Es=None, fs=None, kind=1, atyp=0, Q=0.0, K=-1.0, Rdamp=1.0, pos0=(0.,0.,0.), bSave=False )

    ps_ = ps.copy(); ps_[:,2]+=-2.0;  

    FFout = mmff.sampleSurf_new( ps_, PLQH, mode=mode, Rdamp=1.0 )
    
    Emin = FFout[:,3].min();
    Fmin = FFout[:,2].min();

    #Emin = FF_ref[:,3].min(); 
    #Fmin = FF_ref[:,2].min()

    plt.figure(figsize=(5,10))
    plt.subplot(2,1,1);
    plt.plot( ts, FFout[:,3], "-g", lw=0.5, label="Etot_fit" )
    plt.plot( ts, FF_ref[:,3], ":k", lw=2.0, label="Etot_ref" )
    plt.plot( ts, (FFout[:,3]-FF_ref[:,3])*scErr, "-r", lw=0.5, label=("Etot_err*%.2f" %scErr) )
    #plt.plot( zs, FFout[:,0], "-r", lw=0.5, label="Ftot_x" )
    #plt.plot( zs, FFout[:,1], "-g", lw=0.5, label="Ftot_y" )
    #plt.plot( zs, FFout[:,2], "-b", lw=0.5, label="Etot_z" )
    plt.axhline(0.0, c="k", ls='--', lw=0.5)
    plt.ylim( Emin, -Emin )
    plt.legend()
    plt.subplot(2,1,2);
    plt.plot( ts, FFout[:,iax],  "-g", lw=0.5, label="Ftot_fit" )
    plt.plot( ts, FF_ref[:,iax], ":k", lw=2.0, label="Ftot_ref" )
    plt.plot( ts, (FFout[:,iax]-FF_ref[:,iax])*scErr, "-r", lw=0.5, label=("Ftot_err*%.2f" %scErr) )
    plt.axhline(0.0, c="k", ls='--', lw=0.5)
    plt.ylim( Fmin, -Fmin )
    plt.legend()

    if ( title is not None ): plt.suptitle( title )
    
    #print( "ff.shape ", EFg.shape )
    print( "test_gridFF() DONE" )
    #return EFg


R0 = 3.5
E0 = 1.0
a  = 1.8


#test_1D( lambda x:fu.getLJ(x,R0,E0), lambda x:fu.getLJSplit(x,R0,E0), title="LJ simple",      mode=1 )
#test_1D( lambda x:fu.getLJ(x,R0,E0), lambda x:fu.getLJSplit(x,R0,E0), title="LJ deriv",       mode=2 )
#test_1D( lambda x:fu.getLJ(x,R0,E0), lambda x:fu.getLJSplit(x,R0,E0), title="LJ deriv split", mode=3 )

#test_1D( lambda x:fu.getMorse(x,R0,E0), lambda x:fu.getMorseSplit(x,R0,E0), title="Morse simple",      mode=1 )
#test_1D( lambda x:fu.getMorse(x,R0,E0), lambda x:fu.getMorseSplit(x,R0,E0), title="Morse deriv",       mode=2 )
#test_1D( lambda x:fu.getMorse(x,R0,E0), lambda x:fu.getMorseSplit(x,R0,E0), title="Morse deriv split", mode=3 )

#fu.checkNumDeriv( lambda x: fu.getLJ(x,3.5,1.0),        2.0, 6.0, tol=1e-6, n=400, errSc=100.0, plt=plt, label="LJ"   ,c='r' )
#fu.checkNumDeriv( lambda x: fu.getMorse(x,3.5,1.0,1.8), 2.0, 6.0, tol=1e-6, n=400, errSc=100.0, plt=plt, label="Morse",c='b' )


#test_1D( title="No-Half", mode=1)
#test_1D( title="No-Half", mode=1)


#test_fit_2D( title="test mode=1", mode=1 )
#test_fit_2D( title="test mode=2", mode=3 )

mmff.initParams()
#test_gridFF( mode=1, title="tri-linar force \n(z-cut)"          )
#test_gridFF( mode=4, title="Hybrid Hermite tri-cubic\n(z-cut)" )
test_gridFF( mode=6, title="Bspline (from HH)\n(z-cut)" )

#test_gridFF_lat( mode=1, title="tri-linar force"          )
#test_gridFF_lat( mode=4, title="Hybrid Hermite tri-cubic" )

#Q = 0.0
Q = 0.4
#p0 = [1.0,-5.05,2.0]
p0 = [0.0,0.0,2.0]

#test_gridFF_lat( mode=1, title="tri-linar force \n(y-cut)"         , Q=Q, p0=p0, iax=1, tmin=-10, tmax=10 )
#test_gridFF_lat( mode=4, title="Hybrid Hermite tri-cubic \n(y-cut)", Q=Q, p0=p0, iax=1, tmin=-10, tmax=10. )

#test_gridFF_lat( mode=1, title="tri-linar force \n(x-cut)"         , Q=Q, p0=p0, iax=0 )
#test_gridFF_lat( mode=4, title="Hybrid Hermite tri-cubic \n(x-cut)", Q=Q, p0=p0, iax=0 )

#test_gridFF_lat( mode=4, title="Hybrid Hermite tri-cubic", Q=0.0, p0=[1.0,1.05,2.0], iax=1 )
#test_gridFF_lat( mode=4, title="Hybrid Hermite tri-cubic", Q=0.0, p0=[1.0,1.05,2.0], iax=0 )


test_gridFF_lat( mode=6, title="Bspline tri-cubic", Q=0.0, p0=[1.0,1.05,2.0], iax=0 )

plt.show()