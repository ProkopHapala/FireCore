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
    print( "py======= test_gridFF() START" );
    #print( "test_gridFF() START" )
    #mode = 4
    #mode = 1
    mmff.makeGridFF( name=name, mode=mode )
    #plotGridFF_1D( EFg, ix=20,iy=20 )
    PLQH = getPLQH( R0, E0, a, Q, H )

    zs = np.arange(0.0, 10.0, dsamp)
    ps = np.zeros( (len(zs), 3) )
    #ps[:,0] = 1.05
    #ps[:,1] = 1.05
    ps[:,0] = 0.0
    ps[:,1] = 0.0
    ps[:,2] = zs
    
    FF_ref = mmff.evalGridFFAtPoints( ps, PLQH=PLQH, bSplit=False )
    FF_ref = mmff.evalGridFFAtPoints( ps, PLQH=PLQH, bSplit=True  )
    
    #FFout  = mmff.sample_SplineHermite3D_comb3( ps, EFg, g0=[0.0,0.0,0.0], dg=[0.1,0.1,0.1], fes=None, Cs=PLQH )
    # Es,Fs = sampleSurf( name, zs, Es=None, fs=None, kind=1, atyp=0, Q=0.0, K=-1.0, Rdamp=1.0, pos0=(0.,0.,0.), bSave=False )

    ps_ = ps.copy();

    ps_[:,2]+=-2.0+3.25;
    ps_[:,0]+=2.0;
    ps_[:,1]+=2.0;
    
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
    #print( "test_gridFF() DONE" )
    print( "py======= test_gridFF() DONE" );
    #return EFg

def test_gridFF_lat( name="data/NaCl_1x1_L2", iax=0, tmin=0.0,tmax=10.0, p0=[1.05,1.05,2.0], mode=4, dsamp=0.02,  R0=3.5, E0=0.1, a=1.6, Q=0.4, H=0.0, scErr=100.0, title=None, ):
    print( "py======= test_gridFF_lat() START" );
    #print( "test_gridFF() START" )
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
    ps[:,iax] = ts
    #ps[:,iax] += 0.05
    
    #ps[:,0] = ts
    #ps[:,1] = ts

    FF_ref = mmff.evalGridFFAtPoints( ps, PLQH=PLQH )
    
    #FFout  = mmff.sample_SplineHermite3D_comb3( ps, EFg, g0=[0.0,0.0,0.0], dg=[0.1,0.1,0.1], fes=None, Cs=PLQH )
    # Es,Fs = sampleSurf( name, zs, Es=None, fs=None, kind=1, atyp=0, Q=0.0, K=-1.0, Rdamp=1.0, pos0=(0.,0.,0.), bSave=False )

    ps_ = ps.copy(); 
    ps_[:,2]+=-2.0+3.25;
    ps_[:,0]+=2.0;
    ps_[:,1]+=2.0;

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
    #print( "test_gridFF() DONE" )
    print( "py======= test_gridFF_lat() DONE" );
    #return EFg


R0 = 3.5
E0 = 1.0
a  = 1.8

#Q = 0.0
Q = 0.4
#p0 = [1.0,-5.05,2.0]
#p0 = [0.0,0.0,2.0]
p0 = [-2.0,-2.0,0.0]

mmff.initParams()

test_gridFF    ( mode=6, title="Bspline (from HH)\n(z-cut)" )
test_gridFF_lat( mode=6, title="Bspline tri-cubic", Q=0.0, p0=p0, iax=0 )

#test_gridFF    ( mode=1, title="tri-linar force \n(z-cut)"          )
#test_gridFF_lat( mode=1, title="tri-Linear Force", Q=0.0, p0=p0, iax=0 )

plt.show()