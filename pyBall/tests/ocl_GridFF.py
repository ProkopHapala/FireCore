import os
import numpy as np
import matplotlib.pyplot as plt

from .  import utils as ut
from .. import atomicUtils as au
from .. import MMFF as mmff
from .. import FunctionSampling as fu
from ..OCL.GridFF import GridFF_cl

# =============  Functions

os.environ['PYOPENCL_CTX'] = '0'
ocl_splines = GridFF_cl()

def test_gridFF_vs_ocl( path="data/xyz/NaCl_1x1_L2", mode=6, dsamp=0.02, p0=[0.0,0.0,2.0], R0=3.5, E0=0.1, a=1.6, Q=0.4, H=0.0, scErr=100.0, iax=2, Emax=None, Fmax=None, maxSc=5.0, title=None, bSaveFig=True, bRefine=True, nPBC=[100,100,0], bRealSpace=False ):
    print( "py======= test_gridFF_vs_ocl() START" );
    #print( "test_gridFF() START" )
    #mode = 4
    #mode = 1
    mmff.makeGridFF( name=path, mode=mode, bRefine=bRefine )
    #plotGridFF_1D( EFg, ix=20,iy=20 )
    PLQH = ut.getPLQH( R0, E0, a, Q, H )

    ps,ts = ut.make_sample_points( p0, t0=0.0, tmax=10.0, dsamp=dsamp, iax=iax )    
        
    #FFout  = mmff.sample_SplineHermite3D_comb3( ps, EFg, g0=[0.0,0.0,0.0], dg=[0.1,0.1,0.1], fes=None, Cs=PLQH )
    # Es,Fs = sampleSurf( name, zs, Es=None, fs=None, kind=1, atyp=0, Q=0.0, K=-1.0, Rdamp=1.0, pos0=(0.,0.,0.), bSave=False )

    ps_ = ps.copy();    
    FFout = mmff.sampleSurf_new( ps_, PLQH, mode=mode, Rdamp=1.0 )
    
    #Emin = FFout[:,3].min();
    #Fmin = FFout[:,2].min();

    # ---- OpenCL
    g0   = (-2.0,-2.0,0.0)
    #g0   = (0.0,0.0,0.0)
    #g0   = (2.0,0.0,0.0)
    #g0   = (2.0,2.0,0.0)
    dg   = (0.1,0.1,0.1)
    ps_ocl = ut.points_to_ocl(ps)
    VPLQ = ut.load_potential_comb( path )  #;print( "VPLQ.shape=", VPLQ.shape )
    ocl_splines.prepare_sample3D( g0, dg, VPLQ.shape[:3], VPLQ )
    fe_cl = ocl_splines.sample3D_comb(  ps_ocl, PLQH )

    if bRealSpace:
        #FF_ref = mmff.evalGridFFAtPoints( ps, PLQH=PLQH, bSplit=False, nPBC=nPBC )
        FF_ref = mmff.evalGridFFAtPoints( ps, PLQH=PLQH, bSplit=True, nPBC=nPBC  )
        plt.subplot(2,1,2);
        if Emax is None: Emax = -FF_ref[:,3].min()*maxSc; 
        if Fmax is None: Fmax = -FF_ref[:,2].min()*maxSc;

    plt.figure(figsize=(5,10))
    plt.subplot(2,1,1);
    plt.plot( ts, fe_cl[:,3], '-k' , lw=1.0, label='OpenCL')
    plt.plot( ts, FFout[:,3], "-g" , lw=0.5, label="Etot_fit" )
    if bRealSpace:
        plt.plot( ts, FF_ref[:,3], ":k", lw=2.0, label="Etot_ref" )
        plt.plot( ts, (FFout[:,3]-FF_ref[:,3])*scErr, "-r", lw=0.5, label=("Etot_err*%.2f" %scErr) )
    plt.axhline(0.0, c="k", ls='--', lw=0.5)
    plt.ylim( -Emax, Emax )
    plt.legend()
    plt.subplot(2,1,2);
    plt.plot( ts, fe_cl[:,iax], '-k' , lw=1.0, label='OpenCL')
    plt.plot( ts, FFout[:,iax], "-g" , lw=0.5, label="Ftot_fit" )
    if bRealSpace:
        plt.plot( ts, FF_ref[:,2], ":k", lw=2.0, label="Ftot_ref" )
        plt.plot( ts, (FFout[:,2]-FF_ref[:,2])*scErr, "-r", lw=0.5, label=("Ftot_err*%.2f" %scErr) )
    plt.axhline(0.0, c="k", ls='--', lw=0.5)
    plt.ylim( -Fmax, Fmax )
    plt.legend()

    title_ = "iax="+str(iax)+"\n p0="+str(p0)+"\nnPBC="+str(nPBC)
    if ( title is not None ): title_=title+title_
    plt.suptitle( title_ )
    if ( bSaveFig ): plt.savefig( "test_gridFF_zcut.png", bbox_inches='tight' )
    
    #print( "ff.shape ", EFg.shape )
    #print( "test_gridFF() DONE" )
    print( "py======= test_gridFF() DONE" );
    #return EFg
