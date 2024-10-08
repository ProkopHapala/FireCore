import os
import numpy as np
import matplotlib.pyplot as plt

from . import utils as ut
from .. import atomicUtils as au
from .. import FunctionSampling as fu
from ..OCL.splines import OCLSplines

# =============  Functions

os.environ['PYOPENCL_CTX'] = '0'
ocl_splines = OCLSplines()

def test_gridFF_ocl( fname="./data/xyz/NaCl_1x1_L2.xyz", Element_Types_name="./data/ElementTypes.dat", bSymetrize=True, mode=6, dsamp=0.02, p0=[0.0,0.0,2.0], R0=3.5, E0=0.1, a=1.6, Q=0.4, H=0.0, scErr=100.0, iax=2, Emax=None, Fmax=None, maxSc=5.0, title=None, bSaveFig=True, bRefine=True, nPBC=[100,100,0], bRealSpace=False ):
    print( "py======= test_gridFF() START" );

    #Element_Types_name="/home/prokop/git/FireCore/tests/tMMFF/data/ElementTypes.dat"

    print( os.getcwd() )
    atoms = au.AtomicSystem(  fname=fname )
    if bSymetrize:
        atoms, ws = atoms.symmetrized()
        #print( "ws ",    ws    )
    REvdW = au.getVdWparams( atoms.atypes, fname=Element_Types_name )
    #print( "REvdW ", REvdW )
    if bSymetrize: REvdW[:,1] *= ws
    print( "REvdW ", REvdW )

    na = len(atoms.atypes)
    REQs=np.zeros( (na,4), dtype=np.float32 )
    xyzq=np.zeros( (na,4), dtype=np.float32 )

    xyzq[:,:3] = atoms.apos
    xyzq[:,3]  = atoms.qs
    REQs[:,0]  = REvdW[:,0]
    REQs[:,1]  = REvdW[:,1]
    REQs[:,2]  = atoms.qs
    REQs[:,3]  = 0.0

    V_Paul, V_Lond = ocl_splines.make_MorseFF( xyzq, REQs, nPBC=(4,4,0), dg=(0.1,0.1,0.1), lvec=atoms.lvec, grid_p0=(0.0,0.0,0.0), GFFParams=(0.1,1.5,0.0,0.0)  )

    print( "V_Paul.shape ", V_Paul.shape )
    plt.figure()
    plt.plot( V_Paul[:,0,0], label="V_Paul(0,0,z)" )
    plt.plot( V_Lond[:,0,0], label="V_Lond(0,0,z)" )
    plt.plot( V_Paul[:,20,20], label="V_Paul(20,20,z)" )
    plt.plot( V_Lond[:,20,20], label="V_Lond(20,20,z)" )
    plt.legend()
    #plt.yscale('log')
    plt.grid()

    ix=0
    iy=0
    iz=5
    replicate_count = 3  # Number of replications
    V_Paul_cut = V_Paul[iz,iy,:]; V_Paul_cut_rep = np.tile(V_Paul_cut, replicate_count)
    V_Lond_cut = V_Lond[iz,iy,:]; V_Lond_cut_rep = np.tile(V_Lond_cut, replicate_count)
    
    plt.figure()
    # Plot the replicated 1D cuts
    plt.plot( V_Paul_cut_rep, label=f"V_Paul(x,{iy},{iz}) - replicated")
    plt.plot( V_Lond_cut_rep, label=f"V_Lond(x,{iy},{iz}) - replicated")
    plt.legend()
    plt.grid()
    plt.show()
    plt.show()

    exit()

    #plotGridFF_1D( EFg, ix=20,iy=20 )
    PLQH  = ut.getPLQH( R0, E0, a, Q, H )

    ps,ts = ut.make_sample_points( p0, t0=0.0, tmax=10.0, dsamp=dsamp, iax=iax )    
    
    # ---- OpenCL
    g0   = (-2.0,-2.0,0.0)
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

    # plt.figure(figsize=(5,10))
    # plt.subplot(2,1,1);
    # plt.plot( ts, fe_cl[:,3], '-k' , lw=1.0, label='OpenCL')
    # plt.plot( ts, FFout[:,3], "-g" , lw=0.5, label="Etot_fit" )
    # if bRealSpace:
    #     plt.plot( ts, FF_ref[:,3], ":k", lw=2.0, label="Etot_ref" )
    #     plt.plot( ts, (FFout[:,3]-FF_ref[:,3])*scErr, "-r", lw=0.5, label=("Etot_err*%.2f" %scErr) )
    # plt.axhline(0.0, c="k", ls='--', lw=0.5)
    # plt.ylim( -Emax, Emax )
    # plt.legend()
    # plt.subplot(2,1,2);
    # plt.plot( ts, fe_cl[:,iax], '-k' , lw=1.0, label='OpenCL')
    # plt.plot( ts, FFout[:,iax], "-g" , lw=0.5, label="Ftot_fit" )
    # if bRealSpace:
    #     plt.plot( ts, FF_ref[:,2], ":k", lw=2.0, label="Ftot_ref" )
    #     plt.plot( ts, (FFout[:,2]-FF_ref[:,2])*scErr, "-r", lw=0.5, label=("Ftot_err*%.2f" %scErr) )
    # plt.axhline(0.0, c="k", ls='--', lw=0.5)
    # plt.ylim( -Fmax, Fmax )
    # plt.legend()

    # title_ = "iax="+str(iax)+"\n p0="+str(p0)+"\nnPBC="+str(nPBC)
    # if ( title is not None ): title_=title+title_
    # plt.suptitle( title_ )
    if ( bSaveFig ): plt.savefig( "test_gridFF_zcut.png", bbox_inches='tight' )
    
    #print( "ff.shape ", EFg.shape )
    #print( "test_gridFF() DONE" )
    print( "py======= test_gridFF() DONE" );
    #return EFg
