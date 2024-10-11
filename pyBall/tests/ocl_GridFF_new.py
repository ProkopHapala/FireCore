import os
import numpy as np
import matplotlib.pyplot as plt

#from . import utils as ut
from .. import atomicUtils as au
#from .. import FunctionSampling as fu
from ..OCL.GridFF import GridFF_cl, GridShape

# =============  Functions

os.environ['PYOPENCL_CTX'] = '0'
clgff = GridFF_cl()

def autoPBC(lvec,Rcut=20.0,mask=(1,1,0)):
    nPBC = [0,0,0]
    for i in range(3):
        if mask[i]>0:
            L = np.linalg.norm(lvec[i])  # Length of the cell in each direction
            #print( "L ",L,"[",i,"] Rcut ", Rcut  )
            nPBC[i] = int(Rcut/L)+1
    return tuple(nPBC)

def test_gridFF_ocl( fname="./data/xyz/NaCl_1x1_L1.xyz", Element_Types_name="./data/ElementTypes.dat", bSymetrize=False, mode=6, dsamp=0.02, p0=[0.0,0.0,2.0], R0=3.5, E0=0.1, a=1.6, Q=0.4, H=0.0, scErr=100.0, iax=2, Emax=None, Fmax=None, maxSc=5.0, title=None, bSaveFig=True, bRefine=True, nPBC=[100,100,0], bRealSpace=False ):
    print( "py======= test_gridFF() START" );

    #Element_Types_name="/home/prokop/git/FireCore/tests/tMMFF/data/ElementTypes.dat"

    print( os.getcwd() )
    atoms = au.AtomicSystem(  fname=fname )
    if bSymetrize:
        na_before = len(atoms.atypes)
        atoms, ws = atoms.symmetrized()
        #print( "ws ",    ws    )
    REvdW = au.getVdWparams( atoms.atypes, fname=Element_Types_name )
    #print( "REvdW ", REvdW )
    if bSymetrize: 
        REvdW[:,1] *= ws
        print( "n_atoms (symetrized): ", len(atoms.atypes)," before symmertization: ", na_before )
    #print( "REvdW ", REvdW )

    na = len(atoms.atypes)
    REQs=np.zeros( (na,4), dtype=np.float32 )
    xyzq=np.zeros( (na,4), dtype=np.float32 )

    xyzq[:,:3] = atoms.apos
    xyzq[:,3]  = atoms.qs
    REQs[:,0]  = REvdW[:,0]
    REQs[:,1]  = REvdW[:,1]
    REQs[:,2]  = atoms.qs
    REQs[:,3]  = 0.0


    #---- Test Poisson
    
    grid = GridShape( dg=(0.1,0.1,0.1),  lvec=atoms.lvec)
    clgff.set_grid( grid )
    Vgrid = clgff.makeCoulombEwald( xyzq )

    # xline = Vgrid.sum(axis=(1,2)); plt.plot( xline )
    # yline = Vgrid.sum(axis=(0,2)); plt.plot( yline )
    # zline = Vgrid.sum(axis=(0,1)); plt.plot( zline )

    # plt.figure( figsize=(25,5) )
    # plt.subplot(1,6,1); plt.imshow( Vgrid[166,:,:], cmap='bwr' ); plt.colorbar()
    # plt.subplot(1,6,2); plt.imshow( Vgrid[167,:,:], cmap='bwr' ); plt.colorbar()
    # plt.subplot(1,6,3); plt.imshow( Vgrid[168,:,:], cmap='bwr' ); plt.colorbar()
    # plt.subplot(1,6,4); plt.imshow( Vgrid[169,:,:], cmap='bwr' ); plt.colorbar()
    # plt.subplot(1,6,5); plt.imshow( Vgrid[170,:,:], cmap='bwr' ); plt.colorbar()
    # plt.subplot(1,6,6); plt.imshow( Vgrid[171,:,:], cmap='bwr' ); plt.colorbar()

    plt.figure( figsize=(15,5) )
    plt.subplot(1,3,1); plt.imshow( Vgrid[170,:,:], cmap='bwr' ); plt.colorbar(); plt.title( "Vgrid[170,:,:]" );
    plt.subplot(1,3,2); plt.imshow( Vgrid[:,0,:  ], cmap='bwr' ); plt.colorbar(); plt.title( "Vgrid[:,0,:  ]" );
    plt.subplot(1,3,3); plt.imshow( Vgrid[:,:,0  ], cmap='bwr' ); plt.colorbar(); plt.title( "Vgrid[:,:,0  ]" );

    return
    
    #---- Test Charge-to-grid projection
    """
    #xyzq[:,2]=0.0

    Qgrid = clgff.project_atoms_on_grid_quintic_pbc(xyzq, dg=(0.1,0.1,0.1),  lvec=atoms.lvec )

    plt.figure( figsize=(25,5) )
    plt.subplot(1,6,1); plt.imshow( Qgrid[166,:,:], cmap='bwr' ); plt.colorbar()
    plt.subplot(1,6,2); plt.imshow( Qgrid[167,:,:], cmap='bwr' ); plt.colorbar()
    plt.subplot(1,6,3); plt.imshow( Qgrid[168,:,:], cmap='bwr' ); plt.colorbar()
    plt.subplot(1,6,4); plt.imshow( Qgrid[169,:,:], cmap='bwr' ); plt.colorbar()
    plt.subplot(1,6,5); plt.imshow( Qgrid[170,:,:], cmap='bwr' ); plt.colorbar()
    plt.subplot(1,6,6); plt.imshow( Qgrid[171,:,:], cmap='bwr' ); plt.colorbar()
    #plt.subplot(1,6,6); plt.imshow( Qgrid[172,:,:], cmap='bwr' ); plt.colorbar()

    #print( Qgrid[:,0,0], Qgrid[0,:,0], Qgrid[0,0,:], )
    # for iz in range(10):
    #     plt.figure()
    #     iiz=iz*3
    #     plt.imshow( Qgrid[:,:,iiz] )
    #     plt.title( f"Qgrid[{iz}] {iz*3}" )

    return
    """
    #---- Test Morse
    """
    nPBC = autoPBC(atoms.lvec,Rcut=20.0); print("autoPBC: ", nPBC )

    V_Paul, V_Lond = clgff.make_MorseFF( xyzq, REQs, nPBC=nPBC, lvec=atoms.lvec, g0=(0.0,0.0,0.0), GFFParams=(0.1,1.5,0.0,0.0)  )

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
    """
    print( "py======= test_gridFF() DONE" );

