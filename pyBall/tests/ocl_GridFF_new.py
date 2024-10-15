import os
import numpy as np
import matplotlib.pyplot as plt

#from . import utils as ut
from .. import atomicUtils as au
#from .. import FunctionSampling as fu
from ..OCL.GridFF import GridFF_cl, GridShape
from .Ewald import compute_potential, plot1Dcut

# =============  Functions

mmff = None

os.environ['PYOPENCL_CTX'] = '0'
clgff = GridFF_cl()

def try_load_mmff():
    global mmff
    if mmff is None:
        from .. import MMFF as mmff_
        mmff = mmff_


def autoPBC(lvec,Rcut=20.0,mask=(1,1,0)):
    nPBC = [0,0,0]
    for i in range(3):
        if mask[i]>0:
            L = np.linalg.norm(lvec[i])  # Length of the cell in each direction
            #print( "L ",L,"[",i,"] Rcut ", Rcut  )
            nPBC[i] = int(Rcut/L)+1
    return tuple(nPBC)

def test_Ewald( apos, qs, ns=[100,100,100], dg=(0.1,0.1,0.1), nPBC=[30,30,30], pos0=None, scErr=100.0, order=2, bPython=True, bOCL=True, bPlotPy=False,  bPlotOcl=True, bOMP=False, nBlur=0, cSOR=0.0, cV=0.5, yrange=None, bPlot1D=True , bSlab=False, z_slab=None ):
    try_load_mmff()

    apos = apos.copy()
    print( "apos ", apos )
    Ls = [ ns[0]*dg[0], ns[1]*dg[1], ns[2]*dg[2] ]                                # lenghts along each axis  
    if pos0 is None: pos0=np.array( [ -0.5*Ls[0], -0.5*Ls[1], -0.5*Ls[2] ] )                       # origin of grid - half of the grid size
    mmff.setupEwaldGrid( ns, dg=dg, pos0=pos0 )                                   # setup grid in C++ for projection
    dens = mmff.projectAtomsEwaldGrid( apos, qs, ns=ns, order=order )             # project atoms to grid using C++ 

    # perhaps we should find better reference than direct sum in reals space
    #  * perhaps using [Madelung Constant](https://en.wikipedia.org/wiki/Madelung_constant) for some known crystals

    Qtot = np.abs(dens).sum();        print("Qtot = ", Qtot, np.abs(qs).sum() )   # check if the total charge is conserved in projection ( should be same as sum of qs )
    
    # set voxel through which we cut the 3D grid Vg in debug plots
    ix0=ns[0]//2
    iy0=ns[1]//2
    iz0=ns[2]//2
    i0 = (ix0,iy0,iz0)

    if bPython:
        Vg,(density_fft, Vw, ker_w) = compute_potential(dens, dg, bInternals=True )

        if bPlotPy:
            # dens2d =  np.sum( dens, axis=iax )
            plt.figure(figsize=(15,10))
            extent  =[  -Ls[0]*0.5  , +Ls[0]*0.5,    -Ls[1]*0.5  ,+Ls[1]*0.5      ]
            extent_w=[  -np.pi/dg[0], +np.pi/dg[0],  -np.pi/dg[1],+np.pi/dg[1]    ]
            # --- plot in real-space (cartesian) coordinates
            plt.subplot(2,3,1  ); plt.imshow( dens [iz0,:,:],                        cmap='bwr', extent=extent  , origin="lower" ); plt.colorbar(); plt.title("Charge Density" )
            plt.subplot(2,3,2  ); plt.imshow( ker_w[iz0,:,:],                        cmap='bwr', extent=extent_w, origin="lower" ); plt.colorbar(); plt.title("ker_w" )
            plt.subplot(2,3,3  ); plt.imshow( Vg   [iz0,:,:],   vmin=-1.0,vmax=1.0,  cmap='bwr', extent=extent  , origin="lower" ); plt.colorbar(); plt.title("V ewald " )
            # --- plot in pixel-coordinates
            plt.subplot(2,3,3+1); plt.imshow( dens [iz0,:,:],                        cmap='bwr'               , origin="lower" ); plt.colorbar(); plt.title("Charge Density" )
            plt.subplot(2,3,3+2); plt.imshow( ker_w[iz0,:,:],                        cmap='bwr'               , origin="lower" ); plt.colorbar(); plt.title("ker_w" )
            plt.subplot(2,3,3+3); plt.imshow( Vg   [iz0,:,:],   vmin=-1.0,vmax=1.0,  cmap='bwr'               , origin="lower" ); plt.colorbar(); plt.title("V ewald" )
    
    if bOCL:
        grid = GridShape( dg=dg,  Ls=Ls, g0=tuple(pos0) )
        clgff.set_grid( grid )
        xyzq = np.zeros( (len(apos),4), dtype=np.float32 )
        xyzq[:,0] = apos[:,0]
        xyzq[:,1] = apos[:,1]
        xyzq[:,2] = apos[:,2]
        xyzq[:,3] = qs
        Vocl = clgff.makeCoulombEwald( xyzq )

        print( "Vocl min,max", Vocl.min(), Vocl.max() )

        if bPlotOcl:
            # dens2d =  np.sum( dens, axis=iax )
            plt.figure(figsize=(10,5))
            extent  =[  -Ls[0]*0.5  , +Ls[0]*0.5,    -Ls[1]*0.5  ,+Ls[1]*0.5      ]
            plt.subplot(1,2,1); plt.imshow( Vocl[iz0,:,:],   vmin=-1.0,vmax=1.0,  cmap='bwr', extent=extent  , origin="lower" ); plt.colorbar(); plt.title("V ewald " )
            plt.subplot(1,2,2); plt.imshow( Vocl[:,iy0,:],   vmin=-1.0,vmax=1.0,  cmap='bwr'                 , origin="lower" ); plt.colorbar(); plt.title("V ewald" )

    # ---- 1D debug plots
    if bPlot1D:
        plt.figure(figsize=(15,5))
        plt.subplot(1,3,1); plot1Dcut( apos, qs, Vocl, Vref=Vg, i0=i0, dg=dg, Ls=Ls, iax=0, nPBC=nPBC, scErr=scErr )
        plt.subplot(1,3,2); plot1Dcut( apos, qs, Vocl, Vref=Vg, i0=i0, dg=dg, Ls=Ls, iax=1, nPBC=nPBC, scErr=scErr )
        plt.subplot(1,3,3); plot1Dcut( apos, qs, Vocl, Vref=Vg, i0=i0, dg=dg, Ls=Ls, iax=2, nPBC=nPBC, scErr=scErr )
        #plt.subplot(1,3,3); plot1Dcut( Vg, i0=i0, dg=dg, lvec=lvec, iax=2, nPBC=nPBC, scErr=scErr )
        if yrange is not None:
            plt.subplot(1,3,1); plt.ylim(yrange)
            plt.subplot(1,3,2); plt.ylim(yrange)
            plt.subplot(1,3,3); plt.ylim(yrange)
        #plt.tight_layout()
        #plt.show()
        #super title
        plt.suptitle( "order=" + str(order) + " nBlur=" + str(nBlur) + " cSOR=" + str(cSOR) + " cV=" + str(cV) )


def test_gridFF_ocl( fname="./data/xyz/NaCl_1x1_L1.xyz", Element_Types_name="./data/ElementTypes.dat", bSymetrize=False, bMorse=True, bEwald=False, bFit=True,  mode=6, dsamp=0.02, p0=[0.0,0.0,2.0], R0=3.5, E0=0.1, a=1.6, Q=0.4, H=0.0, scErr=100.0, iax=2, Emax=None, Fmax=None, maxSc=5.0, title=None, bSaveFig=True, bRefine=True, nPBC=[100,100,0], bRealSpace=False, save_name=None ):
    print( "py======= test_gridFF_ocl() START" );

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

    if bEwald:
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
    
    if bMorse:
        nPBC = autoPBC(atoms.lvec,Rcut=20.0); print("autoPBC: ", nPBC )

        if bFit:
            clgff.make_MorseFF( xyzq, REQs, nPBC=nPBC, lvec=atoms.lvec, g0=(0.0,0.0,0.0), GFFParams=(0.1,1.5,0.0,0.0), bReturn=False )
            V_Paul,trj_paul = clgff.fit3D( clgff.V_Paul_buff , bConvTrj=True )
            V_Lond,trj_lond = clgff.fit3D( clgff.V_Lond_buff , bConvTrj=True )

            if save_name is not None:
                path = os.path.basename( fname )
                path = "./data/" + os.path.splitext( path )[0]
                print( "test_gridFF_ocl() path = ", path )
                if not os.path.exists( path ):
                    os.makedirs( path )
                if save_name=='double3':
                    PLQ = np.zeros( V_Paul.shape + (3,) )
                    PLQ[:,:,:,0] = V_Paul
                    PLQ[:,:,:,1] = V_Lond
                    #PLQ[:,:,2] = 0.0 # electrostatic
                    np.save( path+"/Bspline_PLQd.npy", PLQ )

            # --- cdamp scan
            # for damp in [0.05,0.10,0.15,0.20,0.25]:
            #     print(" damp =  ", damp )
            #     #V_Paul,trj_paul = clgff.fit3D( clgff.V_Paul_buff, damp=damp, bConvTrj=True )
            #     #plt.plot( trj_paul[:,0], trj_paul[:,1], label=("Paul |F| %8.4f" %damp) )
            #     V_Lond,trj_lond = clgff.fit3D( clgff.V_Lond_buff, damp=damp, bConvTrj=True )
            #     plt.plot( trj_lond[:,0], trj_lond[:,1], label=("Lond |F| %8.4f" %damp))

            # --- dt scan
            #for dt in [0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45]:
            # for dt in [0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.67]:
            #     print(" dt =  ", dt )
            #     V_Paul,trj_paul = clgff.fit3D( clgff.V_Paul_buff, dt=dt, damp=0.15, bConvTrj=True )
            #     plt.plot( trj_paul[:,0], trj_paul[:,1], label=("Paul |F| dt=%8.4f" %dt) )
            #     #V_Lond,trj_lond = clgff.fit3D( clgff.V_Lond_buff, dt=dt, damp=0.15, bConvTrj=True )
            #     #plt.plot( trj_lond[:,0], trj_lond[:,1], label=("Lond |F| dt=%8.4f" %dt))

            plt.plot( trj_paul[:,0], trj_paul[:,1], label="Paul |F|" )
            plt.plot( trj_paul[:,0], trj_paul[:,2], label="Paul |E|" )
            plt.plot( trj_lond[:,0], trj_lond[:,1], label="Lond |F|" )
            plt.plot( trj_lond[:,0], trj_lond[:,2], label="Lond |E|" )
            plt.legend()
            plt.grid()
            plt.xlabel('iteration')
            plt.yscale('log')
            plt.title( "GridFF Bspline fitting error" )
            plt.show()

        else:

            V_Paul, V_Lond = clgff.make_MorseFF( xyzq, REQs, nPBC=nPBC, lvec=atoms.lvec, g0=(0.0,0.0,0.0), GFFParams=(0.1,1.5,0.0,0.0)  )

            print( "V_Paul.shape ", V_Paul.shape )
            plt.figure() 
            plt.plot( V_Paul[:,0,0],   label="V_Paul(0,0,z)" )
            plt.plot( V_Lond[:,0,0],   label="V_Lond(0,0,z)" )
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
        
    print( "py======= test_gridFF() DONE" );

