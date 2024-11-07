import os
import numpy as np
import matplotlib.pyplot as plt
import time

#exit()
#from . import utils as ut
from .. import atomicUtils as au
#from .. import FunctionSampling as fu
#exit()
from ..OCL.GridFF import GridFF_cl, GridShape
#exit()
#from .Ewald import compute_potential, plot1Dcut
from .utils import compute_potential, plot1Dcut

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

def Bspline_basis5(t):
    ws = np.zeros( (len(t),6) )
    inv6 = 1./6.;
    t2 = t*t;
    t3 = t2*t;
    t4 = t2*t2;
    t5 = t3*t2;                                                 
    ws[:,0]=  -0.008333333333333333*t5  +0.041666666666666666*t4  -0.08333333333333333*t3 +0.08333333333333333*t2  -0.041666666666666666*t   +0.008333333333333333;
    ws[:,1]=   0.041666666666666666*t5  -0.166666666666666666*t4  +0.16666666666666666*t3 +0.16666666666666666*t2  -0.416666666666666666*t   +0.216666666666666666;        
    ws[:,2]=  -0.083333333333333333*t5  +0.250000000000000000*t4                          -0.50000000000000000*t2                            +0.550000000000000000;  
    ws[:,3]=   0.083333333333333333*t5  -0.166666666666666666*t4  -0.16666666666666666*t3 +0.16666666666666666*t2  +0.416666666666666666*t   +0.216666666666666666;
    ws[:,4]=  -0.041666666666666666*t5  +0.041666666666666666*t4  +0.08333333333333333*t3 +0.08333333333333333*t2  +0.041666666666666666*t   +0.008333333333333333; 
    ws[:,5]=   0.008333333333333333*t5;
    return ws

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

def coulomb_brute_2D( atoms, kind='diag_z', bPlot=False, p0=[0.0,0.0,0.1], nPBC=(60,60,0) ):
    
    if kind=='diag_z':
        xs = np.arange( 0.0, 4.0-.0000001,  0.1 )
        zs = np.arange( 0.0, 40.0-.0000001, 0.1 )
        Xs,Zs = np.meshgrid( xs, zs )
        nx=len(xs); ny=len(zs)
        ps = np.zeros( (nx*ny,4), dtype=np.float32 )
        ps[:,0] = Xs.flatten()
        ps[:,1] = Xs.flatten()
        ps[:,2] = Zs.flatten()
        extent=[0.,4.*np.sqrt(2),0.,40.*np.sqrt(2)]
    elif kind=='xy':
        xs = np.arange( 0.0, 4.0-.0000001,  0.1 )
        Xs,Ys = np.meshgrid( xs, xs )
        nx=len(xs); ny=len(xs)
        ps = np.zeros( (nx*ny,4), dtype=np.float32 )
        ps[:,0] = Xs.flatten()
        ps[:,1] = Ys.flatten()
        ps[:,2] = p0[2]
        extent=[0.,4.,0.,4.]
    elif kind=='xz':
        xs = np.arange( 0.0, 4.0-.0000001,  0.1 )
        zs = np.arange( 0.0, 40.0-.0000001, 0.1 )
        Xs,Zs = np.meshgrid( xs, zs )
        nx=len(xs); ny=len(zs)
        ps = np.zeros( (nx*ny,4), dtype=np.float32 )
        ps[:,0] = Xs.flatten()
        ps[:,1] = p0[1]
        ps[:,2] = Zs.flatten()
        extent=[0.,4.,0.,40.]
    FEps = clgff.make_Coulomb_points( atoms, ps, nPBC=nPBC, Ls=[4.0,4.0,40.0], GFFParams=(0.00001, 1.5, 0.0, 0.0), bReturn=True)
    FE = FEps[:,3].reshape( ny,nx )
    #FE-=FE[-1,-1]
    if bPlot:
        vmin=FE.min()
        vmax=FE.max()
        vmax_=max( vmax, -vmin )
        print( "vmin,vmax ", vmin, vmax )
        plt.imshow( FE, origin="lower", cmap='bwr',  extent=extent, vmin=-vmax_, vmax=vmax_ )
        plt.colorbar()
    return FE

def coulomb_brute_1D( atoms, kind='z', p0=[0.0,0.0,2.0], bPlot=False, nPBC=(60,60,0) ):
        if kind=='z':
            ts = np.arange( 0.0, 40.0-.0000001, 0.1 )
            ps = np.zeros( (len(ts),4), dtype=np.float32 )
            ps[:,0] = p0[0]; 
            ps[:,1] = p0[1]
            ps[:,2] = p0[2]
            ps[:,2] = ts
        elif kind=='xy':
            ts = np.arange( 0.0, 4*4.0-.0000001, 0.1 )
            ps = np.zeros( (len(ts),4), dtype=np.float32 )
            ps[:,0] = p0[0]; 
            ps[:,1] = p0[1]
            ps[:,2] = p0[2]
            ps[:,0] = ts
            ps[:,1] = ts
        FEps = clgff.make_Coulomb_points( atoms, ps, nPBC=nPBC, Ls=[4.0,4.0,40.0], GFFParams=(0.00001, 1.5, 0.0, 0.0), bReturn=True)
        vmin=FEps.min()
        vmax=FEps.max()
        print( "coulomb_brute_1D() vmin,vmax ", vmin, vmax )
        if bPlot:
            plt.plot( ts, FEps[:,3], label=(f"V({kind}) "+str(p0)) )
        return FEps
        #plt.show()

def make_atoms_arrays( atoms=None, fname=None, bSymetrize=False, Element_Types_name="./data/ElementTypes.dat", bSqrtEvdw=True ): 
    #print( os.getcwd() )
    if atoms is None:
        atoms = au.AtomicSystem( fname=fname )
    if bSymetrize:
        na_before = len(atoms.atypes)
        atoms, ws = atoms.symmetrized()
        #print( "ws ",    ws    )

    print( "Qtot ", np.sum(atoms.qs)," Qabs ", np.sum(np.abs(atoms.qs)) )
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
    if bSqrtEvdw:
        REQs[:,1]  = np.sqrt(REvdW[:,1])
    else:
        REQs[:,1]  = REvdW[:,1]
    REQs[:,2]  = atoms.qs
    REQs[:,3]  = 0.0

    return xyzq, REQs, atoms


def plotTrjs( trjs, names ):
    plt.figure(figsize=(5,5))
    for i,trj in enumerate(trjs):
        #plt.plot( trj[:,0], trj[:,1], label=names[i] )
        plt.plot( trj[:,0], trj[:,1], label=names[i]+" |F|" )
        plt.plot( trj[:,0], trj[:,2], label=names[i]+" |E|" )
    plt.legend()
    plt.grid()
    plt.xlabel('iteration')
    plt.yscale('log')
    plt.title( "GridFF Bspline fitting error" )
    plt.show()


def test_gridFF_ocl( fname="./data/xyz/NaCl_1x1_L1.xyz", Element_Types_name="./data/ElementTypes.dat", job="PLQ", b2D=False, bSymetrize=False, bFit=True, save_name=None, z0=np.nan ):
    print( "py======= test_gridFF_ocl() START" );

    T00 = time.perf_counter()

    #Element_Types_name="/home/prokop/git/FireCore/tests/tMMFF/data/ElementTypes.dat"

    xyzq, REQs, atoms = make_atoms_arrays( fname=fname, bSymetrize=bSymetrize, Element_Types_name=Element_Types_name )

    if np.isnan(z0): z0 = xyzq[0,2].max()
    print( "test_gridFF_ocl() z0= ", z0 )

    grid = GridShape( dg=(0.1,0.1,0.1),  lvec=atoms.lvec)
    clgff.set_grid( grid )

    if job=="brute":
            
        nPBC=(60,60, 0)
        if b2D:
            plt.figure(figsize=(7,5))
            plt.subplot(1,2,1); coulomb_brute_2D( xyzq, kind='diag_z', bPlot=True )
            plt.subplot(1,2,2); coulomb_brute_2D( xyzq, kind='xy',     bPlot=True )
        else:
            coulomb_brute_1D( xyzq, kind='z',  p0=[0.0,0.0,2.0], bPlot=True )
            coulomb_brute_1D( xyzq, kind='z',  p0=[2.0,2.0,2.0], bPlot=True )
            coulomb_brute_1D( xyzq, kind='xy', p0=[0.0,0.0,2.0], bPlot=True )
        plt.show()
    
    elif job=="PLQ_lin":
        
        path = os.path.basename( fname )
        path = "./data/" + os.path.splitext( path )[0] + "/"
        print( "test_gridFF_ocl() path = ", path )
        
        # ----- Morse -----

        g0 = ( -grid.Ls[0]*0.5, -grid.Ls[1]*0.5, z0 )
        nPBC_mors = autoPBC(atoms.lvec,Rcut=20.0); print("autoPBC(nPBC_mors): ", nPBC_mors )
        FE_Paul, FE_Lond = clgff.make_MorseFF_f4( xyzq, REQs, nPBC=nPBC_mors, lvec=atoms.lvec, g0=g0, GFFParams=(0.1,1.5,0.0,0.0), bReturn=True )
        np.save( path+"FE_Paul.npy", FE_Paul )
        np.save( path+"FE_Lond.npy", FE_Lond )
    
        # ----- Coulomb -----

        Vcoul = clgff.makeCoulombEwald_slab( xyzq, niter=2, bTranspose=True )
        grid3D_shape = (clgff.gcl.g0, clgff.gcl.dg, clgff.gcl.ns)
        
        FE_Coul = clgff.sample3D_grid( clgff.V_Coul_buff, grid3D_shape )
        
        np.save( path+"FE_Coul.npy", FE_Coul )

        print( "SUMMARY:  FE_Paul.shape, FE_Lond.shape, FE_Coul.shape ", FE_Paul.shape, FE_Lond.shape, FE_Coul.shape )

        plt.figure( figsize=(16,8) )
        iChan = 0
        iSlice = 5
        
        plt.subplot(1,3,1); plt.imshow( FE_Paul[iSlice,:,:,iChan].transpose() ); plt.colorbar(); plt.title( f"FE_Paul[{iSlice},:,:,{iChan}]" );
        plt.subplot(1,3,2); plt.imshow( FE_Lond[iSlice,:,:,iChan].transpose() ); plt.colorbar(); plt.title( f"FE_Lond[{iSlice},:,:,{iChan}]" );
        #plt.subplot(1,3,1); plt.imshow( Vcoul[iSlice,:,:]        .transpose() );  plt.colorbar(); plt.title( f"V_Paul[{iSlice},:,:]" );
        #plt.subplot(1,3,2); plt.imshow( FE_Coul[:,:,iSlice,iChan].transpose() );  plt.colorbar(); plt.title( f"FE_Coul[:,:,{iSlice},{iChan}]" );
        plt.subplot(1,3,3); plt.imshow( FE_Coul[iSlice,:,:,iChan].transpose() );  plt.colorbar(); plt.title( f"FE_Coul[{iSlice},:,:,{iChan}]" );
        
    
    elif job=="PLQ":
        
        # --- Morse
        g0 = ( -grid.Ls[0]*0.5, -grid.Ls[1]*0.5, z0 )
        nPBC_mors = autoPBC(atoms.lvec,Rcut=20.0); print("autoPBC(nPBC_mors): ", nPBC_mors )
        clgff.make_MorseFF( xyzq, REQs, nPBC=nPBC_mors, lvec=atoms.lvec, g0=g0, GFFParams=(0.1,1.5,0.0,0.0), bReturn=False )
        V_Paul,trj_paul = clgff.fit3D( clgff.V_Paul_buff, nPerStep=50, nmaxiter=500, damp=0.15, bConvTrj=True ); #T_fit_P = time.perf_counter()
        V_Lond,trj_lond = clgff.fit3D( clgff.V_Lond_buff, nPerStep=50, nmaxiter=500, damp=0.15, bConvTrj=True ); #T_fit_ = time.perf_counter()

        # --- Coulomb
        Vcoul = clgff.makeCoulombEwald_slab( xyzq, niter=2, bSaveQgrid=True, bCheckVin=True, bCheckPoisson=True )

        #VcoulB,trj_coul = clgff.fit3D( clgff.V1_buff, nPerStep=100, nmaxiter=3000, damp=0.05, bConvTrj=True ); #T_fit_ = time.perf_counter()
        VcoulB,trj_coul = clgff.fit3D( clgff.V_Coul_buff, nPerStep=50, nmaxiter=500, damp=0.15, bConvTrj=True ); #T_fit_ = time.perf_counter()

        plotTrjs( [trj_paul,trj_lond,trj_coul], ["Paul", "Lond", "Coul"] )

        # -- total max,min over 3d grid of V_Paul,V_Lond,VcoulB 
        print( "Paul  min,max = ", V_Paul.min(), V_Paul.max() )
        print( "Lond  min,max = ", V_Lond.min(), V_Lond.max() ) 
        print( "Coul  min,max = ", Vcoul.min(),  Vcoul.max()  )  
        print( "CoulB min,max = ", VcoulB.min(), VcoulB.max() )

        axs=(1,2)
        minPaul = V_Paul.max( axis=axs )
        maxPaul = V_Paul.max( axis=axs )
        minLond = V_Lond.min( axis=axs )
        maxLond = V_Lond.max( axis=axs )
        minCoul = VcoulB.min( axis=axs )
        maxCoul = VcoulB.max( axis=axs )
        plt.plot( minPaul, label="V_Paul.min" )
        plt.plot( maxPaul, label="V_Paul.max" )
        plt.plot( minLond, label="V_Lond.min" )
        plt.plot( maxLond, label="V_Lond.max" )
        plt.plot( minCoul, label="V_Coul.min" )
        plt.plot( maxCoul, label="V_Coul.max" )
        plt.show()


        # V_Paul = V_Paul.max(axis=0)
        # V_Lond = V_Lond.max(axis=0)
        # VcoulB = VcoulB.max(axis=0)
        # V_Coul = VcoulB.max(axis=0)
        # V_max = np.maximum( V_Paul, V_Lond )
        # V_max = np.maximum( V_max, VcoulB )
        # print( "V_max.shape ", V_max.shape )

        if save_name=='double3':
            path = os.path.basename( fname )
            path = "./data/" + os.path.splitext( path )[0]
            print( "test_gridFF_ocl() path = ", path )
            if not os.path.exists( path ): os.makedirs( path )
            V_Paul = V_Paul.transpose( (2,1,0) )
            V_Lond = V_Lond.transpose( (2,1,0) )
            V_Coul = VcoulB.transpose( (2,1,0) )
            PLQ = np.zeros( V_Paul.shape + (3,) )
            PLQ[:,:,:,0] = V_Paul
            PLQ[:,:,:,1] = V_Lond
            PLQ[:,:,:,2] = V_Coul
            full_name = path+"/Bspline_PLQd_ocl.npy"; 
            print("test_gridFF_ocl() - save Morse to: ", full_name)
            np.save( full_name, PLQ )

        #cmap='plasma'
        #cmap='inferno'
        #cmap='magma'
        plt.figure( figsize=(16,8) )
        # plt.subplot(1,3,1); plt.imshow( V_Paul[:,:,0] ); plt.colorbar(); plt.title( "V_Paul[:,:,0]" );
        # plt.subplot(1,3,2); plt.imshow( V_Lond[:,:,0] ); plt.colorbar(); plt.title( "V_Lond[:,:,0]" );
        # plt.subplot(1,3,3); plt.imshow( V_Coul[:,:,0] ); plt.colorbar(); plt.title( "V_Coul[:,:,0]" );
        plt.subplot(1,3,1); plt.imshow( V_Paul[0,:,:].transpose() ); plt.colorbar(); plt.title( "V_Paul[0,:,:]" );
        plt.subplot(1,3,2); plt.imshow( V_Lond[0,:,:].transpose() ); plt.colorbar(); plt.title( "V_Lond[0,:,:]" );
        plt.subplot(1,3,3); plt.imshow( V_Coul[0,:,:].transpose() ); plt.colorbar(); plt.title( "V_Coul[0,:,:]" );
        # #plt.subplot(1,2,1); plt.imshow( V_Paul[0,:,:], cmap='bwr' ); plt.colorbar(); plt.title( "V_Paul[0,:,:]" );
        # #plt.subplot(1,2,2); plt.imshow( Vcoul [0,:,:], cmap='bwr' ); plt.colorbar(); plt.title( "Vcoul [0,:,:]" );
        # #plt.subplot(1,2,1); plt.imshow( V_Paul[:,:,0], cmap='bwr' ); plt.colorbar(); plt.title( "V_Paul[0,:,:]" );
        # #plt.subplot(1,2,1); plt.imshow( V_Lond[:,:,0], cmap='bwr' ); plt.colorbar(); plt.title( "V_Lond[0,:,:]" );
        # #plt.subplot(1,2,2); plt.imshow( Vcoul [:,:,0], cmap='bwr' ); plt.colorbar(); plt.title( "Vcoul [0,:,:]" );
        # plt.subplot(1,2,1); plt.imshow( Vcoul [5,:,:] ); plt.colorbar(); plt.title( "V_Lond[0,:,:]" );
        # plt.subplot(1,2,2); plt.imshow( Vcoul [:,:,0] ); plt.colorbar(); plt.title( "Vcoul [0,:,:]" );




        # --- Check the inverted axis of the V grid
        # plt.figure( figsize=(16,8) )
        # V = np.empty( [400,40,40], dtype=np.float32)
        # clgff.cl.enqueue_copy(clgff.queue, V, clgff.V2_buff)
        # plt.subplot(1,2,1); plt.imshow( V [-5,:,:] ); plt.colorbar(); plt.title( "V_[0,:,:]" );
        # plt.subplot(1,2,2); plt.imshow( V [200:,:,0] ); plt.colorbar(); plt.title( "V_ [0,:,:]" );
        # print( "V_Before[:10,0,0] \n", V[::-1,-1,-1][:10] )
        # print( "V_after[:10,0,0] \n", Vcoul[:10,0,0] )

        plt.show()
        #exit(0)



    elif job=="Ewald":
        
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
        #clgff.set_grid( grid )
        #Vcoul = clgff.makeCoulombEwald( xyzq )

        xyzq[:,2] += 0.1*4  # NOTE / TODO : This is strange, not sure why we need to shift the z-coordinate by 4 dg  

        Vcoul = clgff.makeCoulombEwald_slab( xyzq, niter=2 )

        if save_name is not None:
            VcoulB,trj = clgff.fit3D( clgff.V1_buff, nPerStep=10, nmaxiter=5000, damp=0.05, bConvTrj=True ); #T_fit_ = time.perf_counter()

            # NOTE: TODO: We need to cut Potential from oposite side of the cell ( we probably need to make CPU kernel for this )

            Vcoul = VcoulB

            plt.figure(figsize=(5,5))

            #for damp in [0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.50]:
            #for damp in [0.01,0.02,0.03,0.04,0.05,0.06]:
            # for damp in [0.04,0.05,0.06,0.07,0.08,0.09]:
            #     VcoulB,trj = clgff.fit3D( clgff.V1_buff, nPerStep=10, nmaxiter=5000, damp=damp, bConvTrj=True ); #T_fit_ = time.perf_counter()
            #     plt.plot( trj[:,0], trj[:,1], label=("VcoulB |F| damp=%g" %damp) )
            #     #plt.plot( trj[:,0], trj[:,2], label="VcoulB |E|" )

            plt.plot( trj[:,0], trj[:,1], label="VcoulB |F|" )
            plt.plot( trj[:,0], trj[:,2], label="VcoulB |E|" )
            plt.legend()
            plt.grid()
            plt.xlabel('iteration')
            plt.yscale('log')
            plt.title( "GridFF Bspline fitting error" )
            plt.show()

            path = os.path.basename( fname )
            path = "./data/" + os.path.splitext( path )[0]
            print( "test_gridFF_ocl() path = ", path )
            if not os.path.exists( path ):
                os.makedirs( path )
            if save_name=='double3':
                V_coul = VcoulB.transpose( (2,1,0) )
                PLQ = np.zeros( V_coul.shape + (3,) )
                PLQ[:,:,:,0] = 0
                PLQ[:,:,:,1] = 0
                PLQ[:,:,:,2] = V_coul
                full_name = path+"/Bspline_PLQd_ocl.npy"; print("test_gridFF_ocl() - save Morse to: ", full_name)
                np.save( full_name, PLQ )


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

        vmin = Vcoul.min()
        vmax = Vcoul.max()
        print( " Vcoul min, max ", vmin, vmax, )
        vmax = max( vmax, -vmin )
        vmax2=1e-5

        plt.figure( figsize=(20,7) )
        plt.subplot(1,6,1); plt.imshow( Vcoul[:,0,:       ], vmin=-vmax, vmax=vmax,   cmap='bwr' ); plt.colorbar(); plt.title( "Vcoul[:,0,:  ]" );
        plt.subplot(1,6,2); plt.imshow( Vcoul[:,:,20      ], vmin=-vmax, vmax=vmax,   cmap='bwr' ); plt.colorbar(); plt.title( "Vcoul[:,:,0  ]" );
        plt.subplot(1,6,3); plt.imshow( Vcoul[32,:,:      ],                          cmap='bwr' ); plt.colorbar(); plt.title( "Vcoul[32,:,:]" );
        # plt.subplot(1,6,4); plt.imshow( Vcoul[340,:,:     ], vmin=-vmax2, vmax=vmax2, cmap='bwr' ); plt.colorbar(); plt.title( "Vcoul[341,:,:]" );
        # plt.subplot(1,6,5); plt.imshow( Vcoul[100:,20,:],    vmin=-vmax2, vmax=vmax2, cmap='bwr' ); plt.colorbar(); plt.title( "Vcoul[200:,20,:]" );
        # plt.subplot(1,6,6); plt.imshow( Vcoul[100:,:,0 ],    vmin=-vmax2, vmax=vmax2, cmap='bwr' ); plt.colorbar(); plt.title( "Vcoul[200:,:,0]" );
        # plt.subplot(1,6,4); plt.imshow( Vcoul[340,:,:     ],   cmap='bwr' ); plt.colorbar(); plt.title( "Vcoul[341,:,:]"  );
        # plt.subplot(1,6,5); plt.imshow( Vcoul[200:320,20,:],   cmap='bwr' ); plt.colorbar(); plt.title( "Vcoul[200:,20,:]");
        # plt.subplot(1,6,6); plt.imshow( Vcoul[200:320,:,0 ],   cmap='bwr' ); plt.colorbar(); plt.title( "Vcoul[200:,:,0]" );

        # plt.subplot(1,6,5); plt.imshow( Vcoul[:,0,:],   cmap='bwr' ); plt.colorbar(); plt.title( "Vcoul[200:,20,:]" );
    
        # Vbrute = coulomb_brute_2D( xyzq, kind='xz',     bPlot=False )

        # plt.subplot(1,6,6); plt.imshow( Vbrute, cmap='bwr' ); plt.colorbar(); plt.title( "Vbrute[200:300,:,0]" );

        xyzq_sym, _, _ = make_atoms_arrays( atoms=atoms, bSymetrize=True, Element_Types_name=Element_Types_name )
        Vbrute         = coulomb_brute_1D( xyzq_sym, kind='z', p0=[0.0,0.0,0.0], bPlot=False )

        plt.figure( figsize=(7,5) )
        ysc = 0.03
        V_ref = Vbrute[::-1,3] - Vbrute[100,3]
        plt.plot( V_ref   ,':', lw=1.5, label="Vbrute" )
        plt.plot( Vcoul[:,0,0], lw=0.5, label="Vcoul"  )
        plt.ylim(-ysc,ysc)
        plt.legend()
        plt.grid()
        plt.show()

    elif job=="Morse":
    
        nPBC = autoPBC(atoms.lvec,Rcut=20.0); print("autoPBC: ", nPBC )

        if bFit:
            #g0 = ( -grid.Ls[0]*0.5, -grid.Ls[1]*0.5, p0[2] ) 
            g0 = ( -grid.Ls[0]*0.5, -grid.Ls[1]*0.5, 0.0 )
            #g0 = ( -grid.Ls[0]*0.5, -grid.Ls[1]*0.5, -p0[2] ) 
            #T0 = time.perf_counter()
            clgff.make_MorseFF( xyzq, REQs, nPBC=nPBC, lvec=atoms.lvec, g0=g0, GFFParams=(0.1,1.5,0.0,0.0), bReturn=False )
            #T_morse = time.perf_counter()

            # if bDebug:
            #     sh = clgff.gsh.ns[::-1]
            #     VPaul_ref = np.zeros( sh, dtype=np.float32); clgff.cl.enqueue_copy(clgff.queue, VPaul_ref, clgff.V_Paul_buff ); clgff.queue.finish()
            #     VLond_ref = np.zeros( sh, dtype=np.float32); clgff.cl.enqueue_copy(clgff.queue, VLond_ref, clgff.V_Lond_buff ); clgff.queue.finish()
            #     plt.figure(figsize=(10,5))
            #     #plt.subplot(1,2,1); plt.imshow( V_Paul[:,:,1] ); plt.colorbar()
            #     #plt.subplot(1,2,2); plt.imshow( V_Lond[:,:,1] ); plt.colorbar()
            #     plt.subplot(1,2,1); plt.imshow( VPaul_ref[1,:,:] ); plt.colorbar(); plt.title("VPaul_ref")
            #     plt.subplot(1,2,2); plt.imshow( VLond_ref[1,:,:] ); plt.colorbar(); plt.title("VLond_ref")

            nPerStep=50
            #nPerStep=5

            V_Paul,trj_paul = clgff.fit3D( clgff.V_Paul_buff, nPerStep=nPerStep, bConvTrj=True ); #T_fit_P = time.perf_counter()
            V_Lond,trj_lond = clgff.fit3D( clgff.V_Lond_buff, nPerStep=nPerStep, bConvTrj=True ); #T_fit_ = time.perf_counter()

            # if bDebug:
            #         plt.figure(figsize=(10,5))
            #         #plt.subplot(1,2,1); plt.imshow( V_Paul[:,:,1] ); plt.colorbar()
            #         #plt.subplot(1,2,2); plt.imshow( V_Lond[:,:,1] ); plt.colorbar()
            #         plt.subplot(1,2,1); plt.imshow( V_Paul[1,:,:] ); plt.colorbar(); plt.title("V_Paul_fit")
            #         plt.subplot(1,2,2); plt.imshow( V_Lond[1,:,:] ); plt.colorbar(); plt.title("V_Lond_fit")

            if save_name is not None:
                path = os.path.basename( fname )
                path = "./data/" + os.path.splitext( path )[0]
                print( "test_gridFF_ocl() path = ", path )
                if not os.path.exists( path ):
                    os.makedirs( path )
                if save_name=='double3':
                    V_Paul = V_Paul.transpose( (2,1,0) )
                    V_Lond = V_Lond.transpose( (2,1,0) )
                    PLQ = np.zeros( V_Paul.shape + (3,) )
                    PLQ[:,:,:,0] = V_Paul
                    PLQ[:,:,:,1] = V_Lond
                    #PLQ[:,:,2] = 0.0 # electrostatic
                    full_name = path+"/Bspline_PLQd_ocl.npy"; print("test_gridFF_ocl() - save Morse to: ", full_name)
                    np.save( full_name, PLQ )

            # --- cdamp scan
            #plt.figure(figsize=(5,5))
            # for damp in [0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.50]:
            #     print(" damp =  ", damp )
            #     #V_Paul,trj_paul = clgff.fit3D( clgff.V_Paul_buff, damp=damp, bConvTrj=True )
            #     #plt.plot( trj_paul[:,0], trj_paul[:,1], label=("Paul |F| %8.4f" %damp) )
            #     V_Lond,trj_lond = clgff.fit3D( clgff.V_Lond_buff, damp=damp, bConvTrj=True )
            #     plt.plot( trj_lond[:,0], trj_lond[:,1], label=("Lond |F| %8.4f" %damp))

            # --- dt scan
            # plt.figure(figsize=(5,5))
            # #for dt in [0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45]:
            # #for dt in [0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.67]:
            # for dt in [0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.00]:
            #     print(" dt =  ", dt )
            #     # V_Paul,trj_paul = clgff.fit3D( clgff.V_Paul_buff, dt=dt, damp=0.30, bConvTrj=True )
            #     # plt.plot( trj_paul[:,0], trj_paul[:,1], label=("Paul |F| dt=%8.4f" %dt) )
            #     V_Lond,trj_lond = clgff.fit3D( clgff.V_Lond_buff, dt=dt, damp=0.15, bConvTrj=True )
            #     plt.plot( trj_lond[:,0], trj_lond[:,1], label=("Lond |F| dt=%8.4f" %dt))

            plt.figure(figsize=(5,5))
            plt.plot( trj_paul[:,0], trj_paul[:,1], label="Paul |F|" )
            plt.plot( trj_paul[:,0], trj_paul[:,2], label="Paul |E|" )
            plt.plot( trj_lond[:,0], trj_lond[:,1], label="Lond |F|" )
            plt.plot( trj_lond[:,0], trj_lond[:,2], label="Lond |E|" )

            plt.legend()
            plt.grid()
            plt.xlabel('iteration')
            plt.yscale('log')
            plt.title( "GridFF Bspline fitting error" )

            #plt.show()

        # else:
        #     V_Paul, V_Lond = clgff.make_MorseFF( xyzq, REQs, nPBC=nPBC, lvec=atoms.lvec, g0=(0.0,0.0,0.0), GFFParams=(0.1,1.5,0.0,0.0)  )
        #     print( "V_Paul.shape ", V_Paul.shape )
        #     plt.figure() 
        #     plt.plot( V_Paul[:,0,0],   label="V_Paul(0,0,z)" )
        #     plt.plot( V_Lond[:,0,0],   label="V_Lond(0,0,z)" )
        #     plt.plot( V_Paul[:,20,20], label="V_Paul(20,20,z)" )
        #     plt.plot( V_Lond[:,20,20], label="V_Lond(20,20,z)" )
        #     plt.legend()
        #     #plt.yscale('log')
        #     plt.grid()
        #     ix=0
        #     iy=0
        #     iz=5
        #     replicate_count = 3  # Number of replications
        #     V_Paul_cut = V_Paul[iz,iy,:]; V_Paul_cut_rep = np.tile(V_Paul_cut, replicate_count)
        #     V_Lond_cut = V_Lond[iz,iy,:]; V_Lond_cut_rep = np.tile(V_Lond_cut, replicate_count)
        #     plt.figure()
        #     # Plot the replicated 1D cuts
        #     plt.plot( V_Paul_cut_rep, label=f"V_Paul(x,{iy},{iz}) - replicated")
        #     plt.plot( V_Lond_cut_rep, label=f"V_Lond(x,{iy},{iz}) - replicated")
        #     plt.legend()
        #     plt.grid()
        #     #plt.show()
        #     #plt.show()
    
    

        
    print( "py======= test_gridFF() DONE" );

