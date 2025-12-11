import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from . import utils as ut
from .. import MMFF as mmff

COULOMB_CONST  =    14.3996448915 

def test_vs_direct( apos, qs, ns=[100,100,100], dg=[0.1,0.1,0.1], nPBC=[30,30,30], pos0=None, scErr=100.0, order=2, bPython=False, bPlotDebug=True, bOMP=False, nBlur=0, cSOR=0.0, cV=0.5, yrange=None, bPlot1D=True , bSlab=False, z_slab=None ):
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
    iz0=ns[2]//2 + 10
    i0 = (ix0,iy0,iz0)

    if bPython:
        Vg,(density_fft, Vw, ker_w) = ut.compute_potential(dens, dg, bInternals=True )

        if bPlotDebug:
            # ---- 2D debug plots
            # dens2d =  np.sum( dens, axis=iax )
            plt.figure(figsize=(15,10))
            extent  =[  -Ls[0]*0.5  , +Ls[0]*0.5,    -Ls[1]*0.5  ,+Ls[1]*0.5      ]
            extent_w=[  -np.pi/dg[0], +np.pi/dg[0],  -np.pi/dg[1],+np.pi/dg[1]    ]
            # --- plot in real-space (cartesian) coordinates
            plt.subplot(2,3,1  ); plt.imshow( dens [iz0,:,:],                        cmap='bwr', extent=extent  , origin="lower" ); plt.colorbar(); plt.title("Charge Density" )
            plt.subplot(2,3,2  ); plt.imshow( ker_w[iz0,:,:],                        cmap='bwr', extent=extent_w, origin="lower" ); plt.colorbar(); plt.title("ker_w" )
            plt.subplot(2,3,3  ); plt.imshow( Vg   [iz0,:,:],   vmin=-1.0,vmax=1.0,  cmap='bwr', extent=extent  , origin="lower" ); plt.colorbar(); plt.title("V ewald" )
            # --- plot in pixel-coordinates
            plt.subplot(2,3,3+1); plt.imshow( dens [iz0,:,:],                        cmap='bwr'               , origin="lower" ); plt.colorbar(); plt.title("Charge Density" )
            plt.subplot(2,3,3+2); plt.imshow( ker_w[iz0,:,:],                        cmap='bwr'               , origin="lower" ); plt.colorbar(); plt.title("ker_w" )
            plt.subplot(2,3,3+3); plt.imshow( Vg   [iz0,:,:],   vmin=-1.0,vmax=1.0,  cmap='bwr'               , origin="lower" ); plt.colorbar(); plt.title("V ewald" )

    else:
        nz_slab=-1
        if bSlab: 
            if z_slab is None: z_slab = Ls[2] 
            nz_slab = int(z_slab/dg[2])
        Vg = mmff.EwaldGridSolveLaplace(  dens, nz_slab=nz_slab, nBlur=nBlur, cSOR=cSOR, cV=cV, bOMP=bOMP )
        #Vg = mmff.EwaldGridSolveLaplace( dens  )
        
        # dV = dg[0]*dg[1]*dg[2]
        # ntot = ns[0]*ns[1]*ns[2]
        # scEwald = COULOMB_CONST * 4.0* np.pi / (ntot*dV)  
        # print("Ewald(C++) mmff.EwaldGridSolveLaplace() scEwald = ", scEwald)
        # Vg *= scEwald

    # ---- 1D debug plots
    if bPlot1D:
        plt.figure(figsize=(15,5))
        plt.subplot(1,3,1); ut.plot1Dcut( apos, qs, Vg, i0=i0, dg=dg, Ls=Ls, pos0=pos0,iax=0, nPBC=nPBC, scErr=scErr, mmff=mmff )
        plt.subplot(1,3,2); ut.plot1Dcut( apos, qs, Vg, i0=i0, dg=dg, Ls=Ls, pos0=pos0, iax=1, nPBC=nPBC, scErr=scErr, mmff=mmff )
        plt.subplot(1,3,3); ut.plot1Dcut( apos, qs, Vg, i0=i0, dg=dg, Ls=Ls, pos0=pos0,iax=2, nPBC=nPBC, scErr=scErr, mmff=mmff )
        #plt.subplot(1,3,3); plot1Dcut( Vg, i0=i0, dg=dg, lvec=lvec, iax=2, nPBC=nPBC, scErr=scErr )
        if yrange is not None:
            plt.subplot(1,3,1); plt.ylim(yrange)
            plt.subplot(1,3,2); plt.ylim(yrange)
            plt.subplot(1,3,3); plt.ylim(yrange)
        #plt.tight_layout()
        #plt.show()
        #super title
        plt.suptitle( "order=" + str(order) + " nBlur=" + str(nBlur) + " cSOR=" + str(cSOR) + " cV=" + str(cV) )

def test_project_dens( apos, qs, ns=[100,100,100], dg=[0.1,0.1,0.1], pos0=None, order=2, yrange=None ):
    apos = apos.copy()
    #print( "apos ", apos )
    Ls = [ ns[0]*dg[0], ns[1]*dg[1], ns[2]*dg[2] ]                                # lenghts along each axis  
    if pos0 is None: pos0=np.array( [ -0.5*Ls[0], -0.5*Ls[1], -0.5*Ls[2] ] )      # origin of grid - half of the grid size
    mmff.setupEwaldGrid( ns, dg=dg, pos0=pos0 )                                   # setup grid in C++ for projection
    dens = mmff.projectAtomsEwaldGrid( apos, qs, ns=ns, order=order )             # project atoms to grid using C++ 
    Qtot = np.abs(dens).sum();        print("Qtot = ", Qtot, np.abs(qs).sum() )   # check if the total charge is conserved in projection ( should be same as sum of qs )
    ix0=ns[0]//2
    iy0=ns[1]//2
    iz0=ns[2]//2
    plt.figure(figsize=(5,5))
    extent  =[  pos0[0], pos0[0]+Ls[0],   pos0[1], pos0[1]+Ls[1],]
    vmax=+1; vmin=-1; 
    if yrange is not None: vmin,vmax=yrange
    plt.imshow( dens [iz0,:,:],        cmap='bwr', vmin=vmin, vmax=vmax, extent=extent, origin="lower" ); plt.colorbar(); plt.title("Charge Density" )
    plt.gca().add_patch( patches.Rectangle( (extent[0],extent[2]),     (extent[1]-extent[0]), (extent[3]-extent[2]), edgecolor='r', facecolor='none' )  )
    plt.axis('equal')

def test_project2D( xs,  g0=[0.0,0.0], dg=[0.1,0.1], ng=[16,16], order=3, ws=None ):
    Ls = [ng[0]*dg[0], ng[1]*dg[1]  ]
    ys = mmff.projectBspline2D( xs, g0, dg, ng, order=order, ws=ws )
    Qtot = ys.sum(); QtotAbs=np.abs(ys).sum();
    #print( "Qtot ",Qtot," order=",order," xs=", xs );
    print( "Qtot ",QtotAbs," order=",order, " g0 ", g0 );

    extent=(g0[0],g0[0]+Ls[0],g0[1],g0[1]+Ls[1])
    plt.figure()
    #plt.plot( xg, ys, '.-', label=("order=%i" %order) )
    plt.imshow( ys, origin="lower", extent=extent, vmin=-1.0,vmax=1.0,  cmap='bwr' )
    plt.colorbar()
    plt.axis('equal')
    plt.gca().add_patch( patches.Rectangle( (extent[0],extent[2]),     (extent[1]-extent[0]), (extent[3]-extent[2]), edgecolor='r', facecolor='none' )  )

    plt.suptitle("g0="+str(g0)+"order="+str(order) )


if __name__ == "__main__":
    # d=0.6
    # apos=np.array([
    #     [-d,.0,0.],
    #     [+d,.0,0.],
    #     [0.,-d,0.],
    #     [0.,+d,0.],
    # ])
    # qs = [ +1.,+1.,-1.,-1. ]



    d=0.6
    apos=np.array([
        [0.,.0,-d],
        [0.,.0,+d],
    ])
    qs = [ +1.,-1. ]

    # d=0.6
    # apos = []
    # qs   = [] 
    # for ix in range(0,10):
    #     for iy in range(0,10):
    #         apos.append( [(ix-5)+0.5, (iy-5)+0.5, +d] ); qs.append(+0.01)
    #         apos.append( [(ix-5)+0.5, (iy-5)+0.5, -d] ); qs.append(-0.01)


    # ----- testing slab dipole correction
    test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.10,0.10], order=3, bSlab=True,  nPBC=[100,100,0] )
    test_vs_direct( apos, qs,  ns=[100,100,150], dg=[0.10,0.10,0.10], order=3, bSlab=True,  nPBC=[100,100,0] )
    test_vs_direct( apos, qs,  ns=[100,100,200], dg=[0.10,0.10,0.10], order=3, bSlab=True,  nPBC=[100,100,0] )
    test_vs_direct( apos, qs,  ns=[100,100,300], dg=[0.10,0.10,0.10], order=3, bSlab=True,  nPBC=[100,100,0] )
    test_vs_direct( apos, qs,  ns=[100,100,400], dg=[0.10,0.10,0.10], order=3, bSlab=True,  nPBC=[100,100,0] )
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.10,0.10], order=3, bSlab=False, nPBC=[200,200,0] )

    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.10,0.10], order=3, bSlab=True,  nPBC=[100,100,0] )
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.10,0.10], order=3, bSlab=False, nPBC=[100,100,0] )

    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.10,0.10], order=3, bPython=True )  
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.10,0.10], order=3, bPython=True, pos0=[0,0,-5.0] )  
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.10,0.10], order=2, bPython=True, pos0=[0,0,-5.0] )  

    #test_vs_direct( apos, qs,  ns=[16,16,16], dg=[0.10,0.10,0.10], order=2, bPython=True, pos0=[-0.8,-0.8,-0.8], bPlot1D=False )
    #test_vs_direct( apos, qs,  ns=[16,16,16], dg=[0.10,0.10,0.10], order=2, bPython=True, pos0=[ 0.0, 0.0,-0.8], bPlot1D=False )  

    #test_project_dens( apos, qs, ns=[16,16,16], pos0=[0.0,0.0,0.0],     order=2 )
    #test_project_dens( apos, qs, ns=[16,16,16], pos0=[-0.8,-0.8,-0.8], order=2 )
    #test_project_dens( apos, qs, ns=[16,16,16], pos0=[ 0.0, 0.0,-0.8], order=2 )

    #test_project2D( apos[:,:2].copy(), g0=[ 0.0, 0.0], order=3, ws=qs )
    #test_project2D( apos[:,:2].copy(), g0=[-0.8,-0.8], order=3, ws=qs )


    # --- change voxel size  homogeneously in all directions
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.15,0.15,0.15] )   # GOOD, This is perfect
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.07,0.07,0.07] )   # GOOD, This is perfect

    # --- change number of grid points homogeneously in all directions
    #test_vs_direct( apos, qs,  ns=[150,150,150], dg=[0.10,0.10,0.10] )    # GOOD, This is perfect
    #test_vs_direct( apos, qs,  ns=[80,80,80],    dg=[0.10,0.10,0.10] )    # GOOD, This is perfect

    # ========= Changing step size in two directions

    # --- change voxel size  homogeneously in xy-directions
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.15,0.15,0.10] )    # GOOD, This is perfect
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.07,0.07,0.10] )    # GOOD, This is perfect

    # --- change voxel size  homogeneously in xz-directions
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.07,0.10,0.07] )    # Quite good
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.15,0.10,0.15] )    # Quite good

    # --- change voxel size  homogeneously in yz-directions
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.07,0.07] )    # Quite good
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.15,0.15] )    # Quite good


    # ========= Changing step size in one direction


    # --- change voxel size  only in x-directions
    #test_vs_direct( apos, qs, ns=[100,100,100], dg=[0.07,0.10,0.10] )    # Quite good
    #test_vs_direct( apos, qs, ns=[100,100,100], dg=[0.15,0.10,0.10] )    # Quite good
    #test_vs_direct( apos, qs, ns=[100,100,100], dg=[0.20,0.10,0.10] )    # NOT SO PERFECT, iax0 goes to ~1.5 almost

    #test_vs_direct( apos, qs, ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[30,30,30] )    # NOT SO PERFECT, iax1 goes to ~1.5 almost
    #test_vs_direct( apos, qs, ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30] )    # NOT SO PERFECT, iax1 goes to ~1.5 almost

    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[30,30,30], order=3, bPython=True )    # NOT SO PERFECT, iax1 goes to ~1.5 almost

    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[100,100,100], order=3, bPython=True ) 

    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[30,30,30], order=3, bPython=True )    # Almost perfect
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=3 )    # Almost perfect
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[30,60,30], order=3 )    # Almost perfect
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[30,30,60], order=3 )    # Almost perfect
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.20,0.10,0.10], nPBC=[30,60,30], order=3 )    # Almost perfect

    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[120,60,60], order=3, bPython=True )    # Almost perfect

    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[30,30,30], order=3 )    # NOT SO PERFECT, iax1 goes to ~1.5 almost
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=3 )    # Almost perfect


    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=3, nBlur=4, cV=0.85  )    # Almost perfect

    #test_vs_direct( apos, qs, ns=[100,100,100], dg=[0.10,0.20,0.10], order=3 )    # NOT SO PERFECT, iax1 goes to ~1.5 almost
    #test_vs_direct( apos, qs, ns=[100,100,100], dg=[0.10,0.10,0.20] )    # GOOD, This is perfect

    # --- change voxel size  only in y-directions
    #test_vs_direct( apos, qs, ns=[100,100,100], dg=[0.10,0.07,0.10] )    # Quite good
    #test_vs_direct( apos, qs, ns=[100,100,100], dg=[0.10,0.15,0.10] )    # Quite good

    # --- change voxel size  only in z-directions
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.10,0.07] )    # GOOD, This is perfect
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.10,0.15] )    # GOOD, This is perfect




    # ========= Changing number of grid points in different directions

    #test_vs_direct( apos, qs,  ns=[150,100,100], dg=[0.10,0.10,0.10] )    # 
    #test_vs_direct( apos, qs,  ns=[100,150,100], dg=[0.10,0.10,0.10] )    # 
    #test_vs_direct( apos, qs,  ns=[100,100,150], dg=[0.10,0.10,0.10] )    # GOOD, This is perfect


    #test_vs_direct( apos, qs,  ns=[200,100,100], dg=[0.10,0.10,0.10] )    # NOT SO PERFECT - iax=0 goes up to 1.2
    #test_vs_direct( apos, qs,  ns=[100,200,100], dg=[0.10,0.10,0.10] )    # NOT SO PERFECT -  iax=1 goes up to 1.2


    #test_vs_direct( apos, qs,  ns=[50,100,100], dg=[0.10,0.10,0.10] )    # NOT SO PERFECT - iax=1 goes up to 1.2
    #test_vs_direct( apos, qs,  ns=[100,50,100], dg=[0.10,0.10,0.10] )    # NOT SO PERFECT -  iax=0 goes up to 1.2


    #test_vs_direct( apos, qs,  ns=[100,100,200], dg=[0.10,0.10,0.10] )   # GOOD, This is perfect
    #test_vs_direct( apos, qs,  ns=[100,100,50], dg=[0.10,0.10,0.10] )    # GOOD, This is perfect



    # ========= Real space smoothening from Bspline(order=3)
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=3, nBlur=2, cV=0.85, yrange=[0.98,1.02] )    # Almost perfect
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=3, nBlur=2, cV=0.90, yrange=[0.98,1.02] )    # Almost perfect
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=3, nBlur=2, cV=0.93, yrange=[0.98,1.02] )    # Almost perfect  # BEST
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=3, nBlur=2, cV=0.95, yrange=[0.98,1.02] )    # Almost perfect
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=3, nBlur=2, cV=0.97, yrange=[0.98,1.02] )    # Almost perfect
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=3, nBlur=3, cV=0.85, yrange=[0.98,1.02] )    # Almost perfect
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=3, nBlur=3, cV=0.90, yrange=[0.98,1.02] )    # Almost perfect
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=3, nBlur=3, cV=0.94, yrange=[0.98,1.02] )    # Almost perfect
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=3, nBlur=3, cV=0.95, yrange=[0.98,1.02] )    # Almost perfect  # BEST
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=3, nBlur=3, cV=0.96, yrange=[0.98,1.02] )    # Almost perfect
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=3, nBlur=4, cV=0.85, yrange=[0.98,1.02] )    # Almost perfect
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=3, nBlur=4, cV=0.90, yrange=[0.98,1.02] )    # Almost perfect
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=3, nBlur=4, cV=0.95, yrange=[0.98,1.02] )    # Almost perfect
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=3, nBlur=4, cV=0.85, yrange=[0.99,1.01] )    # Almost perfect
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=3, nBlur=4, cV=0.90, yrange=[0.99,1.01] )    # Almost perfect
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=3, nBlur=4, cV=0.95, yrange=[0.99,1.01] )    # Almost perfect  # BEST

    # ========= Real space smoothening from Bspline(order=2)

    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=2, nBlur=4, cV=0.850, yrange=[0.99,1.01] )   # Almost perfect
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=2, nBlur=4, cV=0.900, yrange=[0.99,1.01] )   # Almost perfect
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=2, nBlur=4, cV=0.930, yrange=[0.99,1.01] )   # Almost perfect
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=2, nBlur=4, cV=0.935, yrange=[0.99,1.01] )   # Almost perfect  # BEST
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=2, nBlur=4, cV=0.940, yrange=[0.99,1.01] )   # Almost perfect
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=2, nBlur=4, cV=0.950, yrange=[0.99,1.01] )   # Almost perfect
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=2, nBlur=4, cV=0.960, yrange=[0.99,1.01] )   # Almost perfect
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.15,0.10], nPBC=[45,30,30], order=2, nBlur=4, cV=0.930, yrange=[0.99,1.01] )   # Almost perfect  # BEST
    #test_vs_direct( apos, qs,  ns=[100,100,100], dg=[0.10,0.15,0.10], nPBC=[45,30,30], order=2, nBlur=4, cV=0.935, yrange=[0.99,1.01] )   # Almost perfect  # BEST

    plt.show()