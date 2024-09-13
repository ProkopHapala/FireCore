import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")

from pyBall import MMFF as mmff
#from pyBall import atomicUtils as au
#from pyBall import FunctionSampling as fu

COULOMB_CONST  =    14.3996448915 

FFTW_PRESERVE_INPUT = 1 << 4
FFTW_DESTROY_INPUT  = 1 << 0
FFTW_ESTIMATE       = 1 << 6
FFTW_MEASURE        = 0
FFTW_PATIENT        = 1 << 5
FFTW_EXHAUSTIVE     = 1 << 3

def compute_potential(dens, dg, bInternals=False ):
    density_fft = np.fft.fftn(dens)   
    nx, ny, nz = dens.shape           
    # kx = np.fft.fftfreq(nx, d=dg[0]) * 2 * np.pi
    # ky = np.fft.fftfreq(ny, d=dg[1]) * 2 * np.pi
    # kz = np.fft.fftfreq(nz, d=dg[2]) * 2 * np.pi
    kx = np.fft.fftfreq(nx, d=dg[2]) * 2 * np.pi
    ky = np.fft.fftfreq(ny, d=dg[1]) * 2 * np.pi
    kz = np.fft.fftfreq(nz, d=dg[0]) * 2 * np.pi
    Kx, Ky, Kz = np.meshgrid(kx, ky, kz, indexing='ij')    # this seems give better results, (more independent on swapping dx,dy,dz)
    #Kx, Ky, Kz = np.meshgrid(kx, ky, kz)                  # this seems give worse results,  (more dependent on swapping dx,dy,dz)
    # With 'xy' Indexing: The dimensions of the output arrays are ordered as (len(y), len(x), len(z)).
    # With 'ij' Indexing: The dimensions of the output arrays are ordered as (len(x), len(y), len(z)).
    ker_w = 1./(  Kx**2 + Ky**2 + Kz**2 )             
    ker_w[0, 0, 0] = 0.0          
    Vw = density_fft * ker_w      
    V  = np.fft.ifftn(Vw).real        
    if bInternals:
        internals = ( density_fft, Vw, ker_w )
    else:
        internals = None
    return V, internals

def plot_fft_debug( Vs, nx=4, ny=3, iy=0, label="Python/numpy", iz=50 ):
    plt.subplot( ny, nx, iy*nx+1); plt.imshow( Vs[0][iz,:,:]     , cmap='bwr' ); plt.colorbar(); plt.title("kernel(w)   "+label)
    plt.subplot( ny, nx, iy*nx+2); plt.imshow( Vs[1][iz,:,:].real, cmap='bwr' ); plt.colorbar(); plt.title("FFT(density)"+label)
    plt.subplot( ny, nx, iy*nx+3); plt.imshow( Vs[2][iz,:,:].real, cmap='bwr' ); plt.colorbar(); plt.title("V(w)        "+label)
    plt.subplot( ny, nx, iy*nx+4); plt.imshow( Vs[3][iz,:,:]     , cmap='bwr' ); plt.colorbar(); plt.title("V           "+label)

def plot1Dcut( Vg, i0, dg, lvec, iax=0, nPBC=[10,10,10], Vmax=1.0, scErr=100.0 ):

    ix,iy,iz = i0

    if iax==0:
        Vgl = Vg[ iz, iy, : ]
    elif iax==1:
        Vgl = Vg[ iz, :, ix ]
    elif iax==2:
        Vgl = Vg[ :, iy, ix ]

    ns = Vg.shape[::-1]

    nps = ns[iax]   #;print(  "nps ",nps," ns ", ns, " Vg.shape=", Vg.shape, " iax ", iax );
    ps = np.zeros( (nps,3) )
    ps[:,0] = 0.0 #dg[0]*ix
    ps[:,1] = 0.0 #dg[1]*iy
    ps[:,2] = 0.0 #dg[2]*iz
    ps[:,iax] = np.linspace( ns[iax]*dg[iax]*-0.5, ns[iax]*dg[iax]*0.5, nps, endpoint=False )

    fe  = mmff.sampleCoulombPBC(  ps, apos, qs, lvec=lvec, nPBC=nPBC )
    fe0 = mmff.sampleCoulombPBC(  ps, apos, qs, lvec=lvec, nPBC=[0.0,0.0,0.0] )

    vd_vg= fe[:,3]/Vgl

    #mask  = np.abs(fe[:,3]) < 1.0
    factor = np.average( vd_vg[0:nps//3] ); print( "iax=",  iax ,"factor ", factor )


    #plt.figure()
    plt.plot( ps[:,iax], fe0[:,3],   '-', label="direct0" )
    plt.plot( ps[:,iax], fe [:,3],   ':', label="direct" )
    #plt.plot( ps[:,iax], (fe[:,3]-fe0[:,3])*scErr,   '-', label="direct(0-pbc)", lw=0.5 )
    plt.plot( ps[:,iax], Vgl,'-',   label="ewald"  )
    plt.plot( ps[:,iax], vd_vg,'-', label="ref/ewald"  )
    #plt.ylim( -Vmax, Vmax )
    #plt.ylim( -100.0, 100.0 )
    plt.ylim( -2.0, 2.00  )
    #plt.ylim( 0.995, 1.005  )
    plt.legend()
    plt.grid()


def plot2Dcut( Vg, iaxs, i0, dg, lvec,  nPBC=[0,0,0], nplt=3, iplt=0, vmax=1.00 ):
    iax,iay, iaz = iaxs
    ix,iy,iz = i0

    if iaz==0:
        Vgl = Vg[ :, :, ix ]
    elif iaz==1:
        Vgl = Vg[ :, iy, : ]
    elif iaz==2:
        Vgl = Vg[ iz, :, : ]

    ns = Vg.shape[::-1]

    print("ns ", ns, " dg ", dg, "i0=",ix,iy,iz )

    xs = np.linspace( -ns[iax]*dg[iax]*0.5, ns[iax]*dg[iax]*0.5, ns[iax], endpoint=False )
    ys = np.linspace( -ns[iay]*dg[iay]*0.5, ns[iax]*dg[iax]*0.5, ns[iay], endpoint=False )

    #print("xs: \n", xs, " ys: \n", ys )
    Xs,Ys = np.meshgrid( xs, ys )

    nx = len(xs)
    ny = len(ys)

    ps = np.zeros( (nx*ny,3) )
    ps[:,0] = 0.0 #dg[0]*ix
    ps[:,1] = 0.0 #dg[1]*iy
    ps[:,2] = 0.0 #dg[2]*iz
    ps[:,iax] = Xs.flat
    ps[:,iay] = Ys.flat

    #print("ps:   \n", ps )
    
    fe  = mmff.sampleCoulombPBC(  ps, apos, qs, lvec=lvec, nPBC=nPBC ).reshape( (nx,ny,4) )

    print("fe min,max", fe[:,:,3].min(), fe[:,:,3].max() )
    print("Vgl min,max", Vgl.min(), Vgl.max() )

    plt.subplot(2,nplt,iplt+0*nplt+1); plt.imshow( fe[:,:,3],  vmin=-vmax, vmax=vmax,  cmap='bwr', origin='lower' ); plt.colorbar(); plt.title("V_direct" )
    plt.subplot(2,nplt,iplt+1*nplt+1); plt.imshow( Vgl,        vmin=-vmax, vmax=vmax,  cmap='bwr', origin='lower' ); plt.colorbar(); plt.title("V_grid"   )





def test_vs_direct( apos, qs, ns=[100,120,100], dg=[0.1,0.1/1.2,0.1], nPBC=[20,20,20], iax=0, scErr=100.0, order=1, nBlur=0, cSOR=0.0, cV=1.0, bPython=False ):
    apos = apos.copy()
    print( "apos ", apos )
    lvec=[
        [ ns[0]*dg[0], 0.0, 0.0 ],
        [ 0.0, ns[1]*dg[1], 0.0 ],
        [ 0.0, 0.0, ns[2]*dg[2] ],
    ]
    Ls = [ lvec[0][0], lvec[1][1], lvec[2][2] ]

    #order=1
    order=2
    #order=3  # Quintic spline works badly at the moment

    pos0=np.array( [ -0.5*dg[0]*ns[0],  -0.5*dg[1]*ns[1],  -0.5*dg[2]*ns[2]  ] )
    mmff.setupEwaldGrid( ns, dg=dg, pos0=pos0 )
    dens = mmff.projectAtomsEwaldGrid( apos, qs, ns=ns, order=order )

    Qtot = np.abs(dens).sum(); print("Qtot = ", Qtot)
    
    if bPython:
        #Vg,_ = compute_potential(dens, dg)
        Vg,(density_fft, Vw, ker_w) = compute_potential(dens, dg, bInternals=True )

        sh = Vg.shape
        nx, ny, nz = sh
        print( "test_vs_direct()   bPython = True ", sh, ns, Ls, dg )
        #scEwald = COULOMB_CONST * (nx*ny*ny)**(1./3.)   #* (nx*nx + ny*ny + nz*nz)/3.0
        #scEwald = COULOMB_CONST *  (nx*nx + ny*ny + nz*nz)/3.0
        #scEwald = COULOMB_CONST *  (nx*nx + ny*ny + nz*nz)/3.0    * 1.24

        #V = Ls[0]*Ls[1]*Ls[2]

        #scEwald = COULOMB_CONST *  (nx*nx + ny*ny + nz*nz)/3.0  * 0.37343958303  # ns=[100,100,100], dg=[0.15,0.15,0.15],
        #scEwald = COULOMB_CONST *  (nx*nx + ny*ny + nz*nz)/3.0  * 1.2577692859   # ns=[100,100,100], dg=[0.10,0.10,0.10],
        #scEwald = COULOMB_CONST *  (nx*nx + ny*ny + nz*nz)/3.0  * 3.66329865007  # ns=[100,100,100], dg=[0.07,0.07,0.07],

        #scEwald = COULOMB_CONST *  (nx*nx + ny*ny + nz*nz)/3.0  * 0.87332844942  # ns=[120,120,120], dg=[0.1,0.1,0.1], 
        #scEwald = COULOMB_CONST *  (nx*nx + ny*ny + nz*nz)/3.0  * 1.964433827227 # ns=[80,80,80],    dg=[0.1,0.1,0.1], 

        #scEwald = COULOMB_CONST *  (nx*nx + ny*ny + nz*nz)/3.0  * 0.87332844942  # ns=[120,120,120], dg=[0.1,0.1,0.1], 
        #scEwald = COULOMB_CONST *  (nx*nx + ny*ny + nz*nz)/3.0   * 1.96443382722  # ns=[80,80,80],    dg=[0.1,0.1,0.1], 


        #scEwald = COULOMB_CONST *  (nx*nx + ny*ny + nz*nz)/3.0  * 1.2577692859* 1.198912796701991   # ns=[120,120,120], dg=[0.1/1.2,0.1/1.2,0.1/1.2],

        # 12577.692859040353  / (2*pi)^3  = 50.7062367917
        #   1.2*1.2*1.2 = 1.728     1/1.728 = 0.5787037037
        #   1.5*1.5*1.5 = 3.375     1/3.375 = 0.2961799984

        # Observation : if we change now dx and nx proportionally, so that the size Lx=dx*nx does not change, then the scEwald is the same


        dV = dg[0]*dg[1]*dg[2]  # voxel volume
 
        scEwald = COULOMB_CONST * 4*np.pi / dV      # the normmalization factor 1/(dx*dy*dz) is there because our density is not really density, it is the charge, so we need to divide it by voxel volume  dV
        #scEwald = COULOMB_CONST * 4*np.pi / dV     # ns=[100,100,100], dg=[0.1,0.1,0.1],
        #scEwald = COULOMB_CONST * 12577.692859040353       # ns=[120,120,120], dg=[0.1,0.1,0.1],
        #scEwald = COULOMB_CONST * 12577.692859040353       # ns=[150,150,150], dg=[0.1,0.1,0.1],
        # => The constant scEwald is independent of the number of grid points (nx,ny,nz), but it depends on the grid spacing (dg)

        #scEwald = COULOMB_CONST * 12577.692859040353 *1.728       # ns=[120,120,120], dg=[0.1/1.2,0.1/1.2,0.1/1.2],
    
        #scEwald = COULOMB_CONST * 12577.692859040353 *3.375       # ns=[150,150,150], dg=[0.1/1.5,0.1/1.5,0.1/1.5],


        #scEwald = COULOMB_CONST   * 21714.561967164984         # ns=[120,120,120], dg=[0.1/1.2,0.1/1.2,0.1/1.2],

        #scEwald = COULOMB_CONST *  nx*ny*nz
        #scEwald = COULOMB_CONST *  (nx*nx + ny*ny + nz*nz)/3.0    * 1.24 * 1.0
        #scEwald = COULOMB_CONST #*  V #* 10
        #scEwald = COULOMB_CONST 


    else:
        print( "test_vs_direct()   bPython = False " )
        #Vg, density_fft, Vw, ker_w = compute_potential(dens, dg)
        Vg = mmff.EwaldGridSolveLaplace( dens, nBlur=nBlur, cSOR=cSOR, cV=cV )

        #scEwald = 0.2
        # bulhar1 = 0.9973498924
        # bulhar2 = 0.99**(1./3.)
        # c213 = 2.0**(1.0/3.0)  * bulhar1
        # scEwald = COULOMB_CONST * np.sqrt(2.0)/100.0  ;print("scEwald = ", scEwald)
        #scEwald = COULOMB_CONST*c213/100.0  ;print("scEwald = ", scEwald)
        #scEwald = COULOMB_CONST*np.sqrt(2.0) *(dg[0]*dg[1]*dg[2]) /( ns[0]*ns[1]*ns[2])  ;print("scEwald = ", scEwald)
        #scEwald = COULOMB_CONST*np.sqrt(2.0) /( (dg[0]*dg[1]*dg[2]) * ns[0]*ns[1]*ns[2])  ;print("scEwald = ", scEwald)
        scEwald = COULOMB_CONST*np.sqrt(2.0)/np.power( ns[0]*ns[1]*ns[2], 1./3.)  ;print("scEwald = ", scEwald)
        #scEwald = COULOMB_CONST*np.sqrt(2.0) / ( ( ns[0]*ns[1]*ns[2] ) * (dg[0]*dg[1]*dg[2]) ) ;print("scEwald = ", scEwald)
        #scEwald = COULOMB_CONST*2 *  (dg[0]*dg[1]*dg[2]) *   1.102/np.sqrt( ns[0]**2 + ns[1]**2 + ns[2]**2)   ;print("scEwald = ", scEwald)
        #scEwald1 = np.sqrt(2.0)/np.power( ns[0]*ns[1]*ns[2], 1./3.)  
        #scEwald2 = 2.0         /np.sqrt( ns[0]**2 + ns[1]**2 + ns[2]**2)  
        #scEwald = COULOMB_CONST*np.sqrt( scEwald1*scEwald2 )
        #scEwald1 = np.sqrt(2.0)/np.power( ns[0]*ns[1]*ns[2], 1./3.)  
        #scEwald2 = 2.0         /np.sqrt( ns[0]**2 + ns[1]**2 + ns[2]**2)  
        #scEwald = COULOMB_CONST * (1./3.) * np.sqrt(0.5)/np.sqrt(3.) * ( ns[0]**2 + ns[1]**2 + ns[2]**2 )/( ns[0]*ns[1]*ns[2] ) 
        #scEwald = COULOMB_CONST * (1./3.) * np.sqrt(0.5)/np.sqrt(3.) * ( ns[0]**2 + ns[1]**2 + ns[2]**2 )/( ns[0]*ns[1]*ns[2] ) 

        #scEwald = COULOMB_CONST*( ns[0]**2 + ns[1]**2 + ns[2]**2 )/(ns[0]*ns[1]*ns[2])  ;print("scEwald = ", scEwald)
        #scEwald = COULOMB_CONST*np.power( ns[0]*ns[1]*ns[2], 1./3.)/( ns[0]**2 + ns[1]**2 + ns[2]**2)  ;print("scEwald = ", scEwald)
        #scEwald = COULOMB_CONST/100.0;    print("scEwald = ", scEwald)

    print("scEwald = ", scEwald)
    Vg *= scEwald
    ix0=ns[0]//2
    iy0=ns[1]//2
    iz0=ns[2]//2
    i0 = (ix0,iy0,iz0)


    # ---- Plot 2D
    # dens2d =  np.sum( dens, axis=iax )
    plt.figure(figsize=(15,5))
    plt.subplot(1,3,1); plt.imshow( dens [iz0,:,:],                        cmap='bwr' ); plt.colorbar(); plt.title("Charge Density" )
    plt.subplot(1,3,2); plt.imshow( ker_w[iz0,:,:],                        cmap='bwr' ); plt.colorbar(); plt.title("ker_w" )
    plt.subplot(1,3,3); plt.imshow( Vg   [iz0,:,:],   vmin=-1.0,vmax=1.0,  cmap='bwr' ); plt.colorbar(); plt.title("V ewald C++/FFTW3 " )
    
 
    # # plt.figure(figsize=(15,10))
    # plot2Dcut( Vg, (0,1,2), i0, dg, lvec,  nPBC=[0,0,0], nplt=3, iplt=0 )
    # plot2Dcut( Vg, (0,2,1), i0, dg, lvec,  nPBC=[0,0,0], nplt=3, iplt=1 )
    # plot2Dcut( Vg, (1,2,0), i0, dg, lvec,  nPBC=[0,0,0], nplt=3, iplt=2 )


    plt.figure(figsize=(10,5))
    plt.subplot(1,2,1); plot1Dcut( Vg, i0=i0, dg=dg, lvec=lvec, iax=0, nPBC=nPBC, scErr=scErr )
    plt.subplot(1,2,2); plot1Dcut( Vg, i0=i0, dg=dg, lvec=lvec, iax=1, nPBC=nPBC, scErr=scErr )
    #plt.subplot(1,3,3); plot1Dcut( Vg, i0=i0, dg=dg, lvec=lvec, iax=2, nPBC=nPBC, scErr=scErr )
    
    plt.tight_layout()

    #plt.show()
    



def test_poison( apos, qs, bPlot=True, bDebug=True, iz=50, ns=[90,100,110], dg=[0.1,0.1,0.1], flags=-1, bOMP=False ):

    mmff.setupEwaldGrid( ns, dg=dg )
    dens = mmff.projectAtomsEwaldGrid( apos, qs, ns=ns )

    V, (density_fft, Vw, ker_w) = compute_potential(dens, dg, bInternals=True)

    #print("DEBUG 1 ")
    if bDebug:
        c_V, c_densw, c_kerw, c_VwKer = mmff.EwaldGridSolveLaplaceDebug( dens )
    else:
        c_V = mmff.EwaldGridSolveLaplace( dens, bPrepare=True, bDestroy=True, flags=flags, bOMP=bOMP )
    #print("DEBUG 2 ")
    Vmax= 0.0001
    Qmax= dens.max()
    #print( "Vmax, Qmax ", Vmax, Qmax )
    #print("DEBUG 3 ")
    if bPlot:
        if bDebug:
            #print("DEBUG 4a ")
            plt.figure(figsize=(20,15))
            plt.subplot(3,4,1); plt.imshow( dens[iz,:,:], vmin=-Qmax, vmax=Qmax, cmap='bwr' ); plt.colorbar(); plt.title("Charge Density")
            #plt.subplot(1,3,2); plt.imshow( V   [50,:,:], vmin=-Vmax, vmax=Vmax, cmap='bwr' ); plt.colorbar(); plt.title("Electronstatic Potential (python)")
            plot_fft_debug( [ker_w,  density_fft, Vw , V  ], iy=1, label="Python/numpy" )
            plot_fft_debug( [c_kerw, c_densw, c_VwKer, c_V], iy=2, label="C++/FFTW3   " )
        else:
            #print("DEBUG 4b ")
            plt.figure(figsize=(15,5))
            plt.subplot(1,3,1); plt.imshow( dens[iz,:,:], vmin=-Qmax, vmax=Qmax, cmap='bwr' ); plt.colorbar(); plt.title("Charge Density" )
            #plt.subplot(1,3,2); plt.imshow( V   [iz,:,:], vmin=-Vmax, vmax=Vmax, cmap='bwr' ); plt.colorbar(); plt.title("V Python/numpy" )
            #plt.subplot(1,3,3); plt.imshow( c_V [iz,:,:], vmin=-Vmax, vmax=Vmax, cmap='bwr' ); plt.colorbar(); plt.title("V C++/FFTW3 "   )
            plt.subplot(1,3,2); plt.imshow( V   [iz,:,:], cmap='bwr' ); plt.colorbar(); plt.title("V Python/numpy" )
            plt.subplot(1,3,3); plt.imshow( c_V [iz,:,:], cmap='bwr' ); plt.colorbar(); plt.title("V C++/FFTW3 "   )
        #print("DEBUG 4 ")
        plt.show()

# =========== Main

d=0.6
apos=np.array([
    #[5.0-d,5.0-d,5.0],
    #[5.0+d,5.0-d,5.0],
    #[5.0-d,5.0+d,5.0],
    #[5.0+d,5.0+d,5.0],
    [-d,.0,0.],
    [+d,.0,0.],
    [0.,-d,0.],
    [0.,+d,0.],
])
#qs = [ -1.,+1.,+1.,-1. ]
qs = [ +1.,+1.,-1.,-1. ]



ns=[200,200,200]

#test_poison(  apos, qs, bDebug=True )
#test_poison(  apos, qs, bDebug=False, flags=FFTW_PRESERVE_INPUT | FFTW_ESTIMATE  )
#test_poison(  apos, qs, bDebug=False, flags=FFTW_PRESERVE_INPUT | FFTW_MEASURE   )
#test_poison(  apos, qs, bDebug=False, flags=FFTW_PRESERVE_INPUT | FFTW_PATIENT   )

#test_poison(  apos, qs, bDebug=False, flags=FFTW_DESTROY_INPUT | FFTW_ESTIMATE  )
#test_poison(  apos, qs, bDebug=False, flags=FFTW_DESTROY_INPUT | FFTW_MEASURE   )
#test_poison(  apos, qs, bDebug=False, flags=FFTW_DESTROY_INPUT | FFTW_ESTIMATE, ns=ns, bOMP=False )
#test_poison(  apos, qs, bDebug=False, flags=FFTW_DESTROY_INPUT | FFTW_PATIENT, ns=ns, bOMP=False )
#test_poison(  apos, qs, bDebug=False, flags=FFTW_DESTROY_INPUT | FFTW_PATIENT, ns=ns, bOMP=True  )



#test_vs_direct( apos, qs, order=1, nBlur=8 )
#test_vs_direct( apos, qs, order=2, nBlur=2 )
#test_vs_direct( apos, qs, order=2, nBlur=4 )
#test_vs_direct( apos, qs, order=2, nBlur=8, cSOR=0.0  )
#test_vs_direct( apos, qs, order=2, nBlur=8, cSOR=-0.2 )
#test_vs_direct( apos, qs, order=2, nBlur=8, cSOR=+0.2 )

#test_vs_direct( apos, qs, order=2, nBlur=4, cSOR=0.0  )
#test_vs_direct( apos, qs, order=2, nBlur=4, cSOR=-0.1 )
#test_vs_direct( apos, qs, order=2, nBlur=4, cSOR=-0.2 )
#test_vs_direct( apos, qs, order=2, nBlur=4, cSOR=-0.3 )
#test_vs_direct( apos, qs, order=2, nBlur=4, cSOR=-0.4 )
#test_vs_direct( apos, qs, order=2, nBlur=4, cSOR=-0.5 )

'''
test_vs_direct( apos, qs, order=2, nBlur=4, cV=1.00 )
test_vs_direct( apos, qs, order=2, nBlur=4, cV=0.95 )
test_vs_direct( apos, qs, order=2, nBlur=4, cV=0.90 )
test_vs_direct( apos, qs, order=2, nBlur=4, cV=0.85 ) # The best
test_vs_direct( apos, qs, order=2, nBlur=4, cV=0.80 )
test_vs_direct( apos, qs, order=2, nBlur=4, cV=0.75 )
test_vs_direct( apos, qs, order=2, nBlur=4, cV=0.70 )
test_vs_direct( apos, qs, order=2, nBlur=4, cV=0.65 )
'''

test_vs_direct( apos, qs, order=3, nBlur=0, iax=0, bPython=True )

#test_vs_direct( apos, qs, order=3, nBlur=0, iax=0 )
#test_vs_direct( apos, qs, order=3, nBlur=0, iax=1 )
#test_vs_direct( apos, qs, order=3, nBlur=0, iax=2 )

#test_vs_direct( apos, qs, order=3, nBlur=4, cV=0.85 )
#test_vs_direct( apos, qs, order=3, nBlur=4, cV=0.90 )
#test_vs_direct( apos, qs, order=3, nBlur=4, cV=0.95 )


#test_vs_direct( apos, qs, order=3, nBlur=2, cV=0.94, cSOR=+0.10 )
#test_vs_direct( apos, qs, order=3, nBlur=2, cV=0.95, cSOR=+0.05 )
#test_vs_direct( apos, qs, order=3, nBlur=2, cV=0.95, cSOR=+0.15 )
#test_vs_direct( apos, qs, order=3, nBlur=2, cV=0.94, cSOR=-0.10 )
#test_vs_direct( apos, qs, order=3, nBlur=2, cV=0.94, cSOR=-0.20 )

#test_vs_direct( apos, qs, order=3, nBlur=2, cV=1.00 )
#test_vs_direct( apos, qs, order=3, nBlur=2, cV=0.97 )
#test_vs_direct( apos, qs, order=3, nBlur=2, cV=0.95 )
#test_vs_direct( apos, qs, order=3, nBlur=2, cV=0.94 )
#test_vs_direct( apos, qs, order=3, nBlur=2, cV=0.93 )
#test_vs_direct( apos, qs, order=3, nBlur=2, cV=0.90 )

# test_vs_direct( apos, qs, order=3, nBlur=2, cV=0.90 )
# test_vs_direct( apos, qs, order=3, nBlur=2, cV=0.85 )
# test_vs_direct( apos, qs, order=3, nBlur=2, cV=0.80 )
# test_vs_direct( apos, qs, order=3, nBlur=2, cSOR=-0.10, cV=-2. )
# test_vs_direct( apos, qs, order=3, nBlur=2, cSOR=-0.15, cV=-2. )
# test_vs_direct( apos, qs, order=3, nBlur=2, cSOR=-0.20, cV=-2. )
#test_vs_direct( apos, qs, order=3, nBlur=0 )

#test_vs_direct( apos, qs, order=2, nBlur=4, cV=0.6 )
#test_vs_direct( apos, qs, order=2, nBlur=4, cV=0.4 )
#test_vs_direct( apos, qs, order=2, nBlur=4, cV=0.2 )
#test_vs_direct( apos, qs, order=2, nBlur=0 )
#test_vs_direct( apos, qs, order=3 )
plt.show()

performance_results='''

using Build-opt -Ofast

prepare_laplace() flags=80 n(100,100,100) T(prepare_laplace)= 163.56  [Mticks] T(solve_laplace)= 120.55  [Mticks]      FFTW_PRESERVE_INPUT | FFTW_ESTIMATE
prepare_laplace() flags=16 n(100,100,100) T(prepare_laplace)= 6.46114 [Mticks] T(solve_laplace)= 117.883 [Mticks]      FFTW_PRESERVE_INPUT | FFTW_MEASURE 
prepare_laplace() flags=48 n(100,100,100) T(prepare_laplace)= 6.57928 [Mticks] T(solve_laplace)= 116.977 [Mticks]      FFTW_PRESERVE_INPUT | FFTW_PATIENT
prepare_laplace() flags=65 n(100,100,100) T(prepare_laplace)= 6.64837 [Mticks] T(solve_laplace)= 117.1   [Mticks]      FFTW_DESTROY_INPUT  | FFTW_ESTIMATE
prepare_laplace() flags=1  n(100,100,100) T(prepare_laplace)= 6.57567 [Mticks] T(solve_laplace)= 116.077 [Mticks]      FFTW_DESTROY_INPUT  | FFTW_MEASURE
prepare_laplace() flags=33 n(100,100,100) T(prepare_laplace)= 6.45647 [Mticks] T(solve_laplace)= 116.67  [Mticks]      FFTW_DESTROY_INPUT  | FFTW_PATIENT


prepare_laplace() flags=80 n(100,100,100) T(prepare_laplace)= 17.8152 [Mticks] T(solve_laplace)= 108.935 [Mticks]     FFTW_PRESERVE_INPUT | FFTW_ESTIMATE
prepare_laplace() flags=16 n(100,100,100) T(prepare_laplace)= 16.4503 [Mticks] T(solve_laplace)= 109.146 [Mticks]     FFTW_PRESERVE_INPUT | FFTW_MEASURE 
prepare_laplace() flags=48 n(100,100,100) T(prepare_laplace)= 18.5292 [Mticks] T(solve_laplace)= 107.973 [Mticks]     FFTW_PRESERVE_INPUT | FFTW_PATIENT
prepare_laplace() flags=1 n(100,100,100)  T(prepare_laplace)= 16.908  [Mticks] T(solve_laplace)= 107.678 [Mticks]     FFTW_DESTROY_INPUT  | FFTW_ESTIMATE
prepare_laplace() flags=1 n(100,100,100)  T(prepare_laplace)= 17.912  [Mticks] T(solve_laplace)= 109.982 [Mticks]     FFTW_DESTROY_INPUT  | FFTW_MEASURE
prepare_laplace() flags=33 n(100,100,100) T(prepare_laplace)= 17.9685 [Mticks] T(solve_laplace)= 107.581 [Mticks]     FFTW_DESTROY_INPUT  | FFTW_PATIENT


with OpenMP
prepare_laplace() flags=80 n(100,100,100) T(prepare_laplace)= 17.9227 [Mticks] T(solve_laplace)= 49.495 [Mticks] 

'''