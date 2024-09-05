import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")

from pyBall import MMFF as mmff
#from pyBall import atomicUtils as au
#from pyBall import FunctionSampling as fu


FFTW_PRESERVE_INPUT = 1 << 4
FFTW_DESTROY_INPUT  = 1 << 0
FFTW_ESTIMATE       = 1 << 6
FFTW_MEASURE        = 0
FFTW_PATIENT        = 1 << 5
FFTW_EXHAUSTIVE     = 1 << 3

def compute_potential(dens, dg):
    density_fft = np.fft.fftn(dens)   # Perform the forward 3D FFT of the density
    nx, ny, nz = dens.shape           # Get the dimensions of the density grid
    
    # Create the k-space grid (frequency domain)
    kx = np.fft.fftfreq(nx, d=dg[0]) * 2 * np.pi
    ky = np.fft.fftfreq(ny, d=dg[1]) * 2 * np.pi
    kz = np.fft.fftfreq(nz, d=dg[2]) * 2 * np.pi

    Kx, Ky, Kz = np.meshgrid(kx, ky, kz, indexing='ij')  # Create 3D grids for kx, ky, kz
    ker_w = 1./(  Kx**2 + Ky**2 + Kz**2 ) 
    #K2 = Kx**2 + Ky**2 + Kz**2
    #mask = K2 > 0
    #ker_w[mask] = 1./( K2[mask] )                 # Calculate the squared magnitude of the wave vectors
    ker_w[0, 0, 0] = 0.0             # Avoid division by zero at the zero frequency component (k=0)
    Vw = density_fft * ker_w         # Divide the FFT of the density by k^2 to get the potential in Fourier space
    V  = np.fft.ifftn(Vw).real        # Perform the inverse 3D FFT to get the potential in real space    
    return V, density_fft, Vw, ker_w

def plot_fft_debug( Vs, nx=4, ny=3, iy=0, label="Python/numpy", iz=50 ):
    plt.subplot( ny, nx, iy*nx+1); plt.imshow( Vs[0][iz,:,:]     , cmap='bwr' ); plt.colorbar(); plt.title("kernel(w)   "+label)
    plt.subplot( ny, nx, iy*nx+2); plt.imshow( Vs[1][iz,:,:].real, cmap='bwr' ); plt.colorbar(); plt.title("FFT(density)"+label)
    plt.subplot( ny, nx, iy*nx+3); plt.imshow( Vs[2][iz,:,:].real, cmap='bwr' ); plt.colorbar(); plt.title("V(w)        "+label)
    plt.subplot( ny, nx, iy*nx+4); plt.imshow( Vs[3][iz,:,:]     , cmap='bwr' ); plt.colorbar(); plt.title("V           "+label)


def test_vs_direct( apos, qs, ns=[90,100,110], dg=[0.1,0.1,0.1], nPBC=[50,50,50], iax=0, scErr=100.0, order=2, nBlur=0, cSOR=0.0, cV=1.0 ):
    apos = apos.copy()
    apos[:,0] += dg[0]*ns[0]*0.5
    apos[:,1] += dg[1]*ns[1]*0.5
    apos[:,2] += dg[2]*ns[2]*0.5
    print( "apos ", apos )

    mmff.setupEwaldGrid( ns, dg=dg )
    dens = mmff.projectAtomsEwaldGrid( apos, qs, ns=ns, order=order )
    #Vg, density_fft, Vw, ker_w = compute_potential(dens, dg)
    Vg = mmff.EwaldGridSolveLaplace( dens, nBlur=nBlur, cSOR=cSOR, cV=cV )
    lvec=[
        [ ns[0]*dg[1], 0.0, 0.0 ],
        [ 0.0, ns[1]*dg[1], 0.0 ],
        [ 0.0, 0.0, ns[2]*dg[2] ],
    ]

    #scEwald = 0.2
    COULOMB_CONST  =    14.3996448915 

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

    
    # ---- Plot 2D
    # dens2d =  np.sum( dens, axis=iax )
    # plt.figure(figsize=(15,5))
    # plt.subplot(1,3,1); plt.imshow( dens2d,                      cmap='bwr' ); plt.colorbar(); plt.title("Charge Density" )
    # #plt.subplot(1,3,1); plt.imshow( dens[iz,:,:],                      cmap='bwr' ); plt.colorbar(); plt.title("Charge Density" )
    # plt.subplot(1,3,2); plt.imshow( Vg[iz,:,:],   vmin=-1.0,vmax=1.0,  cmap='bwr' ); plt.colorbar(); plt.title("V ewald C++/FFTW3 " )
    
    # xs = np.linspace( 0.0, ns[0]*dg[0], ns[0], endpoint=False )
    # ys = np.linspace( 0.0, ns[1]*dg[1], ns[1], endpoint=False )
    # Xs,Ys = np.meshgrid( xs, ys )
    # ps = np.zeros( (ns[0],ns[1],3) )
    # ps[:,:,0] = Xs
    # ps[:,:,1] = Ys
    # ps[:,:,2] = dg[2]*iz
    # ps = np.reshape( ps, (ns[0]*ns[1],3) )
    # fe = mmff.sampleCoulombPBC(  ps, apos, qs, lvec=lvec, nPBC=nPBC ).reshape( (ns[0],ns[1],4) )
    # plt.subplot(1,3,3); plt.imshow( fe[:,:,3],    vmin=-1.0,vmax=1.0,  cmap='bwr' ); plt.colorbar(); plt.title("V direct C++      " )
    plt.show()
    

    
    
    ix=ns[0]//2
    iy=ns[1]//2
    iz=ns[2]//2

    if iax==0:
        Vgl = Vg[ iz, iy, : ]
    elif iax==1:
        Vgl = Vg[ iz, :, ix ]
    elif iax==2:
        Vgl = Vg[ :, iy, ix ]

    nps = ns[iax]   ;print(  "nps ",nps, " Vg.shape=", Vg.shape, " iax ", iax );
    ps = np.zeros( (nps,3) )
    ps[:,0] = dg[0]*ix
    ps[:,1] = dg[1]*iy
    ps[:,2] = dg[2]*iz
    ps[:,iax] = np.linspace( 0.0, ns[iax]*dg[iax], nps, endpoint=False )

    fe  = mmff.sampleCoulombPBC(  ps, apos, qs, lvec=lvec, nPBC=nPBC )
    fe0 = mmff.sampleCoulombPBC(  ps, apos, qs, lvec=lvec, nPBC=[0.0,0.0,0.0] )

    # plt.figure(figsize=(10,5))
    # plt.subplot(1,2,1); plt.imshow( dens[iz,:,:], cmap='bwr' ); plt.colorbar(); plt.title("Charge Density" )
    # plt.subplot(1,2,2); plt.imshow( Vg[iz,:,:],   cmap='bwr' ); plt.colorbar(); plt.title("V C++/FFTW3 "   )

    #Vgsc = 0.001
    #Vg*=scEwald
    #Vmin = np.min(Vg); Vmax= np.max(Vg); Vmax= np.maximum( Vmax, -Vmin )
    Vmax = 1.2
    plt.figure()
    plt.plot( ps[:,iax], fe0[:,3],   '-', label="direct0" )
    plt.plot( ps[:,iax], fe [:,3],   ':', label="direct" )
    plt.plot( ps[:,iax], (fe[:,3]-fe0[:,3])*scErr,   '-', label="direct(0-pbc)", lw=0.5 )
    plt.plot( ps[:,iax], Vgl,'-', label="ewald"  )

    plt.plot( ps[:,iax], fe [:,3]/Vgl,'-', label="ref/ewald"  )



    #plt.ylim( -Vmax, Vmax )
    #plt.ylim( -100.0, 100.0 )
    plt.ylim( -1.0, 2.00  )
    #plt.ylim( 0.995, 1.005  )
    plt.legend()
    plt.grid()

    name = f"order_{order}_nBLur_{nBlur}_cV={cV:.2}_cSOR={cSOR:.2}"
    plt.title( name );

    plt.savefig( name+".png", bbox_inches='tight')
    

    #plt.show()
    



def test_poison( apos, qs, bPlot=True, bDebug=True, iz=50, ns=[90,100,110], dg=[0.1,0.1,0.1], flags=-1, bOMP=False ):

    mmff.setupEwaldGrid( ns, dg=dg )
    dens = mmff.projectAtomsEwaldGrid( apos, qs, ns=ns )

    V, density_fft, Vw, ker_w = compute_potential(dens, dg)

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

d=0.4
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

test_vs_direct( apos, qs, order=3, nBlur=0, iax=0 )
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