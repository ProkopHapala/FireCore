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


def test_poison( apos, qs, bPlot=True, bDebug=True, iz=50, ns=[100,100,100], dg=[0.1,0.1,0.1], flags=-1, bOMP=False ):

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
apos=[
    #[5.0-d,5.0-d,5.0],
    #[5.0+d,5.0-d,5.0],
    #[5.0-d,5.0+d,5.0],
    #[5.0+d,5.0+d,5.0],

    [5.0-d,5.0  ,5.0],
    [5.0+d,5.0  ,5.0],
    [5.0  ,5.0-d,5.0],
    [5.0  ,5.0+d,5.0],
]
#qs = [ -1.,+1.,+1.,-1. ]
qs = [ +1.,+1.,-1.,-1. ]



ns=[200,200,200]

#test_poison(  apos, qs, bDebug=True )
#test_poison(  apos, qs, bDebug=False, flags=FFTW_PRESERVE_INPUT | FFTW_ESTIMATE  )
#test_poison(  apos, qs, bDebug=False, flags=FFTW_PRESERVE_INPUT | FFTW_MEASURE   )
#test_poison(  apos, qs, bDebug=False, flags=FFTW_PRESERVE_INPUT | FFTW_PATIENT   )

#test_poison(  apos, qs, bDebug=False, flags=FFTW_DESTROY_INPUT | FFTW_ESTIMATE  )
#test_poison(  apos, qs, bDebug=False, flags=FFTW_DESTROY_INPUT | FFTW_MEASURE   )
test_poison(  apos, qs, bDebug=False, flags=FFTW_DESTROY_INPUT | FFTW_PATIENT, ns=ns, bOMP=False )
test_poison(  apos, qs, bDebug=False, flags=FFTW_DESTROY_INPUT | FFTW_PATIENT, ns=ns, bOMP=True  )

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