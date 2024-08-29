import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")

from pyBall import MMFF as mmff
#from pyBall import atomicUtils as au
#from pyBall import FunctionSampling as fu


def compute_potential(dens, dg):
    density_fft = np.fft.fftn(dens)   # Perform the forward 3D FFT of the density
    nx, ny, nz = dens.shape           # Get the dimensions of the density grid
    
    # Create the k-space grid (frequency domain)
    kx = np.fft.fftfreq(nx, d=dg[0]) * 2 * np.pi
    ky = np.fft.fftfreq(ny, d=dg[1]) * 2 * np.pi
    kz = np.fft.fftfreq(nz, d=dg[2]) * 2 * np.pi

    Kx, Ky, Kz = np.meshgrid(kx, ky, kz, indexing='ij')  # Create 3D grids for kx, ky, kz
    ker_w = 1./( Kx**2 + Ky**2 + Kz**2 )                  # Calculate the squared magnitude of the wave vectors
    ker_w[0, 0, 0] = 0.0             # Avoid division by zero at the zero frequency component (k=0)
    Vw = density_fft * ker_w         # Divide the FFT of the density by k^2 to get the potential in Fourier space
    V  = np.fft.ifftn(Vw).real        # Perform the inverse 3D FFT to get the potential in real space    
    return V

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

ns = [100,100,100]
dg = [0.1,0.1,0.1]
mmff.setupEwaldGrid( ns, dg=dg )
dens = mmff.projectAtomsEwaldGrid( apos, qs, ns=ns )
V = compute_potential(dens, dg)

Vmax= 0.0001
Qmax= dens.max()
print( "Vmax, Qmax ", Vmax, Qmax )

plt.figure(figsize=(10,5))
plt.subplot(1,2,1); plt.imshow( dens[50,:,:], vmin=-Qmax, vmax=Qmax, cmap='bwr' ); plt.colorbar(); plt.title("Charge Density")
plt.subplot(1,2,2); plt.imshow( V   [50,:,:], vmin=-Vmax, vmax=Vmax, cmap='bwr' ); plt.colorbar(); plt.title("Electronstatic Potential")
plt.show()