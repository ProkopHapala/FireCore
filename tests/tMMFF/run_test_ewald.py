import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import MMFF as mmff

COULOMB_CONST  =    14.3996448915 

def compute_potential(dens, dg, bInternals=False):
    """
    Compute the electrostatic potential using Particle Mesh Ewald (PME) method.

    Parameters:
    - dens: 3D numpy array with shape [nz, ny, nx], representing charge density.
    - dg: List or array of grid spacings [dz, dy, dx].
    - bInternals: Boolean flag to return internal variables.

    Returns:
    - V: 3D numpy array of electrostatic potential.
    - internals: Tuple of internal variables if bInternals is True.
    """
    dx, dy, dz = dg                  # Unpack grid spacings (NOTE: grid is in [z, y, x] order) but dg=(dx,dy,dz)
    nz, ny, nx = dens.shape          # we assume 3D arrays are in the order [z, y, x]
    dV         = dx * dy * dz        # voxel volume
    rho        = dens / dV           # from electron count to charge density
    density_fft = np.fft.fftn(rho)   # Perform FFT on the density
    
    # calculate inverse Laplace operator in k-space
    kz = np.fft.fftfreq(nz, d=dz) * 2 * np.pi
    ky = np.fft.fftfreq(ny, d=dy) * 2 * np.pi
    kx = np.fft.fftfreq(nx, d=dx) * 2 * np.pi
    Kz, Ky, Kx = np.meshgrid(kz, ky, kx, indexing='ij')  # Create 3D k-space grids with correct axis alignment
    ker_w = 4*np.pi / (Kx**2 + Ky**2 + Kz**2)            # Compute the Coulomb kernel in k-space
    ker_w[0, 0, 0] = 0.0                                 # Handle the singularity at k=0
    
    Vw = density_fft * ker_w                             # solve Poisson's equation (rho = Laplace * V) in k-space (multiply fourier image of density by inverse Laplace operator)
    V  = np.fft.ifftn(Vw).real                           # Inverse FFT to obtain the potential in real space
    V *= COULOMB_CONST
    if bInternals:
        internals = (density_fft, Vw, ker_w)
    else:
        internals = None
    return V, internals

def plot1Dcut(apos, Vg, i0, dg, Ls, iax=0, nPBC=[10,10,10], Vmax=1.0, scErr=100.0):
    ix, iy, iz = i0  # unpack voxel through which we plot the debug cuts
    lvec=np.array( [[Ls[0],0.0,0.0], [0.0,Ls[1],0.0], [0.0,0.0,Ls[2]] ]) # lattice vectors, cubic grid
    # 1D Cut the Ewald potential along the [iax] axis  
    if iax == 0:
        Vgl = Vg[iz, iy, :]  # Along x-axis
    elif iax == 1:
        Vgl = Vg[iz, :, ix]  # Along y-axis
    elif iax == 2:
        Vgl = Vg[:, iy, ix]  # Along z-axis
    else:
        raise ValueError("Invalid axis index. iax must be 0, 1, or 2.")

    #ns = Vg.shape[::-1]  # grid is ordered  [x,y,z], so ns=(nx,ny,nz)
    nt = Vgl.size

    dt = dg[iax]
    Lt = nt*dt

    # Sampling points at which we calculate potential reference using direct sum in real space
    ts = np.linspace( -0.5*Lt, 0.5*Lt, nt, endpoint=False)
    ps = np.zeros((nt, 3))
    ps[:, iax] = ts  # Only vary along the chosen axis

    # Calculate reference potentials using the direct (real space) sum
    fe  = mmff.sampleCoulombPBC(ps, apos, qs, lvec=lvec, nPBC=nPBC)

    vd_vg = fe[:, 3] / Vgl  # Ratio of direct potential (reference) to Ewald potential
    factor = np.average(vd_vg[0:nt//3])  # Average in a safe distance from charges
    print("plot1Dcut() iax=", iax, " nt ", nt, " dt ", dt, "Lt ", Lt, " factor ", factor)
    # Plotting
    plt.plot(ts, Vgl,   '-k', label="V ewald")
    plt.plot(ts, vd_vg, '-', label="ref/ewald")
    plt.plot(ts, fe[:, 3],  '--', label="V direct " + str(nPBC))
    plt.title("cut1D iax=" + str(iax) + " n=" + str(len(ts)) + " tmax=" + str(ts.max()))
    plt.ylim(-2.0, 2.0)
    plt.legend()
    plt.grid()

def test_vs_direct( apos, qs, ns=[80,120,120], dg=[0.1/0.8,0.1/1.2,0.1/1.2], nPBC=[60,30,30], scErr=100.0, order=2, bPython=False ):
    apos = apos.copy()
    print( "apos ", apos )
    Ls = [ ns[0]*dg[0], ns[1]*dg[1], ns[2]*dg[2] ]                                # lenghts along each axis  
    pos0=np.array( [ -0.5*Ls[0], -0.5*Ls[1], -0.5*Ls[2] ] )                       # origin of grid - half of the grid size
    mmff.setupEwaldGrid( ns, dg=dg, pos0=pos0 )                                   # setup grid in C++ for projection
    dens = mmff.projectAtomsEwaldGrid( apos, qs, ns=ns, order=order )             # project atoms to grid using C++ 

    Qtot = np.abs(dens).sum();        print("Qtot = ", Qtot, np.abs(qs).sum() )   # check if the total charge is conserved in projection ( should be same as sum of qs )
    
    if bPython:
        Vg,(density_fft, Vw, ker_w) = compute_potential(dens, dg, bInternals=True )
    else:
        print( "ERROR not implemented yet !!!" )
        exit()

    # set voxel through which we cut the 3D grid Vg in debug plots
    ix0=ns[0]//2
    iy0=ns[1]//2
    iz0=ns[2]//2
    i0 = (ix0,iy0,iz0)

    # ---- 2D debug plots
    # dens2d =  np.sum( dens, axis=iax )
    plt.figure(figsize=(15,10))
    extent  =[  -Ls[0]*0.5  , +Ls[0]*0.5,    -Ls[1]*0.5  ,+Ls[1]*0.5      ]
    extent_w=[  -np.pi/dg[0], +np.pi/dg[0],  -np.pi/dg[1],+np.pi/dg[1]    ]
    # --- plot in real-space (cartesian) coordinates
    plt.subplot(2,3,1  ); plt.imshow( dens [iz0,:,:],                        cmap='bwr', extent=extent  , origin="lower" ); plt.colorbar(); plt.title("Charge Density" )
    plt.subplot(2,3,2  ); plt.imshow( ker_w[iz0,:,:],                        cmap='bwr', extent=extent_w, origin="lower" ); plt.colorbar(); plt.title("ker_w" )
    plt.subplot(2,3,3  ); plt.imshow( Vg   [iz0,:,:],   vmin=-1.0,vmax=1.0,  cmap='bwr', extent=extent  , origin="lower" ); plt.colorbar(); plt.title("V ewald C++/FFTW3 " )
    # --- plot in pixel-coordinates
    plt.subplot(2,3,3+1); plt.imshow( dens [iz0,:,:],                        cmap='bwr'               , origin="lower" ); plt.colorbar(); plt.title("Charge Density" )
    plt.subplot(2,3,3+2); plt.imshow( ker_w[iz0,:,:],                        cmap='bwr'               , origin="lower" ); plt.colorbar(); plt.title("ker_w" )
    plt.subplot(2,3,3+3); plt.imshow( Vg   [iz0,:,:],   vmin=-1.0,vmax=1.0,  cmap='bwr'               , origin="lower" ); plt.colorbar(); plt.title("V ewald C++/FFTW3 " )
    
    # ---- 1D debug plots
    plt.figure(figsize=(10,5))
    plt.subplot(1,2,1); plot1Dcut( apos, Vg, i0=i0, dg=dg, Ls=Ls, iax=0, nPBC=nPBC, scErr=scErr )
    plt.subplot(1,2,2); plot1Dcut( apos, Vg, i0=i0, dg=dg, Ls=Ls, iax=1, nPBC=nPBC, scErr=scErr )
    #plt.subplot(1,3,3); plot1Dcut( Vg, i0=i0, dg=dg, lvec=lvec, iax=2, nPBC=nPBC, scErr=scErr )
    
    plt.tight_layout()
    #plt.show()
    
d=0.6
apos=np.array([
    [-d,.0,0.],
    [+d,.0,0.],
    [0.,-d,0.],
    [0.,+d,0.],
])
qs = [ +1.,+1.,-1.,-1. ]



#test_vs_direct( apos, qs, bPython=True, ns=[100,100,100], dg=[0.10,0.10,0.10] )   # GOOD, This is perfect

# --- change voxel size  homogeneously in all directions
#test_vs_direct( apos, qs, bPython=True, ns=[100,100,100], dg=[0.15,0.15,0.15] )   # GOOD, This is perfect
#test_vs_direct( apos, qs, bPython=True, ns=[100,100,100], dg=[0.07,0.07,0.07] )   # GOOD, This is perfect

# --- change number of grid points homogeneously in all directions
#test_vs_direct( apos, qs, bPython=True, ns=[150,150,150], dg=[0.10,0.10,0.10] )    # GOOD, This is perfect
#test_vs_direct( apos, qs, bPython=True, ns=[80,80,80],    dg=[0.10,0.10,0.10] )    # GOOD, This is perfect

# ========= Changing step size in two directions

# --- change voxel size  homogeneously in xy-directions
#test_vs_direct( apos, qs, bPython=True, ns=[100,100,100], dg=[0.15,0.15,0.10] )    # GOOD, This is perfect
#test_vs_direct( apos, qs, bPython=True, ns=[100,100,100], dg=[0.07,0.07,0.10] )    # GOOD, This is perfect

# --- change voxel size  homogeneously in xz-directions
#test_vs_direct( apos, qs, bPython=True, ns=[100,100,100], dg=[0.07,0.10,0.07] )    # Quite good
#test_vs_direct( apos, qs, bPython=True, ns=[100,100,100], dg=[0.15,0.10,0.15] )    # Quite good

# --- change voxel size  homogeneously in yz-directions
#test_vs_direct( apos, qs, bPython=True, ns=[100,100,100], dg=[0.10,0.07,0.07] )    # Quite good
#test_vs_direct( apos, qs, bPython=True, ns=[100,100,100], dg=[0.10,0.15,0.15] )    # Quite good


# ========= Changing step size in one direction


# --- change voxel size  only in x-directions
#test_vs_direct( apos, qs, bPython=True, ns=[100,100,100], dg=[0.07,0.10,0.10] )    # Quite good
#test_vs_direct( apos, qs, bPython=True, ns=[100,100,100], dg=[0.15,0.10,0.10] )    # Quite good
#test_vs_direct( apos, qs, bPython=True, ns=[100,100,100], dg=[0.20,0.10,0.10] )    # NOT SO PERFECT, iax0 goes to ~1.5 almost

#test_vs_direct( apos, qs, bPython=True, ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[30,30,30] )    # NOT SO PERFECT, iax1 goes to ~1.5 almost
#test_vs_direct( apos, qs, bPython=True, ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30] )    # NOT SO PERFECT, iax1 goes to ~1.5 almost

test_vs_direct( apos, qs, bPython=True, ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[30,30,30], order=3 )    # NOT SO PERFECT, iax1 goes to ~1.5 almost
test_vs_direct( apos, qs, bPython=True, ns=[100,100,100], dg=[0.10,0.20,0.10], nPBC=[60,30,30], order=3 )    # NOT SO PERFECT, iax1 goes to ~1.5 almost

#test_vs_direct( apos, qs, bPython=True, ns=[100,100,100], dg=[0.10,0.20,0.10], order=3 )    # NOT SO PERFECT, iax1 goes to ~1.5 almost
#test_vs_direct( apos, qs, bPython=True, ns=[100,100,100], dg=[0.10,0.10,0.20] )    # GOOD, This is perfect

# --- change voxel size  only in y-directions
#test_vs_direct( apos, qs, bPython=True, ns=[100,100,100], dg=[0.10,0.07,0.10] )    # Quite good
#test_vs_direct( apos, qs, bPython=True, ns=[100,100,100], dg=[0.10,0.15,0.10] )    # Quite good

# --- change voxel size  only in z-directions
#test_vs_direct( apos, qs, bPython=True, ns=[100,100,100], dg=[0.10,0.10,0.07] )    # GOOD, This is perfect
#test_vs_direct( apos, qs, bPython=True, ns=[100,100,100], dg=[0.10,0.10,0.15] )    # GOOD, This is perfect




# ========= Changing number of grid points in different directions

#test_vs_direct( apos, qs, bPython=True, ns=[150,100,100], dg=[0.10,0.10,0.10] )    # 
#test_vs_direct( apos, qs, bPython=True, ns=[100,150,100], dg=[0.10,0.10,0.10] )    # 
#test_vs_direct( apos, qs, bPython=True, ns=[100,100,150], dg=[0.10,0.10,0.10] )    # GOOD, This is perfect


#test_vs_direct( apos, qs, bPython=True, ns=[200,100,100], dg=[0.10,0.10,0.10] )    # NOT SO PERFECT - iax=0 goes up to 1.2
#test_vs_direct( apos, qs, bPython=True, ns=[100,200,100], dg=[0.10,0.10,0.10] )    # NOT SO PERFECT -  iax=1 goes up to 1.2


#test_vs_direct( apos, qs, bPython=True, ns=[50,100,100], dg=[0.10,0.10,0.10] )    # NOT SO PERFECT - iax=1 goes up to 1.2
#test_vs_direct( apos, qs, bPython=True, ns=[100,50,100], dg=[0.10,0.10,0.10] )    # NOT SO PERFECT -  iax=0 goes up to 1.2


#test_vs_direct( apos, qs, bPython=True, ns=[100,100,200], dg=[0.10,0.10,0.10] )   # GOOD, This is perfect
#test_vs_direct( apos, qs, bPython=True, ns=[100,100,50], dg=[0.10,0.10,0.10] )    # GOOD, This is perfect

plt.show()