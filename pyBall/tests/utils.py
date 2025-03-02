import numpy as np
import os
import matplotlib.pyplot as plt
#import matplotlib.patches as patches

COULOMB_CONST  =    14.3996448915 

def getPLQH( R0, E0, a, Q, H ):
    
    e  = np.exp(a*R0);
    cL = e*E0;
    cP = e*cL;
    cH = e*e*H;
    print( "getPLQH cL,cP,e ", cL, cP, e, "  R0,E0,a ", R0, E0, a )
    return np.array([ cP, cL, Q, cH ])

def make_sample_points( p0, t0=0.0, tmax=10.0, dsamp=0.02, iax=2 ):
    ts = np.arange( t0, tmax, dsamp)
    ps = np.zeros( (len(ts), 3) )
    ps[:,0] = p0[0]
    ps[:,1] = p0[1]
    ps[:,2] = p0[2]
    ps[:,iax] = ts
    return ps, ts

def make_sample_points_f4( p0, t0=0.0, tmax=10.0, dsamp=0.02, iax=2 ):
    ts = np.arange( t0, tmax, dsamp)
    ps = np.zeros( (len(ts), 4), dtype=np.float32 )
    ps[:,0] = p0[0]
    ps[:,1] = p0[1]
    ps[:,2] = p0[2]
    ps[:,3] = 0.0
    ps[:,iax] = ts
    return ps, ts

def points_to_ocl(ps_):
    ps = np.zeros( ( len(ps_), 4), dtype=np.float32 )
    ps[:,0] = ps_[:,0]
    ps[:,1] = ps_[:,1]
    ps[:,2] = ps_[:,2]
    ps[:,3] = 0.0
    return ps

def select_cut_1d(V, iax, i0s):
    ix, iy, iz = i0s  # unpack voxel through which we plot the debug cuts
    # 1D Cut the Ewald potential along the [iax] axis  
    if iax == 0:
        return V[iz, iy, :]  # Along x-axis
    elif iax == 1:
        return V[iz, :, ix]  # Along y-axis
    elif iax == 2:
        return V[:, iy, ix]  # Along z-axis
    else:
        raise ValueError("Invalid axis index. iax must be 0, 1, or 2.")
    
def load_potential_comb( path ):
    VPaul = np.load( os.path.join( path, "debug_BsplinePaul_pbc.npy" ) )
    VLond = np.load( os.path.join( path, "debug_BsplineLond_pbc.npy" ) )
    VCoul = np.load( os.path.join( path, "debug_BsplineCoul_pbc.npy" ) )
    #print("VCoul.shape = ", VCoul.shape)
    #print("VLond.shape = ", VLond.shape)
    #print("VPaul.shape = ", VPaul.shape)
    VPLQ = np.zeros(VCoul.shape+(4,), dtype=np.float32)
    VPLQ[:,:,:,0] = VPaul
    VPLQ[:,:,:,1] = VLond
    VPLQ[:,:,:,2] = VCoul
    VPLQ[:,:,:,3] = 0.0
    #print("VPLQ.shape(BEFORE) = ", VPLQ.shape)
    VPLQ = VPLQ.transpose( (2,1,0,3) ).copy()
    #sh = VPLQ.shape; print("sh = ", sh)
    #print("VPLQ.shape(AFTER) = ", VPLQ.shape)
    return VPLQ

def soft_clamp(y, dy=None, y1=100, y2=200 ):
    """
    Applies a soft clamp to y, smoothly transitioning values above y1 towards y2.
    Also computes the derivative dy accordingly using the chain rule.

    Parameters:
    - y: np.ndarray, input values to clamp
    - dy: np.ndarray, derivatives of y with respect to some variable x
    - y1: float, lower threshold for clamping
    - y2: float, upper threshold for clamping

    Returns:
    - y_new: np.ndarray, clamped y values
    - dy_new: np.ndarray, updated derivatives
    """    
    mask   = y > y1 
    y12    = y2 - y1
    invdy  = 1.0 / y12
    z = (y[mask] - y1) * invdy
    y[mask]   = y1 + y12 * (1 - 1 / (1 + z))
    if dy is not None:
        dy[mask] *= 1.0 / (1.0 + z)**2
    return y, dy

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

def plot1Dcut(apos, qs, Vg, i0, dg, Ls, iax=0, nPBC=[10,10,10], vmax=5.0, scErr=100.0, Vref=None, mmff=None ):
    #ix, iy, iz = i0  # unpack voxel through which we plot the debug cuts
    lvec=np.array( [[Ls[0],0.0,0.0], [0.0,Ls[1],0.0], [0.0,0.0,Ls[2]] ]) # lattice vectors, cubic grid
    Vgl = select_cut_1d( Vg, iax, i0 )

    #ns = Vg.shape[::-1]  # grid is ordered  [x,y,z], so ns=(nx,ny,nz)
    nt = Vgl.size

    dt = dg[iax]
    Lt = nt*dt

    # Sampling points at which we calculate potential reference using direct sum in real space
    ts = np.linspace( -0.5*Lt, 0.5*Lt, nt, endpoint=False)
    ps = np.zeros((nt, 3))
    ps[:, iax] = ts  # Only vary along the chosen axis

    # Calculate reference potentials using the direct (real space) sum
    if Vref is None:
        fe  = mmff.sampleCoulombPBC(ps, apos, qs, lvec=lvec, nPBC=nPBC)
        Vref = fe[:, 3]
    else:
        Vref = select_cut_1d( Vref, iax, i0)

    vd_vg = Vref/ Vgl  # Ratio of direct potential (reference) to Ewald potential
    factor = np.average(vd_vg[0:nt//3])  # Average in a safe distance from charges
    #print("plot1Dcut() iax=", iax, " nt ", nt, " dt ", dt, "Lt ", Lt, " factor ", factor)
    # Plotting
    plt.plot(ts, Vgl,   '-k', label="V ewald")
    plt.plot(ts, vd_vg, '-',  label="ref/ewald")
    plt.plot(ts, Vref,  '--', label="Vref " + str(nPBC))
    plt.title("cut1D iax=" + str(iax) + " n=" + str(len(ts)) + " tmax=" + str(ts.max()))
    plt.ylim(-vmax, vmax )
    plt.legend()
    plt.grid()