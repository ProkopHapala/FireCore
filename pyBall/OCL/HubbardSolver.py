import os
import numpy as np
import pyopencl as cl

from .OpenCLBase import OpenCLBase


kBoltzmann = 8.61733326214511e-5  # eV/K


default_params = {
    "tipDecay"  : 0.3,
    "tipRadius" : 3.0,
    "zTip"      : 5.0,
    "zMirror"   : -0.5,
    "Wcouple"   : 1.0,
}


class HubbardSolver(OpenCLBase):
    """Convenient wrapper around the `solve_minBrute_fly` kernel in `cl/hubbard.cl`.

    Responsibilities:
        1.  Load and build the OpenCL program.
        2.  Handle GPU buffer allocation / re-allocation.
        3.  Provide a single high-level `solve` method that:
                * uploads input data (site positions + tip positions),
                * launches the kernel,
                * downloads the results (Emin, iMin, Itot).

    The class deliberately keeps most operations low-level (explicit buffer
    management & kernel argument setting) to stay transparent while still
    hiding repetitive boiler-plate from the user.
    """

    # ---- names used in buffer_dict / as object attributes
    # buf_pose  = "posE_buff"      #   nSingle  * sizeof(float4)
    # buf_pTips = "pTips_buff"     #   nTips    * sizeof(float4)
    # buf_Emin  = "Emin_buff"      #   nTips    * sizeof(float )
    # buf_iMin  = "iMin_buff"      #   nTips    * sizeof(int  )
    # buf_Itot  = "Itot_buff"      #   nTips    * sizeof(float2)

    def __init__(self, nloc: int = 32, device_index: int = 0):
        super().__init__(nloc=nloc, device_index=device_index)

        # Build the OpenCL program ------------------------------------------------
        base_path = os.path.dirname(os.path.abspath(__file__))
        rel_path  = "cl/hubbard.cl"
        if not self.load_program(rel_path=rel_path, base_path=base_path):
            raise RuntimeError("[HubbardSolver] Failed to load/compile hubbard.cl")

    # ---------------------------------------------------------------------
    # Internal helpers
    # ---------------------------------------------------------------------
    def realloc_buffers(self, nSingle: int, nTips: int):
        """(Re-)allocate GPU buffers as needed for the given problem size."""
        print("HubbardSolver::realloc_buffers() dir(self):")
        sz_f  = 4
        #sz_f4  = 16  # bytes
        #sz_f2  = 8
        buffs = {
            #  name     handle, size 
            "posE":  [ None, nSingle * sz_f*4 ],
            "pTips": [ None, nTips   * sz_f*4 ],
            "Emin":  [ None, nTips   * sz_f*1 ],
            "iMin":  [ None, nTips   * sz_f*1 ],
            "Itot":  [ None, nTips   * sz_f*2 ],
        }
        for name, buff_info  in buffs.items():
            handle=buff_info[0]
            sz=buff_info[1]
            #print("  ",name, sz)
            buffname = name + "_buff"
            if handle is None:
                handle = self.try_make_buff(buffname, sz)
            else:
                if handle.size != sz:
                    handle = self.try_make_buff(buffname, sz)
            buff_info[0] = handle
    # ---------------------------------------------------------------------
    # Public API
    # ---------------------------------------------------------------------
    
    def setup_solve_minBrute_fly(self, nSingle, nTips, tipDecay, tipRadius, zMirror, Wcouple):
        kernel = self.prg.solve_minBrute_fly
        args = (
            np.int32(nSingle),
            self.posE_buff,
            np.int32(nTips),
            self.pTips_buff,
            np.float32(tipDecay),
            np.float32(tipRadius),
            np.float32(zMirror),
            np.float32(Wcouple),
            self.Emin_buff,
            self.iMin_buff,
            self.Itot_buff,
        )
        return args, kernel
    
    def setup_solve_minBrute_boltzmann(self, nSingle, nTips, tipDecay, tipRadius, zMirror, Wcouple, bBoltzmann):
        kernel = self.prg.solve_minBrute_boltzmann
        args = (
            np.int32(nSingle),
            self.posE_buff,
            np.int32(nTips),
            self.pTips_buff,
            np.float32(tipDecay),
            np.float32(tipRadius),
            np.float32(zMirror),
            np.float32(Wcouple),
            np.float32(bBoltzmann),
            self.Emin_buff,
            self.iMin_buff,
            self.Itot_buff,
        )
        return args, kernel

    def solve(
        self,
        posE : np.ndarray = None,  # shape (nSingle,4), float32  (xyz , E0)
        pTips: np.ndarray = None,  # shape (nTips, 4), float32  (xyz , Vbias)
        bRealloc:    bool = True,
        nSingle           = None,
        nTips             = None,
        bBoltzmann:  bool = False,
        params: dict      = default_params,
    ):
        """Run the brute-force Hubbard solver for a batch of tip positions.

        Returns
        -------
        Emin : (nTips,)  float32  ground-state energy for each tip
        iMin : (nTips,)  int32    bit mask of the ground-state occupation
        Itot : (nTips,2) float32  (I_occ, I_unoc) tunnelling currents
        """
        assert posE.dtype  == np.float32 and posE.ndim  == 2 and posE.shape[1]  == 4
        assert pTips.dtype == np.float32 and pTips.ndim == 2 and pTips.shape[1] == 4

        if nSingle is None: nSingle = posE.shape[0]
        if nTips   is None: nTips   = pTips.shape[0]
        if nSingle > 32: raise ValueError("[HubbardSolver] nSingle > 32 exceeds kernel limit")

        # (Re)allocate GPU memory & upload inputs --------------------------------
        if bRealloc: self.realloc_buffers(nSingle, nTips)

        # print("HubbardSolver::solve() dir(self):")
        # for name in dir(self): 
        #     #if name[0] != "_": print("  ",name)
        #     if "_buff" in name: print(name)
        posE_buff  = self.posE_buff
        pTips_buff = self.pTips_buff
        Emin_buff  = self.Emin_buff
        iMin_buff  = self.iMin_buff
        Itot_buff  = self.Itot_buff

        if posE  is not None: self.toGPU_( posE_buff,  posE)
        if pTips is not None: self.toGPU_( pTips_buff, pTips)

        tipDecay  = params["tipDecay"]
        tipRadius = params["tipRadius"]
        zMirror   = params["zMirror"]
        Wcouple   = params["Wcouple"]
        
        # Kernel launch (construct arg list only once per call) ------------------
        global_size = (nTips * self.nloc,)
        local_size  = (self.nloc,)
        if(bBoltzmann):
            bBoltzmann   = 1./(kBoltzmann*params["Temperature"])
            args, kernel = self.setup_solve_minBrute_boltzmann(nSingle, nTips, tipDecay, tipRadius, zMirror, bBoltzmann)
        else:
            args, kernel = self.setup_solve_minBrute_fly(nSingle, nTips, tipDecay, tipRadius, zMirror, Wcouple)
        kernel(self.queue, global_size, local_size, *args)
        
        # Download results -------------------------------------------------------
        Emin = np.empty(nTips, dtype=np.float32)
        iMin = np.empty(nTips, dtype=np.int32)
        Itot = np.empty((nTips, 2), dtype=np.float32)
        self.fromGPU_( Emin_buff, Emin)
        self.fromGPU_( iMin_buff, iMin)
        self.fromGPU_( Itot_buff, Itot)

        self.queue.finish()

        return Emin, iMin, Itot

    def setup_eval_coupling_matrix(self, nSingle, zMirror):
        """Prepare arguments for eval_coupling_matrix kernel."""
        kernel = self.prg.eval_coupling_matrix
        args = (
            np.int32(nSingle),
            self.posE_buff,
            self.Rij_buff,
            self.Wij_buff,
            np.float32(zMirror)
        )
        return args, kernel

    def eval_coupling_matrix(self, posE: np.ndarray, zMirror: float, bRealloc: bool = True):
        """Evaluate the full coupling matrix between all sites.
        
        Args:
            posE: (nSingle,4) array of site positions and energies (x,y,z,E0)
            zMirror: z-coordinate of mirror plane
            bRealloc: whether to reallocate buffers if needed
            
        Returns:
            Rij: (nSingle,nSingle) distance matrix
            Wij: (nSingle,nSingle) coupling matrix
        """
        nSingle = len(posE)
        nTotal = nSingle * nSingle
        
        # Allocate all buffers needed
        if bRealloc:
            self.try_make_buff("posE_buff", posE.nbytes)
            self.try_make_buff("Rij_buff", nTotal * 4)
            self.try_make_buff("Wij_buff", nTotal * 4)
        
        # Upload inputs
        self.toGPU_(self.posE_buff, posE)
        
        # Launch kernel
        global_size = (nTotal,)
        args, kernel = self.setup_eval_coupling_matrix(nSingle, zMirror)
        kernel(self.queue, global_size, None, *args)
        
        # Download results
        Rij = np.empty(nTotal, dtype=np.float32)
        Wij = np.empty(nTotal, dtype=np.float32)
        self.fromGPU_(self.Rij_buff, Rij)
        self.fromGPU_(self.Wij_buff, Wij)
        
        return Rij.reshape(nSingle, nSingle), Wij.reshape(nSingle, nSingle)

# -----------------------------------------------------------------------------
# Utility helpers (stand-alone, not bound to the class) -------------------------
# -----------------------------------------------------------------------------

def make_site_ring( n, R, E0=-0.1, ang0=0.0 ):
    """Make a ring of n sites at radius R and energy E0."""
    pos = np.zeros((n,4), dtype=np.float32)
    ang = np.linspace(0, 2*np.pi, n, endpoint=False, dtype=np.float32) + ang0
    pos[:,0] = R * np.cos(ang)
    pos[:,1] = R * np.sin(ang)
    pos[:,2] = 0.0
    pos[:,3] = E0
    return pos

def save_sites_to_txt(path: str, posE: np.ndarray):
    with open(path, "w") as f:
        for p in posE: 
            l=f"{p[0]:12.6f} {p[1]:12.6f} {p[2]:12.6f} {p[3]:12.6f}\n"
            print(l,end="")
            f.write(l)

def load_sites_from_txt(path: str) -> np.ndarray:
    """Load site parameters (x y z E0) from a whitespace-separated text file.

    The file must have four columns per line.  Returns `float32 (nSites,4)`.
    """
    arr = np.loadtxt(path, dtype=np.float32)
    if arr.ndim == 1:  # single line -> shape (4,) ; make (1,4)
        arr = arr[None, :]
    if arr.shape[1] != 4:
        raise ValueError("File must contain exactly four columns (x y z E0)")
    return arr.astype(np.float32)


def generate_xy_scan(extent, nxy=(10,10), zTip=5.0, Vbias=0.0) -> np.ndarray:
    """Generate a 2-D grid of tip positions at fixed z & Vbias.

    Returns array of shape (nx*ny,4) with columns (x,y,z,Vbias).
    """
    xmin, xmax, ymin, ymax = extent
    nx,ny = nxy
    xs = np.linspace(xmin, xmax, nx, dtype=np.float32)
    ys = np.linspace(ymin, ymax, ny, dtype=np.float32)
    X, Y = np.meshgrid(xs, ys, indexing='ij')
    #pts = np.stack([X.ravel(), Y.ravel(), np.full(X.size, z, np.float32), np.full(X.size, Vbias, np.float32)], axis=1)
    pts = np.empty( (nx*ny,4), dtype=np.float32)
    pts[:,0] = X.flat
    pts[:,1] = Y.flat
    pts[:,2] = zTip
    pts[:,3] = Vbias
    return pts


def generate_xV_scan(p1, p2, Vrange=(0.0,1.0), nxV=(10,10)) -> np.ndarray:
    """Generate tips along a line p1->p2 (included) and a range of Vbias.

    p1, p2 : 3-component sequences giving XYZ.
    Returns array of shape (nx*nV,4).
    """
    p1 = np.asarray(p1, dtype=np.float32)
    p2 = np.asarray(p2, dtype=np.float32)
    V1, V2 = Vrange
    nx, nV = nxV
    ts = np.linspace(0.0, 1.0, nx, dtype=np.float32)[:, None]
    line_pts = p1 + ts * (p2 - p1)  # shape (nx,3)
    Vvals = np.linspace(V1, V2, nV, dtype=np.float32)
    # Cartesian product of points and Vbias
    line_rep = np.repeat(line_pts, nV, axis=0)            # (nx*nV,3)
    V_rep    = np.tile(Vvals,  nx)[:, None]                # (nx*nV,1)
    pts = np.hstack([line_rep, V_rep])                     # (nx*nV,4)
    return pts.astype(np.float32)

# -----------------------------------------------------------------------------
# High-level convenience wrappers ----------------------------------------------
# -----------------------------------------------------------------------------

def solve_xy_scan(solver: HubbardSolver, posE: np.ndarray, extent, nxy=(10,10), params = default_params):
    """Run an XY scan and reshape outputs to (nx,ny)."""
    zTip             = params["zTip"]
    Vbias            = params["Vbias"]
    pTips            = generate_xy_scan(extent, nxy, zTip=zTip, Vbias=Vbias)
    Emin, iMin, Itot = solver.solve(posE, pTips, params=params)
    return Emin.reshape(nxy), iMin.reshape(nxy), Itot.reshape(nxy + (2,))

def solve_xV_scan(solver: HubbardSolver, posE: np.ndarray, p1, p2, Vrange=(0.0,1.0), nxV=(10,10), params=default_params ):
    """Run an X-V scan (position along a line vs bias).  Outputs shaped (nx,nV)."""
    pTips = generate_xV_scan(p1, p2, Vrange, nxV)
    pTips[:,2] = params["zTip"]
    Emin, iMin, Itot = solver.solve(posE, pTips, params=params)
    return Emin.reshape(nxV), iMin.reshape(nxV), Itot.reshape(nxV + (2,))

def solve_pSites(params, posE: np.ndarray, solver: HubbardSolver,  ):
    """
    Debugging function that runs the solver with tip positions matching the sites.
    Returns:
        Emin : (nSites,) ground-state energies
        iMin : (nSites,) occupation patterns
        Itot : (nSites,2) tunneling currents
    """
    # Create tips at site positions with given Vbias
    pTips = posE.copy()
    pTips[:,3] = params["Vbias"]  # Replace E0 with Vbias
    pTips[:,2] = params["zTip"]
    Emin, iMin, Itot = solver.solve(posE, pTips, params=params)
    for i in range(len(Emin)):
        print(f"Site {i}: Emin={Emin[i]:.6f}, iMin={iMin[i]:b}, Itot={Itot[i]} pTips={pTips[i]}")
    return Emin, iMin, Itot

def plot_sites(posE, ax=None, c="r", ms=2.0):
    if ax is None:
        ax = plt.gca()
    ax.plot(posE[:,0], posE[:,1], '.', color=c, markersize=ms)
    ax.set_aspect('equal')
    return ax

def plot2d( data, extent=None, title=None, ax=None, cmap=None, ps=None):
    if ax is None:
        ax = plt.gca()
    im0 = ax.imshow(data, extent=extent, cmap=cmap);        
    fig.colorbar(im0, ax=ax); 
    plot_sites(ps,ax); 
    ax.set_title(title);

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    params={ 
        "tipDecay"  : 0.3,
        "tipRadius" : 3.0,
        "zTip"      : 5.0,
        "zMirror"   : -0.5,
        "Wcouple"   : 1.0,
        "E0"        : -0.1,
        #"Vbias"     : +0.2,
        #"Vbias"     : -0.45,
        "Vbias"     : -0.55,
        "nSingle"   : 3,
        "nTips"     : 1,
    }


    # trimer
    #posE = make_site_ring(3, 5.0, E0=params["E0"] )
    #save_sites_to_txt("trimer.txt", posE)

    # # hexamer
    outer = make_site_ring(3, 5.0, -0.1        )
    inner = make_site_ring(3, 3.0, -0.1, np.pi )
    posE  = np.vstack([outer, inner])
    save_sites_to_txt("hexamer.txt", posE)

    solver = HubbardSolver()

    #solve_pSites(params, posE, solver)

    Rij, Wij = solver.eval_coupling_matrix(posE, params["zMirror"])
    print("Rij:\n", Rij)
    print("Wij:\n", Wij)

    extent=[-10.0,10.0,-10.0,10.0]
    # #posE = load_sites_from_txt("hexamer.txt")
    Emin, iMin, Itot = solve_xy_scan(solver, posE, extent=extent, nxy=(200,200), params=params)
    

    fig,axs = plt.subplots(2,2)
    # im0 = axs[0,0].imshow(Emin,        extent=extent);        fig.colorbar(im0, ax=axs[0,0]); plot_sites(posE,axs[0,0]); axs[0,0].set_title("Emin");
    # im1 = axs[0,1].imshow(iMin,        extent=extent);        fig.colorbar(im1, ax=axs[0,1]); plot_sites(posE,axs[0,1]); axs[0,1].set_title("iMin");
    # im2 = axs[1,0].imshow(Itot[:,:,0], extent=extent); fig.colorbar(im2, ax=axs[1,0]); plot_sites(posE,axs[1,0]); axs[1,0].set_title("I_occ");
    # im3 = axs[1,1].imshow(Itot[:,:,1], extent=extent); fig.colorbar(im3, ax=axs[1,1]); plot_sites(posE,axs[1,1]); axs[1,1].set_title("I_unocc");
    plot2d(Emin.T,        extent=extent, title="Emin", ax=axs[0,0], cmap="viridis", ps=posE);
    plot2d(iMin.T,        extent=extent, title="iMin", ax=axs[0,1], cmap="viridis", ps=posE);
    plot2d(Itot[:,:,0].T, extent=extent, title="I_occ", ax=axs[1,0], cmap="viridis", ps=posE);
    plot2d(Itot[:,:,1].T, extent=extent, title="I_unocc", ax=axs[1,1], cmap="viridis", ps=posE);
    plt.tight_layout()
    plt.show()
