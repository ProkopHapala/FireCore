
# run it like this:
#   python -u -m pyBall.OCL.HubbardSolver | tee OUT

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


        # NEW: A dedicated realloc function for the new local update solver's buffers
    def realloc_local_update_buffers(self, nSingle: int, nTips: int, nMaxNeigh: int):
        """(Re-)allocate GPU buffers specifically for the local update solver."""
        sz_f  = 4
        # Reuse output buffers from brute-force where possible (E_out, Itot_out)
        buffs = {
            "Esite_Thop_buff": (nTips * nSingle * sz_f * 2, None),
            "W_val_buff":      (nSingle * nMaxNeigh * sz_f * 1, None),
            "W_idx_buff":      (nSingle * nMaxNeigh * sz_f * 1, None),
            "nNeigh_buff":     (nSingle * sz_f * 1, None),
            # Re-purpose brute-force outputs
            "E_out_buff":      (nTips * sz_f * 1, "Emin_buff"),
            "state_out_buff":  (nTips * 8, "iMin_buff"), # ulong is 8 bytes. Repurpose iMin
            "Itot_out_buff":   (nTips * sz_f * 2, "Itot_buff"),
        }
        for buffname, (sz, alias) in buffs.items():
            final_name = alias if alias else buffname
            if getattr(self, final_name, None) is None or getattr(self, final_name).size != sz:
                self.try_make_buff(final_name, sz)

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

    def setup_solve_local_updates_parallel(self, nSite, nTips, params):
        """Prepare arguments for the solve_local_updates_parallel kernel."""
        kernel = self.prg.solve_local_updates_parallel
        
        solver_params = np.array([
            params.get("kT", 0.001),
            params.get("nIter", 10000),
            params.get("solverMode", 2),
            params.get("seed", np.random.randint(0, 2**31))
        ], dtype=np.float32)

        args = (
            np.int32(nSite),
            np.int32(nTips),
            solver_params,
            self.Esite_Thop_buff,
            self.W_val_buff,
            self.W_idx_buff,
            self.nNeigh_buff,
            self.iMin_buff, # Re-purposed for state_out (ulong)
            self.Emin_buff, # Re-purposed for E_out (float)
            self.Itot_buff  # Re-purposed for Itot_out (float2)
        )
        return args, kernel


    # NEW: High-level public API for the local update solver
    def solve_local_updates(self, Esite_Thop, W_sparse, params=default_params):
        """
        Run the local-update Monte Carlo solver.

        Args:
            Esite_Thop (np.ndarray): Shape (nTips, nSingle, 2), float32.
                                     Contains {Esite, Thop} for each tip-site pair.
            W_sparse (tuple): A tuple (W_val, W_idx, nNeigh) representing the sparse
                              interaction matrix, typically from `make_sparse_W`.
            params (dict): Dictionary with solver parameters like kT, nIter, etc.

        Returns:
            E_out (np.ndarray): (nTips,) float32 - final energy for each simulation.
            state_out (np.ndarray): (nTips,) uint64 - final state bitmask.
            Itot_out (np.ndarray): (nTips, 2) float32 - final currents {I_occ, I_unoc}.
        """
        nTips, nSite, _ = Esite_Thop.shape
        W_val, W_idx, nNeigh = W_sparse
        nMaxNeigh = W_val.shape[1]

        if nSite > 64:
            raise ValueError("[HubbardSolver] nSite > 64 exceeds ulong bitmask limit.")

        # Reallocate buffers and upload data
        self.realloc_local_update_buffers(nSite, nTips, nMaxNeigh)
        self.toGPU_(self.Esite_Thop_buff, Esite_Thop.reshape(-1, 2)) # Reshape to (nTips*nSingle, 2)
        self.toGPU_(self.W_val_buff, W_val)
        self.toGPU_(self.W_idx_buff, W_idx)
        self.toGPU_(self.nNeigh_buff, nNeigh)

        # Kernel launch
        global_size = (nTips * self.nloc,)
        local_size = (self.nloc,)
        args, kernel = self.setup_solve_local_updates_parallel(nSite, nTips, params)
        kernel(self.queue, global_size, local_size, *args)
        
        # Download results
        state_out = np.empty(nTips, dtype=np.uint64)
        E_out = np.empty(nTips, dtype=np.float32)
        Itot_out = np.empty((nTips, 2), dtype=np.float32)
        
        self.fromGPU_(self.iMin_buff, state_out) # Re-purposed
        self.fromGPU_(self.Emin_buff, E_out)     # Re-purposed
        self.fromGPU_(self.Itot_buff, Itot_out)  # Re-purposed

        self.queue.finish()
        return E_out, state_out, Itot_out
    


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

    def setup_eval_oriented_hopping(self, nSingle, nTips):
        """Prepare arguments for eval_Oriented_Hopping kernel."""
        kernel = self.prg.eval_Oriented_Hopping
        args = (
            np.int32(nSingle),
            np.int32(nTips),
            self.posE_buff,
            self.rots_buff,
            self.orbs_buff,
            self.pTips_buff,
            self.Tout_buff
        )
        return args, kernel

    def eval_oriented_hopping(self, posE: np.ndarray, rots: np.ndarray, orbs: np.ndarray, pTips: np.ndarray, bRealloc: bool = True):
        """Evaluate tunneling amplitudes for each site-tip pair.
        
        Args:
            posE:       (nSingle,4) array of site positions and energies
            rots:       (nSingle,2)  array of unit rotation (cos,sin) per site
            orbs:       (nSingle,4)  array of orbital params (C_s, C_px, C_py, decay)
            pTips:      (nTips,4)    array of tip positions
            bRealloc:   whether to reallocate buffers
        Returns:
            Tout:       (nSingle, nTips) tunneling amplitudes
        """
        nSingle = posE.shape[0]
        nTips   = pTips.shape[0]
        total   = nSingle * nTips
        print("HubbardSolver::eval_oriented_hopping() nSingle: ", nSingle, " nTips: ", nTips, " total: ", total)

        if bRealloc:
            self.try_make_buff("posE_buff",  posE.nbytes)
            self.try_make_buff("rots_buff",  rots.nbytes)
            self.try_make_buff("orbs_buff",  orbs.nbytes)
            self.try_make_buff("pTips_buff", pTips.nbytes)
            self.try_make_buff("Tout_buff",  total * 4)

        # Upload data
        self.toGPU_(self.posE_buff,  posE)
        self.toGPU_(self.rots_buff,  rots)
        self.toGPU_(self.orbs_buff,  orbs)
        self.toGPU_(self.pTips_buff, pTips)

        global_size = (total,)
        args, kernel = self.setup_eval_oriented_hopping(nSingle, nTips)
        kernel(self.queue, global_size, None, *args)
        self.queue.finish()

        Tout = np.empty((nTips,nSingle), dtype=np.float32)
        self.fromGPU_(self.Tout_buff, Tout)
        return Tout

# -----------------------------------------------------------------------------
# Utility helpers (stand-alone, not bound to the class) -------------------------
# -----------------------------------------------------------------------------


# NEW: Helper function to create the sparse interaction matrix from a dense one.
def make_sparse_W(Wij_dense, nMaxNeigh=16):
    """
    Converts a dense (nSite, nSite) interaction matrix to the sparse format
    required by the local_updates kernel.

    Args:
        Wij_dense (np.ndarray): The full, dense interaction matrix.
        nMaxNeigh (int): The maximum number of neighbors to store per site.

    Returns:
        tuple: (W_val, W_idx, nNeigh) ready for GPU upload.
            - W_val (np.ndarray): (nSite, nMaxNeigh) float32, interaction values.
            - W_idx (np.ndarray): (nSite, nMaxNeigh) int32, neighbor indices.
            - nNeigh (np.ndarray): (nSite,) int32, number of actual neighbors.
    """
    nSite = Wij_dense.shape[0]
    
    # Use argsort to find the indices of the largest couplings for each site
    # We negate the matrix because argsort sorts in ascending order.
    sorted_indices = np.argsort(-np.abs(Wij_dense), axis=1)

    # Initialize output arrays
    W_val  = np.zeros((nSite, nMaxNeigh), dtype=np.float32)
    W_idx  = np.zeros((nSite, nMaxNeigh), dtype=np.int32)
    nNeigh = np.zeros(nSite, dtype=np.int32)

    for i in range(nSite):
        count = 0
        for j_idx in sorted_indices[i]:
            if i == j_idx or count >= nMaxNeigh:
                continue # Skip self-interaction and excess neighbors
            
            W_idx[i, count] = j_idx
            W_val[i, count] = Wij_dense[i, j_idx]
            count += 1
        nNeigh[i] = count
        
    return W_val, W_idx, nNeigh

def make_site_ring( n, R, E0=-0.1, ang0=0.0, z=0.0, ang2=0.0 ):
    """Make a ring of n sites at radius R and energy E0."""
    pos = np.zeros((n,4), dtype=np.float32)
    ang = np.linspace(0, 2*np.pi, n, endpoint=False, dtype=np.float32) + ang0
    pos[:,0] = R * np.cos(ang)
    pos[:,1] = R * np.sin(ang)
    pos[:,2] = z
    pos[:,3] = E0
    return pos, ang+ang2

def make_grid_sites( nxy=(4,4),  avec=(1.0,0.0), bvec=(0.0,1.0), z=0.0, E0=-0.1, bCenter=True ):
    """Make a grid of n sites at radius R and energy E0."""
    ntot = np.prod(nxy)
    pos = np.zeros((ntot,4), dtype=np.float32)
    nx,ny = nxy
    for ix in range(ny):
        for jy in range(nx):
            pos[ix*nx+jy,0] = avec[0]*ix + bvec[0]*jy
            pos[ix*nx+jy,1] = avec[1]*ix + bvec[1]*jy
            pos[ix*nx+jy,2] = z
            pos[ix*nx+jy,3] = E0
    if bCenter:
        pos[:,:2] -= np.mean(pos, axis=0)[:2]
    return pos


def save_sites_to_txt(path: str, posE: np.ndarray):
    with open(path, "w") as f:
        for p in posE: 
            #l=f"{p[0]:12.6f} {p[1]:12.6f} {p[2]:12.6f} {p[3]:12.6f}\n"
            l = " ".join([f"{x:12.6f}" for x in p]) + "\n"
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

def plot_sites(posE, ax=None, c="r", ms=2.0,angles=None):
    if ax is None: 
        ax = plt.gca()
    ax.plot(posE[:,0], posE[:,1], '.', color=c, markersize=ms)
    if angles is not None: # with  quiver
        ax.quiver(posE[:,0], posE[:,1], np.cos(angles), np.sin(angles), color=c)
    ax.set_aspect('equal')
    return ax

def plot2d( data, extent=None, title=None, ax=None, cmap=None, ps=None):
    if ax is None:
        ax = plt.gca()
    im0 = ax.imshow(data, extent=extent, cmap=cmap, origin="lower");        
    ax.figure.colorbar(im0, ax=ax) 
    plot_sites(ps,ax); 
    ax.set_title(title);

def solve_hopping_scan(solver: HubbardSolver, posE: np.ndarray, rots: np.ndarray, orbs: np.ndarray, extent, nxy=(10,10), params=default_params):
    """Run oriented hopping scan over XY grid and sum amplitudes over sites."""
    pTips = generate_xy_scan(extent, nxy, zTip=params["zTip"], Vbias=params["Vbias"])
    nSingle = posE.shape[0]
    Tout = solver.eval_oriented_hopping(posE, rots, orbs, pTips, bRealloc=True)
    return Tout.reshape(nxy+(nSingle,))


def test_brute_force_solver():
    # trimer
    #posE = make_site_ring(3, 5.0, E0=params["E0"] )
    #save_sites_to_txt("trimer.txt", posE)

    # hexamer
    params={ 
        "tipDecay"  : 0.3,
        "tipRadius" : 3.0,
        "zTip"      : 5.0,
        "zMirror"   : -0.5,
        "Wcouple"   : 1.0,
        "E0"        : -0.10,
        #"Vbias"     : +0.3,
        #"Vbias"     : +0.2,
        "Vbias"     : +0.1,
        #"Vbias"     : -0.45,    
        #"Vbias"     : -0.55,
    }
    #outer,a1 = make_site_ring(3, 8.0, E0=params["E0"], ang0=np.pi*0.5 )
    #inner,a2 = make_site_ring(3, 5.2, E0=params["E0"], ang0=np.pi*1.5 )

    outer,a1 = make_site_ring(3, 0.86602540378*(2.0/3.0)*20.0, E0=params["E0"], ang0=np.pi*0.5 )
    inner,a2 = make_site_ring(3, 0.86602540378*(2.0/3.0)*10.0, E0=params["E0"], ang0=np.pi*1.5 )

    print( "|p1-p0| ",  np.sqrt( ((inner[0,:]-inner[1,:])**2).sum() ) )

    #outer = make_site_ring(3, 5.0, E0=params["E0"], ang0=0.0 )
    #inner = make_site_ring(3, 3.0, E0=params["E0"], ang0=np.pi )
    posE  = np.vstack([outer, inner])
    angles = np.hstack([a1,a2])
    print("posE.shape", posE.shape)
    #save_sites_to_txt("hexamer.txt", posE)
    #np.savetxt("hexamer.txt", posE)
    xya      = posE.copy()
    xya[:,2] = angles
    save_sites_to_txt("hexamer.txt", xya)
    #np.savetxt("hexamer.txt", xya)

    plot_sites(posE, angles=angles)
    #plt.show(); exit(0)
    



    # # ---- grid
    # params={ 
    #     "tipDecay"  : 0.3,
    #     "tipRadius" : 3.0,
    #     "zTip"      : 1.0,
    #     "zMirror"   : -0.5,
    #     "Wcouple"   : 1.0,
    #     "E0"        : -0.15,
    #     "Vbias"     : 0.2,
    # }
    # posE = make_grid_sites( nxy=(4,4), avec=(7.0,0.0), bvec=(0.0,5.0), E0=params["E0"] )
    # print("posE:\n", posE)
    # save_sites_to_txt("grid.txt", posE)

    solver = HubbardSolver()
    #solve_pSites(params, posE, solver)

    Rij, Wij = solver.eval_coupling_matrix(posE, params["zMirror"])
    #print("Rij:\n", Rij)
    #print("Wij:\n", Wij)

    extent=[-10.0,10.0,-10.0,10.0]
    # #posE = load_sites_from_txt("hexamer.txt")
    # Oriented hopping demo
    # Create random site rotations and orbitals for demonstration
    nSingle = posE.shape[0]
    # Simple rots: all zeros (no rotation)
    phis = np.zeros(nSingle)
    rots = np.vstack([np.cos(phis), np.sin(phis)]).T.astype(np.float32)
    # Simple orbitals: s-wave only with decay=1.0
    orbs = np.zeros((nSingle,4), dtype=np.float32)
    orbs[:,0] = 1.0  # C_s
    orbs[:,3] = 0.2  # decay
    Tsites = solve_hopping_scan(solver, posE, rots, orbs, extent=extent, nxy=(10,10), params=params)
    print("Tsites: shape: ", Tsites.shape, " min: ", Tsites.min(), " max: ", Tsites.max())
    # Plot hopping amplitude scan
    
    fig2, axs = plt.subplots( 1,nSingle, figsize=(nSingle*4,4) )
    for i in range(nSingle):
        plot2d(Tsites[i], extent=extent, title="Total Hopping Amplitude", ax=axs[i], cmap="inferno", ps=posE)
    # Also run XY scan for comparison plots
    
    plt.show(); exit()
    
    Emin, iMin, Itot = solve_xy_scan(solver, posE, extent=extent, nxy=(200,200), params=params)
    fig,axs = plt.subplots(2,2)
    # im0 = axs[0,0].imshow(Emin,        extent=extent);        fig.colorbar(im0, ax=axs[0,0]); plot_sites(posE,axs[0,0]); axs[0,0].set_title("Emin");
    # im1 = axs[0,1].imshow(iMin,        extent=extent);        fig.colorbar(im1, ax=axs[0,1]); plot_sites(posE,axs[0,1]); axs[0,1].set_title("iMin");
    # im2 = axs[1,0].imshow(Itot[:,:,0], extent=extent); fig.colorbar(im2, ax=axs[1,0]); plot_sites(posE,axs[1,0]); axs[1,0].set_title("I_occ");
    # im3 = axs[1,1].imshow(Itot[:,:,1], extent=extent); fig.colorbar(im3, ax=axs[1,1]); plot_sites(posE,axs[1,1]); axs[1,1].set_title("I_unocc");
    plot2d(Emin.T,        extent=extent, title="Emin",    ax=axs[0,0], cmap="viridis", ps=posE);
    plot2d(iMin.T,        extent=extent, title="iMin",    ax=axs[0,1], cmap="viridis", ps=posE);
    plot2d(Itot[:,:,0].T, extent=extent, title="I_occ",   ax=axs[1,0], cmap="viridis", ps=posE);
    plot2d(Itot[:,:,1].T, extent=extent, title="I_unocc", ax=axs[1,1], cmap="viridis", ps=posE);
    plt.tight_layout()
    plt.show()

def test_local_update_solver():
    # --- Example Usage for the new Local Update Solver ---
    
    # 1. Setup a larger system
    posE  = make_grid_sites(nxy=(8, 8), avec=(5.0, 0.0), bvec=(0.0, 5.0), E0=-0.15)
    nSite = posE.shape[0]
    print(f"Created a large system with {nSite} sites.")

    solver = HubbardSolver()

    # 2. Pre-calculate the dense interaction matrix
    _, Wij_dense = solver.eval_coupling_matrix(posE, zMirror=-0.5)

    # 3. Convert the dense matrix to the required sparse format
    nMaxNeigh = 16
    W_sparse = make_sparse_W(Wij_dense, nMaxNeigh=nMaxNeigh)
    print(f"Converted dense Wij to sparse format with nMaxNeigh={nMaxNeigh}.")

    # 4. Pre-calculate Esite and Thop for a few tip positions
    pTips = generate_xy_scan(extent=[-5, 5, -5, 5], nxy=(3, 1), zTip=3.0, Vbias=0.2)
    nTips = pTips.shape[0]

    # This part would typically involve a more complex physical model,
    # but for demonstration, we'll construct it manually.
    Esite_Thop = np.zeros((nTips, nSite, 2), dtype=np.float32)
    for i in range(nTips):
        tip_pos = pTips[i, :3]
        tip_vbias = pTips[i, 3]
        for j in range(nSite):
            dist = np.linalg.norm(tip_pos - posE[j, :3])
            # Simplified model for Esite and Thop
            Esite_Thop[i, j, 0] = posE[j, 3] + tip_vbias * (3.0 / dist) # On-site energy
            Esite_Thop[i, j, 1] = np.exp(-0.3 * dist)                 # Hopping
    
    # 5. Run the local update solver
    print("Running the parallel local update solver...")
    solver_params = {
        "kT": 0.001,
        "nIter": 20000,
        "solverMode": 2, # Annealing
        "seed": 42
    }
    E_out, state_out, Itot_out = solver.solve_local_updates(Esite_Thop, W_sparse, params=solver_params)
    
    # 6. Print results
    for i in range(nTips):
        print(f"--- Tip {i} ---")
        print(f"  Final Energy: {E_out[i]:.6f}")
        print(f"  Final State Mask: {state_out[i]:0{nSite}b}") # Print as binary string
        print(f"  Final Currents: {Itot_out[i]}")

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    #test_brute_force_solver()
    test_local_update_solver()

    

