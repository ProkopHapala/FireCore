# run it like this:
#   python -u -m pyBall.OCL.HubbardSolver | tee OUT

import os
import numpy as np
import pyopencl as cl
import matplotlib.pyplot as plt

from .OpenCLBase import OpenCLBase

COULOMB_CONST = 14.399644730092272   # eV*Angstrom
kBoltzmann    = 8.61733326214511e-5  # eV/K

# unlimited line length
#np.set_printoptions(threshold=np.inf)
np.set_printoptions(linewidth=np.inf)
#np.set_printoptions(suppress=True)


def screened_coulomb(r, screening_length=10.0):
    if r < 1e-6: return 0.0 # Avoid division by zero
    return COULOMB_CONST * np.exp(-r / screening_length) / r

default_params = {
    "tipDecay"  : 0.3,
    "tipRadius" : 3.0,
    "zTip"      : 5.0,
    "zMirror"   : -0.5,
    "Wcouple"   : 1.0,
}


# -----------------------------------------------------------------------------
#                  class  HubbardSolver 
# -----------------------------------------------------------------------------

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

    #occ_bytes = ((nSingle + 7) // 8)  # Ceiling division for number of bytes needed
    occ_bytes  = 32
    max_neighs = 8

    # ---- names used in buffer_dict / as object attributes
    # buf_pose  = "posE_buff"      #   nSingle  * sizeof(float4)
    # buf_pTips = "pTips_buff"     #   nTips    * sizeof(float4)
    # buf_Emin  = "Emin_buff"      #   nTips    * sizeof(float )
    # buf_iMin  = "iMin_buff"      #   nTips    * sizeof(int  )
    # buf_Itot  = "Itot_buff"      #   nTips    * sizeof(float2)

    def __init__(self, nloc: int = 32, device_index: int = 0):
        super().__init__(nloc=nloc, device_index=device_index)

        # Define numpy dtype matching OpenCL struct `cl_Mat3 { float4 a; float4 b; float4 c; };`
        self.cl_mat3_dtype = np.dtype([
            ('a', np.float32, 4),
            ('b', np.float32, 4),
            ('c', np.float32, 4),
        ])

        # Build the OpenCL program ------------------------------------------------
        base_path = os.path.dirname(os.path.abspath(__file__))
        rel_path  = "cl/hubbard.cl"
        if not self.load_program(rel_path=rel_path, base_path=base_path):
            raise RuntimeError("[HubbardSolver] Failed to load/compile hubbard.cl")


    # ---------------------------------------------------------------------
    #  HubbardSolver  :   Buffer Management 
    # ---------------------------------------------------------------------

    def realloc_buffers(self, nSingle: int, nTips: int):
        """(Re-)allocate GPU buffers as needed for the given problem size."""
        print("HubbardSolver::realloc_buffers() dir(self):")
        sz_f  = 4
        buffs = {
            #  name     handle, size 
            "posE":   sz_f*4 * nSingle,
            "pTips":  sz_f*4 * nTips,
            "Emin":   sz_f*1 * nTips,
            "iMin":   sz_f*1 * nTips,
            "Itot":   sz_f*2 * nTips,
        }
        self.try_make_buffers(buffs)

    def realloc_local_update_buffers(self, nSingle: int, nTips: int, nMaxNeigh: int):
        """(Re-)allocate GPU buffers specifically for the local update solver."""
        sz_f  = 4
        # Reuse output buffers from brute-force where possible (E_out, Itot_out)
        buffs = {
            "Esite":   (sz_f*1 * nTips   * nSingle   ),
            "Tsite":   (sz_f*1 * nTips   * nSingle   ),
            "Wval":    (sz_f*1 * nSingle * self.max_neighs ),
            "Widx":    (sz_f*1 * nSingle * self.max_neighs ),
            "WnNeigh": (sz_f*1 * nSingle  ),
            # Re-purpose brute-force outputs
            "Emin":    (sz_f*1 * nTips  ),
            "occ":     (self.occ_bytes * nTips  ), 
            "Itot":    (sz_f*2 * nTips  ),
        }
        self.try_make_buffers(buffs)

    def realloc_precalc_buffers(self, nSingle: int, nTips: int, nMulti: int):
        sz_f  = 4
        """(Re-)allocate GPU buffers needed for the precalc_esite_thop kernel."""
        # Site positions and initial energies
        # self.try_make_buff("posE_buff",   nSingle * sz_f * 4)                 # Tip positions and bias
        # self.try_make_buff("pTips_buff",  nTips * sz_f * 4)                  
        # self.try_make_buff("multipoleCoefs_buff", nSingle * nMulti * sz_f * 4)  # Multipole coefficients per site
        # self.try_make_buff("rotSite_buff", nSingle * self.cl_mat3_dtype.itemsize)  # Rotation matrices per site (cl_Mat3 dtype)
        # self.try_make_buff("Esite_Thop_buff", nTips * nSingle * 2 * 4)   # Output buffer for site energies and hopping amplitudes
        buffs = {
            "posE":      (sz_f*4 * nSingle),
            "pTips":     (sz_f*4 * nTips),
            "mpolCs":    (sz_f*4 * nSingle * nMulti),
            "rotSite":   (self.cl_mat3_dtype.itemsize * nSingle),
            "Esite":     (sz_f*1 * nTips * nSingle),
            "Tsite":     (sz_f*1 * nTips * nSingle),
        }
        self.try_make_buffers(buffs)

    # ---------------------------------------------------------------------
    #  HubbardSolver  :   Setup Kernels
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

    def setup_solve_local_updates(self, nSite, nTips, params={}, initMode=3, T=3.0, nIter=100, solverMode=0 ):
        """Prepare arguments for the solve_local_updates kernel."""
        kernel = self.prg.solve_local_updates
        
        solver_params = np.array([
            params.get("kT",         kBoltzmann*T ),
            params.get("nIter",      nIter),
            params.get("solverMode", solverMode),
            params.get("seed",       np.random.randint(0, 2**31))
        ], dtype=np.float32)

        args = (
            np.int32(nSite),
            np.int32(nTips),
            solver_params,
            self.Esite_buff,
            self.Tsite_buff,
            self.Wval_buff,
            self.Widx_buff,
            self.WnNeigh_buff,
            self.occ_buff,  
            self.Emin_buff, 
            self.Itot_buff, 
            np.int32(self.max_neighs),
            np.int32(initMode),
        )
        return args, kernel

    # MODIFIED: Updated setup method for the pre-calculation kernel
    def setup_precalc_esite_thop(self, nSingle, nTips, nMulti, params, bMirror, bRamp):
        """Prepare arguments for the precalc_Esite_Thop kernel."""
        kernel = self.prg.precalc_Esite_Thop        
        args = (
            np.int32(nSingle),
            np.int32(nTips),
            self.posE_buff,
            self.rotSite_buff,
            self.pTips_buff,
            np.float32(params.get("Rtip", 3.0)),
            np.array([params.get("zMirror", -0.5), params.get("zRampOffset", 0.0)], dtype=np.float32),
            self.mpolCs_buff, # MODIFIED: Pass buffer instead of value
            np.int32(nMulti),         # MODIFIED: Pass nMulti
            np.int32(bMirror),
            np.int32(bRamp),
            np.float32(params.get("Thop_decay", 0.3)),
            np.float32(params.get("Thop_amp", 1.0)),
            self.Esite_buff,
            self.Tsite_buff,
        )
        return args, kernel

    # ---------------------------------------------------------------------
    #  HubbardSolver  :   Run Kernels
    # ---------------------------------------------------------------------


    def precalc_esite_thop(self, posE, pTips, rots=None, multipoleCoefs=None, bMirror=True, bRamp=True, params=default_params):
        """
        Runs the pre-calculation kernel to generate on-site energies and hoppings.

        Args:
            posE (np.ndarray): Shape (nSingle, 4), site positions {x,y,z,E0}.
            pTips (np.ndarray): Shape (nTips, 4), tip positions {x,y,z,Vbias}.
            rots (np.ndarray, optional): Shape (nSingle, 3, 3), rotation matrices.
            multipoleCoefs (np.ndarray, optional): Shape (nSingle, nMulti), site-dependent coefficients.
            bMirror (bool): Enable/disable mirror charge effect.
            bRamp (bool): Enable/disable linear potential ramp.
            params (dict): Dictionary with physical parameters.

        Returns:
            np.ndarray: The calculated Esite_Thop array of shape (nTips, nSingle, 2).
        """
        nSingle = posE.shape[0]
        nTips = pTips.shape[0]

        # Handle multipole coefficients
        if multipoleCoefs is None:
            # Default to a simple monopole (c0=1.0) for all sites
            nMulti = 1
            multipoleCoefs = np.zeros((nSingle, nMulti), dtype=np.float32)
            multipoleCoefs[:, 0] = 1.0
        else:
            if multipoleCoefs.shape[0] != nSingle:
                raise ValueError("multipoleCoefs must have shape (nSingle, nMulti)")
            nMulti = multipoleCoefs.shape[1]

        # Allocate buffers
        self.realloc_precalc_buffers(nSingle, nTips, nMulti)
        
        # Upload data
        self.toGPU_(self.posE_buff, posE)
        self.toGPU_(self.pTips_buff, pTips)
        self.toGPU_(self.mpolCs_buff, multipoleCoefs)
        if rots is not None:
            rot_struct_arr = np.empty(nSingle, dtype=self.cl_mat3_dtype)
            for i in range(nSingle):
                rot_struct_arr[i]['a'][:3] = rots[i, 0, :]
                rot_struct_arr[i]['b'][:3] = rots[i, 1, :]
                rot_struct_arr[i]['c'][:3] = rots[i, 2, :]
            self.toGPU_(self.rotSite_buff, rot_struct_arr)

        # Kernel launch
        global_size = (nTips * nSingle,)
        args, kernel = self.setup_precalc_esite_thop(nSingle, nTips, nMulti, params, bMirror, bRamp)
        kernel(self.queue, global_size, None, *args)
        
        # Download results
        Esite = np.empty((nTips * nSingle), dtype=np.float32)
        Tsite = np.empty((nTips * nSingle), dtype=np.float32)
        self.fromGPU_(self.Esite_buff, Esite)
        self.fromGPU_(self.Tsite_buff, Tsite)

        self.queue.finish()

        return Esite.reshape((nTips, nSingle)), Tsite.reshape((nTips, nSingle))


    # NEW: High-level public API for the local update solver
    def solve_local_updates(self, W_sparse=None, Esite=None, Tsite=None, nTips=None, nSite=None, nMaxNeigh=None,  params=default_params, nWorkGroup=16, bRealloc=True ):
        """
        Run the local-update Monte Carlo solver.

        Args:
            Esite (np.ndarray): Shape (nTips, nSingle, 2), float32. Contains {Esite, Thop} for each tip-site pair.
            W_sparse (tuple): A tuple (W_val, W_idx, nNeigh) representing the sparse interaction matrix, typically from `make_sparse_W`.
            params (dict): Dictionary with solver parameters like kT, nIter, etc.

        Returns:
            E_out (np.ndarray): (nTips,) float32 - final energy for each simulation.
            state_out (np.ndarray): (nTips,) uint64 - final state bitmask.
            Itot_out (np.ndarray): (nTips, 2) float32 - final currents {I_occ, I_unoc}.
        """

        if nSite > 64: raise ValueError("[HubbardSolver] nSite > 64 exceeds ulong bitmask limit.")

        # Reallocate buffers and upload data
        if bRealloc:
            if nSite     is None: nSite = Esite.shape[0]
            if nTips     is None: nTips = Esite.shape[1]
            if nMaxNeigh is None: nMaxNeigh = len(W_sparse[2])
            self.realloc_local_update_buffers(nSite, nTips, nMaxNeigh)
        if Esite is not None: self.toGPU_(self.Esite_buff, Esite)
        if Tsite is not None: self.toGPU_(self.Tsite_buff, Tsite)
        if W_sparse is not None:
            Wval, Widx, WnNeigh = W_sparse
            self.toGPU_(self.Wval_buff,  Wval)
            self.toGPU_(self.Widx_buff,  Widx)
            self.toGPU_(self.WnNeigh_buff, WnNeigh)

        # Kernel launch
        global_size  = (nTips * nWorkGroup,)
        local_size   = (nWorkGroup,)
        args, kernel = self.setup_solve_local_updates(nSite, nTips, params)
        kernel(self.queue, global_size, local_size, *args)
        
        # Download results
        occ_out   = np.empty(nTips*self.occ_bytes, dtype=np.uint8)
        E_out     = np.empty(nTips, dtype=np.float32)
        Itot_out  = np.empty((nTips, 2), dtype=np.float32)
        
        self.fromGPU_(self.occ_buff,  occ_out)  # Re-purposed
        self.fromGPU_(self.Emin_buff, E_out)     # Re-purposed
        self.fromGPU_(self.Itot_buff, Itot_out)  # Re-purposed

        self.queue.finish()
        return E_out, Itot_out, occ_out
    


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
#                             Utility helpers 
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


# NEW: A more physically-motivated way to build the sparse interaction matrix.
def make_sparse_W_pbc(pos, lvecs, Rcut, W_func, nMaxNeigh=16):
    """
    Creates a sparse interaction matrix based on a distance cutoff with 2D
    periodic boundary conditions (PBC).

    For each site `i`, it finds all other sites `j` (including their periodic
    images) within a given `cutoff` distance. If more than `nMaxNeigh` neighbors
    are found, it keeps the ones with the strongest interaction strength.

    Args:
        pos (np.ndarray): Shape (nSite, 3), float. The XYZ positions of the sites.
        lvecs (np.ndarray): Shape (2, 3), float. The two lattice vectors defining
                          the 2D periodic cell (e.g., lvecs[0] is 'a', lvecs[1] is 'b').
        cutoff (float): The maximum distance to consider for interactions.
        W_func (callable): A function that takes a distance `r` (float) and returns
                         the interaction energy `W` (float).
        nMaxNeigh (int): The maximum number of neighbors to store per site.

    Returns:
        tuple: (W_val, W_idx, nNeigh) ready for GPU upload.
            - W_val (np.ndarray): (nSite, nMaxNeigh) float32, interaction values.
            - W_idx (np.ndarray): (nSite, nMaxNeigh) int32, neighbor indices.
            - nNeigh (np.ndarray): (nSite,) int32, number of actual neighbors.
    """
    nSite     = pos.shape[0]
    R2cut     = Rcut * Rcut
    
    # Store found neighbors temporarily in a list of lists.
    # Each entry will be a tuple: (interaction_strength, neighbor_index)
    #all_neighbors = [[] for _ in range(nSite)]

    iNeighs  = [[] for _ in range(nSite)]
    ngVals   = [[] for _ in range(nSite)]
    for i in range(nSite):
        pi = pos[i]
        for j in range(nSite):
            if i == j: continue
            pj_base = pos[j]
            
            # Check the 3x3 grid of periodic cells (including the home cell at [0,0])
            for nx in [-1, 0, 1]:
                for ny in [-1, 0, 1]:                        
                    # Calculate the position of the image of site j
                    shift = nx * lvecs[0, :] + ny * lvecs[1, :]
                    pj = pj_base + shift
                    
                    # Check if the distance is within the cutoff
                    d_vec = pi - pj
                    r2 = np.dot(d_vec, d_vec)
                    
                    if r2 < R2cut:
                        r = np.sqrt(r2)
                        #print(f"Site {i} {j} {r:.3f} < {Rcut:.3f}")
                        wij  = W_func(r)
                        nNeigh = len(iNeighs[i])
                        if nNeigh >= nMaxNeigh: raise ValueError(f"Too many neighbors for site {i} {nNeigh} >= {nMaxNeigh}")
                        iNeighs[i].append(j)
                        ngVals [i].append(wij)

    # Initialize final output arrays
    W_val  = np.zeros((nSite, nMaxNeigh), dtype=np.float32)
    W_idx  = np.zeros((nSite, nMaxNeigh), dtype=np.int32); W_idx[:,:] = -1
    nNeigh = np.zeros(nSite, dtype=np.int32)

    # Process the found neighbors for each site
    for i in range(nSite):
        nng = len(iNeighs[i])
        nNeigh[i]     = nng
        W_val[i,:nng] = ngVals [i]
        W_idx[i,:nng] = iNeighs[i]
        # neighbors_i   = all_neighbors[i]
        # neighbors_i.sort(key=lambda item: abs(item[0]), reverse=True)
        # num_to_store  = min(len(neighbors_i), nMaxNeigh)
        # nNeigh[i]    = num_to_store
        # for k in range(num_to_store):
        #     wij, j_idx  = neighbors_i[k]
        #     W_val[i, k] = wij
        #     W_idx[i, k] = j_idx
            
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
#                     High-level convenience functions 
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

# -----------------------------------------------------------------------------
#                            Tests  and Demos 
# -----------------------------------------------------------------------------

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
    Esite = np.zeros((nTips, nSite), dtype=np.float32)
    Tsite = np.zeros((nTips, nSite), dtype=np.float32)
    for i in range(nTips):
        tip_pos   = pTips[i, :3]
        tip_vbias = pTips[i, 3]
        for j in range(nSite):
            dist = np.linalg.norm(tip_pos - posE[j, :3])
            # Simplified model for Esite and Thop
            Esite[i, j] = posE[j, 3] + tip_vbias * (3.0 / dist) # On-site energy
            Tsite[i, j] = np.exp(-0.3 * dist)                 # Hopping
    
    # 5. Run the local update solver
    print("Running the parallel local update solver...")
    # solver_params = {
    #     "kT": 0.001,
    #     "nIter": 20000,
    #     "solverMode": 2, # Annealing
    #     "seed": 42
    # }
    E_out, state_out, Itot_out = solver.solve_local_updates(Esite, Tsite, W_sparse )
    
    # 6. Print results
    for i in range(nTips):
        print(f"--- Tip {i} ---")
        print(f"  Final Energy: {E_out[i]:.6f}")
        print(f"  Final State Mask: {state_out[i]:0{nSite}b}") # Print as binary string
        print(f"  Final Currents: {Itot_out[i]}")

def make_site_lattice( lvecs, nxy=(6,6), p0=(0.0,0.0,0.0) ):
    nx, ny = nxy
    nSite = nx * ny
    pos = np.zeros((nSite, 3), dtype=np.float32)
    idx = 0
    for iy in range(ny):
        for ix in range(nx):
            pos[idx, 0] = ix * lvecs[0,0] + p0[0]
            pos[idx, 1] = iy * lvecs[1,1] + p0[1]
            pos[idx, 2] = p0[2]
            idx += 1
    return pos

def test_site_coupling( nxy=(6,6), cutoff=8.0, a=5.0, b=5.0, nMaxNeigh=8 ):
    # Lattice vectors define the size of the periodic cell
    lvs = np.array([ [a, 0.0, 0.0], [0.0, b, 0.0] ], dtype=np.float32)
    pos = make_site_lattice(lvs, nxy=nxy, p0=(0.0,0.0,0.0))
    lvs_pbc = lvs.copy(); lvs_pbc[1,:] *= nxy[1]; lvs_pbc[0,:] *= nxy[0] 
    W_val, W_idx, nNeigh = make_sparse_W_pbc( pos, lvs_pbc, cutoff, W_func=screened_coulomb, nMaxNeigh=nMaxNeigh )
    print("nNeigh:", nNeigh)
    for i in range(len(nNeigh)):
        print(f"Site {i}: nNeigh: {nNeigh[i]} ",end="")
        for j in range(nNeigh[i]):
            print(f"  {W_idx[i,j]}: {W_val[i,j]:.2f}", end=" ")
        print()

def demo_precalc_scan(solver: HubbardSolver=None, nxy_sites=(6, 6), nxy_scan=(100, 100), Vbias=0.1):
    """
    Demonstrates the usage of the precalc_esite_thop kernel by running a 2D scan.

    This function visualizes the fundamental interaction landscape between the tip and
    the individual sites before any many-body effects are considered. It plots maps of:
    - The minimum on-site energy induced by the tip.
    - The maximum on-site energy induced by the tip.
    - The minimum hopping amplitude (i.e., the site the tip is closest to).

    Args:
        solver (HubbardSolver): An instance of your solver class.
        nxy_sites (tuple): The (nx, ny) dimensions of the grid of sample sites.
        nxy_scan (tuple): The (nx, ny) resolution of the tip scan.
        Vbias (float): The bias voltage to apply to the tip.
    """
    if solver is None:
        solver = HubbardSolver()

    print("--- Running pre-calculation scan demo ---")

    # 1. Define the system of sites
    posE    = make_grid_sites(nxy=nxy_sites, avec=(5.0, 0.0), bvec=(0.0, 5.0), E0=-0.15)
    nSingle = posE.shape[0]
    print(f"Created a {nxy_sites[0]}x{nxy_sites[1]} grid of {nSingle} sites.")

    # 2. Define the 2D grid of tip positions for the scan
    extent = [-15.0, 15.0, -15.0, 15.0]
    pTips  = generate_xy_scan(extent, nxy=nxy_scan, zTip=3.0, Vbias=Vbias)
    nTips  = pTips.shape[0]
    print(f"Generated a {nxy_scan[0]}x{nxy_scan[1]} grid of {nTips} tip positions.")

    # 3. Define the physical model parameters (e.g., simple monopole)
    multipoleCoefs = np.zeros((nSingle, 1), dtype=np.float32)
    multipoleCoefs[:, 0] = 1.0  # All sites are simple monopoles
    
    params = {
        "Rtip": 3.0,
        "zMirror": -0.5,
        "Thop_decay": 0.5, # A higher decay for more localized hopping
        "Thop_amp": 1.0
    }

    # 4. Call the pre-calculation kernel
    print("Calling precalc_esite_thop kernel...")
    Esite,Tsite = solver.precalc_esite_thop(posE, pTips, multipoleCoefs=multipoleCoefs, params=params)
    # Output shape is (nTips, nSingle, 2)
    print(f"Esite min,max : {np.min(Esite), np.max(Esite)}")
    print(f"Tsite min,max : {np.min(Tsite), np.max(Tsite)}")
    
    # 5. Perform reduction in Python/NumPy to get the desired maps
    # For each tip position (axis 0), find the min/max over all sites (axis 1)
    min_E_map    = np.min(Esite, axis=1) # Min Esite
    max_E_map    = np.max(Esite, axis=1) # Max Esite
    max_Thop_map = np.max(Tsite, axis=1) # Max Thop    # Hopping is always positive, so min is fine. Use max for "strongest hopping".
    
    # 6. Reshape the 1D arrays back into 2D maps for plotting
    nx, ny = nxy_scan
    min_E_map    = min_E_map.reshape((nx, ny))
    max_E_map    = max_E_map.reshape((nx, ny))
    max_Thop_map = max_Thop_map.reshape((nx, ny))

    # 7. Plot the results using the provided helper functions
    fig, axs = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(f"Tip-Site Interaction Maps (Vbias = {Vbias:.2f} V)", fontsize=16)

    plot2d(min_E_map.T,    extent=extent, title="Minimum Site Energy (E_min)", ax=axs[0], cmap="viridis", ps=posE)
    plot2d(max_E_map.T,    extent=extent, title="Maximum Site Energy (E_max)", ax=axs[1], cmap="magma",   ps=posE)
    plot2d(max_Thop_map.T, extent=extent, title="Maximum Hopping (T_max)",     ax=axs[2], cmap="inferno", ps=posE)
           
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()

def plot_occupancy_slice(occupation, nxy_scan, occ_bytes, slice_idx, axis='x'):
    """
    Plots a 1D slice of the occupancy array as a binary image.

    Args:
        occupation (np.ndarray): The 1D array of site occupancies (nTips * occ_bytes).
        nxy_scan (tuple): The (nx, ny) dimensions of the tip scan grid.
        occ_bytes (int): The number of bytes per tip position for occupancy.
        slice_idx (int): The index for the slice (e.g., ix if axis='x', iy if axis='y').
        axis (str): The axis along which to take the slice ('x' or 'y').
    """
    nx, ny = nxy_scan
    nTips  = nx * ny
    nSites = occ_bytes * 8

    # Reshape occupation to (ny, nx, occ_bytes) for easier slicing
    # Note: The order here depends on how the scan grid is flattened to nTips
    # Assuming row-major (y then x) for indexing pTips
    occupancy_2d_bytes = occupation.reshape(ny, nx, occ_bytes)

    if axis == 'x':
        # Slice along y-axis at fixed x_idx
        if not (0 <= slice_idx < nx):
            print(f"Error: slice_idx {slice_idx} out of bounds for nx={nx}")
            return
        slice_data_bytes = occupancy_2d_bytes[:, slice_idx, :]
        title_suffix = f" (x_idx={slice_idx})"
    elif axis == 'y':
        # Slice along x-axis at fixed y_idx
        if not (0 <= slice_idx < ny):
            print(f"Error: slice_idx {slice_idx} out of bounds for ny={ny}")
            return
        slice_data_bytes = occupancy_2d_bytes[slice_idx, :, :]
        title_suffix = f" (y_idx={slice_idx})"
    else:
        print("Error: axis must be 'x' or 'y'")
        return

    # Convert bytes to bits
    # Each row in slice_data_bytes is (occ_bytes,)
    # We want to convert each byte into 8 bits and concatenate them horizontally
    binary_image_data = np.zeros((slice_data_bytes.shape[0], nSites), dtype=np.uint8)
    for i, byte_row in enumerate(slice_data_bytes):
        bits = np.unpackbits(byte_row)
        binary_image_data[i, :] = bits

    plt.figure(figsize=(10, 5))
    plt.imshow(binary_image_data, cmap='gray', aspect='auto', interpolation='nearest')
    plt.title(f"Site Occupancy Slice{title_suffix}")
    plt.xlabel("Site Index (0 to nSites-1)")
    plt.ylabel(f"Tip Position Index along {'Y' if axis=='x' else 'X'} axis")
    plt.colorbar(label="Occupancy (0=empty, 1=occupied)")
    plt.show()

def plot_occupancy_line(occupation, nSingle ):
    #binary_image_data = np.zeros((slice_data_bytes.shape[0], nSites), dtype=np.uint8)
    # for i, byte_row in enumerate(slice_data_bytes):
    #     bits = np.unpackbits(byte_row)
    #     binary_image_data[i, :] = bits
    nbyte = nSingle//8
    bits = np.unpackbits(occupation, axis=1)
    print( "plot_occupancy_line() bits.shape = ", bits.shape, nSingle, nbyte)
    bits = bits[:,:nSingle]
    plt.figure(figsize=(10, 5))
    plt.imshow(bits, cmap='gray', aspect='auto', interpolation='nearest')
    plt.title(f"Site Occupancy Slice")
    plt.xlabel("Site Index (0 to nSites-1)")
    plt.ylabel(f"Tip Position Index ")
    plt.colorbar(label="Occupancy (0=empty, 1=occupied)")
    #plt.show()

def print_occupancy(occupation, nSingle, bHex=False ):
    occ_bytes = nSingle//8
    nTips = occupation.shape[0]
    for i in range(nTips):
        #print(f"Tip {i}: E={energy[i]:.3f}, I_occ={current[i,0]:.3f}, I_unocc={current[i,1]:.3f}")
        #i0 = i*nByteMax
        print(f"CPU iTip {i:3} occ: ", end="")
        for j in range(occ_bytes):
            if bHex:
                print(f"{occupation[i,j]:02x}", end="")
            else:
                print(f"{occupation[i,j]:08b}", end="")
        print()

def demo_local_update(solver: HubbardSolver=None, nxy_sites=(4, 4), nxy_scan=(50, 50), Vbias=0.1, cutoff=8.0, W_amplitude=1.0, T=0.001, nIter=100):
    """
    Demonstrates the usage of the solve_local_updates kernel by running a 2D scan with Monte Carlo optimization.
    
    This function extends demo_precalc_scan by adding the many-body effects through Monte Carlo optimization.
    It calculates site energies and hoppings, then runs the local update solver for each tip position.
    The function plots maps of:
    - The minimum energy after Monte Carlo optimization
    - The maximum hopping amplitude after optimization
    - The total charge (number of occupied sites) for each tip position
    
    Args:
        solver (HubbardSolver): An instance of your solver class.
        nxy_sites (tuple): The (nx, ny) dimensions of the grid of sample sites.
        nxy_scan (tuple): The (nx, ny) resolution of the tip scan.
        Vbias (float): The bias voltage to apply to the tip.
        cutoff (float): Cutoff distance for site-site interactions.
        W_amplitude (float): Strength of the Coulomb interaction between sites.
        T (float): Temperature for the Monte Carlo simulation in Kelvin.
        nIter (int): Number of Monte Carlo iterations.
    """
    if solver is None:
        solver = HubbardSolver()

    print("--- Running local update Monte Carlo solver demo ---")

    # 1. Define the system of sites
    posE    = make_grid_sites(nxy=nxy_sites, avec=(5.0, 0.0), bvec=(0.0, 5.0), E0=-0.15)
    nSingle = posE.shape[0]
    print(f"Created a {nxy_sites[0]}x{nxy_sites[1]} grid of {nSingle} sites.")

    # 2. Define the 2D grid of tip positions for the scan
    extent = [-15.0, 15.0, -15.0, 15.0]
    pTips = generate_xy_scan(extent, nxy=nxy_scan, zTip=3.0, Vbias=Vbias)
    nTips = pTips.shape[0]
    print(f"Generated a {nxy_scan[0]}x{nxy_scan[1]} grid of {nTips} tip positions.")

    # 3. Define the physical model parameters
    # 3a. Multipole coefficients (simple monopole)
    multipoleCoefs = np.zeros((nSingle, 1), dtype=np.float32)
    multipoleCoefs[:, 0] = 1.0  # All sites are simple monopoles
    
    # 3b. Pre-calculation kernel parameters
    params_precalc = {
        "Rtip": 3.0,
        "zMirror": -0.5,
        "Thop_decay": 0.5, # A higher decay for more localized hopping
        "Thop_amp": 1.0
    }
    
    # 3c. Generate site-site interactions (W matrix)
    # Periodic boundary conditions for the lattice
    lvs = np.array([[5.0, 0.0, 0.0], [0.0, 5.0, 0.0]], dtype=np.float32)
    lvs_pbc = lvs.copy()
    lvs_pbc[0,:] *= nxy_sites[0]
    lvs_pbc[1,:] *= nxy_sites[1]
    
    # Function to calculate Coulomb interaction with screening
    def screened_coulomb(r):
        return W_amplitude * np.exp(-r/cutoff) / (r + 0.1)
    
    # Create sparse W matrix with periodic boundary conditions
    nMaxNeigh = 8  # Maximum number of neighbors per site
    # Extract only the position components (x,y,z) from posE for the PBC calculation
    pos_xyz = posE[:, :3]  # Only take x, y, z components, not the energy E0
    W_val, W_idx, nNeigh = make_sparse_W_pbc(pos_xyz, lvs_pbc, cutoff, W_func=screened_coulomb, nMaxNeigh=nMaxNeigh)
    print(f"Created sparse W matrix with cutoff {cutoff} and max {nMaxNeigh} neighbors per site.")
    print("GPU W_val: \n", W_val)
    print("GPU W_idx: \n", W_idx)
    
    # 3d. Local update solver parameters
    params_solver = {
        "kT":    kBoltzmann * T,  # Temperature in energy units
        "nIter": nIter,      # Number of Monte Carlo iterations
        "solverMode": 2,     # Use Metropolis algorithm (0=deterministic, 2=simulated annealing)
        "seed":  np.random.randint(0, 2**31)  # Random seed
    }

    # 4. Call the pre-calculation kernel
    print("Calling precalc_esite_thop kernel...")
    Esite, Tsite = solver.precalc_esite_thop(posE, pTips, multipoleCoefs=multipoleCoefs, params=params_precalc)
    print(f"Kernel finished. Sites: {nSingle}, Tips: {nTips}")
    print(f"Esite shape: {Esite.shape}, min/max: {np.min(Esite):.3f}/{np.max(Esite):.3f}")
    print(f"Tsite shape: {Tsite.shape}, min/max: {np.min(Tsite):.3f}/{np.max(Tsite):.3f}")

    # 5. Run the local update Monte Carlo solver
    print(f"Running local update solver with T={T}K, nIter={nIter}...")
    
    # Allocate buffers for the solver
    solver.realloc_local_update_buffers(nSingle, nTips, nMaxNeigh)
    
    energy, current, occupation = solver.solve_local_updates( W_sparse=(W_val, W_idx, nNeigh), nTips=nTips, nSite=nSingle, nMaxNeigh=nMaxNeigh )

    # occ_bytes = nSingle//8
    # for i in range(nTips):
    #     #print(f"Tip {i}: E={energy[i]:.3f}, I_occ={current[i,0]:.3f}, I_unocc={current[i,1]:.3f}")
    #     print(f"CPU iTip {i:3} occ: ", end="")
    #     for j in range(occ_bytes):
    #         #print(f"{occupation[i*solver.occ_bytes+j]:02x}", end="")
    #         print(f"{occupation[i*solver.occ_bytes+j]:08b}", end="")
    #     print()

    occ_slice = occupation.reshape( nxy_scan[0], nxy_scan[1], -1 )[:, nxy_scan[1]//2,:]

    #print_occupancy(occupation.reshape(nTips,-1), nSingle, bHex=False)

    print_occupancy(occ_slice, nSingle, bHex=False)


    # Plot occupancy slice
    # For demonstration, let's plot a slice along x-axis at y_idx = ny // 2
    # Or along y-axis at x_idx = nx // 2

    plot_occupancy_line( occ_slice, nSingle )
    #plot_occupancy_slice(occupation, nxy_scan, solver.occ_bytes, slice_idx=nxy_scan[0] // 2, axis='x')
    
    # Calculate total charge for each tip position (count bits set to 1)
    total_charge        = np.zeros(nTips, dtype=np.int32)
    occupation_reshaped = occupation.reshape(nTips, solver.occ_bytes)
    bits = np.unpackbits(occupation_reshaped, axis=1)   # shape (nTips, 8*occ_bytes)
    total_charge = bits.sum(axis=1)                     # shape (nTips,)
    
    # Reshape results into 2D maps
    nx, ny = nxy_scan
    energy_map       = energy.reshape((nx, ny))
    current_map_occ  = current[:, 0].reshape((nx, ny))   # Occupied sites current
    current_map_unoc = current[:, 1].reshape((nx, ny))  # Unoccupied sites current
    charge_map       = total_charge.reshape((nx, ny))
    
    print(f"Energy range: {np.min(energy):.3f} to {np.max(energy):.3f}")
    print(f"Charge range: {np.min(total_charge)} to {np.max(total_charge)} sites")
    print(f"Current (occupied) range: {np.min(current[:, 0]):.3f} to {np.max(current[:, 0]):.3f}")
    
    # 6. Plot the results
    fig, axs = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(f"Monte Carlo Optimization Results (T={T}K, nIter={nIter}, W={W_amplitude})", fontsize=16)
    
    plot2d(energy_map.T,       extent=extent, title="Energy (optimized)", ax=axs[0, 0], cmap="viridis", ps=posE)
    plot2d(charge_map.T,       extent=extent, title="Total Charge (# occupied sites)", ax=axs[0, 1], cmap="plasma", ps=posE)
    plot2d(current_map_occ.T,  extent=extent, title="Current (occupied sites)", ax=axs[1, 0], cmap="inferno", ps=posE)
    plot2d(current_map_unoc.T, extent=extent, title="Current (unoccupied sites)", ax=axs[1, 1], cmap="magma", ps=posE)
    
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    #test_brute_force_solver()
    #test_site_coupling()
    #demo_precalc_scan()
    demo_local_update()

    

