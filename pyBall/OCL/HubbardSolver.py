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

    def __init__(self, nloc:int= 32, nloc_MC=16, device_index: int=0 ):
        super().__init__(nloc=nloc, device_index=device_index)

        self.nloc_MC = nloc_MC

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
            "Wval":    (sz_f*1 * nSingle * nMaxNeigh ),
            "Widx":    (sz_f*1 * nSingle * nMaxNeigh ),
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

    # def realloc_mc_buffers(self, nSite: int, nTips: int, nMaxNeigh: int):
    #     """
    #     (Re-)allocate GPU buffers specifically for the solve_MC kernel.
    #     These buffers will hold the best-known state across multiple runs.
    #     """
    #     sz_f = np.nbytes[np.float32]
    #     buffs = {
    #         # --- Input data describing the system ---
    #         "Esite":       (sz_f * 1 * nTips * nSite),
    #         "Tsite":       (sz_f * 1 * nTips * nSite),
    #         "W_val":       (sz_f * 1 * nSite * nMaxNeigh),
    #         "W_idx":       (np.nbytes[np.int32] * 1 * nSite * nMaxNeigh),
    #         "nNeigh":      (np.nbytes[np.int32] * 1 * nSite),
    #         # --- Input/Output buffers for storing the best result ---
    #         "E_best":      (sz_f * 1 * nTips),
    #         "occ_best":    (self.occ_bytes * nTips),
    #         "Itot_out":    (sz_f * 2 * nTips),
    #         "rng_seeds":   (sz_f * nTips * self.nloc_MC ),
    #     }
    #     self.try_make_buffers(buffs)

    def realloc_mc_buffers(self, nSite: int, nTips: int, nMaxNeigh: int):
        """
        (Re-)allocate GPU buffers for the 2-phase MC solver.
        """
        sz_f  = np.nbytes[np.float32]
        sz_i  = np.nbytes[np.int32]
        sz_ui = np.nbytes[np.uint32]

        buffs = {
            # --- System Data ---
            "Esite":       (sz_f * nTips * nSite),
            "Tsite":       (sz_f * nTips * nSite), # Still needed for final calculation
            "W_val":       (sz_f * nSite * nMaxNeigh),
            "W_idx":       (sz_i * nSite * nMaxNeigh),
            "nNeigh":      (sz_i * nSite),
            "rng_seeds":   (sz_ui * nTips * self.nloc_MC ),
            
            # --- Ping-Pong State Buffers ---
            "E_best_A":    (sz_f * nTips), "occ_best_A": (self.occ_bytes * nTips),
            "E_best_B":    (sz_f * nTips), "occ_best_B": (self.occ_bytes * nTips),

            # --- Final Output ---
            "Itot_out":    (sz_f * 2 * nTips),
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

    def setup_solve_local_updates(self, nSite, nTips, params={}, initMode=3, solverMode=0, T=3.0, nIter=100, bNoCoupling=False, max_neighs=None ):
        """Prepare arguments for the solve_local_updates kernel."""
        kernel = self.prg.solve_local_updates
        
        solver_params = np.array([
            params.get("kT",         kBoltzmann*T ),
            params.get("nIter",      nIter),
            params.get("solverMode", solverMode),
            params.get("seed",       np.random.randint(0, 2**31))
        ], dtype=np.float32)

        if max_neighs is None: max_neighs = self.max_neighs
        if bNoCoupling:        max_neighs = 0
        
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
            np.int32(max_neighs),
            np.int32(initMode),

        )
        return args, kernel

    def setup_solve_mc_neigh(self, nSite, nTips, params={}, nLocalIter=10, prob_params=(0.05, 0.45, 0.50, 0.0), nx=None, initMode=3, max_neighs=16):
        """Prepare arguments for the solve_MC kernel."""
        kernel = self.prg.solve_MC_neigh
        
        solver_params = np.array([
            params.get("kT", 300.0 * 8.617333262e-5), # kT is unused, but kept for signature
            params.get("nIter", 1000),
            params.get("solverMode", 1), 
            params.get("seed", np.random.randint(0, 2**31))
        ], dtype=np.float32)

        renorm = np.sum(prob_params)
        prob_params = np.cumsum(prob_params)/renorm
        prob_params_vec = np.array(prob_params, dtype=np.float32)
        print( "prob_params_vec ", prob_params_vec )

        args = (
            np.int32(nSite),
            np.int32(nTips),
            solver_params,
            np.int32(nLocalIter),
            prob_params_vec,                       # Pass the new probability vector
            np.int32(nx),                          # Pass the image width
            self.Esite_buff,
            self.Tsite_buff,
            self.W_val_buff,
            self.W_idx_buff,
            self.nNeigh_buff,
            self.occ_best_A_buff,
            self.E_best_A_buff,
            self.Itot_out_buff,
            self.rng_seeds_buff,
            np.int32(max_neighs),
            np.int32(initMode),
        )
        return args, kernel

    def setup_solve_mc_2phase(self, nSite, nTips, nLocalIter, nGlobalSteps, prob_params, nx, in_buffs, out_buffs ):
        """Pre-assembles arguments for the main optimization kernel."""
        kernel = self.prg.solve_MC_2phase
        prob_params_vec = np.array(np.cumsum(prob_params), dtype=np.float32)
        # // --- Optimization Parameters ---
        # const int nSite,                  // 1 : Number of sites
        # const int nTips,                  // 2 : Number of tips
        # const int nLocalIter,             // 3 : Number of local iterations
        # const float4 prob_params,         // 4 : cummulative mutation probability (load best, load neighbor, random reset)
        # const int nx,                     // 5 : Number of tips per row (for neighbor pixel selection)
        # // --- Ping-Pong State Buffers ---
        # __global const uchar* occ_in,     // 6 : [nTips,OCC_BYTES] Input occupancy buffer
        # __global const float* E_in,       // 7 : [nTips] Input energy buffer
        # __global uchar*       occ_out,    // 8 : [nTips,OCC_BYTES] Output occupancy buffer
        # __global float*       E_out,      // 9 : [nTips] Output energy buffer
        # // --- System Data (Energy-Related Only) ---
        # __global const float* Esite,      // 10 : [nTips,nSite] Site energies
        # const int max_neighs,             // 11 : Maximum number of neighbors per site in couling matrix Wij
        # __global const int*   nNeigh,     // 12 : [nSite] Number of neighbors for each site in couling matrix Wij
        # __global const float* W_val,      // 13 : [nSite,max_neighs] Coulomb interaction strengths     (Wij), sparse matrix, only max_neighs
        # __global const int*   W_idx,      // 14 : [nSite,max_neighs] indexes of neighbor sites  j for  (Wij), sparse matrix
        # __global const uint*  rng_seeds,  // 15 : [nTips] Random number generator seeds for each tip position (image pixel)
        # const uint glob_seed              // 16 : Global random seed ( to xor with rng_seeds)
        # Append buffer handles
        kernel_args = (
            np.int32(nSite),          # 1
            np.int32(nTips),          # 2
            np.int32(nLocalIter),     # 3
            prob_params_vec,          # 4
            np.int32(nx),             # 5
            #---
            in_buffs[0],              # 6 occ_in
            in_buffs[1],              # 7 E_in
            out_buffs[0],             # 8 occ_out
            out_buffs[1],             # 9 E_out
            # ---
            self.Esite_buff,          # 10 Esite
            np.int32(self.max_neighs),# 11 max_neighs
            self.nNeigh_buff,         # 12 nNeigh
            self.W_val_buff,          # 13 W_val
            self.W_idx_buff,          # 14 W_idx
            self.rng_seeds_buff,      # 15 rng_seeds
            # np.int32(gseed),        # 16 --- missing, add when run kernel 
        )
        return kernel_args, kernel

    def setup_calculate_currents(self, nSite, nTips, final_occ_buff):
        """Pre-assembles arguments for the final current calculation kernel."""
        # __kernel void calculate_currents(
        #     const int nSite,                 // 1 : Number of sites
        #     const int nTips,                 // 2 : Number of tips
        #     __global const uchar*  occ,      // 3 : [nTips,OCC_BYTES] The final optimized configurations
        #     __global const float*  Tsite,    // 4 : [nTips,nSite] The tunneling data
        #     __global float2*       Itot_out  // 5 : [nTips,2] The final output buffer for currents
        # ) {
        kernel = self.prg.calculate_currents
        kernel_args = (
            np.int32(nSite),    # 1
            np.int32(nTips),    # 2
            final_occ_buff,     # 3
            self.Tsite_buff,    # 4
            self.Itot_out_buff  # 5
        )
        return kernel_args, kernel

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
    def solve_local_updates(self, W_sparse=None, Esite=None, Tsite=None, nTips=None, nSite=None, nMaxNeigh=None, params=default_params, bRealloc=True, initMode=3, bNoCoupling=False ):
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
            if Esite is not None:
                if nTips is None: nTips = Esite.shape[0]
                if nSite is None: nSite = Esite.shape[1]
            if nMaxNeigh is None:
                if W_sparse is not None: nMaxNeigh = W_sparse[0].shape[1]
                else:                    nMaxNeigh = 0
            self.realloc_local_update_buffers(nSite, nTips, nMaxNeigh)
        if Esite is not None: 
            iDBG = 612
            #print("solve_local_updates() Esite.toGPU_ shape() ", Esite.shape )
            #print(f"solve_local_updates() Esite[{iDBG}] \n", Esite[iDBG] )
            #print("solve_local_updates() Esite.toGPU_ \n", Esite )
            self.toGPU_(self.Esite_buff, Esite)
        if Tsite is not None: self.toGPU_(self.Tsite_buff, Tsite)
        if W_sparse is not None:
            Wval, Widx, WnNeigh = W_sparse
            self.toGPU_(self.Wval_buff,  Wval)
            self.toGPU_(self.Widx_buff,  Widx)
            self.toGPU_(self.WnNeigh_buff, WnNeigh)

        # Kernel launch
        global_size  = (nTips * self.nloc_MC,)
        local_size   = (self.nloc_MC,)
        args, kernel = self.setup_solve_local_updates(nSite, nTips, params, initMode=initMode, bNoCoupling=bNoCoupling, max_neighs=nMaxNeigh)
        kernel(self.queue, global_size, local_size, *args)
        
        # Download results
        occ_out   = np.empty((nTips,self.occ_bytes), dtype=np.uint8)
        E_out     = np.empty(nTips, dtype=np.float32)
        Itot_out  = np.empty((nTips, 2), dtype=np.float32)
        
        self.fromGPU_(self.occ_buff,  occ_out)  # Re-purposed
        self.fromGPU_(self.Emin_buff, E_out)     # Re-purposed
        self.fromGPU_(self.Itot_buff, Itot_out)  # Re-purposed

        self.queue.finish()
        return E_out, Itot_out, occ_out
    
    def solve_mc(self, W_sparse, Esite, Tsite, nTips, nSite,
                 params={}, 
                 nLocalIter=100,
                 prob_params=( 0.1, 0.0, 0.5, 0.0), #  MC_neigh:{LoadBest, LoadNeighbor, RandomReset,  unused}   MC:(keep, reload_best, random_reset, unused)
                 #prob_params=(0.0, 0.5, 0.5, 0.0),
                 nRelaxSteps=50,
                 bRealloc=True,
                 initMode=3,
                 nxy_scan=None
                 ):
        """
        Run the hybrid Monte Carlo solver with global exploration.

        This solver combines fast local updates with periodic global checks
        to avoid local minima and improve ergodicity. The best-found energy
        and configuration are stored on the GPU between calls if bRealloc=False
        and initMode=2.

        Args:
            W_sparse (tuple): (W_val, W_idx, nNeigh) sparse interaction matrix.
            Esite (np.ndarray): (nTips, nSite) array of on-site energies.
            Tsite (np.ndarray): (nTips, nSite) array of tunneling factors.
            nTips (int): Number of tip positions (simulations).
            nSite (int): Number of sites in the system.
            params (dict): Solver parameters (kT, nIter, solverMode, seed).
            nLocalIter (int): Number of local updates before a global check.
            randConfigProb (float): Probability of resetting to a random state.
            nRelaxSteps (int): Initial local relaxation steps.
            nWorkGroup (int): Number of threads per work-group.
            bRealloc (bool): If True, re-allocate all GPU buffers.
            initMode (int): 0=zeros, 1=ones, 2=continue from buffer, 3=random.

        Returns:
            E_out (np.ndarray): (nTips,) best energy for each simulation.
            Itot_out (np.ndarray): (nTips, 2) {I_occ, I_unoc} for the best state.
            occ_out (np.ndarray): (nTips, occ_bytes) bitmask of the best state.
        """
        if nSite > self.occ_bytes * 8:
            raise ValueError(f"nSite ({nSite}) exceeds max supported sites ({self.occ_bytes*8})")
            
        nMaxNeigh = W_sparse[0].shape[1] if W_sparse is not None else 0

        # --- 1. Allocate Buffers and Upload Data ---
        if bRealloc:
            self.realloc_mc_buffers(nSite, nTips, nMaxNeigh)
            # --- Generate and upload the high-quality seeds ---
            #if "seed" in params:
            #    rng = np.random.Generator(np.random.PCG64(params["seed"]))            
            rng = np.random.Generator(np.random.PCG64())
            seeds = rng.integers(0, 2**32, size=nTips * self.nloc_MC, dtype=np.uint32)
            print( f"seeds.shape {seeds.shape}  seeds[:8] ", seeds[:8] )
            # print("self.toGPU(self.rng_seeds_buff, seeds) ", seeds.shape)
            # print( "self.buffer_dict", self.buffer_dict )
            # print( "hasattr(self,`rng_seeds_buff`): ", hasattr(self,"rng_seeds_buff") ) 
            self.toGPU_(self.rng_seeds_buff, seeds)
            #print("self.toGPU(self.rng_seeds_buff, seeds) DONE " )

            E_best_init = np.full(nTips, np.inf, dtype=np.float32)
            self.toGPU_(self.E_best_A_buff, E_best_init)
            occ_init = np.zeros((nTips, self.occ_bytes), dtype=np.uint8)
            self.toGPU_(self.occ_best_A_buff, occ_init)

        self.toGPU_(self.Esite_buff, Esite)
        self.toGPU_(self.Tsite_buff, Tsite)
        if W_sparse is not None:
            W_val, W_idx, nNeigh = W_sparse
            self.toGPU_(self.W_val_buff, W_val)
            self.toGPU_(self.W_idx_buff, W_idx)
            self.toGPU_(self.nNeigh_buff, nNeigh)

        # --- 2. Prepare and Launch Kernel ---
        global_size  = (nTips * self.nloc_MC,)
        local_size   = (self.nloc_MC,)
        
        args, kernel = self.setup_solve_mc_neigh( nSite, nTips, params, nLocalIter=nLocalIter, prob_params=prob_params, nx=nxy_scan[0], initMode=initMode, max_neighs=nMaxNeigh )
        kernel(self.queue, global_size, local_size, *args)
        
        # --- 3. Download Results ---
        occ_out  = np.empty((nTips, self.occ_bytes), dtype=np.uint8)
        E_out    = np.empty(nTips,                   dtype=np.float32)
        Itot_out = np.empty((nTips, 2),              dtype=np.float32)
        
        self.fromGPU_(self.occ_best_A_buff, occ_out  )
        self.fromGPU_(self.E_best_A_buff,   E_out    )
        self.fromGPU_(self.Itot_out_buff, Itot_out )

        self.queue.finish()
        
        return E_out, Itot_out, occ_out

    def solve_mc_2phase(self, 
            W_sparse, 
            Esite, 
            Tsite, 
            nTips, 
            nSite, 
            nx,
            nGlobalSteps=100,
            nLocalIter=100,
            prob_params=(0.1, 0.6, 0.3), # (p_best, p_neighbor, p_random)
            bAlloc=True,
            bFinalize=True,
            printsPerStep=10
        ):
        """
        Drives the 2-phase MC optimization with a fast Python loop.
        Call this once per phase.
        """
        nMaxNeigh = W_sparse[0].shape[1] if W_sparse is not None else 0
        #self.nloc_MC = 16 # Fixed workgroup size

        # --- Phase 1: Allocation and Initial Data Upload ---
        if bAlloc:
            print("solve_mc_2phase().alloc")
            self.realloc_mc_buffers(nSite, nTips, nMaxNeigh)
            
            E_best_init = np.full(nTips, np.inf, dtype=np.float32)
            self.toGPU_(self.E_best_A_buff, E_best_init)
            occ_init = np.zeros((nTips, self.occ_bytes), dtype=np.uint8)
            self.toGPU_(self.occ_best_A_buff, occ_init)

            self.toGPU_(self.Esite_buff, Esite)
            self.toGPU_(self.Tsite_buff, Tsite)
            if W_sparse:
                self.toGPU_(self.W_val_buff, W_sparse[0])
                self.toGPU_(self.W_idx_buff, W_sparse[1])
                self.toGPU_(self.nNeigh_buff, W_sparse[2])
            #rng = np.random.Generator(np.random.PCG64(params.get("seed", 12345)))
            #seeds = rng.integers(0, 2**32, size=nTips * self.nloc_MC, dtype=np.uint32)
            # defualt python random numbers
            seeds = np.random.randint(0, 2**32, size=nTips * self.nloc_MC, dtype=np.uint32)
            self.toGPU_(self.rng_seeds_buff, seeds)

        # --- Pre-computation for the Fast Loop ---
        print(f"solve_mc_2phase().setup_kernels")
        
        # setup_solve_mc_2phase(self, nSite, nTips, nLocalIter, nGlobalSteps, prob_params, nx, in_buffs, out_buffs ):
        # Pre-assemble both sets of arguments (A->B and B->A)
        args1, kernel = self.setup_solve_mc_2phase( nSite, nTips, nLocalIter, nGlobalSteps, prob_params, nx, (self.occ_best_A_buff, self.E_best_A_buff), (self.occ_best_B_buff, self.E_best_B_buff) )
        args2, _      = self.setup_solve_mc_2phase( nSite, nTips, nLocalIter, nGlobalSteps, prob_params, nx, (self.occ_best_B_buff, self.E_best_B_buff), (self.occ_best_A_buff, self.E_best_A_buff) )

        print(f"args from setup_solve_mc_2phase() {len(args1)}:"); 
        for i,arg in enumerate( args1 ): print(f"arg #{i}: {arg}")
        
        loo_seeds = np.random.randint(0, 2**32, size=nGlobalSteps, dtype=np.uint32)

        global_size = (nTips * self.nloc_MC,)
        local_size  = (self.nloc_MC,)
        
        # --- The Fast Python Loop ---
        print("solve_mc_2phase().loop nGlobalSteps", nGlobalSteps)
        for i_global in range(nGlobalSteps):
            gseed = loo_seeds[i_global]                                  # Update the global_seed for this iteration. This is a very fast operation.
            kernel(self.queue, global_size, local_size, *args1, gseed )  # Launch the kernel with the current set of arguments
            args1, args2 = args2, args1
            #print(f"solve_mc_2phase().loop step # {i_global}/{nGlobalSteps}")
            if (i_global + 1) % printsPerStep == 0: 
                self.queue.finish()
                print(f"    ...step {i_global + 1}/{nGlobalSteps} complete.")
        
        # --- Finalization and Data Download ---
        if bFinalize:
            print("solve_mc_2phase().finalize")
            bEven = (nGlobalSteps%2==0)
            if bEven:
                final_occ_buff = self.occ_best_A_buff
                final_E_buff   = self.E_best_A_buff
            else:
                final_occ_buff = self.occ_best_B_buff
                final_E_buff   = self.E_best_B_buff
            # Launch the current calculation kernel
            current_args, current_kernel = self.setup_calculate_currents(nSite, nTips, final_occ_buff)
            print(f"args from setup_calculate_currents() {len(current_args)}:")
            for i,arg in enumerate( current_args ): print(f"arg #{i}: {arg}")
            current_kernel(self.queue, (nTips,), None, *current_args)
            self.queue.finish() # Ensure the very last kernel launch is complete
            # Download all results
            E_out    = np.empty(nTips, dtype=np.float32)
            occ_out  = np.empty((nTips, self.occ_bytes), dtype=np.uint8)
            Itot_out = np.empty((nTips, 2), dtype=np.float32)
            self.fromGPU_(final_E_buff      , E_out)
            self.fromGPU_(final_occ_buff    , occ_out)
            self.fromGPU_(self.Itot_out_buff, Itot_out)
            self.queue.finish()
            return E_out, Itot_out, occ_out
        
        # If not finalizing, just finish the kernel queue and return
        self.queue.finish()
        print("solve_mc_2phase() DONE")
        return None


    def solve(
        self,
        posE :    np.ndarray = None,  # shape (nSingle,4), float32  (xyz , E0)
        pTips:    np.ndarray = None,  # shape (nTips, 4), float32  (xyz , Vbias)
        bRealloc: bool = True,
        nSingle      = None,
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

# New helper functions for demo_local_update
def find_closest_pTip(posE, pTips, nxy_scan):
    """Find the closest tip position index for each site."""
    nx, ny = nxy_scan
    # Reshape tip positions (x,y) into grid
    grid_xy = pTips[:, :2].reshape((nx, ny, 2))
    indices = []
    for x_site, y_site, *_ in posE:
        dx = grid_xy[:, :, 0] - x_site
        dy = grid_xy[:, :, 1] - y_site
        idx_flat = np.argmin(dx*dx + dy*dy)
        ix, iy = np.unravel_index(idx_flat, (nx, ny))
        indices.append(ix*ny + iy)
    return np.array(indices, dtype=int)


def print_occupancy(occupation, nSingle, bHex=False):
    occ_bytes = nSingle//8
    nTips = occupation.shape[0]
    for i in range(nTips):
        print(f"iTip {i:3} occ: ", end="")
        for j in range(occ_bytes):
            if bHex:
                print(f"{occupation[i,j]:02x}", end="")
            else:
                print(f"{occupation[i,j]:08b}", end="")
        print()

def print_site_maps(Esite_map, Tsite_map, occ_map_bits, total_energy, total_current, title=""):
    """
    Prints Esite, Tsite, and occupancy for a single tip configuration in a textual format.
    This provides a quantitative view analogous to plot_site_maps_imshow.

    Args:
        Esite_map (np.ndarray): 2D array of on-site energies, shape (nx, ny).
        Tsite_map (np.ndarray): 2D array of hopping amplitudes, shape (nx, ny).
        occ_map (np.ndarray): 2D array of unpacked occupancies (0 or 1), shape (nx, ny).
        total_energy (float): The total energy for this configuration.
        total_current (np.ndarray): The (I_occ, I_unocc) currents for this configuration.
        title (str, optional): A title for the printout.
    """
    if title: print(f"\n--- {title} ---")
    print(f"Total Energy: {total_energy:.4f} eV, Current (occ, unocc): ({total_current[0]:.4e}, {total_current[1]:.4e})")
    print("-" * 80)
    print(f"{'Site ':<7} | {'E[eV]':<10} | {'T':<10} | {'n'}")
    print("-" * 80)
    #nx, ny = Esite_map.shape
    # The map is stored (nx, ny) but we want to print it as a grid of ny rows and nx columns
    nsite = len(Esite_map)
    #print("occ_map.shape", occ_map.shape())
    #occ_map_ = np.unpackbits(occ_map_bits, axis=1)
    for i in range(nsite):
        e_site, t_site, occ = Esite_map[i], Tsite_map[i], occ_map_bits[i]
        print(f"{i:7} | {e_site:^10.4f} | {t_site:^10.4e} | {occ}")
    print("-" * 80)

