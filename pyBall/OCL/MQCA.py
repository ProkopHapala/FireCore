# run as module:  python -u -m pyBall.OCL.MQCA

import os
import numpy as np
import pyopencl as cl
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from .OpenCLBase import OpenCLBase

MAX_SITES = 16
MAX_NEIGH = 8
WG_SIZE   = 16   # one thread per site

# ---- all 16 possible binary logic functions of 2 inputs ----
# Bit encoding: code = out(A=0,B=0) | out(A=0,B=1)<<1 | out(A=1,B=0)<<2 | out(A=1,B=1)<<3
# combos in order: (0,0),(0,1),(1,0),(1,1)
# Truth tables:
#   FALSE  = [0,0,0,0] = 0
#   AND    = [0,0,0,1] = 8
#   A_GT_B = [0,0,1,0] = 4     A and not B
#   A      = [0,0,1,1] = 12
#   B_GT_A = [0,1,0,0] = 2     B and not A
#   B      = [0,1,0,1] = 10
#   XOR    = [0,1,1,0] = 6
#   OR     = [0,1,1,1] = 14
#   NOR    = [1,0,0,0] = 1
#   XNOR   = [1,0,0,1] = 9
#   NOT_B  = [1,0,1,0] = 5
#   A_GE_B = [1,0,1,1] = 13    A or not B
#   NOT_A  = [1,1,0,0] = 3
#   B_GE_A = [1,1,0,1] = 11    B or not A
#   NAND   = [1,1,1,0] = 7
#   TRUE   = [1,1,1,1] = 15
LOGIC_NAMES = {
     0: 'FALSE',
     8: 'AND',
     4: 'A_GT_B',
    12: 'A',
     2: 'B_GT_A',
    10: 'B',
     6: 'XOR',
    14: 'OR',
     1: 'NOR',
     9: 'XNOR',
     5: 'NOT_B',
    13: 'A_GE_B',
     3: 'NOT_A',
    11: 'B_GE_A',
     7: 'NAND',
    15: 'TRUE',
}

USEFUL_LOGIC = {'AND', 'OR', 'XOR', 'NAND', 'NOR', 'XNOR'}

# =====================================================================
#  Square-lattice helper:  build sparse W matrix from (W1, W2)
# =====================================================================

def sq_lattice_sparse(positions, W1, W2, nSite=None):
    """
    Build sparse coupling arrays for a square-lattice cluster.

    positions : (nSite, 2) integer grid coordinates
    W1 : coupling along Cartesian directions (distance=1)
    W2 : coupling along diagonal directions  (distance=sqrt(2))

    Returns (W_val, W_idx, nNeigh) all shape [nSite, MAX_NEIGH]
    """
    if nSite is None: nSite = len(positions)
    pos = np.array(positions, dtype=int)
    W_val  = np.zeros((nSite, MAX_NEIGH), dtype=np.float32)
    W_idx  = np.zeros((nSite, MAX_NEIGH), dtype=np.int32)
    nNeigh = np.zeros(nSite,              dtype=np.int32)
    for i in range(nSite):
        nn = 0
        for j in range(nSite):
            if i == j: continue
            dx = abs(pos[i,0] - pos[j,0])
            dy = abs(pos[i,1] - pos[j,1])
            if dx <= 1 and dy <= 1:
                if dx + dy == 1:
                    w = W1
                else:                  # dx==dy==1  diagonal
                    w = W2
                if nn < MAX_NEIGH:
                    W_val[i, nn] = w
                    W_idx[i, nn] = j
                    nn += 1
        nNeigh[i] = nn
    return W_val, W_idx, nNeigh


def compute_input_bias(positions, input_positions, input_neighbors, input_vals, W1, W2, shift=0.0):
    """
    Compute bias on active sites from fixed input pads, using proper W1/W2 coupling.

    positions       : (nSite, 2) active site grid coordinates
    input_positions : (2,) list of (x,y) coordinates for input pads A and B
    input_neighbors : list of active-site indices adjacent to each input pad
                      e.g. [[0], [2]] means input A neighbors site 0, input B neighbors site 2
    input_vals      : (2,) int 0 or 1 for each input
    W1, W2          : Cartesian and diagonal coupling strengths
    shift           : float, shift applied to input values (default 0.0)
                      shift=0.0: 0→0, 1→1 (no shift)
                      shift=-0.5: 0→-0.5, 1→+0.5 (centered)
                      shift=-1.0: 0→-1, 1→0 (spin mapping)

    Returns bias array (nSite,) to add to Esite_base.
    """
    bias = np.zeros(len(positions), dtype=np.float32)
    positions = np.array(positions, dtype=int)
    input_positions = np.array(input_positions, dtype=int)

    for inp_idx, sites in enumerate(input_neighbors):
        input_val = input_vals[inp_idx]
        # Apply shift: input_val + shift
        spin_val = float(input_val) + shift
        inp_pos = input_positions[inp_idx]
        for s in sites:
            site_pos = positions[s]
            dx = abs(inp_pos[0] - site_pos[0])
            dy = abs(inp_pos[1] - site_pos[1])
            # Determine coupling based on geometry
            if dx + dy == 1:
                w = W1  # Cartesian neighbor
            elif dx == 1 and dy == 1:
                w = W2  # Diagonal neighbor
            else:
                w = 0.0  # Should not happen for properly placed inputs
            bias[s] += w * spin_val
    return bias


def apply_input_bias(Esite_base, positions, input_positions, input_neighbors, input_vals, W1, W2, shift=0.0):
    """
    Return a copy of Esite_base with input bias applied using proper W coupling.
    """
    Esite = Esite_base.copy()
    bias = compute_input_bias(positions, input_positions, input_neighbors, input_vals, W1, W2, shift)
    Esite += bias
    return Esite


def identify_logic(outputs_4):
    """
    outputs_4 : sequence of 4 ints [n(0,0), n(0,1), n(1,0), n(1,1)]
    Returns (code, name) where code is 0-15.
    """
    code = (int(outputs_4[0]) & 1)       \
         | ((int(outputs_4[1]) & 1) << 1) \
         | ((int(outputs_4[2]) & 1) << 2) \
         | ((int(outputs_4[3]) & 1) << 3)
    return code, LOGIC_NAMES.get(code, '?')


def occ_mask_to_array(occ_mask, nSite):
    """Convert integer bitmask to (nSite,) int array."""
    return np.array([(occ_mask >> i) & 1 for i in range(nSite)], dtype=np.int32)


# =====================================================================
#  MQCASolver  –  thin OpenCL wrapper
# =====================================================================

class MQCASolver(OpenCLBase):
    """
    Brute-force Gray-code MQCA ground-state solver.

    Usage
    -----
    solver = MQCASolver()
    E_min, occ_min = solver.solve(Esite_batch, W_val, W_idx, nNeigh, nSite)
    """

    def __init__(self, device_index=0, preferred_vendor='nvidia', bPrint=True):
        super().__init__(nloc=WG_SIZE, device_index=device_index,
                         preferred_vendor=preferred_vendor, bPrint=bPrint)
        base_path = os.path.dirname(os.path.abspath(__file__))
        if not self.load_program(rel_path='cl/MQCA.cl', base_path=base_path):
            raise RuntimeError('[MQCASolver] Failed to compile cl/MQCA.cl')
        self._nInst_alloc = 0
        self._nSite_alloc = 0
        self._k_gs   = cl.Kernel(self.prg, 'mqca_groundstate')
        self._k_gs_W = cl.Kernel(self.prg, 'mqca_groundstate_batch_W')
        self._k_gs_top8 = None  # Loaded on demand

    # ------------------------------------------------------------------
    def _realloc(self, nInst, nSite):
        if nInst == self._nInst_alloc and nSite == self._nSite_alloc:
            return
        sz_f = np.dtype(np.float32).itemsize
        sz_i = np.dtype(np.int32  ).itemsize
        buffs = {
            'Esite'  : sz_f * nInst * nSite,
            'W_val'  : sz_f * nSite * MAX_NEIGH,
            'W_idx'  : sz_i * nSite * MAX_NEIGH,
            'nNeigh' : sz_i * nSite,
            'E_out'  : sz_f * nInst,
            'occ_out': sz_i * nInst,
        }
        self.try_make_buffers(buffs)
        self._nInst_alloc = nInst
        self._nSite_alloc = nSite

    def _realloc_batch_W(self, nInst, nSite):
        sz_f = np.dtype(np.float32).itemsize
        sz_i = np.dtype(np.int32  ).itemsize
        buffs = {
            'bw_Esite'  : sz_f * nInst * nSite,
            'bw_W_val'  : sz_f * nInst * nSite * MAX_NEIGH,
            'bw_W_idx'  : sz_i * nInst * nSite * MAX_NEIGH,
            'bw_nNeigh' : sz_i * nInst * nSite,
            'bw_E_out'  : sz_f * nInst,
            'bw_occ_out': sz_i * nInst,
        }
        self.try_make_buffers(buffs)

    # ------------------------------------------------------------------
    def solve(self, Esite_batch, W_val, W_idx, nNeigh, nSite):
        """
        Solve ground state for nInstances independent systems
        sharing the same W geometry.

        Esite_batch : (nInst, nSite) float32 – on-site energies per instance
        W_val       : (nSite, MAX_NEIGH) float32
        W_idx       : (nSite, MAX_NEIGH) int32
        nNeigh      : (nSite,) int32
        nSite       : int

        Returns
        -------
        E_min  : (nInst,) float32
        occ_min: (nInst,) int32  – 16-bit bitmask of ground-state occupancy
        """
        Esite_batch = np.ascontiguousarray(Esite_batch, dtype=np.float32)
        W_val  = np.ascontiguousarray(W_val,  dtype=np.float32)
        W_idx  = np.ascontiguousarray(W_idx,  dtype=np.int32)
        nNeigh = np.ascontiguousarray(nNeigh, dtype=np.int32)
        if Esite_batch.ndim == 1:
            Esite_batch = Esite_batch[np.newaxis, :]
        nInst = Esite_batch.shape[0]

        self._realloc(nInst, nSite)
        self.toGPU_( self.Esite_buff,   Esite_batch.ravel())
        self.toGPU_( self.W_val_buff,   W_val.ravel())
        self.toGPU_( self.W_idx_buff,   W_idx.ravel())
        self.toGPU_( self.nNeigh_buff,  nNeigh.ravel())

        kernel = self._k_gs
        kernel.set_args(
            np.int32(nSite),
            np.int32(nInst),
            self.Esite_buff,
            self.W_val_buff,
            self.W_idx_buff,
            self.nNeigh_buff,
            self.E_out_buff,
            self.occ_out_buff,
        )
        global_size = (nInst * WG_SIZE,)
        local_size  = (WG_SIZE,)
        cl.enqueue_nd_range_kernel(self.queue, kernel, global_size, local_size)
        self.queue.finish()

        E_min   = np.empty(nInst, dtype=np.float32)
        occ_min = np.empty(nInst, dtype=np.int32)
        self.fromGPU_(self.E_out_buff,   E_min)
        self.fromGPU_(self.occ_out_buff, occ_min)
        return E_min, occ_min

    # ------------------------------------------------------------------
    def solve_batch_W(self, Esite_batch, W_val_batch, W_idx_batch, nNeigh_batch, nSite):
        """
        Variant where each instance has its own W matrix.

        Esite_batch  : (nInst, nSite)
        W_val_batch  : (nInst, nSite, MAX_NEIGH)
        W_idx_batch  : (nInst, nSite, MAX_NEIGH)
        nNeigh_batch : (nInst, nSite)
        """
        Esite_batch  = np.ascontiguousarray(Esite_batch,  dtype=np.float32)
        W_val_batch  = np.ascontiguousarray(W_val_batch,  dtype=np.float32)
        W_idx_batch  = np.ascontiguousarray(W_idx_batch,  dtype=np.int32)
        nNeigh_batch = np.ascontiguousarray(nNeigh_batch, dtype=np.int32)
        if Esite_batch.ndim == 1:
            Esite_batch = Esite_batch[np.newaxis, :]
        nInst = Esite_batch.shape[0]

        self._realloc_batch_W(nInst, nSite)
        self.toGPU_(self.bw_Esite_buff,   Esite_batch.ravel())
        self.toGPU_(self.bw_W_val_buff,   W_val_batch.ravel())
        self.toGPU_(self.bw_W_idx_buff,   W_idx_batch.ravel())
        self.toGPU_(self.bw_nNeigh_buff,  nNeigh_batch.ravel())

        kernel = self._k_gs_W
        kernel.set_args(
            np.int32(nSite),
            np.int32(nInst),
            self.bw_Esite_buff,
            self.bw_W_val_buff,
            self.bw_W_idx_buff,
            self.bw_nNeigh_buff,
            self.bw_E_out_buff,
            self.bw_occ_out_buff,
        )
        global_size = (nInst * WG_SIZE,)
        local_size  = (WG_SIZE,)
        cl.enqueue_nd_range_kernel(self.queue, kernel, global_size, local_size)
        self.queue.finish()

        E_min   = np.empty(nInst, dtype=np.float32)
        occ_min = np.empty(nInst, dtype=np.int32)
        self.fromGPU_(self.bw_E_out_buff,   E_min)
        self.fromGPU_(self.bw_occ_out_buff, occ_min)
        return E_min, occ_min

    # ------------------------------------------------------------------
    def _ensure_top8_kernel(self):
        """Load top8 kernel on demand (separate program to avoid conflicts)."""
        if self._k_gs_top8 is None:
            base_path = os.path.dirname(os.path.abspath(__file__))
            kernel_src = open(os.path.join(base_path, 'cl', 'MQCA_top8.cl')).read()
            self._prg_top8 = cl.Program(self.ctx, kernel_src).build()
            self._k_gs_top8 = cl.Kernel(self._prg_top8, 'mqca_groundstate_top8')

    def solve_batch_W_top8(self, Esite_batch, W_val_batch, W_idx_batch, nNeigh_batch, nSite):
        """
        Variant tracking top-8 lowest energy states (for degeneracy analysis).

        Returns
        -------
        E_top8   : (nInst, 8) float32 - energies of 8 lowest states (sorted)
        occ_top8 : (nInst, 8) int32   - occupancy masks of 8 lowest states
        """
        self._ensure_top8_kernel()

        Esite_batch  = np.ascontiguousarray(Esite_batch,  dtype=np.float32)
        W_val_batch  = np.ascontiguousarray(W_val_batch,  dtype=np.float32)
        W_idx_batch  = np.ascontiguousarray(W_idx_batch,  dtype=np.int32)
        nNeigh_batch = np.ascontiguousarray(nNeigh_batch, dtype=np.int32)
        if Esite_batch.ndim == 1:
            Esite_batch = Esite_batch[np.newaxis, :]
        nInst = Esite_batch.shape[0]

        self._realloc_batch_W(nInst, nSite)
        # Allocate top8 output buffers
        sz_f = np.dtype(np.float32).itemsize
        sz_i = np.dtype(np.int32).itemsize
        if not hasattr(self, '_top8_buff_alloc') or self._top8_buff_alloc < nInst:
            self.bw_E_top8_buff   = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, sz_f * nInst * 8)
            self.bw_occ_top8_buff = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, sz_i * nInst * 8)
            self._top8_buff_alloc = nInst

        self.toGPU_(self.bw_Esite_buff,   Esite_batch.ravel())
        self.toGPU_(self.bw_W_val_buff,   W_val_batch.ravel())
        self.toGPU_(self.bw_W_idx_buff,   W_idx_batch.ravel())
        self.toGPU_(self.bw_nNeigh_buff,  nNeigh_batch.ravel())

        kernel = self._k_gs_top8
        kernel.set_args(
            np.int32(nSite),
            np.int32(nInst),
            self.bw_Esite_buff,
            self.bw_W_val_buff,
            self.bw_W_idx_buff,
            self.bw_nNeigh_buff,
            self.bw_E_top8_buff,
            self.bw_occ_top8_buff,
        )
        global_size = (nInst * WG_SIZE,)
        local_size  = (WG_SIZE,)
        cl.enqueue_nd_range_kernel(self.queue, kernel, global_size, local_size)
        self.queue.finish()

        E_top8   = np.empty((nInst, 8), dtype=np.float32)
        occ_top8 = np.empty((nInst, 8), dtype=np.int32)
        self.fromGPU_(self.bw_E_top8_buff,   E_top8.ravel())
        self.fromGPU_(self.bw_occ_top8_buff, occ_top8.ravel())
        return E_top8, occ_top8


# =====================================================================
#  High-level: logic-gate evaluation helpers
# =====================================================================

INPUT_COMBOS = [(0,0), (0,1), (1,0), (1,1)]

def check_ground_state_uniqueness(E_top8, threshold=0.01):
    """
    Check if ground state is well-separated from first excited state.

    Parameters
    ----------
    E_top8 : (nInst, 8) float32
        Energies of 8 lowest states from solve_batch_W_top8
    threshold : float
        Minimum energy gap to consider ground state well-defined

    Returns
    -------
    unique : (nInst,) bool
        True if ground state gap >= threshold (ground state is unique)
    ground_gaps : (nInst,) float
        Energy gaps between ground and first excited state
    """
    ground_gaps = E_top8[:, 1] - E_top8[:, 0]
    unique = ground_gaps >= threshold
    return unique, ground_gaps


def eval_logic_table(solver, Esite_base, W_val, W_idx, nNeigh, nSite,
                     positions, input_positions, input_neighbors, output_site, W1, W2, shift=0.0):
    """
    Evaluate logic truth-table for one cluster at one (W1,W2) point.

    positions       : (nSite, 2) for computing input coupling geometry
    input_positions : (2,) list of (x,y) for input pads A and B
    shift           : float, shift applied to input values (default 0.0)

    Returns
    -------
    outputs_4 : (4,) int   output bit for each input combination
    occ_4     : (4, nSite) occupancy arrays for each input combination
    E_4       : (4,) float ground-state energies
    Esite_4   : (4, nSite) on-site energies (for debugging)
    logic_code, logic_name
    """
    Esite_batch = np.zeros((4, nSite), dtype=np.float32)
    for k, (A, B) in enumerate(INPUT_COMBOS):
        Esite_batch[k] = apply_input_bias(Esite_base, positions, input_positions, input_neighbors, [A, B], W1, W2, shift)

    E_4, occ_raw = solver.solve(Esite_batch, W_val, W_idx, nNeigh, nSite)
    outputs_4 = np.array([(occ_raw[k] >> output_site) & 1 for k in range(4)], dtype=np.int32)
    occ_4     = np.array([occ_mask_to_array(occ_raw[k], nSite) for k in range(4)])
    code, name = identify_logic(outputs_4)
    return outputs_4, occ_4, E_4, Esite_batch, code, name


def scan_W1_W2(solver, positions, input_positions, Esite_base, input_neighbors, output_site,
               W1_vals, W2_vals, nSite=None, shift=0.0):
    """
    Scan (W1, W2) parameter space for a fixed cluster geometry.

    shift : float, shift applied to input values (default 0.0)

    Returns
    -------
    logic_map : (nW2, nW1) int   0-15 logic code at each (W1,W2)
    occ_map   : (nW2, nW1, 4, nSite) occupancy for each combo
    """
    if nSite is None: nSite = len(positions)
    nW1, nW2 = len(W1_vals), len(W2_vals)
    nInst = nW1 * nW2 * 4

    # Build Esite batch: shape (nW2, nW1, 4, nSite) → flat (nInst, nSite)
    Esite_batch = np.zeros((nW2, nW1, 4, nSite), dtype=np.float32)
    for iy, W2v in enumerate(W2_vals):
        for ix, W1v in enumerate(W1_vals):
            for k, (A, B) in enumerate(INPUT_COMBOS):
                Esite_batch[iy, ix, k] = apply_input_bias(Esite_base, positions, input_positions, input_neighbors, [A, B], W1v, W2v, shift)

    # W coupling is same per (W1,W2) → replicate into full flat batch
    # Flatten Esite: (nW2*nW1*4, nSite)
    Esite_flat = Esite_batch.reshape(nInst, nSite)

    # Build W batches: (nInst, nSite, MAX_NEIGH)
    W_val_batch  = np.zeros((nInst, nSite, MAX_NEIGH), dtype=np.float32)
    W_idx_batch  = np.zeros((nInst, nSite, MAX_NEIGH), dtype=np.int32)
    nNeigh_batch = np.zeros((nInst, nSite),             dtype=np.int32)
    inst = 0
    for iy, W2 in enumerate(W2_vals):
        for ix, W1 in enumerate(W1_vals):
            Wv, Wi, Wn = sq_lattice_sparse(positions, W1, W2, nSite)
            for k in range(4):
                W_val_batch [inst] = Wv
                W_idx_batch [inst] = Wi
                nNeigh_batch[inst] = Wn
                inst += 1

    E_flat, occ_flat = solver.solve_batch_W(Esite_flat, W_val_batch, W_idx_batch, nNeigh_batch, nSite)

    # Reshape outputs
    E_4d   = E_flat  .reshape(nW2, nW1, 4)
    occ_4d = occ_flat.reshape(nW2, nW1, 4)

    logic_map = np.zeros((nW2, nW1), dtype=np.int32)
    occ_map   = np.zeros((nW2, nW1, 4, nSite), dtype=np.int32)
    for iy in range(nW2):
        for ix in range(nW1):
            outs  = [(occ_4d[iy,ix,k] >> output_site) & 1 for k in range(4)]
            code, _ = identify_logic(outs)
            logic_map[iy, ix] = code
            for k in range(4):
                occ_map[iy, ix, k] = occ_mask_to_array(occ_4d[iy,ix,k], nSite)

    return logic_map, occ_map, E_4d


def scan_W1_W2_top8(solver, positions, input_positions, Esite_base, input_neighbors, output_site,
                   W1_vals, W2_vals, nSite=None, degeneracy_threshold=0.01, shift=0.0):
    """
    Scan (W1, W2) parameter space and check ground state uniqueness using top8 kernel.

    shift : float, shift applied to input values (default 0.0)

    Returns
    -------
    logic_map : (nW2, nW1) int   0-15 logic code at each (W1,W2)
    degenerate_mask : (nW2, nW1) bool  True if ground state is degenerate
    """
    if nSite is None: nSite = len(positions)
    nW1, nW2 = len(W1_vals), len(W2_vals)
    nInst = nW1 * nW2 * 4

    # Build Esite batch: shape (nW2, nW1, 4, nSite) → flat (nInst, nSite)
    Esite_batch = np.zeros((nW2, nW1, 4, nSite), dtype=np.float32)
    for iy, W2v in enumerate(W2_vals):
        for ix, W1v in enumerate(W1_vals):
            for k, (A, B) in enumerate(INPUT_COMBOS):
                Esite_batch[iy, ix, k] = apply_input_bias(Esite_base, positions, input_positions, input_neighbors, [A, B], W1v, W2v, shift)

    # Flatten Esite: (nW2*nW1*4, nSite)
    Esite_flat = Esite_batch.reshape(nInst, nSite)

    # Build W batches: (nInst, nSite, MAX_NEIGH)
    W_val_batch  = np.zeros((nInst, nSite, MAX_NEIGH), dtype=np.float32)
    W_idx_batch  = np.zeros((nInst, nSite, MAX_NEIGH), dtype=np.int32)
    nNeigh_batch = np.zeros((nInst, nSite),             dtype=np.int32)
    inst = 0
    for iy, W2 in enumerate(W2_vals):
        for ix, W1 in enumerate(W1_vals):
            Wv, Wi, Wn = sq_lattice_sparse(positions, W1, W2, nSite)
            for k in range(4):
                W_val_batch [inst] = Wv
                W_idx_batch [inst] = Wi
                nNeigh_batch[inst] = Wn
                inst += 1

    # Use top8 kernel to get ground state uniqueness info
    E_top8, occ_top8 = solver.solve_batch_W_top8(Esite_flat, W_val_batch, W_idx_batch, nNeigh_batch, nSite)

    # Check ground state uniqueness
    unique, ground_gaps = check_ground_state_uniqueness(E_top8, threshold=degeneracy_threshold)

    # Reshape to (nW2, nW1, 4)
    unique_4d = unique.reshape(nW2, nW1, 4)
    occ_4d = occ_top8[:, 0].reshape(nW2, nW1, 4)  # Use ground state occupancy

    # Logic map (use ground state)
    logic_map = np.zeros((nW2, nW1), dtype=np.int32)
    degenerate_mask = np.zeros((nW2, nW1), dtype=bool)

    for iy in range(nW2):
        for ix in range(nW1):
            # Check if any input combo has degenerate ground state
            if not np.all(unique_4d[iy, ix]):
                degenerate_mask[iy, ix] = True
                # Still compute logic from ground state
            outs = [(occ_4d[iy, ix, k] >> output_site) & 1 for k in range(4)]
            code, _ = identify_logic(outs)
            logic_map[iy, ix] = code

    return logic_map, degenerate_mask


# =====================================================================
#  Visualization helpers
# =====================================================================

# Assign a distinct colour to each of the 16 logic codes
_LOGIC_CMAP = matplotlib.colormaps.get_cmap('tab20')
LOGIC_COLORS = {code: _LOGIC_CMAP(code / 15.0) for code in range(16)}


def plot_cluster(ax, positions, occ, input_positions, output_site,
                 W_val=None, W_idx=None, nNeigh=None, title='', tetris_style=False,
                 input_values=None, Esite=None):
    """
    Draw one cluster configuration.

    positions       : (nSite,2) active site grid coords
    occ             : (nSite,) occupancy 0/1
    input_positions : list of (x,y) for input pads (not active)
    output_site     : int index of output site
    tetris_style    : if True, use large adjacent squares like Tetris blocks
    input_values    : (2,) int array [A, B] values 0 or 1 for each input pad
    """
    pos = np.array(positions, dtype=float)
    nSite = len(pos)
    if input_values is None:
        input_values = [0, 0]

    if tetris_style:
        # Tetris-style: large adjacent squares, no lines
        # Determine grid cell size from spacing
        if len(pos) > 1:
            # Use median spacing as cell size
            dxs = np.abs(pos[:,0:1] - pos[:,0:1].T)
            dys = np.abs(pos[:,1:2] - pos[:,1:2].T)
            dxs = dxs[dxs > 0.1]
            dys = dys[dys > 0.1]
            cell = min(np.median(dxs) if len(dxs) > 0 else 1.0,
                       np.median(dys) if len(dys) > 0 else 1.0)
        else:
            cell = 1.0

        # Colors: ON = dark red, OFF = light blue
        color_on = '#c44e4e'   # muted red
        color_off = '#8cb3d9'  # muted blue

        # Draw active sites as squares
        for i, (x, y) in enumerate(pos):
            is_on = occ[i] if occ is not None else False
            facecolor = color_on if is_on else color_off
            is_output = (i == output_site)
            is_input_adj = i in getattr(plot_cluster, '_input_neighbor_sites', [])

            # Main square
            rect = mpatches.Rectangle((x - cell*0.4, y - cell*0.4), cell*0.8, cell*0.8,
                                       linewidth=2, edgecolor='black', facecolor=facecolor, zorder=3)
            ax.add_patch(rect)

            # Markers inside: small dot for input-adjacent, circle for output
            if is_input_adj:
                # Small black dot indicates this site receives input bias
                dot = plt.Circle((x, y), cell*0.12, fill=True, facecolor='black', edgecolor='none', zorder=4)
                ax.add_patch(dot)
            elif is_output:
                # Draw circle
                circ = plt.Circle((x, y), cell*0.25, fill=False, edgecolor='black', linewidth=2, zorder=4)
                ax.add_patch(circ)

            # Display: site number + on-site energy (if provided)
            if Esite is not None:
                label_text = f'{i}\nε={Esite[i]:.2f}'
                fs = 7
            else:
                label_text = str(i)
                fs = 8
            ax.text(x, y, label_text, ha='center', va='center', fontsize=fs, color='white', fontweight='bold', zorder=5)

        # Draw input pads (external) as squares colored by their binary value
        for idx, (x, y) in enumerate(input_positions):
            val = input_values[idx] if idx < len(input_values) else 0
            # Color based on input value: red=1, blue=0
            pad_color = color_on if val else color_off
            rect = mpatches.Rectangle((x - cell*0.35, y - cell*0.35), cell*0.7, cell*0.7,
                                       linewidth=2, edgecolor='black', facecolor=pad_color, zorder=3)
            ax.add_patch(rect)
            # Label A or B
            label = 'A' if idx == 0 else 'B'
            ax.text(x, y, label, ha='center', va='center', fontsize=10, color='white', fontweight='bold', zorder=5)

        # Set limits with padding
        if len(pos) > 0:
            x_min, x_max = pos[:,0].min(), pos[:,0].max()
            y_min, y_max = pos[:,1].min(), pos[:,1].max()
            if input_positions:
                inp = np.array(input_positions)
                x_min, x_max = min(x_min, inp[:,0].min()), max(x_max, inp[:,0].max())
                y_min, y_max = min(y_min, inp[:,1].min()), max(y_max, inp[:,1].max())
            ax.set_xlim(x_min - cell, x_max + cell)
            ax.set_ylim(y_min - cell, y_max + cell)

    else:
        # Original visualization with circles/scatter
        # draw coupling bonds
        if W_val is not None and W_idx is not None and nNeigh is not None:
            for i in range(nSite):
                for k in range(nNeigh[i]):
                    j = W_idx[i, k]
                    if j > i:
                        ax.plot([pos[i,0], pos[j,0]], [pos[i,1], pos[j,1]],
                                'k-', lw=1.0, alpha=0.3, zorder=1)

        # draw active sites
        for i, (x, y) in enumerate(pos):
            c = 'red' if occ[i] else 'steelblue'
            mk = 's' if i == output_site else 'o'
            sz = 220 if i == output_site else 160
            ax.scatter(x, y, c=c, marker=mk, s=sz, edgecolors='black', linewidths=1.2, zorder=3)
            ax.text(x, y, str(i), ha='center', va='center', fontsize=7, color='white', zorder=4)

        # draw input pads
        for x, y in input_positions:
            ax.scatter(x, y, c='limegreen', marker='^', s=260, edgecolors='black', linewidths=1.2, zorder=3)

    ax.set_title(title, fontsize=9)
    ax.set_aspect('equal')
    ax.axis('off')


def plot_ground_states(positions, occ_4, E_4, outputs_4, input_positions,
                       output_site, logic_name,
                       W_val=None, W_idx=None, nNeigh=None,
                       input_neighbors=None, Esite_4=None,
                       W1=None, W2=None, eps0=None,
                       fname=None, show=False, tetris_style=False):
    """
    2×2 subplot showing ground state for each of the 4 input combinations.

    input_neighbors: list of lists, e.g. [[0], [2]] for sites 0 and 2 being input-adjacent
    Esite_4: (4, nSite) array of on-site energies for each input combo
    tetris_style: if True, use square-grid Tetris visualization
    W1, W2, eps0: parameters to display in title
    """
    # Pass input neighbor info to plot_cluster via function attribute
    if input_neighbors is not None:
        input_neighbor_sites = set()
        for sites in input_neighbors:
            input_neighbor_sites.update(sites)
        plot_cluster._input_neighbor_sites = list(input_neighbor_sites)

    fig, axes = plt.subplots(2, 2, figsize=(8, 7))

    # Main title with parameters
    param_str = ""
    if W1 is not None:
        param_str += f"  W1={W1:.2f}"
    if W2 is not None:
        param_str += f"  W2={W2:.2f}"
    if eps0 is not None:
        param_str += f"  ε0={eps0:.2f}"
    fig.suptitle(f'Cluster: {logic_name}{param_str}', fontsize=11)

    for k, ax in enumerate(axes.flat):
        A, B = INPUT_COMBOS[k]
        title = f'In({A},{B}) → Out={int(outputs_4[k])}  E={E_4[k]:.3f}'
        Esite_k = Esite_4[k] if Esite_4 is not None else None
        plot_cluster(ax, positions, occ_4[k], input_positions, output_site,
                     W_val, W_idx, nNeigh, title=title, tetris_style=tetris_style,
                     input_values=[A, B], Esite=Esite_k)

    # Legend
    if tetris_style:
        handles = [
            mpatches.Rectangle((0,0), 1, 1, facecolor='#c44e4e', edgecolor='black', label='ON (n=1)'),
            mpatches.Rectangle((0,0), 1, 1, facecolor='#8cb3d9', edgecolor='black', label='OFF (n=0)'),
            mpatches.Patch(facecolor='white', edgecolor='black', label='Output: ○  In-bias: ●'),
        ]
    else:
        handles = [
            mpatches.Patch(color='red',       label='Occupied (n=1)'),
            mpatches.Patch(color='steelblue', label='Empty    (n=0)'),
            mpatches.Patch(color='limegreen', label='Input pad (fixed)'),
            plt.scatter([], [], marker='s', c='white', edgecolors='black', s=80, label='Output site'),
        ]
    fig.legend(handles=handles, loc='lower center', ncol=4, fontsize=8, frameon=False)
    plt.tight_layout(rect=[0, 0.07, 1, 1])
    if fname: plt.savefig(fname, dpi=150)
    if show:  plt.show()
    plt.close()

    # Clean up
    plot_cluster._input_neighbor_sites = []


def plot_logic_map(W1_vals, W2_vals, logic_map, title='Logic phase diagram',
                   fname=None, show=False):
    """
    imshow of logic_map with annotated logic names.
    """
    nW2, nW1 = logic_map.shape
    # Build RGB image
    img = np.zeros((nW2, nW1, 3))
    for code in range(16):
        mask = logic_map == code
        c = LOGIC_COLORS[code][:3]
        img[mask] = c

    fig, ax = plt.subplots(figsize=(8, 6))
    extent = [W1_vals[0], W1_vals[-1], W2_vals[0], W2_vals[-1]]
    ax.imshow(img, origin='lower', extent=extent, aspect='auto', interpolation='nearest')
    ax.set_xlabel('W1 (Cartesian coupling)', fontsize=11)
    ax.set_ylabel('W2 (Diagonal coupling)',  fontsize=11)
    ax.set_title(title, fontsize=12)

    # Legend patches for codes that actually appear
    present = np.unique(logic_map)
    patches = [mpatches.Patch(color=LOGIC_COLORS[c], label=LOGIC_NAMES[c]) for c in present]
    ax.legend(handles=patches, bbox_to_anchor=(1.01, 1), loc='upper left', fontsize=8, frameon=True)
    plt.tight_layout()
    if fname: plt.savefig(fname, dpi=150, bbox_inches='tight')
    if show:  plt.show()
    plt.close()


def plot_logic_fraction_map(W1_vals, W2_vals, logic_map, target_set=None,
                             title='Useful logic fraction', fname=None, show=False):
    """
    Show which (W1,W2) regions produce any useful logic function.
    target_set: set of logic names to highlight; defaults to USEFUL_LOGIC.
    """
    if target_set is None: target_set = USEFUL_LOGIC
    target_codes = {c for c, n in LOGIC_NAMES.items() if n in target_set}
    useful_mask = np.isin(logic_map, list(target_codes))

    fig, ax = plt.subplots(figsize=(7, 5))
    extent = [W1_vals[0], W1_vals[-1], W2_vals[0], W2_vals[-1]]
    ax.imshow(useful_mask.astype(float), origin='lower', extent=extent,
              aspect='auto', cmap='RdYlGn', vmin=0, vmax=1, interpolation='nearest')
    ax.set_xlabel('W1', fontsize=11)
    ax.set_ylabel('W2', fontsize=11)
    ax.set_title(title, fontsize=12)
    plt.tight_layout()
    if fname: plt.savefig(fname, dpi=150)
    if show:  plt.show()
    plt.close()
