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


def apply_input_bias(Esite_base, input_neighbors, input_vals, W_in):
    """
    Return a copy of Esite_base with input bias applied.

    input_neighbors : list of active-site indices adjacent to each input pad
                      e.g. [[site_A0, site_A1, ...], [site_B0, ...]]
    input_vals      : (2,) int  0 or 1 for each input
    W_in            : coupling strength between input pad and adjacent active site
    """
    Esite = Esite_base.copy()
    for inp_idx, sites in enumerate(input_neighbors):
        bias = W_in * float(input_vals[inp_idx])
        for s in sites:
            Esite[s] += bias
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


# =====================================================================
#  High-level: logic-gate evaluation helpers
# =====================================================================

INPUT_COMBOS = [(0,0), (0,1), (1,0), (1,1)]

def eval_logic_table(solver, Esite_base, W_val, W_idx, nNeigh, nSite,
                     input_neighbors, output_site, W_in=1.0):
    """
    Evaluate logic truth-table for one cluster at one (W1,W2) point.

    Returns
    -------
    outputs_4 : (4,) int   output bit for each input combination
    occ_4     : (4, nSite) occupancy arrays for each input combination
    E_4       : (4,) float ground-state energies
    logic_code, logic_name
    """
    Esite_batch = np.zeros((4, nSite), dtype=np.float32)
    for k, (A, B) in enumerate(INPUT_COMBOS):
        Esite_batch[k] = apply_input_bias(Esite_base, input_neighbors, [A, B], W_in)

    E_4, occ_raw = solver.solve(Esite_batch, W_val, W_idx, nNeigh, nSite)
    outputs_4 = np.array([(occ_raw[k] >> output_site) & 1 for k in range(4)], dtype=np.int32)
    occ_4     = np.array([occ_mask_to_array(occ_raw[k], nSite) for k in range(4)])
    code, name = identify_logic(outputs_4)
    return outputs_4, occ_4, E_4, code, name


def scan_W1_W2(solver, positions, Esite_base, input_neighbors, output_site,
               W1_vals, W2_vals, W_in=1.0, nSite=None):
    """
    Scan (W1, W2) parameter space for a fixed cluster geometry.

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
    for iy, W2 in enumerate(W2_vals):
        for ix, W1 in enumerate(W1_vals):
            for k, (A, B) in enumerate(INPUT_COMBOS):
                Esite_batch[iy, ix, k] = apply_input_bias(Esite_base, input_neighbors, [A, B], W_in)

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


# =====================================================================
#  Visualization helpers
# =====================================================================

# Assign a distinct colour to each of the 16 logic codes
_LOGIC_CMAP = matplotlib.colormaps.get_cmap('tab20')
LOGIC_COLORS = {code: _LOGIC_CMAP(code / 15.0) for code in range(16)}


def plot_cluster(ax, positions, occ, input_positions, output_site,
                 W_val=None, W_idx=None, nNeigh=None, title=''):
    """
    Draw one cluster configuration.

    positions       : (nSite,2) active site grid coords
    occ             : (nSite,) occupancy 0/1
    input_positions : list of (x,y) for input pads (not active)
    output_site     : int index of output site
    """
    pos = np.array(positions, dtype=float)
    nSite = len(pos)

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
                       fname=None, show=False):
    """
    2×2 subplot showing ground state for each of the 4 input combinations.
    """
    fig, axes = plt.subplots(2, 2, figsize=(8, 7))
    fig.suptitle(f'Cluster ground states  |  Logic: {logic_name}', fontsize=12)
    for k, ax in enumerate(axes.flat):
        A, B = INPUT_COMBOS[k]
        title = f'In({A},{B}) → Out={int(outputs_4[k])}  E={E_4[k]:.3f}'
        plot_cluster(ax, positions, occ_4[k], input_positions, output_site,
                     W_val, W_idx, nNeigh, title=title)

    # Legend
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
