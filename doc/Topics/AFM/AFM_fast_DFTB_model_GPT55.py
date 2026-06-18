#!/usr/bin/env python3
"""
afmdftb_edgepair_demo.py

Standalone NumPy demo for a fast rigid-fragment DFTB-like short-range interaction model.

Goal
----
Estimate Pauli repulsion + local donor/acceptor chemistry between two rigid fragments A,B
without SCC/SCF and without diagonalizing the combined system during an AFM scan.

The important design decision is that the accurate P/W Pauli term is NOT simply an
atom-pair loop.  It is an *active inter-fragment edge-pair* contraction:

    e0 = (a0,b0), e1 = (a1,b1)

where a* are tip atoms and b* are surface atoms.  The two cross edges are connected
by internal fragment density blocks P_A, P_B and energy-weighted density blocks W_A,W_B.

This script is intentionally self-contained.  It uses toy P/W matrices and toy
Slater-Koster blocks so it can be run immediately.  In Fireball/DFTB production code:

    - replace make_toy_fragment() by isolated-fragment DFTB precomputation;
    - replace toy_sk_sp_block() by real SK table evaluation for S_AB,H_AB,dS,dH;
    - move build_active_edges(), evaluate_pauli_edgepair(), and evaluate_local_DA()
      to OpenCL/CUDA/WebGPU kernels.

Physics summary
---------------
For isolated fragments in a nonorthogonal AO basis:

    H_X C_X = S_X C_X eps_X,        C_X^T S_X C_X = I

Precompute the occupied density and energy-weighted density:

    P_X[mu,nu] = sum_i f_i C[mu,i] C[nu,i]
    W_X[mu,nu] = sum_i f_i eps_i C[mu,i] C[nu,i]

Do NOT define W blindly as P@H@P unless your spin/idempotency convention has been
checked.  Direct MO summation is safest.

The second-order frozen-density Pauli/orthogonalization expression used here is:

    E = T_H + T_A + T_B

    T_H = - sum_{e0,e1} Tr[ P_A[a0,a1] S[a1,b1] P_B[b1,b0] H[a0,b0]^T ]
    T_A = + 1/2 sum_{e0,e1} Tr[ S[a0,b0] P_B[b0,b1] S[a1,b1]^T W_A[a1,a0] ]
    T_B = + 1/2 sum_{e0,e1} Tr[ S[a0,b0]^T P_A[a0,a1] S[a1,b1] W_B[b1,b0] ]

This combination is invariant to a constant energy-zero shift H -> H + lambda S,
W -> W + lambda P.  The test at the end checks this numerically.

Author: ChatGPT demo for Prokop Hapala / FireCore-AFM design discussion.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional
import math
import time
import numpy as np

# Keep the toy model simple: every atom has an sp basis [s,px,py,pz].
# Fireball/DFTB production code should use per-atom orbital counts and block offsets.
NORB = 4
Block = np.ndarray


def trace4(A: Block, B: Block, C: Block, D: Block) -> float:
    """Trace of a product of four small AO blocks.

    For NORB=4 this is tiny, but in a GPU kernel this should be hand-unrolled.
    Avoid generic matrix multiplication in the innermost production kernel.
    """
    return float(np.trace(A @ B @ C @ D))


@dataclass
class BlockSparse:
    """Simple atom-block sparse matrix.

    blocks[(i,j)] is a dense [NORB,NORB] block.  For P and W we store both
    orientations explicitly: M[j,i] = M[i,j].T.  This makes GPU-like lookup simple.

    Production GPU layout should not be Python dict.  Use ELLPACK-like arrays:

        neigh[i, k]
        block[i, k, row, col]

    with fixed MAX_NEIGH for coalesced loads.
    """
    n_atoms: int
    blocks: Dict[Tuple[int, int], Block] = field(default_factory=dict)
    neigh: List[List[int]] = field(default_factory=list)

    def __post_init__(self):
        if not self.neigh:
            self.neigh = [[] for _ in range(self.n_atoms)]

    def set_block(self, i: int, j: int, M: Block):
        self.blocks[(i, j)] = np.asarray(M, dtype=np.float64)
        if j not in self.neigh[i]:
            self.neigh[i].append(j)

    def get(self, i: int, j: int) -> Optional[Block]:
        return self.blocks.get((i, j), None)


@dataclass
class LocalChannel:
    """Low-rank local donor or acceptor channel.

    This is a local AO hybrid representing a Fukui-like channel.
    Examples:
      - occupied lone-pair donor on O/N;
      - pi donor on aromatic carbon;
      - sigma* or pi* acceptor;
      - local LDOS/Fukui eigenchannel from a small atom cluster.

    eps is a representative one-electron energy.  weight is the Fukui/LDOS weight.
    """
    atom: int
    coeff: np.ndarray     # shape [NORB]
    eps: float
    weight: float = 1.0


@dataclass
class Fragment:
    name: str
    pos: np.ndarray       # [n_atoms,3]
    P: BlockSparse
    W: BlockSparse
    donors: List[LocalChannel] = field(default_factory=list)
    acceptors: List[LocalChannel] = field(default_factory=list)


@dataclass
class Edge:
    """Active inter-fragment edge e=(a,b)."""
    a: int
    b: int
    S: Block              # S_AB block [A_a orbitals, B_b orbitals]
    H: Block              # H_AB block [A_a orbitals, B_b orbitals]
    R: np.ndarray         # vector from B_b to A_a, useful for derivatives/forces


def toy_sk_sp_block(R: np.ndarray, h0: float = -8.0, r0: float = 1.8) -> Tuple[Block, Block]:
    """Toy Slater-Koster-like overlap and Hamiltonian block for an sp basis.

    This is NOT a real DFTB SK table.  It only gives direction-dependent blocks
    so the demo behaves qualitatively like an AO model.

    Production code:
      - interpolate real two-center SK tables for ss,sp,pp_sigma,pp_pi;
      - also return derivatives dS/dR, dH/dR for forces.
    """
    r = float(np.linalg.norm(R))
    if r < 1e-12:
        e = np.array([0.0, 0.0, 1.0])
    else:
        e = R / r
    # Short-range exponential/Gaussian envelope.
    amp = math.exp(-(r / r0) ** 2)

    # Toy SK radial channels.  Signs are chosen only to create angular structure.
    ss = 1.00 * amp
    sp = 0.55 * amp
    pp_sig = 0.80 * amp
    pp_pi = -0.18 * amp

    S = np.zeros((NORB, NORB), dtype=np.float64)
    S[0, 0] = ss
    S[0, 1:4] = sp * e
    S[1:4, 0] = -sp * e
    S[1:4, 1:4] = pp_pi * np.eye(3) + (pp_sig - pp_pi) * np.outer(e, e)

    # DFTB-like two-center H is roughly overlap times an energy scale, but not
    # exactly.  Add a weak anisotropic correction to avoid H being purely eta*S.
    H = h0 * S
    H[1:4, 1:4] += (0.35 * amp) * (np.eye(3) - np.outer(e, e))
    return S, H


def build_active_edges(tip: Fragment, surf: Fragment, tip_shift: np.ndarray,
                       cutoff: float, h_shift: float = 0.0) -> Tuple[List[Edge], Dict[int, List[int]]]:
    """Build sparse A-B edge list for one pixel.

    GPU target: yes.  This is one of the hottest stages.  Use cell lists for the
    surface and one workgroup per pixel.  Store edges in local/shared memory.

    h_shift is used only for the energy-zero invariance test: H_AB -> H_AB + lambda*S_AB.
    """
    edges: List[Edge] = []
    edges_by_b: Dict[int, List[int]] = {}
    for a, ra0 in enumerate(tip.pos):
        ra = ra0 + tip_shift
        # Demo uses brute-force surface scan.  Production must use cell lists / bins.
        for b, rb in enumerate(surf.pos):
            R = ra - rb
            if np.dot(R, R) <= cutoff * cutoff:
                S, H = toy_sk_sp_block(R)
                H = H + h_shift * S
                idx = len(edges)
                edges.append(Edge(a=a, b=b, S=S, H=H, R=R))
                edges_by_b.setdefault(b, []).append(idx)
    return edges, edges_by_b


def evaluate_pauli_pairdiag(tip: Fragment, surf: Fragment, tip_shift: np.ndarray,
                            cutoff: float, h_shift: float = 0.0,
                            W_tip: Optional[BlockSparse] = None,
                            W_surf: Optional[BlockSparse] = None) -> float:
    """Kernel 0: diagonal atom-pair P/W Pauli.

    Very fast and very parallel.  Good first GPU kernel and debugging baseline.
    It neglects off-diagonal density-matrix blocks P[a,a'], W[a,a'], so it misses
    density delocalization, but captures AO directionality and energy weighting.
    """
    W_tip = tip.W if W_tip is None else W_tip
    W_surf = surf.W if W_surf is None else W_surf
    edges, _ = build_active_edges(tip, surf, tip_shift, cutoff, h_shift=h_shift)
    E = 0.0
    for e in edges:
        PA = tip.P.get(e.a, e.a)
        PB = surf.P.get(e.b, e.b)
        WA = W_tip.get(e.a, e.a)
        WB = W_surf.get(e.b, e.b)
        if PA is None or PB is None or WA is None or WB is None:
            continue
        E += -trace4(PA, e.S, PB, e.H.T)
        E += 0.5 * trace4(e.S, PB, e.S.T, WA)
        E += 0.5 * trace4(e.S.T, PA, e.S, WB)
    return E


def evaluate_pauli_edgepair(tip: Fragment, surf: Fragment, tip_shift: np.ndarray,
                            cutoff: float, h_shift: float = 0.0,
                            W_tip: Optional[BlockSparse] = None,
                            W_surf: Optional[BlockSparse] = None) -> Tuple[float, Dict[str, float], int, int]:
    """Kernel 1: sparse edge-pair P/W Pauli.

    This is the main accurate contraction.  It never assembles dense S_AB or dense
    local P_B.  It loops over active cross edges e0=(a0,b0), then only over surface
    P/W neighbors b1 of b0, and active edges e1=(a1,b1).

    GPU target: yes, this should become the main per-pixel workgroup kernel.
    The Python dict/list structure is only for clarity.

    Returns:
        total energy, components, number of active edges, number of edge-pair visits.
    """
    W_tip = tip.W if W_tip is None else W_tip
    W_surf = surf.W if W_surf is None else W_surf
    edges, edges_by_b = build_active_edges(tip, surf, tip_shift, cutoff, h_shift=h_shift)

    EH = 0.0
    EWA = 0.0
    EWB = 0.0
    n_pair = 0

    for e0 in edges:
        a0, b0 = e0.a, e0.b
        S0, H0 = e0.S, e0.H

        # Candidate b1 values: union of surface P and W neighbors of b0.
        # This avoids O(N_edges^2) scans.  In GPU ELLPACK layout this is a small fixed loop.
        cand_b1 = set(surf.P.neigh[b0]) | set(W_surf.neigh[b0])
        for b1 in cand_b1:
            for i_e1 in edges_by_b.get(b1, []):
                e1 = edges[i_e1]
                a1 = e1.a
                S1 = e1.S
                n_pair += 1

                # T_H = -Tr[ P_A[a0,a1] S[a1,b1] P_B[b1,b0] H[a0,b0]^T ]
                PA_01 = tip.P.get(a0, a1)
                PB_10 = surf.P.get(b1, b0)
                if PA_01 is not None and PB_10 is not None:
                    EH += -trace4(PA_01, S1, PB_10, H0.T)

                # T_A = 1/2 Tr[ S[a0,b0] P_B[b0,b1] S[a1,b1]^T W_A[a1,a0] ]
                PB_01 = surf.P.get(b0, b1)
                WA_10 = W_tip.get(a1, a0)
                if PB_01 is not None and WA_10 is not None:
                    EWA += 0.5 * trace4(S0, PB_01, S1.T, WA_10)

                # T_B = 1/2 Tr[ S[a0,b0]^T P_A[a0,a1] S[a1,b1] W_B[b1,b0] ]
                WB_10 = W_surf.get(b1, b0)
                if PA_01 is not None and WB_10 is not None:
                    EWB += 0.5 * trace4(S0.T, PA_01, S1, WB_10)

    comps = {"EH": EH, "EWA": EWA, "EWB": EWB}
    return EH + EWA + EWB, comps, len(edges), n_pair


def evaluate_local_DA(tip: Fragment, surf: Fragment, tip_shift: np.ndarray,
                      cutoff: float, damp: float = 0.5) -> float:
    """Kernel 2: local Fukui / donor-acceptor correction.

    This is a local low-rank alternative to global HOMO/LUMO frontier MOs.
    For every active atom pair (a,b), couple donor channels on one side to acceptor
    channels on the other side.

    Nonorthogonal AO basis correction:
        T = u_D^T [ H_ab - 0.5*(eps_D+eps_A)*S_ab ] v_A

    This avoids spurious CT from pure overlap and is closer to generalized-eigenproblem
    perturbation theory than raw u^T H v.

    GPU target: yes.  This is embarrassingly parallel over active edges and channels.
    """
    # Organize channels by atom for fast lookup.
    tip_D: Dict[int, List[LocalChannel]] = {}
    tip_A: Dict[int, List[LocalChannel]] = {}
    surf_D: Dict[int, List[LocalChannel]] = {}
    surf_A: Dict[int, List[LocalChannel]] = {}
    for ch in tip.donors:    tip_D.setdefault(ch.atom, []).append(ch)
    for ch in tip.acceptors: tip_A.setdefault(ch.atom, []).append(ch)
    for ch in surf.donors:   surf_D.setdefault(ch.atom, []).append(ch)
    for ch in surf.acceptors: surf_A.setdefault(ch.atom, []).append(ch)

    edges, _ = build_active_edges(tip, surf, tip_shift, cutoff)
    E = 0.0
    for e in edges:
        # Tip donor -> surface acceptor
        for d in tip_D.get(e.a, []):
            for acc in surf_A.get(e.b, []):
                T = float(d.coeff @ (e.H - 0.5 * (d.eps + acc.eps) * e.S) @ acc.coeff)
                de = max(acc.eps - d.eps, 0.0)
                denom = math.sqrt(de * de + damp * damp)
                E += -d.weight * acc.weight * (T * T) / denom
        # Surface donor -> tip acceptor.  Use H_ba = H_ab.T, S_ba = S_ab.T.
        for d in surf_D.get(e.b, []):
            for acc in tip_A.get(e.a, []):
                T = float(d.coeff @ (e.H.T - 0.5 * (d.eps + acc.eps) * e.S.T) @ acc.coeff)
                de = max(acc.eps - d.eps, 0.0)
                denom = math.sqrt(de * de + damp * damp)
                E += -d.weight * acc.weight * (T * T) / denom
    return E


def make_symmetric_block(rng: np.random.Generator, scale: float, diag_boost: float = 0.0) -> Block:
    A = rng.normal(size=(NORB, NORB))
    M = scale * 0.5 * (A + A.T)
    M += diag_boost * np.eye(NORB)
    return M


def add_symmetric_block_pair(M: BlockSparse, i: int, j: int, block_ij: Block):
    M.set_block(i, j, block_ij)
    if i != j:
        M.set_block(j, i, block_ij.T)


def make_toy_fragment(name: str, positions: np.ndarray, rng: np.random.Generator,
                      density_cut: float = 4.5) -> Fragment:
    """Generate toy P/W sparse blocks and local donor/acceptor channels.

    Real precompute step should do:
      1. isolated DFTB diagonalization;
      2. P = sum_occ f_i c_i c_i^T;
      3. W = sum_occ f_i eps_i c_i c_i^T;
      4. truncate/taper P,W by atom distance;
      5. construct local donor/acceptor channels from Fukui/LDOS/projected orbitals.
    """
    n = len(positions)
    P = BlockSparse(n)
    W = BlockSparse(n)

    for i in range(n):
        for j in range(i, n):
            rij = np.linalg.norm(positions[i] - positions[j])
            if i == j or rij < density_cut:
                decay = math.exp(-rij / 2.0)
                if i == j:
                    # Diagonal occupied density-like block.  Values are not meant
                    # to be chemically exact; they only make the demo stable.
                    Pij = np.diag([1.5, 0.8, 0.8, 0.8]) + make_symmetric_block(rng, 0.03)
                    Wij = np.diag([-14.0, -8.0, -8.0, -8.0]) + make_symmetric_block(rng, 0.10)
                else:
                    Pij = rng.normal(scale=0.08 * decay, size=(NORB, NORB))
                    Wij = rng.normal(scale=0.25 * decay, size=(NORB, NORB))
                add_symmetric_block_pair(P, i, j, Pij)
                add_symmetric_block_pair(W, i, j, Wij)

    frag = Fragment(name=name, pos=np.asarray(positions, dtype=np.float64), P=P, W=W)

    # Toy local channels: use pz donor/acceptor on every atom, with different energies.
    # In real code, these should come from local Fukui/projected LDOS analysis.
    pz = np.array([0.0, 0.0, 0.0, 1.0])
    s = np.array([1.0, 0.0, 0.0, 0.0])
    for ia in range(n):
        frag.donors.append(LocalChannel(atom=ia, coeff=pz.copy(), eps=-6.0 + 0.2*rng.normal(), weight=0.6))
        frag.acceptors.append(LocalChannel(atom=ia, coeff=s.copy(), eps=-1.5 + 0.2*rng.normal(), weight=0.3))
    return frag


def shifted_W(W: BlockSparse, P: BlockSparse, lam: float) -> BlockSparse:
    """Return W' = W + lam*P with the same sparse support."""
    out = BlockSparse(W.n_atoms)
    keys = set(W.blocks.keys()) | set(P.blocks.keys())
    for key in keys:
        M = np.zeros((NORB, NORB), dtype=np.float64)
        if key in W.blocks:
            M += W.blocks[key]
        if key in P.blocks:
            M += lam * P.blocks[key]
        out.set_block(key[0], key[1], M)
    return out


def build_demo_system(seed: int = 123) -> Tuple[Fragment, Fragment]:
    rng = np.random.default_rng(seed)

    # Small tip: four atoms in a rough tetrahedral/CO-like cluster above the origin.
    tip_pos = np.array([
        [ 0.0,  0.0,  0.0],
        [ 0.0,  0.0, -1.2],
        [ 0.9,  0.0, -1.8],
        [-0.6,  0.7, -1.8],
    ], dtype=np.float64)

    # Surface: 7x7 square grid in z=0 plane.
    surf_pos = []
    for ix in range(-3, 4):
        for iy in range(-3, 4):
            surf_pos.append([1.6 * ix, 1.6 * iy, 0.0])
    surf_pos = np.array(surf_pos, dtype=np.float64)

    tip = make_toy_fragment("tip", tip_pos, rng, density_cut=4.0)
    surf = make_toy_fragment("surface", surf_pos, rng, density_cut=3.3)
    return tip, surf


# -----------------------------------------------------------------------------
# Optional OpenCL sketch.  This is intentionally not executable, because real code
# needs the exact Fireball orbital layout, SK table interpolation, and packed P/W
# buffers.  It shows what should be moved to GPU.
# -----------------------------------------------------------------------------
OPENCL_EDGEPAIR_SKETCH = r'''
// One workgroup per AFM pixel.
// Local memory:
//   Edge edges[MAX_EDGES];
//   int edges_by_B[MAX_LOCAL_B][MAX_EDGES_PER_B];
//
// GPU stages:
//   1) Build active A-B edge list using surface cell list.
//   2) Evaluate S_ab,H_ab,dS_ab,dH_ab from SK tables.
//   3) Edge-pair contraction for T_H,T_A,T_B.
//   4) Optional local donor/acceptor edge loop.
//   5) Workgroup reduction to E,Fx,Fy,Fz.
//
// Hot data layout should be ELLPACK-like, not CSR/dicts:
//   int neighB[Nsurf*MAX_PW_NEIGH];
//   float PB[Nsurf*MAX_PW_NEIGH*NORB*NORB];
//   float WB[Nsurf*MAX_PW_NEIGH*NORB*NORB];
//
// Inner edge-pair loop, schematic only:
__kernel void eval_afm_pixel_edgepair(
    __global const float4* tip_pos,
    __global const float4* surf_pos,
    __global const int* neighB,
    __global const float* PB,
    __global const float* WB,
    __constant const float* PA,
    __constant const float* WA,
    __global float4* out_EF
){
    int pix = get_group_id(0);
    int lid = get_local_id(0);

    // build edges[] in local memory ...
    // barrier(CLK_LOCAL_MEM_FENCE);

    float E = 0.0f;
    for(int ie0=lid; ie0<nEdges; ie0+=get_local_size(0)){
        Edge e0 = edges[ie0];
        for(int kb=0; kb<MAX_PW_NEIGH; kb++){
            int b1 = neighB[e0.b*MAX_PW_NEIGH + kb];
            if(b1<0) continue;
            for(int kk=0; kk<nEdgesByB[b1]; kk++){
                Edge e1 = edges[ edgesByB[b1][kk] ];
                // load PA[a0,a1], WA[a1,a0], PB[b1,b0], WB[b1,b0]
                // E += -trace4(PA, e1.S, PB10, transpose(e0.H));
                // E += 0.5*trace4(e0.S, PB01, transpose(e1.S), WA);
                // E += 0.5*trace4(transpose(e0.S), PA, e1.S, WB);
            }
        }
    }
    // reduce E within workgroup; write out_EF[pix]
}
'''


def main():
    tip, surf = build_demo_system()
    cutoff = 4.2

    # One AFM pixel: tip apex at z=3.0 A above surface, with slight xy offset.
    tip_shift = np.array([0.15, -0.10, 3.0], dtype=np.float64)

    print("=== Single-pixel energy demo ===")
    E0 = evaluate_pauli_pairdiag(tip, surf, tip_shift, cutoff)
    E1, comps, n_edges, n_pairs = evaluate_pauli_edgepair(tip, surf, tip_shift, cutoff)
    Eda = evaluate_local_DA(tip, surf, tip_shift, cutoff)
    print(f"active edges      : {n_edges}")
    print(f"edge-pair visits  : {n_pairs}")
    print(f"pair-diag Pauli   : {E0: .8e}")
    print(f"edge-pair Pauli   : {E1: .8e}  components={comps}")
    print(f"local DA/Fukui    : {Eda: .8e}")
    print(f"total short-range : {E1 + Eda: .8e}")

    # Energy-zero shift invariance test.
    # H_AB -> H_AB + lam*S_AB and W_X -> W_X + lam*P_X should leave E_PW invariant.
    lam = 3.7
    Wt_shift = shifted_W(tip.W, tip.P, lam)
    Ws_shift = shifted_W(surf.W, surf.P, lam)
    E_shift, comps_shift, _, _ = evaluate_pauli_edgepair(
        tip, surf, tip_shift, cutoff, h_shift=lam, W_tip=Wt_shift, W_surf=Ws_shift
    )
    print("\n=== Energy-zero shift test ===")
    print(f"lambda shift      : {lam}")
    print(f"E original        : {E1: .12e}")
    print(f"E shifted         : {E_shift: .12e}")
    print(f"difference        : {E_shift - E1: .3e}")
    print("A small difference means the P/W/H/S contraction and factors are self-consistent.")

    # Tiny scan timing.  This is Python, so do not interpret as production speed.
    # It is useful only for checking scaling with n_edges and n_pairs.
    print("\n=== Tiny 2D scan timing, Python only ===")
    xs = np.linspace(-1.2, 1.2, 9)
    ys = np.linspace(-1.2, 1.2, 9)
    t0 = time.perf_counter()
    vals = []
    max_edges = max_pairs = 0
    for x in xs:
        for y in ys:
            E, _, ne, npair = evaluate_pauli_edgepair(tip, surf, np.array([x, y, 3.0]), cutoff)
            vals.append(E)
            max_edges = max(max_edges, ne)
            max_pairs = max(max_pairs, npair)
    dt = time.perf_counter() - t0
    vals = np.array(vals)
    print(f"pixels            : {len(vals)}")
    print(f"time              : {dt:.3f} s  (NumPy/Python loops, not optimized)")
    print(f"E range           : [{vals.min(): .6e}, {vals.max(): .6e}]")
    print(f"max edges/pixel   : {max_edges}")
    print(f"max edgepairs/pix : {max_pairs}")

    print("\nGPU migration recommendation:")
    print("  CPU/precompute : isolated-fragment DFTB, P/W construction, sparse packing, local Fukui channels")
    print("  GPU/per-pixel  : neighbor search, SK S/H/dS/dH, active edge list, edge-pair Pauli, local DA, reductions")


if __name__ == "__main__":
    main()
