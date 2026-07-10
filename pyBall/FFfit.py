"""Python ctypes wrapper for C++ FFfit force-field parameter fitting library.

Wraps FFfit_lib.cpp → FFfit.h for fitting bond/angle/dihedral stiffness parameters
to reference Hessians via linear least-squares or gradient descent. Also provides
optimized C++ graph algorithms (CSR bond-graph BFS, 1-4 pair detection, dihedral
enumeration) and batch dihedral sensitivity computation.

THEORY OVERVIEW
===============

Given a force field E(r) = Σ_t k_t * f_t(q_t(r)), the Cartesian Hessian is:

    H_αβ = ∂²E/∂r_α ∂r_β = Σ_t k_t * A_t

where A_t = ∂²f_t/∂r_α ∂r_β is the SENSITIVITY MATRIX, linear in k_t at fixed
geometry.  By the chain rule:

    A_t = f''(q) * b⊗b^T + f'(q) * C

where b = ∂q/∂r is the WILSON VECTOR and C = ∂²q/∂r² is the COORDINATE HESSIAN.
The first term (rank-one) dominates at equilibrium; the second (prestress) is
nonzero when geometry ≠ equilibrium.

This module provides three complementary fitting objectives combined into a
single linear least-squares problem (the "hybrid" fit):

1. MODE objective: targets Rayleigh quotients v_s^T D(k) v_s = λ_s for each
   reference vibrational mode, where D = M^{-1/2} H M^{-1/2}.
2. LOCAL objective: fits only graph-local mass-weighted Hessian blocks,
   masked by bounded BFS on the bond graph (local_hessian_mask).
3. INTERNAL objective: projects H into Wilson internal coordinates via
   F = C^{+T} D C^+ (the GF method), fitting the valence force matrix.

All three are linear in k, so the combined problem is a single linear solve.
The direct stacked residual matrix is solved via scipy.optimize.lsq_linear
with column scaling, Tikhonov regularization, and non-negative bounds.

GRAPH ALGORITHMS
================

Bond-graph topology queries use a CSR (Compressed Sparse Row) adjacency built
once per call, with bounded ring-based BFS that only visits atoms within the
requested distance cutoff. This avoids computing the full N×N distance matrix
when only local distances (≤3 hops) are needed. A combined pass
(local_mask_and_14pairs) extracts both the Hessian mask and 1-4 neighbor pairs
in a single BFS traversal per atom.

Batch dihedral sensitivity (dihedral_dHdk_batch_typed_cpp) replaces the Python
compute_dihedral_sensitivity loop — all dihedrals are processed in one C++ call
with in-place FD perturbation and Hessian symmetry exploitation.

See doc/Topics/FFfit/HessianFitting_Theory.md for the full derivation.
Parity tests: tests/tSiNCs/test_parity_graph_cpp.py (14 tests).

Usage:
    from pyBall import FFfit
    fitter = FFfit.FFfit()
    fitter.set_geometry(positions)
    fitter.add_bond(0, 1, 1.48)
    fitter.set_symbols(['Si', 'H', 'H', 'H', 'H'])
    fitter.auto_assign_types()
    k = fitter.fit_lsq(H_ref)
    H_model = fitter.compute_model_hessian()
"""

import numpy as np
import ctypes
from collections import deque
from ctypes import c_int, c_double, c_char_p, c_void_p, c_bool, byref, POINTER
from . import cpp_utils_

c_double_p = ctypes.POINTER(c_double)
c_int_p    = ctypes.POINTER(c_int)
c_char_pp  = ctypes.POINTER(c_char_p)

cpp_utils_.BUILD_PATH = cpp_utils_.BUILD_PATH  # already set by other modules, but ensure
lib = cpp_utils_.loadLib('FFfit_lib', recompile=False)

array1d = np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
array1i = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=1, flags='CONTIGUOUS')

# === Setup argtypes/restype ===

lib.fffit_create.restype  = c_void_p
lib.fffit_create.argtypes = []

lib.fffit_destroy.restype  = None
lib.fffit_destroy.argtypes = [c_void_p]

lib.fffit_set_geometry.restype  = None
lib.fffit_set_geometry.argtypes = [c_void_p, c_int, array1d]

lib.fffit_add_bond.restype  = None
lib.fffit_add_bond.argtypes = [c_void_p, c_int, c_int, c_double]

lib.fffit_add_angle.restype  = None
lib.fffit_add_angle.argtypes = [c_void_p, c_int, c_int, c_int, c_double]

lib.fffit_set_symbols.restype  = None
lib.fffit_set_symbols.argtypes = [c_void_p, c_int, c_char_pp]

lib.fffit_auto_assign_types.restype  = None
lib.fffit_auto_assign_types.argtypes = [c_void_p]

lib.fffit_set_n_free.restype  = None
lib.fffit_set_n_free.argtypes = [c_void_p, c_int]

lib.fffit_set_bond_param.restype  = None
lib.fffit_set_bond_param.argtypes = [c_void_p, c_int, c_int]

lib.fffit_set_angle_param.restype  = None
lib.fffit_set_angle_param.argtypes = [c_void_p, c_int, c_int]

lib.fffit_get_n_free.restype  = c_int
lib.fffit_get_n_free.argtypes = [c_void_p]

lib.fffit_get_natoms.restype  = c_int
lib.fffit_get_natoms.argtypes = [c_void_p]

lib.fffit_get_nbonds.restype  = c_int
lib.fffit_get_nbonds.argtypes = [c_void_p]

lib.fffit_get_nangles.restype  = c_int
lib.fffit_get_nangles.argtypes = [c_void_p]

lib.fffit_fit_lsq.restype  = None
lib.fffit_fit_lsq.argtypes = [c_void_p, array1d, array1d, array1d, array1d]

lib.fffit_fit_gd.restype  = None
lib.fffit_fit_gd.argtypes = [c_void_p, array1d, array1d, array1d, c_double, c_double, c_int, c_double, c_int, array1d]

lib.fffit_compute_model_hessian.restype  = None
lib.fffit_compute_model_hessian.argtypes = [c_void_p, array1d]

lib.fffit_get_params.restype  = None
lib.fffit_get_params.argtypes = [c_void_p, array1d]

lib.fffit_set_params.restype  = None
lib.fffit_set_params.argtypes = [c_void_p, array1d, c_int]

lib.fffit_get_bond_param_idx.restype  = None
lib.fffit_get_bond_param_idx.argtypes = [c_void_p, array1i]

lib.fffit_get_angle_param_idx.restype  = None
lib.fffit_get_angle_param_idx.argtypes = [c_void_p, array1i]

lib.fffit_print_params.restype  = None
lib.fffit_print_params.argtypes = [c_void_p]

lib.fffit_accumulate_normal_equations.restype  = None
lib.fffit_accumulate_normal_equations.argtypes = [c_void_p, array1d, array1d, array1d, array1d, array1d]

lib.fffit_compute_gradient_loss.restype  = None
lib.fffit_compute_gradient_loss.argtypes = [c_void_p, array1d, array1d, array1d, array1d, array1d, POINTER(c_double)]

lib.fffit_solve_normal_equations.restype  = None
lib.fffit_solve_normal_equations.argtypes = [array1d, array1d, array1d, c_int]

# === Graph algorithms ===
lib.fffit_bond_graph_distances.restype  = None
lib.fffit_bond_graph_distances.argtypes = [array1i, c_int, c_int, array1i]

lib.fffit_local_hessian_mask.restype  = None
lib.fffit_local_hessian_mask.argtypes = [array1i, c_int, c_int, c_int, np.ctypeslib.ndpointer(dtype=np.uint8, ndim=1, flags='CONTIGUOUS')]

lib.fffit_find_3rd_neighbor_bonds.restype  = None
lib.fffit_find_3rd_neighbor_bonds.argtypes = [array1i, c_int, c_int, array1d, c_double, array1i, POINTER(c_int)]

lib.fffit_enumerate_dihedrals.restype  = None
lib.fffit_enumerate_dihedrals.argtypes = [array1i, c_int, c_int, array1i, POINTER(c_int)]

# === Wilson B matrix ===
lib.fffit_build_wilson_matrix.restype  = None
lib.fffit_build_wilson_matrix.argtypes = [c_void_p, array1d]

# === Dihedral sensitivity ===
lib.fffit_dihedral_hessian.restype  = None
lib.fffit_dihedral_hessian.argtypes = [c_void_p, c_int, c_int, c_int, c_int, array1d]

lib.fffit_dihedral_dHdk.restype  = None
lib.fffit_dihedral_dHdk.argtypes = [c_void_p, c_int, c_int, c_int, c_int, array1d]

lib.fffit_local_mask_and_14pairs.restype  = None
lib.fffit_local_mask_and_14pairs.argtypes = [array1i, c_int, c_int, c_int, np.ctypeslib.ndpointer(dtype=np.uint8, ndim=1, flags='CONTIGUOUS'), array1d, c_double, array1i, POINTER(c_int)]

lib.fffit_dihedral_dHdk_batch.restype  = None
lib.fffit_dihedral_dHdk_batch.argtypes = [c_void_p, array1i, c_int, array1d]

lib.fffit_dihedral_dHdk_batch_typed.restype  = None
lib.fffit_dihedral_dHdk_batch_typed.argtypes = [c_void_p, array1i, c_int, array1i, c_int, array1d]


class FFfit:
    """Python wrapper for C++ FFfit class.

    Fits bond/angle force-field stiffness parameters to reference Hessians
    using linear least-squares or gradient descent with sparse parameter mapping
    for symmetry constraints.
    """
    def __init__(self):
        self._handle = lib.fffit_create()

    def __del__(self):
        if hasattr(self, '_handle') and self._handle:
            lib.fffit_destroy(self._handle)

    def set_geometry(self, positions):
        """Set atomic positions. positions: (natoms, 3) array in Angstrom."""
        positions = np.ascontiguousarray(positions, dtype=np.float64).ravel()
        lib.fffit_set_geometry(self._handle, len(positions) // 3, positions)

    def add_bond(self, i, j, l0):
        """Add a bond between atoms i,j with equilibrium length l0."""
        lib.fffit_add_bond(self._handle, i, j, l0)

    def add_angle(self, i, j, k, theta0):
        """Add an angle i-j-k (j central) with equilibrium angle theta0 (radians)."""
        lib.fffit_add_angle(self._handle, i, j, k, theta0)

    def set_symbols(self, symbols):
        """Set element symbols for auto type assignment. symbols: list of strings."""
        encoded = [s.encode('utf-8') for s in symbols]
        arr = (c_char_p * len(encoded))(*encoded)
        lib.fffit_set_symbols(self._handle, len(symbols), arr)

    def auto_assign_types(self):
        """Auto-assign parameter indices based on element pairs/triples."""
        lib.fffit_auto_assign_types(self._handle)

    def set_n_free(self, n_free):
        """Explicitly set number of free parameters (for multi-system global fitting)."""
        lib.fffit_set_n_free(self._handle, n_free)

    def set_bond_param(self, ibond, param_idx):
        """Manually set bond parameter index (for custom symmetry schemes)."""
        lib.fffit_set_bond_param(self._handle, ibond, param_idx)

    def set_angle_param(self, iangle, param_idx):
        """Manually set angle parameter index."""
        lib.fffit_set_angle_param(self._handle, iangle, param_idx)

    @property
    def n_free(self):
        return lib.fffit_get_n_free(self._handle)

    @property
    def natoms(self):
        return lib.fffit_get_natoms(self._handle)

    @property
    def nbonds(self):
        return lib.fffit_get_nbonds(self._handle)

    @property
    def nangles(self):
        return lib.fffit_get_nangles(self._handle)

    def fit_lsq(self, H_ref, H_0=None, weight=None):
        """Fit via linear least-squares (normal equations + Gaussian elimination).

        H_ref: (3N x 3N) reference Hessian (flattened or 2D, row-major)
        H_0:   constant part of Hessian (optional, same shape)
        weight: per-element weight (optional, flattened 3N*3N)
        Returns: fitted parameters as (n_free,) array
        """
        H_ref = np.ascontiguousarray(H_ref, dtype=np.float64).ravel()
        H_0_arr = np.zeros_like(H_ref) if H_0 is None else np.ascontiguousarray(H_0, dtype=np.float64).ravel()
        w_arr = np.ones_like(H_ref) if weight is None else np.ascontiguousarray(weight, dtype=np.float64).ravel()
        k = np.zeros(self.n_free, dtype=np.float64)
        lib.fffit_fit_lsq(self._handle, H_ref, H_0_arr, w_arr, k)
        return k

    def fit_gd(self, H_ref, H_0=None, weight=None, lr=1e-4, momentum=0.9,
               max_iter=5000, tol=1e-10, verbose=True):
        """Fit via gradient descent with momentum.

        H_ref: (3N x 3N) reference Hessian
        Returns: fitted parameters as (n_free,) array
        """
        H_ref = np.ascontiguousarray(H_ref, dtype=np.float64).ravel()
        H_0_arr = np.zeros_like(H_ref) if H_0 is None else np.ascontiguousarray(H_0, dtype=np.float64).ravel()
        w_arr = np.ones_like(H_ref) if weight is None else np.ascontiguousarray(weight, dtype=np.float64).ravel()
        k = np.zeros(self.n_free, dtype=np.float64)
        lib.fffit_fit_gd(self._handle, H_ref, H_0_arr, w_arr, lr, momentum, max_iter, tol, int(verbose), k)
        return k

    def compute_model_hessian(self):
        """Compute model Hessian from current fitted parameters.
        Returns: (3N x 3N) array
        """
        n3 = self.natoms * 3
        H = np.zeros(n3 * n3, dtype=np.float64)
        lib.fffit_compute_model_hessian(self._handle, H)
        return H.reshape(n3, n3)

    def get_params(self):
        """Get current fitted parameters as (n_free,) array."""
        k = np.zeros(self.n_free, dtype=np.float64)
        lib.fffit_get_params(self._handle, k)
        return k

    def set_params(self, k):
        """Set parameters directly (for computing model hessian with externally fitted params)."""
        k_arr = np.ascontiguousarray(k, dtype=np.float64)
        lib.fffit_set_params(self._handle, k_arr, len(k_arr))

    def get_bond_param_idx(self):
        """Get bond→parameter index mapping as (nbonds,) array."""
        idx = np.zeros(self.nbonds, dtype=np.int32)
        lib.fffit_get_bond_param_idx(self._handle, idx)
        return idx

    def get_angle_param_idx(self):
        """Get angle→parameter index mapping as (nangles,) array."""
        idx = np.zeros(self.nangles, dtype=np.int32)
        lib.fffit_get_angle_param_idx(self._handle, idx)
        return idx

    def print_params(self):
        """Print fitted parameters to stdout."""
        lib.fffit_print_params(self._handle)

    # === Multi-system accumulation (for global fitting across multiple molecules) ===

    def accumulate_normal_equations(self, G, y, H_ref, H_0=None, weight=None):
        """Accumulate this system's contribution to global normal equations G, y (in-place).

        G: (n_free, n_free) array, y: (n_free,) array — both modified in-place.
        H_ref: (3N x 3N) reference Hessian for this system.
        """
        np_ = self.n_free
        H_ref_arr = np.ascontiguousarray(H_ref, dtype=np.float64).ravel()
        H_0_arr = np.zeros_like(H_ref_arr) if H_0 is None else np.ascontiguousarray(H_0, dtype=np.float64).ravel()
        w_arr = np.ones_like(H_ref_arr) if weight is None else np.ascontiguousarray(weight, dtype=np.float64).ravel()
        G_flat = np.ascontiguousarray(G, dtype=np.float64).ravel()
        y_arr = np.ascontiguousarray(y, dtype=np.float64)
        lib.fffit_accumulate_normal_equations(self._handle, G_flat, y_arr, H_ref_arr, H_0_arr, w_arr)
        G[:] = G_flat.reshape(np_, np_)
        y[:] = y_arr

    def compute_gradient_loss(self, H_ref, k, H_0=None, weight=None, grad_out=None, loss_out=None):
        """Compute gradient and loss for this system with given k. Accumulates into provided arrays.

        H_ref: (3N x 3N) reference Hessian.
        k: (n_free,) current parameters.
        grad_out: (n_free,) array to accumulate gradient into (modified in-place).
        loss_out: list or [float] to accumulate loss into.
        """
        H_ref_arr = np.ascontiguousarray(H_ref, dtype=np.float64).ravel()
        H_0_arr = np.zeros_like(H_ref_arr) if H_0 is None else np.ascontiguousarray(H_0, dtype=np.float64).ravel()
        w_arr = np.ones_like(H_ref_arr) if weight is None else np.ascontiguousarray(weight, dtype=np.float64).ravel()
        k_arr = np.ascontiguousarray(k, dtype=np.float64)
        if grad_out is None:
            grad_out = np.zeros(self.n_free, dtype=np.float64)
        grad_arr = np.ascontiguousarray(grad_out, dtype=np.float64)
        loss = c_double(0.0)
        lib.fffit_compute_gradient_loss(self._handle, H_ref_arr, H_0_arr, w_arr, k_arr, grad_arr, byref(loss))
        grad_out[:] = grad_arr
        if loss_out is not None:
            loss_out[0] += loss.value
        return grad_out, loss.value

    @staticmethod
    def solve_normal_equations(G, y):
        """Solve G k = y via Gaussian elimination. Returns k."""
        np_ = len(y)
        G_flat = np.ascontiguousarray(G, dtype=np.float64).ravel().copy()  # copy: G is modified in-place
        y_arr = np.ascontiguousarray(y, dtype=np.float64).copy()
        k = np.zeros(np_, dtype=np.float64)
        lib.fffit_solve_normal_equations(G_flat, y_arr, k, np_)
        return k

    # === Graph algorithms (C++ accelerated) ===

    @staticmethod
    def bond_graph_distances(bond_pairs, natoms):
        """BFS shortest-path distances between all pairs in the bond graph.
        bond_pairs: (nbonds, 2) or (nbonds*2,) int array.
        Returns: (natoms, natoms) int array, -1 = unreachable.
        """
        bp = np.ascontiguousarray(bond_pairs, dtype=np.int32).ravel()
        dist = np.full(natoms * natoms, -1, dtype=np.int32)
        lib.fffit_bond_graph_distances(bp, len(bp) // 2, natoms, dist)
        return dist.reshape(natoms, natoms)

    @staticmethod
    def local_hessian_mask(bonds, natoms, max_graph_distance=2):
        """Boolean mask for atom pairs within max_graph_distance in bond graph.
        bonds: (nbonds, 2) or list of (i,j) pairs.
        Returns: (3*natoms, 3*natoms) boolean array.
        """
        bp = np.ascontiguousarray([(b[0], b[1]) for b in bonds], dtype=np.int32).ravel()
        mask = np.zeros(natoms * natoms, dtype=np.uint8)
        lib.fffit_local_hessian_mask(bp, len(bp) // 2, natoms, max_graph_distance, mask)
        block = mask.reshape(natoms, natoms).astype(bool)
        return np.repeat(np.repeat(block, 3, axis=0), 3, axis=1)

    @staticmethod
    def find_3rd_neighbor_bonds(bond_pairs, natoms, positions, max_dist=0.0):
        """Find 1-4 pairs (graph distance 3). Returns list of (i, j) pairs.
        positions: (natoms, 3) array. max_dist<=0 means no distance filter.
        """
        bp = np.ascontiguousarray(bond_pairs, dtype=np.int32).ravel()
        pos = np.ascontiguousarray(positions, dtype=np.float64).ravel()
        pairs3 = np.zeros(natoms * natoms * 2, dtype=np.int32)  # over-allocate
        n3 = c_int(0)
        lib.fffit_find_3rd_neighbor_bonds(bp, len(bp) // 2, natoms, pos, max_dist, pairs3, byref(n3))
        return [(pairs3[i*2], pairs3[i*2+1]) for i in range(n3.value)]

    @staticmethod
    def local_mask_and_14pairs(bonds, natoms, max_graph_distance, positions, max_dist=0.0):
        """Combined local Hessian mask + 1-4 pairs in a single BFS pass.
        Returns: (mask_3Nx3N_bool, list_of_14_pairs).
        """
        bp = np.ascontiguousarray([(b[0], b[1]) for b in bonds], dtype=np.int32).ravel()
        pos = np.ascontiguousarray(positions, dtype=np.float64).ravel()
        mask = np.zeros(natoms * natoms, dtype=np.uint8)
        pairs3 = np.zeros(natoms * natoms * 2, dtype=np.int32)
        n3 = c_int(0)
        lib.fffit_local_mask_and_14pairs(bp, len(bp) // 2, natoms, max_graph_distance,
                                          mask, pos, max_dist, pairs3, byref(n3))
        block = mask.reshape(natoms, natoms).astype(bool)
        mask_3N = np.repeat(np.repeat(block, 3, axis=0), 3, axis=1)
        pairs = [(pairs3[i*2], pairs3[i*2+1]) for i in range(n3.value)]
        return mask_3N, pairs

    @staticmethod
    def enumerate_dihedrals(bond_pairs, natoms):
        """Enumerate proper torsions (i-j-k-l) from bond topology.
        bond_pairs: (nbonds, 2) or (nbonds*2,) int array.
        Returns: list of (i, j, k, l) tuples.
        """
        bp = np.ascontiguousarray(bond_pairs, dtype=np.int32).ravel()
        # Max dihedrals: each bond can produce up to (deg-1)^2 dihedrals; over-allocate
        max_dihed = nbonds = len(bp) // 2
        for _ in range(4): max_dihed = max_dihed * max(natoms, 1)
        diheds = np.zeros(max_dihed * 4, dtype=np.int32)
        n_out = c_int(0)
        lib.fffit_enumerate_dihedrals(bp, nbonds, natoms, diheds, byref(n_out))
        return [(diheds[i*4], diheds[i*4+1], diheds[i*4+2], diheds[i*4+3]) for i in range(n_out.value)]

    def build_wilson_matrix_cpp(self):
        """Build Wilson B matrix using C++ implementation.
        Must have geometry, bonds, and angles set on this FFfit instance.
        Returns: (n_internal, 3N) array where n_internal = nbonds + nangles.
        """
        n3 = self.natoms * 3
        n_int = self.nbonds + self.nangles
        B = np.zeros(n_int * n3, dtype=np.float64)
        lib.fffit_build_wilson_matrix(self._handle, B)
        return B.reshape(n_int, n3)

    def dihedral_hessian_cpp(self, i, j, k, l):
        """Compute 12x12 dihedral Hessian via C++ finite differences.
        Returns: (12, 12) array.
        """
        H = np.zeros(144, dtype=np.float64)
        lib.fffit_dihedral_hessian(self._handle, i, j, k, l, H)
        return H.reshape(12, 12)

    def dihedral_dHdk_cpp(self, i, j, k, l):
        """Compute dihedral sensitivity dH/dV scattered into (3N, 3N).
        Returns: (3N, 3N) array.
        """
        n3 = self.natoms * 3
        A = np.zeros(n3 * n3, dtype=np.float64)
        lib.fffit_dihedral_dHdk(self._handle, i, j, k, l, A)
        return A.reshape(n3, n3)

    def dihedral_dHdk_batch_cpp(self, dihedrals):
        """Batch: compute dH/dV for all dihedrals, accumulated into one (3N, 3N) matrix.
        dihedrals: list of (i, j, k, l) tuples or (ndihedrals*4) int array.
        Returns: (3N, 3N) array.
        """
        n3 = self.natoms * 3
        if len(dihedrals) == 0:
            return np.zeros((n3, n3), dtype=np.float64)
        dh = np.ascontiguousarray(dihedrals, dtype=np.int32).ravel()
        A = np.zeros(n3 * n3, dtype=np.float64)
        lib.fffit_dihedral_dHdk_batch(self._handle, dh, len(dh) // 4, A)
        return A.reshape(n3, n3)

    def dihedral_dHdk_batch_typed_cpp(self, dihedrals, type_idx, n_types):
        """Batch with per-type accumulation: each dihedral's Hessian goes into its
        parameter's (3N, 3N) matrix. Equivalent to Python compute_dihedral_sensitivity.
        dihedrals: list of (i,j,k,l) tuples. type_idx: list of param indices per dihedral.
        Returns: (n_types, 3N, 3N) array.
        """
        n3 = self.natoms * 3
        if len(dihedrals) == 0 or n_types == 0:
            return np.zeros((n_types, n3, n3), dtype=np.float64)
        dh = np.ascontiguousarray(dihedrals, dtype=np.int32).ravel()
        ti = np.ascontiguousarray(type_idx, dtype=np.int32)
        A = np.zeros(n_types * n3 * n3, dtype=np.float64)
        lib.fffit_dihedral_dHdk_batch_typed(self._handle, dh, len(dh) // 4, ti, n_types, A)
        return A.reshape(n_types, n3, n3)


def collect_sensitivity_matrices(fitter, extra=None):
    """Return A[p] = dH/dk_p from a linear FFfit model without finite differences."""
    n = fitter.natoms * 3
    npar = fitter.n_free
    old = fitter.get_params()
    A = np.empty((npar, n, n), dtype=np.float64)
    for p in range(npar):
        k = np.zeros(npar)
        k[p] = 1.0
        fitter.set_params(k)
        A[p] = fitter.compute_model_hessian()
    fitter.set_params(old)
    if extra:
        for p, Ap in extra.items():
            if p < 0 or p >= npar:
                raise IndexError(f"extra sensitivity parameter {p} outside [0,{npar})")
            Ap = np.asarray(Ap, dtype=np.float64)
            if Ap.shape != (n, n):
                raise ValueError(f"extra sensitivity {p} has shape {Ap.shape}, expected {(n, n)}")
            A[p] += Ap
    return A


def mass_weight_hessian(H, masses):
    """Transform Cartesian Hessian(s) to D = M^-1/2 H M^-1/2.

    The mass-weighted dynamical matrix D relates to vibrational frequencies by:
        D v_s = λ_s v_s,  where λ_s = ω_s² (eigenvalues are squared frequencies)
    The mass-weighted displacement is u_s = M^{-1/2} v_s, so that:
        v_s = M^{1/2} u_s  (Cartesian ← mass-weighted)
    This transformation makes the kinetic energy diagonal: T = ½ Σ ṁ_i u̇_i².
    """
    H = np.asarray(H, dtype=np.float64)
    masses = np.asarray(masses, dtype=np.float64)
    if masses.ndim != 1 or np.any(~np.isfinite(masses)) or np.any(masses <= 0.0):
        raise ValueError("masses must be a finite positive 1D array")
    n = masses.size * 3
    if H.shape[-2:] != (n, n):
        raise ValueError(f"Hessian trailing shape {H.shape[-2:]} does not match {n} Cartesian coordinates")
    im = 1.0 / np.sqrt(np.repeat(masses, 3))
    lead = (1,) * (H.ndim - 2)
    return H * im.reshape(lead + (n, 1)) * im.reshape(lead + (1, n))


def rigid_and_vibrational_bases(positions, masses, rtol=1e-10):
    """Return orthonormal mass-weighted rigid and complementary vibrational bases.

    PHYSICAL BASIS:
    Rigid-body motion (translation + rotation) produces zero Hessian eigenvalues.
    These 6 (or 5 for linear molecules) null modes must be removed before fitting
    vibrational frequencies, otherwise they dominate the least-squares objective
    with uninformative zero targets.

    CONSTRUCTION:
    In mass-weighted coordinates, the 3 translation vectors are:
        T_α = sqrt(m_i) * e_α  (α = x,y,z; repeated for each atom i)
    The 3 rotation vectors are (about the center of mass):
        R_α = sqrt(m_i) * (e_α × r_i^CM)  (cross product with COM-positioned r_i)
    These 6 vectors are assembled into a (3N × 6) matrix and orthonormalized by SVD.
    The rank of this matrix reveals the number of rigid modes:
        - 6 for a nonlinear molecule (3 trans + 3 rot)
        - 5 for a linear molecule (3 trans + 2 rot; rotation about the axis is zero)
    The SVD left-singular vectors U[:, :rank] span the rigid space; U[:, rank:]
    spans the complementary vibrational space.
    """
    pos = np.asarray(positions, dtype=np.float64)
    masses = np.asarray(masses, dtype=np.float64)
    if pos.shape != (masses.size, 3):
        raise ValueError(f"positions shape {pos.shape} does not match masses shape {masses.shape}")
    if np.any(~np.isfinite(pos)):
        raise ValueError("positions contain NaN or infinity")
    r = pos - np.sum(pos * masses[:, None], axis=0) / np.sum(masses)
    sm = np.sqrt(masses)
    R = np.zeros((3 * masses.size, 6))
    for a in range(3):
        R[a::3, a] = sm
        for i in range(masses.size):
            R[3*i:3*i+3, 3+a] = sm[i] * np.cross(np.eye(3)[a], r[i])
    U, s, _ = np.linalg.svd(R, full_matrices=True)
    tol = rtol * max(R.shape) * s[0]
    rank = int(np.sum(s > tol))
    return U[:, :rank], U[:, rank:]


def reference_vibrational_modes(H_ref, positions, masses):
    """Diagonalize the reference dynamical matrix in the complete non-rigid space.

    The QM dynamical matrix D = M^{-1/2} H M^{-1/2} is projected onto the
    vibrational subspace (rigid modes removed) and diagonalized:
        Q^T D Q C = C Λ
    where Q spans the vibrational space and C are the eigenvectors within it.
    The full-space modes are V = Q C, with eigenvalues λ_s = ω_s².
    Frequencies in cm⁻¹: ν_s = sqrt(|λ_s|) × 521.5.
    """
    D = mass_weight_hessian(H_ref, masses)
    D = 0.5 * (D + D.T)
    _, Q = rigid_and_vibrational_bases(positions, masses)
    lam, C = np.linalg.eigh(Q.T @ D @ Q)
    return Q @ C, lam


def build_wilson_matrix(positions, bonds, angles):
    """Build B=dq/dx for bond lengths and angles in Angstrom and radians.

    The Wilson B matrix (Wilson, Decius, Cross, *Molecular Vibrations*, 1955)
    relates internal coordinate displacements to Cartesian displacements:
        Δq = B Δx
    where q = (r_ij, θ_ijk, ...) are internal coordinates and x = (r_1, ..., r_N)
    are Cartesian positions.  Each row of B is the Wilson vector b_t = ∂q_t/∂x.

    BOND row (i,j):  b = ∂r/∂x, where r = |r_j - r_i|
        b_i = -e,  b_j = +e  (e = unit bond vector),  all other atoms = 0

    ANGLE row (i,j,k):  b = ∂θ/∂x, where θ = arccos(u·v), u = (r_i-r_j)/|r_i-r_j|, v = (r_k-r_j)/|r_k-r_j|
        ∂θ/∂r_i = -(v - ct*u) / (|a| * sin θ)  (from ∂(cos θ)/∂r_i = (v-ct*u)/|a| and ∂θ = -∂(cos θ)/sin θ)
        ∂θ/∂r_k = -(u - ct*v) / (|c| * sin θ)  (by symmetry i↔k)
        ∂θ/∂r_j = -(∂θ/∂r_i + ∂θ/∂r_k)  (translational invariance)

    The B matrix is central to the Wilson GF method: G = B M^{-1} B^T is the
    kinematic matrix, and the vibrational problem becomes G F P = λ P in
    internal coordinates, where F is the valence force constant matrix.
    """
    pos = np.asarray(positions, dtype=np.float64)
    n = pos.shape[0] * 3
    B = np.zeros((len(bonds) + len(angles), n), dtype=np.float64)
    labels = []
    iq = 0
    for bond in bonds:
        i, j = int(bond[0]), int(bond[1])
        d = pos[j] - pos[i]
        r = np.linalg.norm(d)
        if r <= 1e-12:
            raise ValueError(f"zero-length bond ({i},{j})")
        e = d / r
        B[iq, 3*i:3*i+3] = -e
        B[iq, 3*j:3*j+3] = e
        labels.append(("bond", i, j))
        iq += 1
    for angle in angles:
        i, j, k = int(angle[0]), int(angle[1]), int(angle[2])
        a = pos[i] - pos[j]
        c = pos[k] - pos[j]
        ra, rc = np.linalg.norm(a), np.linalg.norm(c)
        if ra <= 1e-12 or rc <= 1e-12:
            raise ValueError(f"zero-length arm in angle ({i},{j},{k})")
        u, v = a / ra, c / rc
        ct = np.clip(np.dot(u, v), -1.0, 1.0)
        st = np.sqrt(max(0.0, 1.0 - ct*ct))
        if st <= 1e-8:
            raise ValueError(f"angle ({i},{j},{k}) is singular at {np.degrees(np.arccos(ct)):.6f} degrees")
        gi = -(v - ct*u) / (ra * st)
        gk = -(u - ct*v) / (rc * st)
        B[iq, 3*i:3*i+3] = gi
        B[iq, 3*k:3*k+3] = gk
        B[iq, 3*j:3*j+3] = -(gi + gk)
        labels.append(("angle", i, j, k))
        iq += 1
    return B, labels


def internal_coordinate_basis(B, masses, rtol=1e-10):
    """Orthonormal basis for the mass-weighted Cartesian space sampled by B.

    The mass-weighted Wilson matrix C = B M^{-1/2} maps mass-weighted Cartesian
    displacements to internal coordinate displacements.  Its row space (obtained
    via SVD: C = U S V^T) spans the subspace of Cartesian space that is physically
    described by the internal coordinates.  The columns of V[:rank]^T form an
    orthonormal basis for this subspace.

    Redundancy: when there are more internal coordinates than independent ones
    (e.g. all bond lengths + all angles in a ring), C has rank < n_internal.
    The SVD truncation handles this automatically.
    """
    B = np.asarray(B, dtype=np.float64)
    masses = np.asarray(masses, dtype=np.float64)
    n = masses.size * 3
    if B.ndim != 2 or B.shape[1] != n:
        raise ValueError(f"Wilson matrix shape {B.shape} is incompatible with {n} coordinates")
    C = B / np.sqrt(np.repeat(masses, 3))[None, :]
    _, s, Vt = np.linalg.svd(C, full_matrices=False)
    if s.size == 0:
        raise ValueError("Wilson matrix has no internal coordinates")
    rank = int(np.sum(s > rtol * max(C.shape) * s[0]))
    if rank == 0:
        raise ValueError("Wilson matrix is numerically rank zero")
    return Vt[:rank].T, s[:rank]


def internal_hessian_projection(H, B, masses, coordinate_scale=None, rtol=1e-10):
    """Project Cartesian Hessian(s) to the least-norm redundant valence force matrix.

    WILSON GF METHOD IMPLEMENTATION:
    In the harmonic approximation, E = ½ q^T F q, giving H = B^T F B in Cartesian
    coordinates.  The mass-weighted dynamical matrix is D = M^{-1/2} H M^{-1/2}.
    Substituting C = B M^{-1/2} (the mass-weighted, scaled Wilson matrix):
        D = C^T F C
    To recover F from D, we use the pseudoinverse C^+ (via SVD):
        F = C^{+T} D C^+
    This is the least-norm solution when coordinates are redundant.

    DIMENSIONLESS COORDINATES:
    Bond displacements are scaled as Δr/r₀ and angles in radians, giving
    dimensionless internal coordinates.  This makes F dimensionless and
    comparable across different coordinate types.

    The resulting F matrix reveals the valence force constants and their couplings.
    For an independent-coordinate model (diagonal F), only the diagonal elements
    are fitted.  Off-diagonal elements (stretch-stretch, stretch-bend) indicate
    where the independent model is incomplete.
    """
    H = np.asarray(H, dtype=np.float64)
    B = np.asarray(B, dtype=np.float64)
    masses = np.asarray(masses, dtype=np.float64)
    if coordinate_scale is None:
        coordinate_scale = np.ones(B.shape[0])
    coordinate_scale = np.asarray(coordinate_scale, dtype=np.float64)
    if coordinate_scale.shape != (B.shape[0],) or np.any(~np.isfinite(coordinate_scale)) or np.any(coordinate_scale <= 0.0):
        raise ValueError("coordinate_scale must contain one finite positive scale per Wilson coordinate")
    Bs = B * coordinate_scale[:, None]
    C = Bs / np.sqrt(np.repeat(masses, 3))[None, :]
    Uq, s, Vt = np.linalg.svd(C, full_matrices=False)
    if s.size == 0:
        raise ValueError("Wilson matrix has no internal coordinates")
    rank = int(np.sum(s > rtol * max(C.shape) * s[0]))
    if rank == 0:
        raise ValueError("Wilson matrix is numerically rank zero")
    Cplus = (Vt[:rank].T / s[:rank][None, :]) @ Uq[:, :rank].T
    D = mass_weight_hessian(H, masses)
    F = (Cplus.T @ D) @ Cplus  # batched BLAS: (npar, n_int, n_int), ~4400x faster than einsum
    return F, {"rank": rank, "singular_values": s[:rank], "scaled_wilson_matrix": Bs, "pseudoinverse": Cplus}


def local_hessian_mask(natoms, bonds, max_graph_distance=2):
    """Cartesian mask for atom blocks separated by at most a bond-graph distance.

    PHYSICAL MOTIVATION:
    A valence force field with only bond and angle terms produces Hessian blocks
    that are nonzero only between atoms sharing a bond or angle (graph distance ≤ 2).
    Atoms farther apart in the bond graph have zero Hessian blocks in the model.

    Fitting the full QM Hessian element-by-element penalizes these zero blocks
    against nonzero QM long-range blocks that the reduced model is structurally
    incapable of generating.  This mask restricts the fit to only the blocks
    the model can represent, preventing long-range QM information from dominating
    the objective and distorting the local parameters.

    The graph distance is computed by BFS on the bond network.  Default d_max=2
    covers 1-2 (bond) and 1-3 (angle) interactions.  d_max=3 would include
    1-4 (torsion) blocks.
    """
    if max_graph_distance < 0:
        raise ValueError("max_graph_distance must be non-negative")
    adj = [[] for _ in range(natoms)]
    for bond in bonds:
        i, j = int(bond[0]), int(bond[1])
        adj[i].append(j)
        adj[j].append(i)
    block = np.zeros((natoms, natoms), dtype=bool)
    for src in range(natoms):
        dist = np.full(natoms, -1, dtype=int)
        dist[src] = 0
        queue = deque([src])
        while queue:
            i = queue.popleft()
            if dist[i] == max_graph_distance:
                continue
            for j in adj[i]:
                if dist[j] < 0:
                    dist[j] = dist[i] + 1
                    queue.append(j)
        block[src, dist >= 0] = dist[dist >= 0] <= max_graph_distance
    return np.repeat(np.repeat(block, 3, axis=0), 3, axis=1)


def _normalized_component(X, b, weight, name):
    if weight < 0.0:
        raise ValueError(f"{name} weight must be non-negative")
    if weight == 0.0:
        return None
    scale = np.linalg.norm(b)
    if not np.isfinite(scale) or scale <= 1e-30:
        raise ValueError(f"{name} reference norm is zero or non-finite")
    fac = np.sqrt(weight) / scale
    return fac * X, fac * b


def assemble_hybrid_hessian_system(A, H_ref, positions, masses, bonds, angles,
                                   mode_weight=1.0, local_weight=1.0, internal_weight=1.0,
                                   mode_balance="frequency", mode_mixing=0.1, frequency_floor_cm1=50.0,
                                   local_graph_distance=2):
    """Assemble direct linear residuals for an all-mode/local/internal hybrid fit.

    HYBRID OBJECTIVE DERIVATION:
    The simple Hessian-element fit asks "which FF Hessian is closest to the QM
    Hessian element-by-element?"  A better question for vibrational spectra is
    "which FF reproduces the action of the QM dynamical matrix on the vibrational
    modes that matter?"  This function combines three complementary objectives:

    1. MODE TERM (Rayleigh quotient targeting):
       For each reference vibrational mode v_s with eigenvalue λ_s = ω_s²:
           Target: v_s^T D(k) v_s = λ_s
       Since D(k) = Σ_p k_p D_p is linear in k, this is a linear equation:
           Σ_p k_p (v_s^T D_p v_s) = λ_s
       Off-diagonal mode mixing v_s^T D(k) v_t (s≠t) is penalized separately
       to preserve eigenvector structure.

       Mode balancing weights:
         - 'equal': w_s = 1 (all modes weighted equally)
         - 'frequency': w_s = 1/max(|λ_s|, λ_floor) — first-order weighting
           appropriate to squared frequency errors (δω² ~ δλ/λ)
         - 'relative': w_s = 1/max(|λ_s|, λ_floor)² — relative eigenvalue accuracy

    2. LOCAL TERM (graph-local Hessian blocks):
       Only mass-weighted Cartesian atom blocks within bond-graph distance
       d_max are fitted: ||D_ij^FF - D_ij^QM||² for d_graph(i,j) ≤ d_max.
       This prevents absent long-range Hessian blocks from dominating the
       reduced local model.

    3. INTERNAL TERM (Wilson GF projection):
       Projects both reference and model Hessians into the Wilson internal-
       coordinate basis: F = C^{+T} D C^+, then fits ||F^FF - F^QM||².
       This targets the valence force constants directly.

    Each term is normalized by its own reference norm before applying user weight,
    so the three components are dimensionless and comparable.  All three are
    linear in k, so the combined problem is a single linear least-squares solve.
    """
    A = np.asarray(A, dtype=np.float64)
    H_ref = np.asarray(H_ref, dtype=np.float64)
    npar, n, n2 = A.shape
    if n != n2 or H_ref.shape != (n, n):
        raise ValueError(f"invalid sensitivity/reference shapes A={A.shape}, H_ref={H_ref.shape}")
    if np.any(~np.isfinite(A)) or np.any(~np.isfinite(H_ref)):
        raise ValueError("sensitivities or reference Hessian contain NaN or infinity")
    Dp = mass_weight_hessian(A, masses)
    Dref = mass_weight_hessian(H_ref, masses)
    V, lam = reference_vibrational_modes(H_ref, positions, masses)
    parts, targets, slices = [], [], {}

    if mode_weight > 0.0:
        if mode_balance == "equal":
            w = np.ones_like(lam)
        elif mode_balance == "frequency":
            lam_floor = (frequency_floor_cm1 / 521.5)**2
            w = 1.0 / np.maximum(np.abs(lam), lam_floor)
        elif mode_balance == "relative":
            lam_floor = (frequency_floor_cm1 / 521.5)**2
            w = 1.0 / np.maximum(np.abs(lam), lam_floor)**2
        else:
            raise ValueError("mode_balance must be 'equal', 'frequency', or 'relative'")
        if mode_mixing < 0.0:
            raise ValueError("mode_mixing must be non-negative")
        Kp = (V.T @ Dp) @ V  # batched BLAS: (npar, nmodes, nmodes), ~2500x faster than einsum
        imode = np.arange(lam.size)
        Xdiag = Kp[:, imode, imode].T * np.sqrt(w)[:, None]
        bdiag = lam * np.sqrt(w)
        Xmode, bmode = Xdiag, bdiag
        if mode_mixing > 0.0 and lam.size > 1:
            off = ~np.eye(lam.size, dtype=bool)
            wmix = (w[:, None] * w[None, :])**0.25
            Xoff = Kp[:, off].T * wmix[off, None]
            Xmode = np.vstack((Xdiag, np.sqrt(mode_mixing) * Xoff))
            bmode = np.concatenate((bdiag, np.zeros(Xoff.shape[0])))
        component = _normalized_component(Xmode, bmode, mode_weight, "mode")
        parts.append(component[0]); targets.append(component[1]); slices["mode"] = parts[-1].shape[0]

    if local_weight > 0.0:
        mask = local_hessian_mask(len(masses), bonds, local_graph_distance)
        component = _normalized_component(Dp[:, mask].T, Dref[mask], local_weight, "local Hessian")
        parts.append(component[0]); targets.append(component[1]); slices["local"] = parts[-1].shape[0]

    B, labels = build_wilson_matrix(positions, bonds, angles)
    coordinate_scale = np.ones(B.shape[0])
    for iq, label in enumerate(labels):
        if label[0] == "bond":
            coordinate_scale[iq] = 1.0 / np.linalg.norm(np.asarray(positions)[label[2]] - np.asarray(positions)[label[1]])
    Fref, internal_info = internal_hessian_projection(H_ref, B, masses, coordinate_scale=coordinate_scale)
    if internal_weight > 0.0:
        Fp, _ = internal_hessian_projection(A, B, masses, coordinate_scale=coordinate_scale)
        component = _normalized_component(Fp.reshape(npar, -1).T, Fref.ravel(), internal_weight, "internal Hessian")
        parts.append(component[0]); targets.append(component[1]); slices["internal"] = parts[-1].shape[0]

    if not parts:
        raise ValueError("at least one hybrid objective weight must be positive")
    return np.vstack(parts), np.concatenate(targets), {
        "n_modes": V.shape[1], "mode_eigenvalues": lam, "internal_rank": internal_info["rank"],
        "n_internal_coordinates": B.shape[0], "wilson_singular_values": internal_info["singular_values"],
        "wilson_matrix": B, "coordinate_scale": coordinate_scale, "coordinate_labels": labels,
        "component_rows": slices,
    }


def solve_regularized_lsq(X, y, prior=None, regularization=0.0, parameter_scale=None,
                          bounds=(0.0, np.inf), rcond=None):
    """Solve a scaled direct least-squares problem, optionally with bounds and a prior.

    MATHEMATICAL FORMULATION:
    The problem is:  min_k ||X k - y||²  subject to lo ≤ k ≤ hi
    with optional Tikhonov regularization:  + ρ Σ_p ((k_p - k_p^0)/σ_p)²

    The regularization pulls poorly constrained parameters toward physical priors
    k_p^0 with strength controlled by ρ and per-parameter scale σ_p.  This is
    equivalent to augmenting the system:
        [X]         [y]
        [R] k  ≈   [R k^0]   where R = diag(sqrt(ρ)/σ_p)

    COLUMN SCALING: Each column of X is normalized by its norm before solving,
    so all parameters have comparable sensitivity.  Bounds are scaled accordingly.
    This prevents parameters with large sensitivity columns from dominating.

    SOLVER CHOICE:
    - With bounds: scipy.optimize.lsq_linear (trust-region reflective method)
      This avoids forming normal equations (which square the condition number:
      κ(G) = κ(X)²) and handles bounds natively.
    - Without bounds: numpy.linalg.lstsq (SVD-based, optimal conditioning)

    UNOBSERVABLE PARAMETERS: Columns with near-zero norm (≤ 1e-30) indicate
    parameters that the data cannot constrain at all — these raise an explicit
    error rather than producing arbitrary values.
    """
    X = np.asarray(X, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if X.ndim != 2 or y.shape != (X.shape[0],):
        raise ValueError(f"incompatible least-squares shapes X={X.shape}, y={y.shape}")
    if np.any(~np.isfinite(X)) or np.any(~np.isfinite(y)):
        raise ValueError("least-squares system contains NaN or infinity")
    npar = X.shape[1]
    prior = np.zeros(npar) if prior is None else np.asarray(prior, dtype=np.float64)
    if prior.shape != (npar,):
        raise ValueError(f"prior shape {prior.shape} does not match {npar} parameters")
    if parameter_scale is None:
        parameter_scale = np.maximum(np.abs(prior), 1.0)
    parameter_scale = np.asarray(parameter_scale, dtype=np.float64)
    if parameter_scale.shape != (npar,) or np.any(~np.isfinite(parameter_scale)) or np.any(parameter_scale <= 0.0):
        raise ValueError("parameter_scale must contain one finite positive value per parameter")
    if regularization < 0.0:
        raise ValueError("regularization must be non-negative")
    if regularization > 0.0:
        R = np.diag(np.sqrt(regularization) / parameter_scale)
        X = np.vstack((X, R))
        y = np.concatenate((y, R @ prior))

    col_norm = np.linalg.norm(X, axis=0)
    if np.any(col_norm <= 1e-30):
        bad = np.flatnonzero(col_norm <= 1e-30)
        raise np.linalg.LinAlgError(f"unobservable force-field parameters: {bad.tolist()}")
    Xs = X / col_norm[None, :]
    lo, hi = bounds
    lo = np.broadcast_to(np.asarray(lo, dtype=np.float64), (npar,)) * col_norm
    hi = np.broadcast_to(np.asarray(hi, dtype=np.float64), (npar,)) * col_norm
    bounded = np.any(np.isfinite(lo)) or np.any(np.isfinite(hi))
    if bounded:
        from scipy.optimize import lsq_linear
        result = lsq_linear(Xs, y, bounds=(lo, hi), method="trf", tol=1e-12, lsmr_tol=1e-12, max_iter=1000, verbose=0)
        if not result.success:
            raise RuntimeError(f"bounded least-squares failed: {result.message}")
        z = result.x
        solver = "scipy.optimize.lsq_linear"
        optimality = result.optimality
    else:
        z, _, _, _ = np.linalg.lstsq(Xs, y, rcond=rcond)
        solver = "numpy.linalg.lstsq"
        optimality = np.linalg.norm(Xs.T @ (Xs @ z - y), ord=np.inf)
    k = z / col_norm
    s = np.linalg.svd(Xs, compute_uv=False)
    residual = X @ k - y
    return k, {
        "solver": solver, "residual_norm": np.linalg.norm(residual),
        "relative_residual": np.linalg.norm(residual) / np.linalg.norm(y),
        "rank": int(np.sum(s > (np.finfo(float).eps * max(Xs.shape) if rcond is None else rcond) * s[0])),
        "condition": s[0] / s[-1] if s[-1] > 0.0 else np.inf,
        "singular_values": s, "column_norms": col_norm, "optimality": optimality,
    }


def fit_hybrid_hessian(systems, mode_weight=1.0, local_weight=1.0, internal_weight=1.0,
                       mode_balance="frequency", mode_mixing=0.1, frequency_floor_cm1=50.0, local_graph_distance=2,
                       prior=None, regularization=1e-2, parameter_scale=None, bounds=(0.0, np.inf)):
    """Fit shared linear FF parameters to one or more hybrid Hessian objectives.

    MULTI-SYSTEM ACCUMULATION:
    When fitting across multiple systems (e.g. SiH4 + Si nanocrystals of different
    sizes), the same global parameters appear in all systems.  Each system's
    hybrid residual (mode + local + internal) is stacked vertically and weighted
    equally by 1/sqrt(n_systems).  This ensures each system contributes equally
    to the combined objective regardless of system size.

    The combined stacked system is then solved by solve_regularized_lsq with
    column scaling, Tikhonov regularization toward physical priors, and
    non-negative bounds (k_p >= 0 for physical stiffnesses).
    """
    Xs, ys, system_info = [], [], []
    npar = None
    for system in systems:
        X, y, info = assemble_hybrid_hessian_system(
            system["A"], system["H_ref"], system["positions"], system["masses"], system["bonds"], system["angles"],
            mode_weight=mode_weight, local_weight=local_weight, internal_weight=internal_weight,
            mode_balance=mode_balance, mode_mixing=mode_mixing, frequency_floor_cm1=frequency_floor_cm1, local_graph_distance=local_graph_distance)
        if npar is None:
            npar = X.shape[1]
        elif X.shape[1] != npar:
            raise ValueError("all systems must use the same global parameter count")
        Xs.append(X / np.sqrt(len(systems)))
        ys.append(y / np.sqrt(len(systems)))
        system_info.append(info)
    if not Xs:
        raise ValueError("no systems supplied")
    k, diagnostics = solve_regularized_lsq(np.vstack(Xs), np.concatenate(ys), prior=prior, regularization=regularization, parameter_scale=parameter_scale, bounds=bounds)
    diagnostics["systems"] = system_info
    return k, diagnostics
