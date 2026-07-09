"""Python ctypes wrapper for C++ FFfit force-field parameter fitting library.

Wraps FFfit_lib.cpp → FFfit.h class for fitting bond/angle stiffness parameters
to reference Hessians via linear least-squares or gradient descent.

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
lib.fffit_destroy.argtypes = []

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
