"""
Python ctypes interface to libwaveplot.so (DFTB+ waveplot library).

Array convention:
  Fortran is column-major, Python/C is row-major.
  - coords[natoms, 3]   -> pass np.asfortranarray(coords) or coords.T
  - eigvecs[nStates, nOrb] -> pass eigvecs.T (Fortran: eigvecs[nOrb, nStates])
  - points[npoints, 3]  -> pass np.asfortranarray(points) or points.T
  - gridVecs[3, 3]      -> step vectors; each column is a grid axis

Usage:
    from pyBall.WavePlot.WavePlot import WavePlot
    wp = WavePlot('/path/to/libwaveplot.so')
    wp.set_geometry(coords_bohr, species)
    wp.set_basis(...)
    wp.set_eigenvectors(eigvecs)
    values = wp.orb2points(istate, points_bohr)
"""

import ctypes
import os
import numpy as np

# ---- default library location (build tree) ----
_DEFAULT_LIB = os.path.join(
    os.path.dirname(__file__),
    '..', '..', '_build', 'app', 'waveplot', 'libwaveplot.so'
)
_DEFAULT_LIB = os.path.normpath(_DEFAULT_LIB)

c_double_p = ctypes.POINTER(ctypes.c_double)
c_int_p    = ctypes.POINTER(ctypes.c_int)

def _dp(arr):
    """Return contiguous float64 array pointer."""
    arr = np.asarray(arr, dtype=np.float64, order='F')
    return arr.ctypes.data_as(c_double_p), arr  # keep ref to avoid GC

def _ip(arr):
    """Return contiguous int32 array pointer."""
    arr = np.asarray(arr, dtype=np.int32, order='F')
    return arr.ctypes.data_as(c_int_p), arr


class WavePlot:
    """Wrapper around libwaveplot.so C-bindable Fortran interface."""

    def __init__(self, libpath=None):
        if libpath is None:
            libpath = _DEFAULT_LIB
        if not os.path.exists(libpath):
            raise FileNotFoundError(f"libwaveplot.so not found at: {libpath}")
        self._lib = ctypes.CDLL(libpath, mode=ctypes.RTLD_LOCAL)
        self._setup_signatures()
        self._lib.waveplot_init()

    def _setup_signatures(self):
        lib = self._lib

        lib.waveplot_init.restype  = None
        lib.waveplot_init.argtypes = []

        lib.waveplot_get_nOrb.restype  = None
        lib.waveplot_get_nOrb.argtypes = [c_int_p]

        lib.waveplot_set_geometry.restype  = None
        lib.waveplot_set_geometry.argtypes = [
            ctypes.c_int,   # natoms
            ctypes.c_int,   # isPeriodic
            c_double_p,     # coords(3,natoms) Fortran
            c_int_p,        # species(natoms)
        ]

        lib.waveplot_set_eigenvectors.restype  = None
        lib.waveplot_set_eigenvectors.argtypes = [
            ctypes.c_int,   # nOrb
            ctypes.c_int,   # nStates
            c_double_p,     # eigvecs(nOrb,nStates) Fortran
        ]

        lib.waveplot_orb2points.restype  = None
        lib.waveplot_orb2points.argtypes = [
            ctypes.c_int,   # iState (1-indexed)
            ctypes.c_int,   # npoints
            c_double_p,     # points(3,npoints) Fortran
            c_double_p,     # out(npoints)
        ]

        lib.waveplot_allorbs2point.restype  = None
        lib.waveplot_allorbs2point.argtypes = [
            c_double_p,     # point(3)
            ctypes.c_int,   # nStates
            c_double_p,     # out(nStates)
        ]

        lib.waveplot_orb2grid.restype  = None
        lib.waveplot_orb2grid.argtypes = [
            ctypes.c_int,   # iState (1-indexed)
            c_double_p,     # origin(3)
            c_double_p,     # gridVecs(3,3) Fortran
            c_int_p,        # nPoints(3)
            c_double_p,     # out(n1*n2*n3)
        ]

        # waveplot_set_basis is called via set_basis() helper below
        lib.waveplot_set_basis.restype  = None
        lib.waveplot_set_basis.argtypes = [
            ctypes.c_int,   # nSpecies
            ctypes.c_int,   # nOrbMax
            ctypes.c_int,   # nPowMax
            ctypes.c_int,   # nAlphaMax
            c_int_p,        # nOrb_arr(nSpecies)
            c_int_p,        # angMom_arr(nOrbMax, nSpecies)  Fortran
            c_double_p,     # cutoff_arr(nOrbMax, nSpecies)  Fortran
            c_double_p,     # occ_arr   (nOrbMax, nSpecies)  Fortran
            c_int_p,        # nAlpha_arr(nOrbMax, nSpecies)  Fortran
            c_int_p,        # nPow_arr  (nOrbMax, nSpecies)  Fortran
            c_double_p,     # alpha_flat(nAlphaMax,nOrbMax,nSpecies) Fortran
            c_double_p,     # aa_flat   (nPowMax,nAlphaMax,nOrbMax,nSpecies) Fortran
            ctypes.c_double,# resolution
        ]

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def get_nOrb(self):
        """Return nOrb (after set_eigenvectors) or -1."""
        n = ctypes.c_int(-1)
        self._lib.waveplot_get_nOrb(ctypes.byref(n))
        return n.value

    def set_geometry(self, coords, species, is_periodic=False):
        """
        coords: array [natoms, 3] in Bohr (float64)
        species: array [natoms] of 1-indexed species indices (int)
        """
        coords  = np.asarray(coords,  dtype=np.float64)
        species = np.asarray(species, dtype=np.int32)
        natoms  = coords.shape[0]
        assert coords.shape == (natoms, 3)
        assert species.shape == (natoms,)
        # Fortran wants (3, natoms) -> pass transpose as F-contiguous
        coords_f = np.asfortranarray(coords.T)
        p_c, _c = _dp(coords_f)
        p_s, _s = _ip(species)
        self._lib.waveplot_set_geometry(
            ctypes.c_int(natoms),
            ctypes.c_int(1 if is_periodic else 0),
            p_c, p_s
        )

    def set_eigenvectors(self, eigvecs):
        """
        eigvecs: array [nStates, nOrb] (float64)
        Fortran expects (nOrb, nStates), so we pass eigvecs.T as F-contiguous.
        """
        eigvecs = np.asarray(eigvecs, dtype=np.float64)
        nStates, nOrb = eigvecs.shape
        evF = np.asfortranarray(eigvecs.T)   # shape (nOrb, nStates) F-order
        p, _r = _dp(evF)
        self._lib.waveplot_set_eigenvectors(
            ctypes.c_int(nOrb), ctypes.c_int(nStates), p
        )

    def set_basis(self, basis_list, resolution=1e-3):
        """
        basis_list: list of dicts (one per species), each with keys:
            angMoms    : list[int]   length nOrb
            cutoffs    : list[float] length nOrb (Bohr)
            occupations: list[float] length nOrb
            stos       : list of dicts with keys:
                alpha  : list[float]  length nAlpha   (STO exponents)
                aa     : 2D list[float] shape (nPow, nAlpha) (STO coefficients)
        resolution: float, radial tabulation step (Bohr)
        """
        nSpecies = len(basis_list)
        nOrbMax  = max(len(b['angMoms']) for b in basis_list)
        nAlphaMax = max(
            max(len(s['alpha']) for s in b['stos'])
            for b in basis_list
        )
        nPowMax = max(
            max(len(s['aa']) for s in b['stos'])
            for b in basis_list
        )

        nOrb_arr   = np.zeros((nSpecies,),             dtype=np.int32)
        angMom_arr = np.zeros((nOrbMax, nSpecies),     dtype=np.int32,   order='F')
        cutoff_arr = np.zeros((nOrbMax, nSpecies),     dtype=np.float64, order='F')
        occ_arr    = np.zeros((nOrbMax, nSpecies),     dtype=np.float64, order='F')
        nAlpha_arr = np.zeros((nOrbMax, nSpecies),     dtype=np.int32,   order='F')
        nPow_arr   = np.zeros((nOrbMax, nSpecies),     dtype=np.int32,   order='F')
        alpha_flat = np.zeros((nAlphaMax, nOrbMax, nSpecies), dtype=np.float64, order='F')
        aa_flat    = np.zeros((nPowMax, nAlphaMax, nOrbMax, nSpecies), dtype=np.float64, order='F')

        for iSp, b in enumerate(basis_list):
            nOrb = len(b['angMoms'])
            nOrb_arr[iSp] = nOrb
            for iOrb in range(nOrb):
                angMom_arr[iOrb, iSp] = b['angMoms'][iOrb]
                cutoff_arr[iOrb, iSp] = b['cutoffs'][iOrb]
                occ_arr   [iOrb, iSp] = b['occupations'][iOrb]
                sto = b['stos'][iOrb]
                alpha = np.asarray(sto['alpha'], dtype=np.float64)
                aa    = np.asarray(sto['aa'],    dtype=np.float64)  # (nPow, nAlpha)
                nA = len(alpha); nP = aa.shape[0]
                nAlpha_arr[iOrb, iSp] = nA
                nPow_arr  [iOrb, iSp] = nP
                alpha_flat[:nA, iOrb, iSp]       = alpha
                aa_flat   [:nP, :nA, iOrb, iSp]  = aa

        p_nOrb,   _n   = _ip(nOrb_arr)
        p_ang,    _ang = _ip(angMom_arr)
        p_cut,    _cut = _dp(cutoff_arr)
        p_occ,    _occ = _dp(occ_arr)
        p_nAlpha, _na  = _ip(nAlpha_arr)
        p_nPow,   _np_ = _ip(nPow_arr)
        p_alpha,  _al  = _dp(alpha_flat)
        p_aa,     _aa  = _dp(aa_flat)

        self._lib.waveplot_set_basis(
            ctypes.c_int(nSpecies),
            ctypes.c_int(nOrbMax),
            ctypes.c_int(nPowMax),
            ctypes.c_int(nAlphaMax),
            p_nOrb, p_ang, p_cut, p_occ,
            p_nAlpha, p_nPow, p_alpha, p_aa,
            ctypes.c_double(resolution)
        )

    def orb2points(self, istate, points):
        """
        Evaluate MO istate (1-indexed) at arbitrary points.
        points: array [npoints, 3] in Bohr (float64)
        Returns: array [npoints] of MO values
        """
        points  = np.asarray(points, dtype=np.float64)
        npoints = points.shape[0]
        assert points.shape == (npoints, 3)
        # Fortran wants (3, npoints)
        pts_f = np.asfortranarray(points.T)
        out   = np.zeros(npoints, dtype=np.float64)
        p_pts, _pts = _dp(pts_f)
        p_out, _out = _dp(out)
        self._lib.waveplot_orb2points(
            ctypes.c_int(istate), ctypes.c_int(npoints), p_pts, p_out
        )
        return out

    def allorbs2point(self, point, nStates):
        """
        Evaluate all nStates MOs at a single point.
        point: array [3] in Bohr
        Returns: array [nStates]
        """
        point = np.asarray(point, dtype=np.float64).ravel()
        assert point.shape == (3,)
        out = np.zeros(nStates, dtype=np.float64)
        p_pt,  _pt  = _dp(point)
        p_out, _out = _dp(out)
        self._lib.waveplot_allorbs2point(p_pt, ctypes.c_int(nStates), p_out)
        return out

    def orb2grid(self, istate, origin, gridVecs, nPoints):
        """
        Evaluate MO istate on a regular 3D grid.
        origin   : array [3] in Bohr
        gridVecs : array [3,3]; each row is a step vector along grid axis
                   (Python row-major; transposed before passing to Fortran)
        nPoints  : tuple/array (n1, n2, n3)
        Returns  : array [n1, n2, n3] in Fortran order (column-major, use .T for C order)
        """
        origin   = np.asarray(origin,   dtype=np.float64).ravel()
        gridVecs = np.asarray(gridVecs, dtype=np.float64)  # [3,3]
        nPoints  = np.asarray(nPoints,  dtype=np.int32).ravel()
        assert origin.shape   == (3,)
        assert gridVecs.shape == (3, 3)
        assert nPoints.shape  == (3,)
        n1, n2, n3 = nPoints
        # gridVecs rows are Python step vectors; Fortran expects columns -> pass .T
        gvF = np.asfortranarray(gridVecs.T)
        out = np.zeros(n1 * n2 * n3, dtype=np.float64)
        p_orig, _o = _dp(origin)
        p_gv,   _g = _dp(gvF)
        p_np,   _n = _ip(nPoints)
        p_out,  _r = _dp(out)
        self._lib.waveplot_orb2grid(
            ctypes.c_int(istate), p_orig, p_gv, p_np, p_out
        )
        return out.reshape((n1, n2, n3), order='F')


# ================================================================
# High-level evaluation functions
# ================================================================

def evaluate_mos_on_grid(wp, mo_indices, origin, gridVecs, nPoints):
    """
    Evaluate multiple MOs on a 3D grid.
    
    Args:
        wp: WavePlot instance (already configured with geometry, basis, eigenvectors)
        mo_indices: list of 1-indexed MO indices to evaluate
        origin: array [3] in Bohr
        gridVecs: array [3,3] step vectors (Bohr)
        nPoints: tuple (n1, n2, n3)
    
    Returns:
        grid_vals: list of arrays, each [n1, n2, n3] in Fortran order
    """
    grid_vals = []
    for istate in mo_indices:
        val = wp.orb2grid(istate, origin, gridVecs, nPoints)
        grid_vals.append(val)
    return grid_vals


def evaluate_mos_on_points(wp, mo_indices, points):
    """
    Evaluate multiple MOs at arbitrary points.
    
    Args:
        wp: WavePlot instance (already configured)
        mo_indices: list of 1-indexed MO indices
        points: array [npoints, 3] in Bohr
    
    Returns:
        point_vals: list of arrays, each [npoints]
    """
    point_vals = []
    for istate in mo_indices:
        val = wp.orb2points(istate, points)
        point_vals.append(val)
    return point_vals


def setup_waveplot_from_dftb(dftb_data, libpath=None):
    """
    Configure WavePlot instance from parsed DFTB+ data.
    
    Args:
        dftb_data: dict with keys:
            - coords_bohr: (natoms, 3) array
            - species_wp: (natoms,) 1-indexed species indices
            - basis: list of basis dicts (for set_basis)
            - resolution: float (Bohr)
            - evecs: (nstates, norb) array
        libpath: path to libwaveplot.so (optional)
    
    Returns:
        wp: configured WavePlot instance
    """
    wp = WavePlot(libpath)
    wp.set_geometry(dftb_data['coords_bohr'], dftb_data['species_wp'])
    wp.set_basis(dftb_data['basis'], dftb_data['resolution'])
    wp.set_eigenvectors(dftb_data['evecs'])
    return wp
