"""
DFTBcore: Simplified Python ctypes interface to DFTB+ core functionality.

Provides direct access to DFTB+ calculations with matrix export.
Pattern follows pyBall/Fireball/FireCore.py and pyBall/WavePlot/WavePlot.py

Usage:
    from pyBall.DFTBcore import DFTBcore
    
    # Initialize from input file
    dftb = DFTBcore(libpath='/path/to/libdftbcore.so')
    dftb.init('h2o.hsd')
    
    # Enable matrix collection
    dftb.enable_matrix_collection(dm=True, h=True, s=True)
    
    # Run SCF
    energy = dftb.run_scf()
    
    # Extract matrices
    dm = dftb.get_dm_dense()
    H = dftb.get_h_dense()
    S = dftb.get_s_dense()
    
    # Cleanup
    dftb.finalize()

Array convention:
  Fortran is column-major, Python/C is row-major.
  - coords[natoms, 3] -> pass as Fortran-contiguous or transpose
"""

import ctypes
import os
import numpy as np

c_int = ctypes.c_int
c_int_p = ctypes.POINTER(c_int)
c_double = ctypes.c_double
c_double_p = ctypes.POINTER(c_double)

# ---- default library location (build tree) ----
# Try multiple possible locations
_POSSIBLE_PATHS = [
    # From pyBall/DFTBcore.py -> ../../../_build/app/dftbcore/libdftbcore.so
    os.path.join(os.path.dirname(__file__), '..', '..', '_build', 'app', 'dftbcore', 'libdftbcore.so'),
    # Alternative: from current working directory
    os.path.join(os.getcwd(), '_build', 'app', 'dftbcore', 'libdftbcore.so'),
    # Installed location
    os.path.expanduser('~/opt/dftbplus/lib/libdftbcore.so'),
    os.path.expanduser('~/git/dftbplus/_build/app/dftbcore/libdftbcore.so'),
]

_DEFAULT_LIB = None
for path in _POSSIBLE_PATHS:
    path = os.path.normpath(os.path.abspath(path))
    if os.path.exists(path):
        _DEFAULT_LIB = path
        break

# Alternative: installed location
_INSTALLED_LIB = os.path.expanduser('~/opt/dftbplus/lib/libdftbcore.so')

c_double_p = ctypes.POINTER(ctypes.c_double)
c_int_p = ctypes.POINTER(ctypes.c_int)


def _dp(arr):
    """Return contiguous float64 array pointer (Fortran order)."""
    arr = np.asarray(arr, dtype=np.float64, order='F')
    return arr.ctypes.data_as(c_double_p), arr


def _ip(arr):
    """Return contiguous int32 array pointer."""
    arr = np.asarray(arr, dtype=np.int32, order='F')
    return arr.ctypes.data_as(c_int_p), arr


class DFTBcore:
    """Wrapper around libdftbcore.so C-bindable Fortran interface."""

    def __init__(self, libpath=None):
        """
        Load the DFTBcore shared library.
        
        Args:
            libpath: Path to libdftbcore.so. If None, searches common locations.
        """
        if libpath is None:
            # Use the first path that exists
            if _DEFAULT_LIB is not None:
                libpath = _DEFAULT_LIB
            elif os.path.exists(_INSTALLED_LIB):
                libpath = _INSTALLED_LIB
            else:
                searched = '\n  '.join(_POSSIBLE_PATHS + [_INSTALLED_LIB])
                raise FileNotFoundError(
                    f"libdftbcore.so not found. Searched:\n  {searched}\n"
                    f"Please build DFTB+ with -DBUILD_SHARED_LIBS=ON or specify libpath."
                )
        
        if not os.path.exists(libpath):
            raise FileNotFoundError(f"libdftbcore.so not found at: {libpath}")
        
        self._lib = ctypes.CDLL(libpath, mode=ctypes.RTLD_LOCAL)
        self._basis_size = 0
        self._setup_signatures()
        
    def _setup_signatures(self):
        """Setup ctypes function signatures."""
        lib = self._lib
        
        # dftbcore_init(input_file, output_file)
        lib.dftbcore_init.restype = None
        lib.dftbcore_init.argtypes = [
            ctypes.c_char_p,  # input file path
            ctypes.c_char_p   # output file path (optional)
        ]
        
        # dftbcore_set_coords(natoms, coords)
        lib.dftbcore_set_coords.restype = None
        lib.dftbcore_set_coords.argtypes = [
            ctypes.c_int,    # natoms
            c_double_p       # coords(3, natoms) in Fortran order
        ]
        
        # dftbcore_set_coords_and_lattice(natoms, coords, lattice)
        lib.dftbcore_set_coords_and_lattice.restype = None
        lib.dftbcore_set_coords_and_lattice.argtypes = [
            ctypes.c_int,    # natoms
            c_double_p,      # coords(3, natoms)
            c_double_p       # lattice(3, 3)
        ]
        
        # dftbcore_enable_matrix_collection(collect_dm, collect_h, collect_s)
        lib.dftbcore_enable_matrix_collection.restype = None
        lib.dftbcore_enable_matrix_collection.argtypes = [
            ctypes.c_int,    # collect_dm (0/1)
            ctypes.c_int,    # collect_h (0/1)
            ctypes.c_int     # collect_s (0/1)
        ]
        
        # dftbcore_set_debug(debug)
        lib.dftbcore_set_debug.restype = None
        lib.dftbcore_set_debug.argtypes = [ctypes.c_int]
        
        # dftbcore_enable_hamiltonian_storage(store)
        lib.dftbcore_enable_hamiltonian_storage.restype = None
        lib.dftbcore_enable_hamiltonian_storage.argtypes = [ctypes.c_int]
        
        # dftbcore_run_scf(energy, ierr)
        lib.dftbcore_run_scf.restype = None
        lib.dftbcore_run_scf.argtypes = [
            c_double_p,      # energy output
            c_int_p          # error code output
        ]
        
        # dftbcore_write_debug_matrices()
        lib.dftbcore_write_debug_matrices.restype = None
        lib.dftbcore_write_debug_matrices.argtypes = []
        
        # dftbcore_get_basis_size(norb)
        lib.dftbcore_get_basis_size.restype = None
        lib.dftbcore_get_basis_size.argtypes = [c_int_p]
        
        # dftbcore_get_energy(energy)
        lib.dftbcore_get_energy.restype = None
        lib.dftbcore_get_energy.argtypes = [c_double_p]
        
        # dftbcore_get_eigvecs_dense(eigvecs(*), eigvals(*), n)  -- flat Fortran column-major buffer
        lib.dftbcore_get_eigvecs_dense.restype = None
        lib.dftbcore_get_eigvecs_dense.argtypes = [c_double_p, c_double_p, c_int]
        
        # dftbcore_get_h_dense(h(*), n)
        lib.dftbcore_get_h_dense.restype = None
        lib.dftbcore_get_h_dense.argtypes = [c_double_p, c_int]
        
        # dftbcore_get_s_dense(s(*), n)
        lib.dftbcore_get_s_dense.restype = None
        lib.dftbcore_get_s_dense.argtypes = [c_double_p, c_int]
        
        # dftbcore_get_dm_dense(dm(*), n)
        lib.dftbcore_get_dm_dense.restype = None
        lib.dftbcore_get_dm_dense.argtypes = [c_double_p, c_int]
        
        # dftbcore_finalize()
        lib.dftbcore_finalize.restype = None
        lib.dftbcore_finalize.argtypes = []

    def init(self, input_file, output_file=None):
        """
        Initialize DFTB+ from an input file.
        
        Args:
            input_file: Path to DFTB+ input file (.hsd or .dftb)
            output_file: Optional path for DFTB+ output log
        """
        input_bytes = input_file.encode('utf-8')
        if output_file:
            output_bytes = output_file.encode('utf-8')
            self._lib.dftbcore_init(input_bytes, output_bytes)
        else:
            self._lib.dftbcore_init(input_bytes, None)

    def set_coords(self, coords):
        """
        Set atomic coordinates (non-periodic system).
        
        Args:
            coords: array [natoms, 3] in Bohr
        """
        coords = np.asarray(coords, dtype=np.float64)
        natoms = coords.shape[0]
        assert coords.shape == (natoms, 3)
        
        # Fortran wants (3, natoms) -> transpose and make Fortran-contiguous
        coords_f = np.asfortranarray(coords.T)
        p_c, _c = _dp(coords_f)
        
        self._lib.dftbcore_set_coords(ctypes.c_int(natoms), p_c)

    def set_coords_and_lattice(self, coords, lattice):
        """
        Set atomic coordinates and lattice vectors (periodic system).
        
        Args:
            coords: array [natoms, 3] in Bohr
            lattice: array [3, 3] lattice vectors in Bohr
        """
        coords = np.asarray(coords, dtype=np.float64)
        lattice = np.asarray(lattice, dtype=np.float64)
        natoms = coords.shape[0]
        
        assert coords.shape == (natoms, 3)
        assert lattice.shape == (3, 3)
        
        coords_f = np.asfortranarray(coords.T)
        lattice_f = np.asfortranarray(lattice.T)
        
        p_c, _c = _dp(coords_f)
        p_l, _l = _dp(lattice_f)
        
        self._lib.dftbcore_set_coords_and_lattice(ctypes.c_int(natoms), p_c, p_l)

    def enable_matrix_collection(self, dm=True, h=True, s=True):
        """
        Enable collection of density matrix, Hamiltonian, and/or overlap matrix.
        Must be called BEFORE run_scf().
        
        Args:
            dm: Collect density matrix
            h: Collect Hamiltonian
            s: Collect overlap matrix
        """
        self._lib.dftbcore_enable_matrix_collection(
            ctypes.c_int(1 if dm else 0),
            ctypes.c_int(1 if h else 0),
            ctypes.c_int(1 if s else 0)
        )

    def set_debug(self, debug=True):
        """
        Enable or disable debug mode.
        
        When enabled, debug matrices are written to files after SCF.
        
        Args:
            debug: Enable debug mode
        """
        self._lib.dftbcore_set_debug(ctypes.c_int(1 if debug else 0))

    def enable_hamiltonian_storage(self, store=True):
        """
        Enable or disable Hamiltonian storage before diagonalization.
        
        When enabled, the Hamiltonian is stored before diagonalization
        (when it contains the actual SCC Hamiltonian), not after
        (when it contains eigenvectors).
        
        Args:
            store: Enable Hamiltonian storage
        """
        self._lib.dftbcore_enable_hamiltonian_storage(ctypes.c_int(1 if store else 0))

    def run_scf(self):
        """
        Run the SCF calculation.
        
        Returns:
            Total energy in Hartree
            
        Raises:
            RuntimeError: If SCF calculation fails
        """
        energy = ctypes.c_double()
        ierr = ctypes.c_int()
        
        self._lib.dftbcore_run_scf(ctypes.byref(energy), ctypes.byref(ierr))
        
        if ierr.value != 0:
            raise RuntimeError(f"DFTB+ SCF calculation failed with error code: {ierr.value}")
        
        # Update basis size
        norb = ctypes.c_int()
        self._lib.dftbcore_get_basis_size(ctypes.byref(norb))
        self._basis_size = norb.value
        
        return energy.value

    def get_basis_size(self):
        """Return the number of basis functions (orbitals)."""
        if self._basis_size == 0:
            norb = ctypes.c_int()
            self._lib.dftbcore_get_basis_size(ctypes.byref(norb))
            self._basis_size = norb.value
        return self._basis_size

    def get_energy(self):
        """
        Get total energy from last SCF calculation.
        
        Returns:
            Total energy in Hartree
        """
        energy = ctypes.c_double()
        self._lib.dftbcore_get_energy(ctypes.byref(energy))
        return energy.value

    def get_eigvecs_dense(self, apply_inverse_lowdin=False):
        """
        Get eigenvectors and eigenvalues.
        
        Args:
            apply_inverse_lowdin: If True, apply S^(-1/2) transform to get atomic basis coefficients
        
        Returns:
            (eigvecs[n,n], eigvals[n]) in C row-major order
        """
        n = self.get_basis_size()
        buf_vecs = np.zeros(n*n, dtype=np.float64)
        buf_vals = np.zeros(n, dtype=np.float64)
        self._lib.dftbcore_get_eigvecs_dense(buf_vecs.ctypes.data_as(c_double_p), buf_vals.ctypes.data_as(c_double_p), c_int(n))
        # Fortran stores column-major: reshape as (n,n) Fortran order then convert to C order
        C = np.asfortranarray(buf_vecs.reshape(n, n, order='F')).T.copy()
        
        if apply_inverse_lowdin:
            S = self.get_s_dense()
            # Compute S^(-1/2) via eigendecomposition: S = X * s * X^T
            eigvals_S, eigvecs_S = np.linalg.eigh(S)
            s_inv_sqrt = 1.0 / np.sqrt(np.maximum(eigvals_S, 1e-12))
            S_inv_sqrt = eigvecs_S @ np.diag(s_inv_sqrt) @ eigvecs_S.T
            # Apply transform: C_atomic = S^(-1/2) * C_lowdin
            C = S_inv_sqrt @ C
        
        return C, buf_vals

    def _get_matrix(self, func_name):
        """Helper: call a Fortran matrix getter (flat Fortran-order buffer, n by value) -> C-order numpy array."""
        n = self.get_basis_size()
        buf = np.zeros(n*n, dtype=np.float64)
        getattr(self._lib, func_name)(buf.ctypes.data_as(c_double_p), c_int(n))
        return buf.reshape(n, n, order='F').T.copy()

    def get_dm_dense(self):   return self._get_matrix('dftbcore_get_dm_dense')
    def get_h_dense(self):    return self._get_matrix('dftbcore_get_h_dense')
    def get_s_dense(self):    return self._get_matrix('dftbcore_get_s_dense')

    def write_debug_matrices(self):
        """
        Write debug matrices to files (debug_H.dat, debug_S.dat, debug_DM.dat).
        
        Only writes if debug mode was enabled via set_debug(True).
        """
        self._lib.dftbcore_write_debug_matrices()

    def print_orbital_coeffs(self, eigvecs, eigenvals, atom_orbital_map=None, max_orbitals=None):
        """
        Print orbital coefficients in tabular format with column headers.
        
        Args:
            eigvecs: (nStates, nOrb) matrix - columns are MOs
            eigenvals: (nStates,) array
            atom_orbital_map: Optional list of (atom_name, orbital_names) tuples.
                            If None, prints without atom labels.
            max_orbitals: Maximum number of orbitals to print (None for all)
        """
        nStates, nOrb = eigvecs.shape
        if max_orbitals is not None:
            nStates = min(nStates, max_orbitals)
        
        print(f"\n{'='*80}")
        print(f"ORBITAL COEFFICIENTS (nStates={nStates}, nOrb={nOrb})")
        print(f"{'='*80}\n")
        
        # Build column headers with zero-padded indices
        if atom_orbital_map is not None:
            col_headers = []
            orb_idx = 0
            for atom_name, orbital_names in atom_orbital_map:
                # Extract element and atom index from "O (atom 0)" -> "O", "0"
                parts = atom_name.split()
                elem = parts[0]
                atom_idx = parts[2].strip(')')
                for orb_name in orbital_names:
                    if orb_idx >= nOrb:
                        break
                    # Create compact header like "O000s", "O000px", etc. with zero-padding
                    header = f"{elem}{int(atom_idx):03d}{orb_name}"
                    col_headers.append(header)
                    orb_idx += 1
        else:
            col_headers = [f"AO{i:03d}" for i in range(nOrb)]
        
        # Calculate column width based on header length
        col_width = max(12, max(len(h) for h in col_headers))
        
        # Print header row
        header_str = f"{'MOs':<6}  {'E[eV]':<12}  |  " + "  ".join(f"{h:>{col_width}s}" for h in col_headers)
        print(header_str)
        print("-" * len(header_str))
        
        # Print each MO row (single line)
        for istate in range(nStates):
            mo_label = f"MO{istate:03d}"
            energy_ev = eigenvals[istate] * 27.2114
            coeffs = [f"{coeff:{col_width}.6f}" for coeff in eigvecs[istate, :]]
            row_str = f"{mo_label:<6}  {energy_ev:12.6f}  |  " + "  ".join(coeffs)
            print(row_str)

    def finalize(self):
        """Finalize DFTB+ and free resources."""
        self._lib.dftbcore_finalize()
        self._basis_size = 0

    def __enter__(self):
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit - auto-finalize."""
        self.finalize()
        return False


# ==============================================================================
# Convenience functions for common workflows
# ==============================================================================

def run_dftb_calculation(input_file, libpath=None, collect_dm=True, collect_h=True, collect_s=True):
    """
    Convenience function: Run DFTB+ calculation and extract all matrices.
    
    Args:
        input_file: Path to DFTB+ input file (.hsd)
        libpath: Path to libdftbcore.so (optional)
        collect_dm: Collect density matrix
        collect_h: Collect Hamiltonian
        collect_s: Collect overlap matrix
        
    Returns:
        Dictionary with:
            - energy: Total energy (Hartree)
            - dm: Density matrix (numpy array)
            - h: Hamiltonian matrix
            - s: Overlap matrix
            - basis_size: Number of orbitals
    """
    with DFTBcore(libpath=libpath) as dftb:
        dftb.init(input_file)
        dftb.enable_matrix_collection(dm=collect_dm, h=collect_h, s=collect_s)
        energy = dftb.run_scf()
        
        result = {
            'energy': energy,
            'basis_size': dftb.get_basis_size()
        }
        
        if collect_dm:
            result['dm'] = dftb.get_dm_dense()
        if collect_h:
            result['h'] = dftb.get_h_dense()
        if collect_s:
            result['s'] = dftb.get_s_dense()
    
    return result
