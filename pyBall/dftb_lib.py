"""
DFTB+ C API wrapper for Python

Provides Hamiltonian, Overlap, and Density matrix access using serial DFTB+ build.
No MPI required - uses dftbp_init() instead of dftbp_init_mpi().

Based on working implementation in test_dftb_c_api.py
"""

import os
import sys
import numpy as np
from ctypes import CDLL, c_int, c_void_p, c_double, POINTER, c_char_p, Structure, CFUNCTYPE, byref, cast


# ============ DFTB+ Structures ============

class DftbPlus(Structure):
    """DFTB+ instance structure"""
    _fields_ = [("pDftbPlus", c_void_p)]


class DftbPlusInput(Structure):
    """DFTB+ input structure"""
    _fields_ = [("pDftbPlusInput", c_void_p)]


# ============ Callback Types ============

DMHSCallBackFunc = CFUNCTYPE(None, c_void_p, c_int, c_int, POINTER(c_int), c_void_p)


# ============ DFTB+ Calculator Class ============

class DftbPlusCalculator:
    """
    DFTB+ calculator with Hamiltonian/Overlap/Density matrix access.
    
    Uses serial DFTB+ C API (dftbp_init) - no MPI required.
    Extracts matrices via callbacks during SCF calculation.
    
    Example:
        >>> calc = DftbPlusCalculator(lib_path="/path/to/libdftbplus.so")
        >>> calc.initialize(input_file="input.dftb")
        >>> calc.register_callbacks()
        >>> energy = calc.calculate()
        >>> H = calc.get_hamiltonian(iK=1, iS=1)
        >>> S = calc.get_overlap(iK=1, iS=1)
        >>> DM = calc.get_density_matrix(iK=1, iS=1)
        >>> calc.finalize()
    """
    
    def __init__(self, lib_path="/home/prokophapala/opt/dftb-asi/lib/libdftbplus.so"):
        """
        Initialize DFTB+ calculator by loading library.
        
        Args:
            lib_path: Path to libdftbplus.so
        """
        self.lib_path = lib_path
        self.dftblib = None
        self.instance = None
        self.input_obj = None
        self.basis_size = 0
        self.nr_atoms = 0
        
        # Matrix storage: (iK, iS) -> numpy array
        self.hamiltonian_storage = {}
        self.overlap_storage = {}
        self.dm_storage = {}
        
        # Callback wrappers (must keep references)
        self.h_cb_wrapper = None
        self.s_cb_wrapper = None
        self.dm_cb_wrapper = None
        
        self._load_library()
        self._setup_api_functions()
    
    def _load_library(self):
        """Load DFTB+ shared library"""
        if not os.path.exists(self.lib_path):
            raise FileNotFoundError(f"DFTB+ library not found: {self.lib_path}")
        
        self.dftblib = CDLL(self.lib_path, mode=0x2)  # RTLD_GLOBAL
    
    def _setup_api_functions(self):
        """Define ctypes function signatures"""
        # Initialization
        self.dftbp_init = self.dftblib.dftbp_init
        self.dftbp_init.argtypes = [POINTER(DftbPlus), c_char_p]
        self.dftbp_init.restype = None
        
        self.dftbp_final = self.dftblib.dftbp_final
        self.dftbp_final.argtypes = [POINTER(DftbPlus)]
        self.dftbp_final.restype = None
        
        # Input handling
        self.dftbp_get_input_from_file = self.dftblib.dftbp_get_input_from_file
        self.dftbp_get_input_from_file.argtypes = [POINTER(DftbPlus), c_char_p, POINTER(DftbPlusInput)]
        self.dftbp_get_input_from_file.restype = None
        
        self.dftbp_process_input = self.dftblib.dftbp_process_input
        self.dftbp_process_input.argtypes = [POINTER(DftbPlus), POINTER(DftbPlusInput)]
        self.dftbp_process_input.restype = None
        
        # System info
        self.dftbp_get_nr_atoms = self.dftblib.dftbp_get_nr_atoms
        self.dftbp_get_nr_atoms.argtypes = [POINTER(DftbPlus)]
        self.dftbp_get_nr_atoms.restype = c_int
        
        self.dftbp_get_basis_size = self.dftblib.dftbp_get_basis_size
        self.dftbp_get_basis_size.argtypes = [POINTER(DftbPlus)]
        self.dftbp_get_basis_size.restype = c_int
        
        self.dftbp_is_hs_real = self.dftblib.dftbp_is_hs_real
        self.dftbp_is_hs_real.argtypes = [POINTER(DftbPlus)]
        self.dftbp_is_hs_real.restype = c_int
        
        # Energy
        self.dftbp_get_energy = self.dftblib.dftbp_get_energy
        self.dftbp_get_energy.argtypes = [POINTER(DftbPlus), POINTER(c_double)]
        self.dftbp_get_energy.restype = None
        
        # Callbacks
        self.dftbp_register_h_callback = self.dftblib.dftbp_register_h_callback
        self.dftbp_register_h_callback.argtypes = [POINTER(DftbPlus), DMHSCallBackFunc, c_void_p]
        self.dftbp_register_h_callback.restype = None
        
        self.dftbp_register_s_callback = self.dftblib.dftbp_register_s_callback
        self.dftbp_register_s_callback.argtypes = [POINTER(DftbPlus), DMHSCallBackFunc, c_void_p]
        self.dftbp_register_s_callback.restype = None
        
        self.dftbp_register_dm_callback = self.dftblib.dftbp_register_dm_callback
        self.dftbp_register_dm_callback.argtypes = [POINTER(DftbPlus), DMHSCallBackFunc, c_void_p]
        self.dftbp_register_dm_callback.restype = None
    
    def initialize(self, input_file, output_file="dftb_output.txt"):
        """
        Initialize DFTB+ instance and load input file.
        
        Args:
            input_file: Path to DFTB+ input file (.dftb or .hsd)
            output_file: Path for DFTB+ output log
        """
        if not os.path.exists(input_file):
            raise FileNotFoundError(f"Input file not found: {input_file}")
        
        # Create instance
        self.instance = DftbPlus()
        self.dftbp_init(byref(self.instance), output_file.encode('utf-8'))
        
        # Load input
        self.input_obj = DftbPlusInput()
        self.dftbp_get_input_from_file(
            byref(self.instance), 
            input_file.encode('utf-8'), 
            byref(self.input_obj)
        )
        
        # Process input
        self.dftbp_process_input(byref(self.instance), byref(self.input_obj))
        
        # Get system info
        self.nr_atoms = self.dftbp_get_nr_atoms(byref(self.instance))
        self.basis_size = self.dftbp_get_basis_size(byref(self.instance))
        self.is_real = bool(self.dftbp_is_hs_real(byref(self.instance)))
    
    def _hamiltonian_callback(self, aux_ptr, iK, iS, blacs_descr, blacs_data):
        """Internal callback for Hamiltonian matrix"""
        if blacs_data and self.basis_size > 0:
            n = self.basis_size
            arr_ptr = cast(blacs_data, POINTER(c_double))
            arr = np.ctypeslib.as_array(arr_ptr, shape=(n, n))
            self.hamiltonian_storage[(iK, iS)] = arr.copy()
    
    def _overlap_callback(self, aux_ptr, iK, iS, blacs_descr, blacs_data):
        """Internal callback for Overlap matrix"""
        if blacs_data and self.basis_size > 0:
            n = self.basis_size
            arr_ptr = cast(blacs_data, POINTER(c_double))
            arr = np.ctypeslib.as_array(arr_ptr, shape=(n, n))
            self.overlap_storage[(iK, iS)] = arr.copy()
    
    def _dm_callback(self, aux_ptr, iK, iS, blacs_descr, blacs_data):
        """Internal callback for Density matrix"""
        if blacs_data and self.basis_size > 0:
            n = self.basis_size
            arr_ptr = cast(blacs_data, POINTER(c_double))
            arr = np.ctypeslib.as_array(arr_ptr, shape=(n, n))
            self.dm_storage[(iK, iS)] = arr.copy()
    
    def register_callbacks(self):
        """Register callbacks for Hamiltonian, Overlap, and Density matrix extraction"""
        # Create callback wrappers (must keep references to prevent garbage collection)
        self.h_cb_wrapper = DMHSCallBackFunc(self._hamiltonian_callback)
        self.s_cb_wrapper = DMHSCallBackFunc(self._overlap_callback)
        self.dm_cb_wrapper = DMHSCallBackFunc(self._dm_callback)
        
        # Register with DFTB+
        self.dftbp_register_h_callback(byref(self.instance), self.h_cb_wrapper, None)
        self.dftbp_register_s_callback(byref(self.instance), self.s_cb_wrapper, None)
        self.dftbp_register_dm_callback(byref(self.instance), self.dm_cb_wrapper, None)
    
    def calculate(self):
        """
        Run DFTB+ calculation and extract energy.
        
        Returns:
            Total energy in Hartree
        """
        energy = c_double()
        self.dftbp_get_energy(byref(self.instance), byref(energy))
        return energy.value
    
    def get_hamiltonian(self, iK=1, iS=1):
        """
        Get Hamiltonian matrix.
        
        Args:
            iK: K-point index (1-based)
            iS: Spin index (1-based)
        
        Returns:
            Hamiltonian matrix as numpy array of shape (basis_size, basis_size)
        """
        if (iK, iS) not in self.hamiltonian_storage:
            raise KeyError(f"Hamiltonian matrix not available for k-point {iK}, spin {iS}")
        return self.hamiltonian_storage[(iK, iS)]
    
    def get_overlap(self, iK=1, iS=1):
        """
        Get Overlap matrix.
        
        Args:
            iK: K-point index (1-based)
            iS: Spin index (1-based)
        
        Returns:
            Overlap matrix as numpy array of shape (basis_size, basis_size)
        """
        if (iK, iS) not in self.overlap_storage:
            raise KeyError(f"Overlap matrix not available for k-point {iK}, spin {iS}")
        return self.overlap_storage[(iK, iS)]
    
    def get_density_matrix(self, iK=1, iS=1):
        """
        Get Density matrix.
        
        Args:
            iK: K-point index (1-based)
            iS: Spin index (1-based)
        
        Returns:
            Density matrix as numpy array of shape (basis_size, basis_size)
        """
        if (iK, iS) not in self.dm_storage:
            raise KeyError(f"Density matrix not available for k-point {iK}, spin {iS}")
        return self.dm_storage[(iK, iS)]
    
    def get_electron_count(self, iK=1, iS=1):
        """
        Calculate electron count from density matrix: Tr(S * DM)
        
        Args:
            iK: K-point index (1-based)
            iS: Spin index (1-based)
        
        Returns:
            Number of electrons
        """
        S = self.get_overlap(iK, iS)
        DM = self.get_density_matrix(iK, iS)
        return np.sum(S * DM.T)
    
    def get_eigenvalues(self, iK=1, iS=1):
        """
        Calculate eigenvalues from Hamiltonian and Overlap matrices.
        
        Solves generalized eigenvalue problem: H * C = E * S * C
        
        Args:
            iK: K-point index (1-based)
            iS: Spin index (1-based)
        
        Returns:
            Eigenvalues as numpy array
        """
        H = self.get_hamiltonian(iK, iS)
        S = self.get_overlap(iK, iS)
        return np.linalg.eigvals(np.linalg.solve(S, H))
    
    def finalize(self):
        """Finalize DFTB+ instance and free resources"""
        if self.instance:
            self.dftbp_final(byref(self.instance))
            self.instance = None
    
    def __enter__(self):
        """Context manager entry"""
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit - auto-finalize"""
        self.finalize()


# ============ Convenience Functions ============

def calculate_with_matrices(input_file, lib_path="/home/prokophapala/opt/dftb-asi/lib/libdftbplus.so"):
    """
    Convenience function: run DFTB+ calculation and extract all matrices.
    
    Args:
        input_file: Path to DFTB+ input file
        lib_path: Path to libdftbplus.so
    
    Returns:
        Dictionary with:
            - energy: Total energy (Hartree)
            - hamiltonian: H matrix (numpy array)
            - overlap: S matrix (numpy array)
            - density_matrix: DM matrix (numpy array)
            - nr_atoms: Number of atoms
            - basis_size: Basis size
            - electron_count: Tr(S*DM)
    """
    with DftbPlusCalculator(lib_path=lib_path) as calc:
        calc.initialize(input_file)
        calc.register_callbacks()
        energy = calc.calculate()
        
        result = {
            'energy': energy,
            'hamiltonian': calc.get_hamiltonian(),
            'overlap': calc.get_overlap(),
            'density_matrix': calc.get_density_matrix(),
            'nr_atoms': calc.nr_atoms,
            'basis_size': calc.basis_size,
            'electron_count': calc.get_electron_count()
        }
    
    return result
