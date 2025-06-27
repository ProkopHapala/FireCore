import os
import numpy as np
import pyopencl as cl
from .OpenCLBase import OpenCLBase
from . import MMparams

# Size constants for better readability
i32sz = 4  # size of int32 in bytes
f32sz = 4  # size of float32 in bytes

class UFF_CL(OpenCLBase):
    """
    PyOpenCL interface for running UFF calculations on GPU,
    following the structure of the provided C++ reference.
    
    This class inherits from OpenCLBase and implements specific functionality
    for Universal Force Field (UFF) calculations.
    """
    
    def __init__(self, nloc=32, kernel_path=None):
        """
        Initialize the UFF OpenCL environment.
        
        Args:
            nloc (int): Local work group size
            kernel_path (str, optional): Path to the UFF kernel file. If None, auto-detect.
        """
        # Initialize the base class
        super().__init__(nloc=nloc)
        
        # Find and load the UFF kernel file
        if kernel_path is None:
            base_path = os.path.dirname(os.path.abspath(__file__))
            rel_path = "../../cpp/common_resources/cl/UFF.cl"
            kernel_path = os.path.join(base_path, rel_path)
        
        # Load the OpenCL program
        if not self.load_program(kernel_path=kernel_path, bPrint=True):
            print(f"Failed to load UFF kernels from {kernel_path}")
            return
        
        # Initialize system parameters
        self.nSystems = 0
        self.natoms = 0
        self.nbonds = 0
        self.nangles = 0
        self.ndihedrals = 0
        self.ninversions = 0
        self.npbc = 0
        
        # Initialize a2f map parameters
        self.a2f_map_size = 0
        
        # Initialize kernel arguments dictionary
        self.kernel_args = {}
        
        # Flag to track if kernel arguments are set up
        self.args_setup = False
    
    def realloc_buffers(self, natoms, nbonds, nangles, ndihedrals, ninversions, npbc, nSystems=1):
        """
        Allocates or reallocates all necessary OpenCL buffers for UFF calculations.
        
        Args:
            natoms (int): Number of atoms per system
            nbonds (int): Number of bonds per system
            nangles (int): Number of angles per system
            ndihedrals (int): Number of dihedrals per system
            ninversions (int): Number of inversions per system
            npbc (int): Number of PBC shift vectors per system
            nSystems (int): Number of systems to simulate
        """
        self.nSystems = nSystems
        self.natoms = natoms
        self.nbonds = nbonds
        self.nangles = nangles
        self.ndihedrals = ndihedrals
        self.ninversions = ninversions
        self.npbc = npbc
        
        # Use proper OpenCL memory flags (READ_WRITE) instead of numpy data types
        self.check_buf("apos",   natoms * nSystems * 4 * f32sz)  
        self.check_buf("fapos",  natoms * nSystems * 4 * f32sz)  
        self.check_buf("fint",   natoms * nSystems * 4 * f32sz)  
        self.check_buf("atype",  natoms * nSystems * i32sz)     
        self.check_buf("aREQ",   natoms * nSystems * 4 * f32sz)  
        
        # Bond parameters
        self.check_buf("bonds",      nbonds * nSystems * 2 * i32sz)       
        self.check_buf("bondParams", nbonds * nSystems * 4 * f32sz)       
        
        # Angle parameters
        self.check_buf("angles",      nangles * nSystems * 3 * i32sz)     
        self.check_buf("angleParams", nangles * nSystems * 4 * f32sz)     
        self.check_buf("angAtoms",    nangles * nSystems * 4 * i32sz)     
        self.check_buf("angNgs",      nangles * nSystems * 4 * i32sz)     
        
        # Dihedral parameters
        self.check_buf("dihedrals",       ndihedrals  * nSystems * 4 * i32sz)     
        self.check_buf("dihedralParams",  ndihedrals  * nSystems * 4 * f32sz)     
        self.check_buf("dihAtoms",        ndihedrals  * nSystems * 4 * i32sz)     
        self.check_buf("dihNgs",          ndihedrals  * nSystems * 4 * i32sz)     
        self.check_buf("inversions",      ninversions * nSystems * 4 * i32sz)     
        self.check_buf("inversionParams", ninversions * nSystems * 4 * f32sz)     
        self.check_buf("invAtoms",        ninversions * nSystems * 4 * i32sz)     
        self.check_buf("invNgs",          ninversions * nSystems * 4 * i32sz)     
        self.check_buf("neighBs",         natoms      * nSystems * 4 * i32sz)     
        self.check_buf("pbcShifts",       npbc        * nSystems * 4 * f32sz)     
        self.check_buf("energies",        5           * nSystems * f32sz)       
        self.check_buf("lvec",            9           * nSystems * f32sz)       
        self.check_buf("params",          10                     * i32sz)       
        self.args_setup = False
        
        print(f"UFF buffers allocated for {nSystems} systems with {natoms} atoms each")
    
    def set_a2f_map_size(self, size):
        """
        Sets the size of the atom-to-force map.
        
        Args:
            size (int): Total number of references in the a2f map
        """
        self.a2f_map_size = size
        
        # Allocate a2f map buffers
        self.check_buf("a2f_offsets", self.natoms * self.nSystems * i32sz)  
        self.check_buf("a2f_counts", self.natoms * self.nSystems * i32sz)   
        self.check_buf("a2f_indices", size * i32sz)                         
        
        print(f"A2F map size set to {size}")
    
    
    def upload_positions(self, positions, iSys=0):
        """
        Uploads atom positions for one system to the GPU.
        
        Args:
            positions (np.ndarray): Array of atom positions (natoms x 3)
            iSys (int): System index
        """
        if positions.shape[0] != self.natoms:
            raise ValueError(f"Expected {self.natoms} atoms, got {positions.shape[0]}")
        
        # Convert to float32 if needed
        if positions.dtype != np.float32:
            positions = positions.astype(np.float32)
        
        # Ensure positions is a 2D array
        if len(positions.shape) == 1:
            positions = positions.reshape(-1, 3)
        
        # Create a padded array with 4 components per position (xyz + padding)
        padded_positions = np.zeros((self.natoms, 4), dtype=np.float32)
        padded_positions[:, :3] = positions
        
        # Calculate offset for this system
        offset = iSys * self.natoms * 4
        
        # Upload to GPU
        cl.enqueue_copy(self.queue, self.buffer_dict["apos"], padded_positions.flatten(), device_offset=offset * f32sz)  
    
    def upload_topology_params(self, uff_data, iSys=0):
        """
        Uploads static topology and parameter data for one system to the GPU buffers.
        
        Args:
            uff_data (dict): Dictionary containing UFF parameter arrays
            iSys (int): System index
        """
        # Calculate offsets for this system
        atom_offset      = iSys * self.natoms
        bond_offset      = iSys * self.nbonds
        angle_offset     = iSys * self.nangles
        dihedral_offset  = iSys * self.ndihedrals
        inversion_offset = iSys * self.ninversions
        
        def safe_copy(buffer_name, data_key, data_type, offset, byte_size):
            #if data_key in uff_data and len(uff_data[data_key]) > 0:
            #cl.enqueue_copy(self.queue, self.buffer_dict[buffer_name], uff_data[data_key].astype(data_type).flatten(), device_offset=offset * byte_size)
            cl.enqueue_copy(self.queue, self.buffer_dict[buffer_name], uff_data[data_key], device_offset=offset * byte_size)

        # Upload atom types and parameters
        safe_copy("atype",          uff_data["atype"]          .astype(np.int32),             device_offset=atom_offset      * i32sz)
        safe_copy("aREQ",           uff_data["aREQ"]           .astype(np.float32).flatten(), device_offset=atom_offset      * 4*f32sz)
        safe_copy("bonds",          uff_data["bonds"]          .astype(np.int32).flatten(),   device_offset=bond_offset      * 2*i32sz)
        safe_copy("bondParams",     uff_data["bondParams"]     .astype(np.float32).flatten(), device_offset=bond_offset      * 4*f32sz)
        safe_copy("angles",         uff_data["angles"]         .astype(np.int32).flatten(),   device_offset=angle_offset     * 3*i32sz)
        safe_copy("angleParams",    uff_data["angleParams"]    .astype(np.float32).flatten(), device_offset=angle_offset     * 4*f32sz)
        

                    
        # Copy dihedrals and parameters - skip if empty
        safe_copy("dihedrals", "dihedrals", np.int32, dihedral_offset, 4*i32sz)
        safe_copy("dihedralParams", "dihedralParams", np.float32, dihedral_offset, 4*f32sz)
        
        # Copy inversions and parameters - skip if empty
        safe_copy("inversions", "inversions", np.int32, inversion_offset, 4*i32sz)
        safe_copy("inversionParams", "inversionParams", np.float32, inversion_offset, 4*f32sz)
        
        # Copy neighbor lists, lattice vectors, and PBC shifts
        safe_copy("neighBs", "neighBs", np.int32, atom_offset, 4*i32sz)
        
        # Copy lattice vectors if present
        if "lvec" in uff_data:
            safe_copy("lvec", "lvec", np.float32, iSys, 9*f32sz)
        
        # Copy PBC shifts if present
        if "pbcShifts" in uff_data:
            safe_copy("pbcShifts", "pbcShifts", np.float32, iSys * self.npbc, 4*f32sz)
        
        # Copy angle atoms and neighbors if angles are present
        if len(uff_data["angles"]) > 0:
            safe_copy("angAtoms", "angAtoms", np.int32, angle_offset, 4*i32sz)
            safe_copy("angNgs", "angNgs", np.int32, angle_offset, 4*i32sz)
        
        # Copy dihedral (torsion) atoms and neighbors if dihedrals are present
        if len(uff_data["dihedrals"]) > 0:
            safe_copy("dihAtoms", "dihAtoms", np.int32, dihedral_offset, 4*i32sz)
            safe_copy("dihNgs", "dihNgs", np.int32, dihedral_offset, 4*i32sz)
        # Copy inversion atoms and neighbors if inversions are present
        if len(uff_data["inversions"]) > 0:
            safe_copy("invAtoms", "invAtoms", np.int32, inversion_offset, 4*i32sz)
            safe_copy("invNgs", "invNgs", np.int32, inversion_offset, 4*i32sz)
        
        # Upload a2f map if present
        # Calculate the indices offset based on a2f map position in the buffer
        indices_offset = iSys * self.a2f_map_size if hasattr(self, 'a2f_map_size') else 0
        
        cl.enqueue_copy(self.queue, self.buffer_dict["a2f_offsets"], uff_data["a2f_offsets"].astype(np.int32), device_offset=atom_offset    * i32sz)  
        cl.enqueue_copy(self.queue, self.buffer_dict["a2f_counts"],  uff_data["a2f_counts"].astype(np.int32),  device_offset=atom_offset    * i32sz)  
        cl.enqueue_copy(self.queue, self.buffer_dict["a2f_indices"], uff_data["a2f_indices"].astype(np.int32), device_offset=indices_offset * i32sz)  
    
    def upload_params(self, params):
        """
        Uploads general parameters to the GPU.
        
        Args:
            params (np.ndarray): Array of parameters
        """
        cl.enqueue_copy(self.queue, self.buffer_dict["params"], params.astype(np.int32))
    
    def prepare_kernel_args(self):
        """
        Prepares kernel arguments for all UFF kernels.
        """
        if self.args_setup:
            return
        
        # Initialize kernel parameters if not already done
        if not hasattr(self, 'kernel_params'):
            self.kernel_params = {}
            # Set basic parameters like natoms, nbonds, etc.
            self.kernel_params['natoms'] = np.int32(self.natoms)
            self.kernel_params['nbonds'] = np.int32(self.nbonds)
            self.kernel_params['nangles'] = np.int32(self.nangles)
            self.kernel_params['ndihedrals'] = np.int32(self.ndihedrals)
            self.kernel_params['ninversions'] = np.int32(self.ninversions)
            self.kernel_params['nSystems'] = np.int32(self.nSystems)
            self.kernel_params['bSubtractVdW'] = np.int32(0) # Default value
        
        # Use OpenCLBase's functionality for generating kernel arguments
        if not hasattr(self, 'kernelheaders') or not self.kernelheaders:
            # If kernel headers are not set, extract them from prg source
            self.kernelheaders = self.extract_kernel_headers(self.prg.get_info(cl.prg_info.SOURCE))
        
        self.kernel_args = {}
        for kernel_name in self.kernelheaders:
            self.kernel_args[kernel_name] = self.generate_kernel_args(kernel_name)
        
        self.args_setup = True
    
    # def _extract_kernel_headers(self, source_code):
    #     """
    #     Extracts kernel function signatures from OpenCL source code.
        
    #     Args:
    #         source_code (str): The entire OpenCL source code as a string.
        
    #     Returns:
    #         dict: Dictionary mapping kernel names to their header signatures
    #     """
    #     import re
        
    #     # Regular expression to match kernel function declarations
    #     kernel_pattern = r'__kernel\s+void\s+(\w+)\s*\((.*?)\)'
        
    #     # Find all kernel declarations in the source code
    #     kernel_matches = re.finditer(kernel_pattern, source_code, re.DOTALL)
        
    #     kernel_headers = {}
    #     for match in kernel_matches:
    #         kernel_name = match.group(1)
    #         kernel_args = match.group(2)
    #         kernel_headers[kernel_name] = kernel_args
        
    #     return kernel_headers
        
    def run_eval_step(self, bClearForce=True):
        """
        Executes one step of UFF evaluation kernels.
        
        Args:
            bClearForce (bool): Whether to clear forces before evaluation
        
        Returns:
            float: Total energy
        """
        if not self.args_setup:
            self.prepare_kernel_args()
        
        # Clear forces if requested
        if bClearForce:
            self.prg.clearForces(self.queue, (self.natoms * self.nSystems,), None,  self.buffer_dict["fapos"])
        self.prg.evalBonds(self.queue, (self.nbonds * self.nSystems,), None,  *self.kernel_args["evalBonds"])
        if self.nangles > 0:      self.prg.evalAngles    (self.queue, (self.nangles * self.nSystems,), None,  *self.kernel_args["evalAngles"])
        if self.ndihedrals > 0:   self.prg.evalDihedrals (self.queue, (self.ndihedrals * self.nSystems,), None,  *self.kernel_args["evalDihedrals"])
        if self.ninversions > 0:  self.prg.evalInversions(self.queue, (self.ninversions * self.nSystems,), None, *self.kernel_args["evalInversions"])
        if self.kernel_params.get("bNonBonded", False): self.prg.evalNonBonded(self.queue, (self.natoms * self.nSystems,), None,  *self.kernel_args["evalNonBonded"])
        self.prg.sumEnergies(self.queue, (self.nSystems,), None,   *self.kernel_args["sumEnergies"])
        energies = np.zeros(5 * self.nSystems, dtype=np.float32)
        cl.enqueue_copy(self.queue, energies, self.buffer_dict["energies"])
        
        # Return total energy (last component for each system)
        return energies[4::5]
    
    def get_forces(self, iSys=None):
        """
        Downloads forces from GPU.
        
        Args:
            iSys (int, optional): System index. If None, download all systems.
        
        Returns:
            np.ndarray: Forces array
        """
        if iSys is None:
            # Download all forces
            forces = np.zeros(self.natoms * self.nSystems * 4, dtype=np.float32)
            cl.enqueue_copy(self.queue, forces, self.buffer_dict["fapos"])
            
            # Reshape to (nSystems, natoms, 4) and remove padding
            forces = forces.reshape(self.nSystems, self.natoms, 4)[:, :, :3]
            return forces
        else:
            # Download forces for one system
            forces = np.zeros(self.natoms * 4, dtype=np.float32)
            offset = iSys * self.natoms * 4
            cl.enqueue_copy(self.queue, forces, self.buffer_dict["fapos"], 
                           device_offset=offset * f32sz)  
    
            # Reshape to (natoms, 4) and remove padding
            forces = forces.reshape(self.natoms, 4)[:, :3]
            return forces
    
    def get_energies(self):
        """
        Downloads energy contributions from GPU.
        
        Returns:
            dict: Dictionary of energy contributions
        """
        energies = np.zeros(5 * self.nSystems, dtype=np.float32)
        cl.enqueue_copy(self.queue, energies, self.buffer_dict["energies"])
        
        # Reshape to (nSystems, 5)
        energies = energies.reshape(self.nSystems, 5)
        
        # Create dictionary of energy components
        energy_dict = {
            "bond": energies[:, 0],
            "angle": energies[:, 1],
            "dihedral": energies[:, 2],
            "inversion": energies[:, 3],
            "total": energies[:, 4]
        }
        
        return energy_dict
    

    def mapAtomInteractions(self, natoms, ndihedrals, ninversions, nangles):
        """
        Maps atom interactions to force pieces using a buckets structure.
        Similar to UFF::mapAtomInteractions in C++.
        
        Args:
            natoms (int): Number of atoms
            ndihedrals (int): Number of dihedrals
            ninversions (int): Number of inversions
            nangles (int): Number of angles
        
        Returns:
            tuple: (a2f_offsets, a2f_counts, a2f_indices) arrays for GPU upload
        """
        # Initialize arrays
        a2f_counts = np.zeros(natoms, dtype=np.int32)
        
        # Count interactions per atom
        # For dihedrals (4 atoms per dihedral)
        for i in range(ndihedrals):
            a2f_counts[self.dihedrals[i, 0]] += 1
            a2f_counts[self.dihedrals[i, 1]] += 1
            a2f_counts[self.dihedrals[i, 2]] += 1
            a2f_counts[self.dihedrals[i, 3]] += 1
        
        # For inversions (4 atoms per inversion)
        for i in range(ninversions):
            a2f_counts[self.inversions[i, 0]] += 1
            a2f_counts[self.inversions[i, 1]] += 1
            a2f_counts[self.inversions[i, 2]] += 1
            a2f_counts[self.inversions[i, 3]] += 1
        
        # For angles (3 atoms per angle)
        for i in range(nangles):
            a2f_counts[self.angles[i, 0]] += 1
            a2f_counts[self.angles[i, 1]] += 1
            a2f_counts[self.angles[i, 2]] += 1
        
        # Calculate total size and offsets
        total_refs = np.sum(a2f_counts)
        a2f_offsets = np.zeros(natoms, dtype=np.int32)
        
        # Calculate offsets
        offset = 0
        for i in range(natoms):
            a2f_offsets[i] = offset
            offset += a2f_counts[i]
        
        # Reset counts for filling indices
        a2f_counts_temp = np.zeros(natoms, dtype=np.int32)
        a2f_indices = np.zeros(total_refs, dtype=np.int32)
        
        # Fill indices for dihedrals
        for i in range(ndihedrals):
            for j in range(4):
                atom_idx                   = self.dihedrals[i, j]
                offset                     = a2f_offsets[atom_idx] + a2f_counts_temp[atom_idx]
                a2f_indices[offset]        = i
                a2f_counts_temp[atom_idx] += 1
        
        # Fill indices for inversions
        for i in range(ninversions):
            for j in range(4):
                atom_idx                   = self.inversions[i, j]
                offset                     = a2f_offsets[atom_idx] + a2f_counts_temp[atom_idx]
                a2f_indices[offset]        = i + ndihedrals  # Offset by ndihedrals
                a2f_counts_temp[atom_idx] += 1
        
        # Fill indices for angles
        for i in range(nangles):
            for j in range(3):
                atom_idx                   = self.angles[i, j]
                offset                     = a2f_offsets[atom_idx] + a2f_counts_temp[atom_idx]
                a2f_indices[offset]        = i + ndihedrals + ninversions  # Offset by ndihedrals + ninversions
                a2f_counts_temp[atom_idx] += 1
        
        return a2f_offsets, a2f_counts, a2f_indices
    
    def bakeDihedralNeighs(self, dihedrals, natoms):
        """
        Prepares dihedral neighbor information.
        Similar to UFF::bakeDihedralNeighs in C++.
        
        Args:
            dihedrals (np.ndarray): Array of dihedral definitions
            natoms (int): Number of atoms
        
        Returns:
            tuple: (dihAtoms, dihNgs) arrays for GPU upload
        """
        ndihedrals = len(dihedrals)
        dihAtoms   = np.zeros((ndihedrals, 4), dtype=np.int32)
        dihNgs     = np.zeros((ndihedrals, 4), dtype=np.int32)
        
        # Initialize neighbor arrays to -1 (no neighbor)
        dihNgs.fill(-1)
        
        # Copy atom indices
        for i in range(ndihedrals):  dihAtoms[i] = dihedrals[i]
        
        # Build neighbor lists
        for i in range(ndihedrals):
            a1, a2, a3, a4 = dihedrals[i]
            # Find other dihedrals sharing atoms with this one
            for j in range(ndihedrals):
                if i == j:  continue
                b1, b2, b3, b4 = dihedrals[j]
                # Check for shared atoms and assign neighbors
                if ((a1 == b1) or (a1 == b2) or (a1 == b3) or (a1 == b4)) and (dihNgs[i, 0] == -1):    dihNgs[i, 0] = j
                if ((a2 == b1) or (a2 == b2) or (a2 == b3) or (a2 == b4)) and (dihNgs[i, 1] == -1):    dihNgs[i, 1] = j
                if ((a3 == b1) or (a3 == b2) or (a3 == b3) or (a3 == b4)) and (dihNgs[i, 2] == -1):    dihNgs[i, 2] = j
                if ((a4 == b1) or (a4 == b2) or (a4 == b3) or (a4 == b4)) and (dihNgs[i, 3] == -1):    dihNgs[i, 3] = j
        
        return dihAtoms, dihNgs
    
    def bakeAngleNeighs(self, angles, natoms):
        """
        Prepares angle neighbor information.
        Similar to UFF::bakeAngleNeighs in C++.
        
        Args:
            angles (np.ndarray): Array of angle definitions
            natoms (int): Number of atoms
        
        Returns:
            tuple: (angAtoms, angNgs) arrays for GPU upload
        """
        nangles = len(angles)
        angAtoms = np.zeros((nangles, 4), dtype=np.int32)  # 4th element is padding
        angNgs = np.zeros((nangles, 4), dtype=np.int32)
        
        # Initialize neighbor arrays to -1 (no neighbor)
        angNgs.fill(-1)
        
        # Copy atom indices and pad with -1
        for i in range(nangles):
            angAtoms[i, :3] = angles[i]
            angAtoms[i, 3] = -1  # Padding
        
        # Build neighbor lists
        for i in range(nangles):
            a1, a2, a3 = angles[i][:3]
            
            # Find other angles sharing atoms with this one
            for j in range(nangles):
                if i == j: continue
                b1, b2, b3 = angles[j][:3]
                # Check for shared atoms and assign neighbors
                if ((a1 == b1) or (a1 == b2) or (a1 == b3)) and (angNgs[i, 0] == -1): angNgs[i, 0] = j
                if ((a2 == b1) or (a2 == b2) or (a2 == b3)) and (angNgs[i, 1] == -1): angNgs[i, 1] = j
                if ((a3 == b1) or (a3 == b2) or (a3 == b3)) and (angNgs[i, 2] == -1): angNgs[i, 2] = j
        
        return angAtoms, angNgs
    
    def bakeInversionNeighs(self, inversions, natoms):
        """
        Prepares inversion neighbor information.
        Similar to UFF::bakeInversionNeighs in C++.
        
        Args:
            inversions (np.ndarray): Array of inversion definitions
            natoms (int): Number of atoms
        
        Returns:
            tuple: (invAtoms, invNgs) arrays for GPU upload
        """
        ninversions = len(inversions)
        invAtoms = np.zeros((ninversions, 4), dtype=np.int32)
        invNgs = np.zeros((ninversions, 4), dtype=np.int32)
        
        # Initialize neighbor arrays to -1 (no neighbor)
        invNgs.fill(-1)
        
        # Copy atom indices
        for i in range(ninversions):
            invAtoms[i] = inversions[i]
        
        # Build neighbor lists
        for i in range(ninversions):
            a1, a2, a3, a4 = inversions[i]
            
            # Find other inversions sharing atoms with this one
            for j in range(ninversions):
                if i == j:
                    continue
                b1, b2, b3, b4 = inversions[j]
                # Check for shared atoms and assign neighbors
                if ((a1 == b1) or (a1 == b2) or (a1 == b3) or (a1 == b4)) and (invNgs[i, 0] == -1): invNgs[i, 0] = j
                if ((a2 == b1) or (a2 == b2) or (a2 == b3) or (a2 == b4)) and (invNgs[i, 1] == -1): invNgs[i, 1] = j
                if ((a3 == b1) or (a3 == b2) or (a3 == b3) or (a3 == b4)) and (invNgs[i, 2] == -1): invNgs[i, 2] = j
                if ((a4 == b1) or (a4 == b2) or (a4 == b3) or (a4 == b4)) and (invNgs[i, 3] == -1): invNgs[i, 3] = j
        return invAtoms, invNgs

    def toUFF(self, mol, bRealloc=True):
        """
        Converts molecular structure to UFF representation.
        Similar to builder.toUFF() in C++.
        
        Args:
            mol (object): Molecular system with atoms, bonds, etc.
            bRealloc (bool): Flag to reallocate buffers
        
        Returns:
            dict: UFF data dictionary for upload
        """
        # Make sure element_types and atom_types are loaded before processing
        if not hasattr(self, 'element_types') or not hasattr(self, 'atom_types'):
            # Load element and atom types from data files
            base_path = os.path.dirname(os.path.abspath(__file__))
            data_path = os.path.join(base_path, "../../cpp/common_resources/")
            self.element_types = MMparams.read_element_types(os.path.join(data_path, 'ElementTypes.dat'))
            self.atom_types = MMparams.read_atom_types(os.path.join(data_path, 'AtomTypes.dat'), self.element_types)
        
        # Generate or retrieve required attributes for UFF calculation
        # Number of pi electrons for each atom
        npi_list = getattr(mol, 'npi_list', [self.atom_types[name].npi if name in self.atom_types else 0 for name in mol.enames])
        # Number of electron pairs for each atom
        nep_list = getattr(mol, 'nep_list', [self.atom_types[name].nepair if name in self.atom_types else 0 for name in mol.enames])
        # Which atoms are nodes in the molecular graph (non-hydrogen atoms usually)
        capping_atoms = ['H', 'F', 'Cl', 'Br', 'I']  # Terminal atoms typically not considered nodes
        isNode = getattr(mol, 'isNode', [0 if self.atom_types[name].element_name in capping_atoms else 1 
                                       for name in mol.enames if name in self.atom_types])
        # REQs (radius, epsilon, charge) parameters
        REQs = getattr(mol, 'REQs', MMparams.generate_REQs_from_atom_types(mol, self.atom_types))
        
        # Extract basic molecular information
        natoms = len(mol.atypes)
        nbonds = len(mol.bonds)
        
        # Count angles, dihedrals, and inversions
        angles = []
        dihedrals = []
        inversions = []
        
        # Build neighbor list
        neighs = [[] for _ in range(natoms)]
        for i, (a, b) in enumerate(mol.bonds):
            neighs[a].append(b)
            neighs[b].append(a)
        
        # Find angles (3 connected atoms)
        for i in range(natoms):
            for j in neighs[i]:
                for k in neighs[j]:
                    if k != i:
                        # Ensure we don't add the same angle twice
                        if i < k:
                            angles.append((i, j, k))
        
        # Find dihedrals (4 connected atoms)
        for i in range(natoms):
            for j in neighs[i]:
                for k in neighs[j]:
                    if k != i:
                        for l in neighs[k]:
                            if l != j:
                                # Ensure we don't add the same dihedral twice
                                if i < l:
                                    dihedrals.append((i, j, k, l))
        
        # Find inversions (4 connected atoms)
        for i in range(natoms):
            for j in neighs[i]:
                for k in neighs[j]:
                    if k != i:
                        for l in neighs[k]:
                            if l != j:
                                for m in neighs[l]:
                                    if m != k:
                                        # Ensure we don't add the same inversion twice
                                        if i < m:
                                            inversions.append((i, j, k, l, m))
        
        # Allocate buffers if needed
        if bRealloc:
            self.realloc_buffers(natoms, nbonds, len(angles), len(dihedrals), len(inversions), 0, 1)
        
        # Convert atom types to UFF types
        atype = np.zeros(natoms, dtype=np.int32)
        aREQ = np.zeros(natoms * 4, dtype=np.float32)
        for i, at in enumerate(mol.atypes):
            uff_type        = self._get_uff_type(at)
            atype[i]        = uff_type
            aREQ[i*4:i*4+4] = self.get_uff_params(uff_type)
        
        # Convert bonds to UFF format
        bonds = np.zeros(nbonds * 2, dtype=np.int32)
        bondParams = np.zeros(nbonds * 4, dtype=np.float32)
        for i, (a, b) in enumerate(mol.bonds):
            bonds[i*2:i*2+2]      = [a, b]
            # get_bond_params returns [bond_length, bond_force] (shape 2), but we need shape 4
            bond_params = self.get_bond_params(uff_type)
            bondParams[i*4]   = bond_params[0]  # bond length
            bondParams[i*4+1] = bond_params[1]  # bond force
            bondParams[i*4+2] = 0.0  # Placeholder for any additional parameters
            bondParams[i*4+3] = 0.0  # Placeholder for any additional parameters
        
        # Convert angles to UFF format
        angles = np.array(angles, dtype=np.int32)
        angleParams = np.zeros(len(angles) * 4, dtype=np.float32)
        for i, (a, b, c) in enumerate(angles):
            # _get_angle_params returns [angle, force] (shape 2), but we need shape 4
            angle_params = self._get_angle_params(uff_type)
            angleParams[i*4]   = angle_params[0]  # equilibrium angle
            angleParams[i*4+1] = angle_params[1]  # force constant
            angleParams[i*4+2] = 0.0  # Placeholder for any additional parameters
            angleParams[i*4+3] = 0.0  # Placeholder for any additional parameters
        
        # Convert dihedrals to UFF format
        dihedrals = np.array(dihedrals, dtype=np.int32)
        dihedralParams = np.zeros(len(dihedrals) * 4, dtype=np.float32)
        for i, (a, b, c, d) in enumerate(dihedrals):
            # _get_dihedral_params returns [angle, force] (shape 2), but we need shape 4
            dihedral_params = self._get_dihedral_params(uff_type)
            dihedralParams[i*4]   = dihedral_params[0]  # equilibrium dihedral angle
            dihedralParams[i*4+1] = dihedral_params[1]  # force constant
            dihedralParams[i*4+2] = 0.0  # Placeholder for any additional parameters
            dihedralParams[i*4+3] = 0.0  # Placeholder for any additional parameters
        
        # Convert inversions to UFF format
        inversions = np.array(inversions, dtype=np.int32)
        inversionParams = np.zeros(len(inversions) * 4, dtype=np.float32)
        for i, (a, b, c, d, e) in enumerate(inversions):
            # _get_inversion_params returns [angle, force] (shape 2), but we need shape 4
            inversion_params = self._get_inversion_params(uff_type)
            inversionParams[i*4]   = inversion_params[0]  # equilibrium angle
            inversionParams[i*4+1] = inversion_params[1]  # force constant
            inversionParams[i*4+2] = 0.0  # Placeholder for any additional parameters
            inversionParams[i*4+3] = 0.0  # Placeholder for any additional parameters
        
        # Create UFF data dictionary
        uff_data = {
            "atype": atype,
            "aREQ": aREQ,
            "bonds": bonds,
            "bondParams": bondParams,
            "angles": angles,
            "angleParams": angleParams,
            "dihedrals": dihedrals,
            "dihedralParams": dihedralParams,
            "inversions": inversions,
            "inversionParams": inversionParams,
            "neighBs": neighs
        }
        
        return uff_data
    

    
    def get_uff_params(self, uff_type):
        """
        Retrieves UFF parameters for a given UFF type.
        
        Args:
            uff_type (int): UFF type index
        
        Returns:
            np.ndarray: Array of UFF parameters [REQ, epsilon, sigma, mass]
        """
        # Use element_types if available
        if hasattr(self, 'element_types'):
            for et_name, et in self.element_types.items():
                if et.iZ == uff_type:
                    return np.array([et.RvdW, et.EvdW, 0.0, et.iZ * 2.0], dtype=np.float32)
        
        # Fallback to basic parameters
        uff_params = {
            1: [1.0, 0.01, 0.0, 1.0],  # Hydrogen
            6: [1.7, 0.05, 0.0, 12.0],  # Carbon
            7: [1.55, 0.07, 0.0, 14.0],  # Nitrogen
            8: [1.52, 0.06, 0.0, 16.0],  # Oxygen
            16: [1.8, 0.05, 0.0, 32.0]  # Sulfur
        }
        return np.array(uff_params.get(uff_type, [0.0, 0.0, 0.0, 0.0]), dtype=np.float32)
    
    def get_bond_params(self, uff_type):
        """
        Retrieves bond parameters for a given UFF type.
        
        Args:
            uff_type (int): UFF type index
        
        Returns:
            np.ndarray: Array of bond parameters [bond length, bond force constant]
        """
        # ElementType should already be loaded before this is called
        for et_name, et in self.element_types.items():
            if et.iZ == uff_type:
                # UFF bond length is related to covalent radius
                bond_length = et.Rcov
                # Bond force constant scales with atomic number
                bond_force = 350.0 + (et.iZ * 2.0)
                return np.array([bond_length, bond_force], dtype=np.float32)
                
        # If we get here, the element type was not found - this is an error
        raise ValueError(f"Element {uff_type} not found in element_types")
    

    
    def _get_dihedral_params(self, uff_type):
        """
        Retrieves dihedral parameters for a given UFF type.
        
        Args:
        
            uff_type (int): UFF type index
        
        Returns:
            np.ndarray: Array of dihedral parameters [dihedral angle, dihedral force constant]
        """
        # AtomType should already be loaded before this is called
        for at_name, at in self.atom_types.items():
            if at.iZ == uff_type:
                # Determine dihedral parameters based on hybridization
                angle = 0.0  # Default barrier angle (pi)
                force = at.iZ * 0.5  # Base force constant
                
                # Stronger barriers for pi-bonded atoms
                if at.npi > 0:
                    force *= 1.5
                    
                return np.array([angle, force], dtype=np.float32)
        
        # If we get here, the atom type was not found - this is an error
        raise ValueError(f"UFF type {uff_type} not found in atom_types")


    def _get_angle_params(self, uff_type):
        """
        Retrieves angle parameters for a given UFF type.
        
        Args:
            uff_type (int): UFF type index
        
        Returns:
            np.ndarray: Array of angle parameters [angle, angle force constant]
        """
        # AtomType should already be loaded before this is called
        for at_name, at in self.atom_types.items():
            if at.iZ == uff_type:
                # Get angle based on hybridization determined by npi
                angle = 109.5  # Default tetrahedral sp3
                if at.npi == 1: angle = 120.0  # sp2
                if at.npi == 2: angle = 180.0  # sp
                
                # Force constant depends on element properties
                force = 50.0 + (10.0 * at.iZ / 8.0)
                return np.array([angle, force], dtype=np.float32)
        
        # If we get here, the atom type was not found - this is an error
        raise ValueError(f"UFF type {uff_type} not found in atom_types")



    def _get_inversion_params(self, uff_type):
        """
        Retrieves inversion parameters for a given UFF type.
        
        Args:
            uff_type (int): UFF type index
    
        Returns:
            np.ndarray: Array of inversion parameters [inversion angle, inversion force constant]
        """
        # AtomType should be already loaded before this is called
        # Lookup atom type by element number
        for at_name, at in self.atom_types.items():
            if at.iZ == uff_type:
                # Inversion parameters depend on element and hybridization
                angle = 0.0  # Default angle
                force = 0.1   # Default force
                
                # Planar elements (sp2 hybridized) have higher barriers
                if at.iZ in [6, 7, 15, 33, 51, 83] and at.npi == 1:
                    force = 4.0 + (at.iZ / 10.0)
                    
                return np.array([angle, force], dtype=np.float32)
                
        # If we get here, the atom type was not found - this is an error
        raise ValueError(f"UFF type {uff_type} not found in atom_types")

    def _get_uff_type(self, atom):
        """
        Get the UFF type for an atom.
        
        Args:
            atom: The atom to get the UFF type for - can be an atom object with Z attribute
                 or directly the atomic number (int or numpy.int32)
            
        Returns:
            int: The UFF atomic type identifier
        """
        # Check if atom is already numeric (int or numpy.int32) or an object with Z attribute
        if hasattr(atom, 'Z'):
            element_num = atom.Z
        else:
            # Assume it's directly the atomic number
            element_num = int(atom)  # Convert to int in case it's numpy.int32
        
        # Verify element exists in element_types
        for et_name, et in self.element_types.items():
            if et.iZ == element_num:
                return element_num
        
        # If we get here, the element type was not found - this is an error
        raise ValueError(f"Element {element_num} not found in element_types")
