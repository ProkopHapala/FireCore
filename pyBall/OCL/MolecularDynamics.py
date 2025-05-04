import sys
import os
import numpy as np

import pyopencl as cl
import pyopencl.array as cl_array
import pyopencl.cltypes as cltypes
import matplotlib.pyplot as plt
import time

from . import clUtils as clu
from .MMFF import MMFF

verbose=False

def pack(iSys, source_array, target_buffer, queue):
    """
    Packs data from a source NumPy array into a target OpenCL buffer at the specified system index.
    """
    # Calculate the offset
    offset = iSys * source_array.size * source_array.dtype.itemsize
    cl.enqueue_copy(queue, target_buffer, source_array, offset=offset)

def copy(source, target, queue, iSys):
    """
    Copies data from source array to target OpenCL buffer.
    """
    pack(iSys, source, target, queue)

def copy_add(source, source_add, target, offset, queue):
    """
    Copies data from source array and adds source_add array to the target buffer.
    """
    # This function needs to be implemented based on specific requirements
    pass

def mat3_to_cl(mat3_np):
    """
    Converts a 3x3 NumPy matrix to a format suitable for OpenCL.
    """
    return mat3_np.flatten().astype(np.float32)

def vec3_to_cl(vec3_np):
    """
    Converts a 3-component vector to a 4-component format for OpenCL.
    """
    return np.append(vec3_np, 0.0).astype(np.float32)

class MolecularDynamics:

    kernelheader_getMMFFf4 = """
__kernel void getMMFFf4(
    const int4 nDOFs,               // 1   (nAtoms,nnode) dimensions of the system
    // Dynamical
    __global float4*  apos,         // 2  [natoms]     positions of atoms (including node atoms [0:nnode] and capping atoms [nnode:natoms] and pi-orbitals [natoms:natoms+nnode] )
    __global float4*  aforce,       // 3  [natoms]     forces on    atoms (just node atoms are evaluated)
    __global float4*  fneigh,       // 4  [nnode*4*2]  recoil forces on neighbors (and pi-orbitals)
    // parameters
    __global int4*    neighs,       // 5  [nnode]  neighboring atoms
    __global int4*    neighCell,    // 5  [nnode]  neighboring atom  cell index
    __global float4*  REQs,         // 6  [natoms] non-boding parametes {R0,E0,Q} i.e. R0: van der Waals radii, E0: well depth and partial charge, Q: partial charge 
    __global float4*  apars,        // 7  [nnode]  per atom forcefield parametrs {c0ss,Kss,c0sp}, i.e. c0ss: cos(equlibrium angle/2) for sigma-sigma; Kss: stiffness of sigma-sigma angle; c0sp: is cos(equlibrium angle) for sigma-pi
    __global float4*  bLs,          // 8  [nnode]  bond length    between node and each neighbor
    __global float4*  bKs,          // 9  [nnode]  bond stiffness between node and each neighbor
    __global float4*  Ksp,          // 10 [nnode]  stiffness of pi-alignment for each neighbor     (only node atoms have pi-pi alignemnt interaction)
    __global float4*  Kpp,          // 11 [nnode]  stiffness of pi-planarization for each neighbor (only node atoms have pi-pi alignemnt interaction)
    __global cl_Mat3* lvecs,        // 12 lattice vectors         for each system
    __global cl_Mat3* ilvecs,       // 13 inverse lattice vectors for each system
    __global float4*  pbc_shifts,
    const int npbc,
    const int bSubtractVdW
){
    """

    kernelheader_updateAtomsMMFFf4 = """
__kernel void updateAtomsMMFFf4(
    const int4        n,            // 1 // (natoms,nnode) dimensions of the system
    __global float4*  apos,         // 2 // positions of atoms  (including node atoms [0:nnode] and capping atoms [nnode:natoms] and pi-orbitals [natoms:natoms+nnode] )
    __global float4*  avel,         // 3 // velocities of atoms 
    __global float4*  aforce,       // 4 // forces on atoms
    __global float4*  cvf,          // 5 // damping coefficients for velocity and force
    __global float4*  fneigh,       // 6 // recoil forces on neighbors (and pi-orbitals)
    __global int4*    bkNeighs,     // 7 // back neighbors indices (for recoil forces)
    __global float4*  constr,       // 8 // constraints (x,y,z,K) for each atom
    __global float4*  constrK,      // 9 // constraints stiffness (kx,ky,kz,?) for each atom
    __global float4*  MDparams,     // 10 // MD parameters (dt,damp,Flimit)
    __global float4*  TDrives,      // 11 // Thermal driving (T,gamma_damp,seed,?)
    __global cl_Mat3* bboxes,       // 12 // bounding box (xmin,ymin,zmin)(xmax,ymax,zmax)(kx,ky,kz)
    __global int*     sysneighs,    // 13 // // for each system contains array int[nMaxSysNeighs] of nearby other systems
    __global float4*  sysbonds      // 14 // // contains parameters of bonds (constrains) with neighbor systems   {Lmin,Lmax,Kpres,Ktens}
){
    """
    
    kernelheader_printOnGPU = """    }
__kernel void printOnGPU(
    const int4        n,            // 1
    const int4        mask,         // 2
    __global float4*  apos,         // 3
    __global float4*  avel,         // 4
    __global float4*  aforce,       // 5
    __global float4*  fneigh,       // 6
    __global int4*    bkNeighs,     // 7
    __global float4*  constr        // 8
){
    """

    kernelheader_cleanForceMMFFf4 = """
__kernel void cleanForceMMFFf4(
    const int4        n,           // 2
    __global float4*  aforce,      // 5
    __global float4*  fneigh       // 6
){
    """

    kernelheader_getNonBond = """
__kernel void getNonBond(
    const int4 ns,                 // 1 // (natoms,nnode) dimensions of the system
    // Dynamical
    __global float4*  apos,        // 2 // positions of atoms  (including node atoms [0:nnode] and capping atoms [nnode:natoms] and pi-orbitals [natoms:natoms+nnode] )
    __global float4*  aforce,      // 3 // forces on atoms
    // Parameters
    __global float4*  REQs,        // 4 // non-bonded parameters (RvdW,EvdW,QvdW,Hbond)
    __global int4*    neighs,      // 5 // neighbors indices      ( to ignore interactions between bonded atoms )
    __global int4*    neighCell,   // 6 // neighbors cell indices ( to know which PBC image should be ignored  due to bond )
    __global cl_Mat3* lvecs,       // 7 // lattice vectors for each system
    const int4        nPBC,        // 8 // number of PBC images in each direction (x,y,z)
    const float4      GFFParams    // 9 // Grid-Force-Field parameters
){
    """

    def __init__(self, nloc=32):
        # Initialize OpenCL context and queue
        self.nloc = nloc
        self.ctx = cl.create_some_context(answers=[0])
        self.queue = cl.CommandQueue(self.ctx)
        
        # Print device info
        clu.get_cl_info(self.ctx.devices[0])
        
        # Grid information placeholders
        self.grid = None  # instance of GridShape, if initialized
        self.gcl = None   # instance of GridCL, if initialized
        
        # Load and compile the OpenCL program - try multiple possible locations
        # Define potential paths where the kernel file might be located
        base_path = os.path.abspath(__file__)
        # /home/prokop/git/FireCore/cpp/common_resources/cl/relax_multi_mini.cl
        rel_path = "../../../cpp/common_resources/cl/relax_multi_mini.cl"
        # form relative path to absolute
        kernel_path = os.path.abspath(os.path.join(base_path, rel_path))
        kernel_found = False
        if os.path.exists(kernel_path):
            with open(kernel_path, 'r') as f:
                self.prg = cl.Program(self.ctx, f.read()).build()
            kernel_found = True
            print(f"Successfully loaded kernel from: {kernel_path}")
        else:
            print(f"Kernel file not found at: {kernel_path}")
            exit(1)

        print("prg.dir()", dir(self.prg) )
        self.getMMFFf4         = self.prg.getMMFFf4
        self.getNonBond        = self.prg.getNonBond
        self.updateAtomsMMFFf4 = self.prg.updateAtomsMMFFf4
        self.cleanForceMMFFf4  = self.prg.cleanForceMMFFf4

        # Initialize other attributes that will be set in realloc
        self.nSystems = 0
        self.mmff_instances = []
        self.buffer_dict = {}



    def realloc(self, nSystems, mmff):
        """
        Reallocate buffers for the given number of systems based on the MMFF template.
        """
        # Store dimensions explicitly to avoid reference issues
        print(f"MMFF dimensions received in realloc: natoms={mmff.natoms}, nvecs={mmff.nvecs}, nnode={mmff.nnode}")
        
        self.nSystems = nSystems
        self.mmff_instances = [mmff] * nSystems  # Assuming all systems use the same MMFF parameters
        
        # Create a new instance for internal use to avoid reference issues
        self.mmff_template = mmff
        
        # Call buffer allocation with explicit dimensions
        self.allocate_cl_buffers(self.mmff_template)

        # Allocate host buffers
        self.allocate_host_buffers(self.mmff_template)
        self.allocate_host_buffers(mmff)

    def allocate_host_buffers(self, mmff ):
        self.atoms  = np.empty(self.nSystems * mmff.nvecs * 4, dtype=np.float32)
        self.forces = np.empty(self.nSystems * mmff.nvecs * 4, dtype=np.float32)

    def allocate_cl_buffers(self, mmff ):
        """
        Allocates OpenCL buffers based on the MMFF template and number of systems.
        """
        nSystems = self.nSystems
        natoms = mmff.natoms
        nvecs  = mmff.nvecs
        nnode  = mmff.nnode
        ncap   = mmff.ncap
        ntors  = mmff.ntors
        nbkng  = nvecs
        nPBC   = mmff.nPBC
        npbc   = mmff.npbc

        self.nDOFs = (natoms,nnode)

        self.natoms = natoms
        self.nvecs  = nvecs
        self.nnode  = nnode
        self.ncap   = ncap
        self.ntors  = ntors
        self.nbkng  = nbkng
        self.nPBC   = nPBC
        self.npbc   = npbc
        
        # Print dimensions for debugging
        print(f"Buffer Allocation Dimensions:\n  nSystems: {nSystems}\n  natoms: {natoms}\n  nvecs: {nvecs}\n  nnode: {nnode}\n  ncap: {ncap}\n  ntors: {ntors}\n  nbkng: {nbkng}")
        
        # Validate dimensions
        if nSystems <= 0 or natoms <= 0 or nvecs <= 0 or nnode <= 0:
            raise ValueError(f"Invalid dimensions for buffer allocation: nSystems={nSystems}, natoms={natoms}, nvecs={nvecs}, nnode={nnode}")
        
        # Calculate buffer sizes and validate
        float_size = np.float32().itemsize
        int_size = np.int32().itemsize
        
        def get_valid_size(dim1, dim2, dim3, elem_size, name):
            size = dim1 * dim2 * dim3 * elem_size
            print(f"  Buffer '{name}' size: {size} bytes ({dim1}*{dim2}*{dim3}*{elem_size})")
            if size <= 0:
                raise ValueError(f"Invalid buffer size for {name}: {size}")
            return size
        
        # Example buffer allocations
        mf = cl.mem_flags
        self.buffer_dict['apos']         = cl.Buffer(self.ctx, mf.READ_WRITE, size=nSystems * nvecs * 3 * np.float32().itemsize)
        self.buffer_dict['aforce']       = cl.Buffer(self.ctx, mf.READ_WRITE, size=nSystems * nvecs * 3 * np.float32().itemsize)
        self.buffer_dict['REQs']         = cl.Buffer(self.ctx, mf.READ_ONLY, size=nSystems * natoms * 4 * np.float32().itemsize)
        self.buffer_dict['neighs']       = cl.Buffer(self.ctx, mf.READ_ONLY, size=nSystems * natoms * 4 * np.int32().itemsize)
        self.buffer_dict['neighCell']    = cl.Buffer(self.ctx, mf.READ_ONLY, size=nSystems * natoms * 4 * np.int32().itemsize)
        self.buffer_dict['bkNeighs']     = cl.Buffer(self.ctx, mf.READ_ONLY, size=nSystems * nvecs * 4 * np.int32().itemsize)
        self.buffer_dict['bkNeighs_new'] = cl.Buffer(self.ctx, mf.READ_ONLY, size=nSystems * nvecs * 4 * np.int32().itemsize)
        self.buffer_dict['avel']         = cl.Buffer(self.ctx, mf.READ_WRITE, size=nSystems * nvecs * 3 * np.float32().itemsize)
        self.buffer_dict['cvf']          = cl.Buffer(self.ctx, mf.READ_WRITE, size=nSystems * nvecs * 3 * np.float32().itemsize)
        self.buffer_dict['fneigh']       = cl.Buffer(self.ctx, mf.READ_WRITE, size=nSystems * nbkng * 3 * np.float32().itemsize)
        self.buffer_dict['fneighpi']     = cl.Buffer(self.ctx, mf.READ_WRITE, size=nSystems * nbkng * 3 * np.float32().itemsize)
        self.buffer_dict['apars']        = cl.Buffer(self.ctx, mf.READ_ONLY, size=nSystems * nnode * 4 * np.float32().itemsize)
        self.buffer_dict['bLs']          = cl.Buffer(self.ctx, mf.READ_ONLY, size=nSystems * nnode * 4 * np.float32().itemsize)
        self.buffer_dict['bKs']          = cl.Buffer(self.ctx, mf.READ_ONLY, size=nSystems * nnode * 4 * np.float32().itemsize)
        self.buffer_dict['Ksp']          = cl.Buffer(self.ctx, mf.READ_ONLY, size=nSystems * nnode * 4 * np.float32().itemsize)
        self.buffer_dict['Kpp']          = cl.Buffer(self.ctx, mf.READ_ONLY, size=nSystems * nnode * 4 * np.float32().itemsize)
        #self.buffer_dict['tors2atom']    = cl.Buffer(self.ctx, mf.READ_ONLY, size=nSystems * ntors * np.int32().itemsize)
        #self.buffer_dict['torsParams']   = cl.Buffer(self.ctx, mf.READ_ONLY, size=nSystems * ntors * 4 * np.float32().itemsize)
        self.buffer_dict['constr']       = cl.Buffer(self.ctx, mf.READ_WRITE, size=nSystems * natoms * 4 * np.float32().itemsize)
        self.buffer_dict['constrK']      = cl.Buffer(self.ctx, mf.READ_WRITE, size=nSystems * natoms * 4 * np.float32().itemsize)
        self.buffer_dict['lvecs']        = cl.Buffer(self.ctx, mf.READ_ONLY, size=nSystems * 9 * np.float32().itemsize)  # 3x3 m
        self.buffer_dict['ilvecs']       = cl.Buffer(self.ctx, mf.READ_ONLY, size=nSystems * 9 * np.float32().itemsize)  # 3x3 matrix
        self.buffer_dict['MDparams']   = cl.Buffer(self.ctx, mf.READ_ONLY, size=nSystems * 4 * np.float32().itemsize)
        self.buffer_dict['TDrives']      = cl.Buffer(self.ctx, mf.READ_ONLY, size=nSystems * 4 * np.float32().itemsize)
        self.buffer_dict['bboxes']       = cl.Buffer(self.ctx, mf.READ_ONLY, size=nSystems * 9 * np.float32().itemsize)
        self.buffer_dict['sysneighs']    = cl.Buffer(self.ctx, mf.READ_ONLY, size=nSystems * np.int32().itemsize)
        self.buffer_dict['sysbonds']     = cl.Buffer(self.ctx, mf.READ_ONLY, size=nSystems * 4 * np.float32().itemsize)
        self.buffer_dict['pbc_shifts']   = cl.Buffer(self.ctx, mf.READ_ONLY, size=nSystems * npbc *4* np.float32().itemsize)
        # Additional buffers as per your C++ code
        # ...

    def pack_system(self, iSys, system, params, bParams=False, bForces=False, bVel=False, blvec=True, l_rnd=-1):
        """
        Packs data from an AtomicSystem into the OpenCL buffers for the specified system index.
        """
        mmff = self.mmff_instances[iSys]
        mmff.assign_force_field_parameters(system, params, bParams=bParams, bEPairs=True, bUFF=False)
        mmff.make_back_neighs()

        # Example: Pack atomic positions
        atoms_flat = mmff.apos.flatten().astype(np.float32)
        cl.enqueue_copy(self.queue, self.buffer_dict['apos'], atoms_flat, byte_offset=iSys * mmff.nvecs * 3 * 4)

        # Pack forces if required
        if bForces:
            forces_flat = mmff.fapos.flatten().astype(np.float32)
            cl.enqueue_copy(self.queue, self.buffer_dict['aforce'], forces_flat, byte_offset=iSys * mmff.nvecs * 3 * 4)

        # Pack parameters if required
        if bParams:
            neighs_flat = mmff.neighs.flatten().astype(np.int32)
            cl.enqueue_copy(self.queue, self.buffer_dict['neighs'], neighs_flat, byte_offset=iSys * mmff.natoms * 4 * 4)

            REQs_flat = mmff.REQs.flatten().astype(np.float32)
            cl.enqueue_copy(self.queue, self.buffer_dict['REQs'], REQs_flat, byte_offset=iSys * mmff.natoms * 4 * 4)

            # Pack other parameters similarly
            # ...

        # Synchronize the queue to ensure data is copied
        self.queue.finish()

    def upload_sys(self, iSys, bParams=False, bForces=False, bVel=True, blvec=True):
        """
        Uploads system data from host to device.
        """
        # Since data is already packed into buffers using enqueue_copy with byte offsets,
        # this function can remain empty or handle additional synchronization if needed
        pass


    def toGPU(self, buf_name, host_data, byte_offset=0):
        """Uploads data to GPU buffer."""
        cl.enqueue_copy(self.queue, self.buffer_dict[buf_name], host_data, device_offset=byte_offset)

    def bufflist(self, names):
        """Returns a list of buffers based on the provided names."""
        return [self.buffer_dict[name] for name in names]

    def pack_system(self, iSys, mmff):
        """Packs data from an MMFF instance into GPU buffers for a specific system index."""
        nvecs = mmff.nvecs
        natoms = mmff.natoms
        nnode = mmff.nnode
        float4_size = 4 * np.float32().itemsize
        int4_size = 4 * np.int32().itemsize

        offset_atoms = iSys * nvecs * float4_size
        self.toGPU('apos', mmff.apos.astype(np.float32).flatten(), byte_offset=offset_atoms)
        self.toGPU('aforce', mmff.fapos.astype(np.float32).flatten(), byte_offset=offset_atoms)

        offset_REQs = iSys * natoms * float4_size
        self.toGPU('REQs', mmff.REQs.astype(np.float32).flatten(), byte_offset=offset_REQs)

        offset_neighs = iSys * natoms * int4_size
        self.toGPU('neighs', mmff.neighs.astype(np.int32).flatten(), byte_offset=offset_neighs)
        self.toGPU('neighCell', mmff.neighCell.astype(np.int32).flatten(), byte_offset=offset_neighs)

        offset_apars = iSys * nnode * float4_size
        self.toGPU('apars', mmff.apars.astype(np.float32).flatten(), byte_offset=offset_apars)
        self.toGPU('bLs', mmff.bLs.astype(np.float32).flatten(), byte_offset=offset_apars)
        self.toGPU('bKs', mmff.bKs.astype(np.float32).flatten(), byte_offset=offset_apars)
        self.toGPU('Ksp', mmff.Ksp.astype(np.float32).flatten(), byte_offset=offset_apars)
        self.toGPU('Kpp', mmff.Kpp.astype(np.float32).flatten(), byte_offset=offset_apars)
        # Continue for other buffers as needed...

    def upload_all_systems(self):
        """Uploads data for all systems to the GPU."""
        for sys_idx in range(self.nSystems):
            self.pack_system(sys_idx, self.mmff_instances[sys_idx])
        #if self.verbose:
        print("All systems uploaded to GPU.")

    def parse_kernel_header(self, header_string):
        """
        Parse a kernel header to extract buffer and parameter information.
        Simple version that assumes arguments are in correct order and one per line.
        
        Args:
            header_string (str): OpenCL kernel function header as a string
            
        Returns:
            tuple: (buffer_names, param_names, param_types)
        """
        # Extract parameter lines (everything between parentheses)
        param_block = header_string[header_string.find('(') + 1:]
        end_idx = param_block.rfind(')')
        if end_idx != -1:
            param_block = param_block[:end_idx]
        
        # Split the parameter block into lines and clean them up
        param_lines = [line.strip() for line in param_block.split('\n')]
        param_lines = [line for line in param_lines if line and not line.startswith('//') and not line.startswith('/*')]
        
        # Extract buffer names and parameter info
        args = [ ]
        
        for line in param_lines:
            # Skip lines that are just comments or empty
            if not line or (line.startswith('//') and '__global' not in line and 'const' not in line):
                continue
                
            # Remove inline comments if present
            if '//' in line:
                line = line.split('//')[0].strip()
            
            # Skip if the line is now empty
            if not line:
                continue
                
            # Strip trailing comma if present
            if line.endswith(','):
                line = line[:-1].strip()
            
            # Check if it's a buffer parameter (__global)
            if '__global' in line:
                # Extract parameter name (removing any pointers)
                parts = line.split()
                param_name = parts[-1].strip()
                # Remove pointer symbol if present
                param_name = param_name.replace('*', '').strip()
                args.append(  (param_name,0)  )
            # Check if it's a direct parameter (const)
            elif 'const' in line:
                # Extract type and name
                parts = line.split()
                if len(parts) >= 3:  # Should have at least 'const', type, and name
                    param_name = parts[2]  # Name is after the type
                    # Clean up the name in case it has trailing characters
                    param_name = param_name.strip(',;')
                    args.append(  (param_name,1)  )
        
        return args

    def init_kernel_params(self):
        """
        Initialize a dictionary of standard kernel parameters.
        This provides default values for common parameters used in kernels.
        """
        # Create a dictionary to store kernel parameters
        self.kernel_params = {
            # Common dimension parameters
            'ns': np.array([self.natoms, self.nnode, 0, 0], dtype=np.int32),
            'n': np.array([self.natoms, self.nnode, 0, 0], dtype=np.int32),
            'mask': np.array([1, 1, 1, 1], dtype=np.int32),
            'nPBC': np.array([1, 1, 1, 0], dtype=np.int32),
            'GFFParams': np.array([0.0, 0.0, 0.0, 0.0], dtype=np.float32),
            # Common scalar parameters
            'npbc': np.int32(1),
            'bSubtractVdW': np.int32(0),
            'nDOFs': np.int32(self.nDOFs),
        }
        
    def get_work_sizes(self):
        """
        Generate standard work sizes based on current system dimensions.
        
        Returns:
            dict: Dictionary containing sz_na, sz_node, sz_nvec, sz_loc
        """
        # Default to the standard work size parameters
        sz_loc = (self.nloc, 1)
        sz_na = (clu.roundup_global_size(self.natoms, self.nloc), self.nSystems)
        sz_node = (clu.roundup_global_size(self.nnode, self.nloc), self.nSystems)
        sz_nvec = (clu.roundup_global_size(self.nvecs, self.nloc), self.nSystems)
        
        # Return all sizes, let the caller decide which to use
        return {
            'sz_na': sz_na,
            'sz_node': sz_node,
            'sz_nvec': sz_nvec,
            'sz_loc': sz_loc
        }
    
    def generate_kernel_args(self, kernel_header):
        """
        Generate argument list for a kernel based on its header definition.
        
        Args:
            kernel_header (str): The kernel header string
            
        Returns:
            list: List of arguments for the kernel call
        """
        # Initialize the kernel parameters if they don't exist yet
        if not hasattr(self, 'kernel_params'):
            self.init_kernel_params()
        args_names = self.parse_kernel_header(kernel_header)
        print("\n=======================\nkernel_header ", kernel_header ,  "\nargs_names: ", args_names)
        args = []
        for name,typ in args_names:
            if typ == 0:
                args.append(self.buffer_dict[name])
            else:
                args.append(self.kernel_params[name])
        return args

    def setup_kernels(self):
        """
        Prepares the kernel arguments for all kernels by parsing their headers.
        Also sets up the work sizes for each kernel.
        """
        # Get all work sizes at once
        work_sizes = self.get_work_sizes()
        sz_loc     = work_sizes['sz_loc']
        sz_na      = work_sizes['sz_na']
        sz_node    = work_sizes['sz_node']
        sz_nvec    = work_sizes['sz_nvec']
        
        # Initialize kernel parameters
        self.init_kernel_params()
        
        # Set global and local work sizes for each kernel
        self.global_size_mmff = (self.nSystems * self.nnode,)
        self.global_size_nonbond = (self.nSystems * self.natoms,)
        self.global_size_update = (self.nSystems * self.nvecs,)
        self.global_size_clean = (self.nSystems * self.natoms,)
        
        # Generate kernel arguments
        # For getMMFFf4 kernel
        args_mmff = [sz_node, sz_loc]  # Start with work sizes
        args_mmff.extend(self.generate_kernel_args(self.kernelheader_getMMFFf4))
        self.kernel_args_getMMFFf4 = args_mmff
        
        # For getNonBond kernel
        args_nonbond = [sz_na, sz_loc]  # Start with work sizes
        args_nonbond.extend(self.generate_kernel_args(self.kernelheader_getNonBond))
        self.kernel_args_getNonBond = args_nonbond
        
        # For updateAtomsMMFFf4 kernel
        args_update = [sz_na, sz_loc]  # Start with work sizes
        args_update.extend(self.generate_kernel_args(self.kernelheader_updateAtomsMMFFf4))
        self.kernel_args_updateAtomsMMFFf4 = args_update
        
        # For cleanForceMMFFf4 kernel
        args_clean = [sz_nvec, sz_loc]  # Start with work sizes
        args_clean.extend(self.generate_kernel_args(self.kernelheader_cleanForceMMFFf4))
        self.kernel_args_cleanForceMMFFf4 = args_clean
        
        print("Kernel arguments prepared for optimization.")


    def run_ocl_opt(self, niter, Fconv=1e-6, nPerVFs=10):
        F2conv = Fconv ** 2
        niterdone = 0
        #mmff = self.mmff_instances[0]  # Assuming all instances are identical
        for _ in range(niter):
            self.getNonBond       (*self.kernel_args_getNonBond)
            self.getMMFFf4        (*self.kernel_args_getMMFFf4)
            self.updateAtomsMMFFf4(*self.kernel_args_updateAtomsMMFFf4)
            niterdone += 1
            if niterdone % nPerVFs == 0:
                # Read back forces to check convergence
                #aforces = np.empty(self.nSystems * mmff.nvecs * 4, dtype=np.float32)
                cl.enqueue_copy(self.queue, self.aforces, self.buffer_dict['aforce'])
                F2max = np.max(np.sum( self.aforces**2, axis=1))
                #if self.verbose: 
                print(f"Iteration {niterdone}: Max |F|^2 = {F2max}")
                if F2max < F2conv:
                    #if self.verbose: 
                    print(f"Converged after {niterdone} iterations.")
                    break
        if self.verbose and F2max >= F2conv: 
            print(f"Did not converge after {niterdone} iterations. Final |F|^2 = {F2max}")
        return niterdone

    def download_results(self):
        mmff   = self.mmff_instances[0]
        cl.enqueue_copy(self.queue, self.atoms,  self.buffer_dict['apos'])
        cl.enqueue_copy(self.queue, self.forces, self.buffer_dict['aforce'])
        self.queue.finish()
        #if self.verbose: print("Downloaded results from GPU.")
        return self.atoms.reshape(self.nSystems, mmff.nvecs, 4), self.forces.reshape(self.nSystems, mmff.nvecs, 4)

    def clean_forces(self):
        self.cleanForceMMFFf4(*self.kernel_args_cleanForceMMFFf4)
        #if self.verbose:  print("Forces cleaned using cleanForceMMFFf4 kernel.")