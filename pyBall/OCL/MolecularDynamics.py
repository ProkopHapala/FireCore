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
        base_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        possible_paths = [
            '../../cpp/common_resources/cl/relax_multi.cl',  # Relative from working dir
            os.path.join(base_path, 'cpp/common_resources/cl/relax_multi.cl'),  # From project root
            os.path.join(os.path.dirname(__file__), '../../cpp/common_resources/cl/relax_multi.cl')  # Relative from module location
        ]
        
        kernel_found = False
        for kernel_path in possible_paths:
            try:
                if os.path.exists(kernel_path):
                    with open(kernel_path, 'r') as f:
                        self.prg = cl.Program(self.ctx, f.read()).build()
                    kernel_found = True
                    print(f"Successfully loaded kernel from: {kernel_path}")
                    break
            except Exception as e:
                # Continue to the next path option
                pass
                
        if not kernel_found:
            print(f"Error: Could not find or compile OpenCL program")
            print(f"MolecularDynamics() called from path: {os.getcwd()}")
            print(f"Tried the following paths:")
            for path in possible_paths:
                print(f"  - {path} (exists: {os.path.exists(path)})")
            exit(0)
        
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
        self.nDOFs = (natoms,nnode)

        self.natoms = natoms
        self.nvecs  = nvecs
        self.nnode  = nnode
        self.ncap   = ncap
        self.ntors  = ntors
        self.nbkng  = nbkng
        
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
        self.buffer_dict['atoms']        = cl.Buffer(self.ctx, mf.READ_WRITE, size=nSystems * nvecs * 3 * np.float32().itemsize)
        self.buffer_dict['aforces']      = cl.Buffer(self.ctx, mf.READ_WRITE, size=nSystems * nvecs * 3 * np.float32().itemsize)
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
        self.buffer_dict['invLvec']      = cl.Buffer(self.ctx, mf.READ_ONLY, size=nSystems * 9 * np.float32().itemsize)  # 3x3 matrix
        self.buffer_dict['MDpars']       = cl.Buffer(self.ctx, mf.READ_ONLY, size=nSystems * 4 * np.float32().itemsize)
        self.buffer_dict['TDrives']      = cl.Buffer(self.ctx, mf.READ_ONLY, size=nSystems * 4 * np.float32().itemsize)
        self.buffer_dict['bboxes']       = cl.Buffer(self.ctx, mf.READ_ONLY, size=nSystems * 9 * np.float32().itemsize)
        self.buffer_dict['sysneighs']    = cl.Buffer(self.ctx, mf.READ_ONLY, size=nSystems * np.int32().itemsize)
        self.buffer_dict['sysbonds']     = cl.Buffer(self.ctx, mf.READ_ONLY, size=nSystems * 4 * np.float32().itemsize)

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
        cl.enqueue_copy(self.queue, self.buffer_dict['atoms'], atoms_flat, byte_offset=iSys * mmff.nvecs * 3 * 4)

        # Pack forces if required
        if bForces:
            forces_flat = mmff.fapos.flatten().astype(np.float32)
            cl.enqueue_copy(self.queue, self.buffer_dict['aforces'], forces_flat, byte_offset=iSys * mmff.nvecs * 3 * 4)

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

    def setup_kernels(self):
        """
        Sets up OpenCL kernels with their respective arguments.
        """
        # Example: Setting up getNonBond kernel
        self.getNonBond_kernel = self.prg.getNonBond
        self.getNonBond_kernel.set_arg(0, self.buffer_dict['atoms'])
        self.getNonBond_kernel.set_arg(1, self.buffer_dict['forces'])
        self.getNonBond_kernel.set_arg(2, self.buffer_dict['REQs'])
        self.getNonBond_kernel.set_arg(3, self.buffer_dict['neighs'])
        self.getNonBond_kernel.set_arg(4, self.buffer_dict['neighCell'])
        self.getNonBond_kernel.set_arg(5, self.buffer_dict['lvecs'])
        self.getNonBond_kernel.set_arg(6, self.buffer_dict['ilvecs'])
        self.getNonBond_kernel.set_arg(7, self.buffer_dict['pbc_shifts'])

        # Similarly set up other kernels like getMMFFf4, updateAtomsMMFFf4, etc.
        # ...

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
        self.toGPU('atoms', mmff.apos.astype(np.float32).flatten(), byte_offset=offset_atoms)
        self.toGPU('aforces', mmff.fapos.astype(np.float32).flatten(), byte_offset=offset_atoms)

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

    def setup_kernels(self):
        self.getMMFFf4         = self.prg.getMMFFf4
        self.getNonBond        = self.prg.getNonBond
        self.updateAtomsMMFFf4 = self.prg.updateAtomsMMFFf4
        self.cleanForceMMFFf4  = self.prg.cleanForceMMFFf4
        #if self.verbose:
        print("Kernels set up.")

    def setup_run_ocl_opt(self):
        sz_loc    = (self.nloc, 1 )
        sz_na   = ( clu.roundup_global_size( self.natoms, self.nloc), self.nSystems )
        sz_node = ( clu.roundup_global_size( self.nnode,  self.nloc), self.nSystems )
        sz_nvec = ( clu.roundup_global_size( self.nvecs,  self.nloc), self.nSystems )

        """
        // ======================================================================
        //                          getMMFFf4()
        // ======================================================================
        // 1.  getMMFFf4() - computes bonding interactions between atoms and nodes and its neighbors (max. 4 neighbors allowed), the resulting forces on atoms are stored "fapos" array and recoil forces on neighbors are stored in "fneigh" array
        //                   kernel run over all atoms and all systems in parallel to exploit GPU parallelism
        //__attribute__((reqd_work_group_size(1,1,1)))
        __kernel void getMMFFf4(
            const int4 nDOFs,               // 1   (nAtoms,nnode) dimensions of the system
            // Dynamical
            __global float4*  apos,         // 2  [natoms]     positions of atoms (including node atoms [0:nnode] and capping atoms [nnode:natoms] and pi-orbitals [natoms:natoms+nnode] )
            __global float4*  fapos,        // 3  [natoms]     forces on    atoms (just node atoms are evaluated)
            __global float4*  fneigh,       // 4  [nnode*4*2]  recoil forces on neighbors (and pi-orbitals)
            // parameters
            __global int4*    neighs,       // 5  [nnode]  neighboring atoms
            __global int4*    neighCell,    // 5  [nnode]  neighboring atom  cell index
            __global float4*  REQKs,        // 6  [natoms] non-boding parametes {R0,E0,Q} i.e. R0: van der Waals radii, E0: well depth and partial charge, Q: partial charge 
            __global float4*  apars,        // 7  [nnode]  per atom forcefield parametrs {c0ss,Kss,c0sp}, i.e. c0ss: cos(equlibrium angle/2) for sigma-sigma; Kss: stiffness of sigma-sigma angle; c0sp: is cos(equlibrium angle) for sigma-pi
            __global float4*  bLs,          // 8  [nnode]  bond length    between node and each neighbor
            __global float4*  bKs,          // 9  [nnode]  bond stiffness between node and each neighbor
            __global float4*  Ksp,          // 10 [nnode]  stiffness of pi-alignment for each neighbor     (only node atoms have pi-pi alignemnt interaction)
            __global float4*  Kpp,          // 11 [nnode]  stiffness of pi-planarization for each neighbor (only node atoms have pi-pi alignemnt interaction)
            __global cl_Mat3* lvecs,        // 12 lattice vectors         for each system
            __global cl_Mat3* ilvecs,       // 13 inverse lattice vectors for each system
            __global float4*  pbc_shifts,   // 14 shifts of pbc images
            const int npbc,
            const int bSubtractVdW
        ){
        """
        global_size_mmff = (self.nSystems * self.nnode,)
        local_size_mmff = None  # Let OpenCL decide
        self.kernel_args_getMMFFf4 = [
            sz_node,
            (1,1), #sz_loc,
            np.int32(self.nDOFs),
            *self.bufflist([
                'atoms', 'aforces', 'fneigh', 'neighs', 'neighCell', 'REQs', 'apars',
                'bLs', 'bKs', 'Ksp', 'Kpp',
                #'angles', 'tors2atom','torsParams', 'constr', 'constrK', 'invLvec', 'pbc_shifts'
            ]),
            np.int32(1),  # npbc
            np.int32(0)   # bSubtractVdW
        ]

        """
        // ======================================================================
        //                           getNonBond()
        // ======================================================================
        // Calculate non-bonded forces on atoms (icluding both node atoms and capping atoms), cosidering periodic boundary conditions
        // It calculate Lenard-Jones, Coulomb and Hydrogen-bond forces between all atoms in the system
        // it can be run in parallel for multiple systems, in order to efficiently use number of GPU cores (for small systems with <100 this is essential to get good performance)
        // This is the most time consuming part of the forcefield evaluation, especially for large systems when nPBC>1
        __attribute__((reqd_work_group_size(32,1,1)))
        __kernel void getNonBond(
            const int4 ns,                  // 1 // (natoms,nnode) dimensions of the system
            // Dynamical
            __global float4*  atoms,        // 2 // positions of atoms  (including node atoms [0:nnode] and capping atoms [nnode:natoms] and pi-orbitals [natoms:natoms+nnode] )
            __global float4*  forces,       // 3 // forces on atoms
            // Parameters
            __global float4*  REQKs,        // 4 // non-bonded parameters (RvdW,EvdW,QvdW,Hbond)
            __global int4*    neighs,       // 5 // neighbors indices      ( to ignore interactions between bonded atoms )
            __global int4*    neighCell,    // 6 // neighbors cell indices ( to know which PBC image should be ignored  due to bond )
            __global cl_Mat3* lvecs,        // 7 // lattice vectors for each system
            const int4        nPBC,         // 8 // number of PBC images in each direction (x,y,z)
            const float4      GFFParams     // 9 // Grid-Force-Field parameters
        );
        """
        global_size_nonbond = (self.nSystems * self.natoms,)
        local_size_nonbond = None
        self.kernel_args_getNonBond = [
            sz_na,
            sz_loc,
            np.array([self.natoms, self.nnode, 0, 0], dtype=np.int32),
            *self.bufflist(['atoms', 'aforces', 'REQs', 'neighs', 'neighCell', 'invLvec'])
        ]

        """
        // Assemble recoil forces from neighbors and  update atoms positions and velocities 
        //__attribute__((reqd_work_group_size(1,1,1)))
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
        )
        """
        global_size_update = (self.nSystems * self.nvecs,)
        local_size_update  = self.nloc
        self.kernel_args_updateAtomsMMFFf4 = [
            sz_na,
            sz_loc,
            *self.bufflist([
                'atoms', 'avel', 'aforces', 'cvf', 'fneigh', 'bkNeighs',
                'constr', 'constrK', 'MDpars', 'TDrives', 'bboxes',
                'sysneighs', 'sysbonds'
            ])
        ]

        """
        __kernel void cleanForceMMFFf4(
            const int4        n,           // 2
            __global float4*  aforce,      // 5
            __global float4*  fneigh       // 6
        ){
        """
        global_size_clean = (self.nSystems * self.natoms,)
        local_size_clean = None
        self.kernel_args_cleanForceMMFFf4 = [
            sz_nvec,
            sz_loc,
            *self.bufflist(['aforces', 'fneigh'])
        ]

        #if self.verbose:
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
                cl.enqueue_copy(self.queue, self.aforces, self.buffer_dict['aforces'])
                F2max = np.max(np.sum( self.aforces**2, axis=1))
                #if self.verbose: 
                print(f"Iteration {niterdone}: Max |F|^2 = {F2max}")
                if F2max < F2conv:
                    #if self.verbose: 
                    print(f"Converged after {niterdone} iterations.")
                    break
        if self.verbose and F2max >= F2conv: print(f"Did not converge after {niterdone} iterations. Final |F|^2 = {F2max}")
        return niterdone

    def download_results(self):
        mmff   = self.mmff_instances[0]
        cl.enqueue_copy(self.queue, self.atoms,  self.buffer_dict['atoms'])
        cl.enqueue_copy(self.queue, self.forces, self.buffer_dict['aforces'])
        self.queue.finish()
        #if self.verbose: print("Downloaded results from GPU.")
        return self.atoms.reshape(self.nSystems, mmff.nvecs, 4), self.forces.reshape(self.nSystems, mmff.nvecs, 4)

    def clean_forces(self):
        self.cleanForceMMFFf4(*self.kernel_args_cleanForceMMFFf4)
        #if self.verbose:  print("Forces cleaned using cleanForceMMFFf4 kernel.")