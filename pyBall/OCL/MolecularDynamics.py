import sys
import os
import numpy as np
import ctypes
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

    def __init__(self, nloc=32 ):
        self.nloc  = nloc
        self.ctx   = cl.create_some_context()
        self.queue = cl.CommandQueue(self.ctx)

        self.grid  = None   # instance of GridShape, if initialized
        self.gcl   = None   # instance of GridCL, if initialized
 
        clu.get_cl_info( self.ctx.devices[0] )

        local_size = 64
        #print( " local_memory_per_workgroup() size=", local_size, " __local []  ", clu.local_memory_per_workgroup( self.ctx.devices[0], local_size=32, sp_per_cu=128 ), " Byte " );
        print( " local_memory_per_workgroup() size=", local_size, " __local []  ", clu.local_memory_float_per_workgroup( self.ctx.devices[0], local_size=32, sp_per_cu=128 ), " float32 " );

        try:
            with open('../../cpp/common_resources/cl/splines.cl', 'r') as f:
                self.prg = cl.Program(self.ctx, f.read()).build()
        except Exception as e:
            print( "MolecularDynamics() called from path=", os.getcwd() )
            print(f"Error compiling OpenCL program: {e}")
            exit(0)

    def __init__(self, nloc=32 ):
        self.ctx   = cl.create_some_context()
        self.queue = cl.CommandQueue(self.ctx)
        clu.get_cl_info( self.ctx.devices[0] )

        # Initialize OpenCL context and queue
        # platforms = cl.get_platforms()
        # if not platforms:                 raise RuntimeError("No OpenCL platforms found.")
        # if platform_id >= len(platforms): raise ValueError("Invalid platform ID.")
        # platform = platforms[platform_id]
        # devices = platform.get_devices()
        # if not devices:                raise RuntimeError("No OpenCL devices found on the selected platform.")
        # if device_id >= len(devices):  raise ValueError("Invalid device ID.")
        # device = devices[device_id]
        # self.context = cl.Context([device])
        # self.queue = cl.CommandQueue(self.context, device)

        try:
            with open('../../cpp/common_resources/cl/relax_multi.cl', 'r') as f:
                self.prg = cl.Program(self.ctx, f.read()).build()
        except Exception as e:
            print(f"Error compiling OpenCL program: {e}")
            print( "MolecularDynamics() called from path=", os.getcwd() )
            exit(0)

        self.nloc  = nloc
        self.grid  = None   # instance of GridShape, if initialized
        self.gcl   = None   # instance of GridCL, if initialized

        # Initialize OpenCL program
        self.program = self.build_program()

        # Initialize other attributes
        self.nSystems = 0
        self.mmff_instances = []
        self.buffer_dict = {}


    def build_program(self):
        """
        Compiles the OpenCL kernels.
        """
        # Load OpenCL kernel code from file or string
        kernel_file = "kernels.cl"  # Replace with your actual kernel file path
        if not os.path.exists(kernel_file):
            raise FileNotFoundError(f"OpenCL kernel file '{kernel_file}' not found.")

        with open(kernel_file, 'r') as f:
            kernel_source = f.read()

        try:
            program = cl.Program(self.context, kernel_source).build()
            if verbose: print("OpenCL program built successfully.")
            return program
        except Exception as e:
            if verbose:
                print("Error building OpenCL program:")
                print(e)
            raise e

    def realloc(self, nSystems, mmff):
        """
        Reallocate buffers for the given number of systems based on the MMFF template.
        """
        self.nSystems = nSystems
        self.mmff_instances = [ MMFF() for _ in range(nSystems)]

        # Example: Assuming mmff_template is an instance of MMFF with predefined parameters
        for mmff in self.mmff_instances:
            mmff.realloc(mmff.nnode, mmff.ncap, mmff.ntors)

        # Allocate OpenCL buffers
        self.allocate_cl_buffers(mmff)
        self.allocate_host_buffers(mmff)

    def allocate_host_buffers(self, mmff ):
        self.atoms  = np.empty(self.nSystems * mmff.nvecs * 4, dtype=np.float32)
        self.forces = np.empty(self.nSystems * mmff.nvecs * 4, dtype=np.float32)

    def allocate_cl_buffers(self, mmff ):
        """
        Allocates OpenCL buffers based on the MMFF template and number of systems.
        """
        # Determine buffer sizes
        nvecs  = mmff.nvecs
        natoms = mmff.natoms
        nnode  = mmff.nnode
        ntors  = mmff.ntors
        nbkng  = nnode * 4 * 2  # As per C++ code
        npbc = 1  # Adjust based on your requirements
        nSystems = self.nSystems

        # Example buffer allocations
        mf = cl.mem_flags
        self.buffer_dict['atoms']        = cl.Buffer(self.context, mf.READ_WRITE, size=nSystems * nvecs * 3 * np.float32().itemsize)
        self.buffer_dict['aforces']      = cl.Buffer(self.context, mf.READ_WRITE, size=nSystems * nvecs * 3 * np.float32().itemsize)
        self.buffer_dict['REQs']         = cl.Buffer(self.context, mf.READ_ONLY, size=nSystems * natoms * 4 * np.float32().itemsize)
        self.buffer_dict['neighs']       = cl.Buffer(self.context, mf.READ_ONLY, size=nSystems * natoms * 4 * np.int32().itemsize)
        self.buffer_dict['neighCell']    = cl.Buffer(self.context, mf.READ_ONLY, size=nSystems * natoms * 4 * np.int32().itemsize)
        self.buffer_dict['bkNeighs']     = cl.Buffer(self.context, mf.READ_ONLY, size=nSystems * nvecs * 4 * np.int32().itemsize)
        self.buffer_dict['bkNeighs_new'] = cl.Buffer(self.context, mf.READ_ONLY, size=nSystems * nvecs * 4 * np.int32().itemsize)
        self.buffer_dict['avel']         = cl.Buffer(self.context, mf.READ_WRITE, size=nSystems * nvecs * 3 * np.float32().itemsize)
        self.buffer_dict['cvf']          = cl.Buffer(self.context, mf.READ_WRITE, size=nSystems * nvecs * 3 * np.float32().itemsize)
        self.buffer_dict['fneigh']       = cl.Buffer(self.context, mf.READ_WRITE, size=nSystems * nbkng * 3 * np.float32().itemsize)
        self.buffer_dict['fneighpi']     = cl.Buffer(self.context, mf.READ_WRITE, size=nSystems * nbkng * 3 * np.float32().itemsize)
        self.buffer_dict['apars']        = cl.Buffer(self.context, mf.READ_ONLY, size=nSystems * nnode * 4 * np.float32().itemsize)
        self.buffer_dict['bLs']          = cl.Buffer(self.context, mf.READ_ONLY, size=nSystems * nnode * 4 * np.float32().itemsize)
        self.buffer_dict['bKs']          = cl.Buffer(self.context, mf.READ_ONLY, size=nSystems * nnode * 4 * np.float32().itemsize)
        self.buffer_dict['Ksp']          = cl.Buffer(self.context, mf.READ_ONLY, size=nSystems * nnode * 4 * np.float32().itemsize)
        self.buffer_dict['Kpp']          = cl.Buffer(self.context, mf.READ_ONLY, size=nSystems * nnode * 4 * np.float32().itemsize)
        #self.buffer_dict['tors2atom']    = cl.Buffer(self.context, mf.READ_ONLY, size=nSystems * ntors * np.int32().itemsize)
        #self.buffer_dict['torsParams']   = cl.Buffer(self.context, mf.READ_ONLY, size=nSystems * ntors * 4 * np.float32().itemsize)
        #self.buffer_dict['constr']       = cl.Buffer(self.context, mf.READ_WRITE, size=nSystems * natoms * 4 * np.float32().itemsize)
        #self.buffer_dict['constrK']      = cl.Buffer(self.context, mf.READ_WRITE, size=nSystems * natoms * 4 * np.float32().itemsize)
        self.buffer_dict['invLvec']      = cl.Buffer(self.context, mf.READ_ONLY, size=nSystems * 9 * np.float32().itemsize)  # 3x3 matrix

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
        self.getNonBond_kernel = self.program.getNonBond
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
        if self.verbose:
            print("All systems uploaded to GPU.")

    def setup_kernels(self):
        self.getMMFFf4         = self.program.getMMFFf4
        self.getNonBond        = self.program.getNonBond
        self.updateAtomsMMFFf4 = self.program.updateAtomsMMFFf4
        self.cleanForceMMFFf4  = self.program.cleanForceMMFFf4
        if self.verbose:
            print("Kernels set up.")

    def setup_run_ocl_opt(self):
        mmff = self.mmff_instances[0]  # Assuming all instances are identical
        sz_loc    = (self.nloc, 1 )
        sz_na   = ( clu.roundup_global_size( self.mmff.natoms, self.nloc), self.nSystems )
        sz_node = ( clu.roundup_global_size( self.mmff.nnode,  self.nloc), self.nSystems )
        sz_nvec = ( clu.roundup_global_size( self.mmff.nvecs,  self.nloc), self.nSystems )

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
        global_size_mmff = (self.nSystems * mmff.nnode,)
        local_size_mmff = None  # Let OpenCL decide
        self.kernel_args_getMMFFf4 = [
            sz_node,
            (1,1), #sz_loc,
            np.int32(mmff.nDOFs),
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
        global_size_nonbond = (self.nSystems * mmff.natoms,)
        local_size_nonbond = None
        self.kernel_args_getNonBond = [
            sz_na,
            sz_loc,
            np.array([mmff.natoms, mmff.nnode, 0, 0], dtype=np.int32),
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
        global_size_update = (self.nSystems * mmff.nvecs,)
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
        global_size_clean = (self.nSystems * mmff.natoms,)
        local_size_clean = None
        self.kernel_args_cleanForceMMFFf4 = [
            sz_nvec,
            sz_loc,
            *self.bufflist(['aforces', 'fneigh'])
        ]

        if self.verbose:
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
                if self.verbose: print(f"Iteration {niterdone}: Max |F|^2 = {F2max}")
                if F2max < F2conv:
                    if self.verbose: print(f"Converged after {niterdone} iterations.")
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