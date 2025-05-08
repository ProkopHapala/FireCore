import sys
import os
import numpy as np
import re

import pyopencl as cl
import pyopencl.array as cl_array
import pyopencl.cltypes as cltypes
import matplotlib.pyplot as plt
import time

from . import clUtils as clu
from .MMFF import MMFF
from .OpenCLBase import OpenCLBase

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

class MolecularDynamics(OpenCLBase):
    """
    Class for molecular dynamics simulations using OpenCL.
    
    This class inherits from OpenCLBase and implements specific functionality
    for molecular dynamics simulations using the relax_multi_mini.cl kernel.
    """
    
    def __init__(self, nloc=32, perBatch=10):
        # Initialize the base class
        super().__init__(nloc=nloc, device_index=0)
        
        # Load the OpenCL program
        base_path = os.path.dirname(os.path.abspath(__file__))
        rel_path = "../../cpp/common_resources/cl/relax_multi_mini.cl"
        if not self.load_program(rel_path=rel_path, base_path=base_path, bPrint=False):
            exit(1)
        
        # Initialize other attributes that will be set in realloc
        self.nSystems = 0
        self.mmff_instances = []
        self.MD_event_batch = None
        self.perBatch = perBatch

    def realloc(self, mmff, nSystems=1 ):
        """
        Reallocate buffers for the given number of systems based on the MMFF template.
        """
        # Store dimensions explicitly to avoid reference issues
        print(f"MolecularDynamics::realloc() natoms={mmff.natoms}, nvecs={mmff.nvecs}, nnode={mmff.nnode}")
        self.nSystems = nSystems
        self.mmff_instances = [mmff] * nSystems  # Assuming all systems use the same MMFF parameters
        self.allocate_cl_buffers(mmff)
        self.allocate_host_buffers()

    def allocate_host_buffers(self):
        self.atoms  = np.zeros((self.nSystems, self.nvecs, 4), dtype=np.float32)
        self.aforce = np.zeros((self.nSystems, self.nvecs, 4), dtype=np.float32)

    def allocate_cl_buffers(self, mmff):
        """
        Allocates OpenCL buffers based on the MMFF template and number of systems.
        Includes all buffers required by the runMD kernel.
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
        
        print(f"MolecularDynamics::allocate_cl_buffers(): nSystems: {nSystems}  natoms: {natoms}  nvecs: {nvecs} nnode: {nnode} ncap: {ncap}  ntors: {ntors}  nbkng: {nbkng}")
        
        if nSystems <= 0 or natoms <= 0 or nvecs <= 0 or nnode <= 0:
            raise ValueError(f"Invalid dimensions for buffer allocation: nSystems={nSystems}, natoms={natoms}, nvecs={nvecs}, nnode={nnode}")
        
        float_size = np.float32().itemsize
        int_size = np.int32().itemsize
        mat3_size = 3 * 4 * float_size  # 3x3 matrix
        
        mf = cl.mem_flags
        
        # Dynamical variables
        self.create_buffer('apos',     nSystems * nvecs * 4 * float_size, mf.READ_WRITE)
        self.create_buffer('aforce',   nSystems * nvecs * 4 * float_size, mf.READ_WRITE)
        self.create_buffer('avel',     nSystems * nvecs * 4 * float_size, mf.READ_WRITE)
        self.create_buffer('fneigh',   nSystems * nnode * 4 * 2 * float_size, mf.READ_WRITE)
        self.create_buffer('cvf',      nSystems * nvecs * 4 * float_size, mf.READ_WRITE)
        
        # Neighbor lists
        self.create_buffer('neighs',    nSystems * natoms * 4 * int_size, mf.READ_ONLY)
        self.create_buffer('neighCell', nSystems * natoms * 4 * int_size, mf.READ_ONLY)
        self.create_buffer('bkNeighs',  nSystems * nbkng * 4 * int_size, mf.READ_ONLY)
        
        # Force field parameters
        self.create_buffer('REQs',     nSystems * natoms * 4 * float_size, mf.READ_ONLY)
        self.create_buffer('apars',    nSystems * nnode * 4 * float_size, mf.READ_ONLY)
        self.create_buffer('bLs',      nSystems * nnode * 4 * float_size, mf.READ_ONLY)
        self.create_buffer('bKs',      nSystems * nnode * 4 * float_size, mf.READ_ONLY)
        self.create_buffer('Ksp',      nSystems * nnode * 4 * float_size, mf.READ_ONLY)
        self.create_buffer('Kpp',      nSystems * nnode * 4 * float_size, mf.READ_ONLY)
        
        # System parameters
        self.create_buffer('lvecs',    nSystems * mat3_size, mf.READ_ONLY)
        self.create_buffer('ilvecs',   nSystems * mat3_size, mf.READ_ONLY)
        self.create_buffer('pbc_shifts', nSystems * npbc * 4 * float_size, mf.READ_ONLY)
        
        # MD parameters and constraints
        self.create_buffer('constr',   nSystems * natoms * 4 * float_size, mf.READ_WRITE)
        self.create_buffer('constrK',  nSystems * natoms * 4 * float_size, mf.READ_WRITE)
        self.create_buffer('MDparams', nSystems * 4 * float_size, mf.READ_ONLY)
        self.create_buffer('TDrives',  nSystems * 4 * float_size, mf.READ_ONLY)
        
        # System interactions
        self.create_buffer('bboxes',   nSystems * mat3_size, mf.READ_ONLY)
        self.create_buffer('sysneighs',nSystems * int_size, mf.READ_ONLY)
        self.create_buffer('sysbonds', nSystems * 4 * float_size, mf.READ_ONLY)
        
        # Grid force field parameters (scalar)
        # self.kernel_params['GFFParams'] = np.zeros(4, dtype=np.float32)
        # self.kernel_params['bSubtractVdW'] = np.int32(0)

    def pack_system(self, iSys, mmff):
        """Packs data from an MMFF instance into GPU buffers for a specific system index."""
        nvecs   = mmff.nvecs
        natoms  = mmff.natoms
        nnode   = mmff.nnode
        float4_size = 4 * np.float32().itemsize
        int4_size   = 4 * np.int32().itemsize

        offset_atoms = iSys * nvecs * float4_size
        self.toGPU('apos',      mmff.apos.astype(np.float32).flatten(), byte_offset=offset_atoms)
        self.toGPU('aforce',    mmff.fapos.astype(np.float32).flatten(), byte_offset=offset_atoms)

        offset_REQs = iSys * natoms * float4_size
        self.toGPU('REQs',      mmff.REQs.astype(np.float32).flatten(), byte_offset=offset_REQs)
        
        offset_neighs = iSys * natoms * int4_size
        self.toGPU('neighs',    mmff.neighs.astype(np.int32).flatten(), byte_offset=offset_neighs)
        self.toGPU('neighCell', mmff.neighCell.astype(np.int32).flatten(), byte_offset=offset_neighs)
        
        offset_apars = iSys * nnode * float4_size
        self.toGPU('apars',     mmff.apars.astype(np.float32).flatten(), byte_offset=offset_apars)
        self.toGPU('bLs',       mmff.bLs.astype(np.float32).flatten(), byte_offset=offset_apars)
        self.toGPU('bKs',       mmff.bKs.astype(np.float32).flatten(), byte_offset=offset_apars)
        self.toGPU('Ksp',       mmff.Ksp.astype(np.float32).flatten(), byte_offset=offset_apars)
        self.toGPU('Kpp',       mmff.Kpp.astype(np.float32).flatten(), byte_offset=offset_apars)

        print("pack_system() iSys=%d" % iSys, "offset_atoms=%d" % offset_atoms, "offset_REQs=%d" % offset_REQs, "offset_neighs=%d" % offset_neighs, "offset_apars=%d" % offset_apars)
        print("pack_system() iSys=%d" % iSys, "atoms.nbytes=%d" % mmff.apos.nbytes, "REQs.nbytes=%d" % mmff.REQs.nbytes, "neighs.nbytes=%d" % mmff.neighs.nbytes, "apars.nbytes=%d" % mmff.apars.nbytes)

        self.toGPU('MDparams',  np.array([mmff.dt, mmff.damp, mmff.Flimit], dtype=np.float32), byte_offset=iSys*float4_size)


    def setup_kernels(self):
        """
        Prepares the kernel arguments for all kernels by parsing their headers.
        Also sets up the work sizes for each kernel.
        """
        # Get all work sizes at once
        work_sizes = self.get_work_sizes()
        sz_loc    = work_sizes['sz_loc']
        self.local_size_opt = sz_loc
        self.init_kernel_params()
        
        # Set global work sizes for each kernel
        self.global_size_mmff    = (self.roundUpGlobalSize(self.nSystems * self.nnode  ),self.nSystems)
        self.global_size_nonbond = (self.roundUpGlobalSize(self.nSystems * self.natoms ),self.nSystems)
        self.global_size_update  = (self.roundUpGlobalSize(self.nSystems * self.nvecs  ),self.nSystems)
        self.global_size_clean   = (self.roundUpGlobalSize(self.nSystems * self.natoms ),self.nSystems)
        
        # Generate kernel arguments
        self.kernel_args_getMMFFf4         = self.generate_kernel_args("getMMFFf4")
        self.kernel_args_getNonBond        = self.generate_kernel_args("getNonBond")
        self.kernel_args_updateAtomsMMFFf4 = self.generate_kernel_args("updateAtomsMMFFf4")
        self.kernel_args_cleanForceMMFFf4  = self.generate_kernel_args("cleanForceMMFFf4")
        self.kernel_args_runMD             = self.generate_kernel_args("runMD")

        # print("MolecularDynamics::setup_kernels() kernel_args_getNonBond:"); 
        # for arg in self.kernel_args_getNonBond: print("    ", arg)
        # print("MolecularDynamics::setup_kernels() kernel_args_getMMFFf4:");  
        # for arg in self.kernel_args_getMMFFf4: print("    ", arg)
        # print("MolecularDynamics::setup_kernels() kernel_args_updateAtomsMMFFf4:")
        # for arg in self.kernel_args_updateAtomsMMFFf4: print("    ", arg)
        # print("MolecularDynamics::setup_kernels() kernel_args_cleanForceMMFFf4:")
        # for arg in self.kernel_args_cleanForceMMFFf4: print("    ", arg)
        # print("MolecularDynamics::setup_kernels() kernel_args_runMD:")
        # for arg in self.kernel_args_runMD: print("    ", arg)

    def init_kernel_params(self):
        """
        Initialize a dictionary of standard kernel parameters.
        This provides default values for common parameters used in kernels.
        """
        # Create a dictionary to store kernel parameters
        self.kernel_params = {
            # Common dimension parameters
            'nDOFs':        np.array([self.natoms, self.nnode, 0, self.perBatch], dtype=np.int32),
            'mask':         np.array([1, 1, 1, 1],         dtype=np.int32),
            'nPBC':         np.array(self.nPBC+(0,),       dtype=np.int32),
            'GFFParams':    np.array([0.0, 0.0, 0.0, 0.0], dtype=np.float32),
            # Common scalar parameters
            'npbc':         np.int32(self.npbc),
            'bSubtractVdW': np.int32(0),
        }
        
    def get_work_sizes(self):
        """
        Generate standard work sizes based on current system dimensions.
        
        Returns:
            dict: Dictionary containing sz_na, sz_node, sz_nvec, sz_loc
        """
        # Default to the standard work size parameters
        sz_loc = (self.nloc, 1)
        sz_na   = (clu.roundup_global_size(self.natoms, self.nloc), self.nSystems)
        sz_node = (clu.roundup_global_size(self.nnode, self.nloc), self.nSystems)
        sz_nvec = (clu.roundup_global_size(self.nvecs, self.nloc), self.nSystems)
        
        # Return all sizes, let the caller decide which to use
        return {
            'sz_na':   sz_na,
            'sz_node': sz_node,
            'sz_nvec': sz_nvec,
            'sz_loc':  sz_loc
        }

    def upload_all_systems(self):
        """Uploads data for all systems to the GPU."""
        for sys_idx in range(self.nSystems):
            self.pack_system(sys_idx, self.mmff_instances[sys_idx])
        print("MolecularDynamics::upload_all_systems() DONE")

    def run_getNonBond(self):
        self.prg.getNonBond(self.queue, self.global_size_nonbond, self.local_size_opt, *self.kernel_args_getNonBond)
        self.queue.finish()
    
    def run_getMMFFf4(self):
        self.prg.getMMFFf4(self.queue, self.global_size_mmff, self.local_size_opt, *self.kernel_args_getMMFFf4)
        self.queue.finish()
    
    def run_updateAtomsMMFFf4(self):
        self.prg.updateAtomsMMFFf4(self.queue, self.global_size_update, self.local_size_opt, *self.kernel_args_updateAtomsMMFFf4)
        self.queue.finish()
    
    def run_cleanForceMMFFf4(self):
        self.prg.cleanForceMMFFf4(self.queue, self.global_size_clean, self.local_size_opt, *self.kernel_args_cleanForceMMFFf4)
        self.queue.finish()
    
    def run_runMD(self):
        self.prg.runMD(self.queue, self.global_size_update, self.local_size_opt, *self.kernel_args_runMD)
        self.queue.finish()
    
    def run_MD_py(self, nsteps):
        #print( "MolecularDynamics::run_MD_py() nsteps=", nsteps)
        #print( "MolecularDynamics::run_MD_py() kernel_args_getNonBond=", self.kernel_args_getNonBond)
        #print( "MolecularDynamics::run_MD_py() kernel_args_getMMFFf4=", self.kernel_args_getMMFFf4)
        #print( "MolecularDynamics::run_MD_py() kernel_args_updateAtomsMMFFf4=", self.kernel_args_updateAtomsMMFFf4)
        for i in range(nsteps):
            self.prg.getNonBond       (self.queue, self.global_size_nonbond, self.local_size_opt, *self.kernel_args_getNonBond)
            self.prg.getMMFFf4        (self.queue, self.global_size_mmff,    self.local_size_opt, *self.kernel_args_getMMFFf4)
            self.prg.updateAtomsMMFFf4(self.queue, self.global_size_update,  self.local_size_opt, *self.kernel_args_updateAtomsMMFFf4)
        self.fromGPU('apos',   self.atoms)
        self.fromGPU('aforce', self.aforce)
        self.queue.finish()
        return self.atoms.reshape(-1, 4), self.aforce.reshape(-1, 4)

    def download_results(self):
        self.fromGPU('apos',   self.atoms)
        self.fromGPU('aforce', self.aforce)
        return self.atoms.reshape(-1, 4), self.aforce.reshape(-1, 4)
