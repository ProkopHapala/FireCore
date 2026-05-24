import sys
import os
import numpy as np
import re
import math

import pyopencl as cl
# from . import clUtils as clu
import pyopencl.array as cl_array
import pyopencl.cltypes as cltypes
import matplotlib.pyplot as plt
import time

from . import clUtils as clu
from .MMFF import MMFF
from .OpenCLBase import OpenCLBase
from .InteractionEnergy import load_xyz_with_REQs

REQ_DEFAULT = np.array([1.7, 0.1, 0.0, 0.0], dtype=np.float32)  # R, E, Q, padding

FOLDED_BASIS_MAX = 128
FOLDED_TYPES_MAX = 8
FOLDED_KERNEL_NAMES = ('orig', 'harmonics', 'workgroup')

verbose=False

def pack(iSys, source_array, target_buffer, queue):
    offset = iSys * source_array.size * source_array.dtype.itemsize
    cl.enqueue_copy(queue, target_buffer, source_array, offset=offset)

def copy(source, target, queue, iSys):
    pack(iSys, source, target, queue)

def copy_add(source, source_add, target, offset, queue):
    pass

def mat3_to_cl(mat3_np):
    return mat3_np.flatten().astype(np.float32)

def vec3_to_cl(vec3_np):
    return np.append(vec3_np, 0.0).astype(np.float32)

COMPONENT_PLQH = {
    'total':   np.array([1.0, 1.0, 1.0, 0.0], dtype=np.float32),
    'pauli':   np.array([1.0, 0.0, 0.0, 0.0], dtype=np.float32),
    'london':  np.array([0.0, 1.0, 0.0, 0.0], dtype=np.float32),
    'coulomb': np.array([0.0, 0.0, 1.0, 0.0], dtype=np.float32),
}

def half_step_from_coords(cs):
    cu = np.unique(np.round(np.asarray(cs, dtype=np.float64), 8))
    if len(cu) < 2:
        return 0.0
    ds = np.diff(np.sort(cu))
    ds = ds[ds > 1e-8]
    if len(ds) == 0:
        return 0.0
    return 0.5 * float(ds.min())

class MolecularDynamics(OpenCLBase):
    """
    Class for molecular dynamics simulations using OpenCL.
    
    This class inherits from OpenCLBase and implements specific functionality
    for molecular dynamics simulations using the relax_multi_mini.cl kernel.
    """
    
    def __init__(self, nloc=32, perBatch=10, debug_build_options=None, enable_nonbond=False):
        # Initialize the base class
        super().__init__(nloc=nloc, device_index=0)
        
        # Load the OpenCL program
        base_path = os.path.dirname(os.path.abspath(__file__))
        rel_path = "../../cpp/common_resources/cl/relax_multi.cl"
        if not self.load_program(rel_path=rel_path, base_path=base_path, bPrint=False, build_options=debug_build_options):
            exit(1)
        
        # Initialize other attributes that will be set in realloc
        self.nSystems         = 0
        self.mmff_list        = []
        self.MD_event_batch   = None
        self.perBatch         = perBatch
        self.nstep            = 1
        self.bPrintPackSystem = False
        self.enable_nonbond   = bool(enable_nonbond)
        self.use_padded_nonbond_small = False
        self._nb_small_ready  = False
        self._nb_prg_small    = None
        self._nb_bufs         = {}
        self.surface_atoms    = None
        self.surface_REQs     = None
        self.surface_lvec     = None
        self.surface_pos0     = np.zeros(4, dtype=np.float32)
        self.rigid_apos0      = None
        self.rigid_REQs0      = None
        self.folded_params    = None
        self.folded_type_ids  = None
        self.folded_fit_info  = None
        self.folded_kernel_kind = 'orig'

    def realloc(self, mmff, nSystems=1 ):
        """
        Reallocate buffers for the given number of systems based on the MMFF template.
        """
        # Store dimensions explicitly to avoid reference issues
        print(f"MolecularDynamics::realloc() natoms={mmff.natoms}, nvecs={mmff.nvecs}, nnode={mmff.nnode}")
        self.nSystems = nSystems
        self.mmff_list = [mmff] * nSystems  # Assuming all systems use the same MMFF parameters
        # Recreate command queue to ensure validity after reallocation/resizing
        self.queue = cl.CommandQueue(self.ctx)
        self.allocate_cl_buffers(mmff)
        self.allocate_host_buffers()
        self._nb_small_ready = False
        self._nb_bufs = {}

    def _setup_small_system_nonbond(self):
        if self._nb_small_ready:
            return
        if self.natoms >= self.nloc:
            return
        # Allocate padded buffers of size nloc per system to satisfy kernels with reqd_work_group_size(32)
        nSystems = self.nSystems
        nloc = self.nloc
        float_size = np.float32().itemsize
        int_size   = np.int32().itemsize
        mf = cl.mem_flags

        self._nb_bufs['apos_nb']      = cl.Buffer(self.ctx, mf.READ_WRITE, nSystems * nloc * 4 * float_size)
        self._nb_bufs['aforce_nb']    = cl.Buffer(self.ctx, mf.READ_WRITE, nSystems * nloc * 4 * float_size)
        self._nb_bufs['REQs_nb']      = cl.Buffer(self.ctx, mf.READ_ONLY,  nSystems * nloc * 4 * float_size)
        self._nb_bufs['neighs_nb']    = cl.Buffer(self.ctx, mf.READ_ONLY,  nSystems * nloc * 4 * int_size)
        self._nb_bufs['neighCell_nb'] = cl.Buffer(self.ctx, mf.READ_ONLY,  nSystems * nloc * 4 * int_size)
        if self.enable_nonbond:
            self._nb_bufs['excl_nb'] = cl.Buffer(self.ctx, mf.READ_ONLY, nSystems * nloc * 16 * int_size)

        src = r'''
        __kernel void copyAposToPadded(const int natoms, const int nvec_main, __global const float4* apos_main, __global float4* apos_nb){
            const int iG = get_global_id(0);
            const int iS = get_global_id(1);
            const int i0m = iS*nvec_main;
            const int i0n = iS*get_global_size(0);
            float4 p = (float4)(0.0f,0.0f,0.0f,0.0f);
            if(iG<natoms){ p = apos_main[i0m + iG]; }
            apos_nb[i0n + iG] = p;
        }

        __kernel void addForcesFromPadded(const int natoms, const int nvec_main, __global float4* aforce_main, __global const float4* aforce_nb){
            const int iG = get_global_id(0);
            const int iS = get_global_id(1);
            if(iG>=natoms) return;
            const int i0m = iS*nvec_main;
            const int i0n = iS*get_global_size(0);
            aforce_main[i0m + iG] += aforce_nb[i0n + iG];
        }
        '''
        self._nb_prg_small = cl.Program(self.ctx, src).build()

        # Prepare padded parameter arrays once (REQs / neighs are constant for the topology)
        natoms = self.natoms
        REQs_pad = np.zeros((nSystems, nloc, 4), dtype=np.float32)
        neighs_pad = np.full((nSystems, nloc, 4), -1, dtype=np.int32)
        neighCell_pad = np.zeros((nSystems, nloc, 4), dtype=np.int32)
        # Download current per-system topology params from mmff_list[0] and broadcast
        mm = self.mmff_list[0]
        REQs_pad[:, :natoms, :] = np.asarray(mm.REQs, dtype=np.float32)[None, :, :]
        neighs_pad[:, :natoms, :] = np.asarray(mm.neighs, dtype=np.int32)[None, :, :]
        neighCell_pad[:, :natoms, :] = np.asarray(mm.neighCell, dtype=np.int32)[None, :, :]
        cl.enqueue_copy(self.queue, self._nb_bufs['REQs_nb'], REQs_pad)
        cl.enqueue_copy(self.queue, self._nb_bufs['neighs_nb'], neighs_pad)
        cl.enqueue_copy(self.queue, self._nb_bufs['neighCell_nb'], neighCell_pad)
        if self.enable_nonbond and (self._nb_bufs.get('excl_nb') is not None):
            if getattr(mm, 'excl', None) is None:
                raise ValueError("_setup_small_system_nonbond(): enable_nonbond=True but mmff.excl is None")
            excl_pad = np.full((nSystems, nloc, 16), -1, dtype=np.int32)
            excl_pad[:, :natoms, :] = np.asarray(mm.excl, dtype=np.int32)[None, :, :]
            cl.enqueue_copy(self.queue, self._nb_bufs['excl_nb'], excl_pad)
        self.queue.finish()

        self._nb_small_ready = True

    def ensure_queue(self):
        if self.queue is None:
            self.queue = cl.CommandQueue(self.ctx)

    def allocate_host_buffers(self):
        self.atoms  = np.zeros((self.nSystems, self.nvecs, 4), dtype=np.float32)
        self.aforce = np.zeros((self.nSystems, self.nvecs, 4), dtype=np.float32)
        self.aforce_old = np.zeros((self.nSystems, self.nvecs, 4), dtype=np.float32)

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
        self.create_buffer('apos',       nSystems * nvecs * 4 * float_size, mf.READ_WRITE)
        self.create_buffer('aforce',     nSystems * nvecs * 4 * float_size, mf.READ_WRITE)
        # relax_multi.cl getMMFFf4() uses argument name 'fapos' for forces; alias it to the same buffer
        self.buffer_dict['fapos'] = self.buffer_dict['aforce']
        # relax_multi.cl nonbond kernels use argument names `atoms` and `forces`
        self.buffer_dict['atoms']  = self.buffer_dict['apos']
        self.buffer_dict['forces'] = self.buffer_dict['aforce']
        self.create_buffer('aforce_old', nSystems * nvecs * 4 * float_size, mf.READ_WRITE)
        self.create_buffer('avel',       nSystems * nvecs * 4 * float_size, mf.READ_WRITE)
        float4_size = 4 * float_size  # sizeof(float4) = 16 bytes
        self.create_buffer('fneigh',     nSystems * nnode * 4 * 2 * float4_size, mf.READ_WRITE)
        self.create_buffer('cvf',        nSystems * nvecs * 4 * float_size, mf.READ_WRITE)
        # Neighbor lists
        self.create_buffer('neighs',    nSystems * natoms * 4 * int_size, mf.READ_ONLY)
        self.create_buffer('neighCell', nSystems * natoms * 4 * int_size, mf.READ_ONLY)
        self.create_buffer('bkNeighs',  nSystems * nvecs  * 4 * int_size, mf.READ_ONLY)
        # Force field parameters
        self.create_buffer('REQs',     nSystems * natoms * 4 * float_size, mf.READ_ONLY)
        # relax_multi.cl getMMFFf4() uses argument name 'REQKs' for REQ parameters; alias it to the same buffer
        self.buffer_dict['REQKs'] = self.buffer_dict['REQs']

        # Non-bonded exclusions (optional; required by getNonBond_ex2)
        # NOTE: kernel uses EXCL_MAX=16 ints per atom
        if self.enable_nonbond:
            self.create_buffer('excl', nSystems * natoms * 16 * int_size, mf.READ_ONLY)
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
        # NOTE: updateAtomsMMFFf4 indexes constr/constrK by iaa=iG+iS*natoms (only atoms have constraints)
        # so the per-system stride must be natoms (not nvecs).
        self.create_buffer('constr',   nSystems * natoms * 4 * float_size, mf.READ_WRITE)
        self.create_buffer('constrK',  nSystems * natoms * 4 * float_size, mf.READ_WRITE)
        self.create_buffer('MDparams', nSystems * 4 * float_size, mf.READ_ONLY)
        self.create_buffer('TDrives',  nSystems * 4 * float_size, mf.READ_ONLY)
        # System interactions
        self.create_buffer('bboxes',   nSystems * mat3_size, mf.READ_ONLY)
        self.create_buffer('sysneighs',nSystems * int_size, mf.READ_ONLY)
        self.create_buffer('sysbonds', nSystems * 4 * float_size, mf.READ_ONLY)
        self.check_buf('surf_mpos', nSystems * 4 * float_size, mf.READ_ONLY)
        self.check_buf('surf_mdip', nSystems * 4 * float_size, mf.READ_ONLY)
        self.check_buf('surf_mQa',  nSystems * 4 * float_size, mf.READ_ONLY)
        self.check_buf('surf_mQb',  nSystems * 4 * float_size, mf.READ_ONLY)
        self.check_buf('surf_mQc',  nSystems * 4 * float_size, mf.READ_ONLY)
        self.check_buf('surf_qQa',  nSystems * 4 * float_size, mf.READ_ONLY)
        self.check_buf('surf_qQb',  nSystems * 4 * float_size, mf.READ_ONLY)
        self.check_buf('surf_qQc',  nSystems * 4 * float_size, mf.READ_ONLY)
        
        # Zero-initialize dynamical state buffers to prevent NaN from uninitialized GPU memory
        zero4 = np.zeros(1, dtype=np.float32)
        for bname in ('avel', 'cvf', 'aforce', 'aforce_old', 'fneigh', 'constr', 'constrK', 'TDrives'):
            buf = self.buffer_dict.get(bname)
            if buf is not None and buf.size > 0:
                cl.enqueue_fill_buffer(self.queue, buf, zero4, 0, buf.size)
        self.queue.finish()

    def set_pack_system_debug(self, enabled=True):
        """Enable/disable detailed parameter printing inside `pack_system()`."""
        self.bPrintPackSystem = bool(enabled)

    def _print_pack_system_params(self, iSys, mmff):
        print(f"[pack_system dbg] system={iSys}")
        print(f"  apos shape={mmff.apos.shape}\n{mmff.apos}")
        print(f"  fapos shape={getattr(mmff, 'fapos', None).shape if hasattr(mmff, 'fapos') else 'N/A'}")
        print(f"  REQs shape={mmff.REQs.shape}\n{mmff.REQs}")
        print(f"  apars shape={mmff.apars.shape}\n{mmff.apars}")
        print(f"  bLs shape={mmff.bLs.shape}\n{mmff.bLs}")
        print(f"  bKs shape={mmff.bKs.shape}\n{mmff.bKs}")
        print(f"  Ksp shape={mmff.Ksp.shape}\n{mmff.Ksp}")
        print(f"  Kpp shape={mmff.Kpp.shape}\n{mmff.Kpp}")
        print(f"  neighs shape={mmff.neighs.shape}\n{mmff.neighs}")
        print(f"  neighCell shape={mmff.neighCell.shape}\n{mmff.neighCell}")
        print(f"  back_neighs shape={mmff.back_neighs.shape}\n{mmff.back_neighs}")

    # NOTE: pack_system is defined further below with lvec/PBC support (the version near update_pbc_shifts).
    # The original version without lvec/PBC was here but is now superseded.

    def init_with_atoms(self, na=None, atoms=None, REQs=None, REQ_default=REQ_DEFAULT):
        """
        Initialize MolecularDynamics directly with atom positions and REQ parameters,
        without using MMFF object.
        
        Args:
            atoms: numpy array of shape (na, 3) with atom positions
            REQs: numpy array of shape (na, 4) with REQ parameters or None for defaults
            nSystems: number of systems
            nloc: local workgroup size
            
        Returns:
            Initialized MolecularDynamics instance
        """
        float_size = np.float32().itemsize

        if na is None:
            na = len(atoms)
        mf = cl.mem_flags

        # Initialize necessary attributes
        self.nSystems = 1
        self.natoms   = na
        self.nvecs    = na   # We don't have nodes or pi-orbitals, just atoms
        self.nnode    = 0         # No node atoms
        self.nDOFs    = (na, 0)  # (natoms, nnode)
        
        # Additional required attributes for initGridFF
        self.nPBC  = (0, 0, 0)  # No periodic boundary conditions
        self.npbc  = 0
        self.nbkng = na
        self.ncap  = 0
        self.ntors = 0
        # Convert atoms to float32 if needed
        
        if atoms is None:
            self.atoms = np.zeros((na, 4), dtype=np.float32)
        else:
            self.atoms = np.asarray(atoms, dtype=np.float32)
        if REQs is None:
            self.REQs    = np.zeros((na, 4), dtype=np.float32)
            self.REQs[:, 0] = REQ_default[0]
            self.REQs[:, 1] = REQ_default[1]
            self.REQs[:, 2] = REQ_default[2]
            self.REQs[:, 3] = REQ_default[3]
        else:
            self.REQs    = np.asarray(REQs, dtype=np.float32)
        self.aforce = np.zeros((na, 4), dtype=np.float32)

        self.create_buffer('apos',   na * 4 * float_size, mf.READ_WRITE)
        self.create_buffer('aforce', na * 4 * float_size, mf.READ_WRITE)
        # relax_multi.cl getMMFFf4() uses argument name 'fapos' for forces; alias it to the same buffer
        self.buffer_dict['fapos'] = self.buffer_dict['aforce']
        self.buffer_dict['atoms']  = self.buffer_dict['apos']
        self.buffer_dict['forces'] = self.buffer_dict['aforce']
        self.create_buffer('REQs',   na * 4 * float_size, mf.READ_ONLY)

        self.toGPU('apos',   self.atoms  )
        self.toGPU('aforce', self.aforce )
        self.toGPU('REQs',   self.REQs   )
        self.queue.finish()

    def init_rigid_molecule_batch(self, apos, REQs, nSystems=1):
        apos = np.asarray(apos, dtype=np.float32)
        REQs = np.asarray(REQs, dtype=np.float32)
        if apos.ndim != 2 or apos.shape[1] < 3:
            raise ValueError(f"init_rigid_molecule_batch(): apos must have shape (natoms,3+) got {apos.shape}")
        if REQs.ndim != 2 or REQs.shape[0] != apos.shape[0] or REQs.shape[1] < 4:
            raise ValueError(f"init_rigid_molecule_batch(): REQs must have shape (natoms,4+) matching apos, got {REQs.shape}")
        self.nSystems = int(nSystems)
        self.natoms = int(apos.shape[0])
        self.nvecs = self.natoms
        self.nnode = 0
        self.nDOFs = (self.natoms, 0)
        self.nPBC = (0, 0, 0)
        self.npbc = 0
        self.nbkng = self.natoms
        self.ncap = 0
        self.ntors = 0
        self.mmff_list = []
        self.rigid_apos0 = np.zeros((self.natoms, 4), dtype=np.float32)
        self.rigid_apos0[:, :3] = apos[:, :3]
        self.rigid_REQs0 = np.zeros((self.natoms, 4), dtype=np.float32)
        self.rigid_REQs0[:, :4] = REQs[:, :4]
        self.allocate_host_buffers()
        float_size = np.float32().itemsize
        mf = cl.mem_flags
        self.check_buf('apos', self.nSystems * self.nvecs * 4 * float_size, mf.READ_WRITE)
        self.check_buf('aforce', self.nSystems * self.nvecs * 4 * float_size, mf.READ_WRITE)
        self.buffer_dict['fapos'] = self.buffer_dict['aforce']
        self.buffer_dict['atoms'] = self.buffer_dict['apos']
        self.buffer_dict['forces'] = self.buffer_dict['aforce']
        self.check_buf('REQs', self.nSystems * self.natoms * 4 * float_size, mf.READ_ONLY)
        self.buffer_dict['REQKs'] = self.buffer_dict['REQs']
        self.check_buf('surf_mpos', self.nSystems * 4 * float_size, mf.READ_ONLY)
        self.check_buf('surf_mdip', self.nSystems * 4 * float_size, mf.READ_ONLY)
        self.check_buf('surf_mQa', self.nSystems * 4 * float_size, mf.READ_ONLY)
        self.check_buf('surf_mQb', self.nSystems * 4 * float_size, mf.READ_ONLY)
        self.check_buf('surf_mQc', self.nSystems * 4 * float_size, mf.READ_ONLY)
        self.check_buf('surf_qQa', self.nSystems * 4 * float_size, mf.READ_ONLY)
        self.check_buf('surf_qQb', self.nSystems * 4 * float_size, mf.READ_ONLY)
        self.check_buf('surf_qQc', self.nSystems * 4 * float_size, mf.READ_ONLY)
        self.check_buf('folded_coeffs', FOLDED_TYPES_MAX * FOLDED_BASIS_MAX * float_size, mf.READ_ONLY)
        self.check_buf('folded_kxyz',   FOLDED_BASIS_MAX * 4 * float_size, mf.READ_ONLY)
        self.check_buf('folded_atom_type', self.natoms * np.int32().itemsize, mf.READ_ONLY)
        reqs_all = np.broadcast_to(self.rigid_REQs0[None, :, :], (self.nSystems, self.natoms, 4)).copy()
        self.toGPU('REQs', reqs_all)
        self.toGPU('folded_atom_type', np.zeros(self.natoms, dtype=np.int32))
        self.toGPU('folded_coeffs', np.zeros((FOLDED_TYPES_MAX, FOLDED_BASIS_MAX), dtype=np.float32))
        self.toGPU('folded_kxyz', np.zeros((FOLDED_BASIS_MAX, 4), dtype=np.float32))
        cl.enqueue_fill_buffer(self.queue, self.buffer_dict['apos'], np.zeros(1, dtype=np.float32), 0, self.buffer_dict['apos'].size)
        cl.enqueue_fill_buffer(self.queue, self.buffer_dict['aforce'], np.zeros(1, dtype=np.float32), 0, self.buffer_dict['aforce'].size)
        self.setup_kernels()
        return self

    def setup_kernels(self):
        """
        Prepares the kernel arguments for all kernels by parsing their headers.
        Also sets up the work sizes for each kernel.
        """
        print("MolecularDynamics::setup_kernels()")
        # Get all work sizes at once
        self.get_work_sizes()
        self.init_kernel_params()

        def can_bind_kernel(kname):
            if kname not in self.kernelheaders:
                return False, [f'kernel:{kname}']
            missing = []
            for aname, typ in self.parse_kernel_header(self.kernelheaders[kname]):
                if typ == 0:
                    if (aname not in self.buffer_dict) and (aname not in self.kernel_params):
                        missing.append(aname)
                else:
                    if aname not in self.kernel_params:
                        missing.append(aname)
            return (len(missing) == 0), missing

        warned = set()
        def warn_skip(kname, missing):
            key = (kname, tuple(missing))
            if key in warned:
                return
            warned.add(key)
            print(f"warning: skipping {kname} because required buffers/params are not initialized: {', '.join(missing)}")
                
        # Generate kernel arguments (only for kernels present in the compiled source)
        # Wrap bonded-force kernels in try/except so init_with_atoms (no MMFF buffers) still works
        self.kernel_args_getMMFFf4 = None
        ok, missing = can_bind_kernel("getMMFFf4")
        if ok:
            self.kernel_args_getMMFFf4 = self.generate_kernel_args("getMMFFf4")
        elif "getMMFFf4" in self.kernelheaders:
            warn_skip("getMMFFf4", missing)
        self.kernel_args_getMMFFf4_rot = None
        if "getMMFFf4_rot" in self.kernelheaders:
            ok, missing = can_bind_kernel("getMMFFf4_rot")
            if ok:
                self.kernel_args_getMMFFf4_rot = self.generate_kernel_args("getMMFFf4_rot")
            else:
                warn_skip("getMMFFf4_rot", missing)
        
        self.kernel_args_getSurfFlat = None
        if "getSurfFlat" in self.kernelheaders:
            ok, missing = can_bind_kernel("getSurfFlat")
            if ok:
                self.kernel_args_getSurfFlat = self.generate_kernel_args("getSurfFlat")
            else:
                warn_skip("getSurfFlat", missing)
        self.kernel_args_getSurfMorse = None
        if "getSurfMorse" in self.kernelheaders:
            ok, missing = can_bind_kernel("getSurfMorse")
            if ok:
                self.kernel_args_getSurfMorse = self.generate_kernel_args("getSurfMorse")
            else:
                warn_skip("getSurfMorse", missing)
        self.kernel_args_getSurfFolded = None
        if "getSurfFolded" in self.kernelheaders:
            ok, missing = can_bind_kernel("getSurfFolded")
            if ok:
                self.kernel_args_getSurfFolded = self.generate_kernel_args("getSurfFolded")
            else:
                warn_skip("getSurfFolded", missing)
        self.kernel_args_getSurfFolded_harmonics = None
        if "getSurfFolded_harmonics" in self.kernelheaders:
            ok, missing = can_bind_kernel("getSurfFolded_harmonics")
            if ok:
                self.kernel_args_getSurfFolded_harmonics = self.generate_kernel_args("getSurfFolded_harmonics")
            else:
                warn_skip("getSurfFolded_harmonics", missing)
        self.kernel_args_getSurfFolded_workgroup = None
        if "getSurfFolded_workgroup" in self.kernelheaders:
            ok, missing = can_bind_kernel("getSurfFolded_workgroup")
            if ok:
                self.kernel_args_getSurfFolded_workgroup = self.generate_kernel_args("getSurfFolded_workgroup")
            else:
                warn_skip("getSurfFolded_workgroup", missing)
        self.kernel_args_sampleGridFF_Bspline_points = None
        if "sampleGridFF_Bspline_points" in self.kernelheaders:
            ok, missing = can_bind_kernel("sampleGridFF_Bspline_points")
            if ok:
                self.kernel_args_sampleGridFF_Bspline_points = self.generate_kernel_args("sampleGridFF_Bspline_points")
            else:
                warn_skip("sampleGridFF_Bspline_points", missing)
        # Non-bonded (optional)
        self.kernel_args_getNonBond = None
        self.kernel_args_getNonBond_ex2 = None
        if self.enable_nonbond:
            if "getNonBond_ex2" in self.kernelheaders:
                ok, missing = can_bind_kernel("getNonBond_ex2")
                if ok:
                    self.kernel_args_getNonBond_ex2 = self.generate_kernel_args("getNonBond_ex2")
                else:
                    warn_skip("getNonBond_ex2", missing)
            if (self.kernel_args_getNonBond_ex2 is None) and ("getNonBond" in self.kernelheaders):
                ok, missing = can_bind_kernel("getNonBond")
                if ok:
                    self.kernel_args_getNonBond = self.generate_kernel_args("getNonBond")
                else:
                    warn_skip("getNonBond", missing)
        # --- NOTE: grid-kernels are intialized in initGridFF()
        #self.kernel_args_getNonBond_GridFF_Bspline = self.generate_kernel_args("getNonBond_GridFF_Bspline")
        #self.kernel_args_getNonBond_GridFF_Bspline_tex = self.generate_kernel_args("getNonBond_GridFF_Bspline_tex")
        self.kernel_args_updateAtomsMMFFf4 = None
        ok, missing = can_bind_kernel("updateAtomsMMFFf4")
        if ok:
            self.kernel_args_updateAtomsMMFFf4 = self.generate_kernel_args("updateAtomsMMFFf4")
        elif "updateAtomsMMFFf4" in self.kernelheaders:
            warn_skip("updateAtomsMMFFf4", missing)
        # New propagator variants (optional)
        self.kernel_args_updateAtomsMMFFf4_rot = None
        if "updateAtomsMMFFf4_rot" in self.kernelheaders:
            ok, missing = can_bind_kernel("updateAtomsMMFFf4_rot")
            if ok:
                self.kernel_args_updateAtomsMMFFf4_rot = self.generate_kernel_args("updateAtomsMMFFf4_rot")
            else:
                warn_skip("updateAtomsMMFFf4_rot", missing)
        #self.kernel_args_updateAtomsMMFFf4_RATTLE = self.generate_kernel_args("updateAtomsMMFFf4_RATTLE")
        self.kernel_args_cleanForceMMFFf4 = None
        ok, missing = can_bind_kernel("cleanForceMMFFf4")
        if ok:
            self.kernel_args_cleanForceMMFFf4 = self.generate_kernel_args("cleanForceMMFFf4")
        elif "cleanForceMMFFf4" in self.kernelheaders:
            warn_skip("cleanForceMMFFf4", missing)
        self.kernel_args_runMD = None
        if "runMD" in self.kernelheaders:
            ok, missing = can_bind_kernel("runMD")
            if ok:
                self.kernel_args_runMD = self.generate_kernel_args("runMD")
            else:
                warn_skip("runMD", missing)

    def init_kernel_params(self):
        """
        Initialize a dictionary of standard kernel parameters.
        This provides default values for common parameters used in kernels.
        """
        print("MolecularDynamics::init_kernel_params()")
        # Create a dictionary to store kernel parameters
        self.kernel_params = {
            # Common dimension parameters
            'nDOFs':        np.array([self.natoms, self.nnode, 0, self.perBatch], dtype=np.int32),
            'mask':         np.array([1, 1, 1, 1],         dtype=np.int32),
            'nPBC':         np.array(self.nPBC+(0,),       dtype=np.int32),
            'GFFParams':    np.array([0.0, 0.0,  0.0, 0.0], dtype=np.float32),
            'MDparams':     np.array([0.1, 0.05, 0.0, 0.0], dtype=np.float32),
            # Common scalar parameters
            'npbc':         np.int32(self.npbc),
            'bSubtractVdW': np.int32(0),
            'grid_ns':      np.array([0,0,0,0], dtype=np.int32),
            'grid_invStep': np.array([0.0,0.0,0.0,0.0], dtype=np.float32),
            'grid_p0':      np.array([0.0,0.0,0.0,0.0], dtype=np.float32),
            # Surface parameters
            'surf_pos0':    np.array([0.0,0.0,0.0,0.0], dtype=np.float32),
            'surf_normal':  np.array([0.0,0.0,1.0,0.0], dtype=np.float32),
            'surf_REQ':     np.array([1.0,1.0,0.0,0.0], dtype=np.float32),
            'surf_param':   np.array([1.0,1.0,0.0,0.0], dtype=np.float32), # K, mode, 0, 0
            'lvec':         np.zeros((3,4), dtype=np.float32),
            'pos0':         np.zeros(4, dtype=np.float32),
            'PLQH':         np.array([1.0, 1.0, 1.0, 0.0], dtype=np.float32),
        }

        # relax_multi.cl uses different arg names for the same dimensions
        # - updateAtomsMMFFf4 uses `n`
        # - getNonBond uses `ns`
        self.kernel_params['n']  = self.kernel_params['nDOFs']
        self.kernel_params['ns'] = self.kernel_params['nDOFs']

    def _surface_supercell_moments(self, apos, qs, lvec, nPBC):
        a = np.array(lvec[0], dtype=np.float64)
        b = np.array(lvec[1], dtype=np.float64)
        c = np.array(lvec[2], dtype=np.float64)
        shifts = []
        for iz in range(-int(nPBC[2]), int(nPBC[2])+1):
            for iy in range(-int(nPBC[1]), int(nPBC[1])+1):
                for ix in range(-int(nPBC[0]), int(nPBC[0])+1):
                    shifts.append(ix*a + iy*b + iz*c)
        pts = np.concatenate([apos + sh[None, :] for sh in shifts], axis=0)
        qq = np.tile(qs, len(shifts))
        qtot = float(qq.sum())
        if abs(qtot) > 1e-10:
            ctr = (qq[:, None] * pts).sum(axis=0) / qtot
        else:
            ctr = 0.5 * (pts.min(axis=0) + pts.max(axis=0))
        d = pts - ctr[None, :]
        mu = (qq[:, None] * d).sum(axis=0)
        Q = np.zeros((3,3), dtype=np.float64)
        r2 = np.sum(d*d, axis=1)
        for i in range(3):
            for j in range(3):
                Q[i, j] = np.sum(qq * (3.0*d[:, i]*d[:, j] - (r2 if i == j else 0.0)))
        return ctr, qtot, mu, Q

    def set_surface(self, surf_xyz, nPBC=(1,1,0), pos0=None, alpha_morse=1.6, r_damp=0.0, bMacro=True, type_map=None):
        apos, REQs, enames, Zs, lvec = load_xyz_with_REQs(surf_xyz, type_map=type_map)
        if lvec is None:
            raise ValueError(f"Surface file {surf_xyz} must contain lattice vectors in comment line")
        self.surface_atoms = np.zeros((len(apos), 4), dtype=np.float32)
        self.surface_atoms[:, :3] = apos.astype(np.float32)
        self.surface_enames = list(enames)
        self.surface_REQs = np.ascontiguousarray(REQs, dtype=np.float32)
        self.surface_lvec = np.zeros((3,4), dtype=np.float32)
        self.surface_lvec[:, :3] = lvec.astype(np.float32)
        self.surface_pos0 = np.array([0.0, 0.0, 0.0, 0.0] if pos0 is None else [pos0[0], pos0[1], pos0[2], 0.0], dtype=np.float32)
        float_size = np.float32().itemsize
        self.check_buf('atoms_s', len(apos) * 4 * float_size, cl.mem_flags.READ_ONLY)
        self.check_buf('REQ_s',   len(apos) * 4 * float_size, cl.mem_flags.READ_ONLY)
        self.check_buf('surf_mpos', self.nSystems * 4 * float_size, cl.mem_flags.READ_ONLY)
        self.check_buf('surf_mdip', self.nSystems * 4 * float_size, cl.mem_flags.READ_ONLY)
        self.check_buf('surf_mQa',  self.nSystems * 4 * float_size, cl.mem_flags.READ_ONLY)
        self.check_buf('surf_mQb',  self.nSystems * 4 * float_size, cl.mem_flags.READ_ONLY)
        self.check_buf('surf_mQc',  self.nSystems * 4 * float_size, cl.mem_flags.READ_ONLY)
        self.check_buf('surf_qQa',  self.nSystems * 4 * float_size, cl.mem_flags.READ_ONLY)
        self.check_buf('surf_qQb',  self.nSystems * 4 * float_size, cl.mem_flags.READ_ONLY)
        self.check_buf('surf_qQc',  self.nSystems * 4 * float_size, cl.mem_flags.READ_ONLY)
        self.toGPU('atoms_s', self.surface_atoms)
        self.toGPU('REQ_s', self.surface_REQs)
        ctr, qtot, mu, Q = self._surface_supercell_moments(apos.astype(np.float64), REQs[:,2].astype(np.float64), lvec.astype(np.float64), nPBC)
        a = np.array(lvec[0], dtype=np.float64)
        b = np.array(lvec[1], dtype=np.float64)
        xmin0, xmax0 = float(apos[:,0].min()), float(apos[:,0].max())
        ymin0, ymax0 = float(apos[:,1].min()), float(apos[:,1].max())
        hx = half_step_from_coords(apos[:,0])
        hy = half_step_from_coords(apos[:,1])
        xmin = xmin0 - hx - float(nPBC[0]) * float(np.linalg.norm(a))
        xmax = xmax0 + hx + float(nPBC[0]) * float(np.linalg.norm(a))
        ymin = ymin0 - hy - float(nPBC[1]) * float(np.linalg.norm(b))
        ymax = ymax0 + hy + float(nPBC[1]) * float(np.linalg.norm(b))
        area = float(np.linalg.norm(np.cross(a, b)))
        if area < 1e-12:
            raise ValueError(f"set_surface(): invalid cell area {area}")
        zs = apos[:,2].astype(np.float64)
        qs = REQs[:,2].astype(np.float64)
        zuniq = []
        layers = []
        sigmas = []
        tol = 1e-4
        for i in np.argsort(zs):
            z = zs[i]
            if (not zuniq) or (abs(z-zuniq[-1]) > tol):
                zuniq.append(z)
        if bMacro and (len(zuniq) > 3):
            raise ValueError(f"set_surface(): getSurfMorse rectangle macro currently supports up to 3 z-layers, got {len(zuniq)}")
        cx = 0.5 * (xmin + xmax)
        cy = 0.5 * (ymin + ymax)
        for z in zuniq:
            m = np.abs(zs-z) < tol
            dx = apos[m,0].astype(np.float64) - cx
            dy = apos[m,1].astype(np.float64) - cy
            dz = apos[m,2].astype(np.float64) - z
            qq = qs[m]
            mu_layer = np.array([np.sum(qq*dx), np.sum(qq*dy), np.sum(qq*dz)], dtype=np.float64) / area
            layers.append(np.array([mu_layer[0], mu_layer[1], mu_layer[2], z], dtype=np.float32))
            sigmas.append(float(np.sum(qq) / area))
        qlayers = []
        for z in zuniq:
            m = np.abs(zs-z) < tol
            dx = apos[m,0].astype(np.float64) - cx
            dy = apos[m,1].astype(np.float64) - cy
            qq = qs[m]
            qxx = float(np.sum(qq * (2.0*dx*dx - dy*dy)) / area)
            qyy = float(np.sum(qq * (2.0*dy*dy - dx*dx)) / area)
            qxy = float(np.sum(qq * (3.0*dx*dy)) / area)
            qlayers.append(np.array([qxx, qxy, qyy, z], dtype=np.float32))
        while len(layers) < 3:
            layers.append(np.zeros(4, dtype=np.float32))
        while len(sigmas) < 3:
            sigmas.append(0.0)
        while len(qlayers) < 3:
            qlayers.append(np.zeros(4, dtype=np.float32))
        mpos = np.tile(np.array([xmin, xmax, ymin, ymax], dtype=np.float32), (self.nSystems, 1))
        mdip = np.tile(layers[0], (self.nSystems, 1))
        mQa  = np.tile(layers[1], (self.nSystems, 1))
        mQb  = np.tile(layers[2], (self.nSystems, 1))
        mQc  = np.tile(np.array([sigmas[0], sigmas[1], sigmas[2], qtot], dtype=np.float32), (self.nSystems, 1))
        qQa  = np.tile(qlayers[0], (self.nSystems, 1))
        qQb  = np.tile(qlayers[1], (self.nSystems, 1))
        qQc  = np.tile(qlayers[2], (self.nSystems, 1))
        self.toGPU('surf_mpos', mpos)
        self.toGPU('surf_mdip', mdip)
        self.toGPU('surf_mQa', mQa)
        self.toGPU('surf_mQb', mQb)
        self.toGPU('surf_mQc', mQc)
        self.toGPU('surf_qQa', qQa)
        self.toGPU('surf_qQb', qQb)
        self.toGPU('surf_qQc', qQc)
        self.kernel_params['ns'] = np.array([self.natoms, self.nnode, len(apos), self.perBatch], dtype=np.int32)
        self.kernel_params['nPBC'] = np.array([int(nPBC[0]), int(nPBC[1]), int(nPBC[2]), 0], dtype=np.int32)
        self.kernel_params['lvec'] = self.surface_lvec.copy()
        self.kernel_params['pos0'] = self.surface_pos0.copy()
        self.kernel_params['GFFParams'] = np.array([r_damp, alpha_morse, 1.0 if bMacro else 0.0, float(len(zuniq))], dtype=np.float32)
        if 'getSurfMorse' in self.kernelheaders:
            self.kernel_args_getSurfMorse = self.generate_kernel_args('getSurfMorse', bPrint=False)
        return {'apos': apos, 'REQs': REQs, 'lvec': lvec, 'center': ctr, 'qtot': qtot, 'mu': mu, 'Q': Q, 'layers': layers, 'sigmas': sigmas, 'bounds': (xmin, xmax, ymin, ymax)}
        
    def get_work_sizes(self):
        """
        Generate standard work sizes based on current system dimensions.
        
        Returns:
            dict: Dictionary containing sz_na, sz_node, sz_nvec, sz_loc
        """
        # Default to the standard work size parameters
        self.sz_loc = (self.nloc, 1)
        self.sz_na   = (clu.roundup_global_size(self.natoms, self.nloc), self.nSystems)
        self.sz_node = (clu.roundup_global_size(self.nnode,  self.nloc), self.nSystems)
        self.sz_nvec = (clu.roundup_global_size(self.nvecs,  self.nloc), self.nSystems)
        
        # Return all sizes, let the caller decide which to use
        return {
            'sz_na':   self.sz_na,
            'sz_node': self.sz_node,
            'sz_nvec': self.sz_nvec,
            'sz_loc':  self.sz_loc
        }

    def upload_all_systems(self):
        """Uploads data for all systems to the GPU."""
        for sys_idx in range(self.nSystems):
            self.pack_system(sys_idx, self.mmff_list[sys_idx])
        #print("MolecularDynamics::upload_all_systems() DONE")

    def run_getNonBond(self):
        self.prg.getNonBond(self.queue, self.sz_na, self.sz_loc, *self.kernel_args_getNonBond)
        self.queue.finish()

    def run_getNonBond_ex2(self):
        self.prg.getNonBond_ex2(self.queue, self.sz_na, self.sz_loc, *self.kernel_args_getNonBond_ex2)
        self.queue.finish()
    
    def run_getNonBond_GridFF_Bspline(self):
        self.ensure_queue()
        if self.use_padded_nonbond_small and (self.natoms < self.nloc):
            self._setup_small_system_nonbond()
            # Copy current positions into padded buffer
            sz = (self.nloc, self.nSystems)
            self._nb_prg_small.copyAposToPadded(self.queue, sz, self.sz_loc, np.int32(self.natoms), np.int32(self.nvecs), self.buffer_dict['apos'], self._nb_bufs['apos_nb'])
            # Zero padded forces
            cl.enqueue_fill_buffer(self.queue, self._nb_bufs['aforce_nb'], np.zeros(1, dtype=np.float32), 0, self._nb_bufs['aforce_nb'].size)
            # Run GridFF kernel on padded buffers
            args = list(self.kernel_args_getNonBond_GridFF_Bspline)
            # Padded buffers are laid out with stride==nloc; set nnode=0 so kernel uses nvec==natoms==nloc
            args[0] = cltypes.make_int4(self.nloc, 0, self.nloc, 0)
            args[1] = self._nb_bufs['apos_nb']
            args[2] = self._nb_bufs['aforce_nb']
            args[3] = self._nb_bufs['REQs_nb']
            args[4] = self._nb_bufs['neighs_nb']
            args[5] = self._nb_bufs['neighCell_nb']
            self.prg.getNonBond_GridFF_Bspline(self.queue, sz, self.sz_loc, *args)
            # Accumulate first natoms forces back to main aforce
            self._nb_prg_small.addForcesFromPadded(self.queue, sz, self.sz_loc, np.int32(self.natoms), np.int32(self.nvecs), self.buffer_dict['aforce'], self._nb_bufs['aforce_nb'])
            self.queue.finish()
        else:
            self.prg.getNonBond_GridFF_Bspline(self.queue, self.sz_na, self.sz_loc, *self.kernel_args_getNonBond_GridFF_Bspline)
            self.queue.finish()

    def run_getNonBond_GridFF_Bspline_tex(self):
        self.ensure_queue()
        if self.use_padded_nonbond_small and (self.natoms < self.nloc):
            self._setup_small_system_nonbond()
            sz = (self.nloc, self.nSystems)
            self._nb_prg_small.copyAposToPadded(self.queue, sz, self.sz_loc, np.int32(self.natoms), np.int32(self.nvecs), self.buffer_dict['apos'], self._nb_bufs['apos_nb'])
            cl.enqueue_fill_buffer(self.queue, self._nb_bufs['aforce_nb'], np.zeros(1, dtype=np.float32), 0, self._nb_bufs['aforce_nb'].size)
            args = list(self.kernel_args_getNonBond_GridFF_Bspline_tex)
            args[0] = cltypes.make_int4(self.nloc, 0, self.nloc, 0)  # ns (natoms padded, nnode, nvec padded)
            args[1] = self._nb_bufs['apos_nb']
            args[2] = self._nb_bufs['aforce_nb']
            args[3] = self._nb_bufs['REQs_nb']
            args[4] = self._nb_bufs['neighs_nb']
            args[5] = self._nb_bufs['neighCell_nb']
            self.prg.getNonBond_GridFF_Bspline_tex(self.queue, sz, self.sz_loc, *args)
            self._nb_prg_small.addForcesFromPadded(self.queue, sz, self.sz_loc, np.int32(self.natoms), np.int32(self.nvecs), self.buffer_dict['aforce'], self._nb_bufs['aforce_nb'])
            self.queue.finish()
        else:
            self.prg.getNonBond_GridFF_Bspline_tex(self.queue, self.sz_na, self.sz_loc, *self.kernel_args_getNonBond_GridFF_Bspline_tex)
            self.queue.finish()

    def run_getNonBond_GridFF_Bspline_ex2(self):
        self.ensure_queue()
        if self.kernel_args_getNonBond_GridFF_Bspline_ex2 is None:
            raise RuntimeError("run_getNonBond_GridFF_Bspline_ex2(): kernel args not initialized; call initGridFF(..., bKernels=True) and ensure enable_nonbond=True")
        if self.use_padded_nonbond_small and (self.natoms < self.nloc):
            self._setup_small_system_nonbond()
            if self._nb_bufs.get('excl_nb') is None:
                raise RuntimeError("run_getNonBond_GridFF_Bspline_ex2(): enable_nonbond=True required for excl buffer")
            sz = (self.nloc, self.nSystems)
            self._nb_prg_small.copyAposToPadded(self.queue, sz, self.sz_loc, np.int32(self.natoms), np.int32(self.nvecs), self.buffer_dict['apos'], self._nb_bufs['apos_nb'])
            cl.enqueue_fill_buffer(self.queue, self._nb_bufs['aforce_nb'], np.zeros(1, dtype=np.float32), 0, self._nb_bufs['aforce_nb'].size)
            args = list(self.kernel_args_getNonBond_GridFF_Bspline_ex2)
            # Padded buffers are laid out with stride==nloc; set nnode=0 so kernel uses nvec==natoms==nloc
            args[0] = cltypes.make_int4(self.nloc, 0, self.nloc, 0)
            args[1] = self._nb_bufs['apos_nb']
            args[2] = self._nb_bufs['aforce_nb']
            args[3] = self._nb_bufs['REQs_nb']
            args[4] = self._nb_bufs['excl_nb']
            self.prg.getNonBond_GridFF_Bspline_ex2(self.queue, sz, self.sz_loc, *args)
            self._nb_prg_small.addForcesFromPadded(self.queue, sz, self.sz_loc, np.int32(self.natoms), np.int32(self.nvecs), self.buffer_dict['aforce'], self._nb_bufs['aforce_nb'])
            self.queue.finish()
        else:
            self.prg.getNonBond_GridFF_Bspline_ex2(self.queue, self.sz_na, self.sz_loc, *self.kernel_args_getNonBond_GridFF_Bspline_ex2)
            self.queue.finish()

    def run_getMMFFf4(self):
        self.ensure_queue()
        self.prg.getMMFFf4(self.queue, self.sz_node, self.sz_loc, *self.kernel_args_getMMFFf4)
        self.queue.finish()
    
    def run_getMMFFf4_rot(self):
        self.prg.getMMFFf4_rot(self.queue, self.sz_node, self.sz_loc, *self.kernel_args_getMMFFf4_rot)
        self.queue.finish()
    
    def run_updateAtomsMMFFf4(self):
        self.ensure_queue()
        # IMPORTANT: updateAtomsMMFFf4 indexes apos/avel/aforce using nvec=natoms+nnode.
        # When nvecs < nloc and global size is padded to nloc, iG can exceed nvec and write out-of-bounds
        # into the next system's memory, corrupting geometry (observed as stretched molecules in multi-system runs).
        sz = (int(self.nvecs), int(self.nSystems))
        self.prg.updateAtomsMMFFf4(self.queue, sz, None, *self.kernel_args_updateAtomsMMFFf4)
        self.queue.finish()
    
    def run_updateAtomsMMFFf4_rot(self):
        self.prg.updateAtomsMMFFf4_rot(self.queue, self.sz_nvec, self.sz_loc, *self.kernel_args_updateAtomsMMFFf4_rot)
        self.queue.finish()
    
    def run_updateAtomsMMFFf4_RATTLE(self):
        self.prg.updateAtomsMMFFf4_RATTLE(self.queue, self.sz_nvec, self.sz_loc, *self.kernel_args_updateAtomsMMFFf4_RATTLE)
        self.queue.finish()
    
    def run_cleanForceMMFFf4(self):
        self.ensure_queue()
        # IMPORTANT: cleanForceMMFFf4 indexes `aforce` using `nvec = natoms+nnode`.
        # Therefore it must be launched with global size covering nvec, not natoms.
        # Also must NOT pad to nloc when nvec < nloc, otherwise iG spans beyond nvec and writes OOB across systems.
        sz = (int(self.nvecs), int(self.nSystems))
        self.prg.cleanForceMMFFf4(self.queue, sz, None, *self.kernel_args_cleanForceMMFFf4)
        # cleanForceMMFFf4 kernel only zeros sigma portion of fneigh (nnode*4 entries).
        # The pi portion (nnode*4 to nnode*8-1) may contain stale data.
        # Zero the full fneigh buffer to prevent accumulation → NaN.
        fneigh_buf = self.buffer_dict.get('fneigh')
        if fneigh_buf is not None:
            import pyopencl as cl
            cl.enqueue_fill_buffer(self.queue, fneigh_buf, np.zeros(1, dtype=np.float32), 0, fneigh_buf.size)
        self.queue.finish()

    def run_getSurfFlat(self):
        # Update scalar arguments from kernel_params
        # args 5,6,7,8 are constants: surf_pos0, surf_normal, surf_REQ, surf_param
        # We need to ensure they are passed as vectors/scalars correctly
        # The kernel arg generator should have picked them up from kernel_params or overrides
        # We can update the values in kernel_params before calling this if needed.
        
        # We need to make sure we are passing the values from self.kernel_params
        # because generate_kernel_args binds the values at generation time if they were present
        # BUT generate_kernel_args uses the LIST/VALUE from kernel_params. 
        # If we updated the array content in place (numpy), it might be fine if the arg list holds reference.
        # But for scalars/vectors in pyopencl, they are usually passed by value or as specific types.
        
        # Safest is to regenerate args or manually override specific args.
        # But let's assume kernel_params values are used.
        
        # Re-generate args to be sure we pick up current kernel_params values
        # Optimization: only do this if params changed. For now, just do it.
        # Or better: check if generate_kernel_args stores references. 
        # It appends self.kernel_params[aname]. 
        # For numpy arrays (used for float4), it passes the array object. 
        # PyOpenCL handles numpy arrays as vector types by value? No, it expects cl types or correct numpy types.
        
        # Let's trust setup_kernels did the job, but we might need to update the args list if we change params.
        # For this specific kernel, we want to allow changing surface params dynamically.
        
        # Let's create a temporary args list with current params
        args = self.kernel_args_getSurfFlat[:]
        # Map of arg name to index is not stored. 
        # We can re-generate.
        args = self.generate_kernel_args("getSurfFlat", bPrint=False)
        
        self.prg.getSurfFlat(self.queue, self.sz_na, self.sz_loc, *args)
        self.queue.finish()

    def run_getSurfMorse(self, nSystems=None):
        if self.kernel_args_getSurfMorse is None:
            self.kernel_args_getSurfMorse = self.generate_kernel_args('getSurfMorse', bPrint=False)
        nSystems = self.nSystems if nSystems is None else int(nSystems)
        lx = min(int(self.nloc), int(self.natoms))
        while (lx > 1) and (int(self.natoms) % lx != 0):
            lx -= 1
        sz = (int(self.natoms), nSystems)
        loc = (int(lx), 1)
        self.prg.getSurfMorse(self.queue, sz, loc, *self.kernel_args_getSurfMorse)
        self.queue.finish()

    def run_getSurfFolded(self, nSystems=None):
        nSystems = self.nSystems if nSystems is None else int(nSystems)
        kind = getattr(self, 'folded_kernel_kind', 'orig')
        if kind == 'orig':
            if self.kernel_args_getSurfFolded is None:
                self.kernel_args_getSurfFolded = self.generate_kernel_args('getSurfFolded', bPrint=False)
            lx = min(int(self.nloc), int(self.natoms))
            while (lx > 1) and (int(self.natoms) % lx != 0):
                lx -= 1
            sz = (int(self.natoms), nSystems)
            loc = (int(lx), 1)
            self.prg.getSurfFolded(self.queue, sz, loc, *self.kernel_args_getSurfFolded)
        elif kind == 'harmonics':
            if self.kernel_args_getSurfFolded_harmonics is None:
                self.kernel_args_getSurfFolded_harmonics = self.generate_kernel_args('getSurfFolded_harmonics', bPrint=False)
            lx = min(int(self.nloc), int(self.natoms))
            while (lx > 1) and (int(self.natoms) % lx != 0):
                lx -= 1
            sz = (int(self.natoms), nSystems)
            loc = (int(lx), 1)
            self.prg.getSurfFolded_harmonics(self.queue, sz, loc, *self.kernel_args_getSurfFolded_harmonics)
        elif kind == 'workgroup':
            if self.kernel_args_getSurfFolded_workgroup is None:
                self.kernel_args_getSurfFolded_workgroup = self.generate_kernel_args('getSurfFolded_workgroup', bPrint=False)
            lx = 64
            sz = (int(math.ceil(float(self.natoms) / lx) * lx), nSystems)
            loc = (lx, 1)
            self.prg.getSurfFolded_workgroup(self.queue, sz, loc, *self.kernel_args_getSurfFolded_workgroup)
        else:
            raise ValueError(f"run_getSurfFolded(): unknown folded kernel '{kind}', expected one of {FOLDED_KERNEL_NAMES}")
        self.queue.finish()

    def run_sampleGridFF_Bspline_points(self, nSystems=None, PLQH=None):
        if self.kernel_args_sampleGridFF_Bspline_points is None:
            self.kernel_args_sampleGridFF_Bspline_points = self.generate_kernel_args('sampleGridFF_Bspline_points', bPrint=False)
        args = self.kernel_args_sampleGridFF_Bspline_points
        if PLQH is not None:
            self.kernel_params['PLQH'] = np.asarray(PLQH, dtype=np.float32).reshape(4)
            self.kernel_args_sampleGridFF_Bspline_points = self.generate_kernel_args('sampleGridFF_Bspline_points', bPrint=False)
            args = self.kernel_args_sampleGridFF_Bspline_points
        nSystems = self.nSystems if nSystems is None else int(nSystems)
        lx = int(self.nloc)
        gx = int(math.ceil(float(self.natoms) / lx) * lx)
        sz = (gx, nSystems)
        loc = (lx, 1)
        self.prg.sampleGridFF_Bspline_points(self.queue, sz, loc, *args)
        self.queue.finish()

    def _folded_unique_types(self, REQs, decimals=6):
        keys = np.round(np.asarray(REQs[:, :4], dtype=np.float64), decimals).astype(np.float32)
        uniq = []
        idxs = np.empty(len(keys), dtype=np.int32)
        mapping = {}
        for i, row in enumerate(keys):
            key = tuple(row.tolist())
            it = mapping.get(key, None)
            if it is None:
                it = len(uniq)
                mapping[key] = it
                uniq.append(np.array(row, dtype=np.float32))
            idxs[i] = it
        return np.array(uniq, dtype=np.float32), idxs

    def _component_plqh(self, component=None):
        if component is None:
            component = 'total'
        key = str(component).lower()
        if key not in COMPONENT_PLQH:
            raise ValueError(f"_component_plqh(): unknown component '{component}', expected one of {sorted(COMPONENT_PLQH.keys())}")
        return COMPONENT_PLQH[key].copy()

    def _resolve_plqh(self, component=None, PLQH=None):
        if PLQH is not None:
            arr = np.asarray(PLQH, dtype=np.float32).reshape(4)
            return arr.copy()
        return self._component_plqh(component)

    def _set_folded_coefficients(self, coeffs):
        coeffs = np.asarray(coeffs, dtype=np.float32)
        if self.folded_params is None:
            raise ValueError('_set_folded_coefficients(): fit_folded_surface_basis() must be called first')
        ntypes = int(coeffs.shape[0])
        kind = getattr(self, 'folded_kernel_kind', 'orig')
        if kind == 'orig':
            nbasis = int(self.folded_params['basis_params'].shape[0])
            coeff_pad = np.zeros((FOLDED_TYPES_MAX, FOLDED_BASIS_MAX), dtype=np.float32)
            coeff_pad[:ntypes, :nbasis] = coeffs[:, :nbasis]
            kxyz_pad = np.zeros((FOLDED_BASIS_MAX, 4), dtype=np.float32)
            kxyz_pad[:nbasis, :] = np.asarray(self.folded_params['basis_params'][:, :4], dtype=np.float32)
            self.toGPU('folded_coeffs', coeff_pad)
            self.toGPU('folded_kxyz', kxyz_pad)
            self.kernel_params['folded_meta'] = np.array([nbasis, ntypes, 0, 0], dtype=np.int32)
            return
        nu = int(self.folded_params['nu']); nv = int(self.folded_params['nv']); nz = int(self.folded_params['nz'])
        if nu > 4 or nv > 4 or nz > 8:
            raise ValueError(f"_set_folded_coefficients(): optimized folded kernels support at most nu<=4 nv<=4 nz<=8, got nu={nu} nv={nv} nz={nz}")
        if kind == 'workgroup' and nu != nv:
            raise ValueError(f"_set_folded_coefficients(): workgroup kernel requires nu==nv, got nu={nu} nv={nv}")
        coeff_tensor = coeffs[:, :nu*nv*nz].reshape(ntypes, nu, nv, nz).transpose(0, 3, 2, 1)
        coeff_pad = np.zeros((FOLDED_TYPES_MAX, FOLDED_BASIS_MAX), dtype=np.float32)
        coeff_pad[:ntypes, :nu*nv*nz] = coeff_tensor.reshape(ntypes, nu*nv*nz)
        self.toGPU('folded_coeffs', coeff_pad)
        basis_params = np.asarray(self.folded_params['basis_params'], dtype=np.float32).reshape(nu, nv, nz, 4)
        u_params = np.zeros((nu, 4), dtype=np.float32); u_params[:, 0] = basis_params[:, 0, 0, 0]
        v_params = np.zeros((nv, 4), dtype=np.float32); v_params[:, 1] = basis_params[0, :, 0, 1]
        z_params = np.zeros((nz, 4), dtype=np.float32); z_params[:, 2:] = basis_params[0, 0, :, 2:4]
        if kind == 'harmonics':
            kxyz = np.vstack([u_params, v_params, z_params])
            self.kernel_params['folded_meta'] = np.array([nu, nv, nz, ntypes], dtype=np.int32)
        elif kind == 'workgroup':
            kxyz = np.vstack([u_params, z_params])
            self.kernel_params['folded_meta'] = np.array([nu, nz, ntypes, 0], dtype=np.int32)
        else:
            raise ValueError(f"_set_folded_coefficients(): unknown folded kernel '{kind}', expected one of {FOLDED_KERNEL_NAMES}")
        kxyz_pad = np.zeros((FOLDED_BASIS_MAX, 4), dtype=np.float32)
        kxyz_pad[:len(kxyz), :] = kxyz
        self.toGPU('folded_kxyz', kxyz_pad)

    def set_folded_kernel_kind(self, kind='orig'):
        kind = str(kind).lower()
        if kind not in FOLDED_KERNEL_NAMES:
            raise ValueError(f"set_folded_kernel_kind(): unknown kind '{kind}', expected one of {FOLDED_KERNEL_NAMES}")
        self.folded_kernel_kind = kind

    def _build_folded_basis_params(self, nu=4, nv=4, nz=4, z0=0.0, z_scale=0.75):
        params = []
        alphas = z_scale * (1.0 + np.arange(int(nz), dtype=np.float32))
        for iu in range(int(nu)):
            for iv in range(int(nv)):
                for az in alphas:
                    params.append((float(iu), float(iv), float(az), float(z0)))
        if len(params) > FOLDED_BASIS_MAX:
            raise ValueError(f'_build_folded_basis_params(): nbasis={len(params)} exceeds FOLDED_BASIS_MAX={FOLDED_BASIS_MAX}')
        return np.array(params, dtype=np.float32)

    def _folded_basis_matrix(self, uvz, basis_params):
        uvz = np.asarray(uvz, dtype=np.float64)
        basis_params = np.asarray(basis_params, dtype=np.float64)
        u = uvz[:, 0][:, None]
        v = uvz[:, 1][:, None]
        z = uvz[:, 2][:, None]
        ku = basis_params[None, :, 0]
        kv = basis_params[None, :, 1]
        az = basis_params[None, :, 2]
        z0 = basis_params[None, :, 3]
        bx = np.cos((2.0 * np.pi) * ku * u)
        by = np.cos((2.0 * np.pi) * kv * v)
        bz = np.exp(-az * np.maximum(0.0, z - z0))
        return (bx * by * bz).astype(np.float64)

    def _folded_fractional_coords(self, xyz):
        if self.surface_lvec is None:
            raise ValueError('_folded_fractional_coords(): surface_lvec is not initialized')
        lvec2d = getattr(self, 'folded_lvec_basis', None)
        if lvec2d is None:
            a = np.array(self.surface_lvec[0, :3], dtype=np.float64)
            b = np.array(self.surface_lvec[1, :3], dtype=np.float64)
        else:
            a = np.array(lvec2d[0, :3], dtype=np.float64)
            b = np.array(lvec2d[1, :3], dtype=np.float64)
        M = np.array([[a[0], b[0]], [a[1], b[1]]], dtype=np.float64)
        det = float(np.linalg.det(M))
        if abs(det) < 1e-12:
            raise ValueError(f'_folded_fractional_coords(): degenerate in-plane lattice det={det}')
        invM = np.linalg.inv(M)
        xy = np.asarray(xyz[:, :2], dtype=np.float64)
        uv = (invM @ xy.T).T
        uv -= np.floor(uv)
        out = np.empty((len(xyz), 3), dtype=np.float64)
        out[:, :2] = uv
        out[:, 2] = np.asarray(xyz[:, 2], dtype=np.float64)
        return out

    def _infer_folded_primitive_lvec2d(self):
        if self.surface_lvec is None:
            raise ValueError('_infer_folded_primitive_lvec2d(): surface_lvec is not initialized')
        a = np.array(self.surface_lvec[0, :3], dtype=np.float64)
        b = np.array(self.surface_lvec[1, :3], dtype=np.float64)
        apos = getattr(self, 'surface_atoms', None)
        enames = getattr(self, 'surface_enames', None)
        if apos is None or enames is None:
            out = np.zeros((3, 4), dtype=np.float32)
            out[0, :3] = a.astype(np.float32)
            out[1, :3] = b.astype(np.float32)
            out[2, :3] = np.array(self.surface_lvec[2, :3], dtype=np.float32)
            return out
        apos = np.asarray(apos[:, :3], dtype=np.float64)
        enames = np.asarray(enames)
        la = np.linalg.norm(a[:2])
        lb = np.linalg.norm(b[:2])
        ah = a[:2] / max(la, 1e-12)
        bh = b[:2] / max(lb, 1e-12)
        ua = apos[:, :2] @ ah
        ub = apos[:, :2] @ bh

        def infer_period(ucoord, vcoord, labels, axis_len):
            best = None
            ztol = 1e-3
            vtol = 1e-3
            for lab in np.unique(labels):
                mask_lab = labels == lab
                if np.count_nonzero(mask_lab) < 2:
                    continue
                zs = apos[mask_lab, 2]
                us = ucoord[mask_lab]
                vs = vcoord[mask_lab]
                for z0 in np.unique(np.round(zs, 3)):
                    mask_z = np.abs(zs - z0) < ztol
                    if np.count_nonzero(mask_z) < 2:
                        continue
                    usz = us[mask_z]
                    vsz = vs[mask_z]
                    vkeys = np.unique(np.round(vsz, 3))
                    for vk in vkeys:
                        mask_v = np.abs(vsz - vk) < vtol
                        vals = np.unique(np.round(usz[mask_v], 6))
                        if len(vals) < 2:
                            continue
                        ds = np.diff(np.sort(vals))
                        ds = ds[ds > 1e-6]
                        if len(ds) <= 0:
                            continue
                        d = float(np.min(ds))
                        if (d > 1e-6) and (d < axis_len - 1e-6):
                            best = d if (best is None or d < best) else best
            return axis_len if best is None else best

        pa = infer_period(ua, ub, enames, la)
        pb = infer_period(ub, ua, enames, lb)
        out = np.zeros((3, 4), dtype=np.float32)
        out[0, :3] = (a * (pa / max(la, 1e-12))).astype(np.float32)
        out[1, :3] = (b * (pb / max(lb, 1e-12))).astype(np.float32)
        out[2, :3] = np.array(self.surface_lvec[2, :3], dtype=np.float32)
        return out

    def fit_folded_surface_basis(self, surf_xyz=None, type_map=None, nPBC=(4,4,0), z_range=(0.5, 8.0), nu=4, nv=4, nz=4, nxy=32, nz_samp=40, r_damp=0.0, alpha_morse=1.8, bMacro=True, components=('total',), fit_mask=None, weight_power=0.0, coulomb_solver='morse', ewald_n_harm=6):
        if self.rigid_REQs0 is None:
            raise ValueError('fit_folded_surface_basis(): call init_rigid_molecule_batch() first')
        if surf_xyz is None:
            surf_xyz = getattr(self, 'surface_source_xyz', None)
        if surf_xyz is None:
            raise ValueError('fit_folded_surface_basis(): surf_xyz must be provided the first time')
        uniq_REQs, atom_type_ids = self._folded_unique_types(self.rigid_REQs0)
        ntypes = len(uniq_REQs)
        if ntypes > FOLDED_TYPES_MAX:
            raise ValueError(f'fit_folded_surface_basis(): ntypes={ntypes} exceeds FOLDED_TYPES_MAX={FOLDED_TYPES_MAX}')
        basis_params = self._build_folded_basis_params(nu=nu, nv=nv, nz=nz, z0=float(z_range[0]), z_scale=1.0/max(1e-6, float(z_range[1]-z_range[0])))
        nbasis = len(basis_params)
        self.folded_lvec_basis = self._infer_folded_primitive_lvec2d()
        a = np.array(self.folded_lvec_basis[0, :3], dtype=np.float32)
        b = np.array(self.folded_lvec_basis[1, :3], dtype=np.float32)
        u_freqs = np.unique(basis_params[:, 0]).astype(int).tolist()
        v_freqs = np.unique(basis_params[:, 1]).astype(int).tolist()
        alphas = np.unique(np.round(basis_params[:, 2].astype(np.float64), 6)).tolist()
        print(f'[folded] basis primitive_lvec a={a.tolist()} b={b.tolist()} |a|={float(np.linalg.norm(a[:2])):.6f} |b|={float(np.linalg.norm(b[:2])):.6f}')
        print(f'[folded] basis u_freqs={u_freqs} v_freqs={v_freqs} z_alphas={alphas} z0={float(z_range[0]):.6f} z1={float(z_range[1]):.6f}')
        us = np.linspace(0.0, 1.0, int(nxy), endpoint=False, dtype=np.float32)
        vs = np.linspace(0.0, 1.0, int(nxy), endpoint=False, dtype=np.float32)
        zs = np.linspace(float(z_range[0]), float(z_range[1]), int(nz_samp), endpoint=True, dtype=np.float32)
        xyz = []
        uvz = []
        for z in zs:
            for v in vs:
                for u in us:
                    p = a*u + b*v
                    xyz.append((float(p[0]), float(p[1]), float(z)))
                    uvz.append((float(u), float(v), float(z)))
        xyz = np.array(xyz, dtype=np.float32)
        uvz = np.array(uvz, dtype=np.float64)
        transforms = np.zeros((len(xyz), 3, 4), dtype=np.float32)
        transforms[:, 0, 0] = 1.0; transforms[:, 1, 1] = 1.0; transforms[:, 2, 2] = 1.0
        transforms[:, :, 3] = xyz
        Phi = self._folded_basis_matrix(uvz, basis_params)
        comp_keys = [str(c).lower() for c in components]
        coeff_sets = {ck: np.zeros((ntypes, nbasis), dtype=np.float32) for ck in comp_keys}
        refs = {ck: [] for ck in comp_keys}
        weights = []
        coulomb_solver = str(coulomb_solver).lower()
        if (coulomb_solver not in ('morse', 'ewald2d', 'none')):
            raise ValueError(f"fit_folded_surface_basis(): unknown coulomb_solver='{coulomb_solver}', expected 'morse','ewald2d','none'")
        want_coulomb = ('coulomb' in comp_keys)
        ew = None
        z_top = None
        if want_coulomb and coulomb_solver == 'ewald2d':
            from .SurfaceEwald import SurfaceEwaldCL
            if surf_xyz is None:
                raise ValueError('fit_folded_surface_basis(): coulomb_solver=ewald2d requires surf_xyz')
            if getattr(self, 'surface_atoms', None) is None:
                raise ValueError('fit_folded_surface_basis(): call set_surface() before coulomb_solver=ewald2d')
            if getattr(self, 'surface_REQs', None) is None:
                raise ValueError('fit_folded_surface_basis(): surface_REQs not available (call set_surface())')
            z_top = float(np.max(self.surface_atoms[:, 2]))
            ion_data = np.zeros((len(self.surface_atoms), 4), dtype=np.float32)
            ion_data[:, :3] = self.surface_atoms[:, :3]
            ion_data[:, 2] -= z_top
            ion_data[:, 3] = np.asarray(self.surface_REQs[:, 2], dtype=np.float32)
            a_vec = np.array(self.surface_lvec[0, :2], dtype=np.float32)
            b_vec = np.array(self.surface_lvec[1, :2], dtype=np.float32)
            ew = SurfaceEwaldCL(platform='nvidia')
            ew.prepare_system(ion_data, a_vec, b_vec, n_harm=int(ewald_n_harm))
            print(f"[folded] coulomb_solver=ewald2d z_top={z_top:.6f} ewald_n_harm={int(ewald_n_harm)}")
        z_mask = np.ones(len(uvz), dtype=bool) if fit_mask is None else np.asarray(fit_mask, dtype=bool).reshape(-1)
        if len(z_mask) != len(uvz):
            raise ValueError(f'fit_folded_surface_basis(): fit_mask length {len(z_mask)} != nsamples {len(uvz)}')
        for it in range(ntypes):
            md1 = MolecularDynamics(nloc=self.nloc, debug_build_options='-DDBG_UFF=0')
            md1.init_rigid_molecule_batch(np.zeros((1,3), dtype=np.float32), uniq_REQs[it:it+1], nSystems=min(max(len(transforms), 1), 8192))
            md1.set_surface(surf_xyz, nPBC=nPBC, pos0=(0.0,0.0,0.0), alpha_morse=alpha_morse, r_damp=r_damp, bMacro=bMacro, type_map=type_map)
            morse_comps = [ck for ck in comp_keys if ck != 'coulomb'] if (want_coulomb and coulomb_solver == 'ewald2d') else comp_keys
            out = md1.eval_rigid_getSurfMorse_components(transforms.reshape(-1,12), chunk_size=md1.nSystems, components=morse_comps)
            for ck in comp_keys:
                if ck == 'coulomb' and want_coulomb and coulomb_solver == 'ewald2d':
                    if ew is None:
                        raise ValueError('fit_folded_surface_basis(): internal error ew is None for coulomb_solver=ewald2d')
                    if z_top is None:
                        raise ValueError('fit_folded_surface_basis(): internal error z_top is None for coulomb_solver=ewald2d')
                    q_probe = float(uniq_REQs[it, 2])
                    X = xyz[:, 0].reshape(-1, 1)
                    Y = xyz[:, 1].reshape(-1, 1)
                    Z = (xyz[:, 2] - z_top).reshape(-1, 1)
                    phi = ew.eval_full(X, Y, Z)
                    y = np.asarray(phi.reshape(-1) * q_probe, dtype=np.float64)
                elif ck == 'coulomb' and want_coulomb and coulomb_solver == 'none':
                    y = np.zeros(len(transforms), dtype=np.float64)
                else:
                    y = np.asarray(out[ck], dtype=np.float64)
                refs[ck].append(y)
                ww = np.ones(len(y), dtype=np.float64)
                ww[~z_mask] = 0.0
                if weight_power != 0.0:
                    ymin = float(np.min(y[z_mask])) if np.any(z_mask) else float(np.min(y))
                    ww *= np.exp(np.clip(weight_power*(ymin - y), -60.0, 60.0))
                Phiw = Phi * ww[:, None]
                yw = y * ww
                S, *_ = np.linalg.lstsq(Phiw, yw, rcond=None)
                coeff_sets[ck][it, :nbasis] = S.astype(np.float32)
                if ck == comp_keys[0]:
                    weights.append(ww.copy())
        self.surface_source_xyz = surf_xyz
        coeffs = coeff_sets['total'] if 'total' in coeff_sets else coeff_sets[comp_keys[0]]
        self.folded_params = {
            'basis_params': basis_params,
            'coeffs': coeffs,
            'coeff_sets': coeff_sets,
            'atom_type_ids': atom_type_ids.astype(np.int32),
            'unique_REQs': uniq_REQs,
            'basis_lvec2d': self.folded_lvec_basis.copy(),
            'nu': int(nu), 'nv': int(nv), 'nz': int(nz),
            'z_range': (float(z_range[0]), float(z_range[1])),
            'nxy': int(nxy), 'nz_samp': int(nz_samp),
        }
        self.folded_fit_info = {
            'uvz': uvz,
            'Phi': Phi,
            'refs': refs,
            'weights': weights,
            'fit_mask': z_mask,
            'weight_power': float(weight_power),
        }
        self.toGPU('folded_atom_type', atom_type_ids.astype(np.int32))
        self.folded_type_ids = atom_type_ids.astype(np.int32)
        self.kernel_params['folded_lvec2d'] = np.array([a[0], b[0], a[1], b[1]], dtype=np.float32)
        self.kernel_params['folded_z0'] = np.array([float(z_range[0]), 0.0, 0.0, 0.0], dtype=np.float32)
        self._set_folded_coefficients(coeffs)
        if 'getSurfFolded' in self.kernelheaders:
            self.kernel_args_getSurfFolded = self.generate_kernel_args('getSurfFolded', bPrint=False)
        if 'getSurfFolded_harmonics' in self.kernelheaders:
            self.kernel_args_getSurfFolded_harmonics = self.generate_kernel_args('getSurfFolded_harmonics', bPrint=False)
        if 'getSurfFolded_workgroup' in self.kernelheaders:
            self.kernel_args_getSurfFolded_workgroup = self.generate_kernel_args('getSurfFolded_workgroup', bPrint=False)
        return self.folded_params

    def eval_rigid_getSurfFolded(self, transforms, chunk_size=None, component='total'):
        if self.folded_params is None:
            raise ValueError('eval_rigid_getSurfFolded(): call fit_folded_surface_basis() first')
        key = str(component).lower()
        coeff_sets = self.folded_params.get('coeff_sets', None) if self.folded_params is not None else None
        if coeff_sets is not None:
            if key not in coeff_sets:
                raise ValueError(f"eval_rigid_getSurfFolded(): component '{component}' not available, have {sorted(coeff_sets.keys())}")
            self._set_folded_coefficients(coeff_sets[key])
        T = self.wrap_rigid_transforms_PBC(transforms).reshape(-1, 3, 4)
        nconf = len(T)
        if nconf <= 0:
            z = np.zeros(0, dtype=np.float32)
            return {'total': z.copy(), 'LJ': z.copy(), 'Coulomb': z.copy(), 'HBond': z.copy()}
        if chunk_size is None:
            chunk_size = self.nSystems
        chunk_size = int(min(chunk_size, self.nSystems))
        if chunk_size <= 0:
            raise ValueError(f'eval_rigid_getSurfFolded(): invalid chunk_size={chunk_size} for nSystems={self.nSystems}')
        out_total = np.empty(nconf, dtype=np.float32)
        float_size = np.float32().itemsize
        sys_bytes = self.nvecs * 4 * float_size
        t_prep = 0.0
        t_kernel = 0.0
        t_download = 0.0
        for i0 in range(0, nconf, chunk_size):
            nch = min(chunk_size, nconf - i0)
            t1 = time.perf_counter()
            self.upload_rigid_transforms(T[i0:i0+nch], iSys0=0)
            cl.enqueue_fill_buffer(self.queue, self.buffer_dict['aforce'], np.zeros(1, dtype=np.float32), 0, nch * sys_bytes)
            self.queue.finish()
            t2 = time.perf_counter()
            self.run_getSurfFolded(nSystems=nch)
            self.queue.finish()
            t3 = time.perf_counter()
            aforce = np.empty((nch, self.nvecs, 4), dtype=np.float32)
            self.fromGPU('aforce', aforce)
            self.queue.finish()
            out_total[i0:i0+nch] = -aforce[:, :self.natoms, 3].sum(axis=1)
            t4 = time.perf_counter()
            t_prep += (t2 - t1)
            t_kernel += (t3 - t2)
            t_download += (t4 - t3)
        z = np.zeros_like(out_total)
        return {'total': out_total, 'LJ': z.copy(), 'Coulomb': z.copy(), 'HBond': z.copy(), 't_prep_s': t_prep, 't_kernel_s': t_kernel, 't_download_s': t_download, 't_total_s': t_prep + t_kernel + t_download}

    def eval_gridff_points(self, points, PLQH=None, chunk_size=None):
        if self.kernel_args_sampleGridFF_Bspline_points is None:
            if 'sampleGridFF_Bspline_points' in self.kernelheaders:
                self.kernel_args_sampleGridFF_Bspline_points = self.generate_kernel_args('sampleGridFF_Bspline_points', bPrint=False)
            else:
                raise ValueError('eval_gridff_points(): sampleGridFF_Bspline_points kernel not available')
        pts = np.asarray(points, dtype=np.float32).reshape(-1, 3)
        nconf = len(pts)
        if nconf <= 0:
            z = np.zeros(0, dtype=np.float32)
            return {'total': z.copy(), 'fx': z.copy(), 'fy': z.copy(), 'fz': z.copy()}
        if chunk_size is None:
            chunk_size = self.nSystems
        chunk_size = int(min(chunk_size, self.nSystems))
        if chunk_size <= 0:
            raise ValueError(f'eval_gridff_points(): invalid chunk_size={chunk_size} for nSystems={self.nSystems}')
        out_total = np.empty(nconf, dtype=np.float32)
        out_fx = np.empty(nconf, dtype=np.float32)
        out_fy = np.empty(nconf, dtype=np.float32)
        out_fz = np.empty(nconf, dtype=np.float32)
        float_size = np.float32().itemsize
        sys_bytes = self.nvecs * 4 * float_size
        plqh = self._resolve_plqh(PLQH=PLQH)
        for i0 in range(0, nconf, chunk_size):
            nch = min(chunk_size, nconf - i0)
            apos = np.zeros((nch, self.nvecs, 4), dtype=np.float32)
            apos[:, 0, :3] = pts[i0:i0+nch, :3]
            self.toGPU('apos', apos, byte_offset=0)
            cl.enqueue_fill_buffer(self.queue, self.buffer_dict['aforce'], np.zeros(1, dtype=np.float32), 0, nch * sys_bytes)
            self.queue.finish()
            self.run_sampleGridFF_Bspline_points(nSystems=nch, PLQH=plqh)
            aforce = np.empty((nch, self.nvecs, 4), dtype=np.float32)
            self.fromGPU('aforce', aforce)
            self.queue.finish()
            out_fx[i0:i0+nch] = aforce[:, 0, 0]
            out_fy[i0:i0+nch] = aforce[:, 0, 1]
            out_fz[i0:i0+nch] = aforce[:, 0, 2]
            out_total[i0:i0+nch] = -aforce[:, 0, 3]
        return {'total': out_total, 'fx': out_fx, 'fy': out_fy, 'fz': out_fz}

    def eval_surface_iso_gridff(self, probe_req, sel_PLQH, col_PLQH, x_range, y_range, z_range, nx=64, ny=64, nz=80, mode='threshold', threshold=0.0, z_top=0.0):
        if 'getSurfaceIsoGridFF' not in self.kernelheaders:
            raise ValueError("eval_surface_iso_gridff(): kernel 'getSurfaceIsoGridFF' not available")
        nx = int(nx); ny = int(ny); nz = int(nz)
        if nx <= 0 or ny <= 0 or nz < 2:
            raise ValueError(f"eval_surface_iso_gridff(): invalid grid ({nx},{ny},{nz})")
        nxy = nx * ny
        float_size = np.float32().itemsize
        self.check_buf('surf_xyzq', nxy * 4 * float_size, cl.mem_flags.READ_WRITE)
        self.check_buf('surf_zc', nxy * 2 * float_size, cl.mem_flags.READ_WRITE)
        x0, x1 = float(x_range[0]), float(x_range[1])
        y0, y1 = float(y_range[0]), float(y_range[1])
        z0, z1 = float(z_range[0]), float(z_range[1])
        dx = 0.0 if nx <= 1 else (x1 - x0) / float(nx - 1)
        dy = 0.0 if ny <= 1 else (y1 - y0) / float(ny - 1)
        dz = 0.0 if nz <= 1 else (z1 - z0) / float(nz - 1)
        imode = 0 if str(mode).lower() == 'threshold' else 1
        kernel = self.prg.getSurfaceIsoGridFF
        kernel.set_args(
            np.asarray(self.kernel_params['grid_ns'], dtype=np.int32),
            self.buffer_dict['BsplinePLQ'],
            np.asarray(self.kernel_params['grid_invStep'], dtype=np.float32),
            np.asarray(self.kernel_params['grid_p0'], dtype=np.float32),
            np.asarray(sel_PLQH, dtype=np.float32).reshape(4),
            np.asarray(col_PLQH, dtype=np.float32).reshape(4),
            np.array([nx, ny, nz, imode], dtype=np.int32),
            np.array([x0, y0, z0, float(threshold)], dtype=np.float32),
            np.array([dx, dy, dz, z1], dtype=np.float32),
            np.array([float(z_top), 0.0, 0.0, 0.0], dtype=np.float32),
            self.buffer_dict['surf_xyzq'],
            self.buffer_dict['surf_zc'],
        )
        local = (min(self.nloc, 16), min(self.nloc, 16))
        global_size = (clu.roundup_global_size(nx, local[0]), clu.roundup_global_size(ny, local[1]))
        cl.enqueue_nd_range_kernel(self.queue, kernel, global_size, local)
        self.queue.finish()
        xyzq = np.empty((nxy, 4), dtype=np.float32)
        zc = np.empty((nxy, 2), dtype=np.float32)
        self.fromGPU('surf_xyzq', xyzq)
        self.fromGPU('surf_zc', zc)
        self.queue.finish()
        xyzq = xyzq.reshape(ny, nx, 4).transpose(1, 0, 2).copy()
        zc = zc.reshape(ny, nx, 2).transpose(1, 0, 2).copy()
        return {'points_world': xyzq[:, :, :3], 'ok_mask': xyzq[:, :, 3] > 0.5, 'z_report': zc[:, :, 0], 'color': zc[:, :, 1]}

    def eval_surface_iso_morse(self, probe_req, sel_PLQH, col_PLQH, x_range, y_range, z_range, nx=64, ny=64, nz=80, mode='threshold', threshold=0.0, z_top=None):
        if 'getSurfaceIsoSurfMorse' not in self.kernelheaders:
            raise ValueError("eval_surface_iso_morse(): kernel 'getSurfaceIsoSurfMorse' not available")
        if getattr(self, 'surface_atoms', None) is None:
            raise ValueError('eval_surface_iso_morse(): call set_surface() first')
        nx = int(nx); ny = int(ny); nz = int(nz)
        if nx <= 0 or ny <= 0 or nz < 2:
            raise ValueError(f"eval_surface_iso_morse(): invalid grid ({nx},{ny},{nz})")
        nxy = nx * ny
        float_size = np.float32().itemsize
        self.check_buf('surf_xyzq', nxy * 4 * float_size, cl.mem_flags.READ_WRITE)
        self.check_buf('surf_zc', nxy * 2 * float_size, cl.mem_flags.READ_WRITE)
        x0, x1 = float(x_range[0]), float(x_range[1])
        y0, y1 = float(y_range[0]), float(y_range[1])
        z0, z1 = float(z_range[0]), float(z_range[1])
        dx = 0.0 if nx <= 1 else (x1 - x0) / float(nx - 1)
        dy = 0.0 if ny <= 1 else (y1 - y0) / float(ny - 1)
        dz = 0.0 if nz <= 1 else (z1 - z0) / float(nz - 1)
        imode = 0 if str(mode).lower() == 'threshold' else 1
        z_top = float(np.max(self.surface_atoms[:, 2])) if z_top is None else float(z_top)
        ns = np.array([1, 0, len(self.surface_atoms), 0], dtype=np.int32)
        kernel = self.prg.getSurfaceIsoSurfMorse
        kernel.set_args(
            ns,
            self.buffer_dict['atoms_s'],
            self.buffer_dict['REQ_s'],
            self.buffer_dict['surf_mpos'],
            self.buffer_dict['surf_mdip'],
            self.buffer_dict['surf_mQa'],
            self.buffer_dict['surf_mQb'],
            self.buffer_dict['surf_mQc'],
            self.buffer_dict['surf_qQa'],
            self.buffer_dict['surf_qQb'],
            self.buffer_dict['surf_qQc'],
            np.asarray(self.kernel_params['nPBC'], dtype=np.int32),
            np.asarray(self.kernel_params['lvec'], dtype=np.float32),
            np.asarray(self.kernel_params['GFFParams'], dtype=np.float32),
            np.asarray(probe_req, dtype=np.float32).reshape(4),
            np.asarray(sel_PLQH, dtype=np.float32).reshape(4),
            np.asarray(col_PLQH, dtype=np.float32).reshape(4),
            np.array([nx, ny, nz, imode], dtype=np.int32),
            np.array([x0, y0, z_top + z0, float(threshold)], dtype=np.float32),
            np.array([dx, dy, dz, z_top + z1], dtype=np.float32),
            self.buffer_dict['surf_xyzq'],
            self.buffer_dict['surf_zc'],
        )
        local = (min(self.nloc, 16), min(self.nloc, 16))
        global_size = (clu.roundup_global_size(nx, local[0]), clu.roundup_global_size(ny, local[1]))
        cl.enqueue_nd_range_kernel(self.queue, kernel, global_size, local)
        self.queue.finish()
        xyzq = np.empty((nxy, 4), dtype=np.float32)
        zc = np.empty((nxy, 2), dtype=np.float32)
        self.fromGPU('surf_xyzq', xyzq)
        self.fromGPU('surf_zc', zc)
        self.queue.finish()
        xyzq = xyzq.reshape(ny, nx, 4).transpose(1, 0, 2).copy()
        zc = zc.reshape(ny, nx, 2).transpose(1, 0, 2).copy()
        zc[:, :, 0] = xyzq[:, :, 2] - z_top
        return {'points_world': xyzq[:, :, :3], 'ok_mask': xyzq[:, :, 3] > 0.5, 'z_report': zc[:, :, 0], 'color': zc[:, :, 1]}

    def run_runMD(self):
        self.prg.runMD(self.queue, self.sz_nvec, self.sz_loc, *self.kernel_args_runMD)
        self.queue.finish()

    def run_sampleGrid_tex(self, apos=None, bUseTexture=False):
        if apos is not None:
            self.toGPU('apos', apos)
        if bUseTexture:
            self.prg.sampleGrid_tex(self.queue, self.sz_na, self.sz_loc, *self.kernel_args_sampleGrid_tex)
        else:
            self.prg.sampleGrid(self.queue, self.sz_na, self.sz_loc, *self.kernel_args_sampleGrid)
        self.fromGPU('aforce', self.aforce)
        self.queue.finish()
        return self.aforce.reshape(-1, 4).copy()

    def run_MD_py(self, nsteps, use_gridff=False):
        """Run molecular dynamics simulation..."""
        for i in range(nsteps):
            if use_gridff and self.has_gridff:
                if self.use_texture:
                    self.prg.getNonBond_GridFF_Bspline_tex(self.queue, self.sz_na, self.sz_loc, *self.kernel_args_getNonBond_GridFF_Bspline_tex)
                else:
                    self.prg.getNonBond_GridFF_Bspline(self.queue, self.sz_na, self.sz_loc, *self.kernel_args_getNonBond_GridFF_Bspline)
            else:
                self.prg.getNonBond       (self.queue, self.sz_na, self.sz_loc,*self.kernel_args_getNonBond)
            self.prg.getMMFFf4        (self.queue, self.sz_node, self.sz_loc, *self.kernel_args_getMMFFf4)
            self.prg.updateAtomsMMFFf4(self.queue, self.sz_nvec, self.sz_loc, *self.kernel_args_updateAtomsMMFFf4)
        self.fromGPU('apos',   self.atoms)
        self.fromGPU('aforce', self.aforce)
        self.queue.finish()
        return self.atoms.reshape(-1, 4), self.aforce.reshape(-1, 4)
    
    def run_step_basic(self, do_nb=False ):
        """Run a single MD step using basic (non-rotational) force kernels."""
        if do_nb: 
            self.prg.getNonBond      (self.queue, self.sz_na,   self.sz_loc, *self.kernel_args_getNonBond)
        else:
            self.prg.cleanForceMMFFf4(self.queue, self.sz_na,   self.sz_loc, *self.kernel_args_cleanForceMMFFf4)
        self.prg.getMMFFf4           (self.queue, self.sz_node, self.sz_loc, *self.kernel_args_getMMFFf4)
        self.prg.updateAtomsMMFFf4   (self.queue, self.sz_nvec, self.sz_loc, *self.kernel_args_updateAtomsMMFFf4)
        self.queue.finish()

    def run_step_rot(self, do_nb=False):
        """Run a single MD step using rotational force and update kernels."""
        if do_nb: 
            self.prg.getNonBond       (self.queue, self.sz_na,   self.sz_loc, *self.kernel_args_getNonBond)
        else:
            self.prg.cleanForceMMFFf4 (self.queue, self.sz_na,   self.sz_loc, *self.kernel_args_cleanForceMMFFf4)
        self.prg.getMMFFf4_rot        (self.queue, self.sz_node, self.sz_loc, *self.kernel_args_getMMFFf4_rot)
        self.prg.updateAtomsMMFFf4_rot(self.queue, self.sz_nvec, self.sz_loc, *self.kernel_args_updateAtomsMMFFf4_rot)
        self.queue.finish()

    def run_MD_step(self, do_clean=True, do_nb=False, do_mmff=True, use_rot=False, force_kernel='basic'):
        """Backward-compatible wrapper for previous API (deprecated)."""
        if use_rot:
            self.run_step_rot(do_clean=do_clean, do_nb=do_nb, do_mmff=do_mmff, integrator=force_kernel)
        else:
            self.run_step_basic(do_clean=do_clean, do_nb=do_nb, do_mmff=do_mmff, integrator=force_kernel)

    def download_results(self):
        self.fromGPU('apos',   self.atoms)
        self.fromGPU('aforce', self.aforce)
        return self.atoms.reshape(-1, 4), self.aforce.reshape(-1, 4)
    
    def initGridFF(self, grid_shape, bspline_data, grid_p0, grid_step, use_texture=False, r_damp=0.0, alpha_morse=0.0, bKernels=True):
        """Initialize GridFF with B-spline data"""
        
        #grid_shape = grid_shape[::-1] #.copy()

        print("MolecularDynamics::initGridFF() grid_shape: ", grid_shape)
        self.has_gridff = True
        self.use_texture = use_texture
        
        # 1. Ensure kernel_params exists
        if not hasattr(self, 'kernel_params'):
            self.init_kernel_params()
            
        # 2. Set grid parameters
        self.kernel_params.update({
            'grid_ns':      np.array([*grid_shape ,0],                   dtype=np.int32),
            'grid_invStep': np.array([1.0/s for s in grid_step] + [0.0], dtype=np.float32),
            'grid_p0':      np.array([*grid_p0   ,0.0],                  dtype=np.float32),
            'GFFParams':    np.array([r_damp, alpha_morse, 0.0, 0.0],    dtype=np.float32),
            'nstep':        np.int32(self.nstep),
        })
        
        # 3. Create buffers BEFORE generating kernel args
        if use_texture:
            print(f"MolecularDynamics::initGridFF() use_texture=True grid_shape={grid_shape} bspline_data.shape={bspline_data.shape} bspline_data.dtype={bspline_data.dtype}")
            print("Original bspline_data dimensions:", bspline_data.ndim)

            #bspline_data = bspline_data.transpose(2,1,0,3).copy(); grid_shape = (40,40,200)   # better 
            #bspline_data = bspline_data.transpose(2,1,0,3).copy(); grid_shape = ( grid_shape[2], grid_shape[1], grid_shape[0] )   # better 
            bspline_data = bspline_data.transpose(2,1,0,3).copy(); grid_shape = ( grid_shape[0], grid_shape[1], grid_shape[2] )   # better 
            print("!!!!! grid_shape: ", grid_shape)

            fmt = cl.ImageFormat(cl.channel_order.RGBA, cl.channel_type.FLOAT)
            tex = cl.Image(self.ctx, cl.mem_flags.READ_ONLY, fmt, shape=tuple(grid_shape))
            
            #tex = cl.Image(self.ctx, cl.mem_flags.READ_ONLY, fmt, shape=tuple(grid_shape[::-1]))

            cl.enqueue_copy(self.queue, tex, bspline_data,  origin=(0, 0, 0), region=grid_shape)
            #cl.enqueue_copy(self.queue, tex, bspline_data,  origin=(0, 0, 0), region=bspline_data.shape[:3])
            self.buffer_dict['BsplinePLQH_tex'] = tex
            if bKernels:
                self.kernel_args_getNonBond_GridFF_Bspline_tex = self.generate_kernel_args("getNonBond_GridFF_Bspline_tex", bPrint=False)            
        else:
            print("MolecularDynamics::initGridFF() use_texture=False")

            #bspline_data = bspline_data.transpose(2,1,0,3).copy(); # grid_shape = ( grid_shape[2], grid_shape[1], grid_shape[0] )   # better 

            print("!!!!! grid_shape: ", grid_shape)

            buf = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY|cl.mem_flags.COPY_HOST_PTR, hostbuf=bspline_data)
            self.buffer_dict['BsplinePLQH'] = buf
            # Kernel expects argument name 'BsplinePLQ'
            self.buffer_dict['BsplinePLQ'] = buf
            if bKernels:
                self.kernel_args_getNonBond_GridFF_Bspline = None
                if "getNonBond_GridFF_Bspline" in self.kernelheaders:
                    try:
                        self.kernel_args_getNonBond_GridFF_Bspline = self.generate_kernel_args("getNonBond_GridFF_Bspline", bPrint=False)
                    except Exception as e:
                        print(f"initGridFF: skipping getNonBond_GridFF_Bspline ({e})")
                if "sampleGridFF_Bspline_points" in self.kernelheaders:
                    try:
                        self.kernel_args_sampleGridFF_Bspline_points = self.generate_kernel_args("sampleGridFF_Bspline_points", bPrint=False)
                    except Exception as e:
                        print(f"initGridFF: skipping sampleGridFF_Bspline_points ({e})")
                self.kernel_args_getNonBond_GridFF_Bspline_ex2 = None
                if self.enable_nonbond and ('excl' in self.buffer_dict) and ("getNonBond_GridFF_Bspline_ex2" in self.kernelheaders):
                    try:
                        self.kernel_args_getNonBond_GridFF_Bspline_ex2 = self.generate_kernel_args("getNonBond_GridFF_Bspline_ex2", bPrint=False)
                    except Exception as e:
                        print(f"initGridFF: skipping getNonBond_GridFF_Bspline_ex2 ({e})")
        print("MolecularDynamics::initGridFF() DONE")

    def scan_1D( self, nsteps=100, d=[0.1,0.0,0.0], p0=[0.0,0.0,0.0], use_texture=False, mmff=None ):
        if mmff is None: mmff = self.mmff_list[0]
        # Initialize position and force grid_dataays
        pos    = np.zeros((nsteps,4), dtype=np.float32)
        forces = np.zeros((nsteps,4), dtype=np.float32)
        print("Running 1D force scan...")
        d  = np.array(d, dtype=np.float32)
        p0 = np.array(p0,  dtype=np.float32)
        for i in range(nsteps):
            # Update atom position along x-axis
            mmff.apos[0,:3] = p0 + d*i
            self.upload_all_systems()
            if use_texture:
                self.run_getNonBond_GridFF_Bspline_tex()
            else:
                self.run_getNonBond_GridFF_Bspline()
            pos_i, force_i = self.download_results()
            forces[i] = force_i[0]  # Get force on first (and only) atom
            pos[i] = pos_i[0]  # Store position for reference
            #if i % 10 == 0:  print(f"Step {i+1}/{nsteps}: x = {x[i]:.2f} Å, F = ({forces[i,0]:.3f}, {forces[i,1]:.3f}, {forces[i,2]:.3f}) kJ/mol/Å")
        return pos, forces


    def update_pbc_shifts(self, lvec, nPBC):
        """
        Generate PBC shifts for neighbor cells.
        lvec: (3,3) matrix (rows are vectors)
        nPBC: (3,) tuple or array
        """
        shifts = []
        nx, ny, nz = nPBC
        # Loop order must match C++: 
        # for(int iz=-nPBC.z; iz<=nPBC.z; iz++){ ... }
        # Note: relax_multi_mini.cl uses pbc_shifts[ipbc0+ic]
        # In C++, the neighbors are mapped to cell indices.
        # We need to ensure the order matches what neighborCell expects?
        # Actually, pbc_shifts in C++ `MolWorld_sp3_multi` seems to be just a list of all shift vectors?
        # No, let's look at `MolecularDynamics.py` usage or existing code.
        # There is no existing pbc shift gen in python.
        # In C++ `MolWorld_sp3_multi.h`, `update_pbc_shifts` generates `(2*nx+1)*(2*ny+1)*(2*nz+1)` shifts.
        # The order is z, y, x loops from -n to +n.
        
        n_shifts = (2*nx+1)*(2*ny+1)*(2*nz+1)
        shifts = np.zeros((n_shifts, 4), dtype=np.float32)
        
        idx = 0
        for iz in range(-nz, nz+1):
            for iy in range(-ny, ny+1):
                for ix in range(-nx, nx+1):
                    # shift = ix*a + iy*b + iz*c
                    sh = ix*lvec[0] + iy*lvec[1] + iz*lvec[2]
                    shifts[idx, :3] = sh
                    shifts[idx, 3]  = 0.0 # padding
                    idx += 1
        return shifts

    def pack_system(self, iSys, mmff):
        """Packs data from an MMFF instance into GPU buffers for a specific system index."""
        #print("pack_system() iSys=%d" % iSys)
        nvecs       = mmff.nvecs
        natoms      = mmff.natoms
        nnode       = mmff.nnode
        float4_size = 4 * np.float32().itemsize
        int4_size   = 4 * np.int32().itemsize
        mat3_size   = 12 * np.float32().itemsize # 3*float4

        if self.bPrintPackSystem:
            self._print_pack_system_params(iSys, mmff)

        # Back-neighbor indices per vector (atoms + pi) for recoil force accumulation; default to -1 if not provided.
        # IMPORTANT: updateAtomsMMFFf4 uses these as direct indices into the global `fneigh` buffer.
        # `fneigh` is laid out per-system with stride (nnode*4*2) float4 entries.
        # Therefore bkNeighs MUST be offset per-system, otherwise replicas read recoil from other systems.
        bk = np.full((nvecs, 4), -1, dtype=np.int32)
        if hasattr(mmff, 'back_neighs') and (mmff.back_neighs is not None):
            ncopy = min(mmff.back_neighs.shape[0], nvecs)
            bk[:ncopy, :] = mmff.back_neighs[:ncopy, :].astype(np.int32)
            if nnode > 0:
                fneigh_stride = int(nnode) * 8  # nnode * (4 neighbors) * (2 channels sigma+pi)
                fneigh_off = int(iSys) * fneigh_stride
                mask = bk[:ncopy, :] >= 0
                bk[:ncopy, :][mask] += np.int32(fneigh_off)
        offset_bk     = iSys * nvecs * int4_size
        offset_atoms  = iSys * nvecs  * float4_size
        offset_REQs   = iSys * natoms * float4_size
        offset_neighs = iSys * natoms * int4_size
        offset_apars  = iSys * nnode  * float4_size
        offset_lvec   = iSys * mat3_size
        
        self.toGPU('apos',      self._flat32(mmff.apos),    byte_offset=offset_atoms)
        self.toGPU('aforce',    self._flat32(mmff.fapos),   byte_offset=offset_atoms)
        self.toGPU('aforce_old', self._flat32(mmff.fapos),  byte_offset=offset_atoms)
        self.toGPU('REQs',      self._flat32(mmff.REQs),    byte_offset=offset_REQs)
        self.toGPU('neighs',    self._int32 (mmff.neighs),    byte_offset=offset_neighs)
        self.toGPU('neighCell', self._int32 (mmff.neighCell), byte_offset=offset_neighs)
        self.toGPU('bkNeighs',  bk,                         byte_offset=offset_bk)

        # Non-bonded exclusions (packed EXCL_MAX=16 ints per atom)
        if self.enable_nonbond and (self.buffer_dict.get('excl') is not None):
            if getattr(mmff, 'excl', None) is None:
                raise ValueError("enable_nonbond=True but mmff.excl is None; build exclusions before pack_system()")
            offset_excl = iSys * natoms * 16 * np.int32().itemsize
            self.toGPU('excl', self._int32(mmff.excl), byte_offset=offset_excl)
        self.toGPU('apars',     self._flat32(mmff.apars),   byte_offset=offset_apars)
        self.toGPU('bLs',       self._flat32(mmff.bLs),     byte_offset=offset_apars)
        self.toGPU('bKs',       self._flat32(mmff.bKs),     byte_offset=offset_apars)
        self.toGPU('Ksp',       self._flat32(mmff.Ksp),     byte_offset=offset_apars)
        self.toGPU('Kpp',       self._flat32(mmff.Kpp),     byte_offset=offset_apars)
        
        # Upload lattice vectors if present
        if hasattr(mmff, 'lvec') and mmff.lvec is not None:
            lvec_padded = np.zeros((3,4), dtype=np.float32)
            lvec_padded[:,:3] = mmff.lvec
            self.toGPU('lvecs', lvec_padded.flatten(), byte_offset=offset_lvec)
            
            # Inverse lattice vectors
            try:
                inv_lvec = np.linalg.inv(mmff.lvec)
                inv_lvec_padded = np.zeros((3,4), dtype=np.float32)
                inv_lvec_padded[:,:3] = inv_lvec
                self.toGPU('ilvecs', inv_lvec_padded.flatten(), byte_offset=offset_lvec)
            except np.linalg.LinAlgError:
                pass
            
            # PBC Shifts
            if hasattr(mmff, 'nPBC') and mmff.nPBC is not None:
                shifts = self.update_pbc_shifts(mmff.lvec, mmff.nPBC)
                # Ensure we don't overflow allocated pbc_shifts buffer
                # The allocation size is nSystems * npbc * 4 * float_size
                # self.npbc should match len(shifts)
                offset_shifts = iSys * self.npbc * float4_size
                # Check size
                if shifts.size <= self.npbc * 4:
                    self.toGPU('pbc_shifts', shifts.flatten(), byte_offset=offset_shifts)
                else:
                    print(f"Warning: Calculated PBC shifts size {shifts.size} exceeds buffer size {self.npbc*4}")

        #print("pack_system() iSys=%d" % iSys, "offset_atoms=%d" % offset_atoms, "offset_REQs=%d" % offset_REQs, "offset_neighs=%d" % offset_neighs, "offset_apars=%d" % offset_apars)
        #print("pack_system() iSys=%d" % iSys, "atoms.nbytes=%d" % mmff.apos.nbytes, "REQs.nbytes=%d" % mmff.REQs.nbytes, "neighs.nbytes=%d" % mmff.neighs.nbytes, "apars.nbytes=%d" % mmff.apars.nbytes)
        # MDparams layout expected by kernels: (dt, damp, friction)
        # Note: kernels currently hardcode Flimit; third component is used as velocity multiplier
        self.toGPU('MDparams',  np.array([mmff.dt, mmff.damp, mmff.damp], dtype=np.float32), byte_offset=iSys*float4_size)




    def scan_2D( self, ns=(50,50), du=[0.1,0.0,0.0], dv=[0.0,0.1,0.0],  p0=[0.0,0.0,0.0], use_texture=False, mmff=None ):
        if mmff is None: mmff = self.mmff_list[0]
        pos    = np.zeros((ns[0],ns[1],4), dtype=np.float32)
        forces = np.zeros((ns[0],ns[1],4), dtype=np.float32)
        du = np.array(du, dtype=np.float32)
        dv = np.array(dv, dtype=np.float32)
        p0 = np.array(p0,  dtype=np.float32)
        print("Running 2D force scan...")
        for ix in range(ns[0]):
            for iy in range(ns[1]):
                p = p0 + du*ix + dv*iy
                mmff.apos[0,:3] = p
                self.upload_all_systems()
                if use_texture:
                    self.run_getNonBond_GridFF_Bspline_tex()
                else:
                    self.run_getNonBond_GridFF_Bspline()
                pos_i, force_i = self.download_results()
                forces[ix, iy,:] = force_i[0,: ]  # Get force on first (and only) atom
                pos   [ix, iy,:] = pos_i  [0,: ]  # Store position for reference
        return pos, forces

    def scanSurfMorse_2D(self, ns=(50,50), du=[0.1,0.0,0.0], dv=[0.0,0.1,0.0], p0=[0.0,0.0,0.0], mmff=None):
        if self.kernel_args_getSurfMorse is None:
            raise ValueError('scanSurfMorse_2D() requires set_surface() and getSurfMorse kernel setup first')
        if mmff is None:
            mmff = self.mmff_list[0]
        pos = np.zeros((ns[0], ns[1], 4), dtype=np.float32)
        forces = np.zeros((ns[0], ns[1], 4), dtype=np.float32)
        du = np.array(du, dtype=np.float32)
        dv = np.array(dv, dtype=np.float32)
        p0 = np.array(p0, dtype=np.float32)
        for ix in range(ns[0]):
            for iy in range(ns[1]):
                p = p0 + du*ix + dv*iy
                mmff.apos[0,:3] = p[:3]
                self.pack_system(0, mmff)
                self.toGPU('aforce', np.zeros((self.nvecs,4), dtype=np.float32))
                self.run_getSurfMorse()
                pos_i, force_i = self.download_results()
                forces[ix, iy,:] = force_i[0,:]
                pos[ix, iy,:] = pos_i[0,:]
        return pos, forces

    def upload_rigid_transforms(self, transforms, iSys0=0):
        if self.rigid_apos0 is None or self.rigid_REQs0 is None:
            raise ValueError('upload_rigid_transforms(): call init_rigid_molecule_batch() first')
        T = self.wrap_rigid_transforms_PBC(transforms).reshape(-1, 3, 4)
        nsys = len(T)
        if (iSys0 < 0) or ((iSys0 + nsys) > self.nSystems):
            raise ValueError(f'upload_rigid_transforms(): systems [{iSys0},{iSys0+nsys}) exceed allocated nSystems={self.nSystems}')
        xyz = np.einsum('aj,nij->nai', self.rigid_apos0[:, :3], T[:, :, :3], optimize=True)
        xyz += T[:, None, :, 3]
        apos = np.zeros((nsys, self.nvecs, 4), dtype=np.float32)
        apos[:, :self.natoms, :3] = xyz
        byte_offset = iSys0 * self.nvecs * 4 * np.float32().itemsize
        self.toGPU('apos', apos, byte_offset=byte_offset)
        return apos

    def wrap_rigid_transforms_PBC(self, transforms):
        if getattr(self, 'surface_lvec', None) is None:
            return np.asarray(transforms, dtype=np.float32).reshape(-1, 3, 4)
        a = np.array(self.surface_lvec[0, :3], dtype=np.float64)
        b = np.array(self.surface_lvec[1, :3], dtype=np.float64)
        M = np.array([[a[0], b[0]], [a[1], b[1]]], dtype=np.float64)
        det = float(np.linalg.det(M))
        if abs(det) < 1e-12:
            raise ValueError(f'wrap_rigid_transforms_PBC(): degenerate lattice vectors det={det} a={a} b={b}')
        invM = np.linalg.inv(M)
        T = np.array(transforms, copy=True, dtype=np.float32).reshape(-1, 3, 4)
        txy = T[:, 0:2, 3].astype(np.float64)
        frac = (invM @ txy.T).T
        frac -= np.round(frac)
        txy2 = (M @ frac.T).T
        T[:, 0, 3] = txy2[:, 0]
        T[:, 1, 3] = txy2[:, 1]
        return T

    def eval_rigid_getSurfMorse(self, transforms, chunk_size=None, component='total', PLQH=None):
        if self.kernel_args_getSurfMorse is None:
            raise ValueError('eval_rigid_getSurfMorse(): call set_surface() before evaluation')
        self.kernel_params['PLQH'] = self._resolve_plqh(component=component, PLQH=PLQH)
        self.kernel_args_getSurfMorse = self.generate_kernel_args('getSurfMorse', bPrint=False)
        T = self.wrap_rigid_transforms_PBC(transforms).reshape(-1, 3, 4)
        nconf = len(T)
        if nconf <= 0:
            z = np.zeros(0, dtype=np.float32)
            return {'total': z.copy(), 'LJ': z.copy(), 'Coulomb': z.copy(), 'HBond': z.copy()}
        if chunk_size is None:
            chunk_size = self.nSystems
        chunk_size = int(min(chunk_size, self.nSystems))
        if chunk_size <= 0:
            raise ValueError(f'eval_rigid_getSurfMorse(): invalid chunk_size={chunk_size} for nSystems={self.nSystems}')
        out_total = np.empty(nconf, dtype=np.float32)
        float_size = np.float32().itemsize
        sys_bytes = self.nvecs * 4 * float_size
        t_prep = 0.0
        t_kernel = 0.0
        t_download = 0.0
        for i0 in range(0, nconf, chunk_size):
            nch = min(chunk_size, nconf - i0)
            t1 = time.perf_counter()
            self.upload_rigid_transforms(T[i0:i0+nch], iSys0=0)
            cl.enqueue_fill_buffer(self.queue, self.buffer_dict['aforce'], np.zeros(1, dtype=np.float32), 0, nch * sys_bytes)
            self.queue.finish()
            t2 = time.perf_counter()
            self.run_getSurfMorse(nSystems=nch)
            self.queue.finish()
            t3 = time.perf_counter()
            aforce = np.empty((nch, self.nvecs, 4), dtype=np.float32)
            self.fromGPU('aforce', aforce)
            self.queue.finish()
            out_total[i0:i0+nch] = -aforce[:, :self.natoms, 3].sum(axis=1)
            t4 = time.perf_counter()
            t_prep += (t2 - t1)
            t_kernel += (t3 - t2)
            t_download += (t4 - t3)
        z = np.zeros_like(out_total)
        return {'total': out_total, 'LJ': z.copy(), 'Coulomb': z.copy(), 'HBond': z.copy(), 't_prep_s': t_prep, 't_kernel_s': t_kernel, 't_download_s': t_download, 't_total_s': t_prep + t_kernel + t_download}

    def eval_rigid_getSurfMorse_components(self, transforms, chunk_size=None, components=('pauli', 'london', 'coulomb')):
        out = {}
        total = None
        timing = None
        for ck in components:
            res = self.eval_rigid_getSurfMorse(transforms, chunk_size=chunk_size, component=ck)
            arr = np.asarray(res['total'], dtype=np.float32)
            out[ck] = arr
            total = arr.copy() if total is None else (total + arr)
            timing = res
        if total is None:
            total = np.zeros(len(np.asarray(transforms).reshape(-1, 12)), dtype=np.float32)
        out['total'] = total
        if timing is not None:
            out['t_prep_s'] = timing['t_prep_s']
            out['t_kernel_s'] = timing['t_kernel_s']
            out['t_download_s'] = timing['t_download_s']
            out['t_total_s'] = timing['t_total_s']
        return out

    def realloc_scan(self, n, na=-1):
        sz_f  = 4
        buffs = {
            "poss":    (sz_f*4 * n),
            "forces":  (sz_f*4 * n),
        }
        if na>0:
            buffs.update({
                "apos":   (sz_f*4 * na),
                "aREQs":  (sz_f*4 * na),
            })
        self.try_make_buffers(buffs)

    def scanNonBond(self, pos, force, REQH, ffpar, bRealloc=True ):
        n = len(pos)
        if bRealloc:  self.realloc_scan(n)
        # Upload data to GPU
        self.toGPU_( self.poss_buff,   pos)
        self.toGPU_( self.forces_buff, force)
        # Get kernel
        kernel = self.prg.scanNonBond
        # Set arguments
        kernel.set_args(
            np.int32(n),
            cl_array.vec.make_float4(*REQH),
            self.poss_buff,
            self.forces_buff,
            cl_array.vec.make_float8(*ffpar)
        )
        # Run kernel
        global_size = (n,)
        local_size = None
        cl.enqueue_nd_range_kernel(self.queue, kernel, global_size, local_size)
        result = self.fromGPU_( self.forces_buff, shape=(n, 4))
        self.queue.finish()
        return result

    def scanNonBond2(self, pos, force, apos, aREQs, REQH0, ffpar, bRealloc=True, nPBC=None, lvec=None, name=""):
        n  = len(pos)
        na = len(apos)
        if bRealloc:  self.realloc_scan(n, na=na)
        # Upload data to GPU
        self.toGPU_( self.poss_buff,   pos)
        self.toGPU_( self.forces_buff, force)
        self.toGPU_( self.apos_buff,   apos)
        self.toGPU_( self.aREQs_buff,  aREQs)  
        # Get kernel

        if nPBC is not None:
            # Convert lvec to cl_Mat3 structure (3 float4 vectors)
            lvec_cl = np.zeros((3,4), dtype=np.float32)
            lvec_cl[:,:3] = lvec
            npbc_cl = np.zeros((4), dtype=np.int32)
            npbc_cl[:3] = nPBC
            #kernel = self.prg.scanNonBond2PBC
            kernel = self.prg.scanNonBond2PBC_2
            kernel.set_args(
                np.int32(n),
                cl_array.vec.make_float4(*REQH0),
                self.poss_buff,
                self.forces_buff,
                np.int32(na),
                self.apos_buff,
                self.aREQs_buff,
                cl_array.vec.make_float8(*ffpar),
                lvec_cl,
                npbc_cl,
            )
        else:
            kernel = self.prg.scanNonBond2        
            kernel.set_args(
                np.int32(n),
                cl_array.vec.make_float4(*REQH0),
                self.poss_buff,
                self.forces_buff,
                np.int32(na),
                self.apos_buff,
                self.aREQs_buff,
                cl_array.vec.make_float8(*ffpar)
            )        
        nloc=32
        local_size  = (nloc,)
        global_size = (clu.roundup_global_size(n, nloc),)
        T0 = time.time()
        cl.enqueue_nd_range_kernel(self.queue, kernel, global_size, local_size)
        self.queue.finish()
        T = time.time() - T0

        if nPBC is not None:
            npbc=(nPBC[0]*2+1)*(nPBC[1]*2+1)*(nPBC[2]*2+1)
            ntot = n*na*npbc
            #print(f"scanNonBond2() {name:<15} | {(T*1.e+9/ntot):>2.6f} [ns/op] {(T*1.e+12/ntot):>4.6f} [TOPS] | ntot: {ntot:<12}  np: {n:<6} na: {na:<6} nPBC({npbc:<6},{nPBC}) time: {T:3.6f} [s]")
            print(f"scanNonBond2PBC() {name:<15} | {T*1.e+9/ntot:>8.4f} [ns/op] {(ntot/(T*1.e+9)):>8.4f} [GOPS] | ntot: {ntot:>12} np: {n:>6} na: {na:>6} nPBC({npbc:>6},{nPBC}) time: {T:>8.4f} [s]")
        else:
            ntot = n*na
            #print(f"scanNonBond2() {name:<15} | {(T*1.e+9/ntot):>2.6f} [ns/op] {(T*1.e+12/ntot):>4.6f} [TOPS] | ntot: {ntot:<12}  np: {n:<6} na: {na:<6} nPBC({npbc:<6},{nPBC}) time: {T:3.6f} [s]")
            print(f"scanNonBond2() {name:<15} | {T*1.e+9/ntot:>8.4f} [ns/op] {(ntot/(T*1.e+9)):>8.4f} [GOPS] | ntot: {ntot:>12} np: {n:>6} na: {na:>6} time: {T:>8.4f} [s]")

        # Download results
        result = self.fromGPU_( self.forces_buff, shape=(n, 4))
        return result