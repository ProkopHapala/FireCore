import os
import numpy as np
import pyopencl as cl
from .OpenCLBase import OpenCLBase
from .UFFbuilder import UFF_Builder
#from .MMFF import MMFF
#from . import MMparams as mmparams    # Do we need it here ? Maybe it is enough to use it in UFFBuilder

# Size constants for better readability
i32sz = 4  # size of int32 in bytes
f32sz = 4  # size of float32 in bytes

class UFF_CL(OpenCLBase):
    """
    PyOpenCL interface for running UFF calculations on GPU.
    This class is responsible for managing OpenCL buffers and running kernels.
    Topology and parameter preparation is delegated to UFF_Builder.
    """

    def __init__(self, nloc=32, kernel_path=None, bPrint=False, debug_build_options=None):
        super().__init__(nloc=nloc)
        if kernel_path is None:
            base_path = os.path.dirname(os.path.abspath(__file__))
            rel_path = "../../cpp/common_resources/cl/UFF.cl"
            kernel_path = os.path.join(base_path, rel_path)
        if not self.load_program(kernel_path=kernel_path, bPrint=bPrint, build_options=debug_build_options):
            print(f"Failed to load UFF kernels from {kernel_path}")
            return

        base_path = os.path.dirname(os.path.abspath(__file__))
        data_path = os.path.join(base_path, "../../cpp/common_resources/")
        #self.params        = read_element_types(os.path.join(data_path, 'ElementTypes.dat'), os.path.join(data_path, 'AtomTypes.dat'))
        #self.element_types = self.params.element_types

        self.nSystems = 0
        self.natoms = 0
        self.nbonds = 0
        self.nangles = 0
        self.ndihedrals = 0
        self.ninversions = 0
        self.npbc = 0
        self.a2f_map_size = 0
        self.kernel_args = {}
        self.args_setup = False
        self.bDoBonds = True
        self.bDoAngles = True
        self.bDoDihedrals = True
        self.bDoInversions = True
        self.bDoNonBonded = False

    def toUFF(self, mol, bRealloc=True, bSimple=True, b141=True, bConj=True, bCumulene=True):
        builder = UFF_Builder(mol, bSimple=bSimple, b141=b141, bConj=bConj, bCumulene=bCumulene)
        uff_data = builder.build()
        if bRealloc:
            self.realloc_buffers(
                natoms=len(mol.apos),
                nbonds=len(uff_data['bonAtoms']),
                nangles=len(uff_data['angAtoms']),
                ndihedrals=len(uff_data['dihAtoms']),
                ninversions=len(uff_data['invAtoms']),
                npbc=0
            )
        a2f_offsets, a2f_counts, a2f_indices = self.mapAtomInteractions(len(mol.apos), uff_data['dihAtoms'], uff_data['invAtoms'], uff_data['angAtoms'])
        self.set_a2f_map_size(len(a2f_indices))
        uff_data['a2f_offsets'] = a2f_offsets
        uff_data['a2f_counts'] = a2f_counts
        uff_data['a2f_indices'] = a2f_indices
        return uff_data

    def realloc_buffers(self, natoms, nbonds, nangles, ndihedrals, ninversions, npbc, nSystems=1):
        self.nSystems = nSystems
        self.natoms = natoms
        self.nbonds = nbonds
        self.nangles = nangles
        self.ndihedrals = ndihedrals
        self.ninversions = ninversions
        self.npbc = npbc
        nA = natoms * nSystems
        nB = nbonds * nSystems
        nAng = nangles * nSystems
        nD = ndihedrals * nSystems
        nInv = ninversions * nSystems
        nf_per_system = (ndihedrals * 4) + (ninversions * 4) + (nangles * 3) + nbonds  # matching C++ UFF.h::realloc
        self.na_Tot, self.nb_Tot, self.nd_Tot, self.ni_Tot = nA, nB, nD, nInv
        self.check_buf("apos", nA * 4 * f32sz)
        self.check_buf("fapos", nA * 4 * f32sz)
        # fint stores force pieces, indexed by interaction-piece id (not by atom)
        self.check_buf("fint", (nf_per_system * nSystems) * 4 * f32sz)
        self.check_buf("atype", nA * i32sz)
        self.check_buf("REQs", nA * 4 * f32sz)
        self.check_buf("bonAtoms", nB * 2 * i32sz)
        self.check_buf("bonParams", nB * 2 * f32sz)
        self.check_buf("angles", nAng * 3 * i32sz)
        self.check_buf("angParams1", nAng * 4 * f32sz)
        self.check_buf("angParams2_w", nAng * f32sz)
        self.check_buf("angAtoms", nAng * 4 * i32sz)
        self.check_buf("angNgs", nAng * 2 * i32sz)  # int2 per angle
        self.check_buf("dihedrals", nD * 4 * i32sz)
        # UFF.cl expects dihParams as float4 per dihedral (xyz used, w ignored)
        self.check_buf("dihParams", nD * 4 * f32sz)
        self.check_buf("dihAtoms", nD * 4 * i32sz)
        self.check_buf("dihNgs", nD * 4 * i32sz)
        self.check_buf("inversions", nInv * 4 * i32sz)
        self.check_buf("invParams", nInv * 4 * f32sz)
        self.check_buf("invAtoms", nInv * 4 * i32sz)
        # UFF.cl uses padded int4 for invNgs on the C++ path; keep int4 for consistent stride
        self.check_buf("invNgs", nInv * 4 * i32sz)
        self.check_buf("neighs", nA * 4 * i32sz)
        self.check_buf("neighCell", nA * 4 * i32sz)
        self.check_buf("neighBs", nA * 4 * i32sz)
        # hneigh is [natoms*4] float4 per system
        self.check_buf("hneigh", (nA * 4) * 4 * f32sz)
        self.check_buf("pbc_shifts", npbc * nSystems * 4 * f32sz)
        self.check_buf("energies", 5 * nSystems * f32sz)
        self.check_buf("lvec", 9 * nSystems * f32sz)
        self.check_buf("params", 10 * i32sz)
        self.check_buf("Ea_contrib", nAng * f32sz)
        self.check_buf("Ed_contrib", nD * f32sz)
        self.check_buf("Ei_contrib", nInv * f32sz)
        self.args_setup = False
        print(f"UFF buffers allocated for {nSystems} systems with {natoms} atoms each")

    def upload_positions(self, positions, iSys=0):
        if positions.shape[0] != self.natoms:
            raise ValueError(f"Expected {self.natoms} atoms, got {positions.shape[0]}")
        if positions.dtype != np.float32:
            positions = positions.astype(np.float32)
        if len(positions.shape) == 1:
            positions = positions.reshape(-1, 3)
        padded_positions = np.zeros((self.natoms, 4), dtype=np.float32)
        padded_positions[:, :3] = positions
        offset = iSys * self.natoms * 4
        cl.enqueue_copy(self.queue, self.buffer_dict["apos"], padded_positions.flatten(), device_offset=offset * f32sz)

    def upload_topology_params(self, uff_data, iSys=0):
        atom_offset = iSys * self.natoms
        bond_offset = iSys * self.nbonds
        angle_offset = iSys * self.nangles
        dihedral_offset = iSys * self.ndihedrals
        inversion_offset = iSys * self.ninversions
        def _upload_if_present(buffer_name, data_key, dtype, offset_elements, element_size):
            if data_key in uff_data and uff_data[data_key] is not None and len(uff_data[data_key]) > 0:
                data = uff_data[data_key].astype(dtype)
                cl.enqueue_copy(self.queue, self.buffer_dict[buffer_name], data, device_offset=(offset_elements * element_size))
        _upload_if_present("atype", "atype", np.int32, atom_offset, i32sz)
        _upload_if_present("REQs", "REQs", np.float32, atom_offset * 4, f32sz)
        _upload_if_present("bonAtoms", "bonAtoms", np.int32, bond_offset * 2, i32sz)
        _upload_if_present("bonParams", "bonParams", np.float32, bond_offset * 2, f32sz)
        _upload_if_present("angAtoms", "angAtoms", np.int32, angle_offset * 4, i32sz)
        if 'angParams' in uff_data and uff_data['angParams'] is not None and len(uff_data['angParams']) > 0:
            ang_params = uff_data['angParams'].astype(np.float32)
            # CPU angParams layout is [K,c0,c1,c2,c3], but kernel expects:
            #   angParams1  = [c0,c1,c2,c3]
            #   angParams2_w= [K]
            ang_params1   = np.ascontiguousarray(ang_params[:, 1:5])
            ang_params2_w = np.ascontiguousarray(ang_params[:, 0])
            cl.enqueue_copy(self.queue, self.buffer_dict["angParams1"],   ang_params1,   device_offset=(angle_offset * 4 * f32sz))
            cl.enqueue_copy(self.queue, self.buffer_dict["angParams2_w"], ang_params2_w, device_offset=(angle_offset * 1 * f32sz))
        _upload_if_present("dihAtoms", "dihAtoms", np.int32, dihedral_offset * 4, i32sz)
        if 'dihParams' in uff_data and uff_data['dihParams'] is not None and len(uff_data['dihParams']) > 0:
            p3 = np.ascontiguousarray(uff_data['dihParams'].astype(np.float32))
            assert p3.ndim == 2 and p3.shape[1] == 3
            p4 = np.zeros((p3.shape[0], 4), dtype=np.float32)
            p4[:, :3] = p3
            cl.enqueue_copy(self.queue, self.buffer_dict["dihParams"], p4, device_offset=(dihedral_offset * 4 * f32sz))
        _upload_if_present("invAtoms", "invAtoms", np.int32, inversion_offset * 4, i32sz)
        _upload_if_present("invParams", "invParams", np.float32, inversion_offset * 4, f32sz)
        _upload_if_present("neighs", "neighs", np.int32, atom_offset * 4, i32sz)
        _upload_if_present("neighBs", "neighBs", np.int32, atom_offset * 4, i32sz)
        # Precomputed hneigh indices for angles/dihedrals/inversions (from bakeNeighs)
        _upload_if_present("angNgs", "angNgs", np.int32, angle_offset * 2, i32sz)  # int2 per angle
        _upload_if_present("dihNgs", "dihNgs", np.int32, dihedral_offset * 4, i32sz)
        if 'invNgs' in uff_data and uff_data['invNgs'] is not None and len(uff_data['invNgs']) > 0:
            ng3 = np.ascontiguousarray(uff_data['invNgs'].astype(np.int32))
            assert ng3.ndim == 2 and ng3.shape[1] == 3
            ng4 = np.full((ng3.shape[0], 4), -1, dtype=np.int32)
            ng4[:, :3] = ng3
            cl.enqueue_copy(self.queue, self.buffer_dict["invNgs"], ng4, device_offset=(inversion_offset * 4 * i32sz))
        _upload_if_present("a2f_offsets", "a2f_offsets", np.int32, atom_offset, i32sz)
        _upload_if_present("a2f_counts",  "a2f_counts",  np.int32, atom_offset, i32sz)
        if 'a2f_indices' in uff_data and uff_data['a2f_indices'] is not None and len(uff_data['a2f_indices']) > 0:
            # a2f_indices is per-system local; for now we upload whole map for iSys==0
            data = np.ascontiguousarray(uff_data['a2f_indices'].astype(np.int32))
            cl.enqueue_copy(self.queue, self.buffer_dict["a2f_indices"], data, device_offset=0)

    def run_eval_step(self, bClearForce=True):
        if not self.args_setup:
            self.prepare_kernel_args()
        queue = self.queue
        if bClearForce:
            cl.enqueue_fill_buffer(queue, self.buffer_dict["fapos"], np.float32(0), 0, self.natoms * self.nSystems * 4 * f32sz)
        # Clear fint when any term writes into it
        if self.bDoAngles or self.bDoDihedrals or self.bDoInversions:
            cl.enqueue_fill_buffer(queue, self.buffer_dict["fint"], np.float32(0), 0, self.buffer_dict["fint"].size)

        if self.bDoBonds:
            self.prg.evalBondsAndHNeigh_UFF(queue, (self.natoms * self.nSystems,), None, *self.kernel_args["evalBondsAndHNeigh_UFF"])
        if self.bDoAngles:
            self.prg.evalAngles_UFF(queue, (self.nangles * self.nSystems,), None, *self.kernel_args["evalAngles_UFF"])
        if self.bDoDihedrals:
            self.prg.evalDihedrals_UFF(queue, (self.ndihedrals * self.nSystems,), None, *self.kernel_args["evalDihedrals_UFF"])
        if self.bDoInversions:
            self.prg.evalInversions_UFF(queue, (self.ninversions * self.nSystems,), None, *self.kernel_args["evalInversions_UFF"])
        # Assemble forces from fint into fapos (bond forces are already accumulated in fapos)
        if self.bDoAngles or self.bDoDihedrals or self.bDoInversions:
            self.prg.assembleForces_UFF(queue, (self.natoms * self.nSystems,), None, *self.kernel_args["assembleForces_UFF"])
        if self.bDoNonBonded: self.prg.evalNonBonded(queue, (self.natoms * self.nSystems,), None, *self.kernel_args["evalNonBonded"])
        queue.finish()
        energies = np.zeros(5 * self.nSystems, dtype=np.float32)
        return energies[4::5]

    def get_forces(self, iSys=None):
        if iSys is None:
            forces = np.zeros(self.natoms * self.nSystems * 4, dtype=np.float32)
            cl.enqueue_copy(self.queue, forces, self.buffer_dict["fapos"])
            return forces.reshape(self.nSystems, self.natoms, 4)[:, :, :3]
        else:
            forces = np.zeros(self.natoms * 4, dtype=np.float32)
            offset = iSys * self.natoms * 4
            cl.enqueue_copy(self.queue, forces, self.buffer_dict["fapos"], device_offset=offset * f32sz)
            return forces.reshape(self.natoms, 4)[:, :3]

    def get_total_energy(self):
        energies = np.zeros(5 * self.nSystems, dtype=np.float32)
        cl.enqueue_copy(self.queue, energies, self.buffer_dict["energies"])
        return energies.reshape(self.nSystems, 5)[:, 4]

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
            self.kernel_params['npbc'] = np.int32(self.npbc)
            # Generic sizes used by clear kernels
            self.kernel_params['n'] = np.int32(self.natoms * self.nSystems)
            # fint layout matching C++ UFF.h::realloc: dih|inv|ang|bon
            nf_per_system = self.ndihedrals * 4 + self.ninversions * 4 + self.nangles * 3 + self.nbonds
            self.kernel_params['nf_per_system'] = np.int32(nf_per_system)
            self.kernel_params['bSubtractVdW'] = np.int32(0) # Default value
            # fint offsets computed by mapAtomInteractions (must be called first)
            self.kernel_params['i0dih'] = np.int32(getattr(self, '_i0dih', 0))
            self.kernel_params['i0inv'] = np.int32(getattr(self, '_i0inv', self.ndihedrals * 4))
            self.kernel_params['i0ang'] = np.int32(getattr(self, '_i0ang', self.ndihedrals * 4 + self.ninversions * 4))
            self.kernel_params['i0bon'] = np.int32(getattr(self, '_i0bon', self.ndihedrals * 4 + self.ninversions * 4 + self.nangles * 3))
            self.kernel_params['SubNBTorsionFactor'] = np.float32(0.0)
            self.kernel_params['Rdamp'] = np.float32(1.0)
            self.kernel_params['FmaxNonBonded'] = np.float32(10.0)
            self.kernel_params['bSubtractBondNonBond'] = np.int32(0)
            self.kernel_params['bSubtractAngleNonBond'] = np.int32(0)
            self.kernel_params['bClearForce'] = np.int32(1)

        # Use OpenCLBase's functionality for generating kernel arguments
        if not hasattr(self, 'kernelheaders') or not self.kernelheaders:
            # If kernel headers are not set, extract them from prg source
            self.kernelheaders = self.extract_kernel_headers(self.prg.get_info(cl.prg_info.SOURCE))

        self.kernel_args = {}

        # UFF.cl may contain additional kernels unrelated to UFF_CL (e.g. MMFF MD kernels).
        # Only prepare argument lists for kernels we actually call from UFF_CL.
        needed = {
            'clear_fapos_UFF',
            'clear_fint_UFF',
            'evalBondsAndHNeigh_UFF',
            'evalAngles_UFF',
            'evalDihedrals_UFF',
            'evalInversions_UFF',
            'assembleForces_UFF',
        }
        if self.bDoNonBonded:
            # Non-bonded kernels are optional; keep only if present
            needed.update({'getNonBond', 'getNonBond_ex2', 'getSurfMorse', 'getNonBond_GridFF_Bspline', 'getNonBond_GridFF_Bspline_ex2'})

        for kernel_name in needed:
            if kernel_name in self.kernelheaders:
                self.kernel_args[kernel_name] = self.generate_kernel_args(kernel_name)

        self.args_setup = True

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

    def mapAtomInteractions(self, natoms, dihedrals, inversions, angles):
        """
        Maps atom interactions to fint force-piece indices using a buckets structure.
        Matches C++ UFF::mapAtomInteractions exactly.
        a2f_indices stores fint piece indices (not interaction indices).

        fint layout (matching C++ UFF.h::realloc):
          i0dih = 0
          i0inv = 4 * ndihedrals
          i0ang = 4 * ndihedrals + 4 * ninversions
          i0bon = i0ang + 3 * nangles  (bonds don't use fint, kept for size compat)
        """
        ndihedrals  = len(dihedrals)
        ninversions = len(inversions)
        nangles     = len(angles)
        # fint offsets matching C++ UFF.h::realloc
        i0dih = 0
        i0inv = 4 * ndihedrals
        i0ang = i0inv + 4 * ninversions
        # Store offsets for later use in prepare_kernel_args
        self._i0dih = i0dih
        self._i0inv = i0inv
        self._i0ang = i0ang
        self._i0bon = i0ang + 3 * nangles

        # --- Phase 1: Count interactions per atom (same as C++) ---
        a2f_counts = np.zeros(natoms, dtype=np.int32)
        for i in range(ndihedrals):
            a2f_counts[dihedrals[i, 0]] += 1
            a2f_counts[dihedrals[i, 1]] += 1
            a2f_counts[dihedrals[i, 2]] += 1
            a2f_counts[dihedrals[i, 3]] += 1
        for i in range(ninversions):
            a2f_counts[inversions[i, 0]] += 1
            a2f_counts[inversions[i, 1]] += 1
            a2f_counts[inversions[i, 2]] += 1
            a2f_counts[inversions[i, 3]] += 1
        for i in range(nangles):
            a2f_counts[angles[i, 0]] += 1
            a2f_counts[angles[i, 1]] += 1
            a2f_counts[angles[i, 2]] += 1

        # --- Phase 2: Compute offsets ---
        total_refs = int(np.sum(a2f_counts))
        a2f_offsets = np.zeros(natoms, dtype=np.int32)
        offset = 0
        for i in range(natoms):
            a2f_offsets[i] = offset
            offset += a2f_counts[i]

        # --- Phase 3: Fill a2f_indices with fint piece indices (matching C++ addToCell) ---
        a2f_counts_temp = np.zeros(natoms, dtype=np.int32)
        a2f_indices = np.zeros(total_refs, dtype=np.int32)

        # Dihedrals: fint piece at i0dih + i*4 + j  (j=0..3 for each atom in dihedral)
        for i in range(ndihedrals):
            i0 = i * 4 + i0dih
            for j in range(4):
                atom_idx                   = dihedrals[i, j]
                offset                     = a2f_offsets[atom_idx] + a2f_counts_temp[atom_idx]
                a2f_indices[offset]        = i0 + j
                a2f_counts_temp[atom_idx] += 1

        # Inversions: fint piece at i0inv + i*4 + j  (j=0..3)
        for i in range(ninversions):
            i0 = i * 4 + i0inv
            for j in range(4):
                atom_idx                   = inversions[i, j]
                offset                     = a2f_offsets[atom_idx] + a2f_counts_temp[atom_idx]
                a2f_indices[offset]        = i0 + j
                a2f_counts_temp[atom_idx] += 1

        # Angles: fint piece at i0ang + i*3 + j  (j=0..2)
        for i in range(nangles):
            i0 = i * 3 + i0ang
            for j in range(3):
                atom_idx                   = angles[i, j]
                offset                     = a2f_offsets[atom_idx] + a2f_counts_temp[atom_idx]
                a2f_indices[offset]        = i0 + j
                a2f_counts_temp[atom_idx] += 1

        return a2f_offsets, a2f_counts, a2f_indices
