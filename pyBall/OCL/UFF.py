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

    def __init__(self, nloc=32, kernel_path=None, bPrint=False):
        super().__init__(nloc=nloc)
        if kernel_path is None:
            base_path = os.path.dirname(os.path.abspath(__file__))
            rel_path = "../../cpp/common_resources/cl/UFF.cl"
            kernel_path = os.path.join(base_path, rel_path)
        if not self.load_program(kernel_path=kernel_path, bPrint=bPrint):
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
        self.na_Tot, self.nb_Tot, self.nd_Tot, self.ni_Tot = nA, nB, nD, nInv
        self.check_buf("apos", nA * 4 * f32sz)
        self.check_buf("fapos", nA * 4 * f32sz)
        self.check_buf("fint", nA * 4 * f32sz)
        self.check_buf("atype", nA * i32sz)
        self.check_buf("REQs", nA * 4 * f32sz)
        self.check_buf("bonAtoms", nB * 2 * i32sz)
        self.check_buf("bonParams", nB * 2 * f32sz)
        self.check_buf("angles", nAng * 3 * i32sz)
        self.check_buf("angParams1", nAng * 4 * f32sz)
        self.check_buf("angParams2_w", nAng * f32sz)
        self.check_buf("angAtoms", nAng * 4 * i32sz)
        self.check_buf("angNgs", nAng * 4 * i32sz)
        self.check_buf("dihedrals", nD * 4 * i32sz)
        self.check_buf("dihParams", nD * 3 * f32sz)
        self.check_buf("dihAtoms", nD * 4 * i32sz)
        self.check_buf("dihNgs", nD * 4 * i32sz)
        self.check_buf("inversions", nInv * 4 * i32sz)
        self.check_buf("invParams", nInv * 4 * f32sz)
        self.check_buf("invAtoms", nInv * 4 * i32sz)
        self.check_buf("invNgs", nInv * 4 * i32sz)
        self.check_buf("neighs", nA * 4 * i32sz)
        self.check_buf("neighCell", nA * 4 * i32sz)
        self.check_buf("neighBs", nA * 4 * i32sz)
        self.check_buf("hneigh", nA * 4 * f32sz)
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
            ang_params1 = np.ascontiguousarray(ang_params[:, :4])
            ang_params2_w = np.ascontiguousarray(ang_params[:, 4])
            cl.enqueue_copy(self.queue, self.buffer_dict["angParams1"], ang_params1, device_offset=(angle_offset * 4 * f32sz))
            cl.enqueue_copy(self.queue, self.buffer_dict["angParams2_w"], ang_params2_w, device_offset=(angle_offset * 1 * f32sz))
        _upload_if_present("dihAtoms", "dihAtoms", np.int32, dihedral_offset * 4, i32sz)
        _upload_if_present("dihParams", "dihParams", np.float32, dihedral_offset * 3, f32sz)
        _upload_if_present("invAtoms", "invAtoms", np.int32, inversion_offset * 4, i32sz)
        _upload_if_present("invParams", "invParams", np.float32, inversion_offset * 4, f32sz)
        _upload_if_present("neighs", "neighs", np.int32, atom_offset * 4, i32sz)
        _upload_if_present("neighBs", "neighBs", np.int32, atom_offset * 4, i32sz)

    def run_eval_step(self, bClearForce=True):
        if not self.args_setup:
            self.prepare_kernel_args()
        queue = self.queue
        if bClearForce: cl.enqueue_fill_buffer(queue, self.buffer_dict["fapos"], np.float32(0), 0, self.natoms * self.nSystems * 4 * f32sz)
        if self.bDoBonds: self.prg.evalBondsAndHNeigh_UFF(queue, (self.natoms * self.nSystems,), None, *self.kernel_args["evalBondsAndHNeigh_UFF"])
        if self.bDoAngles: self.prg.evalAngles_UFF(queue, (self.nangles * self.nSystems,), None, *self.kernel_args["evalAngles_UFF"])
        if self.bDoDihedrals: self.prg.evalDihedrals_UFF(queue, (self.ndihedrals * self.nSystems,), None, *self.kernel_args["evalDihedrals_UFF"])
        if self.bDoInversions: self.prg.evalInversions_UFF(queue, (self.ninversions * self.nSystems,), None, *self.kernel_args["evalInversions_UFF"])
        if self.bDoNonBonded: self.prg.evalNonBonded(queue, (self.natoms * self.nSystems,), None, *self.kernel_args["evalNonBonded"])
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
            self.kernel_params['bSubtractVdW'] = np.int32(0) # Default value
            self.kernel_params['i0bon'] = np.int32(0)
            self.kernel_params['i0ang'] = np.int32(0)
            self.kernel_params['i0dih'] = np.int32(0)
            self.kernel_params['i0inv'] = np.int32(0)
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
        for kernel_name in self.kernelheaders:
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
        Maps atom interactions to force pieces using a buckets structure.
        Similar to UFF::mapAtomInteractions in C++.

        Args:
            natoms (int): Number of atoms
            dihedrals (np.ndarray): Dihedral indices
            inversions (np.ndarray): Inversion indices
            angles (np.ndarray): Angle indices

        Returns:
            tuple: (a2f_offsets, a2f_counts, a2f_indices) arrays for GPU upload
        """
        ndihedrals = len(dihedrals)
        ninversions = len(inversions)
        nangles = len(angles)
        # Initialize arrays
        a2f_counts = np.zeros(natoms, dtype=np.int32)

        # Count interactions per atom
        # For dihedrals (4 atoms per dihedral)
        for i in range(ndihedrals):
            a2f_counts[dihedrals[i, 0]] += 1
            a2f_counts[dihedrals[i, 1]] += 1
            a2f_counts[dihedrals[i, 2]] += 1
            a2f_counts[dihedrals[i, 3]] += 1

        # For inversions (4 atoms per inversion)
        for i in range(ninversions):
            a2f_counts[inversions[i, 0]] += 1
            a2f_counts[inversions[i, 1]] += 1
            a2f_counts[inversions[i, 2]] += 1
            a2f_counts[inversions[i, 3]] += 1

        # For angles (3 atoms per angle)
        for i in range(nangles):
            a2f_counts[angles[i, 0]] += 1
            a2f_counts[angles[i, 1]] += 1
            a2f_counts[angles[i, 2]] += 1

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
                atom_idx                   = dihedrals[i, j]
                offset                     = a2f_offsets[atom_idx] + a2f_counts_temp[atom_idx]
                a2f_indices[offset]        = i
                a2f_counts_temp[atom_idx] += 1

        # Fill indices for inversions
        for i in range(ninversions):
            for j in range(4):
                atom_idx                   = inversions[i, j]
                offset                     = a2f_offsets[atom_idx] + a2f_counts_temp[atom_idx]
                a2f_indices[offset]        = i + ndihedrals  # Offset by ndihedrals
                a2f_counts_temp[atom_idx] += 1

        # Fill indices for angles
        for i in range(nangles):
            for j in range(3):
                atom_idx                   = angles[i, j]
                offset                     = a2f_offsets[atom_idx] + a2f_counts_temp[atom_idx]
                a2f_indices[offset]        = i + ndihedrals + ninversions  # Offset by ndihedrals + ninversions
                a2f_counts_temp[atom_idx] += 1

        return a2f_offsets, a2f_counts, a2f_indices
