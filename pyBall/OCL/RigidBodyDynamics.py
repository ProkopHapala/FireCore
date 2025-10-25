import os
import numpy as np
import pyopencl as cl

from .OpenCLBase import OpenCLBase


DEFAULT_WORKGROUP_SIZE = 32
DEFAULT_MAX_ATOMS_PER_BODY = 128


def _ensure_float4(arr, w_value=0.0):
    arr = np.asarray(arr, dtype=np.float32)
    if arr.ndim != 2:
        raise ValueError(f"Expected 2D array for float4 conversion, got shape {arr.shape}")
    if arr.shape[1] == 3:
        w = np.full((arr.shape[0], 1), np.float32(w_value), dtype=np.float32)
        arr = np.hstack((arr, w))
    if arr.shape[1] != 4:
        raise ValueError(f"Expected array with 4 columns after padding, got shape {arr.shape}")
    return np.ascontiguousarray(arr, dtype=np.float32)


def _ensure_cl_mat3(mat, n_bodies):
    mat = np.asarray(mat, dtype=np.float32)
    if mat.shape[:2] != (n_bodies, 3):
        raise ValueError(f"Expected inertia tensor shape (n_bodies,3,3) or (n_bodies,3,4), got {mat.shape}")
    if mat.shape[2] == 3:
        out = np.zeros((n_bodies, 3, 4), dtype=np.float32)
        out[:, :, :3] = mat
        return out
    if mat.shape[2] == 4:
        return np.ascontiguousarray(mat, dtype=np.float32)
    raise ValueError(f"Unsupported inertia tensor trailing dimension {mat.shape[2]}")


class RigidBodyDynamics(OpenCLBase):
    """
    Simple pyOpenCL wrapper around `rigid_body_dynamics_kernel`.
    Each rigid body is simulated within a single workgroup.
    """

    def __init__(self, nloc=DEFAULT_WORKGROUP_SIZE, max_atoms=DEFAULT_MAX_ATOMS_PER_BODY):
        if nloc != DEFAULT_WORKGROUP_SIZE:
            raise ValueError(f"Kernel expects workgroup size {DEFAULT_WORKGROUP_SIZE}, got {nloc}")
        super().__init__(nloc=nloc, device_index=0)
        base_path = os.path.dirname(os.path.abspath(__file__))
        rel_path = "../../cpp/common_resources/cl/Rigid.cl"
        if not self.load_program(rel_path=rel_path, base_path=base_path, bPrint=False, bMakeHeaders=False):
            raise RuntimeError("Failed to load Rigid.cl kernel")

        self.max_atoms_per_body = max_atoms
        self.n_bodies = 0
        self.num_atoms = 0

        self.kernelheaders = {
            "rigid_body_dynamics_kernel": """__kernel
void rigid_body_dynamics_kernel(
    __global const float4*  pos_in,
    __global const float4*  quats_in,
    __global const float4*  lin_mom_in,
    __global const float4*  ang_mom_in,
    __global const float*   mass,
    __global const float*   inv_mass,
    __global const cl_Mat3* I_body_inv,
    __global const float4*  all_atom_pos_body,
    __global float4* pos_out,
    __global float4* quats_out,
    __global float4* lin_mom_out,
    __global float4* ang_mom_out,
    __global float4* all_atom_pos_world,
    const int num_atoms,
    const int num_integration_steps,
    const float dt
)"""
        }

        self.kernel_params = {}
        self.kernel_args = None

    def realloc(self, n_bodies, num_atoms):
        if num_atoms > self.max_atoms_per_body: raise ValueError(f"num_atoms={num_atoms} exceeds max_atoms_per_body={self.max_atoms_per_body}")
        self.n_bodies = int(n_bodies)
        self.num_atoms = int(num_atoms)

        float_size = np.float32().itemsize
        mat3_size = 3 * 4 * float_size
        atom_block_size = self.max_atoms_per_body * 4 * float_size

        mf = cl.mem_flags
        bytes_per_body = 4 * float_size

        self.create_buffer('pos_in',     self.n_bodies * bytes_per_body, mf.READ_WRITE)
        self.create_buffer('quats_in',   self.n_bodies * bytes_per_body, mf.READ_WRITE)
        self.create_buffer('lin_mom_in', self.n_bodies * bytes_per_body, mf.READ_WRITE)
        self.create_buffer('ang_mom_in', self.n_bodies * bytes_per_body, mf.READ_WRITE)
        self.create_buffer('mass',       self.n_bodies * float_size,     mf.READ_ONLY)
        self.create_buffer('inv_mass',   self.n_bodies * float_size,     mf.READ_ONLY)
        self.create_buffer('I_body_inv', self.n_bodies * mat3_size,      mf.READ_ONLY)
        self.create_buffer('all_atom_pos_body',  self.n_bodies * atom_block_size, mf.READ_ONLY)

        self.create_buffer('pos_out',     self.n_bodies * bytes_per_body, mf.READ_WRITE)
        self.create_buffer('quats_out',   self.n_bodies * bytes_per_body, mf.READ_WRITE)
        self.create_buffer('lin_mom_out', self.n_bodies * bytes_per_body, mf.READ_WRITE)
        self.create_buffer('ang_mom_out', self.n_bodies * bytes_per_body, mf.READ_WRITE)
        self.create_buffer('all_atom_pos_world', self.n_bodies * atom_block_size, mf.READ_WRITE)

        self.kernel_params = {
            'num_atoms':             np.int32(self.num_atoms),
            'num_integration_steps': np.int32(1),
            'dt':                    np.float32(0.01),
        }
        self.kernel_args = self.generate_kernel_args("rigid_body_dynamics_kernel")

    def upload_state(self, pos, quats, lin_mom, ang_mom, mass, inv_mass, inertia_inv, atom_pos_body):
        if self.n_bodies == 0:
            raise RuntimeError("Call realloc() before uploading state")

        pos_in   = _ensure_float4(pos)
        quats_in = _ensure_float4(quats)
        lin_in   = _ensure_float4(lin_mom)
        ang_in   = _ensure_float4(ang_mom)

        mass = np.asarray(mass, dtype=np.float32).reshape(self.n_bodies)
        inv_mass = np.asarray(inv_mass, dtype=np.float32).reshape(self.n_bodies)

        inertia_inv = _ensure_cl_mat3(inertia_inv, self.n_bodies)

        atoms = np.asarray(atom_pos_body, dtype=np.float32)
        if atoms.shape != (self.n_bodies, self.num_atoms, 3) and atoms.shape != (self.n_bodies, self.num_atoms, 4):
            raise ValueError(f"Expected body atom positions shape ({self.n_bodies},{self.num_atoms},3/4), got {atoms.shape}")

        if atoms.shape[2] == 3:
            pad = np.zeros((self.n_bodies, self.num_atoms, 1), dtype=np.float32)
            atoms = np.concatenate((atoms, pad), axis=2)

        atoms_packed = np.zeros((self.n_bodies, self.max_atoms_per_body, 4), dtype=np.float32)
        atoms_packed[:, :self.num_atoms, :] = atoms

        self.toGPU('pos_in', pos_in)
        self.toGPU('quats_in', quats_in)
        self.toGPU('lin_mom_in', lin_in)
        self.toGPU('ang_mom_in', ang_in)
        self.toGPU('mass', mass)
        self.toGPU('inv_mass', inv_mass)
        self.toGPU('I_body_inv', inertia_inv)
        self.toGPU('all_atom_pos_body', atoms_packed)
        self.queue.finish()

    def run(self, num_steps, dt):
        if self.kernel_args is None:
            raise RuntimeError("Kernel arguments not initialized; call realloc() first")

        self.kernel_params['num_integration_steps'] = np.int32(num_steps)
        self.kernel_params['dt'] = np.float32(dt)
        self.kernel_args = self.generate_kernel_args("rigid_body_dynamics_kernel")

        global_size = (self.roundUpGlobalSize(self.n_bodies * self.nloc),)
        local_size = (self.nloc,)

        self.prg.rigid_body_dynamics_kernel(self.queue, global_size, local_size, *self.kernel_args)
        self.queue.finish()

    def download_outputs(self):
        pos         = np.empty((self.n_bodies, 4), dtype=np.float32)
        quats       = np.empty((self.n_bodies, 4), dtype=np.float32)
        lin_mom     = np.empty((self.n_bodies, 4), dtype=np.float32)
        ang_mom     = np.empty((self.n_bodies, 4), dtype=np.float32)
        atoms_world = np.empty((self.n_bodies, self.max_atoms_per_body, 4), dtype=np.float32)
        self.fromGPU('pos_out', pos)
        self.fromGPU('quats_out', quats)
        self.fromGPU('lin_mom_out', lin_mom)
        self.fromGPU('ang_mom_out', ang_mom)
        self.fromGPU('all_atom_pos_world', atoms_world)
        self.queue.finish()
        return {
            'pos':     pos,
            'quats':   quats,
            'lin_mom': lin_mom,
            'ang_mom': ang_mom,
            'atom_positions': atoms_world[:, :self.num_atoms, :],
        }

    def sync_outputs_to_inputs(self):
        cl.enqueue_copy(self.queue, self.buffer_dict['pos_in'],            self.buffer_dict['pos_out'])
        cl.enqueue_copy(self.queue, self.buffer_dict['quats_in'],          self.buffer_dict['quats_out'])
        cl.enqueue_copy(self.queue, self.buffer_dict['lin_mom_in'],        self.buffer_dict['lin_mom_out'])
        cl.enqueue_copy(self.queue, self.buffer_dict['ang_mom_in'],        self.buffer_dict['ang_mom_out'])
        cl.enqueue_copy(self.queue, self.buffer_dict['all_atom_pos_body'], self.buffer_dict['all_atom_pos_world'])
        self.queue.finish()
