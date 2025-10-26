import os
import numpy as np
import pyopencl as cl

from .OpenCLBase import OpenCLBase


DEFAULT_WORKGROUP_SIZE = 32
DEFAULT_MAX_ATOMS_PER_BODY = 128


def _pack_float3(arr):
    vec = np.asarray(arr, dtype=np.float32)
    if vec.shape != (3,):
        raise ValueError(f"Expected array of shape (3,) for float3, got {vec.shape}")
    out = np.zeros(4, dtype=np.float32)
    out[:3] = vec
    return out


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


def _quat_to_matrix_np(q):
    q = np.asarray(q, dtype=np.float32)
    if q.shape != (4,):  raise ValueError(f"Quaternion must have shape (4,), got {q.shape}")
    x, y, z, w = q
    xx, yy, zz = x * x, y * y, z * z
    xy, xz, yz = x * y, x * z, y * z
    wx, wy, wz = w * x, w * y, w * z
    return np.array([
        [1.0 - 2.0 * (yy + zz), 2.0 * (xy - wz),       2.0 * (xz + wy)],
        [2.0 * (xy + wz),       1.0 - 2.0 * (xx + zz), 2.0 * (yz - wx)],
        [2.0 * (xz - wy),       2.0 * (yz + wx),       1.0 - 2.0 * (xx + yy)],
    ], dtype=np.float32)


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
        self.total_atoms = 0
        self.atom_counts = None
        self.mol_offsets = None
        self.max_atoms_body = 0

        self.kernelheaders = {
            "rigid_body_dynamics_kernel": """__kernel
void rigid_body_dynamics_kernel(
    __global const int*      mols,
    __global float4*         poss,
    __global float4*         qrots,
    __global float4*         vposs,
    __global float4*         vrots,
    __global const cl_Mat3*  I_body_inv,
    __global const float4*   apos_body,
    __global float4*         apos_world,
    __global const float4*   anchors,
    const int   natoms,
    const int   niter,
    const float dt,
    const float3  Efield
)"""
        }

        self.kernel_params = {}
        self.kernel_args = None

    def realloc(self, n_bodies, num_atoms):
        if num_atoms > self.max_atoms_per_body:
            raise ValueError(f"num_atoms={num_atoms} exceeds max_atoms_per_body={self.max_atoms_per_body}")
        self.n_bodies = int(n_bodies)
        self.num_atoms = int(num_atoms)
        self.total_atoms = self.n_bodies * self.num_atoms

        float_size = np.float32().itemsize
        int_size = np.int32().itemsize
        mat3_size = 3 * 4 * float_size
        atom_block_size = self.max_atoms_per_body * 4 * float_size  # kept for compatibility, actual total handled per-body
        mf = cl.mem_flags
        bytes_per_body = 4 * float_size

        self.create_buffer('mols', (self.n_bodies + 1) * int_size, mf.READ_ONLY)
        self.create_buffer('poss',   self.n_bodies * bytes_per_body, mf.READ_WRITE)
        self.create_buffer('qrots',  self.n_bodies * bytes_per_body, mf.READ_WRITE)
        self.create_buffer('vposs',  self.n_bodies * bytes_per_body, mf.READ_WRITE)
        self.create_buffer('vrots',  self.n_bodies * bytes_per_body, mf.READ_WRITE)
        self.create_buffer('I_body_inv', self.n_bodies * mat3_size,      mf.READ_ONLY)
        self.create_buffer('anchors', self.total_atoms * 4 * float_size, mf.READ_ONLY)

        total_atom_bytes = self.total_atoms * 4 * float_size
        self.create_buffer('apos_body',  total_atom_bytes, mf.READ_ONLY)
        self.create_buffer('apos_world', total_atom_bytes, mf.READ_WRITE)

        self.kernel_params = {
            'natoms': np.int32(self.total_atoms),
            'niter': np.int32(1),
            'dt': np.float32(0.01),
            'Efield': np.zeros(4, dtype=np.float32),
        }
        self.kernel_args = self.generate_kernel_args("rigid_body_dynamics_kernel")

    def upload_state(self, pos, quats, lin_mom, ang_mom, mass, inv_mass, inertia_inv, atom_pos_body, anchors=None):
        if self.n_bodies == 0:
            raise RuntimeError("Call realloc() before uploading state")

        pos_in   = _ensure_float4(pos)
        quats_in = _ensure_float4(quats)
        lin_in   = _ensure_float4(lin_mom)
        ang_in   = _ensure_float4(ang_mom)

        inertia_inv = _ensure_cl_mat3(inertia_inv, self.n_bodies)

        atoms = np.asarray(atom_pos_body, dtype=np.float32)
        if atoms.shape != (self.n_bodies, self.num_atoms, 3) and atoms.shape != (self.n_bodies, self.num_atoms, 4):
            raise ValueError(f"Expected body atom positions shape ({self.n_bodies},{self.num_atoms},3/4), got {atoms.shape}")

        if atoms.shape[2] == 3:
            pad = np.zeros((self.n_bodies, self.num_atoms, 1), dtype=np.float32)
            atoms = np.concatenate((atoms, pad), axis=2)

        atoms_body = atoms.reshape(self.total_atoms, 4)

        mols = np.arange(0, self.total_atoms + 1, self.num_atoms, dtype=np.int32)

        self.toGPU('mols', mols)
        self.toGPU('poss', pos_in)
        self.toGPU('qrots', quats_in)
        self.toGPU('vposs', lin_in)
        self.toGPU('vrots', ang_in)
        self.toGPU('I_body_inv', inertia_inv)
        self.toGPU('apos_body', atoms_body)
        if anchors is None:
            anchors = np.zeros((self.total_atoms, 4), dtype=np.float32)
        anchors = _ensure_float4(anchors)
        if anchors.shape[0] != self.total_atoms:
            raise ValueError(f"anchors array length {anchors.shape[0]} does not match total atoms {self.total_atoms}")
        self.toGPU('anchors', anchors)

        # initialize world positions consistent with current state
        atoms  = atoms_body.reshape(self.n_bodies, self.num_atoms, 4)
        world_atoms = np.zeros_like(atoms)
        # for ib in range(self.n_bodies):
        #     world_atoms[ib, :, :3] = atoms[ib, :, :3] @ _quat_to_matrix_np(quats_in[ib]).T + pos_in[ib, :3]
        #     world_atoms[ib, :, 3]  = atoms[ib, :, 3]
        # world_atoms = world_atoms.reshape(self.total_atoms, 4)
        self.toGPU('apos_world', world_atoms)
        self.queue.finish()

    def run(self, num_steps, dt, efield=None):
        if self.kernel_args is None:
            raise RuntimeError("Kernel arguments not initialized; call realloc() first")

        self.kernel_params['niter'] = np.int32(num_steps)
        self.kernel_params['dt'] = np.float32(dt)
        if efield is not None:
            self.kernel_params['Efield'] = _pack_float3(efield)
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
        atoms_world = np.empty((self.total_atoms, 4), dtype=np.float32)

        self.fromGPU('poss', pos)
        self.fromGPU('qrots', quats)
        self.fromGPU('vposs', lin_mom)
        self.fromGPU('vrots', ang_mom)
        self.fromGPU('apos_world', atoms_world)
        self.queue.finish()

        atoms_world = atoms_world.reshape(self.n_bodies, self.num_atoms, 4)

        return {
            'pos': pos,
            'quats': quats,
            'lin_mom': lin_mom,
            'ang_mom': ang_mom,
            'atom_positions': atoms_world,
        }

    def sync_outputs_to_inputs(self):
        self.queue.finish()
