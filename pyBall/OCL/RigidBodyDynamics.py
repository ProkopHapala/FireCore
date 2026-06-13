import os
import ast
import numpy as np
import pyopencl as cl

from .OpenCLBase import OpenCLBase
from .InteractionEnergy import load_xyz_with_REQs


DEFAULT_WORKGROUP_SIZE = 32
DEFAULT_MAX_ATOMS_PER_BODY = 128
DEFAULT_ALPHA_MORSE = 1.5


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


def _reqs_to_plq(reqs, alpha=DEFAULT_ALPHA_MORSE):
    """Convert REQ parameters to PLQ coefficients for GridFF sampling.
    
    CRITICAL: This function expects REQ.y to be sqrt(EvdW), NOT raw EvdW.
    If reading from ElementTypes.dat, you MUST sqrt the E value before calling this.
    
    Formula (matching C++ REQ2PLQ in Forces.h):
        e  = exp(alpha * R)
        cL = e * E              # London coefficient
        cP = e * cL = e^2 * E   # Pauli coefficient
        cH = e^2 * H           # H-bond coefficient (usually 0)
    
    The sqrt(E) convention ensures proper mixed interaction:
        Eij = sqrt(Ei * Ej) when GridFF channels contain substrate sqrt(Ej)
    
    Args:
        reqs: (n, 4) array of (R, sqrt(EvdW), Q, H) - E MUST be sqrt(EvdW)
        alpha: alphaMorse parameter (default DEFAULT_ALPHA_MORSE, must match GridFF generation)
    
    Returns:
        (n, 4) array of PLQ coefficients (cP, cL, Q, cH)
    """
    reqs = np.asarray(reqs, dtype=np.float32)
    if reqs.ndim != 2 or reqs.shape[1] != 4:
        raise ValueError(f"Expected REQs shape (n,4), got {reqs.shape}")
    e = np.exp(alpha * reqs[:, 0]).astype(np.float32)
    cL = e * reqs[:, 1]
    cP = e * cL
    cH = e * e * reqs[:, 3]
    out = np.zeros_like(reqs, dtype=np.float32)
    out[:, 0] = cP
    out[:, 1] = cL
    out[:, 2] = reqs[:, 2]
    out[:, 3] = cH
    return out


def _plq_to_coeffs(plq):
    plq = np.asarray(plq, dtype=np.float32)
    if plq.ndim != 2 or plq.shape[1] != 4:
        raise ValueError(f"Expected PLQ shape (n,4), got {plq.shape}")
    return {
        'Pauli': plq[:, 0],
        'London': plq[:, 1],
        'Coulomb': plq[:, 2],
        'Hb': plq[:, 3],
    }


def _guess_mass(enames):
    mass_table = {
        'H': 1.0079, 'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998, 'Na': 22.990,
        'Mg': 24.305, 'Al': 26.982, 'Si': 28.085, 'P': 30.974, 'S': 32.06, 'Cl': 35.45,
        'K': 39.098, 'Ca': 40.078, 'Br': 79.904, 'I': 126.904,
    }
    masses = np.zeros(len(enames), dtype=np.float32)
    for i, e in enumerate(enames):
        if e not in mass_table:
            raise KeyError(f"Missing atomic mass for element '{e}'")
        masses[i] = mass_table[e]
    return masses


def _load_npy_legacy(fname):
    with open(fname, 'rb') as f:
        magic = f.read(6)
        if magic != b'\x93NUMPY':
            raise ValueError(f"Unsupported grid file magic in {fname}")
        major = int.from_bytes(f.read(1), 'little')
        minor = int.from_bytes(f.read(1), 'little')
        if major == 1:
            hlen = int.from_bytes(f.read(2), 'little')
        elif major in (2, 3):
            hlen = int.from_bytes(f.read(4), 'little')
        else:
            raise ValueError(f"Unsupported npy version {(major, minor)} in {fname}")
        header = f.read(hlen).decode('latin1').strip()
        if not header:
            raise ValueError(f"Empty npy header in {fname}")
        meta = ast.literal_eval(header)
        descr = meta.get('descr', None)
        shape = meta.get('shape', None)
        fortran_order = bool(meta.get('fortran_order', False))
        if descr is None or shape is None:
            raise ValueError(f"Incomplete npy header in {fname}: {meta}")
        dtype = np.dtype(descr)
        count = int(np.prod(shape))
        data = np.fromfile(f, dtype=dtype, count=count)
        if data.size != count:
            raise ValueError(f"Unexpected data size in {fname}: expected {count}, got {data.size}")
        arr = data.reshape(shape, order='F' if fortran_order else 'C')
        return arr


def compute_mass_properties(rel_positions, masses):
    rel = np.asarray(rel_positions, dtype=np.float32)
    m = np.asarray(masses, dtype=np.float32)
    if rel.ndim != 2 or rel.shape[1] != 3:
        raise ValueError(f"Expected relative positions shape (n,3), got {rel.shape}")
    if m.shape != (rel.shape[0],):
        raise ValueError(f"Expected masses shape ({rel.shape[0]},), got {m.shape}")
    mtot = float(m.sum())
    if mtot <= 0.0:
        raise ValueError(f"Non-positive total mass {mtot}")
    I = np.zeros((3, 3), dtype=np.float32)
    for mi, r in zip(m, rel):
        rr = np.dot(r, r)
        I += mi * (rr * np.eye(3, dtype=np.float32) - np.outer(r, r).astype(np.float32))
    det = np.linalg.det(I)
    if abs(det) < 1e-10:
        raise ValueError(f"Singular inertia tensor det={det}")
    Iinv = np.linalg.inv(I).astype(np.float32)
    return np.float32(mtot), I.astype(np.float32), Iinv


def _quat_to_matrix_np(q):
    q = np.asarray(q, dtype=np.float32)
    if q.shape == (4,):  # Single quaternion
        x, y, z, w = q
        xx, yy, zz = x * x, y * y, z * z
        xy, xz, yz = x * y, x * z, y * z
        wx, wy, wz = w * x, w * y, w * z
        return np.array([
            [1.0 - 2.0 * (yy + zz), 2.0 * (xy - wz),       2.0 * (xz + wy)],
            [2.0 * (xy + wz),       1.0 - 2.0 * (xx + zz), 2.0 * (yz - wx)],
            [2.0 * (xz - wy),       2.0 * (yz + wx),       1.0 - 2.0 * (xx + yy)],
        ], dtype=np.float32)
    elif q.ndim == 2 and q.shape[1] == 4:  # Multiple quaternions (N, 4)
        x, y, z, w = q[:, 0], q[:, 1], q[:, 2], q[:, 3]
        xx, yy, zz = x * x, y * y, z * z
        xy, xz, yz = x * y, x * z, y * z
        wx, wy, wz = w * x, w * y, w * z
        # Return (N, 3, 3) array of rotation matrices
        return np.stack([
            np.stack([1.0 - 2.0 * (yy + zz), 2.0 * (xy - wz), 2.0 * (xz + wy)], axis=1),
            np.stack([2.0 * (xy + wz), 1.0 - 2.0 * (xx + zz), 2.0 * (yz - wx)], axis=1),
            np.stack([2.0 * (xz - wy), 2.0 * (yz + wx), 1.0 - 2.0 * (xx + yy)], axis=1),
        ], axis=2).astype(np.float32)
    else:
        raise ValueError(f"Quaternion must have shape (4,) or (N,4), got {q.shape}")


class RigidBodyDynamics(OpenCLBase):
    """
    Simple pyOpenCL wrapper around `rigid_body_dynamics_kernel`.
    Each rigid body is simulated within a single workgroup.
    """

    def __init__(self, nloc=DEFAULT_WORKGROUP_SIZE, max_atoms=DEFAULT_MAX_ATOMS_PER_BODY, debug=False):
        if nloc != DEFAULT_WORKGROUP_SIZE:
            raise ValueError(f"Kernel expects workgroup size {DEFAULT_WORKGROUP_SIZE}, got {nloc}")
        super().__init__(nloc=nloc, device_index=0)
        base_path = os.path.dirname(os.path.abspath(__file__))
        rel_path = "../../cpp/common_resources/cl/Rigid.cl"
        build_options = ['-D', f'RIGID_DBG={1 if debug else 0}']
        if not self.load_program(rel_path=rel_path, base_path=base_path, bPrint=False, bMakeHeaders=False, build_options=build_options):
            raise RuntimeError("Failed to load Rigid.cl kernel")

        self.debug = bool(debug)
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
)""",
            "rigid_body_gridff_kernel": """__kernel
void rigid_body_gridff_kernel(
    __global const int*      mols,
    __global float4*         poss,
    __global float4*         qrots,
    __global float4*         vposs,
    __global float4*         vrots,
    __global const cl_Mat3*  I_body_inv,
    __global const float4*   apos_body,
    __global float4*         apos_world,
    __global const float4*   atom_PLQ,
    __global const float4*   BsplinePLQ,
    __global float4*         atom_force,
    __global float4*         body_force,
    __global float4*         body_torque,
    __global const float4*   anchors,
    const int4               grid_ns,
    const float4             grid_invStep,
    const float4             grid_p0,
    const float              dt,
    const float4             md_params,
    const int                niter
)"""
        }

        self.kernel_params = {}
        self.kernel_args = None
        self.gridff_args = None
        self.grid_shape = None
        self.grid_data = None
        self.grid_p0 = None
        self.grid_step = None
        self.atom_PLQ = None
        self.atom_REQ = None
        self.enames = None
        self.atom_masses = None
        self.atom_types_assigned = None
        self.last_atom_force = None
        self.last_body_force = None
        self.last_body_torque = None
        self.atom_body_host = None
        self.mass_total = None
        self.inertia_inv_host = None

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
        self.create_buffer('atom_PLQ',   total_atom_bytes, mf.READ_ONLY)
        self.create_buffer('atom_force', total_atom_bytes, mf.READ_WRITE)
        self.create_buffer('body_force', self.n_bodies * bytes_per_body, mf.READ_WRITE)
        self.create_buffer('body_torque', self.n_bodies * bytes_per_body, mf.READ_WRITE)

        self.kernel_params = {
            'natoms': np.int32(self.total_atoms),
            'niter': np.int32(1),
            'dt': np.float32(0.01),
            'Efield': np.zeros(4, dtype=np.float32),
        }
        self.kernel_args = self.generate_kernel_args("rigid_body_dynamics_kernel")
        self.gridff_args = None

    def upload_state(self, pos, quats, lin_mom, ang_mom, mass, inv_mass, inertia_inv, atom_pos_body, anchors=None, atom_PLQ=None):
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
        self.atom_body_host = atoms_body.copy()
        self.mass_total = float(pos_in[0, 3]) if len(pos_in) else None
        self.inertia_inv_host = inertia_inv.copy()

        mols = np.arange(0, self.total_atoms + 1, self.num_atoms, dtype=np.int32)

        self.toGPU('mols', mols)
        self.toGPU('poss', pos_in)
        self.toGPU('qrots', quats_in)
        self.toGPU('vposs', lin_in)
        self.toGPU('vrots', ang_in)
        self.toGPU('I_body_inv', inertia_inv)
        self.toGPU('apos_body', atoms_body)
        if atom_PLQ is not None:
            plq = _ensure_float4(atom_PLQ)
            if plq.shape[0] != self.total_atoms:
                raise ValueError(f"atom_PLQ length {plq.shape[0]} does not match total atoms {self.total_atoms}")
            self.atom_PLQ = plq.copy()
            self.toGPU('atom_PLQ', self.atom_PLQ)

        # GPU already recomputes apos_world from apos_body+qrots in every kernel step.
        # No need to precompute on CPU - just upload zeros; kernel overwrites on first step.
        world_atoms_flat = np.zeros((self.total_atoms, 4), dtype=np.float32)
        world_atoms_flat[:, 3] = atoms_body[:, 3]  # preserve w (charge/mass)
        # NOTE: CPU backup below (kept for reference/debugging)
        # atoms  = atoms_body.reshape(self.n_bodies, self.num_atoms, 4)
        # rot_mats = _quat_to_matrix_np(quats_in)              # (n_bodies, 3, 3)
        # rotated = np.einsum('bij,bkj->bik', atoms[:, :, :3], rot_mats)
        # world_atoms_flat[:, :3] = (rotated + pos_in[:, :3][:, None, :]).reshape(self.total_atoms, 3)

        if anchors is None:
            anchors = np.zeros_like(world_atoms_flat)
            anchors[:, 3] = -1.0
        else:
            anchors = _ensure_float4(anchors, w_value=-1.0)
            if anchors.shape[0] != self.total_atoms:
                raise ValueError(f"anchors array length {anchors.shape[0]} does not match total atoms {self.total_atoms}")
            anchors = anchors.copy()
        
        self.anchors = anchors
        self.upload_anchors()

        self.toGPU('apos_world', world_atoms_flat)
        self.toGPU('atom_force', np.zeros_like(world_atoms_flat))
        self.toGPU('body_force', np.zeros((self.n_bodies, 4), dtype=np.float32))
        self.toGPU('body_torque', np.zeros((self.n_bodies, 4), dtype=np.float32))
        self.queue.finish()

    def upload_anchors(self):
        self.toGPU('anchors', self.anchors)

    def init_gridff(self, bspline_data, grid_p0, grid_step):
        arr = np.asarray(bspline_data)
        if arr.ndim != 4 or arr.shape[3] not in (3, 4):
            raise ValueError(f"Expected Bspline grid shape (nx,ny,nz,3/4), got {arr.shape}")
        if arr.shape[3] == 3:
            tmp = np.zeros(arr.shape[:3] + (4,), dtype=np.float32)
            tmp[..., :3] = arr.astype(np.float32)
            arr = tmp
        else:
            arr = arr.astype(np.float32, copy=False)
        self.grid_shape = tuple(int(v) for v in arr.shape[:3])
        self.grid_data = np.ascontiguousarray(arr, dtype=np.float32)
        self.grid_p0 = np.array([*grid_p0, 0.0], dtype=np.float32)
        self.grid_step = np.array([*grid_step, 0.0], dtype=np.float32)
        self.buffer_dict['BsplinePLQ'] = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=self.grid_data)
        self.kernel_params['grid_ns'] = np.array([*self.grid_shape, 0], dtype=np.int32)
        self.kernel_params['grid_invStep'] = np.array([1.0 / v for v in grid_step] + [0.0], dtype=np.float32)
        self.kernel_params['grid_p0'] = self.grid_p0
        self.kernel_params['md_params'] = np.array([0.92, 0.88, 1.0, 1.0], dtype=np.float32)
        self.gridff_args = self.generate_kernel_args("rigid_body_gridff_kernel")
        self.krnl_gridff = cl.Kernel(self.prg, "rigid_body_gridff_kernel")

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

    def run_gridff(self, num_steps, dt, lin_damp=0.92, ang_damp=0.88, force_scale=1.0, torque_scale=1.0):
        if self.gridff_args is None:
            raise RuntimeError("GridFF kernel arguments not initialized; call init_gridff(...) first")
        self.kernel_params['dt'] = np.float32(dt)
        self.kernel_params['niter'] = np.int32(num_steps)
        self.kernel_params['md_params'] = np.array([lin_damp, ang_damp, force_scale, torque_scale], dtype=np.float32)
        self.gridff_args = self.generate_kernel_args("rigid_body_gridff_kernel")
        global_size = (self.roundUpGlobalSize(self.n_bodies * self.nloc),)
        local_size = (self.nloc,)
        self.krnl_gridff(self.queue, global_size, local_size, *self.gridff_args)
        self.queue.finish()

    def download_outputs(self):
        pos         = np.empty((self.n_bodies, 4), dtype=np.float32)
        quats       = np.empty((self.n_bodies, 4), dtype=np.float32)
        lin_mom     = np.empty((self.n_bodies, 4), dtype=np.float32)
        ang_mom     = np.empty((self.n_bodies, 4), dtype=np.float32)
        atoms_world = np.empty((self.total_atoms, 4), dtype=np.float32)
        atom_force  = np.empty((self.total_atoms, 4), dtype=np.float32)
        body_force  = np.empty((self.n_bodies, 4), dtype=np.float32)
        body_torque = np.empty((self.n_bodies, 4), dtype=np.float32)

        self.fromGPU('poss', pos)
        self.fromGPU('qrots', quats)
        self.fromGPU('vposs', lin_mom)
        self.fromGPU('vrots', ang_mom)
        self.fromGPU('apos_world', atoms_world)
        self.fromGPU('atom_force', atom_force)
        self.fromGPU('body_force', body_force)
        self.fromGPU('body_torque', body_torque)
        self.queue.finish()

        atoms_world = atoms_world.reshape(self.n_bodies, self.num_atoms, 4)
        atom_force = atom_force.reshape(self.n_bodies, self.num_atoms, 4)
        self.last_atom_force = atom_force
        self.last_body_force = body_force
        self.last_body_torque = body_torque

        return {
            'pos': pos,
            'quats': quats,
            'lin_mom': lin_mom,
            'ang_mom': ang_mom,
            'atom_positions': atoms_world,
            'atom_force': atom_force,
            'body_force': body_force,
            'body_torque': body_torque,
        }

    def download_selected(self, fields):
        req = tuple(fields)
        out = {}
        if 'pos' in req:
            buf = np.empty((self.n_bodies, 4), dtype=np.float32)
            self.fromGPU('poss', buf)
            out['pos'] = buf
        if 'quats' in req:
            buf = np.empty((self.n_bodies, 4), dtype=np.float32)
            self.fromGPU('qrots', buf)
            out['quats'] = buf
        if 'lin_mom' in req:
            buf = np.empty((self.n_bodies, 4), dtype=np.float32)
            self.fromGPU('vposs', buf)
            out['lin_mom'] = buf
        if 'ang_mom' in req:
            buf = np.empty((self.n_bodies, 4), dtype=np.float32)
            self.fromGPU('vrots', buf)
            out['ang_mom'] = buf
        if 'atom_positions' in req:
            buf = np.empty((self.total_atoms, 4), dtype=np.float32)
            self.fromGPU('apos_world', buf)
            out['atom_positions'] = buf.reshape(self.n_bodies, self.num_atoms, 4)
        if 'atom_force' in req:
            buf = np.empty((self.total_atoms, 4), dtype=np.float32)
            self.fromGPU('atom_force', buf)
            out['atom_force'] = buf.reshape(self.n_bodies, self.num_atoms, 4)
        if 'body_force' in req:
            buf = np.empty((self.n_bodies, 4), dtype=np.float32)
            self.fromGPU('body_force', buf)
            self.last_body_force = buf
            out['body_force'] = buf
        if 'body_torque' in req:
            buf = np.empty((self.n_bodies, 4), dtype=np.float32)
            self.fromGPU('body_torque', buf)
            self.last_body_torque = buf
            out['body_torque'] = buf
        self.queue.finish()
        return out

    def sync_outputs_to_inputs(self):
        self.queue.finish()

    def get_debug_dict(self):
        out = {
            'enames': list(self.enames) if self.enames is not None else None,
            'atom_types': list(self.atom_types_assigned) if self.atom_types_assigned is not None else None,
            'REQ': None if self.atom_REQ is None else np.array(self.atom_REQ, copy=True),
            'PLQ': None if self.atom_PLQ is None else np.array(self.atom_PLQ, copy=True),
            'PLQ_coeffs': None if self.atom_PLQ is None else _plq_to_coeffs(self.atom_PLQ),
            'masses': None if self.atom_masses is None else np.array(self.atom_masses, copy=True),
            'grid_shape': self.grid_shape,
            'grid_p0': None if self.grid_p0 is None else np.array(self.grid_p0, copy=True),
            'grid_step': None if self.grid_step is None else np.array(self.grid_step, copy=True),
            'last_atom_force': None if self.last_atom_force is None else np.array(self.last_atom_force, copy=True),
            'last_body_force': None if self.last_body_force is None else np.array(self.last_body_force, copy=True),
            'last_body_torque': None if self.last_body_torque is None else np.array(self.last_body_torque, copy=True),
        }
        return out

    @classmethod
    def from_xyz_and_grid(cls, mol_file, grid_file, substrate_xyz, n_bodies=1, body_positions=None, quats=None, alpha_morse=DEFAULT_ALPHA_MORSE, debug=False, type_map=None, mass_trans=1.0, mass_rot=None):
        apos, reqs, enames, _, _ = load_xyz_with_REQs(mol_file, type_map=type_map)
        masses = _guess_mass(enames)
        apos = np.asarray(apos, dtype=np.float32)
        com0 = (apos * masses[:, None]).sum(axis=0) / masses.sum()
        rel = apos - com0[None, :]
        mtot, _, Iinv = compute_mass_properties(rel, masses)
        mass_trans = float(mass_trans)
        if mass_trans <= 0.0:
            raise ValueError(f"mass_trans must be > 0, got {mass_trans}")
        if mass_rot is None:
            mass_rot = mass_trans
        mass_rot = float(mass_rot)
        if mass_rot <= 0.0:
            raise ValueError(f"mass_rot must be > 0, got {mass_rot}")
        Iinv_relax = Iinv * (mtot / mass_rot)
        if body_positions is None:
            body_positions = np.repeat(com0[None, :], n_bodies, axis=0).astype(np.float32)
        else:
            body_positions = np.asarray(body_positions, dtype=np.float32)
            if body_positions.shape != (n_bodies, 3):
                raise ValueError(f"Expected body_positions shape ({n_bodies},3), got {body_positions.shape}")
        pos4 = np.zeros((n_bodies, 4), dtype=np.float32)
        pos4[:, :3] = body_positions
        pos4[:, 3] = mass_trans
        quat4 = np.zeros((n_bodies, 4), dtype=np.float32)
        quat4[:, 3] = 1.0
        if quats is not None:
            q = _ensure_float4(quats)
            if q.shape[0] != n_bodies:
                raise ValueError(f"Expected quats shape ({n_bodies},4), got {q.shape}")
            quat4[:] = q
        zero4 = np.zeros((n_bodies, 4), dtype=np.float32)
        atom_body = np.repeat(rel[None, :, :], n_bodies, axis=0).astype(np.float32)
        atom_plq_single = _reqs_to_plq(reqs, alpha=alpha_morse)
        atom_plq = np.repeat(atom_plq_single[None, :, :], n_bodies, axis=0).reshape(n_bodies * len(enames), 4)
        try:
            grid = np.load(grid_file)
        except Exception:
            grid = _load_npy_legacy(grid_file)
        
        # Read lattice vectors from comment line manually
        with open(substrate_xyz, 'r') as f:
            lines = f.readlines()
            comment = lines[1].strip()
            lvec = None
            if "lvec:" in comment:
                idx = comment.find("lvec:") + 5
                parts = comment[idx:].split()
            elif "lvs" in comment:
                idx = comment.find("lvs") + 3
                parts = comment[idx:].split()
            else:
                parts = []
            
            try:
                vals = [float(v) for v in parts if v.strip()]
                if len(vals) >= 9:
                    lvec = np.array(vals[:9]).reshape(3,3).astype(np.float32)
            except ValueError:
                pass
            
            if lvec is None:
                raise ValueError(f"Substrate lattice vectors missing in {substrate_xyz}")
                
        ax = float(np.linalg.norm(lvec[0]))
        ay = float(np.linalg.norm(lvec[1]))
        az = float(np.linalg.norm(lvec[2]))
        if abs(lvec[0][1]) > 1e-6 or abs(lvec[1][0]) > 1e-6 or abs(lvec[0][2]) > 1e-6 or abs(lvec[1][2]) > 1e-6:
            raise ValueError(f"Only orthorhombic xy substrate cells supported for now, got lvec={lvec}")
        grid_step = (ax / grid.shape[0], ay / grid.shape[1], az / grid.shape[2])
        grid_p0 = (0.0, 0.0, 0.0)
        rbd = cls(debug=debug)
        rbd.realloc(n_bodies=n_bodies, num_atoms=len(enames))
        rbd.enames = list(enames)
        rbd.atom_types_assigned = [type_map.get(e, e) if type_map is not None else e for e in enames]
        rbd.atom_REQ = reqs.copy()
        rbd.atom_masses = masses.copy()
        rbd.mass_physical = float(mtot)
        rbd.mass_trans = mass_trans
        rbd.mass_rot = mass_rot
        rbd.atom_PLQ = atom_plq.copy()
        rbd.upload_state(pos4, quat4, zero4, zero4, mass_trans, 1.0 / mass_trans, np.repeat(Iinv_relax[None, :, :], n_bodies, axis=0), atom_body, atom_PLQ=atom_plq)
        rbd.init_gridff(grid, grid_p0=grid_p0, grid_step=grid_step)
        return rbd

    def reset_pose(self, pos, quats, lin_mom=None, ang_mom=None):
        if self.atom_body_host is None or self.inertia_inv_host is None or self.mass_total is None:
            raise RuntimeError("RigidBodyDynamics.reset_pose() requires prior upload_state() initialization")
        pos_in = _ensure_float4(pos)
        quats_in = _ensure_float4(quats)
        if lin_mom is None:
            lin_mom = np.zeros((self.n_bodies, 4), dtype=np.float32)
        if ang_mom is None:
            ang_mom = np.zeros((self.n_bodies, 4), dtype=np.float32)
        self.upload_state(
            pos_in,
            quats_in,
            lin_mom,
            ang_mom,
            self.mass_total,
            1.0 / self.mass_total,
            self.inertia_inv_host[:, :, :3],
            self.atom_body_host.reshape(self.n_bodies, self.num_atoms, 4)[:, :, :3],
            atom_PLQ=self.atom_PLQ,
        )

    def update_anchors(self, anchors_world):
        self.anchors = _ensure_float4(anchors_world, w_value=-1.0)
        self.upload_anchors()

    def upload_anchors(self):
        self.toGPU('anchors', self.anchors)
        self.queue.finish()
