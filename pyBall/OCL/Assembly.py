import numpy as np
import pyopencl as cl
from .OpenCLBase import OpenCLBase
from scipy.spatial.transform import Rotation as R

class AssemblyOCL(OpenCLBase):
    def __init__(self, nloc=128, device_index=0):
        super().__init__(nloc=nloc, device_index=device_index)
        # Load kernel
        self.load_program(rel_path="cl/Assembly.cl")
        self.krn_emit_configuration_xyz = cl.Kernel(self.prg, "emit_configuration_xyz")
        self.krn_evaluate_packing_3d    = cl.Kernel(self.prg, "evaluate_packing_3d")
        self.d_base_atoms = None
        self.natoms = 0
        
    def upload_base_atoms(self, atoms):
        # atoms should be float32 (natoms, 4) where w is radius
        self.natoms = np.int32(len(atoms))
        atoms_f32 = self._flat32(atoms)
        self.d_base_atoms = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=atoms_f32)
        
    def evaluate_packing(self, transforms, nmols, max_clash_penalty=50.0):
        # transforms: (n_confs, nmols, 4, 4) float32
        n_confs = np.int32(transforms.shape[0])
        nmols_i32 = np.int32(nmols)
        max_penalty_f32 = np.float32(max_clash_penalty)
        
        transforms_f32 = self._flat32(transforms)
        d_transforms = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=transforms_f32)
        
        results_host = np.zeros(n_confs, dtype=np.float32)
        d_results = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, results_host.nbytes)
        
        local_replica = cl.LocalMemory(self.nloc * 16) # float4
        local_scores = cl.LocalMemory(self.nloc * 4)   # float
        
        global_size = (int(n_confs * self.nloc),)
        local_size = (int(self.nloc),)
        
        self.krn_evaluate_packing_3d.set_args(
            self.d_base_atoms,
            self.natoms,
            d_transforms,
            nmols_i32,
            max_penalty_f32,
            local_replica,
            local_scores,
            d_results,
        )
        cl.enqueue_nd_range_kernel(self.queue, self.krn_evaluate_packing_3d, global_size, local_size).wait()
        
        cl.enqueue_copy(self.queue, results_host, d_results).wait()
        return results_host
        
    def emit_configuration(self, transforms_single_conf, nmols):
        # transforms_single_conf: (nmols, 4, 4) float32
        nmols_i32 = np.int32(nmols)
        transforms_f32 = self._flat32(transforms_single_conf)
        d_transforms = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=transforms_f32)
        
        out_atoms_host = np.zeros((self.natoms * nmols, 4), dtype=np.float32)
        d_out_atoms = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, out_atoms_host.nbytes)
        
        global_size = (int(self.natoms * nmols),)
        local_size = None # let OpenCL decide
        
        self.krn_emit_configuration_xyz.set_args(
            self.d_base_atoms,
            self.natoms,
            d_transforms,
            nmols_i32,
            d_out_atoms,
        )
        cl.enqueue_nd_range_kernel(self.queue, self.krn_emit_configuration_xyz, global_size, local_size).wait()
        
        cl.enqueue_copy(self.queue, out_atoms_host, d_out_atoms).wait()
        return out_atoms_host

def pack_transforms(rotmats, shifts):
    # rotmats: (..., 3, 3), shifts: (..., 3)
    shape = rotmats.shape[:-2]
    out = np.zeros(shape + (4, 4), dtype=np.float32)
    out[..., 0, :3] = rotmats[..., 0, :]
    out[..., 1, :3] = rotmats[..., 1, :]
    out[..., 2, :3] = rotmats[..., 2, :]
    out[..., 3, :3] = shifts
    return out

def super_fibonacci_rotations(N):
    phi = np.sqrt(2)
    psi = (1 + np.sqrt(5)/2 + np.sqrt((5 + 2*np.sqrt(5) + np.sqrt(25 + 20*np.sqrt(5)))/4)) ** 0.25
    quats = np.zeros((N, 4))
    for i in range(N):
        s = i + 1
        t = s / N
        d = 2 * np.pi * s
        r = np.sqrt(t)
        R_val = np.sqrt(1 - t)
        alpha = d * phi
        beta = d * psi
        quats[i] = [r * np.sin(alpha), r * np.cos(alpha), R_val * np.sin(beta), R_val * np.cos(beta)]
    return R.from_quat(quats).as_matrix()

def generate_transform_buffer(lattice_a, lattice_b, n_rot=100, n_shift=6):
    # hexagonal symmetry (6 rotations around Z)
    angles_60 = np.linspace(0, 2 * np.pi, 6, endpoint=False)
    S_sym = R.from_euler('z', angles_60).as_matrix() # (6, 3, 3)
    
    # 3x3 grid shifts, ordered so [0,0] is index 0
    u = np.array([0, -1, 1])
    uu, vv = np.meshgrid(u, u)
    L_lat = np.zeros((9, 3))
    L_lat[:, 0] = uu.flatten() * lattice_a[0] + vv.flatten() * lattice_b[0]
    L_lat[:, 1] = uu.flatten() * lattice_a[1] + vv.flatten() * lattice_b[1]
    L_lat[:, 2] = uu.flatten() * lattice_a[2] + vv.flatten() * lattice_b[2]
    
    # Base search space using Super-Fibonacci Spirals
    R_base = super_fibonacci_rotations(n_rot) # (n_rot, 3, 3)
    
    fa = np.linspace(0, 1.0, n_shift, endpoint=False)
    fb = np.linspace(0, 1.0, n_shift, endpoint=False)
    FA, FB = np.meshgrid(fa, fb, indexing='ij')
    
    T_base = np.zeros((n_shift**2, 3))
    T_base[:, 0] = FA.flatten() * lattice_a[0] + FB.flatten() * lattice_b[0]
    T_base[:, 1] = FA.flatten() * lattice_a[1] + FB.flatten() * lattice_b[1]
    T_base[:, 2] = FA.flatten() * lattice_a[2] + FB.flatten() * lattice_b[2]
    
    N_rot = len(R_base)
    N_shift = len(T_base)
    N_confs = N_rot * N_shift
    
    R_conf = np.repeat(R_base, N_shift, axis=0) # (N_confs, 3, 3)
    T_conf = np.tile(T_base, (N_rot, 1))        # (N_confs, 3)
    
    # Apply 6-fold symmetry: orient and shift
    R_sym = np.einsum('kij,cjl->ckil', S_sym, R_conf) # (N_confs, 6, 3, 3)
    T_sym = np.einsum('kij,cj->cki', S_sym, T_conf)   # (N_confs, 6, 3)
    
    # Tile across 9 lattice cells
    R_all = np.tile(R_sym, (1, 9, 1, 1, 1))           # (N_confs, 9, 6, 3, 3)
    T_all = np.zeros((N_confs, 9, 6, 3))
    for l in range(9):
        T_all[:, l, :, :] = T_sym + L_lat[l].reshape(1, 1, 3)
        
    # Flatten to 54 replicas
    R_all = R_all.reshape(N_confs, 54, 3, 3)
    T_all = T_all.reshape(N_confs, 54, 3)
    
    cl_transforms = pack_transforms(R_all, T_all) # (N_confs, 54, 4, 4)
    return cl_transforms, N_confs, R_conf, T_conf

def generate_transform_buffer_simple(lattice_a, lattice_b, n_rot=4, n_shift=4):
    """Simple version with 1 replica for smoke testing"""
    R_base = np.eye(3)
    T_base = np.zeros(3)
    R_all = np.array([[[R_base]]])
    T_all = np.array([[[T_base]]])
    cl_transforms = pack_transforms(R_all, T_all)
    return cl_transforms, 1, np.array([R_base]), np.array([T_base])
