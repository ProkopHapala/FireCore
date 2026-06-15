"""GPU block-Jacobi solver for vibration dynamic stiffness A(omega) u = f (vib_jacobi.cl)."""
import os
import numpy as np
import pyopencl as cl

from .OpenCLBase import OpenCLBase

VIB_WG_SIZE = 32
VIB_MAX_NEIGH = 32


def pack_block_jacobi_rows(A_csr, n_atoms, max_neigh=VIB_MAX_NEIGH):
    """Pack complex sparse/dense-dynamic-stiffness rows into per-atom neighbor blocks for GPU Jacobi."""
    import scipy.sparse as sp
    if not sp.issparse(A_csr):
        A_csr = sp.csr_matrix(A_csr)
    else:
        A_csr = A_csr.tocsr()
    A_dense = A_csr.toarray()
    row_neigh = np.full((n_atoms, max_neigh), -1, dtype=np.int32)
    row_blk_re = np.zeros((n_atoms, max_neigh, 3, 3), dtype=np.float32)
    row_blk_im = np.zeros((n_atoms, max_neigh, 3, 3), dtype=np.float32)
    row_count = np.zeros(n_atoms, dtype=np.int32)
    diag_inv_re = np.zeros((n_atoms, 3, 3), dtype=np.float32)
    diag_inv_im = np.zeros((n_atoms, 3, 3), dtype=np.float32)
    for i in range(n_atoms):
        si = slice(i * 3, (i + 1) * 3)
        Aii = A_dense[si, si]
        Aii_inv = np.linalg.inv(Aii)
        diag_inv_re[i] = Aii_inv.real.astype(np.float32)
        diag_inv_im[i] = Aii_inv.imag.astype(np.float32)
        cnt = 0
        for j in range(n_atoms):
            if i == j:
                continue
            sj = slice(j * 3, (j + 1) * 3)
            blk = A_dense[si, sj]
            if np.max(np.abs(blk)) < 1e-14:
                continue
            if cnt >= max_neigh:
                raise ValueError(f"atom {i} has >{max_neigh} off-diagonal blocks; raise VIB_MAX_NEIGH")
            row_neigh[i, cnt] = j
            row_blk_re[i, cnt] = blk.real.astype(np.float32)
            row_blk_im[i, cnt] = blk.imag.astype(np.float32)
            cnt += 1
        row_count[i] = cnt
    return row_neigh, row_blk_re, row_blk_im, row_count, diag_inv_re, diag_inv_im


class VibSolver(OpenCLBase):
    """pyOpenCL wrapper for vib_jacobi_solve (complex block Jacobi + heavy-ball)."""

    def __init__(self, *, workgroup_size: int = VIB_WG_SIZE, max_neigh: int = VIB_MAX_NEIGH):
        if workgroup_size <= 0 or (workgroup_size & (workgroup_size - 1)) != 0:
            raise ValueError("workgroup_size must be a positive power of two")
        if max_neigh != VIB_MAX_NEIGH:
            raise ValueError(f"vib_jacobi.cl compiled with VIB_MAX_NEIGH={VIB_MAX_NEIGH}, got {max_neigh}")
        super().__init__(nloc=workgroup_size, device_index=0)
        base_path = os.path.dirname(os.path.abspath(__file__))
        rel_path = "../../cpp/common_resources/cl/vib_jacobi.cl"
        if not self.load_program(rel_path=rel_path, base_path=base_path, bPrint=False, bMakeHeaders=False):
            raise RuntimeError("Failed to load vib_jacobi.cl")
        self.max_neigh = max_neigh
        self.n_atoms = 0
        self.kernelheaders = {
            "vib_jacobi_solve": """__kernel void vib_jacobi_solve(
    __global const int* row_count,
    __global const int* row_neigh,
    __global const float* row_blk_re,
    __global const float* row_blk_im,
    __global const float* diag_inv_re,
    __global const float* diag_inv_im,
    __global const float* rhs_re,
    __global const float* rhs_im,
    __global float* u_re,
    __global float* u_im,
    const int natoms,
    const int n_iter,
    const float b_mix
)""",
        }
        self.kernel_params = {
            'natoms': np.int32(0),
            'n_iter': np.int32(0),
            'b_mix': np.float32(0.0),
        }
        self.kernel_args: dict[str, tuple] = {}

    def realloc(self, n_atoms: int):
        self.n_atoms = int(n_atoms)
        if self.n_atoms <= 0:
            raise ValueError("n_atoms must be positive")
        fs = np.float32().itemsize
        isz = np.int32().itemsize
        mf = cl.mem_flags
        nn = self.n_atoms
        mn = self.max_neigh
        self.create_buffer('row_count', nn * isz, mf.READ_ONLY)
        self.create_buffer('row_neigh', nn * mn * isz, mf.READ_ONLY)
        self.create_buffer('row_blk_re', nn * mn * 9 * fs, mf.READ_ONLY)
        self.create_buffer('row_blk_im', nn * mn * 9 * fs, mf.READ_ONLY)
        self.create_buffer('diag_inv_re', nn * 9 * fs, mf.READ_ONLY)
        self.create_buffer('diag_inv_im', nn * 9 * fs, mf.READ_ONLY)
        self.create_buffer('rhs_re', nn * 3 * fs, mf.READ_ONLY)
        self.create_buffer('rhs_im', nn * 3 * fs, mf.READ_ONLY)
        self.create_buffer('u_re', nn * 3 * fs, mf.READ_WRITE)
        self.create_buffer('u_im', nn * 3 * fs, mf.READ_WRITE)

    def upload_system(self, row_neigh, row_blk_re, row_blk_im, row_count, diag_inv_re, diag_inv_im):
        if self.n_atoms == 0:
            raise RuntimeError("Call realloc(n_atoms) before upload_system")
        self.toGPU('row_count', np.ascontiguousarray(row_count, dtype=np.int32))
        self.toGPU('row_neigh', np.ascontiguousarray(row_neigh, dtype=np.int32))
        self.toGPU('row_blk_re', np.ascontiguousarray(row_blk_re, dtype=np.float32).reshape(-1))
        self.toGPU('row_blk_im', np.ascontiguousarray(row_blk_im, dtype=np.float32).reshape(-1))
        self.toGPU('diag_inv_re', np.ascontiguousarray(diag_inv_re, dtype=np.float32).reshape(-1))
        self.toGPU('diag_inv_im', np.ascontiguousarray(diag_inv_im, dtype=np.float32).reshape(-1))
        self.queue.finish()

    def solve(self, rhs, *, n_iter: int = 200, b_mix: float = 0.2, u0=None):
        rhs = np.asarray(rhs, dtype=np.complex128).ravel()
        if rhs.size != self.n_atoms * 3:
            raise ValueError(f"rhs size {rhs.size} != {self.n_atoms * 3}")
        rhs_re = rhs.real.astype(np.float32).reshape(self.n_atoms, 3)
        rhs_im = rhs.imag.astype(np.float32).reshape(self.n_atoms, 3)
        if u0 is None:
            u_re = np.zeros((self.n_atoms, 3), dtype=np.float32)
            u_im = np.zeros((self.n_atoms, 3), dtype=np.float32)
        else:
            u0 = np.asarray(u0, dtype=np.complex128).ravel()
            u_re = u0.real.astype(np.float32).reshape(self.n_atoms, 3)
            u_im = u0.imag.astype(np.float32).reshape(self.n_atoms, 3)
        self.toGPU('rhs_re', rhs_re.reshape(-1))
        self.toGPU('rhs_im', rhs_im.reshape(-1))
        self.toGPU('u_re', u_re.reshape(-1))
        self.toGPU('u_im', u_im.reshape(-1))
        overrides = {
            'natoms': np.int32(self.n_atoms),
            'n_iter': np.int32(n_iter),
            'b_mix': np.float32(b_mix),
        }
        args = self.generate_kernel_args('vib_jacobi_solve', overrides=overrides)
        local_size = (self.nloc,)
        global_size = (self.roundUpGlobalSize(self.nloc),)
        self.prg.vib_jacobi_solve(self.queue, global_size, local_size, *args)
        self.queue.finish()
        u_re_out = np.empty(self.n_atoms * 3, dtype=np.float32)
        u_im_out = np.empty(self.n_atoms * 3, dtype=np.float32)
        self.fromGPU('u_re', u_re_out)
        self.fromGPU('u_im', u_im_out)
        return u_re_out + 1j * u_im_out


def solve_dynamic_stiffness_gpu(H_sparse, M, omega, eta, rhs, *, n_iter=200, b_mix=0.2, stabilize=1e-6, solver=None):
    """Solve A(omega) u = rhs on GPU; returns complex displacement vector."""
    from pyBall import FTIR
    n_atoms = len(rhs) // 3
    A_csr = FTIR._dynamic_stiffness_sparse(H_sparse, M, omega, eta=eta, stabilize=stabilize)
    row_neigh, row_blk_re, row_blk_im, row_count, diag_inv_re, diag_inv_im = pack_block_jacobi_rows(A_csr, n_atoms)
    if solver is None:
        solver = VibSolver()
    if solver.n_atoms != n_atoms:
        solver.realloc(n_atoms)
    solver.upload_system(row_neigh, row_blk_re, row_blk_im, row_count, diag_inv_re, diag_inv_im)
    return solver.solve(rhs, n_iter=n_iter, b_mix=b_mix)
