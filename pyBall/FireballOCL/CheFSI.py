import numpy as np
import pyopencl as cl
import pyopencl.array as cl_array
import time
import os
from ..OCL.OpenCLBase import select_device

class FrontierSolver:
    def __init__(self, ctx=None, queue=None, n_vecs=8, tile_size=32, verbose=True):
        # Prefer NVIDIA device selection if context not provided
        if ctx is None:
            self.ctx = select_device(preferred_vendor='nvidia', bPrint=False)
        else:
            self.ctx = ctx
        if queue is None:
            self.queue = cl.CommandQueue(self.ctx)
        else:
            self.queue = queue
        self.n_vecs = int(n_vecs)
        self.tile_size = int(tile_size)
        self.verbose = verbose
        
        # Load Kernel Source
        kernel_path = os.path.join(os.path.dirname(__file__), 'cl/CheFSI.cl')
        with open(kernel_path, 'r') as f:
            src = f.read()

        build_options = [
            f"-DM={self.n_vecs}",
            f"-DK={self.tile_size}",
            f"-DWG_RED_SIZE=256",
            "-DBLOCK_SIZE=64" 
        ]
        
        if self.verbose:
            print(f"[Solver] Compiling kernels with M={self.n_vecs}, K={self.tile_size}...")
            
        # Build program on chosen context
        self.prg = cl.Program(self.ctx, src).build(options=build_options)
        self.k_spmm = cl.Kernel(self.prg, "spmm_ellpack_chebyshev")
        self.k_tsg = cl.Kernel(self.prg, "tall_skinny_gram")
        self.k_sum_blocks = cl.Kernel(self.prg, "sum_partial_blocks_parallel")
        self.k_subspace_rotate = cl.Kernel(self.prg, "subspace_rotate")
        self.k_axpy = cl.Kernel(self.prg, "axpy_block")
        
        self.matvec_count = 0
        self.d_partial_sums = None
        self.d_final_gram = None

    @staticmethod
    def estimate_spectral_bounds(indices, values, safety_factor=1.1):
        n_atoms = indices.shape[0]
        row_indices = np.arange(n_atoms).reshape(-1, 1)
        valid_mask = (indices != -1)
        diag_mask = (indices == row_indices) & valid_mask
        off_diag_mask = (~diag_mask) & valid_mask
        
        h_ii = np.sum(values * diag_mask, axis=1)
        r_i = np.sum(np.abs(values) * off_diag_mask, axis=1)
        
        center = (np.max(h_ii + r_i) + np.min(h_ii - r_i)) / 2.0
        radius = (np.max(h_ii + r_i) - np.min(h_ii - r_i)) / 2.0
        return (center - radius * safety_factor, center + radius * safety_factor)

    def prepare(self, indices, H_values, S_values=None):
        self.n_atoms = indices.shape[0]
        self.n_max_neighs = indices.shape[1]
        
        self.h_indices = indices.astype(np.int32)
        self.h_H_val = H_values.astype(np.float32)
        self.d_indices = cl_array.to_device(self.queue, self.h_indices)
        self.d_H_val = cl_array.to_device(self.queue, self.h_H_val)
        
        if S_values is not None:
            self.h_S_val = S_values.astype(np.float32)
            self.d_S_val = cl_array.to_device(self.queue, self.h_S_val)
        else:
            s_host = np.zeros_like(H_values, dtype=np.float32)
            for i in range(self.n_atoms):
                for c in range(self.n_max_neighs):
                    if indices[i, c] == i: s_host[i, c] = 1.0
            self.h_S_val = s_host
            self.d_S_val = cl_array.to_device(self.queue, s_host)

        shape = (self.n_atoms, self.n_vecs)
        self.d_X = cl_array.zeros(self.queue, shape, dtype=np.float32)
        self.d_X_prev = cl_array.zeros(self.queue, shape, dtype=np.float32)
        self.d_Work = cl_array.zeros(self.queue, shape, dtype=np.float32)
        self.d_Work2 = cl_array.zeros(self.queue, shape, dtype=np.float32)
        
        self.n_groups = (self.n_atoms + self.tile_size - 1) // self.tile_size
        self.d_partial_sums = cl_array.zeros(self.queue, self.n_groups * self.n_vecs**2, dtype=np.float32)
        self.d_final_gram = cl_array.zeros(self.queue, self.n_vecs**2, dtype=np.float32)
        
        self.spmm_global = ((self.n_atoms + 63)//64 * 64, )
        self.spmm_local = (64, )

    def _spmm(self, d_matrix_val, d_in, d_out):
        self.matvec_count += 1
        self.k_spmm(
            self.queue, self.spmm_global, self.spmm_local,
            np.int32(self.n_atoms), np.int32(self.n_max_neighs), np.int32(self.n_vecs),
            self.d_indices.data, d_matrix_val.data,
            d_in.data, d_in.data, d_out.data,
            np.float32(0), np.float32(1), np.int32(0)
        )

    def _calc_gram_matrix(self, d_A, d_B):
        global_size_1 = (self.n_groups * self.n_vecs, self.n_vecs)
        local_size_1 = (self.n_vecs, self.n_vecs)
        
        evt1 = self.k_tsg(
            self.queue, global_size_1, local_size_1,
            np.int32(self.n_atoms),
            d_A.data, d_B.data, self.d_partial_sums.data
        )
        
        global_size_2 = (self.n_vecs * self.n_vecs * 256, )
        local_size_2 = (256, )
        
        evt2 = self.k_sum_blocks(
            self.queue, global_size_2, local_size_2,
            np.int32(self.n_groups),
            self.d_partial_sums.data,
            self.d_final_gram.data,
            wait_for=[evt1]
        )
        
        gram_flat = self.d_final_gram.get()
        return gram_flat.reshape(self.n_vecs, self.n_vecs)

    def check_gram_matrix_accuracy(self, d_A, d_B, gpu_result):
        """Sanity check: Compare GPU reduction with CPU numpy dot product"""
        h_A = d_A.get()
        h_B = d_B.get()
        cpu_result = h_A.T @ h_B
        diff = np.max(np.abs(cpu_result - gpu_result))
        return diff

    def _orthogonalize(self, d_V):
        # 1. Compute Metric O = V^T * S * V
        self._spmm(self.d_S_val, d_V, self.d_Work) # Work = S*V
        O = self._calc_gram_matrix(d_V, self.d_Work) # V^T * Work
        
        # Diagnostics
        if self.verbose:
            # Check kernel accuracy (expensive, only for debug)
            # diff = self.check_gram_matrix_accuracy(d_V, self.d_Work, O)
            # print(f"    [Gram Check] GPU vs CPU diff: {diff:.2e}")
            pass

        # Check Condition Number to detect Mode Collapse
        cond_num = np.linalg.cond(O)
        
        # 2. CPU Inverse Sqrt
        evals, evecs = np.linalg.eigh(O)
        
        # Guard against collapse
        evals = np.maximum(evals, 1e-12)
        inv_sqrt = evecs @ np.diag(1.0/np.sqrt(evals)) @ evecs.T
        
        # 3. Apply Rotation
        d_mat = cl_array.to_device(self.queue, inv_sqrt.astype(np.float32))
        self.k_subspace_rotate(
            self.queue, self.spmm_global, None,
            np.int32(self.n_atoms), np.int32(self.n_vecs),
            d_V.data, d_mat.data, self.d_Work.data
        )
        cl.enqueue_copy(self.queue, d_V.data, self.d_Work.data)

        # Report orthogonality BEFORE correction (how bad was it?)
        off = O - np.eye(self.n_vecs, dtype=np.float32)
        max_non_orth = np.max(np.abs(off))
        
        return max_non_orth, cond_num

    def solve(self, target_energy, bounds=None, n_iter=10, filter_order=20, tol=1e-4):
        if bounds is None:
            bounds = self.estimate_spectral_bounds(self.h_indices, self.h_H_val)
        e_min, e_max = bounds
        R_val = max(abs(e_max - target_energy), abs(e_min - target_energy))
        R2 = (R_val * 1.1)**2
        alpha = 1.8 / R2 
        
        if self.verbose: 
            print(f"[Solver] Bounds={bounds}, Target={target_energy}, alpha={alpha:.2e}")

        # Init
        rand_host = np.random.randn(self.n_atoms, self.n_vecs).astype(np.float32)
        cl.enqueue_copy(self.queue, self.d_X.data, rand_host)
        
        # Initial Ortho
        self._orthogonalize(self.d_X)
        
        logs = {'max_non_orth': [], 'cond': []}
        
        for it in range(n_iter):
            t_start = time.time()
            
            # Gradient Descent Loop
            for k in range(filter_order):
                # 1. R = (H - eS) X
                self._spmm(self.d_H_val, self.d_X, self.d_X_prev) # H X
                self._spmm(self.d_S_val, self.d_X, self.d_Work)   # S X
                self.k_axpy(self.queue, (self.n_atoms * self.n_vecs,), None,
                    np.int32(self.n_atoms * self.n_vecs), np.float32(1.0), self.d_X_prev.data, 
                    np.float32(-target_energy), self.d_Work.data, self.d_Work2.data) # Work2 = R
                
                # 2. G = (H - eS) R
                self._spmm(self.d_H_val, self.d_Work2, self.d_X_prev) # H R
                self._spmm(self.d_S_val, self.d_Work2, self.d_Work)   # S R
                self.k_axpy(self.queue, (self.n_atoms * self.n_vecs,), None,
                    np.int32(self.n_atoms * self.n_vecs), np.float32(1.0), self.d_X_prev.data, 
                    np.float32(-target_energy), self.d_Work.data, self.d_Work2.data) # Work2 = G
                
                # 3. X = X - alpha * G
                self.k_axpy(self.queue, (self.n_atoms * self.n_vecs,), None,
                    np.int32(self.n_atoms * self.n_vecs), np.float32(1.0), self.d_X.data, 
                    np.float32(-alpha), self.d_Work2.data, self.d_X.data)

            # Check Condition and Orthogonalize
            max_badness, cond = self._orthogonalize(self.d_X)
            logs['max_non_orth'].append(max_badness)
            logs['cond'].append(cond)
            
            # Rayleigh-Ritz
            self._spmm(self.d_H_val, self.d_X, self.d_Work) # HX
            H_sub = self._calc_gram_matrix(self.d_X, self.d_Work)
            
            evals, evecs = np.linalg.eigh(H_sub)
            
            # Sort
            dist = np.abs(evals - target_energy)
            idx = np.argsort(dist)
            evals = evals[idx]
            evecs = evecs[:, idx]
            
            # Rotate
            # Ensure C-contiguous layout; LAPACK returns Fortran-ordered by default
            d_rot = cl_array.to_device(self.queue, np.ascontiguousarray(evecs.astype(np.float32)))
            self.k_subspace_rotate(
                self.queue, self.spmm_global, None,
                np.int32(self.n_atoms), np.int32(self.n_vecs),
                self.d_X.data, d_rot.data, self.d_Work.data
            )
            cl.enqueue_copy(self.queue, self.d_X.data, self.d_Work.data)
            
            if self.verbose:
                print(f"Iter {it+1}: E0={evals[0]:.4f} Gap={dist[0]:.2e} NonOrth={max_badness:.2e} Cond={cond:.1e}")
                
            if max_badness < tol and cond < 100.0:
                print("Converged.")
                break

        logs['matvec_count'] = self.matvec_count
        return evals, self.d_X, logs