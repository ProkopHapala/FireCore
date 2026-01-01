import numpy as np
import pyopencl as cl
import os
from ..OCL.OpenCLBase import select_device

"""
# OMM_OCL: Orbital Minimization Method using PyOpenCL

## Mathematical Background
The Orbital Minimization Method (OMM) finds the electronic ground state by minimizing an unconstrained 
energy functional with respect to localized molecular orbitals (LMOs). For a given Hamiltonian H 
and overlap matrix S, the energy functional is:
    E = Tr[ (2S - S^2) H ]  (approximate)
or more precisely, it involves orthogonalization via a penalty or a direct inverse:
    \\tilde{C} = C (C^T S C)^{-1/2}

This implementation focuses on calculating gradients of the OMM functional for sparse, localized orbitals.
The gradient calculation is split into 3 kernels to avoid atomic operations and maximize throughput:
1. **apply_operator**: Computes dual vectors |\\phi'> = Op |\\phi> on expanded supports.
2. **project_orbitals**: Computes scalar overlaps <\\phi_i | \\phi'_j> between orbitals.
3. **assemble_gradient**: Combines scalars and dual vectors to form the final gradient on the original support.

## Data Layout
- **Supports**: Discrete sets of atom indices where an orbital has non-zero coefficients.
- **Expanded Supports**: The union of neighbors of atoms in the original support, ensuring all 
  interactions (H, S) are captured when applying operators.
- **Dual Vectors**: Intermediate vectors living on expanded supports.
"""

class OMM_OCL:
    def __init__(self, ctx=None, queue=None):
        """
        Initialize the OMM solver. 
        Selects an OpenCL device (preferring NVIDIA) and builds kernels.
        """
        if ctx is None:
            # Prefer NVIDIA silently; fallback to default device 0 without interactive prompt
            self.ctx = select_device(preferred_vendor='nvidia', bPrint=False)
        else:
            self.ctx = ctx
        
        if queue is None:
            self.queue = cl.CommandQueue(self.ctx)
        else:
            self.queue = queue

        self.load_kernels()
        self.buffers = {}

    def load_kernels(self):
        """Builds OpenCL kernels from OMM.cl."""
        cl_path = os.path.join(os.path.dirname(__file__), 'cl/OMM.cl')
        with open(cl_path, 'r') as f:
            source = f.read()
        self.program = cl.Program(self.ctx, source).build()
        self.k_apply_operator = self.program.apply_operator
        self.k_project_orbitals = self.program.project_orbitals
        self.k_assemble_gradient = self.program.assemble_gradient

    def allocate_buffers(self, nOrbs, nPairs, nMaxSupport, nMaxSupportExt, nAtoms, nMaxNeigh):
        """
        Preallocates GPU buffers to avoid overhead during iterations.
        """
        mf = cl.mem_flags
        
        # Orbital data
        self.buffers['orbSupport']    = cl.Buffer(self.ctx, mf.READ_WRITE, nOrbs * nMaxSupport * 4)
        self.buffers['orbCoefs']      = cl.Buffer(self.ctx, mf.READ_WRITE, nOrbs * nMaxSupport * 4)
        self.buffers['nSupport']      = cl.Buffer(self.ctx, mf.READ_WRITE, nOrbs * 4)
        
        # Expanded orbital data
        self.buffers['orbSupportExt'] = cl.Buffer(self.ctx, mf.READ_WRITE, nOrbs * nMaxSupportExt * 4)
        self.buffers['nSupportExt']   = cl.Buffer(self.ctx, mf.READ_WRITE, nOrbs * 4)
        
        # Intermediate 'twiddle' (dual) vectors
        self.buffers['twiddleH']      = cl.Buffer(self.ctx, mf.READ_WRITE, nOrbs * nMaxSupportExt * 4)
        self.buffers['twiddleS']      = cl.Buffer(self.ctx, mf.READ_WRITE, nOrbs * nMaxSupportExt * 4)
        
        # Pairs and scalars
        self.buffers['orbPairs']      = cl.Buffer(self.ctx, mf.READ_WRITE, nPairs * 2 * 4)
        self.buffers['pairRowPtr']    = cl.Buffer(self.ctx, mf.READ_WRITE, (nOrbs + 1) * 4)
        # scalarLambda should be float32 (4 bytes), but lambda_mat.flatten() might be float64?
        # nPairs=9, 9*4 = 36 bytes. Error says 72 bytes. 72/9 = 8 bytes -> float64.
        self.buffers['scalarLambda']  = cl.Buffer(self.ctx, mf.READ_WRITE, nPairs * 4)
        
        # System geometry/matrices (usually constant during OMM)
        self.buffers['atomNeigh']     = cl.Buffer(self.ctx, mf.READ_WRITE, nAtoms * nMaxNeigh * 4)
        self.buffers['atomH']         = cl.Buffer(self.ctx, mf.READ_WRITE, nAtoms * nMaxNeigh * 4) 
        self.buffers['atomS']         = cl.Buffer(self.ctx, mf.READ_WRITE, nAtoms * nMaxNeigh * 4)
        
        # Output gradient
        self.buffers['outGrads']      = cl.Buffer(self.ctx, mf.READ_WRITE, nOrbs * nMaxSupport * 4)

        # DEBUG: Print buffer sizes
        #for name, buf in self.buffers.items():
        #    print(f"DEBUG: Buffer {name} size={buf.size}")

    def upload_data(self, **kwargs):
        """Uploads data from host to preallocated GPU buffers."""
        for name, data in kwargs.items():
            if name in self.buffers:
                # Force float32 for lambda and coefs if they are float64 from numpy
                if name in ['scalarLambda', 'orbCoefs', 'atomH', 'atomS', 'outGrads', 'twiddleH', 'twiddleS']:
                    data_np = np.require(data, dtype=np.float32, requirements='C')
                elif name in ['orbSupport', 'orbSupportExt', 'nSupport', 'nSupportExt', 'orbPairs', 'pairRowPtr', 'atomNeigh']:
                    data_np = np.require(data, dtype=np.int32, requirements='C')
                else:
                    data_np = np.require(data, requirements='C')
                
                buf = self.buffers[name]
                if data_np.nbytes > buf.size:
                    raise ValueError(f"Data for buffer '{name}' is too large: {data_np.nbytes} > {buf.size} (dtype={data_np.dtype})")
                cl.enqueue_copy(self.queue, buf, data_np)
            else:
                print(f"[WARNING] Buffer {name} not preallocated.")

    def download_data(self, name, shape, dtype=np.float32):
        """Downloads data from GPU buffer for debugging or final results."""
        res = np.empty(shape, dtype=dtype)
        cl.enqueue_copy(self.queue, res, self.buffers[name])
        return res

    def build_expanded_supports(self, nOrbs, orbSupport, nSupport, atomNeigh):
        """
        Mathematical purpose: Ensure dual vectors have enough space to store the result 
        of H|phi> or S|phi> where the operator has a finite range.
        Algorithm: For each orbital, the expanded support is the union of neighbors 
        of all atoms in the original support.
        """
        expanded_supports = []
        for i in range(nOrbs):
            supp = orbSupport[i, :nSupport[i]]
            ext_supp = set(supp)
            for iat in supp:
                if iat < 0: continue
                neighs = atomNeigh[iat]
                for jnat in neighs:
                    if jnat >= 0:
                        ext_supp.add(jnat)
            expanded_supports.append(sorted(list(ext_supp)))
        
        nMaxSupportExt = max(len(s) for s in expanded_supports)
        orbSupportExt = np.full((nOrbs, nMaxSupportExt), -1, dtype=np.int32)
        nSupportExt = np.zeros(nOrbs, dtype=np.int32)
        
        for i, s in enumerate(expanded_supports):
            nSupportExt[i] = len(s)
            orbSupportExt[i, :len(s)] = s
            
        return orbSupportExt, nSupportExt, nMaxSupportExt

    def apply_operator(self, nOrbs, nMaxSupport, nMaxSupportExt, nMaxNeigh, bDownload=False):
        """
        Kernel 1: Computes |\phi'> = Op |\phi>
        Purpose: Apply sparse H or S matrix to localized orbitals.
        """
        local_size = (min(max(nMaxSupport, nMaxSupportExt), 256),)
        global_size = (nOrbs * local_size[0],)

        # Apply H
        self.k_apply_operator(
            self.queue, global_size, local_size,
            np.int32(nOrbs),
            self.buffers['orbSupport'], self.buffers['orbCoefs'], self.buffers['nSupport'],
            self.buffers['orbSupportExt'], self.buffers['nSupportExt'],
            self.buffers['atomNeigh'], self.buffers['atomH'],
            self.buffers['twiddleH'],
            np.int32(nMaxSupport), np.int32(nMaxSupportExt), np.int32(nMaxNeigh)
        )
        # Apply S
        self.k_apply_operator(
            self.queue, global_size, local_size,
            np.int32(nOrbs),
            self.buffers['orbSupport'], self.buffers['orbCoefs'], self.buffers['nSupport'],
            self.buffers['orbSupportExt'], self.buffers['nSupportExt'],
            self.buffers['atomNeigh'], self.buffers['atomS'],
            self.buffers['twiddleS'],
            np.int32(nMaxSupport), np.int32(nMaxSupportExt), np.int32(nMaxNeigh)
        )
        
        if bDownload:
            return (self.download_data('twiddleH', (nOrbs, nMaxSupportExt)),
                    self.download_data('twiddleS', (nOrbs, nMaxSupportExt)))
        return None

    def project_orbitals(self, nPairs, nMaxSupport, nMaxSupportExt, bH=True, bDownload=False, orbPairs=None):
        """
        Kernel 2: Computes scalar overlaps s_ij = <\phi_i | \phi'_j>
        Purpose: Project dual vectors back onto original supports.
        """
        local_size = (256,)
        global_size = (nPairs * local_size[0],)

        # Select which twiddle buffer to project
        d_twiddle = self.buffers['twiddleH'] if bH else self.buffers['twiddleS']
        
        # Select which pairs to use (default to pre-allocated all-pairs)
        if orbPairs is not None:
            d_orbPairs = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=orbPairs.astype(np.int32))
        else:
            d_orbPairs = self.buffers['orbPairs']

        out = np.zeros(nPairs, dtype=np.float32)
        d_out = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, out.nbytes)
        
        self.k_project_orbitals(
            self.queue, global_size, local_size,
            np.int32(nPairs),
            d_orbPairs,
            self.buffers['orbSupport'], self.buffers['orbCoefs'], self.buffers['nSupport'],
            self.buffers['orbSupportExt'], d_twiddle, self.buffers['nSupportExt'],
            d_out, np.int32(nMaxSupport), np.int32(nMaxSupportExt)
        )
        
        if bDownload:
            cl.enqueue_copy(self.queue, out, d_out)
            return out
        return d_out

    def assemble_gradient(self, nOrbs, nMaxSupport, nMaxSupportExt, bDownload=False):
        """
        Kernel 3: Computes final gradient vector.
        Purpose: Sum up contributions from neighbors using Lambda mixing matrix.
        """
        local_size = (256,)
        global_size = (nOrbs * local_size[0],)

        self.k_assemble_gradient(
            self.queue, global_size, local_size,
            np.int32(nOrbs),
            self.buffers['pairRowPtr'], self.buffers['orbPairs'], self.buffers['scalarLambda'],
            self.buffers['orbSupport'], self.buffers['nSupport'],
            self.buffers['orbSupportExt'], self.buffers['twiddleH'], self.buffers['twiddleS'],
            self.buffers['nSupportExt'],
            self.buffers['outGrads'],
            np.int32(nMaxSupport), np.int32(nMaxSupportExt)
        )
        
        if bDownload:
            return self.download_data('outGrads', (nOrbs, nMaxSupport))
        return None

    def get_fireball_data(self):
        """
        Retrieves sparse Hamiltonian and Overlap matrices from Fireball.
        """
        from pyBall import FireCore as fc
        dims = fc.get_HS_dims()
        sparse_data = fc.get_HS_neighs(dims) 
        fc.get_HS_sparse(dims, data=sparse_data)
        
        # Extract s-orbitals for now (simplification for testing)
        h_mat = sparse_data.h_mat[:, :, 0, 0].astype(np.float32)
        s_mat = sparse_data.s_mat[:, :, 0, 0].astype(np.float32)
        atomNeigh = sparse_data.neigh_j.astype(np.int32) - 1
        
        return h_mat, s_mat, atomNeigh, dims.natoms, dims.neigh_max

    def setup_from_fireball(self, nOrbs, orbSupport, nSupport, K_stiff=100.0, orbPairs=None):
        """
        One-liner big-picture: Initializes GPU solver using current Fireball state.
        Algorithm: Fetches sparse H/S, builds expanded supports, and preallocates GPU memory.
        """
        h_mat, s_mat, atomNeigh, nAtoms, nMaxNeigh = self.get_fireball_data()
        
        # Build expanded supports
        orbSupportExt, nSupportExt, nMaxSupportExt = self.build_expanded_supports(
            nOrbs, orbSupport, nSupport, atomNeigh
        )
        
        # Pairs for testing/calculation
        if orbPairs is None:
            # Default to all-to-all pairs for testing
            orbPairs = []
            for i in range(nOrbs):
                for j in range(nOrbs):
                    orbPairs.append([i, j])
            orbPairs = np.array(orbPairs, dtype=np.int32)
        
        nPairs = orbPairs.shape[0]
            
        pairRowPtr = np.zeros(nOrbs + 1, dtype=np.int32)
        # Build CSR-like row pointers from the orbPairs list
        # Assuming orbPairs is sorted by the first element (orbital i)
        curr_p = 0
        for i in range(nOrbs):
            pairRowPtr[i] = curr_p
            while curr_p < nPairs and orbPairs[curr_p, 0] == i:
                curr_p += 1
        pairRowPtr[nOrbs] = nPairs
            
        nMaxSupport = orbSupport.shape[1]
        
        self.allocate_buffers(nOrbs, nPairs, nMaxSupport, nMaxSupportExt, nAtoms, nMaxNeigh)
        
        # Initialize orbCoefs with something if they weren't passed (handled by test script)
        # We need to make sure we upload orbCoefs here too
        
        self.upload_data(
            orbSupport=orbSupport,
            nSupport=nSupport,
            orbSupportExt=orbSupportExt,
            nSupportExt=nSupportExt,
            atomNeigh=atomNeigh,
            atomH=h_mat,
            atomS=s_mat,
            orbPairs=orbPairs,
            pairRowPtr=pairRowPtr
        )
        
        # Meta-data for kernel launches
        self.meta = {
            'nOrbs': nOrbs,
            'nPairs': nPairs,
            'nMaxSupport': nMaxSupport,
            'nMaxSupportExt': nMaxSupportExt,
            'nMaxNeigh': nMaxNeigh,
            'K_stiff': K_stiff
        }

    def run_omm_step(self, bDownload=True, return_overlaps=False):
        """
        Executes a full 3-kernel sequence to compute gradients.
        One-liner: Core OMM pipeline (Apply Op -> Project -> Mix).
        Algorithm:
        1. |\phi^H'> = H|\phi>, |\phi^S'> = S|\phi> (apply_operator)
        2. S_ij = <\phi_i | \phi^S'_j> (project_orbitals)
        3. \lambda_{ij} = K*(S_ij - \delta_{ij}) (penalty mixing on CPU)
        4. g_i = 2 * ( \phi^H'_i - \sum_j \lambda_{ij} \phi^S'_j ) (assemble_gradient)
        """
        m = self.meta
        # 1. Apply H and S to get twiddle vectors
        self.apply_operator(m['nOrbs'], m['nMaxSupport'], m['nMaxSupportExt'], m['nMaxNeigh'])
        
        # 2. Project to get overlaps S_ij (required for Lambda/orthogonalization)
        s_ij_flat = self.project_orbitals(m['nPairs'], m['nMaxSupport'], m['nMaxSupportExt'], bH=False, bDownload=True)
        s_ij = s_ij_flat.reshape(m['nOrbs'], m['nOrbs'])
        
        # 3. CPU: Compute Lambda mixing matrix
        # For penalty-based orthogonalization:
        # Lambda_ij = K * (S_ij - delta_ij)
        lambda_mat = m['K_stiff'] * (s_ij - np.eye(m['nOrbs']))
        self.upload_data(scalarLambda=lambda_mat.flatten())
        
        h_ij = None
        if return_overlaps:
            h_ij_flat = self.project_orbitals(m['nPairs'], m['nMaxSupport'], m['nMaxSupportExt'], bH=True, bDownload=True)
            h_ij = h_ij_flat.reshape(m['nOrbs'], m['nOrbs'])
        
        # 4. Assemble final gradient
        grads = self.assemble_gradient(m['nOrbs'], m['nMaxSupport'], m['nMaxSupportExt'], bDownload=bDownload)
        
        # The functional derivative for penalty-based OMM:
        # dE/dc_i = 2 * ( H|phi_i> - K * sum_j (S_ij - delta_ij) S|phi_j> )
        # Our kernel computes ( H|phi_i> - sum_j Lambda_ij S|phi_j> )
        # So we just multiply by 2.0.
        grads = 2.0 * grads

        # IMPORTANT: Fireball Hamiltonian often has different sign convention or needs to be minimized.
        # Check if we should minimize or maximize.
        # Actually, in OMM we minimize E = Tr[C^T H C]. 
        # If diag(H) is negative, minimizing E makes it more negative (blow up) 
        # unless we have orthogonality constraints.

        if return_overlaps:
            return grads, s_ij, h_ij
        return grads

    def optimize(self, orbCoefs, steps=10, lr=1e-3, log=True):
        """
        Simple gradient-descent loop over coefficients.
        Returns history dict and final overlaps.
        """
        # Ensure coefficients uploaded
        self.upload_data(orbCoefs=orbCoefs)

        m = self.meta
        hist_maxE, hist_totE, hist_offdiag, hist_diagdev, hist_gnorm = [], [], [], [], []
        last_s = None
        last_h = None
        grads = None

        for step in range(steps):
            # Gradient descent update (negative gradient to minimize)
            # 1. Run step to get grads and overlaps
            grads, s_ij, h_ij = self.run_omm_step(bDownload=True, return_overlaps=True)
            last_s, last_h = s_ij, h_ij

            # 2. Record metrics BEFORE update (for current state)
            diag_h = np.diag(h_ij)
            maxE = np.max(diag_h)
            totE = np.sum(diag_h)
            offdiag_mask = ~np.eye(m['nOrbs'], dtype=bool)
            max_offdiag = np.max(np.abs(s_ij[offdiag_mask]))
            max_diagdev = np.max(np.abs(np.diag(s_ij) - 1.0))
            gnorm = np.sqrt(np.sum(grads**2))

            hist_maxE.append(maxE)
            hist_totE.append(totE)
            hist_offdiag.append(max_offdiag)
            hist_diagdev.append(max_diagdev)
            hist_gnorm.append(gnorm)

            if log:
                print(f"[ITER {step}] maxE={maxE:.6f} totE={totE:.6f} max|S_off|={max_offdiag:.6e} max|Sii-1|={max_diagdev:.6e} |g|={gnorm:.6e}")

            # Gradient descent update (minimize energy)
            # functional derivative g = dE/dc
            # To minimize E, we move c_new = c - lr * g
            orbCoefs = orbCoefs - lr * grads

            # 4. Renormalization: C_i = C_i / sqrt(S_ii)
            # Essential for stability in OMM
            # We re-compute S_ii after update to ensure we stay on the manifold
            # For efficiency we might reuse old S_ii, but for stability we re-project
            # push updated coeffs to GPU first to get correct S_ii
            self.upload_data(orbCoefs=orbCoefs)
            diag_pairs = np.array([[i,i] for i in range(m['nOrbs'])], dtype=np.int32)
            diag_s = self.project_orbitals(m['nOrbs'], m['nMaxSupport'], m['nMaxSupportExt'], bH=False, bDownload=True, orbPairs=diag_pairs)
            
            norm_factors = np.sqrt(np.maximum(diag_s, 1e-12))
            orbCoefs = orbCoefs / norm_factors[:, np.newaxis]
            
            # 5. Push final normalized coeffs to GPU
            self.upload_data(orbCoefs=orbCoefs)

        history = {
            'hist_maxE': hist_maxE,
            'hist_totE': hist_totE,
            'hist_offdiag': hist_offdiag,
            'hist_diagdev': hist_diagdev,
            'hist_gnorm': hist_gnorm,
        }
        return history, last_s, last_h, grads, orbCoefs

