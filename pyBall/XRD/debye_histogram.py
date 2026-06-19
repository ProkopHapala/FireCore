"""XRDDebye: pyOpenCL engine for pair-distance histogram + Debye transform."""
import os
import numpy as np
import pyopencl as cl
from pyBall.OCL.OpenCLBase import OpenCLBase
from .form_factors import get_form_factor_table


def load_xyz_atoms(path: str):
    """Load XYZ → (atoms_array, unique_Z) where atoms_array is (N,4) float32 [x,y,z,elem_idx]."""
    with open(path, 'r') as f:
        lines = f.readlines()
    n = int(lines[0].strip())
    elems = []
    pos = []
    for line in lines[2:2 + n]:
        parts = line.split()
        elems.append(parts[0])
        pos.append([float(parts[1]), float(parts[2]), float(parts[3])])
    pos = np.array(pos, dtype=np.float32)
    symbol_to_Z = {'H': 1, 'C': 6, 'Si': 14, 'N': 7, 'O': 8}
    Z = np.array([symbol_to_Z.get(e, 0) for e in elems], dtype=np.int32)
    if np.any(Z == 0):
        missing = {e for e, z in zip(elems, Z) if z == 0}
        raise ValueError(f"Unknown elements in XYZ: {missing}")
    unique_Z = np.unique(Z)
    Z_to_idx = {int(z): i for i, z in enumerate(unique_Z)}
    elem_idx = np.array([Z_to_idx[int(z)] for z in Z], dtype=np.float32)
    atoms = np.column_stack([pos, elem_idx]).astype(np.float32)
    return atoms, unique_Z


class XRDDebye(OpenCLBase):
    """GPU-accelerated powder XRD via pair-distance histogram + Debye transform.

    Usage:
        engine = XRDDebye()
        atoms, unique_Z = load_xyz_atoms('crystal.xyz')
        Q = np.linspace(0.5, 15.0, 500, dtype=np.float32)
        I = engine.powder_spectrum(atoms, unique_Z, Q, r_max=20.0, dr=0.02)
    """

    def __init__(self, cl_src_dir=None, nloc=32, preferred_vendor='nvidia'):
        super().__init__(nloc=nloc, preferred_vendor=preferred_vendor, bPrint=True)
        if cl_src_dir is None:
            d = os.path.dirname(os.path.abspath(__file__))
            cl_src_dir = os.path.realpath(os.path.join(d, '..', '..', 'cpp', 'common_resources', 'cl'))
        self.cl_src_dir = cl_src_dir
        kernel_path = os.path.join(cl_src_dir, 'XRDDebye.cl')
        print(f"XRDDebye: compiling {kernel_path}")
        # Reset prg so we can compile this kernel even if another was loaded by a parent
        self.prg = None
        self.load_program(kernel_path=kernel_path, bPrint=False)
        print("XRDDebye: kernel loaded OK")

    # ── public API ────────────────────────────────────────────────────────────

    def powder_spectrum(self, atoms: np.ndarray, unique_Z: np.ndarray, Q: np.ndarray,
                        r_max: float = 20.0, dr: float = 0.02, bTryAllocate: bool = True) -> np.ndarray:
        """One-shot powder spectrum: atoms -> I(Q)."""
        hist, bin_edges = self.compute_histogram(atoms, r_max=r_max, dr=dr, bTryAllocate=bTryAllocate)
        ff_table = get_form_factor_table(unique_Z, Q)
        I = self.compute_spectrum(hist, bin_edges, Q, ff_table, unique_Z, bTryAllocate=bTryAllocate)
        return I

    def compute_histogram(self, atoms: np.ndarray, r_min: float = 0.0, r_max: float = 20.0,
                          dr: float = 0.02, bTryAllocate: bool = True) -> tuple:
        """Build pair-distance histogram on GPU.

        Args:
            atoms: (N,4) float32 [x,y,z,elem_idx].
            r_min: lower histogram edge (Å).
            r_max: upper histogram edge (Å).
            dr: bin width (Å).

        Returns:
            hist: (n_pair_types, n_bins) float32 counts.
            bin_edges: (n_bins+1,) float32.
        """
        n_atoms = atoms.shape[0]
        n_species = int(atoms[:, 3].max()) + 1
        n_bins = int(np.ceil((r_max - r_min) / dr))
        n_pair_types = n_species * (n_species + 1) // 2

        if bTryAllocate:
            buffs = {
                'atoms': 4 * 4 * n_atoms,
                'hist': 4 * n_pair_types * n_bins,
            }
            self.try_make_buffers(buffs, suffix='_cl')

        # upload atoms
        self.toGPU('atoms_cl', atoms)
        # zero histogram
        hist_zero = np.zeros((n_pair_types, n_bins), dtype=np.float32)
        self.toGPU('hist_cl', hist_zero)

        gs = self.roundUpGlobalSize(n_atoms)
        self.prg.pair_histogram(self.queue, (gs,), (self.nloc,),
                                self.buffer_dict['atoms_cl'],
                                self.buffer_dict['hist_cl'],
                                np.int32(n_atoms),
                                np.int32(n_species),
                                np.int32(n_bins),
                                np.float32(r_min),
                                np.float32(dr),
                                np.int32(n_pair_types))
        cl.enqueue_copy(self.queue, hist_zero, self.buffer_dict['hist_cl'])
        self.queue.finish()

        bin_edges = r_min + np.arange(n_bins + 1, dtype=np.float32) * dr
        return hist_zero, bin_edges

    def compute_spectrum(self, hist: np.ndarray, bin_edges: np.ndarray, Q: np.ndarray,
                         ff_table: np.ndarray, unique_Z: np.ndarray, bTryAllocate: bool = True) -> np.ndarray:
        """Debye-transform a histogram into I(Q).

        Args:
            hist: (n_pair_types, n_bins) float32 from compute_histogram.
            bin_edges: (n_bins+1,) float32.
            Q: (n_Q,) float32 scattering vectors (Å⁻¹).
            ff_table: (n_species, n_Q) float32 precomputed form factors.
            unique_Z: (n_species,) atomic numbers (for self-term counting).

        Returns:
            I: (n_Q,) float32 intensity.
        """
        n_pair_types, n_bins = hist.shape
        n_Q = len(Q)
        n_species = ff_table.shape[0]
        assert ff_table.shape[1] == n_Q
        assert n_pair_types == n_species * (n_species + 1) // 2

        # Build pair-type -> species mapping and ff_prod[pair_type, Q]
        species_a = np.empty(n_pair_types, dtype=np.int32)
        species_b = np.empty(n_pair_types, dtype=np.int32)
        idx = 0
        for a in range(n_species):
            for b in range(a, n_species):
                species_a[idx] = a
                species_b[idx] = b
                idx += 1

        ff_prod = ff_table[species_a, :] * ff_table[species_b, :]  # (n_pair_types, n_Q)

        # Self-scattering term: Σ_i f_i² = Σ_s n_s * f_s²
        # We compute this on CPU and add after GPU transform
        # (we don't have n_s here; caller should provide or we count from atoms)
        # For now skip self-term; user can add if they have atom counts.

        r_centers = bin_edges[:-1] + 0.5 * (bin_edges[1] - bin_edges[0])
        r_centers = r_centers.astype(np.float32)
        Q = Q.astype(np.float32)
        ff_prod = ff_prod.astype(np.float32)
        hist = hist.astype(np.float32)

        if bTryAllocate:
            buffs = {
                'hist': 4 * n_pair_types * n_bins,
                'ff_prod': 4 * n_pair_types * n_Q,
                'Q_vals': 4 * n_Q,
                'r_centers': 4 * n_bins,
                'I_Q': 4 * n_Q,
            }
            self.try_make_buffers(buffs, suffix='_cl')

        self.toGPU('hist_cl', hist)
        self.toGPU('ff_prod_cl', ff_prod)
        self.toGPU('Q_vals_cl', Q)
        self.toGPU('r_centers_cl', r_centers)

        gs = self.roundUpGlobalSize(n_Q)
        self.prg.debye_transform_q(self.queue, (gs,), (self.nloc,),
                                   self.buffer_dict['hist_cl'],
                                   self.buffer_dict['ff_prod_cl'],
                                   self.buffer_dict['Q_vals_cl'],
                                   self.buffer_dict['r_centers_cl'],
                                   self.buffer_dict['I_Q_cl'],
                                   np.int32(n_pair_types),
                                   np.int32(n_bins),
                                   np.int32(n_Q))
        I = np.empty(n_Q, dtype=np.float32)
        cl.enqueue_copy(self.queue, I, self.buffer_dict['I_Q_cl'])
        self.queue.finish()
        return I

    def compute_histogram_gaussian(self, atoms: np.ndarray, sigma: np.ndarray,
                                   r_min: float = 0.0, r_max: float = 20.0,
                                   dr: float = 0.02, bTryAllocate: bool = True) -> tuple:
        """Build Gaussian-broadened pair-distance histogram on GPU.

        Args:
            atoms: (N,4) float32 [x,y,z,elem_idx].
            sigma: (N,N) float32 per-pair Gaussian width (Å), or scalar broadcast.
            r_min, r_max, dr: histogram range and bin width.

        Returns:
            hist, bin_edges
        """
        n_atoms = atoms.shape[0]
        n_species = int(atoms[:, 3].max()) + 1
        n_bins = int(np.ceil((r_max - r_min) / dr))
        n_pair_types = n_species * (n_species + 1) // 2

        if np.isscalar(sigma):
            sigma_arr = np.full((n_atoms, n_atoms), float(sigma), dtype=np.float32)
        else:
            sigma_arr = np.asarray(sigma, dtype=np.float32)
            if sigma_arr.shape != (n_atoms, n_atoms):
                raise ValueError(f"sigma shape {sigma_arr.shape} != ({n_atoms},{n_atoms})")

        if bTryAllocate:
            buffs = {
                'atoms': 4 * 4 * n_atoms,
                'sigma': 4 * n_atoms * n_atoms,
                'hist': 4 * n_pair_types * n_bins,
            }
            self.try_make_buffers(buffs, suffix='_cl')

        self.toGPU('atoms_cl', atoms)
        self.toGPU('sigma_cl', sigma_arr)
        hist_zero = np.zeros((n_pair_types, n_bins), dtype=np.float32)
        self.toGPU('hist_cl', hist_zero)

        gs = self.roundUpGlobalSize(n_atoms)
        self.prg.pair_histogram_gaussian(self.queue, (gs,), (self.nloc,),
                                         self.buffer_dict['atoms_cl'],
                                         self.buffer_dict['sigma_cl'],
                                         self.buffer_dict['hist_cl'],
                                         np.int32(n_atoms),
                                         np.int32(n_species),
                                         np.int32(n_bins),
                                         np.float32(r_min),
                                         np.float32(dr),
                                         np.int32(n_pair_types),
                                         np.int32(0))
        cl.enqueue_copy(self.queue, hist_zero, self.buffer_dict['hist_cl'])
        self.queue.finish()

        bin_edges = r_min + np.arange(n_bins + 1, dtype=np.float32) * dr
        return hist_zero, bin_edges

    def powder_spectrum_with_self(self, atoms: np.ndarray, unique_Z: np.ndarray, Q: np.ndarray,
                                  r_max: float = 20.0, dr: float = 0.02) -> np.ndarray:
        """Compute I(Q) including self-scattering term Σ_i f_i²."""
        n_species = len(unique_Z)
        elem_idx = atoms[:, 3].astype(np.int32)
        counts = np.bincount(elem_idx, minlength=n_species).astype(np.float32)

        ff_table = get_form_factor_table(unique_Z, Q)
        I_pair = self.powder_spectrum(atoms, unique_Z, Q, r_max=r_max, dr=dr)

        # Self term
        I_self = np.sum(counts[:, None] * (ff_table ** 2), axis=0)
        return I_pair + I_self

    def ensemble_histogram(self, atoms_list, r_min=0.0, r_max=20.0, dr=0.02):
        """Accumulate pair-distance histograms from multiple crystal structures.

        All structures must have the same n_species (element set).
        Histograms are summed in float64 on host for numerical stability.

        Args:
            atoms_list: list of (N_i, 4) float32 arrays [x,y,z,elem_idx].
                Each can have different N.
            r_min, r_max, dr: histogram grid.

        Returns:
            hist_sum: (n_pair_types, n_bins) float64 accumulated counts.
            bin_edges: (n_bins+1,) float32.
            total_atoms: int, sum of all N_i.
        """
        n_species = 0
        for atoms in atoms_list:
            ns = int(atoms[:, 3].max()) + 1
            n_species = max(n_species, ns)
        n_bins = int(np.ceil((r_max - r_min) / dr))
        n_pair_types = n_species * (n_species + 1) // 2
        hist_sum = np.zeros((n_pair_types, n_bins), dtype=np.float64)
        total_atoms = 0
        for i, atoms in enumerate(atoms_list):
            hist_i, bin_edges = self.compute_histogram(atoms, r_min=r_min, r_max=r_max, dr=dr)
            hist_sum += hist_i.astype(np.float64)
            total_atoms += atoms.shape[0]
        return hist_sum, bin_edges, total_atoms

    def ensemble_histogram_gaussian(self, atoms_list, sigma_list, r_min=0.0, r_max=20.0, dr=0.02):
        """Accumulate Gaussian-broadened histograms from multiple structures.

        Args:
            atoms_list: list of (N_i, 4) float32 arrays.
            sigma_list: list of sigma values — either scalar or (N_i, N_i) per-pair.
            r_min, r_max, dr: histogram grid.

        Returns:
            hist_sum: (n_pair_types, n_bins) float64.
            bin_edges: (n_bins+1,) float32.
            total_atoms: int.
        """
        n_species = 0
        for atoms in atoms_list:
            ns = int(atoms[:, 3].max()) + 1
            n_species = max(n_species, ns)
        n_bins = int(np.ceil((r_max - r_min) / dr))
        n_pair_types = n_species * (n_species + 1) // 2
        hist_sum = np.zeros((n_pair_types, n_bins), dtype=np.float64)
        total_atoms = 0
        for i, (atoms, sigma) in enumerate(zip(atoms_list, sigma_list)):
            hist_i, bin_edges = self.compute_histogram_gaussian(atoms, sigma, r_min=r_min, r_max=r_max, dr=dr)
            hist_sum += hist_i.astype(np.float64)
            total_atoms += atoms.shape[0]
        return hist_sum, bin_edges, total_atoms

    def ensemble_spectrum(self, hist_sum, bin_edges, Q, unique_Z, total_atoms=None):
        """Debye-transform an accumulated ensemble histogram.

        The self-scattering term uses total_atoms * <f²> if total_atoms is provided.

        Args:
            hist_sum: (n_pair_types, n_bins) float64 or float32 accumulated histogram.
            bin_edges: (n_bins+1,) float32.
            Q: (n_Q,) float32.
            unique_Z: (n_species,) atomic numbers.
            total_atoms: total atom count across ensemble (for self-term).

        Returns:
            I: (n_Q,) float32 intensity.
        """
        ff_table = get_form_factor_table(unique_Z, Q)
        I_pair = self.compute_spectrum(hist_sum.astype(np.float32), bin_edges, Q, ff_table, unique_Z)
        if total_atoms is not None:
            n_species = len(unique_Z)
            # Approximate self-term: total_atoms * mean(f²) — exact if all same species
            # For mixed species, we'd need per-structure counts; caller can add exact self-term
            I_self = total_atoms * np.mean(ff_table ** 2, axis=0) if n_species == 1 else np.zeros(len(Q), dtype=np.float32)
            return I_pair + I_self
        return I_pair


def compute_sigma_from_sparse_blocks(pos: np.ndarray, neigh_idx: np.ndarray,
                                     neigh_count: np.ndarray, blocks: np.ndarray,
                                     kBT: float = 0.02585, default_sigma: float = 0.02) -> np.ndarray:
    """Compute per-pair thermal broadening σ_ij from local Hessian blocks.

    Uses the frozen-stiffness approximation:
        stiffness_ij = u^T (H_ii + H_jj - 2*H_ij) u
        σ_ij = sqrt(kBT / stiffness_ij)

    For pairs without a sparse block (beyond neighbor cutoff), returns default_sigma.

    Args:
        pos: (N,3) float64 atomic positions.
        neigh_idx, neigh_count, blocks: sparse Hessian from MMFF.getHessianSparseBlocks.
        kBT: thermal energy in eV (default 0.02585 eV ≈ 300 K).
        default_sigma: fallback for distant pairs without sparse block.

    Returns:
        sigma: (N,N) float32 symmetric matrix.
    """
    n_atoms = pos.shape[0]
    sigma = np.full((n_atoms, n_atoms), default_sigma, dtype=np.float32)

    # Precompute self-neighbor indices for each atom
    self_idx = np.empty(n_atoms, dtype=np.int32)
    for i in range(n_atoms):
        mask = neigh_idx[i, :neigh_count[i]] == i
        if not np.any(mask):
            continue
        self_idx[i] = np.where(mask)[0][0]

    # Build a lookup: for each atom i, map neighbor atom j -> block index k
    neigh_map = [{} for _ in range(n_atoms)]
    for i in range(n_atoms):
        for k in range(neigh_count[i]):
            j = neigh_idx[i, k]
            if j >= 0:
                neigh_map[i][j] = k

    for i in range(n_atoms):
        Hii = blocks[i, self_idx[i]]
        for j in range(i + 1, n_atoms):
            k_ij = neigh_map[i].get(j, -1)
            if k_ij < 0:
                continue  # no sparse block -> default_sigma already set
            k_ji = neigh_map[j].get(i, -1)
            if k_ji < 0:
                continue

            Hjj = blocks[j, self_idx[j]]
            # blocks[i,k_ij] is placed at H[j,i] in dense reconstruction,
            # so it represents H_ji = H_ij^T
            Hij = blocks[i, k_ij].T

            d = pos[j] - pos[i]
            r = np.linalg.norm(d)
            if r < 1e-8:
                continue
            u = d / r

            stiff = float(u @ Hii @ u + u @ Hjj @ u - 2.0 * (u @ Hij @ u))
            if stiff > 1e-12:
                sig = np.sqrt(kBT / stiff)
                sigma[i, j] = sig
                sigma[j, i] = sig

    return sigma


def compute_sigma_exact(pos: np.ndarray, neigh_idx: np.ndarray, neigh_count: np.ndarray,
                        blocks: np.ndarray, masses: np.ndarray = None,
                        kBT: float = 0.02585, default_sigma: float = 0.02,
                        rigid_proj: bool = True) -> np.ndarray:
    """Compute per-pair thermal broadening σ_ij using exact Hessian pseudoinverse.

    σ²_ij = kBT · b^T H_mw⁻¹ b

    where H_mw is the mass-weighted Hessian and b is the bond-direction force probe
    in mass-weighted coordinates. Rigid-body zero modes are projected out before
    inversion.

    Reuses FTIR.build_sparse_hessian_from_blocks() for Hessian assembly.

    Args:
        pos: (N,3) float64 atomic positions.
        neigh_idx, neigh_count, blocks: sparse Hessian blocks from MMFF.
        masses: (N,) float64 atomic masses in amu. If None, use C=12.
        kBT: thermal energy in eV (0.02585 eV ≈ 300 K).
        default_sigma: fallback for pairs where computation fails.
        rigid_proj: if True, project out 6 rigid-body zero modes before inversion.

    Returns:
        sigma: (N,N) float32 symmetric matrix.
    """
    from scipy.sparse import csc_matrix
    from scipy.sparse.linalg import splu
    from pyBall.FTIR import build_sparse_hessian_from_blocks

    n_atoms = pos.shape[0]
    ndof = 3 * n_atoms

    if masses is None:
        masses = np.full(n_atoms, 12.0, dtype=np.float64)

    # Reuse FTIR.build_sparse_hessian_from_blocks — do NOT reimplement
    H_sparse = build_sparse_hessian_from_blocks(neigh_idx, neigh_count, blocks, symmetrize=True)
    H = H_sparse.toarray()

    # Mass-weight: H_mw = M⁻¹ᐟ² H M⁻¹ᐟ²
    inv_sqrt_m = np.zeros(ndof, dtype=np.float64)
    for i in range(n_atoms):
        inv_sqrt_m[3*i:3*i+3] = 1.0 / np.sqrt(masses[i])
    H_mw = (inv_sqrt_m[:, None] * H) * inv_sqrt_m[None, :]

    # Project out rigid-body modes (6 zero modes: 3 translation + 3 rotation)
    if rigid_proj:
        sqrt_m = np.zeros(ndof, dtype=np.float64)
        for i in range(n_atoms):
            sqrt_m[3*i:3*i+3] = np.sqrt(masses[i])
        rb = np.zeros((ndof, 6), dtype=np.float64)
        for i in range(n_atoms):
            for k in range(3):
                rb[3*i+k, k] = sqrt_m[3*i+k]
            r = pos[i]
            rb[3*i+1, 3] = -sqrt_m[3*i+1] * r[2]; rb[3*i+2, 3] =  sqrt_m[3*i+2] * r[1]
            rb[3*i+0, 4] =  sqrt_m[3*i+0] * r[2]; rb[3*i+2, 4] = -sqrt_m[3*i+2] * r[0]
            rb[3*i+0, 5] = -sqrt_m[3*i+0] * r[1]; rb[3*i+1, 5] =  sqrt_m[3*i+1] * r[0]
        from scipy.linalg import qr
        Q, _ = qr(rb, mode='economic')
        I_minus_QQt = np.eye(ndof) - Q @ Q.T
        H_mw = I_minus_QQt @ H_mw @ I_minus_QQt

    # Regularize: add small diagonal to make non-singular
    eps = 1e-8 * np.max(np.abs(np.diag(H_mw)))
    H_mw += eps * np.eye(ndof)

    # Factorize once, then batch-solve all pair RHS vectors at once
    print(f"  compute_sigma_exact: factorizing {ndof}×{ndof} mass-weighted Hessian ...")
    lu = splu(csc_matrix(H_mw))

    # Build all pair RHS vectors: for each pair (i,j), b has u/sqrt(m_i) at i, -u/sqrt(m_j) at j
    # Vectorized: no Python loop over pairs for the solve
    inv_sqrt_m_atom = 1.0 / np.sqrt(masses)  # (N,)
    # Pair list
    tri_i, tri_j = np.triu_indices(n_atoms, k=1)
    d_ij = pos[tri_j] - pos[tri_i]  # (npairs, 3)
    r_ij = np.linalg.norm(d_ij, axis=1)  # (npairs,)
    valid = r_ij > 1e-8
    tri_i = tri_i[valid]; tri_j = tri_j[valid]
    d_ij = d_ij[valid]; r_ij = r_ij[valid]
    u_ij = d_ij / r_ij[:, None]  # (npairs, 3) unit bond directions
    npairs = len(tri_i)
    print(f"  compute_sigma_exact: {npairs} pairs, batched solve ...")

    # Build RHS matrix B (ndof, npairs) — vectorized sparse fill
    B = np.zeros((ndof, npairs), dtype=np.float64)
    pair_idx = np.arange(npairs)
    # b[3*i:3*i+3] = u / sqrt(m_i)
    B[3*tri_i + 0, pair_idx] = u_ij[:, 0] * inv_sqrt_m_atom[tri_i]
    B[3*tri_i + 1, pair_idx] = u_ij[:, 1] * inv_sqrt_m_atom[tri_i]
    B[3*tri_i + 2, pair_idx] = u_ij[:, 2] * inv_sqrt_m_atom[tri_i]
    # b[3*j:3*j+3] = -u / sqrt(m_j)
    B[3*tri_j + 0, pair_idx] -= u_ij[:, 0] * inv_sqrt_m_atom[tri_j]
    B[3*tri_j + 1, pair_idx] -= u_ij[:, 1] * inv_sqrt_m_atom[tri_j]
    B[3*tri_j + 2, pair_idx] -= u_ij[:, 2] * inv_sqrt_m_atom[tri_j]

    # Batched solve: X = H_mw⁻¹ B  (single scipy call, internally loops in C)
    X = lu.solve(B)
    # σ² = kBT * sum(B * X, axis=0)  — diagonal of B^T X
    var_ij = kBT * np.sum(B * X, axis=0)

    sigma = np.full((n_atoms, n_atoms), default_sigma, dtype=np.float32)
    pos_mask = var_ij > 0
    sig_vals = np.sqrt(var_ij[pos_mask])
    si = tri_i[pos_mask]; sj = tri_j[pos_mask]
    sigma[si, sj] = sig_vals.astype(np.float32)
    sigma[sj, si] = sig_vals.astype(np.float32)

    return sigma
