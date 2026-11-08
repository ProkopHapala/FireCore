import numpy as np

# Sparse matrix support (optional - only imported when needed)
_HAS_SCIPY = False
try:
    import scipy.sparse as sp
    from scipy.sparse.linalg import spsolve
    _HAS_SCIPY = True
except ImportError:
    pass

def dynamic_stiffness(K, M, omega, eta=1e-3, stabilize=1e-6):
    """
    Dynamic stiffness: A = K - (omega + i*eta)^2 M.
    """
    z = omega + 1j * eta
    A = K - (z * z) * M
    if stabilize > 0:
        A = A + stabilize * np.eye(K.shape[0])
    return A

def solve_response(K, M, omega, eta, direction_vec, charges, dim=3, stabilize=1e-6):
    """Solve for displacement under dipole/custom force at a single omega."""
    ndof = K.shape[0]
    n_nodes = ndof // dim
    A = dynamic_stiffness(K, M, omega, eta=eta, stabilize=stabilize)
    rhs = np.zeros((ndof, 1), dtype=np.complex128)
    for n in range(n_nodes):
        rhs[n * dim : n * dim + dim, 0] = charges[n] * direction_vec[n] if direction_vec.ndim == 2 else charges[n] * direction_vec[:dim]
    try:
        U = np.linalg.solve(A, rhs)
        return U[:, 0].reshape(n_nodes, dim)
    except np.linalg.LinAlgError:
        return np.zeros((n_nodes, dim), dtype=np.complex128)

def mechanical_greens_probing(K, M, omegas, eta=1e-3, direction_vec=None, charges=None, dim=3, stabilize=1e-6):
    """
    Dipole-driven probing of mechanical Green's function.
    Returns dict with spectra and dipole couplings.
    """
    ndof = K.shape[0]
    if ndof % dim != 0:
        raise ValueError("DOF count not divisible by dim")
    n_nodes = ndof // dim
    
    if direction_vec is None:
        direction_vec = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    direction_vec = np.asarray(direction_vec, dtype=np.float64)
    
    if charges is None:
        charges = np.ones(n_nodes)
    charges = np.asarray(charges)
    if charges.shape[0] != n_nodes:
        raise ValueError("charges size mismatch")

    spectrum_energy = np.zeros(len(omegas))
    spectrum_dipole = np.zeros((len(omegas), dim), dtype=np.complex128)

    for io, omega in enumerate(omegas):
        A = dynamic_stiffness(K, M, omega, eta=eta, stabilize=stabilize)
        # Dipole RHS: f_i = q_i * e_dir
        rhs = np.zeros((ndof, 1), dtype=np.complex128)
        for n in range(n_nodes):
            dir_v = direction_vec[n] if direction_vec.ndim == 2 else direction_vec[:dim]
            rhs[n * dim : n * dim + dim, 0] = charges[n] * dir_v
        
        try:
            U = np.linalg.solve(A, rhs)
            spectrum_energy[io] = np.sum(np.abs(U) ** 2) / n_nodes
            disp_nodes = U[:, 0].reshape(n_nodes, dim)
            dip = (charges[:, None] * disp_nodes).sum(axis=0)
            spectrum_dipole[io] = dip
        except np.linalg.LinAlgError:
            print(f"Warning: Singular matrix at omega={omega}")
            spectrum_energy[io] = 0.0
            spectrum_dipole[io] = np.zeros(dim)

    return {
        "omega": np.asarray(omegas),
        "energy": spectrum_energy,
        "dipole": spectrum_dipole,
        "n_probes": len(omegas),
    }


def build_sparse_hessian_from_blocks(neigh_idx, neigh_count, blocks, symmetrize=True):
    """
    Build scipy.sparse BSR matrix from sparse Hessian blocks.

    Parameters:
        neigh_idx: (natoms, max_neigh) neighbor indices (-1 for invalid)
        neigh_count: (natoms,) actual neighbor count per atom
        blocks: (natoms, max_neigh, 3, 3) force constant blocks
        symmetrize: ensure H is symmetric (default True)

    Returns:
        scipy.sparse.bsr_matrix of shape (3*natoms, 3*natoms)
    """
    if not _HAS_SCIPY:
        raise ImportError("scipy is required for sparse matrix operations")
    natoms, max_neigh = neigh_idx.shape
    dim = 3
    ndof = natoms * dim

    # Collect data for BSR format
    row_blocks = []
    col_blocks = []
    data_blocks = []

    for p in range(natoms):
        n = int(neigh_count[p])
        for j in range(n):
            o = int(neigh_idx[p, j])
            if o < 0:
                continue
            # Block H[o,p] = -dF_o/du_p
            blk = blocks[p, j]
            row_blocks.append(o)
            col_blocks.append(p)
            data_blocks.append(blk)
            # Symmetrize: also add H[p,o] = H[o,p].T
            if symmetrize and o != p:
                row_blocks.append(p)
                col_blocks.append(o)
                data_blocks.append(blk.T)

    # Build BSR matrix
    data = np.array(data_blocks)  # (nnz, 3, 3)
    indices = np.array(col_blocks, dtype=np.int32)  # BSR uses per-block column indices
    indptr = np.zeros(natoms + 1, dtype=np.int32)

    # Build indptr from row block counts
    row_counts = np.zeros(natoms, dtype=np.int32)
    for r in row_blocks:
        row_counts[r] += 1
    indptr[1:] = np.cumsum(row_counts)

    # Sort blocks by row for BSR format
    order = np.argsort(row_blocks, kind='stable')
    data = data[order]
    indices = indices[order]

    H_sparse = sp.bsr_matrix((data, indices, indptr), shape=(ndof, ndof), blocksize=(3, 3))
    return H_sparse


def harmonic_stick_hessian_blocks(ri, rj, k, l0):
    """3×3 Hessian blocks for E=½k(|rj−ri|−l0)² (MMFFL linear stick / K₁₂/K₁₃/K₁₄)."""
    d = np.asarray(rj, dtype=np.float64) - np.asarray(ri, dtype=np.float64)
    r = float(np.linalg.norm(d))
    if r < 1e-12:
        raise ValueError(f"harmonic_stick_hessian_blocks: zero-length stick ri={ri} rj={rj}")
    n = d / r
    I = np.eye(3)
    nn = np.outer(n, n)
    fac = k * (nn + (r - l0) / r * (I - nn))
    return fac, -fac


def build_hessian_from_linear_topology(pos, neigh_idx, bond_k, bond_l0, stick_class=None, neigh_count=None, symmetrize=True):
    """Dense 3N×3N Hessian from exported 03_topology.npz springs (bond/angle/dihedral sticks)."""
    pos = np.asarray(pos, dtype=np.float64)
    natoms = int(pos.shape[0])
    neigh_idx = np.asarray(neigh_idx, dtype=np.int32)
    bond_k = np.asarray(bond_k, dtype=np.float64)
    bond_l0 = np.asarray(bond_l0, dtype=np.float64)
    if neigh_idx.shape != bond_k.shape or neigh_idx.shape != bond_l0.shape:
        raise ValueError(f"build_hessian_from_linear_topology: shape mismatch neigh={neigh_idx.shape} k={bond_k.shape} l0={bond_l0.shape}")
    N, M = neigh_idx.shape
    if natoms != N:
        raise ValueError(f"build_hessian_from_linear_topology: pos natoms={natoms} neigh_idx rows={N}")
    sc = None if stick_class is None else np.asarray(stick_class, dtype=np.int32)
    nc = None if neigh_count is None else np.asarray(neigh_count, dtype=np.int32).reshape(-1)
    ndof = 3 * natoms
    H = np.zeros((ndof, ndof), dtype=np.float64)
    seen = set()
    for i in range(N):
        m = int(nc[i]) if nc is not None else M
        for s in range(m):
            j = int(neigh_idx[i, s])
            if j < 0 or j <= i:
                continue
            if sc is not None and int(sc[i, s]) == 0:
                continue
            key = (i, j)
            if key in seen:
                continue
            seen.add(key)
            kij, l0ij = float(bond_k[i, s]), float(bond_l0[i, s])
            Hii, Hij = harmonic_stick_hessian_blocks(pos[i], pos[j], kij, l0ij)
            ia, ib = i * 3, j * 3
            H[ia:ia + 3, ia:ia + 3] += Hii
            H[ib:ib + 3, ib:ib + 3] += Hii
            H[ia:ia + 3, ib:ib + 3] += Hij
            H[ib:ib + 3, ia:ia + 3] += Hij
    if symmetrize:
        H = 0.5 * (H + H.T)
    if np.isnan(H).any() or np.isinf(H).any():
        raise ValueError("build_hessian_from_linear_topology: NaN/Inf in Hessian")
    return H


def mass_matrix_from_Z(Z):
    """Diagonal 3N×3N mass matrix from atomic numbers (amu)."""
    from pyBall import elements
    from pyBall.io.crystal_npz import enames_from_Z
    Z = np.asarray(Z, dtype=np.int32).reshape(-1)
    natoms = Z.shape[0]
    masses = []
    for sym in enames_from_Z(Z):
        if sym not in elements.ELEMENT_DICT:
            raise ValueError(f"mass_matrix_from_Z: unknown element {sym}")
        masses.append(float(elements.ELEMENT_DICT[sym][10]))
    M = np.zeros(natoms * 3, dtype=np.float64)
    for i, m in enumerate(masses):
        M[i * 3:i * 3 + 3] = m
    return np.diag(M)


def rigid_mode_shift_vectors(M, pos):
    """Return 6 rigid-body shift vectors v_j = M^{1/2} q_j (length 3N) for Hessian projection."""
    n_atoms = pos.shape[0]
    dim = 3 * n_atoms
    m = np.diag(M).reshape(n_atoms, 3)[:, 0]
    total_mass = np.sum(m)
    com = np.sum(pos * m[:, None], axis=0) / total_mass
    pos_c = pos - com
    V = np.zeros((6, dim))
    for i in range(n_atoms):
        V[0, i*3 + 0] = 1.0; V[1, i*3 + 1] = 1.0; V[2, i*3 + 2] = 1.0
    for i in range(n_atoms):
        x, y, z = pos_c[i]
        V[3, i*3 : i*3+3] = [0.0, -z, y]
        V[4, i*3 : i*3+3] = [z, 0.0, -x]
        V[5, i*3 : i*3+3] = [-y, x, 0.0]
    Vm = V.copy()
    for j in range(6):
        for i in range(n_atoms):
            Vm[j, i*3 : i*3+3] *= np.sqrt(m[i])
    Q, _ = np.linalg.qr(Vm.T)
    Q = Q.T
    sqrt_M = np.sqrt(M)
    return [sqrt_M @ Q[j] for j in range(6)]


def apply_rigid_mode_shift(H, M, pos, shift=1e6):
    """Add rank-6 rigid-body penalty to H (dense or sparse).

    Free clusters have 6 zero modes of K: 3 translations + 3 rotations. In a frequency scan
    A(omega)=K-(omega+i*eta)^2 M, those modes cause a huge spurious response at omega->0.

    Algorithm:
      1. Build 6 rigid-body vectors in Cartesian coords (translations + rotations about COM).
      2. Mass-orthonormalize: q_j = M^{1/2} v_j, Gram-Schmidt -> Q.
      3. Add shift * (M^{1/2} q_j)(M^{1/2} q_j)^T to H for each j.

    Effect: eigenvalues of the 6 rigid directions are pushed to ~shift (very stiff), so they
    do not appear as low-frequency spectral peaks. Does NOT change internal vibrational modes.
    """
    vecs = rigid_mode_shift_vectors(M, pos)
    H_out = H.copy() if not sp.issparse(H) else H.tocsr(copy=True)
    for v in vecs:
        upd = shift * np.outer(v, v)
        H_out = H_out + (sp.csr_matrix(upd) if sp.issparse(H_out) else upd)
    return H_out


def prepare_sparse_hessian(neigh_idx, neigh_count, blocks, M, pos, shift=1e6, symmetrize=True):
    """Build scipy BSR Hessian from MMFF blocks and apply rigid-mode projection."""
    H_sparse = build_sparse_hessian_from_blocks(neigh_idx, neigh_count, blocks, symmetrize=symmetrize)
    if shift > 0:
        H_sparse = apply_rigid_mode_shift(H_sparse, M, pos, shift=shift)
    return H_sparse


def _build_probe_rhs(n_nodes, direction_vec, charges, dim=3):
    if direction_vec is None:
        direction_vec = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    direction_vec = np.asarray(direction_vec, dtype=np.float64)
    if charges is None:
        charges = np.ones(n_nodes)
    charges = np.asarray(charges)
    f = np.zeros(n_nodes * dim, dtype=np.float64)
    for n in range(n_nodes):
        dir_v = direction_vec[n] if direction_vec.ndim == 2 else direction_vec[:dim]
        f[n*dim:(n+1)*dim] = charges[n] * dir_v
    return f, charges


def _mass_diag_csr(M, n_nodes, dim=3):
    m = np.diag(M).reshape(n_nodes, dim)[:, 0]
    return sp.diags(np.repeat(m, dim), format='csr')


def _dynamic_stiffness_sparse(H_sparse, M, omega, eta=1e-3, stabilize=1e-6, dim=3):
    n_nodes = H_sparse.shape[0] // dim
    z = omega + 1j * eta
    A = H_sparse + (-(z * z)) * _mass_diag_csr(M, n_nodes, dim=dim)
    if stabilize > 0:
        A = A + stabilize * sp.eye(H_sparse.shape[0], format='csr')
    return A.tocsr()


def _momentum_bmix(k, max_iters, b_start=0.2, b_end=0.2, b_istart=3, b_iend=None, b_last=0.0):
    if b_iend is None:
        b_iend = max(b_istart, max_iters - 1)
    if k == 0 or k >= max_iters - 1:
        return b_last
    if k < b_istart:
        return 0.0
    if k >= b_iend:
        return b_end
    span = b_iend - b_istart
    if span <= 0:
        return b_end
    return b_start + (k - b_istart) / span * (b_end - b_start)


def pack_block_rows_from_sparse(A_csr, n_atoms, max_neigh=64):
    """Pack per-atom row blocks from CSR dynamic stiffness for block Jacobi."""
    A_dense = A_csr.toarray()
    row_neigh = np.full((n_atoms, max_neigh), -1, dtype=np.int32)
    row_blk = np.zeros((n_atoms, max_neigh, 3, 3), dtype=np.complex128)
    row_count = np.zeros(n_atoms, dtype=np.int32)
    diag_inv = np.zeros((n_atoms, 3, 3), dtype=np.complex128)
    for i in range(n_atoms):
        si = slice(i * 3, (i + 1) * 3)
        Aii = A_dense[si, si]
        diag_inv[i] = np.linalg.inv(Aii)
        cnt = 0
        for j in range(n_atoms):
            if i == j:
                continue
            sj = slice(j * 3, (j + 1) * 3)
            blk = A_dense[si, sj]
            if np.max(np.abs(blk)) < 1e-14:
                continue
            if cnt >= max_neigh:
                raise ValueError(f"atom {i} exceeds max_neigh={max_neigh}")
            row_neigh[i, cnt] = j
            row_blk[i, cnt] = blk
            cnt += 1
        row_count[i] = cnt
    return row_neigh, row_blk, row_count, diag_inv


def solve_block_jacobi_momentum(A, b, max_iter=500, tol=1e-8, b_start=0.2, b_end=0.2, b_istart=3, verbose=False):
    """Block Jacobi (3x3 atom blocks) + heavy-ball for complex sparse/dense A."""
    if not _HAS_SCIPY:
        raise ImportError("scipy is required for block Jacobi solver")
    b = np.asarray(b, dtype=np.complex128).ravel()
    n_atoms = b.size // 3
    A_csr = A.tocsr() if sp.issparse(A) else sp.csr_matrix(A)
    row_neigh, row_blk, row_count, diag_inv = pack_block_rows_from_sparse(A_csr, n_atoms)
    u = np.zeros_like(b)
    u_old = np.zeros_like(b)
    for k in range(max_iter):
        u_new = np.zeros_like(b)
        for i in range(n_atoms):
            si = slice(i * 3, (i + 1) * 3)
            acc = b[si].copy()
            for t in range(int(row_count[i])):
                j = int(row_neigh[i, t])
                sj = slice(j * 3, (j + 1) * 3)
                acc -= row_blk[i, t] @ u[sj]
            u_new[si] = diag_inv[i] @ acc
        bmix = _momentum_bmix(k, max_iter, b_start=b_start, b_end=b_end, b_istart=b_istart)
        u_out = u_new + bmix * (u - u_old)
        res = np.linalg.norm(b - A_csr @ u_out) / max(np.linalg.norm(b), 1e-12)
        if verbose and (k < 5 or k % 50 == 0 or res < tol):
            print(f"  block-jacobi iter {k}: rel_res={res:.3e} bmix={bmix:.3f}")
        if res < tol:
            return u_out, k + 1, res
        u_old, u = u, u_out
    raise np.linalg.LinAlgError(f"Block Jacobi did not converge in {max_iter} iterations (rel_res={res:.3e})")


def block_diag_dominance_ratio(A_csr, n_atoms):
    """Mean ratio |A_ii| / sum_{j!=i}|A_ij| per atom block (scalar diagnostic)."""
    A_dense = A_csr.toarray() if sp.issparse(A_csr) else np.asarray(A_csr)
    ratios = []
    for i in range(n_atoms):
        si = slice(i * 3, (i + 1) * 3)
        d_norm = np.linalg.norm(A_dense[si, si])
        off = sum(np.linalg.norm(A_dense[si, j * 3:(j + 1) * 3]) for j in range(n_atoms) if j != i)
        if off > 1e-30:
            ratios.append(d_norm / off)
    return float(np.mean(ratios)) if ratios else np.inf


def inertial_shift_matrix(M, n_atoms, mu, dim=3):
    """PD-style inertial shift mu*M. A' = A + mu*M improves Jacobi conditioning when mu >> 0."""
    m = np.diag(M).reshape(n_atoms, dim)[:, 0]
    return sp.diags(np.repeat(m * mu, dim), format='csr', dtype=np.complex128)


def solve_block_jacobi_inertial(A, M, b, mu=0.0, max_iter=500, tol=1e-8, b_start=0.0, b_end=0.0, b_istart=0, verbose=False):
    """Block Jacobi on (A + mu*M)u = b. mu>0 adds PD-like inertia; approximates A u=b when mu is small."""
    if not _HAS_SCIPY:
        raise ImportError("scipy is required")
    n_atoms = len(b) // 3
    A_csr = A.tocsr() if sp.issparse(A) else sp.csr_matrix(A)
    if mu > 0:
        A_csr = A_csr + inertial_shift_matrix(M, n_atoms, mu)
    return solve_block_jacobi_momentum(A_csr, b, max_iter=max_iter, tol=tol, b_start=b_start, b_end=b_end, b_istart=b_istart, verbose=verbose)


def solve_preconditioned_richardson(A, b, max_iter=2000, tol=1e-8, tau=0.1, use_block_prec=True, verbose=False):
    """Preconditioned Richardson: u <- u + tau*P^{-1}(b - A u). Needs tau < 2/lambda_max(P^{-1}A); often too slow for bonded K."""
    if not _HAS_SCIPY:
        raise ImportError("scipy is required")
    b = np.asarray(b, dtype=np.complex128).ravel()
    n_atoms = b.size // 3
    A_csr = A.tocsr() if sp.issparse(A) else sp.csr_matrix(A)
    if tau == 'auto':
        tau = estimate_richardson_tau(A_csr, n_atoms, use_block_prec=use_block_prec)
    if use_block_prec:
        _, _, _, diag_inv = pack_block_rows_from_sparse(A_csr, n_atoms)
        def apply_prec(r):
            out = np.zeros_like(r)
            for i in range(n_atoms):
                si = slice(i * 3, (i + 1) * 3)
                out[si] = diag_inv[i] @ r[si]
            return out
    else:
        dinv = 1.0 / A_csr.diagonal()
        def apply_prec(r):
            return dinv * r
    u = np.zeros_like(b)
    bnorm = max(np.linalg.norm(b), 1e-30)
    for k in range(max_iter):
        r = b - A_csr @ u
        u = u + tau * apply_prec(r)
        res = np.linalg.norm(r) / bnorm
        if verbose and (k < 5 or k % 200 == 0 or res < tol):
            print(f"  richardson iter {k}: rel_res={res:.3e} tau={tau}")
        if res < tol:
            return u, k + 1, res
    raise np.linalg.LinAlgError(f"Richardson did not converge in {max_iter} iterations (rel_res={res:.3e})")


def estimate_richardson_tau(A_csr, n_atoms, use_block_prec=True, n_power=25):
    """Estimate tau ~ 1.9/lambda_max(P^{-1}A) via power iteration."""
    n = A_csr.shape[0]
    if use_block_prec:
        _, _, _, diag_inv = pack_block_rows_from_sparse(A_csr, n_atoms)
        def apply_PinvA(x):
            out = np.zeros_like(x)
            Ax = A_csr @ x
            for i in range(n_atoms):
                si = slice(i * 3, (i + 1) * 3)
                out[si] = diag_inv[i] @ Ax[si]
            return out
    else:
        dinv = 1.0 / A_csr.diagonal()
        def apply_PinvA(x):
            return dinv * (A_csr @ x)
    v = np.random.randn(n) + 1j * np.random.randn(n)
    v = v / np.linalg.norm(v)
    for _ in range(n_power):
        v = apply_PinvA(v)
        nrm = np.linalg.norm(v)
        if nrm < 1e-30:
            break
        v = v / nrm
    lam = np.linalg.norm(A_csr @ v) / max(np.linalg.norm(apply_PinvA(A_csr @ v)), 1e-30)
    return 1.9 / max(lam, 1e-12)


def solve_cg_spd(A, b, max_iter=None, tol=1e-9, verbose=False):
    """Conjugate gradient for real symmetric positive-definite A. No diagonal-dominance requirement."""
    if not _HAS_SCIPY:
        raise ImportError("scipy is required")
    A_csr = A.tocsr() if sp.issparse(A) else sp.csr_matrix(A)
    if np.iscomplexobj(A_csr.data):
        A_csr = sp.csr_matrix(np.real(A_csr.toarray()), dtype=np.float64)
    b = np.asarray(b, dtype=np.float64).ravel()
    n = b.size
    if max_iter is None:
        max_iter = 10 * n
    x = np.zeros(n)
    r = b - A_csr @ x
    p = r.copy()
    rs = float(r @ r)
    bnorm = max(np.linalg.norm(b), 1e-30)
    for k in range(max_iter):
        Ap = A_csr @ p
        denom = float(p @ Ap)
        if abs(denom) < 1e-30:
            raise np.linalg.LinAlgError("CG breakdown: p^T A p ~ 0 (A may not be SPD)")
        alpha = rs / denom
        x += alpha * p
        r -= alpha * Ap
        rsnew = float(r @ r)
        res = np.sqrt(rsnew) / bnorm
        if verbose and (k < 5 or k % 50 == 0 or res < tol):
            print(f"  CG iter {k}: rel_res={res:.3e}")
        if res < tol:
            return x, k + 1, res
        p = r + (rsnew / rs) * p
        rs = rsnew
    raise np.linalg.LinAlgError(f"CG did not converge in {max_iter} iterations (rel_res={res:.3e})")


def solve_cgne(A, b, max_iter=None, tol=1e-10, verbose=False):
    """CG on normal equations A^H A u = A^H b for general complex A. Converges but squares condition number."""
    if not _HAS_SCIPY:
        raise ImportError("scipy is required")
    A_csr = A.tocsr() if sp.issparse(A) else sp.csr_matrix(A)
    b = np.asarray(b, dtype=np.complex128).ravel()
    n = b.size
    if max_iter is None:
        max_iter = 10 * n
    AH = A_csr.conj().T
    AHA = AH @ A_csr
    b2 = AH @ b
    x = np.zeros(n, dtype=np.complex128)
    r = b2 - AHA @ x
    p = r.copy()
    rs = np.vdot(r, r).real
    bnorm = max(np.linalg.norm(b2), 1e-30)
    for k in range(max_iter):
        Ap = AHA @ p
        denom = np.vdot(p, Ap).real
        if abs(denom) < 1e-30:
            raise np.linalg.LinAlgError("CGNE breakdown")
        alpha = rs / denom
        x += alpha * p
        r -= alpha * Ap
        rsnew = np.vdot(r, r).real
        res = np.sqrt(rsnew) / bnorm
        if verbose and (k < 5 or k % 100 == 0 or res < tol):
            print(f"  CGNE iter {k}: rel_res={res:.3e}")
        if res < tol:
            return x, k + 1, res
        p = r + (rsnew / rs) * p
        rs = rsnew
    raise np.linalg.LinAlgError(f"CGNE did not converge in {max_iter} iterations (rel_res={res:.3e})")


def solve_gmres(A, b, restart=30, max_cycles=50, tol=1e-10, verbose=False):
    """Restarted GMRES for general complex nonsymmetric A. Robust; Arnoldi orthogonalization is GPU-heavy."""
    if not _HAS_SCIPY:
        raise ImportError("scipy is required")
    A_csr = A.tocsr() if sp.issparse(A) else sp.csr_matrix(A)
    b = np.asarray(b, dtype=np.complex128).ravel()
    n = b.size
    x = np.zeros(n, dtype=np.complex128)
    bnorm = max(np.linalg.norm(b), 1e-30)
    total_iters = 0
    for cycle in range(max_cycles):
        r = b - A_csr @ x
        beta = np.linalg.norm(r)
        res = beta / bnorm
        if res < tol:
            return x, total_iters, res
        if verbose and cycle < 3:
            print(f"  GMRES cycle {cycle}: rel_res={res:.3e}")
        V = np.zeros((n, restart + 1), dtype=np.complex128)
        H = np.zeros((restart + 1, restart), dtype=np.complex128)
        V[:, 0] = r / beta
        m_eff = restart
        for j in range(restart):
            total_iters += 1
            w = A_csr @ V[:, j]
            for i in range(j + 1):
                H[i, j] = np.vdot(V[:, i], w)
                w = w - H[i, j] * V[:, i]
            H[j + 1, j] = np.linalg.norm(w)
            if abs(H[j + 1, j]) < 1e-14:
                m_eff = j + 1
                break
            V[:, j + 1] = w / H[j + 1, j]
        Hm = H[:m_eff + 1, :m_eff]
        e1 = np.zeros(m_eff + 1, dtype=np.complex128)
        e1[0] = beta
        y, _, _, _ = np.linalg.lstsq(Hm, e1, rcond=None)
        x = x + V[:, :m_eff] @ y
        r = b - A_csr @ x
        res = np.linalg.norm(r) / bnorm
        if res < tol:
            return x, total_iters, res
    raise np.linalg.LinAlgError(f"GMRES did not converge in {max_cycles} cycles (rel_res={res:.3e})")


def solve_linear_system(A, b, method='auto', M=None, omega=None, eta=0.0, **kwargs):
    """Unified entry: 'spsolve', 'cg' (real SPD), 'cgne', 'gmres', 'richardson', 'jacobi'. method='auto' picks by matrix structure."""
    if not _HAS_SCIPY:
        raise ImportError("scipy is required")
    b = np.asarray(b)
    if method == 'auto':
        A_csr = A.tocsr() if sp.issparse(A) else sp.csr_matrix(A)
        if not np.iscomplexobj(A_csr.data) and not np.iscomplexobj(b):
            w = np.linalg.eigvalsh(A_csr.toarray()) if A_csr.shape[0] <= 512 else None
            if w is not None and w[0] > 0:
                method = 'cg'
            else:
                method = 'gmres'
        else:
            method = 'gmres'
    if method == 'spsolve':
        u = spsolve(A.tocsr() if sp.issparse(A) else sp.csr_matrix(A), b.astype(np.complex128))
        return u, 1, 0.0
    if method == 'cg':
        u, nit, res = solve_cg_spd(A, b, **kwargs)
        return u, nit, res
    if method == 'cgne':
        u, nit, res = solve_cgne(A, b, **kwargs)
        return u, nit, res
    if method == 'gmres':
        u, nit, res = solve_gmres(A, b, **kwargs)
        return u, nit, res
    if method == 'richardson':
        u, nit, res = solve_preconditioned_richardson(A, b, **kwargs)
        return u, nit, res
    if method == 'jacobi':
        u, nit, res = solve_linear_jacobi_momentum(A, b, **kwargs)
        return u, nit, res
    raise ValueError(f"Unknown solver method '{method}'")


def solve_linear_jacobi_momentum(A, b, max_iter=500, tol=1e-8, b_start=0.2, b_end=0.2, b_istart=3, stabilize_diag=0.0, verbose=False, method='block'):
    """Iterative solve for complex A. method='block' (3x3 atom blocks, default) or 'point' (diagonal Jacobi)."""
    if method == 'block':
        return solve_block_jacobi_momentum(A, b, max_iter=max_iter, tol=tol, b_start=b_start, b_end=b_end, b_istart=b_istart, verbose=verbose)
    if not _HAS_SCIPY:
        raise ImportError("scipy is required for iterative sparse solvers")
    b = np.asarray(b, dtype=np.complex128).ravel()
    if sp.issparse(A):
        A_csr = A.tocsr()
        diag = A_csr.diagonal().copy()
    else:
        A_csr = sp.csr_matrix(A)
        diag = np.diag(A)
    if stabilize_diag > 0:
        diag = diag + stabilize_diag
    if np.any(np.abs(diag) < 1e-14):
        raise np.linalg.LinAlgError("Jacobi failed: near-zero diagonal entry in dynamic stiffness")
    D_inv = 1.0 / diag
    x = np.zeros_like(b)
    x_old = np.zeros_like(b)
    for k in range(max_iter):
        Ax = A_csr @ x
        x_jac = D_inv * (b - Ax + diag * x)
        bmix = _momentum_bmix(k, max_iter, b_start=b_start, b_end=b_end, b_istart=b_istart)
        x_new = x_jac + bmix * (x - x_old)
        res = np.linalg.norm(b - A_csr @ x_new) / max(np.linalg.norm(b), 1e-12)
        if verbose and (k < 5 or k % 50 == 0 or res < tol):
            print(f"  jacobi iter {k}: rel_res={res:.3e} bmix={bmix:.3f}")
        if res < tol:
            return x_new, k + 1, res
        x_old, x = x, x_new
    raise np.linalg.LinAlgError(f"Jacobi did not converge in {max_iter} iterations (rel_res={res:.3e})")


def mechanical_greens_probing_sparse(H_sparse, M, omegas, eta=1e-3, direction_vec=None, charges=None, dim=3, stabilize=1e-6,
                                     pos=None, shift_rigid=0.0, fail_loud=False, solver='spsolve', jacobi_max_iter=500, jacobi_tol=1e-8,
                                     gmres_restart=30, gmres_cycles=50):
    """
    Sparse mechanical Green's function probing at each frequency.
    solver: 'spsolve' | 'gmres' | 'cgne' | 'cg' (real SPD only) | 'richardson' | 'jacobi'
    If pos is given and shift_rigid>0, apply rigid-mode projection before solving.
    """
    if not _HAS_SCIPY:
        raise ImportError("scipy is required for sparse matrix operations")
    _valid = ('spsolve', 'jacobi', 'gmres', 'cgne', 'cg', 'richardson')
    if solver not in _valid:
        raise ValueError(f"Unknown solver '{solver}', expected one of {_valid}")

    ndof = H_sparse.shape[0]
    n_nodes = ndof // dim
    f, charges = _build_probe_rhs(n_nodes, direction_vec, charges, dim=dim)
    if pos is not None and shift_rigid > 0:
        H_sparse = apply_rigid_mode_shift(H_sparse, M, pos, shift=shift_rigid)

    spectrum_energy = np.zeros(len(omegas))
    spectrum_dipole = np.zeros((len(omegas), dim), dtype=np.complex128)

    for io, omega in enumerate(omegas):
        A_csr = _dynamic_stiffness_sparse(H_sparse, M, omega, eta=eta, stabilize=stabilize, dim=dim)
        try:
            if solver == 'spsolve':
                u = spsolve(A_csr, f.astype(np.complex128))
            elif solver == 'gmres':
                u, _, _ = solve_gmres(A_csr, f.astype(np.complex128), restart=gmres_restart, max_cycles=gmres_cycles, tol=jacobi_tol)
            elif solver == 'cgne':
                u, _, _ = solve_cgne(A_csr, f.astype(np.complex128), max_iter=jacobi_max_iter, tol=jacobi_tol)
            elif solver == 'cg':
                u, _, _ = solve_cg_spd(A_csr, f, max_iter=jacobi_max_iter, tol=jacobi_tol)
            elif solver == 'richardson':
                u, _, _ = solve_preconditioned_richardson(A_csr, f.astype(np.complex128), max_iter=jacobi_max_iter, tol=jacobi_tol, tau='auto')
            else:
                u, _, _ = solve_linear_jacobi_momentum(A_csr, f.astype(np.complex128), max_iter=jacobi_max_iter, tol=jacobi_tol)
            spectrum_energy[io] = np.sum(np.abs(u)**2) / n_nodes
            disp_nodes = u.reshape(n_nodes, dim)
            dip = (charges[:, None] * disp_nodes).sum(axis=0)
            spectrum_dipole[io] = dip
        except Exception as e:
            if fail_loud:
                raise RuntimeError(f"sparse solve failed at omega={omega}") from e
            print(f"Warning: sparse solve failed at omega={omega}: {e}")
            spectrum_energy[io] = 0.0
            spectrum_dipole[io] = np.zeros(dim)

    return {
        "omega": np.asarray(omegas),
        "energy": spectrum_energy,
        "dipole": spectrum_dipole,
        "n_probes": len(omegas),
        "solver": solver,
    }


def vibration_spectrum_from_modes(K, M, omegas, eta=1e-3, direction_vec=None, charges=None, dim=3):
    """
    Fast vibration spectrum via dense diagonalization + eigenmode summation.
    Avoids per-frequency linear solves; suitable for systems up to ~3000 DOF.

    Steps:
        1. Mass-weight: K_mw = M^{-1/2} K M^{-1/2}
        2. Diagonalize: K_mw v_i = w_i v_i
        3. Frequencies: omega_i = sqrt(w_i)
        4. Response at omega: sum_i (v_i^T f)^2 / (omega_i^2 - omega^2 - i*eta*omega)
    """
    ndof = K.shape[0]
    n_nodes = ndof // dim
    if direction_vec is None:
        direction_vec = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    direction_vec = np.asarray(direction_vec, dtype=np.float64)
    if charges is None:
        charges = np.ones(n_nodes)
    charges = np.asarray(charges)

    # Build RHS force vector f[n*dim:(n+1)*dim] = charges[n] * direction_vec
    f = np.zeros(ndof, dtype=np.float64)
    for n in range(n_nodes):
        dir_v = direction_vec[n] if direction_vec.ndim == 2 else direction_vec[:dim]
        f[n*dim:(n+1)*dim] = charges[n] * dir_v

    # Mass-weighted Hessian
    m = np.diag(M).reshape(n_nodes, dim)[:, 0]
    m_sqrt = np.repeat(np.sqrt(m), dim)
    m_inv_sqrt = 1.0 / m_sqrt
    K_mw = (m_inv_sqrt[:, None] * K) * m_inv_sqrt[None, :]

    # Diagonalize
    w, v = np.linalg.eigh(K_mw)
    from pyBall.FFfit_utils import assert_harmonic_spectrum_at_minimum
    om_cm = np.sign(w) * 521.5 * np.sqrt(np.abs(w))
    assert_harmonic_spectrum_at_minimum(om_cm, ctx="vibration_spectrum_from_modes: ")
    f_tilde = v.T @ (f * m_inv_sqrt)
    omegas_modes = np.sqrt(np.maximum(w, 0.0))

    spectrum_energy = np.zeros(len(omegas), dtype=np.float64)
    spectrum_dipole = np.zeros((len(omegas), dim), dtype=np.complex128)

    for io, omega in enumerate(omegas):
        denom = omegas_modes**2 - omega**2 - 1j * eta * omega
        # Avoid division by zero at exact resonance
        denom = np.where(np.abs(denom) < 1e-12, 1e-12, denom)
        coeffs = f_tilde / denom
        # Displacement in mass-weighted coords
        u_mw = v @ coeffs
        # Un-mass-weight
        u = u_mw * m_inv_sqrt
        spectrum_energy[io] = np.sum(np.abs(u)**2) / n_nodes
        disp_nodes = u.reshape(n_nodes, dim)
        dip = (charges[:, None] * disp_nodes).sum(axis=0)
        spectrum_dipole[io] = dip

    return {
        "omega": np.asarray(omegas),
        "energy": spectrum_energy,
        "dipole": spectrum_dipole,
        "omegas_modes": omegas_modes,
        "n_modes": len(omegas_modes),
    }


def project_rigid_modes(H, M, pos, shift=1e6):
    """
    Project out the 6 rigid body modes (translation + rotation) from the Hessian
    by shifting their eigenvalues to a high frequency, preventing response at low energy.
    """
    return apply_rigid_mode_shift(H, M, pos, shift=shift)

def get_mass_matrix(mmff, n_atoms):
    """
    Create a diagonal mass matrix by querying atom types from MMFF.
    """
    from pyBall import elements
    masses = []
    for i in range(n_atoms):
        type_name = mmff.getTypeName(i)
        # Extract letters for the element symbol (e.g. 'C_3' -> 'C', 'H_CH3' -> 'H')
        elem_str = type_name.split('_')[0]
        # Sometimes typeName can be tricky, fallback just in case
        if len(elem_str) == 0:
            elem_str = "C"
            
        if elem_str in elements.ELEMENT_DICT:
            masses.append(elements.ELEMENT_DICT[elem_str][10])
        elif elem_str.capitalize() in elements.ELEMENT_DICT:
            masses.append(elements.ELEMENT_DICT[elem_str.capitalize()][10])
        else:
            print(f"Warning: Unknown element type '{elem_str}' for atom {i}, defaulting to mass 12.0")
            masses.append(12.0)
            
    M = np.zeros(n_atoms*3)
    for i, m in enumerate(masses):
        M[i*3:i*3+3] = m
    return np.diag(M)

def save_xyz_vib(fname, elements, pos, vecs, comments=None):
    """
    Save atoms and vibration vectors into an XYZ file readable by Jmol.
    
    elements: list of strings (e.g., ['C', 'H', ...])
    pos: (N, 3) array of atomic positions
    vecs: (N, 3) array of displacement vectors, or a list/array of shape (N_frames, N, 3)
    comments: list of comment strings for each frame
    """
    if vecs.ndim == 2:
        vecs = [vecs]
    if comments is None:
        comments = [f"Frame {i}" for i in range(len(vecs))]
    elif isinstance(comments, str):
        comments = [comments]
        
    n_atoms = len(elements)
    with open(fname, 'w') as f:
        for i, vec in enumerate(vecs):
            f.write(f"{n_atoms}\n")
            f.write(f"{comments[i]}\n")
            for j in range(n_atoms):
                f.write(f"{elements[j]:3s} {pos[j,0]:10.5f} {pos[j,1]:10.5f} {pos[j,2]:10.5f} "
                        f"{vec[j,0]:10.5f} {vec[j,1]:10.5f} {vec[j,2]:10.5f}\n")


# =========================== Hessian Parameter Fitting

def build_design_matrix_from_basis(bond_bases, bond_atoms, angle_bases, angle_atoms, natoms):
    """
    Build design matrix A for linear least squares Hessian fitting.
    
    The linear model is: H_target = sum_p k_p * B_p
    where B_p are the geometric basis matrices (computed at equilibrium geometry).
    
    Parameters:
    -----------
    bond_bases : ndarray (nbonds, 4, 3, 3)
        Basis matrices for bonds [Bii, Bij, Bji, Bjj]
    bond_atoms : ndarray (nbonds, 2)
        Atom indices [i, j] for each bond
    angle_bases : ndarray (nangles, 9, 9)
        Basis matrices for angles (9x9 block for 3 atoms)
    angle_atoms : ndarray (nangles, 3)
        Atom indices [i, j, k] for each angle
    natoms : int
        Total number of atoms
        
    Returns:
    --------
    A : ndarray (n_observations, n_params)
        Design matrix where each column is a flattened basis Hessian
    atom_pairs : list
        List of (param_type, atom_indices) for each parameter
    """
    nb = len(bond_bases) if bond_bases is not None else 0
    na = len(angle_bases) if angle_bases is not None else 0
    n_params = nb + na
    dim = natoms * 3
    
    # Upper triangle indices (including diagonal)
    # H is symmetric, so we only need unique elements
    triu_idx = np.triu_indices(dim)
    n_obs = len(triu_idx[0])
    
    A = np.zeros((n_obs, n_params))
    atom_pairs = []
    
    # Add bond contributions
    for ib in range(nb):
        B_blocks = bond_bases[ib]  # (4, 3, 3) = [Bii, Bij, Bji, Bjj]
        i, j = bond_atoms[ib]
        
        # Assemble sparse block into full 3N x 3N, then extract upper triangle
        B_full = np.zeros((dim, dim))
        # Bii at (i,i)
        B_full[i*3:(i+1)*3, i*3:(i+1)*3] += B_blocks[0]
        # Bij at (i,j)
        B_full[i*3:(i+1)*3, j*3:(j+1)*3] += B_blocks[1]
        # Bji at (j,i)
        B_full[j*3:(j+1)*3, i*3:(i+1)*3] += B_blocks[2]
        # Bjj at (j,j)
        B_full[j*3:(j+1)*3, j*3:(j+1)*3] += B_blocks[3]
        
        A[:, ib] = B_full[triu_idx]
        atom_pairs.append(('bond', (i, j)))
    
    # Add angle contributions
    for ia in range(na):
        B = angle_bases[ia]  # (9, 9) block
        i, j, k = angle_atoms[ia]
        atoms = [i, j, k]
        
        # Assemble 9x9 block into full 3N x 3N
        B_full = np.zeros((dim, dim))
        for r in range(9):
            for c in range(9):
                atom_r = atoms[r // 3]
                atom_c = atoms[c // 3]
                coord_r = r % 3
                coord_c = c % 3
                B_full[atom_r*3+coord_r, atom_c*3+coord_c] += B[r, c]
        
        A[:, nb + ia] = B_full[triu_idx]
        atom_pairs.append(('angle', (i, j, k)))
    
    return A, atom_pairs, triu_idx


def fit_hessian_parameters(H_target, bond_bases, bond_atoms, angle_bases, angle_atoms, 
                            natoms, lam=1e-6, method='ridge'):
    """
    Fit force field stiffness parameters to match target Hessian.
    
    Solves: min || sum_p k_p * B_p - H_target ||²_F + λ * ||k||²
    
    Parameters:
    -----------
    H_target : ndarray (3N, 3N)
        Target Hessian matrix (e.g., from DFT)
    bond_bases : ndarray (nbonds, 4, 3, 3)
        Bond basis matrices from getBondHessianBases()
    bond_atoms : ndarray (nbonds, 2)
        Bond atom indices from getBondAtomsHess()
    angle_bases : ndarray (nangles, 9, 9)
        Angle basis matrices from getAngleHessianBases()
    angle_atoms : ndarray (nangles, 3)
        Angle atom indices from getAngleAtomsHess()
    natoms : int
        Number of atoms
    lam : float
        Ridge regularization parameter (default 1e-6)
    method : str
        'ridge' for Tikhonov regularization, 'lstsq' for simple least squares
        
    Returns:
    --------
    result : dict with keys:
        'k_opt': optimal parameters [k_bonds..., k_angles...]
        'k_bonds': bond parameters only
        'k_angles': angle parameters only
        'H_fit': fitted Hessian = sum_p k_p * B_p
        'residual': Frobenius norm of (H_fit - H_target)
        'r2': R² coefficient of determination
        'A': design matrix
        'condition_number': condition number of A.T @ A
    """
    # Build design matrix
    A, atom_pairs, triu_idx = build_design_matrix_from_basis(
        bond_bases, bond_atoms, angle_bases, angle_atoms, natoms
    )
    
    # Target vector (upper triangle of H_target)
    b = H_target[triu_idx]
    
    n_params = A.shape[1]
    nb = len(bond_bases) if bond_bases is not None else 0
    
    # Solve linear least squares with regularization
    if method == 'ridge':
        # Tikhonov: (A.T @ A + λI) k = A.T @ b
        ATA = A.T @ A
        ATb = A.T @ b
        k_opt = np.linalg.solve(ATA + lam * np.eye(n_params), ATb)
    elif method == 'lstsq':
        k_opt, residuals, rank, s = np.linalg.lstsq(A, b, rcond=None)
    else:
        raise ValueError(f"Unknown method: {method}")
    
    # Assemble fitted Hessian
    H_fit = np.zeros_like(H_target)
    
    # Add bond contributions
    for ib in range(nb):
        k = k_opt[ib]
        B_blocks = bond_bases[ib]
        i, j = bond_atoms[ib]
        H_fit[i*3:(i+1)*3, i*3:(i+1)*3] += k * B_blocks[0]
        H_fit[i*3:(i+1)*3, j*3:(j+1)*3] += k * B_blocks[1]
        H_fit[j*3:(j+1)*3, i*3:(i+1)*3] += k * B_blocks[2]
        H_fit[j*3:(j+1)*3, j*3:(j+1)*3] += k * B_blocks[3]
    
    # Add angle contributions
    na = len(angle_bases) if angle_bases is not None else 0
    for ia in range(na):
        k = k_opt[nb + ia]
        B = angle_bases[ia]
        i, j, k_atoms = angle_atoms[ia]
        atoms = [i, j, k_atoms]
        
        for r in range(9):
            for c in range(9):
                atom_r = atoms[r // 3]
                atom_c = atoms[c // 3]
                coord_r = r % 3
                coord_c = c % 3
                H_fit[atom_r*3+coord_r, atom_c*3+coord_c] += k * B[r, c]
    
    # Compute metrics
    residual = np.linalg.norm(H_fit - H_target, 'fro')
    H_norm = np.linalg.norm(H_target, 'fro')
    r2 = 1.0 - (residual**2) / (H_norm**2) if H_norm > 0 else 0.0
    
    # Condition number for diagnostics
    cond_num = np.linalg.cond(A.T @ A) if A.shape[1] <= A.shape[0] else np.inf
    
    return {
        'k_opt': k_opt,
        'k_bonds': k_opt[:nb],
        'k_angles': k_opt[nb:] if na > 0 else np.array([]),
        'H_fit': H_fit,
        'residual': residual,
        'r2': r2,
        'A': A,
        'condition_number': cond_num,
        'atom_pairs': atom_pairs
    }


def fit_hessian_parameters_ctx(H_target, ctx, lam=1e-6, method='ridge'):
    """
    Fit force field stiffness parameters using context from getHessianContext().
    
    This is a convenience wrapper that extracts data from the context dict
    and calls fit_hessian_parameters().
    
    Parameters:
    -----------
    H_target : ndarray (3N, 3N)
        Target Hessian matrix (e.g., from DFT)
    ctx : dict
        Context dictionary from MMFF.getHessianContext()
    lam : float
        Ridge regularization parameter (default 1e-6)
    method : str
        'ridge' for Tikhonov regularization, 'lstsq' for simple least squares
        
    Returns:
    --------
    result : dict with fitting results (see fit_hessian_parameters)
    """
    return fit_hessian_parameters(
        H_target,
        ctx['bond_bases'], ctx['bond_atoms'],
        ctx['angle_bases'], ctx['angle_atoms'],
        ctx['n_atoms'], lam=lam, method=method
    )


def validate_hessian_fit(mmff, H_target, k_opt, inds, dx=1e-4):
    """
    Validate fitted parameters by computing Hessian with new parameters
    and comparing to target.
    
    Parameters:
    -----------
    mmff : MMFF module
        Initialized MMFF with UFF forcefield
    H_target : ndarray (3N, 3N)
        Target Hessian
    k_opt : ndarray
        Fitted parameters [k_bonds..., k_angles...]
    inds : ndarray
        Atom indices passed to getHessian3Nx3N
    dx : float
        Finite difference step for validation Hessian
        
    Returns:
    --------
    dict with validation metrics
    """
    # Set fitted parameters
    nb = mmff.getNBondsHess()
    na = mmff.getNAnglesHess()
    
    mmff.setBondParamsHess(k_opt[:nb])
    if na > 0:
        mmff.setAngleParamsHess(k_opt[nb:])
    
    # Compute Hessian with fitted parameters
    H_fitted = mmff.getHessian3Nx3N(inds, dx=dx)
    H_fitted = 0.5 * (H_fitted + H_fitted.T)  # Symmetrize
    
    # Compare
    diff = H_fitted - H_target
    rel_error = np.linalg.norm(diff, 'fro') / np.linalg.norm(H_target, 'fro')
    
    max_diff = np.max(np.abs(diff))
    mean_diff = np.mean(np.abs(diff))
    
    return {
        'H_fitted': H_fitted,
        'relative_error': rel_error,
        'max_absolute_diff': max_diff,
        'mean_absolute_diff': mean_diff,
        'diff_matrix': diff
    }

