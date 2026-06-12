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


def mechanical_greens_probing_sparse(H_sparse, M, omegas, eta=1e-3, direction_vec=None, charges=None, dim=3, stabilize=1e-6):
    """
    Sparse version of mechanical Green's function probing.
    Uses scipy.sparse.linalg.spsolve at each frequency.
    Suitable for medium systems (up to ~10000 DOF) where dense solve would be too slow.
    """
    if not _HAS_SCIPY:
        raise ImportError("scipy is required for sparse matrix operations")

    ndof = H_sparse.shape[0]
    n_nodes = ndof // dim

    if direction_vec is None:
        direction_vec = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    direction_vec = np.asarray(direction_vec, dtype=np.float64)

    if charges is None:
        charges = np.ones(n_nodes)
    charges = np.asarray(charges)

    # Build RHS force vector once
    f = np.zeros(ndof, dtype=np.float64)
    for n in range(n_nodes):
        dir_v = direction_vec[n] if direction_vec.ndim == 2 else direction_vec[:dim]
        f[n*dim:(n+1)*dim] = charges[n] * dir_v

    # Convert mass to vector
    m = np.diag(M).reshape(n_nodes, dim)[:, 0]
    m_repeated = np.repeat(m, dim)
    M_sparse = sp.diags(m_repeated, format='csr')

    spectrum_energy = np.zeros(len(omegas))
    spectrum_dipole = np.zeros((len(omegas), dim), dtype=np.complex128)

    for io, omega in enumerate(omegas):
        # Dynamic stiffness: A = H - (omega + i*eta)^2 M
        z = omega + 1j * eta
        coeff = -(z * z)
        A = H_sparse + coeff * M_sparse
        if stabilize > 0:
            A = A + stabilize * sp.eye(ndof, format='csr')

        try:
            # Convert to CSR for spsolve (more efficient)
            A_csr = A.tocsr()
            # Solve A u = f for u
            u = spsolve(A_csr, f.astype(np.complex128))
            spectrum_energy[io] = np.sum(np.abs(u)**2) / n_nodes
            disp_nodes = u.reshape(n_nodes, dim)
            dip = (charges[:, None] * disp_nodes).sum(axis=0)
            spectrum_dipole[io] = dip
        except Exception as e:
            print(f"Warning: sparse solve failed at omega={omega}: {e}")
            spectrum_energy[io] = 0.0
            spectrum_dipole[io] = np.zeros(dim)

    return {
        "omega": np.asarray(omegas),
        "energy": spectrum_energy,
        "dipole": spectrum_dipole,
        "n_probes": len(omegas),
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
    # Project force onto modes
    f_tilde = v.T @ (f * m_inv_sqrt)

    # Frequencies
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
    
    H: (3N, 3N) unweighted Hessian
    M: (3N, 3N) diagonal mass matrix
    pos: (N, 3) positions of atoms
    shift: eigenvalue shift (force constant penalty) applied to rigid body modes
    """
    n_atoms = pos.shape[0]
    dim = 3 * n_atoms
    
    # Get mass array
    m = np.diag(M).reshape(n_atoms, 3)[:, 0] # mass of each atom
    
    # Center of mass
    total_mass = np.sum(m)
    com = np.sum(pos * m[:, None], axis=0) / total_mass
    
    pos_c = pos - com
    
    # Construct the 6 rigid body vectors in unweighted coordinates
    V = np.zeros((6, dim))
    
    # Translations
    for i in range(n_atoms):
        V[0, i*3 + 0] = 1.0
        V[1, i*3 + 1] = 1.0
        V[2, i*3 + 2] = 1.0
        
    # Rotations: cross product (pos_c x v)
    for i in range(n_atoms):
        x, y, z = pos_c[i]
        V[3, i*3 : i*3+3] = [0.0, -z, y]
        V[4, i*3 : i*3+3] = [z, 0.0, -x]
        V[5, i*3 : i*3+3] = [-y, x, 0.0]
        
    # Mass-weight the basis vectors
    Vm = V.copy()
    for j in range(6):
        for i in range(n_atoms):
            Vm[j, i*3 : i*3+3] *= np.sqrt(m[i])
            
    # Orthonormalize the mass-weighted vectors
    Q, _ = np.linalg.qr(Vm.T)
    Q = Q.T # shape (6, 3N)
    
    sqrt_M = np.sqrt(M)
    H_shifted = H.copy()
    for j in range(6):
        # Delta H = shift * (M^{1/2} q) (M^{1/2} q)^T
        vec = sqrt_M @ Q[j]
        H_shifted += shift * np.outer(vec, vec)
        
    return H_shifted

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

