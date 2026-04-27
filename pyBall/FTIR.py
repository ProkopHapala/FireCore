import numpy as np

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

