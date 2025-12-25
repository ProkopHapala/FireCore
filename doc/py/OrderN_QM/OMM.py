import numpy as np

# DEBUG/Dev notes:
# - Pure Python with explicit loops (C-style) for clarity and later porting.
# - Dense matrices used as provided by hydrogen_chain_1d (small demo systems).
# - No dictionaries; intersections are done by linear scans on sorted index arrays.


def _unique_sorted(arr):
    """Return sorted unique copy of 1D array (small sizes expected)."""
    if len(arr) == 0:
        return np.zeros((0,), dtype=int)
    arr_sorted = np.sort(np.asarray(arr, dtype=int))
    unique = [arr_sorted[0]]
    for i in range(1, len(arr_sorted)):
        if arr_sorted[i] != arr_sorted[i - 1]:
            unique.append(arr_sorted[i])
    return np.array(unique, dtype=int)


def build_orbital_masks(cell_ids, support_cells=1, wrap_ends=False):
    """
    Build masks for localized orbitals: one orbital per cell, support = cell +/- support_cells.
    Returns list of index arrays (sorted unique site indices per orbital).
    """
    cell_ids = np.asarray(cell_ids, dtype=int)
    n_atoms = len(cell_ids)
    max_cell = cell_ids.max()
    # Assume cells are contiguous starting at 0
    masks = []
    for c in range(max_cell + 1):
        members = []
        for ic in range(-support_cells, support_cells + 1):
            cc = c + ic
            if wrap_ends:
                cc = cc % (max_cell + 1)
            if cc < 0 or cc > max_cell:
                continue
            for idx in range(n_atoms):
                if cell_ids[idx] == cc:
                    members.append(idx)
        masks.append(_unique_sorted(members))
    print(f"#DEBUG build_orbital_masks n_cells={max_cell+1} support_cells={support_cells} wrap={wrap_ends}")
    return masks


def build_orbital_neighbors(masks):
    """Neighbor list: orbitals overlap if their masks intersect (non-empty)."""
    n_orb = len(masks)
    neighbors = [[] for _ in range(n_orb)]
    for i in range(n_orb):
        for j in range(n_orb):
            if i == j:
                neighbors[i].append(j)
                continue
            mi = masks[i]
            mj = masks[j]
            # linear intersection test
            pi = pj = 0
            overlap = False
            while pi < len(mi) and pj < len(mj):
                a = mi[pi]
                b = mj[pj]
                if a == b:
                    overlap = True
                    break
                if a < b:
                    pi += 1
                else:
                    pj += 1
            if overlap:
                neighbors[i].append(j)
    print(f"#DEBUG build_orbital_neighbors n_orb={n_orb}")
    return neighbors


def init_coeffs_constant(masks, value=1.0):
    """Initialize coefficients list with constant value on each mask."""
    C_list = []
    for m in masks:
        if len(m) == 0:
            C_list.append(np.zeros((0,), dtype=float))
        else:
            C_list.append(np.full(len(m), value, dtype=float))
    print(f"#DEBUG init_coeffs_constant n_orb={len(masks)} value={value}")
    return C_list


def init_coeffs_random(masks, scale=1.0, seed=None):
    """Initialize coefficients with small random values on each mask."""
    rng = np.random.default_rng(seed)
    C_list = []
    for m in masks:
        if len(m) == 0:
            C_list.append(np.zeros((0,), dtype=float))
        else:
            C_list.append(scale * rng.standard_normal(len(m)))
    print(f"#DEBUG init_coeffs_random n_orb={len(masks)} scale={scale} seed={seed}")
    return C_list


def spmv_expand_S(S, masks, C_list):
    """
    Compute S * psi for each orbital psi (defined by mask + coeffs).
    Returns list of (indices, values) arrays with possible one-shell expansion.
    Implementation: dense S, explicit loops, accumulation into dense temp then compressed.
    """
    n_basis = S.shape[0]
    n_orb = len(masks)
    Spsi_idx = []
    Spsi_val = []
    for i in range(n_orb):
        temp = np.zeros(n_basis, dtype=float)
        mask_i = masks[i]
        coeffs = C_list[i]
        for k in range(len(mask_i)):
            row = mask_i[k]
            cval = coeffs[k]
            # dense row multiply
            for col in range(n_basis):
                s_val = S[row, col]
                if s_val != 0.0:
                    temp[col] += s_val * cval
        nz = np.nonzero(temp)[0]
        Spsi_idx.append(nz.astype(int))
        Spsi_val.append(temp[nz])
    print(f"#DEBUG spmv_expand_S n_orb={n_orb} n_basis={n_basis}")
    return Spsi_idx, Spsi_val


def _lookup(indices, values, target):
    """Find value for target index in sorted indices array via linear scan; return 0.0 if absent."""
    for k in range(len(indices)):
        if indices[k] == target:
            return values[k]
        if indices[k] > target:
            return 0.0
    return 0.0


def compute_overlap_errors(C_list, masks, Spsi_idx, Spsi_val, neighbors):
    """
    Compute sparse overlap errors R_ij = <psi_i|S|psi_j> - delta_ij for neighbor pairs.
    Returns list of (i,j,error) and aggregated stats.
    """
    n_orb = len(masks)
    errors = []
    total_abs = 0.0
    max_abs = 0.0
    for i in range(n_orb):
        mi = masks[i]
        ci = C_list[i]
        for j in neighbors[i]:
            # compute dot product over mask of i intersect Spsi_j
            idx_j = Spsi_idx[j]
            val_j = Spsi_val[j]
            dot = 0.0
            for k in range(len(mi)):
                mu = mi[k]
                sj_mu = _lookup(idx_j, val_j, mu)
                if sj_mu != 0.0:
                    dot += ci[k] * sj_mu
            target = 1.0 if i == j else 0.0
            err = dot - target
            errors.append((i, j, err))
            aerr = abs(err)
            total_abs += aerr
            if aerr > max_abs:
                max_abs = aerr
    return errors, total_abs, max_abs


def _errors_to_matrix(errors, n_orb):
    """Convert sparse error list to dense matrix (small n_orb)."""
    R = np.zeros((n_orb, n_orb), dtype=float)
    for (i, j, err) in errors:
        R[i, j] = err
    return R


def coeffs_to_dense(C_list, masks, n_basis):
    """Expand ragged coeffs to dense (n_orb, n_basis) array."""
    n_orb = len(masks)
    dense = np.zeros((n_orb, n_basis), dtype=float)
    for i in range(n_orb):
        mi = masks[i]
        ci = C_list[i]
        for k in range(len(mi)):
            dense[i, mi[k]] = ci[k]
    return dense


def normalize_orbitals(C_list, masks, S):
    """Return new coefficient list normalized w.r.t. overlap S."""
    n_basis = S.shape[0]
    normed = []
    for i, (mi, ci) in enumerate(zip(masks, C_list)):
        vec = np.zeros(n_basis, dtype=float)
        for k in range(len(mi)):
            vec[mi[k]] = ci[k]
        # S-norm
        nrm2 = vec @ (S @ vec)
        if nrm2 <= 1e-16:
            normed.append(ci.copy())
            continue
        scale = 1.0 / np.sqrt(nrm2)
        normed_ci = ci * scale
        normed.append(normed_ci)
    return normed


def _debug_print_coeffs(C_list, masks, tag=""):
    """Print coefficients on finite support for debugging."""
    print(f"#DEBUG OMM coeffs {tag}")
    for i, (mi, ci) in enumerate(zip(masks, C_list)):
        line = " ".join(f"{mi[k]}:{ci[k]:+.6e}" for k in range(len(mi)))
        print(f"#DEBUG orb {i:02d} | {line}")


def jacobi_orthogonalize(
    S,
    masks,
    C_list,
    neighbors,
    n_iter=8,
    damping=0.05,
    eps=1e-8,
    inertia=1.0,
    max_step=0.05,
    verbosity=1,
):
    """
    Orthogonalize localized orbitals using Jacobi updates with finite support.
    S: overlap matrix (dense)
    masks: list of index arrays
    C_list: list of coefficient arrays (in/out)
    neighbors: neighbor list per orbital (indices)
    """
    n_orb = len(masks)
    for it in range(n_iter):
        Spsi_idx, Spsi_val = spmv_expand_S(S, masks, C_list)
        errors, total_abs, max_abs = compute_overlap_errors(C_list, masks, Spsi_idx, Spsi_val, neighbors)
        if verbosity >= 1:
            print(f"#DEBUG OMM iter={it} total_abs_err={total_abs:.6e} max_abs_err={max_abs:.6e}")
        if verbosity >= 3:
            Rmat = _errors_to_matrix(errors, n_orb)
            print("#DEBUG OMM R-matrix\n", Rmat)
            _debug_print_coeffs(C_list, masks, tag=f"iter={it}")
        new_C = [c.copy() for c in C_list] if max_step is not None or True else C_list
        for i in range(n_orb):
            mi = masks[i]
            ci = C_list[i]
            ci_new = new_C[i]
            for k in range(len(mi)):
                mu = mi[k]
                b = 0.0
                a = 0.0
                for j in neighbors[i]:
                    sj_mu = _lookup(Spsi_idx[j], Spsi_val[j], mu)
                    if sj_mu == 0.0:
                        continue
                    # find error(i,j)
                    # linear search in errors list (small)
                    err_ij = 0.0
                    for (ei, ej, ev) in errors:
                        if ei == i and ej == j:
                            err_ij = ev
                            break
                    weight = 2.0 if i == j else 1.0  # self-overlap derivative is doubled
                    b -= weight * err_ij * sj_mu
                    a += weight * sj_mu * sj_mu
                a += inertia
                if a < eps:
                    continue
                delta = damping * b / a
                if (max_step is not None) and (delta > max_step):
                    delta = max_step
                if (max_step is not None) and (delta < -max_step):
                    delta = -max_step
                ci_new[k] = ci[k] + delta
        C_list = new_C
        # Renormalize each iteration to enforce orthonormality trend
        C_list = normalize_orbitals(C_list, masks, S)
    # final diagnostics
    Spsi_idx, Spsi_val = spmv_expand_S(S, masks, C_list)
    errors, total_abs, max_abs = compute_overlap_errors(C_list, masks, Spsi_idx, Spsi_val, neighbors)
    if verbosity >= 1:
        print(f"#DEBUG OMM final total_abs_err={total_abs:.6e} max_abs_err={max_abs:.6e}")
    if verbosity >= 3:
        Rmat = _errors_to_matrix(errors, n_orb)
        print("#DEBUG OMM R-matrix\n", Rmat)
        _debug_print_coeffs(C_list, masks, tag="final")
    return C_list, errors, total_abs, max_abs


def run_omm(
    S,
    cell_ids,
    support_cells=1,
    wrap_ends=False,
    n_iter=8,
    damping=0.05,
    inertia=1.0,
    max_step=0.05,
    init_mode="const",
    rand_scale=0.1,
    rand_seed=None,
    verbosity=1,
):
    """
    Convenience driver: build masks, neighbors, init coefficients, and orthogonalize.
    Returns (C_list, errors, total_abs, max_abs, masks).
    """
    masks = build_orbital_masks(cell_ids, support_cells=support_cells, wrap_ends=wrap_ends)
    neighbors = build_orbital_neighbors(masks)
    if init_mode == "const":
        C_list = init_coeffs_constant(masks, value=1.0)
    elif init_mode == "rand":
        C_list = init_coeffs_random(masks, scale=rand_scale, seed=rand_seed)
    else:
        raise ValueError(f"Unknown init_mode {init_mode}")
    C_list, errors, total_abs, max_abs = jacobi_orthogonalize(
        S,
        masks,
        C_list,
        neighbors,
        n_iter=n_iter,
        damping=damping,
        inertia=inertia,
        max_step=max_step,
        verbosity=verbosity,
    )
    return C_list, errors, total_abs, max_abs, masks
