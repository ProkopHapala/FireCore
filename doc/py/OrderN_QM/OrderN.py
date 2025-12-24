import numpy as np
import matplotlib.pyplot as plt


def hopping_from_distance(r, prefactor=-1.0, decay=1.0):
    """Generic exponential coupling."""
    return prefactor * np.exp(-r / decay)


def find_neighbor_pairs(pos, cell_ids=None, cell_ix=None, cell_iy=None, cell_delta_max=1, tol=1.05, max_distance=None):
    """
    Identify nearest-neighbor pairs limited by cell proximity.
    If max_distance is provided, use it as cutoff; otherwise r_cut = tol * r_min.
    Returns list of (i,j,r) where r <= r_cut among considered pairs.
    """
    n = pos.shape[0]
    use_1d = cell_ids is not None
    if not use_1d and (cell_ix is None or cell_iy is None):
        raise ValueError("Provide either cell_ids (1D) or both cell_ix, cell_iy (2D)")

    def within_cells(i, j):
        if use_1d:
            return abs(cell_ids[i] - cell_ids[j]) <= cell_delta_max
        dx = abs(cell_ix[i] - cell_ix[j])
        dy = abs(cell_iy[i] - cell_iy[j])
        return max(dx, dy) <= cell_delta_max

    r_min = None
    for i in range(n):
        for j in range(i + 1, n):
            if not within_cells(i, j):
                continue
            r = np.linalg.norm(pos[i] - pos[j])
            if r_min is None or r < r_min:
                r_min = r
    if r_min is None:
        print("#DEBUG find_neighbor_pairs no neighbor candidates")
        return []
    r_cut = max_distance if max_distance is not None else tol * r_min
    pairs = []
    for i in range(n):
        for j in range(i + 1, n):
            if not within_cells(i, j):
                continue
            r = np.linalg.norm(pos[i] - pos[j])
            if r <= r_cut:
                pairs.append((i, j, r))
    print(f"#DEBUG find_neighbor_pairs n={n} r_min={r_min:.3f} r_cut={r_cut:.3f} n_pairs={len(pairs)}")
    return pairs


def build_matrix_from_pairs(n, pairs, diag_value=0.0, prefactor=-1.0, decay=1.0):
    """Fill symmetric matrix from pair list (i,j,r) using exponential decay."""
    M = np.zeros((n, n))
    np.fill_diagonal(M, diag_value)
    for i, j, r in pairs:
        v = hopping_from_distance(r, prefactor=prefactor, decay=decay)
        M[i, j] = v
        M[j, i] = v
    print(f"#DEBUG build_matrix_from_pairs n={n} diag={diag_value} prefactor={prefactor} decay={decay}")
    return M


def extract_edges(M):
    """Return list of (i,j,val) for non-zero off-diagonal entries (upper triangle)."""
    i_idx, j_idx = np.nonzero(np.triu(M, k=1))
    vals = M[i_idx, j_idx]
    edges = list(zip(i_idx, j_idx, vals))
    print(f"#DEBUG extract_edges n_edges={len(edges)}")
    return edges


def plot_matrix(M, title="Matrix", symmetric=True, cmap="bwr"):
    fig, ax = plt.subplots(figsize=(7, 6))
    vmax = np.abs(M).max() if symmetric else None
    im = ax.imshow(
        M,
        cmap=cmap,
        origin="lower",
        vmin=-vmax if symmetric else None,
        vmax=vmax if symmetric else None,
    )
    ax.set_title(title)
    ax.set_xlabel("j")
    ax.set_ylabel("i")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label=title)
    plt.tight_layout()
    return fig, ax


def solve_generalized(H, S=None):
    """Solve generalized eigenproblem H C = S C e (or standard if S is None)."""
    if S is None:
        e, C = np.linalg.eigh(H)
    else:
        # symmetric orthogonalization: S = U s U^T, H' = s^-1/2 U^T H U s^-1/2
        s, U = np.linalg.eigh(S)
        eps = 1e-12
        s_inv_sqrt = np.diag(1.0 / np.sqrt(np.clip(s, eps, None)))
        X = U @ s_inv_sqrt @ U.T
        H_ortho = X @ H @ X
        e, Ctilde = np.linalg.eigh(H_ortho)
        C = X @ Ctilde
    print(f"#DEBUG solve_generalized n={H.shape[0]} use_S={S is not None}")
    return e, C


def mulliken_density(C, S=None, n_occ=None):
    """
    Mulliken population: rho_i = sum_occ sum_j C_i C_j S_ij.
    If S is None, reduces to sum_occ |C_i|^2.
    """
    n_basis = C.shape[0]
    if n_occ is None:
        n_occ = n_basis // 2
    Cocc = C[:, :n_occ]
    if S is None:
        rho = (Cocc**2).sum(axis=1)
    else:
        rho = np.sum(Cocc * (S @ Cocc), axis=1)
    print(f"#DEBUG mulliken_density n_basis={n_basis} n_occ={n_occ}")
    return rho


def canonical_reference(H, S=None, n_occ=None):
    """Convenience: solve eigenproblem and return (e, C, rho)."""
    e, C = solve_generalized(H, S)
    rho = mulliken_density(C, S, n_occ=n_occ)
    return e, C, rho


def plot_spectrum_and_map(e, C, n_occ=None, positions=None, nbins=200, padding=0.2, cmap="bwr", plot_spectrum=True, info_text=None):
    """
    Plot eigen-spectrum (occupied vs unoccupied colors + Fermi) and 2D map of MO coefficients vs energy.
    positions: optional x-axis values for sites; if None use site indices.
    Set plot_spectrum=False to skip the line spectrum figure.
    """
    n_basis = C.shape[0]
    if n_occ is None:
        n_occ = n_basis // 2
    e = np.asarray(e)
    C = np.asarray(C)
    e_f = 0.5 * (e[n_occ - 1] + e[n_occ]) if n_occ < len(e) else e[-1]
    xvals = np.arange(n_basis) if positions is None else np.asarray(positions)

    # spectrum
    fig_spec = ax_spec = None
    if plot_spectrum:
        fig_spec, ax_spec = plt.subplots(figsize=(5, 6))
        for idx, ei in enumerate(e):
            color = "tab:blue" if idx < n_occ else "tab:red"
            ax_spec.plot([-0.4, 0.4], [ei, ei], "-", c=color, lw=1.5)
        ax_spec.axhline(e_f, color="k", ls="--", lw=1, label="Fermi")
        ax_spec.set_xlabel("orbital")
        ax_spec.set_ylabel("energy")
        ax_spec.set_title("Eigen-spectrum")
        ax_spec.legend()
        plt.tight_layout()

    # energy-site map
    emin, emax = e.min() - padding, e.max() + padding
    bins = np.linspace(emin, emax, nbins)
    grid = np.zeros((nbins, n_basis))
    for idx, ei in enumerate(e):
        b = np.clip(np.searchsorted(bins, ei) - 1, 0, nbins - 1)
        grid[b, :] += C[:, idx]
    vmax = np.abs(grid).max()
    fig_map, ax_map = plt.subplots(figsize=(6, 6))
    im = ax_map.imshow(
        grid,
        origin="lower",
        aspect="auto",
        cmap=cmap,
        vmin=-vmax,
        vmax=vmax,
        extent=[xvals.min(), xvals.max(), emin, emax],
    )
    ax_map.axhline(e_f, color="k", ls="--", lw=1, label="Fermi")
    ax_map.set_ylabel("energy")
    ax_map.set_xlabel("site index" if positions is None else "x")
    ax_map.set_title("MO coefficients vs energy")
    ax_map.set_yticks(e)
    ax_map.set_yticklabels([f"{ei:.2f}" for ei in e])
    ax_map.legend(loc="upper right")
    if info_text:
        ax_map.text(0.02, 0.98, info_text, ha="left", va="top", transform=ax_map.transAxes, fontsize=9, bbox=dict(facecolor="white", alpha=0.6, edgecolor="none"))
    fig_map.colorbar(im, ax=ax_map, fraction=0.046, pad=0.04, label="C_i")
    plt.tight_layout()
    return (fig_spec, ax_spec), (fig_map, ax_map)
