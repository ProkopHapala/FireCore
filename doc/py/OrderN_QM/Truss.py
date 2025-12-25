import numpy as np

# DEBUG: utilities for building and assembling truss systems


def node_index(ix, iy, nx):
    return iy * nx + ix


def build_triangular_grid(nx, ny, a=1.0, jitter=0.0, seed=None):
    """
    Build 2D grid of points (z=0) arranged on square lattice with optional jitter.
    Geometry is 3D-compatible (z=0) to keep solver general.
    """
    if seed is not None:
        np.random.seed(seed)
    pos = np.zeros((nx * ny, 3))
    idx = 0
    for iy in range(ny):
        for ix in range(nx):
            x = a * ix
            y = a * iy
            if jitter > 0:
                x += jitter * (np.random.rand() - 0.5)
                y += jitter * (np.random.rand() - 0.5)
            pos[idx, 0] = x
            pos[idx, 1] = y
            idx += 1
    print(f"#DEBUG build_triangular_grid nx={nx} ny={ny} a={a} jitter={jitter}")
    return pos


def grid_edges(nx, ny, include_diag=True):
    """
    Create edge list (i,j) connecting nearest neighbors on rectangular grid.
    include_diag adds both diagonals to form triangles.
    """
    edges = []
    for iy in range(ny):
        for ix in range(nx):
            i = node_index(ix, iy, nx)
            if ix + 1 < nx:
                edges.append((i, node_index(ix + 1, iy, nx)))
            if iy + 1 < ny:
                edges.append((i, node_index(ix, iy + 1, nx)))
            if include_diag and ix + 1 < nx and iy + 1 < ny:
                edges.append((i, node_index(ix + 1, iy + 1, nx)))
            if include_diag and ix - 1 >= 0 and iy + 1 < ny:
                edges.append((i, node_index(ix - 1, iy + 1, nx)))
    print(f"#DEBUG grid_edges nx={nx} ny={ny} n_edges={len(edges)} include_diag={include_diag}")
    return edges


def edge_lengths(pos, edges):
    lens = np.zeros(len(edges))
    for idx, (i, j) in enumerate(edges):
        lens[idx] = np.linalg.norm(pos[j] - pos[i])
    print(f"#DEBUG edge_lengths n_edges={len(edges)} l_min={lens.min():.4f} l_max={lens.max():.4f}")
    return lens


def assemble_stiffness_dense(pos, edges, k_spring=1.0, dim=3):
    """
    Assemble dense stiffness matrix for linear springs between nodes.
    K size is (dim*N, dim*N). Contributions are k * (u u^T).
    """
    n_nodes = pos.shape[0]
    ndof = dim * n_nodes
    K = np.zeros((ndof, ndof), dtype=np.float64)
    for (i, j) in edges:
        d = pos[j] - pos[i]
        L = np.linalg.norm(d)
        if L <= 1e-12:
            print(f"#DEBUG assemble_stiffness_dense zero length edge i={i} j={j}")
            continue
        u = d / L
        k_fac = k_spring / (L * L)
        # block indices
        ia = i * dim
        ja = j * dim
        outer = k_fac * np.outer(u, u)
        K[ia:ia+dim, ia:ia+dim] += outer
        K[ja:ja+dim, ja:ja+dim] += outer
        K[ia:ia+dim, ja:ja+dim] -= outer
        K[ja:ja+dim, ia:ia+dim] -= outer
    print(f"#DEBUG assemble_stiffness_dense ndof={ndof} n_edges={len(edges)} k={k_spring}")
    return K


def mass_matrix(masses, dim=3):
    masses = np.asarray(masses)
    n_nodes = masses.size
    diag = np.repeat(masses, dim)
    M = np.diag(diag)
    print(f"#DEBUG mass_matrix n_nodes={n_nodes} dim={dim} m_min={masses.min():.4f} m_max={masses.max():.4f}")
    return M


def boundary_nodes(nx, ny, which="bottom"):
    nodes = []
    if which == "bottom":
        nodes = [node_index(ix, 0, nx) for ix in range(nx)]
    elif which == "top":
        nodes = [node_index(ix, ny - 1, nx) for ix in range(nx)]
    elif which == "left":
        nodes = [node_index(0, iy, nx) for iy in range(ny)]
    elif which == "right":
        nodes = [node_index(nx - 1, iy, nx) for iy in range(ny)]
    elif which == "none":
        nodes = []
    else:
        raise ValueError(f"Unknown boundary selector {which}")
    print(f"#DEBUG boundary_nodes which={which} n={len(nodes)}")
    return nodes


def apply_dirichlet(K, M, fixed_nodes, dim=3):
    """
    Remove DOFs corresponding to fixed nodes. Returns reduced matrices and mask.
    """
    n_nodes = K.shape[0] // dim
    mask = np.ones(n_nodes * dim, dtype=bool)
    for n in fixed_nodes:
        start = n * dim
        mask[start:start+dim] = False
    K_red = K[np.ix_(mask, mask)]
    M_red = M[np.ix_(mask, mask)]
    print(f"#DEBUG apply_dirichlet fixed={len(fixed_nodes)} ndof_full={n_nodes*dim} ndof_red={K_red.shape[0]}")
    return K_red, M_red, mask


def build_test_truss(
    nx=6,
    ny=6,
    a=1.0,
    jitter=0.0,
    k_spring=1.0,
    mass_value=1.0,
    fixed_boundary="none",
    dim=3,
):
    pos = build_triangular_grid(nx, ny, a=a, jitter=jitter)
    edges = grid_edges(nx, ny, include_diag=True)
    K = assemble_stiffness_dense(pos, edges, k_spring=k_spring, dim=dim)
    masses = np.full(nx * ny, mass_value)
    M = mass_matrix(masses, dim=dim)
    fixed = boundary_nodes(nx, ny, which=fixed_boundary)
    if len(fixed) == 0:
        mask = np.ones(K.shape[0], dtype=bool)
        K_red, M_red = K, M
    else:
        K_red, M_red, mask = apply_dirichlet(K, M, fixed, dim=dim)
    print(f"#DEBUG build_test_truss ndof_full={K.shape[0]} ndof_red={K_red.shape[0]} fixed={len(fixed)}")
    return pos, edges, K_red, M_red, mask
