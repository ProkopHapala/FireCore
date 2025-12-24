import numpy as np
import matplotlib.pyplot as plt

from OrderN import (
    find_neighbor_pairs,
    build_matrix_from_pairs,
    extract_edges,
    plot_matrix,
)


def build_lattice_vectors(a_len=2.46, b_len=2.46, gamma_deg=60.0):
    """Return lattice vectors A, B (3D) for planar system."""
    gamma = np.deg2rad(gamma_deg)
    A = np.array([a_len, 0.0, 0.0])
    B = np.array([b_len * np.cos(gamma), b_len * np.sin(gamma), 0.0])
    print(f"#DEBUG build_lattice_vectors a={a_len} b={b_len} gamma_deg={gamma_deg} A={A} B={B}")
    return A, B


def build_unit_positions(d_atom=1.42, phi_deg=0.0):
    """Two-atom basis; phi is angle of displacement vs lattice vector A."""
    phi = np.deg2rad(phi_deg)
    disp = np.array([d_atom * np.cos(phi), d_atom * np.sin(phi), 0.0])
    atoms = [np.zeros(3), disp]
    print(f"#DEBUG build_unit_positions d_atom={d_atom} phi_deg={phi_deg} disp={disp}")
    return np.array(atoms)


def replicate_cells(A, B, unit_atoms, nx=6, ny=4):
    """Tile the unit cell positions into an (nx, ny) finite patch."""
    atoms = []
    cell_ix = []
    cell_iy = []
    for ix in range(nx):
        for iy in range(ny):
            origin = ix * A + iy * B
            for a in unit_atoms:
                atoms.append(origin + a)
                cell_ix.append(ix)
                cell_iy.append(iy)
    pos = np.array(atoms)
    cell_ix = np.array(cell_ix, dtype=int)
    cell_iy = np.array(cell_iy, dtype=int)
    print(f"#DEBUG replicate_cells nx={nx} ny={ny} n_atoms={len(atoms)}")
    return pos, cell_ix, cell_iy


def build_hamiltonian(pos, cell_ix, cell_iy, onsite=0.0, decay=1.0, t0=-1.0):
    pairs = find_neighbor_pairs(pos, cell_ix=cell_ix, cell_iy=cell_iy, tol=1.05)
    H = build_matrix_from_pairs(len(pos), pairs, diag_value=onsite, prefactor=t0, decay=decay)
    return H


def build_overlap(pos, cell_ix, cell_iy, diag=1.0, decay=1.0, s0=1.0):
    pairs = find_neighbor_pairs(pos, cell_ix=cell_ix, cell_iy=cell_iy, tol=1.05)
    S = build_matrix_from_pairs(len(pos), pairs, diag_value=diag, prefactor=s0, decay=decay)
    return S


def plot_hamiltonian(H, title="2D lattice Hamiltonian"):
    plot_matrix(H, title=title, symmetric=True)


def plot_overlap(S, title="2D overlap matrix"):
    plot_matrix(S, title=title, symmetric=True)

def plot_geometry(pos, A, B, nx, ny, edges=None, title="2D geometry"):
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(pos[:, 0], pos[:, 1], c="k", s=18, label="H")
    if edges:
        for (i, j, t) in edges:
            ax.plot([pos[i, 0], pos[j, 0]], [pos[i, 1], pos[j, 1]], "-", c="tab:blue", lw=1.5, alpha=0.8)
    # draw cell grid
    for ix in range(nx + 1):
        p0 = ix * A
        p1 = p0 + ny * B
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]], "--", c="gray", lw=0.7)
    for iy in range(ny + 1):
        p0 = iy * B
        p1 = p0 + nx * A
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]], "--", c="gray", lw=0.7)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x [Å]")
    ax.set_ylabel("y [Å]")
    ax.set_title(title)
    ax.legend(loc="upper right")
    plt.tight_layout()
    

if __name__ == "__main__":
    params = dict(
        d_atom    =1.42,       # (1) distance between the two atoms
        a_len     =2.46,       # (2) |A|
        b_len     =2.46,       # (2) |B|
        gamma_deg =60.0,       # (3) angle between A, B
        phi_deg   =90.0,        # (4) angle of basis displacement vs A
        nx        =6, 
        ny        =4,          # (5) replication counts
        onsite    =-5.0,
        decay     =1.2,
        t0        =-1.0,
        s0        =0.3,
        TS_decay_factor=1.5,
    )
    A, B = build_lattice_vectors(
        a_len=params["a_len"],
        b_len=params["b_len"],
        gamma_deg=params["gamma_deg"],
    )
    unit_atoms = build_unit_positions(d_atom=params["d_atom"], phi_deg=params["phi_deg"])
    pos, cell_ix, cell_iy = replicate_cells(A, B, unit_atoms, nx=params["nx"], ny=params["ny"])
    H = build_hamiltonian(pos, cell_ix, cell_iy, onsite=params["onsite"], decay=params["decay"], t0=params["t0"])
    S = build_overlap(pos, cell_ix, cell_iy, diag=1.0, decay=params["decay"]*params["TS_decay_factor"], s0=params["s0"])
    edges = extract_edges(H)
    print(f"#DEBUG |H| max={np.abs(H).max():.3f} min={np.abs(H).min():.3f}")
    print(f"#DEBUG |S| max={np.abs(S).max():.3f} min={np.abs(S).min():.3f}")
    plot_geometry(pos, A, B, params["nx"], params["ny"], edges=edges, title="2D graphene-like patch geometry")
    plot_hamiltonian(H, title="2D graphene-like patch (tight-binding)")
    plot_overlap(S, title="2D overlap matrix")
    plt.show()