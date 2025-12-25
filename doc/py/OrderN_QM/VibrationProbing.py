"""
Mechanical Green-function probing for vibration spectra of a spring-mass truss.


Problem setup
-------------
- Nodes are mass points in 3D (we keep z=0 for planar demos).
- Springs (edges) define stiffness matrix K (assembled from axial springs).
- Masses define diagonal mass matrix M.
- After removing fixed DOFs (Dirichlet), we study frequency response to external forces.

Dynamic stiffness
-----------------
For real driving frequency ω and small damping η>0 we use the shifted operator
    A(ω) = K - (ω + i η)^2 M
Poles of A^{-1} are near the real axis; small η lifts them slightly to keep the
factorization stable while preserving sharp modal peaks.

Probing strategy (dipole-driven)
--------------------------------
We care about how a homogeneous electric field E couples to point charges q_i:
  - The field exerts force f_i = q_i E (same direction for all nodes).
  - This is exactly a dipole probe; no need for deterministic combs.
  - For each ω we form the dipole RHS once (or per Cartesian direction),
    factor A(ω) once (Cholesky), and solve A(ω) U = f.
  - We scan ω near the real axis with small η to sample the vibration spectrum
    and resolve modes that couple strongly to the dipole operator.

Dipole response with point charges
----------------------------------
Given charges q_i and displacement u_i(ω) from the dipole force,
the induced dipole is
    Δp(ω) = Σ_i q_i u_i(ω)
which is the linear-response quantity multiplying the external field E via
interaction energy -E·Δp. Peaks of |Δp(ω)| highlight vibrational modes that
are IR-active under the chosen field direction.

Outputs
-------
mechanical_greens_probing returns:
  - omega: sampled frequencies
  - energy: proxy amplitude per ω (||U||^2 average)
  - dipole: complex dipole response per ω (vector of length 3)
  - n_probes: total probes used

Usage
-----
Call demo_run() for a triangular grid example; tune nx, ny, k_spring, mass,
and η to trade resolution vs speed.

Didactic notes / how to read results
------------------------------------
- Zero/near-zero modes: A free 3D body has 6 rigid modes (3 translations, 3 rotations) ⇒ ω≈0. Our planar springs also leave z-DOFs loose, so additional tiny eigenvalues (~1e-8) appear unless you add out-of-plane stiffness or run with dim=2.
- Spectrum plot: black/gray vertical lines mark eigenfrequencies; selected ones are highlighted (red). The blue curve is |Δp(ω)| from dipole forcing; peaks align with IR-active modes.
- Mode panels: each panel overlays the eigenvector (blue) and the forced response at that eigenfrequency (orange). Agreement indicates the dynamic solve matches the eigensolver at ω=ω_mode.
- Charge pattern: corner quadrupole (+,+,-,-) keeps net charge and dipole zero, so responses reflect quadrupolar coupling instead of rigid drift.
- Damping/stabilize: η shifts poles off the real axis to keep solves stable; stabilize adds a tiny diagonal to prevent singularity for free modes.


"""

import argparse
import numpy as np
import matplotlib.pyplot as plt

from Truss import (
    build_triangular_grid,
    grid_edges,
    assemble_stiffness_dense,
    mass_matrix,
    boundary_nodes,
    apply_dirichlet,
)

# DEBUG: mechanical probing inspired by GF.py (deterministic probing + Cholesky reuse)

# unlimited line lenght in numpy when printing
np.set_printoptions(linewidth=np.inf)

def dynamic_stiffness(K, M, omega, eta=1e-3, stabilize=1e-6):
    """
    Dynamic stiffness: A = K - (omega + i*eta)^2 M.
    eta>0 pushes poles slightly into upper half-plane so we see modes sharply.
    stabilize>0 adds small diagonal real shift to keep system invertible.
    """
    z = omega + 1j * eta
    A = K - (z * z) * M
    if stabilize > 0:
        A = A + stabilize * np.eye(K.shape[0])
    print(f"#DEBUG dynamic_stiffness omega={omega:.6f} eta={eta:.2e} stab={stabilize:.2e}")
    return A


def cholesky_factor(A):
    # Dense Cholesky; may fail if A not HPD -> raise loudly (debug preference)
    L = np.linalg.cholesky(A)
    print(f"#DEBUG cholesky_factor shape={A.shape}")
    return L


def cholesky_solve(L, B):
    # Solve A X = B given A = L L^H.
    Y = np.linalg.solve(L, B)
    X = np.linalg.solve(L.T.conj(), Y)
    return X


def solve_response(K, M, omega, eta, charges, direction_vec, dim=3, stabilize=1e-6):
    """Solve for displacement under dipole force at a single omega."""
    ndof = K.shape[0]
    n_nodes = ndof // dim
    A = dynamic_stiffness(K, M, omega, eta=eta, stabilize=stabilize)
    rhs = np.zeros((ndof, 1), dtype=np.complex128)
    for n in range(n_nodes):
        rhs[n * dim : n * dim + dim, 0] = charges[n] * direction_vec[:dim]
    U = np.linalg.solve(A, rhs)
    return U[:, 0].reshape(n_nodes, dim)


def mechanical_greens_probing(K, M, omegas, eta=1e-3, direction_vec=None, charges=None, dim=3, stabilize=1e-6):
    """
    Dipole-driven probing of mechanical Green's function.
    - K, M: reduced (Dirichlet) matrices
    - omegas: array of target frequencies (real)
    - eta: small damping to push poles above real axis
    - direction_vec: 3-vector direction of the homogeneous field
    - charges: per-node charges for dipole coupling (len = n_nodes)
    Returns dict with spectra and dipole couplings.
    The RHS is the uniform dipole force f_i = q_i * e_dir (unit field).
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
            rhs[n * dim : n * dim + dim, 0] = charges[n] * direction_vec[:dim]
        U = np.linalg.solve(A, rhs)
        spectrum_energy[io] = np.sum(np.abs(U) ** 2) / n_nodes
        disp_nodes = U[:, 0].reshape(n_nodes, dim)
        dip = (charges[:, None] * disp_nodes).sum(axis=0)
        spectrum_dipole[io] = dip
        print(f"#DEBUG omega={omega:.4f} dipole_probe dir={direction_vec} |dip|={np.linalg.norm(dip):.3e}")
    return {
        "omega": np.asarray(omegas),
        "energy": spectrum_energy,
        "dipole": spectrum_dipole,
        "n_probes": len(omegas),  # one dipole probe per frequency
    }


def expand_displacement(disp_reduced, mask, dim=3):
    """Map reduced displacement (after Dirichlet) back to full node array."""
    ndof_full = mask.size
    n_nodes_full = ndof_full // dim
    disp_full = np.zeros((n_nodes_full, dim), dtype=disp_reduced.dtype)
    node_mask = mask.reshape(n_nodes_full, dim)[:, 0]
    disp_full[node_mask, :] = disp_reduced
    return disp_full


def plot_truss(pos, edges, charges, disp=None, scale=0.2, title="truss"):
    """Plot nodes colored by charge and optional displacement vectors."""
    fig, ax = plt.subplots(figsize=(6, 6))
    charges = np.asarray(charges)
    if disp is not None:
        disp = np.real(disp)
    sc = ax.scatter(pos[:, 0], pos[:, 1], c=charges, cmap="coolwarm", s=40, edgecolors="k")
    for i, j in edges:
        ax.plot([pos[i, 0], pos[j, 0]], [pos[i, 1], pos[j, 1]], "k-", lw=0.5, alpha=0.6)
    if disp is not None:
        ax.quiver(pos[:, 0], pos[:, 1], disp[:, 0], disp[:, 1], angles="xy", scale_units="xy", scale=1.0 / scale, color="tab:green", width=0.005, alpha=0.8)
    ax.set_aspect("equal")
    ax.set_title(title)
    fig.colorbar(sc, ax=ax, label="charge")
    plt.tight_layout()
    return fig, ax


def plot_spectrum(omegas, res, eigfreq=None, sel=None):
    """Plot response magnitude vs frequency with eigenvalue markers."""
    dip = res["dipole"]
    mag = np.linalg.norm(dip, axis=1)
    fig, ax = plt.subplots(figsize=(5, 3))
    ax.plot(omegas, mag, "-", ms=3, lw=1)
    if eigfreq is not None:
        eigfreq = np.asarray(eigfreq)
        sel = set(sel) if sel is not None else set()
        for i, f in enumerate(eigfreq):
            color = "r" if i in sel else "0.6"
            ax.axvline(f, color=color, lw=1, alpha=0.8)
    ax.set_xlabel("omega")
    ax.set_ylabel("|dipole response|")
    ax.set_title("Dipole-coupled spectrum")
    plt.tight_layout()
    return fig, ax


def build_test_truss( nx=6, ny=6, a=1.0, jitter=0.0, k_spring=1.0, mass_value=1.0, fixed_boundary="bottom", dim=3,):
    pos = build_triangular_grid(nx, ny, a=a, jitter=jitter)
    edges = grid_edges(nx, ny, include_diag=True)
    K = assemble_stiffness_dense(pos, edges, k_spring=k_spring, dim=dim)
    masses = np.full(nx * ny, mass_value)
    M = mass_matrix(masses, dim=dim)
    fixed = boundary_nodes(nx, ny, which=fixed_boundary)
    K_red, M_red, mask = apply_dirichlet(K, M, fixed, dim=dim)
    print(f"#DEBUG build_test_truss ndof_full={K.shape[0]} ndof_red={K_red.shape[0]}")
    return pos, edges, K_red, M_red, mask


def corner_quadrupole_charges(nx, ny, charge_val=1.0):
    """
    Neutral quadrupole: (+,+,-,-) on corners to kill net charge and dipole.
    Layout:
        (0,ny-1) : +q
        (nx-1,0) : +q
        (0,0)    : -q
        (nx-1,ny-1): -q
    """
    charges = np.zeros(nx * ny)
    idx = lambda ix, iy: iy * nx + ix
    charges[idx(0   ,    0)] = -charge_val
    charges[idx(0   , ny-1)] = charge_val
    charges[idx(nx-1, 0   )] = -charge_val
    charges[idx(nx-1, ny-1)] = charge_val
    return charges

def assemble_weighted_stiffness(pos, edges, k_edges, dim=3):
    """Assemble stiffness with per-edge spring constants k_edges (len==edges)."""
    n_nodes = pos.shape[0]
    ndof = dim * n_nodes
    K = np.zeros((ndof, ndof), dtype=np.float64)
    for (i, j), k_spring in zip(edges, k_edges):
        d = pos[j] - pos[i]
        L = np.linalg.norm(d)
        if L <= 1e-12:
            continue
        u = d / L
        k_fac = k_spring / (L * L)
        ia = i * dim
        ja = j * dim
        outer = k_fac * np.outer(u, u)
        K[ia:ia+dim, ia:ia+dim] += outer
        K[ja:ja+dim, ja:ja+dim] += outer
        K[ia:ia+dim, ja:ja+dim] -= outer
        K[ja:ja+dim, ia:ia+dim] -= outer
    return K


def classify_edge(pos, i, j, tol=1e-6):
    d = pos[j] - pos[i]
    dx, dy = abs(d[0]), abs(d[1])
    if dy < tol:
        return "x"
    if dx < tol:
        return "y"
    return "diag"


def parse_indices(s):
    if s is None or s == "":
        return []
    return [int(x) for x in s.split(",") if x.strip() != ""]


def plot_modes_with_response(pos, edges, charges, eigvecs_full, respvecs_full, freqs, sel_idx, scale=0.15):
    nsel = len(sel_idx)
    if nsel == 0:
        return None
    cols = min(3, nsel)
    rows = (nsel + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(4 * cols, 4 * rows))
    axes = np.atleast_1d(axes).ravel()
    charges = np.asarray(charges)
    for ax, idx, eig, resp in zip(axes, sel_idx, eigvecs_full, respvecs_full):
        eig = np.real(eig)
        resp = np.real(resp)
        sc = ax.scatter(pos[:, 0], pos[:, 1], c=charges, cmap="coolwarm", s=30, edgecolors="k")
        for i, j in edges:
            ax.plot([pos[i, 0], pos[j, 0]], [pos[i, 1], pos[j, 1]], "k-", lw=0.5, alpha=0.5)
        ax.quiver(pos[:, 0], pos[:, 1], eig[:, 0], eig[:, 1], angles="xy", scale_units="xy", scale=1.0 / scale, color="tab:blue", width=0.004, alpha=0.8, label="eig")
        ax.quiver(pos[:, 0], pos[:, 1], resp[:, 0], resp[:, 1], angles="xy", scale_units="xy", scale=1.0 / scale, color="tab:orange", width=0.004, alpha=0.8, label="resp")
        ax.set_aspect("equal")
        ax.set_title(f"mode {idx} f={freqs[idx]:.3f}")
    for ax in axes[nsel:]:
        ax.axis("off")
    fig.colorbar(sc, ax=axes.tolist(), fraction=0.025, pad=0.02, label="charge")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right")
    fig.tight_layout()
    return fig


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Dipole-driven vibration probing of a spring-mass grid.")
    p.add_argument("--nx", type=int, default=3, help="Grid size in x")
    p.add_argument("--ny", type=int, default=3, help="Grid size in y")
    p.add_argument("--a", type=float, default=1.0, help="Lattice spacing")
    p.add_argument("--jitter", type=float, default=0.02, help="Position jitter")
    p.add_argument("--k1", type=float, default=5.0, help="Spring stiffness along x")
    p.add_argument("--k2", type=float, default=5.0, help="Spring stiffness along y")
    p.add_argument("--kdiag", type=float, default=5.0, help="Spring stiffness on diagonals")
    p.add_argument("--mass", type=float, default=1.0, help="Mass per node")
    p.add_argument("--charge", type=float, default=1.0, help="Corner charge magnitude (quadrupole)")
    p.add_argument("--theta", type=float, default=90.0, help="Dipole angle in degrees (0=x, 90=y)")
    p.add_argument("--fmin", type=float, default=0.1, help="Min frequency")
    p.add_argument("--fmax", type=float, default=6.0, help="Max frequency")
    p.add_argument("--nfreq", type=int, default=1000, help="Number of frequency samples")
    p.add_argument("--eta", type=float, default=0.001, help="Damping shift")
    p.add_argument("--stabilize", type=float, default=1e-6, help="Diagonal stabilizer")
    p.add_argument("--eig_idx", type=str, default="15,16,17,18,19", help="Comma indices of eigenmodes to plot")
    p.add_argument("--verbosity", type=int, default=2, help="Verbosity level; >2 prints eigenvalues")
    args = p.parse_args()

    nx, ny = args.nx, args.ny
    pos = build_triangular_grid(nx, ny, a=args.a, jitter=args.jitter)
    edges = grid_edges(nx, ny, include_diag=True)
    k_edges = []
    for (i, j) in edges:
        cls = classify_edge(pos, i, j)
        if cls == "x":
            k_edges.append(args.k1)
        elif cls == "y":
            k_edges.append(args.k2)
        else:
            k_edges.append(args.kdiag)
    K = assemble_weighted_stiffness(pos, edges, k_edges, dim=3)
    masses = np.full(nx * ny, args.mass)
    M = mass_matrix(masses, dim=3)
    fixed = boundary_nodes(nx, ny, which="none")
    K_red, M_red, mask = apply_dirichlet(K, M, fixed, dim=3)

    omegas = np.linspace(args.fmin, args.fmax, args.nfreq)
    charges_full = corner_quadrupole_charges(nx, ny, charge_val=args.charge)
    node_mask = mask.reshape(pos.shape[0], 3)[:, 0]
    charges = charges_full[node_mask]
    theta = np.deg2rad(args.theta)
    direction_vec = np.array([np.cos(theta), np.sin(theta), 0.0], dtype=np.float64)

    res = mechanical_greens_probing(K_red, M_red, omegas, eta=args.eta, direction_vec=direction_vec, charges=charges, stabilize=args.stabilize)
    omega_idx = min(len(omegas) - 1, max(0, args.nfreq // 4))
    disp_red = solve_response(K_red, M_red, omegas[omega_idx], eta=args.eta, charges=charges, direction_vec=direction_vec, dim=3, stabilize=args.stabilize)
    disp_full = expand_displacement(disp_red, mask, dim=3)
    plot_truss(pos, edges, charges_full, disp=disp_full, scale=0.15, title=f"omega={omegas[omega_idx]:.2f} dipole response")
    # eigen-spectrum for comparison
    A = np.linalg.inv(M_red) @ K_red
    w, V = np.linalg.eigh(A)
    freq = np.sqrt(np.clip(w, 0.0, None))
    if args.verbosity > 1:
        print("#DEBUG eig freq:", freq)
    sel_idx = [i for i in parse_indices(args.eig_idx) if i >= 0 and i < len(freq)]
    resp_red_list = []
    eig_full_list = []
    resp_full_list = []
    for m in sel_idx:
        vec_red = V[:, m].reshape(-1, 3)
        eig_full_list.append(expand_displacement(vec_red, mask, dim=3))
        resp_red = solve_response(K_red, M_red, freq[m], eta=args.eta, charges=charges, direction_vec=direction_vec, dim=3, stabilize=args.stabilize)
        resp_full_list.append(expand_displacement(resp_red, mask, dim=3))
    plot_spectrum(omegas, res, eigfreq=freq, sel=sel_idx)
    plot_modes_with_response(pos, edges, charges_full, eig_full_list, resp_full_list, freq, sel_idx, scale=0.15)
    plt.show()
    print("#DEBUG spectra sample", res["energy"][:5])
