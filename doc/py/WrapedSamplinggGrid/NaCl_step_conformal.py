import argparse
import numpy as np
import matplotlib.pyplot as plt

Q0 = 0.7


def parse_args():
    p = argparse.ArgumentParser(description="NaCl step-edge slab generator")
    p.add_argument('--a',       type=float, default=2.82065, help='lattice constant (Å)')
    p.add_argument('--nx_sub',  type=int,   default=12, help='number of units along x for substrate (lower full layers)')
    p.add_argument('--nx_step', type=int,   default=6, help='number of units along x for step top (shorter upper terrace)')
    p.add_argument('--ny',      type=int,   default=8, help='number of units along y')
    p.add_argument('--nz_sub',  type=int,   default=2, help='number of full substrate layers (thick bottom)')
    p.add_argument('--nz_step', type=int,   default=1, help='number of layers in step region (upper terrace height)')
    p.add_argument('--outfile', type=str,   default='NaCl_step_edge.xyz', help='output XYZ filename')
    p.add_argument('--plot',    type=int,   default=1, help='show side-view plot (y~0 slice)')
    p.add_argument('--yslice',  type=float, default=0.0, help='y-slice (Å) for plotting')
    p.add_argument('--tol',     type=float, default=0.2, help='tolerance (fraction of a) for including atoms in slice plot')
    p.add_argument('--edge_span_x',     type=int, default=4, help='number of lattice sites to keep on each side of the step edge (along x)')
    p.add_argument('--edge_top_layers', type=int, default=2, help='number of top layers (from surface) eligible as boundary')
    p.add_argument('--neighbor_mode',   type=str, default='8', choices=['4','8'], help='boundary test connectivity (4 or 8 neighbors in x–z)')
    p.add_argument('--u_density', type=float,  default=1, help='conformal sampling density for U isolines')
    p.add_argument('--v_density', type=int,    default=10,  help='conformal sampling density for V streamlines')
    p.add_argument('--u_min', type=float, default=24.0,   help='explicit U minimum (log|W|)')
    p.add_argument('--u_max', type=float, default=30.0,  help='explicit U maximum (log|W|)')
    p.add_argument('--z_max', type=float, default=None,  help='optional absolute z max for sampling/plot (Å)')
    return p.parse_args()


def generate_nacl_step(nx_sub, nx_step, ny, nz_sub, nz_step, a):
    """
    Build an orthogonal NaCl slab with a single straight step edge.
    - Bottom (substrate) spans nx_sub in x and nz_sub layers.
    - Top terrace spans nx_step in x and sits above the substrate by nz_sub layers.
    Na/Cl alternate by parity (ix + iy + iz).
    Returns list of (symbol, x, y, z, charge).
    """
    positions = []
    # Substrate (full) block
    for iz in range(nz_sub):
        for iy in range(ny):
            for ix in range(nx_sub):
                x = ix * a
                y = iy * a
                z = iz * a
                parity = ix + iy + iz
                if parity % 2 == 0:
                    S, Q = 'Na', +Q0
                else:
                    S, Q = 'Cl', -Q0
                positions.append((S, x, y, z, Q, ix, iy, iz))

    # Step (upper) block shifted up by nz_sub layers
    z_base = nz_sub * a
    for iz in range(nz_step):
        for iy in range(ny):
            for ix in range(nx_step):
                x = ix * a
                y = iy * a
                z = z_base + iz * a
                parity = ix + iy + (iz + nz_sub)  # continuous parity in height
                if parity % 2 == 0:
                    S, Q = 'Na', +Q0
                else:
                    S, Q = 'Cl', -Q0
                positions.append((S, x, y, z, Q, ix, iy, iz + nz_sub))
    return positions


def save_xyz(filename, positions, nx_sub, nx_step, ny, nz_sub, nz_step, a):
    nx_cell = max(nx_sub, nx_step)
    nz_cell = nz_sub + nz_step
    cell = [
        [nx_cell * a, 0, 0],
        [0, ny * a, 0],
        [0, 0, nz_cell * a]
    ]
    with open(filename, 'w') as fout:
        fout.write(f"{len(positions)}\n")
        fout.write(f"lvs {cell[0][0]} {cell[0][1]} {cell[0][2]}   {cell[1][0]} {cell[1][1]} {cell[1][2]}   {cell[2][0]} {cell[2][1]} {cell[2][2]}\n")
        for S, x, y, z, Q, _, _, _ in positions:
            fout.write(f"{S} {x:10.5f} {y:10.5f} {z:10.5f} {Q}\n")


def select_boundary_atoms(positions, nx_sub, nx_step, ny, nz_sub, nz_step, a, edge_span_x=2, edge_top_layers=2, y_slice=0.0, tol=0.2, neighbor_mode='8'):
    """
    Select atoms on the step boundary within a y-slice, limited x-span near the edge, and top layers.
    Boundary test: any missing neighbor in 8-connected (dx,dz) grid for the same iy.
    """
    nx_cell = max(nx_sub, nx_step)
    nz_cell = nz_sub + nz_step
    occ = {(ix, iy, iz) for (_, _, _, _, _, ix, iy, iz) in positions}

    edge_x_center = nx_step  # interface between terrace (0..nx_step-1) and lower step continuation
    x_min = max(0, edge_x_center - edge_span_x)
    x_max = min(nx_cell - 1, edge_x_center + edge_span_x)
    z_min_sel = max(0, nz_cell - edge_top_layers)

    selected = []
    # neighbor stencil
    if neighbor_mode == '4':
        neighbor_offsets = [(1,0), (-1,0), (0,1), (0,-1)]
    else:
        neighbor_offsets = [(dx,dz) for dx in (-1,0,1) for dz in (-1,0,1) if not (dx==0 and dz==0)]

    for (S, x, y, z, Q, ix, iy, iz) in positions:
        if not (x_min <= ix <= x_max):
            continue
        if iz < z_min_sel:
            continue
        if abs(y - y_slice) > tol * a:
            continue
        boundary = False
        for dx, dz in neighbor_offsets:
            n_ix, n_iz = ix + dx, iz + dz
            if n_ix < 0 or n_ix >= nx_cell or n_iz < 0 or n_iz >= nz_cell:
                boundary = True
                break
            if (n_ix, iy, n_iz) not in occ:
                boundary = True
                break
        if boundary:
            selected.append((S, x, y, z, Q, ix, iy, iz))
    return selected


def plot_side_view(positions, boundary_atoms=None, y_slice=0.0, a=2.82065, tol=0.1):
    """
    Plot x–z side view for atoms whose y is within tol*a of y_slice.
    """
    xs_Na, zs_Na = [], []
    xs_Cl, zs_Cl = [], []
    for S, x, y, z, Q, ix, iy, iz in positions:
        if abs(y - y_slice) <= tol * a:
            if S == 'Na':
                xs_Na.append(x)
                zs_Na.append(z)
            else:
                xs_Cl.append(x)
                zs_Cl.append(z)

    plt.figure(figsize=(8, 4))
    plt.scatter(xs_Cl, zs_Cl, marker='o', color='orange', label='Cl')
    plt.scatter(xs_Na, zs_Na, marker='*', color='steelblue', label='Na')

    if boundary_atoms:
        bx = [p[1] for p in boundary_atoms]
        bz = [p[3] for p in boundary_atoms]
        plt.scatter(bx, bz, marker='s', color='red', s=60, facecolors='none', linewidths=1.5, label='Boundary')
    plt.xlabel('x (Å)')
    plt.ylabel('z (Å)')
    plt.title('NaCl step-edge side view (x–z slice)')
    plt.legend()
    plt.gca().set_aspect('equal', adjustable='box')
    plt.tight_layout()



def run_conformal_sampling(boundary_atoms, u_density, v_density, u_min, u_max, x_min, x_max, z_min, z_max):
    """
    Lightweight ConformalFlow6 sampler restricted to given atoms.
    boundary_atoms: list of (S, x, y, z, Q, ix, iy, iz) ; S in {'Na','Cl'}
    """
    # Map to (x,y,r,multiplicity)
    r_Na = 0.2
    r_Cl = 0.25
    atoms = []  # (x, z, r, multiplicity)
    for S, x, _, z, _, _, _, _ in boundary_atoms:
        if S == 'Cl':
            atoms.append((x, z, r_Cl, 2))
        else:
            atoms.append((x, z, r_Na, 1))

    if not atoms:
        print("No boundary atoms provided; skipping conformal sampling")
        return []

    com_x = np.mean([a[0] for a in atoms])
    com_z = np.mean([a[1] for a in atoms])

    base_poly = [1.0]
    for (cx, cz, r, weight) in atoms:
        weight_int = int(round(weight))
        center = cx + 1j * cz
        for _ in range(weight_int):
            base_poly = np.convolve(base_poly, [1.0, -center])

    u_surface_vals = []
    for (cx, cz, r, _) in atoms:
        for angle in [0, np.pi/2, np.pi, 3*np.pi/2]:
            test_z = (cx + r * np.cos(angle)) + 1j * (cz + r * np.sin(angle))
            val = np.log(np.abs(np.polyval(base_poly, test_z)))
            u_surface_vals.append(val)
    step_U = 1.0 / u_density
    freq_U = u_density
    freq_V = v_density / (2 * np.pi)

    U_targets = np.arange(u_min, u_max, step_U)
    V_targets = np.linspace(-np.pi, np.pi, v_density, endpoint=False)

    sample_points = []
    point_id_counter = 0
    def match_roots(prev_roots, curr_roots):
        ordered = []
        pool = list(curr_roots)
        for p in prev_roots:
            dists = [abs(c - p) for c in pool]
            best = np.argmin(dists)
            ordered.append(pool.pop(best))
        return ordered

    for v_idx, v in enumerate(V_targets):
        prev_roots = None
        for u_idx, u in enumerate(U_targets):
            w_target = np.exp(u + 1j * v)
            poly = np.copy(base_poly)
            poly[-1] -= w_target
            roots = np.roots(poly)
            if prev_roots is None:
                ordered_roots = roots
            else:
                ordered_roots = match_roots(prev_roots, roots)
            prev_roots = ordered_roots
            for branch_idx, root in enumerate(ordered_roots):
                rx, rz = np.real(root), np.imag(root)
                if rx < x_min or rx > x_max:
                    continue
                if rz < z_min or (z_max is not None and rz > z_max):
                    continue
                point_data = {
                    'id': point_id_counter,
                    'x': rx, 'z': rz,
                    'u_idx': u_idx, 'v_idx': v_idx,
                    'u_val': u, 'v_val': v,
                    'branch': branch_idx
                }
                sample_points.append(point_data)
                point_id_counter += 1
    print(f"Conformal sample points (clipped): {len(sample_points)} | U_targets={len(U_targets)} | U_range=({u_min},{u_max})")
    return sample_points, atoms, base_poly, U_targets, V_targets


def plot_sampling(atoms, sample_points, base_poly, u_density, v_density, z_max=None):
    # Checkerboard background in x–z plane
    xs_all = [a[0] for a in atoms]
    zs_all = [a[1] for a in atoms]
    if sample_points:
        xs_all += [p['x'] for p in sample_points]
        zs_all += [p['z'] for p in sample_points]
    pad_x = 2.0
    pad_z = 3.0
    x_min, x_max = min(xs_all)-pad_x, max(xs_all)+pad_x
    z_min_plot, z_max_plot = min(zs_all)-pad_z, max(zs_all)+pad_z
    if z_max is not None:
        z_max_plot = z_max
    res = 500
    x_lin = np.linspace(x_min, x_max, res)
    z_lin = np.linspace(z_min_plot, z_max_plot, res)
    X, Z = np.meshgrid(x_lin, z_lin)
    Zc = X + 1j * Z
    W = np.polyval(base_poly, Zc)
    U = np.log(np.abs(W))
    V = np.angle(W)
    freq_U = u_density
    freq_V = v_density / (2 * np.pi)

    U_targets = np.arange(args.u_min, args.u_max, 1.0 / u_density)
    V_targets = np.linspace(-np.pi, np.pi, v_density, endpoint=False)

    plt.figure(figsize=(7, 6))
    ax = plt.gca()
    ax.imshow((np.floor(U * freq_U) + np.floor(V * freq_V)) % 2, extent=[x_min, x_max, z_min_plot, z_max_plot], origin='lower', cmap='Greys', alpha=0.3)

    has_atoms = False
    for (x, z, r, w) in atoms:
        circle = plt.Circle((x, z), r, fill=False, linestyle='--', color='gray', alpha=0.6)
        ax.add_patch(circle)
        ax.plot(x, z, 'k+', markersize=6)
        has_atoms = True
    if sample_points:
        xs = [p['x'] for p in sample_points]
        zs = [p['z'] for p in sample_points]
        ax.plot(xs, zs, 'r.', alpha=0.5, markersize=4, label='Samples')
    ax.set_aspect('equal')
    ax.set_xlabel('x (Å)')
    ax.set_ylabel('z (Å)')
    ax.set_title('Conformal sampling around boundary atoms (x–z plane)')
    handles, labels = ax.get_legend_handles_labels()
    if labels:
        ax.legend()
    plt.tight_layout()

if __name__ == '__main__':
    args = parse_args()
    positions = generate_nacl_step(args.nx_sub, args.nx_step, args.ny, args.nz_sub, args.nz_step, args.a)
    save_xyz(args.outfile, positions, args.nx_sub, args.nx_step, args.ny, args.nz_sub, args.nz_step, args.a)
    boundary_atoms = select_boundary_atoms(
        positions,
        args.nx_sub, args.nx_step, args.ny,
        args.nz_sub, args.nz_step,
        args.a,
        edge_span_x=args.edge_span_x,
        edge_top_layers=args.edge_top_layers,
        y_slice=args.yslice,
        tol=args.tol,
        neighbor_mode=args.neighbor_mode,
    )
    print(f"Selected boundary atoms: {len(boundary_atoms)} (edge_span_x={args.edge_span_x}, edge_top_layers={args.edge_top_layers})")
    # Run conformal sampling around boundary atoms (Na/Cl multiplicities)
    atoms_for_sampling = []
    for S, x, y, z, _, _, _, _ in boundary_atoms:
        if S == 'Cl':
            atoms_for_sampling.append((x, z, 0.25, 2))
        else:
            atoms_for_sampling.append((x, z, 0.2, 1))
    # Spatial clamps from boundary atoms
    x_min = min(a[0] for a in atoms_for_sampling)
    x_max = max(a[0] for a in atoms_for_sampling)
    z_min = min(a[1] for a in atoms_for_sampling)
    z_max = args.z_max

    sample_points, atoms_for_sampling, base_poly, U_targets, V_targets = run_conformal_sampling(
        boundary_atoms, args.u_density, args.v_density, args.u_min, args.u_max,
        x_min, x_max, z_min, z_max,
    )
    if args.plot:
        plot_side_view(positions, boundary_atoms=boundary_atoms, y_slice=args.yslice, a=args.a, tol=args.tol)
        plot_sampling(atoms_for_sampling, sample_points, base_poly, args.u_density, args.v_density, z_max=z_max)
    plt.show()