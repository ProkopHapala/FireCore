import argparse
import os
import numpy as np
import matplotlib.pyplot as plt

from sampling_lib import (
    AtomicSystem,
    blend_uv_grids,
    build_grid_uv,
    compute_sdf_2d,
    count_5d,
    estimate_scan_total,
    generate_tilts,
    local_up_from_streamline,
    make_atoms_nacl_step,
    make_u_levels,
    make_y_levels,
    make_pose_coords_local,
    parse_int_list,
    parse_u_list,
    save_xyz_frames,
)


AXES = ['u', 'v', 'y', 'tilt', 'phi']


def _local_path(default_name):
    return os.path.join(os.path.dirname(__file__), default_name)


def parse_args(argv=None):
    ap = argparse.ArgumentParser()

    # --- SDF grid args (mostly from SDF_Edge_Sonet3) ---
    ap.add_argument('--a',           type=float, default=2.82065)
    ap.add_argument('--r_Na',        type=float, default=1.4)
    ap.add_argument('--r_Cl',        type=float, default=2.0)
    ap.add_argument('--corner',      type=str,   default='Na', choices=['Na', 'Cl'])
    ap.add_argument('--n_sub',       type=int,   default=9, help='samples per isoline between L and R (n_v)')
    ap.add_argument('--u_levels',    type=int,   default=20)
    ap.add_argument('--u_min',       type=float, default=0.05)
    ap.add_argument('--u_max',       type=float, default=5.0)
    ap.add_argument('--u_schedule',  type=str,   default='sinh', choices=['linear', 'sinh', 'geom'])
    ap.add_argument('--u_sinh_a',    type=float, default=2.5)
    ap.add_argument('--u_list',      type=str,   default=None)
    ap.add_argument('--combine_p',   type=float, default=4.0)
    ap.add_argument('--walk_res',    type=float, default=0.005)
    ap.add_argument('--sdf_res',     type=int,   default=1000)

    # --- y levels ---
    # TODO: y-sample should be in unit of lattice constant, and because Na/Cl atoms swap position in neighbor unit cell, we should interpolated between the two cases
    ap.add_argument('--y_min',       type=float, default=-2.83)
    ap.add_argument('--y_max',       type=float, default=+2.83)
    ap.add_argument('--y_levels',    type=int,   default=5)

    # --- tilt/phi levels ---
    ap.add_argument('--fans',        type=int,   default=4)
    ap.add_argument('--subdiv',      type=int,   default=2)
    ap.add_argument('--theta_max',   type=float, default=45.0)
    ap.add_argument('--radius',      type=float, default=1.0)
    ap.add_argument('--phi_samples', type=int,   default=8)
    ap.add_argument('--symmetry_m',  type=int,   default=1)
    ap.add_argument('--phi_max',     type=float, default=180.0)

    # --- selection / mode ---
    ap.add_argument('--mode',        type=str, default='counts', choices=['counts', 'random', 'scan'])

    ap.add_argument('--coarse_u',    type=str, default='0,10,19', help='comma-separated u indices into valid u grid')
    ap.add_argument('--coarse_v',    type=str, default='0,3,5,8', help='comma-separated v indices into grid columns')
    ap.add_argument('--coarse_y',    type=str, default='0,4', help='comma-separated y indices')
    ap.add_argument('--coarse_tilt', type=str, default='0', help='comma-separated tilt indices')
    ap.add_argument('--coarse_phi',  type=str, default='0,4', help='comma-separated phi indices')

    ap.add_argument('--scan_axes',   type=str, default='u,v,y,tilt,phi', help='which axes to do 1D scan; comma-separated')

    ap.add_argument('--Nrandom',     type=int, default=0)
    ap.add_argument('--seed',        type=int, default=123)

    ap.add_argument('--Ncap',        type=int, default=200000, help='cap for total number of frames emitted in scan mode')

    # --- outputs ---
    ap.add_argument('--png',         type=str, default=None, help='scan visualization PNG path')
    ap.add_argument('--xyz',         type=str, default='../../../cpp/common_resources/xyz/H2O_O.xyz', help='input molecule xyz')
    ap.add_argument('--out_prefix',  type=str, default=_local_path('samples'), help='output prefix; files are <prefix>_<tag>.xyz')

    try:
        return ap.parse_args(argv)
    except SystemExit:
        return ap.parse_args([])


def build_spatial_grid(args):
    a = args.a
    x_left = 0
    x_right = 2 * a
    pad = args.u_max + a
    x_range = (x_left - pad, x_right + pad)
    z_range = (-a - 0.3, a + args.u_max + 1.5)

    if args.u_list:
        u_levels = parse_u_list(args.u_list)
    else:
        u_levels = make_u_levels(args.u_min, args.u_max, args.u_levels, args.u_schedule, args.u_sinh_a)

    atoms_na = make_atoms_nacl_step(a, args.r_Na, args.r_Cl, corner_species='Na')
    X_grid_na, Z_grid_na, sdf_na, x_lin_na, z_lin_na = compute_sdf_2d(atoms_na, x_range, z_range, args.sdf_res, args.combine_p)
    grid_uv_na, u_vals_na, u_idx_valid_na, arc_lengths_na = build_grid_uv(sdf_na, x_lin_na, z_lin_na, u_levels, x_left, x_right, args.n_sub, args.walk_res)
    if grid_uv_na is None:
        raise RuntimeError('No valid isolines found for corner=Na; grid_uv_na is None')

    atoms_cl = make_atoms_nacl_step(a, args.r_Na, args.r_Cl, corner_species='Cl')
    X_grid_cl, Z_grid_cl, sdf_cl, x_lin_cl, z_lin_cl = compute_sdf_2d(atoms_cl, x_range, z_range, args.sdf_res, args.combine_p)
    grid_uv_cl, u_vals_cl, u_idx_valid_cl, arc_lengths_cl = build_grid_uv(sdf_cl, x_lin_cl, z_lin_cl, u_levels, x_left, x_right, args.n_sub, args.walk_res)
    if grid_uv_cl is None:
        raise RuntimeError('No valid isolines found for corner=Cl; grid_uv_cl is None')

    if grid_uv_na.shape != grid_uv_cl.shape:
        raise RuntimeError(f"Na/Cl grids have different shapes: {grid_uv_na.shape} vs {grid_uv_cl.shape}")

    # keep u_vals from Na (should match Cl if u_levels same and both valid)
    return (atoms_na, atoms_cl), (X_grid_na, Z_grid_na, sdf_na, x_lin_na, z_lin_na), (grid_uv_na, grid_uv_cl, u_vals_na, u_idx_valid_na, arc_lengths_na, x_left, x_right)


def plot_up_vectors(grid_uv_y, out_png, u_sel=None, v_sel=None, scale=0.4):
    # grid_uv_y: (n_y,n_u,n_v,2)
    n_y, n_u, n_v = grid_uv_y.shape[:3]
    if u_sel is None:
        u_sel = [0, n_u//2, n_u-1]
    if v_sel is None:
        v_sel = [0, n_v//2, n_v-1]
    y_idx = n_y // 2
    grid_uv = grid_uv_y[y_idx]

    fig, ax = plt.subplots(1, 1, figsize=(7, 6))
    ax.set_title(f"Local up-vectors from streamline connector (y_idx={y_idx})")
    ax.plot(grid_uv[:, :, 0].ravel(), grid_uv[:, :, 1].ravel(), '.', ms=2, alpha=0.2, color='#888')
    for iv in v_sel:
        for iu in u_sel:
            xz = grid_uv[iu, iv]
            up = local_up_from_streamline(grid_uv, iu, iv)
            dx, dz = up[0]*scale, up[2]*scale
            ax.arrow(xz[0], xz[1], dx, dz, head_width=0.05, head_length=0.08, fc='r', ec='r', alpha=0.9)
            ax.text(xz[0], xz[1], f"({iu},{iv})", fontsize=7, color='k')
    ax.set_aspect('equal')
    ax.set_xlabel('x (A)'); ax.set_ylabel('z (A)')
    fig.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches='tight')
    plt.close(fig)


def build_angular_grid(args):
    verts_base, wedges, sub_points, disk_pts, ring_ids, dome_pts, tilts, phi_grid = generate_tilts(
        args.fans, args.subdiv, args.theta_max, args.radius, args.phi_samples, args.symmetry_m, phi_max_deg=args.phi_max,
    )
    return (verts_base, wedges, sub_points, disk_pts, ring_ids, dome_pts, tilts, phi_grid)


def _parse_coarse(args, n_u, n_v, n_y, n_tilt, n_phi):
    coarse_raw = {
        'u': parse_int_list(args.coarse_u),
        'v': parse_int_list(args.coarse_v),
        'y': parse_int_list(args.coarse_y),
        'tilt': parse_int_list(args.coarse_tilt),
        'phi': parse_int_list(args.coarse_phi),
    }

    coarse = {}
    for ax, n in [('u', n_u), ('v', n_v), ('y', n_y), ('tilt', n_tilt), ('phi', n_phi)]:
        vals = coarse_raw[ax] if coarse_raw[ax] is not None else [0]
        if n <= 0:
            vals = []
        fixed = []
        for i in vals:
            if 0 <= i < n:
                fixed.append(i)
            else:
                print(f"WARNING: coarse_{ax} index {i} out of range [0,{n}); dropping")
        if not fixed and n > 0:
            fixed = [0]
            print(f"WARNING: coarse_{ax} became empty; defaulting to [0]")
        coarse[ax] = sorted(set(fixed))

    return coarse


def _parse_scan_axes(args):
    scan_axes = [s.strip() for s in args.scan_axes.split(',') if s.strip() != '']
    for ax in scan_axes:
        if ax not in AXES:
            raise ValueError(f'Unknown axis in --scan_axes: {ax}')
    return scan_axes


def _save_scan_png(args, grid_uv, u_vals, y_vals, tilts, phi_grid, coarse, scan_axes, out_png):
    fig = plt.figure(figsize=(12, 10))
    gs = fig.add_gridspec(2, 2)

    ax0 = fig.add_subplot(gs[0, 0])
    ax0.set_title('Spatial grid (x,z): scan lines overlay')
    y_idx = len(y_vals)//2
    g = grid_uv[y_idx]
    ax0.scatter(g[:, :, 0].ravel(), g[:, :, 1].ravel(), s=3, c='#bbb', alpha=0.25, edgecolors='none')

    n_u, n_v = g.shape[:2]

    if 'u' in scan_axes:
        for j in coarse['v']:
            col = g[:, j]
            ax0.plot(col[:, 0], col[:, 1], lw=1.4, alpha=0.9, label=f'scan u @ v={j}')

    if 'v' in scan_axes:
        for i in coarse['u']:
            row = g[i, :]
            ax0.plot(row[:, 0], row[:, 1], lw=1.4, alpha=0.9, label=f'scan v @ u={i}')

    ax0.set_aspect('equal')
    ax0.set_xlabel('x (A)'); ax0.set_ylabel('z (A)')
    ax0.legend(fontsize=7)

    ax1 = fig.add_subplot(gs[0, 1])
    ax1.set_title('Y levels: scan / coarse')
    ax1.plot(np.zeros_like(y_vals), y_vals, 'k.', ms=8, alpha=0.6, label='all y')
    ax1.plot(np.zeros(len(coarse['y'])), y_vals[coarse['y']], 'ro', ms=8, label='coarse y')
    if 'y' in scan_axes:
        ax1.plot(np.zeros_like(y_vals), y_vals, 'g-', lw=2.0, alpha=0.5, label='scan y')
    ax1.set_xlabel(''); ax1.set_ylabel('y (A)')
    ax1.legend(fontsize=7)

    ax2 = fig.add_subplot(gs[1, 0], projection='3d')
    ax2.set_title('Tilt unit vectors: scan / coarse')
    ax2.scatter(tilts[:, 0], tilts[:, 1], tilts[:, 2], s=10, c='#888', alpha=0.35)
    ax2.scatter(tilts[coarse['tilt'], 0], tilts[coarse['tilt'], 1], tilts[coarse['tilt'], 2], s=60, c='r', alpha=0.9)
    if 'tilt' in scan_axes:
        ax2.plot(tilts[:, 0], tilts[:, 1], tilts[:, 2], color='g', lw=1.5, alpha=0.5)
    max_range = (tilts.max(axis=0) - tilts.min(axis=0)).max()
    mid = tilts.mean(axis=0)
    for setter, m in zip([ax2.set_xlim, ax2.set_ylim, ax2.set_zlim], mid):
        setter(m - 0.5 * max_range, m + 0.5 * max_range)

    ax3 = fig.add_subplot(gs[1, 1])
    ax3.set_title('Phi samples: scan / coarse')
    ax3.plot(np.cos(phi_grid), np.sin(phi_grid), 'k.', ms=8, alpha=0.6, label='all phi')
    ax3.plot(np.cos(phi_grid[coarse['phi']]), np.sin(phi_grid[coarse['phi']]), 'ro', ms=8, label='coarse phi')
    if 'phi' in scan_axes:
        ax3.plot(np.cos(phi_grid), np.sin(phi_grid), 'g-', lw=2.0, alpha=0.5, label='scan phi')
    ax3.axis('equal')
    ax3.legend(fontsize=7)

    fig.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches='tight')
    plt.close(fig)


def _iter_random_indices(rng, n_u, n_v, n_y, n_tilt, n_phi, N):
    for _ in range(N):
        yield (
            int(rng.integers(0, n_u)),
            int(rng.integers(0, n_v)),
            int(rng.integers(0, n_y)),
            int(rng.integers(0, n_tilt)),
            int(rng.integers(0, n_phi)),
        )


def _iter_coarse_product(coarse):
    for iu in coarse['u']:
        for iv in coarse['v']:
            for iy in coarse['y']:
                for it in coarse['tilt']:
                    for ip in coarse['phi']:
                        yield (iu, iv, iy, it, ip)


def _write_xyz_for_indices(tag, indices_iter, sys, grid_uv, y_vals, tilts, phi_grid, out_prefix):
    coords_frames = []
    comments = []
    for (iu, iv, iy, it, ip) in indices_iter:
        xz = grid_uv[iy, iu, iv]
        y = y_vals[iy]
        up = local_up_from_streamline(grid_uv[iy], iu, iv)
        tvec = tilts[it]
        phi = phi_grid[ip]
        pos = np.array([xz[0], y, xz[1]])
        coords = make_pose_coords_local(sys, pos, up, tvec, phi, center=True)
        coords_frames.append(coords)
        comments.append(f"iu={iu} iv={iv} iy={iy} it={it} ip={ip}  x={pos[0]:.6f} y={pos[1]:.6f} z={pos[2]:.6f}  phi_deg={np.rad2deg(phi):.3f}  tilt=({tvec[0]:.6f},{tvec[1]:.6f},{tvec[2]:.6f})")

    out_xyz = f"{out_prefix}_{tag}.xyz"
    save_xyz_frames(sys, coords_frames, comments, out_xyz)
    print(f"Saved {len(coords_frames)} frames -> {out_xyz}")


def _write_xyz_line(axis_name, sweep_values, fixed, sys, grid_uv, y_vals, tilts, phi_grid, out_prefix):
    # fixed: dict with keys {'u','v','y','tilt','phi'} giving fixed indices for non-swept axes
    label_parts = []
    for k in ['u', 'v', 'y', 'tilt', 'phi']:
        if k == axis_name:
            continue
        label_parts.append(f"{k}={fixed[k]}")
    fixed_label = '_'.join(label_parts)
    out_xyz = f"{out_prefix}_{axis_name}_line_{fixed_label}.xyz"

    coords_frames = []
    comments = []
    for idx in sweep_values:
        iu = fixed['u'] if axis_name != 'u' else idx
        iv = fixed['v'] if axis_name != 'v' else idx
        iy = fixed['y'] if axis_name != 'y' else idx
        it = fixed['tilt'] if axis_name != 'tilt' else idx
        ip = fixed['phi'] if axis_name != 'phi' else idx
        xz = grid_uv[iy, iu, iv]
        y = y_vals[iy]
        up = local_up_from_streamline(grid_uv[iy], iu, iv)
        tloc = tilts[it]
        phi = phi_grid[ip]
        pos = np.array([xz[0], y, xz[1]])
        coords = make_pose_coords_local(sys, pos, up, tloc, phi, center=True)
        coords_frames.append(coords)
        comments.append(f"iu={iu} iv={iv} iy={iy} it={it} ip={ip}  x={pos[0]:.6f} y={pos[1]:.6f} z={pos[2]:.6f}  phi_deg={np.rad2deg(phi):.3f}  tilt=({tloc[0]:.6f},{tloc[1]:.6f},{tloc[2]:.6f})  up=({up[0]:.6f},{up[1]:.6f},{up[2]:.6f})")

    save_xyz_frames(sys, coords_frames, comments, out_xyz)
    print(f"Saved {len(coords_frames)} frames -> {out_xyz}")


def do_counts(args):
    atoms_pair, sdfpack, gridpack = build_spatial_grid(args)
    grid_uv_na, grid_uv_cl, u_vals, u_idx_valid, arc_lengths, x_left, x_right = gridpack
    y_vals = make_y_levels(args.y_min, args.y_max, args.y_levels)
    grid_uv_y, w_na, w_cl = blend_uv_grids(grid_uv_na, grid_uv_cl, y_vals, mode='cos')
    angpack = build_angular_grid(args)
    tilts = angpack[6]
    phi_grid = angpack[7]

    n_y, n_u, n_v = grid_uv_y.shape[:3]
    n_y = len(y_vals)
    n_tilt = len(tilts)
    n_phi = len(phi_grid)

    coarse = _parse_coarse(args, n_u, n_v, n_y, n_tilt, n_phi)
    coarse_counts = {ax: len(coarse[ax]) for ax in AXES}
    coarse_prod = np.prod(list(coarse_counts.values()))

    axis_sizes = {'u': n_u, 'v': n_v, 'y': n_y, 'tilt': n_tilt, 'phi': n_phi}
    scan_totals = {}
    scan_breakdown = {}
    for ax in AXES:
        mult = axis_sizes[ax]
        parts = [f"{ax}={axis_sizes[ax]}"]
        for other in AXES:
            if other == ax:
                continue
            mult *= coarse_counts[other]
            parts.append(f"{other}_coarse={coarse_counts[other]}")
        scan_totals[ax] = mult
        scan_breakdown[ax] = " * ".join(parts)

    print('--- sizes ---')
    print(f"n_u(valid)  = {n_u}")
    print(f"n_v         = {n_v}")
    print(f"n_y         = {n_y}")
    print(f"n_tilt      = {n_tilt}")
    print(f"n_phi       = {n_phi}")
    print('--- coarse grid (defaults or provided) ---')
    print(f"coarse_u={coarse['u']} ({coarse_counts['u']})")
    print(f"coarse_v={coarse['v']} ({coarse_counts['v']})")
    print(f"coarse_y={coarse['y']} ({coarse_counts['y']})")
    print(f"coarse_tilt={coarse['tilt']} ({coarse_counts['tilt']})")
    print(f"coarse_phi={coarse['phi']} ({coarse_counts['phi']})")
    print(f"coarse_product = {coarse_prod}")
    print('--- totals ---')
    print(f"N_total_5D  = {count_5d(n_u, n_v, n_y, n_tilt, n_phi)}")
    print('--- per-axis 1D scan totals (full axis × coarse others) ---')
    for ax in AXES:
        print(f"scan_{ax}: {scan_totals[ax]}   ( {scan_breakdown[ax]} )")
    print(f"scan_all_sum = {sum(scan_totals.values())}")


def do_random(args):
    if args.Nrandom <= 0:
        raise ValueError('--Nrandom must be > 0 for mode=random')

    atoms_pair, sdfpack, gridpack = build_spatial_grid(args)
    grid_uv_na, grid_uv_cl, u_vals, u_idx_valid, arc_lengths, x_left, x_right = gridpack
    y_vals = make_y_levels(args.y_min, args.y_max, args.y_levels)
    grid_uv = blend_uv_grids(grid_uv_na, grid_uv_cl, y_vals, mode='cos')[0]
    angpack = build_angular_grid(args)
    tilts = angpack[6]
    phi_grid = angpack[7]

    n_y, n_u, n_v = grid_uv.shape[:3]
    n_tilt = len(tilts)
    n_phi = len(phi_grid)

    sys = AtomicSystem(fname=args.xyz)

    rng = np.random.default_rng(args.seed)
    indices_iter = _iter_random_indices(rng, n_u, n_v, n_y, n_tilt, n_phi, args.Nrandom)
    _write_xyz_for_indices(f"random_N{args.Nrandom}", indices_iter, sys, grid_uv, y_vals, tilts, phi_grid, args.out_prefix)


def do_scan(args):
    scan_axes = _parse_scan_axes(args)

    atoms_pair, sdfpack, gridpack = build_spatial_grid(args)
    grid_uv_na, grid_uv_cl, u_vals, u_idx_valid, arc_lengths, x_left, x_right = gridpack
    y_vals = make_y_levels(args.y_min, args.y_max, args.y_levels)
    grid_uv = blend_uv_grids(grid_uv_na, grid_uv_cl, y_vals, mode='cos')[0]
    angpack = build_angular_grid(args)
    tilts = angpack[6]
    phi_grid = angpack[7]

    n_y, n_u, n_v = grid_uv.shape[:3]
    n_tilt = len(tilts)
    n_phi = len(phi_grid)

    coarse = _parse_coarse(args, n_u, n_v, n_y, n_tilt, n_phi)

    axis_sizes = {'u': n_u, 'v': n_v, 'y': n_y, 'tilt': n_tilt, 'phi': n_phi}
    coarse_counts = {ax: len(coarse[ax]) for ax in AXES}

    estimate_scan_total(axis_sizes, scan_axes, coarse_counts, args.Ncap)

    out_png = args.png or f"{args.out_prefix}_scan_lines.png"
    _save_scan_png(args, grid_uv, u_vals, y_vals, tilts, phi_grid, coarse, scan_axes, out_png)
    print(f"Saved scan visualization -> {out_png}")

    out_up = f"{args.out_prefix}_up_vectors.png"
    plot_up_vectors(grid_uv, out_up)
    print(f"Saved up-vector debug plot -> {out_up}")

    sys = AtomicSystem(fname=args.xyz)

    def fixed_dict(u=None, v=None, y=None, tilt=None, phi=None):
        return {'u': u, 'v': v, 'y': y, 'tilt': tilt, 'phi': phi}

    if 'u' in scan_axes:
        for iv in coarse['v']:
            for iy in coarse['y']:
                for it in coarse['tilt']:
                    for ip in coarse['phi']:
                        _write_xyz_line('u', range(n_u), fixed_dict(u=None, v=iv, y=iy, tilt=it, phi=ip), sys, grid_uv, y_vals, tilts, phi_grid, args.out_prefix)

    if 'v' in scan_axes:
        for iu in coarse['u']:
            for iy in coarse['y']:
                for it in coarse['tilt']:
                    for ip in coarse['phi']:
                        _write_xyz_line('v', range(n_v), fixed_dict(u=iu, v=None, y=iy, tilt=it, phi=ip), sys, grid_uv, y_vals, tilts, phi_grid, args.out_prefix)

    if 'y' in scan_axes:
        for iu in coarse['u']:
            for iv in coarse['v']:
                for it in coarse['tilt']:
                    for ip in coarse['phi']:
                        _write_xyz_line('y', range(n_y), fixed_dict(u=iu, v=iv, y=None, tilt=it, phi=ip), sys, grid_uv, y_vals, tilts, phi_grid, args.out_prefix)

    if 'tilt' in scan_axes:
        for iu in coarse['u']:
            for iv in coarse['v']:
                for iy in coarse['y']:
                    for ip in coarse['phi']:
                        _write_xyz_line('tilt', range(n_tilt), fixed_dict(u=iu, v=iv, y=iy, tilt=None, phi=ip), sys, grid_uv, y_vals, tilts, phi_grid, args.out_prefix)

    if 'phi' in scan_axes:
        for iu in coarse['u']:
            for iv in coarse['v']:
                for iy in coarse['y']:
                    for it in coarse['tilt']:
                        _write_xyz_line('phi', range(n_phi), fixed_dict(u=iu, v=iv, y=iy, tilt=it, phi=None), sys, grid_uv, y_vals, tilts, phi_grid, args.out_prefix)


def main(argv=None):
    args = parse_args(argv)

    if args.mode == 'counts':
        do_counts(args)
    elif args.mode == 'random':
        do_random(args)
    elif args.mode == 'scan':
        do_scan(args)
    else:
        raise ValueError(f'Unknown mode {args.mode}')


if __name__ == '__main__':
    main()
