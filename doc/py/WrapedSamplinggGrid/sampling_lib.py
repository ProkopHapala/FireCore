import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.interpolate import interp1d


# Make sure repo root is on sys.path for pyBall import
HERE = os.path.dirname(__file__)
REPO_ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
for p in (REPO_ROOT, os.path.join(REPO_ROOT, 'pyBall')):
    if p not in sys.path:
        sys.path.insert(0, p)

from pyBall.AtomicSystem import AtomicSystem


# ============================== SDF sampling ==============================

def make_atoms_nacl_step(a, r_Na, r_Cl, corner_species='Na'):
    def sp(ix, iz):
        p = (ix + iz + 1) % 2
        if corner_species == 'Na':
            return ('Na', r_Na) if p == 0 else ('Cl', r_Cl)
        else:
            return ('Cl', r_Cl) if p == 0 else ('Na', r_Na)

    atoms = []
    for ix, iz in [(0, 1), (0, 0), (1, 0), (2, 0)]:
        S, r = sp(ix, iz)
        atoms.append({'x': ix*a, 'z': iz*a, 'r': r, 'species': S, 'ix': ix, 'iz': iz})
    return atoms


def compute_sdf_2d(atoms, x_range, z_range, res, combine_p):
    x_lin = np.linspace(x_range[0], x_range[1], res)
    z_lin = np.linspace(z_range[0], z_range[1], res)
    X, Z = np.meshgrid(x_lin, z_lin)

    dists = []
    for at in atoms:
        di = np.sqrt((X - at['x'])**2 + (Z - at['z'])**2) - at['r']
        dists.append(di)

    stack = np.stack(dists, axis=0)

    if np.isinf(combine_p):
        sdf = np.min(stack, axis=0)
    else:
        k = max(float(combine_p), 1e-6)
        m = np.min(stack, axis=0)
        soft = -np.log(np.sum(np.exp(-k * (stack - m)), axis=0)) / k + m
        sdf = soft

    return X, Z, sdf, x_lin, z_lin


def make_u_levels(u_min, u_max, n, schedule, sinh_a):
    t = np.linspace(0, 1, n)
    if schedule == 'linear':
        return u_min + t * (u_max - u_min)
    elif schedule == 'sinh':
        return u_min + (u_max - u_min) * np.sinh(sinh_a * t) / np.sinh(sinh_a)
    elif schedule == 'geom':
        return np.geomspace(u_min, u_max, n)
    raise ValueError(f"Unknown schedule '{schedule}'")


def parse_u_list(u_list_str):
    vals = [float(s) for s in u_list_str.split(',') if s.strip() != '']
    vals = np.array(vals, dtype=float)
    vals = np.unique(vals)
    vals.sort()
    return vals


def walk_isoline_segment(sdf, x_lin, z_lin, u_level, x_left, x_right, walk_res=0.005):
    X, Z = np.meshgrid(x_lin, z_lin)
    fig_tmp, ax_tmp = plt.subplots()
    cs = ax_tmp.contour(X, Z, sdf, levels=[u_level])
    plt.close(fig_tmp)

    best_seg = None
    best_n = 0

    for path in cs.get_paths():
        verts = path.vertices
        if len(verts) < 4:
            continue
        xv, zv = verts[:, 0], verts[:, 1]
        if xv.max() < x_right - 0.05 or xv.min() > x_left + 0.05:
            continue
        if np.median(zv) < -0.2:
            continue
        if len(verts) > best_n:
            best_n = len(verts)
            best_seg = verts

    if best_seg is None:
        return None

    raw_s = np.concatenate([[0], np.cumsum(np.sqrt(np.diff(best_seg[:, 0])**2 + np.diff(best_seg[:, 1])**2))])
    total_raw = raw_s[-1]
    if total_raw < 0.01:
        return None

    n_fine = max(int(total_raw / walk_res) + 1, 4)
    s_fine = np.linspace(0, total_raw, n_fine)

    fx = interp1d(raw_s, best_seg[:, 0], kind='linear', bounds_error=False, fill_value='extrapolate')
    fz = interp1d(raw_s, best_seg[:, 1], kind='linear', bounds_error=False, fill_value='extrapolate')

    x_fine = fx(s_fine)
    z_fine = fz(s_fine)

    def find_crossing_idx(x_arr, x_target, direction='right'):
        crossings = []
        for i in range(len(x_arr) - 1):
            if direction == 'right':
                if x_arr[i] <= x_target <= x_arr[i + 1]:
                    crossings.append(i)
            else:
                if x_arr[i] >= x_target >= x_arr[i + 1]:
                    crossings.append(i)
        return crossings

    cross_L = find_crossing_idx(x_fine, x_left, 'right')
    cross_R = find_crossing_idx(x_fine, x_right, 'right')

    if (not cross_L) or (not cross_R):
        cross_L = find_crossing_idx(x_fine, x_left, 'left')
        cross_R = find_crossing_idx(x_fine, x_right, 'left')

    if (not cross_L) or (not cross_R):
        return None

    best_arc = None
    best_z = -np.inf
    for iL in cross_L:
        for iR in cross_R:
            if iR <= iL:
                continue
            arc_z = z_fine[iL:iR + 1]
            med_z = np.median(arc_z)
            if med_z > best_z:
                best_z = med_z
                best_arc = (iL, iR)

    if best_arc is None:
        return None

    iL, iR = best_arc

    def interp_crossing(i, x_arr, z_arr, x_target):
        if i + 1 >= len(x_arr):
            return x_arr[i], z_arr[i]
        dx = x_arr[i + 1] - x_arr[i]
        if abs(dx) < 1e-12:
            return x_arr[i], z_arr[i]
        t = (x_target - x_arr[i]) / dx
        return x_target, z_arr[i] + t * (z_arr[i + 1] - z_arr[i])

    xL_exact, zL_exact = interp_crossing(iL, x_fine, z_fine, x_left)
    xR_exact, zR_exact = interp_crossing(iR, x_fine, z_fine, x_right)

    x_seg = np.concatenate([[xL_exact], x_fine[iL + 1:iR + 1], [xR_exact]])
    z_seg = np.concatenate([[zL_exact], z_fine[iL + 1:iR + 1], [zR_exact]])
    pts = np.column_stack([x_seg, z_seg])

    diffs = np.diff(pts, axis=0)
    ds = np.sqrt((diffs**2).sum(axis=1))
    s = np.concatenate([[0.0], np.cumsum(ds)])
    return pts, s, s[-1]


def sample_isoline_uniform(pts, s, n_v):
    s_total = s[-1]
    s_samples = np.linspace(0, s_total, n_v)
    fx = interp1d(s, pts[:, 0], kind='linear', bounds_error=False, fill_value=(pts[0, 0], pts[-1, 0]))
    fz = interp1d(s, pts[:, 1], kind='linear', bounds_error=False, fill_value=(pts[0, 1], pts[-1, 1]))
    return np.column_stack([fx(s_samples), fz(s_samples)])


def build_grid_uv(sdf, x_lin, z_lin, u_levels, x_left, x_right, n_v, walk_res):
    rows = []
    u_vals_valid = []
    u_idx_valid = []
    arc_lengths = []

    for i, u in enumerate(u_levels):
        result = walk_isoline_segment(sdf, x_lin, z_lin, u, x_left, x_right, walk_res)
        if result is None:
            print(f"    u_idx={i:2d}  u={u:.4f}Å  — isoline not found")
            continue
        pts, s, s_total = result
        samples = sample_isoline_uniform(pts, s, n_v)
        rows.append(samples)
        u_vals_valid.append(u)
        u_idx_valid.append(i)
        arc_lengths.append(s_total)
        print(f"    u_idx={i:2d}  u={u:.4f}Å  arc={s_total:.3f}Å  ✓")

    if not rows:
        return None, np.array([]), [], []

    grid = np.array(rows)
    return grid, np.array(u_vals_valid), u_idx_valid, arc_lengths


# ============================== y sampling ==============================

def make_y_levels(y_min, y_max, n_y):
    if n_y <= 1:
        return np.array([0.5 * (y_min + y_max)], dtype=float)
    return np.linspace(y_min, y_max, n_y)


def y_blend_weights(y_vals, mode='cos'):
    """Weights for blending Na and Cl configurations across y.

    We interpret y_vals as sampling within one unit cell; endpoints correspond
    to the two extreme configurations.

    Returns:
      w_Na, w_Cl arrays of shape (n_y,)
    """
    y_vals = np.asarray(y_vals, dtype=float)
    if len(y_vals) == 1:
        return np.array([0.5]), np.array([0.5])
    t = (y_vals - y_vals[0]) / (y_vals[-1] - y_vals[0])
    t = np.clip(t, 0.0, 1.0)
    if mode == 'linear':
        w_na = 1.0 - t
    elif mode == 'cos':
        # smooth periodic-ish blend: w=cos^2(pi/2 * t)
        w_na = np.cos(0.5 * np.pi * t) ** 2
    else:
        raise ValueError(f"Unknown y_blend mode '{mode}'")
    w_cl = 1.0 - w_na
    return w_na, w_cl


def blend_uv_grids(grid_uv_na, grid_uv_cl, y_vals, mode='cos'):
    """Blend two (n_u,n_v,2) uv grids into (n_y,n_u,n_v,2) using y-dependent weights."""
    w_na, w_cl = y_blend_weights(y_vals, mode=mode)
    out = []
    for wa, wc in zip(w_na, w_cl):
        out.append(wa * grid_uv_na + wc * grid_uv_cl)
    return np.stack(out, axis=0), w_na, w_cl


# ============================== angular sampling ==============================

def build_fan_points(M, n_sub):
    center = np.array([[0.0, 0.0]])
    angles = np.linspace(0, 2 * np.pi, M + 1)
    outer = np.stack([np.cos(angles), np.sin(angles)], axis=1)[:-1]
    verts_base = np.vstack([center, outer])

    wedges = [(2 * np.pi * i / M, 2 * np.pi * (i + 1) / M) for i in range(M)]

    pts = [center[0]]
    for k in range(1, n_sub + 1):
        r = k / n_sub
        angs = np.linspace(0, 2 * np.pi, M * k, endpoint=False)
        ring = np.stack([r * np.cos(angs), r * np.sin(angs)], axis=1)
        pts.extend(ring)
    sub_points = np.array(pts)
    return verts_base, wedges, sub_points


def map_to_disk(sub_points, wedges, M, n_sub):
    disk_pts = [np.array([0.0, 0.0])]
    ring_ids = [0.0]
    for k in range(1, n_sub + 1):
        r = k / n_sub
        angs = np.linspace(0, 2 * np.pi, M * k, endpoint=False)
        ring = np.stack([r * np.cos(angs), r * np.sin(angs)], axis=1)
        disk_pts.extend(ring)
        ring_ids.extend([r] * len(ring))
    return np.array(disk_pts), np.array(ring_ids)


def lift_to_dome(disk_pts, theta_max_deg, radius):
    theta_max = np.deg2rad(theta_max_deg)
    if theta_max >= np.pi / 2:
        theta_max = min(theta_max, np.pi / 2 - 1e-6)
    r_xy = np.linalg.norm(disk_pts, axis=1)
    r_xy = np.clip(r_xy, 0.0, 1.0)
    theta = r_xy * theta_max
    phi = np.arctan2(disk_pts[:, 1], disk_pts[:, 0])
    x = radius * np.sin(theta) * np.cos(phi)
    y = radius * np.sin(theta) * np.sin(phi)
    z = radius * np.cos(theta)
    return np.stack([x, y, z], axis=1)


def sample_phi(phi_samples, symmetry_m):
    phi_max = 2 * np.pi / symmetry_m if symmetry_m > 0 else 2 * np.pi
    return np.linspace(0, phi_max, phi_samples, endpoint=False)


def generate_tilts(fans, subdiv, theta_max, radius, phi_samples, symmetry_m, phi_max_deg=None):
    verts_base, wedges, sub_points = build_fan_points(fans, subdiv)
    disk_pts, ring_ids = map_to_disk(sub_points, wedges, fans, subdiv)
    dome_pts = lift_to_dome(disk_pts, theta_max, radius)
    tilts = dome_pts / np.linalg.norm(dome_pts, axis=1, keepdims=True)
    if phi_max_deg is not None:
        phi_range = np.deg2rad(phi_max_deg)
        phi_grid = np.linspace(0, phi_range, phi_samples, endpoint=False)
    else:
        phi_grid = sample_phi(phi_samples, symmetry_m)
    return verts_base, wedges, sub_points, disk_pts, ring_ids, dome_pts, tilts, phi_grid


# ============================== orientation + XYZ ==============================

def rotation_matrix(axis, angle):
    axis = np.asarray(axis, dtype=float)
    norm = np.linalg.norm(axis)
    if norm < 1e-12:
        return np.eye(3)
    axis = axis / norm
    x, y, z = axis
    c = np.cos(angle)
    s = np.sin(angle)
    C = 1 - c
    return np.array([
        [c + x * x * C, x * y * C - z * s, x * z * C + y * s],
        [y * x * C + z * s, c + y * y * C, y * z * C - x * s],
        [z * x * C - y * s, z * y * C + x * s, c + z * z * C],
    ])


def align_z_to(target):
    target = np.asarray(target, dtype=float)
    n = np.linalg.norm(target)
    if n < 1e-12:
        return np.eye(3)
    t = target / n
    z = np.array([0.0, 0.0, 1.0])
    if np.linalg.norm(t - z) < 1e-12:
        return np.eye(3)
    if np.linalg.norm(t + z) < 1e-12:
        return rotation_matrix([1, 0, 0], np.pi)
    axis = np.cross(z, t)
    angle = np.arccos(np.clip(np.dot(z, t), -1.0, 1.0))
    return rotation_matrix(axis, angle)


def _ensure_dir_for_file(path):
    d = os.path.dirname(os.path.abspath(path))
    if d and (not os.path.exists(d)):
        os.makedirs(d, exist_ok=True)


def save_xyz_frames(sys: AtomicSystem, coords_frames, comments, path):
    sys.ensure_numpy_arrays()
    enames = sys.enames
    _ensure_dir_for_file(path)
    with open(path, 'w') as f:
        for coords, cmt in zip(coords_frames, comments):
            f.write(f"{len(enames)}\n")
            f.write(str(cmt) + "\n")
            for e, p in zip(enames, coords):
                f.write(f"{e}  {p[0]: .6f}  {p[1]: .6f}  {p[2]: .6f}\n")


def make_oriented_coords(sys: AtomicSystem, tilt_vec, phi, center=True):
    sys.ensure_numpy_arrays()
    coords0 = sys.apos.copy()
    if center:
        coords0 = coords0 - coords0.mean(axis=0)
    tvec = np.asarray(tilt_vec, dtype=float)
    R_align = align_z_to(tvec)
    R_phi = rotation_matrix(tvec, phi)
    R = R_phi @ R_align
    rot_coords = (R @ coords0.T).T
    return rot_coords


def make_pose_coords(sys: AtomicSystem, pos_xyz, tilt_vec, phi, center=True):
    rot_coords = make_oriented_coords(sys, tilt_vec, phi, center=center)
    return rot_coords + np.asarray(pos_xyz, dtype=float)[None, :]


def local_up_from_streamline(grid_uv, iu, iv):
    """Local 'up' direction from streamline connector at fixed v.

    We use forward/backward difference along u-index.
    Returns normalized 3D vector (x,0,z).
    """
    n_u = grid_uv.shape[0]
    if n_u < 2:
        return np.array([0.0, 0.0, 1.0])
    if iu <= 0:
        d = grid_uv[1, iv] - grid_uv[0, iv]
    elif iu >= n_u - 1:
        d = grid_uv[n_u - 1, iv] - grid_uv[n_u - 2, iv]
    else:
        d = grid_uv[iu + 1, iv] - grid_uv[iu - 1, iv]
    v = np.array([d[0], 0.0, d[1]], dtype=float)
    n = np.linalg.norm(v)
    if n < 1e-12:
        return np.array([0.0, 0.0, 1.0])
    return v / n


def _orthonormal_frame_from_up(up, ref=np.array([0.0, 1.0, 0.0])):
    up = np.asarray(up, dtype=float)
    up = up / max(np.linalg.norm(up), 1e-12)
    e1 = np.cross(ref, up)
    n1 = np.linalg.norm(e1)
    if n1 < 1e-12:
        ref2 = np.array([1.0, 0.0, 0.0])
        e1 = np.cross(ref2, up)
        n1 = np.linalg.norm(e1)
    e1 = e1 / max(n1, 1e-12)
    e2 = np.cross(up, e1)
    e2 = e2 / max(np.linalg.norm(e2), 1e-12)
    # columns are basis vectors
    return np.stack([e1, e2, up], axis=1)


def make_oriented_coords_local(sys: AtomicSystem, up_vec, tilt_vec_local, phi, center=True):
    """Orient molecule with tilt measured from local up_vec.

    tilt_vec_local: unit vector defined in the local frame where +Z is up_vec.
    We map it to global by frame @ tilt_vec_local, then rotate around global tilt axis by phi.
    """
    sys.ensure_numpy_arrays()
    coords0 = sys.apos.copy()
    if center:
        coords0 = coords0 - coords0.mean(axis=0)

    frame = _orthonormal_frame_from_up(up_vec)
    tloc = np.asarray(tilt_vec_local, dtype=float)
    tloc = tloc / max(np.linalg.norm(tloc), 1e-12)
    tglob = frame @ tloc
    # rotate +Z to tglob, then apply phi around tglob
    R_align = align_z_to(tglob)
    R_phi = rotation_matrix(tglob, phi)
    R = R_phi @ R_align
    return (R @ coords0.T).T


def make_pose_coords_local(sys: AtomicSystem, pos_xyz, up_vec, tilt_vec_local, phi, center=True):
    rot_coords = make_oriented_coords_local(sys, up_vec, tilt_vec_local, phi, center=center)
    return rot_coords + np.asarray(pos_xyz, dtype=float)[None, :]


# ============================== index parsing + counting ==============================

def parse_int_list(s):
    if s is None:
        return None
    ss = str(s).strip()
    if ss == '':
        return None
    out = []
    for tok in ss.split(','):
        tok = tok.strip()
        if tok == '':
            continue
        out.append(int(tok))
    return out


def count_5d(n_u, n_v, n_y, n_tilt, n_phi):
    return int(n_u) * int(n_v) * int(n_y) * int(n_tilt) * int(n_phi)


def estimate_scan_total(axis_sizes, scan_axes, coarse_counts_other_axes, ncap_total):
    """Estimate total frames across all requested 1D scans.

    axis_sizes: dict {axis: full_count} for axes in {'u','v','y','tilt','phi'}
    scan_axes: list[str] which axes are scanned (1D)
    coarse_counts_other_axes: dict {axis: coarse_count} for non-scanned axes

    total = Σ_{axis in scan_axes} full_count(axis) * Π_{other!=axis} coarse_count(other)
    """
    total = 0
    for ax in scan_axes:
        mult = axis_sizes[ax]
        for other, cc in coarse_counts_other_axes.items():
            if other == ax:
                continue
            mult *= cc
        total += mult
    if (ncap_total is not None) and (total > ncap_total):
        raise RuntimeError(f"Ncap exceeded: estimated scan total {total} > Ncap {ncap_total}")
    return total
