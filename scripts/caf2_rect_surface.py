#!/usr/bin/env python3

import os, sys, argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if BASE not in sys.path:
    sys.path.append(BASE)

from pyBall.AtomicSystem import AtomicSystem


COLORS = {
    'Ca': '#1f77b4',
    'F':  '#2ca02c',
}


def norm(v):
    v = np.array(v, dtype=float)
    n = np.linalg.norm(v)
    if n <= 0.0:
        raise ValueError(f'norm(): zero-length vector {v}')
    return v / n, n


def frac_to_cart(frac, lvec):
    return np.asarray(frac, dtype=float) @ np.asarray(lvec, dtype=float)


def cart_to_frac(cart, lvec):
    return np.asarray(cart, dtype=float) @ np.linalg.inv(np.asarray(lvec, dtype=float))


def build_rotation(a, b, c):
    ex, _ = norm(a)
    bz = b - ex * np.dot(b, ex)
    ey, by = norm(bz)
    if by <= 1e-8:
        raise ValueError(f'build_rotation(): in-plane vectors are collinear a={a} b={b}')
    ez = np.cross(ex, ey)
    ez, _ = norm(ez)
    if np.dot(ez, c) < 0.0:
        ey *= -1.0
        ez *= -1.0
    return np.vstack([ex, ey, ez])


def find_rectangular_supercell(a, b, max_coeff=4, ortho_tol=1e-6):
    a = np.array(a, dtype=float)
    b = np.array(b, dtype=float)
    ua, la = norm(a)
    best = None
    for m in range(-max_coeff, max_coeff + 1):
        for n in range(-max_coeff, max_coeff + 1):
            if m == 0 and n == 0:
                continue
            if abs(n) < 1:
                continue
            v = m * a + n * b
            lv = np.linalg.norm(v)
            if lv < 1e-8:
                continue
            area = np.linalg.norm(np.cross(a, v))
            if area < 1e-8:
                continue
            perp = abs(np.dot(ua, v)) / lv
            score = (perp, abs(m) + abs(n), lv)
            if (best is None) or (score < best[0]):
                best = (score, m, n, v)
    if best is None:
        raise ValueError('find_rectangular_supercell(): failed to find non-collinear integer combination of a,b')
    perp, m, n, v = best[0][0], best[1], best[2], best[3]
    if perp > ortho_tol:
        raise ValueError(f'find_rectangular_supercell(): best candidate not orthogonal enough perp={perp} coeff=({m},{n}) vec={v}')
    M = np.array([[1.0, 0.0, 0.0], [float(m), float(n), 0.0], [0.0, 0.0, 1.0]], dtype=float)
    return {'A': a.copy(), 'B': v.copy(), 'coeff_B': (m, n), 'area_scale': abs(int(n)), 'M': M}


def deduplicate_atoms(apos, enames, atypes, qs=None, tol=1e-5):
    keys = {}
    keep = []
    for i, (p, e) in enumerate(zip(apos, enames)):
        k = (e, int(np.round(p[0] / tol)), int(np.round(p[1] / tol)), int(np.round(p[2] / tol)))
        if k in keys:
            continue
        keys[k] = i
        keep.append(i)
    keep = np.array(keep, dtype=np.int32)
    out = {
        'apos': apos[keep].copy(),
        'enames': [enames[i] for i in keep],
        'atypes': atypes[keep].copy() if atypes is not None else None,
        'qs': qs[keep].copy() if qs is not None else None,
    }
    return out


def largest_circular_gap_shift(u):
    u = np.asarray(u, dtype=float)
    if len(u) == 0:
        return 0.0, 1.0
    u = np.mod(u, 1.0)
    u = np.sort(u)
    uu = np.concatenate([u, u[:1] + 1.0])
    gaps = uu[1:] - uu[:-1]
    i = int(np.argmax(gaps))
    mid = uu[i] + 0.5 * gaps[i]
    return float(mid % 1.0), float(gaps[i])


def choose_supercell_origin(gfrac):
    sx, gx = largest_circular_gap_shift(gfrac[:, 0])
    sy, gy = largest_circular_gap_shift(gfrac[:, 1])
    return np.array([sx, sy, 0.0], dtype=float), {'gap_x': gx, 'gap_y': gy}


def build_supercell_from_transform(sys0, M, selection_tol=1e-10, dedup_tol=1e-5, bShiftToGap=True):
    if sys0.lvec is None:
        raise ValueError(f'Input file {getattr(sys0, "fname", None)} has no lattice vectors in lvs comment')
    L0 = np.array(sys0.lvec, dtype=float)
    M = np.array(M, dtype=float)
    detM = np.linalg.det(M)
    if abs(detM) < 1e-12:
        raise ValueError(f'build_supercell_from_transform(): singular transform M={M}')
    if detM <= 0.0:
        raise ValueError(f'build_supercell_from_transform(): transform must preserve orientation, det(M)={detM}')
    Ls = M @ L0
    Minv = np.linalg.inv(M)
    apos0 = np.array(sys0.apos, dtype=float)
    f0 = cart_to_frac(apos0, L0)
    enames0 = list(sys0.enames)
    atypes0 = np.array(sys0.atypes, copy=True) if sys0.atypes is not None else None
    qs0 = np.array(sys0.qs, copy=True) if sys0.qs is not None else None

    nxy = max(2, int(np.max(np.abs(M[:2, :2]))) + 2)
    g_all = []
    enames_all = []
    atypes_all = []
    qs_all = []
    for ia, fi in enumerate(f0):
        for ix in range(-nxy, nxy + 1):
            for iy in range(-nxy, nxy + 1):
                sh = np.array([ix, iy, 0.0], dtype=float)
                g = (fi + sh) @ Minv
                if (g[0] >= -selection_tol) and (g[0] < 1.0 - selection_tol) and (g[1] >= -selection_tol) and (g[1] < 1.0 - selection_tol):
                    g_all.append(g)
                    enames_all.append(enames0[ia])
                    if atypes0 is not None:
                        atypes_all.append(atypes0[ia])
                    if qs0 is not None:
                        qs_all.append(qs0[ia])
    if not g_all:
        raise ValueError('build_supercell_from_transform(): no atoms selected into transformed supercell')

    g_all = np.array(g_all, dtype=float)
    apos_all = frac_to_cart(g_all, Ls)
    atypes_all = np.array(atypes_all, dtype=atypes0.dtype) if atypes0 is not None else None
    qs_all = np.array(qs_all, dtype=qs0.dtype) if qs0 is not None else None
    dd = deduplicate_atoms(apos_all, enames_all, atypes_all, qs_all, tol=dedup_tol)

    g = cart_to_frac(dd['apos'], Ls)
    g[:, 0] -= np.floor(g[:, 0])
    g[:, 1] -= np.floor(g[:, 1])
    if bShiftToGap:
        shift, shift_info = choose_supercell_origin(g)
        g[:, 0] = np.mod(g[:, 0] - shift[0], 1.0)
        g[:, 1] = np.mod(g[:, 1] - shift[1], 1.0)
    else:
        shift = np.zeros(3, dtype=float)
        shift_info = {'gap_x': None, 'gap_y': None}
    dd['apos'] = frac_to_cart(g, Ls)
    dd = deduplicate_atoms(dd['apos'], dd['enames'], dd['atypes'], dd['qs'], tol=dedup_tol)
    return dd, Ls, {'M': M, 'detM': detM, 'origin_shift_frac': shift, 'shift_info': shift_info}


def build_rectangular_base(sys0, coeff_search=4, selection_tol=1e-8, dedup_tol=1e-5):
    if sys0.lvec is None:
        raise ValueError(f'Input file {getattr(sys0, "fname", None)} has no lattice vectors in lvs comment')
    lvec0 = np.array(sys0.lvec, dtype=float)
    a, b, c = lvec0[0], lvec0[1], lvec0[2]
    cell = find_rectangular_supercell(a, b, max_coeff=coeff_search)
    dd, lvec_super, info_super = build_supercell_from_transform(sys0, cell['M'], selection_tol=selection_tol, dedup_tol=dedup_tol, bShiftToGap=True)
    A = lvec_super[0].copy()
    B = lvec_super[1].copy()
    C = lvec_super[2].copy()
    R = build_rotation(A, B, C)
    LA = np.linalg.norm(A)
    LB = np.linalg.norm(B)
    LC = abs(np.dot(C, R[2]))
    if LC <= 1e-8:
        raise ValueError(f'build_rectangular_base(): invalid c projection LC={LC}')
    lvec_rect = np.array([[LA, 0.0, 0.0], [0.0, LB, 0.0], [0.0, 0.0, LC]], dtype=float)
    apos_rect = dd['apos'] @ R.T
    out = deduplicate_atoms(apos_rect, dd['enames'], dd['atypes'], dd['qs'], tol=dedup_tol)
    return out, lvec_rect, {'A': A, 'B': B, 'R': R, 'coeff_B': cell['coeff_B'], 'area_scale': cell['area_scale'], 'M': cell['M'], 'detM': info_super['detM'], 'origin_shift_frac': info_super['origin_shift_frac'], 'shift_info': info_super['shift_info']}


def cluster_z_planes(z, tol=0.12):
    order = np.argsort(z)
    planes = []
    cur = []
    last = None
    for i in order:
        zi = float(z[i])
        if (last is None) or (abs(zi - last) <= tol):
            cur.append(i)
        else:
            planes.append(cur)
            cur = [i]
        last = zi
    if cur:
        planes.append(cur)
    centers = np.array([np.mean(z[idx]) for idx in planes], dtype=float)
    counts = np.array([len(idx) for idx in planes], dtype=np.int32)
    return planes, centers, counts


def group_planes_to_layers(centers, gap_scale=1.5):
    if len(centers) == 0:
        raise ValueError('group_planes_to_layers(): no z planes found')
    if len(centers) == 1:
        return [np.array([0], dtype=np.int32)], np.array([], dtype=float)
    gaps = np.diff(centers)
    ref = np.median(gaps)
    if ref <= 0.0:
        raise ValueError(f'group_planes_to_layers(): non-positive reference gap ref={ref} gaps={gaps}')
    breaks = np.where(gaps > (gap_scale * ref))[0]
    layers = []
    i0 = 0
    for ib in breaks:
        layers.append(np.arange(i0, ib + 1, dtype=np.int32))
        i0 = ib + 1
    layers.append(np.arange(i0, len(centers), dtype=np.int32))
    return layers, gaps


def select_layers(apos, enames, atypes, qs, n_keep=None, keep_from='top', ztol=0.12, gap_scale=1.5):
    planes, centers, counts = cluster_z_planes(apos[:, 2], tol=ztol)
    layers, gaps = group_planes_to_layers(centers, gap_scale=gap_scale)
    nlayers = len(layers)
    if n_keep is None:
        i_layers = np.arange(nlayers, dtype=np.int32)
    else:
        if n_keep <= 0 or n_keep > nlayers:
            raise ValueError(f'select_layers(): requested n_keep={n_keep} but detected nlayers={nlayers}')
        if keep_from == 'top':
            i_layers = np.arange(nlayers - n_keep, nlayers, dtype=np.int32)
        elif keep_from == 'bottom':
            i_layers = np.arange(0, n_keep, dtype=np.int32)
        else:
            raise ValueError(f'select_layers(): unsupported keep_from={keep_from}')
    plane_ids = np.concatenate([layers[i] for i in i_layers])
    atom_ids = np.concatenate([np.array(planes[ip], dtype=np.int32) for ip in plane_ids])
    atom_ids.sort()
    out = {
        'apos': apos[atom_ids].copy(),
        'enames': [enames[i] for i in atom_ids],
        'atypes': atypes[atom_ids].copy() if atypes is not None else None,
        'qs': qs[atom_ids].copy() if qs is not None else None,
        'planes': planes,
        'centers': centers,
        'counts': counts,
        'layers': layers,
        'gaps': gaps,
        'kept_layer_ids': i_layers,
        'kept_plane_ids': plane_ids,
    }
    return out


def replicate_rectangular(apos, enames, atypes, qs, lvec, nx=1, ny=1):
    if nx <= 0 or ny <= 0:
        raise ValueError(f'replicate_rectangular(): invalid nx,ny=({nx},{ny})')
    A = np.array(lvec[0], dtype=float)
    B = np.array(lvec[1], dtype=float)
    apos_out = []
    enames_out = []
    atypes_out = []
    qs_out = []
    for ix in range(nx):
        for iy in range(ny):
            sh = ix * A + iy * B
            apos_out.append(apos + sh[None, :])
            enames_out.extend(enames)
            if atypes is not None:
                atypes_out.append(atypes)
            if qs is not None:
                qs_out.append(qs)
    apos_out = np.concatenate(apos_out, axis=0)
    atypes_out = np.concatenate(atypes_out, axis=0) if atypes_out else None
    qs_out = np.concatenate(qs_out, axis=0) if qs_out else None
    lvec_out = np.array(lvec, dtype=float)
    lvec_out[0] *= nx
    lvec_out[1] *= ny
    return apos_out, enames_out, atypes_out, qs_out, lvec_out


def make_output_lvec(lvec_xy, apos, vacuum=0.0):
    zmin = float(np.min(apos[:, 2]))
    zmax = float(np.max(apos[:, 2]))
    apos2 = np.array(apos, copy=True)
    apos2[:, 2] -= zmin
    thickness = float(np.max(apos2[:, 2]) - np.min(apos2[:, 2]))
    lvec_out = np.array(lvec_xy, copy=True)
    lvec_out[2] = np.array([0.0, 0.0, thickness + float(vacuum)], dtype=float)
    return apos2, lvec_out, thickness


def plot_structure(apos, enames, lvec, png_path, title=''):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    ax0, ax1 = axes
    for e in sorted(set(enames)):
        m = np.array([x == e for x in enames], dtype=bool)
        c = COLORS.get(e, '#444444')
        ax0.scatter(apos[m, 0], apos[m, 1], s=20, c=c, label=e, alpha=0.8)
        ax1.scatter(apos[m, 0], apos[m, 2], s=20, c=c, label=e, alpha=0.8)
    Lx = lvec[0, 0]
    Ly = lvec[1, 1]
    Lz = lvec[2, 2]
    ax0.plot([0, Lx, Lx, 0, 0], [0, 0, Ly, Ly, 0], 'k-', lw=1.0)
    ax1.plot([0, Lx, Lx, 0, 0], [0, 0, Lz, Lz, 0], 'k-', lw=1.0)
    ax0.set_xlabel('x [A]')
    ax0.set_ylabel('y [A]')
    ax1.set_xlabel('x [A]')
    ax1.set_ylabel('z [A]')
    ax0.set_aspect('equal', adjustable='box')
    ax1.set_aspect('equal', adjustable='box')
    ax0.legend(loc='best')
    if title:
        fig.suptitle(title)
    fig.savefig(png_path, dpi=150)
    plt.close(fig)


def parse_args():
    ap = argparse.ArgumentParser(description='Build rectangular CaF2 slab from skewed xyz+lvs input, replicate it, cut layers, and export xyz/png.')
    ap.add_argument('--in_xyz', required=True, help='Input xyz file with lvs comment line')
    ap.add_argument('--nx', type=int, default=1, help='Number of repeats along rectangular a direction')
    ap.add_argument('--ny', type=int, default=None, help='Number of repeats along rectangular b direction')
    ap.add_argument('--nz', type=int, default=None, help='Alias for repeats along rectangular second in-plane direction')
    ap.add_argument('--layers', type=int, default=None, help='How many detected z-layer groups to keep; default keeps all')
    ap.add_argument('--layers_from', choices=['top', 'bottom'], default='top', help='Which side of the slab to preserve when trimming layers')
    ap.add_argument('--z_tol', type=float, default=0.12, help='Tolerance for clustering atoms into z-planes [A]')
    ap.add_argument('--gap_scale', type=float, default=1.5, help='Gap threshold multiplier for splitting planes into layer groups')
    ap.add_argument('--vacuum', type=float, default=None, help='Vacuum to add above the slab in output lvec[2,2]; default preserves original vacuum amount')
    ap.add_argument('--out_dir', default=None, help='Directory for generated files; default <input_dir>/generated_rect')
    ap.add_argument('--out_prefix', default=None, help='Prefix for output files; default derived from input and options')
    ap.add_argument('--coeff_search', type=int, default=4, help='Max integer coefficient when searching rectangular in-plane supercell')
    ap.add_argument('--no_charges', action='store_true', help='Do not write charges into the xyz output')
    return ap.parse_args()


def main():
    args = parse_args()
    if (args.ny is not None) and (args.nz is not None) and (args.ny != args.nz):
        raise ValueError(f'Use either --ny or --nz for the second in-plane repeat, or set them equal; got ny={args.ny} nz={args.nz}')
    ny = args.nz if args.nz is not None else (args.ny if args.ny is not None else 1)
    sys0 = AtomicSystem(fname=args.in_xyz, bPreinit=False)
    if sys0.lvec is None:
        raise ValueError(f'Input xyz has no lvs comment line: {args.in_xyz}')

    base, lvec_rect0, info = build_rectangular_base(sys0, coeff_search=args.coeff_search)
    layer_sel = select_layers(base['apos'], base['enames'], base['atypes'], base['qs'], n_keep=args.layers, keep_from=args.layers_from, ztol=args.z_tol, gap_scale=args.gap_scale)

    vac0 = float(np.linalg.norm(sys0.lvec[2]) - (np.max(sys0.apos[:, 2]) - np.min(sys0.apos[:, 2])))
    if args.vacuum is None:
        vacuum = max(0.0, vac0)
    else:
        vacuum = float(args.vacuum)
    if vacuum < 0.0:
        raise ValueError(f'vacuum must be >=0, got {vacuum}')

    apos_trim, lvec_trim, thickness = make_output_lvec(lvec_rect0, layer_sel['apos'], vacuum=vacuum)
    apos_out, enames_out, atypes_out, qs_out, lvec_out = replicate_rectangular(apos_trim, layer_sel['enames'], layer_sel['atypes'], layer_sel['qs'], lvec_trim, nx=args.nx, ny=ny)

    out_dir = args.out_dir or os.path.join(os.path.dirname(os.path.abspath(args.in_xyz)), 'generated_rect')
    os.makedirs(out_dir, exist_ok=True)
    stem = os.path.splitext(os.path.basename(args.in_xyz))[0]
    prefix = args.out_prefix or f'{stem}_rect_nx{args.nx}_nz{ny}_L{args.layers if args.layers is not None else len(layer_sel["layers"])}_{args.layers_from}'
    xyz_path = os.path.join(out_dir, prefix + '.xyz')
    png_path = os.path.join(out_dir, prefix + '.png')

    sys_out = AtomicSystem(apos=apos_out, atypes=atypes_out, enames=enames_out, lvec=lvec_out, qs=qs_out, bPreinit=False)
    sys_out.saveXYZ(xyz_path, blvec=True, bQs=(not args.no_charges))
    plot_structure(apos_out, enames_out, lvec_out, png_path, title=prefix)

    print('Input        :', args.in_xyz)
    print('Rect coeff_B :', info['coeff_B'])
    print('Base atoms   :', len(base['apos']))
    print('Planes       :', len(layer_sel['planes']), 'counts=', layer_sel['counts'].tolist())
    print('Layer groups :', len(layer_sel['layers']), 'kept=', layer_sel['kept_layer_ids'].tolist())
    print('Supercell    :', f'nx={args.nx} nz={ny}')
    print('Vacuum [A]   :', vacuum)
    print('Thickness[A] :', thickness)
    print('Output xyz   :', xyz_path)
    print('Output png   :', png_path)


if __name__ == '__main__':
    main()
