#!/usr/bin/env python3

import os
import sys
import json
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
if ROOT not in sys.path:
    sys.path.append(ROOT)

from pyBall.SurfaceSampling import ProbeParticle, build_backend, sample_surface_map
from pyBall.VispyUtils import make_surface_mesh, render_surface_png
from pyBall.OCL.InteractionEnergy import load_xyz_with_REQs


def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)


def parse_type_map(txt):
    if not txt:
        return {}
    out = {}
    for item in txt.split(','):
        item = item.strip()
        if not item:
            continue
        if '=' not in item:
            raise ValueError(f"parse_type_map(): item '{item}' missing '='")
        k, v = item.split('=', 1)
        out[k.strip()] = v.strip()
    return out


def finite_range(arr):
    a = np.asarray(arr, dtype=np.float64)
    m = np.isfinite(a)
    if not np.any(m):
        raise ValueError('finite_range(): no finite values')
    return float(np.min(a[m])), float(np.max(a[m]))


def seam_metrics(arr):
    a = np.asarray(arr, dtype=np.float64)
    if a.ndim != 2:
        raise ValueError(f'seam_metrics(): expected 2D array, got {a.shape}')
    if not np.isfinite(a).all():
        raise ValueError('seam_metrics(): array contains non-finite values')
    return {
        'seam_x_max_abs': float(np.max(np.abs(a[0, :] - a[-1, :]))),
        'seam_y_max_abs': float(np.max(np.abs(a[:, 0] - a[:, -1]))),
    }


def plot_heightmaps(out_png, result, title):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    z = result.zs.T
    c = result.color.T
    ext = [float(result.xs[0]), float(result.xs[-1]), float(result.ys[0]), float(result.ys[-1])]
    im0 = axes[0].imshow(z, origin='lower', extent=ext, aspect='equal')
    axes[0].set_title('z_node(x,y) [A]')
    axes[0].set_xlabel('x [A]')
    axes[0].set_ylabel('y [A]')
    fig.colorbar(im0, ax=axes[0], fraction=0.046)
    im1 = axes[1].imshow(c, origin='lower', extent=ext, aspect='equal', cmap='coolwarm')
    axes[1].set_title(f'{result.color_component}(x,y,z_node) [eV]')
    axes[1].set_xlabel('x [A]')
    axes[1].set_ylabel('y [A]')
    fig.colorbar(im1, ax=axes[1], fraction=0.046)
    fig.suptitle(title)
    fig.savefig(out_png, dpi=160)
    plt.close(fig)


# keep separate for future JSON export stability

def build_title(args, probe, result):
    src = f"gridff:{os.path.basename(args.gridff)}" if args.source == 'gridff' else f"xyz:{os.path.basename(args.sub)}"
    return f"{src} | probe={probe.name} q={probe.charge:+.3f} | mode={result.diagnostics['mode']} | selector={result.selector} | color={result.color_component}"


def main():
    here = os.path.abspath(os.path.dirname(__file__))
    if os.getcwd() != here:
        raise RuntimeError(f"Run this script from tests/tMMFF, current cwd={os.getcwd()}, expected={here}")
    ap = argparse.ArgumentParser(description='Render isoheight surface z_node(x,y) for probe/substrate interaction and save static PNGs.')
    ap.add_argument('--source', choices=['gridff', 'xyz'], required=True, help='Potential source type')
    ap.add_argument('--gridff', type=str, default=None, help='Path to Bspline_PLQd.npy when --source=gridff')
    ap.add_argument('--gridff_xyz', type=str, default=None, help='Optional xyz corresponding to GridFF (if not inferable)')
    ap.add_argument('--sub', type=str, default=None, help='Substrate xyz path when --source=xyz')
    ap.add_argument('--xyz_backend', choices=['reference', 'fast_gpu', 'folded_gpu'], default='reference', help='XYZ substrate evaluation backend')
    ap.add_argument('--sub_types', type=str, default='', help='Type map e.g. Ca=Ca+2')
    ap.add_argument('--probe', type=str, default='H', help='Probe atom type from AtomTypes.dat')
    ap.add_argument('--charge', type=float, default=0.2, help='Probe charge [e]')
    ap.add_argument('--alpha', type=float, default=2.0, help='Probe Morse alpha')
    ap.add_argument('--E0', type=float, default=None, help='Optional probe EvdW override')
    ap.add_argument('--R0', type=float, default=None, help='Optional probe RvdW override')
    ap.add_argument('--selector', choices=['total', 'nonel', 'coulomb', 'pauli', 'london'], default='total', help='Component used for z-node search')
    ap.add_argument('--color', choices=['total', 'nonel', 'coulomb', 'pauli', 'london'], default='coulomb', help='Component used for coloring')
    ap.add_argument('--mode', choices=['threshold', 'first_minimum'], default='threshold', help='How z_node is selected')
    ap.add_argument('--threshold', type=float, default=0.0, help='Threshold for mode=threshold [eV]')
    ap.add_argument('--xrange', type=str, default='-1.0,5.0', help='x min,max [A]')
    ap.add_argument('--yrange', type=str, default='-1.0,5.0', help='y min,max [A]')
    ap.add_argument('--zrange', type=str, default='1.0,8.0', help='z min,max [A]')
    ap.add_argument('--nx', type=int, default=64, help='Number of x samples')
    ap.add_argument('--ny', type=int, default=64, help='Number of y samples')
    ap.add_argument('--nz', type=int, default=80, help='Number of z samples per xy point')
    ap.add_argument('--outdir', type=str, default='output_surface_iso', help='Output directory under tests/tMMFF')
    ap.add_argument('--name', type=str, default=None, help='Optional output basename prefix')
    ap.add_argument('--strict', action='store_true', help='Fail loudly if any xy point has no valid z node')
    args = ap.parse_args()

    def parse_range(txt, name):
        vals = [float(x.strip()) for x in txt.split(',') if x.strip()]
        if len(vals) != 2:
            raise ValueError(f"{name}: expected 'min,max', got '{txt}'")
        if vals[1] <= vals[0]:
            raise ValueError(f"{name}: max must be > min, got {vals}")
        return tuple(vals)

    xr = parse_range(args.xrange, 'xrange')
    yr = parse_range(args.yrange, 'yrange')
    zr = parse_range(args.zrange, 'zrange')
    outdir = os.path.join(here, args.outdir)
    ensure_dir(outdir)

    probe = ProbeParticle(name=args.probe, charge=args.charge, alpha=args.alpha, E0=args.E0, R0=args.R0)
    backend = build_backend(
        source=args.source,
        substrate_xyz=(os.path.abspath(args.sub) if args.sub else None),
        gridff_path=(os.path.abspath(args.gridff) if args.gridff else None),
        probe=probe,
        xyz_type_map=parse_type_map(args.sub_types),
        backend=args.xyz_backend,
        xyz_path=(os.path.abspath(args.gridff_xyz) if args.gridff_xyz else None),
    )
    result = sample_surface_map(backend, x_range=xr, y_range=yr, z_range=zr, nx=args.nx, ny=args.ny, nz=args.nz, selector=args.selector, color_component=args.color, mode=args.mode, threshold=args.threshold, bFailOnMiss=args.strict)
    if not np.isfinite(result.zs[result.ok_mask]).all():
        raise RuntimeError('Non-finite z values in valid surface points')
    if not np.isfinite(result.color[result.ok_mask]).all():
        raise RuntimeError('Non-finite color values in valid surface points')
    if np.any(~result.ok_mask):
        raise RuntimeError(f'Incomplete surface: missing {int(np.size(result.ok_mask) - np.sum(result.ok_mask))} of {int(result.ok_mask.size)} points')
    seam_z = seam_metrics(result.zs)
    seam_c = seam_metrics(result.color)
    base = args.name
    if base is None:
        src = os.path.splitext(os.path.basename(args.gridff if args.source == 'gridff' else args.sub))[0]
        src_tag = f"{args.source}_{args.xyz_backend}" if args.source == 'xyz' else 'gridff'
        base = f"{src}_{src_tag}_{args.probe}_q{args.charge:+.2f}_{args.mode}_{args.selector}_color_{args.color}".replace('+', 'p').replace('-', 'm')
    title = build_title(args, probe, result)
    mesh = make_surface_mesh(result.xs, result.ys, result.zs, scalar=result.color, cmap='coolwarm', symmetric=False, mask=result.ok_mask)

    png_surface = os.path.join(outdir, base + '_surface.png')
    png_maps = os.path.join(outdir, base + '_maps.png')
    npz_out = os.path.join(outdir, base + '.npz')
    json_out = os.path.join(outdir, base + '.json')

    atom_points = None
    if hasattr(backend, 'sub_apos'):
        atom_points = np.asarray(backend.sub_apos, dtype=np.float32)
    render_surface_png(png_surface, mesh, atom_points=atom_points, title=title)
    plot_heightmaps(png_maps, result, title)
    np.savez(npz_out, xs=result.xs, ys=result.ys, zs=result.zs, color=result.color, ok_mask=result.ok_mask.astype(np.uint8), points=result.points)
    zmin, zmax = finite_range(result.zs)
    cmin, cmax = finite_range(result.color)
    meta = {
        'title': title,
        'probe': probe.as_dict(),
        'source': args.source,
        'xyz_backend': args.xyz_backend,
        'selector': args.selector,
        'color': args.color,
        'mode': args.mode,
        'threshold': float(args.threshold),
        'xrange': xr,
        'yrange': yr,
        'zrange': zr,
        'shape': [int(args.nx), int(args.ny), int(args.nz)],
        'n_ok': int(np.sum(result.ok_mask)),
        'n_total': int(result.ok_mask.size),
        'z_node_range': [zmin, zmax],
        'color_range': [cmin, cmax],
        'seam_z': seam_z,
        'seam_color': seam_c,
        'fail_points': result.diagnostics['fail_points'][:50],
        'surface_png': png_surface,
        'maps_png': png_maps,
        'npz': npz_out,
    }
    with open(json_out, 'w') as f:
        json.dump(meta, f, indent=2)
    print('surface_png:', png_surface)
    print('maps_png   :', png_maps)
    print('npz        :', npz_out)
    print('json       :', json_out)
    print('n_ok/n_tot :', meta['n_ok'], '/', meta['n_total'])
    print('z range    :', meta['z_node_range'])
    print('color range:', meta['color_range'])
    print('seam z     :', seam_z)
    print('seam color :', seam_c)


if __name__ == '__main__':
    main()
