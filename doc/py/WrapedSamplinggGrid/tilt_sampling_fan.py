"""
Tilt + Phi Sampling via Triangle Fan → Disk → Spherical Dome

Constructs a tilt grid by subdividing a triangle fan (center pole + M fan edges),
then maps the fan to a disk (concentric rings), then lifts to a spherical dome
up to a max tilt theta_max. Also samples phi rotation around the tilt axis with
symmetry reduction (m-fold).

Visualizations (four panels):
1) Raw fan (2D) with M wedges
2) Subdivided fan (2D) with n_sub per edge (rings shown)
3) Disk mapping with concentric rings
4) Dome (3D) points on sphere cap

CLI:
  --fans          M fan edges (e.g., 6 for hex)
  --subdiv        subdivisions per wedge edge (n>=1)
  --theta_max     max tilt in degrees
  --phi_samples   samples around axis
  --symmetry_m    m-fold symmetry (phi range 0..2π/m)
  --radius        sphere radius
  --save          filename to save (PNG); default tilt_sampling.png next to script

"""

import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 3D projection

from sampling_lib import (
    AtomicSystem,
    align_z_to,
    build_fan_points,
    generate_tilts,
    lift_to_dome,
    map_to_disk,
    rotation_matrix,
    sample_phi,
)


def plot_all(verts_base, sub_points, disk_pts, dome_pts, phi_grid, ring_ids, args):
    fig = plt.figure(figsize=(12,9))
    gs = fig.add_gridspec(2,2)

    ax0 = fig.add_subplot(gs[0,0])
    ax0.set_title(f"Fan (M={args.fans})")
    ax0.plot(verts_base[:,0], verts_base[:,1], 'k-o', lw=1)
    ax0.scatter(sub_points[:,0], sub_points[:,1], s=10, c='C1', alpha=0.7)
    for i, p in enumerate(sub_points):
        ax0.text(p[0], p[1], str(i), fontsize=7, color='k')
    # ring guides in fan space: lines for constant (a+b)/n_sub mapped as concentric polygons
    for i in range(1, args.subdiv+1):
        rfrac = i/args.subdiv
        angs = np.linspace(0, 2*np.pi, args.fans+1)
        ax0.plot(rfrac*np.cos(angs), rfrac*np.sin(angs), color='#888888', lw=0.8, ls=':')
    ax0.axis('equal')

    ax1 = fig.add_subplot(gs[0,1])
    ax1.set_title(f"Disk mapping (rings)")
    ax1.scatter(disk_pts[:,0], disk_pts[:,1], s=12, c='C2', alpha=0.7)
    for i, p in enumerate(disk_pts):
        ax1.text(p[0], p[1], str(i), fontsize=7, color='k')
    unique_r = np.unique(ring_ids)
    angs = np.linspace(0, 2*np.pi, 361)
    ring_loops = []
    for r in unique_r:
        ax1.plot(r*np.cos(angs), r*np.sin(angs), color='#999', lw=0.8, ls=':')
        mask = np.isclose(ring_ids, r)
        pts_r = disk_pts[mask]
        if len(pts_r) > 2:
            ang_r = np.arctan2(pts_r[:,1], pts_r[:,0])
            order = np.argsort(ang_r)
            loop = pts_r[order]
            ring_loops.append((r, loop))
            loop_closed = np.vstack([loop, loop[:1]])
            ax1.plot(loop_closed[:,0], loop_closed[:,1], color='#444', lw=1.0, alpha=0.8)
    ax1.axis('equal')

    ax2 = fig.add_subplot(gs[1,0], projection='3d')
    ax2.set_title(f"Dome up to θ_max={args.theta_max}°")
    ax2.scatter(dome_pts[:,0], dome_pts[:,1], dome_pts[:,2], s=12, c=dome_pts[:,2], cmap='viridis')
    ax2.set_xlabel('x'); ax2.set_ylabel('y'); ax2.set_zlabel('z')
    # ring wireframe on dome
    for r, loop in ring_loops:
        mask = np.isclose(ring_ids, r)
        pts_d = dome_pts[mask]
        ang_d = np.arctan2(pts_d[:,1], pts_d[:,0])
        order = np.argsort(ang_d)
        loop_d = pts_d[order]
        loop_d = np.vstack([loop_d, loop_d[:1]])
        ax2.plot(loop_d[:,0], loop_d[:,1], loop_d[:,2], color='#444', lw=1.0, alpha=0.8)
    tilts = dome_pts / np.linalg.norm(dome_pts, axis=1, keepdims=True)
    for i, t in enumerate(tilts):
        ax2.text(t[0], t[1], t[2], str(i), color='k', fontsize=7)
    # equal aspect
    xyz = dome_pts
    max_range = (xyz.max(axis=0) - xyz.min(axis=0)).max()
    mid = xyz.mean(axis=0)
    for i, axis in enumerate([ax2.set_xlim, ax2.set_ylim, ax2.set_zlim]):
        axis(mid[i] - 0.5*max_range, mid[i] + 0.5*max_range)

    ax3 = fig.add_subplot(gs[1,1])
    ax3.set_title(f"Phi samples (m={args.symmetry_m})")
    ax3.scatter(np.cos(phi_grid), np.sin(phi_grid), s=40, c='C3')
    ax3.axis('equal')

    fig.tight_layout()
    out_png = args.save or os.path.join(os.path.dirname(__file__), 'tilt_sampling.png')
    fig.savefig(out_png, dpi=150, bbox_inches='tight')
    print(f"Saved plot to {out_png}")
    plt.show()


# ---------- Movie generation using AtomicSystem ----------
def generate_tilts_args(args):
    return generate_tilts(args.fans, args.subdiv, args.theta_max, args.radius, args.phi_samples, args.symmetry_m, phi_max_deg=args.phi_max)


def save_movie_xyz(sys: AtomicSystem, tilts, phi_grid, args, path):
    sys.ensure_numpy_arrays()
    coords0 = sys.apos - sys.apos.mean(axis=0)  # center
    enames = sys.enames
    tilt_indices = range(len(tilts))
    # ordering: phi_fast -> tilt outer, phi inner; else phi outer, tilt inner
    if args.phi_fast:
        order = [(ti, phi) for ti in tilt_indices for phi in phi_grid]
    else:
        order = [(ti, phi) for phi in phi_grid for ti in tilt_indices]
    with open(path, 'w') as f:
        for ti, phi in order:
            tvec = tilts[ti]
            R_align = align_z_to(tvec)
            R_phi = rotation_matrix(tvec, phi)
            R = R_phi @ R_align
            rot_coords = (R @ coords0.T).T
            f.write(f"{len(enames)}\n")
            f.write(f"tilt_idx={ti} phi_deg={np.rad2deg(phi):.2f} tilt=({tvec[0]:.3f},{tvec[1]:.3f},{tvec[2]:.3f})\n")
            for e, p in zip(enames, rot_coords):
                f.write(f"{e}  {p[0]: .6f}  {p[1]: .6f}  {p[2]: .6f}\n")
    print(f"Saved XYZ movie with {len(order)} frames → {path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fans',        type=int,   default=4, help='number of fan edges (wedges)')
    ap.add_argument('--subdiv',      type=int,   default=2, help='subdivisions per wedge edge (>=1)')
    ap.add_argument('--theta_max',   type=float, default=45.0, help='max tilt angle (deg)')
    ap.add_argument('--radius',      type=float, default=1.0, help='sphere radius')
    ap.add_argument('--phi_samples', type=int,   default=1, help='samples around axis')
    ap.add_argument('--symmetry_m',  type=int,   default=1, help='m-fold symmetry; phi range 0..2π/m')
    ap.add_argument('--save',        type=str,   default=None, help='output PNG path')
    ap.add_argument('--movie_xyz',   type=str,   default='orient_movie.xyz', help='output XYZ movie path; if set, generates orientations')
    ap.add_argument('--xyz',         type=str,   default='../../../cpp/common_resources/xyz/H2O_O.xyz', help='input XYZ for movie (uses AtomicSystem)')
    ap.add_argument('--phi_max',     type=float, default=360.0, help='phi range in degrees (e.g., 180 for 2-fold symmetry)')
    ap.add_argument('--phi_fast',    action='store_true', help='iterate phi fastest (tilt outer loop)')
    ap.add_argument('--plot_tilts',  type=str,   default=None, help='optional tilt plot PNG with indices')
    args = ap.parse_args()

    verts_base, wedges, sub_points, disk_pts, ring_ids, dome_pts, tilts, phi_grid = generate_tilts_args(args)
    if args.plot_tilts:
        # simple tilt plot with indices
        fig = plt.figure(figsize=(5,5))
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(tilts[:,0], tilts[:,1], tilts[:,2], s=30, c='C0')
        for i, t in enumerate(tilts):
            ax.text(t[0], t[1], t[2], str(i), color='k', fontsize=8)
        xyz = tilts
        max_range = (xyz.max(axis=0) - xyz.min(axis=0)).max()
        mid = xyz.mean(axis=0)
        for setter, m in zip([ax.set_xlim, ax.set_ylim, ax.set_zlim], mid):
            setter(m-0.5*max_range, m+0.5*max_range)
        fig.tight_layout()
        fig.savefig(args.plot_tilts, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"Saved tilt plot to {args.plot_tilts}")

    # Visualization of grids
    plot_all(verts_base, sub_points, disk_pts, dome_pts, phi_grid, ring_ids, args)

    # Movie generation
    if args.movie_xyz and args.xyz:
        sys = AtomicSystem(fname=args.xyz)
        # adjust phi range if requested
        if args.phi_max is not None:
            phi_grid = np.linspace(0, np.deg2rad(args.phi_max), args.phi_samples, endpoint=False)
        out_path = args.movie_xyz
        save_movie_xyz(sys, tilts, phi_grid, args, out_path)


if __name__ == '__main__':
    main()
