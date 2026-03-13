#!/usr/bin/env python3

import os, sys, argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
if ROOT not in sys.path:
    sys.path.append(ROOT)

from pyBall.AtomicSystem import AtomicSystem
from pyBall import MMFF as mmff
from pyBall.tests import ocl_GridFF_new as gff_ocl
from pyBall.tests import utils as ut


def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)


def export_plq_xsf(prefix, plq, lvec):
    cell = np.array(lvec, dtype=np.float64)
    mmff.saveXSF(prefix + '_VPaul.xsf', np.array(plq[..., 0], dtype=np.float64), cell=cell)
    mmff.saveXSF(prefix + '_VLond.xsf', np.array(plq[..., 1], dtype=np.float64), cell=cell)
    mmff.saveXSF(prefix + '_VCoul.xsf', np.array(plq[..., 2], dtype=np.float64), cell=cell)


def load_plq(path):
    arr = np.load(path)
    if arr.ndim != 4 or arr.shape[-1] < 3:
        raise ValueError(f'Unexpected PLQ shape in {path}: {arr.shape}')
    return arr


def find_plq_file(out_dir):
    cands = [
        os.path.join(out_dir, 'Bspline_PLQd_ocl.npy'),
        os.path.join(out_dir, 'Bspline_PLQd.npy'),
    ]
    for p in cands:
        if os.path.exists(p):
            return p
    raise FileNotFoundError(f'GPU GridFF output missing, tried: {cands}')


def combine_plq(gpu, plqh):
    return gpu[..., 0] * plqh[0] + gpu[..., 1] * plqh[1] + gpu[..., 2] * plqh[2]


def print_stats(name, arr):
    arr = np.asarray(arr)
    print(f'{name:24s} min={arr.min(): .6e} max={arr.max(): .6e} range={(arr.max()-arr.min()): .6e}')


def choose_reasonable_cut_z(gpu, dg, z_top, plqh, zmin_above=0.8, zmax_above=3.0):
    total = combine_plq(gpu, plqh)
    nz = gpu.shape[2]
    iz0 = max(0, min(int(round((z_top + zmin_above) / dg[2])), nz - 1))
    iz1 = max(iz0, min(int(round((z_top + zmax_above) / dg[2])), nz - 1))
    best_iz = iz0
    best_score = -1.0
    for iz in range(iz0, iz1 + 1):
        g = total[:, :, iz]
        score = g.max() - g.min()
        if score > best_score:
            best_score = score
            best_iz = iz
    z_sel = best_iz * dg[2]
    print(f'choose_reasonable_cut_z() z_top={z_top:.6f} iz_range=({iz0},{iz1}) z_range=({iz0*dg[2]:.6f},{iz1*dg[2]:.6f}) best_iz={best_iz} best_z={z_sel:.6f} best_range={best_score:.6e}')
    return best_iz, z_sel


def plot_plq_components(gpu, out_png, slice_idx=None):
    nz = gpu.shape[2]
    if slice_idx is None:
        slice_idx = nz // 2
    slice_idx = max(0, min(int(slice_idx), nz - 1))
    names = ['VPaul', 'VLond', 'VCoul']
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
    for i, name in enumerate(names):
        g = gpu[:, :, slice_idx, i].T
        print_stats(f'{name} xy iz={slice_idx}', g)
        im = axes[i].imshow(g, origin='lower')
        axes[i].set_title(f'{name} GPU z={slice_idx}')
        fig.colorbar(im, ax=axes[i], fraction=0.046)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def plot_plq_xy(gpu, out_png, lvec, z_cut_idx, plqh):
    nx, ny, nz, _ = gpu.shape
    z_cut_idx = max(0, min(int(z_cut_idx), nz - 1))
    cell = np.array(lvec, dtype=float)
    dx = np.linalg.norm(cell[0])
    dy = np.linalg.norm(cell[1])
    names = ['VPaul', 'VLond', 'VCoul', 'VTotal']
    total = combine_plq(gpu, plqh)
    fig, axes = plt.subplots(1, 4, figsize=(16, 4), constrained_layout=True)
    for i, name in enumerate(names):
        gxy = gpu[:, :, z_cut_idx, i].T if name != 'VTotal' else total[:, :, z_cut_idx].T
        print_stats(f'{name} xy iz={z_cut_idx}', gxy)
        im = axes[i].imshow(gxy, origin='lower', extent=[0.0, nx*dx, 0.0, ny*dy], aspect='equal')
        axes[i].set_title(f'{name} xy @ z={z_cut_idx}')
        axes[i].set_xlabel('x [A]')
        axes[i].set_ylabel('y [A]')
        fig.colorbar(im, ax=axes[i], fraction=0.046)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def plot_plq_xz(gpu, out_png, lvec, iy_idx, plqh):
    nx, ny, nz, _ = gpu.shape
    iy_idx = max(0, min(int(iy_idx), ny - 1))
    cell = np.array(lvec, dtype=float)
    dx = np.linalg.norm(cell[0])
    dz = np.linalg.norm(cell[2])
    names = ['VPaul', 'VLond', 'VCoul', 'VTotal']
    total = combine_plq(gpu, plqh)
    fig, axes = plt.subplots(1, 4, figsize=(16, 4), constrained_layout=True)
    for i, name in enumerate(names):
        gxz = gpu[:, iy_idx, :, i].T if name != 'VTotal' else total[:, iy_idx, :].T
        print_stats(f'{name} xz iy={iy_idx}', gxz)
        im = axes[i].imshow(gxz, origin='lower', extent=[0.0, nx*dx, 0.0, nz*dz], aspect='equal')
        axes[i].set_title(f'{name} xz @ y={iy_idx}')
        axes[i].set_xlabel('x [A]')
        axes[i].set_ylabel('z [A]')
        fig.colorbar(im, ax=axes[i], fraction=0.046)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def plot_plq_surface_cuts(gpu, out_png, lvec, z_cut_idx, plqh, iy_idx=None):
    nx, ny, nz, _ = gpu.shape
    z_cut_idx = max(0, min(int(z_cut_idx), nz - 1))
    if iy_idx is None:
        iy_idx = ny // 2
    iy_idx = max(0, min(int(iy_idx), ny - 1))
    cell = np.array(lvec, dtype=float)
    dx = np.linalg.norm(cell[0])
    dy = np.linalg.norm(cell[1])
    dz = np.linalg.norm(cell[2])
    names = ['VPaul', 'VLond', 'VCoul']
    total = combine_plq(gpu, plqh)
    fig, axes = plt.subplots(2, 4, figsize=(17, 8), constrained_layout=True)
    for i, name in enumerate(names):
        gxy = gpu[:, :, z_cut_idx, i].T
        gxz = gpu[:, iy_idx, :, i].T
        print_stats(f'{name} xy iz={z_cut_idx}', gxy)
        print_stats(f'{name} xz iy={iy_idx}', gxz)
        im0 = axes[0, i].imshow(gxy, origin='lower', extent=[0.0, nx*dx, 0.0, ny*dy], aspect='auto')
        axes[0, i].set_title(f'{name} xy @ iz={z_cut_idx}')
        axes[0, i].set_xlabel('x [A]')
        axes[0, i].set_ylabel('y [A]')
        fig.colorbar(im0, ax=axes[0, i], fraction=0.046)
        im1 = axes[1, i].imshow(gxz, origin='lower', extent=[0.0, nx*dx, 0.0, nz*dz], aspect='auto')
        axes[1, i].set_title(f'{name} xz @ iy={iy_idx}')
        axes[1, i].set_xlabel('x [A]')
        axes[1, i].set_ylabel('z [A]')
        fig.colorbar(im1, ax=axes[1, i], fraction=0.046)
    gxy = total[:, :, z_cut_idx].T
    gxz = total[:, iy_idx, :].T
    print_stats(f'VTotal xy iz={z_cut_idx}', gxy)
    print_stats(f'VTotal xz iy={iy_idx}', gxz)
    im0 = axes[0, 3].imshow(gxy, origin='lower', extent=[0.0, nx*dx, 0.0, ny*dy], aspect='auto')
    axes[0, 3].set_title(f'VTotal xy @ iz={z_cut_idx}')
    axes[0, 3].set_xlabel('x [A]')
    axes[0, 3].set_ylabel('y [A]')
    fig.colorbar(im0, ax=axes[0, 3], fraction=0.046)
    im1 = axes[1, 3].imshow(gxz, origin='lower', extent=[0.0, nx*dx, 0.0, nz*dz], aspect='auto')
    axes[1, 3].set_title(f'VTotal xz @ iy={iy_idx}')
    axes[1, 3].set_xlabel('x [A]')
    axes[1, 3].set_ylabel('z [A]')
    fig.colorbar(im1, ax=axes[1, 3], fraction=0.046)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def plot_plq_linecuts(gpu, out_png, plqh, ix=None, iy=None, iz=None):
    nz = gpu.shape[2]
    nx = gpu.shape[0]
    ny = gpu.shape[1]
    if iz is None:
        iz = nz // 2
    if ix is None:
        ix = nx // 2
    if iy is None:
        iy = ny // 2
    ix = max(0, min(int(ix), nx - 1))
    iy = max(0, min(int(iy), ny - 1))
    iz = max(0, min(int(iz), nz - 1))
    egz = gpu[ix, iy, :, 0] * plqh[0] + gpu[ix, iy, :, 1] * plqh[1] + gpu[ix, iy, :, 2] * plqh[2]
    egx = gpu[:, iy, iz, 0] * plqh[0] + gpu[:, iy, iz, 1] * plqh[1] + gpu[:, iy, iz, 2] * plqh[2]
    egy = gpu[ix, :, iz, 0] * plqh[0] + gpu[ix, :, iz, 1] * plqh[1] + gpu[ix, :, iz, 2] * plqh[2]
    print_stats(f'VTotal line z ix={ix} iy={iy}', egz)
    print_stats(f'VTotal line x iy={iy} iz={iz}', egx)
    print_stats(f'VTotal line y ix={ix} iz={iz}', egy)
    fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)
    ax.plot(egz, '-', lw=0.8, label='GPU z')
    ax.plot(egx, '-', lw=0.8, label='GPU x')
    ax.plot(egy, '-', lw=0.8, label='GPU y')
    ax.grid(True)
    ax.legend()
    ax.set_xlabel('grid index')
    ax.set_ylabel('combined PLQ potential')
    ax.set_title(f'Combined GPU PLQ linecuts ix={ix} iy={iy} iz={iz}')
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def parse_args():
    ap = argparse.ArgumentParser(description='Generate GPU GridFF for rectangular charged CaF2 slab, save npy/xsf, and plot debug figures.')
    ap.add_argument('--src_xyz', default=os.path.join(ROOT, 'cpp', 'common_resources', 'Substrates', 'generated_rect', 'CaF2_6L_Ni3_rect_nx2_nz1_L2_top.xyz'), help='Source CaF2 xyz with lvs')
    ap.add_argument('--save_name', default='double3', help='Output format selector for GPU helper')
    ap.add_argument('--job', default='PLQ', choices=['PLQ', 'MorseFit', 'PLQ_lin', 'Ewald', 'Morse'], help='GPU GridFF job')
    ap.add_argument('--use_CG', type=int, default=1, help='Use conjugate gradient in GPU fit when supported')
    ap.add_argument('--nmaxiter', type=int, default=1000, help='GPU fit iteration limit')
    ap.add_argument('--nPerStep', type=int, default=25, help='GPU fit steps between convergence checks')
    ap.add_argument('--damp', type=float, default=0.15, help='GPU MD fit damping if MD fit is used')
    ap.add_argument('--dgx', type=float, default=23.175/225.0, help='Grid spacing along x [A]')
    ap.add_argument('--dgy', type=float, default=20.070/200.0, help='Grid spacing along y [A]')
    ap.add_argument('--dgz', type=float, default=48.472/382.0, help='Grid spacing along z [A]')
    ap.add_argument('--R0', type=float, default=1.443, help='Probe R0 for combined linecut plot')
    ap.add_argument('--E0', type=float, default=0.00190802, help='Probe E0 for combined linecut plot')
    ap.add_argument('--a', type=float, default=1.5, help='Probe a for combined linecut plot')
    ap.add_argument('--Q', type=float, default=0.4, help='Probe charge for combined linecut plot')
    ap.add_argument('--slice_z', type=int, default=None, help='z-slice index for 2D component plots')
    ap.add_argument('--z_above_top', type=float, default=-1.0, help='Height above topmost slab atom for xy cut [A], negative means auto-select')
    return ap.parse_args()


def main():
    args = parse_args()
    cwd = os.getcwd()
    here = os.path.abspath(os.path.dirname(__file__))
    if cwd != here:
        raise RuntimeError(f'Run this script from tests/tMMFF, current cwd={cwd}, expected={here}')

    src_xyz = os.path.abspath(args.src_xyz)
    if not os.path.exists(src_xyz):
        raise FileNotFoundError(f'src_xyz not found: {src_xyz}')

    name = os.path.splitext(os.path.basename(src_xyz))[0]
    out_dir = os.path.join(here, 'data', name)
    ensure_dir(out_dir)

    atoms = AtomicSystem(fname=src_xyz, bPreinit=False)
    if atoms.lvec is None:
        raise ValueError(f'Source xyz has no lvs: {src_xyz}')
    if atoms.qs is None:
        raise ValueError(f'Source xyz has no charges: {src_xyz}')
    if np.max(np.abs(atoms.qs)) <= 1e-8:
        raise ValueError(f'Source xyz charges are all near zero: {src_xyz}')

    dg = (args.dgx, args.dgy, args.dgz)
    atom_types = os.path.join(ROOT, 'cpp', 'common_resources', 'AtomTypes.dat')
    element_types = os.path.join(ROOT, 'cpp', 'common_resources', 'ElementTypes.dat')
    gff_ocl.test_gridFF_ocl(fname=src_xyz, Atom_Types_name=atom_types, Element_Types_name=element_types, job=args.job, save_name=args.save_name, use_CG=bool(args.use_CG), nmaxiter=args.nmaxiter, nPerStep=args.nPerStep, damp=args.damp, save_fig=True, fig_path=os.path.join(out_dir, 'convergence_gpu.png'), dg=dg)

    gpu_npy = find_plq_file(out_dir)

    gpu = load_plq(gpu_npy)
    export_plq_xsf(os.path.join(out_dir, 'gpu'), gpu, atoms.lvec)

    plqh = ut.getPLQH(R0=args.R0, E0=args.E0, a=args.a, Q=args.Q, H=0.0)
    print_stats('VPaul global', gpu[..., 0])
    print_stats('VLond global', gpu[..., 1])
    print_stats('VCoul global', gpu[..., 2])
    print_stats('VTotal global', combine_plq(gpu, plqh))
    z_top = float(np.max(atoms.apos[:,2]))
    if args.z_above_top >= 0.0:
        z_cut = z_top + args.z_above_top
        z_cut_idx = int(round(z_cut / dg[2]))
    else:
        z_cut_idx, z_cut = choose_reasonable_cut_z(gpu, dg, z_top, plqh)
    plot_plq_components(gpu, os.path.join(out_dir, 'plq_components_gpu.png'), slice_idx=(args.slice_z if args.slice_z is not None else z_cut_idx))
    plot_plq_xy(gpu, os.path.join(out_dir, 'plq_xy_gpu.png'), atoms.lvec, z_cut_idx, plqh)
    plot_plq_xz(gpu, os.path.join(out_dir, 'plq_xz_gpu.png'), atoms.lvec, iy_idx=gpu.shape[1]//2, plqh=plqh)
    plot_plq_linecuts(gpu, os.path.join(out_dir, 'plq_linecuts_gpu.png'), plqh, ix=gpu.shape[0]//2, iy=gpu.shape[1]//2, iz=z_cut_idx)

    print('Source xyz   :', src_xyz)
    print('Output dir   :', out_dir)
    print('Grid spacing :', dg)
    print('Top z / cut z:', z_top, z_cut, ' -> iz=', z_cut_idx)
    print('GPU npy      :', gpu_npy)
    print('GPU XSF pref :', os.path.join(out_dir, 'gpu_*'))
    print('Plots        :', os.path.join(out_dir, 'plq_components_gpu.png'), os.path.join(out_dir, 'plq_xy_gpu.png'), os.path.join(out_dir, 'plq_xz_gpu.png'), os.path.join(out_dir, 'plq_linecuts_gpu.png'))


if __name__ == '__main__':
    main()
