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
from pyBall.OCL.MMparams import read_atom_types, read_element_types


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

def sigma_tag(sigma):
    return f'sigma_{sigma:.3f}'.replace('.', 'p')

def parse_float_list(txt):
    return [float(s.strip()) for s in txt.split(',') if len(s.strip()) > 0]

def parse_probe_specs(txt):
    out = []
    for item in txt.split(','):
        item = item.strip()
        if len(item) == 0:
            continue
        name, q = item.split(':')
        out.append((name.strip(), float(q)))
    return out

def load_atom_type_table(atom_types_path, element_types_path):
    etypes = read_element_types(element_types_path)
    return read_atom_types(atom_types_path, etypes)

def get_probe_plqh(probe_name, q, atom_types, a_default):
    if probe_name not in atom_types:
        raise KeyError(f'Probe atom type {probe_name} not found in AtomTypes.dat')
    at = atom_types[probe_name]
    return ut.getPLQH(R0=at.RvdW, E0=at.EvdW, a=a_default, Q=q, H=0.0), at

def get_probe_plqh_custom(probe_name, q, atom_types, a, e0_scale=1.0, e0_override=None, r0_override=None):
    if probe_name not in atom_types:
        raise KeyError(f'Probe atom type {probe_name} not found in AtomTypes.dat')
    at = atom_types[probe_name]
    R0 = at.RvdW if r0_override is None else float(r0_override)
    E0 = at.EvdW if e0_override is None else float(e0_override)
    E0 *= float(e0_scale)
    plqh = ut.getPLQH(R0=R0, E0=E0, a=a, Q=q, H=0.0)
    info = {'name': probe_name, 'R0_src': at.RvdW, 'E0_src': at.EvdW, 'R0': R0, 'E0': E0, 'a': float(a), 'q': float(q), 'cP': float(plqh[0]), 'cL': float(plqh[1]), 'cQ': float(plqh[2]), 'cH': float(plqh[3])}
    return plqh, at, info

def find_top_atom_indices(atoms, species=('F', 'Ca')):
    out = {}
    enames = np.array(atoms.enames, dtype=object)
    zs = np.array(atoms.apos[:, 2], dtype=float)
    for name in species:
        sel = np.where(enames == name)[0]
        if len(sel) == 0:
            raise ValueError(f'No atoms of species {name} found in structure')
        out[name] = int(sel[np.argmax(zs[sel])])
    return out

def coord_to_index(x, dg, n, g0=0.0):
    return max(0, min(int(round((float(x) - float(g0)) / float(dg))), n - 1))

def extract_vertical_scan(gpu, dg, ix, iy, plqh):
    z = np.arange(gpu.shape[2], dtype=float) * dg[2]
    VPaul = gpu[ix, iy, :, 0].copy()
    VLond = gpu[ix, iy, :, 1].copy()
    VCoul = gpu[ix, iy, :, 2].copy()
    VNonElec = VPaul * plqh[0] + VLond * plqh[1]
    VTotal = VNonElec + VCoul * plqh[2]
    return {'z': z, 'VPaul': VPaul, 'VLond': VLond, 'VCoul': VCoul, 'VNonElec': VNonElec, 'VTotal': VTotal}

def pair_prefactors(surface_at, probe_info):
    ai = probe_info['a']
    Ri = float(surface_at.RvdW)
    Ei = float(surface_at.EvdW)
    Rj = float(probe_info['R0'])
    Ej = float(probe_info['E0'])
    return {
        'alpha': ai,
        'surface_R0': Ri,
        'surface_E0': Ei,
        'probe_R0': Rj,
        'probe_E0': Ej,
        'pair_cP': Ei * Ej * np.exp(2.0 * ai * (Ri + Rj)),
        'pair_cL': Ei * Ej * np.exp(ai * (Ri + Rj)),
        'r0_sum': Ri + Rj,
    }

def print_pair_diagnostics(target_name, probe_name, q_probe, surface_at, probe_info):
    p = pair_prefactors(surface_at, probe_info)
    print(f'PAIR {target_name:>2s}-{probe_name:<2s} q={q_probe:+.3f} alpha={p["alpha"]:.6f} Rsum={p["r0_sum"]:.6f} surface(R0={p["surface_R0"]:.6f},E0={p["surface_E0"]:.6e}) probe(R0={p["probe_R0"]:.6f},E0={p["probe_E0"]:.6e}) pair(cP={p["pair_cP"]:.6e}, cL={p["pair_cL"]:.6e})')
    for dr in (-0.5, 0.0, 0.5, 1.0):
        r = max(0.1, p['r0_sum'] + dr)
        e = np.exp(-p['alpha'] * (r - p['surface_R0']))
        Vp = surface_at.EvdW * e * e * probe_info['cP']
        Vl = surface_at.EvdW * e * (-2.0) * probe_info['cL']
        Vc = ut.COULOMB_CONST / r
        print(f'  r={r:7.3f}A d={dr:+5.2f}  VPauli={Vp: .6e}  VLond={Vl: .6e}  VCoul*q(+1,+1)={Vc: .6e}  VCoul*q(+1,-1)={-Vc: .6e}')

def find_zero_crossing(h, y):
    s = np.sign(y)
    for i in range(1, len(y)):
        if s[i] == 0:
            return float(h[i])
        if s[i-1] == 0:
            return float(h[i-1])
        if s[i] != s[i-1]:
            t = y[i-1] / (y[i-1] - y[i])
            return float(h[i-1] + t * (h[i] - h[i-1]))
    return np.nan

def summarize_scan_case(target_name, probe_name, sigma, scan):
    itot_min = int(np.argmin(scan['VTotal']))
    itot_max = int(np.argmax(scan['VTotal']))
    icoul_near = int(np.argmin(np.abs(scan['h'] - max(0.5, scan['h'][0]))))
    zcross = find_zero_crossing(scan['h'], scan['VTotal'])
    print(f'SCAN {target_name:>2s}-{probe_name:<2s} sigma={sigma:.3f} h[min]={scan["h"][itot_min]:.6f} Vtot_min={scan["VTotal"][itot_min]: .6e} h[max]={scan["h"][itot_max]:.6f} Vtot_max={scan["VTotal"][itot_max]: .6e} Vcoul@near={scan["VCoul"][icoul_near]*scan["probe_q"]: .6e} Vnonel@near={scan["VNonElec"][icoul_near]: .6e} zero_cross={zcross: .6f}')

def save_scan_npz(path, data):
    np.savez(path, **data)

def plot_scan_compare(out_png, scans, target_name, probe_name, q):
    fig, axes = plt.subplots(1, 4, figsize=(19, 4), constrained_layout=True)
    for sigma, data in scans:
        lbl = f'sigma={sigma:.2f} A'
        axes[0].plot(data['h'], data['VTotal'], lw=1.5, label=lbl)
        axes[1].plot(data['h'], data['VCoul'] * q, lw=1.5, label=lbl)
        axes[2].plot(data['h'], data['VPaul'] * data['plqh'][0], lw=1.5, label=lbl)
        axes[3].plot(data['h'], data['VNonElec'], lw=1.5, label=lbl)
    axes[0].set_title(f'{target_name} / {probe_name} q={q:+.2f} Total')
    axes[1].set_title(f'{target_name} / {probe_name} q={q:+.2f} Coulomb')
    axes[2].set_title(f'{target_name} / {probe_name} q={q:+.2f} Pauli')
    axes[3].set_title(f'{target_name} / {probe_name} q={q:+.2f} Non-electrostatic')
    for ax in axes:
        ax.set_xlabel('z - z_atom [A]')
        ax.set_ylabel('Energy [eV]')
        ax.grid(True)
        ax.legend()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


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
    ap.add_argument('--a', type=float, default=2.0, help='Probe a for combined linecut plot')
    ap.add_argument('--Q', type=float, default=0.4, help='Probe charge for combined linecut plot')
    ap.add_argument('--slice_z', type=int, default=None, help='z-slice index for 2D component plots')
    ap.add_argument('--z_above_top', type=float, default=-1.0, help='Height above topmost slab atom for xy cut [A], negative means auto-select')
    ap.add_argument('--sigma', type=float, default=0.0, help='Gaussian smearing sigma [A] applied in reciprocal-space Poisson solve')
    ap.add_argument('--sigma_scan', default='0.0,1.0', help='Comma-separated Gaussian sigma values [A] for comparative z-line scans')
    ap.add_argument('--probe_specs', default='H:0.2,O:-0.4', help='Comma-separated probe definitions name:q for z-line scans')
    ap.add_argument('--scan_zmin', type=float, default=0.0, help='Minimum height above target atom for z-line scan [A]')
    ap.add_argument('--scan_zmax', type=float, default=6.0, help='Maximum height above target atom for z-line scan [A]')
    ap.add_argument('--probe_a_scan', type=float, default=None, help='Override Morse alpha for scan probes; default uses --a')
    ap.add_argument('--probe_e0_scale', type=float, default=1.0, help='Scale scan-probe EvdW to strengthen/soften Pauli+London')
    ap.add_argument('--probe_h_e0', type=float, default=None, help='Explicit EvdW override for H probe in scans')
    ap.add_argument('--probe_o_e0', type=float, default=None, help='Explicit EvdW override for O probe in scans')
    ap.add_argument('--probe_h_r0', type=float, default=None, help='Explicit RvdW override for H probe in scans')
    ap.add_argument('--probe_o_r0', type=float, default=None, help='Explicit RvdW override for O probe in scans')
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
    sigma_scan = parse_float_list(args.sigma_scan)
    if args.sigma not in sigma_scan:
        sigma_scan = [args.sigma] + [s for s in sigma_scan if abs(s - args.sigma) > 1e-12]
    probe_specs = parse_probe_specs(args.probe_specs)
    atom_type_table = load_atom_type_table(atom_types, element_types)
    scan_alpha = args.a if args.probe_a_scan is None else args.probe_a_scan
    probe_overrides = {'H': {'E0': args.probe_h_e0, 'R0': args.probe_h_r0}, 'O': {'E0': args.probe_o_e0, 'R0': args.probe_o_r0}}
    target_inds = find_top_atom_indices(atoms)
    scans_dir = os.path.join(out_dir, 'z_scans')
    ensure_dir(scans_dir)
    gpu_by_sigma = {}
    gpu_npy_by_sigma = {}
    z_top = float(np.max(atoms.apos[:,2]))
    z_cut = np.nan
    z_cut_idx = None

    for sigma in sigma_scan:
        stag = sigma_tag(sigma)
        gff_ocl.test_gridFF_ocl(fname=src_xyz, Atom_Types_name=atom_types, Element_Types_name=element_types, job=args.job, save_name=args.save_name, use_CG=bool(args.use_CG), nmaxiter=args.nmaxiter, nPerStep=args.nPerStep, damp=args.damp, save_fig=True, fig_path=os.path.join(out_dir, f'convergence_gpu_{stag}.png'), dg=dg, sigma=sigma, alpha_morse=scan_alpha)
        gpu_npy = find_plq_file(out_dir)
        gpu = load_plq(gpu_npy)
        gpu_by_sigma[sigma] = gpu
        sigma_npy = os.path.join(out_dir, f'Bspline_PLQd_{stag}.npy')
        np.save(sigma_npy, gpu)
        gpu_npy_by_sigma[sigma] = sigma_npy
        export_plq_xsf(os.path.join(out_dir, f'gpu_{stag}'), gpu, atoms.lvec)
        plqh = ut.getPLQH(R0=args.R0, E0=args.E0, a=args.a, Q=args.Q, H=0.0)
        print_stats(f'VPaul global {stag}', gpu[..., 0])
        print_stats(f'VLond global {stag}', gpu[..., 1])
        print_stats(f'VCoul global {stag}', gpu[..., 2])
        print_stats(f'VTotal global {stag}', combine_plq(gpu, plqh))
        if z_cut_idx is None:
            if args.z_above_top >= 0.0:
                z_cut = z_top + args.z_above_top
                z_cut_idx = int(round(z_cut / dg[2]))
            else:
                z_cut_idx, z_cut = choose_reasonable_cut_z(gpu, dg, z_top, plqh)
        plot_plq_components(gpu, os.path.join(out_dir, f'plq_components_gpu_{stag}.png'), slice_idx=(args.slice_z if args.slice_z is not None else z_cut_idx))
        plot_plq_xy(gpu, os.path.join(out_dir, f'plq_xy_gpu_{stag}.png'), atoms.lvec, z_cut_idx, plqh)
        plot_plq_xz(gpu, os.path.join(out_dir, f'plq_xz_gpu_{stag}.png'), atoms.lvec, iy_idx=gpu.shape[1]//2, plqh=plqh)
        plot_plq_linecuts(gpu, os.path.join(out_dir, f'plq_linecuts_gpu_{stag}.png'), plqh, ix=gpu.shape[0]//2, iy=gpu.shape[1]//2, iz=z_cut_idx)

    report_lines = [
        '# CaF2 Gaussian-smearing z-scan report',
        '',
        f'- source_xyz: `{src_xyz}`',
        f'- output_dir: `{out_dir}`',
        f'- grid_spacing: `{dg}`',
        f'- sigma_scan: `{sigma_scan}`',
        f'- probe_specs: `{probe_specs}`',
        f'- scan_probe_alpha: `{scan_alpha}`',
        f'- scan_probe_e0_scale: `{args.probe_e0_scale}`',
        f'- z_window_above_atom: `{(args.scan_zmin, args.scan_zmax)}`',
        '',
    ]

    gpu0 = list(gpu_by_sigma.values())[0]
    for target_name, ia in target_inds.items():
        p = atoms.apos[ia]
        surface_at = atom_type_table[target_name]
        g0x = -0.5 * float(atoms.lvec[0][0])
        g0y = -0.5 * float(atoms.lvec[1][1])
        ix = coord_to_index(p[0], dg[0], gpu0.shape[0], g0=g0x)
        iy = coord_to_index(p[1], dg[1], gpu0.shape[1], g0=g0y)
        z_atom = float(p[2])
        for probe_name, q_probe in probe_specs:
            ov = probe_overrides.get(probe_name, {})
            probe_plqh, at, probe_info = get_probe_plqh_custom(probe_name, q_probe, atom_type_table, scan_alpha, e0_scale=args.probe_e0_scale, e0_override=ov.get('E0', None), r0_override=ov.get('R0', None))
            print_pair_diagnostics(target_name, probe_name, q_probe, surface_at, probe_info)
            scans = []
            for sigma in sigma_scan:
                scan = extract_vertical_scan(gpu_by_sigma[sigma], dg, ix, iy, probe_plqh)
                h = scan['z'] - z_atom
                mask = (h >= args.scan_zmin) & (h <= args.scan_zmax)
                scan_cut = {k: (v[mask].copy() if isinstance(v, np.ndarray) else v) for k, v in scan.items()}
                scan_cut['h'] = h[mask].copy()
                scan_cut['ix'] = ix
                scan_cut['iy'] = iy
                scan_cut['target_z'] = z_atom
                scan_cut['probe_name'] = probe_name
                scan_cut['probe_q'] = q_probe
                scan_cut['probe_R0'] = probe_info['R0']
                scan_cut['probe_E0'] = probe_info['E0']
                scan_cut['probe_a'] = scan_alpha
                scan_cut['plqh'] = probe_plqh.copy()
                tag = f'zscan_{target_name}_{probe_name}_q{q_probe:+.2f}_{sigma_tag(sigma)}'.replace('+', 'p').replace('-', 'm')
                save_scan_npz(os.path.join(scans_dir, tag + '.npz'), scan_cut)
                summarize_scan_case(target_name, probe_name, sigma, scan_cut)
                scans.append((sigma, scan_cut))
            out_png = os.path.join(scans_dir, f'zscan_compare_{target_name}_{probe_name}_q{q_probe:+.2f}.png'.replace('+', 'p').replace('-', 'm'))
            plot_scan_compare(out_png, scans, target_name, probe_name, q_probe)
            report_lines += [
                f'## {target_name} target / {probe_name} probe q={q_probe:+.2f}',
                f'- target_index: `{ia}`',
                f'- target_position: `{tuple(float(x) for x in p)}`',
                f'- grid_indices_xy: `{(ix, iy)}`',
                f'- probe_R0_E0_a: `{(probe_info["R0"], probe_info["E0"], scan_alpha)}`',
                f'- pair_prefactors: `{pair_prefactors(surface_at, probe_info)}`',
                f'- plot: `{out_png}`',
            ]
            for sigma, scan_cut in scans:
                report_lines.append(f'- sigma={sigma:.3f}: VTotal range=`{float(scan_cut["VTotal"].min()):.6f} .. {float(scan_cut["VTotal"].max()):.6f}` VCoul*q range=`{float((scan_cut["VCoul"]*q_probe).min()):.6f} .. {float((scan_cut["VCoul"]*q_probe).max()):.6f}`')
            report_lines.append('')

    report_path = os.path.join(out_dir, 'gaussian_smearing_report.md')
    with open(report_path, 'w') as f:
        f.write('\n'.join(report_lines) + '\n')

    print('Source xyz   :', src_xyz)
    print('Output dir   :', out_dir)
    print('Grid spacing :', dg)
    print('Sigma main   :', args.sigma)
    print('Sigma scan   :', sigma_scan)
    print('Top z / cut z:', z_top, z_cut, ' -> iz=', z_cut_idx)
    for sigma in sigma_scan:
        stag = sigma_tag(sigma)
        print(f'GPU npy {stag:>12s}:', gpu_npy_by_sigma[sigma])
        print(f'Plots  {stag:>12s}:', os.path.join(out_dir, f'plq_components_gpu_{stag}.png'), os.path.join(out_dir, f'plq_xy_gpu_{stag}.png'), os.path.join(out_dir, f'plq_xz_gpu_{stag}.png'), os.path.join(out_dir, f'plq_linecuts_gpu_{stag}.png'))
    print('Z scans      :', scans_dir)
    print('Report       :', report_path)


if __name__ == '__main__':
    main()
