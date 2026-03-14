#!/usr/bin/env python3

import os, sys, argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
if ROOT not in sys.path:
    sys.path.append(ROOT)

from pyBall.OCL.GridFF import GridFF_cl, GridShape


def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)


def make_system(case='single', sep=1.0, q=1.0):
    if case == 'single':
        apos = np.array([[0.0, 0.0, 0.0]], dtype=np.float32)
        qs = np.array([q], dtype=np.float32)
    elif case == 'pair':
        apos = np.array([[-0.5 * sep, 0.0, 0.0], [+0.5 * sep, 0.0, 0.0]], dtype=np.float32)
        qs = np.array([q, -q], dtype=np.float32)
    else:
        raise ValueError(f'Unknown case {case}')
    xyzq = np.zeros((len(apos), 4), dtype=np.float32)
    xyzq[:, :3] = apos
    xyzq[:, 3] = qs
    return xyzq, apos, qs


def get_grid(Lx=16.0, Ly=16.0, Lz=16.0, dg=(0.1, 0.1, 0.1)):
    lvec = np.array([[Lx, 0.0, 0.0], [0.0, Ly, 0.0], [0.0, 0.0, Lz]], dtype=np.float32)
    g0 = (-0.5 * Lx, -0.5 * Ly, -0.5 * Lz)
    grid = GridShape(dg=dg, lvec=lvec, g0=g0)
    return grid, lvec, g0


def make_solver(grid):
    clgff = GridFF_cl()
    clgff.set_grid(grid)
    return clgff


def project_density(grid, xyzq):
    clgff = make_solver(grid)
    return clgff.project_atoms_on_grid_quintic_pbc(xyzq, ng=tuple(clgff.gsh.ns), dg=tuple(clgff.gsh.dg[:3]), lvec=clgff.gsh.lvec, g0=clgff.gsh.g0, bReturn=True)


def apply_filter(grid, xyzq, sigma=0.0, bDivideByK2=False):
    clgff = make_solver(grid)
    clgff.project_atoms_on_grid_quintic_pbc(xyzq, ng=tuple(clgff.gsh.ns), dg=tuple(clgff.gsh.dg[:3]), lvec=clgff.gsh.lvec, g0=clgff.gsh.g0, bReturn=False)
    out = clgff.fft_filter(bReturn=True, sigma=sigma, bDivideByK2=bDivideByK2, bApplyEwaldScale=bDivideByK2)
    return out[:, :, :, 0].copy()


def linecut_x(arr):
    iz = arr.shape[0] // 2
    iy = arr.shape[1] // 2
    return arr[iz, iy, :].copy()


def fit_sigma_from_density(xs, rho, x0=0.0):
    w = np.maximum(rho, 0.0)
    s = w.sum()
    if s <= 0.0:
        return np.nan
    mu = np.sum(xs * w) / s
    var = np.sum(((xs - mu) ** 2) * w) / s
    return np.sqrt(max(var, 0.0))

def measure_fwhm(xs, rho):
    y = np.maximum(rho, 0.0)
    if len(y) < 3:
        return np.nan
    imax = int(np.argmax(y))
    ymax = y[imax]
    if ymax <= 0.0:
        return np.nan
    yhm = 0.5 * ymax
    xl = np.nan
    xr = np.nan
    for i in range(imax, 0, -1):
        y0 = y[i-1]
        y1 = y[i]
        if (y0 <= yhm <= y1) or (y1 <= yhm <= y0):
            t = 0.0 if abs(y1 - y0) < 1e-30 else (yhm - y0) / (y1 - y0)
            xl = xs[i-1] + t * (xs[i] - xs[i-1])
            break
    for i in range(imax, len(y) - 1):
        y0 = y[i]
        y1 = y[i+1]
        if (y0 >= yhm >= y1) or (y1 >= yhm >= y0):
            t = 0.0 if abs(y1 - y0) < 1e-30 else (yhm - y0) / (y1 - y0)
            xr = xs[i] + t * (xs[i+1] - xs[i])
            break
    if np.isnan(xl) or np.isnan(xr):
        return np.nan
    return xr - xl

def gaussian_fwhm_from_sigma(sigma):
    return 2.0 * np.sqrt(2.0 * np.log(2.0)) * sigma

def sigma_from_fwhm(fwhm):
    if not np.isfinite(fwhm):
        return np.nan
    return fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))


def print_stats(name, arr):
    arr = np.asarray(arr)
    print(f'{name:24s} min={arr.min(): .6e} max={arr.max(): .6e} range={(arr.max()-arr.min()): .6e}')


def plot_results(xs, rho0, rhoS, v0, vS, out_png, title, sigma_in, sigma_fit, sigma_fwhm):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4), constrained_layout=True)
    axes[0].plot(xs, rho0, label='density raw', lw=1.5)
    axes[0].plot(xs, rhoS, label=f'density smear sigma={sigma_in:.3f} A', lw=1.5)
    axes[0].set_xlabel('x [A]')
    axes[0].set_ylabel('density [arb.]')
    axes[0].set_title(f'Density cut, sigma_m2={sigma_fit:.3f} A, sigma_FWHM={sigma_fwhm:.3f} A')
    axes[0].grid(True)
    axes[0].legend()
    axes[1].plot(xs, v0, label='potential raw', lw=1.5)
    axes[1].plot(xs, vS, label=f'potential smear sigma={sigma_in:.3f} A', lw=1.5)
    axes[1].set_xlabel('x [A]')
    axes[1].set_ylabel('potential [eV]')
    axes[1].set_title('Potential cut')
    axes[1].grid(True)
    axes[1].legend()
    fig.suptitle(title)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def run_case(case='single', sigma=0.7, sep=1.0, q=1.0, dg=(0.1, 0.1, 0.1), Ls=(16.0, 16.0, 16.0), out_dir='data/gauss_smear'):
    ensure_dir(out_dir)
    grid, lvec, g0 = get_grid(Lx=Ls[0], Ly=Ls[1], Lz=Ls[2], dg=dg)
    xyzq, apos, qs = make_system(case=case, sep=sep, q=q)
    rho0 = project_density(grid, xyzq)
    rhoS = apply_filter(grid, xyzq, sigma=sigma, bDivideByK2=False)
    v0 = make_solver(grid).makeCoulombEwald(xyzq, sigma=0.0, bDensityOnly=False)
    vS = make_solver(grid).makeCoulombEwald(xyzq, sigma=sigma, bDensityOnly=False)
    xs = np.arange(rho0.shape[2], dtype=float) * dg[0] + g0[0]
    c_rho0 = linecut_x(rho0)
    c_rhoS = linecut_x(rhoS)
    c_v0 = linecut_x(v0)
    c_vS = linecut_x(vS)
    sigma_fit = fit_sigma_from_density(xs, c_rhoS)
    fwhm_fit = measure_fwhm(xs, c_rhoS)
    sigma_fwhm = sigma_from_fwhm(fwhm_fit)
    fwhm_ref = gaussian_fwhm_from_sigma(sigma)
    print(f'run_case() case={case} sigma_in={sigma:.6f} sigma_fit={sigma_fit:.6f} sigma_fwhm={sigma_fwhm:.6f} fwhm_fit={fwhm_fit:.6f} fwhm_ref={fwhm_ref:.6f} sep={sep:.6f}')
    print_stats('rho raw cut', c_rho0)
    print_stats('rho smear cut', c_rhoS)
    print_stats('V raw cut', c_v0)
    print_stats('V smear cut', c_vS)
    base = os.path.join(out_dir, f'{case}_sigma_{sigma:.3f}'.replace('.', 'p'))
    np.save(base + '_rho_raw.npy', rho0)
    np.save(base + '_rho_smear.npy', rhoS)
    np.save(base + '_V_raw.npy', v0)
    np.save(base + '_V_smear.npy', vS)
    np.savez(base + '_linecuts.npz', xs=xs, rho_raw=c_rho0, rho_smear=c_rhoS, V_raw=c_v0, V_smear=c_vS, sigma_in=sigma, sigma_fit=sigma_fit, sigma_fwhm=sigma_fwhm, fwhm_fit=fwhm_fit, fwhm_ref=fwhm_ref)
    plot_results(xs, c_rho0, c_rhoS, c_v0, c_vS, base + '_cuts.png', f'{case} q={q} sep={sep}', sigma, sigma_fit, sigma_fwhm)
    return {'sigma_fit': sigma_fit, 'sigma_fwhm': sigma_fwhm, 'fwhm_fit': fwhm_fit, 'fwhm_ref': fwhm_ref, 'base': base, 'rho0': rho0, 'rhoS': rhoS, 'v0': v0, 'vS': vS}


def main():
    ap = argparse.ArgumentParser(description='Validate Gaussian smearing in the GridFF pyOpenCL Fourier solver.')
    ap.add_argument('--case', default='single', choices=['single', 'pair'])
    ap.add_argument('--sigma', type=float, default=0.7, help='Gaussian sigma in Angstrom')
    ap.add_argument('--sep', type=float, default=1.0, help='Pair separation in Angstrom for case=pair')
    ap.add_argument('--q', type=float, default=1.0, help='Charge magnitude in e')
    ap.add_argument('--dgx', type=float, default=0.1)
    ap.add_argument('--dgy', type=float, default=0.1)
    ap.add_argument('--dgz', type=float, default=0.1)
    ap.add_argument('--Lx', type=float, default=16.0)
    ap.add_argument('--Ly', type=float, default=16.0)
    ap.add_argument('--Lz', type=float, default=16.0)
    ap.add_argument('--out_dir', default=os.path.join(os.path.dirname(__file__), 'data', 'gauss_smear'))
    args = ap.parse_args()
    run_case(case=args.case, sigma=args.sigma, sep=args.sep, q=args.q, dg=(args.dgx, args.dgy, args.dgz), Ls=(args.Lx, args.Ly, args.Lz), out_dir=args.out_dir)


if __name__ == '__main__':
    main()
