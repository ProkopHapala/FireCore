#!/usr/bin/env python3
"""
test_interaction_scan.py - Test molecule-substrate interaction scanning
Test cases: H2O on NaCl and PTCDA on NaCl

Uses MMparams/AtomTypes.dat for proper MMFF REQ parameter initialization.

Demonstrates:
1. Z-approach curve (1D) with H2O
2. Lateral (x,y) energy map (2D)
3. Rotation scan (1D)
4. SLERP path scan
5. Constrained relaxation scan
6. PTCDA on NaCl (larger molecule with proper atom typing)
"""

import sys, os, argparse
import numpy as np
import time

# Ensure pyBall is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from pyBall.OCL.InteractionEnergy import InteractionScanner, load_xyz_with_REQs
from pyBall.OCL.MolecularDynamics import MolecularDynamics
from pyBall.OCL import ScanUtils

# ======== Paths ========
XYZ_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'cpp', 'common_resources', 'xyz')
OUT_DIR = os.path.join(os.path.dirname(__file__), 'output_interaction_scan')

MOL_PRESETS = {
    'h2o': ('H2O_O.xyz', {}, 1.0),
    'ptcda': ('PTCDA.xyz', {'C': 'C_R', 'O': 'O_2'}, 3.5),
}

SUB_PRESETS = {
    'nacl_1x1': 'NaCl_1x1_L3.xyz',
    'nacl_8x8': 'NaCl_8x8_L3.xyz',
}


def resolve_xyz_path(path_or_name):
    if path_or_name is None:
        return None
    if os.path.isabs(path_or_name):
        return path_or_name
    if os.path.exists(path_or_name):
        return os.path.abspath(path_or_name)
    return os.path.join(XYZ_DIR, path_or_name)


def scenario_tag(mol_file, sub_file):
    return f"{os.path.splitext(os.path.basename(mol_file))[0]}__{os.path.splitext(os.path.basename(sub_file))[0]}"


def default_xy_range(sub_file, mol_file):
    s = os.path.basename(sub_file).lower()
    m = os.path.basename(mol_file).lower()
    if '8x8' in s:
        if 'ptcda' in m:
            return (-10.0, 30.0)
        return (-5.0, 35.0)
    if 'ptcda' in m:
        return (-4.0, 4.0)
    return (0.0, 4.0)


def setup_reference_scanner(mol_file, sub_file, mol_type_map=None, nPBC=(4, 4, 0), enable_macro=True):
    scanner = InteractionScanner(nloc=32)
    scanner.load_molecule_xyz(mol_file, type_map=mol_type_map)
    scanner.load_substrate_xyz(sub_file)
    scanner.enable_LJ = False
    scanner.enable_Coulomb = True
    scanner.enable_HBond = False
    scanner.enable_Morse = True
    scanner.enable_macro = bool(enable_macro)
    scanner.nPBC[:] = tuple(int(v) for v in nPBC)
    if enable_macro:
        scanner._update_macro_from_substrate()
    return scanner


def _primitive_steps(scanner):
    if scanner.lvec is None:
        raise ValueError('Scanner substrate lattice vectors are not available')
    a = np.array(scanner.lvec[0], dtype=float)
    b = np.array(scanner.lvec[1], dtype=float)
    return a, b


def run_reference_case(mol_file, sub_file, out_dir=OUT_DIR, mol_type_map=None, nPBC=(4, 4, 0), z0=None, nx=121, ny=121, ns_line=9, xy_range=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    os.makedirs(out_dir, exist_ok=True)
    scanner = setup_reference_scanner(mol_file, sub_file, mol_type_map=mol_type_map, nPBC=nPBC, enable_macro=True)
    tag = scenario_tag(mol_file, sub_file)
    if z0 is None:
        z0 = 3.5 if 'ptcda' in os.path.basename(mol_file).lower() else 1.0
    if xy_range is None:
        xy_range = default_xy_range(sub_file, mol_file)
    p0 = np.array([0.0, 0.0, float(z0)], dtype=float)
    step_x, step_y = _primitive_steps(scanner)
    tx = [p0 + i * step_x for i in range(ns_line)]
    ty = [p0 + i * step_y for i in range(ns_line)]
    transforms_x = ScanUtils.pack_transforms([np.eye(3)] * len(tx), tx)
    transforms_y = ScanUtils.pack_transforms([np.eye(3)] * len(ty), ty)
    Ex = np.array(scanner.evaluate(transforms_x)['total'], dtype=float)
    Ey = np.array(scanner.evaluate(transforms_y)['total'], dtype=float)
    dEx = Ex - Ex[0]
    dEy = Ey - Ey[0]
    labs = np.arange(ns_line)
    fig, ax = plt.subplots(1, 1, figsize=(8, 4.5))
    ax.plot(labs, dEx, 'o-', lw=1.5, label=f'step a = ({step_x[0]:.2f},{step_x[1]:.2f}) Å')
    ax.plot(labs, dEy, 's-', lw=1.5, label=f'step b = ({step_y[0]:.2f},{step_y[1]:.2f}) Å')
    ax.axhline(0.0, c='k', lw=0.8)
    ax.set_xlabel('Equivalent-site index')
    ax.set_ylabel('ΔE [eV]')
    ax.set_title(f'{os.path.basename(mol_file)} on {os.path.basename(sub_file)} equivalent-site error, z={z0:.2f} Å')
    ax.grid(True, alpha=0.3)
    ax.legend(frameon=False, fontsize=8)
    line_png = os.path.join(out_dir, f'{tag}_equiv_line_macro.png')
    fig.tight_layout()
    fig.savefig(line_png, dpi=180)
    plt.close(fig)

    t0 = time.perf_counter()
    res_xy = scanner.scan_lateral(z=float(z0), x_range=xy_range, y_range=xy_range, nx=int(nx), ny=int(ny))
    wall_s = time.perf_counter() - t0
    E2d = np.array(res_xy['total'], dtype=float).reshape(res_xy['scan_info']['shape'])
    xs = res_xy['scan_info']['x']
    ys = res_xy['scan_info']['y']
    vabs = np.nanmax(np.abs(E2d - np.nanmean(E2d)))
    if not np.isfinite(vabs) or (vabs <= 0.0):
        raise ValueError(f'Invalid XY energy span {vabs} for {tag}')
    fig, ax = plt.subplots(1, 1, figsize=(8, 7))
    im = ax.imshow((E2d - np.mean(E2d)).T, extent=[xs[0], xs[-1], ys[0], ys[-1]], origin='lower', aspect='equal', cmap='bwr', vmin=-vabs, vmax=vabs)
    ax.set_xlabel('x [Å]')
    ax.set_ylabel('y [Å]')
    ax.set_title(f'{os.path.basename(mol_file)} on {os.path.basename(sub_file)} XY scan at z={z0:.2f} Å')
    plt.colorbar(im, ax=ax, label='E - <E> [eV]')
    xy_png = os.path.join(out_dir, f'{tag}_XYscan_macro.png')
    fig.tight_layout()
    fig.savefig(xy_png, dpi=180)
    plt.close(fig)
    print(f'[reference] {tag} line max|ΔE| x={np.max(np.abs(dEx)):.6e} y={np.max(np.abs(dEy)):.6e} 2D wall={wall_s:.3f}s')
    print(f'[reference] saved {line_png}')
    print(f'[reference] saved {xy_png}')
    return {'tag': tag, 'line_png': line_png, 'xy_png': xy_png, 'dEx': dEx, 'dEy': dEy, 'wall_s': wall_s, 'xy_range': xy_range, 'z0': z0}


def run_fast_case(mol_file, sub_file, out_dir=OUT_DIR, mol_type_map=None, nPBC=(4, 4, 0), z0=None, nx=121, ny=121, ns_line=9, xy_range=None, nSystems=8192, chunk_size=8192):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    os.makedirs(out_dir, exist_ok=True)
    scanner = setup_reference_scanner(mol_file, sub_file, mol_type_map=mol_type_map, nPBC=nPBC, enable_macro=True)
    md, mol_pos, mol_REQs, mol_names = setup_fast_gpu_scanner(mol_file, sub_file, nSystems=nSystems, mol_type_map=mol_type_map)
    tag = scenario_tag(mol_file, sub_file)
    if z0 is None:
        z0 = 3.5 if 'ptcda' in os.path.basename(mol_file).lower() else 1.0
    if xy_range is None:
        xy_range = default_xy_range(sub_file, mol_file)
    p0 = np.array([0.0, 0.0, float(z0)], dtype=float)
    step_x, step_y = _primitive_steps(scanner)
    tx = np.array([p0 + i * step_x for i in range(ns_line)], dtype=np.float32)
    ty = np.array([p0 + i * step_y for i in range(ns_line)], dtype=np.float32)
    transforms_x = ScanUtils.pack_transforms([np.eye(3)] * len(tx), tx)
    transforms_y = ScanUtils.pack_transforms([np.eye(3)] * len(ty), ty)
    cpu_x = np.array(scanner.evaluate(transforms_x)['total'], dtype=float)
    cpu_y = np.array(scanner.evaluate(transforms_y)['total'], dtype=float)
    gpu_x = np.array(eval_fast_gpu(md, transforms_x, chunk_size=chunk_size)['total'], dtype=float)
    gpu_y = np.array(eval_fast_gpu(md, transforms_y, chunk_size=chunk_size)['total'], dtype=float)
    dEx_cpu = cpu_x - cpu_x[0]
    dEy_cpu = cpu_y - cpu_y[0]
    dEx_gpu = gpu_x - gpu_x[0]
    dEy_gpu = gpu_y - gpu_y[0]
    labs = np.arange(ns_line)
    fig, ax = plt.subplots(1, 1, figsize=(8, 4.5))
    ax.plot(labs, dEx_cpu, 'o-', lw=1.5, label='CPU ref a')
    ax.plot(labs, dEy_cpu, 's-', lw=1.5, label='CPU ref b')
    ax.plot(labs, dEx_gpu, 'o--', lw=1.2, label='GPU fast a')
    ax.plot(labs, dEy_gpu, 's--', lw=1.2, label='GPU fast b')
    ax.axhline(0.0, c='k', lw=0.8)
    ax.set_xlabel('Equivalent-site index')
    ax.set_ylabel('ΔE [eV]')
    ax.set_title(f'Fast parity {os.path.basename(mol_file)} on {os.path.basename(sub_file)}')
    ax.grid(True, alpha=0.3)
    ax.legend(frameon=False, fontsize=8)
    line_png = os.path.join(out_dir, f'{tag}_equiv_line_fast_gpu.png')
    fig.tight_layout()
    fig.savefig(line_png, dpi=180)
    plt.close(fig)

    transforms_xy, info_xy = ScanUtils.scan_lateral_2d(float(z0), xy_range, xy_range, R=np.eye(3), nx=int(nx), ny=int(ny))
    cpu_xy = np.array(scanner.evaluate(transforms_xy)['total'], dtype=float)
    gpu_xy = eval_fast_gpu(md, transforms_xy, chunk_size=chunk_size)
    Egpu = np.array(gpu_xy['total'], dtype=float).reshape(info_xy['shape'])
    Ecpu = cpu_xy.reshape(info_xy['shape'])
    dE2d = Egpu - Ecpu
    xs = info_xy['x']
    ys = info_xy['y']
    vmax_map = np.nanmax(np.abs(Egpu - np.nanmean(Egpu)))
    vmax_err = np.nanmax(np.abs(dE2d))
    fig, axs = plt.subplots(1, 2, figsize=(14, 6))
    im0 = axs[0].imshow((Egpu - np.mean(Egpu)).T, extent=[xs[0], xs[-1], ys[0], ys[-1]], origin='lower', aspect='equal', cmap='bwr', vmin=-vmax_map, vmax=vmax_map)
    axs[0].set_xlabel('x [Å]')
    axs[0].set_ylabel('y [Å]')
    axs[0].set_title('Fast GPU XY scan')
    plt.colorbar(im0, ax=axs[0], label='E - <E> [eV]')
    im1 = axs[1].imshow(dE2d.T, extent=[xs[0], xs[-1], ys[0], ys[-1]], origin='lower', aspect='equal', cmap='bwr', vmin=-vmax_err, vmax=vmax_err)
    axs[1].set_xlabel('x [Å]')
    axs[1].set_ylabel('y [Å]')
    axs[1].set_title('Fast GPU - CPU reference')
    plt.colorbar(im1, ax=axs[1], label='ΔE [eV]')
    xy_png = os.path.join(out_dir, f'{tag}_XYscan_fast_gpu.png')
    fig.tight_layout()
    fig.savefig(xy_png, dpi=180)
    plt.close(fig)
    print(f'[fast] {tag} line max|ΔE| x={np.max(np.abs(gpu_x-cpu_x)):.6e} y={np.max(np.abs(gpu_y-cpu_y)):.6e} 2D max|ΔE|={np.max(np.abs(dE2d)):.6e} wall={gpu_xy["wall_s"]:.3f}s prep={gpu_xy.get("t_prep_s",0.0):.3f}s kernel={gpu_xy.get("t_kernel_s",0.0):.3f}s download={gpu_xy.get("t_download_s",0.0):.3f}s')
    print(f'[fast] saved {line_png}')
    print(f'[fast] saved {xy_png}')
    return {'tag': tag, 'line_png': line_png, 'xy_png': xy_png, 'wall_s': gpu_xy['wall_s'], 't_prep_s': gpu_xy.get('t_prep_s', 0.0), 't_kernel_s': gpu_xy.get('t_kernel_s', 0.0), 't_download_s': gpu_xy.get('t_download_s', 0.0), 'dE2d_max': float(np.max(np.abs(dE2d)))}


def default_scenarios():
    scs = []
    for mol_key, (mol_name, mol_type_map, z0) in MOL_PRESETS.items():
        for sub_key, sub_name in SUB_PRESETS.items():
            scs.append({
                'mol_file': resolve_xyz_path(mol_name),
                'sub_file': resolve_xyz_path(sub_name),
                'mol_type_map': dict(mol_type_map),
                'z0': z0,
            })
    return scs


def run_batch_cases(args):
    results = []
    cases = default_scenarios() if args.batch_presets else [{
        'mol_file': resolve_xyz_path(args.mol or MOL_PRESETS['h2o'][0]),
        'sub_file': resolve_xyz_path(args.sub or SUB_PRESETS['nacl_8x8']),
        'mol_type_map': dict(MOL_PRESETS.get(args.mol_preset, ('', {}, 0.0))[1]) if args.mol_preset else {},
        'z0': args.z0,
    }]
    for case in cases:
        if args.mode in ('reference', 'all'):
            results.append(run_reference_case(case['mol_file'], case['sub_file'], out_dir=args.out_dir, mol_type_map=case['mol_type_map'], nPBC=tuple(args.npbc), z0=case['z0'], nx=args.nx, ny=args.ny, ns_line=args.ns_line))
        if args.mode in ('fast', 'all'):
            results.append(run_fast_case(case['mol_file'], case['sub_file'], out_dir=args.out_dir, mol_type_map=case['mol_type_map'], nPBC=tuple(args.npbc), z0=case['z0'], nx=args.nx, ny=args.ny, ns_line=args.ns_line, nSystems=args.fast_nsystems, chunk_size=args.fast_chunk))
    return results


def setup_fast_gpu_scanner(mol_file, sub_file, nSystems=8192, mol_type_map=None):
    mol_pos, mol_REQs, mol_names, mol_Zs, mol_lvec = load_xyz_with_REQs(mol_file, type_map=mol_type_map)
    md = MolecularDynamics(nloc=32, debug_build_options='-DDBG_UFF=0')
    md.init_rigid_molecule_batch(mol_pos, mol_REQs, nSystems=nSystems)
    md.set_surface(sub_file, nPBC=(4, 4, 0), pos0=(0.0, 0.0, 0.0), alpha_morse=1.8, r_damp=0.0, bMacro=True)
    return md, mol_pos, mol_REQs, mol_names


def eval_fast_gpu(md, transforms, chunk_size=None):
    t0 = time.perf_counter()
    out = md.eval_rigid_getSurfMorse(transforms, chunk_size=chunk_size)
    out['wall_s'] = time.perf_counter() - t0
    return out


def run_macro_fast_gpu_scan():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    os.makedirs(OUT_DIR, exist_ok=True)
    mol_file = os.path.join(XYZ_DIR, 'H2O_O.xyz')
    sub_file = os.path.join(XYZ_DIR, 'NaCl_8x8_L3.xyz')

    scanner = InteractionScanner(nloc=32)
    scanner.load_molecule_xyz(mol_file)
    scanner.load_substrate_xyz(sub_file)
    scanner.enable_LJ = True
    scanner.enable_Coulomb = True
    scanner.enable_HBond = False
    scanner.enable_Morse = True
    scanner.enable_macro = True
    scanner.nPBC[:] = (4, 4, 0)
    scanner._update_macro_from_substrate()

    md, mol_pos, mol_REQs, mol_names = setup_fast_gpu_scanner(mol_file, sub_file, nSystems=8192, mol_type_map=None)

    p0 = np.array([1.0, 1.0, 1.0], dtype=float)
    step_x = np.array([4.0, 0.0, 0.0], dtype=float)
    step_y = np.array([0.0, 4.0, 0.0], dtype=float)
    ns = 9
    tx = np.array([p0 + i * step_x for i in range(ns)], dtype=np.float32)
    ty = np.array([p0 + i * step_y for i in range(ns)], dtype=np.float32)
    transforms_x = ScanUtils.pack_transforms([np.eye(3)] * len(tx), tx)
    transforms_y = ScanUtils.pack_transforms([np.eye(3)] * len(ty), ty)

    cpu_x = np.array(scanner.evaluate(transforms_x)['total'], dtype=float)
    cpu_y = np.array(scanner.evaluate(transforms_y)['total'], dtype=float)
    gpu_x = np.array(eval_fast_gpu(md, transforms_x)['total'], dtype=float)
    gpu_y = np.array(eval_fast_gpu(md, transforms_y)['total'], dtype=float)
    dEx_cpu = cpu_x - cpu_x[0]
    dEy_cpu = cpu_y - cpu_y[0]
    dEx_gpu = gpu_x - gpu_x[0]
    dEy_gpu = gpu_y - gpu_y[0]
    err_x = gpu_x - cpu_x
    err_y = gpu_y - cpu_y
    labs = np.arange(ns)
    print('Fast GPU parity along x [eV]:')
    for i, ec, eg, de in zip(labs, cpu_x, gpu_x, err_x):
        print(f'  ix={i:2d}  E_cpu={ec:+.6e}  E_gpu={eg:+.6e}  dE={de:+.6e}')
    print('Fast GPU parity along y [eV]:')
    for i, ec, eg, de in zip(labs, cpu_y, gpu_y, err_y):
        print(f'  iy={i:2d}  E_cpu={ec:+.6e}  E_gpu={eg:+.6e}  dE={de:+.6e}')

    fig, ax = plt.subplots(1, 1, figsize=(8, 4.5))
    ax.plot(labs, dEx_cpu, 'o-', lw=1.5, label='CPU ref step (4,0) Å')
    ax.plot(labs, dEy_cpu, 's-', lw=1.5, label='CPU ref step (0,4) Å')
    ax.plot(labs, dEx_gpu, 'o--', lw=1.2, label='GPU fast step (4,0) Å')
    ax.plot(labs, dEy_gpu, 's--', lw=1.2, label='GPU fast step (0,4) Å')
    ax.axhline(0.0, c='k', lw=0.8)
    ax.set_xlabel('Equivalent-site index')
    ax.set_ylabel('ΔE [eV]')
    ax.set_title('Equivalent-site parity: CPU reference vs fast GPU getSurfMorse')
    ax.grid(True, alpha=0.3)
    ax.legend(frameon=False, fontsize=8)
    line_png = os.path.join(OUT_DIR, 'H2O_NaCl8x8_equiv_line_gpu_macro.png')
    fig.tight_layout()
    fig.savefig(line_png, dpi=180)
    plt.close(fig)

    transforms_xy, info_xy = ScanUtils.scan_lateral_2d(1.0, (-5.0, 35.0), (-5.0, 35.0), R=np.eye(3), nx=161, ny=161)
    t0 = time.perf_counter()
    cpu_xy = np.array(scanner.evaluate(transforms_xy)['total'], dtype=float)
    cpu_xy_s = time.perf_counter() - t0
    gpu_xy = eval_fast_gpu(md, transforms_xy, chunk_size=8192)
    gpu_xy_e = np.array(gpu_xy['total'], dtype=float)
    Ecpu = cpu_xy.reshape(info_xy['shape'])
    Egpu = gpu_xy_e.reshape(info_xy['shape'])
    dE2d = Egpu - Ecpu
    xs = info_xy['x']
    ys = info_xy['y']
    vmax_map = np.nanmax(np.abs(Egpu - np.nanmean(Egpu)))
    vmax_err = np.nanmax(np.abs(dE2d))
    if not np.isfinite(vmax_map) or vmax_map <= 0.0:
        raise ValueError(f'Invalid fast GPU XY energy span {vmax_map}')
    if not np.isfinite(vmax_err):
        raise ValueError(f'Invalid fast GPU XY parity error span {vmax_err}')

    fig, axs = plt.subplots(1, 2, figsize=(14, 6))
    im0 = axs[0].imshow((Egpu - np.mean(Egpu)).T, extent=[xs[0], xs[-1], ys[0], ys[-1]], origin='lower', aspect='equal', cmap='bwr', vmin=-vmax_map, vmax=vmax_map)
    axs[0].set_xlabel('x [Å]')
    axs[0].set_ylabel('y [Å]')
    axs[0].set_title('Fast GPU XY scan at z=1.0 Å (mean-shifted)')
    plt.colorbar(im0, ax=axs[0], label='E - <E> [eV]')
    im1 = axs[1].imshow(dE2d.T, extent=[xs[0], xs[-1], ys[0], ys[-1]], origin='lower', aspect='equal', cmap='bwr', vmin=-vmax_err, vmax=vmax_err)
    axs[1].set_xlabel('x [Å]')
    axs[1].set_ylabel('y [Å]')
    axs[1].set_title('Fast GPU - CPU reference parity error')
    plt.colorbar(im1, ax=axs[1], label='ΔE [eV]')
    xy_png = os.path.join(OUT_DIR, 'H2O_NaCl8x8_XYscan_-5_35_npbc4_macro_gpu.png')
    fig.tight_layout()
    fig.savefig(xy_png, dpi=180)
    plt.close(fig)

    transforms_xy_speed, info_speed = ScanUtils.scan_lateral_2d(1.0, (-5.0, 35.0), (-5.0, 35.0), R=np.eye(3), nx=200, ny=200)
    gpu_speed = eval_fast_gpu(md, transforms_xy_speed, chunk_size=8192)
    nxy_speed = len(transforms_xy_speed)
    print(f'Fast GPU timing 200x200 ({nxy_speed} points): {gpu_speed["wall_s"]:.6f} s = {1e3*gpu_speed["wall_s"]/nxy_speed:.6f} ms/point')
    print(f'CPU reference timing 161x161 ({len(transforms_xy)} points): {cpu_xy_s:.6f} s')
    print(f'Fast GPU timing 161x161 ({len(transforms_xy)} points): {gpu_xy["wall_s"]:.6f} s')
    print(f'Fast GPU parity max|ΔE| line-x = {np.max(np.abs(err_x)):.6e} eV')
    print(f'Fast GPU parity max|ΔE| line-y = {np.max(np.abs(err_y)):.6e} eV')
    print(f'Fast GPU parity max|ΔE| 2D = {np.max(np.abs(dE2d)):.6e} eV')

    return {
        'line_png': line_png,
        'xy_png': xy_png,
        'cpu_line_x': cpu_x,
        'cpu_line_y': cpu_y,
        'gpu_line_x': gpu_x,
        'gpu_line_y': gpu_y,
        'cpu_xy': Ecpu,
        'gpu_xy': Egpu,
        'dE2d': dE2d,
        'gpu_xy_wall_s': gpu_xy['wall_s'],
        'cpu_xy_wall_s': cpu_xy_s,
        'gpu_speed_wall_s': gpu_speed['wall_s'],
    }


def run_macro_reference_scan():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    os.makedirs(OUT_DIR, exist_ok=True)
    mol_file = os.path.join(XYZ_DIR, 'H2O_O.xyz')
    sub_file = os.path.join(XYZ_DIR, 'NaCl_8x8_L3.xyz')
    scanner = InteractionScanner(nloc=32)
    scanner.load_molecule_xyz(mol_file)
    scanner.load_substrate_xyz(sub_file)
    scanner.enable_LJ = True
    scanner.enable_Coulomb = True
    scanner.enable_HBond = False
    scanner.enable_Morse = True
    scanner.enable_macro = True
    scanner.nPBC[:] = (4, 4, 0)
    scanner._update_macro_from_substrate()
    p0 = np.array([1.0, 1.0, 1.0], dtype=float)
    step_x = np.array([4.0, 0.0, 0.0], dtype=float)
    step_y = np.array([0.0, 4.0, 0.0], dtype=float)
    ns = 9
    tx = [p0 + i * step_x for i in range(ns)]
    ty = [p0 + i * step_y for i in range(ns)]
    transforms_x = ScanUtils.pack_transforms([np.eye(3)] * len(tx), tx)
    transforms_y = ScanUtils.pack_transforms([np.eye(3)] * len(ty), ty)
    Ex = np.array(scanner.evaluate(transforms_x)['total'], dtype=float)
    Ey = np.array(scanner.evaluate(transforms_y)['total'], dtype=float)
    dEx = Ex - Ex[0]
    dEy = Ey - Ey[0]
    labs = np.arange(ns)
    print('Equivalent-site line errors along x [eV]:')
    for i, de in zip(labs, dEx):
        print(f'  ix={i:2d}  dE={de:+.6e}')
    print('Equivalent-site line errors along y [eV]:')
    for i, de in zip(labs, dEy):
        print(f'  iy={i:2d}  dE={de:+.6e}')
    fig, ax = plt.subplots(1, 1, figsize=(8, 4.5))
    ax.plot(labs, dEx, 'o-', lw=1.5, label='step (4,0) Å')
    ax.plot(labs, dEy, 's-', lw=1.5, label='step (0,4) Å')
    ax.axhline(0.0, c='k', lw=0.8)
    ax.set_xlabel('Equivalent-site index along a+b')
    ax.set_ylabel('ΔE [eV]')
    ax.set_title('Equivalent-site energy error, H2O_O on NaCl8x8, z=1.0 Å')
    ax.grid(True, alpha=0.3)
    ax.legend(frameon=False)
    line_png = os.path.join(OUT_DIR, 'H2O_NaCl8x8_equiv_line_macro.png')
    fig.tight_layout()
    fig.savefig(line_png, dpi=180)
    plt.close(fig)

    res_xy = scanner.scan_lateral(z=1.0, x_range=(-5.0, 35.0), y_range=(-5.0, 35.0), nx=161, ny=161)
    E2d = np.array(res_xy['total'], dtype=float).reshape(res_xy['scan_info']['shape'])
    xs = res_xy['scan_info']['x']
    ys = res_xy['scan_info']['y']
    vabs = np.nanmax(np.abs(E2d - np.nanmean(E2d)))
    if not np.isfinite(vabs) or (vabs <= 0.0):
        raise ValueError(f'Invalid XY energy span {vabs}')
    fig, ax = plt.subplots(1, 1, figsize=(8, 7))
    im = ax.imshow((E2d - np.mean(E2d)).T, extent=[xs[0], xs[-1], ys[0], ys[-1]], origin='lower', aspect='equal', cmap='bwr', vmin=-vabs, vmax=vabs)
    ax.set_xlabel('x [Å]')
    ax.set_ylabel('y [Å]')
    ax.set_title('H2O_O on NaCl8x8 XY scan at z=1.0 Å (macro corrected, mean-shifted)')
    plt.colorbar(im, ax=ax, label='E - <E> [eV]')
    xy_png = os.path.join(OUT_DIR, 'H2O_NaCl8x8_XYscan_-5_35_npbc4_macro.png')
    fig.tight_layout()
    fig.savefig(xy_png, dpi=180)
    plt.close(fig)
    return {'line_png': line_png, 'xy_png': xy_png, 'dEx': dEx, 'dEy': dEy, 'scanner': scanner}


def test_h2o_on_nacl():
    """Test H2O on NaCl - basic small system."""
    print("=" * 60)
    print("Test 1: H2O on NaCl (MMFF REQ params from AtomTypes.dat)")
    print("=" * 60)

    mol_file = os.path.join(XYZ_DIR, 'H2O_O.xyz')
    sub_file = os.path.join(XYZ_DIR, 'NaCl_1x1_L3.xyz')

    scanner = InteractionScanner(nloc=32)
    mol_pos, mol_REQs, mol_names = scanner.load_molecule_xyz(mol_file)
    sub_pos, sub_REQs, sub_names = scanner.load_substrate_xyz(sub_file)

    print(f"Molecule: {len(mol_names)} atoms: {mol_names}")
    print(f"  positions:\n{mol_pos}")
    print(f"  REQs (R, sqrt(E), Q, H):\n{mol_REQs}")
    print(f"Substrate: {len(sub_names)} atoms: {sub_names}")
    print(f"  REQs (R, sqrt(E), Q, H):\n{sub_REQs}")

    scanner.enable_LJ      = True
    scanner.enable_Coulomb  = True
    scanner.enable_HBond    = True
    scanner.enable_Morse    = False

    # -------- PBC invariance sanity check --------
    # With finite image sums (nPBC), invariance under translation by lattice vectors is only guaranteed
    # if we wrap the molecule translation into the primary cell (InteractionScanner.wrap_PBC=True).
    print("\n--- 0. PBC invariance check (E(t) == E(t+a) == E(t+b)) ---")
    assert scanner.lvec is not None, "NaCl_1x1_L3.xyz must provide lattice vectors in comment line"
    a = np.array(scanner.lvec[0], dtype=float)
    b = np.array(scanner.lvec[1], dtype=float)
    pos = np.array([0.1, 0.2, 3.0], dtype=float)
    E0 = scanner.evaluate_single(pos=pos)['total']
    Ea = scanner.evaluate_single(pos=pos + a)['total']
    Eb = scanner.evaluate_single(pos=pos + b)['total']
    dEa = float(Ea - E0)
    dEb = float(Eb - E0)
    print(f"  E0={E0:+.6e}  Ea={Ea:+.6e}  dEa={dEa:+.3e}  Eb={Eb:+.6e}  dEb={dEb:+.3e}")
    assert abs(dEa) < 1e-5 and abs(dEb) < 1e-5, "PBC invariance failed; check wrap_PBC or kernel PBC loops"

    # -------- 1. Z-approach curve --------
    print("\n--- 1. Z-approach scan ---")
    res_z = scanner.scan_z(pos_xy=(0.0, 0.0), z_range=(1.5, 8.0), nz=60)
    zs = res_z['scan_info']['z']
    imin = np.argmin(res_z['total'])
    print(f"  Min total energy: {res_z['total'][imin]:.4f} eV at z = {zs[imin]:.3f} Ang")
    print(f"    LJ:      {res_z['LJ'][imin]:.4f}")
    print(f"    Coulomb: {res_z['Coulomb'][imin]:.4f}")
    print(f"    HBond:   {res_z['HBond'][imin]:.4f}")

    # -------- 2. Lateral 2D scan --------
    print("\n--- 2. Lateral (x,y) scan at z=3.0 ---")
    res_xy = scanner.scan_lateral(z=3.0, x_range=(0, 4), y_range=(0, 4), nx=30, ny=30)
    E2d = res_xy['total'].reshape(res_xy['scan_info']['shape'])
    imin2d = np.unravel_index(np.argmin(E2d), E2d.shape)
    xs, ys = res_xy['scan_info']['x'], res_xy['scan_info']['y']
    print(f"  Min energy: {E2d[imin2d]:.4f} eV at x={xs[imin2d[0]]:.2f}, y={ys[imin2d[1]]:.2f}")

    # -------- 3. Rotation scan --------
    print("\n--- 3. Rotation scan (z-axis) at pos=(0,0,3) ---")
    res_rot = scanner.scan_rotation(pos=(0.0, 0.0, 3.0), axis=(0,0,1), nrot=36)
    angles = res_rot['scan_info']['angles']
    imin_rot = np.argmin(res_rot['total'])
    print(f"  Min energy: {res_rot['total'][imin_rot]:.4f} eV at angle={angles[imin_rot]:.1f} deg")

    # -------- 4. SLERP path --------
    print("\n--- 4. SLERP path scan ---")
    q0 = ScanUtils.quat_from_axis_angle([0,0,1], 0)
    q1 = ScanUtils.quat_from_axis_angle([0,0,1], np.pi)
    res_sl = scanner.scan_slerp(q0, q1, t0=[0,0,3], t1=[4,0,3], npts=40)
    imin_sl = np.argmin(res_sl['total'])
    print(f"  Min energy along SLERP path: {res_sl['total'][imin_sl]:.4f} eV at t={res_sl['scan_info']['t'][imin_sl]:.3f}")

    # -------- 5. Constrained relaxation --------
    print("\n--- 5. Z-approach with constrained relaxation ---")
    scanner.spring_k    = np.float32(5.0)
    scanner.relax_dt    = np.float32(0.005)
    scanner.relax_nsteps = 100
    res_relax = scanner.scan_z(pos_xy=(0.0, 0.0), z_range=(2.0, 6.0), nz=40, relax=True)
    zs_r = res_relax['scan_info']['z']
    imin_r = np.argmin(res_relax['total'])
    print(f"  Min relaxed energy: {res_relax['total'][imin_r]:.4f} eV at z = {zs_r[imin_r]:.3f}")

    return scanner, res_z, res_xy, res_rot, res_sl, res_relax


def test_ptcda_on_nacl():
    """Test PTCDA on NaCl - larger molecule with proper atom typing."""
    print("\n" + "=" * 60)
    print("Test 2: PTCDA on NaCl (38 atoms, proper atom typing)")
    print("=" * 60)

    mol_file = os.path.join(XYZ_DIR, 'PTCDA.xyz')
    sub_file = os.path.join(XYZ_DIR, 'NaCl_8x8_L3.xyz')  # 384 atoms

    # PTCDA atom type mapping:
    # - All C are aromatic/conjugated -> C_R
    # - Bridge O (ether in anhydride) -> O_3
    # - Carbonyl O -> O_2
    # - H bonded to aromatic C -> H
    # But since we load from xyz (only element names), we need a type_map
    # The xyz has 24 C, 6 O, 8 H
    # We map C->C_R for all carbons (simplification; some are C_2 carbonyl)
    ptcda_type_map = {'C': 'C_R', 'O': 'O_2'}  # O_2 works for both ether and carbonyl vdW

    scanner = InteractionScanner(nloc=32)
    mol_pos, mol_REQs, mol_names = scanner.load_molecule_xyz(mol_file, type_map=ptcda_type_map)
    sub_pos, sub_REQs, sub_names = scanner.load_substrate_xyz(sub_file)

    print(f"Molecule: {len(mol_names)} atoms")
    print(f"  Elements: {set(mol_names)}")
    print(f"  REQs sample (first 5):\n{mol_REQs[:5]}")
    print(f"Substrate: {len(sub_names)} atoms")
    print(f"  Elements: {set(sub_names)}")

    scanner.enable_LJ      = True
    scanner.enable_Coulomb  = True
    scanner.enable_HBond    = False
    scanner.enable_Morse    = False

    # -------- Z-approach --------
    print("\n--- PTCDA Z-approach scan ---")
    res_z = scanner.scan_z(pos_xy=(0.0, 0.0), z_range=(2.5, 10.0), nz=80)
    zs = res_z['scan_info']['z']
    imin = np.argmin(res_z['total'])
    print(f"  Min total energy: {res_z['total'][imin]:.4f} eV at z = {zs[imin]:.3f} Ang")
    print(f"    LJ:      {res_z['LJ'][imin]:.4f}")
    print(f"    Coulomb: {res_z['Coulomb'][imin]:.4f}")

    # -------- Lateral 2D scan --------
    print("\n--- PTCDA Lateral scan at z=3.5 ---")
    res_xy = scanner.scan_lateral(z=3.5, x_range=(-5, 5), y_range=(-5, 5), nx=40, ny=40)
    E2d = res_xy['total'].reshape(res_xy['scan_info']['shape'])
    imin2d = np.unravel_index(np.argmin(E2d), E2d.shape)
    xs, ys = res_xy['scan_info']['x'], res_xy['scan_info']['y']
    print(f"  Min energy: {E2d[imin2d]:.4f} eV at x={xs[imin2d[0]]:.2f}, y={ys[imin2d[1]]:.2f}")

    # -------- Rotation scan --------
    print("\n--- PTCDA Rotation scan at z=3.5 ---")
    res_rot = scanner.scan_rotation(pos=(0, 0, 3.5), axis=(0,0,1), nrot=72)
    angles = res_rot['scan_info']['angles']
    imin_rot = np.argmin(res_rot['total'])
    print(f"  Min energy: {res_rot['total'][imin_rot]:.4f} eV at angle={angles[imin_rot]:.1f} deg")
    print(f"  Energy range: [{res_rot['total'].min():.4f}, {res_rot['total'].max():.4f}] eV")

    # -------- Single pose evaluation --------
    print("\n--- PTCDA single pose ---")
    e_single = scanner.evaluate_single(pos=(0, 0, 3.5))
    print(f"  E_total={e_single['total']:.4f}  E_LJ={e_single['LJ']:.4f}  E_Coul={e_single['Coulomb']:.4f}")

    return scanner, res_z, res_xy, res_rot


def make_argparser():
    parser = argparse.ArgumentParser(description='Headless interaction scan validation and fast-path parity runner')
    parser.add_argument('--mol', type=str, default=None, help='Molecule xyz file or filename in common_resources/xyz')
    parser.add_argument('--sub', type=str, default=None, help='Substrate xyz file or filename in common_resources/xyz')
    parser.add_argument('--mol-preset', choices=sorted(MOL_PRESETS.keys()), default='h2o', help='Preset molecule when --mol is not given')
    parser.add_argument('--sub-preset', choices=sorted(SUB_PRESETS.keys()), default='nacl_8x8', help='Preset substrate when --sub is not given')
    parser.add_argument('--batch-presets', action='store_true', help='Run H2O/PTCDA on NaCl 1x1 and 8x8')
    parser.add_argument('--mode', choices=['reference', 'fast', 'all'], default='all', help='Which scan backend(s) to run')
    parser.add_argument('--nx', type=int, default=121, help='2D scan resolution in x')
    parser.add_argument('--ny', type=int, default=121, help='2D scan resolution in y')
    parser.add_argument('--ns-line', type=int, default=9, help='Number of equivalent-site points in 1D line check')
    parser.add_argument('--z0', type=float, default=None, help='Scan height override')
    parser.add_argument('--npbc', type=int, nargs=3, default=(4, 4, 0), help='PBC replication used for macro correction')
    parser.add_argument('--fast-nsystems', type=int, default=8192, help='Number of replicas allocated for fast GPU batch')
    parser.add_argument('--fast-chunk', type=int, default=8192, help='Chunk size used during fast GPU evaluation')
    parser.add_argument('--out-dir', type=str, default=OUT_DIR, help='Output directory for PNGs')
    parser.add_argument('--legacy-suite', action='store_true', help='Also run the older H2O/PTCDA demo suite')
    return parser


def main(argv=None):
    parser = make_argparser()
    args = parser.parse_args(argv)
    if not args.batch_presets:
        if args.mol is None:
            args.mol = MOL_PRESETS[args.mol_preset][0]
            if args.z0 is None:
                args.z0 = MOL_PRESETS[args.mol_preset][2]
        if args.sub is None:
            args.sub = SUB_PRESETS[args.sub_preset]
    print('\n' + '=' * 60)
    print('Headless interaction scan runner')
    print('=' * 60)
    print(f'mode={args.mode} batch_presets={args.batch_presets} nx={args.nx} ny={args.ny} npbc={tuple(args.npbc)} out_dir={args.out_dir}')
    results = run_batch_cases(args)
    if args.legacy_suite:
        print('\nRunning legacy demonstration suite ...')
        test_h2o_on_nacl()
        test_ptcda_on_nacl()
    print('\n' + '=' * 60)
    print(f'Completed {len(results)} scan scenario runs successfully')
    print('=' * 60)
    return results


if __name__ == '__main__':
    main()
