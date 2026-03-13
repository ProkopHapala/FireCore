#!/usr/bin/env python3
import os, sys, json, argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from pyBall.OCL import ScanUtils
from pyBall.OCL.MolecularDynamics import MolecularDynamics
from pyBall.OCL.InteractionEnergy import InteractionScanner, make_REQs_from_enames, load_xyz_with_REQs
from pyBall.OCL import MMparams

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
XYZ_DIR = os.path.join(ROOT, 'cpp', 'common_resources', 'xyz')
SUB_FILE = os.path.join(XYZ_DIR, 'NaCl_8x8_L3.xyz')
OUT_DIR = os.path.join(os.path.dirname(__file__), 'output_flat_folded_components_periodic')
COMPONENTS = ('pauli', 'london', 'coulomb')
COMP_LABELS = {'pauli': 'Pauli', 'london': 'London', 'coulomb': 'Coulomb', 'total': 'Total'}
COLORS = {'reference': 'k', 'embedding': 'tab:blue', 'folded': 'tab:red'}
FOLDED_KERNELS = ('orig', 'harmonics', 'workgroup')

def ensure_dir(path): os.makedirs(path, exist_ok=True)

def atom_types_data():
    d = os.path.join(ROOT, 'cpp', 'common_resources')
    et = MMparams.read_element_types(os.path.join(d, 'ElementTypes.dat'))
    return MMparams.read_atom_types(os.path.join(d, 'AtomTypes.dat'), et)

def make_probe_REQ(atom_name, q, atom_types):
    return np.array(make_REQs_from_enames([atom_name], np.array([q], dtype=np.float32), atom_types)[0], dtype=np.float32)

def rmse(a, b):
    d = np.asarray(a, dtype=np.float64) - np.asarray(b, dtype=np.float64)
    return float(np.sqrt(np.mean(d*d)))

def arr_stats(arr):
    arr = np.asarray(arr, dtype=np.float64)
    return {'min': float(np.min(arr)), 'max': float(np.max(arr))}

def err_stats(ref, test):
    ref = np.asarray(ref, dtype=np.float64); test = np.asarray(test, dtype=np.float64); d = test - ref
    return {'min': float(np.min(d)), 'max': float(np.max(d)), 'max_abs': float(np.max(np.abs(d))), 'rmse': rmse(ref, test)}

def surface_info(sub_file):
    apos, REQs, enames, Zs, lvec = load_xyz_with_REQs(sub_file)
    z_top = float(np.max(apos[:, 2])); a = np.array(lvec[0, :3], dtype=np.float64); b = np.array(lvec[1, :3], dtype=np.float64)
    sites = {'na': np.array([0.0, 0.0], dtype=np.float32), 'cl': np.array([2.0, 2.0], dtype=np.float32)}
    return {'z_top': z_top, 'lvec': np.array(lvec, dtype=np.float32), 'a': a, 'b': b, 'sites': sites, 'apos': apos, 'REQs': REQs, 'enames': enames}

def setup_md(req, nSystems, sub_file, nPBC, bMacro=False):
    md = MolecularDynamics(nloc=32, debug_build_options='-DDBG_UFF=0')
    md.init_rigid_molecule_batch(np.zeros((1, 3), dtype=np.float32), req[None, :], nSystems=nSystems)
    md.set_surface(sub_file, nPBC=nPBC, pos0=(0.0, 0.0, 0.0), alpha_morse=1.8, r_damp=0.0, bMacro=bMacro)
    return md

def setup_macro_scanner(req, sub_file, nPBC):
    sc = InteractionScanner(nloc=32)
    sc.set_molecule(np.zeros((1, 3), dtype=np.float32), req[None, :])
    sc.load_substrate_xyz(sub_file)
    sc.enable_LJ = False; sc.enable_Coulomb = True; sc.enable_HBond = False; sc.enable_Morse = False
    sc.enable_macro = True; sc.nPBC[:] = tuple(int(v) for v in nPBC); sc.wrap_PBC = True; sc._update_macro_from_substrate()
    return sc

def pack_transforms(xyz):
    xyz = np.asarray(xyz, dtype=np.float32).reshape(-1, 3)
    R = [np.eye(3, dtype=np.float32)] * len(xyz)
    return ScanUtils.pack_transforms(R, [p for p in xyz])

def make_z_transforms(xy, z_top, z_heights):
    zs = np.asarray(z_top + z_heights, dtype=np.float32)
    xyz = np.column_stack([np.full_like(zs, xy[0]), np.full_like(zs, xy[1]), zs])
    return pack_transforms(xyz), zs

def make_lateral_line(info, z_abs, nxy=240, ncell_span=3.0, axis='a'):
    ts = np.linspace(-0.5, -0.5 + ncell_span, int(nxy), dtype=np.float32)
    if axis == 'a':
        step = info['a'][:2]
    elif axis == 'b':
        step = info['b'][:2]
    else:
        raise ValueError(f"make_lateral_line(): unknown axis='{axis}', expected 'a' or 'b'")
    xy = np.outer(ts, step)
    xyz = np.column_stack([xy[:, 0], xy[:, 1], np.full(len(ts), z_abs, dtype=np.float32)])
    return pack_transforms(xyz), ts

def apply_macro(scanner, transforms, local_components):
    T = scanner._wrap_transforms_PBC(np.asarray(transforms, dtype=np.float32))
    out = {'total': np.array(local_components['total'], dtype=np.float64).copy(), 'LJ': np.zeros(len(local_components['total']), dtype=np.float64), 'Coulomb': np.array(local_components['coulomb'], dtype=np.float64).copy(), 'HBond': np.zeros(len(local_components['total']), dtype=np.float64)}
    out2 = scanner._apply_macro_correction(T, out)
    corrected = {'pauli': np.array(local_components['pauli'], dtype=np.float64), 'london': np.array(local_components['london'], dtype=np.float64), 'coulomb': np.array(out2['Coulomb'], dtype=np.float64)}
    corrected['total'] = corrected['pauli'] + corrected['london'] + corrected['coulomb']
    corrected['macro_delta'] = corrected['coulomb'] - np.array(local_components['coulomb'], dtype=np.float64)
    corrected['timing'] = dict(local_components.get('timing', {}))
    return corrected

def eval_md_components(md, transforms, chunk_size):
    out = md.eval_rigid_getSurfMorse_components(transforms, chunk_size=chunk_size, components=COMPONENTS)
    res = {k: np.array(out[k], dtype=np.float64) for k in COMPONENTS}
    res['total'] = np.array(out['total'], dtype=np.float64)
    res['timing'] = {k: float(out.get(k, 0.0)) for k in ('t_prep_s', 't_kernel_s', 't_download_s', 't_total_s')}
    return res

def eval_folded_components(md, transforms, chunk_size):
    out = {}; timing = None
    for ck in COMPONENTS:
        res = md.eval_rigid_getSurfFolded(transforms, chunk_size=chunk_size, component=ck)
        out[ck] = np.array(res['total'], dtype=np.float64); timing = res
    baselines = getattr(md, 'folded_baselines', {})
    for ck in COMPONENTS:
        if ck in baselines:
            out[ck] -= float(baselines[ck])
    out['total'] = out['pauli'] + out['london'] + out['coulomb']
    out['timing'] = {k: float(timing.get(k, 0.0)) for k in ('t_prep_s', 't_kernel_s', 't_download_s', 't_total_s')}
    out['macro_delta'] = np.zeros(len(out['total']), dtype=np.float64)
    return out

def build_fit_mask_and_weights(z_sample, z_min, E_sample, repel_cut, weight_k):
    z_sample = np.asarray(z_sample, dtype=np.float64)
    E_sample = np.asarray(E_sample, dtype=np.float64)
    mask = np.asarray(z_sample >= float(z_min), dtype=bool)
    Emin = float(np.min(E_sample)) if len(E_sample) else 0.0
    w = np.exp(np.clip(-float(weight_k) * (E_sample - Emin), -60.0, 60.0))
    if repel_cut is not None:
        w[E_sample > float(repel_cut)] = 0.0
        mask &= np.asarray(E_sample <= float(repel_cut), dtype=bool)
    return mask, w

def build_folded(req, sub_file, nSystems, z_abs_range, nPBC_small, z_fit_min=1.5, repel_cut=0.2, weight_power=12.0, folded_nxy=24, folded_nz=96, folded_nu=2, folded_nv=2, folded_nzbasis=8, kernel_kind='orig'):
    md = setup_md(req, nSystems=nSystems, sub_file=sub_file, nPBC=nPBC_small, bMacro=False)
    md.set_folded_kernel_kind(kernel_kind)
    ref_scan = setup_macro_scanner(req, sub_file, nPBC_small)
    info = surface_info(sub_file)
    z_top = float(info['z_top'])
    md.fit_folded_surface_basis(sub_file, nPBC=nPBC_small, z_range=z_abs_range, nu=folded_nu, nv=folded_nv, nz=folded_nzbasis, nxy=folded_nxy, nz_samp=folded_nz, r_damp=0.0, alpha_morse=1.8, bMacro=False, components=COMPONENTS, fit_mask=None, weight_power=0.0)
    uvz = np.asarray(md.folded_fit_info['uvz'], dtype=np.float64)
    a = md.folded_lvec_basis[0, :3].astype(np.float64)
    b = md.folded_lvec_basis[1, :3].astype(np.float64)
    xyz = np.empty((len(uvz), 3), dtype=np.float32)
    xyz[:, 0] = uvz[:, 0]*a[0] + uvz[:, 1]*b[0]
    xyz[:, 1] = uvz[:, 0]*a[1] + uvz[:, 1]*b[1]
    xyz[:, 2] = uvz[:, 2]
    transforms = pack_transforms(xyz)
    tref = apply_macro(ref_scan, transforms, eval_md_components(md, transforms, min(len(xyz), md.nSystems)))
    zrel = xyz[:, 2] - float(z_abs_range[0])
    z_thresh = np.percentile(zrel, 95.0)
    far_mask = zrel >= z_thresh
    baselines = {}
    for ck in COMPONENTS:
        baselines[ck] = float(np.mean(tref[ck][far_mask])) if np.any(far_mask) else 0.0
        tref[ck] = np.asarray(tref[ck], dtype=np.float64) - baselines[ck]
    tref['total'] = tref['pauli'] + tref['london'] + tref['coulomb']
    baselines['total'] = baselines['pauli'] + baselines['london'] + baselines['coulomb']
    z_above_top = xyz[:, 2] - z_top
    fit_mask, wmask = build_fit_mask_and_weights(z_above_top, z_fit_min, tref['total'], repel_cut, weight_power)
    Phi = np.asarray(md.folded_fit_info['Phi'], dtype=np.float64)
    ww = np.asarray(wmask, dtype=np.float64) * np.asarray(fit_mask, dtype=np.float64)
    if not np.any(ww > 0.0):
        raise ValueError('build_folded(): all fit weights are zero')
    coeff_sets = {}
    nbasis = Phi.shape[1]
    ntypes = int(md.folded_params['unique_REQs'].shape[0])
    if ntypes != 1:
        raise ValueError(f'build_folded(): expected single probe type for this harness, got ntypes={ntypes}')
    ws = np.sqrt(ww)
    Phiw = Phi * ws[:, None]
    for ck in COMPONENTS:
        yw = np.asarray(tref[ck], dtype=np.float64) * ws
        S, *_ = np.linalg.lstsq(Phiw, yw, rcond=None)
        coeff = np.zeros((ntypes, nbasis), dtype=np.float32)
        coeff[0, :] = S.astype(np.float32)
        coeff_sets[ck] = coeff
    md.folded_params['coeff_sets'] = coeff_sets
    md.folded_params['coeffs'] = coeff_sets[COMPONENTS[0]].copy()
    md._set_folded_coefficients(coeff_sets[COMPONENTS[0]])
    md.folded_fit_info['fit_ref_total'] = np.asarray(tref['total'], dtype=np.float64)
    md.folded_fit_info['repel_cut'] = float(repel_cut)
    md.folded_fit_info['z_fit_min'] = float(z_fit_min)
    md.folded_fit_info['z_top'] = z_top
    md.folded_fit_info['kernel_kind'] = str(kernel_kind)
    md.folded_fit_info['weights'] = [wmask]
    md.folded_fit_info['fit_mask'] = fit_mask
    md.folded_fit_info['baselines'] = baselines
    md.folded_baselines = baselines
    return md

def make_weight_plot(path_png, md, title):
    uvz = md.folded_fit_info['uvz']; ww = np.asarray(md.folded_fit_info['weights'][0], dtype=np.float64); mask = np.asarray(md.folded_fit_info['fit_mask'], dtype=bool); E = np.asarray(md.folded_fit_info.get('fit_ref_total', np.zeros(len(uvz))), dtype=np.float64)
    z = uvz[:, 2]
    zuniq = np.unique(np.round(z, 6))
    Ez = np.zeros(len(zuniq)); wz = np.zeros(len(zuniq)); mz = np.zeros(len(zuniq))
    for i, zz in enumerate(zuniq):
        m = np.abs(z - zz) < 1e-6
        Ez[i] = np.mean(E[m]) if np.any(m) else 0.0
        wz[i] = np.mean(ww[m]) if np.any(m) else 0.0
        mz[i] = np.any(mask[m]).astype(float)
    fig, ax = plt.subplots(1, 1, figsize=(8, 4.5))
    ax.plot(zuniq, Ez, '-', lw=1.8, c='k', label='Reference total (macro-corrected)')
    ax.axhline(md.folded_fit_info.get('repel_cut', 0.0), ls='--', lw=1.0, c='tab:red', label='Repulsion cutoff')
    ax.axvline(float(md.folded_fit_info.get('z_top', float(zuniq.min()))) + md.folded_fit_info.get('z_fit_min', 0.0), ls='--', lw=1.0, c='tab:green', label='z fit min')
    ax.set_xlabel('absolute z [Å]'); ax.set_ylabel('Energy [eV]'); ax.grid(True, alpha=0.25)
    ax2 = ax.twinx()
    ax2.plot(zuniq, wz, '-', lw=1.4, c='tab:purple', label='Fit weight')
    for axx in (ax, ax2):
        axx.legend(loc='upper right', frameon=False)
    ax.set_title(title)
    fig.tight_layout(); fig.savefig(path_png, dpi=180); plt.close(fig)

def compute_weights_only(req, sub_file, z_abs_range, nPBC_small, z_fit_min, repel_cut, weight_power, folded_nxy=16, folded_nz=80):
    # Build sample grid first
    xyz_tmp = []
    z_samples = np.linspace(float(z_abs_range[0]), float(z_abs_range[1]), int(folded_nz), endpoint=True, dtype=np.float32)
    for z in z_samples:
        for iy in range(int(folded_nxy)):
            v = (iy + 0.5) / float(folded_nxy) - 0.5
            for ix in range(int(folded_nxy)):
                u = (ix + 0.5) / float(folded_nxy) - 0.5
                xyz_tmp.append((u, v, float(z)))
    xyz_tmp = np.array(xyz_tmp, dtype=np.float32)
    # Instantiate MD with enough systems to batch upload
    md = setup_md(req, nSystems=min(len(xyz_tmp), 4096), sub_file=sub_file, nPBC=nPBC_small, bMacro=False)
    ref_scan = setup_macro_scanner(req, sub_file, nPBC_small)
    info = surface_info(sub_file)
    z_top = float(info['z_top'])
    # Map fractional (u,v) to Cartesian using surface_lvec
    a = md.surface_lvec[0, :3].astype(np.float64); b = md.surface_lvec[1, :3].astype(np.float64)
    xyz = np.empty_like(xyz_tmp)
    xyz[:, 0] = xyz_tmp[:, 0]*a[0] + xyz_tmp[:, 1]*b[0]
    xyz[:, 1] = xyz_tmp[:, 0]*a[1] + xyz_tmp[:, 1]*b[1]
    xyz[:, 2] = xyz_tmp[:, 2]
    T = pack_transforms(xyz)
    chunk = md.nSystems
    tref = apply_macro(ref_scan, T, eval_md_components(md, T, chunk))
    fit_mask, wmask = build_fit_mask_and_weights(xyz[:, 2] - z_top, z_fit_min, tref['total'], repel_cut, weight_power)
    md.folded_fit_info = {
        'uvz': np.column_stack([np.zeros(len(xyz)), np.zeros(len(xyz)), xyz[:, 2]]),
        'fit_ref_total': np.asarray(tref['total'], dtype=np.float64),
        'repel_cut': float(repel_cut),
        'z_fit_min': float(z_fit_min),
        'z_top': z_top,
        'weights': [wmask],
        'fit_mask': fit_mask,
    }
    return md

def save_case_numeric(path_csv, x, xname, ref, emb, fold):
    cols = [np.asarray(x, dtype=np.float64)]; names = [xname]
    for tag, data in [('ref', ref), ('emb', emb), ('fold', fold)]:
        for ck in ('pauli', 'london', 'coulomb', 'total'):
            cols.append(np.asarray(data[ck], dtype=np.float64)); names.append(f'{tag}_{ck}')
    for tag, data in [('emb', emb), ('fold', fold)]:
        for ck in ('pauli', 'london', 'coulomb', 'total'):
            cols.append(np.asarray(data[ck], dtype=np.float64) - np.asarray(ref[ck], dtype=np.float64)); names.append(f'err_{tag}_{ck}')
    np.savetxt(path_csv, np.column_stack(cols), delimiter=',', header=','.join(names), comments='')

def make_scan_plot(path_png, title, x, xlabel, ref, emb, fold, zoom=None):
    fig, axs = plt.subplots(2, 2, figsize=(14, 10), sharex=True); comps = ('pauli', 'london', 'coulomb', 'total')
    for ax, ck in zip(axs.flat, comps):
        ax.plot(x, ref[ck], '-', c=COLORS['reference'], lw=2.0, label='Reference')
        ax.plot(x, emb[ck], '-', c=COLORS['embedding'], lw=1.4, label='Embedding')
        ax.plot(x, fold[ck], '-', c=COLORS['folded'], lw=1.4, label='Folded')
        if zoom is not None: ax.set_ylim(*zoom)
        ax2 = ax.twinx()
        ax2.plot(x, emb[ck] - ref[ck], '--', c=COLORS['embedding'], lw=0.9, alpha=0.85, label='Emb-Ref')
        ax2.plot(x, fold[ck] - ref[ck], '--', c=COLORS['folded'], lw=0.9, alpha=0.85, label='Fold-Ref')
        ax.set_title(COMP_LABELS[ck]); ax.grid(True, alpha=0.25); ax.set_ylabel('Energy [eV]'); ax2.set_ylabel('Error [eV]')
    axs[1, 0].set_xlabel(xlabel); axs[1, 1].set_xlabel(xlabel)
    h1, l1 = axs[0, 0].get_legend_handles_labels(); h2, l2 = axs[0, 0].twinx().get_legend_handles_labels()
    fig.suptitle(title); fig.legend(h1 + h2, l1 + l2, loc='upper center', ncol=5, frameon=False); fig.tight_layout(rect=(0, 0, 1, 0.95)); fig.savefig(path_png, dpi=180); plt.close(fig)

def print_case_stats(prefix, ref, emb, fold):
    for ck in ('pauli', 'london', 'coulomb', 'total'):
        ee = err_stats(ref[ck], emb[ck]); ef = err_stats(ref[ck], fold[ck]); sr = arr_stats(ref[ck]); sf = arr_stats(fold[ck])
        print(f"[{prefix}] {ck:7s} ref[min,max]=({sr['min']:.6e},{sr['max']:.6e}) fold[min,max]=({sf['min']:.6e},{sf['max']:.6e}) err_emb[maxabs,rmse]=({ee['max_abs']:.6e},{ee['rmse']:.6e}) err_fold[maxabs,rmse]=({ef['max_abs']:.6e},{ef['rmse']:.6e})")

def make_record(kind, probe_name, q, site_name, ref, emb, fold, timing_ref, timing_emb, timing_fold, files, focus_mask=None, kernel_kind='orig'):
    rec = {'kind': kind, 'probe': probe_name, 'charge': float(q), 'site': site_name, 'kernel': str(kernel_kind), 'files': files, 'timing_s': {'reference': timing_ref, 'embedding': timing_emb, 'folded': timing_fold}, 'components': {}}
    for ck in ('pauli', 'london', 'coulomb', 'total'):
        rec['components'][ck] = {'reference': arr_stats(ref[ck]), 'embedding': arr_stats(emb[ck]), 'folded': arr_stats(fold[ck]), 'embedding_error': err_stats(ref[ck], emb[ck]), 'folded_error': err_stats(ref[ck], fold[ck])}
        if focus_mask is not None:
            m = np.asarray(focus_mask, dtype=bool)
            rec['components'][ck]['embedding_error_focus'] = err_stats(np.asarray(ref[ck])[m], np.asarray(emb[ck])[m])
            rec['components'][ck]['folded_error_focus'] = err_stats(np.asarray(ref[ck])[m], np.asarray(fold[ck])[m])
    return rec

def make_benchmark_record(probe_name, q, kernel_kind, timing, nconf):
    tt = timing.get('timing', timing)
    return {'probe': probe_name, 'charge': float(q), 'kernel': str(kernel_kind), 'nconf': int(nconf), 'timing_s': {k: float(tt.get(k, 0.0)) for k in ('t_prep_s', 't_kernel_s', 't_download_s', 't_total_s')}}

def make_benchmark_plot(path_png, bench_records, title):
    kernels = [r['kernel'] for r in bench_records]
    vals = [r['timing_s']['t_total_s'] for r in bench_records]
    fig, ax = plt.subplots(1, 1, figsize=(7, 4.5))
    ax.bar(np.arange(len(kernels)), vals, color=['tab:red', 'tab:orange', 'tab:green'][:len(kernels)])
    ax.set_xticks(np.arange(len(kernels)), kernels)
    ax.set_ylabel('Wall time [s]')
    ax.set_title(title)
    ax.grid(True, axis='y', alpha=0.25)
    fig.tight_layout(); fig.savefig(path_png, dpi=180); plt.close(fig)

def run_kernel_benchmarks(req, info, z_abs_range, nPBC_small, kernel_kinds, probe_name, q, nconf):
    xyz = []
    nu = 2
    nv = 2
    nzb = 8
    nxy = int(np.sqrt(max(4, int(nconf // max(1, len(np.linspace(z_abs_range[0], z_abs_range[1], 8)))))))
    zs = np.linspace(float(z_abs_range[0]), float(z_abs_range[1]), 8, endpoint=True, dtype=np.float32)
    for z in zs:
        for iy in range(nxy):
            v = (iy + 0.5) / float(nxy) - 0.5
            for ix in range(nxy):
                u = (ix + 0.5) / float(nxy) - 0.5
                p = u*info['a'] + v*info['b']
                xyz.append((p[0], p[1], float(z)))
                if len(xyz) >= int(nconf):
                    break
            if len(xyz) >= int(nconf):
                break
        if len(xyz) >= int(nconf):
            break
    transforms = pack_transforms(np.asarray(xyz, dtype=np.float32))
    benches = []
    for kernel_kind in kernel_kinds:
        md = build_folded(req, SUB_FILE, min(len(xyz), 4096), z_abs_range, nPBC_small, folded_nu=nu, folded_nv=nv, folded_nzbasis=nzb, kernel_kind=kernel_kind)
        out = eval_folded_components(md, transforms, chunk_size=min(len(xyz), md.nSystems))
        benches.append(make_benchmark_record(probe_name, q, kernel_kind, out, len(xyz)))
    return benches

def save_summary_markdown(path_md, records, config):
    lines = ['# Rigorous periodic flat-surface folded-component test summary', '', '## Configuration', '']
    for k, v in config.items(): lines.append(f'- **{k}**: `{v}`')
    lines += ['', '## Main observations', '']
    total_max_fold = max(r['components']['total']['folded_error']['max_abs'] for r in records)
    total_max_emb = max(r['components']['total']['embedding_error']['max_abs'] for r in records)
    lines.append(f'- **Worst total folded max|ΔE|**: `{total_max_fold:.6e} eV`')
    lines.append(f'- **Worst total embedding max|ΔE|**: `{total_max_emb:.6e} eV`')
    for r in records:
        lines += ['', f"## {r['kind']} : {r['probe']} q={r['charge']:+.3f} site={r['site']} kernel={r.get('kernel', 'orig')}", '']
        lines.append(f"- **PNG**: `{r['files']['png']}`"); lines.append(f"- **CSV**: `{r['files']['csv']}`"); lines.append(f"- **JSON**: `{r['files']['json']}`")
        for ck in ('pauli', 'london', 'coulomb', 'total'):
            ce = r['components'][ck]['embedding_error']; cf = r['components'][ck]['folded_error']
            lines.append(f"- **{ck}**: emb max|Δ|=`{ce['max_abs']:.6e}` rmse=`{ce['rmse']:.6e}` ; fold max|Δ|=`{cf['max_abs']:.6e}` rmse=`{cf['rmse']:.6e}`")
            cef = r['components'][ck].get('embedding_error_focus', None); cff = r['components'][ck].get('folded_error_focus', None)
            if (cef is not None) and (cff is not None):
                lines.append(f"  - focus: emb max|Δ|=`{cef['max_abs']:.6e}` rmse=`{cef['rmse']:.6e}` ; fold max|Δ|=`{cff['max_abs']:.6e}` rmse=`{cff['rmse']:.6e}`")
    if 'benchmarks' in config:
        lines += ['', '## Folded kernel benchmark', '']
        for b in config['benchmarks']:
            lines.append(f"- **{b['probe']} q={b['charge']:+.3f} kernel={b['kernel']} n={b['nconf']}**: total=`{b['timing_s']['t_total_s']:.6e}` prep=`{b['timing_s']['t_prep_s']:.6e}` kernel=`{b['timing_s']['t_kernel_s']:.6e}` download=`{b['timing_s']['t_download_s']:.6e}`")
    open(path_md, 'w').write('\n'.join(lines) + '\n')

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument('--zmin', type=float, default=1.6)
    ap.add_argument('--zmax', type=float, default=8.0)
    ap.add_argument('--nz', type=int, default=120)
    ap.add_argument('--zfit-min', type=float, default=1.5)
    ap.add_argument('--weight-k', type=float, default=12.0)
    ap.add_argument('--repel-cut', type=float, default=0.2)
    ap.add_argument('--lateral-z', type=float, default=2.0)
    ap.add_argument('--nxy-line', type=int, default=320)
    ap.add_argument('--folded-kernels', type=str, default='orig,harmonics,workgroup')
    ap.add_argument('--bench-n', type=int, default=32768)
    ap.add_argument('--weights-only', action='store_true', help='Only compute reference total + weights and save weight plot; skip fitting/scans')
    ap.add_argument('--z-only', action='store_true', help='Run only z-scans; skip lateral scans and benchmarks')
    return ap.parse_args()

def main():
    args = parse_args()
    ensure_dir(OUT_DIR)
    kernel_kinds = [s.strip().lower() for s in str(args.folded_kernels).split(',') if s.strip()]
    for k in kernel_kinds:
        if k not in FOLDED_KERNELS:
            raise ValueError(f"unknown folded kernel '{k}', expected one of {FOLDED_KERNELS}")
    info = surface_info(SUB_FILE); atom_types = atom_types_data(); z_heights = np.linspace(args.zmin, args.zmax, args.nz, dtype=np.float32)
    z_abs_range = (float(info['z_top'] + z_heights.min()), float(info['z_top'] + z_heights.max())); nPBC_ref = (12, 12, 0); nPBC_small = (4, 4, 0)
    probes = [('H', 0.2), ('O', 0.2)]; records = []; bench_records = []
    config = {'substrate': os.path.relpath(SUB_FILE, ROOT), 'periodic_reference': 'wrapped PBC + macro', 'z_scan_above_top_A': f'{float(z_heights.min()):.3f}..{float(z_heights.max()):.3f}', 'lateral_scan_height_A': float(args.lateral_z), 'nz': len(z_heights), 'nPBC_reference': nPBC_ref, 'nPBC_embedding_local': nPBC_small, 'folded_components': COMPONENTS, 'folded_basis': 'nu=2 nv=2 nzbasis=8 nxy=24 nzsamp=96', 'zfit_min_A': float(args.zfit_min), 'repel_cut_eV': float(args.repel_cut), 'weight_k': float(args.weight_k), 'folded_kernels': kernel_kinds, 'bench_n': int(args.bench_n)}
    for probe_name, q in probes:
        req = make_probe_REQ(probe_name, q, atom_types)
        probe_root = os.path.join(OUT_DIR, f'{probe_name}_q{q:+.2f}'.replace('+', 'p').replace('-', 'm')); ensure_dir(probe_root)
        if args.weights_only:
            mdw = compute_weights_only(req, SUB_FILE, z_abs_range, nPBC_small, args.zfit_min, args.repel_cut, args.weight_k, folded_nxy=24, folded_nz=96)
            make_weight_plot(os.path.join(probe_root, f'{probe_name}_weights_only.png'), mdw, f'{probe_name} reference + weights')
            continue
        md_ref = setup_md(req, len(z_heights), SUB_FILE, nPBC_ref, bMacro=False)
        md_emb = setup_md(req, len(z_heights), SUB_FILE, nPBC_small, bMacro=False)
        macro_ref = setup_macro_scanner(req, SUB_FILE, nPBC_ref)
        macro_emb = setup_macro_scanner(req, SUB_FILE, nPBC_small)
        for kernel_kind in kernel_kinds:
            probe_dir = os.path.join(probe_root, kernel_kind); ensure_dir(probe_dir)
            md_fold = build_folded(req, SUB_FILE, len(z_heights), z_abs_range, nPBC_small, z_fit_min=args.zfit_min, repel_cut=args.repel_cut, weight_power=args.weight_k, kernel_kind=kernel_kind)
            make_weight_plot(os.path.join(probe_dir, f'{probe_name}_{kernel_kind}_weights.png'), md_fold, f'{probe_name} folded fit weights [{kernel_kind}]')
            for site_name, xy in info['sites'].items():
                transforms, _ = make_z_transforms(xy, info['z_top'], z_heights)
                ref = apply_macro(macro_ref, transforms, eval_md_components(md_ref, transforms, len(z_heights)))
                emb = apply_macro(macro_emb, transforms, eval_md_components(md_emb, transforms, len(z_heights)))
                fold = eval_folded_components(md_fold, transforms, len(z_heights))
                prefix = f'{probe_name}_q{q:+.2f}_{kernel_kind}_{site_name}_zscan'.replace('+', 'p').replace('-', 'm')
                print_case_stats(prefix, ref, emb, fold)
                png_file = os.path.join(probe_dir, f'{prefix}.png'); csv_file = os.path.join(probe_dir, f'{prefix}.csv'); json_file = os.path.join(probe_dir, f'{prefix}.json')
                make_scan_plot(png_file, f'{probe_name} q={q:+.2f} periodic z-scan @ {site_name} [{kernel_kind}]', z_heights, 'z above top layer [Å]', ref, emb, fold, zoom=(-0.2, 0.2))
                save_case_numeric(csv_file, z_heights, 'z', ref, emb, fold)
                focus_mask = np.asarray(z_heights >= args.zfit_min, dtype=bool)
                rec = make_record('z_scan', probe_name, q, site_name, ref, emb, fold, ref['timing'], emb['timing'], fold['timing'], {'png': png_file, 'csv': csv_file, 'json': json_file}, focus_mask=focus_mask, kernel_kind=kernel_kind)
                open(json_file, 'w').write(json.dumps(rec, indent=2)); records.append(rec)
            if args.z_only:
                continue
            for axis in ('a', 'b'):
                lateral_trans, ts = make_lateral_line(info, info['z_top'] + args.lateral_z, nxy=args.nxy_line, ncell_span=3.0, axis=axis)
                md_refL = setup_md(req, len(ts), SUB_FILE, nPBC_ref, bMacro=False)
                md_embL = setup_md(req, len(ts), SUB_FILE, nPBC_small, bMacro=False)
                md_foldL = build_folded(req, SUB_FILE, len(ts), z_abs_range, nPBC_small, z_fit_min=args.zfit_min, repel_cut=args.repel_cut, weight_power=args.weight_k, kernel_kind=kernel_kind)
                refL = apply_macro(macro_ref, lateral_trans, eval_md_components(md_refL, lateral_trans, len(ts)))
                embL = apply_macro(macro_emb, lateral_trans, eval_md_components(md_embL, lateral_trans, len(ts)))
                foldL = eval_folded_components(md_foldL, lateral_trans, len(ts))
                prefix = f'{probe_name}_q{q:+.2f}_{kernel_kind}_{axis}_lateral'.replace('+', 'p').replace('-', 'm')
                print_case_stats(prefix, refL, embL, foldL)
                png_file = os.path.join(probe_dir, f'{prefix}.png'); csv_file = os.path.join(probe_dir, f'{prefix}.csv'); json_file = os.path.join(probe_dir, f'{prefix}.json')
                make_scan_plot(png_file, f'{probe_name} q={q:+.2f} periodic lateral scan along {axis} z={args.lateral_z:.2f} Å [{kernel_kind}]', ts, f'path coordinate along {axis} [primitive cells]', refL, embL, foldL, zoom=(-0.2, 0.2))
                save_case_numeric(csv_file, ts, 's_cell', refL, embL, foldL)
                rec = make_record('lateral_scan', probe_name, q, f'{axis}_periodic', refL, embL, foldL, refL['timing'], embL['timing'], foldL['timing'], {'png': png_file, 'csv': csv_file, 'json': json_file}, focus_mask=np.ones(len(ts), dtype=bool), kernel_kind=kernel_kind)
                open(json_file, 'w').write(json.dumps(rec, indent=2)); records.append(rec)
        if not args.z_only:
            probe_bench = run_kernel_benchmarks(req, info, z_abs_range, nPBC_small, kernel_kinds, probe_name, q, args.bench_n)
            bench_records.extend(probe_bench)
            make_benchmark_plot(os.path.join(probe_root, f'{probe_name}_kernel_bench.png'), probe_bench, f'{probe_name} folded kernel benchmark')
    if args.weights_only:
        print('[weights-only] completed weight plots; skipping fits/scans/summary')
        return
    config['benchmarks'] = bench_records
    summary_json = os.path.join(OUT_DIR, 'summary.json'); summary_md = os.path.join(OUT_DIR, 'summary.md')
    open(summary_json, 'w').write(json.dumps({'config': config, 'records': records, 'benchmarks': bench_records}, indent=2)); save_summary_markdown(summary_md, records, config)
    print(f'[flat-folded-periodic] saved summary {summary_json}'); print(f'[flat-folded-periodic] saved summary {summary_md}')

if __name__ == '__main__': main()
