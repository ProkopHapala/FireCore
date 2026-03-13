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
    out['total'] = out['pauli'] + out['london'] + out['coulomb']
    out['timing'] = {k: float(timing.get(k, 0.0)) for k in ('t_prep_s', 't_kernel_s', 't_download_s', 't_total_s')}
    return out

def build_fit_mask_and_weights(z_sample, z_min, E_sample=None, repel_cut=None):
    mask = np.asarray(z_sample >= float(z_min), dtype=bool)
    if E_sample is not None and repel_cut is not None:
        mask &= np.asarray(E_sample <= float(repel_cut), dtype=bool)
    return mask

def build_folded(req, sub_file, nSystems, z_abs_range, nPBC_small, z_fit_min=1.5, repel_cut=0.2, weight_power=12.0, folded_nxy=24, folded_nz=96, folded_nu=2, folded_nv=2, folded_nzbasis=12):
    md = setup_md(req, nSystems=nSystems, sub_file=sub_file, nPBC=nPBC_small, bMacro=False)
    ref_scan = setup_macro_scanner(req, sub_file, nPBC_small)
    xyz = []
    for iz, z in enumerate(np.linspace(float(z_abs_range[0]), float(z_abs_range[1]), int(folded_nz), endpoint=True, dtype=np.float32)):
        for iy in range(int(folded_nxy)):
            v = (iy + 0.5) / float(folded_nxy) - 0.5
            for ix in range(int(folded_nxy)):
                u = (ix + 0.5) / float(folded_nxy) - 0.5
                p = u*md.surface_lvec[0, :3].astype(np.float64) + v*md.surface_lvec[1, :3].astype(np.float64)
                xyz.append((p[0], p[1], float(z)))
    xyz = np.array(xyz, dtype=np.float32)
    tref = apply_macro(ref_scan, pack_transforms(xyz), eval_md_components(md, pack_transforms(xyz), min(len(xyz), md.nSystems)))
    fit_mask = build_fit_mask_and_weights(xyz[:, 2] - float(z_abs_range[0]), z_fit_min, tref['total'], repel_cut)
    md.fit_folded_surface_basis(sub_file, nPBC=nPBC_small, z_range=z_abs_range, nu=folded_nu, nv=folded_nv, nz=folded_nzbasis, nxy=folded_nxy, nz_samp=folded_nz, r_damp=0.0, alpha_morse=1.8, bMacro=False, components=COMPONENTS, fit_mask=fit_mask, weight_power=weight_power)
    md.folded_fit_info['fit_ref_total'] = np.asarray(tref['total'], dtype=np.float64)
    md.folded_fit_info['repel_cut'] = float(repel_cut)
    md.folded_fit_info['z_fit_min'] = float(z_fit_min)
    return md

def make_weight_plot(path_png, md, title):
    uvz = md.folded_fit_info['uvz']; ww = np.asarray(md.folded_fit_info['weights'][0], dtype=np.float64); mask = np.asarray(md.folded_fit_info['fit_mask'], dtype=bool); E = np.asarray(md.folded_fit_info.get('fit_ref_total', np.zeros(len(uvz))), dtype=np.float64)
    z = uvz[:, 2]
    zuniq = np.unique(np.round(z, 6))
    Ez = np.zeros(len(zuniq)); wz = np.zeros(len(zuniq)); mz = np.zeros(len(zuniq))
    for i, zz in enumerate(zuniq):
        m = np.abs(z - zz) < 1e-6
        Ez[i] = np.min(E[m])
        wz[i] = np.max(ww[m])
        mz[i] = np.any(mask[m]).astype(float)
    fig, ax = plt.subplots(1, 1, figsize=(8, 4.5))
    ax.plot(zuniq, Ez, '-', lw=1.8, c='k', label='Reference total (macro-corrected)')
    ax.axhline(md.folded_fit_info.get('repel_cut', 0.0), ls='--', lw=1.0, c='tab:red', label='Repulsion cutoff')
    ax.axvline(float(zuniq.min()) + md.folded_fit_info.get('z_fit_min', 0.0), ls='--', lw=1.0, c='tab:green', label='z fit min')
    ax.set_xlabel('absolute z [Å]'); ax.set_ylabel('Energy [eV]'); ax.grid(True, alpha=0.25)
    ax2 = ax.twinx()
    ax2.plot(zuniq, wz, '-', lw=1.4, c='tab:purple', label='Fit weight')
    ax2.plot(zuniq, mz, '--', lw=1.0, c='tab:blue', label='Fit mask')
    ax2.set_ylabel('Weight / mask')
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1+h2, l1+l2, loc='best', frameon=False, fontsize=8)
    ax.set_title(title)
    fig.tight_layout(); fig.savefig(path_png, dpi=180); plt.close(fig)

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

def make_record(kind, probe_name, q, site_name, ref, emb, fold, timing_ref, timing_emb, timing_fold, files, focus_mask=None):
    rec = {'kind': kind, 'probe': probe_name, 'charge': float(q), 'site': site_name, 'files': files, 'timing_s': {'reference': timing_ref, 'embedding': timing_emb, 'folded': timing_fold}, 'components': {}}
    for ck in ('pauli', 'london', 'coulomb', 'total'):
        rec['components'][ck] = {'reference': arr_stats(ref[ck]), 'embedding': arr_stats(emb[ck]), 'folded': arr_stats(fold[ck]), 'embedding_error': err_stats(ref[ck], emb[ck]), 'folded_error': err_stats(ref[ck], fold[ck])}
        if focus_mask is not None:
            m = np.asarray(focus_mask, dtype=bool)
            rec['components'][ck]['embedding_error_focus'] = err_stats(np.asarray(ref[ck])[m], np.asarray(emb[ck])[m])
            rec['components'][ck]['folded_error_focus'] = err_stats(np.asarray(ref[ck])[m], np.asarray(fold[ck])[m])
    return rec

def save_summary_markdown(path_md, records, config):
    lines = ['# Rigorous periodic flat-surface folded-component test summary', '', '## Configuration', '']
    for k, v in config.items(): lines.append(f'- **{k}**: `{v}`')
    lines += ['', '## Main observations', '']
    total_max_fold = max(r['components']['total']['folded_error']['max_abs'] for r in records)
    total_max_emb = max(r['components']['total']['embedding_error']['max_abs'] for r in records)
    lines.append(f'- **Worst total folded max|ΔE|**: `{total_max_fold:.6e} eV`')
    lines.append(f'- **Worst total embedding max|ΔE|**: `{total_max_emb:.6e} eV`')
    for r in records:
        lines += ['', f"## {r['kind']} : {r['probe']} q={r['charge']:+.3f} site={r['site']}", '']
        lines.append(f"- **PNG**: `{r['files']['png']}`"); lines.append(f"- **CSV**: `{r['files']['csv']}`"); lines.append(f"- **JSON**: `{r['files']['json']}`")
        for ck in ('pauli', 'london', 'coulomb', 'total'):
            ce = r['components'][ck]['embedding_error']; cf = r['components'][ck]['folded_error']
            lines.append(f"- **{ck}**: emb max|Δ|=`{ce['max_abs']:.6e}` rmse=`{ce['rmse']:.6e}` ; fold max|Δ|=`{cf['max_abs']:.6e}` rmse=`{cf['rmse']:.6e}`")
            cef = r['components'][ck].get('embedding_error_focus', None); cff = r['components'][ck].get('folded_error_focus', None)
            if (cef is not None) and (cff is not None):
                lines.append(f"  - focus: emb max|Δ|=`{cef['max_abs']:.6e}` rmse=`{cef['rmse']:.6e}` ; fold max|Δ|=`{cff['max_abs']:.6e}` rmse=`{cff['rmse']:.6e}`")
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
    return ap.parse_args()

def main():
    args = parse_args()
    ensure_dir(OUT_DIR)
    info = surface_info(SUB_FILE); atom_types = atom_types_data(); z_heights = np.linspace(args.zmin, args.zmax, args.nz, dtype=np.float32)
    z_abs_range = (float(info['z_top'] + z_heights.min()), float(info['z_top'] + z_heights.max())); nPBC_ref = (12, 12, 0); nPBC_small = (4, 4, 0)
    probes = [('H', 0.2), ('O', 0.2)]; records = []
    config = {'substrate': os.path.relpath(SUB_FILE, ROOT), 'periodic_reference': 'wrapped PBC + macro', 'z_scan_above_top_A': f'{float(z_heights.min()):.3f}..{float(z_heights.max()):.3f}', 'lateral_scan_height_A': float(args.lateral_z), 'nz': len(z_heights), 'nPBC_reference': nPBC_ref, 'nPBC_embedding_local': nPBC_small, 'folded_components': COMPONENTS, 'folded_basis': 'nu=2 nv=2 nzbasis=12 nxy=24 nzsamp=96', 'zfit_min_A': float(args.zfit_min), 'repel_cut_eV': float(args.repel_cut), 'weight_k': float(args.weight_k)}
    for probe_name, q in probes:
        req = make_probe_REQ(probe_name, q, atom_types)
        md_ref = setup_md(req, len(z_heights), SUB_FILE, nPBC_ref, bMacro=False)
        md_emb = setup_md(req, len(z_heights), SUB_FILE, nPBC_small, bMacro=False)
        md_fold = build_folded(req, SUB_FILE, len(z_heights), z_abs_range, nPBC_small, z_fit_min=args.zfit_min, repel_cut=args.repel_cut, weight_power=args.weight_k)
        macro_ref = setup_macro_scanner(req, SUB_FILE, nPBC_ref)
        macro_emb = setup_macro_scanner(req, SUB_FILE, nPBC_small)
        probe_dir = os.path.join(OUT_DIR, f'{probe_name}_q{q:+.2f}'.replace('+', 'p').replace('-', 'm')); ensure_dir(probe_dir)
        make_weight_plot(os.path.join(probe_dir, f'{probe_name}_weights.png'), md_fold, f'{probe_name} folded fit weights')
        for site_name, xy in info['sites'].items():
            transforms, _ = make_z_transforms(xy, info['z_top'], z_heights)
            ref = apply_macro(macro_ref, transforms, eval_md_components(md_ref, transforms, len(z_heights)))
            emb = apply_macro(macro_emb, transforms, eval_md_components(md_emb, transforms, len(z_heights)))
            fold = apply_macro(macro_emb, transforms, eval_folded_components(md_fold, transforms, len(z_heights)))
            prefix = f'{probe_name}_q{q:+.2f}_{site_name}_zscan'.replace('+', 'p').replace('-', 'm')
            print_case_stats(prefix, ref, emb, fold)
            png_file = os.path.join(probe_dir, f'{prefix}.png'); csv_file = os.path.join(probe_dir, f'{prefix}.csv'); json_file = os.path.join(probe_dir, f'{prefix}.json')
            make_scan_plot(png_file, f'{probe_name} q={q:+.2f} periodic z-scan @ {site_name}', z_heights, 'z above top layer [Å]', ref, emb, fold, zoom=(-0.2, 0.2))
            save_case_numeric(csv_file, z_heights, 'z', ref, emb, fold)
            focus_mask = np.asarray(z_heights >= args.zfit_min, dtype=bool)
            rec = make_record('z_scan', probe_name, q, site_name, ref, emb, fold, ref['macro_delta'].size*0.0, emb['macro_delta'].size*0.0, fold['macro_delta'].size*0.0, {'png': png_file, 'csv': csv_file, 'json': json_file}, focus_mask=focus_mask)
            open(json_file, 'w').write(json.dumps(rec, indent=2)); records.append(rec)
        for axis in ('a', 'b'):
            lateral_trans, ts = make_lateral_line(info, info['z_top'] + args.lateral_z, nxy=args.nxy_line, ncell_span=3.0, axis=axis)
            md_refL = setup_md(req, len(ts), SUB_FILE, nPBC_ref, bMacro=False)
            md_embL = setup_md(req, len(ts), SUB_FILE, nPBC_small, bMacro=False)
            md_foldL = build_folded(req, SUB_FILE, len(ts), z_abs_range, nPBC_small, z_fit_min=args.zfit_min, repel_cut=args.repel_cut, weight_power=args.weight_k)
            refL = apply_macro(macro_ref, lateral_trans, eval_md_components(md_refL, lateral_trans, len(ts)))
            embL = apply_macro(macro_emb, lateral_trans, eval_md_components(md_embL, lateral_trans, len(ts)))
            foldL = apply_macro(macro_emb, lateral_trans, eval_folded_components(md_foldL, lateral_trans, len(ts)))
            prefix = f'{probe_name}_q{q:+.2f}_{axis}_lateral'.replace('+', 'p').replace('-', 'm')
            print_case_stats(prefix, refL, embL, foldL)
            png_file = os.path.join(probe_dir, f'{prefix}.png'); csv_file = os.path.join(probe_dir, f'{prefix}.csv'); json_file = os.path.join(probe_dir, f'{prefix}.json')
            make_scan_plot(png_file, f'{probe_name} q={q:+.2f} periodic lateral scan along {axis} z={args.lateral_z:.2f} Å', ts, f'path coordinate along {axis} [primitive cells]', refL, embL, foldL, zoom=(-0.2, 0.2))
            save_case_numeric(csv_file, ts, 's_cell', refL, embL, foldL)
            rec = make_record('lateral_scan', probe_name, q, f'{axis}_periodic', refL, embL, foldL, 0.0, 0.0, 0.0, {'png': png_file, 'csv': csv_file, 'json': json_file}, focus_mask=np.ones(len(ts), dtype=bool))
            open(json_file, 'w').write(json.dumps(rec, indent=2)); records.append(rec)
    summary_json = os.path.join(OUT_DIR, 'summary.json'); summary_md = os.path.join(OUT_DIR, 'summary.md')
    open(summary_json, 'w').write(json.dumps({'config': config, 'records': records}, indent=2)); save_summary_markdown(summary_md, records, config)
    print(f'[flat-folded-periodic] saved summary {summary_json}'); print(f'[flat-folded-periodic] saved summary {summary_md}')

if __name__ == '__main__': main()
