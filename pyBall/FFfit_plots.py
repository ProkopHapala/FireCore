"""Plotting and visualization utilities for FF fitting.

Separated from FFfit_utils to keep matplotlib imports isolated.
All functions use Agg backend by default (no GUI dependency).
"""

import os, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from pyBall.FFfit_utils import angle_type_key, dihedral_angle

# === Spectrum plotting ===

def broaden(freqs, x, width=20.0):
    """Sum of normalized Gaussians centered at freqs, evaluated on x grid."""
    spec = np.zeros_like(x)
    if len(freqs) == 0:
        return spec
    for f in freqs:
        spec += np.exp(-((x - f)**2) / (2 * width**2))
    return spec

def plot_spectrum(freqs_ref, freqs_model, label, outdir=None, exp_spectrum=None,
                  width=20.0, xmax=2500, noshow=True):
    """Plot reference vs model vibrational spectrum (stick + broadened).

    If exp_spectrum is a (2, n) array or a path to a CSV with two columns
    (frequency, intensity), it is overlaid after normalizing to the model peak.
    """
    x = np.arange(0, xmax + 1, 0.5)
    spec_ref = broaden(freqs_ref, x, width)
    spec_model = broaden(freqs_model, x, width)
    if spec_model.max() > 0:
        spec_model /= spec_model.max()
    if spec_ref.max() > 0:
        spec_ref /= spec_ref.max()
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    fig.suptitle(f'Vibrational spectrum — {label}', fontsize=13)
    ax1.vlines(freqs_ref, 0, 1, colors='#2563eb', linewidth=0.6, alpha=0.75, label='reference (PySCF)')
    ax1.vlines(freqs_model, 0, 0.7, colors='#dc2626', linewidth=0.6, alpha=0.75, label='model (FF)')
    ax1.set_ylim(0, 1.1)
    ax1.set_ylabel('Stick (arb.)')
    ax1.legend(fontsize=9, loc='upper right')
    ax1.grid(True, alpha=0.25, ls=':')
    ax1.tick_params(labelbottom=False)
    ax2.plot(x, spec_ref, color='#1d4ed8', lw=1.4, label=f'Reference (σ={width:.0f} cm⁻¹)')
    ax2.fill_between(x, 0, spec_ref, color='#3b82f6', alpha=0.15)
    ax2.plot(x, spec_model, color='#b91c1c', lw=1.4, label=f'Model (σ={width:.0f} cm⁻¹)')
    ax2.fill_between(x, 0, spec_model, color='#f87171', alpha=0.12)
    exp_x = exp_y = None
    if exp_spectrum is not None:
        try:
            if isinstance(exp_spectrum, str):
                exp = np.loadtxt(exp_spectrum)
            else:
                exp = np.asarray(exp_spectrum)
            exp_x, exp_y = exp[:, 0], exp[:, 1]
            if exp_y.max() > 0:
                exp_y = exp_y / exp_y.max()
            ax2.plot(exp_x, exp_y, 'k-', lw=1.0, alpha=0.7, label='Experiment')
        except Exception as e:
            print(f"  WARNING: could not load experimental spectrum: {e}")
    ax2.set_xlim(0, xmax)
    ax2.set_xlabel('Frequency (cm⁻¹)')
    ax2.set_ylabel('Intensity (normalized)')
    ax2.legend(fontsize=9, loc='upper right')
    ax2.grid(True, alpha=0.25, ls=':')
    plt.tight_layout()
    if outdir:
        os.makedirs(outdir, exist_ok=True)
        out = os.path.join(outdir, f'{label}_spectrum.png')
        fig.savefig(out, dpi=150)
        print(f"  Saved spectrum plot: {out}")
    if not noshow:
        plt.show()
    plt.close(fig)
    return out if outdir else None

def plot_spectra_overlay(system_data, outdir=None, width=20.0, xmax=2500, noshow=True):
    """Overlay broadened reference/model spectra across multiple systems.

    system_data: list of (label, freqs_ref, freqs_model) tuples.
    """
    x = np.arange(0, xmax + 1, 0.5)
    fig, ax = plt.subplots(figsize=(13, 7))
    colors = plt.cm.get_cmap('tab10')
    for i, (label, freqs_ref, freqs_model) in enumerate(system_data):
        spec_ref = broaden(freqs_ref, x, width)
        spec_model = broaden(freqs_model, x, width)
        if spec_ref.max() > 0:
            spec_ref /= spec_ref.max()
        if spec_model.max() > 0:
            spec_model /= spec_model.max()
        c = colors(i % 10)
        ax.plot(x, spec_ref + i, color=c, lw=1.2, alpha=0.6, linestyle='-', label=f'{label} ref')
        ax.plot(x, spec_model + i, color=c, lw=1.2, alpha=0.9, linestyle='--', label=f'{label} model')
    ax.set_xlim(0, xmax)
    ax.set_xlabel('Frequency (cm⁻¹)')
    ax.set_ylabel('Normalized intensity (offset by system)')
    ax.set_title('Overlay of reference vs model vibrational spectra')
    ax.legend(fontsize=8, loc='upper right', ncol=2)
    ax.grid(True, alpha=0.25, ls=':')
    if outdir:
        os.makedirs(outdir, exist_ok=True)
        out = os.path.join(outdir, 'spectrum_overlay.png')
        fig.savefig(out, dpi=150)
        print(f"  Saved overlay plot: {out}")
    if not noshow:
        plt.show()
    plt.close(fig)

def plot_comparison_spectra(case_name, freq_ref, model_freqs, outdir, width=20.0, xmax=2500.0):
    """Overlay normalized DFT and fitted-model vibrational spectra."""
    os.makedirs(outdir, exist_ok=True)
    x = np.arange(0.0, xmax + 0.5, 0.5)
    fig, ax = plt.subplots(figsize=(13, 7))
    curves = [('DFT reference', freq_ref, '#111827', 2.5, '-')]
    styles = [('#0072B2', 1.8, '--'), ('#009E73', 1.8, '-.'), ('#D55E00', 1.7, ':'), ('#CC79A7', 1.8, (0, (5, 2)))]
    for (label, freqs), style in zip(model_freqs.items(), styles): curves.append((label, freqs, *style))
    for label, freqs, color, lw, ls in curves:
        y = broaden(freqs, x, width)
        if np.max(y) > 0.0: y /= np.max(y)
        ax.plot(x, y, color=color, lw=lw, ls=ls, label=label)
    ax.set(xlim=(0, xmax), ylim=(0, 1.08), xlabel='Frequency (cm$^{-1}$)', ylabel='Normalized density', title=f'{case_name}: progressive force-field vibration fit')
    ax.grid(True, alpha=0.2, ls=':')
    ax.legend(frameon=False, fontsize=9, ncol=2)
    out = os.path.join(outdir, f'{case_name}_model_comparison_spectra.png')
    fig.savefig(out, dpi=180, bbox_inches='tight')
    plt.close(fig)
    return out

def plot_hessian_comparison(H_ref, H_model, label, outdir):
    """Plot signed values and log magnitudes of reference, model, and residual Hessians."""
    os.makedirs(outdir, exist_ok=True)
    matrices = [H_ref, H_model, H_model - H_ref]
    names = ['DFT reference', 'FF model', 'model - DFT']
    vmax = max(np.max(np.abs(H_ref)), np.max(np.abs(H_model)))
    floor = max(vmax * 1e-10, np.finfo(float).tiny)
    log_matrices = [np.log10(np.maximum(np.abs(H), floor)) for H in matrices]
    log_min = np.log10(floor)
    log_max = np.log10(vmax)
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), constrained_layout=True)
    for i, (H, name) in enumerate(zip(matrices, names)):
        im = axes[0, i].imshow(H, cmap='RdBu_r', vmin=-vmax, vmax=vmax, interpolation='nearest')
        axes[0, i].set_title(name)
        axes[0, i].set_xlabel('Cartesian coordinate')
        axes[0, i].set_ylabel('Cartesian coordinate')
        fig.colorbar(im, ax=axes[0, i], shrink=0.8, label='eV/A^2')
        imlog = axes[1, i].imshow(log_matrices[i], cmap='magma', vmin=log_min, vmax=log_max, interpolation='nearest')
        axes[1, i].set_title(f'log10 |{name}|')
        axes[1, i].set_xlabel('Cartesian coordinate')
        axes[1, i].set_ylabel('Cartesian coordinate')
        fig.colorbar(imlog, ax=axes[1, i], shrink=0.8, label='log10(eV/A^2)')
    out = os.path.join(outdir, f'{label}_hessian.png')
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"  Saved Hessian comparison: {out}")
    return out

# === Equilibrium distributions ===

def equilibrium_distributions(systems, include_third=True, include_dihedrals=True):
    """Collect relaxed internal-coordinate values by transferable interaction type."""
    values = {'bond': {}, 'angle': {}, '1-4': {}, 'dihedral': {}}
    for sys in systems:
        symbols, positions = sys.get('atom_types', sys['symbols']), sys['positions']
        for i, j, _ in sys['bonds']:
            key = '-'.join(sorted((symbols[i], symbols[j])))
            values['bond'].setdefault(key, []).append(np.linalg.norm(positions[j] - positions[i]))
        for i, j, k, _ in sys['angles']:
            key = '-'.join(angle_type_key(symbols, i, j, k, elements=sys['symbols'], central_only=sys.get('angle_central_only', False)))
            u = positions[i] - positions[j]; u /= np.linalg.norm(u)
            v = positions[k] - positions[j]; v /= np.linalg.norm(v)
            values['angle'].setdefault(key, []).append(np.degrees(np.arccos(np.clip(np.dot(u, v), -1.0, 1.0))))
        if include_third:
            for i, j, _ in sys.get('bonds3', []):
                key = '-'.join(sorted((symbols[i], symbols[j])))
                values['1-4'].setdefault(key, []).append(np.linalg.norm(positions[j] - positions[i]))
        if include_dihedrals:
            for i, j, k, l, _, _ in sys.get('dihedrals', []):
                key = '-'.join((symbols[i], symbols[j], symbols[k], symbols[l]))
                values['dihedral'].setdefault(key, []).append(np.degrees(dihedral_angle(positions[[i, j, k, l]])))
    return values

_SI_SUBTYPE_PREFIXES = ('SiH3', 'SiH2', 'SiH', 'Si')

def _element_family_key(subtype_key):
    """Map a subtype key string (e.g. 'H-SiH2') to element-only family (e.g. 'H-Si')."""
    parts = subtype_key.split('-')
    elem_parts = []
    for p in parts:
        matched = False
        for pref in _SI_SUBTYPE_PREFIXES:
            if p == pref:
                elem_parts.append('Si')
                matched = True
                break
        if not matched:
            elem_parts.append(p)
    return '-'.join(elem_parts)

def _plot_stacked_by_family(values, kind, title, unit, outdir, filename):
    """Plot one subplot per element-family, with subtypes as stacked histograms."""
    raw = values[kind]
    if not raw:
        fig, ax = plt.subplots(1, 1, figsize=(6, 4), constrained_layout=True)
        ax.text(0.5, 0.5, 'not enabled', ha='center', va='center', transform=ax.transAxes)
        ax.set_title(title)
        fig.savefig(os.path.join(outdir, filename), dpi=170)
        plt.close(fig)
        return
    families = {}
    for sk in sorted(raw.keys()):
        fk = _element_family_key(sk)
        families.setdefault(fk, []).append(sk)
    fam_keys = sorted(families.keys())
    ncol = min(3, len(fam_keys))
    nrow = int(np.ceil(len(fam_keys) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.5*ncol, 4.0*nrow), constrained_layout=True, squeeze=False)
    rows = []
    for fi, fk in enumerate(fam_keys):
        ax = axes[fi // ncol, fi % ncol]
        subtypes = families[fk]
        all_vals = [np.asarray(raw[sk]) for sk in subtypes]
        vmin = min(v.min() for v in all_vals); vmax = max(v.max() for v in all_vals)
        span = vmax - vmin
        pad = 0.04*span if span > 1e-10 else max(abs(vmin)*1e-6, 1e-6)
        bins = np.linspace(vmin - pad, vmax + pad, 31)
        plot_vals = []
        labels = []
        for si, sk in enumerate(subtypes):
            v = np.asarray(raw[sk])
            plot_vals.append(v)
            labels.append(f'{sk} (n={v.size})')
            rows.append((kind, sk, v.size, np.mean(v), np.std(v), np.min(v), np.max(v), unit))
        ax.hist(plot_vals, bins=bins, stacked=True, histtype='stepfilled', alpha=0.7, linewidth=0.8, label=labels)
        ax.set_title(fk, fontsize=11)
        ax.set_xlabel(unit)
        ax.set_ylabel('count')
        ax.ticklabel_format(axis='x', style='plain', useOffset=False)
        ax.grid(True, alpha=0.2, ls=':')
        ax.legend(fontsize=7, loc='best')
    for fi in range(len(fam_keys), nrow*ncol):
        axes[fi // ncol, fi % ncol].axis('off')
    fig.suptitle(title, fontsize=13)
    out = os.path.join(outdir, filename)
    fig.savefig(out, dpi=170)
    plt.close(fig)
    return rows, out

def plot_equilibrium_distributions(systems, outdir):
    """Plot relaxed bond, angle, 1-4, and signed torsion distributions."""
    values = equilibrium_distributions(systems)
    os.makedirs(outdir, exist_ok=True)
    all_rows = []
    specs = [
        ('bond',     'Bond length distributions by element family',   'Angstrom', 'equilibrium_bond_distributions.png'),
        ('angle',    'Bond angle distributions by element family',    'degree',   'equilibrium_angle_distributions.png'),
        ('1-4',      '1-4 endpoint distance by element family',       'Angstrom', 'equilibrium_1-4_distributions.png'),
        ('dihedral', 'Signed torsion distributions by element family','degree',   'equilibrium_dihedral_distributions.png'),
    ]
    for kind, title, unit, filename in specs:
        result = _plot_stacked_by_family(values, kind, title, unit, outdir, filename)
        if result is not None:
            all_rows.extend(result[0])
            print(f"  Saved {kind} distributions: {result[1]}")
    table = os.path.join(outdir, 'equilibrium_parameter_distributions.csv')
    with open(table, 'w') as f:
        f.write('coordinate,type,count,mean,std,min,max,unit\n')
        for row in all_rows: f.write(','.join(map(str, row)) + '\n')
    print(f"  Saved equilibrium statistics: {table}")
    return values, all_rows, outdir

# === DFT stiffness distributions ===

def dft_stiffness_distributions(systems):
    """Collect diagonal indicators from the least-norm redundant Wilson projection."""
    from pyBall import FFfit as FFfit_cpp
    from pyBall.FFfit_utils import angle_type_key
    values = {'bond': {}, 'angle': {}}
    for sys in systems:
        symbols = sys.get('atom_types', sys['symbols'])
        elements = sys['symbols']
        positions = sys['positions']
        masses = sys['data']['masses']
        H_ref = sys['H_ref']
        bonds = sys['bonds']
        angles = sys['angles']
        angle_central_only = sys.get('angle_central_only', False)
        B, labels = FFfit_cpp.build_wilson_matrix(positions, bonds, angles)
        coord_scale = FFfit_cpp.dimensionless_wilson_scale(positions, labels)
        F, info = FFfit_cpp.internal_hessian_projection(H_ref, B, masses, coordinate_scale=coord_scale)
        Fdiag = np.diag(F)
        for ib, (i, j, _) in enumerate(bonds):
            key = '-'.join(sorted((symbols[i], symbols[j])))
            r0 = 1.0 / coord_scale[ib]
            values['bond'].setdefault(key, []).append(Fdiag[ib] / (r0 * r0))
        for ia, (i, j, k, _) in enumerate(angles):
            key = '-'.join(FFfit_cpp.angle_type_key(symbols, i, j, k, elements=elements, central_only=angle_central_only)) if hasattr(FFfit_cpp, 'angle_type_key') else '-'.join(angle_type_key(symbols, i, j, k, elements=elements, central_only=angle_central_only))
            values['angle'].setdefault(key, []).append(Fdiag[len(bonds) + ia])
    for kind in values:
        for key in values[kind]:
            values[kind][key] = np.asarray(values[kind][key])
    return values

def plot_dft_stiffness_distributions(systems, outdir):
    """Plot least-norm Wilson diagonal indicators by element family."""
    values = dft_stiffness_distributions(systems)
    os.makedirs(outdir, exist_ok=True)
    all_rows = []
    specs = [
        ('bond',  'DFT projected bond indicator (least-norm Wilson diagonal)',  'eV/Å²',  'dft_bond_stiffness_distributions.png'),
        ('angle', 'DFT projected angle indicator (least-norm Wilson diagonal)', 'eV/rad²', 'dft_angle_stiffness_distributions.png'),
    ]
    for kind, title, unit, filename in specs:
        result = _plot_stacked_by_family(values, kind, title, unit, outdir, filename)
        if result is not None:
            all_rows.extend(result[0])
            print(f"  Saved {kind} stiffness distributions: {result[1]}")
    table = os.path.join(outdir, 'dft_stiffness_distributions.csv')
    with open(table, 'w') as f:
        f.write('coordinate,type,count,mean,std,min,max,unit\n')
        for row in all_rows: f.write(','.join(map(str, row)) + '\n')
    print(f"  Saved DFT stiffness statistics: {table}")
    return values, all_rows, outdir

# === Stiffness visualization: clustering + interactive HTML ===

def cluster_1d(values, gap_threshold=2.0):
    """Identify disjunct clusters in 1D data via gap-based splitting."""
    values = np.asarray(values, dtype=float)
    order = np.argsort(values)
    v = values[order]
    if len(v) < 2:
        return [(float(v.mean()), float(v.std()) if len(v) > 1 else 0.0, order)]
    gaps = np.diff(v)
    median_gap = np.median(gaps)
    if median_gap < 1e-12:
        return [(float(v.mean()), float(v.std()), order)]
    is_gap = gaps > gap_threshold * median_gap
    split_idx = np.where(is_gap)[0] + 1
    clusters = []
    prev = 0
    for si in split_idx:
        idx = np.arange(prev, si)
        clusters.append((float(v[idx].mean()), float(v[idx].std()), order[idx]))
        prev = si
    idx = np.arange(prev, len(v))
    clusters.append((float(v[idx].mean()), float(v[idx].std()), order[idx]))
    return clusters

def prepare_stiffness_viz_data(systems):
    """Compute least-norm Wilson diagonal indicators and prepare their visualization."""
    from pyBall import FFfit as FFfit_cpp
    all_bond_stiff = {}
    all_angle_stiff = {}
    sys_data = []
    for si, sys in enumerate(systems):
        symbols = sys.get('atom_types', sys['symbols'])
        elements = sys['symbols']
        positions = sys['positions']
        masses = sys['data']['masses']
        H_ref = sys['H_ref']
        bonds = sys['bonds']
        angles = sys['angles']
        angle_central_only = sys.get('angle_central_only', False)
        B, labels = FFfit_cpp.build_wilson_matrix(positions, bonds, angles)
        coord_scale = FFfit_cpp.dimensionless_wilson_scale(positions, labels)
        F, info = FFfit_cpp.internal_hessian_projection(H_ref, B, masses, coordinate_scale=coord_scale)
        Fdiag = np.diag(F)
        bond_records = []
        for ib, (i, j, _) in enumerate(bonds):
            r0 = 1.0 / coord_scale[ib]
            k_val = Fdiag[ib] / (r0 * r0)
            sk = '-'.join(sorted((symbols[i], symbols[j])))
            fk = _element_family_key(sk)
            bond_records.append({'i': i, 'j': j, 'r0': r0, 'stiff': k_val, 'subtype': sk, 'family': fk})
            all_bond_stiff.setdefault(sk, []).append((si, ib, k_val))
        angle_records = []
        for ia, (i, j, k, _) in enumerate(angles):
            k_val = Fdiag[len(bonds) + ia]
            u = positions[i] - positions[j]; u /= np.linalg.norm(u)
            v = positions[k] - positions[j]; v /= np.linalg.norm(v)
            theta = np.degrees(np.arccos(np.clip(np.dot(u, v), -1.0, 1.0)))
            ak = angle_type_key(symbols, i, j, k, elements=elements, central_only=angle_central_only)
            sk = '-'.join(ak)
            fk = _element_family_key(sk)
            angle_records.append({'i': i, 'j': j, 'k': k, 'theta': theta, 'stiff': k_val, 'subtype': sk, 'family': fk})
            all_angle_stiff.setdefault(sk, []).append((si, ia, k_val))
        sys_data.append({
            'name': sys['name'], 'positions': positions, 'symbols': sys['symbols'],
            'atom_types': symbols, 'bonds': bond_records, 'angles': angle_records,
        })
    bond_clusters = {}
    for sk, items in all_bond_stiff.items():
        vals = np.array([x[2] for x in items])
        clusters = cluster_1d(vals)
        for ci, (center, std, indices) in enumerate(clusters):
            for idx in indices:
                si, bi, _ = items[idx]
                bond_clusters[(si, bi)] = ci
    angle_clusters = {}
    for sk, items in all_angle_stiff.items():
        vals = np.array([x[2] for x in items])
        clusters = cluster_1d(vals)
        for ci, (center, std, indices) in enumerate(clusters):
            for idx in indices:
                si, ai, _ = items[idx]
                angle_clusters[(si, ai)] = ci
    bond_family_vals = {}
    angle_family_vals = {}
    for si, sd in enumerate(sys_data):
        for bi, br in enumerate(sd['bonds']):
            br['cluster'] = bond_clusters.get((si, bi), 0)
            bond_family_vals.setdefault(br['family'], []).append(br['stiff'])
        for ai, ar in enumerate(sd['angles']):
            ar['cluster'] = angle_clusters.get((si, ai), 0)
            angle_family_vals.setdefault(ar['family'], []).append(ar['stiff'])
    family_ranges = {'bond': {}, 'angle': {}}
    for fk, vals in bond_family_vals.items():
        family_ranges['bond'][fk] = (float(np.min(vals)), float(np.max(vals)))
    for fk, vals in angle_family_vals.items():
        family_ranges['angle'][fk] = (float(np.min(vals)), float(np.max(vals)))
    for sd in sys_data:
        sd['family_ranges'] = family_ranges
    return sys_data

def _json_safe(obj):
    """Recursively convert numpy types to native Python for JSON serialization."""
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, dict):
        return {k: _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    return obj

def build_stiffness_html(sys_data, outdir):
    """Generate a self-contained interactive p5.js HTML for a single nanocrystal.

    Features:
      - 3D ball-and-stick rendering with mouse rotation
      - Bonds colored by stiffness (rainbow, vmin/vmax from family)
      - Angles shown as arcs colored by stiffness
      - Dropdown to select: all bonds, specific family, specific cluster
      - Toggle bonds/angles visibility
    """
    name = sys_data['name']
    positions = sys_data['positions']
    symbols = sys_data['symbols']
    atom_types = sys_data['atom_types']
    bonds = sys_data['bonds']
    angles = sys_data['angles']
    family_ranges = sys_data['family_ranges']
    data = {
        'name': name,
        'positions': positions.tolist(),
        'symbols': symbols,
        'atom_types': atom_types,
        'bonds': [{'i': b['i'], 'j': b['j'], 'r0': b['r0'], 'stiff': b['stiff'], 'subtype': b['subtype'], 'family': b['family'], 'cluster': b['cluster']} for b in bonds],
        'angles': [{'i': a['i'], 'j': a['j'], 'k': a['k'], 'theta': a['theta'], 'stiff': a['stiff'], 'subtype': a['subtype'], 'family': a['family'], 'cluster': a['cluster']} for a in angles],
        'family_ranges': family_ranges,
    }
    data_json = json.dumps(_json_safe(data))
    bond_families = sorted(set(b['family'] for b in bonds))
    angle_families = sorted(set(a['family'] for a in angles))
    bond_subtypes = sorted(set(b['subtype'] for b in bonds))
    angle_subtypes = sorted(set(a['subtype'] for a in angles))
    html = f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Stiffness Map — {name}</title>
<script src="https://cdnjs.cloudflare.com/ajax/libs/p5.js/1.9.0/p5.min.js"></script>
<style>
body {{ margin: 0; padding: 0; background: #1a1a2e; color: #eee; font-family: sans-serif; overflow: hidden; }}
#controls {{ position: absolute; top: 10px; left: 10px; z-index: 10; background: rgba(20,20,40,0.85); padding: 12px; border-radius: 8px; max-width: 320px; }}
#controls h3 {{ margin: 0 0 8px 0; font-size: 14px; }}
#controls label {{ font-size: 12px; display: block; margin: 6px 0 2px 0; }}
#controls select, #controls button {{ font-size: 12px; width: 100%; margin: 2px 0; }}
#info {{ position: absolute; bottom: 10px; left: 10px; z-index: 10; background: rgba(20,20,40,0.85); padding: 8px; border-radius: 8px; font-size: 11px; max-width: 400px; }}
canvas {{ display: block; }}
.legend-bar {{ width: 200px; height: 12px; border-radius: 4px; margin: 4px 0; }}
#histPanel {{ position: absolute; top: 10px; right: 10px; z-index: 10; background: rgba(20,20,40,0.9); padding: 10px; border-radius: 8px; width: 340px; }}
#histPanel h3 {{ margin: 0 0 6px 0; font-size: 13px; }}
#histCanvas {{ background: #111; border-radius: 4px; }}
#histPanel label {{ font-size: 11px; display: inline-block; margin: 4px 2px 0 0; }}
#histPanel input[type=range] {{ width: 120px; vertical-align: middle; }}
</style>
</head>
<body>
<div id="controls">
  <h3>{name} — Stiffness Map</h3>
  <label>Show:</label>
  <select id="modeSelect">
    <option value="bonds_all">All bonds</option>
    <option value="angles_all">All angles</option>
  </select>
  <label>Filter by family:</label>
  <select id="familySelect"><option value="all">All families</option></select>
  <label>Filter by cluster:</label>
  <select id="clusterSelect"><option value="all">All clusters</option></select>
  <label>Bond color mode:</label>
  <select id="colorMode">
    <option value="stiffness">Stiffness (rainbow)</option>
    <option value="length">Bond length (rainbow)</option>
    <option value="cluster">Cluster (discrete)</option>
    <option value="subtype">Subtype (discrete)</option>
  </select>
  <label>Atom size:</label>
  <input type="range" id="atomSize" min="5" max="30" value="14" style="width:100%">
  <div id="legendContainer"></div>
</div>
<div id="info"></div>
<div id="histPanel">
  <h3>Stiffness Histogram</h3>
  <canvas id="histCanvas" width="320" height="180"></canvas>
  <div id="histControls">
    <label>vmin: <input type="range" id="vminSlider" min="0" max="100" value="0" step="0.1"></label>
    <label>vmax: <input type="range" id="vmaxSlider" min="0" max="100" value="100" step="0.1"></label>
    <span id="rangeLabel" style="font-size:11px;color:#8af;"></span>
  </div>
</div>
<script>
const DATA = {data_json};
let rotX = -0.3, rotY = 0.5, camZoom = 1.0;
let dragging = false, lastMX = 0, lastMY = 0;
let mode = 'bonds_all', familyFilter = 'all', clusterFilter = 'all', colorMode = 'stiffness';
let atomSize = 14;
let center = [0, 0, 0];
let scaleF = 1.0;
let vminUser = null, vmaxUser = null;

function setup() {{
  createCanvas(windowWidth, windowHeight, WEBGL);
  let pos = DATA.positions;
  for (let i = 0; i < pos.length; i++) {{
    center[0] += pos[i][0]; center[1] += pos[i][1]; center[2] += pos[i][2];
  }}
  center = [center[0]/pos.length, center[1]/pos.length, center[2]/pos.length];
  let maxR = 0;
  for (let i = 0; i < pos.length; i++) {{
    let d = dist(pos[i][0], pos[i][1], pos[i][2], center[0], center[1], center[2]);
    if (d > maxR) maxR = d;
  }}
  scaleF = 250 / maxR;
  let famSel = document.getElementById('familySelect');
  let bondFams = {json.dumps(bond_families)};
  let angFams = {json.dumps(angle_families)};
  for (let f of bondFams) {{ let o = document.createElement('option'); o.value = 'bond_' + f; o.text = 'Bond: ' + f; famSel.appendChild(o); }}
  for (let f of angFams) {{ let o = document.createElement('option'); o.value = 'angle_' + f; o.text = 'Angle: ' + f; famSel.appendChild(o); }}
  let modeSel = document.getElementById('modeSelect');
  modeSel.addEventListener('change', e => {{ mode = e.target.value; updateClusterDropdown(); }});
  famSel.addEventListener('change', e => {{ familyFilter = e.target.value; updateClusterDropdown(); }});
  document.getElementById('clusterSelect').addEventListener('change', e => {{ clusterFilter = e.target.value; updateHistSliders(); }});
  document.getElementById('colorMode').addEventListener('change', e => {{ colorMode = e.target.value; updateLegend(); }});
  document.getElementById('atomSize').addEventListener('input', e => {{ atomSize = parseFloat(e.target.value); }});
  let vminSlider = document.getElementById('vminSlider');
  let vmaxSlider = document.getElementById('vmaxSlider');
  vminSlider.addEventListener('input', e => {{
    let v = parseFloat(e.target.value);
    let range = getStiffRange();
    vminUser = range[0] + (v/100) * (range[1] - range[0]);
    updateRangeLabel();
  }});
  vmaxSlider.addEventListener('input', e => {{
    let v = parseFloat(e.target.value);
    let range = getStiffRange();
    vmaxUser = range[0] + (v/100) * (range[1] - range[0]);
    updateRangeLabel();
  }});
  updateLegend();
  updateHistSliders();
}}

function getStiffRange() {{
  let isAngle = mode.startsWith('angle');
  let items = isAngle ? DATA.angles : DATA.bonds;
  if (familyFilter !== 'all') {{
    let prefix = isAngle ? 'angle_' : 'bond_';
    if (familyFilter.startsWith(prefix)) {{
      let fam = familyFilter.substring(prefix.length);
      items = items.filter(it => it.family === fam);
    }} else return [0, 1];
  }}
  if (items.length === 0) return [0, 1];
  let vals = items.map(it => it.stiff);
  return [Math.min(...vals), Math.max(...vals)];
}}

function updateRangeLabel() {{
  let lbl = document.getElementById('rangeLabel');
  let isAngle = mode.startsWith('angle');
  let unit = isAngle ? 'eV/rad²' : 'eV/Å²';
  if (vminUser !== null && vmaxUser !== null)
    lbl.textContent = ' [' + vminUser.toFixed(3) + ' — ' + vmaxUser.toFixed(3) + ' ' + unit + ']';
  else
    lbl.textContent = ' [auto]';
}}

function updateHistSliders() {{
  let r = getStiffRange();
  vminUser = r[0]; vmaxUser = r[1];
  document.getElementById('vminSlider').value = 0;
  document.getElementById('vmaxSlider').value = 100;
  updateRangeLabel();
}}

function getFamilyFilteredItems() {{
  let isAngle = mode.startsWith('angle');
  let items = isAngle ? DATA.angles : DATA.bonds;
  if (familyFilter !== 'all') {{
    let prefix = isAngle ? 'angle_' : 'bond_';
    if (familyFilter.startsWith(prefix)) {{
      let fam = familyFilter.substring(prefix.length);
      items = items.filter(it => it.family === fam);
    }} else return [];
  }}
  return items;
}}

function drawHistogram() {{
  let canvas = document.getElementById('histCanvas');
  let ctx = canvas.getContext('2d');
  let W = canvas.width, H = canvas.height;
  ctx.clearRect(0, 0, W, H);
  ctx.fillStyle = '#111';
  ctx.fillRect(0, 0, W, H);
  let isAngle = mode.startsWith('angle');
  let items = getFamilyFilteredItems();
  if (items.length === 0) {{ ctx.fillStyle = '#666'; ctx.font = '12px sans-serif'; ctx.fillText('No data', W/2-30, H/2); return; }}
  let vals = items.map(it => it.stiff);
  let vMin = Math.min(...vals), vMax = Math.max(...vals);
  let pad = (vMax - vMin) * 0.05;
  vMin -= pad; vMax += pad;
  let nBins = 40;
  let binW = (vMax - vMin) / nBins;
  let bins = new Array(nBins).fill(0);
  let binClusters = Array.from({{length: nBins}}, () => new Set());
  for (let it of items) {{
    let bi = Math.min(nBins - 1, Math.max(0, Math.floor((it.stiff - vMin) / binW)));
    bins[bi]++;
    binClusters[bi].add(it.cluster);
  }}
  let maxCount = Math.max(...bins);
  if (maxCount === 0) maxCount = 1;
  let plotW = W - 50, plotH = H - 30, ox = 40, oy = 10;
  ctx.strokeStyle = '#444'; ctx.lineWidth = 1;
  ctx.beginPath(); ctx.moveTo(ox, oy); ctx.lineTo(ox, oy + plotH); ctx.lineTo(ox + plotW, oy + plotH); ctx.stroke();
  for (let i = 0; i < nBins; i++) {{
    if (bins[i] === 0) continue;
    let x = ox + (i / nBins) * plotW;
    let bw = plotW / nBins - 1;
    let bh = (bins[i] / maxCount) * plotH;
    let clusters = [...binClusters[i]].sort((a,b)=>a-b);
    if (clusterFilter !== 'all') {{
      let c = parseInt(clusterFilter);
      if (clusters.includes(c)) {{
        let col = clusterColor(c);
        ctx.fillStyle = 'rgb(' + col.join(',') + ')';
      }} else {{
        ctx.fillStyle = 'rgba(80,80,80,0.3)';
      }}
    }} else {{
      let midVal = vMin + (i + 0.5) * binW;
      let fr = getFamilyRange(items[0].family, isAngle);
      let col = rainbow(midVal, fr[0], fr[1]);
      ctx.fillStyle = 'rgb(' + Math.round(col[0]) + ',' + Math.round(col[1]) + ',' + Math.round(col[2]) + ')';
    }}
    ctx.fillRect(x, oy + plotH - bh, bw, bh);
  }}
  if (vminUser !== null && vmaxUser !== null) {{
    let xvmin = ox + ((vminUser - vMin) / (vMax - vMin)) * plotW;
    let xvmax = ox + ((vmaxUser - vMin) / (vMax - vMin)) * plotW;
    ctx.strokeStyle = '#ff0'; ctx.lineWidth = 1.5;
    ctx.beginPath(); ctx.moveTo(xvmin, oy); ctx.lineTo(xvmin, oy + plotH); ctx.stroke();
    ctx.beginPath(); ctx.moveTo(xvmax, oy); ctx.lineTo(xvmax, oy + plotH); ctx.stroke();
    ctx.fillStyle = 'rgba(255,255,0,0.08)';
    ctx.fillRect(xvmin, oy, xvmax - xvmin, plotH);
  }}
  ctx.fillStyle = '#aaa'; ctx.font = '10px sans-serif';
  ctx.fillText(vMin.toFixed(2), ox, H - 5);
  ctx.fillText(vMax.toFixed(2), ox + plotW - 30, H - 5);
  let unit = isAngle ? 'eV/rad²' : 'eV/Å²';
  ctx.fillText(unit, W - 50, H - 5);
  ctx.fillText('count', 5, oy + 10);
}}

function getFamilyRange(fk, isAngle) {{
  let ranges = isAngle ? DATA.family_ranges.angle : DATA.family_ranges.bond;
  return ranges[fk] ? ranges[fk] : [0, 1];
}}

function mouseWheel(event) {{
  camZoom *= event.delta > 0 ? 1.1 : 0.9;
  camZoom = Math.max(0.1, Math.min(10, camZoom));
  return false;
}}

function mousePressed() {{
  if (mouseButton === LEFT) {{ dragging = true; lastMX = mouseX; lastMY = mouseY; }}
}}
function mouseReleased() {{ dragging = false; }}
function mouseDragged() {{
  if (dragging) {{
    rotY += (mouseX - lastMX) * 0.01;
    rotX += (mouseY - lastMY) * 0.01;
    rotX = Math.max(-Math.PI/2, Math.min(Math.PI/2, rotX));
    lastMX = mouseX; lastMY = mouseY;
  }}
}}

function updateClusterDropdown() {{
  let sel = document.getElementById('clusterSelect');
  sel.innerHTML = '<option value="all">All clusters</option>';
  let items = getFamilyFilteredItems();
  let clusters = [...new Set(items.map(it => it.cluster))].sort((a,b)=>a-b);
  for (let c of clusters) {{
    let count = items.filter(it => it.cluster === c).length;
    let o = document.createElement('option'); o.value = c; o.text = 'Cluster ' + c + ' (' + count + ')'; sel.appendChild(o);
  }}
  clusterFilter = 'all';
  updateHistSliders();
}}

function getFilteredItems() {{
  let isAngle = mode.startsWith('angle');
  let items = isAngle ? DATA.angles : DATA.bonds;
  if (familyFilter !== 'all') {{
    let prefix = isAngle ? 'angle_' : 'bond_';
    if (familyFilter.startsWith(prefix)) {{
      let fam = familyFilter.substring(prefix.length);
      items = items.filter(it => it.family === fam);
    }} else return [];
  }}
  if (clusterFilter !== 'all') {{
    let c = parseInt(clusterFilter);
    items = items.filter(it => it.cluster === c);
  }}
  return items;
}}

function rainbow(t, vmin, vmax) {{
  let u = (t - vmin) / (vmax - vmin + 1e-12);
  u = Math.max(0, Math.min(1, u));
  let r, g, b;
  if (u < 0.25) {{ r = 0; g = 4*u; b = 1; }}
  else if (u < 0.5) {{ r = 0; g = 1; b = 1 - 4*(u-0.25); }}
  else if (u < 0.75) {{ r = 4*(u-0.5); g = 1; b = 0; }}
  else {{ r = 1; g = 1 - 4*(u-0.75); b = 0; }}
  return [r*255, g*255, b*255];
}}

function clusterColor(c) {{
  let colors = [[255,100,100],[100,255,100],[100,100,255],[255,255,100],[255,100,255],[100,255,255],[200,200,100],[200,100,200]];
  return colors[c % colors.length];
}}

function subtypeColor(sk, allSubtypes) {{
  let idx = allSubtypes.indexOf(sk);
  let colors = [[255,100,100],[100,255,100],[100,100,255],[255,255,100],[255,100,255],[100,255,255],[200,200,100],[200,100,200],[150,200,50],[50,200,150]];
  return colors[idx % colors.length];
}}

function updateLegend() {{
  let container = document.getElementById('legendContainer');
  let isAngle = mode.startsWith('angle');
  let items = getFilteredItems();
  if (items.length === 0) {{ container.innerHTML = ''; return; }}
  if (colorMode === 'stiffness' || colorMode === 'length') {{
    let ranges = isAngle ? DATA.family_ranges.angle : DATA.family_ranges.bond;
    let unit = colorMode === 'stiffness' ? (isAngle ? 'eV/rad²' : 'eV/Å²') : 'Å';
    let label = colorMode === 'stiffness' ? 'Stiffness' : 'Bond length';
    let fams = [...new Set(items.map(it => it.family))].sort();
    let html = '';
    for (let fk of fams) {{
      let vmin = ranges[fk] ? ranges[fk][0] : 0;
      let vmax = ranges[fk] ? ranges[fk][1] : 1;
      html += '<label>' + label + ' [' + fk + ']: ' + vmin.toFixed(3) + ' — ' + vmax.toFixed(3) + ' ' + unit + '</label>' +
        '<div class="legend-bar" style="background: linear-gradient(to right, rgb(0,0,255), rgb(0,255,255), rgb(0,255,0), rgb(255,255,0), rgb(255,0,0));"></div>';
    }}
    container.innerHTML = html;
  }} else if (colorMode === 'cluster') {{
    let clusters = [...new Set(items.map(it => it.cluster))].sort((a,b)=>a-b);
    let html = '<label>Clusters:</label>';
    for (let c of clusters) {{
      let col = clusterColor(c);
      html += '<div style="display:inline-block;margin:2px 6px;"><span style="display:inline-block;width:10px;height:10px;background:rgb(' + col.join(',') + ');border-radius:50%;"></span> C' + c + '</div>';
    }}
    container.innerHTML = html;
  }} else {{
    let subs = isAngle ? {json.dumps(angle_subtypes)} : {json.dumps(bond_subtypes)};
    let html = '<label>Subtypes:</label>';
    for (let s of subs) {{
      let col = subtypeColor(s, subs);
      html += '<div style="display:inline-block;margin:2px 6px;"><span style="display:inline-block;width:10px;height:10px;background:rgb(' + col.join(',') + ');border-radius:50%;"></span> ' + s + '</div>';
    }}
    container.innerHTML = html;
  }}
}}

function draw() {{
  background(26, 26, 46);
  let hw = width / (2 * camZoom), hh = height / (2 * camZoom);
  ortho(-hw, hw, -hh, hh, -2000, 2000);
  rotateX(rotX);
  rotateY(rotY);
  ambientLight(80);
  directionalLight(200, 200, 200, 0.5, 0.5, -1);
  let pos = DATA.positions;
  let atomColors = {{'Si': [100,180,255], 'H': [255,255,255], 'SiH': [120,200,255], 'SiH2': [140,220,255], 'SiH3': [160,240,255], 'C': [100,100,100], 'O': [255,100,100], 'N': [100,100,255]}};
  for (let i = 0; i < pos.length; i++) {{
    push();
    let at = DATA.atom_types[i];
    let col = atomColors[at] || atomColors[DATA.symbols[i]] || [200,200,200];
    fill(col[0], col[1], col[2]);
    noStroke();
    translate((pos[i][0]-center[0])*scaleF, (pos[i][1]-center[1])*scaleF, (pos[i][2]-center[2])*scaleF);
    sphere(atomSize * (DATA.symbols[i] === 'Si' ? 1.0 : 0.6));
    pop();
  }}
  let isAngle = mode.startsWith('angle');
  if (!isAngle) {{
    for (let b of DATA.bonds) {{
      let p1 = pos[b.i], p2 = pos[b.j];
      push();
      stroke(60, 60, 70);
      strokeWeight(1.5);
      let x1 = (p1[0]-center[0])*scaleF, y1 = (p1[1]-center[1])*scaleF, z1 = (p1[2]-center[2])*scaleF;
      let x2 = (p2[0]-center[0])*scaleF, y2 = (p2[1]-center[1])*scaleF, z2 = (p2[2]-center[2])*scaleF;
      line(x1, y1, z1, x2, y2, z2);
      pop();
    }}
  }}
  let items = getFilteredItems();
  let ranges = isAngle ? DATA.family_ranges.angle : DATA.family_ranges.bond;
  let allSubtypes = isAngle ? {json.dumps(angle_subtypes)} : {json.dumps(bond_subtypes)};
  for (let it of items) {{
    let col;
    let inRange = true;
    if (vminUser !== null && vmaxUser !== null) {{
      inRange = (it.stiff >= vminUser && it.stiff <= vmaxUser);
    }}
    if (colorMode === 'stiffness') {{
      let fr = ranges[it.family];
      let vmin = fr ? fr[0] : 0, vmax = fr ? fr[1] : 1;
      col = rainbow(it.stiff, vmin, vmax);
    }} else if (colorMode === 'length' && !isAngle) {{
      let vmin = Math.min(...items.map(x => x.r0)); vmax = Math.max(...items.map(x => x.r0));
      col = rainbow(it.r0, vmin, vmax);
    }} else if (colorMode === 'cluster') {{
      col = clusterColor(it.cluster);
    }} else {{
      col = subtypeColor(it.subtype, allSubtypes);
    }}
    if (!isAngle) {{
      let p1 = pos[it.i], p2 = pos[it.j];
      push();
      if (inRange) {{
        stroke(col[0], col[1], col[2]);
        strokeWeight(5);
      }} else {{
        stroke(col[0]*0.3, col[1]*0.3, col[2]*0.3);
        strokeWeight(2);
      }}
      let x1 = (p1[0]-center[0])*scaleF, y1 = (p1[1]-center[1])*scaleF, z1 = (p1[2]-center[2])*scaleF;
      let x2 = (p2[0]-center[0])*scaleF, y2 = (p2[1]-center[1])*scaleF, z2 = (p2[2]-center[2])*scaleF;
      line(x1, y1, z1, x2, y2, z2);
      pop();
    }} else {{
      let p1 = pos[it.i], p2 = pos[it.j], p3 = pos[it.k];
      push();
      if (inRange) {{
        stroke(col[0], col[1], col[2]);
        strokeWeight(5);
      }} else {{
        stroke(col[0]*0.3, col[1]*0.3, col[2]*0.3);
        strokeWeight(2);
      }}
      let x2 = (p2[0]-center[0])*scaleF, y2 = (p2[1]-center[1])*scaleF, z2 = (p2[2]-center[2])*scaleF;
      let x1 = (p1[0]-center[0])*scaleF, y1 = (p1[1]-center[1])*scaleF, z1 = (p1[2]-center[2])*scaleF;
      let x3 = (p3[0]-center[0])*scaleF, y3 = (p3[1]-center[1])*scaleF, z3 = (p3[2]-center[2])*scaleF;
      line(x2, y2, z2, x1, y1, z1);
      line(x2, y2, z2, x3, y3, z3);
      pop();
    }}
  }}
  drawHistogram();
  let info = document.getElementById('info');
  let inCount = items.filter(it => (vminUser === null || (it.stiff >= vminUser && it.stiff <= vmaxUser))).length;
  info.innerHTML = 'System: ' + DATA.name + ' | Atoms: ' + pos.length + ' | Showing: ' + inCount + '/' + items.length + ' ' + (isAngle ? 'angles' : 'bonds') + ' | Drag to rotate, scroll to zoom';
}}

function windowResized() {{
  resizeCanvas(windowWidth, windowHeight);
}}
</script>
</body>
</html>'''
    outpath = os.path.join(outdir, f'stiffness_map_{name}.html')
    with open(outpath, 'w') as f:
        f.write(html)
    return outpath

def generate_stiffness_html(systems, outdir):
    """Generate interactive stiffness HTML maps for each nanocrystal."""
    os.makedirs(outdir, exist_ok=True)
    sys_data_list = prepare_stiffness_viz_data(systems)
    paths = []
    for sd in sys_data_list:
        p = build_stiffness_html(sd, outdir)
        paths.append(p)
        print(f"  Generated {sd['name']}: {p}")
    index_html = '<!DOCTYPE html>\n<html><head><meta charset="utf-8"><title>Stiffness Maps Index</title>\n'
    index_html += '<style>body{background:#1a1a2e;color:#eee;font-family:sans-serif;padding:20px;}a{color:#6cf;text-decoration:none;font-size:16px;display:block;margin:8px 0;}a:hover{text-decoration:underline;}</style>\n'
    index_html += '</head><body><h1>Stiffness Maps</h1>\n'
    for p in paths:
        fname = os.path.basename(p)
        name = fname.replace('stiffness_map_', '').replace('.html', '')
        index_html += f'<a href="{fname}">{name}</a>\n'
    index_html += '</body></html>'
    index_path = os.path.join(outdir, 'index.html')
    with open(index_path, 'w') as f:
        f.write(index_html)
    print(f"  Index: {index_path}")
    return paths
