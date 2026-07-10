#!/usr/bin/env python3
"""Fit bond/angle force-field parameters to PySCF reference Hessians.

Loads PySCF hessian .npy data from results directory, builds bond/angle topology
from relaxed .xyz geometry, and solves linear least-squares for stiffness parameters.

Usage:
    python3 tests/tSiNCs/test_FFfit.py
    python3 tests/tSiNCs/test_FFfit.py --case SiH4
    python3 tests/tSiNCs/test_FFfit.py --case Si_R5p0 --mass-weight
"""

import os, sys, json, argparse
import numpy as np
from collections import deque
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pyBall import FFfit as FFfit_cpp
from pyBall.FFfit_utils import (
    BOHR_TO_ANG, HARTREE_TO_EV, BOHR_TO_ANG_INV, HARTREE_PER_BOHR2_TO_EV_PER_ANG2,
    hessian_ha_bohr_to_ev_ang2,
    angle_type_key, dihedral_type_key, assign_si_environment_types,
    interaction_type_counts, print_interaction_type_counts, parent_parameter_name, subtype_shrinkage_rows, family_mean_prior_rows,
    build_cross_param_maps, compute_cross_sensitivity, add_linear_hessian,
    shortest_path_distances, build_3rd_neighbor_bonds, build_topology, build_dihedrals,
    dihedral_energy_gradient, dihedral_angle, dihedral_hessian,
    compute_dihedral_sensitivity, add_dihedral_hessian,
    assign_types, ParamMap, build_global_param_map, make_param_map_for_system,
    build_sensitivity_matrices, accumulate_normal_equations, fit_hessian, compute_model_hessian,
    compute_gradient_term_blocks, fit_gradient_descent, fit_gradient_descent_multi,
    compute_averaged_equilibrium, rebuild_topology_with_averaged,
    get_reference_modes_and_freqs, get_acoustic_projector, get_frequencies_cm1,
    compare_frequencies, one_to_one_frequency_metrics,
    make_cpp_fitter, fit_gradient_descent_multi_cpp,
    get_sensitivity_action, accumulate_mode_normal_equations, fit_modes_multi,
)
from pyBall.FFfit_plots import (
    broaden, plot_spectrum, plot_spectra_overlay, plot_hessian_comparison,
    plot_comparison_spectra, equilibrium_distributions, plot_equilibrium_distributions,
    dft_stiffness_distributions, plot_dft_stiffness_distributions,
    cluster_1d, prepare_stiffness_viz_data, build_stiffness_html, generate_stiffness_html,
)

RESULTS_DIR = "/home/prokop/SIMULATIONS/jobs_pyscf_vib_OUT_small_nc/results"

def load_pyscf_case(case_dir, geometry_only=False):
    """Load PySCF results for a case.

    If geometry_only=True, skip loading hessian.npy, modes.npy, and frequencies
    (the large/expensive files) — only load relaxed geometry and masses.
    """
    d = {}
    if not geometry_only:
        d['hessian'] = np.load(os.path.join(case_dir, 'hessian.npy'))      # (natoms,3,natoms,3) in Ha/Bohr^2
        d['modes']   = np.load(os.path.join(case_dir, 'modes.npy'))         # (nmodes, natoms, 3)
        d['freqs']   = np.load(os.path.join(case_dir, 'frequencies_cm1.npy'))
    d['masses']  = np.load(os.path.join(case_dir, 'masses.npy'))        # (natoms,)
    # Load relaxed geometry
    xyz_path = os.path.join(case_dir, 'relaxed.xyz')
    with open(xyz_path) as f:
        lines = f.readlines()
    natoms = int(lines[0].strip())
    symbols = []
    positions = np.zeros((natoms, 3))
    for i in range(natoms):
        parts = lines[2 + i].split()
        symbols.append(parts[0])
        positions[i] = [float(parts[1]), float(parts[2]), float(parts[3])]
    d['symbols'] = symbols
    d['positions'] = positions  # in Angstrom
    d['natoms'] = natoms
    # Load status
    with open(os.path.join(case_dir, 'status.json')) as f:
        d['status'] = json.load(f)
    return d


def fit_comparison_variant(base_systems, use_third, use_dihedrals, args, parent_prior=None):
    """Fit one interaction-set variant with local DFT equilibrium coordinates."""
    systems = []
    for base in base_systems:
        sys = dict(base)
        sys['bonds'] = list(base['bonds'])
        sys['angles'] = list(base['angles'])
        sys['bonds3'] = list(base['bonds3']) if use_third else []
        sys['dihedrals'] = list(base['dihedrals']) if use_dihedrals else []
        systems.append(sys)
    all_bonds = [s['bonds'] for s in systems]; all_angles = [s['angles'] for s in systems]
    all_bonds3 = [s['bonds3'] for s in systems]; all_dihedrals = [s['dihedrals'] for s in systems]
    all_symbols = [s.get('atom_types', s['symbols']) for s in systems]
    all_elements = [s['symbols'] for s in systems]
    angle_central_only = any(s.get('angle_central_only', False) for s in systems)
    if angle_central_only and not all(s.get('angle_central_only', False) for s in systems):
        raise ValueError("all systems in one fit must use the same angle typing scope")
    if (args.stretch_stretch or args.stretch_bend) and args.equilibrium != 'local':
        raise ValueError("explicit valence cross terms require --equilibrium local")
    bmap, b3map, amap, dmap, nbase = build_global_param_map(all_bonds, all_angles, all_symbols, all_bonds3, all_dihedrals, all_elements, angle_central_only)
    ssmap, sbmap, npar = build_cross_param_maps(all_angles, all_symbols, nbase, all_elements, angle_central_only, args.stretch_stretch, args.stretch_bend)
    print(f"  parameter types: bond={len(bmap)} 1-4={len(b3map)} angle={len(amap)} torsion={len(dmap)} stretch-stretch={len(ssmap)} stretch-bend={len(sbmap)} total={npar}")
    if args.equilibrium == 'type-average':
        avg_l0, avg_l0_3, avg_t0 = compute_averaged_equilibrium(all_bonds, all_angles, all_symbols, [s['positions'] for s in systems], bmap, amap, b3map, all_bonds3, all_elements, angle_central_only)
        for sys in systems:
            sys['bonds'], sys['angles'], sys['bonds3'] = rebuild_topology_with_averaged(sys['bonds'], sys['angles'], sys['bonds3'], sys.get('atom_types', sys['symbols']), sys['positions'], avg_l0, avg_t0, avg_l0_3, sys['symbols'], angle_central_only)
    for sys in systems:
        sys['dihedral_A'] = compute_dihedral_sensitivity(sys['positions'], sys.get('atom_types', sys['symbols']), sys['dihedrals'], dmap, h=args.dihedral_h) if use_dihedrals else {}
        sys['cross_A'] = compute_cross_sensitivity(sys['positions'], sys['bonds'], sys['angles'], sys.get('atom_types', sys['symbols']), ssmap, sbmap, elements=sys['symbols'], angle_central_only=angle_central_only)
    print("  sensitivity stage: torsions complete")
    fitters = [make_cpp_fitter(sys, (bmap, b3map, amap), npar) for sys in systems]
    hybrid_systems = []
    for f, sys in zip(fitters, systems):
        extra_A = dict(sys['dihedral_A']); extra_A.update(sys['cross_A'])
        hybrid_systems.append({'A': FFfit_cpp.collect_sensitivity_matrices(f, extra=extra_A), 'H_ref': sys['H_ref'],
                               'positions': sys['positions'], 'masses': sys['data']['masses'], 'bonds': sys['bonds'], 'angles': sys['angles']})
    print("  sensitivity stage: stacked model complete")
    names = [''] * npar
    for key, t in bmap.items(): names[t] = 'bond:' + '-'.join(key)
    for key, t in b3map.items(): names[t] = '1-4:' + '-'.join(key)
    for key, t in amap.items(): names[t] = 'angle:' + '-'.join(key)
    for key, t in dmap.items(): names[t] = 'torsion:' + '/'.join(('-'.join(key[0]), '-'.join(key[1])))
    for key, t in ssmap.items(): names[t] = 'stretch-stretch:' + '-'.join(key)
    for key, t in sbmap.items(): names[t] = 'stretch-bend:' + '-'.join(key)
    prior = np.ones(npar)
    for t in bmap.values(): prior[t] = 5.0
    for t in b3map.values(): prior[t] = 0.1
    for t in amap.values(): prior[t] = 1.0
    for t in dmap.values(): prior[t] = 0.1
    for t in list(ssmap.values()) + list(sbmap.values()): prior[t] = 0.0
    scale = np.maximum(np.abs(prior), 0.1)
    for t in list(ssmap.values()) + list(sbmap.values()): scale[t] = args.cross_scale
    coupling_matrix = None
    coupling_target = None
    if parent_prior is not None:
        family_scale = np.empty(nbase)
        for i, name in enumerate(names[:nbase]):
            parent = parent_parameter_name(name)
            if parent not in parent_prior:
                raise KeyError(f"missing elemental parent scale {parent} for subtype parameter {name}")
            family_scale[i] = max(abs(parent_prior[parent]), 0.1)
        R_dev, y_dev, groups = subtype_shrinkage_rows(names[:nbase], args.subtype_shrinkage, family_scale)
        R_mean, y_mean = family_mean_prior_rows(names[:nbase], parent_prior, args.reg)
        R_dev = np.pad(R_dev, ((0, 0), (0, npar - nbase)))
        R_mean = np.pad(R_mean, ((0, 0), (0, npar - nbase)))
        coupling_matrix = np.vstack((R_dev, R_mean))
        coupling_target = np.concatenate((y_dev, y_mean))
        print(f"  subtype hierarchy: deviation_strength={args.subtype_shrinkage:g}, mean_prior_strength={args.reg:g}, rows={len(coupling_target)}, families={sum(len(v) > 1 for v in groups.values())}")
    cross_indices = np.array(list(ssmap.values()) + list(sbmap.values()), dtype=int)
    if cross_indices.size:
        R_cross = np.zeros((cross_indices.size, npar))
        R_cross[np.arange(cross_indices.size), cross_indices] = np.sqrt(args.cross_regularization) / args.cross_scale
        coupling_matrix = R_cross if coupling_matrix is None else np.vstack((coupling_matrix, R_cross))
        y_cross = np.zeros(cross_indices.size)
        coupling_target = y_cross if coupling_target is None else np.concatenate((coupling_target, y_cross))
    k, diag = FFfit_cpp.fit_hybrid_hessian(
        hybrid_systems, mode_weight=args.hybrid_mode, local_weight=args.hybrid_local, internal_weight=args.hybrid_internal,
        mode_balance=args.mode_weight, mode_mixing=args.mode_mixing, frequency_floor_cm1=args.frequency_floor,
        local_graph_distance=args.local_graph_distance, prior=prior, regularization=args.reg if parent_prior is None else 0.0, parameter_scale=scale,
        coupling_matrix=coupling_matrix, coupling_target=coupling_target,
        bounds=(np.where(np.isin(np.arange(npar), cross_indices), -args.cross_bound, 0.0), np.where(np.isin(np.arange(npar), cross_indices), args.cross_bound, np.inf)))
    diag['subtype_shrinkage'] = args.subtype_shrinkage if parent_prior is not None else 0.0
    print(f"  solve complete: condition={diag['condition']:.3g} residual={diag['relative_residual']:.3g}")
    models = []
    for f, sys in zip(fitters, systems):
        f.set_params(k)
        H = f.compute_model_hessian()
        add_linear_hessian(H, k, sys['dihedral_A'])
        add_linear_hessian(H, k, sys['cross_A'])
        models.append(H)
    return systems, models, k, names, diag


def run_model_comparison(base_systems, args):
    """Fit and compare bond/angle, 1-4, torsion, and combined models."""
    variants = [
        ('Bond + angle', False, False),
        ('Bond + angle + 1-4', True, False),
        ('Bond + angle + UFF torsion', False, True),
        ('Bond + angle + UFF torsion + 1-4', True, True),
    ]
    fitted = {}
    for label, use_third, use_dihedrals in variants:
        print(f"\n=== Comparison model: {label} ===")
        fitted[label] = fit_comparison_variant(base_systems, use_third, use_dihedrals, args)
    rows = []
    os.makedirs(args.plot_dir, exist_ok=True)
    for isys, base in enumerate(base_systems):
        freq_ref = get_frequencies_cm1(base['H_ref'], base['data']['masses'], data=base['data'], freq_floor=10.0)
        model_freqs = {}
        rows.append((base['name'], 'DFT reference', 0, 0.0, 0.0, 0.0, 0, 0, 1.0, ''))
        for label, _, _ in variants:
            systems, models, k, names, diag = fitted[label]
            sys, H = systems[isys], models[isys]
            freqs = get_frequencies_cm1(H, sys['data']['masses'], freq_floor=10.0, positions=sys['positions'], project_acoustic=True)
            model_freqs[label] = freqs
            rmse, mae, unmatched = one_to_one_frequency_metrics(freq_ref, freqs)
            lam = FFfit_cpp.reference_vibrational_modes(H, sys['positions'], sys['data']['masses'])[1]
            negative = int(np.sum(lam < -(10.0/521.5)**2))
            rel_frob = 100.0 * np.linalg.norm(H - sys['H_ref']) / np.linalg.norm(sys['H_ref'])
            ptxt = '; '.join(f'{name}={value:.4g}' for name, value in zip(names, k))
            rows.append((base['name'], label, len(k), rmse, mae, rel_frob, unmatched, negative, diag['condition'], ptxt))
        out = plot_comparison_spectra(base['name'], freq_ref, model_freqs, args.plot_dir, args.spectrum_width, args.spectrum_xmax)
        print(f"  Saved model spectrum comparison: {out}")
    header = ('case', 'model', 'npar', 'freq_RMSE_cm-1', 'freq_MAE_cm-1', 'relFrob_percent', 'unmatched_modes', 'negative_modes', 'condition', 'parameters')
    csv_path = os.path.join(args.plot_dir, 'model_comparison.csv')
    md_path = os.path.join(args.plot_dir, 'model_comparison.md')
    with open(csv_path, 'w') as f:
        f.write(','.join(header) + '\n')
        for row in rows: f.write(','.join(f'"{x}"' if isinstance(x, str) else str(x) for x in row) + '\n')
    with open(md_path, 'w') as f:
        f.write('| Case | Model | Npar | RMSE cm^-1 | MAE cm^-1 | relFrob % | Unmatched | Negative | Condition |\n')
        f.write('|---|---|---:|---:|---:|---:|---:|---:|---:|\n')
        for row in rows:
            f.write(f'| {row[0]} | {row[1]} | {row[2]} | {row[3]:.3f} | {row[4]:.3f} | {row[5]:.3f} | {row[6]} | {row[7]} | {row[8]:.3g} |\n')
    print(f"\n{'Model':43s} {'N':>3s} {'RMSE':>9s} {'MAE':>9s} {'Frob%':>8s} {'miss':>5s} {'neg':>4s} {'cond':>8s}")
    for row in rows:
        if row[1] == 'DFT reference': continue
        print(f"{row[1]:43s} {row[2]:3d} {row[3]:9.2f} {row[4]:9.2f} {row[5]:8.2f} {row[6]:5d} {row[7]:4d} {row[8]:8.2f}")
    print(f"  Saved comparison table: {md_path}")
    return rows


def run_typing_comparison(base_systems, args):
    """Compare elemental and H-coordination Si typing using bond+angle models only."""
    elemental = []
    subtyped = []
    for base in base_systems:
        se = dict(base); se['atom_types'] = list(base['symbols']); se['angle_central_only'] = False; elemental.append(se)
        ss = dict(base); ss['atom_types'] = assign_si_environment_types(base['symbols'], base['bonds'], enabled=True); ss['angle_central_only'] = args.si_subtype_scope == 'central'; subtyped.append(ss)
    print("\n=== Elemental typing support ===")
    print_interaction_type_counts(elemental, args.min_type_count)
    print("\n=== Si environment typing support ===")
    print_interaction_type_counts(subtyped, args.min_type_count)
    elemental_fit = fit_comparison_variant(elemental, False, False, args)
    elemental_prior = dict(zip(elemental_fit[3], elemental_fit[2]))
    fits = {
        'Elemental Si typing': elemental_fit,
        'SiH3/SiH2/SiH/Si typing': fit_comparison_variant(subtyped, False, False, args, parent_prior=elemental_prior),
    }
    rows = []
    for isys, base in enumerate(base_systems):
        freq_ref = get_frequencies_cm1(base['H_ref'], base['data']['masses'], data=base['data'], freq_floor=10.0)
        model_freqs = {}
        for label, (systems, models, k, names, diag) in fits.items():
            sys, H = systems[isys], models[isys]
            freqs = get_frequencies_cm1(H, sys['data']['masses'], freq_floor=10.0, positions=sys['positions'], project_acoustic=True)
            model_freqs[label] = freqs
            rmse, mae, unmatched = one_to_one_frequency_metrics(freq_ref, freqs)
            rel_frob = 100.0*np.linalg.norm(H - sys['H_ref'])/np.linalg.norm(sys['H_ref'])
            rows.append((base['name'], label, len(k), rmse, mae, rel_frob, unmatched, diag['condition']))
        out = plot_comparison_spectra(base['name'] + '_typing', freq_ref, model_freqs, args.plot_dir, args.spectrum_width, args.spectrum_xmax)
        print(f"  Saved typing spectrum comparison: {out}")
    csv_path = os.path.join(args.plot_dir, 'typing_comparison.csv')
    md_path = os.path.join(args.plot_dir, 'typing_comparison.md')
    param_path = os.path.join(args.plot_dir, 'typing_parameters.csv')
    with open(csv_path, 'w') as f:
        f.write('case,typing,npar,freq_RMSE_cm-1,freq_MAE_cm-1,relFrob_percent,unmatched_modes,condition\n')
        for row in rows: f.write(','.join(f'"{x}"' if isinstance(x, str) else str(x) for x in row) + '\n')
    with open(md_path, 'w') as f:
        f.write('| Case | Typing | Npar | RMSE cm^-1 | MAE cm^-1 | relFrob % | Unmatched | Condition |\n')
        f.write('|---|---|---:|---:|---:|---:|---:|---:|\n')
        for row in rows: f.write(f'| {row[0]} | {row[1]} | {row[2]} | {row[3]:.3f} | {row[4]:.3f} | {row[5]:.3f} | {row[6]} | {row[7]:.3g} |\n')
        f.write('| **Mean** | **Elemental Si typing** | | {:.3f} | {:.3f} | {:.3f} | | |\n'.format(*[np.mean([r[i] for r in rows if r[1] == 'Elemental Si typing']) for i in (3, 4, 5)]))
        f.write('| **Mean** | **SiH3/SiH2/SiH/Si typing** | | {:.3f} | {:.3f} | {:.3f} | | |\n'.format(*[np.mean([r[i] for r in rows if r[1].startswith('SiH3')]) for i in (3, 4, 5)]))
    with open(param_path, 'w') as f:
        f.write('typing,parameter,value\n')
        for label, (_, _, k, names, _) in fits.items():
            for name, value in zip(names, k): f.write(f'"{label}","{name}",{value}\n')
    print(f"  Saved typing comparison table: {md_path}")
    print(f"  Saved typing parameters: {param_path}")
    return rows

def main():
    parser = argparse.ArgumentParser(description='Fit FF parameters to PySCF Hessians via C++ FFfit library')
    parser.add_argument('--cases', default='SiH4', help='Comma-separated case names or "all_Si" or "all"')
    parser.add_argument('--mass-weight', action='store_true', help='Use mass-weighted Hessian')
    parser.add_argument('--bond-cutoff', type=float, default=2.5, help='Bond distance cutoff (Angstrom)')
    parser.add_argument('--third-bonds', action='store_true', help='Add 3rd-neighbor (1-4) distance springs to mimic dihedral/cross-terms')
    parser.add_argument('--third-bond-cutoff', type=float, default=4.0, help='Distance cutoff for 3rd-neighbor springs (Angstrom)')
    parser.add_argument('--dihedrals', action='store_true', help='Add UFF/Prokop torsion (dihedral) terms to the fit')
    parser.add_argument('--dihedral-n', type=int, default=3, help='Dihedral periodicity n (default 3 for sp3)')
    parser.add_argument('--dihedral-d', type=float, default=1.0, help='Dihedral phase sign d (default 1.0)')
    parser.add_argument('--dihedral-h', type=float, default=1e-5, help='Finite-difference step for dihedral Hessian (Angstrom)')
    parser.add_argument('--equilibrium', choices=['local', 'type-average'], default='local', help='Use each relaxed internal coordinate (default) or one transferable equilibrium value per type')
    parser.add_argument('--plot-distributions', action='store_true', help='Plot relaxed bond/angle/1-4/torsion distributions by type')
    parser.add_argument('--plot-only', action='store_true', help='Only plot equilibrium distributions (skip Hessian loading and fitting)')
    parser.add_argument('--plot-stiffness', action='store_true', help='Plot DFT-extracted bond/angle stiffnesses via Wilson GF projection (requires Hessian)')
    parser.add_argument('--stiffness-html', action='store_true', help='Generate interactive p5.js HTML stiffness maps per nanocrystal (requires Hessian)')
    parser.add_argument('--compare-models', action='store_true', help='Fit the four progressive interaction sets and write a comparison table/spectrum overlay')
    parser.add_argument('--si-subtypes', action='store_true', help='Type Si atoms as SiH3 (including SiH4), SiH2, SiH, or bulk Si by bonded-H count')
    parser.add_argument('--si-subtype-scope', choices=['central', 'full'], default='central', help='Subtype angle centers only (default) or all three angle atoms')
    parser.add_argument('--subtype-shrinkage', type=float, default=0.001, help='Dimensionless within-element family variance penalty for --compare-typing (0 disables subtype offsets)')
    parser.add_argument('--stretch-stretch', action='store_true', help='Optionally fit symmetric dimensionless bond-stretch cross terms per angle type')
    parser.add_argument('--stretch-bend', action='store_true', help='Optionally fit symmetric stretch-bend cross terms per angle type')
    parser.add_argument('--cross-regularization', type=float, default=0.02, help='Zero-centered regularization strength for signed valence cross terms')
    parser.add_argument('--cross-scale', type=float, default=5.0, help='Physical scale in eV for dimensionless valence cross coefficients')
    parser.add_argument('--cross-bound', type=float, default=20.0, help='Absolute bound in eV for each signed valence cross coefficient')
    parser.add_argument('--compare-typing', action='store_true', help='Compare elemental and Si H-coordination typing using bond+angle models')
    parser.add_argument('--min-type-count', type=int, default=3, help='Flag interaction types represented by fewer observations')
    parser.add_argument('--mode-weight', type=str, default='frequency', choices=['equal', 'frequency', 'relative'], help='All-mode balance: equal Hessian accuracy, frequency-error balance (default), or relative eigenvalue accuracy')
    parser.add_argument('--hybrid-mode', type=float, default=1.0, help='Weight of the complete vibrational-mode objective')
    parser.add_argument('--hybrid-local', type=float, default=1.0, help='Weight of graph-local mass-weighted Hessian blocks')
    parser.add_argument('--hybrid-internal', type=float, default=1.0, help='Weight of the orthonormal Wilson row-space curvature objective')
    parser.add_argument('--mode-mixing', type=float, default=0.1, help='Relative penalty for off-diagonal mixing of reference modes')
    parser.add_argument('--local-graph-distance', type=int, default=2, help='Maximum bond-graph distance included in local Hessian blocks')
    parser.add_argument('--frequency-floor', type=float, default=50.0, help='Frequency floor (cm^-1) for stable mode weighting')
    parser.add_argument('--reg', type=float, default=2e-3, help='Dimensionless per-parameter prior regularization weight')
    parser.add_argument('--use-nnls', action='store_true', help='Deprecated compatibility flag; hybrid fitting is non-negative by default')
    parser.add_argument('--allow-negative', action='store_true', help='Allow unphysical negative stiffnesses for diagnostic comparison only')
    parser.add_argument('--plot', action='store_true', help='Plot reference vs model vibrational spectra')
    parser.add_argument('--plot-dir', type=str, default='tests/tSiNCs/OUT_FFfit_plots', help='Output directory for spectrum plots')
    parser.add_argument('--exp-spectrum', type=str, default=None, help='Optional experimental IR/FTIR spectrum file (two columns: frequency cm^-1, intensity)')
    parser.add_argument('--spectrum-width', type=float, default=20.0, help='Gaussian broadening width (cm^-1)')
    parser.add_argument('--spectrum-xmax', type=float, default=2500.0, help='Maximum frequency in spectrum plots (cm^-1)')
    args = parser.parse_args()
    if args.stretch_stretch or args.stretch_bend:
        if not (args.compare_typing or args.compare_models):
            parser.error("--stretch-stretch/--stretch-bend are currently supported by --compare-typing or --compare-models")
        if args.cross_regularization < 0.0 or args.cross_scale <= 0.0 or args.cross_bound <= 0.0:
            parser.error("cross regularization must be non-negative; cross scale and bound must be positive")
    
    # Determine case list
    if args.cases == 'all_Si':
        case_names = ['SiH4', 'Si_R3p8', 'Si_R4p5', 'Si_R5p0', 'Si_R5p5', 'Si_R6p0']
    elif args.cases == 'all':
        case_names = ['SiH4', 'Si_R3p8', 'Si_R4p5', 'Si_R5p0', 'Si_R5p5', 'Si_R6p0',
                      'C_diamond_R2p8', 'C_diamond_R3p2', 'C_diamond_R3p8', 'C_diamond_R4p2']
    else:
        case_names = [c.strip() for c in args.cases.split(',')]
    
    # === Load all systems ===
    geo_only = args.plot_only
    print(f"=== Loading {len(case_names)} systems: {case_names} ===" + (' (geometry only)' if geo_only else ''))
    systems = []
    for name in case_names:
        case_dir = os.path.join(RESULTS_DIR, name)
        if not os.path.isdir(case_dir):
            print(f"  WARNING: skipping {name} (directory not found)")
            continue
        data = load_pyscf_case(case_dir, geometry_only=geo_only)
        natoms = data['natoms']
        if geo_only:
            H_ref = np.zeros((natoms*3, natoms*3))  # placeholder, not used
            H_ref_w = H_ref
            weight = None
        else:
            H_ref = hessian_ha_bohr_to_ev_ang2(
                data['hessian'].transpose(0, 2, 1, 3).reshape(natoms*3, natoms*3))
            if args.mass_weight:
                masses = data['masses']
                inv_sqrt_m = 1.0 / np.sqrt(np.repeat(masses, 3))
                weight = (inv_sqrt_m[:, None] * inv_sqrt_m[None, :]).flatten()
                H_ref_w = (inv_sqrt_m[:, None] * inv_sqrt_m[None, :]) * H_ref
            else:
                weight = None
                H_ref_w = H_ref
        bonds, angles, bonds3 = build_topology(data['symbols'], data['positions'], args.bond_cutoff,
                                                third_bonds=args.third_bonds or args.compare_models, third_bond_cutoff=args.third_bond_cutoff)
        atom_types = assign_si_environment_types(data['symbols'], bonds, enabled=args.si_subtypes or args.compare_typing)
        dihedrals = build_dihedrals(data['symbols'], data['positions'], bonds,
                                    d=args.dihedral_d, n=args.dihedral_n, dihedral=args.dihedrals or args.compare_models)
        sys = {
            'name': name, 'data': data, 'natoms': natoms, 'H_ref': H_ref, 'H_ref_w': H_ref_w,
            'bonds': bonds, 'angles': angles, 'bonds3': bonds3, 'dihedrals': dihedrals, 'weight': weight,
            'symbols': data['symbols'], 'atom_types': atom_types, 'angle_central_only': (args.si_subtypes or args.compare_typing) and args.si_subtype_scope == 'central', 'positions': data['positions'],
        }
        systems.append(sys)
        print(f"  {name}: natoms={natoms}, bonds={len(bonds)}, angles={len(angles)}, 3rd_bonds={len(bonds3)}, dihedrals={len(dihedrals)}, "
              f"H_ref [{H_ref.min():.2e}, {H_ref.max():.2e}]")
        if args.si_subtypes or args.compare_typing:
            unique, count = np.unique(atom_types, return_counts=True)
            print(f"    Si environment populations: {dict(zip(unique.tolist(), count.tolist()))}")
    
    if not systems:
        print("ERROR: no valid systems found"); sys.exit(1)

    print_interaction_type_counts(systems, min_count=args.min_type_count)

    if args.plot_only:
        plot_equilibrium_distributions(systems, args.plot_dir)
        return

    if args.plot_stiffness:
        plot_dft_stiffness_distributions(systems, args.plot_dir)
        return

    if args.stiffness_html:
        generate_stiffness_html(systems, args.plot_dir)
        return

    if args.compare_typing:
        plot_equilibrium_distributions(systems, args.plot_dir)
        run_typing_comparison(systems, args)
        return

    if args.compare_models:
        plot_equilibrium_distributions(systems, args.plot_dir)
        run_model_comparison(systems, args)
        return
    
    # === Build global parameter mapping ===
    all_bonds = [s['bonds'] for s in systems]
    all_angles = [s['angles'] for s in systems]
    all_bonds3 = [s['bonds3'] for s in systems]
    all_dihedrals = [s['dihedrals'] for s in systems]
    all_symbols = [s.get('atom_types', s['symbols']) for s in systems]
    all_elements = [s['symbols'] for s in systems]
    all_positions = [s['positions'] for s in systems]
    global_bond_map, global_bond3_map, global_angle_map, global_dihedral_map, n_total = build_global_param_map(
        all_bonds, all_angles, all_symbols, all_bonds3=all_bonds3, all_dihedrals=all_dihedrals,
        all_elements=all_elements, angle_central_only=args.si_subtypes and args.si_subtype_scope == 'central')
    n_cpp = len(global_bond_map) + len(global_bond3_map) + len(global_angle_map)
    print(f"\n=== Global Parameter Mapping ===")
    print(f"  bond types ({len(global_bond_map)}): {global_bond_map}")
    if global_bond3_map:
        print(f"  3rd-neighbor bond types ({len(global_bond3_map)}): {global_bond3_map}")
    print(f"  angle types ({len(global_angle_map)}): {global_angle_map}")
    if global_dihedral_map:
        print(f"  dihedral types ({len(global_dihedral_map)}): {global_dihedral_map}")
    print(f"  total free params: {n_total} (C++ params: {n_cpp}, dihedral params: {len(global_dihedral_map)})")
    
    # === Compute averaged equilibrium l0/theta0 across all systems ===
    print(f"\n=== Averaged Equilibrium Parameters (multi-system) ===")
    avg_l0, avg_l0_3, avg_theta0 = compute_averaged_equilibrium(
        all_bonds, all_angles, all_symbols, all_positions,
        global_bond_map, global_angle_map, global_bond3_map=global_bond3_map, all_bonds3=all_bonds3,
        all_elements=all_elements, angle_central_only=args.si_subtypes and args.si_subtype_scope == 'central')
    
    if args.plot_distributions or args.plot:
        plot_equilibrium_distributions(systems, args.plot_dir)

    # Local Hessian fitting shares stiffnesses but expands each interaction at its
    # actual relaxed DFT coordinate. Type-averaged minima are a separate energy-model choice.
    if args.equilibrium == 'type-average':
        for sys in systems:
            sys['bonds'], sys['angles'], sys['bonds3'] = rebuild_topology_with_averaged(
                sys['bonds'], sys['angles'], sys['bonds3'], sys.get('atom_types', sys['symbols']), sys['positions'],
                avg_l0, avg_theta0, avg_l0_3=avg_l0_3 if args.third_bonds else None, elements=sys['symbols'],
                angle_central_only=args.si_subtypes and args.si_subtype_scope == 'central')

    # === Precompute per-system dihedral sensitivity matrices (Python side) ===
    if args.dihedrals:
        print(f"\n=== Precomputing Dihedral Sensitivities ===")
        for sys in systems:
            sys['dihedral_A'] = compute_dihedral_sensitivity(
                sys['positions'], sys.get('atom_types', sys['symbols']), sys['dihedrals'], global_dihedral_map, h=args.dihedral_h)
            print(f"  {sys['name']}: computed {len(sys['dihedral_A'])} dihedral sensitivity matrices")
    else:
        for sys in systems:
            sys['dihedral_A'] = {}

    # === Create C++ FFfit instances with global param indices ===
    fitters = [make_cpp_fitter(sys, (global_bond_map, global_bond3_map, global_angle_map), n_total) for sys in systems]
    
    # C++ fitting knows only bond/angle/3rd-bond parameters; switch to n_cpp for Methods 1 and 2
    for f in fitters: f.set_n_free(n_cpp)

    # === Method 1: Multi-system linear least-squares (C++) ===
    print(f"\n=== Method 1: Multi-System Linear Least-Squares (C++) ===")
    G = np.zeros((n_cpp, n_cpp))
    y = np.zeros(n_cpp)
    for f, sys in zip(fitters, systems):
        f.accumulate_normal_equations(G, y, sys['H_ref_w'].ravel(), weight=sys.get('weight'))
    k_lsq = FFfit_cpp.FFfit.solve_normal_equations(G, y)
    for key, t in global_bond_map.items():
        print(f"  k_bond[{key}] = {k_lsq[t]:.6e} eV/A^2")
    for key, t in global_bond3_map.items():
        print(f"  k_bond3[{key}] = {k_lsq[t]:.6e} eV/A^2")
    for key, t in global_angle_map.items():
        print(f"  k_angle[{key}] = {k_lsq[t]:.6e} eV/rad^2")

    # === Method 2: Multi-system gradient descent (C++) ===
    print(f"\n=== Method 2: Multi-System Gradient Descent (C++) ===")
    k_gd = fit_gradient_descent_multi_cpp(fitters, systems, n_cpp, lr=1e-4, momentum=0.9, max_iter=5000, tol=1e-10, verbose=True)
    for key, t in global_bond_map.items():
        print(f"  k_bond[{key}] = {k_gd[t]:.6e} eV/A^2")
    for key, t in global_bond3_map.items():
        print(f"  k_bond3[{key}] = {k_gd[t]:.6e} eV/A^2")
    for key, t in global_angle_map.items():
        print(f"  k_angle[{key}] = {k_gd[t]:.6e} eV/rad^2")

    print(f"\n=== Method Comparison ===")
    print(f"  LSQ vs GD max |diff|: {np.max(np.abs(k_lsq - k_gd)):.6e}")
    
    # === Per-system Hessian comparison (C++) ===
    print(f"\n=== Per-System Hessian Comparison ===")
    k = k_lsq
    for f, sys in zip(fitters, systems):
        f.set_params(k)
        H_model = f.compute_model_hessian()
        H_diff = H_model - sys['H_ref']
        rmsd = np.sqrt(np.mean(H_diff**2))
        nref = np.linalg.norm(sys['H_ref'])
        ndiff = np.linalg.norm(H_diff)
        nmodel = np.linalg.norm(H_model)
        rel_frob = ndiff / nref * 100
        print(f"  {sys['name']:12s}: RMSD={rmsd:.4e} eV/A^2  relFrob={rel_frob:.2f}%  "
              f"||H_ref||={nref:.4e}  ||H_model||={nmodel:.4e}  ||diff||={ndiff:.4e}")
    
    # === Method 3: robust hybrid all-mode/local/internal fit (Python) ===
    print(f"\n=== Method 3: Hybrid All-Mode + Local + Internal-Coordinate Fit ===")
    for f in fitters: f.set_n_free(n_total)
    hybrid_systems = []
    for f, sys in zip(fitters, systems):
        A = FFfit_cpp.collect_sensitivity_matrices(f, extra=sys.get('dihedral_A', {}))
        hybrid_systems.append({'A': A, 'H_ref': sys['H_ref'], 'positions': sys['positions'], 'masses': sys['data']['masses'],
                               'bonds': sys['bonds'], 'angles': sys['angles']})
    k_prior = np.ones(n_total)
    for t in global_bond_map.values(): k_prior[t] = 5.0
    for t in global_bond3_map.values(): k_prior[t] = 0.1
    for t in global_angle_map.values(): k_prior[t] = 1.0
    for t in global_dihedral_map.values(): k_prior[t] = 0.1
    bounds = (-np.inf, np.inf) if args.allow_negative else (0.0, np.inf)
    k_mode, fit_diag = FFfit_cpp.fit_hybrid_hessian(
        hybrid_systems, mode_weight=args.hybrid_mode, local_weight=args.hybrid_local, internal_weight=args.hybrid_internal,
        mode_balance=args.mode_weight, mode_mixing=args.mode_mixing, frequency_floor_cm1=args.frequency_floor, local_graph_distance=args.local_graph_distance,
        prior=k_prior, regularization=args.reg, parameter_scale=k_prior, bounds=bounds)
    print(f"  solver={fit_diag['solver']} residual={fit_diag['relative_residual']:.6e} rank={fit_diag['rank']}/{n_total} condition={fit_diag['condition']:.6e}")
    for sys, info in zip(systems, fit_diag['systems']):
        print(f"  {sys['name']}: all_modes={info['n_modes']} internal_rank={info['internal_rank']}/{info['n_internal_coordinates']} rows={info['component_rows']}")
    for key, t in global_bond_map.items():
        print(f"  k_bond[{key}] = {k_mode[t]:.6e} eV/A^2")
    for key, t in global_bond3_map.items():
        print(f"  k_bond3[{key}] = {k_mode[t]:.6e} eV/A^2")
    for key, t in global_angle_map.items():
        print(f"  k_angle[{key}] = {k_mode[t]:.6e} eV/rad^2")
    for key, t in global_dihedral_map.items():
        print(f"  k_dihedral[{key}] = {k_mode[t]:.6e} eV")

    # === Mode-basis fit diagnostics ===
    print(f"\n=== Mode-Basis Fit Diagnostics ===")
    k = k_mode
    for f, sys in zip(fitters, systems):
        f.set_params(k)
        H_model = f.compute_model_hessian()
        add_dihedral_hessian(H_model, k, sys.get('dihedral_A', {}))
        H_diff = H_model - sys['H_ref']
        rmsd = np.sqrt(np.mean(H_diff**2))
        nref = np.linalg.norm(sys['H_ref'])
        ndiff = np.linalg.norm(H_diff)
        nmodel = np.linalg.norm(H_model)
        rel_frob = ndiff / nref * 100
        print(f"  {sys['name']:12s}: RMSD={rmsd:.4e} eV/A^2  relFrob={rel_frob:.2f}%  "
              f"||H_ref||={nref:.4e}  ||H_model||={nmodel:.4e}  ||diff||={ndiff:.4e}")
    
    # === Per-system frequency comparison ===
    print(f"\n=== Per-System Frequency Comparison (Mode-Basis Fit) ===")
    system_spectra_data = []
    for f, sys in zip(fitters, systems):
        f.set_params(k_mode)
        H_model = f.compute_model_hessian()
        add_dihedral_hessian(H_model, k_mode, sys.get('dihedral_A', {}))
        nz_ref, nz_model = compare_frequencies(sys['H_ref'], H_model, sys['data']['masses'], data=sys['data'], label=sys['name'], positions=sys['positions'])
        if args.plot:
            plot_spectrum(nz_ref, nz_model, sys['name'], outdir=args.plot_dir,
                          exp_spectrum=args.exp_spectrum, width=args.spectrum_width, xmax=args.spectrum_xmax)
            plot_hessian_comparison(sys['H_ref'], H_model, sys['name'], args.plot_dir)
        system_spectra_data.append((sys['name'], nz_ref, nz_model))

    if args.plot:
        plot_spectra_overlay(system_spectra_data, outdir=args.plot_dir, width=args.spectrum_width, xmax=args.spectrum_xmax)

    print(f"\n=== Done ===")

if __name__ == '__main__':
    main()
