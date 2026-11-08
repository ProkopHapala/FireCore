#!/usr/bin/env python3
"""Systematic study of FIRE optimizer settings for eFF relaxation.

DOFs have very different stiffness (electron size ~1e5 vs positions ~1e2),
so we test: dt_max, damping, mass schemes (unit/physical/Hessian-preconditioned),
and FIRE internals (finc/fdec/falpha/minLastNeg/dtmin_frac).

Outputs to fire_study_out/.
"""
import os, sys, time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))
import fire_study_lib as fs
from pyBall import eFF as eff

OUT = os.path.join(os.path.dirname(__file__), 'fire_study_out')
os.makedirs(OUT, exist_ok=True)

XYZ_H2O = os.path.join(os.path.dirname(__file__), 'H2O_spins.xyz')
MARKS = (1e-3, 1e-4, 1e-5)
NMAX = 30000


def iters_to(tr, fc):
    i = np.nonzero(tr['Fu'] < fc)[0]
    return int(i[0]) if len(i) else -1


def reload(xyz):
    na, ne, fixed = fs.load_e(xyz)
    return na, ne, fixed


def phys_masses(na, ne):
    return fs.phys_masses(na, ne)


def prec_masses(na, ne, fixed, k_target=30.0):
    return fs.prec_masses(na, ne, fixed, k_target=k_target)


# ============ Stage 1: baseline + Hessian spectrum ============
def stage1_baseline(xyz):
    print("\n=== Stage 1: baseline (dt=0.001, damp=0.1, unit mass) ===")
    na, ne, fixed = reload(xyz)
    t0 = time.time()
    tr = fs.fire(NMAX, dt=0.001, damping=0.1, f_limit=100.0, Fconv=1e-3, fixed=fixed)
    print(f"wall {time.time()-t0:.1f}s")
    fs.summarize(tr, 'baseline H2O_spins')
    np.savez(os.path.join(OUT, 'baseline_trace.npz'), **{k: v for k, v in tr.items() if isinstance(v, np.ndarray)})

    # Hessian at relaxed geometry
    H = fs.hessian(h=1e-4, fixed=fixed)
    ev, U = np.linalg.eigh(H)
    np.save(os.path.join(OUT, 'hessian_eigvals.npy'), ev)
    g = fs.groups(na, ne)
    print(f"Hessian eigenvalues: min {ev.min():.3f} max {ev.max():.3f}  cond(|max|/|min nonzero|) {ev.max()/max(ev[np.abs(ev)>1e-6].min(),1e-9):.2e}")
    kd = np.diag(H)
    for name in ('apos', 'epos', 'esize'):
        d = kd[g[name]]
        print(f"  diag k[{name:5s}]: min {d.min():12.3e} max {d.max():12.3e} med {np.median(d):12.3e}")
    return tr


# ============ Stage 2: dt x damping grid, unit masses ============
def stage2_grid(xyz):
    print("\n=== Stage 2: dt_max x damping grid (unit masses) ===")
    dts    = [5e-4, 1e-3, 2e-3, 3e-3, 5e-3, 1e-2, 2e-2]
    damps  = [0.02, 0.05, 0.1, 0.2, 0.4]
    conv = np.full((len(damps), len(dts)), -1)
    Efin = np.full((len(damps), len(dts)), np.nan)
    for iy, dm in enumerate(damps):
        for ix, dt in enumerate(dts):
            na, ne, fixed = reload(xyz)
            tr = fs.fire(NMAX, dt=dt, damping=dm, f_limit=100.0, Fconv=1e-6, fixed=fixed)
            n = iters_to(tr, 1e-3)
            conv[iy, ix] = n if n > 0 else NMAX
            Efin[iy, ix] = tr['E_final']
            print(f"  dt={dt:g} damp={dm:g} -> iters(1e-3)={n:6d} E={tr['E_final']:.3f} {'' if np.isfinite(tr['E_final']) else 'DIVERGED'}")
    np.savez(os.path.join(OUT, 'grid_dt_damp.npz'), dts=dts, damps=damps, conv=conv, Efin=Efin)
    fig, ax = plt.subplots(figsize=(8, 5))
    im = ax.imshow(np.log10(conv), aspect='auto', origin='lower', cmap='viridis_r')
    ax.set_xticks(range(len(dts)), [f'{d:g}' for d in dts]); ax.set_yticks(range(len(damps)), [f'{d:g}' for d in damps])
    ax.set_xlabel('dt_max'); ax.set_ylabel('damp_max'); ax.set_title('iters to |F|<1e-3 (log10, unit masses)')
    for iy in range(len(damps)):
        for ix in range(len(dts)):
            ax.text(ix, iy, f'{conv[iy,ix]/1000:.1f}k', ha='center', va='center', fontsize=7, color='w')
    fig.colorbar(im); fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'grid_dt_damp.png'), dpi=130)
    print("REVIEW: fire_study_out/grid_dt_damp.png")
    return dts, damps, conv


# ============ Stage 3: mass schemes ============
# each entry: (tag, mass_mode, k_target, dt, damp)
# mass_mode: 'unit' | 'phys' | 'prec'
STAGE3 = [
    ('unit   dt=0.002', 'unit', 0,      0.002, 0.1),
    ('phys   dt=0.005', 'phys', 0,      0.005, 0.1),
    ('phys   dt=0.010', 'phys', 0,      0.010, 0.1),
    ('phys   dt=0.020', 'phys', 0,      0.020, 0.1),
    ('prec1e2 dt=0.050','prec', 100.0,  0.050, 0.1),
    ('prec1e2 dt=0.100','prec', 100.0,  0.100, 0.1),
    ('prec1e2 dt=0.200','prec', 100.0,  0.200, 0.1),
    ('prec1e3 dt=0.020','prec', 1000.0, 0.020, 0.1),
    ('prec1e3 dt=0.050','prec', 1000.0, 0.050, 0.1),
    ('prec1e3 dt=0.100','prec', 1000.0, 0.100, 0.1),
    ('prec1e4 dt=0.005','prec', 1e4,    0.005, 0.1),
    ('prec1e4 dt=0.010','prec', 1e4,    0.010, 0.1),
    ('prec1e4 dt=0.020','prec', 1e4,    0.020, 0.1),
    ('prec1e4 dt=0.050','prec', 1e4,    0.050, 0.1),
]

def stage3_masses(xyz, entries=STAGE3):
    print("\n=== Stage 3: mass schemes (iters to |F|_unweighted < marks) ===")
    results = {}
    for tag, mode, kt, dt, dm in entries:
        na, ne, fixed = reload(xyz)
        if mode == 'unit':
            invM = np.ones(eff.nDOFs)
        elif mode == 'phys':
            invM = phys_masses(na, ne); invM[fixed] = 0.0
        else:
            invM, kd = prec_masses(na, ne, fixed, k_target=kt)
        tr = fs.fire(NMAX, dt=dt, damping=dm, invM=invM, f_limit=100.0, Fconv=1e-6, fixed=fixed)
        n = iters_to(tr, 1e-3); n4 = iters_to(tr, 1e-4); n5 = iters_to(tr, 1e-5)
        results[tag] = tr
        print(f"  {tag:22s} iters(1e-3)={n:6d} (1e-4)={n4:6d} (1e-5)={n5:6d} E={tr['E_final']:.3f}")
    return results


# ============ Stage 4: FIRE internals ============
def stage4_internals(xyz, dt, dm, invM_mode='prec', k_target=100.0):
    print("\n=== Stage 4: FIRE internals ===")
    base = dict(finc=1.1, fdec=0.5, falpha=0.8, minLastNeg=3, dtmin_frac=0.1)
    sweeps = [
        ('finc',        [1.02, 1.05, 1.1, 1.2, 1.5]),
        ('fdec',        [0.3, 0.5, 0.7, 0.9]),
        ('falpha',      [0.5, 0.7, 0.8, 0.9, 0.95, 0.99]),
        ('minLastNeg',  [1, 3, 5, 10, 20]),
        ('dtmin_frac',  [0.01, 0.05, 0.1, 0.3]),
    ]
    results = {}
    for pname, vals in sweeps:
        for v_ in vals:
            kw = dict(base); kw[pname] = v_
            na, ne, fixed = reload(xyz)
            invM = np.ones(eff.nDOFs)
            if invM_mode == 'prec':
                invM, _ = prec_masses(na, ne, fixed, k_target=k_target)
            tr = fs.fire(NMAX, dt=dt, damping=dm, invM=invM, f_limit=100.0, Fconv=1e-6, fixed=fixed, **kw)
            n = iters_to(tr, 1e-3)
            print(f"  {pname}={v_:g} -> iters(1e-3)={n:6d} E={tr['E_final']:.3f}")
            results[f'{pname}={v_}'] = (n, tr)
    # smooth variant
    for cmn, cmx in [(-0.1, 0.1), (-0.3, 0.05), (-0.5, 0.5)]:
        na, ne, fixed = reload(xyz)
        invM, _ = prec_masses(na, ne, fixed, k_target=k_target) if invM_mode == 'prec' else (np.ones(eff.nDOFs), None)
        tr = fs.fire(NMAX, dt=dt, damping=dm, invM=invM, f_limit=100.0, Fconv=1e-6, fixed=fixed, smooth=True, cvf_min=cmn, cvf_max=cmx, **base)
        n = iters_to(tr, 1e-3)
        print(f"  smooth cvf=({cmn},{cmx}) -> iters(1e-3)={n:6d} E={tr['E_final']:.3f}")
        results[f'smooth{cmn},{cmx}'] = (n, tr)
    return results


# ============ Stage 5: trace plots for selected configs ============
def stage5_traces(xyz, configs):
    print("\n=== Stage 5: |F| traces ===")
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.5), sharex=True)
    for tag, kw in configs:
        na, ne, fixed = reload(xyz)
        invM = kw.pop('invM', None)
        tr = fs.fire(NMAX, fixed=fixed, Fconv=1e-6, invM=invM, **kw)
        axs[0].semilogy(tr['Fu'], label=tag, lw=1)
        axs[1].semilogy(tr['Fpos_e'], lw=1); axs[2].semilogy(tr['Fsize'], lw=1)
        n = iters_to(tr, 1e-3)
        print(f"  {tag:40s} iters(1e-3)={n:6d}")
    axs[0].set_title('|F| unweighted'); axs[1].set_title('|F| epos'); axs[2].set_title('|F| esize')
    for ax in axs: ax.set_xlabel('iter'); ax.grid(alpha=0.3)
    axs[0].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, 'traces_compare.png'), dpi=130)
    print("REVIEW: fire_study_out/traces_compare.png")


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    fs.set_atom_params()

    if which in ('all', '1'):
        stage1_baseline(XYZ_H2O)
    if which in ('all', '2'):
        stage2_grid(XYZ_H2O)
    if which in ('all', '3'):
        stage3_masses(XYZ_H2O)
    if which in ('all', '4'):
        stage4_internals(XYZ_H2O, 0.01, 0.1)
    if which in ('all', '5'):
        na, ne, fixed = reload(XYZ_H2O)
        invMp, _ = prec_masses(na, ne, fixed, k_target=100.0)
        stage5_traces(XYZ_H2O, [
            ('default dt=0.001 damp=0.1 unit', dict(dt=0.001, damping=0.1, f_limit=100.0)),
            ('dt=0.01 damp=0.1 unit',          dict(dt=0.01,  damping=0.1, f_limit=100.0)),
            ('dt=0.05 damp=0.1 prec-k100',     dict(dt=0.05,  damping=0.1, f_limit=100.0, invM=invMp)),
        ])
