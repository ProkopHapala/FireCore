#!/usr/bin/env python3
"""
MD with MMFFsp3 on GPU (pyOpenCL) for formic acid


- Relevant Files:
- cpp/common_resources/cl/relax_multi_mini.cl
- pyBall/OCL/MolecularDynamics.py
- pyBall/OCL/MMFF.py
- pyBall/AtomicSystem.py
- Relevant kernels: `cleanForceMMFFf4`, `getNonBond`, `getMMFFf4`, `updateAtomsMMFFf4`, `updateAtomsMMFFf4_rot`, `updateAtomsMMFFf4_RATTLE` in `cpp/common_resources/cl/relax_multi_mini.cl`
- Python wrappers: `MolecularDynamics.setup_kernels()`, `run_*()` in `pyBall/OCL/MolecularDynamics.py`
- MMFF builder: `MMFF.toMMFFsp3_loc()` populates `apos[0:natoms]` (atoms) and `apos[natoms:natoms+nnode]` (pi vectors), neighbors and params in `pyBall/OCL/MMFF.py`
- Notes: `bkNeighs` currently uploaded as -1 placeholders; recoil force accumulation is therefore disabled

More details in - `FireCore/doc/DevNotes/MD_MMFF_OCL_notes.md`
"""

import os, sys, argparse
import numpy as np
import math
import matplotlib.pyplot as plt

# Add FireCore root to PYTHONPATH
BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(BASE)

from pyBall.AtomicSystem import AtomicSystem
from pyBall.OCL.MMFF import MMFF
from pyBall.OCL.MolecularDynamics import MolecularDynamics
from pyBall.plotUtils import plotAtoms, plotBonds
from pyBall import elements


# infinite line length for numpy print options
np.set_printoptions(linewidth=np.inf)

DATA_MOL = os.path.join(BASE, 'cpp/common_resources/mol/formic_acid.mol2')

MONITOR_PROPERTY_TYPES = {
    'F': 'vector',
    'T': 'vector',
    'P': 'vector',
    'L': 'vector',
    'Rcm': 'vector',
    'Vcm': 'vector',
    'Fcm': 'vector',
    'Tcm': 'vector',
    'Ekin': 'scalar',
    'Epot': 'scalar',
    'Etotal': 'scalar',
}
DEFAULT_MONITOR_PROPS = ['F', 'P', 'L', 'Rcm', 'Vcm', 'Ekin', 'Epot', 'Etotal']


def parse_monitor_props(spec):
    if spec is None:
        return list(DEFAULT_MONITOR_PROPS)
    s = spec.strip()
    if s == '':
        return list(DEFAULT_MONITOR_PROPS)
    if s.lower() in {'none', 'off', 'false', '0'}:
        return []
    if s.lower() == 'all':
        return list(MONITOR_PROPERTY_TYPES.keys())
    props = [p.strip() for p in s.split(',') if p.strip()]
    invalid = [p for p in props if p not in MONITOR_PROPERTY_TYPES]
    if invalid:
        raise ValueError(f"Unknown monitor properties: {', '.join(invalid)}. Valid options: {sorted(MONITOR_PROPERTY_TYPES.keys())}")
    return props


def get_atom_masses(mol):
    masses = np.zeros(len(mol.enames), dtype=np.float32)
    for i, sym in enumerate(mol.enames):
        entry = elements.ELEMENT_DICT.get(sym)
        if entry is None:
            # Try common suffix patterns (e.g. "C_2", "O1") used in some MOL2 exports
            candidates = []
            if '_' in sym:
                candidates.append(sym.split('_')[0])
            stripped_digits = sym.rstrip('0123456789')
            if stripped_digits != sym:
                candidates.append(stripped_digits)
            alpha_only = ''.join(ch for ch in sym if ch.isalpha())
            if alpha_only and alpha_only != sym:
                candidates.append(alpha_only)
            for cand in candidates:
                if cand in elements.ELEMENT_DICT:
                    entry = elements.ELEMENT_DICT[cand]
                    break
        if entry is None:
            raise ValueError(f"No element data for symbol '{sym}' while building masses")
        masses[i] = float(entry[elements.index_mass])
    return masses


def compute_energies(avel_atoms, masses, aforce_atoms_full):
    v = np.asarray(avel_atoms, dtype=np.float32)
    m = np.asarray(masses, dtype=np.float32).reshape(-1, 1)
    kin = 0.5 * m * (v * v).sum(axis=1, keepdims=True)
    Ekin = float(kin.sum())
    af = np.asarray(aforce_atoms_full, dtype=np.float32)
    Epot = float(af[:, 3].sum())
    return {
        'Ekin': Ekin,
        'Epot': Epot,
        'Etotal': Ekin + Epot,
    }


def finalize_monitor_series(monitor_data):
    series = {}
    for key, values in monitor_data.items():
        if not values:
            continue
        arrays = [np.array(v, dtype=np.float32) for v in values]
        if arrays[0].ndim == 0:
            series[key] = np.asarray(arrays, dtype=np.float32).reshape(-1)
        else:
            series[key] = np.stack(arrays, axis=0)
    return series


def plot_monitor_series(series, props, show=True):
    if not props:
        return None, None
    fig, axes = plt.subplots(len(props), 1, figsize=(8, 2.6 * len(props)), sharex=True)
    if isinstance(axes, np.ndarray):
        axes_list = list(axes.ravel())
    elif isinstance(axes, (list, tuple)):
        axes_list = list(axes)
    else:
        axes_list = [axes]
    steps = None
    for ax, name in zip(axes_list, props):
        data = series.get(name)
        if data is None:
            ax.set_visible(False)
            continue
        if steps is None:
            steps = np.arange(data.shape[0])
        if data.ndim == 1:
            ax.plot(steps, data, label=name)
        else:
            comps = ['x', 'y', 'z']
            ncomp = data.shape[1]
            for j in range(ncomp):
                comp_label = comps[j] if j < len(comps) else str(j)
                ax.plot(steps, data[:, j], label=f"{name}_{comp_label}")
            ax.legend(loc='best', fontsize=8)
        ax.set_ylabel(name)
        ax.grid(True, alpha=0.3)
    if axes_list:
        axes_list[-1].set_xlabel('step')
    fig.tight_layout()
    if show:
        plt.show()
    return fig, axes_list


def build_mmff_from_mol(mol2_path):
    mol = AtomicSystem(fname=mol2_path)
    # Ensure neighbors present
    if getattr(mol, 'ngs', None) is None:
        mol.neighs()
    mm = MMFF(bTorsion=False, verbosity=1)
    mm.toMMFFsp3_loc(mol=mol, atom_types=mm.atom_types, bRealloc=True, bEPairs=False, bUFF=False)
    return mol, mm


def _stats(name, arr):
    import numpy as _np
    A = _np.asarray(arr)
    finite = _np.isfinite(A)
    n_nan = int(_np.isnan(A).sum())
    n_inf = int(_np.isinf(A).sum())
    mn = float(_np.min(A)) if A.size else 0.0
    mx = float(_np.max(A)) if A.size else 0.0
    l2 = float(_np.linalg.norm(A)) if A.size else 0.0
    print(f"{name}: shape={A.shape} min={mn:.6e} max={mx:.6e} ||.||_2={l2:.6e} NaN={n_nan} Inf={n_inf} finite={finite.all()}")


def zero_dynamic_buffers(md):
    """Zero avel, aforce, cvf, fneigh for all systems."""
    import numpy as _np
    z_vec = _np.zeros(md.nSystems * md.nvecs * 4, dtype=_np.float32)
    z_fng = _np.zeros(md.nSystems * md.nnode * 8, dtype=_np.float32)
    md.toGPU('aforce', z_vec)
    md.toGPU('avel',   z_vec)
    md.toGPU('cvf',    z_vec)
    md.toGPU('fneigh', z_fng)


def fetch_arrays(md, mm):
    """Download common arrays. Returns dict with atoms-only slices for convenience."""
    apos, aforce = md.download_results()
    avel  = md.download_buf('avel').reshape(md.nSystems, md.nvecs, 4)
    fng   = md.download_buf('fneigh').reshape(md.nSystems, md.nnode, 8)
    cvf   = md.download_buf('cvf').reshape(md.nSystems, md.nvecs, 4)
    # Atoms-only views from system 0
    atoms = {
        'apos_atoms':      apos[:mm.natoms, :3],
        'avel_atoms':      avel[0, :mm.natoms, :3],
        'afor_atoms':      aforce[:mm.natoms, :3],
        'afor_atoms_full': aforce[:mm.natoms, :4],
    }
    return {'apos': apos, 'aforce': aforce, 'avel': avel, 'fneigh': fng, 'cvf': cvf, **atoms}


def run_step(md, do_clean=True, do_nb=True, do_mmff=True, mode='basic'):
    """Run one force/integrator step according to flags."""
    if do_clean:
        md.run_cleanForceMMFFf4()
    if do_nb:
        md.run_getNonBond()
    if do_mmff:
        md.run_getMMFFf4()
    if mode == 'rot':
        md.run_updateAtomsMMFFf4_rot()
    elif mode == 'rattle':
        md.run_updateAtomsMMFFf4_RATTLE()
    elif mode == 'basic':
        md.run_updateAtomsMMFFf4()
    elif mode == 'none':
        pass


def compute_totals(apos_atoms, avel_atoms, aforce_atoms, masses=None, origin=None):
    """
    Compute global conservation diagnostics for a single frame.
    Returns dict with keys: F, T, P, L, Rcm, Vcm, Fcm, Tcm.
    """
    r = np.asarray(apos_atoms, dtype=np.float32)
    v = np.asarray(avel_atoms, dtype=np.float32)
    f = np.asarray(aforce_atoms, dtype=np.float32)
    n = r.shape[0]
    if masses is None:
        m = np.ones((n, 1), dtype=np.float32)
    else:
        m = np.asarray(masses, dtype=np.float32).reshape(n, 1)

    mtot = float(m.sum())
    if mtot <= 0.0:
        raise ValueError("Total mass must be positive to compute COM diagnostics")

    rcm = (m * r).sum(axis=0) / mtot
    vcm = (m * v).sum(axis=0) / mtot

    Ftot = f.sum(axis=0)
    P    = (m * v).sum(axis=0)

    r0 = np.zeros(3, dtype=np.float32) if origin is None else np.asarray(origin, dtype=np.float32)
    rrel = r - r0
    Ttot = np.cross(rrel, f).sum(axis=0)
    L    = np.cross(rrel, (m * v)).sum(axis=0)

    Fcm = Ftot.copy()
    Tcm = np.cross(rcm - r0, Ftot)

    return {
        'F': Ftot,
        'T': Ttot,
        'P': P,
        'L': L,
        'Rcm': rcm,
        'Vcm': vcm,
        'Fcm': Fcm,
        'Tcm': Tcm,
    }


def print_totals(label, totals):
    F, T, P, L = totals['F'], totals['T'], totals['P'], totals['L']
    def fmt(x): return f"({x[0]: .3e}, {x[1]: .3e}, {x[2]: .3e})"
    extras = []
    for key in ('Ekin', 'Epot', 'Etotal'):
        if key in totals:
            extras.append(f"{key}={totals[key]:.6e}")
    extra_str = " " + " ".join(extras) if extras else ""
    print(f"[{label}] F={fmt(F)} | T={fmt(T)} | P={fmt(P)} | L={fmt(L)}{extra_str}")


def plot_trajectories(trj_atoms, dim='xy', labels=None, title='Trajectories', save_path=None, show=True):
    """
    Plot atom trajectories from trj_atoms with shape (nSteps, natoms, 3).
    dim: 'xy' | 'xz' | 'yz'
    """
    proj = {'xy': (0,1), 'xz': (0,2), 'yz': (1,2)}.get(dim, (0,1))
    ii, jj = proj
    nsteps, natoms, _ = trj_atoms.shape
    fig, ax = plt.subplots(figsize=(6,6))
    for a in range(natoms):
        xy = trj_atoms[:, a, :]
        ax.plot(xy[:,ii], xy[:,jj], '-', lw=1, label=(labels[a] if labels is not None else None), zorder=1)
        ax.plot(xy[-1,ii], xy[-1,jj], 'o', ms=3, zorder=2)
    ax.set_xlabel(['x','y','z'][ii])
    ax.set_ylabel(['x','y','z'][jj])
    ax.set_title(title)
    if labels is not None and len(labels) <= 12:
        ax.legend(loc='best', fontsize=8)
    ax.set_aspect('equal', adjustable='box')
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150)
        print(f"Saved plot to {save_path}")
    if show:
        plt.show()
    return fig, ax


def run_md(steps=100, mode='basic', dt=0.01, damp=0.98, flim=10.0, do_clean=True, do_nb=True, do_mmff=True,
           print_stats=False, print_arrays=False, record=False, plot=False, plot_dim='xy', save_plot=None, save_trj=None,
           plot_bonds=False, subtract_vdw=False, plot_label_mode='none', print_params=False,
           monitor=True, monitor_props=None, monitor_plot=False, monitor_save=None, monitor_save_data=None):
    mol, mm = build_mmff_from_mol(DATA_MOL)

    md = MolecularDynamics(nloc=32)
    md.realloc(mm, nSystems=1)

    # MD params
    mm.dt = dt
    mm.damp = damp
    mm.Flimit = flim

    md.upload_all_systems()

    masses = get_atom_masses(mol)

    # Zero dynamic buffers to avoid uninitialized data (helps prevent NaNs)
    zero_dynamic_buffers(md)
    md.setup_kernels()
    # Update bonded/non-bonded subtraction flag before regenerating kernel argument lists
    md.kernel_params['bSubtractVdW'] = np.int32(1 if subtract_vdw else 0)
    md.kernel_args_getMMFFf4 = md.generate_kernel_args('getMMFFf4')
    md.kernel_args_runMD     = md.generate_kernel_args('runMD')

    # Optional recorders for diagnostics
    if monitor_props is None:
        monitor_props = list(DEFAULT_MONITOR_PROPS) if monitor else []
    monitor_enabled = monitor and len(monitor_props) > 0
    trj = [] if record else None
    need_history = record or monitor_enabled
    totals_hist = [] if need_history else None
    monitor_data = {name: [] for name in monitor_props} if monitor_enabled else None

    # Single-step MD loop using chosen integrator
    for i in range(steps):
        run_step(md, do_clean=do_clean, do_nb=do_nb, do_mmff=do_mmff, mode=mode)
        if record or monitor_enabled:
            buf = fetch_arrays(md, mm)
            apos_step = buf['apos_atoms'].copy()
            avel_step = buf['avel_atoms'].copy()
            afor_step = buf['afor_atoms'].copy()
            afor_full_step = buf['afor_atoms_full'].copy()
            totals = compute_totals(apos_step, avel_step, afor_step, masses=masses)
            totals.update(compute_energies(avel_step, masses, afor_full_step))
            if totals_hist is not None:
                totals_hist.append(totals)
            if record:
                trj.append(apos_step)
            if monitor_enabled:
                for name in monitor_props:
                    val = totals.get(name)
                    if val is None:
                        continue
                    if np.isscalar(val):
                        monitor_data[name].append(float(val))
                    else:
                        monitor_data[name].append(np.array(val, dtype=np.float32))

    buf = fetch_arrays(md, mm)
    apos, aforce = buf['apos'], buf['aforce']
    if print_stats or print_arrays:
        avel, fng, cvf = buf['avel'], buf['fneigh'], buf['cvf']
        if print_stats:
            _stats('apos', apos)
            _stats('avel', avel)
            _stats('aforce', aforce)
            _stats('fneigh', fng)
            _stats('cvf', cvf)
        if print_arrays:
            print('apos:\n', apos)
            print('avel:\n', avel)
            print('aforce:\n', aforce)
            print('fneigh:\n', fng)
            print('cvf:\n', cvf)

    # Conservation diagnostics over trajectory
    if totals_hist:
        F = np.array([x['F'] for x in totals_hist])
        T = np.array([x['T'] for x in totals_hist])
        P = np.array([x['P'] for x in totals_hist])
        L = np.array([x['L'] for x in totals_hist])
        def dmax(A):
            d = A - A[0]
            return float(np.max(np.linalg.norm(d, axis=1)))
        energy_keys = ('Ekin', 'Epot', 'Etotal')
        energy_drifts = []
        for key in energy_keys:
            if key in totals_hist[0]:
                arr = np.array([x[key] for x in totals_hist], dtype=np.float64)
                energy_drifts.append(f"|Δ{key}|max={np.max(np.abs(arr - arr[0])):.3e}")
        drift_msg = f"Conservation drift over {steps} steps: |ΔF|max={dmax(F):.3e} |ΔT|max={dmax(T):.3e} |ΔP|max={dmax(P):.3e} |ΔL|max={dmax(L):.3e}"
        if energy_drifts:
            drift_msg += " " + " ".join(energy_drifts)
        print(drift_msg)
        if record:
            print_totals('totals[0]', totals_hist[0])
            print_totals('totals[-1]', totals_hist[-1])

    # Save/plot trajectory
    if record and trj is not None:
        trj = np.stack(trj, axis=0)
        if save_trj:
            np.save(save_trj, trj)
            print(f"Saved trajectory to {save_trj} with shape {trj.shape}")
        if plot:
            fig, ax = plot_trajectories(trj, dim=plot_dim, title=f'Trajectories ({plot_dim})', save_path=None, show=False)
            axes = {'xy': (0,1), 'xz': (0,2), 'yz': (1,2)}.get(plot_dim, (0,1))
            if plot_bonds:
                try:
                    if getattr(mol, 'bonds', None) is None or len(mol.bonds) == 0:
                        mol.findBonds(Rcut=3.0, RvdwCut=0.5)
                    bonds_arr = np.asarray(mol.bonds, dtype=np.int32)
                    if bonds_arr.ndim == 1:
                        bonds_arr = bonds_arr.reshape(1, -1)
                    links = bonds_arr[:, :2]
                    last = trj[-1]
                    labels = None
                    if plot_label_mode == 'number':
                        labels = [str(i) for i in range(last.shape[0])]
                    elif plot_label_mode == 'type' and getattr(mol, 'enames', None) is not None:
                        labels = mol.enames
                    plotAtoms(apos=last, es=getattr(mol, 'enames', None), axes=axes, sizes=60., colors='#404040', marker='o', labels=labels, bNumbers=False)
                    plotBonds(links=links, ps=last, axes=axes, colors='k', lws=1.5)
                except Exception as e:
                    print(f"WARNING: plot_bonds failed: {e}")
            if save_plot:
                fig.savefig(save_plot, dpi=150)
                print(f"Saved plot to {save_plot}")
            else:
                plt.show()

    if monitor_enabled and monitor_data:
        series = finalize_monitor_series(monitor_data)
        if monitor_save_data:
            np.savez(monitor_save_data, **series)
            print(f"Saved monitor data to {monitor_save_data}")
        if monitor_plot:
            fig, _ = plot_monitor_series(series, monitor_props, show=not monitor_save)
            if monitor_save and fig is not None:
                fig.savefig(monitor_save, dpi=150)
                print(f"Saved monitor plot to {monitor_save}")
        elif monitor_save:
            fig, _ = plot_monitor_series(series, monitor_props, show=False)
            if fig is not None:
                fig.savefig(monitor_save, dpi=150)
                print(f"Saved monitor plot to {monitor_save}")

    na = mm.natoms
    print(f"Done. Steps={steps}  natoms={na}  first-atom pos={apos[0,:3]}  force={aforce[0,:3]}")
    return apos, aforce

"""
Examples (omit defaults):
- Minimal MD with diagnostics
  python test_MMFFsp3_pyOCL.py --steps 50 --print-stats 1
- Use rotational integrator
  python test_MMFFsp3_pyOCL.py --steps 50 --mode rot
- Forces only (no motion)
  python test_MMFFsp3_pyOCL.py --steps 10 --mode none --print-stats 1
- Record and plot XY trajectories
  python test_MMFFsp3_pyOCL.py --steps 100 --record 1 --plot 1 --plot-dim xy 
- Monotor total invariants like momentum, energy, torque etc.
  python test_MMFFsp3_pyOCL.py --steps 200 --monitor 1 --monitor-plot 1
"""

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--mode',         choices=['basic','rot','rattle','none'], default='basic')
    ap.add_argument('--molecule',                 default=DATA_MOL)
    ap.add_argument('--steps',        type=int,   default=100)
    ap.add_argument('--dt',           type=float, default=0.01)
    ap.add_argument('--damp',         type=float, default=0.98)
    ap.add_argument('--flim',         type=float, default=10.0)
    ap.add_argument('--do-clean',     type=int,   default=1, help='Run cleanForceMMFFf4 before forces')
    ap.add_argument('--do-nb',        type=int,   default=0, help='Run getNonBond each step')
    ap.add_argument('--do-mmff',      type=int,   default=1, help='Run getMMFFf4 (bonded) each step')
    ap.add_argument('--print-stats',  type=int,   default=1, help='Print min/max/norm/NaN counts for buffers at end')
    ap.add_argument('--print-arrays', type=int,   default=1, help='Print full arrays at end')
    ap.add_argument('--record',       type=int,   default=1, help='Record trajectory and totals each step')
    ap.add_argument('--plot',         type=int,   default=0, help='Plot trajectories')
    ap.add_argument('--plot-dim',     type=str,   default='xy', choices=['xy','xz','yz'], help='Projection for plotting')
    ap.add_argument('--save-plot',    type=str,   default=None, help='Path to save trajectory plot (png)')
    ap.add_argument('--save-trj',     type=str,   default=None, help='Path to save trajectory (npy)')
    ap.add_argument('--plot-bonds',   type=int,   default=1, help='Overlay molecule bonds on trajectory plot')
    ap.add_argument('--subtract-vdw', type=int,   default=0, help='Subtract bonded LJ contributions inside getMMFFf4 (requires non-bonded forces)')
    ap.add_argument('--plot-labels',  type=str,   default='number', choices=['none','number','type'], help='Annotate atoms when plotting trajectories')
    ap.add_argument('--print-params', type=int,   default=0, help='Dump neighbor and bonded parameter arrays before GPU upload')
    ap.add_argument('--monitor',      type=int,   default=1, help='Collect invariant diagnostics each step')
    ap.add_argument('--monitor-props', type=str,  default=','.join(DEFAULT_MONITOR_PROPS), help='Comma separated invariants to track (or all/none)')
    ap.add_argument('--monitor-plot', type=int,   default=1, help='Plot monitored invariants at the end')
    ap.add_argument('--save-monitor', type=str,   default=None, help='Path to save monitor plot (png)')
    ap.add_argument('--save-monitor-data', type=str, default=None, help='Path to save monitor data (npz)')
    args = ap.parse_args()

    if os.path.abspath(args.molecule) != os.path.abspath(DATA_MOL):
        DATA_MOL = args.molecule

    monitor_props = parse_monitor_props(args.monitor_props)

    run_md(
        steps=args.steps, mode=args.mode, dt=args.dt, damp=args.damp, flim=args.flim,
        do_clean=bool(args.do_clean), do_nb=bool(args.do_nb), do_mmff=bool(args.do_mmff),
        print_stats=bool(args.print_stats), print_arrays=bool(args.print_arrays),
        record=bool(args.record), plot=bool(args.plot), plot_dim=args.plot_dim,
        save_plot=args.save_plot, save_trj=args.save_trj,
        plot_bonds=bool(args.plot_bonds), subtract_vdw=bool(args.subtract_vdw),
        plot_label_mode=args.plot_labels, print_params=bool(args.print_params),
        monitor=bool(args.monitor), monitor_props=monitor_props,
        monitor_plot=bool(args.monitor_plot), monitor_save=args.save_monitor,
        monitor_save_data=args.save_monitor_data
    )
