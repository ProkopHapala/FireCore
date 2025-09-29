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

Examples (omit defaults):
- Minimal MD with diagnostics
  python test_MD_OCL_formic_acid.py --steps 50 --print-stats 1
- Use rotational integrator
  python test_MD_OCL_formic_acid.py --steps 50 --mode rot
- Forces only (no motion)
  python test_MD_OCL_formic_acid.py --steps 10 --mode none --print-stats 1
- Record and plot XY trajectories
  python test_MD_OCL_formic_acid.py --steps 100 --record 1 --plot 1 --plot-dim xy --save-plot traj.png
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
from pyBall.plotUtils import plotBonds

DATA_MOL = os.path.join(BASE, 'cpp/common_resources/mol/formic_acid.mol2')


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
        'apos_atoms':  apos[:mm.natoms, :3],
        'avel_atoms':  avel[0, :mm.natoms, :3],
        'afor_atoms':  aforce[:mm.natoms, :3],
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
    - Net force:    Ftot = sum_i f_i
    - Net torque:   Ttot = sum_i r_i x f_i        (about 'origin' or center of mass if provided separately)
    - Momentum:     P    = sum_i m_i v_i
    - Ang. mom.:    L    = sum_i r_i x (m_i v_i)
    If masses is None, unit masses are assumed.
    """
    r = np.asarray(apos_atoms, dtype=np.float32)[:, :3]
    v = np.asarray(avel_atoms, dtype=np.float32)[:, :3]
    f = np.asarray(aforce_atoms, dtype=np.float32)[:, :3]
    n = r.shape[0]
    if masses is None:
        m = np.ones((n, 1), dtype=np.float32)
    else:
        m = np.asarray(masses, dtype=np.float32).reshape(n, 1)

    Ftot = f.sum(axis=0)
    P    = (m * v).sum(axis=0)

    r0 = np.zeros(3, dtype=np.float32) if origin is None else np.asarray(origin, dtype=np.float32)
    rrel = r - r0
    Ttot = np.cross(rrel, f).sum(axis=0)
    L    = np.cross(rrel, (m * v)).sum(axis=0)

    return {
        'F': Ftot,
        'T': Ttot,
        'P': P,
        'L': L,
    }


def print_totals(label, totals):
    F, T, P, L = totals['F'], totals['T'], totals['P'], totals['L']
    def fmt(x): return f"({x[0]: .3e}, {x[1]: .3e}, {x[2]: .3e})"
    print(f"[{label}] F={fmt(F)} | T={fmt(T)} | P={fmt(P)} | L={fmt(L)}")


def plot_trajectories(trj_atoms, dim='xy', labels=None, title='Trajectories', save_path=None):
    """
    Plot atom trajectories from trj_atoms with shape (nSteps, natoms, 3).
    dim: 'xy' | 'xz' | 'yz'
    """
    proj = {'xy': (0,1), 'xz': (0,2), 'yz': (1,2)}.get(dim, (0,1))
    ii, jj = proj
    nsteps, natoms, _ = trj_atoms.shape
    plt.figure(figsize=(6,6))
    for a in range(natoms):
        xy = trj_atoms[:, a, :]
        plt.plot(xy[:,ii], xy[:,jj], '-', lw=1, label=(labels[a] if labels is not None else None))
        plt.plot(xy[-1,ii], xy[-1,jj], 'o', ms=3)
    plt.xlabel(['x','y','z'][ii]); plt.ylabel(['x','y','z'][jj])
    plt.title(title)
    if labels is not None and len(labels) <= 12:
        plt.legend(loc='best', fontsize=8)
    plt.axis('equal'); plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150)
        print(f"Saved plot to {save_path}")
    else:
        plt.show()


def run_md(steps=100, mode='basic', dt=0.01, damp=0.98, flim=10.0, do_clean=True, do_nb=True, do_mmff=True,
           print_stats=False, print_arrays=False, record=False, plot=False, plot_dim='xy', save_plot=None, save_trj=None,
           plot_bonds=False):
    mol, mm = build_mmff_from_mol(DATA_MOL)

    md = MolecularDynamics(nloc=32)
    md.realloc(mm, nSystems=1)

    # MD params
    mm.dt = dt
    mm.damp = damp
    mm.Flimit = flim

    md.upload_all_systems()

    # Zero dynamic buffers to avoid uninitialized data (helps prevent NaNs)
    zero_dynamic_buffers(md)
    md.setup_kernels()

    # Optional recorders for diagnostics
    trj = [] if record else None
    totals_hist = [] if record else None

    # Single-step MD loop using chosen integrator
    for i in range(steps):
        run_step(md, do_clean=do_clean, do_nb=do_nb, do_mmff=do_mmff, mode=mode)
        if record:
            # Download current state for diagnostics (atoms only)
            buf = fetch_arrays(md, mm)
            apos_step = buf['apos_atoms'].copy()
            avel_step = buf['avel_atoms'].copy()
            afor_step = buf['afor_atoms'].copy()
            trj.append(apos_step)
            totals_hist.append(compute_totals(apos_step, avel_step, afor_step))

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
    if record and totals_hist:
        F = np.array([x['F'] for x in totals_hist])
        T = np.array([x['T'] for x in totals_hist])
        P = np.array([x['P'] for x in totals_hist])
        L = np.array([x['L'] for x in totals_hist])
        def dmax(A):
            d = A - A[0]
            return float(np.max(np.linalg.norm(d, axis=1)))
        print(f"Conservation drift over {steps} steps: |ΔF|max={dmax(F):.3e} |ΔT|max={dmax(T):.3e} |ΔP|max={dmax(P):.3e} |ΔL|max={dmax(L):.3e}")
        # Print first/last for quick inspection
        print_totals('totals[0]', totals_hist[0])
        print_totals('totals[-1]', totals_hist[-1])

    # Save/plot trajectory
    if record and trj is not None:
        trj = np.stack(trj, axis=0)
        if save_trj:
            np.save(save_trj, trj)
            print(f"Saved trajectory to {save_trj} with shape {trj.shape}")
        if plot:
            plot_trajectories(trj, dim=plot_dim, title=f'Trajectories ({plot_dim})', save_path=None)
            # Overlay molecular skeleton (bonds) from last frame using plotUtils
            try:
                # Collect bonds as links (nb,2)
                if hasattr(mol, 'bonds') and plot_bonds:
                    links = np.array([[b.i, b.j] for b in mol.bonds], dtype=np.int32)
                    last = trj[-1]
                    axmap = {'xy': (0,1), 'xz': (0,2), 'yz': (1,2)}
                    axes = axmap.get(plot_dim, (0,1))
                    plotBonds(links=links, ps=last, axes=axes, colors='k', lws=1.5)
                    if save_plot:
                        plt.savefig(save_plot, dpi=150)
                        print(f"Saved plot with bonds to {save_plot}")
                    else:
                        plt.show()
                else:
                    if save_plot:
                        plt.savefig(save_plot, dpi=150)
                        print(f"Saved plot to {save_plot}")
            except Exception as e:
                print(f"WARNING: plot_bonds failed: {e}")

    na = mm.natoms
    print(f"Done. Steps={steps}  natoms={na}  first-atom pos={apos[0,:3]}  force={aforce[0,:3]}")
    return apos, aforce


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--molecule', default=DATA_MOL)
    ap.add_argument('--steps', type=int, default=100)
    ap.add_argument('--mode', choices=['basic','rot','rattle','none'], default='basic')
    ap.add_argument('--dt', type=float, default=0.01)
    ap.add_argument('--damp', type=float, default=0.98)
    ap.add_argument('--flim', type=float, default=10.0)
    ap.add_argument('--do-clean', type=int, default=1, help='Run cleanForceMMFFf4 before forces')
    ap.add_argument('--do-nb', type=int, default=1, help='Run getNonBond each step')
    ap.add_argument('--do-mmff', type=int, default=1, help='Run getMMFFf4 (bonded) each step')
    ap.add_argument('--print-stats', type=int, default=1, help='Print min/max/norm/NaN counts for buffers at end')
    ap.add_argument('--print-arrays', type=int, default=1, help='Print full arrays at end')
    ap.add_argument('--record', type=int, default=1, help='Record trajectory and totals each step')
    ap.add_argument('--plot', type=int, default=0, help='Plot trajectories')
    ap.add_argument('--plot-dim', type=str, default='xy', choices=['xy','xz','yz'], help='Projection for plotting')
    ap.add_argument('--save-plot', type=str, default=None, help='Path to save trajectory plot (png)')
    ap.add_argument('--save-trj', type=str, default=None, help='Path to save trajectory (npy)')
    ap.add_argument('--plot-bonds', type=int, default=1, help='Overlay molecule bonds on trajectory plot')
    args = ap.parse_args()

    if os.path.abspath(args.molecule) != os.path.abspath(DATA_MOL):
        DATA_MOL = args.molecule

    run_md(
        steps=args.steps, mode=args.mode, dt=args.dt, damp=args.damp, flim=args.flim,
        do_clean=bool(args.do_clean), do_nb=bool(args.do_nb), do_mmff=bool(args.do_mmff),
        print_stats=bool(args.print_stats), print_arrays=bool(args.print_arrays),
        record=bool(args.record), plot=bool(args.plot), plot_dim=args.plot_dim,
        save_plot=args.save_plot, save_trj=args.save_trj,
        plot_bonds=bool(args.plot_bonds)
    )
