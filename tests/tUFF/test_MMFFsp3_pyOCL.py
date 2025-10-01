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
import matplotlib.pyplot as plt

# Add FireCore root to PYTHONPATH
BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(BASE)
from pyBall.MD_test_utils import *

# infinite line length for numpy print options
np.set_printoptions(linewidth=np.inf)

DATA_MOL = os.path.join(BASE, 'cpp/common_resources/mol/formic_acid.mol2')
#DATA_MOL = os.path.join(BASE, 'cpp/common_resources/mol/CH.mol2')
#DATA_MOL = os.path.join(BASE, 'cpp/common_resources/mol/methanol.mol2')

"""
Examples (omit defaults):
- Minimal MD with diagnostics
  python test_MMFFsp3_pyOCL.py --steps 50 --print-stats 1
- Use rotational integrator (updateAtomsMMFFf4_rot)
  python test_MMFFsp3_pyOCL.py --steps 50 --mode rot
- Use rotational force kernel (getMMFFf4_rot) with rotational integrator
  python test_MMFFsp3_pyOCL.py --steps 50 --mode rot --use-rot-force 1
- Compare basic vs rotational dynamics on formic acid
  python test_MMFFsp3_pyOCL.py --steps 200 --mode basic --monitor 1 --monitor-plot 1
  python test_MMFFsp3_pyOCL.py --steps 200 --mode rot --use-rot-force 1 --monitor 1 --monitor-plot 1
- Forces only (no motion)
  python test_MMFFsp3_pyOCL.py --steps 10 --mode none --print-stats 1
- Record and plot XY trajectories
  python test_MMFFsp3_pyOCL.py --steps 100 --record 1 --plot 1 --plot-dim xy 
- Monitor total invariants like momentum, energy, torque etc.
  python test_MMFFsp3_pyOCL.py --steps 200 --monitor 1 --monitor-plot 1
- Scan energy/force
  python test_MMFFsp3_pyOCL.py --scan 1 --scan-atom 0 --scan-axis x --scan-nsamp 21
  python test_MMFFsp3_pyOCL.py --scan 1
"""

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--rot-dyn',      type=int,   default=1, help='Use rotational dynamics for pi-orbitals (updateAtomsMMFFf4_rot and getMMFFf4_rot)')
    ap.add_argument('--molecule',                 default=DATA_MOL)
    ap.add_argument('--steps',        type=int,   default=100)
    ap.add_argument('--dt',           type=float, default=0.01)
    ap.add_argument('--damp',         type=float, default=1.0)
    ap.add_argument('--flim',         type=float, default=10.0)
    ap.add_argument('--drive-temp',   type=float, default=0.0, help='Langevin thermostat target temperature (0 disables)')
    ap.add_argument('--drive-gamma',  type=float, default=0.0, help='Langevin damping coefficient gamma (0 disables)')
    ap.add_argument('--drive-seed',   type=float, default=0.0, help='Seed offset for Langevin random numbers')
    ap.add_argument('--do-nb',        type=int,   default=0, help='Run getNonBond each step')
    ap.add_argument('--print-stats',  type=int,   default=1, help='Print min/max/norm/NaN counts for buffers at end')
    ap.add_argument('--print-arrays', type=int,   default=1, help='Print full arrays at end')
    ap.add_argument('--record',       type=int,   default=1, help='Record trajectory and totals each step')
    ap.add_argument('--plot',         type=int,   default=0, help='Plot trajectories')
    ap.add_argument('--plot-dim',     type=str,   default='xy', choices=['xy','xz','yz'], help='Projection for plotting')
    ap.add_argument('--save-plot',    type=str,   default="invarants.png", help='Path to save trajectory plot (png)')
    ap.add_argument('--save-trj',     type=str,   default=None, help='Path to save trajectory (npy)')
    ap.add_argument('--save-xyz',     type=str,   default="trj.xyz", help='Path to save trajectory (xyz)')
    ap.add_argument('--plot-bonds',   type=int,   default=1, help='Overlay molecule bonds on trajectory plot')
    ap.add_argument('--subtract-vdw', type=int,   default=0, help='Subtract bonded LJ contributions inside getMMFFf4 (requires non-bonded forces)')
    ap.add_argument('--plot-labels',  type=str,   default='number', choices=['none','number','type'], help='Annotate atoms when plotting trajectories')
    ap.add_argument('--print-params', type=int,   default=1, help='Dump neighbor and bonded parameter arrays before GPU upload')
    ap.add_argument('--monitor',      type=int,   default=1, help='Collect invariant diagnostics each step')
    ap.add_argument('--monitor-props', type=str,  default=','.join(DEFAULT_MONITOR_PROPS), help='Comma separated invariants to track (or all/none)')
    ap.add_argument('--monitor-plot',  type=int,  default=1, help='Plot monitored invariants at the end')
    ap.add_argument('--use-real-mass', type=int,  default=0, help='Use tabulated atomic masses (default: uniform mass=1)')
    ap.add_argument('--save-monitor',  type=str,  default=None, help='Path to save monitor plot (png)')
    ap.add_argument('--save-monitor-data', type=str, default=None, help='Path to save monitor data (npz)')
    ap.add_argument('--scan',         type=int,   default=0, help='Enable energy derivative scan instead of MD')
    ap.add_argument('--scan-atom',    type=int,   default=0, help='Atom index to displace during scan')
    ap.add_argument('--scan-axis',    type=str,   default='x', choices=['x','y','z'], help='Axis to displace (x/y/z)')
    ap.add_argument('--scan-dx',      type=float, default=1e-3, help='Displacement step for scan')
    ap.add_argument('--scan-nsamp',   type=int,   default=1001, help='Number of displacement samples (odd)')
    ap.add_argument('--scan-show',    type=int,   default=1, help='Show scan plot if enabled')
    ap.add_argument('--scan-save',    type=str,   default=None, help='Optional path to save scan plot (png)')
    args = ap.parse_args()

    mol_path = args.molecule

    monitor_props = parse_monitor_props(args.monitor_props)

    axis_map = {'x':0, 'y':1, 'z':2}

    if args.scan:
        mol, mm, md = configure_md( mol_path, args.dt,args.damp, args.flim,  bool(args.subtract_vdw),  args.drive_temp, args.drive_gamma, args.drive_seed, print_params=bool(args.print_params), )
        scan = scan_energy_force(mm, md, atom_index=args.scan_atom, axis=axis_map[args.scan_axis], dx=args.scan_dx, nsamp=args.scan_nsamp, restore=True, do_nb=bool(args.do_nb))
        stats = scan['diff_stats']
        print(f"Scan stats: energy[{stats['energy_min']:.6e}, {stats['energy_max']:.6e}] force[{stats['force_min']:.6e}, {stats['force_max']:.6e}] diff_min={stats['diff_min']:.6e} diff_max={stats['diff_max']:.6e} diff_rms={stats['diff_rms']:.6e}")
        plot_energy_force_scan(scan, axis_label=f"{args.scan_axis}-displacement", show=bool(args.scan_show), save_path=args.scan_save)
        sys.exit(0)

    mol, mm, md = configure_md( mol_path,args.dt, args.damp, args.flim, bool(args.subtract_vdw), args.drive_temp,args.drive_gamma,args.drive_seed, print_params=bool(args.print_params), )
    masses = get_atom_masses(mol, use_real=args.use_real_mass)
    monitor_props, monitor_enabled, totals_hist, monitor_data, trj = init_observers(args.record, args.monitor, monitor_props)

    for _ in range(args.steps):
        if (args.rot_dyn):
            md.run_step_rot  (do_nb=args.do_nb )
        else:
            md.run_step_basic( do_nb=args.do_nb )
        if (args.record or args.monitor):
            collect_diagnostics(md, mm, masses, args.record, args.monitor, monitor_props, totals_hist, monitor_data, trj)

    buf = fetch_arrays(md, mm)
    dump_buffers(buf, print_stats=args.print_stats, print_arrays=args.print_arrays)
    report_conservation(totals_hist, args.steps, args.record)
    mode_label = "rotational" if args.rot_dyn else "basic"
    mol_label = os.path.basename(mol_path)
    plot_title = f"{mol_label} – {mode_label} dynamics"
    trj_arr = None
    if args.record and trj:
        trj_arr = np.stack(trj, axis=0)
        if args.save_trj:
            np.save(args.save_trj, trj_arr)
            print(f"Saved trajectory to {args.save_trj} with shape {trj_arr.shape}")
        want_plot = (args.plot or bool(args.save_plot))
        if want_plot:
            labels = None
            if   args.plot_labels == 'number':  labels = [str(i) for i in range(trj_arr.shape[1])]
            elif args.plot_labels == 'type' and getattr(mol, 'enames', None) is not None:   labels = mol.enames
            fig, _ = plot_trajectories(trj_arr, dim=args.plot_dim, labels=labels, title=f'Trajectories ({args.plot_dim})', save_path=None, show=False)
            if fig is not None and plot_title: fig.suptitle(plot_title, fontsize=12)
            if args.plot_bonds: overlay_bonds(mol, trj_arr[-1], args.plot_dim, labels=labels)
            if args.save_plot:
                fig.savefig(args.save_plot, dpi=150)
                print(f"Saved plot to {args.save_plot}")
            if args.plot: plt.show()
            else:        plt.close(fig)
    elif args.plot or args.save_trj or args.save_plot:   print("No trajectory recorded; skipping plot/trajectory export.")

    if args.save_xyz and trj_arr is not None:
        symbols = mol.enames[:mm.natoms]
        write_xyz_trajectory(args.save_xyz, trj_arr[:, :mm.natoms, :], symbols)
    finalize_monitoring(monitor_enabled, monitor_data, monitor_props, args.save_monitor_data, args.monitor_plot, args.save_monitor)

    apos, aforce = buf['apos'], buf['aforce']
    print(f"Done. Steps={args.steps}  natoms={mm.natoms}  first-atom pos={apos[0,:3]}  force={aforce[0,:3]}")