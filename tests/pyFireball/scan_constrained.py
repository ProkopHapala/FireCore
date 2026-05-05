#!/usr/bin/env python3
"""Constrained scan: move H atoms between ribbons while fixing N atoms and relaxing all others.

Scans the H-bond coordinate (position of H along the N...N axis) as a reaction coordinate.
Fixes N atoms in place, relaxes all C and remaining H atoms.
Outputs: XYZ movie, energy profile PNG, per-step geometry plots with force arrows.
"""
import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))
from pyBall import dftb_utils as dftbu
import pyBall.plotUtils as plotUtils
from pyBall.KekuleBackend import build_two_ribbon_cell

DFTB_PATH    = '/home/prokophapala/miniconda3/bin/dftb+'
BASIS_PATH   = '/home/prokophapala/SIMULATIONS/dftbplus/slakos/3ob-3-1/'

def run_scan(width_chains=4, Lx=2.4, a_CC=1.42, L_Hb_eq=2.0,
             scan_range=(1.0, 3.2), n_steps=23,
             do_relax=True, nk_x=16, outdir='.'):
    """Run constrained H-transfer scan."""
    print("=" * 60)
    print(f"Constrained H-transfer scan  w={width_chains}  Lx={Lx}")
    print(f"  L_H range: {scan_range[0]:.2f} – {scan_range[1]:.2f} Å  ({n_steps} steps)")
    print(f"  Relaxed: {do_relax}")
    print("=" * 60)

    os.makedirs(outdir, exist_ok=True)

    # Build equilibrium geometry
    apos, atypes, enames, lvs = build_two_ribbon_cell(width_chains=width_chains, length_cells=1, Lx=Lx, a_CC=a_CC, L_Hb=L_Hb_eq)

    # Identify atoms (generic helper)
    h_scan_idx, n_donor_idx, n_acceptor_idx, fixed_idx = dftbu.identify_hbond_transfer(enames, apos, h_symbol='H', heavy_symbol='N', h_select_axis=1, h_select_mode='abs_min', h_select_value=0.0, verbose=True)

    L_H_vals = np.linspace(scan_range[0], scan_range[1], n_steps)
    tag = f'w{width_chains}_Lx{Lx}_{"relax" if do_relax else "rigid"}'

    p_nd = apos[n_donor_idx]
    p_na = apos[n_acceptor_idx]
    path = np.zeros((1, n_steps, 3), dtype=float)
    path[0, :, :] = dftbu.make_axis_path(p_nd, p_na, L_H_vals)

    if do_relax:
        relax_params = dftbu.default_params.copy()
        relax_params["Optimizer"] = "FIRE{StepSize = 0.1}"
        relax_params["MaxSteps"] = 500
        relax_params["GradElem"] = 1e-3
        params = relax_params
        Temperature = 300
        MixingParameter = 0.02
        MaxScc = 500
        SCCTolerance = 1e-6
    else:
        params = dftbu.default_params
        Temperature = 600
        MixingParameter = 0.1
        MaxScc = 500
        SCCTolerance = 1e-4

    def _plot_step(istep, apos_out, forces, rmeta):
        L_H = rmeta['L_H']
        f_fixed_max = np.max([np.linalg.norm(forces[j]) for j in fixed_idx]) if forces is not None else float('nan')
        pfname = os.path.join(outdir, f'geom_{tag}_step{istep:03d}.png')
        plotUtils.plotGeometryWithForces(apos_out, enames, lvs=lvs, forces=forces, fixed_idx=fixed_idx, highlight={n_donor_idx:{'edgecolor':'darkgreen','linewidth':4,'size':150,'label':f'N{n_donor_idx+1}\\nDONOR'}, n_acceptor_idx:{'edgecolor':'orange','linewidth':4,'size':150,'label':f'N{n_acceptor_idx+1}\\nACCEPT'}, h_scan_idx:{'edgecolor':'magenta','linewidth':4,'size':200,'label':f'H{h_scan_idx+1}\\nSCAN'}}, scan_path=(n_donor_idx, n_acceptor_idx), title=f'Step {istep+1}: L_H={L_H:.3f} Å   E={rmeta["E"]:.4f} eV   |F_N|={f_fixed_max:.3f} eV/Å', fname=pfname)

    results = dftbu.constrained_scan(
        apos0=apos, enames=enames, lvs=lvs,
        moved_idx=[h_scan_idx], path=path,
        fixed_idx=fixed_idx,
        outdir=outdir,
        do_relax=do_relax,
        use_prev_relaxed=True,
        basis_path=BASIS_PATH,
        dftb_exe=DFTB_PATH,
        nk=(nk_x, 1, 1),
        Temperature=Temperature, MixingParameter=MixingParameter, MaxScc=MaxScc, SCCTolerance=SCCTolerance,
        params=params,
        step_prefix='tmp_scan', results_prefix=f'scan_{tag}',
        save_xyz=False,
        key_func=lambda istep, apos_out, forces, E: float(L_H_vals[istep]),
        key_name='L_H',
        plot_step_func=_plot_step,
    )

    if False:
        results  = []
        for i, L_H in enumerate(L_H_vals):
            print(f"\nStep {i+1}/{n_steps}  L_H={L_H:.3f} Å")

            # Move scanned H to target position; keep previous relaxed geometry if available
            apos_step = set_h_position(
                results[-1]['apos'] if (results and do_relax) else apos,
                enames, h_scan_idx, n_donor_idx, n_acceptor_idx, L_H)

            temp_dir = os.path.join(outdir, f'tmp_scan_{i:03d}')
            # Use stable settings for relaxed scans: FIRE optimizer, tight SCC, low mixing
            if do_relax:
                relax_params = dftbu.default_params.copy()
                relax_params["Optimizer"] = "FIRE{StepSize = 0.1}"  # More stable than LBFGS for constraints
                relax_params["MaxSteps"] = 500  # Moderate steps
                relax_params["GradElem"] = 1e-3  # Moderate convergence
                E, apos_out, forces = run_dftb(
                    apos_step, enames, lvs, temp_dir,
                    do_relax=True, fixed_atoms=fixed_idx, nk_x=nk_x,
                    Temperature=300, MixingParameter=0.02, MaxScc=500, SCCTolerance=1e-6,
                    params=relax_params)
            else:
                E, apos_out, forces = run_dftb(
                    apos_step, enames, lvs, temp_dir,
                    do_relax=False, fixed_atoms=fixed_idx, nk_x=nk_x)

            if E is None:
                print(f"  FAILED at step {i+1}")
                continue

            f_fixed_max = np.max([np.linalg.norm(forces[j]) for j in fixed_idx]) if forces is not None else float('nan')
            print(f"  E={E:.4f} eV  |F_N|_max={f_fixed_max:.3f} eV/Å")

            r = {'L_H': L_H, 'E': E, 'apos': apos_out, 'enames': enames,
                 'forces': forces, 'fixed_idx': fixed_idx}
            results.append(r)

            # Per-step geometry plot
            pfname = os.path.join(outdir, f'geom_{tag}_step{i:03d}.png')
            plotUtils.plotGeometryWithForces(apos_out, enames, lvs=lvs, forces=forces, fixed_idx=fixed_idx, highlight={n_donor_idx:{'edgecolor':'darkgreen','linewidth':4,'size':150,'label':f'N{n_donor_idx+1}\\nDONOR'}, n_acceptor_idx:{'edgecolor':'orange','linewidth':4,'size':150,'label':f'N{n_acceptor_idx+1}\\nACCEPT'}, h_scan_idx:{'edgecolor':'magenta','linewidth':4,'size':200,'label':f'H{h_scan_idx+1}\\nSCAN'}}, scan_path=(n_donor_idx, n_acceptor_idx), title=f'Step {i+1}: L_H={L_H:.3f} Å   E={E:.4f} eV   |F_N|={f_fixed_max:.3f} eV/Å', fname=pfname)

    if not results:
        print("No converged steps — nothing to save.")
        return

    # ─── Save XYZ movie ──────────────────────────────────────────────────────
    xyz_fname = os.path.join(outdir, f'scan_{tag}.xyz')
    dftbu.save_xyz_movie(results, xyz_fname, lvs=lvs, key_order=['L_H', 'E'])

    # ─── Energy profile plot ─────────────────────────────────────────────────
    L_H_res = [r['L_H'] for r in results]
    E_res   = [r['E']   for r in results]
    E_rel   = np.array(E_res) - min(E_res)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    axes[0].plot(L_H_res, E_rel, 'o-', color='blue', lw=2, ms=6)
    axes[0].set_xlabel('H position L_H (Å)');  axes[0].set_ylabel('Relative Energy (eV)')
    axes[0].set_title(f'Energy barrier — w{width_chains} Lx={Lx} {"relaxed" if do_relax else "rigid"}')
    axes[0].grid(True, alpha=0.3)
    i_max = np.argmax(E_rel)
    axes[0].annotate(f'Barrier\n{E_rel[i_max]:.3f} eV',
                     xy=(L_H_res[i_max], E_rel[i_max]),
                     xytext=(L_H_res[i_max]+0.2, E_rel[i_max]*0.8),
                     arrowprops=dict(arrowstyle='->', color='red'), color='red')

    if any(r['forces'] is not None for r in results):
        F_fixed = []
        for r in results:
            if r['forces'] is not None:
                F_fixed.append(max(np.linalg.norm(r['forces'][j]) for j in r['fixed_idx']))
            else:
                F_fixed.append(float('nan'))
        axes[1].plot(L_H_res, F_fixed, 's-', color='red', lw=2, ms=6)
        axes[1].set_xlabel('H position L_H (Å)')
        axes[1].set_ylabel('Max residual force on N atoms (eV/Å)')
        axes[1].set_title('Reaction force on fixed N atoms')
        axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    png_fname = os.path.join(outdir, f'scan_{tag}.png')
    plt.savefig(png_fname, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\nSaved energy profile: {png_fname}")
    print("Done!")

# ─── main ────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description='Constrained H-transfer scan')
    parser.add_argument('--width',   type=int,   default=4,   help='Ribbon width in chains')
    parser.add_argument('--Lx',      type=float, default=2.4, help='Lattice constant Lx (Å)')
    parser.add_argument('--n_steps', type=int,   default=23,  help='Number of scan steps')
    parser.add_argument('--L_min',   type=float, default=1.0, help='Min H position (Å)')
    parser.add_argument('--L_max',   type=float, default=3.2, help='Max H position (Å)')
    parser.add_argument('--rigid',   action='store_true',     help='Rigid scan (no relaxation)')
    parser.add_argument('--outdir',  type=str,   default='.',  help='Output directory')
    args = parser.parse_args()

    run_scan(
        width_chains=args.width,
        Lx=args.Lx,
        scan_range=(args.L_min, args.L_max),
        n_steps=args.n_steps,
        do_relax=not args.rigid,
        outdir=args.outdir,
    )

if __name__ == '__main__':
    main()
