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

sys.path.append("../../")
from pyBall import dftb_utils as dftbu
import pyBall.plotUtils as plotUtils

# Reuse geometry builder from scan_LHb.py
ELEM_MAP     = {'H': 1, 'C': 6, 'N': 7, 'O': 8}
ELEM_MAP_INV = {v: k for k, v in ELEM_MAP.items()}
DFTB_PATH    = '/home/prokophapala/miniconda3/bin/dftb+'
BASIS_PATH   = '/home/prokophapala/SIMULATIONS/dftbplus/slakos/3ob-3-1/'
ELEM_COLORS  = {'H': 'gray', 'C': 'black', 'N': 'blue', 'O': 'red'}
ELEM_SIZES   = {'H': 80,     'C': 120,     'N': 120,    'O': 120}

# ─── geometry builders (copied from scan_LHb.py) ────────────────────────────

def build_ribbon(passivation, width_chains, length_cells, Lx, a_CC=1.42):
    from doc.Topics.Kekule_Topology.GrapheneRibbonBuilder import GrapheneRibbonBuilder
    xa_nom = a_CC * np.cos(np.pi / 6)
    b = GrapheneRibbonBuilder(a_CC=a_CC)
    scale_x = Lx / (2.0 * xa_nom)
    pos2d, elems, bonds = b.build_zigzag_ribbon(
        width_chains=width_chains, length_cells=length_cells,
        passivation=passivation, scale_x=scale_x)
    atypes = np.array([ELEM_MAP[e] for e in elems], dtype=np.int32)
    return np.array(pos2d), atypes, elems

def build_two_ribbon_cell(width_chains=4, length_cells=1, Lx=2.4, a_CC=1.42, L_Hb=2.0, shift_x=0.0):
    pos2d_N,  atypes_N,  elems_N  = build_ribbon('N',  width_chains, length_cells, Lx, a_CC)
    pos2d_NH, atypes_NH, elems_NH = build_ribbon('NH', width_chains, length_cells, Lx, a_CC)
    apos_N  = np.zeros((len(atypes_N),  3));  apos_N[:,  0:2] = pos2d_N
    apos_NH = np.zeros((len(atypes_NH), 3));  apos_NH[:, 0:2] = pos2d_NH
    apos_N[:,  1] -= np.mean(apos_N[:,  1])
    apos_NH[:, 1] -= np.mean(apos_NH[:, 1])
    apos_NH[:, 0] += shift_x * Lx
    y_max_N  = np.max(apos_N[:,  1])
    y_min_NH = np.min(apos_NH[:, 1])
    apos_NH[:, 1] += (y_max_N + L_Hb) - y_min_NH
    y_span_N  = np.max(apos_N[:,  1]) - np.min(apos_N[:,  1])
    y_span_NH = np.max(apos_NH[:, 1]) - np.min(apos_NH[:, 1])
    Ly = y_span_N + y_span_NH + 2 * L_Hb
    apos  = np.vstack([apos_N, apos_NH])
    atypes = np.concatenate([atypes_N, atypes_NH])
    elems  = list(elems_N) + list(elems_NH)
    apos[:, 2]  = 0.0
    apos[:, 1] -= np.mean(apos[:, 1])
    Lz = 20.0
    apos[:, 2] += 0.5 * Lz
    lvs = np.array([[Lx, 0.0, 0.0], [0.0, Ly, 0.0], [0.0, 0.0, Lz]])
    return apos, atypes, elems, lvs

# ─── DFTB+ I/O helpers ──────────────────────────────────────────────────────

def parse_forces(fname, natoms):
    """Parse forces from DFTB+ detailed.out. Returns (natoms,3) array [eV/Ang]."""
    forces = np.zeros((natoms, 3))
    HAU2EVA = 51.422067  # Hartree/Bohr -> eV/Ang
    with open(fname, 'r') as f:
        lines = f.readlines()
    in_forces = False
    idx = 0
    for line in lines:
        if 'Total Forces' in line:
            in_forces = True;  continue
        if in_forces and idx < natoms:
            parts = line.split()
            if len(parts) >= 3:
                try:
                    forces[idx] = [float(parts[0]), float(parts[1]), float(parts[2])]
                    forces[idx] *= HAU2EVA
                    idx += 1
                except ValueError:
                    continue
    return forces

def run_dftb(apos, enames, lvs, temp_dir, do_relax=False, fixed_atoms=None,
             nk_x=16, Temperature=600, MixingParameter=0.1, MaxScc=500, SCCTolerance=1e-4,
             params=None):
    """Run one DFTB+ calculation. Returns (E, apos_out, forces)."""
    if params is None:
        params = dftbu.default_params
    cwd = os.getcwd()
    os.makedirs(temp_dir, exist_ok=True)
    os.chdir(temp_dir)
    try:
        dftbu.makeDFTBjob_pbc(
            enames=enames, apos=apos, lvs=lvs, fname='dftb_in.hsd',
            basis_path=BASIS_PATH,
            nk=(nk_x, 1, 1), k_shift=(0.5, 0.0, 0.0),
            opt=do_relax, params=params,
            Temperature=Temperature, MixingParameter=MixingParameter,
            MaxScc=MaxScc, SCCTolerance=SCCTolerance,
            fixed_atoms=fixed_atoms
        )
        os.system(f'{DFTB_PATH} > OUT')

        # Parse energy
        Estr = os.popen('grep "Total Energy" OUT | tail -1 | cut -b 52-70').read().strip()
        assert Estr, "Could not parse energy from DFTB+ output"
        E = float(Estr)

        # Read relaxed geometry (prefer GenFormat over corrupted XYZ)
        apos_out = apos.copy()
        if do_relax:
            if os.path.exists('geom.out.gen'):
                # Parse GenFormat
                with open('geom.out.gen') as fgen:
                    lines = fgen.readlines()
                n = int(lines[0].split()[0])
                apos_out = np.zeros((n, 3))
                for j in range(n):
                    parts = lines[2+j].split()
                    apos_out[j] = [float(parts[2]), float(parts[3]), float(parts[4])]
            elif os.path.exists('geom.out.xyz'):
                # Fallback to XYZ (may be corrupted in some DFTB+ versions)
                with open('geom.out.xyz') as fxyz:
                    n = int(fxyz.readline())
                    fxyz.readline()
                    apos_out = np.zeros((n, 3))
                    for j in range(n):
                        parts = fxyz.readline().split()
                        if len(parts) >= 4:
                            apos_out[j] = [float(parts[1]), float(parts[2]), float(parts[3])]

        # Parse forces
        forces = None
        if os.path.exists('detailed.out'):
            forces = parse_forces('detailed.out', len(enames))

        return E, apos_out, forces
    finally:
        os.chdir(cwd)

# ─── plotting ────────────────────────────────────────────────────────────────

def plot_step(apos, atypes, lvs, forces, fixed_idx, h_scan_idx, n_donor_idx, n_acceptor_idx, step_label, fname, axes=(0,1)):
    """Plot geometry with bond lengths, highlighted N-H-N path and force arrows."""
    ax1, ax2 = axes
    natoms = len(atypes)

    fig, ax = plt.subplots(figsize=(12, 10))

    # --- scan path line (donor N -> acceptor N) ---
    p_nd = apos[n_donor_idx]
    p_na = apos[n_acceptor_idx]
    ax.plot([p_nd[ax1], p_na[ax1]], [p_nd[ax2], p_na[ax2]], 
            'g--', lw=3, alpha=0.5, zorder=1, label='Scan path (N-N axis)')
    
    # --- atoms ---
    for i, (e, pos) in enumerate(zip(atypes, apos)):
        c = ELEM_COLORS.get(e, 'green')
        s = ELEM_SIZES.get(e, 80)
        
        # Special highlighting for scan atoms
        if i == n_donor_idx:
            ec = 'darkgreen'; lw = 4; s = 150  # Donor N - thick green
            label = f'Donor N{i+1}'
        elif i == n_acceptor_idx:
            ec = 'orange'; lw = 4; s = 150     # Acceptor N - thick orange  
            label = f'Acceptor N{i+1}'
        elif i == h_scan_idx:
            ec = 'magenta'; lw = 4; s = 200    # Scanning H - thick magenta
            label = f'Scanning H{i+1}'
        elif i in fixed_idx:
            ec = 'red'; lw = 3; label = 'Fixed N'  # Other fixed N
        else:
            ec = 'black'; lw = 1; label = None
            
        ax.scatter(pos[ax1], pos[ax2], c=c, s=s, edgecolors=ec, linewidths=lw, zorder=10)
        
        # Labels with special markers for key atoms
        if i == n_donor_idx:
            ax.text(pos[ax1], pos[ax2], f'N{i+1}\nDONOR', color='white',
                    ha='center', va='center', fontsize=6, fontweight='bold', zorder=11)
        elif i == n_acceptor_idx:
            ax.text(pos[ax1], pos[ax2], f'N{i+1}\nACCEPT', color='white',
                    ha='center', va='center', fontsize=6, fontweight='bold', zorder=11)
        elif i == h_scan_idx:
            ax.text(pos[ax1], pos[ax2], f'H{i+1}\nSCAN', color='white',
                    ha='center', va='center', fontsize=6, fontweight='bold', zorder=11)
        else:
            ax.text(pos[ax1], pos[ax2], f'{e}{i+1}', color='white',
                    ha='center', va='center', fontsize=7, zorder=11)

    # Legend
    ax.scatter([], [], c='blue', s=150, edgecolors='darkgreen', linewidths=4, label=f'Donor N{n_donor_idx+1} (y={apos[n_donor_idx,1]:.2f})')
    ax.scatter([], [], c='blue', s=150, edgecolors='orange', linewidths=4, label=f'Acceptor N{n_acceptor_idx+1} (y={apos[n_acceptor_idx,1]:.2f})')
    ax.scatter([], [], c='gray', s=200, edgecolors='magenta', linewidths=4, label=f'Scanning H{h_scan_idx+1} (y={apos[h_scan_idx,1]:.2f})')

    # --- bonds (within cutoff 1.8 Å, periodic) ---
    bond_dist = 2.0
    for i in range(natoms):
        for j in range(i+1, natoms):
            for six in range(-1, 2):
                for siy in range(-1, 2):
                    shift = np.array([six, siy, 0]) @ lvs
                    d = np.linalg.norm(apos[i] - (apos[j] + shift))
                    if d < bond_dist:
                        pi = apos[i];  pj = apos[j] + shift
                        ax.plot([pi[ax1], pj[ax1]], [pi[ax2], pj[ax2]], 'k-', lw=2, zorder=5)
                        mx, my = 0.5*(pi[ax1]+pj[ax1]), 0.5*(pi[ax2]+pj[ax2])
                        ax.text(mx, my, f'{d:.2f}', color='darkred', fontsize=7,
                                ha='center', va='center',
                                bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.7), zorder=12)

    # --- force arrows ---
    if forces is not None:
        fmax = np.max(np.linalg.norm(forces, axis=1))
        scale = 1.5 / max(fmax, 0.1)  # 1.5 Å per max-force
        for i, (pos, frc) in enumerate(zip(apos, forces)):
            fmag = np.linalg.norm(frc)
            if fmag < 0.005: continue
            ax.annotate('', xy=(pos[ax1] + frc[ax1]*scale, pos[ax2] + frc[ax2]*scale),
                        xytext=(pos[ax1], pos[ax2]),
                        arrowprops=dict(arrowstyle='->', color='red', lw=1.5))
        ax.scatter([], [], c='red', marker=r'$\rightarrow$', s=100,
                   label=f'Force (max={fmax:.3f} eV/Å)')

    # --- cell box ---
    corners = np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,0]]) @ lvs
    ax.plot(corners[:, ax1], corners[:, ax2], 'b--', lw=1.5, alpha=0.6, label='Unit cell')

    ax.set_aspect('equal')
    ax.set_xlabel(f"{'xyz'[ax1]} (Å)");  ax.set_ylabel(f"{'xyz'[ax2]} (Å)")
    ax.set_title(step_label)
    ax.legend(loc='upper right', fontsize=8)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(fname, dpi=120, bbox_inches='tight')
    plt.close()
    print(f"  Saved plot: {fname}")

def save_xyz_movie(results, fname, lvs):
    """Save XYZ movie; comment line includes E, L_H, forces on fixed atoms."""
    with open(fname, 'w') as f:
        for r in results:
            apos   = r['apos'];    enames = r['enames']
            E      = r['E'];       L_H    = r['L_H']
            forces = r.get('forces')
            f_fixed_str = ''
            if forces is not None and r.get('fixed_idx'):
                ff = [np.linalg.norm(forces[i]) for i in r['fixed_idx']]
                f_fixed_str = f' F_fixed={np.mean(ff):.3f}'
            f.write(f"{len(enames)}\n")
            f.write(f"L_H={L_H:.4f} E={E:.6f} eV{f_fixed_str} "
                    f"Lx={lvs[0,0]:.2f} Ly={lvs[1,1]:.2f} Lz={lvs[2,2]:.2f}\n")
            for e, pos in zip(enames, apos):
                f.write(f"{e} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}\n")
    print(f"Saved XYZ movie: {fname}")

# ─── scan logic ──────────────────────────────────────────────────────────────

def set_h_position(apos, enames, h_scan_idx, n_donor_idx, n_acceptor_idx, L_H):
    """Move the scanned H atom to distance L_H from n_donor along the N...N axis."""
    apos_new = apos.copy()
    p_nd = apos[n_donor_idx]
    p_na = apos[n_acceptor_idx]
    axis = p_na - p_nd
    axis_len = np.linalg.norm(axis)
    axis_hat = axis / axis_len
    apos_new[h_scan_idx] = p_nd + axis_hat * L_H
    return apos_new

def identify_scan_atoms(apos, enames, verbose=True):
    """Find H and N atom indices and identify the H-bond pair to scan.

    Returns: h_scan_idx, n_donor_idx, n_acceptor_idx, fixed_idx
      - h_scan_idx   : index of the H to move (H on NH ribbon, highest y)
      - n_donor_idx  : N it currently belongs to (NH ribbon side)
      - n_acceptor_idx: N it will transfer to (N ribbon side)
      - fixed_idx    : all N atom indices (to fix during relaxation)
    """
    h_idx = [i for i, e in enumerate(enames) if e == 'H']
    n_idx = [i for i, e in enumerate(enames) if e == 'N']

    if verbose:
        print(f"  H atoms: {h_idx}  positions y={[f'{apos[i,1]:.2f}' for i in h_idx]}")
        print(f"  N atoms: {n_idx}  positions y={[f'{apos[i,1]:.2f}' for i in n_idx]}")

    # H to scan: H closest to y=0 (at the synaptic junction between ribbons)
    h_scan_idx = h_idx[np.argmin([abs(apos[i, 1]) for i in h_idx])]
    
    # Find the two closest N atoms to this H
    # Closest N = donor (covalent N-H bond)
    # Second closest N = acceptor (hydrogen bond H...N)
    dists_n_to_h = [(n, np.linalg.norm(apos[h_scan_idx] - apos[n])) for n in n_idx]
    dists_n_to_h.sort(key=lambda x: x[1])  # Sort by distance
    
    n_donor_idx = dists_n_to_h[0][0]      # Closest N
    n_acceptor_idx = dists_n_to_h[1][0]   # Second closest N
    d_donor = dists_n_to_h[0][1]
    d_acceptor = dists_n_to_h[1][1]

    if verbose:
        print(f"  Scanning H: idx={h_scan_idx}, y={apos[h_scan_idx,1]:.3f}")
        print(f"  Donor N:    idx={n_donor_idx},  y={apos[n_donor_idx,1]:.2f}, d_HN={d_donor:.2f}Å")
        print(f"  Acceptor N: idx={n_acceptor_idx}, y={apos[n_acceptor_idx,1]:.2f}, d_HN={d_acceptor:.2f}Å")
        dNN = np.linalg.norm(apos[n_donor_idx] - apos[n_acceptor_idx])
        print(f"  N-N distance: {dNN:.3f} Å")
        print(f"  H-bond: N{ n_donor_idx+1 }-H{ h_scan_idx+1 }...N{ n_acceptor_idx+1 }  -->  N{ n_donor_idx+1 }...H{ h_scan_idx+1 }-N{ n_acceptor_idx+1 }")

    # Atoms to fix: all N atoms + the scanning H atom
    fixed_idx = n_idx + [h_scan_idx]
    return h_scan_idx, n_donor_idx, n_acceptor_idx, fixed_idx

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
    apos, atypes, enames, lvs = build_two_ribbon_cell(
        width_chains=width_chains, length_cells=1, Lx=Lx, a_CC=a_CC, L_Hb=L_Hb_eq)

    # Identify atoms
    h_scan_idx, n_donor_idx, n_acceptor_idx, fixed_idx = identify_scan_atoms(apos, enames)

    L_H_vals = np.linspace(scan_range[0], scan_range[1], n_steps)
    results  = []
    tag = f'w{width_chains}_Lx{Lx}_{"relax" if do_relax else "rigid"}'

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
        plot_step(apos_out, enames, lvs, forces, fixed_idx, h_scan_idx, n_donor_idx, n_acceptor_idx,
                  f'Step {i+1}: L_H={L_H:.3f} Å   E={E:.4f} eV   |F_N|={f_fixed_max:.3f} eV/Å',
                  pfname)

    if not results:
        print("No converged steps — nothing to save.")
        return

    # ─── Save XYZ movie ──────────────────────────────────────────────────────
    xyz_fname = os.path.join(outdir, f'scan_{tag}.xyz')
    save_xyz_movie(results, xyz_fname, lvs)

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
