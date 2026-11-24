#!/usr/bin/env python3
"""
Reproducible tutorial experiments: eFF relaxation of H2O and H-bonded systems.

Run from tests/tEFF/ (needs pyBall on the path and cpp/common_resources/*.dat):

    python3 tutorial_H2O_relax.py dimer        # 90deg-misoriented dimer -> correct H-bond
    python3 tutorial_H2O_relax.py dimer_bad    # reversed acceptor -> anti-aligned trap (path dependence)
    python3 tutorial_H2O_relax.py forces       # force decomposition: who converged, who jitters
    python3 tutorial_H2O_relax.py water15      # 15-water lattice relaxation + H-bond inventory
    python3 tutorial_H2O_relax.py dtscan       # same start, different dt -> different basins

Figures and relaxed structures go to artifacts/. Numbers printed are quoted in tutorial_eFF_H2O.md.
"""
import sys, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.append('../../')
from pyBall import eFF as eff

FIG = 'artifacts'
os.makedirs(FIG, exist_ok=True)


def load_mol2(fname):
    """Load .mol2 via MM builder, no relaxation yet (nstepMax=0)."""
    eff.setVerbosity(0)
    eff.processMol2(fname, nstepMax=0, dt=0.1, xyz_out=None, fgo_out=None, bOutputs=(0, 0, 0, 0))
    eff.initOpt(dt=0.1, damping=0.1, f_limit=100.0)


def force_blocks():
    """Split the force vector into atom / electron-position / electron-size blocks."""
    eff.bVel = True; eff.getBuffs()
    na3 = 3 * eff.na
    Fa = np.linalg.norm(eff.fDOFs[:na3])
    Fe_pos = np.linalg.norm(eff.fDOFs[na3:na3 + 3 * eff.ne])
    Fe_siz = np.linalg.norm(eff.fDOFs[na3 + 3 * eff.ne:])
    return Fa, Fe_pos, Fe_siz


def plot_mol(ax, title='', elev=20, azim=-60):
    """3D view: O red, H light blue, e-pair centers magenta, H-bonds dashed."""
    ax.scatter(eff.apos[:, 0], eff.apos[:, 1], eff.apos[:, 2],
               c=['r' if z > 4 else '#7fbfff' for z in eff.aPars[:,0]], s=120, edgecolors='k')
    ax.scatter(eff.epos[:, 0], eff.epos[:, 1], eff.epos[:, 2], c='m', s=12, alpha=0.7)
    for i in range(eff.na):                      # covalent O-H sticks
        for j in range(i + 1, eff.na):
            if (eff.aPars[i,0] > 4) != (eff.aPars[j,0] > 4):
                d = np.linalg.norm(eff.apos[i] - eff.apos[j])
                if d < 1.3:
                    ax.plot(*zip(eff.apos[i], eff.apos[j]), 'b-', lw=1.5)
    ax.set_title(title); ax.view_init(elev, azim)
    ax.set_xlabel('x [A]'); ax.set_ylabel('y [A]'); ax.set_zlabel('z [A]')
    set_equal(ax)


def set_equal(ax):
    r = np.vstack([eff.apos, eff.epos]).ptp(axis=0).max() / 2
    c = np.vstack([eff.apos, eff.epos]).mean(axis=0)
    ax.set_xlim(c[0] - r, c[0] + r); ax.set_ylim(c[1] - r, c[1] + r); ax.set_zlim(c[2] - r, c[2] + r)


def lp_alignment():
    """Dimer analysis: identify donor O (whose own H is closest to the other O) vs acceptor O,
    then measure how each acceptor e-pair aligns with the incoming donor-H direction.
    Returns (O...O distance, H...O distance, [(e_idx, cos, r_Oe) sorted by cos desc])."""
    iO = [i for i in range(eff.na) if eff.aPars[i, 0] > 4]
    iH = [i for i in range(eff.na) if eff.aPars[i, 0] < 4]
    O = [eff.apos[i] for i in iO]
    ownH = lambda i: [eff.apos[h] for h in iH
                      if np.argmin([np.linalg.norm(eff.apos[h] - O[j]) for j in range(2)]) == i]
    dmin = [min(np.linalg.norm(h - O[1 - i]) for h in ownH(i)) for i in range(2)]
    iD, iA = (0, 1) if dmin[0] < dmin[1] else (1, 0)
    O_d, O_a = O[iD], O[iA]
    Hb = min(ownH(iD), key=lambda h: np.linalg.norm(h - O_a))      # bridging donor H
    hb = (Hb - O_a) / np.linalg.norm(Hb - O_a)                     # H-bond direction at acceptor
    eA = [i for i in range(eff.ne)
          if np.linalg.norm(eff.epos[i] - O_a) < np.linalg.norm(eff.epos[i] - O_d)]
    out = []
    for i in eA:
        v = eff.epos[i] - O_a
        out.append((i, v @ hb / np.linalg.norm(v), np.linalg.norm(v)))
    out.sort(key=lambda t: -t[1])
    return np.linalg.norm(O_d - O_a), np.linalg.norm(Hb - O_a), out


def run_dimer(fname, tag, nstepMax=60000):
    tag = os.environ.get('TAG', tag)          # e.g. TAG=dimer_90_origZ for A/B runs
    load_mol2(fname)
    E0 = eff.eval()
    ntot, E, Fa, dt, hE, hF = eff.relax_auto(nstepMax=nstepMax, damping=0.3)
    Fa, Fp, Fs = force_blocks()
    dOO, dHO, lps = lp_alignment()
    print(f"\n=== {tag} ===")
    print(f"E: {E0:.2f} -> {E:.3f} eV in {ntot} iters (dt_final={dt:.3g})")
    print(f"|Fa|={Fa:.5f}  |Fepos|={Fp:.4f}  |Fesize|={Fs:.4f}")
    print(f"O...O = {dOO:.3f} A,  H(donor)...O(acc) = {dHO:.3f} A")
    for i, c, r in lps:
        print(f"  e-pair {i:2d} on acceptor: r(O-e)={r:.2f} A, cos(lp, H-bond dir)={c:+.2f}")
    eff.save_xyz(f'{FIG}/{tag}_relaxed.xyz', bAppend=False)
    fig = plt.figure(figsize=(5, 4)); ax = fig.add_subplot(111, projection='3d')
    plot_mol(ax, f'{tag}: E={E:.2f} eV')
    fig.savefig(f'{FIG}/{tag}_relaxed.png', dpi=110, bbox_inches='tight')
    print(f"saved {FIG}/{tag}_relaxed.png/.xyz")
    return ntot, E, Fa


def cmd_dimer():     run_dimer('EMolecules/H2O/dimer_90.mol2',  'dimer_90')
def cmd_dimer_bad(): run_dimer('EMolecules/H2O/dimer_bad.mol2', 'dimer_bad')


def cmd_forces():
    """Show that atoms converge while electron DOFs keep jittering -> why |Fa| is the
    convergence criterion, not total |F|."""
    load_mol2('EMolecules/H2O/dimer_90.mol2')
    eff.setAutoMass(1, 100.0)
    eff.bVel = True; eff.getBuffs()
    oE = np.full(40000, np.nan); oF = np.full(40000, np.nan)
    print("iters      E [eV]     |Fa|     |Fepos|  |Fesize|")
    for i in range(40):
        eff.run(500, dt=0.05, Fconv=1e-12, ialg=2, outE=oE, outF=oF)
        Fa, Fp, Fs = force_blocks()
        if i % 4 == 0 or (i > 30 and i % 2 == 0):
            print(f"{(i + 1) * 500:6d}  {eff.Es[0]:10.3f}  {Fa:8.5f}  {Fp:7.3f}  {Fs:8.4f}")


def cmd_water15():
    load_mol2('EMolecules/H2O/15.mol2')
    E0 = eff.eval()
    ntot, E, Fa, dt, hE, hF = eff.relax_auto(nstepMax=60000, damping=0.3)
    Fa, Fp, Fs = force_blocks()
    print(f"\n15 waters: E {E0:.1f} -> {E:.3f} eV in {ntot} iters")
    print(f"|Fa|={Fa:.5f} |Fepos|={Fp:.4f} |Fesize|={Fs:.4f}  (na={eff.na}, ne={eff.ne})")
    iO = [i for i in range(eff.na) if eff.aPars[i,0] > 4]
    iH = [i for i in range(eff.na) if eff.aPars[i,0] < 4]
    nhb = 0
    for i in iO:
        for h in iH:
            d = np.linalg.norm(eff.apos[h] - eff.apos[i])
            ownO = min(iO, key=lambda j: np.linalg.norm(eff.apos[h] - eff.apos[j]))
            if i != ownO and d < 2.6:
                nhb += 1
    print(f"H-bonds (H...O < 2.6 A, not covalent): {nhb}")
    eff.save_xyz(f'{FIG}/water15_relaxed.xyz', bAppend=False)
    fig = plt.figure(figsize=(10, 4))
    ax = fig.add_subplot(121, projection='3d'); plot_mol(ax, 'top-ish view', elev=80, azim=-90)
    ax = fig.add_subplot(122, projection='3d'); plot_mol(ax, 'side view', elev=10, azim=-70)
    fig.savefig(f'{FIG}/water15_relaxed.png', dpi=110, bbox_inches='tight')
    print(f"saved {FIG}/water15_relaxed.png/.xyz")


def cmd_dtscan():
    """Same start, different dt: shows dt is not just a step size - it picks the basin."""
    print(f"{'dt':>6} {'E [eV]':>10} {'|Fa|':>9} {'O...O':>7}")
    for dt in (0.02, 0.05, 0.1, 0.3):
        load_mol2('EMolecules/H2O/dimer_90.mol2')
        eff.setAutoMass(1, 100.0)
        oE = np.full(20000, np.nan); oF = np.full(20000, np.nan)
        eff.run(20000, dt=dt, Fconv=1e-12, ialg=2, outE=oE, outF=oF)
        Fa, Fp, Fs = force_blocks()
        dOO, _, _ = lp_alignment()
        print(f"{dt:6.2f} {eff.Es[0]:10.3f} {Fa:9.5f} {dOO:7.3f}")


if __name__ == '__main__':
    cmds = {'dimer': cmd_dimer, 'dimer_bad': cmd_dimer_bad, 'forces': cmd_forces,
            'water15': cmd_water15, 'dtscan': cmd_dtscan}
    if len(sys.argv) < 2 or sys.argv[1] not in cmds:
        print(__doc__); sys.exit(1)
    for c in sys.argv[1:]:
        cmds[c]()
