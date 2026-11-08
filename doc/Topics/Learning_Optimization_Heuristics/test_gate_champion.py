"""Decisive test (§15.8): does the gated two-timescale momentum (gate_cos/gate_g/vs_ema)
stack with the full champion FIRE configs (autoMass/massMetric), or is it subsumed?
Replays each saved champion via meta_opt's own eval path, then re-evals with the gate on.
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import meta_opt as M

_D = os.path.dirname(os.path.abspath(__file__))
CASES = [
    ('mmff_propandiol', f'{_D}/meta_results/fire_mmff_prop.json'),
    ('uff_propandiol',  f'{_D}/meta_results/fire_uff_prop.json'),
    ('eff_H2O',         f'{_D}/meta_results/fire_eff.json'),
]
GATES = [(0.0, 'off'), (0.25, 'g0.25'), (0.5, 'g0.5')]   # gate_g values (gate_cos=0.7 fixed)

def main():
    print('gate-vs-champion replay: evals to |F|<1e-6 (meta_opt._cost semantics)', flush=True)
    only = sys.argv[1] if len(sys.argv) > 1 else None
    for sysname, jf in CASES:
        if only and sysname != only: continue
        th = json.load(open(jf))['theta']
        if sysname.startswith('eff'):
            os.chdir(os.path.join(M._REPO, 'tests/tEFF'))
            import fire_study_lib as fs
            from pyBall import eFF as eff
            fs.set_atom_params(); fs.load_e('H2O_pairs_fc.xyz', verbosity=0)
            p0 = fs.snapshot()
            for g, tag in GATES:
                eff.setOptParam('gate_cos', 0.7); eff.setOptParam('gate_g', g)
                c = M._eval_eff(fs, eff, p0, 'fire', th)
                print(f'  {sysname:16s} gate_g={tag:6s} -> {c}', flush=True)
            eff.setOptParam('gate_g', 0.0)
        else:
            os.chdir(M._REPO)
            mmff = M._mmff_mod() if hasattr(M, '_mmff_mod') else None
            if mmff is None:
                from pyBall import MMFF as mmff
                bUFF = sysname.startswith('uff'); mol = sysname.split('_', 1)[1]
                kw = dict(xyz_name=f'cpp/common_resources/xyz/{mol}.xyz' if bUFF else f'tests/tMMFF/data/xyz/{mol}.xyz',
                          sElementTypes='cpp/common_resources/ElementTypes.dat', sAtomTypes='cpp/common_resources/AtomTypes.dat',
                          sBondTypes='cpp/common_resources/BondTypes.dat', sAngleTypes='cpp/common_resources/AngleTypes.dat',
                          sDihedralTypes='cpp/common_resources/DihedralTypes.dat',
                          bMMFF=not bUFF, bUFF=bUFF, b141=True, bSimple=True, bConj=True, bCumulene=True, nPBC=(0, 0, 0))
                mmff.setVerbosity(verbosity=0, idebug=0); mmff.init(**kw)
                if bUFF:
                    mmff.getBuffs_UFF(); mmff.setSwitches(PBC=-1, NonBonded=1, NonBondNeighs=-1, SurfAtoms=-1, GridFF=-1)
                    mmff.setSwitchesUFF_NB(NonBonded=1, NonBondNeighs=-1)
                else:
                    mmff.getBuffs()
            DOFs = mmff.DOFs; p0 = DOFs.copy()
            for g, tag in GATES:
                mmff.setOptParam('gate_cos', 0.7); mmff.setOptParam('gate_g', g)
                c = M._eval_mmff(mmff, DOFs, p0, 'fire', th)
                print(f'  {sysname:16s} gate_g={tag:6s} -> {c}', flush=True)
            mmff.setOptParam('gate_g', 0.0)

if __name__ == '__main__':
    main()
