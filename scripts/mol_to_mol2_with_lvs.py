#!/usr/bin/env python3

import os, sys, argparse
import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(BASE)

from pyBall.AtomicSystem import AtomicSystem


def read_lvs_from_mol2(mol2_path):
    lvs = None
    with open(mol2_path, 'r') as f:
        for line in f:
            s = line.strip()
            if s.startswith('#lvs') or s.startswith('#LVS'):
                parts = s[4:].split()
                nums = []
                for w in parts:
                    try:
                        nums.append(float(w))
                    except Exception:
                        pass
                if len(nums) < 9:
                    raise ValueError(f"Failed to parse 9 floats from lvs line in {mol2_path}: '{s}'")
                lvs = np.array(nums[:9], dtype=float).reshape(3, 3)
                break
    if lvs is None:
        raise ValueError(f"No #lvs line found in {mol2_path}")
    return lvs


def lvs_to_comment(lvs):
    lvs = np.array(lvs, dtype=float).reshape(3, 3)
    a = lvs
    return f"#lvs  {a[0,0]} {a[0,1]} {a[0,2]}     {a[1,0]} {a[1,1]} {a[1,2]}           {a[2,0]} {a[2,1]} {a[2,2]} "


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mol', required=True, help='Input .mol (V2000) file')
    ap.add_argument('--lvs_from_mol2', required=True, help='Reference .mol2 file containing #lvs line')
    ap.add_argument('--out', required=True, help='Output .mol2 file')
    args = ap.parse_args()

    lvs = read_lvs_from_mol2(args.lvs_from_mol2)
    mol = AtomicSystem(fname=args.mol)
    mol.lvec = lvs

    # Basic sanity
    es = list(mol.enames)
    # allow typed names like C_R, C_3, O_2, ... (we want to preserve them)
    bad = [e for e in es if (e.split('_')[0] not in ('H','C','N','O','F','Cl','Br','I','S','P','Si'))]
    if bad:
        raise ValueError(f"Unexpected element symbols in {args.mol}: sample={bad[:10]}")

    mol.save_mol2(args.out, comment=lvs_to_comment(lvs), simple_names=False)

    # Post-check
    with open(args.out, 'r') as f:
        txt = f.read()
    if 'Du' in txt:
        raise ValueError(f"Output mol2 still contains 'Du': {args.out}")
    if '#lvs' not in txt.split('\n', 1)[0]:
        raise ValueError(f"Output mol2 does not start with #lvs header: {args.out}")

    print(f"Wrote {args.out} from {args.mol} with lvs from {args.lvs_from_mol2}")


if __name__ == '__main__':
    main()
