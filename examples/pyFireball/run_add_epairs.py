#!/usr/bin/env python3

import sys
import os
#sys.path.append(os.path.join(os.path.dirname(__file__), '../../pyBall'))
sys.path.append("../../")
#from pyBall import atomicUtils as au
from pyBall.AtomicSystem import AtomicSystem

#fname='H2O'
#fname='NH3'
#fname='CH2O'
fname='CH2NH'
fname=[
    'H2O',
    'NH3',
    'CH2O',
    'CH2NH',
]

for fname in fname:
    molecule_path =  './molecules/%s.xyz' %fname
    asys = AtomicSystem(fname=molecule_path)
    asys.findBonds()
    asys.add_electron_pairs()
    output_path = './molecules/%s_ep.xyz' %fname
    asys.saveXYZ(output_path)
    print(f"Saved molecule with electron pairs to {output_path}")
