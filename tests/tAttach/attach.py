

import sys
import os

sys.path.append('../../')
from pyBall             import atomicUtils as au
from pyBall.atomicUtils import AtomicSystem

import numpy as np
import matplotlib.pyplot as plt
#from functools import partial

# ============ MAIN
# A - H-N(12-1)    up=(1,10)
# T - H-N(11-3)    up=(12,8)
# backbone  F-C(1-2)    Cl-C(17-9)

'''
# -------- A-T pair
BB = AtomicSystem(fname='backbone.xyz' )
G1 = AtomicSystem(fname='endgroups/adenine.xyz'  )
G2 = AtomicSystem(fname='endgroups/thymine.xyz'  )
BB.attach_group( G1,  1, 12, (1,10),   (1 ,2),  up=(0.,0.,1.),  _0=1  )
BB.attach_group( G2,  3, 11, (12,8),   (17,9),  up=(0.,0.,1.),  _0=1  )
BB.delete_atoms( [1-1,17-1] )
BB.saveXYZ( "BB_A_T.xyz" )
'''

# -------- A-T pair
BB = AtomicSystem(fname='backbone.xyz' )
G1 = AtomicSystem(fname='endgroups/penta_hb3_acceptor2.xyz'  )
G2 = AtomicSystem(fname='endgroups/penta_hb3_donor.xyz'  )
BB.attach_group( G1, 10, 11, (11,8),   (1 ,2),  up=(0.,0.,1.),  _0=1  )
BB.attach_group( G2, 13, 16, (13,8),   (17,9),  up=(0.,0.,1.),  _0=1  )
BB.delete_atoms( [1-1,17-1] )
BB.saveXYZ( "BB_G1_G2.xyz" )


