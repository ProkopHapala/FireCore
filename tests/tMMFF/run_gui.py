import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF as mmff
from pyBall import MolGUI as gui

#======== Body

mmff.setVerbosity( verbosity=1, idebug=0 )
#mmff.init( xyz_name="data/pyridine", surf_name="data/NaCl_sum-center" )                             # all
W=mmff.init( xyz_name="data/pyridine", surf_name="data/NaCl_sym-center", bMMFF=False  )              # without MMFF

#mmff.init( xyz_name="data/pyridine", surf_name="data/NaCl_sym-center", bMMFF=False, gridStep=-1 )  # without gridFF
mmff.getBuffs()
mmff.eval()

gui.init(W)
gui.run( 1000000 )