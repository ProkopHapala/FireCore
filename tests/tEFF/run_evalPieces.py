import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

sys.path.append("../../")
from pyBall import eFF as eff
const_bohr = 0.5291772105638411

#eff.setVerbosity(3)

#outE=eff.eval_mol("Li_eFF"          ,fUnits=const_bohr, outE=True)     ;exit()
#outE=eff.eval_mol("Li_eFF_fixcore"  ,fUnits=const_bohr, outE=True)     ;exit()

#outE=eff.eval_mol("CH4_lmps"        ,fUnits=const_bohr, outE=True)    ;exit()
outE=eff.eval_mol("CH4_lmps_fixcore",fUnits=const_bohr, outE=True)     ;exit()