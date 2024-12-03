import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")

from pyBall import MMFF as mmff
#from pyBall import atomicUtils as au
#from pyBall import FunctionSampling as fu


d = 0.5
atoms1=np.array([
    [ -d, 0.0, 0.0,   -1.0 ],
    [ +d, 0.0, 0.0,    1.0 ],
    #[5.0  ,5.0-d,5.0],
    #[5.0  ,5.0+d,5.0],
])

atoms2=np.array([
    [ -d, 5.0, 0.0,   -1.0 ],
    [ +d, 5.0, 0.0,    1.0 ],
    #[5.0  ,5.0-d,5.0],
    #[5.0  ,5.0+d,5.0],
])

ps1 = atoms1[:,0:3].copy(); qs1 = atoms1[:,3].copy()
ps2 = atoms1[:,0:3].copy(); qs2 = atoms2[:,3].copy()

p1,cs1 = mmff.projectMultiPole( ps1, qs1 );
p2,cs2 = mmff.projectMultiPole( ps1, qs1 );

fe2m = mmff.sampleMultipole( ps2, p1, cs1, order=2 )
fe1m = mmff.sampleMultipole( ps1, p2, cs2, order=2 )

fe1 = mmff.sampleCoulombPBC(  ps1, ps2, qs2, nPBC=[0,0,0] )
fe2 = mmff.sampleCoulombPBC(  ps2, ps1, qs1, nPBC=[0,0,0] )