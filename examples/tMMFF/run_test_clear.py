import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF as mmff


#======== Body

mmff.setVerbosity( verbosity=1, idebug=0 )

nstepMax = 1000; dt=0.05; Fconv=1e-2
#fname = "H2O"
fname = "hydropentacene_cross"

print("\n\n#============= 1st run =============\n\n")

mmff.setTrjName( "opt_run_1.xyz", savePerNsteps=1 );
mmff.init( xyz_name="data/"+fname );
nitr = mmff.run(nstepMax=nstepMax, dt=dt, Fconv=Fconv )
mmff.clear()

print("\n\n#============= 2nd run =============\n\n")

mmff.setTrjName( "opt_run_2.xyz", savePerNsteps=1 );
mmff.init( xyz_name="data/"+fname );
nitr = mmff.run(nstepMax=nstepMax, dt=dt, Fconv=Fconv )


print("\n\n#============= ALL DONE =============\n\n")