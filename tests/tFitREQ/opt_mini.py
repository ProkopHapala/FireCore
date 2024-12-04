#!/usr/bin/python3 -u

import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

#sys.path.append("/home/niko/work/FIRECORE/FireCore/")
sys.path.append("/home/prokop/git/FireCore-fitREQH")
from pyBall import FitREQ as fit
from pyBall import atomicUtils as au

# ============== Setup
imodel = 2        #  0=LJQ     1=LJQH1     2=LJQH2     3=LJQH1H2
                  #  4=BuckQ   5=BuckQH1   6=BuckQH2   7=BuckQH1H2
                  #  8=MorseQ  9=MorseQH1 10=MorseQH2 11=MorseQH1H2
                  # 12=LJx2Q  13=LJx2QH1  14=LJx2QH2  15=LJx2QH1H2
                  # 16=LJr8Q  17=LJr8QH1  18=LJr8QH2  19=LJr8QH1H2
                  # 20=LJr9Q  21=LJr9QH1  22=LJr9QH2  23=LJr9QH1H2
isampmode   = 2    # do not change it
ialg        = 2         # 0=GD 1=MD 2=GD_BB_short 3=GD_BB_long
#nstep       = 10
nstep       = 1
dt          = 0.01
ErrMax      = 1e-8
bRegularize = True
bClamp      = False
max_step    = 0.01
bEpairs     = True
bAddEpairs  = bEpairs
bOutXYZ     = False
verbosity   = 3    # Added to enable debug printing

# ------ load stuff
#fit.setVerbosity(1)
fit.setVerbosity(verbosity)
fit.loadTypes_new( )     # load atom types
fit.loadTypeSelection_walls( fname="typeSelection.dat" )     # load atom types
nbatch = fit.loadXYZ_new( "input_all.xyz", bAddEpairs, bOutXYZ )     # load reference geometry
ws     = np.genfromtxt( "weights_all.dat" )
fit.setWeights(ws)

# ------ write unoptimized results
Es = fit.getEs( imodel=imodel, isampmode=isampmode, bEpairs=bEpairs )
np.savetxt("firecore0.dat", Es )

# ------ optimize parameters (fit)
Err = fit.run( nstep=nstep, ErrMax=ErrMax, dt=dt, imodel=imodel, isampmode=isampmode, ialg=ialg, bRegularize=bRegularize, bClamp=bClamp, max_step=max_step, bEpairs=bEpairs )

# ------ write optimized results
Es = fit.getEs( imodel=imodel, isampmode=isampmode, bEpairs=bEpairs )
np.savetxt("firecore.dat", Es )
