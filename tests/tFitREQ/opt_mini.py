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

fit.plt = plt

np.set_printoptions(linewidth=200)

# ============== Setup
imodel = 1        #  0=LJQ     1=LJQH1     2=LJQH2     3=LJQH1H2
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
verbosity   = 2    # Added to enable debug printing

#bMorse = False   # Lenard-Jones
bMorse = True   # Morse


# ============== functions

# ============== Setup

# ------ load stuff
#fit.setVerbosity(1)
fit.setVerbosity(verbosity, PrintDOFs=1, PrintfDOFs=1, PrintBeforReg=-1, PrintAfterReg=1 )
#fit.setVerbosity(verbosity, PrintDOFs=1, PrintfDOFs=1, PrintBeforReg=-1, PrintAfterReg=1 )
fit.loadTypes( )     # load atom types



#fname = "input_single.xyz"
fname = "input_all.xyz"
#fname = "input_2CH2NH.xyz"

if bMorse:
    fit.loadDOFSelection( fname="dofSelection_Morse.dat" )
else:
    fit.loadDOFSelection( fname="dofSelection_LJ.dat" )   

#fname = "input_2CH2NH.xyz"
#fit.loadDOFSelection( fname="dofSelection_N2.dat" )          

nbatch = fit.loadXYZ( fname, bAddEpairs, bOutXYZ )     # load reference geometry
#nbatch = fit.loadXYZ( "input_small.xyz", bAddEpairs, bOutXYZ )     # load reference geometry
#nbatch = fit.loadXYZ( "input_single.xyz", bAddEpairs, bOutXYZ )     # load reference geometry

Erefs, x0s = fit.read_xyz_data(fname)  #;print( "x0s:\n", x0s )
#weights = split_and_weight_curves(Erefs, x0s, n_before_min=4)

EminPlot = np.min(Erefs)
EminPlot = -0.5

weights0 = np.ones( len(Erefs) )*0.5

if bMorse:
    fit.setup( imodel=2, EvalJ=1, WriteJ=1, Regularize=1 )
    weights0 = fit.split_and_weight_curves( Erefs, x0s, n_before_min=100, weight_func=lambda E: fit.exp_weight_func(E,a=0.5, alpha=4.0) )
else:
    fit.setup( imodel=1, EvalJ=1, WriteJ=1, Regularize=1 )
    weights0 = fit.split_and_weight_curves( Erefs, x0s, n_before_min=2, weight_func=lambda E: fit.exp_weight_func(E,a=0.5, alpha=4.0) )
# plotEWs( Erefs=Erefs, weights0=weights0, Emin=-1.5 ); plt.title( "Weighting" )
# plt.show(); exit()

fit.setWeights( weights0 )
fit.getBuffs()

#print( "fit.weights ", fit.weights )
#ws     = np.genfromtxt( "weights_all.dat" )
ev2kcal = 23.060548
#Erefs   = fit.export_Erefs()*ev2kcal  #;print( "Erefs:\n", Erefs )
#weights = genWeights( Erefs, Ecut=-2.0 )
#weights = split_and_weight_curves( Erefs, n_before_min=4, jump_threshold=1.0 )
#plotWeights( Erefs, weights ); 
#plt.plot(x0s)
#plt.show(); exit()

#fit.setFilter( EmodelCutStart=0.0, EmodelCut=0.5, iWeightModel=2, PrintOverRepulsive=1, DiscardOverRepulsive=1, SaveOverRepulsive=-1, ListOverRepulsive=1 )
#fit.setFilter( EmodelCutStart=0.0, EmodelCut=0.5, iWeightModel=2, PrintOverRepulsive=-1, DiscardOverRepulsive=1, SaveOverRepulsive=1, ListOverRepulsive=-1 )
fit.setFilter( EmodelCutStart=0.0, EmodelCut=0.5, PrintOverRepulsive=-1, DiscardOverRepulsive=-1, SaveOverRepulsive=-1, ListOverRepulsive=-1 )
#fit.setFilter( EmodelCutStart=0.0, EmodelCut=0.5, iWeightModel=2, PrintOverRepulsive=-1, DiscardOverRepulsive=1, SaveOverRepulsive=1, ListOverRepulsive=-1 )
#fit.setFilter( EmodelCutStart=0.0, EmodelCut=0.5, PrintOverRepulsive=-1, DiscardOverRepulsive=-1, SaveOverRepulsive=-1, ListOverRepulsive=-1 )




E,Es,Fs = fit.getEs( bOmp=False, bDOFtoTypes=False, bEs=True, bFs=False )
fit.plotEWs( Erefs=Erefs, Emodel=Es, weights=fit.weights, weights0=weights0,  Emin=EminPlot ); plt.title( "BEFORE OPTIMIZATION" )
#plt.show(); exit()

if bMorse:
    #Err = fit.run( iparallel=0, ialg=0, nstep=1000, Fmax=1e-4, dt=0.1, max_step=-1,  bClamp=True )
    Err = fit.run( iparallel=0, ialg=1, nstep=1000, Fmax=1e-4, dt=0.5, damping=0.1,   max_step=-1,  bClamp=True )
else:
    #Err = fit.run( iparallel=0, ialg=0, nstep=1000, Fmax=1e-4, dt=0.01, max_step=-1,  bClamp=True )
    Err = fit.run( iparallel=0, ialg=1, nstep=1000, Fmax=1e-4, dt=0.1, damping=0.1,   max_step=-1,  bClamp=True )

# ----- Combined hybrid optimization ( start with gradient descent, continue with dynamical descent) )
#Err = fit.run( iparallel=0, ialg=0, nstep=20,  Fmax=1e-2, dt=0.005, max_step=-1,  bClamp=False )
#Err = fit.run( iparallel=0, ialg=1, nstep=100, Fmax=1e-8, dt=0.01, damping=0.1,   max_step=-1,  bClamp=True )

print( "fit.fDOFmin ", fit.fDOFbounds[:,0] )
print( "fit.fDOFmax ", fit.fDOFbounds[:,1] )

E,Es,Fs = fit.getEs( bOmp=False, bDOFtoTypes=False, bEs=True, bFs=False );
fit.plotEWs( Erefs=Erefs, Emodel=Es, weights=fit.weights, Emin=EminPlot );   plt.title( "AFTER OPTIMIZATION" )


plt.show(); # exit()



#test_getEs_openmp()

#print( "fit.nDOFs ", fit.nDOFs )
DOFnames = [
"E_N3.Q", # 0
"E_NR.Q", # 1
"E_N2.Q", # 2 
"E_O3.Q", # 3
"E_O2.Q", # 4
"N_3.H",  # 5
"N_R.H",  # 6
"N_2.H",  # 7
"O_3.H",  # 8
"O_2.H",  # 9
"H_N.H",  # 10
"H_O.H"   # 11
]

# ---- Plot 2D parameter scan
# iDOFx = 11 ; xs = np.linspace(   0.0,  1.0,  30 ) # "H_O.H"  # 11
# #iDOFy = 8  ; ys = np.linspace(  -1.0,  0.0,  30 ) # "O_3.H"  # 8
# iDOFy = 9  ; ys = np.linspace(  -1.0,  0.0,  30 ) # "O_3.H"  # 9
# Es,Fx,Fy = fit.scanParam2D( iDOFx, iDOFy, xs, ys, imodel=2, bRegularize=False)
# extent = [xs[0],xs[-1],ys[0],ys[-1]]
# plt.imshow(Es, origin='lower', extent=extent)
# plt.colorbar()
# plt.xlabel(DOFnames[iDOFx])
# plt.ylabel(DOFnames[iDOFy])
# plt.show()

# ------ Plot 1D parameter scans

#fit.plotDOFscans( [0,1,2,3,4], np.linspace(  -1.0,  0.0,  100 ), DOFnames, label="Q Epairs"  )
#fit.plotDOFscans( [5,6,7,8,9], np.linspace(   1.0,  0.0,  100 ), DOFnames, label="H2  X=O,N" )
#fit.plotDOFscans( [10,11]    , np.linspace(   0.0,  1.0,  100 ), DOFnames, label="H2  H-"   )
#fit.plotDOFscans( [10,11]    , np.linspace(   0.0+1e-6,  1.0-1e-6,  100 ), DOFnames, label="H2  H-" , bFs=True , bEvalSamples=False  )
#plt.show()


# ------ write unoptimized results
#def getEs(imodel=0, Es=None, Fs=None, bOmp=False, bDOFtoTypes=False, bEs=True, bFs=False ):

# # ------ optimize parameters (fit)
#fit.setSwitches(EvalJ=0, WriteJ=0, CheckRepulsion=0, Regularize=0, Epairs=0)
# nstep = 5000
# fit.setSwitches(EvalJ=1, WriteJ=1  ); Err = fit.run( nstep=nstep, Fmax=1e-300, dt=0.0, imodel=imodel, iparallel=0, ialg=0, bClamp=bClamp )
# fit.setSwitches(EvalJ=1, WriteJ=-1 ); Err = fit.run( nstep=nstep, Fmax=1e-300, dt=0.0, imodel=imodel, iparallel=0, ialg=0, bClamp=bClamp )






# ------ Test diffirent write-J options with different parallelization
# nstep = 5000
# fit.setSwitches(EvalJ=-1, WriteJ=-1 )
# Err = fit.run( nstep=nstep, Fmax=1e-300, dt=0.0, imodel=imodel, iparallel=0, ialg=0,  bClamp=bClamp )
# Err = fit.run( nstep=nstep, Fmax=1e-300, dt=0.0, imodel=imodel, iparallel=1, ialg=0,  bClamp=bClamp )
# Err = fit.run( nstep=nstep, Fmax=1e-300, dt=0.0, imodel=imodel, iparallel=2, ialg=0,  bClamp=bClamp )
# fit.setSwitches(EvalJ=1, WriteJ=1 )
# Err = fit.run( nstep=nstep, Fmax=1e-300, dt=0.0, imodel=imodel, iparallel=0, ialg=0,  bClamp=bClamp )
# Err = fit.run( nstep=nstep, Fmax=1e-300, dt=0.0, imodel=imodel, iparallel=1, ialg=0,  bClamp=bClamp )
# Err = fit.run( nstep=nstep, Fmax=1e-300, dt=0.0, imodel=imodel, iparallel=2, ialg=0,  bClamp=bClamp )
# fit.setSwitches(EvalJ=1, WriteJ=-1 )
# Err = fit.run( nstep=nstep, Fmax=1e-300, dt=0.0, imodel=imodel, iparallel=0, ialg=0,  bClamp=bClamp )
# Err = fit.run( nstep=nstep, Fmax=1e-300, dt=0.0, imodel=imodel, iparallel=1, ialg=0,  bClamp=bClamp )
# Err = fit.run( nstep=nstep, Fmax=1e-300, dt=0.0, imodel=imodel, iparallel=2, ialg=0,  bClamp=bClamp )

# # ------ write optimized results
# Es = fit.getEs( imodel=imodel, isampmode=isampmode, bEpairs=bEpairs )
# np.savetxt("firecore.dat", Es )

#plot_dofs_fdofs(output_file="OUT", figsize=(12, 8))
