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
verbosity   = 0    # Added to enable debug printing


# ============== functions

def check_array_difference(arr1, arr2, name, max_error=1e-8, err_message="arrays differs" ):
    dmax = (arr1-arr2).max()
    print(f"{name} dmax={dmax}")
    if not np.allclose(arr1, arr2, atol=max_error):
        print(f"{name} arrays differ:")
        for i, (v1, v2) in enumerate(zip(arr1, arr2)):
            if not np.isclose(v1, v2): print(f"{i}\t{v1}\t{v2}")
        assert False, f"{name} "+err_message

def test_getEs_openmp():
    E1, Es1, Fs1 = fit.getEs(imodel=imodel, bOmp=False, bEs=True, bFs=True)
    E2, Es2, Fs2 = fit.getEs(imodel=imodel, bOmp=True,  bEs=True, bFs=True)
    check_array_difference(Es1, Es2, "Es w/o OpenMP")
    check_array_difference(Fs1, Fs2, "Fs w/o OpenMP")
    print( "test_getEs_openmp() E1,E2", E1, E2 )
    assert np.isclose(E1, E2), f"E values differ: E1={E1}, E2={E2}"
    print( "test_getEs_openmp() passed " )

def read_xyz_data(fname="input_all.xyz"):
    """Read XYZ file and extract Etot and x0 values from comment lines"""
    #print("read_xyz_data()\n")
    #print("Reading XYZ file:", fname)
    Etots = []
    x0s = []
    with open(fname, 'r') as f:
        while True:
            line = f.readline()
            #print(line)
            if not line: break
            if line.startswith('# n0'):
                #print(line)
                # Parse line like "# n0 5 Etot .70501356708840164618 x0 1.40"
                parts = line.split()
                Etot  = float(parts[4])
                x0    = float(parts[6])
                Etots.append(Etot)
                x0s.append(x0)
            # Skip the rest of the xyz structure
            #natoms = int(line) if line[0].isdigit() else 0
            #for _ in range(natoms):
            #    f.readline()
    return np.array(Etots), np.array(x0s)

def split_and_weight_curves(Etots, x0s, n_before_min=4):
    """
    Split energy curves based on x0 discontinuities and assign weights.
    
    Args:
        Etots: numpy array of total energies
        x0s: numpy array of x0 values (monotonic within each curve)
        n_before_min: number of points before minimum to keep with positive weight
    
    Returns:
        weights: numpy array of weights (0.0 or 1.0)
    """
    weights = np.zeros_like(Etots)
    
    # Find where x0 values reset (non-monotonic changes)
    dx0 = np.diff(x0s)
    curve_starts = np.where(dx0 < 0)[0] + 1
    
    # Add start and end indices to process all segments
    all_splits = np.concatenate(([0], curve_starts, [len(x0s)]))
    
    # Process each curve segment
    for start, end in zip(all_splits[:-1], all_splits[1:]):
        segment = Etots[start:end]
        if len(segment) == 0:
            continue
            
        # Find minimum in this segment
        min_idx = np.argmin(segment) + start
        icut = min_idx-n_before_min
        weight_start = max(icut, start)
        weights[weight_start:end] = 1.0
    
    return weights

def genWeights(Erefs, Ecut ):
    mask = Erefs<Ecut
    weights = np.zeros( len(Erefs) )
    weights[mask] = 1.0
    return weights

def plotWeights(Erefs, weights):
    plt.plot( Erefs  ,'.-', lw=0.5, ms=1.0, label="E_ref")
    plt.plot( weights, lw=1.0, label="weights")
    plt.legend()
    plt.xlabel("#sample(conf)")
    plt.ylabel("E [kcal/mol]")


def plotDOFscans( iDOFs, xs, label ):
    plt.figure()
    for iDOF in iDOFs:
        y = fit.DOFs[iDOF]    # store backup value of this DOF
        Es,Fs = fit.scanParam( iDOF, xs, imodel=imodel )   # do 1D scan
        #print( "iDOF", iDOF, DOFnames[iDOF], "Es", Es )
        plt.plot(xs,Es, '-', label=DOFnames[iDOF] )       # plot 1D scan
        fit.DOFs[iDOF] = y    # restore
    plt.legend()
    plt.xlabel("DOF value")
    plt.ylabel("E [kcal/mol]")    
    plt.title( label )
    plt.grid()
plt.show()

# ============== Setup

# ------ load stuff
#fit.setVerbosity(1)
fit.setVerbosity(verbosity)
fit.loadTypes( )     # load atom types

fname = "input_all.xyz"
fit.loadTypeSelection( fname="typeSelection.dat" )     # load atom types
nbatch = fit.loadXYZ( fname, bAddEpairs, bOutXYZ )     # load reference geometry
#nbatch = fit.loadXYZ( "input_small.xyz", bAddEpairs, bOutXYZ )     # load reference geometry
#nbatch = fit.loadXYZ( "input_single.xyz", bAddEpairs, bOutXYZ )     # load reference geometry

fit.getBuffs()

#ws     = np.genfromtxt( "weights_all.dat" )
ev2kcal = 23.060548
#Erefs   = fit.export_Erefs()*ev2kcal  #;print( "Erefs:\n", Erefs )
#weights = genWeights( Erefs, Ecut=-2.0 )
#weights = split_and_weight_curves( Erefs, n_before_min=4, jump_threshold=1.0 )
Etots, x0s = read_xyz_data(fname)  #;print( "x0s:\n", x0s )
weights = split_and_weight_curves(Etots, x0s, n_before_min=4)
#plotWeights( Etots, weights ); 
#plt.plot(x0s)
#plt.show(); exit()
fit.setWeights( weights )

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
# plotDOFscans( [0,1,2,3,4], np.linspace(  -1.0,  0.0,  30 ), label="Q Epairs"  )
# plotDOFscans( [5,6,7,8,9], np.linspace(  0.99, 0.0,  30 ), label="H2  X=O,N" )
# plotDOFscans( [10,11]    , np.linspace(  0.0,  0.99, 30 ), label="H2  H-"   )
# plt.show()


# ------ write unoptimized results
#def getEs(imodel=0, Es=None, Fs=None, bOmp=False, bDOFtoTypes=False, bEs=True, bFs=False ):

# # ------ optimize parameters (fit)
#fit.setSwitches(EvalJ=0, WriteJ=0, CheckRepulsion=0, Regularize=0, Epairs=0)
nstep = 5000
fit.setSwitches(EvalJ=1, WriteJ=1  ); Err = fit.run( nstep=nstep, Fmax=1e-300, dt=0.0, imodel=imodel, iparallel=0, ialg=0, bClamp=bClamp )
fit.setSwitches(EvalJ=1, WriteJ=-1 ); Err = fit.run( nstep=nstep, Fmax=1e-300, dt=0.0, imodel=imodel, iparallel=0, ialg=0, bClamp=bClamp )

# fit.setSwitches(EvalJ=1, WriteJ=1 )
# fit.setSwitches(EvalJ=1, WriteJ=-1 )
# nstep = 5000
# Err = fit.run( nstep=nstep, Fmax=1e-300, dt=0.0, imodel=imodel, iparallel=0, ialg=0,  bClamp=bClamp )
# Err = fit.run( nstep=nstep, Fmax=1e-300, dt=0.0, imodel=imodel, iparallel=1, ialg=0,  bClamp=bClamp )
# Err = fit.run( nstep=nstep, Fmax=1e-300, dt=0.0, imodel=imodel, iparallel=2, ialg=0,  bClamp=bClamp )

# # ------ write optimized results
# Es = fit.getEs( imodel=imodel, isampmode=isampmode, bEpairs=bEpairs )
# np.savetxt("firecore.dat", Es )
