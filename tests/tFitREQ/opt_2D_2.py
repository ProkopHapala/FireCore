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

ref_path = "/home/prokop/Desktop/CARBSIS/PEOPLE/Paolo/FitREQ/DFT_2D/"
#name = "H2O-D1_H2O-A1"
#name = "H2O-D1_H2O-A1"
name = "HCOOH-D1_HCOOH-A1"

# ============== Setup
imodel      = 1    
isampmode   = 2    
ialg        = 2         
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

# ============== Setup

#ref_path = "/home/prokop/Desktop/CARBSIS/PEOPLE/Paolo/FitREQ/DFT_2D/"

#ref_dirs = fit.find_all_dirs( ref_path ); #print(ref_dirs)
#frags1, frags2 = fit.extract_fragment_names( base_path=ref_path ); print( "donors:\n", frags1, "\nacceptors:\n", frags2 )

#donors    = ['C4H3NO2-D1', 'C4H5N-D1', 'CH2NH-D1',            'H2O-D1', 'HCONH2-D1', 'HCOOH-D1',             'NH3-D1'] 
#acceptors = ['C4H3NO2-A1', 'C5H5N-A1', 'CH2NH-A1', 'CH2O-A1', 'H2O-A1', 'HCONH2-A1', 'HCOOH-A1', 'HCOOH-A2', 'NH3-A1']

donors    = [
# 'C4H3NO2-D1', 
# 'C4H5N-D1', 
# 'CH2NH-D1',            
'H2O-D1', 
# 'HCONH2-D1', 
# 'HCOOH-D1',             
# 'NH3-D1'
] 
acceptors = [
# 'C4H3NO2-A1', 
# 'C5H5N-A1', 
# 'CH2NH-A1', 
#'CH2O-A1', 
'H2O-A1', 
# 'HCONH2-A1', 
# 'HCOOH-A1', 
# 'HCOOH-A2', 
#'NH3-A1', 
]

ref_dirs = fit.combine_fragments( donors, acceptors )  ;print( "ref_dirs:\n", ref_dirs )

marks    = fit.concatenate_xyz_files( directories=ref_dirs, base_path=ref_path, fname='all.xyz', output_file='all.xyz' )

#fname = "input_single.xyz"
#fname = ref_path +"/"+ name + "/all.xyz"
#fname = ref_path +"/"+ "/concatenated_all.xyz"
fname = 'all.xyz'
#fname = "input_2CH2NH.xyz"

# comments          = fit.read_file_comments(fname) #;print( "comments:\n", comments )
type_names,comments = fit.extract_comments_and_types(fname)
marks, angle_data   = fit.mark_molecule_blocks( comments )
# print( "marks:\n", marks )
# print( "angle_data:\n", angle_data )

# Erefs, x0s = fit.read_xyz_data(fname)  #;print( "x0s:\n", x0s )
# Eplots = fit.slice_and_reshape(Erefs, marks, angle_data)

# fit.plot_Epanels(Eplots, ref_dirs, bColorbar=True)

# plt.show()

# exit()


# ------ load stuff
#fit.setVerbosity(1)
fit.setVerbosity(verbosity, PrintDOFs=1, PrintfDOFs=1, PrintBeforReg=-1, PrintAfterReg=1 )
#fit.setVerbosity(verbosity, PrintDOFs=1, PrintfDOFs=1, PrintBeforReg=-1, PrintAfterReg=1 )
fit.loadTypes( )     # load atom types


if bMorse:
    #fit.loadDOFSelection( fname="dofSelection_Morse.dat" )
    fit.loadDOFSelection( fname="dofSelection_H2O_Morse.dat" )
    #fit.comment_non_matching_lines( fname_in="dofSelection_Morse.dat"); fit.loadDOFSelection()
    #fit.loadDOFSelection( fname="dofSelection_HCOOH_Morse.dat" )
    #fit.loadDOFSelection( fname="dofSelection_HCOOH_Morse.dat" )
else:
    #fit.loadDOFSelection( fname="dofSelection_LJ.dat" )   
    #fit.loadDOFSelection( fname="dofSelection_H2O_LJ.dat" )  
    fit.loadDOFSelection( fname="dofSelection_H2O_LJ.dat" )  
    #fit.loadDOFSelection( fname="dofSelection_CH2NH_LJ.dat" )   
    #fit.loadDOFSelection( fname="dofSelection_HCOOH_LJ.dat" ) 
    #fit.loadDOFSelection( fname="dofSelection_HCOOH_LJ.dat" ) 

#fname = "input_2CH2NH.xyz"
#fit.loadDOFSelection( fname="dofSelection_N2.dat" )          

nbatch = fit.loadXYZ( fname, bAddEpairs, bOutXYZ )     # load reference geometry
#nbatch = fit.loadXYZ( "input_small.xyz", bAddEpairs, bOutXYZ )     # load reference geometry
#nbatch = fit.loadXYZ( "input_single.xyz", bAddEpairs, bOutXYZ )     # load reference geometry
#exit(0)

Erefs, x0s = fit.read_xyz_data(fname)  #;print( "x0s:\n", x0s )
#weights = split_and_weight_curves(Erefs, x0s, n_before_min=4)

EminPlot = np.min(Erefs)*fit.ev2kcal
EminRef = np.min(Erefs)
#EminPlot = -0.5

weights0 = np.ones( len(Erefs) )*0.5


fit.setGlobalParams( kMorse=1.8, Lepairs=0.7 )
if bMorse:
    fit.setup( imodel=2, EvalJ=1, WriteJ=1, Regularize=1 )
    weights0, lens = fit.split_and_weight_curves( Erefs, x0s, n_before_min=100, weight_func=lambda E: fit.exp_weight_func(E,a=1.0, alpha=4.0) )
else:
    fit.setup( imodel=1, EvalJ=1, WriteJ=1, Regularize=1 )
    #fit.setup( imodel=3, EvalJ=1, WriteJ=1, Regularize=1 )
    weights0, lens = fit.split_and_weight_curves( Erefs, x0s, n_before_min=2, weight_func=lambda E: fit.exp_weight_func(E,a=1.0, alpha=4.0) )
# plotEWs( Erefs=Erefs, weights0=weights0, Emin=-1.5 ); plt.title( "Weighting" )
# plt.show(); exit()

fit.setWeights( weights0 )
fit.getBuffs()

#print( "fit.weights ", fit.weights )
#ws     = np.genfromtxt( "weights_all.dat" )
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

# Eplots_ref = fit.slice_and_reshape(Erefs, marks, angle_data)
# Eplots_mod = fit.slice_and_reshape(Es,    marks, angle_data)
# fig = fit.plot_Epanels_diff(Eplots_mod, Eplots_ref, ref_dirs, Emin=EminRef*fit.ev2kcal, bColorbar=True, bKcal=True )
# plt.savefig( "opt_2D.png" )

#plt.show(); exit()
if bMorse:
    #Err = fit.run( iparallel=0, ialg=0, nstep=1000, Fmax=1e-4, dt=0.1, max_step=-1,  bClamp=True )
    Err = fit.run( iparallel=0, ialg=1, nstep=100, Fmax=1e-8, dt=0.5, damping=0.1,   max_step=-1,  bClamp=True )
else:
    #Err = fit.run( iparallel=0, ialg=0, nstep=1000, Fmax=1e-4, dt=0.01, max_step=-1,  bClamp=True )
    Err = fit.run( iparallel=0, ialg=1, nstep=100, Fmax=1e-4, dt=0.1, damping=0.1,   max_step=-1,  bClamp=True )

# ----- Combined hybrid optimization ( start with gradient descent, continue with dynamical descent) )
#Err = fit.run( iparallel=0, ialg=0, nstep=20,  Fmax=1e-2, dt=0.005, max_step=-1,  bClamp=False )
#Err = fit.run( iparallel=0, ialg=1, nstep=100, Fmax=1e-8, dt=0.01, damping=0.1,   max_step=-1,  bClamp=True )

# print( "fit.fDOFmin ", fit.fDOFbounds[:,0] )
# print( "fit.fDOFmax ", fit.fDOFbounds[:,1] )

#E,Es,Fs = fit.getEs( bOmp=False, bDOFtoTypes=False, bEs=True, bFs=False, xyz_name="opt_2D_out.xyz");
E,Es,Fs = fit.getEs( bOmp=False, bDOFtoTypes=False, bEs=True, bFs=False);
fit.plotEWs( Erefs=Erefs, Emodel=Es, weights=fit.weights, Emin=EminPlot );   plt.title( "AFTER OPTIMIZATION" )

# Eplot     = reformat_and_pad_data(Es   , lens)  # Reformat and pad data
# Eplot_ref = reformat_and_pad_data(Erefs, lens)
#plot_data(Eplot, ref_dirs)  # Plot the data

#plot_data_panels(Eplot, Eplot_ref, ref_dirs, bColorbar=True)


Eplots_ref = fit.slice_and_reshape(Erefs, marks, angle_data)
Eplots_mod = fit.slice_and_reshape(Es,    marks, angle_data)

fig = fit.plot_Epanels_diff(Eplots_mod, Eplots_ref, ref_dirs, Emin=EminRef*fit.ev2kcal, bColorbar=True, bKcal=True )
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
