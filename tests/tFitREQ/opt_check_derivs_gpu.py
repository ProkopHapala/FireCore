#!/usr/bin/python3 -u

import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

#sys.path.append("/home/niko/work/FIRECORE/FireCore/")
#sys.path.append("/home/prokop/git/FireCore-fitREQH")

sys.path.append("../../")
from pyBall import FitREQ as fit
from pyBall.OCL.NonBondFitting import FittingDriver, extract_macro_block
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
#fname = 'all.xyz'
#fname = "input_2CH2NH.xyz"
#fname="H2O_1D.xyz"
fname="H2O_single.xyz"
#fname="just_Epair_2x2.xyz"
#fname="just_Epair_1x1_ee.xyz"
#fname="just_Epair_1x1_eh.xyz"

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
fit.setVerbosity(verbosity, PrintDOFs=1, PrintfDOFs=1, PrintBeforReg=-1, PrintAfterReg=-1 )
#fit.setVerbosity(verbosity, PrintDOFs=1, PrintfDOFs=1, PrintBeforReg=-1, PrintAfterReg=1 )
fit.loadTypes( )     # load atom types

dof_fname = "dofSelection.dat"
if bMorse:
    #dof_fname="dofSelection_Morse.dat"
    dof_fname = "dofSelection_H2O_Morse.dat" 
else:
    #dof_fname="dofSelection_LJ.dat" 
    #dof_fname="dofSelection_H2O_LJ.dat" 
    dof_fname="dofSelection_H2O_LJSR.dat" 
fit.loadDOFSelection( dof_fname)
dof_names, dof_specs = fit.loadDOFnames( dof_fname, return_specs=True )

nbatch = fit.loadXYZ( fname, bAddEpairs, bOutXYZ )     # load reference geometry
#nbatch = fit.loadXYZ( fname, bAddEpairs=False, bOutXYZ=bOutXYZ )     # load reference geometry

Erefs, x0s = fit.read_xyz_data(fname)  #;print( "x0s:\n", x0s )
#weights = split_and_weight_curves(Erefs, x0s, n_before_min=4)

EminPlot = np.min(Erefs)*fit.ev2kcal
EminRef = np.min(Erefs)
#EminPlot = -0.5

weights0 = np.ones( len(Erefs) )*0.5

fit.setGlobalParams( kMorse=1.8, Lepairs=0.7 )
if bMorse:
    imodel = 2
    weights0, lens = fit.split_and_weight_curves( Erefs, x0s, n_before_min=100, weight_func=lambda E: fit.exp_weight_func(E,a=1.0, alpha=4.0) )
else:
    #imodel = 1
    imodel = 3
    weights0, lens = fit.split_and_weight_curves( Erefs, x0s, n_before_min=2, weight_func=lambda E: fit.exp_weight_func(E,a=1.0, alpha=4.0) )
#fit.setup( imodel=imodel, EvalJ=1, WriteJ=1, Regularize=1 )
fit.setup( imodel=imodel, EvalJ=1, WriteJ=1, Regularize=-1 )
# plotEWs( Erefs=Erefs, weights0=weights0, Emin=-1.5 ); plt.title( "Weighting" )
# plt.show(); exit()

fit.setWeights( weights0 )
fit.getBuffs()


#fit.setFilter( EmodelCutStart=0.0, EmodelCut=0.5, iWeightModel=2, PrintOverRepulsive=1, DiscardOverRepulsive=1, SaveOverRepulsive=-1, ListOverRepulsive=1 )
#fit.setFilter( EmodelCutStart=0.0, EmodelCut=0.5, iWeightModel=2, PrintOverRepulsive=-1, DiscardOverRepulsive=1, SaveOverRepulsive=1, ListOverRepulsive=-1 )
fit.setFilter( EmodelCutStart=0.0, EmodelCut=0.5, PrintOverRepulsive=-1, DiscardOverRepulsive=-1, SaveOverRepulsive=-1, ListOverRepulsive=-1 )
#fit.setFilter( EmodelCutStart=0.0, EmodelCut=0.5, iWeightModel=2, PrintOverRepulsive=-1, DiscardOverRepulsive=1, SaveOverRepulsive=1, ListOverRepulsive=-1 )
#fit.setFilter( EmodelCutStart=0.0, EmodelCut=0.5, PrintOverRepulsive=-1, DiscardOverRepulsive=-1, SaveOverRepulsive=-1, ListOverRepulsive=-1 )

#E,Es,Fs = fit.getEs( bOmp=False, bDOFtoTypes=False, bEs=True, bFs=False, xyz_name="all_out_debug.xyz" )
#fit.plotEWs( Erefs=Erefs, Emodel=Es, weights=fit.weights, weights0=weights0,  Emin=EminPlot ); plt.title( "BEFORE OPTIMIZATION" )
#plt.show(); #exit()

#exclude = set(["E_O3.Q", "E_HO.Q", "E_O3.H", "E_HO.H", "O_3.H", "H_O.H"])
#exclude = set(["E_O3.R", "E_HO.R", "O_3.H", "H_O.H"])
#exclude = set(["E_O3.R", "E_HO.R" ])
exclude = set([])
iDOFs = list(range(len(dof_names)))
dof_names_ = [ dof_names[i] for i in iDOFs if dof_names[i] not in exclude ]
iDOFs_     = [ iDOFs[i]     for i in iDOFs if dof_names[i] not in exclude ]
print( "dof_fnames  ", dof_names )
print( "dof_fnames_ ", dof_names_ )
print( "iDOFs_      ", iDOFs_ )
#exit()

# ------ Plot 1D parameter scans
print( "len(dof_names)", len(dof_names), dof_names )
fit.setup( imodel=imodel, Regularize=-1 )
#fit.plotDOFscans( list(range(len(dof_names))), np.linspace( -1.0+1e-6,  1.0-1e-6,  100 ), dof_names, title="DOF scan 1D" , bFs=True , bEvalSamples=True  )
# ---- plot on common range
#fit.checkDOFderiv( 1, x0=0.5, d=0.001, bEvalSamples=True )

# =====================
# GPU comparison section
# =====================

def setup_gpu_driver(xyz_file, dof_file, model_macro, hb_gate=1, regularize=False, verbose=0):
    """Create and initialize the OpenCL fitting driver with a selected model macro."""
    this_dir   = os.path.dirname(os.path.abspath(__file__))
    repo_root  = os.path.abspath(os.path.join(this_dir, '..', '..'))
    atom_types_file = os.path.join(repo_root, 'cpp', 'common_resources', 'AtomTypes.dat')
    forces_path     = os.path.join(repo_root, 'cpp', 'common_resources', 'cl', 'Forces.cl')

    fit_ocl = FittingDriver(verbose=verbose)
    fit_ocl.load_atom_types(atom_types_file)
    fit_ocl.load_data(xyz_file if os.path.isabs(xyz_file) else os.path.join(this_dir, xyz_file))
    fit_ocl.load_dofs(dof_file if os.path.isabs(dof_file) else os.path.join(this_dir, dof_file))
    fit_ocl.init_and_upload()

    macro_der = extract_macro_block(forces_path, model_macro)
    macros = {
        'MODEL_PAIR_ACCUMULATION': macro_der,
        'HBOND_GATE_DEFINE': f"#define HBOND_GATE {int(hb_gate)}",
    }
    fit_ocl.compile_with_model(macros=macros, bPrint=False)
    fit_ocl.set_regularization_enabled(enabled=bool(regularize))
    return fit_ocl

# TODO/DEBUG: Old GPU plotting helper kept for reference; superseded by unified plotting via fit.plotDOFscan_one(data=...)
# def gpu_plotDOFscan_one(...): pass

def gpu_scan_dof(fit_ocl, iDOF, xs):
    """Compute (Es, Fs) on GPU for a given DOF over grid xs without plotting.
    Returns a dict {'xs': xs, 'Es': Es, 'Fs': Fs} suitable for fit.plotDOFscan_one(data=...).
    """
    x0 = np.array([d['xstart'] for d in fit_ocl.dof_definitions], dtype=np.float32)
    Es = np.zeros_like(xs, dtype=np.float64)
    Fs = np.zeros_like(xs, dtype=np.float64)
    for j, x in enumerate(xs):
        xv = x0.copy(); xv[iDOF] = float(x)
        J, g = fit_ocl.getErrorDerivs(xv)
        Es[j] = J
        Fs[j] = g[iDOF]
    return {'xs': xs, 'Es': Es, 'Fs': Fs}

# Build GPU driver once
this_dir = os.path.dirname(os.path.abspath(__file__))
xyz_abs = fname if os.path.isabs(fname) else os.path.join(this_dir, fname)
dof_abs = dof_fname if os.path.isabs(dof_fname) else os.path.join(this_dir, dof_fname)
model_macro = 'MODEL_MorseQ_PAIR' if bMorse else 'MODEL_LJQH2_PAIR'
fit_ocl = setup_gpu_driver(xyz_abs, dof_abs, model_macro=model_macro, hb_gate=1, regularize=False, verbose=0)

############################################
# Precompute scans: first CPU for all DOFs #
############################################

def cpu_scan_dof(iDOF, xs, bEvalSamples=True):
    y_backup = fit.DOFs[iDOF]
    Es, Fs = fit.scanParam(iDOF, xs, bEvalSamples=bEvalSamples)
    fit.DOFs[iDOF] = y_backup
    return {'xs': xs, 'Es': Es, 'Fs': Fs}

fit.setVerbosity(0)
npts = 100
scans_cpu = {}
for i in iDOFs_:
    fit.loadDOFSelection(dof_fname)
    xmin = dof_specs[i]['min'] if dof_specs[i]['min'] is not None else 0.0
    xmax = dof_specs[i]['max'] if dof_specs[i]['max'] is not None else 1.0
    xs = np.linspace(xmin+1e-6, xmax-1e-6, npts)
    scans_cpu[i] = cpu_scan_dof(i, xs, bEvalSamples=True)

###########################################
# Then GPU scans for all DOFs (same order) #
###########################################

scans_gpu = {}
for i in iDOFs_:
    xs = scans_cpu[i]['xs']
    scans_gpu[i] = gpu_scan_dof(fit_ocl, iDOF=i, xs=xs)

########################
# Plot overlays per DOF #
########################

for i in iDOFs_:
    fig = plt.figure(figsize=(8,10.0))
    axE = plt.subplot(2,1,1)
    axF = plt.subplot(2,1,2)
    fit.plotDOFscan_one(i, DOFname=f"{dof_names[i]} [CPU]", bEs=True, bFs=True, verb=1, axE=axE, axF=axF, color='C0', data=scans_cpu[i])
    fit.plotDOFscan_one(i, DOFname=f"{dof_names[i]} [GPU]", bEs=True, bFs=True, verb=1, axE=axE, axF=axF, color='C1', data=scans_gpu[i])
    axE.legend(); axE.set_xlabel("DOF value"); axE.set_ylabel("E [kcal/mol]");   axE.grid( alpha=0.2)
    axF.legend(); axF.set_xlabel("DOF value"); axF.set_ylabel("F [kcal/mol/A]"); axF.grid( alpha=0.2)
    plt.suptitle(f"DOF scan: {dof_names[i]}")

plt.show()

# ------ write unoptimized results
#def getEs(imodel=0, Es=None, Fs=None, bOmp=False, bDOFtoTypes=False, bEs=True, bFs=False ):

# # ------ optimize parameters (fit)
#fit.setSwitches(EvalJ=0, WriteJ=0, CheckRepulsion=0, Regularize=0, Epairs=0)
# nstep = 5000
# fit.setSwitches(EvalJ=1, WriteJ=1  ); Err = fit.run( nstep=nstep, Fmax=1e-300, dt=0.0, imodel=imodel, iparallel=0, ialg=0, bClamp=bClamp )
# fit.setSwitches(EvalJ=1, WriteJ=-1 ); Err = fit.run( nstep=nstep, Fmax=1e-300, dt=0.0, imodel=imodel, iparallel=0, ialg=0, bClamp=bClamp )




# # ------ write optimized results
# Es = fit.getEs( imodel=imodel, isampmode=isampmode, bEpairs=bEpairs )
# np.savetxt("firecore.dat", Es )

#plot_dofs_fdofs(output_file="OUT", figsize=(12, 8))
