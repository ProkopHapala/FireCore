#!/usr/bin/python


import os
from pyBall import  atomicUtils as au
from pyBall import dftb_utils as dftbu
from pyBall import Process as pc

names=[
"2_ammonia_ammonia",
"2_formamide_formamide",
"2_formicAcid_formicAcid", 
"2_HCN_HCN",
"2_water_ammonia", 
"2_water_water",
"ammonia",
"formamide",
"formicAcid",  
"HCN",
"water",
]

# ========= Setup

workdir="/home/prokop/git/FireCore/tests/dftb/input_relaxed_tight-bak/"

opt_frag = True

params = dftbu.default_params.copy()
params.update({
    "Optimizer"   : "Rational{}",
    #"Optimizer"   : "LBFGS{  Memory = 20 }",
    #"GradElem"    : 1E-4,
    #"DispElem"    : 1E-3,
    #"EConv"       : 1E-7,
    #'Temperature' : 50,
    'Temperature' : 300,
    #'Mixer': 'Broyden{ MixingParameter = 0.02 }',
    'Mixer': 'Anderson{ MixingParameter = 0.05 }',
    #'SCCTolerance' : 1e-7,
    #'MaxSccIterations' : 200,
})

nproc_max = 10


# ========= Functions

            
# ========= Main

all_procs = []
procs     = []

workdir = cwd = os.getcwd()

bFrags = True
out_prefix = ""
in_dir     = "input_xyz" 

for name in names:
    atoms = au.AtomicSystem( in_dir + '/'+name+".xyz" )

    '''    
    if bFrags:
        dirname1 = name+"_1"
        dirname2 = name+"_2"   

        ins,outs = atoms.selectBondedCluster( {0} )

        os.mkdir(dirname1)
        A = atoms.selectSubset( ins  )
        A.saveXYZ(dirname1+"/input.xyz", bQs=False)
        dftbu.makeDFTBjob( enames=A.enames, fname=dirname1+"/dftb_in.hsd", params=params, opt=opt_frag )
        os.chdir( dirname1 )
        procs = pc.addJobWhenFree( dirname1, ['dftb+',"psi.in"], procs, nproc_max=nproc_max, wait_sec=0.5 )
        os.chdir( workdir )

        os.mkdir(dirname2)
        B = atoms.selectSubset( outs )
        B.saveXYZ(dirname2+"/input.xyz", bQs=False)
        dftbu.makeDFTBjob( enames=B.enames, fname=dirname2+"/dftb_in.hsd", params=params, opt=opt_frag )
        os.chdir( dirname2 )
        procs = pc.addJobWhenFree( dirname2, ['dftb+',"psi.in"], procs, nproc_max=nproc_max, wait_sec=0.5 )
        os.chdir( workdir )
    '''
        
    dirname = out_prefix+"/"+name
    os.mkdir( dirname )
    dftbu.makeDFTBjob( enames=atoms.enames, fname=dirname+"/dftb_in.hsd", params=params )
    os.chdir( dirname )
    procs = pc.addJobWhenFree( dirname, ['psi4',"psi.in"], procs, nproc_max=nproc_max, wait_sec=0.5 )
    os.chdir( workdir )
            
pc.wait_proc_finish( procs )
print( "DONE " )

