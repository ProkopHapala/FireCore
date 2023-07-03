

import sys
import os

sys.path.append('../../')
from pyBall import atomicUtils as au
from pyBall import FFFit as fff
from pyBall import dftb_utils as dftbu

import numpy as np
import matplotlib.pyplot as plt
#from functools import partial

# ============= SETUP 

params = dftbu.default_params.copy()
params.update({
    'own_dir': True,
    'method':'D3H5',
    'cell':None,
    'basis':"/home/prokop/SIMULATIONS/dftbplus/slakos/3ob-3-1/",
    'opt':False,
    'Temperature' : 300
})


params_ = dftbu.default_params.copy()
params_.update({
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


# ============ MAIN

#xs = np.arange(0.3,4.0,0.2)  #;print("shifts ", xs)

xs = fff.scan_ranges(  )
#print(len(xs))
#plt.plot( xs, np.ones(len(xs)), '.-' )
#plt.show()

names=[
"2_ammonia_ammonia",
"2_formamide_formamide",
"2_formicAcid_formicAcid", 
"2_HCN_HCN",
"2_water_water",

"2_carbonicAcid_carbonicAcid",
"2_cyanidenitrate_urea",
"2_formicAcid_formamide",
"2_nitro_diamine",
"2_urea_urea",
"2_water_ammonia",

"2_water_carbonicAcid",
"2_water_urea",
]

in_dir   = "/home/prokophapala/git/FireCore/tests/dftb/HBsmall_relax"
out_dir  = "/home/prokophapala/git/FireCore/tests/dftb/HBsmall_relax/scans"


wd = os.getcwd()
for name in names:
    print( name )
    dirname=in_dir+'/'+name  
    atoms = au.AtomicSystem( fname=dirname+"/geom.out.xyz" )
    ins,outs = atoms.selectBondedCluster( {0} )
    Es = fff.linearScan_1mol( (atoms.apos,atoms.enames), outs, xs, dir=(1.,0,0), xyz_file=out_dir+"/"+name+"_scan.xyz" )

    # os.mkdir( name )
    # os.chdir( name )
    #os.system("jmol_ scan_in.xyz")
    # Es, xs = fff.scan_xyz( "scan_in.xyz", fxyzout="scan_out_.xyz", Eout=None, params=params, callback=dftbu.run, xs=xs )
    # plt.plot( xs, Es, '.-' )
    # np.savetxt( "scan_.E.dat", np.array([xs,Es]).transpose() )
    os.chdir( wd )

plt.show()


