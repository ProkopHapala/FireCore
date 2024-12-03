

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

'''
params = dftbu.default_params.copy()
params.update({
    'own_dir': True,
    'method':'D3H5',
    'cell':None,
    'basis':"/home/prokop/SIMULATIONS/dftbplus/slakos/3ob-3-1/",
    'opt':False,
    'Temperature' : 300
})
'''


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

names=[
# "2_ammonia_ammonia",
# "2_formamide_formamide",
# "2_formicAcid_formicAcid", 
 "2_HCN_HCN",
# "2_water_water",

# "2_carbonicAcid_carbonicAcid",
# "2_cyanidenitrate_urea",
# "2_formicAcid_formamide",
# "2_nitro_diamine",
# "2_urea_urea",
# "2_water_ammonia",

# "2_water_carbonicAcid",
# "2_water_urea",
]

in_dir   = "/home/prokophapala/git/FireCore/tests/dftb/HBsmall_relax/scans"
#out_dir  = "/home/prokophapala/git/FireCore/tests/dftb/HBsmall_scan"
out_dir  = "/home/prokophapala/git/FireCore/tests/dftb/HBsmall_scan_3ob"
basis_path="/home/prokophapala/SIMULATIONS/dftbplus/slakos/3ob-3-1/"

# ============ Functions

def make_job( geom, id, params, basis_path, out_path, comment="#comment" ):
    apos,es = geom
    dirname_ = out_path+("/%03i" %id)
    os.mkdir( dirname_ )
    au.saveXYZ       ( es, apos,        dirname_+"/input.xyz", comment=comment )
    dftbu.makeDFTBjob( enames=es, fname=dirname_+"/dftb_in.hsd", params=params, basis_path=basis_path, method='3ob', opt=False )

# ============ MAIN

wd = os.getcwd()
for name in names:
    out_path = out_dir+"/"+name    ;print(out_path)
    os.mkdir( out_path )
    au.scan_xyz( in_dir+"/"+name+"_scan.xyz", callback=make_job, kwargs={"out_path":out_path, "params":params, "basis_path":basis_path } )

plt.show()


