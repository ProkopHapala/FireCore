

import sys
import os

sys.path.append('../../')
from pyBall import atomicUtils as au
from pyBall import FFFit as fff
from pyBall import dftb_utils as dftbu
from pyBall import psi4_utils as psi4u

#import numpy as np
#import matplotlib.pyplot as plt
#from functools import partial

# ============= SETUP 


params_block='''
    scf_type df 
''' 

# params_block='''
#     scf_type df 
    
#     opt_type MIN
#     geom_maxiter 1000
#     g_convergence INTERFRAG_TIGHT
#     print_trajectory_xyz_file true
#     DYNAMIC_LEVEL 1
# ''' 

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

in_dir   = "/home/prokophapala/git/FireCore/tests/tPsi4resp/scans"
out_dir  = "/home/prokophapala/git/FireCore/tests/tPsi4resp/HBsmall_scan"

bOpt = False
mem  = '15GB'
method="b3lyp-d3"
basis="cc-pVDZ"

# ============ Functions

def make_job( geom, id, params, out_path, comment="#comment" ):
    apos,es = geom
    dirname_ = out_path+("/%03i" %id)
    os.mkdir( dirname_ )
    #au.saveXYZ       ( es, apos, dirname_+"/input.xyz", comment=comment )
    atoms = au.AtomicSystem( enames=es, apos=apos )
    ins,outs = atoms.selectBondedCluster( {0} )
    nhyphen=len(ins)-1
    lines  = au.geomLines( apos, es )
    #dftbu.makeDFTBjob( enames=es, fname=dirname_+"/dftb_in.hsd", params=params, basis_path=basis_path, opt=False )
    psi4u.write_psi4_in( lines, fname=dirname_+"/psi.in", nhyphen=nhyphen,  mem=mem, method=method, basis=basis, bsse="['cp','nocp']", params_block=params, opt=bOpt )

# ============ MAIN

wd = os.getcwd()
for name in names:
    out_path = out_dir+"/"+name    ;print(out_path)
    os.mkdir( out_path )
    au.scan_xyz( in_dir+"/"+name+"_scan.xyz", callback=make_job, kwargs={"out_path":out_path, "params":params_block } )



