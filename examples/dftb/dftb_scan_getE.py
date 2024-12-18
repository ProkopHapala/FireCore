

import sys
import os

sys.path.append('../../')
from pyBall import atomicUtils as au
#from pyBall import FFFit as fff
#from pyBall import dftb_utils as dftbu
import subprocess

import numpy as np
import matplotlib.pyplot as plt
#from functools import partial

# ============= SETUP 

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

in_dir   = "/home/prokophapala/git/FireCore/tests/dftb/HBsmall_relax/scans"
#out_dir  = "/home/prokophapala/git/FireCore/tests/dftb/HBsmall_scan"
out_dir  = "/home/prokophapala/git/FireCore/tests/dftb/HBsmall_scan_3ob"
basis_path="/home/prokophapala/SIMULATIONS/dftbplus/slakos/3ob-3-1/"

# ============ Functions

def post_proc( geom, id, out_path, comment="#comment" ):
    apos,es = geom
    fname = out_path+("/%03i/stdout.log" %id)
    hosts = subprocess.check_output( 'grep "Total Energy" %s | tail -1' %fname, shell=True )
    ws = hosts.split()
    x = float(comment.split()[1])
    E = float(ws[4])
    print( id, x, E )
    return ( id, x, E )

# ============ MAIN

wd = os.getcwd()
for name in names:
    print( name )
    out_path = out_dir+"/"+name    #;print(out_path)
    #os.mkdir( out_path )
    res = au.scan_xyz( in_dir+"/"+name+"_scan.xyz", callback=post_proc, kwargs={"out_path":out_path } )
    dat = np.array(res)
    dat[:,2] *= 23.0609;
    #print( "dat ", dat )
    np.savetxt( out_dir+"/"+name+".dat", dat[:,1:] )

    plt.plot( dat[:,1], dat[:,2]-dat[-1,2], label=name )

plt.xlabel( "x [A]" )
plt.ylabel( "E [kcal/mol]" )
plt.legend()
plt.grid()
plt.show()


