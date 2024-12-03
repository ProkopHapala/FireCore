

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

# ============ MAIN

#xs = np.arange(0.3,4.0,0.2)  #;print("shifts ", xs)

xs = fff.scan_ranges(  )
#print(len(xs))
#plt.plot( xs, np.ones(len(xs)), '.-' )
#plt.show()

#names=['AT','CG']
names=['CG']

wd = os.getcwd()
for name in names:
    os.mkdir( name )
    os.chdir( name )
    fname = "../"+name+".xyz"
    print( fname )
    atoms = au.AtomicSystem( fname=fname )
    ins,outs = atoms.selectBondedCluster( {0} )
    Es = fff.linearScan_1mol( (atoms.apos,atoms.enames), outs, xs, dir=(1.,0,0), xyz_file="scan_in.xyz" )
    #os.system("jmol_ scan_in.xyz")
    Es, xs = fff.scan_xyz( "scan_in.xyz", fxyzout="scan_out_.xyz", Eout=None, params=params, callback=dftbu.run, xs=xs )
    plt.plot( xs, Es, '.-' )
    np.savetxt( "scan_.E.dat", np.array([xs,Es]).transpose() )
    os.chdir( wd )

plt.show()


