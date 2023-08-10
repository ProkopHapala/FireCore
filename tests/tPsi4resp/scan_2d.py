

import sys
import os


sys.path.append('../../')
from pyBall import atomicUtils as au
from pyBall import FFFit as fff
from pyBall import psi4_utils as psi4u

import numpy as np
import matplotlib.pyplot as plt
#from functools import partial

import psi4;   psi4u.psi4=psi4


# https://psicode.org/psi4manual/master/api/psi4.core.Molecule.html
# https://github.com/psi4/psi4numpy/blob/master/Tutorials/01_Psi4NumPy-Basics/1b_molecule.ipynb    
#    * Example: Fitting Lennard-Jones Parameters from Potential Energy Scan


# ============= SETUP 

params = {  # These are the default settings
"geom_maxiter": 1000,                # increase limit for geometry relaxation
"intrafrag_step_limit"    : 0.1,    # this helps with geometry relaxation convergence
"intrafrag_step_limit_min": 0.1,
"intrafrag_step_limit_max": 0.1,
"opt_coordinates" : "cartesian",
"step_type":  "nr",
"print_trajectory_xyz_file": True,

#'method' : 'scf',
#'method' : 'pbe',
#'method' : 'b3lyp',
'method' : 'b3lyp-d3',
#'method' : 'mp2',

#'basis'  : 'sto-3g',
# 'basis'  : '6-31+G',
# 'basis'  : '6-311+G*',
# 'basis'  : '6-311++G**',
# 'basis'  : '6-311++G(3df,3pd)',
'basis'  : 'cc-pvdz',
# 'basis'  : 'aug-cc-pvtz',
# 'basis'  : 'def2-QZVPPD',
'bsse'     : 'cp'

}

# ============ MAIN


'''
#mol = au.AtomicSystem( fname='HBond_OCH2_vs_H2O.xyz' )
#mol = au.AtomicSystem( fname='HBond_OCH2_vs_H2O-x-2.xyz' )
mol = au.AtomicSystem( fname='HBond_H2O_vs_H2O.xyz' )
#mol.orient( 2, (5,2), (5,6), trans=(2,1,0)  )
#mol.orient( 1, (1,2), (2,3), trans=(2,1,0)  )
#mol.saveXYZ('HBond_OCH2_vs_H2O-.xyz' )

#os.system("jmol_ name.xyz")


#xs = np.arange(-0.3,5.0,0.1)  #;print("shifts ", xs)
rs   = np.arange(-0.6,2.0,0.1)             ;print("rs   ", rs   )
#rs    = np.array([0.]) 
#ang0  = np.pi*(0.825/5. - 0.5 )
#angs  = np.array([ang0]) 

#angs  = np.arange(-0.5,0.5+1e-8,0.1) * np.pi + ang0 
angs  = np.arange(0.0,0.5+1e-8,0.1) * np.pi #+ ang0 

#angs = np.arange(-np.pi*3/4.,np.pi/4,0.1)      ;print("angs ", angs )
#angs = np.arange(-np.pi*3/4.,0.0,0.1)      ;print("angs ", angs )

#fff.linearScan_1mol( (mol.apos,mol.enames), [0,1,2], xs, dir=(1.,0,0), xyz_file="scan_in.xyz" )

#fff.angularScan_1mol(  (mol.apos,mol.enames), [4,5,6], rs, angs, ax1=0,ax2=1, xyz_file="scan_in.xyz" )

#fff.angularScan_1mol_vecs(  (mol.apos,mol.enames), [4,5,6], rs, angs, dir=(0.0,0.0,1.0),up=(0.0,1.0,0.0), xyz_file="scan_in.xyz" )
#fff.angularScan_1mol_vecs(  (mol.apos,mol.enames),  [4,5,6], rs, angs, dir=(0.0,0.0,-1.0),up=(0.0,1.0,0.0), xyz_file="scan_in.xyz" )
fff.angularScan_1mol(  (mol.apos,mol.enames),  [3,4,5], rs, angs, ax1=0,ax2=1, dir=(-1.0,0.0,0.0), xyz_file="scan_in.xyz" )

os.system("jmol_ scan_in.xyz")
'''


params['ifrag_line'] = 3
Es, xs = fff.scan_xyz( "scan_in.xyz", fxyzout="scan_out_.xyz", Eout=None, params=params, callback=psi4u.eval )
#Es, xs = fff.scan_xyz( "scan_in.xyz", fxyzout="scan_out_.xyz", Eout=None, params=params )


exit()

