

import sys
import os
import psi4

sys.path.append('../../')
from pyBall import atomicUtils as au
from pyBall import FFFit as fff
from pyBall import psi4_utils as psi4u

import psi4
import numpy as np
import matplotlib.pyplot as plt
#from functools import partial


# https://psicode.org/psi4manual/master/api/psi4.core.Molecule.html
# https://github.com/psi4/psi4numpy/blob/master/Tutorials/01_Psi4NumPy-Basics/1b_molecule.ipynb    
#    * Example: Fitting Lennard-Jones Parameters from Potential Energy Scan


# ============= SETUP 

params = {  # These are the default settings
#    'gradientmax': 0.45e-6,  # Eh/[Bohr|rad]
#    'gradientrms': 0.15e-6,  # Eh/[Bohr|rad]
#    'stepmax': 1.8e-3,       # [Bohr|rad]
#    'steprms': 1.2e-3,       # [Bohr|rad]

#'method' : 'scf',
#'method' : 'pbe',
#'method' : 'b3lyp',
'method' : 'mp2',

'basis'  : 'sto-3g',
# 'basis'  : '6-31+G',
# 'basis'  : '6-311+G*',
# 'basis'  : '6-311++G**',
# 'basis'  : '6-311++G(3df,3pd)',
#'basis'  : 'cc-pvdz',
# 'basis'  : 'aug-cc-pvtz',
# 'basis'  : 'def2-QZVPPD',

}

# ============ MAIN

# ------ H2O Scan
# apos,Zs,es,qs = au.loadAtomsNP( fname='input/H2O.xyz' )   ;print("apos ", apos) ;print("es ", es)
# apos1=au.orient( 0,(1,2),(1,0), apos, _0=0, trans=(2,1,0) )
# apos2=au.orient( 0,(0,1),(0,2), apos, _0=0, trans=(0,2,1) )     ;apos2[:,1]*=-1.0 ;apos2[:,1]+=4.0
# xs = np.arange(-1.5,2.0,0.1) # ;print("shifts ", xs)
# Es = fff.linearScan( (apos1,es), (apos2,es), xs, dir=(0.,1.,0.), xyz_file="scan_H2O-H2O_lin.xyz"  )

# ------ HCOOH dimer scan
# apos,Zs,es,qs = au.loadAtomsNP( fname='input/HCOOH.xyz' )   ;print("apos ", apos) ;print("es ", es)
# apos1=au.orient( 0,(0,1),(2,3), apos, _0=0, trans=(2,1,0) )
# apos2=au.orient( 0,(1,0),(3,2), apos, _0=0, trans=(2,1,0) )    ;apos2[:,0]+=7.5
# #apos2=au.orient( 0,(0,1),(0,2), apos, _0=0, trans=(0,2,1) )   # ;apos2[:,1]*=-1.0 ;apos2[:,1]+=4.0
# xs = np.arange(-1.5,2.0,0.1)  #;print("shifts ", xs)
# Es = fff.linearScan( (apos1,es), (apos2,es), xs, dir=(1.,0.,0.), xyz_file="scan_HCOOH-dimer_lin.xyz"  )

apos,Zs,es,qs = au.loadAtomsNP( fname='input/HCOOH_dimer.xyz' )   #;print("apos ", apos) ;print("es ", es)

'''
psi4.core.set_output_file('relax.log', False)
mol = psi4u.relax( (apos,es), params=params )
mol.save_xyz_file("HCOOH_relaxed.xyz",1)
#apos, es = psi4u.
'''

params['ifrag_line']=5
#xs         = np.arange(0,4.0,0.1) #;print()
xs         = fff.exprange( 4.0, n=20 ) #;print(len(xs))
#Es = np.ones(len(xs))
selection  = list(range(5,10))    #;print("selection ",selection)
Es = fff.linearScan_1mol( (apos,es), selection, xs, dir=(1.,0,0), xyz_file="scan_in.xyz" )

# ---- from xyz movie
psi4.core.set_output_file('scan.log', False)
fff.verbosity=1
Es, xs = fff.scan_xyz( "scan_in.xyz", fxyzout="scan_out.xyz", Eout=None, params=params, callback=psi4u.eval, xs=xs )
#Es, xs = fff.scan_xyz( "scan_H2O-H2O_lin.xyz", fxyzout="scan_out.xyz", Eout=None, params=params, callback=psi4u.eval )



plt.plot( xs, Es ,'.-' );   plt.grid(); 
plt.show()
