# https://pyscf.org/quickstart.html
# https://pyscf.org/pyscf_api_docs/pyscf.gto.html#module-pyscf.gto.mole

'''
Molecules:
H2O
NH3
HCOOH
CH2=O
'''

import sys
import os
#from tkinter import UNITS
import pyscf
sys.path.append('../../')
from pyBall import atomicUtils as au
from pyBall import FFFit as fff
from pyBall import pyscf_utils as scfu

import numpy as np
import matplotlib.pyplot as plt
#from functools import partial


# ============= Setup

# ============ MAIN

mol = scfu.preparemol(fname='relax.xyz')

apos,es = scfu.unpack_mol( mol )   #;print("apos ",apos) ;print("es ",es)
apos1=au.orient( 0,(1,2),(1,0), apos, _0=0, trans=(2,1,0) )
apos2=au.orient( 0,(0,1),(0,2), apos, _0=0, trans=(0,2,1) )     ;apos2[:,1]*=-1.0 ;apos2[:,1]+=4.0
#mol=pyscf.M(atom= pack_mol( apos1, es ), unit='B')   #;print( " mol._atom \n", mol.atom )
#mol.atom.extend( pack_mol( apos2, es ) )             #;print( " mol._atom \n", mol.atom )
#mol.build()
#print( " mol._atom \n", mol._atom )
#mol_1=pyscf.M(atom=mol, unit='B')

# ------- Linear scan
#xs = np.arange(-1.5,4.0,0.2)  ;print("shifts ", xs)
#Es = linearScan( (apos1,es), (apos2,es), xs, dir=(0.,1.,0), bEcho=True, xyz_file="lin_scan.xyz",  Eout_file="lin_scan_Es.dat" )

# -------- Angular scan
#nrot = 16
#dang= (np.pi*(3./2.))/nrot  ;print("dang ", dang )
#dang= 2.0*np.pi/nrot  ;print("dang ", dang )
#rot0 = au.makeRotMatAng( dang*nrot*0.5, ax1=1, ax2=2 ).transpose()
#rot  = au.makeRotMatAng( dang , ax1=1, ax2=2 ).transpose()
#print( rot ); # exit()
#mulpos( apos2, rot0 )
#Es,xs = angularScan( nrot, dang, (apos1,es), (apos2,es), rot=rot, bEcho=True, xyz_file="ang_scan.xyz", Eout_file="ang_scan_Es.dat" )


# -------- Angular scan
#nrot = 16
#dang = -1.5*np.pi/nrot  ;print("dang ", dang )
#Es,xs = fff.angularScan( nrot,dang, (apos1,es),(apos2,es), ang0=dang*-0.5*nrot, ax1=1, ax2=2, xyz_file="scan_H2O-H2O_rot.xyz", Eout_file=None, callback=None )
#Es,xs = fff.angularScan( nrot, dang, (apos1,es), (apos2,es), ang0=dang*0.5*nrot, ax1=1, ax2=2, xyz_file="scan_H2O-H2O_rot.xyz", Eout_file=None, callback=scfu.evalHf )

# ---- Linear
#s = np.arange(-1.5,2.0,0.1)  ;print("shifts ", xs)
#Es = fff.linearScan( (apos1,es), (apos2,es), xs, dir=(0.,1.,0.), xyz_file="scan.xyz", Eout_file=None, callback=None )
#Es = fff.linearScan( (apos1,es), (apos2,es), xs, dir=(0.,1.,0.), xyz_file="scan_H2O-H2O_lin.xyz", Eout_file=None, callback=scfu.evalHf )

# ---- from xyz movie
fff.verbosity=1
Es, xs = fff.scan_xyz( "scan_H2O-H2O_lin.xyz", fxyzout="scan_out.xyz", Eout=None, callback=scfu.evalHf )


#np.savetxt( "E_scan.dat", np.array( [xs, Es*hartree2eV] ).transpose() )
plt.plot( xs, Es );   plt.grid(); 
plt.show()