import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF as mmff
from pyBall import MolGUI as gui
from pyBall import Lattice2D as lat

#======== Body

Rmax = 20.0
#lat0 = np.genfromtxt("data/NaCl_sym-center.lvs")   ;print("lat1\n",lat0)
lat0 = np.genfromtxt("data/NaCl_1x1.lvs")   ;lat0=lat0[:2,:2].copy()  ;print("lat1\n",lat0)
lat1 = np.genfromtxt("data/polymer-2.lvs")  ;lat1=lat1[:2,:2].copy()  ;lat1[0,0]=16.0; lat1[0,1]=0     ;print("lat2\n",lat1)

print("DEBUG 1")
n = lat.match(lat0, lat1, Rmax=Rmax, dRmax=0.05, dAngMax=0.1  )    #;print("n", n)
print("DEBUG 2")
inds,ns,errs,Ks = lat.getMatches(inds=None, errs=None, bSort=True, Ks=(1.,1.,1.,0.) )


lat.plotLattice(lat0, plt, ls='k'   ,label='NaCL'    ) 
lat.plotLattice(lat1, plt, ls='gray',label='polymer' )
lat.plotSuperLattices( lat1, inds, plt, n=2 )
plt.legend(); plt.grid(); plt.axis('equal')
plt.show()

exit()

mmff.setVerbosity( verbosity=1, idebug=0 )
#mmff.init( xyz_name="data/pyridine", surf_name="data/NaCl_sum-center" )                             # all
W=mmff.init( xyz_name="data/pyridine", surf_name="data/NaCl_sym-center", bMMFF=False  )              # without MMFF

#mmff.init( xyz_name="data/pyridine", surf_name="data/NaCl_sym-center", bMMFF=False, gridStep=-1 )  # without gridFF
mmff.getBuffs()
mmff.eval()

gui.init(W)
gui.run( 1000000 )