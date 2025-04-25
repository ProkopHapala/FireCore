import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import plotUtils   as plu
from pyBall import MMFF as mmff

#======== Body

#mmff.setVerbosity( verbosity=1, idebug=0 )

#xs = np.linspace( -2., 2., 100 );
#y_sin = mmff.sample_func( xs*np.pi, kind=0); plt.plot( xs, y_sin, label="y_sin" )
#y_px2 = mmff.sample_func( xs,       kind=1); plt.plot( xs, y_px2, label="y_px2" )


xs = np.linspace( 0., 5.0, 100 );
#y_px2 = mmff.sample_func( xs,       kind=2, params=[-1.0,2.0,2.5,3.5],); plt.plot( xs, y_px2, label="y_px2" )

plt.figure(figsize=(5,10))

Emin=-1.0; Rmin=2.0; Rcut=2.5; Rcut2=3.5
EFs = mmff.sample_funcEF( xs, kind=2, params=[Emin, Rmin, Rcut, Rcut2]); plu.plotEF( xs, EFs, label='SR_x2' )
plt.subplot(2,1,1)
plt.axvline(Rmin, ls='--', c='r', label='E_min')
plt.axvline(Rcut, ls='--', c='g', label='E_cut')
plt.axvline(Rcut2,ls='--', c='b', label='E_cut2')
plt.subplot(2,1,2)
plt.axvline(Rmin, ls='--', c='r', label='E_min')
plt.axvline(Rcut, ls='--', c='g', label='E_cut')
plt.axvline(Rcut2,ls='--', c='b', label='E_cut2')

plt.legend()
#plt.grid()
plt.show()

print("\n\n#============= ALL DONE =============\n\n")