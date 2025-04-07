import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

sys.path.append("../../")
from pyBall import eFF as eff

from pyBall import eFF_terms as pyeff
const_bohr = 0.5291772105638411
eff.setVerbosity(0) # If 1: it will  write more stuff to console; If 0: it will wirte less stuff to console


nconf = 10

eff.load_fgo("data/H2_eFF.fgo" )   
eff.getBuffs()

# --- setup constraints
fixed_inds = np.array([
#    index, binary mask    
(0, 0b111), # 1st atom, x,y,z
(1, 0b111), # 2nd atom, x,y,z
])

# --- setup positions of fixed atoms (1st atoms at (0,0,0), 2nd at (x,0,0) )
xs = np.linspace( 0.5, 3.0,nconf,endpoint=False);      
fixed_poss = np.zeros((nconf, eff.na, 4 ))
fixed_poss[:,1,0] = xs   # set x coordinate of 2nd atom

apos, epos, Es = eff.relaxed_scan( fixed_poss, fixed_inds, nstepMax=1000, dt=1e-2, Fconv=1e-6, ialg=0 )





# plt.plot(radiusSpace, Etot_value,"o--", label='Etot', color='blue')
plt.plot( xs, Es[:,0],".-", label='Etot', color='blue')
plt.plot( xs, epos[:,0,3],".--", label='e1 size', color='blue')
plt.plot( xs, epos[:,1,3],".--", label='e2 size', color='blue')
plt.xlabel('Radius [A]')
plt.ylabel('Energy [eV]')
plt.legend()
plt.grid()
plt.show()