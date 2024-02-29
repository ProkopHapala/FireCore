import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF as mmff

def numDeriv( xs, Es): 
    dx = xs[1]-xs[0]
    Fs = (Es[2:]-Es[:-2])/(2*dx)
    return Fs

# ----- SplineConstr
dx    = 1.5
x0    = 0.5 
Eps   = np.array( [1.0, 0.0,-1.0,-0.5,-0.2,-0.1] )
xp    = (np.array(range(len(Eps)))-1)*dx + x0
xs    = np.linspace(0.0,6.0,100)
Es,Fs = mmff.sample_SplineHermite( xs, Eps, x0=x0, dx=dx )  # ;print("Fs",Fs)
plt.figure(); 
plt.plot(xp, Eps, 'o-k', lw=0.2, label="Eps"); 
plt.plot(xs, Es,'.-',    label="E"); 
plt.plot(xs, Fs,'-',     label="F_ana");  
plt.plot(xs[1:-1],numDeriv(xs,Es),':', label="F_num"); 
plt.grid(); plt.legend()

plt.show()