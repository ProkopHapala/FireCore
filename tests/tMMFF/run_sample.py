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

def cos_half( a ):
    E = (1-np.cos(a*0.5))*0.5
    F = np.sin(a*0.5)*0.25
    return E,F

#======== Body

'''
#mmff.sample_evalAngleCos( xs, lmin=1, lmax=1, kmin=1, kmax=1, flim=1e+300, Es=None, Fs=None)
xs    = np.linspace(-np.pi,np.pi,100)
#Es,Fs = mmff.sample_evalAngleCos( xs, c0=0.5 )
#Es,Fs = mmff.sample_evalAngleCos( xs, c0=-0.5 )
Es,Fs = mmff.sample_evalAngleCos( xs, c0=-0.99 )
#Es,Fs = cos_half( xs + np.pi*0.5 )
plt.figure(); plt.plot(xs, Es, label="E"); plt.plot(xs, Fs, label="F_ana");  plt.plot(xs[1:-1], numDeriv(xs,Es), label="F_num"); plt.grid(); plt.legend()
'''

#mmff.sample_DistConstr( xs, lmin=1, lmax=1, kmin=1, kmax=1, flim=1e+300, Es=None, Fs=None)
xs    = np.linspace(0.0,3.0,100)
Es,Fs = mmff.sample_DistConstr( xs, lmin=0.9, lmax=1.2, flim=0.5 )   ;print("Fs",Fs)
plt.figure(); plt.plot(xs, Es,'.-', label="E"); plt.plot(xs, Fs, label="F_ana");  plt.plot(xs[1:-1], numDeriv(xs,Es), label="F_num"); plt.grid(); plt.legend()


plt.show()