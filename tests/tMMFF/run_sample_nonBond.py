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

#======== Body

REQi = (1.487,np.sqrt(0.0006808),0.0)
REQj = (1.487,np.sqrt(0.0006808),0.0)
eij = np.sqrt(REQi[1]*REQi[1])

#mmff.sample_evalAngleCos( xs, lmin=1, lmax=1, kmin=1, kmax=1, flim=1e+300, Es=None, Fs=None)
xs    = np.linspace(0.0,6.0,1000)
#xs    = np.linspace(0.0,6.0,100)


#Es,Fs = mmff.sample_evalAngleCosHalf( xs, ang0=np.pi*0.1 )
EsSR,FsSR = mmff.sampleNonBond( xs, kind=4, REQi=REQi, REQj=REQj, K=eij*0.2, Rdamp=0.5 )   # short-range repulsion R4func
#EsMo,FsMo = mmff.sampleNonBond( xs, kind=1, REQi=REQi, REQj=REQj, K=1.0, Rdamp=0.1 )       # Morse
EsLJ,FsLJ = mmff.sampleNonBond( xs, kind=2, REQi=REQi, REQj=REQj )       # L-J

#FnumSR = numDeriv(xs,EsSR)     #ratios = Fs[1:-1]/Fnum    ;print("ratios", ratios)
#Es,Fs = cos_half( xs + np.pi*0.5 )
xs/=np.pi
plt.figure(); 
plt.plot(xs, EsSR-eij, label="E_SR"); 
plt.plot(xs, FsSR, label="F_SR");  
#plt.plot(xs[1:-1], -FnumSR, label="F_num"); 
#plt.plot(xs, EsMo, label="E_Mo"); 
#plt.plot(xs, FsMo, label="F_Mo"); 
plt.plot(xs, EsLJ, label="E_LJ"); 
plt.plot(xs, FsLJ, label="F_LJ"); 
plt.grid(); 
plt.legend()




plt.ylim(-eij*2,eij*2)


plt.show()