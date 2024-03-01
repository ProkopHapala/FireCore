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

# # ----- Test 1: sample_SplineHermite 1D 
# dx    =  1.5
# x0    = -0.1 
# Eps   = np.array( [1.0, 0.0,-1.0,-0.5,-0.4,+0.1] )
# xp    = (np.array(range(len(Eps)))-1)*dx + x0
# xs    = np.linspace(0.0,4.2,100)
# Es,Fs = mmff.sample_SplineHermite( xs, Eps, x0=x0, dx=dx )  # ;print("Fs",Fs)
# plt.figure(); 
# plt.plot(xp, Eps, 'o-k', lw=0.2, label="Eps"); 
# plt.plot(xs, Es,'.-',    label="E"); 
# plt.plot(xs, Fs,'-',     label="F_ana");  
# plt.plot(xs[1:-1],numDeriv(xs,Es),':', label="F_num"); 
# plt.grid(); plt.legend()


# ----- Test 2: sample_SplineHermite 2D

Eps   = np.array([ 
    [1.0, 0.0,-1.0,-0.5,-0.4,+0.1], 
    [1.0, 0.0,-1.0,-0.5,-0.4,+0.1],
    [1.0, 0.0,-2.0,-1.5,-0.4,+0.1], 
    [1.0, 0.0,-3.0,-2.5,-1.4,+0.1], 
    [1.0, 0.0,-1.0,-0.5,-0.4,+0.1],
    [1.0, 0.0,-1.0,-0.5,-0.4,+0.1],
])


nsamp=100
xs = np.linspace(0.0,4.2,nsamp)
ps = np.zeros((nsamp,2))

# # ----- Test 2.1: sample_SplineHermite2D --- line along x
# ps[:,1] = 1.0
# ps[:,0] = xs
# FEs = mmff.sample_SplineHermite2D( ps, Eps, g0=[-0.1,-0.1], dg=[1.5,1.5] )  # ;print("Fs",Fs)
# plt.figure();
# plt.plot( xs, FEs[:,2], '.-', label="E" )
# plt.plot( xs, FEs[:,0], '.-', label="Fx" )
# plt.plot( xs, FEs[:,1], '.-', label="Fy" )
# plt.plot( xs[1:-1],numDeriv(xs,FEs[:,2]), ':',  label="Fx_num" )
# plt.grid(); plt.legend()

# # ----- Test 2.1: sample_SplineHermite2D --- line along y
# ps[:,1] = xs
# ps[:,0] = 1.0
# FEs = mmff.sample_SplineHermite2D( ps, Eps, g0=[-0.1,-0.1], dg=[1.5,1.5] )  # ;print("Fs",Fs)
# plt.figure();
# plt.plot( xs, FEs[:,2], '.-', label="E" )
# plt.plot( xs, FEs[:,0], '.-', label="Fx" )
# plt.plot( xs, FEs[:,1], '.-', label="Fy" )
# plt.plot( xs[1:-1],numDeriv(xs,FEs[:,2]), ':',  label="Fy_num" )
# plt.grid(); plt.legend()


#nsamp=200
nsamp=10
extent = [0.0,4.2,0.0,4.2]
xs,ys = np.meshgrid( np.linspace(extent[0],extent[1],nsamp), np.linspace(extent[2],extent[3],nsamp) )
ps = np.vstack([xs.flatten(), ys.flatten()]).T.copy()
print("ps",ps)
print("ps.shape: ",ps.shape)

FEs = mmff.sample_SplineHermite2D( ps, Eps, g0=[-0.1,-0.1], dg=[1.5,1.5] )  # ;print("Fs",Fs)
FEs = FEs.reshape(nsamp,nsamp,3)

cmap='plasma'
plt.figure(figsize=(15,5));
plt.subplot(1,3,1); plt.imshow( FEs[:,:,2], extent=extent, origin='lower', cmap=cmap ); plt.colorbar(); plt.title("E")
plt.subplot(1,3,2); plt.imshow( FEs[:,:,0], extent=extent, origin='lower', cmap=cmap ); plt.colorbar(); plt.title("Fx")
plt.subplot(1,3,3); plt.imshow( FEs[:,:,1], extent=extent, origin='lower', cmap=cmap ); plt.colorbar(); plt.title("Fy")



plt.show()