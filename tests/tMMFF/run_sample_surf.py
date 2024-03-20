import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF as mmff

mmff.init( xyz_name="data/H2O", surf_name="data/NaCl_1x1_L2", bMMFF=False )


# ======= Setup

kind = 12
RvdW = 1.487
#EvdW = 1.0
EvdW = 0.0006808
Q    = +0.2 

# ========= Main

# ========== Test 1: 1D plot
'''
n  = 200
xs = np.linspace( -10., 10.0,  n )    #;print(xs)
ps = np.zeros( (n,3) )
ps[:,2] = xs
FEs = mmff.sampleSurf_vecs( ps, kind=kind, RvdW=RvdW, EvdW=EvdW, Q=Q, npbc=5 )
Emin = FEs[:,3].min()
plt.plot( xs, FEs[:,3] )
plt.ylim(Emin*1.2,-2*Emin)
plt.grid()
'''

# ========== Test 2: 2D plot
'''
nsamp=200
extent = [ -10.0,10.0,  -10.0,10.0 ]
xs,ys = np.meshgrid( np.linspace(extent[0],extent[1],nsamp), np.linspace(extent[2],extent[3],nsamp) )
#ps = np.vstack([xs.flatten(), ys.flatten(), ys.flatten()*0.0 ]).T.copy()

ps = np.vstack([xs.flatten(), ys.flatten()*0.0, ys.flatten() ]).T.copy()

# #print("ps",ps)
# print("ps.shape: ",ps.shape)

FEs = mmff.sampleSurf_vecs( ps, kind=kind, RvdW=RvdW, EvdW=EvdW, Q=Q, npbc=5 )
Emin = FEs[:,3].min();   Fmin=Emin*1.0
FEs = FEs.reshape(nsamp,nsamp,4)

#cmap='plasma'
cmap='bwr'
plt.figure(figsize=(20,5));
plt.subplot(1,4,1); plt.imshow( FEs[:,:,3], vmin=Emin*1.2,vmax=-1.2*Emin, extent=extent, origin='lower', cmap=cmap ); plt.colorbar(); plt.title("E")
plt.subplot(1,4,2); plt.imshow( FEs[:,:,0], vmin=Fmin,vmax=-Fmin, extent=extent, origin='lower', cmap=cmap ); plt.colorbar(); plt.title("Fx")
plt.subplot(1,4,3); plt.imshow( FEs[:,:,1], vmin=Fmin,vmax=-Fmin, extent=extent, origin='lower', cmap=cmap ); plt.colorbar(); plt.title("Fy")
plt.subplot(1,4,4); plt.imshow( FEs[:,:,2], vmin=Fmin,vmax=-Fmin, extent=extent, origin='lower', cmap=cmap ); plt.colorbar(); plt.title("Fz")
'''

# ========== Test 2: 2D different FFs

nsamp=200
extent = [ -10.0,10.0,  -10.0,10.0 ]
xs,ys = np.meshgrid( np.linspace(extent[0],extent[1],nsamp), np.linspace(extent[2],extent[3],nsamp) )

ps = np.vstack([xs.flatten(), ys.flatten()*0.0, ys.flatten() ]).T.copy()

FEs_f   = mmff.sampleSurf_vecs( ps, kind=12, RvdW=RvdW, EvdW=EvdW, Q=Q, npbc=5 )
FEs_d   = mmff.sampleSurf_vecs( ps, kind=13, RvdW=RvdW, EvdW=EvdW, Q=Q, npbc=5 )
FEs_cub = mmff.sampleSurf_vecs( ps, kind=14, RvdW=RvdW, EvdW=EvdW, Q=Q, npbc=5 )

Emin = FEs_f[:,3].min();   Fmin=Emin*1.0
FEs_f   = FEs_f.reshape(nsamp,nsamp,4)
FEs_d   = FEs_d.reshape(nsamp,nsamp,4)
FEs_cub = FEs_cub.reshape(nsamp,nsamp,4)

#cmap='plasma'
cmap='bwr'
plt.figure(figsize=(15,5));
# plt.subplot(1,3,1); plt.imshow( FEs_f  [:,:,3], vmin=Emin*1.2,vmax=-1.2*Emin, extent=extent, origin='lower', cmap=cmap ); plt.colorbar(); plt.title("E_float")
# plt.subplot(1,3,2); plt.imshow( FEs_d  [:,:,3], vmin=Emin*1.2,vmax=-1.2*Emin, extent=extent, origin='lower', cmap=cmap ); plt.colorbar(); plt.title("E_double")
#plt.subplot(1,3,3); plt.imshow( FEs_cub[:,:,3], vmin=Emin*1.2,vmax=-1.2*Emin, extent=extent, origin='lower', cmap=cmap ); plt.colorbar(); plt.title("E_tricubic")


plt.subplot(1,3,1); plt.imshow( FEs_f  [:,:,3], extent=extent, origin='lower', cmap=cmap ); plt.colorbar(); plt.title("E_float")
plt.subplot(1,3,2); plt.imshow( FEs_d  [:,:,3], extent=extent, origin='lower', cmap=cmap ); plt.colorbar(); plt.title("E_double")
plt.subplot(1,3,3); plt.imshow( FEs_cub[:,:,3], extent=extent, origin='lower', cmap=cmap ); plt.colorbar(); plt.title("E_tricubic")
plt.show()