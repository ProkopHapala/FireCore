import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append("../../")
from pyBall import grid_utils as gu

# ========= Functions


# ========= MAIN

path="data/NaCl_1x1_L2/"

# files=[
# "debug_BsplineCoul_pbc.npy",
# "debug_BsplineLond_pbc.npy",
# "debug_BsplinePaul_pbc.npy",
# "debug_VCoul_pbc.npy",
# "debug_VLond_pbc.npy",
# "debug_VPaul_pbc.npy",    
# ]
# #f =[0.10,0.5,0.5]
# #f2=[0.15,0.5,0.5]
# #f3=[0.20,0.5,0.5]
# f =[0.10,0.50,0.50]
# f2=[0.10,0.25,0.25]
# f3=[0.20,0.00,0.00]
# iax = 1
# fig1 = plt.figure(1,figsize=(15,10))
# #fig2 = plt.figure(2,figsize=(15,10))
# for i,name in enumerate(files):
#     dat = np.load(path+name)
#     vmax = np.abs(dat).max()
#     print( "name ", name, dat.shape," vmax=", vmax )
#     #plt.figure(fig1.number); plt.subplot(2,3,i+1);  cut1d_xyz(dat, f=f, n=2 )
#     plt.figure(fig1.number); plt.subplot(2,3,i+1);  
#     gu.cut1d_xyz(dat, f=f , n=2, iaxs=[1,2] )
#     gu.cut1d_xyz(dat, f=f2, n=2, iaxs=[1,2] )
#     gu.cut1d_xyz(dat, f=f3, n=2, iaxs=[1,2] )
#     #plt.figure(fig2.number); plt.subplot(2,3,i+1);  gu.cut2(dat,f=f[iax],iax=iax, n=[0,0], bPlot=True); plt.title(name)  

fname= "debug_VCoul_pbc.npy"
dat = np.load(path+fname)

#plt.plot( dat[:,:,1] );
#plt.plot( dat[:,1,:] );
#plt.plot( dat[1,:,:] );
nz,ny,nx = dat.shape
# iy = int( 0.5*ny )
# ix = int( 0.5*nx )
# iz = int( 0.5*nz )
# iy = int( 0.5*ny )
# ix = int( 0.5*nx )


#plt.plot( dat[iz,iy,:], label='x' );
#plt.plot( dat[iz,:,ix], label='y' );

plt.plot( dat[:,int(0.0*ny),int(0.0*nx)], label='z (0.0,0.0)' );
plt.plot( dat[:,int(0.5*ny),int(0.0*nx)], label='z (0.5,0.0)' );
plt.plot( dat[:,int(0.0*ny),int(0.5*nx)], label='z (0.0,0.5)' );
plt.plot( dat[:,int(0.5*ny),int(0.5*nx)], label='z (0.5,0.5)' );
plt.legend()
plt.grid()
plt.title(fname)
plt.show()