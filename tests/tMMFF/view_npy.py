import numpy as np
import matplotlib.pyplot as plt

dat = np.load("data/NaCl_1x1_L2/debug_BsplineCoul_pbc.npy")

print( "npy.shape ", dat.shape )

(nx,ny,nz) = dat.shape

print( "min,max: ", dat.min(), dat.max() )

plt.figure(figsize=(15,5))
#plt.subplot(1,3,1); plt.imshow( dat[nx//2,:,:], vmin=-1.0, vmax=1.0 ); plt.colorbar()
#plt.subplot(1,3,2); plt.imshow( dat[:,ny//2,:], vmin=-1.0, vmax=1.0 ); plt.colorbar()
#plt.subplot(1,3,3); plt.imshow( dat[:,:,nz//2], vmin=-1.0, vmax=1.0 ); plt.colorbar()

plt.subplot(1,3,1); plt.imshow( dat[int(nx*0.8),:,:], vmin=-1.0, vmax=1.0 ); plt.colorbar()
plt.subplot(1,3,2); plt.imshow( dat[:,int(nx*0.0),:], vmin=-1.0, vmax=1.0 ); plt.colorbar()
plt.subplot(1,3,3); plt.imshow( dat[:,:,int(nx*0.0)], vmin=-1.0, vmax=1.0 ); plt.colorbar()

plt.show()