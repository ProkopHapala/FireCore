#!/usr/bin/python
import sys
import numpy as np
import matplotlib.pyplot as plt


fname1=sys.argv[1]
fname2=sys.argv[2]

dat1 = np.genfromtxt(fname1); print("dat1.shape ", dat1.shape); dat1=dat1.reshape((-1,dat1.shape[1]//2,2))
dat2 = np.genfromtxt(fname2); dat2=dat2.reshape((-1,dat2.shape[1]//2,2))

plt.figure(figsize=(15,5))
plt.subplot(1,3,1); plt.imshow( dat1[:,:,0]             ); plt.colorbar(); plt.title(fname1)
plt.subplot(1,3,2); plt.imshow( dat2[:,:,0]             ); plt.colorbar(); plt.title(fname2)
plt.subplot(1,3,3); plt.imshow( dat1[:,:,0]-dat2[:,:,0] ); plt.colorbar(); plt.title("difference")
plt.show()
