import numpy as np
import matplotlib.pyplot as plt

import os
import sys
sys.path.append("../../")
from pyBall import grid_utils as gu

# ========= Functions

def normalpath(path):
    return os.path.abspath(os.path.normpath(os.path.expanduser(path)))

# # normalize path (to remove ~)
# def normalize(path):
#     if path[0]=='~':
#         path=path.replace('~', os.environ['HOME'])
#     return path

fname1 = normalpath( "~/git/FireCore/tests/tMMFF/data/NaCl_1x1_L3/Bspline_PLQd.npy" )
fname2 = normalpath( "~/git/FireCore/tests/tMMFF/data/NaCl_1x1_L3/Bspline_PLQd_ocl.npy" )
#fname2 = normalpath( "~/git/FireCore/tests/tMMFF/data/NaCl_8x8_L3/Bspline_PLQd.npy" )


dat1 = np.load(fname1)  ;print("dat1.shape", dat1.shape)
dat2 = np.load(fname2)  ;print("dat2.shape", dat2.shape)

nz1,ny1,nx1,_ = dat1.shape
nz2,ny2,nx2,_ = dat2.shape

iChan=0
cx=0.0
cy=0.0

# ----- 1D plot
plt.figure(figsize=(15,5))
channame=['Pauli', 'London', 'Coulomb']
for iChan in range(3):
    print(iChan)
    plt.subplot(1,3,iChan+1);
    plt.plot( dat1[int(cx*nx1),int(cy*ny1),:,iChan], label=channame[iChan]+' NaCl_1x1_L3' );
    plt.plot( dat2[int(cx*nx2),int(cy*ny2),:,iChan], label=channame[iChan]+' NaCl_1x1_L3_ocl' );
    #plt.plot( dat2[int(cx*nx2),int(cy*ny2),:,iChan], label=channame[iChan]+' NaCl_8x8_L3' );

    #plt.plot( dat1[int(cx*nx1),:,5,iChan], label=channame[iChan]+' NaCl_1x1_L3' );
    #plt.plot( dat2[int(cx*nx2),:,5,iChan], label=channame[iChan]+' NaCl_1x1_L3 ocl' );
    #plt.plot( dat2[int(cx*nx2),:,5,iChan], label=channame[iChan]+' NaCl_8x8_L3' );
    plt.legend()
    plt.grid()
plt.title(fname1+"\n"+fname2)
plt.show()

# ----- 2D plot
# plt.figure(figsize=(10,5))
# plt.subplot(1,2,1); plt.imshow( dat1[int(cx*nx1),:,:,iChan], ); plt.title( 'NaCl_1x1_L3' );
# plt.subplot(1,2,2); plt.imshow( dat2[int(cx*nx2),:,:,iChan], ); plt.title( 'NaCl_8x8_L3' );
# plt.show()