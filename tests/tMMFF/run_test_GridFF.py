import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import MMFF as mmff
from pyBall.tests import GridFF as gff
from pyBall.tests import utils as ut

# =============  Functions

R0 = 3.5
E0 = 1.0
a  = 1.8

#Q = 0.0
Q = 0.4
#p0 = [1.0,-5.05,2.0]
#p0 = [0.0,0.0,2.0]
p0 = [-2.0,-2.0,0.0]

mmff.initParams()

#gff.test_gridFF    ( mode=6, title="Bspline_o3 \n(z-cut)" ,    Q=0.0, )
#gff.test_gridFF_lat( mode=6, title="Bspline_o3 \n(lat iax=0)", Q=0.0, )
#gff.test_gridFF_lat( mode=6, title="Bspline_o3 \n(lat iax=0)", Q=0.0, p0=p0, iax=0 )

#gff.test_gridFF    ( mode=1, title="tri-linar force \n(z-cut)"          )
#gff.test_gridFF_lat( mode=1, title="tri-Linear Force", Q=0.0, p0=p0, iax=0 )

#Emax=0.01 
#Fmax=0.01

Emax=0.00001 
Fmax=0.00001

#gff.test_gridFF    ( mode=6, title="Bspline_o3 \n(z-cut)" , p0=[0.0,0.0,2.0],  Q=0.4, E0=0, bRefine=False, nPBC=[50 ,50 ,0], Emax=Emax, Fmax=Fmax )
#gff.test_gridFF    ( mode=6, title="Bspline_o3 \n(z-cut)" , p0=[0.0,0.0,2.0],  Q=0.4, E0=0, bRefine=False, nPBC=[100,100,0], Emax=Emax, Fmax=Fmax )
#gff.test_gridFF    ( mode=6, title="Bspline_o3 \n(z-cut)" , p0=[0.0,0.0,2.0],  Q=0.4, E0=0, bRefine=False, nPBC=[150,150,0], Emax=Emax, Fmax=Fmax )

#gff.test_gridFF    ( mode=6, title="Bspline_o3 \n(z-cut)" , p0=[0.0,0.0,2.0],  Q=0.4, E0=0, bRefine=False, nPBC=[100,100,0], Emax=Emax, Fmax=Fmax )
#gff.test_gridFF    ( mode=6, title="Bspline_o3 \n(z-cut)" , p0=[0.0,0.0,2.0],  Q=0.4, E0=0, bRefine=False, nPBC=[300,300,0], Emax=Emax, Fmax=Fmax )
#gff.test_gridFF    ( mode=6, title="Bspline_o3 \n(z-cut)" , p0=[0.0,0.0,2.0],  Q=0.4, E0=0, bRefine=False, nPBC=[400,400,0], Emax=Emax, Fmax=Fmax )
#gff.test_gridFF    ( mode=6, title="Bspline_o3 \n(z-cut)" , p0=[2.0,2.0,2.0],  Q=0.4, E0=0, bRefine=False, nPBC=[400,400,0], Emax=Emax, Fmax=Fmax )

#gff.test_gridFF    ( name="data/xyz/NaCl_8x8_L3", mode=6, title="Bspline_o3 \n(z-cut)" , p0=[0.0,0.0,2.0],  Q=0.0, E0=1.0, bRefine=False, nPBC=[400,400,0], Emax=Emax, Fmax=Fmax )

Emax=0.01 
Fmax=0.01

#gff.test_gridFF_lat( mode=6, title="Bspline_o3 \n(lat iax=0)", p0=[0.0,0.0,2.0], nPBC=[400,400,0],  Q=0.4, E0=0.0, Emax=Emax, Fmax=Fmax  )
#gff.test_gridFF_lat( mode=6, title="Bspline_o3 \n(lat iax=0)", p0=[2.0,2.0,2.0], nPBC=[400,400,0],  Q=0.4, E0=0.0, Emax=Emax, Fmax=Fmax  )

Emax=0.1 
Fmax=0.1

# plt.figure(figsize=(10,5))
# PLQ = np.load("./data/NaCl_1x1_L2/Bspline_PLQd.npy"); print("PLQ.shape ", PLQ.shape )
# VPaul = PLQ[:,:,:,0]
# VLond = PLQ[:,:,:,1]
# VCoul = PLQ[:,:,:,2]
# print( "Bspline_PLQd Min,Max VPaul ", VPaul.min(), VPaul.max(), " VLond ", VLond.min(), VLond.max(), " VCoul ", VCoul.min(), VCoul.max() )
# plt.subplot(1,3,1); plt.imshow( VPaul[:,:,1] ); plt.colorbar(); plt.title("VPaul CPU ")
# #plt.subplot(1,3,1); plt.imshow( VPaul[1,:,:] ); plt.colorbar(); plt.title("VPaul CPU ")
# #plt.show()
# PLQ_ = np.load("./data/NaCl_1x1_L2/Bspline_PLQd_ocl.npy"); print("PLQ_.shape ", PLQ_.shape )
# VPaul_ = PLQ_[:,:,:,0]
# VLond_ = PLQ_[:,:,:,1]
# VCoul_ = PLQ_[:,:,:,2]
# print( "Bspline_PLQd_ocl Min,Max VPaul ", VPaul_.min(), VPaul_.max(), " VLond ", VLond_.min(), VLond_.max(), " VCoul ", VCoul_.min(), VCoul_.max() )
# plt.subplot(1,3,2); plt.imshow( VPaul_[:,:,1] ); plt.colorbar(); plt.title("VPaul GPU ")
# #plt.subplot(1,3,2); plt.imshow( VPaul_[1,:,:] ); plt.colorbar(); plt.title("VPaul GPU ")
# #plt.subplot(1,3,3); plt.imshow( VPaul_[:,:,1]-VPaul[:,:,1]  ); plt.colorbar(); plt.title("VPaul CPU-GPU ")
# plt.show()

# PLQ_cpu = np.load("./data/NaCl_1x1_L2/Bspline_PLQd.npy"); print("PLQ_cpu.shape ", PLQ_cpu.shape )
# PLQ_gpu = np.load("./data/NaCl_1x1_L2/Bspline_PLQd_ocl.npy"); print("PLQ_gpu.shape ", PLQ_gpu.shape )
# ix,iy=0,0
# # plt.plot( PLQ_cpu[ix,iy,:, 0],  label="V_Paul CPU" )
# # plt.plot( PLQ_gpu[ix,iy,:, 0],  label="V_Paul GPU" )
# # plt.plot( PLQ_cpu[ix,iy,:, 1],  label="V_Lond CPU" )
# # plt.plot( PLQ_gpu[ix,iy,:, 1],  label="V_Lond GPU" )
# PLQH = ut.getPLQH( R0, E0, a, Q, 0 )
# plt.plot( PLQ_cpu[ix,iy,:, 0]*PLQH[0] + PLQ_cpu[ix,iy,:, 1]*PLQH[1],    ':',    lw=1.5, label="V_Morse CPU" )
# plt.plot( PLQ_gpu[ix,iy,:, 0]*PLQH[0] + PLQ_gpu[ix,iy,:, 1]*PLQH[1],    '-',    lw=0.5, label="V_Morse GPU" )
# plt.ylim(-1.0,1.0)
# plt.grid()
# plt.legend()
# plt.show()

gff.test_gridFF    ( mode=6, name="data/xyz/NaCl_1x1_L2", p0=[0.0,0.0,2.0],  Q=0.0, E0=0.1, bRefine=False, nPBC=[5,5,0], Emax=Emax, Fmax=Fmax )
#gff.test_gridFF    ( mode=6, name="data/xyz/NaCl_8x8_L3", p0=[2.0,2.0,2.0],  Q=0.0, E0=0.1, bRefine=False, nPBC=[5,5,0], Emax=Emax, Fmax=Fmax )

#gff.test_gridFF    ( mode=6, title="Bspline_o3 \n(z-cut)" , p0=[0.0,0.0,2.0],  Q=0.0, E0=0.1, bRefine=False, nPBC=[5,5,0], Emax=Emax, Fmax=Fmax )
#gff.test_gridFF    ( mode=6, title="Bspline_o3 \n(z-cut)" , p0=[2.0,2.0,2.0],  Q=0.0, E0=0.1, bRefine=False, nPBC=[5,5,0], Emax=Emax, Fmax=Fmax )
#gff.test_gridFF_lat( mode=6, title="Bspline_o3 \n(lat iax=0)", p0=[0.0,0.0,2.0], Q=0.0, E0=0.1 )
#gff.test_gridFF_lat( mode=6, title="Bspline_o3 \n(lat iax=0)", p0=[2.0,2.0,2.0], Q=0.0, E0=0.1 )

#gff.test_gridFF_npy( ps_xy=[(0.0,0.0),(0.0,0.5),(0.5,0.0),(0.5,0.5)], mode=6, title="" )

#gff.test_gridFF_npy( ps_xy=[(0.0,0.0)], mode=6, title="" )
#gff.test_gridFF_npy_lat( ps_zy=[(0.0,0.0)], mode=6, title="" )
#gff.test_gridFF_npy_lat( ps_zy=[(0.1,0.1)], mode=6, title="" )

plt.show()