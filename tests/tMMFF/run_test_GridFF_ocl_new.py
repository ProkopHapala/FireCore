import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")

from pyBall.tests import ocl_GridFF_new as gff


# ==== Just plot basis
# t = np.linspace(0.0, 1.0, 10)
# ws = gff.Bspline_basis5(t)
# plt.plot(ws)
# plt.figure()
# plt.imshow( ws )
# plt.show()
# exit()


# =============  Functions

R0 = 3.5
E0 = 0.1
a  = 1.8

Q = 0.4
p0 = [-2.0,-2.0,0.0]


#name="NaCl_1x1_L2"

#name="NaCl_1x1_L3"
#name="NaCl_8x8_L3"
#name="NaCl_8x8_L3_NaHole"
#name="NaCl_8x8_L3_ClHole"
#name="NaCl_8x8_L3_NaClHole"
name="NaCl_8x8_L3_step"

#mol_name="PTCDA.xyz"
#gff.test_gridFF_ocl( fname="data/xyz/NaCl_1x1_L2.xyz" )
#gff.test_gridFF_ocl( fname="data/xyz/"+name+".xyz", save_name="double3", bMorse=True, bEwald=False  )
#gff.test_gridFF_ocl( fname="/home/prokop/git/FireCore/tests/pyutils/NaCl_8x8_L3.xyz" )

gff.test_gridFF_ocl( fname="data/xyz/"+name+".xyz", save_name="double3", job="PLQ" )




# name="NaCl_1x1_L3"
# print("# =================== \n\n\n\ NOW ", name, " =================== \n")
# gff.test_gridFF_ocl( fname="data/xyz/"+name+".xyz", save_name="double3", job="PLQ" )

#gff.test_gridFF_ocl( fname="data/xyz/"+name+".xyz", save_name="double3", job="PLQ_lin" )

# PLQ = np.load("./data/"+name+"/Bspline_PLQd_ocl.npy")
# VPaul = PLQ[:,:,:,0]
# VLond = PLQ[:,:,:,1]
# VCoul = PLQ[:,:,:,2]
# print( "Min,Max VPaul ", VPaul.min(), VPaul.max(), " VLond ", VLond.min(), VLond.max(), " VCoul ", VCoul.min(), VCoul.max() )
# plt.figure( figsize=(15,5)); 
# plt.subplot(1,3,1); plt.imshow( VPaul[:,:,1] ); plt.colorbar(); plt.title("VPaul fit GPU")
# plt.subplot(1,3,2); plt.imshow( VLond[:,:,1] ); plt.colorbar(); plt.title("VLond fit GPU")
# plt.subplot(1,3,3); plt.imshow( VCoul[:,:,1] ); plt.colorbar(); plt.title("VCoul fit GPU")
# plt.show()

# ======== Ewald

d=0.6
apos=np.array([
    [-d,.0,0.],
    [+d,.0,0.],
    [0.,-d,0.],
    [0.,+d,0.],
])
qs = [ +1.,+1.,-1.,-1. ]

# d=0.6
# apos=np.array([
#     [0.,.0,-d],
#     [0.,.0,+d],
# ])
# qs = [ +1.,-1. ]

#gff.test_Ewald( apos, qs,  ns=(100,100,100), dg=(0.10,0.10,0.10), order=3, bSlab=True,  nPBC=(100,100,0), bOld=True )
#gff.test_Ewald( apos, qs,  ns=(100,100,100), dg=(0.10,0.10,0.10), order=3, bSlab=True,  nPBC=(100,100,0) )
# gff.test_Ewald( apos, qs,  ns=[100,100,150], dg=[0.10,0.10,0.10], order=3, bSlab=True,  nPBC=[100,100,0] )
# gff.test_Ewald( apos, qs,  ns=[100,100,200], dg=[0.10,0.10,0.10], order=3, bSlab=True,  nPBC=[100,100,0] )
# gff.test_Ewald( apos, qs,  ns=[100,100,300], dg=[0.10,0.10,0.10], order=3, bSlab=True,  nPBC=[100,100,0] )
# gff.test_Ewald( apos, qs,  ns=[100,100,400], dg=[0.10,0.10,0.10], order=3, bSlab=True,  nPBC=[100,100,0] )

plt.show()