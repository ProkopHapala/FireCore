import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import MMFF as mmff
from pyBall.tests import GridFF as gff
from pyBall.tests import utils as ut

# =============  Functions

def load_potential(filename):
    """Load potential components from a .npy file."""
    data = np.load(filename)
    VPaul, VLond, VCoul = data[..., 0], data[..., 1], data[..., 2]
    return VPaul, VLond, VCoul

def plot_potentials(V_cpu, V_gpu, slice_idx=4):
    """Plot CPU vs GPU potentials for each component."""
    components = ['VPaul', 'VLond', 'VCoul']
    plt.figure(figsize=(15, 5))
    for i, comp in enumerate(components, 1):

        Vcpui = V_cpu[i-1][:, :, slice_idx]
        Vgpui = V_gpu[i-1][:, :, slice_idx]

        plt.subplot(3, 3, i)
        plt.imshow( Vcpui, origin='lower')
        plt.colorbar()
        plt.title(f"{comp} CPU")
        
        plt.subplot(3, 3, i+3)
        plt.imshow( Vgpui, origin='lower')
        plt.colorbar()
        plt.title(f"{comp} GPU")

        plt.subplot(3, 3, i+6)
        plt.imshow( Vgpui-Vcpui, origin='lower')
        plt.colorbar()
        plt.title(f"{comp} GPU")
    plt.tight_layout()
    #plt.show()

def plot_1D_comparison(V_cpu, V_gpu, PLQH, ix=0, iy=0, iz=4, ylim=(-1.0, 1.0)):
    """Plot 1D combined potentials for CPU and GPU."""
    Ecpu_z = V_cpu[0][ix, iy, :] * PLQH[0] + V_cpu[1][ix, iy, :] * PLQH[1] + V_cpu[2][ix, iy, :] * PLQH[2]
    Egpu_z = V_gpu[0][ix, iy, :] * PLQH[0] + V_gpu[1][ix, iy, :] * PLQH[1] + V_gpu[2][ix, iy, :] * PLQH[2]

    Ecpu_x = V_cpu[0][:, iy, iz] * PLQH[0] + V_cpu[1][:, iy, iz] * PLQH[1] + V_cpu[2][:, iy, iz] * PLQH[2]
    Egpu_x = V_gpu[0][:, iy, iz] * PLQH[0] + V_gpu[1][:, iy, iz] * PLQH[1] + V_gpu[2][:, iy, iz] * PLQH[2]

    Ecpu_y = V_cpu[0][ix, :, iz] * PLQH[0] + V_cpu[1][ix, :, iz] * PLQH[1] + V_cpu[2][ix, :, iz] * PLQH[2]
    Egpu_y = V_gpu[0][ix, :, iz] * PLQH[0] + V_gpu[1][ix, :, iz] * PLQH[1] + V_gpu[2][ix, :, iz] * PLQH[2]

    plt.figure(figsize=(10, 5))
    #plt.plot(combined_cpu[4:], ':', lw=1.5, label="V_PLQ CPU "+str(PLQH))   # NOTE: TODO: We see that GPU is shifted by 4 voxels with respect to GPU
    plt.plot(Ecpu_z, ':', lw=1.5, label="V_PLQ(z) CPU "+str(PLQH)) 
    plt.plot(Egpu_z, '-', lw=0.5, label="V_PLQ(z) GPU "+str(PLQH))

    plt.plot(Ecpu_x, ':', lw=1.5, label="V_PLQ(x) CPU "+str(PLQH)) 
    plt.plot(Egpu_x, '-', lw=0.5, label="V_PLQ(x) GPU "+str(PLQH))

    plt.plot(Ecpu_y, ':', lw=1.5, label="V_PLQ(y) CPU "+str(PLQH)) 
    plt.plot(Egpu_y, '-', lw=0.5, label="V_PLQ(y) GPU "+str(PLQH))

    #plt.ylim(ylim)
    plt.grid(True)
    plt.legend()
    plt.xlabel('Index (Z-axis)')
    plt.ylabel('Combined Potential')
    plt.title('1D Potential Slice Comparison')
    #plt.show()

def compare_potentials( name="NaCl_1x1_L3", R0=3.5, E0=1.0, a=-1.5, Q=0.0 ):
    # File paths
    cpu_file = "./data/"+name+"/Bspline_PLQd.npy"
    gpu_file = "./data/"+name+"/Bspline_PLQd_ocl.npy"
    print( f"Comparing CPU and GPU potentials for {name}." )
    # Load potentials
    V_cpu = load_potential(cpu_file)
    V_gpu = load_potential(gpu_file)

    plot_potentials(V_cpu, V_gpu, slice_idx=1)

    PLQH = ut.getPLQH(R0=R0, E0=E0, Q=Q, H=0, a=a )  # Replace with actual parameters

    plot_1D_comparison(V_cpu, V_gpu, PLQH, ix=0, iy=0)


# =============  Main

#R0 = 3.5; E0 = 0.1
R0 = 1.443; E0=0.00190802
a  = 1.5
#Q = 0.0
Q = 0.4

#p0 = [1.0,-5.05,2.0]
#p0 = [0.0,0.0,2.0]
p0 = [-2.0,-2.0,0.0]


# Qcpu = np.load("./data/NaCl_1x1_L3/Qgrid.npy")
# Qgpu = np.load("./data/NaCl_1x1_L3/Qgrid_ocl.npy")
# plt.figure(figsize=(15,5))
# print( "Qcpu min,max ", Qcpu.min(), Qcpu.max() )
# print( "Qgpu min,max ", Qgpu.min(), Qgpu.max() )
# # where is the maximum, in 3d ?
# argmax_cpu = np.unravel_index(np.argmax(Qcpu), Qcpu.shape)
# argmax_gpu = np.unravel_index(np.argmax(Qgpu), Qgpu.shape)
# print( "argmax(Qcpu) ", argmax_cpu, "\nargmax(Qgpu) ", argmax_gpu )
# plt.subplot(1,3,1); plt.imshow( Qcpu[300:,:,20] ); plt.colorbar(); plt.title("CPU Qgrid")
# plt.subplot(1,3,2); plt.imshow( Qgpu[300:,:,20] ); plt.colorbar(); plt.title("GPU Qgrid")
# #plt.subplot(1,3,3); plt.imshow( Qgpu[:,:,0]-Qcpu[:,:,0] ); plt.colorbar(); plt.title("CPU-GPU Qgrid")
# plt.show()
# exit()

#compare_potentials( name="NaCl_1x1_L3", R0=R0, E0=E0, a=a, Q=Q ); plt.show(); exit()

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

tmax=10.0
tmin=-5.0
Q    = +0.4
R0 = 1.443; E0=np.sqrt(0.00190802)
a  = 1.5


#gff.test_gridFF    ( mode=6, name="data/xyz/NaCl_1x1_L2", p0=[0.0,0.0,2.0],  Q=0.0, E0=0.1, bRefine=False, nPBC=[5,5,0], Emax=Emax, Fmax=Fmax )
#gff.test_gridFF    ( mode=6, name="data/xyz/NaCl_1x1_L3", p0=[0.0,0.0,2.0],  Q=0.0, E0=0.1, bRefine=False, nPBC=[5,5,0], Emax=Emax, Fmax=Fmax )

# -- above Cl
# gff.test_gridFF    ( mode=6, name="data/xyz/Cl_surf", p0=[0.0,0.0,2.0],  Q=0.0, E0=E0, R0=R0, tmin=tmin,tmax=tmax, bRefine=False, nPBC=[0,0,0], Emax=Emax, Fmax=Fmax, title="just Morse" )
# gff.test_gridFF    ( mode=6, name="data/xyz/Cl_surf", p0=[0.0,0.0,2.0],  Q=Q,   E0=E0, R0=R0, tmin=tmin,tmax=tmax, bRefine=False, nPBC=[0,0,0], Emax=Emax, Fmax=Fmax, title="Morse+Coulomb" )
# gff.test_gridFF    ( mode=6, name="data/xyz/Cl_surf", p0=[0.0,0.0,2.0],  Q=Q,   E0=0.0, R0=R0,tmin=tmin,tmax=tmax, bRefine=False, nPBC=[0,0,0], Emax=Emax, Fmax=Fmax, title="just Coulomb" )

# -- above Cl
gff.test_gridFF    ( mode=6, name="data/xyz/NaCl_1x1_L3", p0=[2.0,2.0,2.0],  Q=0.0, E0=E0, R0=R0, tmin=tmin,tmax=tmax, bRefine=False, nPBC=[5,5,0], Emax=Emax, Fmax=Fmax, title="just Morse" )
gff.test_gridFF    ( mode=6, name="data/xyz/NaCl_1x1_L3", p0=[2.0,2.0,2.0],  Q=Q,   E0=E0, R0=R0, tmin=tmin,tmax=tmax, bRefine=False, nPBC=[2000,10000,0], Emax=Emax, Fmax=Fmax, title="Morse+Coulomb" )
gff.test_gridFF    ( mode=6, name="data/xyz/NaCl_1x1_L3", p0=[2.0,2.0,2.0],  Q=Q,   E0=0.0, R0=R0,tmin=tmin,tmax=tmax, bRefine=False, nPBC=[2000,10000,0], Emax=Emax, Fmax=Fmax, title="just Coulomb" )

# -- above Na
#gff.test_gridFF    ( mode=6, name="data/xyz/NaCl_1x1_L3", p0=[0.0,0.0,2.0],  Q=0.0, E0=E0, R0=R0, tmin=tmin,tmax=tmax, bRefine=False, nPBC=[5,5,0], Emax=Emax, Fmax=Fmax, title="just Morse" )
#gff.test_gridFF    ( mode=6, name="data/xyz/NaCl_1x1_L3", p0=[0.0,0.0,2.0],  Q=Q,   E0=E0, R0=R0, tmin=tmin,tmax=tmax, bRefine=False, nPBC=[200,200,0], Emax=Emax, Fmax=Fmax, title="Morse+Coulomb" )
#gff.test_gridFF    ( mode=6, name="data/xyz/NaCl_1x1_L3", p0=[0.0,0.0,2.0],  Q=Q,   E0=0.0, R0=R0,tmin=tmin,tmax=tmax, bRefine=False, nPBC=[200,200,0], Emax=Emax, Fmax=Fmax, title="just Coulomb" )
#gff.test_gridFF    ( mode=6, name="data/xyz/NaCl_1x1_L3", p0=[0.0,0.0,2.0],  Q=Q, E0=E0, R0=R0, bRefine=False, nPBC=[5,5,0], Emax=Emax, Fmax=Fmax )
#gff.test_gridFF    ( mode=6, name="data/xyz/NaCl_8x8_L3", p0=[2.0,2.0,2.0],  Q=0.0, E0=0.1, bRefine=False, nPBC=[5,5,0], Emax=Emax, Fmax=Fmax )

#gff.test_gridFF    ( mode=6, title="Bspline_o3 \n(z-cut)" , p0=[0.0,0.0,2.0],  Q=0.0, E0=0.1, bRefine=False, nPBC=[5,5,0], Emax=Emax, Fmax=Fmax )
#gff.test_gridFF    ( mode=6, title="Bspline_o3 \n(z-cut)" , p0=[2.0,2.0,2.0],  Q=0.0, E0=0.1, bRefine=False, nPBC=[5,5,0], Emax=Emax, Fmax=Fmax )
#gff.test_gridFF_lat( mode=6, title="Bspline_o3 \n(lat iax=0)", p0=[0.0,0.0,2.0], Q=0.0, E0=0.1 )
#gff.test_gridFF_lat( mode=6, title="Bspline_o3 \n(lat iax=0)", p0=[2.0,2.0,2.0], Q=0.0, E0=0.1 )

#gff.test_gridFF_npy( ps_xy=[(0.0,0.0),(0.0,0.5),(0.5,0.0),(0.5,0.5)], mode=6, title="" )

#gff.test_gridFF_lat( mode=6, name="data/xyz/NaCl_1x1_L3", p0=[0.0,0.0,0.5], Q=0.0, E0=0.1, bRefine=False, nPBC=[5,5,0] )

#gff.test_gridFF_npy( ps_xy=[(0.0,0.0)], mode=6, title="" )
#gff.test_gridFF_npy_lat( ps_zy=[(0.0,0.0)], mode=6, title="" )
#gff.test_gridFF_npy_lat( ps_zy=[(0.1,0.1)], mode=6, title="" )

plt.show()