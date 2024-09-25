import numpy as np
import matplotlib.pyplot as plt

# ========= Functions

def cut(dat,f=[0.8,0.5,0.5],iax=0, n=1, bPlot=True, label="", ls='-', lw='0.5'):
    sh = dat.shape
    if (   iax==0 ):
        ys = dat[:,int(sh[1]*f[1]),int(sh[2]*f[2])]
    elif ( iax==1 ):
        ys = dat[int(sh[0]*f[0]),:,int(sh[2]*f[2])]
    elif( iax==2 ):
        ys = dat[int(sh[0]*f[0]),int(sh[1]*f[1]),:]
    if n>1: ys = np.tile( ys, n )
    if bPlot:  plt.plot( ys, ls, lw=lw, label=label+str(iax)+" f="+str(f) )
    return ys

def cut2(dat,f=0.5,iax=0, n=[1,1], bPlot=True, vmax=0.1 ):
    ni = dat.shape[iax]
    if ( iax==0 ):
        ys = dat[int(ni*f),:,:]
    elif ( iax==1 ):
        ys = dat[:,int(ni*f),:]
    elif( iax==2 ):
        ys = dat[:,:,int(ni*f)]
    #print("ys.shape ", iax, ys.shape)
    if (n[0]>1) or (n[1]>1): ys = np.tile( ys, n )
    if bPlot: plt.imshow( ys, origin='lower', vmin=-vmax,vmax=vmax, cmap='bwr' )
    return ys

def cut1d_xyz( dat, f=[0.8,0.5,0.5],n=2, iaxs=[0,1,2], label=""):
    for iax in iaxs: cut(dat,f=f,iax=iax, n=n, bPlot=True, label=label)
    # cut(dat,f=f,iax=0, n=n, bPlot=True )
    # cut(dat,f=f,iax=1, n=n, bPlot=True)
    # cut(dat,f=f,iax=2, n=n, bPlot=True)
    plt.legend()
    plt.grid()
    plt.title(name)    

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
#     cut1d_xyz(dat, f=f , n=2, iaxs=[1,2] )
#     cut1d_xyz(dat, f=f2, n=2, iaxs=[1,2] )
#     cut1d_xyz(dat, f=f3, n=2, iaxs=[1,2] )
#     #plt.figure(fig2.number); plt.subplot(2,3,i+1);  cut2(dat,f=f[iax],iax=iax, n=[0,0], bPlot=True); plt.title(name)  



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