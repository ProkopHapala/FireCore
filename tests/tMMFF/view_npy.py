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

def cut1d_xyz( dat, f=[0.8,0.5,0.5],n=2):
    cut(dat,f=f,iax=0, n=n, bPlot=True )
    cut(dat,f=f,iax=1, n=n, bPlot=True)
    cut(dat,f=f,iax=2, n=n, bPlot=True)
    plt.legend()
    plt.grid()
    plt.title(name)    

# ========= MAIN

path="data/NaCl_1x1_L2/"

files=[
"debug_BsplineCoul_pbc.npy",
"debug_BsplineLond_pbc.npy",
"debug_BsplinePaul_pbc.npy",
"debug_VCoul_pbc.npy",
"debug_VLond_pbc.npy",
"debug_VPaul_pbc.npy",    
]


f=[0.8,0.5,0.5]
iax = 1

fig1 = plt.figure(1,figsize=(15,10))
fig2 = plt.figure(2,figsize=(15,10))

for i,name in enumerate(files):
    dat = np.load(path+name)
    print( "name ", name, dat.shape )
    plt.figure(fig1.number); plt.subplot(2,3,i+1);  cut1d_xyz(dat, f=f, n=2 )
    plt.figure(fig2.number); plt.subplot(2,3,i+1);  cut2(dat,f=f[iax],iax=iax, n=[0,0], bPlot=True); plt.title(name)  


plt.show()