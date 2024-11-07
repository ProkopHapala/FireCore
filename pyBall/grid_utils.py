import numpy as np
import matplotlib.pyplot as plt

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
    plt.title(label)    
