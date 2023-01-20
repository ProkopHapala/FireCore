import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

sys.path.append("../../")
from pyBall import SchroedingerGreen2D as sch

#sch.init(10,10)
sch.init(100,100)
#sch.init(1000,1000)
sch.getBuffs()

nx,ny = sch.V.shape

def prepareFucntions( extent=[ -5,5, -5,5 ], width =0.5, strenght=0.1, x0=0.,y0=0. ):
    xs    = np.linspace( extent[0], extent[1], nx ) #- nx*0.5
    ys    = np.linspace( extent[2], extent[3], ny ) #- ny*0.5
    Xs,Ys = np.meshgrid(xs,ys)
    #----- source
    sch.source[:,:] = np.exp( -((Xs-x0)**2 + (Ys-y0)**2)/(width*width)   ) * strenght
    # ---- Potential
    #sch.V[:,:] = (Xs**2 + Ys**2)*(16.0/100.0)
    wV = 3.0
    sch.V[:,:]  = np.exp( -((Xs-0)**2 + (Ys-0)**2)/(wV*wV)   ) * -1.0  + 1.0
    sch.V[:,:] += np.sin(Xs)**2*np.cos(Ys*3)**2* -0.1
    sch.V[:,:] *= 10.0
    return sch.V, sch.source

prepareFucntions( x0=1., y0=1. )      #;print( "V\n", V )


#sch.V[:,:] = 0

#sch.psi[:,:] = 1
sch.psi[:,:] = 0
#sch.psi[:,:] = 0.01; sch.psi[nx//2,ny//2] = 1.
#sch.source[:,:] = 0;   sch.source[nx//2,ny//2] = -1.

nstep=20
plt.figure(figsize=(3*nstep,6))
Es=[]
Fs=[]

#E0=0;
E0=4;
#sch.EQF[3] = 2.0   #   set Energy
sch.EQF[3] = E0   #   set Energy

'''
perView = 50
#perView = 1
for i in range(nstep):
    for j in range(perView):
        F2=sch.step( E0 = E0, dt=0.3 )               #;print( "fpsi \n", sch.fpsi ); 
    ij = i*perView 
    Es.append(sch.EQF[0])
    Fs.append(np.sqrt(F2))
    plt.subplot(2,nstep,i+1      ); plt.imshow(sch.psi ); plt.title("Psi_%i" %(ij+1)     ); plt.colorbar()
    plt.subplot(2,nstep,i+nstep+1); plt.imshow(sch.fpsi); plt.title("dPsi/dt_%i" %(ij+1) ); plt.colorbar()
    plt.figure(); plt.plot(Es,label="E"); plt.plot(Fs,label="|Ferr|"); plt.axhline(0,ls='--',c='k'); plt.legend(); plt.grid()
'''

perView = 100
for i in range(nstep):
    for j in range(perView):
        F2=sch.step_Green( )   #;print(F2)  #;print( "fpsi \n", sch.fpsi ); 
        print(i,j,F2)
        f = np.sqrt(F2)
        Fs.append( f )
        if f<1e-4: break
    ij = i*perView 
    plt.subplot(2,nstep,i+1      ); plt.imshow(sch.psi ); plt.title("Psi_%i"     %(ij+1) ); plt.colorbar()
    if f<1e-4: break
#print(Fs)
plt.figure();plt.plot(Fs,label="|Ferr|"); plt.axhline(0,ls='--',c='k'); plt.yscale('log'); plt.legend(); plt.grid()

vmin=sch.psi.min()
vmax=sch.psi.max()
vmax=max(-vmin,vmax)

plt.figure(figsize=(20,5))
plt.subplot(1,4,1); plt.imshow(sch.V);      plt.title("V "     );  plt.colorbar()
plt.subplot(1,4,2); plt.imshow(sch.source); plt.title("source" );  plt.colorbar()
plt.subplot(1,4,3); plt.imshow(sch.psi, cmap='seismic', vmin=-vmax,vmax=vmax);    plt.title("Psi"  );  plt.colorbar()
plt.subplot(1,4,4); plt.imshow(sch.psi**2); plt.title("|Psi|^2"  );  plt.colorbar()

#plt.axis('equal')
plt.show(  )


