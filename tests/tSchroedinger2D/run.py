import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

hbar_SI = 1.054571817e-34;    # [J*s]
Me_SI   = 9.1093837e-31;      # [kg ]
eV_SI   = 1.602176634e-19;    # [J  ]

pref      = (hbar_SI*hbar_SI)/(2*Me_SI)    ;print(pref)
pref_eV   = pref/eV_SI                     ;print(pref_eV)
pref_eVA2 = pref/(eV_SI*(1e-10**2))        ;print(pref_eVA2)

#exit()


sys.path.append("../../")
from pyBall import SchroedingerGreen2D as sch

extent=[ -100,100, -100,100 ]
L=(extent[1]-extent[0])

#sch.init(10,10, L=L )
sch.init(50,50, L=L)
#sch.init(100,100, L=L)
#sch.init(1000,1000, L=L)
sch.getBuffs()

#setStep(dstep=0.1,m_Me=1.0, L=None)

nx,ny = sch.V.shape


def gauss( w, x0=0., y0=0. ):
    return np.exp( -((Xs-x0)**2 + (Ys-y0)**2)/(w*w) )

def prepareFucntions( extent=[ -5,5, -5,5 ], width =10.0, strenght=0.1, x0=0.,y0=0. ):
    global Xs,Ys
    #sch.setStep(L=extent[1]-extent[0])
    xs    = np.linspace( extent[0], extent[1], nx ) #- nx*0.5
    ys    = np.linspace( extent[2], extent[3], ny ) #- ny*0.5
    Xs,Ys = np.meshgrid(xs,ys)
    #----- source
    sch.source[:,:] = gauss( width, x0, y0 )
    # ---- Potential
    #sch.V[:,:] = (Xs**2 + Ys**2)*(16.0/100.0)
    wV = 40.0
    sch.V[:,:]  = (1. - gauss( 40.0 )   ) * 0.05
    sch.V[:,:] += (1. - np.sin(Xs*(2*np.pi/50.0))**2*np.cos(Ys*(2*np.pi/100.0))**2 ) * 0.01
    print( "V min,max: ", sch.V.min(), sch.V.max() )
    return sch.V, sch.source

prepareFucntions( extent=extent, x0=10., y0=10. )      #;print( "V\n", V )

#sch.V[:,:] = 0

sch.psi[:,:] = 1
#sch.psi[:,:] = 0
#sch.psi[:,:] = 0.01; sch.psi[nx//2,ny//2] = 1.
#sch.source[:,:] = 0;   sch.source[nx//2,ny//2] = -1.


#E0=0;
E0=0.1
#sch.EQF[3] = 2.0   #   set Energy
sch.EQF[3] = E0   #   set Energy



nstep=20
plt.figure(figsize=(5*nstep,3*5))
Es=[]
Fs=[]
#perView = 50
perView = 5
for i in range(nstep):
    for j in range(perView):
        F2=sch.stepResonance( E0=E0, dt=0.5 )        #;print("F2 ",np.sqrt(F2) )       #;print( "fpsi \n", sch.fpsi ); 
        Es.append(sch.EQF[0])
        Fs.append(np.sqrt(F2))
    ij = i*perView 
    plt.subplot(3,nstep,i        +1); plt.imshow(sch.psi ); plt.title("Psi_%i" %(ij+1)     ); plt.colorbar()
    plt.subplot(3,nstep,i+ nstep +1); plt.imshow(sch.Apsi); plt.title("APsi_%i" %(ij+1)    ); plt.colorbar()
    plt.subplot(3,nstep,i+2*nstep+1); plt.imshow(sch.fpsi); plt.title("dPsi/dt_%i" %(ij+1) ); plt.colorbar()

plt.figure(); plt.plot(Es,label="E");      plt.plot(Fs,label="|Ferr|"); plt.axhline(0,ls='--',c='k'); plt.legend(); plt.grid()
#plt.figure(); plt.plot(Fs,label="|Ferr|"); plt.axhline(0,ls='--',c='k'); plt.yscale('log'); plt.legend(); plt.grid()
plt.show()
exit()



'''
nstep=20
plt.figure(figsize=(3*nstep,6))
Es=[]
Fs=[]
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


'''
nstep=20
plt.figure(figsize=(3*nstep,6))
Es=[]
Fs=[]
perView = 100
nend=0
for i in range(nstep):
    for j in range(perView):
        F2=sch.step_Green( )   #;print(F2)  #;print( "fpsi \n", sch.fpsi ); 
        #print(i,j,F2)
        f = np.sqrt(F2)
        Fs.append( f )
        if f<1e-4: break
    ij = i*perView 
    plt.subplot(2,nstep,i+1 ); plt.imshow(sch.psi, extent=extent ); plt.title("Psi_%i"     %(ij+1) ); plt.colorbar()
    if f<1e-4: 
        break
#print(Fs)
plt.figure();plt.plot(Fs,label="|Ferr|"); plt.axhline(0,ls='--',c='k'); plt.yscale('log'); plt.legend(); plt.grid()
'''

xs = np.linspace( extent[0]*0.9, extent[1]*0.9, 40 )
ys = np.linspace( extent[2]*0.9, extent[3]*0.9, 40 )

#sources = []
#for x in xs:
#    for y in ys:
#       sources.appned( (x,y) ) 

Qs = np.zeros((len(xs),len(ys)))
i=0
for ix,x in enumerate(xs):
    for iy,y in enumerate(ys):
        #for i,p in enumerate(sources):
        sch.source[:,:] = gauss( 10.0, x, y )
        nitr = sch.solve_Green( maxErr=1e-2, maxIters=1000 ); 
        Q    = np.sum(sch.psi**2 ) 
        Qs[ix,iy] = Q
        print( i,x,y, Q, nitr )
        i+=1


plt.figure(figsize=(20,5))
plt.imshow( Qs )


vmin=sch.psi.min()
vmax=sch.psi.max()
vmax=max(-vmin,vmax)

plt.figure(figsize=(20,5))
plt.subplot(1,4,1); plt.imshow(sch.V,      extent=extent); plt.title("V "     );  plt.colorbar()
plt.subplot(1,4,2); plt.imshow(sch.source, extent=extent); plt.title("source" );  plt.colorbar()
plt.subplot(1,4,3); plt.imshow(sch.psi**2, extent=extent); plt.title("|Psi|^2"  );  plt.colorbar()
plt.subplot(1,4,4); plt.imshow(sch.psi,    extent=extent, cmap='seismic', vmin=-vmax,vmax=vmax);    plt.title("Psi"  );  plt.colorbar()

plt.savefig("solution.png",bbox_inches='tight' )

#plt.axis('equal')
plt.show(  )


