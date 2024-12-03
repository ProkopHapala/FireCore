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
from pyBall import SchroedingerGreen1D as sch

extent=[ -100,100, -100,100 ]
L=(extent[1]-extent[0])

#sch.init(10,10, L=L )
sch.init(50,50, L=L)
#sch.init(100,100, L=L)
#sch.init(1000,1000, L=L)
sch.getBuffs()

#setStep(dstep=0.1,m_Me=1.0, L=None)

nx,ny = sch.V.shape

def gauss( w, x0=0. ):
    return np.exp( -( ( xs-x0)**2 )/(w*w) )

def prepareFucntions( extent=[ -5,5 ], width =10.0, strenght=0.1, x0=0., ):
    global xs
    xs    = np.linspace( extent[0], extent[1], nx ) #- nx*0.5
    #----- source
    sch.source[:,:] = gauss( width, x0 )
    # ---- Potential
    wV = 40.0
    sch.V[:,:]  = (1. - gauss( 40.0 )   ) * 0.05
    #sch.V[:,:] += (1. - np.sin(xs*(2*np.pi/50.0))**2 ) * 0.01
    print( "V min,max: ", sch.V.min(), sch.V.max() )
    return sch.V, sch.source

prepareFucntions( extent=extent, x0=10. )      #;print( "V\n", V )

#sch.V[:,:] = 0

sch.psi[:,:] = 1
#sch.psi[:,:] = 0
#sch.psi[:,:] = 0.01; sch.psi[nx//2,ny//2] = 1.
#sch.source[:,:] = 0;   sch.source[nx//2,ny//2] = -1.


#E0=0.0;
E0=0.05
#sch.EQF[3] = 2.0   #   set Energy
sch.EQF[3] = E0   #   set Energy

nstep=20
plt.figure(figsize=(5*nstep,3*5))
Es=[]
Fs=[]
#perView = 50
perView = 1
for i in range(nstep):
    for j in range(perView):
        F2=sch.stepResonance( E0=E0, dt=10.1 )        #;print("F2 ",np.sqrt(F2) )       #;print( "fpsi \n", sch.fpsi ); 
        Es.append(sch.EQF[0])
        Fs.append(np.sqrt(F2))
    ij = i*perView 

    plt.subplot(1,3,1); plt.plot(sch.psi ,label=str(i));
    plt.subplot(1,3,2); plt.plot(sch.Apsi,label=str(i));
    plt.subplot(1,3,3); plt.plot(sch.fpsi,label=str(i));


plt.subplot(1,3,1); plt.title("Psi_%i" %(ij+1)     );
plt.subplot(1,3,2); plt.title("APsi_%i" %(ij+1)    );
plt.subplot(1,3,3); plt.title("dPsi/dt_%i" %(ij+1) ); 

plt.figure(); plt.plot(Es,label="E");      plt.plot(Fs,label="|Ferr|"); plt.axhline(0,ls='--',c='k'); plt.legend(); plt.grid()
#plt.figure(); plt.plot(Fs,label="|Ferr|"); plt.axhline(0,ls='--',c='k'); plt.yscale('log'); plt.legend(); plt.grid()
plt.show()
exit()
