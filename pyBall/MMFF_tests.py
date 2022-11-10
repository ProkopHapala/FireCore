import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF as mmff

def deriv( xs,Es ):
    xs_ =(xs[2:]+xs[:-2])*0.5
    dxs =(xs[2:]-xs[:-2])
    dEs =(Es[2:]-Es[:-2])
    return dEs/dxs,xs_

def GridFF( surf_name="data/NaCl_sym-center", atyp=2, Q =-0.01, K=-1.6, kind=1, Rdamp=2.0, pos0=(0.0,0.0,0.0) ):
    mmff.tryInit()
    # ------------------ 1D
    #rs = np.linspace(-5.0,5.0,150)
    #rs = np.linspace(0.0,20.0,100)
    rs = np.linspace(2.2,10.0,78)
    #EsM,fsM  = mmff.sampleNonBond( rs, kind=kind, REQi=(R,e  ,Q*0), K=K, Rdamp=Rdamp)  
    #EsC,fsC  = mmff.sampleNonBond( rs, kind=kind, REQi=(R,e*0,Q  ), K=K, Rdamp=Rdamp)   
    Eg,fg   = mmff.sampleSurf( surf_name,   rs, kind=12, atyp=atyp, Q=Q, K=K, Rdamp=Rdamp, pos0=pos0, bSave=True)  
    Ea,fa   = mmff.sampleSurf( None,        rs, kind=1 , atyp=atyp, Q=Q, K=K, Rdamp=Rdamp, pos0=pos0, bSave=True)  
    '''
    Eg,fg   = mmff.sampleSurf( "data/H_atom",  rs, kind=11, atyp=atyp, Q=Q, K=K, Rdamp=Rdamp, pos0=pos0, bSave=True) 
    #Eg,fg   = mmff.sampleSurf( "data/H_atom",  rs, kind=5, atyp=atyp, Q=Q, K=K, Rdamp=Rdamp, pos0=pos0, bSave=True) 
    #Ea,fa   = mmff.sampleSurf( None,           rs, kind=1, atyp=atyp, Q=Q, K=K, Rdamp=Rdamp, pos0=pos0, bSave=True)  
    Ea,fa   = mmff.sampleSurf( None,  rs, kind=5, atyp=atyp, Q=Q, K=K, Rdamp=Rdamp, pos0=pos0, bSave=True) 
    #Eg*=1.25
    '''
    #print("Es \n", Es);
    fg_,xfs = deriv( rs,Eg )
    fa_,xfs = deriv( rs,Ea )
    #Emin=Es.min()
    #Emax=Es.max()
    #Emin=Es.min()
    plt.figure(figsize=(5,10))
    plt.subplot(2,1,1); 
    plt.plot(rs,Eg,label='Grid');  plt.plot(rs,Ea,"--",label='Atoms'); plt.grid(); plt.legend(); plt.axhline(0,c='k',ls='--');                 #plt.ylim(Emin*sc,Emax*sc)  
    plt.subplot(2,1,2); 
    plt.plot(rs,-fg,'-r', label='Fg'); #plt.plot(xfs,fg_, 'r--',label='Fg_num')    
    plt.plot(rs,-fa,'-b', label='Fa'); plt.plot(xfs,fa_, 'b:',label='Fa_num')   
    plt.grid(); plt.axhline(0,c='k',ls='--'); plt.legend(); #plt.ylim(Emin*sc,Emax*sc)  
    plt.show()
    #exit()
