import numpy             as np
import matplotlib.pyplot as plt

# ======= Setup

Ei = 0.1
Ej = 0.15
Ri = 1.4
Rj = 1.9

# ======= Functions

def numDeriv( xs, Es ):
    Fs  = ( Es[2:]-Es[:-2] )/( xs[2:]-xs[:-2] )
    xs_ = xs[1:-1]
    return Fs,xs_

def getMorse( x, Eij=0.01, Rij=1.6+1.6, b=1.0 ):
    expr = np.exp( -b*( x - Rij ) )
    V    =      Eij*( expr*expr - 2*expr )
    F    = -2*b*Eij*( expr*expr -   expr ) 
    return V,F

'''
def getMorse_Half( x, Ei=Ei, Ri=Ri, b=1.0 ):
    expr = np.exp( -b*( x - Ri ) )
    EP    =      Ei*expr*expr
    FP    = -2*b*Ei*expr*expr 
    EL    =  -2*Ei*expr
    FL    = 2*b*Ei*expr
    return EP,FP,EL,FL
'''

'''
def getMorse_Half( x, Ei=Ei, Ri=Ri, b=1.0 ):
    expr = np.exp( -b*( x - Ri ) )
    eM    = Ei*expr
    EP    =      eM*expr
    FP    = -2*b*eM*expr 
    EL    = -2*  eM
    FL    =  2*b*eM
    return EP,FP,EL,FL
'''

def getMorse_Half( x, Ei=Ei, Ri=Ri, b=1.0 ):
    expr = np.exp( -b*( x - Ri ) )
    eM    = Ei*expr
    de    = -2*b*eM
    EP    =    eM*expr
    FP    =    de*expr 
    EL    = -2*eM
    FL    =   -de
    return EP,FP,EL,FL


def getMorse_Coefs( x, Ej=Ej, Rj=Rj, b=1.0 ):
    expr = np.exp( b*Rj )
    cP    = Ej*expr*expr
    cL    = Ej*expr
    return cP,cL

# ========== Main

x   = np.linspace(2.0,5.0,300)
E,F        = getMorse( x, Eij=Ei*Ej, Rij=Ri+Rj, b=1.0 )
F_num,xnum = numDeriv( x, E )


EP,FP,EL,FL = getMorse_Half( x,  Ei=Ei, Ri=Ri, b=1.0 )
cP,cL       = getMorse_Coefs( x, Ej=Ej, Rj=Rj, b=1.0 )

E_ = EP*cP + EL*cL
F_ = FP*cP + FL*cL

plt.figure(figsize=(5,10))
plt.subplot(2,1,1); plt.plot( x,E,label='E_Morse'); plt.plot( x,E_,':', label='E_Mod' );                                              plt.grid(); plt.legend()
plt.subplot(2,1,2); plt.plot( x,F,label='F_Morse'); plt.plot( x,F_,':', label='F_Mod' );  plt.plot( xnum,F_num, ':', label='F_num' ); plt.grid(); plt.legend()



plt.show()