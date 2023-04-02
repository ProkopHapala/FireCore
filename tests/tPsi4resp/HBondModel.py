

import sys
import os

import numpy as np
import matplotlib.pyplot as plt
#from functools import partial
sys.path.append('../../')
from pyBall import Forces as FF

Hartree2eV   = 27.2114079527
Hartree2kcal = 627.503 

# ================ Setup

r0 = 1.91

REQs={
    'O': ( 1.661, 0.0006808, -0.3 ),
    'H': ( 1.487, 0.0091063, +0.3 )
}

Q     = 0.5
a     = 0.5
b     = 0.5

E0Hb = 0.015

R0 = REQs['O'][0]+REQs['H'][0] 
E0 = np.sqrt( REQs['O'][1]*REQs['H'][1] )
QQ = REQs['O'][2]*REQs['H'][2]

# ================ Main

data = np.genfromtxt('H2O_scan.E.dat')

xs=data[:,0] + r0
Es=data[:,1]

ELJ,FLJ  = FF.LenardJones(xs, E0,   R0 )
EHb,FHb  = FF.Vdw        (xs, E0Hb, R0 )
Eel1,Fel = FF.Coulomb    (xs-a  , -Q*Q )
Eel2,Fel = FF.Coulomb    (xs    , +Q*Q )
Eel3,Fel = FF.Coulomb    (xs+b  , -Q*Q )
Eel4,Fel = FF.Coulomb    (xs-a+b, +Q*Q )
Eel = Eel1 + Eel2 + Eel3 + Eel4



plt.plot( xs, Es*Hartree2eV,'.-k', label="Ref"   ) 
plt.plot( xs, ELJ+Eel+EHb  ,'--c', label="LJ+El+Hb" ) 
plt.plot( xs, ELJ+Eel      ,'--m', label="LJ+El" ) 
plt.plot( xs, Eel          ,':r' , label="El   "    )
plt.plot( xs, ELJ          ,':b' , label="LJ   " ) 
plt.plot( xs, EHb          ,':g' , label="Hb   " ) 
plt.ylim(-0.5,0.5)
plt.legend()
plt.ylabel("x[A]");  plt.ylabel("E[eV]  ");  plt.grid(); 


#plt.plot( xs, Es*Hartree2kcal,'.-' ); plt.ylabel("x[A]");  plt.ylabel("E[kcal]");  plt.grid(); 
plt.show()
