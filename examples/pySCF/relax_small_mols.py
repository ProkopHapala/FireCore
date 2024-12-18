# https://pyscf.org/quickstart.html
# https://pyscf.org/pyscf_api_docs/pyscf.gto.html#module-pyscf.gto.mole

'''
Molecules:
H2O
NH3
HCOOH
CH2=O
'''

import sys
import os
#from tkinter import UNITS
import pyscf
sys.path.append('../../')
from pyBall import atomicUtils as au
from pyBall import FFFit       as fff
from pyBall import pyscf_utils as scfu
from pyBall import plotUtils   as plu
import numpy as np
import matplotlib.pyplot as plt


geoms=[
    # --- Linear
    ("HCCH",   'H -1 0 0; C 0 0 0; C 1 0 0; H 2 0 0'),  # HCCN
    ("HCN",    'H -1 0 0; C 0 0 0; N 1 0 0'),           # HCN
    ("CO",     'O -1 0 0; C 0 0 0; '),                  # CO
    ("HF",     'H -1 0 0; F 0 0 0'),            # HF
    ("HCl",    'H -1 0 0; Cl 1.0 0 0'),         # HCl
    ("NaCl",   'Na 0 0 0; Cl 2.5 0 0'),         # NaCl
    # --- pyramidal
    ("H2O",    'H -1 0 0; O 0 0 0; H 0 1 0'),             # H2O
    ("NH3",    'H -1 0 0; N 0 0 0; H 1 +1 1; H 1 -1 1'),  # NH3
    # --- carbonyl
    ("OCH2",  'O -1 0 0; C 0 0 0; H 1 +1 0; H 1 -1 0; '),       # HCOH
    ("HCOOH", 'H -1 0 0; O 0 0 0; C 1 1 0; O 0 2 0; H 2 1 0; '),  # HCOOH
]

#geoms=[
#    ("H2O",    'H -1 0 0; O 0 0 0; H 0 1 0'),             # H2O
#    ("NH3",    'H -1 0 0; N 0 0 0; H 1 +1 1; H 1 -1 1'),  # NH3
#] 

def plot( es, apos):
    plt.figure(figsize=(10,5))
    plt.subplot(1,2,1); plu.plotAtoms( es=es, apos=apos, ax1=0, ax2=1 )  ;plt.xlabel("x[A]");plt.ylabel("y[A]");
    plt.subplot(1,2,2); plu.plotAtoms( es=es, apos=apos, ax1=2, ax2=1 )  ;plt.xlabel("z[A]");plt.ylabel("y[A]");
    plt.title(geom[0])

scfu.verbosity = 4

for i,geom in enumerate(geoms):
    print( i, " ", geom[0], geom[1] )
    mol = pyscf.M(atom=geom[1])
    mol.build()

    #opt=pyscf.scf.UHF(mol).Gradients().optimizer(solver='berny')
    opt = pyscf.scf.RHF(mol).Gradients().optimizer(solver='berny')
    mol = opt.kernel()

    apos, es = scfu.unpack_mol( mol ) #;print(apos)
    
    cog=np.sum(apos,axis=0)/len(apos)
    
    # H-pointing orientialtion
    p1 =apos[0,:]
    fw1=apos[1,:]-apos[0,:]
    up1=cog-apos[1,:]
    if(len(es)>2): up1[:]+=apos[2,:]-apos[1,:]
    apos = au.orient_vs( p1, fw1, up1, apos.copy()-p1[None,:], trans=None ); 
    apos-=apos[0,:][None,:]
    au.saveXYZ( es, apos, geom[0]+"_acid.xyz" )
    
    # by center
    #p2 =apos[1,:]
    #fw2=apos[1,:]-cog        ;print("fw2 ", fw2)
    #up2=apos[1,:]-apos[0,:]  ;print("up2 ", up2)
    #apos = au.orient_vs( p2, fw2, up2, apos.copy()-p1[None,:], trans=None );
    #au.saveXYZ( es, apos, geom[0]+"_base.xyz" )

    plot( es, apos) ; plt.show()

    #h2o.verbose = verbosity
    #calc = pyscf.scf.RHF(h2o)
    plt.show()

#plt.plot( xs, Es );   plt.grid(); 
plt.show()