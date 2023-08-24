import sys
import os
import matplotlib.pyplot as plt
from   matplotlib import collections  as mc
sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import plotUtils   as plu


# ============ Setup

sz      =5
dirname ='.'
if len(sys.argv)>1: dirname=sys.argv[1]

# ============ Functions

def plotBondLenghts( mol, axes=(0,1), HBs=None, Bs=None ):
    plu.plotSystem( mol,                    axes=axes, bBonds=False, bLabels=False )
  
    # Bonds
    if Bs is None:
        bs,rbs = mol.findBonds ( )                 #;print( bs, rbs )
    else:
        bs,rbs = Bs
    rb_labs = [ ("%3.2f" %r) for r in rbs ] 
    plu.plotBonds ( ps=mol.apos, links=bs,  axes=axes, colors="#808080", labels=rb_labs )
    
    # H-Bonds
    if HBs is None:
        hbs,rhbs = mol.findHBonds( bPrint=True )     #;print( hbs, rhbs )
    else:
        hbs,rhbs = HBs
    if len(hbs)>0:
        rh_labs = [ ("%3.2f" %r) for r in rhbs ] 
        plu.plotBonds ( ps=mol.apos, links=hbs, axes=axes, colors="g",       labels=rh_labs )
    return (bs,rbs),(hbs,rhbs)
        
# ============ Main

ax_names =['x','y','z']
ax_names2=['PCA.a','PCA.b','PBC.x']
axes=(0,1)

names = os.listdir( dirname )    #;print(names)

for name in names:
    name_,ext = os.path.splitext(name)
    #print(name, ext)
    if ext != '.xyz': continue
    
    mol = au.AtomicSystem( name ) 

    plt.figure(figsize=(2*sz,sz))

    plt.subplot(1,2,1); Bs,HBs = plotBondLenghts( mol, axes=axes );  plt.xlabel(ax_names[axes[0]]);plt.ylabel(ax_names[axes[1]])
    mol.orientPCA()
    #mol.saveXYZ("post_PCA."+name  ,   other_lines=other_lines)
    plt.subplot(1,2,2); plotBondLenghts( mol, axes=axes, Bs=Bs, HBs=HBs );  plt.xlabel(ax_names2[axes[0]]);plt.ylabel(ax_names2[axes[1]])
    
    plt.savefig(name+".png", bbox_inches='tight' )
    plt.close(plt.gcf())

#plt.show()
