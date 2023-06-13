import sys,os
import numpy as np
import matplotlib.pyplot as plt
#sys.path.append("../../")
sys.path.append("/home/prokop/git/FireCore/")
from pyBall import atomicUtils as au
from pyBall import plotUtils   as plu

# ============ Functions

def plotBondLenghts( mol, axes=(0,1) ):
    hbs,rhbs = mol.findHBonds( bPrint=True )     #;print( hbs, rhbs )
    bs,  rbs = mol.findBonds ( )                 #;print( bs, rbs )
    rh_labs = [ ("%3.2f" %r) for r in rhbs ] 
    rb_labs = [ ("%3.2f" %r) for r in rbs ] 
    plu.plotSystem( mol,                    axes=axes, bBonds=False, bLabels=False )
    plu.plotBonds ( ps=mol.apos, links=hbs, axes=axes, colors="g",       labels=rh_labs )
    plu.plotBonds ( ps=mol.apos, links=bs,  axes=axes, colors="#808080", labels=rb_labs )
    plt.savefig(fname+".png", bbox_inches='tight' )

# ============= Main

fname = sys.argv[1]
mol   = au.AtomicSystem( fname  )
mol.orientPCA()
mol.saveXYZ("CG_PCA.xyz")


#plt.figure(figsize=(5,5)) 
plt.figure(figsize=(8,8)) 
plotBondLenghts( mol, axes=(1,0) )
#plotBondLenghts( mol, axes=(0,1) )
plt.show()
