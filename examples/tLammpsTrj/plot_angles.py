import sys,os
import numpy as np
import matplotlib.pyplot as plt
#sys.path.append("../../")
sys.path.append("/home/prokop/git/FireCore/")
from pyBall import atomicUtils as au
from pyBall import plotUtils   as plu


# ========= Functions

def plotAngles( mol, typ='N', bDih=False, bB=False, bHb=True, neighTyps={'H':(1,2)}, axes=(0,1) ):
    bs,  rbs = mol.findBonds ( )                 #;print( bs, rbs )
    ngs = mol.neighs()
    plu.plotSystem( mol, axes=axes, bBonds=True, bLabels=False )
    if bB:
        rb_labs = [ ("%3.2f" %r) for r in rbs ] 
        plu.plotBonds ( ps=mol.apos, links=bs,  axes=axes, colors="#808080", labels=rb_labs )
    if bHb:
        hbs,rhbs = mol.findHBonds( bPrint=True )     #;print( hbs, rhbs )
        rh_labs = [ ("%3.2f" %r) for r in rhbs ] 
        plu.plotBonds ( ps=mol.apos, links=hbs, axes=axes, colors="g",       labels=rh_labs )
    if bDih:
        sel = mol.select_by_neighType( ngs, typ=typ, neighTyps=neighTyps )     #;print("selection: ", sel)
        angs, iangs = mol.findDihedral( ngs=ngs, select=sel, neighTyp={'H'} )
    else:   
        sel = mol.select_by_ename( typ )   #;print("selection: ", sel)
        angs, iangs = mol.findAngles( select=sel, ngs=ngs )    # passing ngs is just performance-optimization
    plu.plotAngles( iangs, angs, mol.apos,  axes=axes, colors='b',  bPoly=True )



# ========= Main

fnames=[
"HHH-hhS1_NNO-hh_1S-3-final.xyz",
"HH-hh-pS1_HNO-pS3-final.xyz",
]

for fname in fnames:
    plt.figure(figsize=(5,5))
    mol   = au.AtomicSystem( fname  )
    #plotAngles( mol, typ='N', bDih=False, bB=False, bHb=True, axes=(0,1) )
    plotAngles( mol, typ='N', bDih=True, bB=False, bHb=True, axes=(0,1) )
    plt.savefig(fname+".png", bbox_inches='tight' )

plt.show()