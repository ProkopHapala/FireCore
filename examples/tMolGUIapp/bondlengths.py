#!/usr/bin/python

import sys
import numpy      as np
import matplotlib.pyplot as plt

sys.path.append('../../')
from pyBall import atomicUtils as au
from pyBall import plotUtils   as plu

def addToHist( hist, u, w=1.0 ):
    #u = (x-x0)/dx
    i = int(u);
    f = u-i
    if (i>0) and (i<len(hist)-1):
        hist[i  ] += w*(1-f)
        hist[i+1] += w*(f)

def makeHists( typeParis, bonds, apos, types, histgrid=(0.5,2.5,0.1) ):
    x0,xmax,dx = histgrid
    #nhist = int( (xmax-x0)/dx )+1
    xs    = np.arange(x0,xmax+1e-9,dx) 
    nhist = len(xs)
    hists = { ts : np.zeros(nhist) for ts in typeParis }
    print( "types ",  types )
    for ib,b in enumerate(bonds):
        i,j = b
        d = apos[j] - apos[i]
        l = np.sqrt( np.dot(d,d)  )
        ti=types[i]; 
        tj=types[j]; 
        if(ti>tj): 
            ts=(ti,tj)
        else:
            ts=(tj,ti)
        #print( ib, ts, l )
        if ts in hists: addToHist( hists[ts], (l-x0)/dx )
    return hists, xs

#mol = au.AtomicSystem("Si4930_110-H-relaxed.xyz")
#mol = au.AtomicSystem("Si2647_100-H-relaxed.xyz")
mol = au.AtomicSystem("Si2505_111-H-relaxed.xyz")

mols = [
"Si2505_111-H-relaxed.xyz",
"Si2505_111-noH-SiH3-relaxed.xyz",
"Si2505_111-H-brak-110-relaxed.xyz",
]

ts = (14,14)

for imol,molname in enumerate(mols):
    mol = au.AtomicSystem(molname)
    mol.findBonds()
    #print( "mol.bonds ", mol.bonds )
    #print( "len(mol.bonds) ", len(mol.bonds) )
    #typeParis=[(14,14),(14,1)]
    typeParis=[ ts ]
    hists, xs = makeHists( typeParis, mol.bonds, mol.apos, mol.atypes, histgrid=(2.2,2.5,0.001) )
    #for ts in typeParis:
    #    plt.plot( xs, hists[ts], label=str(ts) )
    plt.plot( xs, hists[ts], label=str(molname) )
    #plt.figure( figsize=(15,15) )
    #plu.plotSystem( mol, bLabels=False )


plt.legend()
plt.savefig("histogram.png")# box_inches='tight' )
plt.show()
