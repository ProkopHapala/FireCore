import sys
import os
import matplotlib.pyplot as plt
sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import plotUtils   as plu

# ============ Functions

mol = au.AtomicSystem( "/home/prokop/Desktop/CARBSIS/Mithun/resp-Hbond_close_small/cc-opt-CHONH2-e1_vs_CHONH2-H0.xyz" ) 

hbs,rbs = mol.findHBonds( bPrint=True )
print( hbs, rbs )

r_labesl = [ ("%3.3fA" %r) for r in rbs ] 
plu.plotSystem( mol, axes=(1,2) )
plu.plotBonds( ps=mol.apos, links=hbs, axes=(1,2), colors="g", labels=r_labesl )

plt.show()
