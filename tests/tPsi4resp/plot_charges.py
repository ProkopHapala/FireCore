import sys
import os
import time
import psi4
import resp
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall  import atomicUtils as au
from pyBall  import plotUtils   as pltu

# ========= Setup

#indir="./input/"
indir="./pbe/cc-pvdz/"
# ======== Main

names = [ f.split('.')[0] for f in os.listdir(indir) ]
#names =["hexa_hb3_donor"]
print(names)

for name in names:
    print( "plotting %s " %name )
    fig = plt.figure(figsize=(5,5))
    apos,Zs,es,qs = au.loadAtomsNP( fname=indir+name+".xyz" ) 
    bonds = au.findBondsNP( apos, atypes=Zs )
    #labels = [ "%7.2f" %q for q in qs ]
    labels = [ "%s%7.2f" %(es[i],qs[i]) for i in range(len(es)) ]
    pltu.plotBonds( links=bonds, ps=apos, lws=None )
    pltu.plotAtoms( apos=apos, es=es, ax1=0, ax2=1, labels=labels, marker='o', colors=qs, sizes=150 )
    plt.savefig( indir+name+".png", bbox_inches='tight' )
    #plt.close(fig)

#plt.show()