import sys
import matplotlib.pyplot as plt
sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import plotUtils   as plu

#selection=[1,2,3,4]

selection=range(100)
trj = au.readLammpsTrj( fname="traj.lammpstrj", nmax=100, selection=set(selection) )
print( len(trj) )

for i,sys in enumerate(trj):
    print("plot # ", i)
    fig = plt.figure(figsize=(5,5))
    #plu.plotAtoms( apos=sys.apos, es=sys.enames, sizes=100., colors='#808080', marker='o', axes=(0,1) )
    #sys.findBonds( Rcut=3.0, RvdwCut=0.8 )
    sys.findBonds( Rcut=3.0, RvdwCut=0.9 )
    colors = [ au.elements.ELEMENT_DICT[e][8]     for e in sys.enames ]
    sizes  = [ au.elements.ELEMENT_DICT[e][6]*50. for e in sys.enames ]
    plu.plotBonds( links=sys.bonds, ps=sys.apos, axes=(0,1) )
    plu.plotAtoms( apos=sys.apos, es=sys.enames, sizes=sizes, colors=colors, marker='o', axes=(0,1) )

    plt.xlim(-10.0,10.0)
    plt.ylim( 0.0,20.0)
    plt.savefig( "mol_%03i.png" %selection[i], bbox_inches='tight' )
    plt.close(fig)

#plt.show()