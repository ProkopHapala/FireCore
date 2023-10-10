import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall  import atomicUtils as au
from pyBall.atomicUtils import AtomicSystem
from pyBall  import plotUtils   as plu

#from pyBall import MMFF as mmff

import polymerUtils as pu

sys.path.append("../tLammpsTrj")
import BasePairUtils as bpu

# ======== Setup

folder    = "/home/prokop/Desktop/CARBSIS/Paolo/endgroups/"
dir_meta  = folder+"endgroups/"
dir_relax = folder+"endgroups_relaxed/mols/"

'''
sequence = [
("HNH-h","OHO-h_1"),
("OHO-h_1","HNH-h"),
("HNH-h","OHO-h_1"),
]
'''

'''
sequence = [
("HH-hh-p","NN-hh"),
("NN-hh","HH-hh-p"),
("HH-hh-p","NN-hh"),
]
'''


sequence = [
("HHH-h-p","NNN-hhh"),
("NNN-hhh","HHH-h-p"),
("HHH-h-p","NNN-hhh"),
]


out_name = "sequence_1"
amargin = 5.0
#======== Functions

#======== Body

#print(dir_relax)
#print(dir_meta)

names = os.listdir( dir_relax )
#print( names )

group_dict, typs = pu.load_groups( names )

for typ in typs: print( typ, typs[typ]  )

B = AtomicSystem(fname='backbone.xyz' )
B.lvec = np.array( [[25.,0.,0.],[0.,5.,0.],[0.,0.,20.0]  ] )

#monomers = [ ]

fig = plt.figure(figsize=(16,4))

BBB = None

inds1_ = []
inds2_ = []
inds3_ = [] 

for i,pair in enumerate(sequence):
    name1,name2 = pair
    BB, inds1, inds2 = pu.attachPair( B, name1, name2, group_dict, amargin=amargin )

    # find the indices of nitrogen atoms which are hydrogen bond donors
    cls1   = bpu.name_to_class( name1 )
    inds_  = [ inds1[i] if cls1[i]=='H' else inds2[i] for i in range(len(inds1)) ]    ;print(  "inds_ ",  inds_ )
    inds3 = BB.getNeighsOfType( inds_, typ='N')
    inds3 = [ a[0] for a in inds3 ]

    BB.apos[:,:] += B.lvec[1,:][None,:] * i 

    i0=0
    if BBB is None:
        BBB = BB.clonePBC()
    else:
        i0 = len(BBB.apos)
        BBB.append_atoms( BB )
        #BBB.lvec[1,:] += B.lvec[1,:] 
    print(  "inds1, inds2, inds3 ", len(inds1), len(inds2), len(inds3),  inds1, inds2, inds3 )
    inds1_ += [ i0+j for j in inds1 ]
    inds2_ += [ i0+j for j in inds2 ]
    inds3_ += [ i0+j for j in inds3 ]

    #monomers.append( BB )

    '''
    plt.subplot(1,2,1)
    axes=(0,1)
    plu.plotSystem( BB, bBonds=True, colors=None, sizes=None, extent=None, sz=50., RvdwCut=0.5, axes=axes, bLabels=True, labels=None, _0=1, HBs=None, bHBlabels=True, bBLabels=False )
    plu.plotAtoms( BB.apos, colors='#FF0000', marker='o', axes=axes, selection = inds1 )
    #plu.plotAtoms( BB.apos, colors='#FFFF00', marker='o', axes=axes, selection = inds1b )
    plu.plotAtoms( BB.apos, colors='#0000FF', marker='o', axes=axes, selection = inds2 )
    axes=(0,2)
    plt.subplot(1,2,2)
    plu.plotSystem( BB, bBonds=True, colors=None, sizes=None, extent=None, sz=50., RvdwCut=0.5, axes=axes, bLabels=True, labels=None, _0=1, HBs=None, bHBlabels=True, bBLabels=False )
    plu.plotAtoms( BB.apos, colors='#FF0000', marker='o', axes=axes, selection = inds1 )
    #plu.plotAtoms( BB.apos, colors='#FFFF00', marker='o', axes=axes, selection = inds1b )
    plu.plotAtoms( BB.apos, colors='#0000FF', marker='o', axes=axes, selection = inds2 )
    '''
    
print(  "inds1, inds2, inds3 ", inds1_, inds2_, inds3_ )
BBB.lvec[1,:] = B.lvec[1,:]*3
BBB.findBonds( )

BBB.saveXYZ( "./out/"+out_name+".xyz", comment=" Hbonds={'X':"+str(inds1_)+",'Y':"+str(inds2_)+"'N':"+str(inds3_)+"}" )
pu.saveMolGUIscript( out_name, (inds1_,inds2_,inds3_), path="./out/", amargin=amargin-3.0 )

plt.subplot(1,2,1)
axes=(0,1)
plu.plotSystem( BBB, bBonds=True,  axes=axes, bLabels=False )
plu.plotAtoms( BBB.apos, colors='#FF0000', marker='o', axes=axes, selection = inds1_ )
plu.plotAtoms( BBB.apos, colors='#0000FF', marker='o', axes=axes, selection = inds2_ )
plu.plotAtoms( BBB.apos, colors='#FFFF00', marker='o', axes=axes, selection = inds3_ )
axes=(0,2)
plt.subplot(1,2,2)
plu.plotSystem( BBB, bBonds=True,axes=axes, bLabels=False )
plu.plotAtoms( BBB.apos, colors='#FF0000', marker='o', axes=axes, selection = inds1_ )
plu.plotAtoms( BBB.apos, colors='#0000FF', marker='o', axes=axes, selection = inds2_ )
plu.plotAtoms( BBB.apos, colors='#FFFF00', marker='o', axes=axes, selection = inds3_ )

plt.tight_layout()
#plt.savefig( odir+name+".png", bbox_inches='tight' )
#plt.close(fig)

plt.show()