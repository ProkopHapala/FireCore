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

# ======== Setup

folder    = "/home/prokop/Desktop/CARBSIS/Paolo/endgroups/"
dir_meta  = folder+"endgroups/"
dir_relax = folder+"endgroups_relaxed/mols/"

pairs = [
("HNH-h","OHO-h_1")
]

pairTypes = [
("HeH","eHe"),
#("HHe","Hee"),
#("HHH","eee"),
#("HH","ee"),
#("HH","eee"),
#("HHH","ee"),
]

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
#for pair in pairs:
#    name1, name2 = pair
#    attachPair( name1, name2, group_dict )

amargin = 5.0
for pairTyp in pairTypes:
    names1 = typs[pairTyp[0]]
    names2 = typs[pairTyp[1]]

    bPure = all( 'H' == c for c in pairTyp[0] )

    odir = "out_"+pairTyp[0]+"_"+pairTyp[1]+"/"
    try:
        os.mkdir( odir )
        os.mkdir( odir+"/2x2/" )
    except:
        pass
    for name1 in names1:
        for name2 in names2:
            #print( name1, name2 )
            BB, inds1, inds2 = pu.attachPair( B, name1, name2, group_dict, amargin=amargin )
            
            # if pairTyp[0]='H' take from inds1 else from inds2
            inds_ = [ inds1[i] if pairTyp[0][i]=='H' else inds2[i] for i in range(len(inds1)) ]    ;print(  "inds_ ",  inds_ )
            inds1b = BB.getNeighsOfType( inds_, typ='N')
            
            print( "inds1, inds2, inds1b ", inds1, inds2, inds_, inds1b )
            inds1b = [ a[0] for a in inds1b ] 

            comment = " Hbonds={'X':"+str(inds1)+",'Y':"+str(inds2)+"}"

            name = "BB."+name1+"."+name2

            print( name+".Qtot: ", BB.qs.sum() )

            BB.saveXYZ( odir+name+".xyz", comment=comment )

            BB_ = BB.clonePBC( (2,2,1) )
            BB_.saveXYZ( odir+"/2x2/"+name+"_2x2.xyz", comment=comment )

            if len(inds1)==len(inds2):
                pu.saveMolGUIscript( name, (inds1,inds2,inds1b), path="./out/", amargin=amargin-3.0 )

            
            fig = plt.figure(figsize=(16,4))
            axes=(0,1)
            plt.subplot(1,2,1)
            plu.plotSystem( BB, bBonds=True, colors=None, sizes=None, extent=None, sz=50., RvdwCut=0.5, axes=axes, bLabels=True, labels=None, _0=1, HBs=None, bHBlabels=True, bBLabels=False )
            plu.plotAtoms( BB.apos, colors='#FF0000', marker='o', axes=axes, selection = inds1 )
            plu.plotAtoms( BB.apos, colors='#FFFF00', marker='o', axes=axes, selection = inds1b )
            plu.plotAtoms( BB.apos, colors='#0000FF', marker='o', axes=axes, selection = inds2 )
            axes=(0,2)
            plt.subplot(1,2,2)
            plu.plotSystem( BB, bBonds=True, colors=None, sizes=None, extent=None, sz=50., RvdwCut=0.5, axes=axes, bLabels=True, labels=None, _0=1, HBs=None, bHBlabels=True, bBLabels=False )
            plu.plotAtoms( BB.apos, colors='#FF0000', marker='o', axes=axes, selection = inds1 )
            plu.plotAtoms( BB.apos, colors='#FFFF00', marker='o', axes=axes, selection = inds1b )
            plu.plotAtoms( BB.apos, colors='#0000FF', marker='o', axes=axes, selection = inds2 )

            plt.tight_layout()
            plt.savefig( odir+name+".png", bbox_inches='tight' )
            plt.savefig( odir+name+".svg", bbox_inches='tight' )
            plt.close(fig)