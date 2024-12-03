import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall  import atomicUtils as au
from pyBall.atomicUtils import AtomicSystem
#from pyBall  import plotUtils   as plu
#from pyBall import MMFF as mmff

import polymerUtils as pu

# ======== Setup

folder="/home/prokop/Desktop/CARBSIS/Paolo/endgroups/"
dir_meta  = folder+"endgroups/"
dir_relax = folder+"endgroups_relaxed/mols/"
#print(dir_relax)
#print(dir_meta)

pairs = [
("HNH-h","OHO-h_1")
]

pairTypes = [
("HeH","eHe"),
#("HHH","eee"),
#("HHe","Hee"),
#("HHH","eee"),
#("HH","ee"),
]

#======== Functions

#======== Body

names = os.listdir( dir_relax )
#print( names )

group_dict, typs = pu.load_groups( names, folder )

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

            #print( inds )
            comment = " Hbonds={'X':"+str(inds1)+",'Y':"+str(inds2)+"}"
            #BB.enames[ inds1 ] = 'As'
            #BB.enames[ inds2 ] = 'P'

            #BB.enames[ inds1[0] ] = 'F';BB.enames[ inds1[1] ] = 'Cl';BB.enames[ inds1[2] ] = 'Br';
            #BB.enames[ inds2[0] ] = 'Ne';BB.enames[ inds2[1] ] = 'Ar'; BB.enames[ inds2[2] ] = 'Kr';

            name = "BB."+name1+"."+name2

            print( name+".Qtot: ", BB.qs.sum() )

            BB.saveXYZ( odir+name+".xyz", comment=comment )

            BB_ = BB.clonePBC( (2,2,1) )
            BB_.saveXYZ( odir+"/2x2/"+name+"_2x2.xyz", comment=comment )

            pu.saveMolGUIscript( name, (inds1,inds2), path="./out/", amargin=amargin-3.0 )







'''
mmff.setVerbosity( verbosity=1, idebug=0 )
mmff.init( xyz_name="data/pyridine", surf_name="data/NaCl_1x1_L2", bMMFF=False  )              # without MMFF
mmff.getBuffs()
mmff.eval()
print( "FORCES:\n mmff.fapos:\n ", mmff.fapos )
mmff.plot(bForce=True, Fscale=10.0 )
plt.show()
exit(0)
'''




