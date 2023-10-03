import sys
import os
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import plotUtils   as plu

import BasePairUtils as bpu

# ============ Setup 

#Econtrast_min  = 10.0
#Econtrast_min  = 5.0

dEcontrast_min = 10.0





#f_geoms="/home/prokop/Desktop/CARBSIS/Paolo/pairs_compare_DFTB_vs_B3LYP/2-plot_confs/dirs/"
#f_egs="/home/prokop/Desktop/CARBSIS/Paolo/endgroups/endgroups/"
f_egs="/home/prokop/Desktop/CARBSIS/Paolo/endgroups/B3LYP_oriented/"
dir_geom="/home/prokop/Desktop/CARBSIS/Paolo/B3LYP_finished/pairs_b3lyp/confs/"

'''
# ----- DD
sel_names = [
("HH-hh-p", 'r'),
("HH-h-p",   'b'),
#("HH-h",   'b'),
("HH-h_1",   'g'),
]
conuters = [
    "OO-h", 
    "NO-h",
    "NO-p",
    "NN-hh",
    "NN-hp",
    "NN-pp",
    "NNO-hh_1","NNO-hh_2","NNN-hhh",
]
'''


# ----- DDD
sel_names = [
("HHH-h-p", 'r'),
("HHH-h"  , 'b'),
]
conuters = [
    "OO-h",#"NO-h",
    "NO-p","NN-hh","NN-hp","NN-pp",
    "NNO-hh_1","NNO-hh_2","NNN-hhh",
]


'''
# ----- Guanine
sel_names = [
("HHO-h-p_1", 'r'),
("HHO-h-p_2", 'b'),
("HHO-h"  , 'g'),
]
conuters = [
    "HNO-h", "HNO-p", "HNO-h-p", "HNO-hh", "HNN-hh", "HNN-hp", "HNO-hp",
]
'''

'''
# ----- Cytosine
sel_names = [
("HNO-h-p", 'r'),
("HNO-h",   'b'),
]
conuters = [
    "HHO-h-p_1", "HHO-h-p_2", "HHO-h", "HHO-p", "HHN-hh", "HHO-hh", "HHO-hp", "HNN-hp",
]
'''

'''
# ----- DAD
sel_names = [
("HNH-h-p", 'r'),
("HNH-h",   'b'),
]
conuters = [
    "OHO-h_1", "OHO-h_2", "OHO-h-p", "NHO-hh", "NHO-hp", "OHO-p",
]
'''

'''
# ----- ADA
sel_names = [
("OHO-h-p", 'r'),
("OHO-h_1", 'b'),
]
conuters = [
    "HNH-hh", "HNH-hp", "HNH-h-p", "HNH-h", "HNH-p"
]
'''

sel_names_dct = dict( sel_names )
conuters_set = set( conuters )


lss=['-', '--', ':']


# ============ Main

#_,Es,names,pair_names = bpu.read_dat( "/home/prokop/Desktop/CARBSIS/Paolo/pairs_compare_DFTB_vs_B3LYP/4-dimers_energy_sorted.dat", ni=1, nf=3, iname=0, toRemove=bases_to_remove )
_,Es,names,pair_names           = bpu.read_dat( "/home/prokop/Desktop/CARBSIS/Paolo/B3LYP_finished/1-binding_energy_all.dat",   ni=4, nf=2, iname=0 )
_, Angs,   names_angs,_         = bpu.read_dat( "/home/prokop/Desktop/CARBSIS/Paolo/B3LYP_finished/dihedral_angle_ave.dat",   ni=0, nf=1, iname=0 )
_, Hbonds, names_hbonds,_       = bpu.read_dat( "/home/prokop/Desktop/CARBSIS/Paolo/B3LYP_finished/bond_length_min.dat",      ni=0, nf=1, iname=0 )
_, Hbonds_av, names_hbonds_av,_ = bpu.read_dat( "/home/prokop/Desktop/CARBSIS/Paolo/B3LYP_finished/bond_length_ave.dat",      ni=0, nf=1, iname=0 )

old_inds = np.array(  range(len(Es)) )
Angs     = np.array( Angs ); print( "Angs ", Angs ); #exit()

#my_set = set( ['NO-h', "HH-h-p" ])
# my_set = set(["HHH-h-p", "NNO-hh_1"])
# for p in pair_names: 
#     if (p[0] in my_set) and (p[1] in my_set):
#         print( "found ", p )

Emins, uique_pnames =   bpu.find_minim_energy_confs( Es, pair_names, Emax=1000, ipivot=0 )                          # find minimum energy conformers for each pair
#print("uique_pnames ", uique_pnames)

imins = np.array([ Emins[k][1]                                                    for k in  uique_pnames ])     # index of minimum energy conformer for each pair as numpy array
EBs   = np.array([-Emins[k][0]                                                    for k in  uique_pnames ])     # compute binding energy  for each pair as numpy array
ECs   = np.array([-Emins[k][0]+0.5*( Emins[(k[0],k[0])][0]+Emins[(k[1],k[1])][0]) for k in  uique_pnames ])     # compute energy contrast for each pair as numpy array
Amins = Angs[imins] 

print( "imins ", imins )
print( "Amins ", Amins )

selections=[]
for n,c in sel_names:
    sel = [ (i,n2 if(n1==n) else n1, EBs[i]) for i,(n1,n2) in enumerate(uique_pnames) if (n1==n)or(n2==n) ] 
    selections.append( sel )

    print( "##### ", n )
    for t in sel: print( t, Amins[t[0]] )

sel0 = [ (i,n,e) for i,n,e in selections[0] if n in conuters_set ]
sel0.sort( key=lambda k: k[2] )

print( "sel0 ",  sel0 )


xs = np.arange( len(sel0) )

plt.figure( figsize=( len(sel0)*0.20+0.5,5) )

for isel,sel in enumerate(selections):
    dct = { n:(i,n,e) for i,n,e in sel }
    sel = [ dct[n] for i,n,e in sel0 ] 
    selections[isel] = sel
    inds,names,Es = zip(*sel)
    inds=np.array(inds)
    plt.plot( xs,EBs[inds], '.b', ls=lss[isel], label=sel_names[isel][0] )
    plt.plot( xs,ECs[inds], '.r', ls=lss[isel] )
    plt.plot( xs,Amins[inds], '.g', ls=lss[isel] )

font_prop = FontProperties(family='monospace', size=12 , weight='bold')
plt.xticks( xs, [ n  for i,n,e in sel0 ], rotation='vertical', fontproperties=font_prop )  # set xticks to names of pairs in this class

plt.grid()
plt.ylabel( "Energy [kcal/mol]", fontproperties=font_prop ) 
plt.ylim(0.,+60)
plt.xlim(-0.5,len(sel0)-0.5)
plt.yticks(np.arange( 0.0, 45.0+0.001, 5.0 ))

plt.legend()
plt.tight_layout()
plt.savefig( "Econtrast_2endgroup.png", bbox_inches='tight' )
plt.savefig( "Econtrast_2endgroup.svg", bbox_inches='tight' )




plt.show()