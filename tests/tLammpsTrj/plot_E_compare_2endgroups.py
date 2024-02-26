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
#f_egs="/home/prokop/Desktop/CARBSIS/Paolo/endgroups/B3LYP_oriented/"
#dir_geom="/home/prokop/Desktop/CARBSIS/Paolo/B3LYP_finished/pairs_b3lyp/confs/"
path = "/home/prokop/Desktop/CARBSIS/PEOPLE/Paolo/"

f_egs   =path+"endgroups/B3LYP_oriented/"
dir_geom=path+"B3LYP_finished/pairs_b3lyp/confs/"

sel_names = {
"HHH-h-p": 'r',
"HHH-h"  : 'b',
}



# ============ Main


#_,Es,names,pair_names          = bpu.read_dat( "/home/prokop/Desktop/CARBSIS/Paolo/pairs_compare_DFTB_vs_B3LYP/4-dimers_energy_sorted.dat", ni=1, nf=3, iname=0, toRemove=bases_to_remove )
_,Es,names,pair_names           = bpu.read_dat( path+"B3LYP_finished/1-binding_energy_all.dat",   ni=4, nf=2, iname=0 )
_, Angs,   names_angs,_         = bpu.read_dat( path+"B3LYP_finished/dihedral_angle_ave.dat",   ni=0, nf=1, iname=0 )
_, Hbonds, names_hbonds,_       = bpu.read_dat( path+"B3LYP_finished/bond_length_min.dat",      ni=0, nf=1, iname=0 )
_, Hbonds_av, names_hbonds_av,_ = bpu.read_dat( path+"B3LYP_finished/bond_length_ave.dat",      ni=0, nf=1, iname=0 )

old_inds = np.array(  range(len(Es)) )
Angs = np.array( Angs ); print( "Angs ", Angs ); #exit()
#intervals = np.arange( np.min(Es), np.max(Es), 5.0 ); 
intervals = np.arange( 0, 50+1, 5.0 ); print( "intervals ", intervals   ); #exit()


Emins =   bpu.find_minim_energy_confs( Es, pair_names, Emax=1000, ipivot=0 )                          # find minimum energy conformers for each pair

EBs   = np.array([-Emins[k][0]                                                    for k in  pair_names ])     # compute binding energy  for each pair as numpy array
ECs   = np.array([-Emins[k][0]+0.5*( Emins[(k[0],k[0])][0]+Emins[(k[1],k[1])][0]) for k in  pair_names ])     # compute energy contrast for each pair as numpy array


pair_sel = [( i,(n1,n2), sel_names[n1] if n1 in sel_names else sel_names[n2]) for i,(n1,n2) in enumerate(pair_names) if n1 in sel_names or n2 in sel_names]
pair_sel = sorted( pair_sel, key=lambda k: EBs[k[0]] )           # sort pairs by energy of the minimum energy conformer    (k[1] is the name of the pair)

inds, names, colors = zip(*pair_sel)

names = [ (n1,n2) if n1 in sel_names else (n2,n1) for n1,n2 in names  ]


inds= np.array(inds)

#inds     = [ i for i,n in pair_sel ]                                      # get indices of the pairs in the list pair_names

print( pair_sel )
print( inds )

EBs_  = EBs[inds]
ECs_  = ECs[inds]
Angs_ = Angs[inds]
old_inds_ = old_inds[inds]

#best_in_group =  find_in_interval( intervals, EBs_, ECs_, nbest=3 )   ;print( "best_in_group ", best_in_group )


for ii,i in enumerate(inds):
    print( i, names[ii], EBs[i], ECs[i], Angs[i] )

# ====== Figure 2 Plotting: 1D plot of energies for each pair in this class
font_prop = FontProperties(family='monospace', size=12 , weight='bold')

#plt.figure( figsize=(len(ns)*0.2,4) )
plt.figure( figsize=(len(pair_sel)*0.18,5) )
xs=range(len(pair_sel))
plt.xticks( xs, [ f"{n1:<10} {n2:<8}" for (n1,n2) in names ], rotation='vertical', fontproperties=font_prop )  # set xticks to names of pairs in this class
for l, c in zip(plt.gca().get_xticklabels(), colors):
    l.set_color(c)

plt.plot( EBs_, '-b',  lw=4, label='E_bind'     )   # plot energies of the minimum energy conformer for each pair in this class
plt.plot( ECs_, '.-r', lw=2, label='E_Contrast' )   # plot energies of the minimum energy conformer for each pair in this class

plt.plot( Angs_, 'o-g', lw=1, label='Angs' )  

plt.grid()
plt.ylabel( "Energy [kcal/mol]", fontproperties=font_prop ) 
plt.ylim(0.,+45)
plt.xlim(-0.5,len(pair_sel)-0.5)
plt.yticks(np.arange( 0.0, 45.0+0.001, 5.0 ))



plt.tight_layout()
#plt.savefig( "Econtrast_trashold.png", bbox_inches='tight' )
#plt.savefig( "Econtrast_trashold.svg", bbox_inches='tight' )




plt.show()