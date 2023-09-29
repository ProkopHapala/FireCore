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
Econtrast_min  = 5.0

dEcontrast_min = 10.0

class_prop__ = [ 
    ('D ... A'     ,'y'),
    ('DD ... AA'   ,'brown'),
    ('DD ... AAA'  ,'m'),
    ('DDD ... AA'  ,'c'),
    ('DDD ... AAA' ,'r'),
    ('DAD ... ADA' ,'b'),
    ('DDA ... AAD' ,'g'), 
] 

class_prop_ = [ 
    ('H_e'     ,'y'),
    ('HH_ee'   ,'brown'),
    ('HH_eee'  ,'m'),
    ('ee_HHH'  ,'c'),
    ('HHH_eee' ,'r'),
    ('HeH_eHe' ,'b'),
    ('HHe_Hee' ,'g'), 
] 
class_prop = dict( class_prop_ )



#f_geoms="/home/prokop/Desktop/CARBSIS/Paolo/pairs_compare_DFTB_vs_B3LYP/2-plot_confs/dirs/"
#f_egs="/home/prokop/Desktop/CARBSIS/Paolo/endgroups/endgroups/"
f_egs="/home/prokop/Desktop/CARBSIS/Paolo/endgroups/B3LYP_oriented/"

dir_geom="/home/prokop/Desktop/CARBSIS/Paolo/B3LYP_finished/pairs_b3lyp/confs/"




# ============ Functions

def find_in_interval( intervals, EBs, ECs, nbest=1 ):
    isels = []
    print( "EBs ", EBs )
    print( "ECs ", ECs )
    for i in range(len(intervals)-1):
        emin=intervals[i]
        emax=intervals[i+1]
        isel = [ j for j in range(len(EBs)) if EBs[j]>emin and EBs[j]<emax ]
        isel.sort( key=lambda k: EBs[k]-ECs[k] )
        isels.append( isel[:nbest] )
    return isels



# ============ Main

bases_to_remove = set( ['NNO-hp','ONO-p','NO-h-p'] )

#_,Es,names,pair_names = bpu.read_dat( "/home/prokop/Desktop/CARBSIS/Paolo/pairs_compare_DFTB_vs_B3LYP/4-dimers_energy_sorted.dat", ni=1, nf=3, iname=0, toRemove=bases_to_remove )
_,Es,names,pair_names = bpu.read_dat( "/home/prokop/Desktop/CARBSIS/Paolo/B3LYP_finished/1-binding_energy_all.dat",   ni=4, nf=2, iname=0, toRemove=bases_to_remove )

old_inds = np.array(  range(len(Es)) )

_, Angs,   names_angs,_         = bpu.read_dat( "/home/prokop/Desktop/CARBSIS/Paolo/B3LYP_finished/dihedral_angle_ave.dat",   ni=0, nf=1, iname=0, toRemove=bases_to_remove )
_, Hbonds, names_hbonds,_       = bpu.read_dat( "/home/prokop/Desktop/CARBSIS/Paolo/B3LYP_finished/bond_length_min.dat",      ni=0, nf=1, iname=0, toRemove=bases_to_remove )
_, Hbonds_av, names_hbonds_av,_ = bpu.read_dat( "/home/prokop/Desktop/CARBSIS/Paolo/B3LYP_finished/bond_length_ave.dat",      ni=0, nf=1, iname=0, toRemove=bases_to_remove )


Angs = np.array( Angs ); print( "Angs ", Angs ); #exit()

#intervals = np.arange( np.min(Es), np.max(Es), 5.0 ); 
intervals = np.arange( 0, 50+1, 5.0 ); print( "intervals ", intervals   ); #exit()


Emins =   bpu.find_minim_energy_confs( Es, pair_names, Emax=1000, ipivot=0 )                          # find minimum energy conformers for each pair

EBs   = np.array([-Emins[k][0]                                                    for k in  pair_names ])     # compute binding energy  for each pair as numpy array
ECs   = np.array([-Emins[k][0]+0.5*( Emins[(k[0],k[0])][0]+Emins[(k[1],k[1])][0]) for k in  pair_names ])     # compute energy contrast for each pair as numpy array




#for n,e in Emins.items(): print( n, e )

pair_sel = [ (i,n) for i,n in enumerate(pair_names) if ECs[i]>Econtrast_min ]    # select pairs with energy contrast above Econtrast_min   (i is the index of the pair in the list pair_names)

pair_sel = sorted( pair_sel, key=lambda k: Emins[k[1]][0] )           # sort pairs by energy of the minimum energy conformer    (k[1] is the name of the pair)
inds     = [ i for i,n in pair_sel ]                                      # get indices of the pairs in the list pair_names
colors   = [  class_prop.get(  bpu.name_to_class( n[0] ) + "_" + bpu.name_to_class( n[1] ), 'k' ) for i,n in pair_sel ]   # assign colors to each pair according to the class of the pair


EBs_  = EBs[inds]
ECs_  = ECs[inds]
Angs_ = Angs[inds]
old_inds_ = old_inds[inds]

#best_in_group =  find_in_interval( intervals, EBs_, ECs_, nbest=3 )   ;print( "best_in_group ", best_in_group )

explocit_sel = set( [('HHH-h-p','NNN-hhh')] )

dECmin=2.0
best_in_group =  [ i for i in range(len(inds)) if (((EBs_[i]-ECs_[i])<dECmin) or (pair_sel[i][1] in explocit_sel ))  ] ; print( "best_in_group ", best_in_group )

old_inds_best = old_inds_[best_in_group]

pair_names_best = [ pair_names[i] for i in old_inds_best ]




# ====== Figure 2 Plotting: 1D plot of energies for each pair in this class
font_prop = FontProperties(family='monospace', size=12 , weight='bold')

#plt.figure( figsize=(len(ns)*0.2,4) )
plt.figure( figsize=(len(pair_sel)*0.18,5) )
xs=range(len(pair_sel))
plt.xticks( xs, [ f"{n1:<10} {n2:<8}"  for i,(n1,n2) in pair_sel ], rotation='vertical', fontproperties=font_prop )  # set xticks to names of pairs in this class
for l, c in zip(plt.gca().get_xticklabels(), colors):
    l.set_color(c)

plt.plot( EBs_, '-b',  lw=4, label='E_bind'     )   # plot energies of the minimum energy conformer for each pair in this class
plt.plot( ECs_, '.-r', lw=2, label='E_Contrast' )   # plot energies of the minimum energy conformer for each pair in this class

plt.plot( Angs_, 'o-g', lw=1, label='Angs' )  

#plt.scatter( xs, EBs[inds],  c=colors   ,zorder=5  )              # plot energies of the minimum energy conformer for each pair in this class
#plt.scatter( xs, ECs[inds],  c=colors   ,zorder=5, marker="o"  )   # plot energies of the minimum energy conformer for each pair in this class


plt.plot( best_in_group, EBs_[best_in_group], 'o', ms=10, markerfacecolor="None", markeredgecolor='k' ) 

#for ibest in best_in_group:
#    plt.plot( ibest, EBs_[ibest], 'o', ms=10, markerfacecolor="None", markeredgecolor='k' ) 
#    plt.plot( ibest, ECs_[ibest], 'o', ms=10, markerfacecolor="None", markeredgecolor='k' ) 

plt.grid()
plt.ylabel( "Energy [kcal/mol]", fontproperties=font_prop ) 
plt.ylim(0.,+45)
plt.xlim(-0.5,len(pair_sel)-0.5)
plt.yticks(np.arange( 0.0, 45.0+0.001, 5.0 ))



plt.tight_layout()
plt.savefig( "Econtrast_trashold.png", bbox_inches='tight' )
plt.savefig( "Econtrast_trashold.svg", bbox_inches='tight' )

plt.figure();
for cls in class_prop__:
    plt.plot( [0,1], [0,1], 'o', c=cls[1], label=cls[0] )
#plt.legend(fontsize='large')
font = {
    #'family': 'serif',
    #'color':  'darkred',
    'weight': 'bold',
    'size': 16,
}
plt.legend(prop=font)
plt.savefig( "Legend_classes.png", bbox_inches='tight' )
plt.savefig( "Legend_classes.svg", bbox_inches='tight' )
plt.close()

'''
for i in range( len(best_in_group) ):
    print( i, old_inds_best[i],  names[old_inds_best[i]] , pair_names_best[i] )
    name = names[old_inds_best[i]]
    mol = au.AtomicSystem( dir_geom + name+".xyz" )
    mol.findBonds()
    hbs,rhbs = mol.findHBonds( bPrint=True )     #;print( hbs, rhbs )
    n1,n2 = pair_names_best[i]
    fig = plt.figure( figsize=(3,4) )
    plu.plotSystem ( mol, bLabels=False, HBs=(hbs,rhbs) )
    #plt.title( n1,n2, EBs_[best_in_group[i]], ECs_[best_in_group[i]]  )
    fig.patch.set_visible(False)
    sEB = ("\nE_XY=%6.1f kcal/mol" %EBs_[best_in_group[i]] )
    sEC = ("\nE_C=%6.1f kcal/mol"  %ECs_[best_in_group[i]] )
    plt.title( n1+" "+n2+sEB+sEC )
    plt.tight_layout()
    plt.gca().axis('off')
    plt.savefig( ("%03i_" %i )+name+".svg", bbox_inches='tight' )
    plt.close()
'''

plt.show()