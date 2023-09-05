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


# ============ Functions


f_geoms="/home/prokop/Desktop/CARBSIS/Paolo/pairs_compare_DFTB_vs_B3LYP/2-plot_confs/dirs/"
f_energy="/home/prokop/Desktop/CARBSIS/Paolo/pairs_compare_DFTB_vs_B3LYP/4-dimers_energy_sorted.dat"

#f_egs="/home/prokop/Desktop/CARBSIS/Paolo/endgroups/endgroups/"

f_egs="/home/prokop/Desktop/CARBSIS/Paolo/endgroups/B3LYP_oriented/"

# ============ Main

# Read energy file
f = open( f_energy, 'r' )


_,Es,names = bpu.read_dat( f_energy, ni=1, nf=3, nn=2 )

names__  = [ bpu.split_pair_with_S( n1 ) for n1,n2 in names ]   # convert names format from HHH-hhS1_NNO-hpS1  to (HHH-hh,NNO-hp)
names_   = [  n1+"_"+n2 for n1,n2 in names__ ]                  # convert names format from (HHH-hh,NNO-hp) to HHH-hh_NNO-hp


for i,(n1,n2) in enumerate(names__): print( i,n1,n2 )

#exit()

#ng_set = set( s for pair in names__ for s in pair ) 
#print(ng_set)

#pair_classes = set( bpu.name_to_class( n1 ) + "_" + bpu.name_to_class( n2 ) for n1,n2 in names__ )    ;print( pair_classes )

names__set = set(names__)

pair_classes ={
    'H_e'      : [],
    'HH_ee'    : [],
    'ee_HHH'   : [],
    'HH_eee'   : [],
    'HHH_eee'  : [],
    'HeH_eHe'  : [], 
    'HHe_Hee'  : [],
}

for n1,n2 in names__set:
    n = bpu.name_to_class( n1 ) + "_" + bpu.name_to_class( n2 )
    pair_classes[n].append( n1+"_"+n2 )
for n in pair_classes: print( n,"\n", pair_classes[n] )

Emins = bpu.find_minim_energy_confs( Es, names_, Emax=1000, ipivot=0 )
for n,e in Emins.items():  print(n) # print( n, e )

sel_classes =[
    #'H_e',
    'HH_ee',
    #'HH_eee',  
    #'ee_HHH',
    'HHH_eee', 
    'HeH_eHe',  
    'HHe_Hee', 
]


font_prop = FontProperties(family='monospace', size=10)

# Plot energies for selected classes of pairs ordered by energy of the minimum energy conformer 
for n in sel_classes:        # for each class of pairs
    ns = pair_classes[n]     # list of pairs in this class
    ns.sort( key=lambda k: Emins[k][0] )   # sort names of pairs in this class by energy of the minimum energy conformer

    print( "#=========", n )

    nps = [ (names__[Emins[k][1]], Emins[k][0]) for k in ns ]

    for (n1,n2),e in nps: print( n1, n2, e )

    donors   = { n1:[0,0] for (n1,n2),e in nps }
    acceptors= { n2:[0,0] for (n1,n2),e in nps }
    #comps = { nn:[0,0] for p,e in nps for nn in p }   #;print( comps )
    #for (n1,n2),e in nps: 
    #    l1=comps[n1]; l1[0]+=e; l1[1]+=1
    #    l2=comps[n2]; l2[0]+=e; l2[1]+=1

    # sort donors and acceptors by energy (average over all pairs)
    for (n1,n2),e in nps: 
        l1=donors   [n1]; l1[0]+=e; l1[1]+=1
        l2=acceptors[n2]; l2[0]+=e; l2[1]+=1
    for nn in donors:    donors   [nn] = donors   [nn][0]/donors   [nn][1]
    for nn in acceptors: acceptors[nn] = acceptors[nn][0]/acceptors[nn][1]
    donors_sorted    = sorted( donors   .items(), key=lambda kv: kv[1] )
    acceptors_sorted = sorted( acceptors.items(), key=lambda kv: kv[1] )
    for i,(nn,e) in enumerate(donors_sorted):       print( "D#%i " %i, nn, e )
    for i,(nn,e) in enumerate(acceptors_sorted):    print( "A#%i " %i,nn, e )


    donors_d    = { k[0]:i for i,k in enumerate( donors_sorted)    }
    acceptors_d = { k[0]:i for i,k in enumerate( acceptors_sorted) }
    Eg = np.zeros( (len(donors_d),len(acceptors_d)) ); Eg.fill( np.nan )
    
    print( donors_d )
    print( acceptors_d )

    plt.figure( figsize=(10,10) )

    # Plot geometries of constituents of each pair in this class
    mol_size = 10
    donor_geoms    = [ au.AtomicSystem( f_egs+n+".xyz" ) for n in donors_d    ] 
    acceptor_geoms = [ au.AtomicSystem( f_egs+n+".xyz" ) for n in acceptors_d ]
    for i in range(len(donors_d)):    # donors along x axis
        g = donor_geoms[i]
        g.apos[:,1] += (i+0.5)*mol_size
        g.apos[:,0] += -mol_size*0.5
        plu.plotSystem( g, axes=(0,1), bLabels=False )
    for i in range(len(acceptors_d)): # acceptors along y axis
        g = acceptor_geoms[i]
        g.apos[:,0] += (i+0.5)*mol_size
        g.apos[:,1] += -mol_size*0.5
        plu.plotSystem( g, axes=(0,1), bLabels=False )   

    extent = (0, mol_size*(len(acceptors_d)), 0, mol_size*(len(donors_d)) )
    #extent = (0, mol_size*(len(donors_d)),    0,   mol_size*(len(acceptors_d)) )

    plt.title( n )                                       # set title to name of this class
    for (n1,n2),e in nps: 
        Eg[ donors_d[n1], acceptors_d[n2] ] = e
    plt.imshow( Eg, cmap='jet', interpolation='nearest', origin='lower', extent=extent )
    #plt.imshow( Eg.transpose(), cmap='jet', interpolation='nearest', origin='lower', extent=extent )
    plt.xticks( [(i+0.5)*mol_size for i in range(len(acceptors_d))], [ n for n,e in acceptors_sorted ], rotation='vertical'  , fontproperties=font_prop   )  # set xticks to names of pairs in this class
    plt.yticks( [(i+0.5)*mol_size for i in range(len(donors_d   ))], [ n for n,e in donors_sorted    ], rotation='horizontal', fontproperties=font_prop )  # set xticks to names of pairs in this class

    plt.axis('equal')
    #plt.xlim(-mol_size, mol_size*(len(donors_d   )+1) )
    #plt.ylim(-mol_size, mol_size*(len(acceptors_d)+1) )
    plt.colorbar()
    plt.tight_layout()

    plt.savefig( n+"_2D.png", bbox_inches='tight' )

    
    #   1D plot of minimum energy conformer for each pair in this class
    plt.figure( figsize=(len(ns)*0.2+1,4) )
    plt.title( n )                                       # set title to name of this class
    #plt.xticks(range(len(ns)), [ f"{n1:>8} {n2:>8}"  for n1,n2 in nps ], rotation='vertical', fontproperties=font_prop )  # set xticks to names of pairs in this class
    plt.xticks(range(len(ns)), [ f"{n1:<10} {n2:<8}"  for (n1,n2),e in nps ], rotation='vertical', fontproperties=font_prop )  # set xticks to names of pairs in this class
    plt.plot( [Emins[k][0] for k in ns], 'o-' )          # plot energies of the minimum energy conformer for each pair in this class
    plt.grid()
    plt.ylabel( "Energy [kcal/mol]" ) 
    plt.ylim(-60.,0.)
    plt.tight_layout()
    plt.savefig( n+".png", bbox_inches='tight' )
    
    


    plt.plot

    for k in pair_classes[n]:
        print( k, Emins[k] )
        #print( k )

plt.show()