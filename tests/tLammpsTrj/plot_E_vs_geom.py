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


nicknames = {
    'HHO-h-p_1' : 'guanine',
    'HNO-h'     : 'cytosine',
    'HN-h-p'    : 'adenine',
    'OHO-h_2'   : 'thymine',
    'OHO-h_1'   : 'uracil',
}



f_geoms="/home/prokop/Desktop/CARBSIS/Paolo/pairs_compare_DFTB_vs_B3LYP/2-plot_confs/dirs/"

#f_egs="/home/prokop/Desktop/CARBSIS/Paolo/endgroups/endgroups/"

f_egs="/home/prokop/Desktop/CARBSIS/Paolo/endgroups/B3LYP_oriented/"

# ============ Main

bases_to_remove = set( ['NNO-hp','ONO-p','NO-h-p'] )

#_,Es,names,pair_names = bpu.read_dat( "/home/prokop/Desktop/CARBSIS/Paolo/pairs_compare_DFTB_vs_B3LYP/4-dimers_energy_sorted.dat", ni=1, nf=3, iname=0, toRemove=bases_to_remove )
_,Es,names,pair_names = bpu.read_dat( "/home/prokop/Desktop/CARBSIS/Paolo/B3LYP_finished/1-binding_energy_all.dat",   ni=4, nf=2, iname=0, toRemove=bases_to_remove )

_, Angs,   names_angs,_   = bpu.read_dat( "/home/prokop/Desktop/CARBSIS/Paolo/B3LYP_finished/dihedral_angle_ave.dat",   ni=0, nf=1, iname=0, toRemove=bases_to_remove )
_, Hbonds, names_hbonds,_ = bpu.read_dat( "/home/prokop/Desktop/CARBSIS/Paolo/B3LYP_finished/bond_length_min.dat",      ni=0, nf=1, iname=0, toRemove=bases_to_remove )



#names__  = [ bpu.split_pair_with_S( n1 ) for n1,n2 in names ]   # convert names format from HHH-hhS1_NNO-hpS1  to (HHH-hh,NNO-hp)
#names__  = [ (n1,n2) for n1,n2 in names__ if not ((n1 in bases_to_remove) or (n2 in bases_to_remove)) ]   # filter out pairs from bases_to_remove

#pair_names_  = [  n1+"_"+n2 for n1,n2 in pair_names ]                  # convert names format from (HHH-hh,NNO-hp) to HHH-hh_NNO-hp

for i,(n1,n2) in enumerate(pair_names): print( i,n1,n2 )

#exit()
#ng_set = set( s for pair in names__ for s in pair ) 
#print(ng_set)

pair_classes = {  ( bpu.name_to_class( n1 ) + "_" + bpu.name_to_class( n2 ) ):[] for n1,n2 in pair_names }    ;print( "pair_classes \n", pair_classes )

pair_name_set = set(pair_names)

for n1,n2 in pair_name_set:
    n = bpu.name_to_class( n1 ) + "_" + bpu.name_to_class( n2 )
    #pair_classes[n].append( n1+"_"+n2 )
    pair_classes[n].append( (n1,n2) )
for n in pair_classes: print( n,"\n", pair_classes[n] )

Emins = bpu.find_minim_energy_confs( Es, pair_names, Emax=1000, ipivot=0 )
for n,e in Emins.items():  print(n) # print( n, e )

sel_classes =[
    #'H_e',
    'HH_ee',
    #'HH_eee',  
    #'ee_HHH',
    #'HHH_eee', 
    #'HeH_eHe',  
    #'HHe_Hee', 
]


font_prop = FontProperties(family='monospace', size=10)

# Plot energies for selected classes of pairs ordered by energy of the minimum energy conformer 
for n in sel_classes:        # for each class of pairs
    ns = pair_classes[n]     # list of pairs in this class
    ns.sort( key=lambda k: Emins[k][0] )   # sort names of pairs in this class by energy of the minimum energy conformer

    print( "#=========", n )
    nps = [ (pair_names[Emins[k][1]], Emins[k][0], Emins[k][0] - 0.5*( Emins[(k[0],k[0])][0]+Emins[(k[1],k[1])][0]) ) for k in ns ]

    #for (n1,n2),e in nps: print( n1, n2, e )

    # ====== Figure 1 Evluation: 2D plot of energies for each pair in this class

    donors   = { n1:[1000.0,0.0,0] for (n1,n2),e,ec in nps }
    acceptors= { n2:[1000.0,0.0,0] for (n1,n2),e,ec in nps }
    #comps = { nn:[0,0] for p,e in nps for nn in p }   #;print( comps )
    #for (n1,n2),e in nps: 
    #    l1=comps[n1]; l1[0]+=e; l1[1]+=1
    #    l2=comps[n2]; l2[0]+=e; l2[1]+=1

    # sort donors and acceptors by energy (average over all pairs)
    for (n1,n2),e,ec in nps: 
        l1=donors   [n1]; l1[0]=min(l1[0],e); l1[1]+=e; l1[2]+=1    # find minimum energy and average energy for each donor
        l2=acceptors[n2]; l2[0]=min(l2[0],e); l2[1]+=e; l2[2]+=1    # find minimum energy and average energy for each acceptor
    for nn in donors:    donors   [nn][1] /= donors   [nn][2]  # finish average energy for each donor
    for nn in acceptors: acceptors[nn][1] /= acceptors[nn][2]  # finish average energy for each acceptor
    #donors_sorted    = sorted( donors   .items(), key=lambda kv: kv[1][0] ) # sort donors by minimum energy 
    #acceptors_sorted = sorted( acceptors.items(), key=lambda kv: kv[1][0] )
    donors_sorted    = sorted( donors   .items(), key=lambda kv: kv[1][1] ) # sort donors by average energy
    acceptors_sorted = sorted( acceptors.items(), key=lambda kv: kv[1][1] )
    #for i,(nn,e) in enumerate(donors_sorted):       print( "D#%i " %i, nn, e )
    #for i,(nn,e) in enumerate(acceptors_sorted):    print( "A#%i " %i, nn, e )


    donors_d    = { k[0]:i for i,k in enumerate( donors_sorted)    }
    acceptors_d = { k[0]:i for i,k in enumerate( acceptors_sorted) }
    Eg = np.zeros( (len(donors_d),len(acceptors_d)) ); Eg.fill( np.nan )
    Ec = np.zeros( (len(donors_d),len(acceptors_d)) ); Eg.fill( np.nan )
    
    Eg_xx = np.expand_dims( np.array( [  Emins[n,n][0] for n,e in donors_sorted    ] ), axis=1)
    Eg_yy = np.expand_dims( np.array( [  Emins[n,n][0] for n,e in acceptors_sorted ] ), axis=0)

    #print( donors_d )
    #print( acceptors_d )

    # ====== Figure 1 Plotting: 2D plot of energies for each pair in this class

    #plt.figure( figsize=(10,10) )

    sz=0.75
    plt.figure( figsize=(len(acceptors_d)*sz+3*sz, len(donors_d)*sz+2*sz) )
    # Plot geometries of constituents of each pair in this class
    mol_size = 10
    donor_geoms    = [ au.AtomicSystem( f_egs+n+".xyz" ) for n in donors_d    ] 
    acceptor_geoms = [ au.AtomicSystem( f_egs+n+".xyz" ) for n in acceptors_d ]
    for i in range(len(donors_d)):    # donors along y axis
        g = donor_geoms[i]
        g.apos[:,0] += (i+0.5)*mol_size
        g.apos[:,1] += -mol_size*0.5
        plu.plotSystem( g, axes=(1,0), bLabels=False )
    for i in range(len(acceptors_d)): # acceptors along x axis
        g = acceptor_geoms[i]
        g.apos[:,1] += (i+0.5)*mol_size
        g.apos[:,0] += -mol_size*0.5
        plu.plotSystem( g, axes=(1,0), bLabels=False )   
    extent = (0, mol_size*(len(acceptors_d)), 0, mol_size*(len(donors_d)) )
    #extent = (0, mol_size*(len(donors_d)),    0,   mol_size*(len(acceptors_d)) )

    plt.title( n )                                       # set title to name of this class
    for (n1,n2),e,ec in nps: 
        Eg[ donors_d[n1], acceptors_d[n2] ] = e
        Ec[ donors_d[n1], acceptors_d[n2] ] = ec

    vmin=-np.nanmax(Eg)
    vmax=-np.nanmin(Eg)
    print( "vmin, vmax ",  vmin, vmax )

    plt.imshow( -Eg_xx, vmin=vmin, vmax=vmax, cmap='jet', interpolation='nearest', origin='lower', extent=[ -mol_size, 0, extent[2],extent[3] ] )
    plt.imshow( -Eg_yy, vmin=vmin, vmax=vmax, cmap='jet', interpolation='nearest', origin='lower', extent=[ extent[0],extent[1], -mol_size, 0 ] )
    plt.imshow( -Eg   , vmin=vmin, vmax=vmax, cmap='jet', interpolation='nearest', origin='lower', extent=extent )
    #plt.imshow( Eg.transpose(), cmap='jet', interpolation='nearest', origin='lower', extent=extent )
    plt.xticks( [(i+0.5)*mol_size for i in range(len(acceptors_d))], [ n+bpu.try_nickname(n,nicknames) for n,e in acceptors_sorted ], rotation='vertical'  , fontproperties=font_prop   )  # set xticks to names of pairs in this class
    plt.yticks( [(i+0.5)*mol_size for i in range(len(donors_d   ))], [ n+bpu.try_nickname(n,nicknames) for n,e in donors_sorted    ], rotation='horizontal', fontproperties=font_prop )  # set xticks to names of pairs in this class

    #plt.axhline(0,lw=2,c='k',ls='--')
    #plt.axvline(0,lw=2,c='k',ls='--')
    plt.axhline(0,lw=2,c='k',ls='-')
    plt.axvline(0,lw=2,c='k',ls='-')
    #plt.axhline(0,lw=2,c='gray',ls='-')
    #plt.axvline(0,lw=2,c='gray',ls='-')

    for i in range(len(donors_d)): 
        for j in range(len(acceptors_d)):
            plt.text( (j+0.5)*mol_size, (i+0.5)*mol_size, "%3.0f" %-Eg[i,j], ha="center", va="top",    color="gray" )
            plt.text( (j+0.5)*mol_size, (i+0.5)*mol_size, "%3.0f" %-Ec[i,j], ha="center", va="bottom", color="gray", weight='bold' )

    for i in range(len(donors_d)): 
        plt.text( -mol_size, (i+1)*mol_size, "%3.0f" %-Eg_xx[i,0], ha="left", va="top", color="gray" )
    
    for j in range(len(acceptors_d)):
        plt.text( (j)*mol_size, 0, "%3.0f" %-Eg_yy[0,j], ha="left", va="top", color="gray" )


    #plt.axis('equal')
    plt.xlim(-mol_size, mol_size*(len(acceptors_d)) )
    plt.ylim(-mol_size, mol_size*(len(donors_d   )) )
    plt.colorbar()
    plt.tight_layout()

    plt.savefig( n+"_2D.png", bbox_inches='tight' )

    
    # ====== Figure 2 Plotting: 1D plot of energies for each pair in this class
    
    plt.figure( figsize=(len(ns)*0.2+1,4) )
    plt.title( n )                                       # set title to name of this class
    #plt.xticks(range(len(ns)), [ f"{n1:>8} {n2:>8}"  for n1,n2 in nps ], rotation='vertical', fontproperties=font_prop )  # set xticks to names of pairs in this class
    plt.xticks(range(len(ns)), [ f"{n1:<10} {n2:<8}"  for (n1,n2),e,ec in nps ], rotation='vertical', fontproperties=font_prop )  # set xticks to names of pairs in this class
    plt.plot( [-Emins[k][0] for k in ns], '.-', label='Ebind' )          # plot energies of the minimum energy conformer for each pair in this class
    plt.plot( [-Emins[k][0]+0.5*( Emins[(k[0],k[0])][0]+Emins[(k[1],k[1])][0]) for k in ns], '.-', label='Contrast' )          # plot energies of the minimum energy conformer for each pair in this class
    
    plt.plot( [ Angs  [ Emins[k][1] ] for k in ns], '.-', label='Angles' )   
    plt.plot( [ Hbonds[ Emins[k][1] ] for k in ns], '.-', label='HBonds' )   
    
    plt.grid()
    plt.ylabel( "Energy [kcal/mol]" ) 
    #plt.ylim(-60.,0.)
    #plt.ylim(-45.,0.)
    plt.ylim(0.,+45)
    #plt.yticks(np.arange( -45.0, 0.00001, 5.0 ))
    plt.yticks(np.arange( 0.0, 45.0+0.001, 5.0 ))
    plt.tight_layout()
    plt.savefig( n+".png", bbox_inches='tight' )
    
    '''
    for k in pair_classes[n]:
        Exy = Emins[k][0]
        Exx = Emins[ (k[0],k[0]) ][0]
        Eyy = Emins[ (k[1],k[1]) ][0]
        print( k, Exx, Exx - 0.5*( Exx + Eyy ), "    |   ",  Exx, Eyy  )  
    '''

plt.show()