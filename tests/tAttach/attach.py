

import sys
import os

sys.path.append('../../')
from pyBall             import plotUtils   as plu
from pyBall             import atomicUtils as au
from pyBall.atomicUtils import AtomicSystem

import numpy as np
import matplotlib.pyplot as plt
#from functools import partial

# ============ MAIN
# A - H-N(12-1)    up=(1,10)
# T - H-N(11-3)    up=(12,8)
# backbone  F-C(1-2)    Cl-C(17-9)

# =============== SETUP ==================

#       name                   attachment         H-Bonds           
#                             C/N  H    Up
groups = [                        
( "penta_hb3_acceptor2",  ( ( 10, 11, (14, 8) ),  [ 8, 6,14 ] ) ),
( "penta_hb3_acceptor" ,  ( ( 10, 11, ( 7, 8) ),  [ 8, 6, 7 ] ) ),
( "penta_hb3_donor"    ,  ( ( 13, 16, (13, 8) ),  [12, 9,14 ] ) ),

( "adenine_mod"        ,  ( (  4,  7, (13, 8) ),  [11, 9    ] ) ),
( "uracil_mod"         ,  ( (  3, 12, ( 5, 1) ),  [ 8,10    ] ) ), 

( "adenine"            ,  ( (  1, 12, (10,11) ),  [15, 7    ] ) ),
( "thymine"            ,  ( (  3, 11, (12, 8) ),  [ 9,12    ] ) ),
( "uracil"             ,  ( (  3, 10, ( 7, 8) ),  [ 8,11    ] ) ),     # = hexa_hb3_acceptor

( "hexa_hb3_donor"     ,  ( (  3, 15, (12, 8) ),  [11, 9,14 ] ) ),
( "uracil"             ,  ( (  3, 10, ( 7, 8) ),  [ 8,11, 7 ] ) ),     # = hexa_hb3_acceptor

( "guanine"            ,  ( (  3, 13, (11,10) ),  [10,14,16 ] ) ),
( "citosine"           ,  ( (  3, 13, ( 7, 8) ),  [ 9, 5, 7 ] ) ),

( "penta_hb2_acceptor" ,  ( (  3,  9, ( 6, 8) ),  [ 8, 6    ] ) ),     # other attachments  (3,9), (5,7), (2,4)
( "penta_hb2_donor"    ,  ( (  3,  4, ( 7, 8) ),  [12, 9    ] ) ),     # other attachments  (3,4), (5,7), (2,10)

( "naphta_hb2_acceptor",  ( ( 16, 10, (15, 1) ),  [ 9,13    ] ) ),     # other attachments  (16,10), (14,8), (4,7), (2,6), 
( "naphta_hb2_donor"   ,  ( (  2,  6, ( 1,12) ),  [14, 7    ] ) ),     # other attachments  (2,6),   (4,16), (10,13), 
]

pairs = [
("penta_hb3_acceptor2","penta_hb3_donor"),
("penta_hb3_acceptor" ,"penta_hb3_donor"),
("adenine_mod"        ,"uracil_mod"),
("adenine"            ,"uracil"),
("hexa_hb3_donor"     ,"uracil"),
("guanine"            ,"citosine"),
("penta_hb2_acceptor" ,"penta_hb2_donor"),
("naphta_hb2_acceptor","naphta_hb2_donor"),
]

# =============== Functions ==================

def plotGroup( G, inds, bFlip=False ):
    iis,his = inds
    his     = np.array(his)-1
    G.orient( iis[0], iis[2], (iis[0],iis[1]), trans=(1,2,0) )  
    if(bFlip): G.apos[:,0]*=-1
    plu.plotSystem( G, sz=1000. );  
    ps=G.apos[his];                 plt.scatter(ps[:,0],ps[:,1],color='k',zorder=5);  
    ps=G.apos[[iis[0]-1,iis[1]-1]]; plt.scatter(ps[:,0],ps[:,1],color=['b','g'],zorder=5);  
    #plt.title(name1)

# =============== MAIN ==================

group_dict = dict(groups)


# ------- Group attachent
#B = AtomicSystem(fname='backbone.xyz' )
B = AtomicSystem(fname='backbone.mol2' )


B.lvec = np.array( [[25.,0.,0.],[0.,5.,0.],[0.,0.,20.0]  ] )
for pair in pairs:
    print(pair)
    BB = B.clonePBC()
    name1,name2 = pair
    G1 = AtomicSystem(fname="endgroups/"+name1+".xyz"  )
    G2 = AtomicSystem(fname="endgroups/"+name2+".xyz"  )
    
    BB.preinitialize_atomic_properties()
    G1.preinitialize_atomic_properties()
    G2.preinitialize_atomic_properties()
    
    inds1, Hs1 = group_dict[name1]
    inds2, Hs2 = group_dict[name2]
    #BB.attach_group( G1, inds1[0], inds1[1], inds1[2], (1 ,2), up=(0.,0.,1.), _0=1  )
    #BB.attach_group( G2, inds2[0], inds2[1], inds2[2], (17,9), up=(0.,0.,1.), _0=1  )
    #BB.delete_atoms( [1-1,17-1] )
    BB.attach_group( G1, inds1[0], inds1[1], inds1[2], (18,8), up=(0.,0.,1.), _0=1 , pre="X" )
    BB.attach_group( G2, inds2[0], inds2[1], inds2[2], (17,6), up=(0.,0.,1.), _0=1 , pre="Y" )
    BB.delete_atoms( [17-1,18-1] )

    inds1 = BB.remap( [ "X"+str(i-1) for i in Hs1 ] )
    inds2 = BB.remap( [ "Y"+str(i-1) for i in Hs2 ] )
    #print( inds )
    comment = " Hbonds:X"+str(inds1)+"Y"+str(inds2)

    BB.print()
    BB.saveXYZ( "BB."+name1+"."+name2+".xyz", comment=comment )




'''
# ------- Plotting
for pair in pairs:
    name1,name2 = pair
    G1 = AtomicSystem(fname="endgroups/"+name1+".xyz"  )
    G2 = AtomicSystem(fname="endgroups/"+name2+".xyz"  )
    fig = plt.figure(figsize=(10.0,5.0))
    plt.subplot(1,2,1); plotGroup( G1, group_dict[name1], bFlip=True  ); plt.title(name1) 
    plt.subplot(1,2,2); plotGroup( G2, group_dict[name2], bFlip=False ); plt.title(name2)
    
    plt.savefig( name1+"."+name2+".png", bbox_inches='tight', )
    plt.close(fig)
    #G2 = AtomicSystem(fname="endgroups/"+name2+".xyz"  )
    #inds1 = group_dict[name1][0]
'''










'''
# -------- A-T pair
BB = AtomicSystem(fname='backbone.xyz' )
G1 = AtomicSystem(fname='endgroups/adenine.xyz'  )
G2 = AtomicSystem(fname='endgroups/thymine.xyz'  )
BB.attach_group( G1,  1, 12, (1,10),   (1 ,2),  up=(0.,0.,1.),  _0=1  )
BB.attach_group( G2,  3, 11, (12,8),   (17,9),  up=(0.,0.,1.),  _0=1  )
BB.delete_atoms( [1-1,17-1] )
BB.saveXYZ( "BB_A_T.xyz" )
'''

'''
# -------- A-T pair
BB = AtomicSystem(fname='backbone.xyz' )
G1 = AtomicSystem(fname='endgroups/penta_hb3_acceptor2.xyz'  )
G2 = AtomicSystem(fname='endgroups/penta_hb3_donor.xyz'  )
BB.attach_group( G1, 10, 11, (11,8),   (1 ,2),  up=(0.,0.,1.),  _0=1  )
BB.attach_group( G2, 13, 16, (13,8),   (17,9),  up=(0.,0.,1.),  _0=1  )
BB.delete_atoms( [1-1,17-1] )
BB.saveXYZ( "BB_G1_G2.xyz" )
'''




