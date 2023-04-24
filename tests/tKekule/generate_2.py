#! /bin/python3

import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

sys.path.append("../../")
#from pyBall import Kekule as kek
#from pyBall import atomicUtils as au
from pyBall import plotUtils   as pu




'''
      Donor         Acceptor
-----------------------------
1    -OH,-NH2             
2        -NH-         -O-,=O
3                     =N-
'''


groups={
'-CH3':1,'=CH2':2,'-CH2-':2,'=CH-':3,'C':4,'C/N':(3,4),
'-NH2':1,'=NH' :2,'-NH-' :2,'=N-' :3,
'-OH' :1,'=O'  :2,'-O-'  :2,
}

non=[None]
cdon =[ '-NH-',      ]    # central donors
cacc =[ '=N-', '-O-' ]
sdon =[ '-NH2'       ]
sacc =[ '=O'         ]
Cs=['=CH-']
Cc=['C']
CN=['C/N']

central = cdon+cacc
side    = sdon+sacc 
don     = cdon+sdon
acc     = cacc+sacc 




# ------- Sites
sites={
"A1": [(-2.4, 0.0),  side+central+non   ],
"A2":[( 0.0, 0.0),       central       ], 
"A3":[( 2.4, 0.0),  side+central+non   ],
"B1":[(-3.6,-0.7), Cs ],
"B2":[(-1.2,-0.7), Cc ], 
"B3":[( 1.2,-0.7), Cc ], 
"B4":[( 3.6,-0.7), Cs ],
"C1":[(-3.6,-2.1), Cs ], 
"C2":[(-1.2,-2.1), Cc ], 
"C3":[( 1.2,-2.1), Cc ], 
"C4":[( 3.6,-2.1), Cs ],
"D1":[(-2.4,-2.8), CN+non ], 
"D2":[( 0.0,-2.8), CN+non ], 
"D3":[( 2.4,-2.8), CN+non ],
}


sites_can={
"A1":[(-2.4, 0.0),  0     ],
"A2":[( 0.0, 0.0),  1     ], 
"A3":[( 2.4, 0.0),  2     ],
"B1":[(-3.6,-0.7), "=CH-" ],
"B2":[(-1.2,-0.7), "C"    ], 
"B3":[( 1.2,-0.7), "C"    ], 
"B4":[( 3.6,-0.7), "=CH-" ],
"C1":[(-3.6,-2.1), "=CH-" ], 
"C2":[(-1.2,-2.1), "C"    ], 
"C3":[( 1.2,-2.1), "C"    ], 
"C4":[( 3.6,-2.1), "=CH-" ],
"D1":[(-2.4,-2.8), "C/N"  ], 
"D2":[( 0.0,-2.8), "C/N"  ], 
"D3":[( 2.4,-2.8), "C/N"  ],
}


# ====== Generate skeletons (no-types assignment)

atoms0 = ["A2","B2","B3","C2","C3"]
bonds0 = [ ("A2","B2"),("A2","B3"), ("B2","C2"), ("B3","C3") ]

nmol=0
skeletons = [] 
for h1 in [0,1,5,6]:
    atoms1 = []
    bonds1 = []
    if   h1==0: # missing
        pass
    elif h1==1: # terminating
        atoms1.append( "A1" )
        bonds1.append( ("A1","B2") )
    elif h1==5: # pentagon
        atoms1 += [ "A1","B1","C1" ]
        bonds1 += [ ("A1","B2"),("A1","B1"),("B1","C1"),("C1","C2"), ]
    elif h1==6: # hexagon
        atoms1 += [ "A1","B1","C1","D1" ]
        bonds1 += [ ("A1","B2"),("A1","B1"),("B1","C1"),("C1","D1"),("D1","C2"), ]
            
    for h3 in [0,1,5,6]:

        if(h3>h1): continue    # prevent symmetric variants

        atoms3 = []
        bonds3 = []
        if   h3==0: # missing
            pass
        elif h3==1: # end
            atoms3.append( "A3" )
            bonds3.append( ("A3","B3") )
        elif h3==5: # pentagon
            atoms3 += [ "A3","B4","C4" ]
            bonds3 += [ ("A3","B3"),("A3","B4"),("B4","C4"),("C4","C3"), ]
        elif h3==6: # hexagon
            atoms3 += [ "A3","B4","C4","D3" ]
            bonds3 += [ ("A3","B3"),("A3","B4"),("B4","C4"),("C4","D3"),("D3","C3"), ]

        for h2 in [5,6]:
            atoms2 = []
            bonds2 = []
            if   h2==5:  # pentagon
                bonds2 += [ ("C2","C3") ]
            elif h2==6: # hexagon
                atoms2 += [ "D2" ]
                bonds2 += [ ("D2","C2"),("D2","C3") ]

            nmol+=1
            print( nmol, h1,h3,h2 )
            
            molecule = ( atoms0+atoms1+atoms2+atoms3,  bonds0+bonds1+bonds2+bonds3, (h1,h2,h3),  )
            skeletons.append( molecule )

print( len(skeletons) )

def assing_type( val, key, aaa ):
    if isinstance(val, int):
        return aaa[val]
    else:
        return val

molecules = []
for i,mol in enumerate(skeletons):
    atoms,bonds,variant = mol
    v1,v2,v3=variant
    dct={ k:i for i,k in enumerate(atoms) }
    ibonds = [ (dct[b[0]],dct[b[1]]) for b in bonds ]
    apos   = np.array([ sites[k][0] for k in atoms ])

    if   v1==0:  # missing 
        groups1=[None] 
    elif v1==1:  # end
        groups1=side
    else:  # cycle
        groups1=central

    if   v3==0:  # missing 
        groups3=[None] 
    elif v3==1:  # end
        groups3=side
    else:  # cycle
        groups3=central

    group2 = central

    for a1 in groups1:
        for a3 in groups3:
            for a2 in group2:
                aaa = (a1,a2,a3)
                atypes= [ assing_type( sites_can[k][1], k, aaa ) for k in atoms ] 
                mol = ( atypes, apos, ibonds )
                molecules.append(mol)


print(len(molecules))

sz=20
plt.figure( figsize=(sz,sz) )
iimg = 0
for i,mol in enumerate(molecules):
    atypes, apos, ibonds = mol
    print( "#### Molecule[%i] " %i, atypes )

    plt.subplot(5,5,iimg+1)
    pu.plotAtoms( apos, labels=atypes, axes=(0,1) )
    pu.plotBonds( links=ibonds, ps=apos, axes=(0,1) )
    plt.axis('equal')
    plt.xlim(-4.,4.)
    plt.ylim(-7.,1.)
    plt.title("Mol[%i]" %i )
    iimg+=1
    if(iimg >= 25):
        plt.savefig( "endgroups_%03i_%03i.png" %(i-25+1,i+1) , bbox_inches='tight')  
        iimg=0  
        plt.figure( figsize=(sz,sz) )


'''
plt.figure( figsize=(10,10) )
for i,mol in enumerate(skeletons):
    atoms,bonds,variant = mol
    print( "#### Molecule[%i] " %i, variant, atoms )
    #print( "atoms: ", atoms )
    #print( "bonds"  , bonds )

    plt.subplot(5,5,i+1)
    
    na   = len(atoms)
    apos = np.array([ sites[k][0] for k in atoms ])
    es   = ['C']*na
    dct={ k:i for i,k in enumerate(atoms) }
    #print( "dct: ", dct )

    ibonds = [ (dct[b[0]],dct[b[1]]) for b in bonds ]
    
    #print("ibonds", ibonds)
    pu.plotAtoms( apos, es=es, axes=(0,1) )
    pu.plotBonds( links=ibonds, ps=apos, axes=(0,1) )
    plt.axis('equal')
    plt.xlim(-4.,4.)
    plt.ylim(-7.,1.)
plt.show()
'''
