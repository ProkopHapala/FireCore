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


# ====== Generate skeletons (no-types assignment)

atoms0 = ["A2","B2","B3","C2","C3"]
bonds0 = [ ("A2","B2"),("A2","B3"), ("B2","C2"), ("B3","C3") ]

nmol=0
molecules = [] 
for h1 in [0,1,5,6]:
    atoms1 = []
    bonds1 = []
    if   h1==0: # nothing
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
        if   h3==0: # nothing
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
            molecules.append( molecule )

print( len(molecules) )


plt.figure( figsize=(10,10) )

for i,mol in enumerate(molecules):
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

#pu.plotAtoms( apos=apos, es=None, atoms=None, bNumbers=False, labels=None, sizes=100., colors='#808080', marker='o', axes=(0,1) ):
#pu.plotBonds( lps=None, links=None, ps=None, lws=None, axes=(0,1) ):



exit()




'''
center =([2,3,3,3,2], [(0,1),(1,2),(2,3),(3,4),(4,5)]  )
bridge2=(    [2,3,2], [(0,1),(1,2),(2,3)]              )
bridge1=(      [2,2], [(0,1),(1,2)]                    ) 
'''

'''
kek.setVerbosity(2,1)
kek.setVerbosity(0,0)
kek.init( len(AOs), len(b2a), b2a, AOs, seed=np.random.randint(15454) )
kek.init_buffers()
kek.getBuffs()
kek.setDefaultBondOrders();
print("kek.Ks ", kek.Ks)
#kek.Ks[0]=0;   # Katom
#kek.Ks[1]=0;   # Kbond
#kek.Ks[2]=0;   # KatomInt
#kek.Ks[3]=0;   # KbondInt
#Eout, xs = kek.testBond( n=40, d=0.1, ib=0, val0=0.0 )
#plt.plot( xs,Eout ); plt.show()
T0 = time.time_ns()
#kek.relax(maxIter=100, dt=0.1, ialg=0)
kek.relax(maxIter=100, dt=0.5,  ialg=1)
#kek.relaxComb()
#print("#================= Relaxed in %g[ms]" %((time.time_ns()-T0)*1e-6) )
#print( "b2a", b2a )
#print( "kek.bond2atom", kek.bond2atom )
#print( kek.Es )
print( "kek.atomValenceMin", kek.atomValenceMin )
print( "kek.atomValenceMax", kek.atomValenceMax )
print( "kek.bondOrderMin", kek.bondOrderMin )
print( "kek.bondOrderMax", kek.bondOrderMax )
print( "kek.atomValence", kek.atomValence )
print( "kek.bondOrder", kek.bondOrder )
'''


# np.angs=np.arange(kek.natoms)*2.0*np.pi/kek.natoms
# apos = np.zeros( (kek.natoms,2) )
# apos[:,0] = np.cos(np.angs)
# apos[:,1] = np.sin(np.angs)


removes = [  [0,7,8], [6,13,12], [9,10,11]   ]

#= au.addBond( base, [], bNew=True )

#(As,Bs),old_i = au.removeGroup( (AOs,b2a), removes[0] )    ;As[0]-=1
#(As,Bs),old_i = au.removeGroup( (AOs,b2a), removes[0] )    ;As[2]-=1
(As,Bs),old_i = au.removeGroup( (AOs,b2a), removes[0] )    ;As[4]-=1
#(As,Bs),old_i = au.disolveAtom( (AOs,b2a), 8 )            ;As[1] -=1 
apos_ = apos[old_i,:]                            
#pu.plotAtoms( apos_, bNumbers=True )
pu.plotAtoms( apos_, labels=np.array(As), marker='.' )
#pu.plotBonds( links=Bs, ps=apos_ )

kek.runSystem( As, Bs )
pu.plotBonds( links=Bs, ps=apos_, lws=(kek.bondOrder-0.8)*5 )

plt.figure(); pu.plotAtoms( apos_, bNumbers=True )



'''
(As,Bs),old_i = au.removeGroup( (AOs,b2a), removes[0]+removes[1] )
#(As,Bs),old_i = au.removeGroup( (AOs,b2a), removes[0] )
print( "old_i ", old_i)
print( "Bs ", Bs )
apos_ = apos[old_i,:]

#pu.plotAtoms( apos_, bNumbers=True )
pu.plotAtoms( apos_, labels=As )
pu.plotBonds( links=Bs, ps=apos_ )
'''

#plt.figure()
#pu.plotAtoms( apos, bNumbers=True )
#pu.plotAtoms( apos, labels=AOs )
#pu.plotBonds( links=b2a, ps=apos )

plt.axis('equal')
plt.show(  )


