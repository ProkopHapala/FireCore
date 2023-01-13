import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

sys.path.append("../../")
from pyBall import Kekule as kek
from pyBall import atomicUtils as au
from pyBall import plotUtils   as pu


groups={
'-CH3':1,'=CH2':2,'-CH2-':2,'=CH-':3,'C':4,
'-NH2':1,'=NH' :2,'-NH-' :2,'=N-' :3,
'-OH' :1,'=O'  :2,'-O-'  :2,
}

donors    =[ '-OH','-NH2','-NH-' ]
acceptors =[ '=O' ,'-O-' ,'=N-'  ]

'''
      Donor         Acceptor
-----------------------------
1    -OH,-NH2             
2        -NH-         -O-,=O
3                     =N-
'''


#======== Body

# ========== Ethylene
#b2a = np.array([(0,1)],dtype=np.int32)
#AOs = np.array([2,2],dtype=np.float64)

# ========== Benzene
#b2a = np.array([(0,1),(1,2),(2,3),(3,4),(4,5),(5,0)],dtype=np.int32)
#AOs = np.array([3,3,3,3,3,3],dtype=np.float64)

# ========== Antracene
b2a = np.array([
    (0,1),(1,2),(2,3),(3,4),(4,5),(5,6),
    (7,8),(8,9),(9,10),(10,11),(11,12),(12,13),
    (0,7),      (2,9),         (4,11),   (6,13),     
],dtype=np.int32)
AOs = np.array( [3,3,4,3,4,3,3, 3,3,4,3,4,3,3],dtype=np.float64)
apos=np.array([
    (-3.,-0.5),(-2.,-1.),(-1.,-0.5),(0.,-1.),(+1.,-0.5),(+2.,-1.),(+3.,-0.5),
    (-3.,+0.5),(-2.,+1.),(-1.,+0.5),(0.,+1.),(+1.,+0.5),(+2.,+1.),(+3.,+0.5),
])



'''
center =([2,3,3,3,2], [(0,1),(1,2),(2,3),(3,4),(4,5)]  )
bridge2=(    [2,3,2], [(0,1),(1,2),(2,3)]              )
bridge1=(      [2,2], [(0,1),(1,2)]                    ) 
'''




kek.setVerbosity(2,0)

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


