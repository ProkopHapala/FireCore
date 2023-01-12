import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import Kekule as kek

#======== Body

b2a = np.array([(0,1),(1,2),(2,3),(3,4),(4,5),(5,0)],dtype=np.int32)
AOs = np.array([3,3,3,3,3,3],dtype=np.float64)

#b2a = np.array([(0,1)],dtype=np.int32)
#AOs = np.array([2,2],dtype=np.float64)

kek.setVerbosity(2,1)
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

#kek.relax(maxIter=100, dt=0.1, ialg=0)  ;print("#================= Relaxed ")
#kek.relax(maxIter=100, dt=0.5,  ialg=1)  ;print("#================= Relaxed ")
kek.relaxComb()  ;print("#================= Relaxed ")

#print( "b2a", b2a )
#print( "kek.bond2atom", kek.bond2atom )
#print( kek.Es )
print( "kek.atomValenceMin", kek.atomValenceMin )
print( "kek.atomValenceMax", kek.atomValenceMax )
print( "kek.bondOrderMin", kek.bondOrderMin )
print( "kek.bondOrderMax", kek.bondOrderMax )
print( "kek.atomValence", kek.atomValence )
print( "kek.bondOrder", kek.bondOrder )



np.angs=np.arange(kek.natoms)*2.0*np.pi/kek.natoms
apos = np.zeros( (kek.natoms,2) )
apos[:,0] = np.cos(np.angs)
apos[:,1] = np.sin(np.angs)


from matplotlib import collections  as mc
lines=np.zeros( (kek.nbonds,2,2) )
lines[:,0,:] = apos[ b2a[:,0],: ]
lines[:,1,:] = apos[ b2a[:,1],: ]
lc = mc.LineCollection(lines, linewidths=(kek.bondOrder-0.8)*5 )
fig, ax = plt.subplots()
ax.add_collection(lc)
#ax.autoscale()
#ax.margins(0.1)


plt.plot(apos[:,0],apos[:,1],"o")
plt.axis('equal')
plt.show(  )


