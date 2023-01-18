import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

sys.path.append("../../")
from pyBall import Lattice2D as lat

#    LM.lat0[0].set(  1.3,  0.0 );
#    LM.lat0[1].set( -0.5,  0.6 );
#    LM.lat1[0].set(  0.3,  0.2 );
#    LM.lat1[1].set( -0.13, 0.9 );

lat0=np.array([
    [ 1.3, 0.0 ],
    [-0.5, 0.6 ]
])

lat1=np.array([
    [ 0.30, 0.2 ],
    [-0.13, 0.9 ]
])

print("lat0\n", lat0)
print("lat1\n", lat1)

#Rmax = 10.0
Rmax = 5.0

t0=time.time_ns()
n = lat.match(lat0, lat1, Rmax=Rmax, dRmax=0.1, dAngMax=0.1  )    #;print("n", n)
print( "Time{lat.match} [ms]",  (time.time_ns()-t0)*1.e-6  ) 

t0=time.time_ns()
inds,errs,Ks = lat.getMatches(inds=None, errs=None, bSort=True, Ks=(1.,1.,1.,0.) )
print( "Time{lat.getMatches} [ms]",  (time.time_ns()-t0)*1.e-6  ) 


E = errs[:,0]**2 + errs[:,1]**2 + errs[:,2]**2

print( inds )
#print( errs )
#print( E )

#plt.plot( lat0 )

#plt.plot( lat0[0] )
#plt.plot( lat0[1] )


ps = np.zeros( (30,40,2) )
ns = np.arange(ps.shape[1]) - 20
#print( ns )
for i in range( len(ps) ):
    ia=i-15
    ps[i,:,0] = lat0[0,0]*ia  + lat0[1,0]*ns
    ps[i,:,1] = lat0[0,1]*ia  + lat0[1,1]*ns
#print(ps)
ps = ps.reshape( (ps.shape[0]*ps.shape[1],2) )
#print(ps.shape)

plt.plot(ps[:,0],ps[:,1], 'k.', markersize=0.5 )
plt.plot( [lat0[0,0],0.0,lat0[1,0], lat0[1,0]+lat0[0,0], lat0[0,0] ],  [lat0[0,1],0.0,lat0[1,1],lat0[1,1]+lat0[0,1], lat0[0,1]],  'k' )      # plot latticle 1 
plt.plot( [lat1[0,0],0.0,lat1[1,0], lat1[1,0]+lat1[0,0], lat1[0,0] ],  [lat1[0,1],0.0,lat1[1,1],lat1[1,1]+lat1[0,1], lat1[0,1]],  'gray' )     # plot latticle 2 


# ----- Plot best found latticle matches
for i in range(5):
    u=lat0[0,:]*inds[i,0] + lat0[1,:]*inds[i,1] 
    v=lat0[0,:]*inds[i,2] + lat0[1,:]*inds[i,3]
    plt.plot( [u[0],0.,v[0]], [u[1],0.,v[1]]  )

# ----- Plot limit circle
angs= np.linspace(0,2.0*np.pi,100)
plt.plot( np.cos(angs)*Rmax, np.sin(angs)*Rmax , 'k-' )

plt.axis('equal')
plt.ylim(-Rmax,Rmax)
plt.xlim(-Rmax,Rmax)
plt.show(  )


