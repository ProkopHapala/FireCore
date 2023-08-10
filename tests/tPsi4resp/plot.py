#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt

Hartree2kcal = 627.503
rad2deg      = 180.0/np.pi 


fname='Emap.dat'
if len(sys.argv)>1: fname=sys.argv[1]

data = np.genfromtxt(fname, comments='?')
#data = np.loadtxt('Emap.dat', comments='?')
#print(data)
print(data.shape)
rs=data[:,2]    #;print("rs   ", rs  )
angs=data[:,4]  #;print("angs ", angs)
Es=data[:,6]    #;print("Es   ", Es  )

Es*=Hartree2kcal

cmap='plasma'

vmin=Es.min()
Es = np.reshape( Es, (-1,6) )

r0=2.0
extent=[rs.min()+r0,rs.max()+r0, angs.min()*rad2deg,angs.max()*rad2deg ]

print(extent)

#plt.imshow( Es.T, cmap, vmax=vmin+2.0, vmin=vmin, origin='lower', extent=extent )
plt.imshow( Es.T, cmap, vmax=0, vmin=vmin, origin='lower', extent=extent )
plt.colorbar()
plt.axis('tight')
plt.xlabel('Distance [A]')
plt.ylabel('Angle [deg.]')
plt.savefig( 'Emap.png', bbox_inches='tight' )
plt.show()