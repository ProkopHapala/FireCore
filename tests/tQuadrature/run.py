import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

sys.path.append("../../")
from pyBall import quadrature as quad


# ==================================== Functions

def plot_points( ps, ws ):
    # equal axis in 3d plot
    ax = plt.figure( figsize=(5,5) ).add_subplot(projection='3d', proj_type = 'ortho' )
    #plt.axes('equal')
    #ax.scatter( ps[:,0], ps[:,1], zs=ps[:,2], c=ws, zdir='z', )
    ax.scatter( ps[:,0], ps[:,1], zs=ps[:,2], s=ws*300, zdir='z', )
    ax.plot( [0,1, 0,0,0, 0,0,0,   1,], [0,0, 1,0,1, 0,0,0,  0,], zs=[0,0, 0,0,0, 1,0,1,  0], zdir='z', color='k' )
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    ax.set_zlim(0,1)

# ==================================== Main



'''
#quad.setParamsSize(   np.array([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15] )   )

#quad.distributPointsTetrahedron( 5, bRef=True, bAlloc=True, bOpen=False )

quad.distributPointsTetrahedron( 4, bRef=True, bAlloc=True, imode=2 )
quad.getBuffs()

print( "npqs_ref", quad.npqs_ref )
print( "pqs_ref\n", quad.qps_ref  )

plot_points( quad.qps_ref, quad.qws_ref )
'''

dat = np.genfromtxt( "quadrature_rules/tet/6-24.txt" )
ps = dat[:,:3]*0.5 + 0.5
ws = dat[:,3]

plot_points( ps, ws )
#print( dat )



plt.show()