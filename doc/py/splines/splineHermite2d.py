import numpy             as np
import matplotlib.pyplot as plt
#       x3   x2   x   1
#  ----------------------
#  y0   2   -3        1
#  y1  -2   +3
# dy0   1   -2    1
# dy1   1   -1


H=np.array([
##  x3   x2   x     1
 [  2., -3.,  0.0,  1. ],  # p0   
 [ -2., +3.,  0.0,  0. ],  # p1
 [  1., -2.,  1.0,  0. ],  # d0
 [  1., -1.,  0.0,  0. ],  # d1
])

Hp = np.zeros((4,4))

Hp[0,:] = H[2,:]*-0.5
Hp[1,:] = H[3,:]*-0.5 + H[0,:]
Hp[2,:] = H[2,:]* 0.5 + H[1,:]
Hp[3,:] = H[3,:]* 0.5
print( "Hp:\n", Hp )

dHp = Hp.copy();   dHp[:,0]*=3; dHp[:,1]*=2; dHp[:,3]*=0
print( "dHp:\n", dHp )

ddHp = dHp.copy();   ddHp[:,0]*=2; ddHp[:,2]*=0
print( "ddHp:\n", ddHp )

'''
def spline_hermite( t, f0, f1, d0, d1 ):
    print( "t.shape ",  t.shape )
    print( "f0.shape ", f0.shape )
    print( "f1.shape ", f1.shape )
    print( "d0.shape ", d0.shape )
    print( "d1.shape ", d1.shape )
    t2  = t*t
    t3  = t2*t
    h00 =  2*t3 - 3*t2 + 1
    h10 =    t3 - 2*t2 + t
    h01 = -2*t3 + 3*t2
    h11 =    t3 -   t2
    if isinstance(f0, np.ndarray):
        #print( "f0 is an array" )
        #return h00*f0 + h01*f1 + h10*d0  + h11*d1
        #return h00*f0[None,:] + h01*f1[None,:] + h10*d0[None,:]  + h11*d1[None,:]
        return h00[None,:]*f0[:,None] + h01[None,:]*f1[:,None] + h10[None,:]*d0[:,None] + h11[None,:]*d1[:,None]
    else:
        return h00*f0 + h01*f1 + h10*d0  + h11*d1


def spline_hermite_2d( p, v00, v01, v10, v11 ):
    v0y = spline_hermite( p[1], np.array(v00[:-1]), np.array(v01[:-1]), np.array([v00[2],0]), np.array([v01[2],0]) )
    print( "v0y.shape ", v0y.shape )
    #exit()
    v1y = spline_hermite( p[1], np.array(v10[:-1]), np.array(v11[:-1]), np.array([v10[2],0]), np.array([v11[2],0]) )
    return spline_hermite( p[0], v0y[0], v1y[0], v0y[1], v1y[1] )


quad = np.array([
[1.0,1.0,1.0,1.0],
[1.0,-1.0,3.0,1.0],
[1.0,1.0,2.0,1.0],
[1.0,1.0,1.0,1.0],
])

xs,ys = np.meshgrid( np.linspace(0,1,10), np.linspace(0,1,10) )

E = spline_hermite_2d( [xs.flat.copy(),ys.flat.copy()], 
    [quad[1,1], (quad[2,1]-quad[0,1])*0.5,  (quad[1,2]-quad[1,0])*0.5 ], 
    [quad[2,1], (quad[3,1]-quad[1,1])*0.5,  (quad[2,2]-quad[2,0])*0.5 ], 
    [quad[1,2], (quad[2,2]-quad[0,2])*0.5,  (quad[1,3]-quad[1,1])*0.5 ], 
    [quad[2,2], (quad[3,2]-quad[1,2])*0.5,  (quad[2,3]-quad[2,1])*0.5 ],     
 )

E=  E.reshape(10,10)

plt.imshow( E, extent=(0,1,0,1), origin='lower' )

plt.legend()
plt.show()
'''