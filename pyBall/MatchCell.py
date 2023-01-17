#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

# ==================================



'''

y = sqrt( R**2 - x**2 )      #    Circle
k = ay/ax
y = k*(x-x0)

k*(x-x0)                               = sqrt( R**2 - x**2 ) 
(k**2)*( x**2  + (- 2*x0)*x + x0**2 )  =   R**2 - x**2

(k**2 +1)*(x**2) + (k**2)*(-2*x0)*x +    (k**2)*(x0**2) - R**2

A = (k**2 +1)
B = (k**2)*(-2*x0)
C = (k**2)*(x0**2) - R**2



deter =0
0 =      B**2           - 4 *       A           C         
0 = ((k**2)*(-2*x0))**2 - 4 * (k**2 +1)*( (k**2)*(x0**2) - R**2 )
solve for x0

x0**2 ( 4* (k**2)**2   - 4*(k**2 +1)*(k**2 )  = 4*(k**2 +1) * R**2
x0**2 = 4*(k**2 +1)(R**2) / ( -4*k**2 )
x0 = +/-   ( R/(2*k) )

'''

def determinant(A,B,C):
    return B*B - 4*A*C

def solveQuadratic(A,B,C):
    sD = np.sqrt( determinant(A,B,C) )
    denom = 1/(2*A)
    x1 = (-B - sD)*denom
    x2 = (-B + sD)*denom
    return x1,x2 


def find_y(  ax, ay, R, x0 ):
    k = ay/ax
    k2=k*k
    A = (k2 +1)
    B = (k2)*(-2*x0)
    C = (k2)*(x0**2) - R**2
    x1,x2 = solveQuadratic(A,B,C)
    y1 = (x1-x0)*k
    y2 = (x2-x0)*k
    return (x1,y1), (x2,y2)

R = 13.0
angs = np.linspace(0.0,2*np.pi,100);

A = np.array((0.5, 1.0))
B = np.array((2.5, 0.0))

nB = B*3.0


'''
for i in range(-5,5):
    p1,p2 = find_y(  A[0], A[1], R, B[0]*i )   ;print(i,p1,p2)
    plt.plot( [p1[0], p2[0]], [p1[1], p2[1]], '.-' )
'''

x0s = np.arange(0, (int(R/B[0])*2 +2) ) * B[0]
p1,p2 = find_y(  A[0], A[1], R, x0s )   ;print("p1s",p1,"\n p2s ",p2)

plt.plot( R*np.cos(angs), R*np.sin(angs), )
plt.plot( p1[0], p1[1], 'o' )
plt.plot( p2[0], p2[1], 'o' )
plt.axis('equal'); plt.grid(); plt.show()



'''



plt.plot( [0., B[0]], [0., B[1]], '.-' )
plt.plot( [0.,nB[0]], [0.,nB[1]], '.-')
plt.plot( [nB[0],nB[0]+A[0]], [nB[1],nB[1]+A[1]], '.-')
plt.plot( R*np.cos(angs), R*np.sin(angs), )

plt.axis('equal'); plt.grid()

plt.show()
'''