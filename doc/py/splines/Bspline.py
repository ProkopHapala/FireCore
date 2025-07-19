import numpy             as np
import matplotlib.pyplot as plt

'''
# Cubic B-spline 
# See:  https://www.studocu.com/en-us/document/university-of-michigan/computer-graphics/lecture-notes-lecture-37/839972
#       https://web.eecs.umich.edu/~sugih/courses/eecs487/lectures/37-B-splines+NURBS.pdf
#       https://dokumen.tips/documents/eecs-487-interactive-computer-graphics-sugihcourseseecs487-eecs-487-interactive.html?page=1


Basis:

(1/6) *  u^3
(1/6) * -3u^3 + 3u^2 + 3u + 1   = 3*u*  ( 1 + u - u2 )   +  1
(1/6) * +3u^3 - 6u^2 + 4        = 3*u^2*( u - 2      )   +  4
(1/6) * (1-u)^3

dBasis:

(1/2)   u^2
(1/2)  -3u^2 + 2u + 1 
(1/2)   3u^2  - 4u
(1/2)  -(1-u)^2

ddBasis:

a1 = u
a2 =  -3*u + 1
a3 =   3*u - 2
a4 = 1-u

Points:

  Basis: [ 0.0, 1./6.,   2./3.,   1./6. ]
 dBasis: [ 0.0, 0.5,     0.0,    -0.5   ]
ddBasis: [ 0.0, 1.0,    -2.0,     1.0   ]

'''

def Bbasis( u ):
    inv6 = 1./6.
    b1 = inv6 *  u**3
    b2 = inv6 *( -3*u**3 + 3*u**2 + 3*u + 1 )  
    b3 = inv6 *( +3*u**3 - 6*u**2 + 4       )
    b4 = inv6 * (1-u)**3
    return b1,b2,b3,b4

def dBbasis( u ):
    d1 = 0.5* u**2
    d2 = 0.5*( -3*u**2 + 2*u + 1  )
    d3 = 0.5*(  3*u**2 - 4*u     )
    d4 = 0.5* -(1-u)**2
    return d1,d2,d3,d4

def ddBbasis( u ):
    a1 = u
    a2 =  -3*u + 1
    a3 =   3*u - 2
    a4 = 1-u
    return a1,a2,a3,a4


us = np.linspace( 0.0,1.0,100 )

b1,b2,b3,b4 = Bbasis  ( us );   print( b1[0],b2[0],b3[0],b4[0] )
d1,d2,d3,d4 = dBbasis ( us );   print( d1[0],d2[0],d3[0],d4[0] )
a1,a2,a3,a4 = ddBbasis( us );   print( a1[0],a2[0],a3[0],a4[0] )

plt.figure(figsize=(5,15))
plt.subplot(3,1,1)
plt.plot( us,  b1, label="B1" )
plt.plot( us,  b2, label="B2" )
plt.plot( us,  b3, label="B3" )
plt.plot( us,  b4, label="B4" )
plt.grid()
plt.legend()

plt.subplot(3,1,2)
plt.plot( us,  d1, label="dB1" )
plt.plot( us,  d2, label="dB2" )
plt.plot( us,  d3, label="dB3" )
plt.plot( us,  d4, label="dB4" )
plt.grid()
plt.legend()

plt.subplot(3,1,3)
plt.plot( us,  a1, label="ddB1" )
plt.plot( us,  a2, label="ddB2" )
plt.plot( us,  a3, label="ddB3" )
plt.plot( us,  a4, label="ddB4" )
plt.grid()
plt.legend()


plt.show()