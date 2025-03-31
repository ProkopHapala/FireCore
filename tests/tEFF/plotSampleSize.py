import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

sys.path.append("../../")
from pyBall import eFF as eff
from pyBall import eFF_terms as pyeff
# const_bohr = 0.5291772105638411



# Define a vector function (example: E, fx, fy, fz)
def example_vector_function(points):
    Fout = eff.sample_ee(points, 1, bEvalCoulomb=True, bEvalPauli=True)
    return Fout

nx = 100
r = 10
size = np.linspace(0, r, nx)
sizeConst = 1.0

points1 = np.zeros( (len(size), 3) )
points1[:,0] = 3 # vzdaleonst elektronu od sebe
points1[:,1] = size
points1[:,2] = sizeConst

FEout1 = example_vector_function( points1)

plt.plot(size, FEout1[:,0], label="fx1")
plt.plot(size, FEout1[:,1], label="fs11")
plt.plot(size, FEout1[:,2], label="fs21")
plt.plot(size, FEout1[:,3], label="E1")

print(points1)
print(FEout1)

# f' = (f(x+d) - f(x))/d
# f' = (f(x+d) - f(x-d))/(2*d)
# fx = FEout1[:,3]

# df_dx =  ( fx[2:] - fx[:-2] )/ ( X[2] - X[0] )

plt.legend()
plt.title("Sample ee")
plt.show()