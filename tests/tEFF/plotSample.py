import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

sys.path.append("../../")
from pyBall import eFF as eff
from pyBall import eFF_terms as pyeff

size11 = 1.0
size12 = 2.0
size21 = size11*2
size22 = size12*2

def example_vector_function(points):
    Fout = eff.sample_ee(points, 1, bEvalCoulomb=True, bEvalPauli=False)
    return Fout

nx = 100
r = 10
X = np.linspace(-r, r, nx)

points1 = np.zeros((len(X), 3))
points1[:,0] = X
points1[:,1] = size11
points1[:,2] = size12

FEout1 = example_vector_function(points1)

points2 = np.zeros((len(X), 3))
points2[:,0] = X
points2[:,1] = size21
points2[:,2] = size22

FEout2 = example_vector_function(points2)

plt.plot(X, FEout1[:,0], label="fx1")
plt.plot(X, FEout1[:,1], label="fs11")
plt.plot(X, FEout1[:,2], label="fs21")
plt.plot(X, FEout1[:,3], label="E1")

fx = FEout1[:,3]
df_dx = (fx[2:] - fx[:-2])/(X[2] - X[0])

plt.plot(X, FEout2[:,0], label="fx2")
plt.plot(X, FEout2[:,1], label="fs12")
plt.plot(X, FEout2[:,2], label="fs22")
plt.plot(X, FEout2[:,3], label="E2")

plt.legend()
plt.title("Sample ee")
plt.show()
