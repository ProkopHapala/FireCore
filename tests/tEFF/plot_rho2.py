from ctypes import wstring_at
import sys
import numpy as np
import os
import time
import argparse
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import eFF as eff



elementPath_e = "export/scan_data/single_H2O_ee2_relaxed.xyz"


if __name__ == "__main__":
    distance = 3.0
    iterations = 100
    eff.setVerbosity(1,0)
    differentspin = eff.FindMinrho2( iterations, elementPath_e, nstepMax=6000000,optAlg=2, dt=0.0001, Fconv=1e-7)
    H1 = np.split(differentspin, 3)
    x = np.linspace(0, -0.3, iterations)
    print("rho2", x)
    H1[2] = H1[2] * 180 / np.pi
    print("H1", H1)
    plt.plot( x, H1[0], label="Hydrogen1" )
    plt.plot( x, H1[1], label="Hydrogen2" )
    plt.plot( x, H1[2], label = "Angle" )
    plt.axhline(y=0.96, color='r', linestyle='-')
    plt.xlabel("rho2")
    plt.ylabel("distance")
    plt.grid() 
    plt.legend()
    plt.show()
