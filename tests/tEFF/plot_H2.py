from ctypes import wstring_at
import sys
import numpy as np
import os
import time
import argparse
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import eFF as eff



#evalfile = "Emolecules/CH3/.mol2"


if __name__ == "__main__":
    distance = -0.120350
    _start = -0.120325
    iterations = 200
    if os.path.exists("processXYZ.xyz"):
        os.remove("processXYZ.xyz") # deleting useless information
    H2 = np.split(eff.EvalH2(_start, iterations, distance), 1)
    #print("H1", H2[1])
    x = np.linspace(_start, distance, iterations)
    plt.plot( x, H2[0], label="Distance" )
    #plt.plot( x, H2[1], label="e1" )
    #plt.plot( x, H2[2], label = "e2" )
    plt.xlabel("Aparam")
    plt.ylabel("Force")
    plt.axhline(y=0, color='r', linestyle='--')
    plt.grid() 
    plt.legend()
    plt.show()