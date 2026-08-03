from ctypes import wstring_at
import sys
import numpy as np
import os
import time
import argparse
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import eFF as eff






if __name__ == "__main__":
    distance = -0.3
    iterations = 100
    if os.path.exists("processXYZ.xyz"):
        os.remove("processXYZ.xyz") # deleting useless information
    H2 = np.split(eff.EvalH2( iterations, distance), 3)
    print("works?")
    #print("H1", H2[1])
    x = np.linspace(0, distance, iterations)
    plt.plot( x, H2[0], label="F" )
    #plt.plot( x, H2[1], label="e1" )
    #plt.plot( x, H2[2], label = "e2" )
    plt.xlabel("Aparam")
    plt.ylabel("Energy")
    plt.grid() 
    plt.legend()
    plt.show()