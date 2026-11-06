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
    distance = 20
    iterations = 1000
    kinetic = eff.EvalHE( iterations, distance, False, True)
    both = eff.EvalHE( iterations, distance, True, True )
    coulomb = eff.EvalHE( iterations, distance, True, False)
    x = np.linspace(0.4, distance, iterations)
    plt.plot( x, kinetic, label="Kinetic" )
    plt.plot( x, both, label="Both" )
    plt.plot( x, coulomb, label="Coulomb" )
    plt.xlabel("Size")
    plt.ylabel("Energy")
    plt.grid() 
    plt.legend()
    plt.show()