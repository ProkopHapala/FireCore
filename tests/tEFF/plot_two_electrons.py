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
    distance = 3.0
    iterations = 1000
    differentspin = eff.EvalTwoElectrons( iterations, distance, False, False, True )
    samespin = eff.EvalTwoElectrons( iterations, distance, True, False, True )
    x = np.linspace(0, distance, iterations)
    plt.plot( x, differentspin, label="Different Spin" )
    plt.plot( x, samespin, label="Same Spin" )
    plt.xlabel("Distance")
    plt.ylabel("Energy")
    plt.grid() 
    plt.legend()
    plt.show()