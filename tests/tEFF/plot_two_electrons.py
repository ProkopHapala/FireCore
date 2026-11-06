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
    distance = 1.8
    iterations = 1000
    evToAu = 0.0367502
    differentspin = np.split(eff.EvalTwoElectrons( iterations, distance, False, False, True )*evToAu, 2)
    samespin = np.split(eff.EvalTwoElectrons( iterations, distance, True, False, True )*evToAu, 2)
    coulomb = np.split(eff.EvalTwoElectrons( iterations, distance, True, True, False)*evToAu, 2)
    x = np.linspace(0, distance, iterations)/0.529
    fig, axs = plt.subplots(1,2)
    axs[0].plot( x, differentspin[0], label="Different Spin" )
    axs[0].plot( x, samespin[0], label="Same Spin" )
    axs[0].plot( x, coulomb[0], label="Coulomb" )
    axs[0].set_xlabel("Distance (Bohr)")
    axs[0].set_ylabel("Energy (au)")
    axs[0].grid() 
    axs[0].legend()
    axs[1].plot( x, differentspin[1], label="Different Spin" )
    axs[1].plot( x, samespin[1], label="Same Spin" )
    axs[1].plot( x, coulomb[1], label="Coulomb" )
    axs[1].set_xlabel("Distance (Bohr)")
    axs[1].set_ylabel("Energy (au)")
    axs[1].grid() 
    axs[1].legend()
    plt.show()