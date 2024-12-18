#!/usr/bin/python

import numpy      as np
import matplotlib.pyplot as plt

dat = np.genfromtxt('lattice_scan_1d_all.dat').transpose()

plt.plot( dat[2], dat[11] )

plt.show()
