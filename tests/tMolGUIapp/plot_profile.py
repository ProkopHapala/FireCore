#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

dat = np.genfromtxt('gridFF_EFprofile.log', skip_header=1 )
plt.plot( dat[:,3], dat[:,7], label="E" )
plt.plot( dat[:,3], dat[:,6], label="F_z" )
plt.legend()
plt.grid()
plt.show()