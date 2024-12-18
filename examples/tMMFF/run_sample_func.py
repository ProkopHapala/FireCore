import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF as mmff

#======== Body

#mmff.setVerbosity( verbosity=1, idebug=0 )

xs = np.linspace( -2., 2., 100 );

y_sin = mmff.sample_func( xs*np.pi, kind=0); plt.plot( xs, y_sin, label="y_sin" )
y_px2 = mmff.sample_func( xs,       kind=1); plt.plot( xs, y_px2, label="y_px2" )

plt.legend()
plt.grid()
plt.show()

print("\n\n#============= ALL DONE =============\n\n")