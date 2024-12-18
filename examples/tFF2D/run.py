import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

sys.path.append("../../")
from pyBall import FF2D as ff
from pyBall import atomicUtils as au
from pyBall import plotUtils   as pu

s='''
6 0 2 3 4;
6 0 1;
6 0 1;
6 0 1; 
'''

s='''
6 1 6 2 7; 
6 1 1 3; 
6 1 2 4 8; 
6 1 3 5;
6 1 4 6 9; 
6 1 5 1;
6 0 1;
6 0 3;
6 0 5;
'''


n = ff.init(s )
apos = np.zeros( (n,2) )

'''
for i in range(5):
    ff.step()
    ff.toArrays( apos=apos )
    print(i,apos)
'''

ff.run()

ff.toArrays( apos=apos )
plt.plot( apos[:,0], apos[:,1], 'o'); plt.axis('equal'); plt.grid()

plt.show()