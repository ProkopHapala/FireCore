import matplotlib.pyplot as plt
import numpy as np

def plotAtoms( es=None, apos=None, atoms=None, ax1=0, ax2=1 ):
    if apos is None: apos = np.array([ a[1] for a in atoms ])  #;print(apos)
    plt.plot( apos[:,ax1],apos[:,ax2], 'o' ); plt.axis('equal'); plt.grid()