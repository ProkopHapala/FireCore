import numpy as np
import matplotlib.pyplot as plt
from matplotlib import collections  as mc

def plotAtoms( apos=None, es=None, atoms=None, ax1=0, ax2=1, bNumbers=False, labels=None, marker='o' ):
    if apos is None: apos = np.array([ a[1] for a in atoms ])  #;print(apos)
    plt.plot( apos[:,ax1],apos[:,ax2], marker ); plt.axis('equal'); #plt.grid()
    bLabels = labels is not None
    if bNumbers or bLabels:
        na = len(apos)
        if not bLabels:
            labels = range(na)
        ax= plt.gca()
        for i in range(na):
            ax.annotate( str(labels[i]), (apos[i,0], apos[i,1]))


def plotBonds( lps=None, links=None, ps=None, lws=None ):
    if lps is None:
        links=np.array(links,dtype=np.int32)
        lps=np.zeros( (len(links),2,2) )
        lps[:,0,:] = ps[ links[:,0],: ]
        lps[:,1,:] = ps[ links[:,1],: ]
    #lws = (kek.bondOrder-0.8)*5
    lc = mc.LineCollection(lps, linewidths=lws )
    ax= plt.gca()
    ax.add_collection(lc)
    #ax.autoscale()
    #ax.margins(0.1)