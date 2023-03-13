import numpy as np
import matplotlib.pyplot as plt
from matplotlib import collections  as mc

def plotAtoms( apos=None, es=None, atoms=None, bNumbers=False, labels=None, sizes=100., colors='#808080', marker='o', axes=(0,1) ):
    ax1,ax2=axes
    if apos is None: apos = np.array([ a[1] for a in atoms ])  #;print(apos)
    plt.scatter( apos[:,ax1],apos[:,ax2], marker=marker, c=colors, s=sizes, cmap='seismic' ); plt.axis('equal'); #plt.grid()
    bLabels = labels is not None
    if bNumbers or bLabels:
        na = len(apos)
        if not bLabels:
            labels = range(na)
        ax= plt.gca()
        for i in range(na):
            ax.annotate( str(labels[i]), (apos[i,ax1], apos[i,ax2]))


def plotBonds( lps=None, links=None, ps=None, lws=None, axes=(0,1) ):
    ax1,ax2=axes
    ax_inds=[ax1,ax2]
    if lps is None:
        ps_=ps
        if ps_.shape[1]!=2:
            ps_ = ps[:,ax_inds]
        links=np.array(links,dtype=np.int32)
        lps=np.zeros( (len(links),2,2) )
        #print( "lps.shape, ps_.shape ", lps.shape, ps_.shape )
        lps[:,0,:] = ps_[ links[:,0],: ]
        lps[:,1,:] = ps_[ links[:,1],: ]
    #lws = (kek.bondOrder-0.8)*5
    lc = mc.LineCollection(lps, linewidths=lws, colors='k' )
    ax= plt.gca()
    ax.add_collection(lc)
    #ax.autoscale()
    #ax.margins(0.1)