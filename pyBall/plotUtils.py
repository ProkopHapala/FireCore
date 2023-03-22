import numpy as np
import matplotlib.pyplot as plt
from   matplotlib import collections  as mc
from . import elements

def plotAtoms( apos=None, es=None, atoms=None, bNumbers=False, labels=None, sizes=100., colors='#808080', marker='o', axes=(0,1) ):
    ax1,ax2=axes
    if apos is None: apos = np.array([ a[1] for a in atoms ])  #;print(apos)
    plt.scatter( apos[:,ax1],apos[:,ax2], marker=marker, c=colors, s=sizes, cmap='seismic', zorder=2 ); plt.axis('equal'); #plt.grid()
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


def plotSystem( sys , bBonds=True, colors=None, sizes=None, extent=None, sz=50., RvdwCut=0.5, axes=(0,1), bLabels=True, labels=None, _0=1 ):    
    if( bBonds ):
        if sys.bonds is None:
            sys.findBonds( Rcut=3.0, RvdwCut=RvdwCut )
        plotBonds( links=sys.bonds, ps=sys.apos, axes=axes )
    
    if(colors is None): colors = [ elements.ELEMENT_DICT[e][8]    for e in sys.enames ]
    if(sizes  is None): sizes  = [ elements.ELEMENT_DICT[e][6]*sz for e in sys.enames ]
    if((labels is None) and bLabels): labels=[ "%s%i" %(e,i+_0) for i,e in enumerate(sys.enames) ]
    plotAtoms( apos=sys.apos, es=sys.enames, sizes=sizes, colors=colors, marker='o', axes=axes, labels=labels )

    if extent is not None:
        plt.xlim(extent[0],extent[1])
        plt.ylim(extent[2],extent[3]) 


def plotTrj( trj, bBonds=True, sz=50., numbers=None, axes=(0,1), extent=None, prefix="mol_", RvdwCut=0.5, figsize=(5,5) ):
    if numbers is None: numbers=range(len(trj))
    for i,sys in enumerate(trj):
        print("plot # ", i)
        fig = plt.figure(figsize=figsize)
        
        if( bBonds ):
            if sys.bonds is None:
                sys.findBonds( Rcut=3.0, RvdwCut=RvdwCut )
                plotBonds( links=sys.bonds, ps=sys.apos, axes=(0,1) )
        
        colors = [ elements.ELEMENT_DICT[e][8]    for e in sys.enames ]
        sizes  = [ elements.ELEMENT_DICT[e][6]*sz for e in sys.enames ]
        plotAtoms( apos=sys.apos, es=sys.enames, sizes=sizes, colors=colors, marker='o', axes=axes )

        if extent is not None:
            plt.xlim(extent[0],extent[1])
            plt.ylim(extent[2],extent[3])
        
        plt.savefig( prefix+("%03i.png" %numbers[i]), bbox_inches='tight' )
        plt.close(fig)