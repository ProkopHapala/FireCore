import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

sys.path.append('../../')
from pyBall import atomicUtils as au
from pyBall import plotUtils   as plu

# ================== Functions

def skelet2triangles( bonds, bond2center, npoint=None ):
    #print( "bonds ", bonds )
    if npoint is None: npoint =  np.array(bonds,dtype=np.int32).max() + 1 
    tris = []
    nb = len(bonds)
    for ib,b in enumerate(bonds):
        tris.append( [b[0],b[1],npoint+bond2center[ib]] )
    return tris

def skelet2triangles2( c2p, binds, bonds, npoint=None ):
    #print( "bonds ", bonds )
    if npoint is None: npoint =  np.array(bonds,dtype=np.int32).max() + 1 
    tris = []
    nb = len(bonds)
    for ic,ci in enumerate(c2p):
        for ip in ci:
            ib = binds[ip]
            b = bonds[ib]
            tris.append( [b[0],b[1],npoint+ic ] )
    return tris

def drawTriangles(tris, points, ax=None, ec='gray', fc=None  ):
    if ax is None: ax = plt.gca()
    triangles       = np.array( [points[tri,:2] for tri in tris] )
    poly_collection = PolyCollection(triangles, edgecolors=ec, facecolors=fc )
    ax.add_collection(poly_collection)

def bonder_edges(polygons, bonds, border_polys=None, border_atoms=None, neighs=None, np0=0 ):
    if border_atoms is None:
        is_capping   = [ (len(ngs)==1) for ngs in neighs ]
        is_kink      = [ (len(ngs)==2) for ngs in neighs ]
        border_atoms = np.logical_or( is_capping, is_kink )
    if border_polys is None:
        is_capping_b = [ (is_capping[b[0]] or is_capping[b[1]]) for b in bonds ]
        border_polys = [ any( is_capping_b[ binds[edge] ] for edge in polygon) for polygon in polygons ]
    edges = []
    for ip,poly in enumerate(polygons):
        if not border_polys[ip]: continue
        for ib in poly:
            (i,j) = bonds[ib]
            if   border_atoms[i]:
                edges.append( (i,ip+np0) )
            elif border_atoms[j]:
                edges.append( (j,ip+np0) )
    return edges

def order_contour(bonds, start=None ):
    # Create a dictionary to store adjacency list
    adjacency = {}
    for i, j in bonds:
        if i not in adjacency: adjacency[i] = []
        if j not in adjacency: adjacency[j] = []
        adjacency[i].append(j)
        adjacency[j].append(i)
    # Start from an arbitrary point
    bClosed = True
    if start is None: 
        endpoints = [point for point, neighbors in adjacency.items() if len(neighbors) == 1]
        if len(endpoints) == 0:
            start = bonds[0][0]
        else:
            bClosed = False
            start = endpoints[0]
    if bClosed:
        contour = []
    else:
        contour = [start]
    current = start
    prev    = None 
    while True:
        next_point = None
        ngs = adjacency[current]
        bFound = False
        for ng in ngs:
            if ng != prev:
                prev       = current
                current    = ng
                bFound = True
                break
        if not bFound: break 
        contour.append(current)
        if current == start: break
    return contour, bClosed


def polygonContours( polys, bonds ):
    controus    = []
    cont_closed = []
    bonds = np.array( bonds, dtype=np.int32 )
    #print( ">>>>> bonds", bonds )
    #print( ">>>>> polys", polys )
    for poly in polys: 
        ps, closed = order_contour( bonds[ list(poly) ] )
        controus.append( ps )
        cont_closed.append( closed )
    return controus, cont_closed

def offset_contour( controur, ps, l=2.0, cAve=0.5):
    print("========")
    n = len(controur)
    p0 = np.average( ps, axis=0 )
    ps2 = np.zeros( (n,2) ) 
    pout = None
    for i in range(0, n ):
        pp  = ps[ controur[i-2] ]
        pi  = ps[ controur[i-1] ]
        pm  = ps[ controur[i  ] ]
        d1  = pp-pi;   d1/=np.sqrt( np.dot(d1,d1) )
        d2  = pm-pi;   d2/=np.sqrt( np.dot(d2,d2) )
        d   = d1 + d2; d /=np.sqrt( np.dot(d ,d) )
        d0  = pi - p0; c=np.dot(d0,d)
        if c<0: d*=-1
        #print( i, pi )
        ps2[i-1,:] = pi + d*l
    if cAve > 0:
        ps2 = (ps2[:,:] + cAve*np.roll(ps2,1,axis=0) + cAve*np.roll(ps2,-1,axis=0) )/(1. + 2.*cAve)
    return ps2

def contours2bonds( contours, bonds ):
    bdict = { b:i for i,b in enumerate(bonds) }
    print( "bdict ", bdict )
    bss = []
    for c in contours:
        bs = []
        for ic in range( len(c) ):
            i,j = c[ic],c[ic-1]
            if i<j:
                b=(i,j)
            else:
                b=(j,i)
            ib = bdict.get(b,-1)
            if(ib>=0):
                bs.append( ib )
        bss.append(bs)
    return bss

def controusPoints( contours, ps, centers, beta=0.2 ):
    ps2 = [] 
    for ic,ips in enumerate(contours):
        ps2.append( (1-beta)*ps[ips,:] +  beta*centers[ic][None,:] )
        #plt.plot( (1-beta)*ps[ips,0] + beta*centers[ic][None,0], (1-beta)*ps[ips,1] + beta*centers[ic][None,1], 'o:', c='k' )
    return ps2

def small_trinagles( conts, cont_closed, bss, cont0, bss0 ):
    tris = []
    nc = 0
    for i,bs in enumerate(bss):
        c      = conts[i]
        closed = cont_closed[i]
        nci = len(c)
        if closed:
            tris.append(  [ bss0, nc+cont0, nci-1+nc+cont0 ] )
            pass
        for j,ib in enumerate(bs):
            j1 = j
            if closed: 
                j1-=1
            else:
                j1+=1
            if( j1>=0 ): # j1 = nci
                tris.append(  [ ib+bss0, j+nc+cont0, j1+nc+cont0 ] )
                pass
        nc +=  nci
    return tris
