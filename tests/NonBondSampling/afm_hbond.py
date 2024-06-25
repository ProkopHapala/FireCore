import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection


sys.path.append('../../')
from pyBall import atomicUtils as au
from pyBall import plotUtils   as plu


# ================== Setup

mol = au.AtomicSystem( "../../cpp/common_resources/PTCDA-.xyz" )  #; mol.print()

#mol = au.AtomicSystem( "../../cpp/common_resources/formic_dimer.xyz" )  #; mol.print()

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

def drawTriangles(tris, points, ax=None):
    if ax is None: ax = plt.gca()
    # Create a list of triangle vertices from point indexes
    #triangles = []
    #for tri in tris:
    #    triangles.append([points[tri[0]], points[tri[1]], points[tri[2]]])
    triangles = np.array( [points[tri,:2] for tri in tris] )
    #print( "triangles ", triangles )
    #print( "triangles.shape ", triangles.shape )
    #poly_collection = PolyCollection(triangles, edgecolors='k', facecolors='none')
    poly_collection = PolyCollection(triangles, edgecolors='gray' )
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
    if start is None: start = bonds[0][0]
    #contour = [start]
    contour = []
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
    return contour

def offset_contour( controur, ps, l=2.0, cAve=1.0 ):
    print("========")
    n = len(controur)
    p0 = np.average( ps, axis=0 )
    ps2 = np.zeros( (n,2) ) 
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




# ================== Body

mol.findBonds()                                             #; print( "mol.bonds ", mol.bonds )
#bonds        = au.filterBonds( mol.bonds, mol.enames, set('H') )   #; print( "bonds2 ", bonds )
bonds        = mol.bonds
binds        = np.repeat(np.arange(len(bonds)), 2)   ;print("binds ", binds )  # indexes of bonds for each point
#bsamp       = au.makeBondSamples( bonds, mol.apos, where=[-0.4,0.0,0.4] )
bsamp        = au.makeBondSamples( bonds, mol.apos, where=None )
centers, polygons = au.colapse_to_means( bsamp, R=0.7, binds=binds )

polys = [  set( binds[p] for p in poly )  for poly in polygons ]

print( "polys ", polys )

#tris = skelet2triangles( bonds, b2c, npoint=len(mol.apos)+1 );   print("tris ", tris)
tris = skelet2triangles2( polygons, binds, bonds, npoint=len(mol.apos) ); #  print("tris ", tris)

#print( "mol.apos.shape, centers.shape ", mol.apos.shape, centers.shape  )
points = np.concatenate( (mol.apos[:,:2], centers[:,:2]), axis=0 )

plt.show()

#exit()
neighs = au.neigh_atoms( len(mol.enames), mol.bonds )

is_capping   = [ (len(neighbors)==1) for neighbors in neighs ]
is_kink      = [ (len(neighbors)==2) for neighbors in neighs ]
is_capping_b = [ (is_capping[b[0]] or is_capping[b[1]]) for b in bonds ]
border_atoms = np.logical_or( is_capping, is_kink )
#is_peripheral = [ ((len(neighbors)==1) or (len(neighbors)==2))  for neighbors in neighs ]
#is_per = [  for ci in  cis ]
border_polys = [ any( is_capping_b[ binds[edge] ] for edge in polygon) for polygon in polygons ]


edges = bonder_edges( polys, bonds, border_polys=border_polys, border_atoms=border_atoms, np0=len(mol.apos),  );  print( "edges ", edges )
#edges = bonder_edges( polys, bonds, neighs=neighs, np0=len(mol.apos) );  print( "edges ", edges )

contour = order_contour(edges) #;print( "contour ", contour )

plt.plot( points[contour,0], points[contour,1], 'o:r', lw=3  )

contour2  = offset_contour( contour, points, l=2.0 )     #;print( "contour2 ", contour2.shape, contour2  )

plt.plot( contour2[:,0], contour2[:,1], 'o-b', lw=3  )

#plu.plotBonds( links=edges, ps=points, lws=5, colors='r' )
     



asamp = au.makeAtomSamples    ( neighs, mol.apos, mol.enames,  ignore=set(['H']), where=[-0.6] )
esamp = au.makeEndAtomSamples ( neighs, mol.apos, mol.enames, ignore=set(['H']),  whereX=[-0.6, +0.6 ], whereY=[0.0,+0.6] )
ksamp = au.makeKinkAtomSamples( neighs, mol.apos, where=[-0.6, +0.6 ] )

nbsamp = len(bsamp)
nasamp = len(asamp)
nesamp = len(esamp)
nksamp = len(ksamp)
                                       
print( " tot npoint= ", nbsamp+nasamp+nesamp+nksamp ," nbsamp=", nbsamp," nasamp=", nasamp," nasamp=", nesamp," nksamp=", nksamp )


plt.plot( mol.apos[:,0],mol.apos[:,1], "o", ms=10, )

drawTriangles(tris, points, ax=None)

plt.plot( centers[:,0],centers[:,1], "o", ms=15. )
#plt.plot( bsamp[:,0],bsamp[:,1], "." )


plt.plot( mol.apos[border_atoms,0],mol.apos[border_atoms,1], "xm", ms=15. )
plt.plot( centers[border_polys,0],centers[border_polys,1], "*m", ms=15. )

#plt.plot( esamp[:,0],esamp[:,1], "." )
#plt.plot( ksamp[:,0],ksamp[:,1], "." )
#plt.plot( asamp[:,0],asamp[:,1], "." )





ax=plt.gca()
ax.scatter(points[:, 0], points[:, 1], c='b', marker='o')    
for i, point in enumerate(points):  ax.annotate(str(i), (point[0], point[1]), textcoords="offset points", xytext=(5,5), ha='center')


plu.plotSystem( mol, bLabels=False )

#print( "neighs ", neighs )
#cycles = au.find_cycles( neighs, max_length=7)
#print( cycles )


'''
# ------------ Cycles
adj_list           = au.convert_to_adjacency_list(neighs)
preprocessed_graph = au.preprocess_graph(adj_list.copy())     ;#print( "preprocessed_graph ", preprocessed_graph  )
cycles             = au.find_cycles(preprocessed_graph, max_length=7)
print( "cycles ", cycles )


centers = np.zeros( (len(cycles),3) )
for i,cycle in enumerate(cycles):
    ps = mol.apos[cycle,:]
    centers[i] = ps.sum(axis=0)/len(ps)
    plt.plot( ps[:,0], ps[:,1] )

plt.plot( centers[:,0], centers[:,1], '*', )
'''

plt.axis('equal')
#plt.grid()

plt.show()