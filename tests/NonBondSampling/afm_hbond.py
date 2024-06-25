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

# ================== Body

mol.findBonds()                                             #; print( "mol.bonds ", mol.bonds )
#bonds        = au.filterBonds( mol.bonds, mol.enames, set('H') )   #; print( "bonds2 ", bonds )
bonds        = mol.bonds
binds        = np.repeat(np.arange(len(bonds)), 2)   ;print("binds ", binds )  # indexes of bonds for each point
#bsamp       = au.makeBondSamples( bonds, mol.apos, where=[-0.4,0.0,0.4] )
bsamp        = au.makeBondSamples( bonds, mol.apos, where=None )
centers, cis = au.colapse_to_means( bsamp, R=0.7, binds=binds )

#tris = skelet2triangles( bonds, b2c, npoint=len(mol.apos)+1 );   print("tris ", tris)
tris = skelet2triangles2( cis, binds, bonds, npoint=len(mol.apos) );   print("tris ", tris)

print( "mol.apos.shape, centers.shape ", mol.apos.shape, centers.shape  )
points = np.concatenate( (mol.apos[:,:2], centers[:,:2]), axis=0 )

plt.show()

#exit()
neighs = au.neigh_atoms( len(mol.enames), mol.bonds )

asamp = au.makeAtomSamples    ( neighs, mol.apos, mol.enames,  ignore=set(['H']), where=[-0.6] )
esamp = au.makeEndAtomSamples ( neighs, mol.apos, mol.enames, ignore=set(['H']),  whereX=[-0.6, +0.6 ], whereY=[0.0,+0.6] )
ksamp = au.makeKinkAtomSamples( neighs, mol.apos, where=[-0.6, +0.6 ] )

nbsamp = len(bsamp)
nasamp = len(asamp)
nesamp = len(esamp)
nksamp = len(ksamp)
                                       
print( " tot npoint= ", nbsamp+nasamp+nesamp+nksamp ," nbsamp=", nbsamp," nasamp=", nasamp," nasamp=", nesamp," nksamp=", nksamp )


plt.plot( mol.apos[:,0],mol.apos[:,1], "o", ms=10, )
plt.plot( centers[:,0],centers[:,1], "o", ms=15. )
#plt.plot( bsamp[:,0],bsamp[:,1], "." )

#plt.plot( esamp[:,0],esamp[:,1], "." )
#plt.plot( ksamp[:,0],ksamp[:,1], "." )
#plt.plot( asamp[:,0],asamp[:,1], "." )



drawTriangles(tris, points, ax=None)

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