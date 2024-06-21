import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append('../../')
from pyBall import atomicUtils as au
from pyBall import plotUtils   as plu



mol = au.AtomicSystem( "../../cpp/common_resources/PTCDA-.xyz" )  #; mol.print()

#mol = au.AtomicSystem( "../../cpp/common_resources/formic_dimer.xyz" )  #; mol.print()


mol.findBonds()                                             #; print( "mol.bonds ", mol.bonds )
bonds = au.filterBonds( mol.bonds, mol.enames, set('H') )   #; print( "bonds2 ", bonds )
#bsamp       = au.makeBondSamples( bonds, mol.apos, where=[-0.4,0.0,0.4] )
bsamp        = au.makeBondSamples( bonds, mol.apos, where=None )
centers, cis = au.colapse_to_means( bsamp, R=0.7 )

#plt.show()

#exit()
neighs = au.neigh_atoms( len(mol.enames), mol.bonds )

asamp = au.makeAtomSamples( neighs, mol.apos,mol.enames, ignore=set(['H']), where=[-0.6] )
esamp = au.makeEndAtomSamples( neighs, mol.apos, mol.enames, ignore=set(['H']),  whereX=[-0.6, +0.6 ], whereY=[0.0,+0.6] )
ksamp = au.makeKinkAtomSamples( neighs, mol.apos, where=[-0.6, +0.6 ] )

nbsamp = len(bsamp)
nasamp = len(asamp)
nesamp = len(esamp)
nksamp = len(ksamp)
                                       
print( " tot npoint= ", nbsamp+nasamp+nesamp+nksamp ," nbsamp=", nbsamp," nasamp=", nasamp," nasamp=", nesamp," nksamp=", nksamp )


plt.plot( mol.apos[:,0],mol.apos[:,1], "o", ms=10, )
plt.plot( centers[:,0],centers[:,1], "o", ms=15. )
plt.plot( bsamp[:,0],bsamp[:,1], "." )

#plt.plot( esamp[:,0],esamp[:,1], "." )
#plt.plot( ksamp[:,0],ksamp[:,1], "." )
#plt.plot( asamp[:,0],asamp[:,1], "." )





#print( "neighs ", neighs )
#cycles = au.find_cycles( neighs, max_length=7)
#print( cycles )



# # Example usage:
# graph = [
#     {2, 13, 7}, {9, 3, 12}, {0, 8, 4}, {1, 20, 6}, {2, 5, 14}, {17, 4, 6},
#     {10, 3, 5}, {0, 16, 15}, {19, 2, 10}, {24, 1, 27}, {8, 6, 22}, {24, 26, 20},
#     {1, 17, 31}, {0, 18, 21}, {4, 21, 30}, {34, 19, 7}, {25, 29, 7}, {32, 12, 5},
#     {25, 28, 13}, {8, 33, 15}, {11, 3, 23}, {35, 13, 14}, {10, 36, 23}, {20, 37, 22},
#     {9, 11}, {16, 18}, {11}, {9}, {18}, {16}, {14}, {12}, {17}, {19}, {15}, {21}, {22}, {23}
# ]

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

plt.axis('equal')
plt.grid()

plt.show()