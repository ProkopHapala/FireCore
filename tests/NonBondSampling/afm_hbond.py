import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

sys.path.append('../../')
from pyBall import atomicUtils as au
from pyBall import plotUtils   as plu
import molmesh2d as mm 

# ================== Setup


mol = au.AtomicSystem( "./AFM/Samples/Nucleobases/OHO-h_1-uracil.xyz")


#mol = au.AtomicSystem( "./AFM/Samples/Nucleobases/HNO-h-cytosine.xyz")
##mol.apos[:,2] *= -1
##mol.apos[[7,8,9],2] -= mol.apos[7,2] 
##mol.saveXYZ( "./AFM/Samples/Nucleobases/HNO-h-cytosine-.xyz")

#mol = au.AtomicSystem( "./AFM/Samples/Nucleobases/HN-h-p-adenine.xyz")
##mol.apos[[10,9,13],2] -= mol.apos[10,2] 
##mol.saveXYZ( "./AFM/Samples/Nucleobases/HN-h-p-adenine-.xyz")


#mol = au.AtomicSystem( "./AFM/Samples/Nucleobases/HHO-h-p_1-guanine.xyz")
##mol.apos[[10,14,15],2] -= mol.apos[10,2] 
##mol.saveXYZ( "./AFM/Samples/Nucleobases/HHO-h-p_1-guanine-.xyz")




#mol = au.AtomicSystem( "../../cpp/common_resources/formic_dimer.xyz" )  #; mol.print()
#mol = au.AtomicSystem( "../../cpp/common_resources/PTCDA-.xyz" )  #; mol.print()


#mol = au.AtomicSystem( "./AFM/Tips/H2O_H.xyz")
#mol.orientPCA(perm=(2,0,1))
#mol.saveXYZ( "./AFM/Tips/H2O_O.xyz")

#mol = au.AtomicSystem( "./AFM/Tips/NH3_N.xyz")
#mol.orient( 1, (1,0), (0,[2,3]) ) 
#mol.saveXYZ( "./AFM/Tips/NH3_H.xyz" )

tip = au.AtomicSystem( "./AFM/Tips/H2O_H.xyz" )
tip.apos[:,2] -= tip.apos[:,2].min()

# plt.figure(figsize=(5,10))
# plt.subplot(2,1,1); plu.plotSystem( mol, axes=(0,1) ); plt.xlabel("x");plt.ylabel("y"); plt.grid()
# plt.subplot(2,1,2); plu.plotSystem( mol, axes=(2,1) ); plt.xlabel("z");plt.ylabel("y"); plt.grid()
# plt.show(); plt.exit()

# ================= Functions

def printPointInfo( ps, ptypes, ngs ):
    np = len(ps)
    nt = len(ptypes)
    ng = len(ngs)
    if np!= nt or np!= ng: 
         print( "Error: printPointInfo() size mismatch ps(%i) ptypes(%i) ngs(%i)" %(np,nt,ng) )
    #     exit()
    n = len(ps)
    for i in range(n):
        print( i, ptypes[i], ps[i], ngs[i] )

def dictToString( dict ):
    lst = [ "%s %s" for k,v in dict.items() ]
    return "\n".join(lst)


params={
    "basis"    : "cc-pvdz",
    "scf_type" : "df",
}

def toPsi4( fname, mol1, mol2, params ):
    fout = open( fname, "w" )
    #fout.write( au.psi4frags2string( mol1.enames, mol1.apos ) )
    #fout.write( au.psi4frags2string( mol1.enames, mol1.apos) )
    fout.write(
    '''memory 16GB
molecule dimer {
0 1
''')
    fout.write( "".join(mol1.toLines( )) )
    fout.write( "--\n" )
    fout.write( "".join(mol2.toLines( )) )
    fout.write(
    '''units angstrom
no_reorient
symmetry c1
}
set {
''');
    fout.write(  dictToString( params ) )
    fout.write('''
}
gradient('b3lyp-d3', bsse_type='cp' )
''');
    fout.close()

def toXYZ( fout, mol1, mol2, comment="#comment" ):
    fout.write( "%i\n" % (len(mol1.enames)+len(mol2.enames)) )
    fout.write( "%s\n" %comment )
    mol1.toXYZ(fout)
    mol2.toXYZ(fout)

def make_scan_geoms( dirname, sample, tip, points, zs = [0.0,0.1,0.3,0.4,0.5], params=params ):
    n  = len(points)
    nz = len(zs)
    apos0 = tip.apos.copy()
    try:
        os.mkdir( dirname )
    except:
        pass
    fxyz = open( dirname+"/movie.xyz", "w" )
    for i in range(n):
        pi = points[i]
        for iz in range(nz):
            p = np.array( [ pi[0],pi[1],zs[iz] ] )
            tip.apos[:,:] = apos0[:,:] + p[None,:]
            path = dirname+"/"+("p%03i_z%02i" %(i,iz))
            os.mkdir( path )
            #toPsi4( "psi_%03i_%02i.in" %(i,iz), tip, sample, params )
            toXYZ( fxyz, tip, sample, comment="ip %3i iz %3i " %(i,iz) )
            toPsi4( path+"/psi4.in", tip, sample, params )
    fxyz.close()



# ================== Body
#plu.plotSystem( mol, bLabels=False )

mol.findBonds()                                             #; print( "mol.bonds ", mol.bonds )
mol.neighs()
anew, bnew = mm.makeKinkDummy( mol.apos, mol.ngs, angMin=10.0, l=1.0 )  #;print("anew ", anew, "bnew ", bnew )  ; exit()

bonds = [(i, j) if i < j else (j, i) for i, j in mol.bonds ]    # order bonds ( i < j )

bonds_bak = bonds.copy()

bonds = bonds + bnew

nao = len(mol.apos)
if len(anew) > 0:
    apos = np.concatenate( (mol.apos[:,:2], anew), axis=0 ) 
else:
    apos = mol.apos[:,:2].copy()

binds        = np.repeat(np.arange(len(bonds)), 2)   ;print("binds ", binds )  # indexes of bonds for each point
#bsamp       = au.makeBondSamples( bonds, apos, where=[-0.4,0.0,0.4] )
bsamp        = au.makeBondSamples( bonds, apos, where=None )
centers, polygons = au.colapse_to_means( bsamp, R=0.7, binds=binds )     # polygons are lists of bond indices

#print( "mol.apos.shape, centers.shape ", mol.apos.shape, centers.shape  )
points = np.concatenate( (apos, centers[:,:2]), axis=0 )     # points = atoms + centers

na  = len(apos)
nc  = len(centers)
nac = na + nc

polys = [  set( binds[p] for p in poly )  for poly in polygons ]       # for each polygon center list of bonds adjecent to it

conts, cont_closed = mm.polygonContours( polys, bonds )                      ;print( "conts 1 ", conts  )   # order contours ( ordered list of points for each polygon )
cps_,cpis = mm.controusPoints( conts, points, centers, beta=0.3 )   #;print( "cps 1 ", cps  )       # orderd list of points for each contour
cps       = np.concatenate( cps_, axis=0 )                           #;print( "cps ", cps  )

bss = mm.contours2bonds( conts, bonds )                    #;print("bss ", bss )     # bond-centers to bonds
bcs = au.makeBondSamples( bonds, apos, where=[0.0] )       #;print("bcs ", bcs )     # bond-centers

bss0   = len(points)         # start bond-centers points
cont0  = bss0 + len(bcs) 
points = np.concatenate( (points, bcs, cps ), axis=0 )        # points = atoms + bond-centers + contours

ptyps = mol.enames  + [ "kink" ]*len(anew) + [ "center" ]*len(conts) + ['bond']*len(bcs) + ['cp']*len(cps)
pngs  = [ [] ]*nao  + [ [j] for i,j in bnew ]*len(anew) + conts      + bonds               + cpis


print( "na(%i) nac(%i) bss0(%i) cont0(%i)" %( na, nac, bss0, cont0) )
print( "len: mol.apos(%i) anew(%i) centers(%i) bonds(%i) bcs(%i) cps(%i) cpis(%i) tot(%i)" %( len(mol.apos), len(anew), len(centers), len(bonds), len(bcs), len(cps), len(cpis), len(mol.apos)+len(anew)+len(centers)+len(bss)+len(cps) ) )

#printPointInfo( points, ptyps, pngs )

zs = np.array( [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2] ) + 1.6

make_scan_geoms( "guanine_H2O", mol, tip, mol.apos, zs = zs, params=params )


exit()

# ---------- PLOT points
#plu.plotSystem( mol, bLabels=False )
plu.plotSystem( mol, bLabels=False )
#plt.plot( points[:,0],          points[:,1],           'ok' );    # all points
plt.plot( points[0:na,0],       points[0:na,1],        'ok' );    # bond centers
plt.plot( points[na:nac,0],     points[na:nac,1],      'og' );    # bond centers
plt.plot( points[bss0:cont0,0], points[bss0:cont0,1],  'or' );    # bond centers
plt.plot( points[cont0:,0],     points[cont0:,1],      'ob' );    # contours centers
ax = plt.gca()
for i, point in enumerate(points):  ax.annotate(str(i), (point[0], point[1]), textcoords="offset points", xytext=(5,5), ha='center')


#for ic,centers in enumerate(centers): print( "center[%i]=p[%i]" %(ic,na+ic), polygons[ic] )
for ic,centers in enumerate(centers): print( "center[%i]=p[%i]" %(ic,na+ic), conts[ic] )


# --------- Plot Triangles

#tris2 = mm.small_trinagles( conts, cont_closed, bss, cont0, bss0 )   # ;print( "tris2 ", tris2 )    # triangles from bond-centers to atom-contours-offsets
#mm.drawTriangles(tris2, points, fc='y', ec='r' )    


#plt.show() ;exit()

#print( "polys ", polys )

#tris = skelet2triangles( bonds, b2c, npoint=len(mol.apos)+1 )   ;print("tris ", tris)
#tris = mm.skelet2triangles2( polygons, binds, bonds, npoint=len(mol.apos) )  ;print("tris ", tris)

#exit()
neighs = au.neigh_atoms( len(mol.enames), mol.bonds )

is_capping   = [ (len(neighbors)==1) for neighbors in neighs ]
is_kink      = [ (len(neighbors)==2) for neighbors in neighs ]
is_capping_b = [ (is_capping[b[0]] or is_capping[b[1]]) for b in bonds_bak ]
border_atoms = np.logical_or( is_capping, is_kink )
#is_peripheral = [ ((len(neighbors)==1) or (len(neighbors)==2))  for neighbors in neighs ]
#is_per = [  for ci in  cis ]
border_polys = [ any( is_capping_b[ binds[edge] ] for edge in polygon) for polygon in polygons ]
edges = mm.bonder_edges( polys, bonds, border_polys=border_polys, border_atoms=border_atoms, np0=len(mol.apos),  );  print( "edges ", edges )
#edges = bonder_edges( polys, bonds, neighs=neighs, np0=len(mol.apos) );  print( "edges ", edges )
contour, closed = mm.order_contour(edges) #;print( "contour ", contour )


plt.plot( mol.apos[is_capping,0],        mol.apos[is_capping,1],        '*', ms=25 );    # bond centers


plt.show() ;exit()


#plt.plot( points[contour,0], points[contour,1], 'o:r', lw=3  )

# beta=0.2
# for ic,ips in enumerate(conts):
#     plt.plot( (1-beta)*points[ips,0] + beta*centers[ic][None,0], (1-beta)*points[ips,1] + beta*centers[ic][None,1], 'o:', c='k' ) 

contour2  = mm.offset_contour( contour, points, l=2.0 )     #;print( "contour2 ", contour2.shape, contour2  )
plt.plot( contour2[:,0], contour2[:,1], 'o-b', lw=3  )

plt.plot( bcs[:,0], bcs[:,1], 'ob'  )


plt.show() ;exit()

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

#drawTriangles(tris, points, ax=None)

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