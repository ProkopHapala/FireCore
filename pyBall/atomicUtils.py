#!/usr/bin/python

from random import random
import numpy as np
from . import elements
#import elements

neg_types_set = { "O", "N" }

def findAllBonds( atoms, Rcut=3.0, RvdwCut=0.7 ):
    bonds     = []
    bondsVecs = []
    ps     = atoms[:,1:]
    iatoms = np.arange( len(atoms), dtype=int )
    Rcut2 = Rcut*Rcut
    for i,atom in enumerate(atoms):
        p    = atom[1:]
        dp   = ps - p
        rs   = np.sum( dp**2, axis=1 )
        for j in iatoms[:i][ rs[:i] < Rcut2 ]:
            ei = int( atoms[i,0] )
            ej = int( atoms[j,0] )
            Rcut_ij =  elements.ELEMENTS[ ei ][7] + elements.ELEMENTS[ ej ][7]
            #print ( i, j, ei, ej, Rcut_ij )
            rij =  np.sqrt( rs[j] )
            if ( rij < ( RvdwCut * Rcut_ij ) ):
                bonds.append( (i,j) )
                bondsVecs.append( ( rij, dp[j]/rij ) )
    return bonds, bondsVecs

def convert_to_adjacency_list(graph):
    adj_list = {i: list(neighbors) for i, neighbors in enumerate(graph)}
    return adj_list

def preprocess_graph(graph):
    changed = True
    while changed:
        changed = False
        to_remove = [node for node in graph if len(graph[node]) == 1]
        if to_remove:
            changed = True
            for node in to_remove:
                neighbor = graph[node][0]
                graph[neighbor].remove(node)
                del graph[node]
    return graph

def find_cycles(graph, max_length=7):
    def unblock(node, blocked, blocked_nodes):
        stack = [node]
        while stack:
            n = stack.pop()
            if n in blocked:
                blocked.remove(n)
                stack.extend(blocked_nodes[n])
                blocked_nodes[n].clear()

    def circuit(node, start, blocked, blocked_nodes, stack):
        found_cycle = False
        stack.append(node)
        blocked.add(node)
        #print( graph )
        gnd = graph.get( node, None )
        if gnd is not None:
            for neighbor in gnd:
                if neighbor == start and len(stack) <= max_length:
                    cycles.append(stack[:])
                    found_cycle = True
                elif neighbor not in blocked and len(stack) < max_length:
                    if circuit(neighbor, start, blocked, blocked_nodes, stack):
                        found_cycle = True
            if found_cycle:
                unblock(node, blocked, blocked_nodes)
            else:
                for neighbor in gnd:
                    if node not in blocked_nodes[neighbor]:
                        blocked_nodes[neighbor].append(node)
        stack.pop()
        return found_cycle

    def find_all_cycles():
        blocked = set()
        blocked_nodes = {node: [] for node in graph}
        stack = []
        nodes = list(graph.keys())
        for start in nodes:
            circuit(start, start, blocked, blocked_nodes, stack)
            while nodes and nodes[0] != start:
                nodes.pop(0)
            graph.pop(start, None)

    cycles = []
    find_all_cycles()
    return [cycle for cycle in cycles if 3 <= len(cycle) <= max_length]

def filterBonds( bonds, enames, ignore ):
    return [ (i,j) for (i,j) in bonds if not ( ( enames[i] in ignore ) or ( enames[j] in ignore ) ) ]

def colapse_to_means( bsamp, R=0.7, binds=None ):
    R2 = R*R
    cs = [ ]
    ns = [ ]
    ci = [ ]
    # n = len(bsamp)
    # bMap = False; b2c=None
    # if binds is not None: 
    #     b2c = np.full(n,-1, dtype=np.int32 )
    #     bMap = True
    for ip,p in enumerate(bsamp):
        imin  = -1
        for ic,c in enumerate(cs):
            r2 = np.sum( (p-c)**2 )
            if( r2<R2 ):
                n = ns[ic]
                cs[ic] = (c*n + p)/(n+1.0)
                ci[ic].append(ip)
                #if(bMap): b2c[ip] = ic
                ns[ic] +=1
                imin    = ic
                break
        if imin<0:
            #if(bMap): b2c[ip] = len(cs)
            cs.append( p  )
            ns.append( 1. )
            ci.append( [ip] )
    print( "centers ", cs)
    cs = np.array( cs )
    return cs, ci


def makeBondSamples( bonds, apos, where=[-0.2,0.0,0.2] ):
    bsamp = []
    for ib,(i,j) in enumerate(bonds):
        p1 = apos[i,:2]
        p2 = apos[j,:2]
        c  = 0.5*(p1+p2)
        d  = (p2-p1)
        l = np.sqrt(np.sum(d*d))
        d /=l
        q  = d[[1,0]]; q[0]*=-1
        if where is None:
            ws = [-l,l]
        else:
            ws = where    
        for w in ws: 
            bsamp.append( c + q*w )
    return np.array(bsamp)

def makeAtomSamples( neighs, apos, enames, ignore=set(['H']), where=[-0.5,0.0,0.5] ):
    samps = []
    nw = len(where)
    for ia,ngs in enumerate(neighs):
        if enames[ia] in ignore: continue
        p1 = apos[ia,:2]
        samps.append( p1 )
        if nw > 0:
            for j in ngs:
                p2 = apos[j,:2]
                d  = (p2-p1)
                d /= np.sqrt(np.sum(d*d))
                for w in where: 
                    samps.append( p1 + d*w )
    return np.array(samps)

def makeEndAtomSamples( neighs, apos,enames, ignore=set(['H']),  whereX=[-0.6, +0.6 ], whereY=[0.0,+0.6] ):
    samps = []
    nw = len(whereY)*len(whereX)
    for ia,ngs in enumerate(neighs):
        if len(ngs) != 1:        continue
        if enames[ia] in ignore: continue
        p1 = apos[ia,:2]
        (j,)  = ngs
        p2 = apos[j,:2]
        d  = (p1-p2)
        d /= np.sqrt(np.sum(d*d))
        q  = d[[1,0]]; q[0]*=-1
        for x in whereX:
            for y in whereY: 
                #print( x,y, d, q )
                samps.append( p1 + d*y + q*x )
    return np.array(samps)

def makeKinkAtomSamples( neighs, apos, where=[-0.6, +0.6 ] ):
    samps = []
    for ia,ngs in enumerate(neighs):
        if len(ngs) != 2:        continue
        p0 = apos[ia,:2]
        (i,j)  = ngs
        d1  = (apos[i,:2]-p0);  d1/=np.sqrt(np.sum(d1*d1))
        d2  = (apos[j,:2]-p0);  d2/=np.sqrt(np.sum(d2*d2))
        d = d1  + d2;           d2/=np.sqrt(np.sum(d2*d2))
        for x in where:
            #print( x,y, d, q )
            samps.append( p0 + d*x )
    return np.array(samps)

def getAtomRadius( atypes, eparams=elements.ELEMENTS, icol=6 ):
    # icol=7 RvdW, icol=6 covalent radius
    #print( eparams[ 6 ][7], eparams[ 6 ] )
    return [ eparams[ ei ][icol] for ei in atypes ]

def getAtomRadiusNP( atypes, eparams=elements.ELEMENTS ):
    return np.array( getAtomRadius( atypes, eparams ) ) 

def findBondsNP( apos, atypes=None, Rcut=3.0, RvdwCut=1.5, RvdWs=None, byRvdW=True ):
    bonds  = []
    rbs    = []
    iatoms = np.arange( len(apos), dtype=int )
    if byRvdW:
        if  RvdWs is None:
            RvdWs = getAtomRadiusNP( atypes, eparams=elements.ELEMENTS )
            #print( "findBondsNP() RvdWs=", RvdWs, RvdwCut  )
    else:
        RvdWs = np.ones(len(apos))*Rcut
    for i,pi in enumerate(apos):
        j0=i+1
        rs   = np.sqrt( np.sum( (apos[j0:,:] - pi[None,:] )**2, axis=1 ) )
        mask = rs[:] < ( RvdWs[j0:]+RvdWs[i] )*RvdwCut
        dbonds = [ (i,j) for j in iatoms[j0:][mask]  ]
        bonds += dbonds 
        rbs   += [ rs[b[1]-j0] for b in dbonds ]
    #print( "bonds=", len(bonds), bonds )
    #print( "rbs=",   len(rbs),   rbs )
    return np.array( bonds, dtype=np.int32 ), np.array( rbs )

def findHBondsNP( apos, atypes=None, Rb=1.5, Rh=2.5, angMax=60.0, typs1={"H"}, typs2=neg_types_set, bPrint=False, bHbase=False ):
    bonds  = []
    rbs    = []
    iatoms = np.arange( len(apos), dtype=int )
    cos_min = np.cos( angMax*np.pi/180.0 )
    #print( "cos_min ", cos_min )
    for i,pi in enumerate(apos):
        if ( atypes[ i ] not in typs1) : continue
        
        ds   = apos[:,:] - pi[None,:]
        rs   = np.sqrt( np.sum( ds**2, axis=1 ) )
        rs[i]=100.0
        jmin = np.argmin(rs)   # nearest neighbor i.e. bond
        
        cs   = np.dot( ds, ds[jmin] ) / ( rs*rs[jmin] )
        mask = np.logical_and( cs<-cos_min, rs<Rh )
        # print( i, atypes[i], jmin, rs[jmin], atypes[jmin] )
        # if( atypes[jmin]=="N" ):
        #     print(" ==== NEIGHS of ", i )
        #     for j in range(len(apos)):
        #         print( atypes[j], j, cs[j], rs[j], mask[j], cs[j]<-cos_min,  rs[j]<Rh )
        if(bHbase):
            dbonds = [ (i,j,jmin) for j in iatoms[:][mask] if atypes[j] in typs2 ]
        else:
            dbonds = [ (i,j) for j in iatoms[:][mask] if atypes[j] in typs2 ]
        bonds += dbonds
        rbs   += [ rs[b[1]] for b in dbonds ]
    return np.array( bonds, dtype=np.int32 ), np.array( rbs )

def neigh_bonds( natoms, bonds ):
    neighs = [{} for i in range(natoms) ]
    #print( "neigh_bonds() bonds=", bonds )
    for ib, b in enumerate(bonds):
        i = b[0]; j = b[1]; 
        neighs[i][j] = ib
        neighs[j][i] = ib
    return neighs

def neigh_atoms( natoms, bonds ):
    neighs = [ set() for i in range(natoms) ]
    for b in bonds:
        i = b[0]; j = b[1]; 
        neighs[i].add(j)
        neighs[j].add(i)
    return neighs

def findNeighOfType( ia, atypes, neighs, typ='N' ):
    result=[]
    for ja in neighs[ia]:
        #print( ia, atypes[ia], ja, atypes[ja], typ )
        if atypes[ja]==typ:
            result.append( ja )
    return result

def findNeighsOfType( selection, atypes, neighs, typ='N' ):
    found =  []
    for ia in selection:
        out = findNeighOfType( ia, atypes, neighs, typ )
        found.append( out )
    return found

def findTypeNeigh( atoms, neighs, typ, neighTyps=[(1,2,2)] ):
    '''
    find atoms of type 'typ' that have neighbors of types in 'neighTyps'
    '''
    typ_mask = ( atoms[:,0] == typ )   # boolean mask of atoms of type 'typ'
    satoms   = atoms[typ_mask]         # selected atoms of type 'typ'
    iatoms   = np.arange(len(atoms),dtype=int)[typ_mask]   # indices of selected atoms of type 'typ'
    selected = []
    for i,atom in enumerate(satoms):
        iatom = iatoms[i]
        #for jatom in neighs[ iatom ]:
        #    jtyp = atoms[jatom,0]
        count = {}   # count number of neighbors of certain type
        for jatom in neighs[ iatom ]:
            jtyp = atoms[jatom,0]
            count[jtyp] = count.get(jtyp, 0) + 1
        for jtyp, (nmin,nmax) in list(neighTyps.items()):   
            n = count.get(jtyp,0)
            if( (n>=nmin)and(n<=nmax) ):
                selected.append( iatom )
    return selected

def findTypeNeigh_( types, neighs, typ='N', neighTyps={'H':(1,2)} ):
    '''
    find atoms of type 'typ' that have neighbors of types in 'neighTyps'
    '''
    #typ_mask = ( types == typ )
    #satoms   = atoms[typ_mask]
    select = [ i for i,t in enumerate(types) if (t==typ) ]
    selected = [] 
    for ia in select:
        count = {}
        for jatom in neighs[ ia ].keys():  # count number of neighbors of certain type
            jtyp = types[jatom]
            count[jtyp] = count.get(jtyp, 0) + 1
        for jtyp, (nmin,nmax) in neighTyps.items():
            n = count.get(jtyp,0)
            if( (n>=nmin)and(n<=nmax) ):
                selected.append( ia )
    return selected


def findAngles( apos, neighs, select=None ):
    if select is None:
        select = range(len(apos))
    iang=[]
    angs=[]
    for ia in select:
        ngs=neighs[ia]
        for ja in ngs.keys():
            a = apos[ja] - apos[ia]
            for jb in ngs.keys():
                if jb<ja:
                    b = apos[jb] - apos[ia] 
                    angs.append( np.arccos( np.dot(a,b)/np.sqrt( np.dot(a,a)*np.dot(b,b) ) ) )
                    iang.append( (ja,ia,jb) )
    return angs, iang

def findDihedral( apos, enames, neighs, select, neighTyp={'H'} ):
    if select is None:
        select = range(len(apos))
    iang=[]
    angs=[]
    for ia in select:
        ngs=neighs[ia]
        if len(ngs)<3:
            print( f"ERROR in findDihedral: atom {ia} has <3 neighbors" )
        js = list(ngs.keys())
        if(len(js))<3:
            continue
        #print("------- ", js )
        for ja in js:
            if enames[ja] in neighTyp:
                js.remove(ja)
                a = apos[ja] - apos[ia]
                for jb in js:
                    for jc in js:
                        if jc>jb:
                            b = apos[jb] - apos[ia]
                            c = apos[jc] - apos[ia] 
                            n = np.cross(b,c)
                            #print( ja,jb,jc," b: ", b, " c: ",c, " n: ", n )
                            #print( " n: ", n, " a: ", a )
                            sa = np.dot(n,a) / np.sqrt(np.dot(n,n)*np.dot(a,a))
                            ca = np.sqrt( 1. - sa**2 )
                            angs.append(  np.arccos( ca ) )
                            #iang.append( (ja,ia,jb,jc) )
                            #iang.append( (ia,ja,jb,jc) )
                            #iang.append( (jb,jc,ia,ja) )
                            #iang.append( (jb,jc,ja,ia) )
                            iang.append( (ja,jb,ia,jc) )
                break
    return angs, iang 

def getAllNeighsOfSelected( selected, neighs, atoms, typs={1} ):
    result = {}
    for iatom in selected:
        for jatom in neighs[ iatom ]:
            if( atoms[jatom,0] in typs ):
                if jatom in result:
                    result[jatom].append( iatom )
                else:
                    result[jatom] = [iatom]
    return result 
    
def findPairs( select1, select2, atoms, Rcut=2.0 ):    
    ps = atoms[select2,1:]
    Rcut2 = Rcut*Rcut
    pairs = []
    select2 = np.array( select2 )
    for iatom in select1:
        p = atoms[iatom,1:]
        rs = np.sum( (ps - p)**2, axis=1 )
        for jatom in select2[ rs < Rcut2 ]:
            pairs.append( (iatom,jatom) )
    return pairs

def findPairs_one( select1, atoms, Rcut=2.0 ):    
    ps = atoms[select1,1:]
    Rcut2 = Rcut*Rcut
    pairs = []
    select1 = np.array( select1 )
    for i,iatom in enumerate(select1):
        p = atoms[iatom,1:]
        rs = np.sum( (ps - p)**2, axis=1 )
        #print ( i, iatom, rs )
        for jatom in select1[:i][ rs[:i] < Rcut2 ]:
            pairs.append( (iatom,jatom) )
    return pairs  
    
def pairsNotShareNeigh( pairs, neighs ):
    pairs_ = []
    for pair in pairs:
        ngis = neighs[ pair[0] ]
        ngjs = neighs[ pair[1] ]
        share_ng = False
        for ngi in ngis:
            if ngi in ngjs:
                share_ng = True
                break
        if not share_ng:
            pairs_.append( pair )
    return pairs_

def rotMatPCA( ps, bUnBorm=False ):
    cog = ps.sum(axis=0)/len(ps)
    ps = ps - cog[None,:]
    M = np.dot( ps.T, ps )    #;print( M )
    es, vs = np.linalg.eigh(M)
    inds = np.argsort(es)[::-1]
    #inds = np.argsort(es)
    es   = es[inds]
    vs   = vs.T[inds]
    #for ii,i in enumerate(inds):
    #    print( ['x','y','z'][ii], es[ii], vs[ii] )
    #print(es)
    #print(vs)
    if(bUnBorm):
        emax=es.max()
        for i in range(3): vs[i]*=(es[i]/emax)
    return( vs )
    
def makeRotMatAng( ang, ax=(0,1) ):
    ax1,ax2=ax
    ca=np.cos(ang)
    sa=np.sin(ang)
    rot=np.eye(3)
    rot[ax1,ax1]=ca;  rot[ax2,ax2]=ca
    rot[ax1,ax2]=-sa; rot[ax2,ax1]=sa
    return rot 

def makeRotMat( fw, up ):    
    fw   = fw/np.linalg.norm(fw)
    up   = up - fw*np.dot(up,fw)
    ru   = np.linalg.norm(up)
    if ru<1e-4:  # if colinear
        print("WARRNING: makeRotMat() up,fw are colinear => randomize up ")
        up = np.random.rand(3)
        up = up - fw*np.dot(up,fw)
        ru = np.linalg.norm(up)
    up   = up/ru
    left = np.cross(fw,up)
    left = left/np.linalg.norm(left) 
    return np.array([left,up,fw])

def makeRotMatAng2( fw, up, ang ):    
    fw   = fw/np.linalg.norm(fw)
    up   = up - fw*np.dot(up,fw)
    ru   = np.linalg.norm(up)
    if ru<1e-4:  # if colinear
        print("WARRNING: makeRotMat() up,fw are colinear => randomize up ")
        up = np.random.rand(3)
        up = up - fw*np.dot(up,fw)
        ru = np.linalg.norm(up)
    up   = up/ru
    left = np.cross(fw,up)
    left = left/np.linalg.norm(left)
    ca = np.cos(ang)
    sa = np.sin(ang)
    fw_ = fw*ca + up*-sa
    up_ = fw*sa + up*ca
    return np.array([left,up_,fw_])

def rotmat_from_points( ps, ifw=None, iup=None,  fw=None, up=None, _0=1 ):
    if fw is None: fw  = ps[ifw[1]-_0]-ps[ifw[0]-_0]
    if up is None: up  = ps[iup[1]-_0]-ps[iup[0]-_0]
    return makeRotMat( fw, up )

def mulpos( ps, rot ):
    for ia in range(len(ps)): ps[ia,:] = np.dot(rot, ps[ia,:])

def orient_vs( p0, fw, up, apos, trans=None ):
    rot = makeRotMat( fw, up )
    if trans is not None: rot=rot[trans,:]   #;print( "orient() rot=\n", rot ) 
    apos[:,:]-=p0[None,:]
    mulpos( apos, rot )
    return apos

def orient( i0, ip1, ip2, apos, _0=1, trans=None, bCopy=True ):
    if(bCopy): apos=apos.copy()
    p0  = apos[i0-_0]
    fw  = apos[ip1[1]-_0]-apos[ip1[0]-_0]
    up  = apos[ip2[1]-_0]-apos[ip2[0]-_0]
    return orient_vs( p0, fw, up, apos, trans=trans )

def orientPCA(ps, perm=None):
    M = rotMatPCA( ps )
    if perm is not None: M = M[perm,:]
    mulpos( ps, M )

def groupToPair( p1, p2, group, up, up_by_cog=False ):
    center = (p1+p2)*0.5
    fw  = p2-p1;    
    if up_by_cog:
        up  = center - up
    rotmat = makeRotMat( fw, up )
    ps  = group[:,1:]
    #ps_ = ps
    #ps_ = np.transpose( np.dot( rotmat, np.transpose(ps) ) )
    ps_ = np.dot( ps, rotmat ) 
    group[:,1:] = ps_ + center
    return group
    
def replacePairs( pairs, atoms, group, up_vec=(np.array((0.0,0.0,0.0)),1) ):
    replaceDict = {}
    for ipair,pair in enumerate(pairs):
        for iatom in pair:
            replaceDict[iatom] = 1
            #if( iatom in replaceDict ):
            #    replaceDict[iatom].append(ipair)
            #else:
            #    replaceDict[iatom] = [ipair]
    atoms_ = []
    for iatom,atom in enumerate( atoms ):
        if(iatom in replaceDict): continue
        atoms_.append(atom)
    for pair in pairs:
        group_ = groupToPair( atoms[pair[0],1:], atoms[pair[1],1:], group.copy(), up_vec[0], up_vec[1] )
        for atom in group_:
            atoms_.append( atom )
        #break
    return atoms_

def findNearest( p, ps, rcut=1e+9 ):
    rs = np.sum( (ps - p)**2, axis=1 )
    imin = np.argmin(rs)
    if rs[imin]<(rcut**2):
        return imin
    else: 
        return -1 

def countTypeBonds( atoms, ofAtoms, rcut ):
    bond_counts = np.zeros(len(atoms), dtype=int )
    ps = ofAtoms[:,1:]
    for i,atom in enumerate(atoms):
        p = atom[1:]
        rs = np.sum( (ps - p)**2, axis=1 )
        bond_counts[i] = np.sum( rs < (rcut**2) )
    return bond_counts
	
def findBondsTo( atoms, typ, ofAtoms, rcut ):
    found     = []
    foundDict = {}
    ps = ofAtoms[:,1:] 
    for i,atom in enumerate(atoms):
        if atom[0]==typ:
            p = atom[1:]
            ineigh = findNearest( p, ps, rcut )
            if ineigh >= 0:
                foundDict[i] = len(found)
                found.append( (i, p - ps[ineigh]) )
    return found, foundDict
	
def replace( atoms, found, to=17, bond_length=2.0, radial=0.0, prob=0.75 ):
    replace_mask = np.random.rand(len(found)) < prob
    for i,foundi in enumerate(found):
        if replace_mask[i]:
            iatom = foundi[0]
            bvec  = foundi[1]
            rb    = np.linalg.norm(bvec)
            bvec *= ( bond_length - rb )/rb
            #if radial > 0:
            #    brad = atoms[iatom,1:]
            #    brad = brad/np.linalg.norm(brad)
            #    cdot = np.dot( brad, bvec )
            #    bvec = (1-radial)*bvec + brad*radial/cdot
            atoms[iatom,0]   = to
            atoms[iatom,1:] += bvec  
    return atoms

def saveAtoms( atoms, fname, xyz=True ):
    fout = open(fname,'w')
    fout.write("%i\n"  %len(atoms) )
    if xyz==True : fout.write("\n") 
    for i,atom in enumerate( atoms ):
        if isinstance( atom[0], str ):
            fout.write("%s %f %f %f\n"  %( atom[0], atom[1], atom[2], atom[3] ) )
        else:
            fout.write("%i %f %f %f\n"  %( atom[0], atom[1], atom[2], atom[3] ) )
    fout.close() 


def psi4frags2string( enames, apos, frags=None ):
    n = len(enames)
    s = []
    if frags is None: frags = [ range(n) ]
    for i,frag in enumerate(frags):
        #print( i, frag )
        if i>0: s.append( "--" )
        for ia in frag:
            xyz = apos[ia]
            s.append( "%s %f %f %f"  %( enames[ia], xyz[0], xyz[1], xyz[2]) )
    #print("s = ", s)
    return "\n".join(s)

def writeToXYZ( fout, es, xyzs, qs=None, Rs=None, comment="#comment", bHeader=True, ignore_es=None, other_lines=None ):
    na=len(xyzs)
    if(bHeader):
        if other_lines is not None:
            na += len(other_lines)
        fout.write("%i\n"  %na )
        fout.write(comment+"\n")
    if ignore_es is not None:
        mask = [ (  e not in ignore_es) for e in es ]
        na = sum(mask)
    else:
        mask = [True]*na
    #print( "writeToXYZ len(es,xyzs,qs,mask) ", len(es), len(xyzs),  len(qs), len(mask) )
    # print( "writeToXYZ es ", es )
    # print( "writeToXYZ qs ", qs )
    # print( "writeToXYZ xyzs ", xyzs )
    if   (Rs is not None):
        for i,xyz in enumerate( xyzs ):
            if mask[i]: fout.write("%s %f %f %f %f %f \n"  %( es[i], xyz[0], xyz[1], xyz[2], qs[i], Rs[i] ) )
    elif (qs is not None):
        for i,xyz in enumerate( xyzs ):
            if mask[i]: fout.write("%s %f %f %f %f\n"  %( es[i], xyz[0], xyz[1], xyz[2], qs[i] ) )
    else:
        for i,xyz in enumerate( xyzs ):
            if mask[i]: fout.write("%s %f %f %f\n"  %( es[i], xyz[0], xyz[1], xyz[2] ) )
    if other_lines is not None:
        for l in other_lines:
            fout.write(l)

def saveXYZ( es, xyzs, fname, qs=None, Rs=None, mode="w", comment="#comment", ignore_es=None, other_lines=None ):
    fout = open(fname, mode )
    writeToXYZ( fout, es, xyzs, qs, Rs=Rs, comment=comment, ignore_es=ignore_es, other_lines=other_lines )
    fout.close() 

'''
def savePDB( es, xyzs, bonds, fname, mode="w" ):
    fout = open( fname, mode )
    fout.write("COMPND    UNNAMED\n")
    fout.write("AUTHOR    GENERATED BY FireCore\n")
    for i,xyz in enumerate( xyzs ):
        fout.write("HETATM %i %s  UNL   1    %f %f %f  1.0 0.0     %s \n"  %( es[i], xyz[0], xyz[1], xyz[2], qs[i], Rs[i] ) )
    for i,b in enumerate( xyzs ):
        fout.write("HETATM %i %s  UNL   1    %f %f %f  1.0 0.0     %s \n"  %( i, b[0], b[1], xyz[2], qs[i], Rs[i] ) )
'''

def makeMovie( fname, n, es, func ):
    fout = open(fname, "w")
    for i in range(n):
        xyzs, qs = func(i)
        writeToXYZ( fout, es, xyzs, qs, comment=("frame %i " %i) )
    fout.close() 

def loadAtomsNP(fname=None, fin=None, bReadN=False, nmax=10000, comments=None ):
    #print(" HELLO !!!!!!!", fname)
    bClose=False
    if fin is None: 
        fin=open(fname, 'r')
        bClose=True
    xyzs   = [] 
    Zs     = []
    enames = []
    qs     = []
    ia=0
    for line in fin:
        if comments is not None:
            if line[0]=='#':
                comments.append(line)
                continue        
        wds = line.split()
        try:
            #print( "line", line )
            #print( "wds", wds )
            xyz = ( float(wds[1]), float(wds[2]), float(wds[3]) )
            try:
                q = float(wds[4])
            except:
                q = 0
            try:
                iz    = int(wds[0]) 
                ename = elements.ELEMENTS[iz-1][1]
            except:
                ename = wds[0]
                iz    = elements.ELEMENT_DICT[ename][0]
            enames.append( ename )
            Zs    .append( iz    )
            qs    .append( q     )
            xyzs  .append( xyz   )
            ia+=1
        except:
            #print("loadAtomsNP("+fname+")cannot interpet line: ", line)
            if bReadN and (ia==0):
                try:
                    nmax=int(wds[0])
                    #print("nmax: ", nmax)
                except:
                    pass
        if(ia>=nmax): break
    if(bClose): fin.close()
    xyzs = np.array( xyzs )
    Zs   = np.array( Zs, dtype=np.int32 )
    qs   = np.array( qs )
    #print( len(enames), enames )
    #print( len(Zs), Zs )
    #print( len(xyzs), xyzs )
    #print( "loadAtomsNP ", fname, len(xyzs)  ,len(Zs),len(enames),len(qs) )
    return xyzs,Zs,enames,qs


def load_xyz(fname=None, fin=None, bReadN=False, bReadComment=True, nmax=10000 ):
    bClose=False
    if fin is None: 
        fin=open(fname, 'r')
        bClose=True
    xyzs   = [] 
    Zs     = []
    enames = []
    qs     = []
    comment = None
    ia=0
    for line in fin:
        wds = line.split()
        try:
            xyz = ( float(wds[1]), float(wds[2]), float(wds[3]) )
            try:
                q = float(wds[4])
            except:
                q = 0
            try:
                iz    = int(wds[0]) 
                ename = elements.ELEMENTS[iz-1][1]
            except:
                ename = wds[0]
                iz    = elements.ELEMENT_DICT[ename][0]
            enames.append( ename )
            Zs    .append( iz    )
            qs    .append( q     )
            xyzs  .append( xyz   )
            ia+=1
        except:
            #print("cannot interpet line: ", line)
            if bReadN and (ia==0):
                try:
                    nmax=int(wds[0])
                except:
                    comment=line.strip()    
        if(ia>=nmax): break
    if(bClose): fin.close()
    xyzs = np.array( xyzs )
    Zs   = np.array( Zs, dtype=np.int32 )
    qs   = np.array( qs )
    return xyzs,Zs,enames,qs,comment

def string_to_matrix( s, nx=3,ny=3, bExactSize=False ):
    elements = []
    for item in s.split():
        try:  # Try to convert each element to a float
            elements.append(float(item))
        except ValueError:  # If conversion fails, ignore (e.g., "lvs")
            continue
    n=len(elements)
    nxy=nx*ny
    if (n<nxy):
        print( f"string_to_matrix(): n({n})<nx({nx})*ny(ny)" )
        exit()
    elif (n>nxy) and (bExactSize):
        print( f"string_to_matrix(): n({n})>nx({nx})*ny(ny)" )
        exit()
    else:
        elements=elements[:nxy]
    matrix = np.array(elements).reshape(nx,ny)   # Convert the list of elements into a 3x3 NumPy array
    return matrix


def loadMol(fname=None, fin=None, bReadN=False, nmax=10000 ):
    bClose=False
    if fin is None: 
        fin=open(fname, 'r')
        bClose=True
    xyzs   = [] 
    Zs     = []
    enames = []
    qs     = []
    bonds  = []
    ia=0
    na=0
    nb=0
    for il,line in enumerate(fin):
        wds = line.split()
        if (il>4) and ((il-4)<na):
            xyzs.append( ( float(wds[0]), float(wds[1]), float(wds[2]) ) )
            ename = wds[0]
            enames.append( ename )
            Zs    .append( elements.ELEMENT_DICT[ename][0] )
            ia+=1
            if (bReadN and ia>=na): break
        if (il>(4+na) ) and ((il-4-na)<nb):
            bonds.append( wds[1],wds[2] )  # ignoring bond order
        elif(il==4):
            na=int(wds[0])
            nb=int(wds[2])
    if(bClose): fin.close()
    xyzs  = np.array( xyzs )
    Zs    = np.array( Zs, dtype=np.int32 )
    qs    = np.array(qs)
    bonds = np.array( bonds, dtype=np.int32 )
    return xyzs,Zs,enames,qs,bonds


def loadMol2(fname, bReadN=True):
    """
    Load an AtomicSystem from a .mol2 file.
    
    The mol2 file is expected to contain at least the following sections:
      - @<TRIPOS>MOLECULE
      - @<TRIPOS>ATOM
      - @<TRIPOS>BOND   (optional: if bonds exist)
    
    In the MOLECULE section, if a comment line starting with '#' is present
    and contains "lvs", the 9 numbers following "lvs" will be parsed as the
    lattice vectors (lvec). They are arranged row‐wise in a 3×3 numpy array.
    
    The ATOM section is assumed to have lines in the format:
    
          atom_id  atom_name  x  y  z  atom_type  substructure_id  residue_name  charge
    
    The element symbol is taken as follows:
      - If the atom_name (the second token) is a valid chemical symbol (i.e.
        found in elements.ELEMENT_DICT), that is used.
      - Otherwise, the part before the period in the atom_type token is used.
    
    The atomic number is determined via elements.ELEMENT_DICT.
    
    The BOND section (if present) is assumed to have lines in the format:
    
          bond_id  origin_atom_id  target_atom_id  bond_type
    
    and bonds are returned as zero-based index tuples.
    
    Parameters:
      fname (str): The name of the mol2 file.
      bReadN (bool): (Unused here; provided for compatibility.)
    
    Returns:
      An AtomicSystem instance with fields:
         - apos: numpy array (N,3) of atomic coordinates.
         - atypes: numpy array (N,) of atomic numbers.
         - enames: numpy array (N,) of element symbols (strings).
         - qs: numpy array (N,) of charges.
         - bonds: numpy array of shape (nB,2) with bonds (zero-based).
         - lvec: if a lattice comment is found, a (3,3) numpy array.
    """
    apos    = []
    atypes  = []
    enames  = []
    qs      = []
    bonds   = []
    lvec    = None

    with open(fname, 'r') as fin:
        lines = fin.readlines()

    # --- First, search for the MOLECULE section and a lattice comment if present.
    in_molecule = False
    in_atom     = False
    in_bond     = False

    # The counts (number of atoms, bonds, etc.) can be parsed from the MOLECULE section;
    # however we won’t depend on that for reading.
    for i, line in enumerate(lines):
        line = line.strip()
        if line.upper().startswith("@<TRIPOS>MOLECULE"):
            in_molecule = True
            in_atom = False
            in_bond = False
            continue
        if in_molecule:
            # If the line starts with '#' and contains 'lvs', parse lattice vectors.
            if line.startswith('#') and 'lvs' in line:
                # For example: "# lvs   20.0 0.0 0.0   0.0 5.0 0.0    0.0 0.0 20.0"
                parts = line.split()
                try:
                    idx = parts.index('lvs')
                except ValueError:
                    idx = -1
                if idx >= 0 and len(parts) >= idx+10:
                    # Read next 9 numbers:
                    nums = [ float(parts[idx + j]) for j in range(1, 10) ]
                    lvec = np.array(nums).reshape(3,3)
            # We assume the MOLECULE block ends when we hit the next section header.
            if line.upper().startswith("@<TRIPOS>ATOM"):
                in_molecule = False
                in_atom = True
                continue
        # --- Now, process the ATOM section.
        if line.upper().startswith("@<TRIPOS>ATOM"):
            in_atom = True
            in_bond = False
            continue
        if line.upper().startswith("@<TRIPOS>BOND"):
            in_atom = False
            in_bond = True
            continue

        if in_atom:
            if line == "" or line.startswith("@<TRIPOS>"):
                continue
            # Split the line.
            # Expected tokens:
            # 0: atom_id (integer)
            # 1: atom_name (string)
            # 2,3,4: x, y, z coordinates
            # 5: atom_type (string, e.g., "C.1", "O.3", etc.)
            # 6: substructure_id (can be ignored)
            # 7: residue_name (can be ignored)
            # 8: charge (optional)
            tokens = line.split()
            if len(tokens) < 6:
                continue  # skip malformed lines
            try:
                x = float(tokens[2])
                y = float(tokens[3])
                z = float(tokens[4])
            except:
                continue
            apos.append( [x, y, z] )
            
            # Determine the element symbol.
            # First try the atom_name (token[1]).
            cand = tokens[1]
            # Check if cand (or its first letter) is in the ELEMENT_DICT.
            if cand in elements.ELEMENT_DICT:
                elem = cand
            else:
                # Otherwise, use the atom_type (token[5]) and take the part before '.'
                elem = tokens[5].split('.')[0]
            enames.append( elem )
            
            # Determine atomic number.
            try:
                # Using the dictionary that maps element symbol to a tuple.
                # For example: elements.ELEMENT_DICT["C"][0] should be 6.
                znum = elements.ELEMENT_DICT[elem][0]
            except KeyError:
                # If not found, try to convert cand to an integer.
                try:
                    znum = int(cand)
                except:
                    znum = 0
            atypes.append( znum )
            
            # Charge is the last token if present (some mol2 files provide it).
            if len(tokens) >= 9:
                try:
                    charge = float(tokens[8])
                except:
                    charge = 0.0
            else:
                charge = 0.0
            qs.append( charge )
            
        if in_bond:
            if line == "" or line.startswith("@<TRIPOS>"):
                continue
            # Expected tokens for bonds:
            # 0: bond_id
            # 1: origin_atom_id
            # 2: target_atom_id
            # 3: bond_type (ignored here)
            tokens = line.split()
            if len(tokens) < 3:
                continue
            try:
                iatom = int(tokens[1]) - 1  # convert to zero-based index
                jatom = int(tokens[2]) - 1
                bonds.append( (iatom, jatom) )
            except:
                continue

    # Convert lists to numpy arrays:
    apos_np   = np.array(apos, dtype=float)
    atypes_np = np.array(atypes, dtype=int)
    qs_np     = np.array(qs, dtype=float)
    enames_np = np.array(enames)
    #bonds_np  = np.array(bonds, dtype=int)

    # Create an AtomicSystem instance.
    #system = AtomicSystem(apos=apos_np, atypes=atypes_np, enames=enames_np, qs=qs_np, bonds=bonds, lvec=lvec)
    
    return apos_np, atypes_np, enames_np, qs_np, bonds


def readAtomsXYZ( fin, na ):
    apos=[]
    es  =[] 
    for i in range(na):
        ws = fin.readline().split()
        apos.append( [ float(ws[1]),float(ws[2]),float(ws[3]) ] )
        es  .append( ws[0] )
    return np.array(apos), es

def read_lammps_lvec( fin ):
    '''
    https://docs.lammps.org/Howto_triclinic.html
    a = (xhi-xlo,0,0); b = (xy,yhi-ylo,0); c = (xz,yz,zhi-zlo). 
    '''
    xlo_bound,xhi_bound,xy = ( float(w) for w in fin.readline().split() )
    ylo_bound,yhi_bound,xz = ( float(w) for w in fin.readline().split() )
    zlo,zhi,yz             = ( float(w) for w in fin.readline().split() )
    xlo = xlo_bound - np.min( [0.0,xy,xz,xy+xz] )
    xhi = xhi_bound - np.max( [0.0,xy,xz,xy+xz] )
    ylo = ylo_bound - min(0.0,yz)
    yhi = yhi_bound - max(0.0,yz)
    return np.array( ( (xhi-xlo,0.,0.), (xy,yhi-ylo,0.), (xz,yz,zhi-zlo) ) )

def readLammpsTrj(fname=None, fin=None, bReadN=False, nmax=100, selection=None ):

    bClose=False
    if fin is None: 
        fin=open(fname, 'r')
        bClose=True
    trj  = []
    isys = 0
    while True:
        line = fin.readline()
        if not line: break
        if( len(line)>=5 ):
            if(line[:5]=="ITEM:"):
                wds = line.split()
                if wds[1]=='NUMBER':
                    na = int(fin.readline())
                    isys+=1
                elif( isys in selection ):
                    if wds[1]=='BOX':
                        lvec = read_lammps_lvec( fin )
                    elif wds[1]=='ATOMS':
                        apos,es = readAtomsXYZ( fin, na )
                        S = AtomicSystem( lvec=lvec, enames=es, apos=apos)
                        trj.append( S )
    return trj

def loadAtoms( name ):
    f = open(name,"r")
    n=0;
    l = f.readline()
    try:
        n=int(l)
    except:
        raise ValueError("First line of a xyz file should contain the number of atoms. Aborting...")
    line = f.readline() 
    if (n>0):
        n=int(l)
        e=[];x=[];y=[]; z=[]; q=[]
        i = 0;
        for line in f:
            words=line.split()
            nw = len( words)
            ie = None
            if( nw >=4 ):
                e.append( words[0] )
                x.append( float(words[1]) )
                y.append( float(words[2]) )
                z.append( float(words[3]) )
                if ( nw >=5 ):
                    q.append( float(words[4]) )
                else:
                    q.append( 0.0 )
                i+=1
            else:
                print(" skipped line : ", line)
    f.close()
    return [ e,x,y,z,q ]


def load_xyz_movie( fname ):
    f = open(fname,"r")
    il=0
    imgs=[]
    while True:
        line = f.readline()
        if il==0:
            n=int(line.split()[0]) 
        else:
            apos=np.array((n,3))
            es  =[]
            for i in range(n):
                line = f.readline()
                words=line.split()
                es.append( words[0] )
                apos[i,0] = float(words[1])
                apos[i,1] = float(words[1])
                apos[i,2] = float(words[1])
                il+=1
            imgs.append( (apos,es) )
    return imgs

#def loadCoefs( characters=['s','px','py','pz'] ):
def loadCoefs( characters=['s'] ):
    dens = None
    coefs = []
    for char in characters:
        fname  = 'phi_0000_%s.dat' %char
        raw = np.genfromtxt(fname,skip_header=1)
        Es  = raw[:,0]
        cs  = raw[:,1:]
        sh  = cs.shape
        cs  = cs.reshape(sh[0],sh[1]//2,2)
        d   = cs[:,:,0]**2 + cs[:,:,1]**2
        coefs.append( cs[:,:,0] + 1j*cs[:,:,1] )
        if dens is None:
            dens  = d 
        else:
            dens += d
    return dens, coefs, Es

def findCOG( ps, byBox=False ):
    if(byBox):
        xmin=ps[:,0].min(); xmax=ps[:,0].max();
        ymin=ps[:,1].min(); ymax=ps[:,1].max();
        zmin=ps[:,2].min(); zmax=ps[:,2].max();
        return np.array( (xmin+xmax, ymin+ymax, zmin+zmax) ) * 0.5
    else:
        cog = np.sum( ps, axis=0 )
        cog *=(1.0/len(ps))
        return cog

def projectAlongBondDir( apos, i0, i1 ):
    dir = (apos[i1]-apos[i0])
    dir*=(1/np.sqrt(np.dot(dir,dir)))   # normalize
    prjs = np.dot( apos, dir[:,None] )  #;print(prjs)
    return prjs

def histR( ps, dbin=None, Rmax=None, weights=None ):
    rs = np.sqrt(np.sum((ps*ps),axis=1))
    bins=100
    if dbin is not None:
        if Rmax is None:
            Rmax = rs.max()+0.5
        bins = np.linspace( 0,Rmax, int(Rmax/(dbin))+1 )
    return np.histogram(rs, bins, weights=weights)

# ================= Topology Builder

def addGroup( base, group, links ):
    A0s,B0s = base
    A1s,B1s = group
    n0 = len(A0s)
    n1 = len(A1s)
    As = A0s + A1s
    Bs = B0s + []
    for b in links:
        b_ = (  b[0], b[1]+n0 )
        Bs.append(b_)
        As[b_[0]] += 1
        As[b_[1]] += 1
    for b in B1s:
        b_ = (  b[0]+n0, b[1]+n0 )
        Bs.append(b_)

def addBond( base, link, bNew=True ):
    As,Bs = base
    if bNew:
        As=As+[]
        Bs=Bs+[]
    As[link[0]] += 1
    As[link[1]] += 1
    Bs.append( link )
    return (As,Bs)


def disolveAtom( base, ia ):
    A0s,B0s = base
    Bs     = []
    neighs = []
    for b in B0s:
        i=b[0]
        j=b[1]
        if(i>ia):
            i-=1
        if(j>ia):
            j-=1
        
        if b[0]==ia:
            neighs.append(j)
            continue
        if b[1]==ia:
            neighs.append(i)
            continue
        print( "add B ", (i,j), b )
        Bs.append( (i,j) )
    As = list(A0s) + []
    nng = len(neighs)
    
    if( nng ==2 ):
        Bs.append( (neighs[0],neighs[1]) )
    else:
        print("ERROR: disolveAtoms applicable only for atom with 2 neighbors ( not %i )" %nng )
        print( neighs )
        exit()
    
    #for i in neighs:
    #    As[i]-=1
    old_i = list( range(len(As)) )
    old_i.pop(ia)
    As.pop(ia)
    len( As )
    print("old_i", old_i)
    print( len(As), len(old_i) )
    return (As,Bs), old_i

def removeGroup( base, remove ):
    remove = set(remove)
    A0s,B0s = base
    #left = []
    As = []
    new_i = [-1]*len(A0s)
    old_i = []
    j = 0
    for i,a in enumerate(A0s):
        if not (i in remove):
            As   .append( a )
            old_i.append( i )
            new_i[i] = j
            j+=1
    Bs = []
    for b in B0s:
        ia=b[0]
        ja=b[1]
        bi=( ia in remove )
        bj=( ja in remove )
        ia_ = new_i[ia]
        ja_ = new_i[ja]
        if not( bi or bj ):
            Bs.append( (ia_,ja_) )
        elif bi != bj:
            if bi:  # ja is in, ia not
                As[ ja_ ] -=1
            else :  # ia is in, ja not
                As[ ia_ ] -=1
    return (As,Bs), old_i

def selectBondedCluster( s, bonds ):
    #s = { i0 }
    for i in range( len(bonds) ):
        n = len(s)
        for b in bonds:
            if   b[0] in s: s.add( b[1] )
            elif b[1] in s: s.add( b[0] )
        if( len(s) <= n ): break
    return s

def scan_xyz( fxyzin, callback=None, kwargs=None ):
    fin =open(fxyzin,'r')
    i=0
    results = []
    while True:
        apos,Zs,es,qs,comment = load_xyz( fin=fin, bReadN=True )
        if(len(es)==0): break
        if callback is not None: 
            if( kwargs is not None ): 
                res = callback( (apos,es), id=i, **kwargs, comment=comment )
            else:
                res = callback( (apos,es), id=i, comment=comment )
            results.append( res )
        i+=1
    return results

def geomLines( apos, enames ):
    lines = []
    for i,pos in enumerate(apos):
        lines.append(  "%s %3.5f %3.5f %3.5f\n" %(enames[i], pos[0],pos[1],pos[2]) )
    return lines

def tryAverage( ip, apos, _0=1 ):
    if hasattr(ip, '__iter__'):
        ip = np.array(ip) - _0
        p0  = apos[ip].mean(axis=0)
    else:
        p0 = apos[ip-_0]
    return p0

def makeVectros( apos, ip0, b1, b2, _0=1 ):
    p0 = tryAverage( ip0, apos, _0=_0 )
    if ( b1==None ):
        return p0, None, None
    fw = tryAverage( b1[1], apos, _0=_0 ) - tryAverage( b1[0], apos, _0=_0 )
    if ( ( b2==None ) or (len(apos)<3) ):
        up = np.cross( fw, np.random.rand(3) )
        up = up/np.linalg.norm(up)
    else:
        up = tryAverage( b2[1], apos, _0=_0 ) - tryAverage( b2[0], apos, _0=_0 )
    return p0, fw, up

def loadElementTypes( fname='ElementTypes.dat', bDict=False ):
    lst = []
    with open(fname,'r') as fin:
        lines = fin.readlines()
        for line in lines:
            if( line[0]=='#' ): continue
            wds = line.split()
            # He        2   2   0   0   0xFFC0CB  0.849     1.1810    0.00242838984   0.098   0.00000000000   0.00000000000
            name = wds[0]
            rec = [ name ] + [ int(w) for w in wds[1:4] ] + [ wds[5] ] + [ float(w) for w in wds[6:12] ]
            lst.append( rec )
    if bDict: return { rec[0]:rec for rec in lst }
    return lst

def getVdWparams( iZs, etypes=None, fname='ElementTypes.dat' ):
    if etypes is None: etypes = loadElementTypes( fname=fname, bDict=False )
    return np.array( [( etypes[i][6],etypes[i][7] ) for i in iZs  ] )

def iz2enames( iZs ):
    return [ elements.ELEMENTS[iz-1][1] for iz in iZs ]

def atoms_symmetrized( atypes, apos, lvec, qs=None, REQs=None, d=0.1):
    """
    Symmetrize atoms in a unit cell by replicating atoms near the cell boundaries.

    Parameters:
    - n (int): Number of atoms.
    - atypes (np.ndarray): Array of atom types with shape (n,).
    - apos (np.ndarray): Array of atom positions with shape (n, 3).
    - REQs (np.ndarray): Array of quaternions with shape (n, 4).
    - grid_cell (np.ndarray): 3x3 matrix representing the unit cell vectors as columns.
    - d (float): Threshold distance from the cell boundaries (default is 0.1).

    Returns:
    - new_atypes (np.ndarray): Array of symmetrized atom types.
    - new_apos (np.ndarray): Array of symmetrized atom positions.
    - new_REQs (np.ndarray): Array of symmetrized quaternions.
    """
    n = len(atypes)
    # Compute inverse transformation matrix M
    M = np.linalg.inv(lvec)

    # Define boundary thresholds
    cmax = -0.5 + d
    cmin =  0.5 - d

    # Extract lattice vectors a and b from grid_cell
    a = lvec[:, 0]  # First column
    b = lvec[:, 1]  # Second column

    # Transform atom positions using the inverse matrix M
    p_transformed = apos @ M.T  # Shape: (n, 3)
    p_a = p_transformed[:, 0]
    p_b = p_transformed[:, 1]

    # Determine if atoms are near the boundaries in a and b directions
    alo = p_a < cmax
    ahi = p_a > cmin
    blo = p_b < cmax
    bhi = p_b > cmin

    aa = alo | ahi  # Atoms near the a-direction boundaries
    bb = blo | bhi  # Atoms near the b-direction boundaries

    # Calculate weighting factor based on replica count
    ws = 1.0 / ((1 + aa.astype(float)) * (1 + bb.astype(float)))

    bREQs = REQs is not None
    bQs   = qs   is not None

    new_REQs = None
    if bREQs:
        REQs_adj = REQs.copy()
        REQs_adj[:, 2] *= ws  # Adjust Q
        REQs_adj[:, 1] *= ws  # Adjust E0
        new_REQs = list(REQs_adj)

    new_qs = None
    if bQs:
        qs_adj = qs.copy()
        qs_adj *= ws
        new_qs = list(qs_adj)

    # Initialize lists with original atoms
    new_atypes = list(atypes)
    new_apos   = list(apos)
    new_ws     = list(ws)
    
    # Determine shifts based on boundary conditions
    shift_a = np.where(alo[:, np.newaxis], a, -a)  # Shape: (n, 3)
    shift_b = np.where(blo[:, np.newaxis], b, -b)  # Shape: (n, 3)

    # Replicate atoms shifted by a
    if np.any(aa):
        indices_a = np.where(aa)[0]
        new_atypes.extend(atypes[indices_a])
        new_apos.extend(apos[indices_a] + shift_a[indices_a])
        new_ws.extend( ws[indices_a] )
        if bREQs: new_REQs.extend(REQs_adj[indices_a])
        if bQs:   new_qs  .extend(qs_adj[indices_a])

    # Replicate atoms shifted by b
    if np.any(bb):
        indices_b = np.where(bb)[0]
        new_atypes.extend(atypes[indices_b])
        new_apos  .extend(apos[indices_b] + shift_b[indices_b])
        new_ws    .extend( ws[indices_b] )
        if bREQs: new_REQs.extend(REQs_adj[indices_b])
        if bQs:   new_qs.  extend(qs_adj[indices_b])

        # Replicate atoms shifted by both a and b
        indices_ab = np.where(aa & bb)[0]
        if len(indices_ab) > 0:
            new_atypes.extend(atypes[indices_ab])
            new_apos  .extend(apos[indices_ab] + shift_a[indices_ab] + shift_b[indices_ab])
            new_ws    .extend( ws[indices_ab] )
            if bREQs:  new_REQs.extend(REQs_adj[indices_ab])
            if bQs:    new_qs.extend(qs_adj[indices_ab])

    # Convert lists back to NumPy arrays
    new_atypes = np.array(new_atypes, dtype=atypes.dtype )
    new_apos   = np.array(new_apos,   dtype=apos.dtype   )
    new_ws     = np.array(new_ws,     dtype=ws.dtype   )
    if bREQs: new_REQs   = np.array(new_REQs, dtype=REQs.dtype   )
    if bQs:   new_qs     = np.array(new_qs,   dtype=qs.dtype     )

    return new_atypes, new_apos, new_qs, new_REQs, new_ws

# ========================== Class Geom

class AtomicSystem( ):

    def __init__(self,fname=None, apos=None, atypes=None, enames=None, lvec=None, qs=None, Rs=None, bonds=None, ngs=None, bReadN=True ) -> None:
        self.apos    = apos
        self.atypes  = atypes
        self.enames  = enames
        self.qs      = qs
        self.Rs      = Rs
        self.bonds   = bonds
        self.ngs     = ngs 
        self.lvec    = lvec
        self.aux_labels = None
        if fname is not None:
            ext = fname.split('.')[-1]
            #print( f"AtomicSystem.__init__({fname}) ext=", ext  )
            if( 'mol' == ext ):
                self.apos,self.atypes,self.enames,self.qs,self.bonds = loadMol(fname=fname, bReadN=bReadN )
            elif ( 'mol2' == ext ):
                self.apos,self.atypes,self.enames,self.qs,self.bonds = loadMol2(fname=fname, bReadN=bReadN )
            elif ( 'xyz' == ext ):
                self.apos,self.atypes,self.enames,self.qs, comment = load_xyz(fname=fname, bReadN=bReadN )
                if comment is not None:
                    if comment[:3] == 'lvs':      
                        self.lvec = string_to_matrix( comment, nx=3,ny=3, bExactSize=False )
                        #print( f"AtomicSystem.__init__({fname}) lvec=\n", self.lvec   )
                #print( f"AtomicSystem.__init__({fname}) comment=", comment  )
            else:
                self.apos,self.atypes,self.enames,self.qs = loadAtomsNP(fname=fname , bReadN=bReadN )

    def saveXYZ(self, fname, mode="w", blvec=True, comment="", ignore_es=None, bQs=True, other_lines=None ):
        if blvec and (self.lvec is not None):
            #print( self.lvec )
            comment= ( "lvs %6.3f %6.3f %6.3f   %6.3f %6.3f %6.3f   %6.3f %6.3f %6.3f" %(self.lvec[0,0],self.lvec[0,1],self.lvec[0,2],  self.lvec[1,0],self.lvec[1,1],self.lvec[1,2],  self.lvec[2,0],self.lvec[2,1],self.lvec[2,2]   ) ) + comment
        qs = self.qs
        if(not bQs): qs=None
        saveXYZ( self.enames, self.apos, fname, qs=qs, Rs=self.Rs, mode=mode, comment=comment, ignore_es=ignore_es, other_lines=other_lines )

    def toLines(self):
        #lines = []
        #for i,pos in enumerate(self.apos):
        #    lines.append(  "%s %3.5f %3.5f %3.5f\n" %(self.enames[i], pos[0],pos[1],pos[2]) )
        return geomLines( self.apos, self.enames )

    def toXYZ(self, fout, comment="#comment", ignore_es=None, other_lines=None, bHeader=False ):
        writeToXYZ( fout, self.enames, self.apos, qs=self.qs, Rs=self.Rs, bHeader=bHeader, comment=comment, ignore_es=ignore_es, other_lines=other_lines )

    def print(self):
        print( len(self.atypes), len(self.enames), len(self.apos) )
        for i in range(len(self.apos)):
            print( "[%i] %i=%s p(%10.5f,%10.5f,%10.5f)" %( i, self.atypes[i],self.enames[i], self.apos[i,0], self.apos[i,1], self.apos[i,2] ), end =" " )
            if(self.aux_labels is not None): print(self.aux_labels[i], end =" ")
            print("")

    def getValenceElectrons( self ):
        return  np.array( [ elements.ELEMENT_DICT[e][9] for e in self.enames ] )

    def subtractValenceE(self, f0=-1.0, f=+1.0 ):
        self.qs[:] = self.qs[:]*f0 + self.getValenceElectrons()*f       

    def printBonds(self):
        for i in range(len(self.bonds)):
            print( "[%i] (%i,%i) (%s,%s)" %( i, self.bonds[i,0],self.bonds[i,1],  self.enames[self.bonds[i,0]], self.enames[self.bonds[i,1]] ) )

    def findBonds(self, Rcut=3.0, RvdwCut=1.5, RvdWs=None, byRvdW=True ):
        if self.atypes is None:
            self.atypes = [ elements.ELEMENT_DICT[e][0] for e in self.enames ]
        self.bonds, rs = findBondsNP( self.apos, self.atypes, Rcut=Rcut, RvdwCut=RvdwCut, RvdWs=RvdWs, byRvdW=byRvdW )
        return self.bonds, rs

    def findHBonds(self, Rb=1.5, Rh=2.5, angMax=60.0, typs1={"H"}, typs2=neg_types_set, bPrint=False, bHbase=False ):
        return findHBondsNP( self.apos, atypes=self.enames, Rb=Rb, Rh=Rh, angMax=angMax, typs1=typs1, typs2=typs2, bPrint=bPrint,  bHbase=bHbase )

    def findBondsOfAtom(self, ia, bAtom=False ):
        if bAtom: 
            return [ b[1] for b in self.bonds if(b[0]==ia) ] + [ b[0] for b in self.bonds if(b[1]==ia) ] 
        else:
            return [i for i,b in enumerate(self.bonds) if (b[0]==ia) or (b[1]==ia) ]

    def neighs( self, bBond=True ):
        if(self.bonds is None):
            self.findBonds()
        self.ngs = neigh_bonds( len(self.apos), self.bonds )
        return self.ngs

    def find_groups(self):
        if self.ngs is None: self.neighs()
        ngs = self.ngs
        #print( ngs )
        groups = { }
        for inod in range(len(self.apos)):
            if len(ngs[inod]) > 1: groups[inod] = [inod]
        for inod,g in groups.items():
            inod = g[0] 
            g += [ ia for ia in ngs[inod].keys() if ia not in groups ] 
        return groups

    def select_by_ename( self, elist ):
        return [ i for i,e in enumerate(self.enames) if e in elist ]

    def getNeighsOfType( self, selection, typ='N'):
        if self.ngs is None: self.neighs()
        return findNeighsOfType( selection, self.enames, self.ngs, typ=typ ) 

    def select_by_neighType( self, neighs, typ='N', neighTyps={'H':(1,2)} ):
        return findTypeNeigh_( self.enames, neighs, typ=typ, neighTyps=neighTyps )

    # def findTypeNeigh( atoms, neighs=None, typ, neighTyps=[(1,2,2)] ):
    #     if 
    #     def findTypeNeigh( atoms, neighs, typ, neighTyps=[(1,2,2)] ):

    def findAngles(self, select=None, ngs=None, ):
        if ngs is None:
            ngs = self.neighs()
        return findAngles( self.apos, select=select, neighs=ngs )

    def findDihedral( self, select=None, ngs=None, neighTyp={'H'} ):
        if ngs is None:
            ngs = self.neighs()
        return findDihedral( self.apos, self.enames, ngs, select=select, neighTyp=neighTyp ) 

    def findCOG(self, apos, byBox=False ):
        return findCOG( apos, byBox=byBox )
    
    def projectAlongBondDir( self, i0, i1 ):
        return projectAlongBondDir( self.apos, i0, i1 )

    def store_bond_lengths(self):
        bond_lengths = {}
        bonds = self.findBonds()  # Get all bonds in the system
        for bond in bonds:
            ia,ja = bond
            if ia>ja: ia,ja = ja,ia
            length = np.linalg.norm(self.apos[ia]-self.apos[ja])
            bond_lengths[(ia,ja)] = length
        self.bond_legths = bond_lengths
        return bond_lengths

    def restore_bond_length(self, ij, L=None ):
        ia,ja= ij
        d = self.apos[ja] - self.apos[ia]
        Lnow = np.linalg.norm(d)
        if L is None:
            if ia>ja: i,j = ja,ia
            else:     i,j = ia,ja
            L = self.bond_lengths[(i,j)]
        f = L / Lnow
        self.apos[ia] = self.apos[ja] + d * f


    def clonePBC(self,nPBC=(1,1,1) ):
        nx,ny,nz= nPBC
        nxyz=nx*ny*nz
        na = len(self.apos)
        apos   = np.zeros((na*nxyz,3))
        #print( "clonePBC ", na, len(self.atypes) )
        if self.atypes is not None: 
            atypes = np.zeros(na*nxyz,np.int32)
        else:
            atypes = None

        if self.enames is not None: 
            enames = []
        else:
            enames = None

        if self.qs is not None: 
            qs = np.zeros(na*nxyz) 
        else:
            qs = None

        #print( nxyz, na, apos.shape, atypes.shape )
        if( nxyz > 1 ):
            lvec   = np.array([ self.lvec[0,:]*nx,self.lvec[1,:]*ny,self.lvec[2,:]*nz ]) 
            i0=0
            for iz in range(nz):
                for iy in range(ny):
                    for ix in range(nx):
                        shift = self.lvec[0,:]*ix  + self.lvec[1,:]*iy + self.lvec[2,:]*iz
                        apos  [i0:i0+na,:] = self.apos[:,:] + shift[None,:]
                        if atypes is not None: atypes[i0:i0+na  ] = self.atypes
                        if qs     is not None: qs    [i0:i0+na  ] = self.qs    
                        if enames is not None: enames[i0:i0+na  ] = self.enames
                        #if enames is not None: enames += self.enames
                        i0+=na
        else:
            lvec=self.lvec
            apos  [:,:] = self.apos[:,:]
            if atypes is not None: atypes[:] = self.atypes[:]
            if qs     is not None: qs    [:] = self.qs    [:]  
            if enames is not None: enames[:] = self.enames[:]

        return AtomicSystem(apos=apos, atypes=atypes, enames=enames, lvec=lvec, qs=qs ) 

    def symmetrized(self, d=0.1 ):
        # def atoms_symmetrized( atypes, apos, lvec, qs=None, REQs=None, d=0.1):
        atypes, apos, qs, REQs, ws = atoms_symmetrized( self.atypes, self.apos, self.lvec, qs=self.qs, d=d );
        enames = iz2enames( atypes )
        return AtomicSystem( apos=apos, atypes=atypes, enames=enames, lvec=self.lvec.copy(), qs=qs ), ws 

    def selectSubset(self, inds ):
        if self.atypes is not None: 
                atypes = self.atypes[inds]
        else:
            atypes = None

        if self.enames is not None: 
            enames = [ self.enames[i] for i in inds ]
        else:
            enames = None

        if self.qs is not None: 
            qs = self.qs[inds]
        else:
            qs = None

        lvec=self.lvec
        apos  = self.apos[inds,:]

        return AtomicSystem(apos=apos, atypes=atypes, enames=enames, lvec=lvec, qs=qs ) 

    def selectBondedCluster( self, s ):
        na = len(self.apos)
        if self.bonds is None: self.findBonds()
        s     = selectBondedCluster( s, self.bonds )
        ins  = [ i for i in range(na) if (i in s) ]
        outs = [ i for i in range(na) if (i not in s) ] 
        return ins,outs

    def makeRotMat( self, ip1, ip2, _0=1 ):
        fw  = self.apos[ip1[1]-_0]-self.apos[ip1[0]-_0]
        up  = self.apos[ip2[1]-_0]-self.apos[ip2[0]-_0]
        return makeRotMat( fw, up )

    def orient_mat(self, rot, p0=None, bCopy=False ):
        apos=self.apos  
        if(bCopy): apos=apos.copy()
        if p0  is not None: apos[:,:]-=p0[None,:]
        if rot is not None: mulpos( apos, rot )
        return apos

    def orient_vs(self, fw, up, p0=None, trans=None, bCopy=False ):
        if fw is None:
            rot = None
        else:
            rot = makeRotMat( fw, up )
            if trans is not None: rot=rot[trans,:]
        return self.orient_mat( rot, p0, bCopy )

    def orient( self, i0, b1, b2, _0=1, trans=None, bCopy=False ):
        #print( "orient i0 ", i0, " ip1 ", ip1, " ip2 ",ip2 )
        # p0  = self.apos[i0-_0]
        # fw  = self.apos[ip1[1]-_0]-self.apos[ip1[0]-_0]
        # up  = self.apos[ip2[1]-_0]-self.apos[ip2[0]-_0]
        p0, fw, up = makeVectros( self.apos, i0, b1, b2, _0=_0 )
        return self.orient_vs( fw, up, p0, trans=trans, bCopy=bCopy )
    
    def orientPCA(self, perm=None):
        orientPCA(self.apos, perm=perm )

    def shift(self, vec, sel=None ):
        if sel is None: 
            self.apos[:,0] += vec[0]
            self.apos[:,1] += vec[1]
            self.apos[:,2] += vec[2]
        else:
            self.apos[sel,0] += vec[0]
            self.apos[sel,1] += vec[1]
            self.apos[sel,2] += vec[2]

    def rotate_ax(self, ang, ax=(0,1), p0=None ):
        rot = makeRotMatAng( ang, ax=ax )
        if p0  is not None: self.apos[:,:]-=p0[None,:]
        mulpos( self.apos, rot )
        if p0  is not None: self.apos[:,:]+=p0[None,:]

    def delete_atoms(self, lst ):
        st = set(lst)
        if( self.apos   is not None ): self.apos   =  np.delete( self.apos,   lst, axis=0 )
        if( self.atypes is not None ): self.atypes =  np.delete( self.atypes, lst )
        if( self.qs     is not None ): self.qs     =  np.delete( self.qs,     lst )
        if( self.Rs     is not None ): self.Rs     =  np.delete( self.Rs,     lst )
        if( self.enames is not None ): self.enames =  np.delete( self.enames, lst )
        if( self.aux_labels is not None ): self.aux_labels = [ v for i,v in enumerate(self.aux_labels) if i not in st ] 


    def preinitialize_atomic_properties(self):
        """
        Preinitialize per-atom arrays for an AtomicSystem.
        
        This function assumes that the system’s atypes (or enames) have been set.
        It uses the global 'elements.ELEMENTS' (a list of lists) to set default values:
        - qs: set to the element’s default valence electron count (column index 9)
        - Rs: set to the element’s van der Waals radius (column index 7)
        - aux_labels: set to a default label (simply the atom’s index as a string)
        
        Parameters:
        atomicSystem (AtomicSystem): an instance of AtomicSystem.
        
        Raises:
        ValueError: if atomicSystem.atypes is not defined.
        """
        natoms = len(self.apos)
        
        if self.atypes is None:
            raise ValueError("The system does not have atypes defined. "
                            "Please initialize the system’s atypes (or enames) first.")
        
        # Initialize qs if not defined.
        if self.qs is None:
            # Assume atypes is an array of atomic numbers (e.g. 6 for carbon, etc.)
            qs = []
            for z in self.atypes:
                # our ELEMENTS list is zero-based: for atomic number z, use ELEMENTS[z-1]
                qs.append(elements.ELEMENTS[z-1][9])
            self.qs = np.array(qs)
        
        # Initialize Rs if not defined.
        if self.Rs is None:
            # For each atom, use the vdW radius (column index 7)
            Rs = []
            for z in self.atypes:
                Rs.append(elements.ELEMENTS[z-1][7])
            self.Rs = np.array(Rs)
        
        # Initialize aux_labels if not defined.
        if self.aux_labels is None:
            self.aux_labels = [str(i) for i in range(natoms)]
            
        # (If you have other arrays you want to preinitialize, do it here.)
        
        print(f"Pre-initialized atomic properties for {natoms} atoms.")

        
    def check_atomic_properties(atomicSystem):
        """
        Check that the per-atom arrays (qs, Rs, aux_labels) are defined and
        have the correct length. If not, raise an error telling the user
        to run preinitialize_atomic_properties().
        """
        natoms = len(atomicSystem.apos)
        if (atomicSystem.qs is None or len(atomicSystem.qs) != natoms or
            atomicSystem.Rs is None or len(atomicSystem.Rs) != natoms or
            atomicSystem.aux_labels is None or len(atomicSystem.aux_labels) != natoms):
            raise ValueError("Not all per-atom arrays are initialized correctly. "
                            "Please call preinitialize_atomic_properties() on your system.")
                            
                            
    # Example: modify append_atoms() to check rather than auto-initialize
    def append_atoms(self, B, pre="A"):
        # Ensure both systems have been pre-initialized:
        self.check_atomic_properties()
        B.check_atomic_properties()
        
        # Number of atoms in self and in B
        nA = len(self.apos)
        nB = len(B.apos)
        
        self.apos   = np.append(self.apos,   B.apos, axis=0)
        self.atypes = np.append(self.atypes, B.atypes)
        self.qs     = np.append(self.qs,     B.qs)
        self.Rs     = np.append(self.Rs,     B.Rs)
        self.enames = np.append(self.enames, B.enames)
        
        # Extend the aux_labels list:
        self.aux_labels.extend(B.aux_labels)


    def remap( self, lst ):
        dct = {   key:value for (value,key) in enumerate(self.aux_labels) }
        return [ dct.get(key,-1) for key in lst ]


    def attach_group( self, G,  i0, i1, iup,   bond,  up=(0.,0.,1.),  _0=1, pre="A"  ): 
        """
        Attach an end–group (G) to the backbone (self) at a specified bond.
        
        The attachment is done in two steps:
        1. **Internal Orientation of the Group:**  
            The group is reoriented in its own frame by calling:
                G.orient(i0, (i0, i1), iup, _0=_0)
            - *i0*: the index (or indices) for the pivot atom in the group. This atom
                    is moved to the attachment position.
            - *(i0, i1)*: a tuple defining a bond in the group that determines the
                        forward (direction) vector. The forward vector is computed as
                        the difference between the positions of the atom at i1 and i0.
                        The atom corresponding to i1 is then deleted (replaced) in the group.
            - *iup*: a tuple (or list) of two indices that defines the up vector in the group.
                    The up vector is computed (typically as the difference between the
                    positions of the atoms provided) and is used to fix the rotation about
                    the forward axis.
        
        2. **Alignment to the Backbone:**  
            The backbone provides the attachment bond (given by `bond`) and a backbone
            up vector (given by `up`). A rotation matrix is computed with:
                rot = rotmat_from_points(self.apos, ifw=bond, up=up, _0=_0)
            This matrix aligns the backbone’s forward vector (computed from the bond) with
            the group’s forward vector. The group is then rotated by this matrix (via
                G.orient_mat(rot)
            ) and translated so that the pivot atom of the group coincides with the backbone’s
            attachment position.
        
        Parameters:
        G      : AtomicSystem
                The end–group to attach. It must have its atoms pre‐oriented as per the
                expected coordinate system.
        i0     : int or iterable
                The index (or indices) of the pivot atom in G (1-based indexing if _0=1).
        i1     : int
                The index (1-based) of the atom in G used to define the forward vector.
                This atom will be removed after orientation.
        iup    : tuple (i_up0, i_up1)
                A pair of indices (1-based) in G whose difference defines the up vector.
        bond   : tuple (i_backbone1, i_backbone2)
                A pair of atom indices (1-based) in the backbone that define the bond where
                the end–group is attached. The forward vector on the backbone is computed as
                the vector from i_backbone1 to i_backbone2.
        up     : 3-tuple or array, optional (default=(0.,0.,1.))
                The up vector for the backbone. This is used to fix the rotation about the
                forward axis.
        _0     : int, optional (default=1)
                An offset to account for whether the provided indices are 0-based or 1-based.
        pre    : str, optional (default="A")
                A prefix for labeling the atoms that come from the group.
        
        After executing, the group G is reoriented, rotated, and translated so that its
        pivot atom is placed at the backbone’s attachment site. The atom used for forward
        definition (i1) is deleted, and the group’s atoms are appended to the backbone.
        """
        up  = np.array( up )
        rot = rotmat_from_points( self.apos, ifw=bond, up=up, _0=1 );   
        rot = rot.transpose()
        p0  = self.apos[bond[0]-_0]
        
        if( G.aux_labels is None ): G.aux_labels = [ pre+str(i) for  i in range(len(G.apos)) ]

        G.orient( i0,(i0,i1),iup, _0=_0 )
        G.orient_mat( rot ); 
        G.apos[:,:]+=p0[None,:]
        G.delete_atoms( [i1-_0] )

        self.append_atoms( G, pre=pre )

    #def orient_vs( p0, fw, up, apos, trans=None, bool bCopy ):
    #def orient( i0, ip1, ip2, apos, _0=1, trans=None, bCopy=True ):