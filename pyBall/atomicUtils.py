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

def getRvdWs( atypes, eparams=elements.ELEMENTS ):
    #print( eparams[ 6 ][7], eparams[ 6 ] )
    return [ eparams[ ei ][7] for ei in atypes ]

def getRvdWsNP( atypes, eparams=elements.ELEMENTS ):
    return np.array( getRvdWs( atypes, eparams ) ) 

def findBondsNP( apos, atypes=None, Rcut=3.0, RvdwCut=0.5, RvdWs=None, byRvdW=True ):
    bonds  = []
    iatoms = np.arange( len(apos), dtype=int )
    Rcut2  = Rcut*Rcut
    if byRvdW:
        if  RvdWs is None:
            RvdWs = getRvdWsNP( atypes, eparams=elements.ELEMENTS )
            #print( RvdWs )
    else:
        RvdWs = np.ones(len(apos))*Rcut
    for i,pi in enumerate(apos):
        j0=i+1
        rs   = np.sqrt( np.sum( (apos[j0:,:] - pi[None,:] )**2, axis=1 ) )
        mask = rs[:] < ( RvdWs[j0:]+RvdWs[i] )*RvdwCut
        bonds += [ (i,j) for j in iatoms[j0:][mask] ]
    return np.array( bonds, dtype=np.int32 )

def findHBondsNP( apos, atypes=None, Rb=1.5, Rh=2.2, angMax=30.0, typs1={"H"}, typs2=neg_types_set, bPrint=False ):
    bonds  = []
    rbs    = []
    iatoms = np.arange( len(apos), dtype=int )
    cos_min = np.cos( angMax*np.pi/180.0 )
    print( "cos_min ", cos_min )
    for i,pi in enumerate(apos):
        if ( atypes[ i ] not in typs1) : continue
        ds   = apos[:,:] - pi[None,:]
        rs   = np.sqrt( np.sum( ds**2, axis=1 ) )
        rs[i]=100.0
        jmin = np.argmin(rs)   # nearest neighbor i.e. bond
        cs   = np.dot( ds, ds[jmin] ) / ( rs*rs[jmin] )
        mask = np.logical_and( cs<-cos_min, rs<Rh )
        dbonds = [ (i,j) for j in iatoms[:][mask] if atypes[j] in typs2 ]
        bonds += dbonds
        rbs   += [ rs[b[1]] for b in dbonds ]

    return np.array( bonds, dtype=np.int32 ), np.array( rbs )

def neighs( natoms, bonds ):
    neighs = [{} for i in range(natoms) ]
    for ib, b in enumerate(bonds):
        i = b[0]; j = b[1]; 
        neighs[i][j] = ib
        neighs[j][i] = ib
    return neighs

def findTypeNeigh( atoms, neighs, typ, neighTyps=[(1,2,2)] ):
    typ_mask = ( atoms[:,0] == typ )
    satoms   = atoms[typ_mask]
    iatoms   = np.arange(len(atoms),dtype=int)[typ_mask]
    selected = []
    for i,atom in enumerate(satoms):
        iatom = iatoms[i]
        #for jatom in neighs[ iatom ]:
        #    jtyp = atoms[jatom,0]
        count = {}
        for jatom in neighs[ iatom ]:
            jtyp = atoms[jatom,0]
            count[jtyp] = count.get(jtyp, 0) + 1
        for jtyp, (nmin,nmax) in list(neighTyps.items()):
            n = count.get(jtyp,0)
            if( (n>=nmin)and(n<=nmax) ):
                selected.append( iatom )
    return selected
    
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

def makeRotMatAng( ang, ax1=0, ax2=1 ):
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

def writeToXYZ( fout, es, xyzs, qs=None, Rs=None, comment="#comment", bHeader=True, ignore_es=None ):
    if(bHeader):
        na=len(xyzs)
        if ignore_es is not None:
            mask = [ (  e not in ignore_es) for e in es ]
            na = sum(mask)
        else:
            mask = [True]*na
        fout.write("%i\n"  %na )
        fout.write(comment+"\n")
    if   (Rs is not None):
        for i,xyz in enumerate( xyzs ):
            if mask[i]: fout.write("%s %f %f %f %f %f \n"  %( es[i], xyz[0], xyz[1], xyz[2], qs[i], Rs[i] ) )
    elif (qs is not None):
        for i,xyz in enumerate( xyzs ):
            if mask[i]: fout.write("%s %f %f %f %f\n"  %( es[i], xyz[0], xyz[1], xyz[2], qs[i] ) )
    else:
        for i,xyz in enumerate( xyzs ):
            if mask[i]: fout.write("%s %f %f %f\n"  %( es[i], xyz[0], xyz[1], xyz[2] ) )

def saveXYZ( es, xyzs, fname, qs=None, Rs=None, mode="w", comment="#comment", ignore_es=None ):
    fout = open(fname, mode )
    writeToXYZ( fout, es, xyzs, qs, Rs=Rs, comment=comment, ignore_es=ignore_es )
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

def loadAtomsNP(fname=None, fin=None, bReadN=False, nmax=10000 ):
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
        wds = line.split()
        try:
            xyzs.append( ( float(wds[1]), float(wds[2]), float(wds[3]) ) )
            try:
                iz    = int(wds[0]) 
                Zs    .append(iz)
                enames.append( elements.ELEMENTS[iz] )
            except:
                ename = wds[0]
                enames.append( ename )
                Zs    .append( elements.ELEMENT_DICT[ename][0] )
            try:
                q = float(wds[4])
            except:
                q = 0
            qs.append(q)
            ia+=1
        except:
            #print("cannot interpet line: ", line)
            if bReadN and (ia==0):
                try:
                    nmax=int(wds[0])
                    print("nmax: ", nmax)
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
    return xyzs,Zs,enames,qs

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
                        S = Atoms( lvec=lvec, enames=es, apos=apos)
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
            n=int(l.split()[0])
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

# ========================== Class Geom

class AtomicSystem():

    def __init__(self,fname=None, apos=None, atypes=None, enames=None, lvec=None, qs=None, Rs=None, bonds=None ) -> None:
        self.apos    = apos
        self.atypes  = atypes
        self.enames  = enames
        self.qs      = qs
        self.Rs      = Rs
        self.bonds   = bonds
        self.lvec    = lvec
        self.aux_labels = None
        if fname is not None:
            if( '.mol' == fname.split('.')[0] ):
                self.apos,self.atypes,self.enames,self.qs,self.bonds = loadMol(fname=fname, bReadN=True )
            else:
                self.apos,self.atypes,self.enames,self.qs = loadAtomsNP(fname=fname , bReadN=True )

    def saveXYZ(self, fname, mode="w", blvec=True, comment="", ignore_es=None ):
        if blvec and (self.lvec is not None):
            print( self.lvec )
            comment= ( "lvs %6.3f %6.3f %6.3f   %6.3f %6.3f %6.3f   %6.3f %6.3f %6.3f" %(self.lvec[0,0],self.lvec[0,1],self.lvec[0,2],  self.lvec[1,0],self.lvec[1,1],self.lvec[1,2],  self.lvec[2,0],self.lvec[2,1],self.lvec[2,2]   ) ) + comment
        saveXYZ( self.enames, self.apos, fname, qs=self.qs, Rs=self.Rs, mode=mode, comment=comment, ignore_es=ignore_es )
    
    def toXYZ(self, fout ):
        writeToXYZ( fout, self.enames, self.apos, qs=self.qs, Rs=self.Rs, bHeader=False )

    def print(self):
        #print( len(self.atypes), len(self.enames), len(self.apos) )
        for i in range(len(self.apos)):
            print( "[%i] %i=%s p(%10.5f,%10.5f,%10.5f)" %( i, self.atypes[i],self.enames[i], self.apos[i,0], self.apos[i,1], self.apos[i,2] ), end =" " )
            if(self.aux_labels is not None): print(self.aux_labels[i], end =" ")
            print("")
            

    def printBonds(self):
        for i in range(len(self.bonds)):
            print( "[%i] (%i,%i) (%s,%s)" %( i, self.bonds[i,0],self.bonds[i,1],  self.enames[self.bonds[i,0]], self.enames[self.bonds[i,1]] ) )

    def findBonds(self, Rcut=3.0, RvdwCut=0.7, RvdWs=None, byRvdW=True ):
        if self.atypes is None:
            self.atypes = [ elements.ELEMENT_DICT[e][0] for e in self.enames ]
        self.bonds = findBondsNP( self.apos, self.atypes, Rcut=Rcut, RvdwCut=RvdwCut, RvdWs=RvdWs, byRvdW=byRvdW )

    def findHBonds(self, Rb=1.5, Rh=2.2, angMax=30.0, typs1={"H"}, typs2=neg_types_set, bPrint=False ):
        return findHBondsNP( self.apos, atypes=self.enames, Rb=Rb, Rh=Rh, angMax=angMax, typs1=typs1, typs2=typs2, bPrint=True )

    def findCOG(self, apos, byBox=False ):
        return findCOG( apos, byBox=byBox )
        
    def clonePBC(self,nPBC=(1,1,1) ):
        nx,ny,nz= nPBC
        nxyz=nx*ny*nz
        na = len(self.apos)
        apos   = np.zeros((na*nxyz,3))

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
                        if enames is not None: enames[i0:i0+na  ] += self.enames
                        #if enames is not None: enames += self.enames
                        i0+=na
        else:
            lvec=self.lvec
            apos  [:,:] = self.apos[:,:]
            if atypes is not None: atypes[:] = self.atypes[:]
            if qs     is not None: qs    [:] = self.qs    [:]  
            if enames is not None: enames[:] = self.enames[:]

        return AtomicSystem(apos=apos, atypes=atypes, enames=enames, lvec=lvec, qs=qs ) 

    def orient_mat(self, rot, p0=None, bCopy=False ):
        apos=self.apos  
        if(bCopy): apos=apos.copy()
        if p0 is not None: apos[:,:]-=p0[None,:]
        mulpos( apos, rot )
        return apos

    def orient_vs(self, fw, up, p0=None, trans=None, bCopy=False ):
        rot = makeRotMat( fw, up )
        if trans is not None: rot=rot[trans,:]
        return self.orient_mat( rot, p0, bCopy )

    def orient( self, i0, ip1, ip2, _0=1, trans=None, bCopy=False ):
        #print( "orient i0 ", i0, " ip1 ", ip1, " ip2 ",ip2 )
        p0  = self.apos[i0-_0]
        fw  = self.apos[ip1[1]-_0]-self.apos[ip1[0]-_0]
        up  = self.apos[ip2[1]-_0]-self.apos[ip2[0]-_0]
        return self.orient_vs( fw, up, p0, trans=trans, bCopy=bCopy )

    def delete_atoms(self, lst ):
        st = set(lst)
        if( self.apos   is not None ): self.apos   =  np.delete( self.apos,   lst, axis=0 )
        if( self.atypes is not None ): self.atypes =  np.delete( self.atypes, lst )
        if( self.qs     is not None ): self.qs     =  np.delete( self.qs,     lst )
        if( self.Rs     is not None ): self.Rs     =  np.delete( self.Rs,     lst )
        if( self.enames is not None ): self.enames =  np.delete( self.enames, lst )
        if( self.aux_labels is not None ): self.aux_labels = [ v for i,v in enumerate(self.aux_labels) if i not in st ] 

    def append_atoms(self, B, pre="A" ):
        if( self.aux_labels is None ):
            print( 'self.aux_labels is None', pre ) 
            self.aux_labels = [ str(i) for  i in range(len(self.apos)) ]
        
        #if( B.auxl is None ): B.auxl = [ pre+str(i) for  i in range(len(B.apos)) ]
        self.aux_labels += B.aux_labels

        if( self.apos   is not None ): self.apos   =  np.append( self.apos,   B.apos, axis=0 )
        if( self.atypes is not None ): self.atypes =  np.append( self.atypes, B.atypes )
        if( self.qs     is not None ): self.qs     =  np.append( self.qs,     B.qs )
        if( self.Rs     is not None ): self.Rs     =  np.append( self.Rs,     B.Rs )
        if( self.enames is not None ): self.enames =  np.append( self.enames, B.enames )
        #print( type(self.enames),   type(B.enames),   )
        #print( "self.enames ", self.enames, "B.enames ", B.enames )
        #if( self.enames is not None ): self.enames += B.enames
       
        #print(auxl)
        #self.aux_labels += auxl
        #print( self.aux_labels )
        #print( "len( self.aux_labels) ", len( self.aux_labels), "len( self.apos) ", len( self.apos)  )


    def remap( self, lst ):
        dct = {   key:value for (value,key) in enumerate(self.aux_labels) }
        return [ dct.get(key,-1) for key in lst ]


    def attach_group( self, G,  i0, i1, iup,   bond,  up=(0.,0.,1.),  _0=1, pre="A"  ): 
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