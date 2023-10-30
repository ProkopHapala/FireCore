import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall             import atomicUtils  as au
#from pyBall             import plotUtils   as plu

#======== Functions

def read_Endgroup_metadata(fname):
    fin = open(fname, 'r')
    ws = fin.readline().split(); ifw=(int(ws[0]),int(ws[1]))
    ws = fin.readline().split(); ilf=(int(ws[0]),int(ws[1]))
    n = int(fin.readline().split()[0])
    sites = []
    for i in range(n):
        sites.append( int(fin.readline().split()[1]) )
    #print( ifw, ilf, sites )
    return ifw,ilf, sites

def attachPair( B, name1, name2, group_dict, _0=1, bSaveDebugG=False, amargin=5.0 ):
    
    BB = B.clonePBC()
    #G1 = AtomicSystem(fname="endgroups/"+name1+".xyz"  )
    #G2 = AtomicSystem(fname="endgroups/"+name2+".xyz"  )
    inds1, Hs1, G1 = group_dict[name1]
    inds2, Hs2, G2 = group_dict[name2]

    G1 = G1.clonePBC()
    G2 = G2.clonePBC()
    
    if bSaveDebugG:
        G1.enames[  inds1[1]-_0 ]='Cl'
        G2.enames[  inds2[1]-_0 ]='Cl'
        G1.saveXYZ( "G1."+name1+".xyz" )
        G2.saveXYZ( "G2."+name1+".xyz" )
    
    #BB.attach_group( G1, inds1[0], inds1[1], inds1[2], (1 ,2), up=(0.,0.,1.), _0=1  )
    #BB.attach_group( G2, inds2[0], inds2[1], inds2[2], (17,9), up=(0.,0.,1.), _0=1  )
    #BB.delete_atoms( [1-1,17-1] )
    BB.attach_group( G1, inds1[0], inds1[1], inds1[2], (18,8), up=(0.,0.,1.), _0=1 , pre="X" )
    BB.attach_group( G2, inds2[0], inds2[1], inds2[2], (17,6), up=(0.,0.,1.), _0=1 , pre="Y" )
    BB.delete_atoms( [17-_0, 18-_0 ] )
    
    inds1 = BB.remap( [ "X"+str(i-1) for i in Hs1 ] )
    inds2 = BB.remap( [ "Y"+str(i-1) for i in Hs2 ] )

    dir1 = BB.apos[ inds1[-1] ] - BB.apos[ inds1[-0] ]; 
    if(dir1[2]<0.): 
        inds1=inds1[::-1] 
    dir2 = BB.apos[ inds2[-1] ] - BB.apos[ inds2[-0] ]; 
    if(dir2[2]<0.): 
        inds2=inds2[::-1]

    avec = BB.apos[ inds1[0] ] - BB.apos[ inds2[0] ]
    r = np.sqrt(np.dot(avec,avec))
    print( r, avec )
    avec*=( (r+2.0)/(r) )
    avec[2] = 0.0
    avec[0] += amargin
    
    BB.lvec[0] = avec
    
    return BB, inds1, inds2

def findLastHydrogen( atoms, ifw, ilf ):
    rot = atoms.makeRotMat( ifw, ilf, _0=1 )
    dir = rot[2,:]
    rs  = np.dot( atoms.apos, dir[:,None] )  #;print(prjs)
    rs[atoms.atypes != 1 ] = 1000
    iH = np.argmin( rs )
    return iH

def nameToTyp( s ):
    return s.split('-')[0].replace("O", "e" ).replace("N", "e" )

def saveMolGUIscript(name, inds, cls='HHH', path="./", amargin=5.0, bBoncConstr=True, bAngConstr=True ):
    fsh = open( path+name+".sh", 'w')
    fsh.write(f"./MolGUIapp -x {name} -b {name}.hbonds -dlvec {-amargin},0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0\n")
    fsh.close()
    fhb = open( path+name+".hbonds", 'w')

    print( "inds[0], inds[1], inds[2] ", inds[0], inds[1], inds[2], len(inds[0]), len(inds[1]), len(inds[2]) )
    nmin = min( len(inds[0]), len(inds[1]) )
    # write constrains for /home/prokop/git/FireCore/cpp/common/molecular/constrains.h
    if bBoncConstr:
        # sscanf( line, "b %i %i      %lf %lf       %lf %lf   %lf    %lf %lf %lf",   &cons.ias.a,&cons.ias.b,  &cons.ls.a,&cons.ls.b,  &cons.ks.a,&cons.ks.b,  &cons.flim, &cons.shift.a,&cons.shift.b,&cons.shift.c );
        for i in range(nmin): 
            fhb.write( "b %4i %4i 1.7 1.5    5.0 0.0  10.0   1.0 0.0 0.0 \n" %(inds[0][i],inds[1][i]) ) 
    if bAngConstr:
        # sscanf( line, "g %i %i %i   %lf %lf %lf   %i %i %i  %i %i %i",    &cons.ias.a,&cons.ias.b,&cons.ias.c,   &ang,   &cons.k,   &cons.flim,   &cons.acell.a,&cons.acell.b,&cons.acell.c,   &cons.bcell.a,&cons.bcell.b,&cons.bcell.c );
        #for i in range(len(inds[0])): fhb.write( "g %4i %4i %4i  180.0 1.5 10.0  0 0 0    -1 0 0  \n" %(inds[0][i],inds[1][i],inds[2][i]) )  
        for i in range(nmin):
            if cls[i]=='H':
                fhb.write( "g %4i %4i %4i  180.0 1.5 10.0  1 0 0    0 0 0  \n" %(inds[0][i],inds[1][i],inds[2][i]) ) 
            else: 
                fhb.write( "g %4i %4i %4i  180.0 1.5 10.0  -1 0 0    0 0 0  \n" %(inds[1][i],inds[0][i],inds[2][i]) )
    fhb.close()

def load_groups( names, folder="/home/prokop/Desktop/CARBSIS/Paolo/endgroups/", dir_meta=None, dir_relax=None ):
    if dir_meta is None:  dir_meta  = folder+"endgroups/"
    if dir_relax is None: dir_relax = folder+"endgroups_relaxed/mols/"
    typs={}
    group_dict = {}
    for name in names:
        typ = nameToTyp( name )
        if typ in typs: 
            typs[typ].append(name) 
        else: 
            typs[typ] = [name] 
        
        ifw,ilf, sites = read_Endgroup_metadata(dir_meta+name+".txt" )
        #atoms = au.AtomicSystem( dir_meta+name+".xyz"  )
        atoms  = au.AtomicSystem( dir_relax+name+"/final.xyz"  )
        atoms.subtractValenceE()

        #print( name+".qs:   ", atoms.qs );
        #print( name+".Qtot: ", atoms.qs.sum() );
        
        iH = findLastHydrogen( atoms, ifw, ilf  )
        atoms.findBonds()
        iC = atoms.findBondsOfAtom( iH, bAtom=True )[0]
        #print( iC, len(atoms.apos) )

        # --- Debug - repleace H-bond atoms with Cl and Si and save
        # atoms.enames[iH] = 'Cl'
        # atoms.enames[iC] = 'Si'
        # atoms.saveXYZ("out/"+name+".xyz")

        group_dict[name] = ( ( iC+1, iH+1, ilf ), sites, atoms )
        #       name                   attachment         H-Bonds           
        #                                C/N  H    Up
        ## ( "penta_hb3_acceptor2",  ( ( 10, 11, (14, 8) ),  [ 8, 6,14 ] ) ),
    return group_dict, typs

