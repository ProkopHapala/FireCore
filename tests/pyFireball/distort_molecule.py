import sys
import os
import numpy as np

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall.AtomicSystem import AtomicSystem


# ======== Body

path = './relaxed_mols/'
#mol = au.AtomicSystem( path+'CH2NH.xyz' )
#mol = au.AtomicSystem( path+'NH2NH2.xyz' )
#mol = au.AtomicSystem( path+'NH2OH.xyz' )
#mol = au.AtomicSystem( path+'CH3NH2.xyz' )



def move_group( mol, gi, gj, l=0.1 ):
    #i,j=ij
    #print( "move_group : groups", groups )
    #gi=groups[i]
    #gj=groups[j]
    #print( "move_group i,j: ", i,j )
    ia = gi[0]
    ja = gj[0]
    #print( "move_group ia,ja: ", ia,ja )
    d  = mol.apos[ja] - mol.apos[ia]
    mol.apos[gj,:] += d*l
    #mol.apos[gi,:] -= d*l

def scan_group_distance( mol, groups, l0=-0.4, dl=0.1 ):
    fname = path+"scan_"+name
    print( "scan_group_distance: ", fname )
    gi,gj = groups[0],groups[1]
    move_group( mol, gi, gj, l0 )
    mol.saveXYZ( fname, mode='w' )
    for i in range(10):
        move_group( mol, gi, gj, dl )
        mol.saveXYZ( fname, mode='a' )

def scan_angle( mol, angs, ia=0, j=1, jbs=[2] ):
    fname = path+"angscan_"+name
    #r1 = mol.apos[j1] - mol.apos[i]
    #r2 = mol.apos[j2] - mol.apos[i]
    ax = au.normalize(mol.apos[j] - mol.apos[ia]) 
    vs = []
    rs = []
    for j in jbs:
        v = mol.apos[j] - mol.apos[ia]
        r = np.linalg.norm(v)
        vs.append(v/r)
        rs.append(r)
    mol.saveXYZ( fname, mode='w' )
    for a in angs:
        ca = np.cos(a)
        sa = np.sin(a)
        for jj,j in enumerate(jbs):
            r = rs[jj]
            mol.apos[j,:] = mol.apos[ia,:] + ax*(r*ca) + vs[jj]*(r*sa)
        mol.saveXYZ( fname, mode='a' )


# find all molecules in path
#molecules = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
#print(molecules)

molecules = ["H2O.xyz"]

#mol = AtomicSystem( path+molecules[0] )
for name in molecules:
    mol    = AtomicSystem( path+name )
    
    # scan group distance
    # groups = mol.find_groups()
    # groups = list(groups.values())
    # groups.sort( key=lambda g: g[0] )
    # print(  name + ".groups ", groups )
    # if len(groups)<2: continue
    # scan_group_distance( mol, groups, l0=-0.4, dl=0.1 )

    # scan angle
    angs = np.linspace(np.pi/4, np.pi*7./8., 10)
    scan_angle( mol, angs, ia=0, j=1, jbs=[2] )