import sys
import os

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import FireCore as fc

# ======== Body

path = './relaxed_mols/'
#mol = au.AtomicSystem( path+'CH2NH.xyz' )
#mol = au.AtomicSystem( path+'NH2NH2.xyz' )
#mol = au.AtomicSystem( path+'NH2OH.xyz' )
#mol = au.AtomicSystem( path+'CH3NH2.xyz' )

# find all molecules in path
molecules = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
print(molecules)

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
    gi,gj = groups[0],groups[1]
    move_group( mol, gi, gj, l0 )
    mol.saveXYZ( fname, mode='w' )
    for i in range(10):
        move_group( mol, gi, gj, dl )
        mol.saveXYZ( fname, mode='a' )

#mol = au.AtomicSystem( path+molecules[0] )
for name in molecules:
    mol    = au.AtomicSystem( path+name )
    groups = mol.find_groups()
    groups = list(groups.values())
    groups.sort( key=lambda g: g[0] )
    print(  name + ".groups ", groups )
    if len(groups)<2: continue
    scan_group_distance( mol, groups, l0=-0.4, dl=0.1 )