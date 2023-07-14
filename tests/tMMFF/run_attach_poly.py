import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall  import atomicUtils as au
from pyBall.atomicUtils import AtomicSystem
from pyBall  import plotUtils   as plu


from pyBall import MMFF as mmff

#======== Functions

def read_Endgroup_metadata(fname):
    fin = open(fname, 'r')
    ws = fin.readline().split(); ifw=(int(ws[0]),int(ws[1]))
    ws = fin.readline().split(); ilf=(int(ws[0]),int(ws[1]))
    n = int(fin.readline().split()[0])
    sites = []
    for i in range(n):
        sites.append( int(fin.readline().split()[1]) )
    print( ifw, ilf, sites )
    return ifw,ilf, sites

def plotGroup( G, inds, bFlip=False ):
    iis,his = inds
    his     = np.array(his)-1
    G.orient( iis[0], iis[2], (iis[0],iis[1]), trans=(1,2,0) )  
    if(bFlip): G.apos[:,0]*=-1
    plu.plotSystem( G, sz=1000. );  
    ps=G.apos[his];                 plt.scatter(ps[:,0],ps[:,1],color='k',zorder=5);  
    ps=G.apos[[iis[0]-1,iis[1]-1]]; plt.scatter(ps[:,0],ps[:,1],color=['b','g'],zorder=5);  
    #plt.title(name1)

def attachPair( name1, name2, group_dict, _0=1, bSaveDebugG=False ):
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
    #print( inds )
    comment = " Hbonds:X"+str(inds1)+"Y"+str(inds2)

    #BB.print()
    BB.saveXYZ( "BB."+name1+"."+name2+".xyz", comment=comment )

def findLastHydrogen( atoms, ifw, ilf ):
    rot = atoms.makeRotMat( ifw, ilf, _0=1 )
    dir = rot[2,:]
    rs  = np.dot( atoms.apos, dir[:,None] )  #;print(prjs)
    rs[atoms.atypes != 1 ] = 1000
    iH = np.argmin( rs )
    return iH

def nameToTyp( s ):
    return s.split('-')[0].replace("O", "e" ).replace("N", "e" )

#def selectPartners():

#======== Body

folder="/home/prokop/Desktop/CARBSIS/Paolo/endgroups/"
dir_meta  = folder+"endgroups/"
dir_relax = folder+"endgroups_relaxed/mols/"
#print(dir_relax)
#print(dir_meta)

names = os.listdir( dir_relax )
#print( names )



typs={}
group_dict = {}
for name in names:
    typ = nameToTyp( name )
    if typ in typs: 
        typs[typ].append(name) 
    else: 
        typs[typ] = [name] 
    

    ifw,ilf, sites = read_Endgroup_metadata(dir_meta+name+".txt" )
    atoms = AtomicSystem( dir_meta+name+".xyz"  )
    
    iH = findLastHydrogen( atoms, ifw, ilf  )
    atoms.findBonds()
    iC = atoms.findBondsOfAtom( iH, bAtom=True )[0]
    #print( iC, len(atoms.apos) )

    # atoms.enames[iH] = 'Cl'
    # atoms.enames[iC] = 'Si'
    # atoms.saveXYZ("out/"+name+".xyz")

    group_dict[name] = ( ( iC+1, iH+1, ilf ), sites, atoms )
    #       name                   attachment         H-Bonds           
    #                                C/N  H    Up
    ## ( "penta_hb3_acceptor2",  ( ( 10, 11, (14, 8) ),  [ 8, 6,14 ] ) ),

for typ in typs: print( typ, typs[typ]  )

pairs = [
("HNH-h","OHO-h_1")
]

pairTypes = [
("HeH","eHe"),
#("HHH","eee"),
#("HHe","Hee"),
#("HH","ee"),
]

B = AtomicSystem(fname='backbone.xyz' )
B.lvec = np.array( [[25.,0.,0.],[0.,5.,0.],[0.,0.,20.0]  ] )
#for pair in pairs:
#    name1, name2 = pair
#    attachPair( name1, name2, group_dict )

for pairTyp in pairTypes:
    names1 = typs[pairTyp[0]]
    names2 = typs[pairTyp[1]]
    for name1 in names1:
        for name2 in names2:
            print( name1, name2 )
            attachPair( name1, name2, group_dict )







'''
mmff.setVerbosity( verbosity=1, idebug=0 )
mmff.init( xyz_name="data/pyridine", surf_name="data/NaCl_1x1_L2", bMMFF=False  )              # without MMFF
mmff.getBuffs()
mmff.eval()
print( "FORCES:\n mmff.fapos:\n ", mmff.fapos )
mmff.plot(bForce=True, Fscale=10.0 )
plt.show()
exit(0)
'''




