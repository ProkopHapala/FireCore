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
    #print( ifw, ilf, sites )
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

def attachPair( name1, name2, group_dict, _0=1, bSaveDebugG=False, amargin=5.0 ):
    
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

def saveMolGUIscript(name, inds, path="./", amargin=5.0 ):
    fsh = open( path+name+".sh", 'w')
    fsh.write(f"../../cpp/Build/apps/MolecularEditor/MolGUIapp -x {name} -b {name}.hbonds -dlvec {-amargin},0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0\n")
    fsh.close()
    fhb = open( path+name+".hbonds", 'w')
    for i in range(len(inds[0])):
        fhb.write( "%i %i 1.7 1.5    5.0 0.0     1.0 0.0 0.0   10.0 \n" %(inds[0][i],inds[1][i]) ) 
    fhb.close()




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
    #atoms = AtomicSystem( dir_meta+name+".xyz"  )
    atoms = AtomicSystem( dir_relax+name+"/final.xyz"  )
    atoms.subtractValenceE()

    #print( name+".qs: ", atoms.qs );
    print( name+".Qtot: ", atoms.qs.sum() );
    
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
#("HeH","eHe"),
#("HHH","eee"),
#("HHe","Hee"),
("HHH","eee"),
#("HH","ee"),
]

B = AtomicSystem(fname='backbone.xyz' )
B.lvec = np.array( [[25.,0.,0.],[0.,5.,0.],[0.,0.,20.0]  ] )
#for pair in pairs:
#    name1, name2 = pair
#    attachPair( name1, name2, group_dict )

amargin = 5.0
for pairTyp in pairTypes:
    names1 = typs[pairTyp[0]]
    names2 = typs[pairTyp[1]]

    odir = "out_"+pairTyp[0]+"_"+pairTyp[1]+"/"
    try:
        os.mkdir( odir )
        os.mkdir( odir+"/2x2/" )
    except:
        pass
    for name1 in names1:
        for name2 in names2:
            #print( name1, name2 )
            BB, inds1, inds2 = attachPair( name1, name2, group_dict, amargin=amargin )
            inds1b = BB.getNeighsOfType( inds1, typ='N')
            inds1b = [ a[0] for a in inds1b ] 
            print( "inds1, inds2, inds1b ", inds1, inds2, inds1b )

            comment = " Hbonds={'X':"+str(inds1)+",'Y':"+str(inds2)+"}"
            #BB.enames[ inds1 ] = 'As'
            #BB.enames[ inds2 ] = 'P'

            #BB.enames[ inds1[0] ] = 'F';BB.enames[ inds1[1] ] = 'Cl';BB.enames[ inds1[2] ] = 'Br';
            #BB.enames[ inds2[0] ] = 'Ne';BB.enames[ inds2[1] ] = 'Ar'; BB.enames[ inds2[2] ] = 'Kr';

            name = "BB."+name1+"."+name2

            print( name+".Qtot: ", BB.qs.sum() )

            BB.saveXYZ( odir+name+".xyz", comment=comment )

            BB_ = BB.clonePBC( (2,2,1) )
            BB_.saveXYZ( odir+"/2x2/"+name+"_2x2.xyz", comment=comment )

            saveMolGUIscript( name, (inds1,inds2), path="./out/", amargin=amargin-3.0 )

            
            fig = plt.figure(figsize=(16,4))
            axes=(0,1)
            plt.subplot(1,2,1)
            plu.plotSystem( BB, bBonds=True, colors=None, sizes=None, extent=None, sz=50., RvdwCut=0.5, axes=axes, bLabels=True, labels=None, _0=1, HBs=None, bHBlabels=True, bBLabels=False )
            plu.plotAtoms( BB.apos, colors='#FF0000', marker='o', axes=axes, selection = inds1 )
            plu.plotAtoms( BB.apos, colors='#FFFF00', marker='o', axes=axes, selection = inds1b )
            plu.plotAtoms( BB.apos, colors='#0000FF', marker='o', axes=axes, selection = inds2 )
            axes=(0,2)
            plt.subplot(1,2,2)
            plu.plotSystem( BB, bBonds=True, colors=None, sizes=None, extent=None, sz=50., RvdwCut=0.5, axes=axes, bLabels=True, labels=None, _0=1, HBs=None, bHBlabels=True, bBLabels=False )
            plu.plotAtoms( BB.apos, colors='#FF0000', marker='o', axes=axes, selection = inds1 )
            plu.plotAtoms( BB.apos, colors='#FFFF00', marker='o', axes=axes, selection = inds1b )
            plu.plotAtoms( BB.apos, colors='#0000FF', marker='o', axes=axes, selection = inds2 )

            plt.tight_layout()
            plt.savefig( odir+name+".png", bbox_inches='tight' )
            plt.close(fig)




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




