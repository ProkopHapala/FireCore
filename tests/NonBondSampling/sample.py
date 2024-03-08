import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append('../../')
from pyBall import atomicUtils as au
from pyBall import plotUtils   as plu

# ================== Settings ==================

Rvdws = {
'H'  : 1.40,
'C'  : 1.80,
'N'  : 1.70,
'O'  : 1.60,
'F'  : 1.60,
'Cl' : 1.90,
}

mols=[
["C2H2",    0.5  ],
["C2H4",    0.5  ],
["CH2NH",   1.00 ],
["CH2O",    0.50 ],
["H2O",     0.50 ],
["HCl",     0.50 ],
["HCN",     0.50 ],
["HCONH2",  1.00 ],
["HCOOH",   1.00 ], 
["HF",      0.50 ],
#["NH3",     0.50 ],
["NH3",     1.00 ],
]

# ================== Functions ==================

def distFromMol( p, mol, Rs ):
    ds = np.sqrt( np.sum( (mol.apos[:,:] - p[None,:])**2, axis=1 ) ) -  Rs
    #print ("ds ",ds)
    return ds.min()

def rayMarch( mol, p0, dir, Rs=None, tol=1e-3 ):
    if Rs is None:
        Rs = np.array( [ Rvdws[ mol.atypes[i] ] for i in range( len(mol.atypes) ) ] )
    p = p0
    for i in range(100):
        r = distFromMol(p, mol, Rs)
        if r < tol:
            break
        p += dir*r
    return p

def makeDirs( phis, z0=0.0 ):
    dirs = np.zeros( (len(phis),3) )
    dirs[:,0] = np.cos(phis)
    dirs[:,1] = np.sin(phis)
    dirs[:,2] = z0
    return dirs

def pathAroundMolecule( mol, phis, z0=0.0, R0=0.0, r_start=10.0 ):
    aRs  = np.array( [ Rvdws[ k ] for k in mol.enames ] ) + R0   # ;print("aRs ",aRs)
    p0s  = np.zeros( (len(phis),3) )
    dirs = makeDirs( phis, z0=z0 )
    for ip, ph in enumerate( phis ):
        dir     = dirs[ip]
        p0s[ip] = rayMarch( mol, dir*r_start , dir*-1, Rs=aRs, tol=1e-3 )
    return p0s, dirs

def radialPoints(  mol, phis, Rs, z0=0.0, R0=0.0, r_start=10.0 ):
    ps = []
    p0s, dirs = pathAroundMolecule(  mol, phis, z0=z0, R0=R0, r_start=r_start )
    ps = np.zeros( (len(phis),len(Rs),3) )
    #print( "p0s.shape ", p0s.shape )
    #print( "dirs.shape ", dirs.shape )
    for ir,r in enumerate( Rs ):
        ps[:,ir,:] = p0s + dirs*r 
    return ps, p0s, dirs

def plotMolAllSides(mol):
    plt.figure(figsize=(15,5))
    plt.subplot(1,3,1); plu.plotSystem( mol, bBonds=False, bLabels=False, axes=(0,1) ); plt.grid(); plt.xlabel("x [A]"); plt.ylabel("y [A]");
    plt.subplot(1,3,2); plu.plotSystem( mol, bBonds=False, bLabels=False, axes=(0,2) ); plt.grid(); plt.xlabel("x [A]"); plt.ylabel("z [A]");
    plt.subplot(1,3,3); plu.plotSystem( mol, bBonds=False, bLabels=False, axes=(1,2) ); plt.grid(); plt.xlabel("y [A]"); plt.ylabel("z [A]");

def makeRotXY( phi ):
    c = np.cos(phi)
    s = np.sin(phi)
    return np.array( [ [c, -s, 0], [s, c, 0], [0, 0, 1] ] )

def makeGeoms( ps, phis, Rs, A, B_, dirs=None, fname="out.xyz" ):
    fout = open( fname, "w" )
    for ip, phi in enumerate( phis ):
        rot  = makeRotXY( phis[ip]+np.pi )
        for ir,r in enumerate( Rs ):
            #print( "ir,ip ", ir, ip )
            B  = B_.clonePBC()
            p  = ps  [ip,ir,:]
            #d  = dirs[ip] 
            au.mulpos( B.apos, rot ) 
            B.apos[:,:] += p[None,:]
            B.append_atoms( A )
            B.toXYZ(fout, comment="%i %i R %g phi %g" %(ir,ip, r, phi), bHeader=True )
    fout.close()

# ================== Main ==================

#name1 = "HF"
#name1 = "H2O"
#name1 = "NH3"
#name1 = "C2H2"
#name1 = "C2H4"
#name1 = "CH2NH"
#name1 = "CH2O"
#name1 = "HCONH2"
#name1 = "HCOOH"

#name2 = "HF"
#name2 = "H2O"
#name2 = "NH3"
    
mol2s = [ "HF", "H2O", "NH3" ]

#mol2s = [ "HF" ]

Rs   = [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 2.0, 2.5, 3.0, 3.5, 4.5, 6.0, 8.0, 10.0, 14.0, 20.0 ]
 
for im2,mol2 in enumerate(mol2s):
    name2 = mol2
    for im1,mol1 in enumerate(mols):
        name1 = mol1[0]
        phi_max =    mol1[1]*np.pi*2.0
        nphi    = round(36*mol1[1]) + 1 
        print( "name1, name2 ", name1, name2, " phi_max ", phi_max)

        name = name1+"_"+name2

        path = './inputs/'
        mol1   = au.AtomicSystem( path+name1+".xyz", bReadN=False )  # ;print( "mol.enames ",  mol1.enames )
        mol2   = au.AtomicSystem( path+name2+".xyz", bReadN=False )

        #print( "mol.apos.shape ", mol1.apos.shape )
        # na=len( mol.apos )

        mol1.orient( 1, None, None )   # just to center the molecule around atom #0
        #mol2.orient( -1, (0,-1), (0, range(1,len(mol2.apos)-1)), _0=0, trans=(2,1,0) ) # orient the molecule by the last hydrogen atom
        mol2.orient( -1, (0,-1), (0, range(1,len(mol2.apos)-1)), _0=0, trans=(2,0,1) ) # orient the molecule by the last hydrogen atom

        #phi_max = np.pi
        #phi_max = np.pi*2
        #phis = np.arange( 0.0, phi_max+0.0001, 0.1 )
        phis = np.linspace( 0.0, phi_max, nphi )
        ps, p0s, dirs  = radialPoints(  mol1, phis, Rs )

        makeGeoms( ps, phis, Rs, mol1, mol2, fname=name+".xyz" )

        ps = ps.reshape( (-1,3) )

        #print( "ps.shape ", ps.shape )
        print( im2,im1, name + " len(ps) ", len(ps) )

        plt.figure(figsize=(5,5))
        plu.plotSystem( mol1, bBonds=False, bLabels=False, axes=(0,1) ); plt.xlabel("x [A]"); plt.ylabel("y [A]");
        plt.plot( p0s[:,0], p0s[:,1], '.-g', lw=0.2, label="Samples" )
        plt.plot( ps[:,0], ps[:,1], '.k', lw=0.2, ms=1.0, label="Samples" )
        plt.grid()
        plt.axis('equal')
        plt.savefig(name1+".png")
        #plt.show()