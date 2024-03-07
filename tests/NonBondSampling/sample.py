import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append('../../')
from pyBall import atomicUtils as au
from pyBall import plotUtils   as plu

Rvdws = {
'H'  : 1.2,
'C'  : 1.90,
'N'  : 1.78,
'O'  : 1.66,
'F'  : 1.66,
'Cl' : 1.948,
}

mols=[
"C2H2",
"C2H4",
"CH2NH",
"CH2O",
"H2O",
"HCl",
"HCN",
"HCONH2",
"HCOOH",
"HF",
"NH3",
]

def distFromMol( p, mol, Rs ):
    ds = np.sqrt( np.sum( (mol.apos[:,:] - p[None,:])**2, axis=1 ) ) -  Rs
    #print ("ds ",ds)
    return ds.min()

def rayMarch( p0, dir, Rs=None, tol=1e-3 ):
    if Rs is None:
        Rs = np.array( [ Rvdws[ mol.atypes[i] ] for i in range( len(mol.atypes) ) ] )
    p = p0
    for i in range(100):
        r = distFromMol(p, mol, Rs)
        if r < tol:
            break
        p += dir*r
    return p

#name = "HCOOH"
#name = "H2O"
#name = "NH3"
#name = "HF"
#name = "HCONH2"
name="CH2NH"

path = './inputs/'
mol   = au.AtomicSystem( path+name+".xyz" )

print( "mol.apos.shape ", mol.apos.shape )

# na=len( mol.apos )
# mol.orient( 0, (0,-1), (0, range(1,na-1)), _0=0, trans=(2,1,0) )

# plt.figure(figsize=(15,5))
# plt.subplot(1,3,1); plu.plotSystem( mol, bBonds=False, bLabels=False, axes=(0,1) ); plt.grid(); plt.xlabel("x [A]"); plt.ylabel("y [A]");
# plt.subplot(1,3,2); plu.plotSystem( mol, bBonds=False, bLabels=False, axes=(0,2) ); plt.grid(); plt.xlabel("x [A]"); plt.ylabel("z [A]");
# plt.subplot(1,3,3); plu.plotSystem( mol, bBonds=False, bLabels=False, axes=(1,2) ); plt.grid(); plt.xlabel("y [A]"); plt.ylabel("z [A]");

#phi_max = np.pi
phi_max = np.pi*2
#phis = np.arange( 0.0, phi_max+0.0001, 0.1 )
phis = np.linspace( 0.0, phi_max, 30 )
Rs   = [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 2.0, 2.5, 3.0, 3.5, 4.5, 6.0, 8.0, 10.0, 14.0, 20.0 ]

R0 = 0.0

print( "mol.enames ",  mol.enames )

aRs = np.array( [ Rvdws[ k ] for k in mol.enames ] ) + R0    ;print("aRs ",aRs)
p0s = []
for ip, ph in enumerate( phis ):
    dir = np.array( [ np.cos(ph), np.sin(ph), 0.0 ] )
    ls = mol.apos[:,0] * dir[0] + mol.apos[:,1] * dir[1] + mol.apos[:,2] * dir[2]
    r0  = 10.0
    p0 = rayMarch( dir*r0, dir*-1, Rs=aRs, tol=1e-3 )
    p0s.append( p0 )
    #exit()
p0s = np.array( p0s )

ps = []
for ip, ph in enumerate( phis ):
    dir = np.array( [ np.cos(ph), np.sin(ph), 0.0 ] )
    p0  = p0s[ip]
    for ir,r in enumerate( Rs ):
        p = p0 + dir*r 
        ps.append( p )
ps = np.array( ps )


plt.figure(figsize=(5,5))
plu.plotSystem( mol, bBonds=False, bLabels=False, axes=(0,1) ); plt.xlabel("x [A]"); plt.ylabel("y [A]");
plt.plot( p0s[:,0], p0s[:,1], '.-g', lw=0.2, label="Samples" )

plt.plot( ps[:,0], ps[:,1], '.k', lw=0.2, label="Samples" )

print( "len(ps) ", len(ps) )



plt.grid()
plt.axis('equal')
plt.show()