# https://pyscf.org/quickstart.html
# https://pyscf.org/pyscf_api_docs/pyscf.gto.html#module-pyscf.gto.mole

'''
Molecules:
H2O
NH3
HCOOH
CH2=O
'''

import sys
import os
from tkinter import UNITS
import pyscf
from pyscf.geomopt.berny_solver import optimize
sys.path.append('../../')
from pyBall import atomicUtils as au
import numpy as np
import matplotlib.pyplot as plt

# ============= Setup
verbosity = 0
conv_params = {  # These are the default settings
    'gradientmax': 0.45e-6,  # Eh/[Bohr|rad]
    'gradientrms': 0.15e-6,  # Eh/[Bohr|rad]
    'stepmax': 1.8e-3,       # [Bohr|rad]
    'steprms': 1.2e-3,       # [Bohr|rad]
}
hartree2eV = 27.211396641308
bohr2A     = 0.52917721090380 
# ============ Functions

def unpack_mol( mol ):
    apos= np.array([ a[1] for a in mol._atom ])
    es  = np.array([ a[0] for a in mol._atom ])
    return apos, es

def pack_mol( apos, es ):
    return [ (es[i],apos[i]) for i in range(len(es))]

def plotAtoms( es=None, apos=None, atoms=None, ax1=0, ax2=1 ):
    if apos is None: apos = np.array([ a[1] for a in atoms ])  #;print(apos)
    plt.plot( apos[:,ax1],apos[:,ax2], 'o' ); plt.axis('equal'); plt.grid()

def printlist(lst):
    for item in lst: print( item )

def printObj(obj):
    printlist(dir(obj))

def saveAtoms(fname,atoms, unit=bohr2A ):
    apos = np.array([ a[1] for a in atoms ])*unit
    es   = [ a[0] for a in atoms ]
    au.saveXYZ( es, apos, fname )

def mulpos( ps, rot ):
    for ia in range(len(ps)): ps[ia,:] = np.dot(rot, ps[ia,:])

def preparemol(fname='relaxed.xyz'):
    if os.path.isfile(fname):
        print("found(%s) => no need for relaxation " %fname )
        mol = pyscf.M(atom=fname)
    else:
        h2o = pyscf.M(atom = 'O 0 0 0; H 1 0 0; H 0 1 0'          )
        #nh3 = pyscf.M(atom = 'N 0 0 0; H 1 0 0; H 0 1 0; H 0 0 1;')
        h2o.verbose = verbosity
        calc = pyscf.scf.RHF(h2o)
        #calc = pyscf.dft.RKS(h2o); calc.xc = 'b3lyp'; #calc.init_guess='atom'
        mol = optimize(calc, maxsteps=1000, **conv_params)
        saveAtoms(fname,mol._atom)
    return mol 

def linearScan( molFix, molMov, xs, dir=(1.,0,0), bEcho=True, xyz_file="scan.xyz", Eout_file=None ):
    dir=np.array(dir); dir/=np.linalg.norm(dir)
    if xyz_file is not None: fout = open(xyz_file,'w')
    Es = np.zeros(xs.shape)
    apos1,es1=molFix
    apos2,es2=molMov
    for i,x in enumerate(xs):
        apos_ = apos2.copy()
        apos_[:,:] += dir[None,:]*x
        mol=pyscf.M(atom= pack_mol( apos1, es1 ), unit='B')   #;print( " mol._atom \n", mol.atom )
        mol.atom.extend( pack_mol( apos_, es2 ) ) 
        mol.build()
        #print( mol._atom )
        mol.verbose = verbosity 
        out = pyscf.scf.UHF(mol).run()
        #printObj(out);exit()
        #printObj( out.energy_elec() );
        #printObj( out.energy_nuc()  ); exit()
        if(bEcho):print( x, out.e_tot, out.energy_elec(), out.energy_nuc() )
        Es[i] = out.e_tot
        #print( " out.mol._atom  \n", out.mol._atom )
        if xyz_file is not None:
            apos, es_ = unpack_mol(out.mol)  ;apos*=bohr2A
            au.writeToXYZ( fout, es_, apos, commet=(" x %g E_tot %g " %(x, out.e_tot)) )
    if xyz_file is not None: fout.close()
    if Eout_file is not None: np.savetxt( Eout_file, np.array( [xs, Es*hartree2eV] ).transpose() )
    return Es

def angularScan( n, dang, molFix, molMov, rot=None, ax1=0,ax2=1, bEcho=True, xyz_file="scan.xyz", Eout_file=None ):
    if xyz_file is not None: fout = open(xyz_file,'w')
    if rot is None: rot = au.makeRotMatAng( dang, ax1=ax1, ax2=ax2 )
    Es   = np.zeros(n)
    angs = np.zeros(n)
    apos1,es1=molFix
    apos2,es2=molMov
    apos_ = apos2.copy()
    for i in range(n):
        angs[i]=dang*i
        mol=pyscf.M(atom= pack_mol( apos1, es1 ), unit='B')   #;print( " mol._atom \n", mol.atom )
        mol.atom.extend( pack_mol( apos_, es2 ) ) 
        mol.build()
        #print( mol._atom )
        mol.verbose = verbosity 
        if True:   # solver callback
            out = pyscf.scf.UHF(mol).run()
            if(bEcho):print( angs[i], out.e_tot, out.energy_elec(), out.energy_nuc() )
            Es[i] = out.e_tot
        #mol = out.mol
        #print( " out.mol._atom  \n", out.mol._atom )
        if xyz_file is not None:
            apos, es_ = unpack_mol(mol)  ;apos*=bohr2A
            au.writeToXYZ( fout, es_, apos, commet=(" x %g E_tot %g " %(dang,Es[i])) )
        mulpos( apos_, rot )
    if xyz_file is not None: fout.close()
    if Eout_file is not None: np.savetxt( Eout_file, np.array( [angs, Es*hartree2eV] ).transpose() )
    return Es,angs

# ============ MAIN

mol = preparemol(fname='relax.xyz')

apos,es = unpack_mol( mol )   #;print("apos ",apos) ;print("es ",es)
apos1=au.orient( 0,(1,2),(1,0), apos, _0=0, trans=(2,1,0) )
apos2=au.orient( 0,(0,1),(0,2), apos, _0=0, trans=(0,2,1) )     ;apos2[:,1]*=-1.0 ;apos2[:,1]+=6.0
#mol=pyscf.M(atom= pack_mol( apos1, es ), unit='B')   #;print( " mol._atom \n", mol.atom )
#mol.atom.extend( pack_mol( apos2, es ) )             #;print( " mol._atom \n", mol.atom )
#mol.build()
#print( " mol._atom \n", mol._atom )
#mol_1=pyscf.M(atom=mol, unit='B')

# ------- Linear scan
#xs = np.arange(-1.5,4.0,0.2)  ;print("shifts ", xs)
#Es = linearScan( (apos1,es), (apos2,es), xs, dir=(0.,1.,0), bEcho=True, xyz_file="lin_scan.xyz",  Eout_file="lin_scan_Es.dat" )

# -------- Angular scan
nrot = 16
#dang= (np.pi*(3./2.))/nrot  ;print("dang ", dang )
dang= 2.0*np.pi/nrot  ;print("dang ", dang )
rot0 = au.makeRotMatAng( dang*nrot*0.5, ax1=1, ax2=2 ).transpose()
rot  = au.makeRotMatAng( dang , ax1=1, ax2=2 ).transpose()

print( rot ); # exit()
mulpos( apos2, rot0 )

Es,xs = angularScan( nrot, dang, (apos1,es), (apos2,es), rot=rot, bEcho=True, xyz_file="ang_scan.xyz", Eout_file="ang_scan_Es.dat" )

plt.plot( xs, Es*hartree2eV );   plt.grid(); 


#np.savetxt( "E_scan.dat", np.array( [xs, Es*hartree2eV] ).transpose() )
#plt.plot( )




#print("---------")
#plotAtoms( atoms=mol._atom )
#plotAtoms( atoms=mol_1._atom )
#plotAtoms( atoms=out.mol._atom )
#plotAtoms( atoms=h2o._atom )
plt.show()


#au.makeRotMat( fw, up )

