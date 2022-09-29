
import sys
import os
import numpy as np
import pyscf
from pyscf.geomopt.berny_solver import optimize
from . import atomicUtils as au

hartree2eV = 27.211396641308
bohr2A     = 0.52917721090380 
verbosity=1

verbosity = 0
default_conv_params = {  # These are the default settings
    'gradientmax': 0.45e-6,  # Eh/[Bohr|rad]
    'gradientrms': 0.15e-6,  # Eh/[Bohr|rad]
    'stepmax': 1.8e-3,       # [Bohr|rad]
    'steprms': 1.2e-3,       # [Bohr|rad]
}

# ============ Functions

def unpack_mol( mol, units=bohr2A ):
    apos= np.array([ a[1] for a in mol._atom ]) * units
    es  = np.array([ a[0] for a in mol._atom ])
    return apos, es

def pack_mol( apos, es ):
    return [ (es[i],apos[i]) for i in range(len(es))]

def printlist(lst):
    for item in lst: print( item )

def printObj(obj):
    printlist(dir(obj))

def saveAtoms(fname,atoms, unit=bohr2A ):
    apos = np.array([ a[1] for a in atoms ])*unit
    es   = [ a[0] for a in atoms ]
    au.saveXYZ( es, apos, fname )

def preparemol(fname='relaxed.xyz', conv_params=None, atoms='O 0 0 0; H 1 0 0; H 0 1 0' ):
    if conv_params is None: conv_params=default_conv_params
    if os.path.isfile(fname):
        print("found(%s) => no need for relaxation " %fname )
        mol = pyscf.M(atom=fname)
    else:
        h2o = pyscf.M(atom=atoms)
        h2o.verbose = verbosity
        calc = pyscf.scf.RHF(h2o)
        #calc = pyscf.dft.RKS(h2o); calc.xc = 'b3lyp'; #calc.init_guess='atom'
        mol = optimize(calc, maxsteps=1000, **conv_params)
        saveAtoms(fname,mol._atom)
    return mol 

def evalHf(inp):
    apos,es = inp            #;print( apos, es )
    m = pack_mol( apos, es)  #;print( m )
    mol=pyscf.M( atom=pack_mol( apos, es) )   #;print( " mol._atom \n", mol.atom )
    #mol.build()             #;print( mol._atom )
    mol.verbose = verbosity 
    out = pyscf.scf.UHF(mol).run()
    #if(verbosity>0): print( x, out.e_tot, out.energy_elec(), out.energy_nuc() )
    return out.e_tot*hartree2eV

def optHf(atoms, conv_params=None ):
    if conv_params is None: conv_params=default_conv_params
    print(atoms)
    job = pyscf.M(atom=atoms)
    job.SCF_max_cycle=100
    job.verbose = verbosity
    calc = pyscf.scf.RHF(job)
    mol = optimize(calc, maxsteps=1000, **conv_params)
    printlist(mol)
    return mol