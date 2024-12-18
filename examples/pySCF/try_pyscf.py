# https://pyscf.org/quickstart.html
# https://pyscf.org/pyscf_api_docs/pyscf.gto.html#module-pyscf.gto.mole

import pyscf
#from pyscf import scf
#from pyscf import gto
#from pyscf.geomopt.geometric_solver import optimize
from pyscf.geomopt.berny_solver import optimize

# geometric
#from pyscf.geomopt.geometric_solver import optimize
#mol_eq = optimize(mf, maxsteps=100)
#print(mol_eq.atom_coords())

# pyberny
#from pyscf.geomopt.berny_solver import optimize
#mol_eq = optimize(mf, maxsteps=100)
#print(mol_eq.atom_coords())

#mol = pyscf.gto.M(atom = 'O 0 0 0; H 0 1 0; H 0 0 1', basis = 'ccpvdz')
h2o = pyscf.M(atom = 'O 0 0 0; H 0 1 0; H 0 0 1'          )
nh3 = pyscf.M(atom = 'N 0 0 0; H 1 0 0; H 0 1 0; H 0 0 1;')
#mol  = gto.M(atom = 'C 0 0 .625; C 0 0 -.625'  , symmetry = 'd2h')

mol=pyscf.M(atom = h2o._atom,basis = 'ccpvdz')
#mol._atom.extend(h2o._atom)
mol._atom.extend(nh3._atom)

#mol.build()
#print( dir(mol) )
#print( mol._atom )
exit()

mol.verbose = 0
mol.output = 'pyscf.log'
mol.max_memory = 1000 # MB

hartree2eV = 27.211396641308


# ------ Hartree Focks
rhf_h2o = pyscf.scf.RHF(mol)
#e_h2o   = rhf_h2o.kernel()
mol_eq = optimize(rhf_h2o, maxsteps=100)
print(mol_eq.atom_coords())

# ----- DFT
#rks_h2o = pyscf.dft.RKS(mol) # likewise for UKS and ROKS
#rks_h2o.xc = 'b3lyp'

# ------ MP2
#mp2_h2o  = pyscf.mp.MP2(rhf_h2o)
#e_h2o    = mp2_c2.kernel()[0]