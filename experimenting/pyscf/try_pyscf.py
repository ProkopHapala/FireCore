# https://pyscf.org/quickstart.html

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

mol = pyscf.gto.M(atom = 'O 0 0 0; H 0 1 0; H 0 0 1', basis = 'ccpvdz')
#mol  = gto.M(atom = 'C 0 0 .625; C 0 0 -.625'  , symmetry = 'd2h')

mol.verbose = 0
mol.output = 'pyscf.log'
mol.max_memory = 1000 # MB

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