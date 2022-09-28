# see https://psicode.org/psi4manual/master/psiapi.html
# installation : https://psicode.org/installs/v161/

# This is the important part
import psi4

psi4.set_memory('500 MB')
psi4.core.set_output_file('psi4.log', False)
#psi4.core.be_quiet()


h2o = psi4.geometry("""
O
H 1 0.96
H 1 0.96 2 104.5
""")

#psi4.energy('scf/cc-pvdz')

psi4.set_options({'reference': 'rhf'})
E = psi4.optimize('scf/cc-pvdz', molecule=h2o)

hartree2eV = 27.211396641308
print( "================================" )
print( "E = %g [Hartree]", E*hartree2eV )