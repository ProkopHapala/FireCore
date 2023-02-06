import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

sys.path.append("../../")
from pyBall import eFF as eff

from pyBall import eFF_terms as pyeff
const_bohr = 0.5291772105638411

'''
rho=-0.2; r=1.2; si=0.3; sj=0.8
rho=-0.2; r=0.7; si=1.3; sj=0.15

Euuref,DTref, Sref = pyeff.pyeff_E_up_up   ( r, si, sj, rho,)
Eudref             = pyeff.pyeff_E_up_down ( r, si, sj, rho,)
Euu,Eud, DT, S     = pyeff.pyeff_EPaul     ( r, si, sj, rho,)

print( "Eupdown ", Eudref,  Eud  )
print( "Eupup   ", Euuref,  Euu  )
print( "S       ", Sref,  S  )
print( "DT      ", DTref, DT )

eff.eval_ee( r, si, sj )
'''
# ================== Test Derivatives of various terms (Analytical vs Numerical Derivatives)

#eff.check_DerivsPauli( r0=0.0,r1=2.5,   s0=0.5,s1=0.5, n=100 , spin=1 )
#eff.check_DerivsPauli( r0=0.0,r1=2.5,   s0=0.5,s1=0.5, n=100 , spin=-1 )
#eff.check_DerivsPauli( r0=0.7,r1=0.7,   s0=0.2,s1=2.5,  sj=0.5, n=100 , spin=1)
#eff.check_DerivsPauli( r0=0.7,r1=0.7,   s0=0.2,s1=2.5,  sj=0.5, n=100 , spin=-1)
#exit()

#eff.setSwitches( pauli=1,   kinetic=-1, coulomb=-1, AA=-1, AE=-1, AECoulomb=-1, AEPauli=-1  )
#eff.setSwitches( pauli=-1,   kinetic=1, coulomb=1, AA=1, AE=1, AECoulomb=1, AEPauli=1  )
#eff.setSwitches( pauli=1,   kinetic=-1, coulomb=1, AA=1, AE=1, AECoulomb=-1, AEPauli=1  )

#eff.check_Derivs_ie( "e2_1g_2o_triplet", ie=0, r0=0.5,r1=1.5, s0=0.5,s1=0.5, n=100 )
#eff.check_Derivs_ie( "H2_eFF_asym", ie=0, r0=0.5,r1=1.5, s0=0.5,s1=0.5, n=100 )
#eff.check_Derivs_ie( "H2_eFF_asym", ie=0, r0=1.0,r1=1.0, s0=0.5,s1=1.5, n=100 )
#eff.checkNumDerivs( "H2_eFF_asym"); exit()

eff.setVerbosity(1)


# ================ Check Pauli terms vs python (pyeff) reference   (see reference implementation https://github.com/ProkopHapala/pyeff fork of  https://github.com/pyflosic/pyeff )

#Euu,_,_,_    = pyeff.E_up_up__( 1.125e-7, 1.7007535132532354, 1.7007535132532354, rho=-0.2 )  
#Euu,DT,S     = pyeff.pyeff_E_up_up(  1.125e-7, 1.7007535132532354, 1.7007535132532354, rho=-0.2 ); print( "!!!!! pyeff_E_up_up(): Euu",Euu, "DT",DT, "S",S )
#Euu,DT,S     = pyeff.pyeff_E_up_up(  1.125e-7, 0.9/const_bohr, 0.9/const_bohr, rho=-0.2 );     print( "!!!!! pyeff_E_up_up(): Euu",Euu, "DT",DT, "S",S )
#Euu,Eud,DT,S = pyeff.pyeff_EPaul( 1e-5, 0.9/const_bohr, 0.9/const_bohr, rho=-0.2 );            print( "!!!!! pyeff_EPaul():   Euu",Euu, "Eud",Eud, "DT",DT, "S",S )

# ================ Check basic objects ( Hydrogen atom and H2 molecule  )

#eff.test_Hatom()        ;exit()
#eff.test_Hatom(True)    ;exit()
#eff.check_H2(False)     ;exit()
#eff.check_H2()          ;exit()

#eff.check_H2( False, "H2_eFF_upup")     ;exit()
#eff.check_H2( False, "H2_eFF_asym")     ;exit()
#eff.check_H2( False, "H2_eFF_sym")     ;exit()
#eff.check_H2( True, "H2_eFF")     ;exit()

# ================= Relax Molecules to compare with pyeff (see reference implementation https://github.com/ProkopHapala/pyeff fork of  https://github.com/pyflosic/pyeff )

#outE=eff.eval_mol("Li_eFF"          ,fUnits=const_bohr, outE=True)     ;exit()
#outE=eff.eval_mol("Li_eFF_fixcore"  ,fUnits=const_bohr, outE=True)     ;exit()

#outE=eff.eval_mol("CH4_lmps"        ,fUnits=const_bohr, outE=True)    ;exit()
#outE=eff.eval_mol("CH4_lmps_fixcore",fUnits=const_bohr, outE=True)    ;exit()

# ================= Relax Molecules

#outE=eff.relax_mol("H_eFF"  ,outE=True)
#outE=eff.relax_mol("H2_eFF" ,outE=True)

#outE=eff.relax_mol("O-7",         dt=0.03,  outE=True)   #;plt.plot(outE); plt.show(); exit()                                #   Oxygen with one electrons (no fixed core) - does it converge ?
#outE=eff.relax_mol("O-6",         dt=0.02,  outE=True)   #;plt.plot(outE); plt.show(); exit()                                #   Oxygen with 2   electrons (no fixed core) - does it converge ?
#outE=eff.relax_mol("H2O",         dt=0.005, outE=True, nMaxIter=50000, perN=10 )   #;plt.plot(outE); plt.show(); exit()      #   H2O   (no fixed core)                     - does it converge ?
#outE=eff.relax_mol("H2O_fixcore", dt=0.03, outE=True)   
#outE=eff.relax_mol("CH4", outE=True)
#outE=eff.relax_mol("CH4_fixcore", outE=True)
#outE=eff.relax_mol("CH4_lmps_fixcore", fUnits=const_bohr, outE=True)

#outE=eff.relax_mol("H2O_fixcore", dt=0.03, outE=True)  ;exit()
#outE=eff.relax_mol("H2O_noe1_", dt=0.03, outE=True, bFixNuclei=True)   ;exit()
#outE=eff.relax_mol("H2O_noe2_", dt=0.03, outE=True,  bFixNuclei=True)   ;exit()

outE=eff.relax_mol("H2O_noe1_", dt=0.03, outE=True, bFixNuclei=False)   ;exit()
#outE=eff.relax_mol("H2O_noe2_", dt=0.03, outE=True,  bFixNuclei=True)   ;exit()

#outE=eff.relax_mol("H2O_noe1", dt=0.03, outE=True)   ;exit()
#outE=eff.relax_mol("H2O_noe2", dt=0.03, outE=True)  ;exit()

#plt.plot(outE-outE[-1]); plt.yscale('log'); plt.grid(); plt.show(); exit()

# ==================== Scan Dependence of bond-lenght in CH4, NH3, H2O on core size (with fized core electrons) 

name="H2O_fixcore"
#name="NH3_fixcore"
#name="CH4_fixcore"
#core_sizes = [0.5,0.4,0.3,0.2,0.1]
#core_sizes = np.arange(0.05,0.50,0.01) #[::-1]
core_sizes  = np.arange(0.1,0.40,0.01) #[::-1]
bondLengths = eff.scan_core_size( name, core_sizes, dt=0.01,nMaxIter=100000, fUnits=1. ); 
plt.plot(core_sizes,bondLengths, '.-' ); plt.xlabel("core_size(X)[A]"); plt.ylabel("L(X-H) [A]"); plt.grid(); plt.savefig(name+"_CoreSizeScan.png", bbox_inches='tight');  plt.show(); exit()