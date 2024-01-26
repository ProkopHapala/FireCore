
import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("/home/niko/work/FIRECORE/FireCore/")
from pyBall import atomicUtils as au
from pyBall import MMFF as mmff


#======== Body
mmff.setVerbosity( verbosity=1, idebug=1 )
#mmff.setVerbosity( verbosity=2, idebug=1 )

# MolWorld_sp3::init
mmff.init( xyz_name=str(sys.argv[1]),
           surf_name =None, 
           smile_name=None, 
           sElementTypes  = "/home/niko/work/FIRECORE/FireCore/tests/tUFF/data_UFF/ElementTypes.dat",
           sAtomTypes     = "/home/niko/work/FIRECORE/FireCore/tests/tUFF/data_UFF/AtomTypes.dat", 
           sBondTypes     = "/home/niko/work/FIRECORE/FireCore/tests/tUFF/data_UFF/BondTypes.dat", 
           sAngleTypes    = "/home/niko/work/FIRECORE/FireCore/tests/tUFF/data_UFF/AngleTypes.dat",
           sDihedralTypes = "/home/niko/work/FIRECORE/FireCore/tests/tUFF/data_UFF/DihedralTypes.dat",
           bMMFF=True,
           bEpairs=False,
           nPBC=(1,1,1),
           gridStep=0.1,
           bUFF=True,
           b141=True,
           bSimple=MYbSimple,
           bConj=MYbConj,
           bCumulene=True)
#exit(0)

# MolWorld_sp3::getBuffs
#mmff.getBuffs()
#exit(0)

# MolWorld_sp3::eval
mmff.eval()
exit(0)

#mmff.relax(1000)
print( "FORCES:\n mmff.fapos:\n ", mmff.fapos )
mmff.plot(bForce=True, Fscale=10.0 )
plt.show()
exit(0)


'''
mmff.init()
#mmff.init_params( "data/AtomTypes.dat", "data/BondTypes.dat", "data/AngleTypes.dat" )
#mmff.insertSMILES("CC");
#mmff.insertSMILES("C=C");
#mmff.insertSMILES("C#C");
#mmff.insertSMILES("C#CCN=C", True );
#mmff.insertSMILES("C1#CCN=C1", True );
#mmff.insertSMILES("C=C1NC#CC1CO", True, True );

#mmff.initWithSMILES( "C=C1NC#CC1CO" )
mmff.initWithSMILES( "C=C" )
mmff.getBuffs()
mmff.relax(1000)
mmff.plot()
plt.show()

exit()
'''

'''
# ======== Oritent Molecule
xyzs,Zs,enames,qs = au.loadAtomsNP( "data/Benzene_deriv.xyz" )
au.orient( 2, (5,2), (1,3), xyzs, bFlipXZ=True )
au.saveXYZ( enames, xyzs, "data/Benzene_deriv_.xyz", qs=qs, Rs=None )
plt.plot( xyzs[:,0],xyzs[:,1], "o" )
plt.axis('equal')
plt.show()
exit()
'''

'''
# ============== C2H4,xyz
#mmff.initWithMolFile( "C2H4.xyz", bNonBonded=False, bOptimizer=True)
#mmff.printBuffNames()
#mmff.getBuffs() #;print( mmff.ndims )
#mmff.eval()
#mmff.relax(1000, bWriteTrj=True )
#Es=mmff.scanRotation( [1,4,5], 0, 0,1, np.pi*2, 100, bWriteTrj=True)   ;print("Es=", Es)
#plt.plot(Es)
#print( "Es(Etot,Eb,Ea,Eps,EppT,EppI):", mmff.Es )
#nsel = mmff.splitAtBond(6-1)  ;print( "split to:\n", mmff.selection[:nsel],"\n", mmff.selection[nsel:] )
'''

# ============== Benzene_deriv.xyz
mmff.initWithMolFile( "data/Benzene_deriv.xyz", bNonBonded=False, bOptimizer=True)
mmff.getBuffs() #;print( mmff.ndims )

#nsel = mmff.splitAtBond(5)  ;print( "split to:\n", mmff.selection[:nsel],"\n", mmff.selection[nsel:] )
#nsel = mmff.splitAtBond(6)  ;print( "split to:\n", mmff.selection[:nsel],"\n", mmff.selection[nsel:] )
#nsel = mmff.splitAtBond(10)  ;print( "split to:\n", mmff.selection[:nsel],"\n", mmff.selection[nsel:] )
#nsel = mmff.splitAtBond(2)  ;print( "split to:\n", mmff.selection[:nsel],"\n", mmff.selection[nsel:] )
#nsel = mmff.splitAtBond(4)  ;print( "split to:\n", mmff.selection[:nsel],"\n", mmff.selection[nsel:] )
#print( "nsel ", nsel, len(mmff.selection)-nsel )
#Es=mmff.scanRotation( 1, 1,11, np.pi*2, 100, bWriteTrj=True, _0=1) ;plt.plot(Es) ;print("Es=", Es)

#Es = mmff.scanBondRotation( 6, np.pi*2, 100, bWriteTrj=True );  plt.plot(Es)
Es = mmff.scanBondRotation( 2, np.pi*2, 100, bWriteTrj=True );  plt.plot(Es); plt.grid()


#mmff.eval()
#mmff.relax(1000, Ftol=1e-4, bWriteTrj=True )
#Es=mmff.scanRotation( 1, 1,11, np.pi*2, 100, sel=[11,13,14,20]+[29,30,31,32], bWriteTrj=True, _0=1) ;plt.plot(Es) ;print("Es=", Es)

plt.figure()
mmff.plot()
#mmff.plot_selection( mmff.selection[:nsel] )
#mmff.plot_selection( [1,2,3] )

plt.show()
