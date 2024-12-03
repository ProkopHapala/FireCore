import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF as mmff



def scanPlot( nscan = 1000, span=(0.0,8.0), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="E" ):
    ts = np.linspace( span[0],span[1], nscan, endpoint=False)
    poss  = np.zeros( (nscan,3) )
    poss[:,0] = p0[0] + ts*dir[0]
    poss[:,1] = p0[1] + ts*dir[1]
    poss[:,2] = p0[2] + ts*dir[2]

    Es,Fs,Ps = mmff.scan( poss, bF=True, bP=True )
    #print( "Es.shape ", Es.shape )
    plt.plot( ts, Es, '-', lw=0.5, label=label  )


#======== Body

mmff.setVerbosity( verbosity=1, idebug=1 )

#mmff.init( xyz_name="data/xyz/pyridine", surf_name="data/NaCl_1x1_L2" )    
#mmff.init( xyz_name="data/xyz/nHexadecan_dicarboxylic", bMMFF=True  )     
mmff.init( xyz_name="data/xyz/O", surf_name="data/xyz/NaCl_1x1_L3" )  
#mmff.init( xyz_name="data/xyz/H2O", surf_name="data/xyz/NaCl_1x1_L3" )    
#mmff.init( xyz_name="data/xyz/PTCDA", surf_name="data/xyz/NaCl_1x1_L3" )    
mmff.getBuffs()

#print( "ffflags ", mmff.ffflags )

mmff.setSwitches( NonBonded=-1, MMFF=-1, SurfAtoms=0, GridFF=1 )

mmff.PLQs[:,2 ] = 0.0 # delete Coulomb (charges)
#mmff.PLQs[:,:2] = 0.0 # delete Morse (EvdW)
scanPlot( nscan=1000, span=(0.0,8.0), dir=(1.0,0.0,0.0), p0=(0.0,0.0,0.0),  label="E_x" )
scanPlot( nscan=1000, span=(0.0,8.0), dir=(0.0,1.0,0.0), p0=(0.0,0.0,0.0),  label="E_y" )
#scanPlot( nscan=1000, span=(-5.0,5.0), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="E_z" )

plt.legend()
plt.grid()
plt.show()


#mmff.run()


# E = 0.0
# mmff.run(omp=True, nstepMax=20)
# print( "E=", E )
#for i in range(200):
#mmff.addSnapshot()
# with open("gopt_trajectory.xyz", "w") as file:pass
# mmff.printDatabase()

#mmff.init( xyz_name="dataxyz//pyridine", surf_name="dataxyz/NaCl_1x1_L2", bMMFF=False, gridStep=-1 )  # without gridFF
#mmff.getBuffs()
#mmff.eval()
#mmff.relax(1000)
#print( "FORCES:\n mmff.fapos:\n ", mmff.fapos )
#mmff.plot(bForce=True, Fscale=10.0 )
#plt.show()
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
xyzs,Zs,enames,qs = au.loadAtomsNP( "data/xyz/Benzene_deriv.xyz" )
au.orient( 2, (5,2), (1,3), xyzs, bFlipXZ=True )
au.saveXYZ( enames, xyzs, "data/xyz/Benzene_deriv_.xyz", qs=qs, Rs=None )
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
#mmff.initWithMolFile( "data/xyz/Benzene_deriv.xyz", bNonBonded=False, bOptimizer=True)


# mmff.initWithMolFile( "data/xyz/PTCDA.xyz", bNonBonded=False, bOptimizer=True)
# mmff.getBuffs() #;print( mmff.ndims )

#nsel = mmff.splitAtBond(5)  ;print( "split to:\n", mmff.selection[:nsel],"\n", mmff.selection[nsel:] )
#nsel = mmff.splitAtBond(6)  ;print( "split to:\n", mmff.selection[:nsel],"\n", mmff.selection[nsel:] )
#nsel = mmff.splitAtBond(10)  ;print( "split to:\n", mmff.selection[:nsel],"\n", mmff.selection[nsel:] )
#nsel = mmff.splitAtBond(2)  ;print( "split to:\n", mmff.selection[:nsel],"\n", mmff.selection[nsel:] )
#nsel = mmff.splitAtBond(4)  ;print( "split to:\n", mmff.selection[:nsel],"\n", mmff.selection[nsel:] )
#print( "nsel ", nsel, len(mmff.selection)-nsel )
#Es=mmff.scanRotation( 1, 1,11, np.pi*2, 100, bWriteTrj=True, _0=1) ;plt.plot(Es) ;print("Es=", Es)

#Es = mmff.scanBondRotation( 6, np.pi*2, 100, bWriteTrj=True );  plt.plot(Es)
#Es = mmff.scanBondRotation( 2, np.pi*2, 100, bWriteTrj=True );  plt.plot(Es); plt.grid()


#mmff.eval()
#mmff.relax(1000, Ftol=1e-4, bWriteTrj=True )
#Es=mmff.scanRotation( 1, 1,11, np.pi*2, 100, sel=[11,13,14,20]+[29,30,31,32], bWriteTrj=True, _0=1) ;plt.plot(Es) ;print("Es=", Es)

# plt.figure()
# mmff.plot()
# #mmff.plot_selection( mmff.selection[:nsel] )
# #mmff.plot_selection( [1,2,3] )

# plt.show()
