import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF as mmff


#======== Body
'''
###### Mexican hat potential
mmff.init( xyz_name="data/H", bMMFF=True )
collectiveVariable = np.array([0], dtype=np.int32)
E = mmff.compute_Free_energy(0.5, 2.0, collectiveVariable)
print("E=", E)
print("Konec Milane")
'''

'''
mmff.setVerbosity( verbosity=1, idebug=1 )
mmff.init( xyz_name="data/enthropic_spring_10", constr_name="data/enthropic_spring_10.cons" ,bMMFF=True)
mmff.setSwitches(CheckInvariants=-1, PBC_nonBond=-1, PBC_evalAtom=-1, NonBonded=-1, MMFF=1, doBonds=1, Angles=-1, PiSigma=-1, PiPiI=-1, bNonBondNeighs=-1, bSurfAtoms=-1, bGridFF=-1, bTricubic=-1, bConstrZ=-1, bConstrains=-1)
colectiveVariable = np.array([0], dtype=np.int32)
E = mmff.compute_Free_energy(1.0, 7.0, colectiveVariable)
# print("E=", E)
'''

mmff.init( xyz_name="data/three_atoms", bMMFF=True, nPBC=(2, 2, 2))#, constr_name="data/three_atom.cons" ) #nPBC=(2, 2, 2)
colectiveVariable = np.array([0], dtype=np.int32)
nbStep = 100
nMDsteps = 1000000
t_damp = 20
T = 300
dt = 0.05
E = mmff.compute_Free_energy(0.0, 5.0, colectiveVariable, nbStep=nbStep, nMDsteps=nMDsteps, t_damp=t_damp, T=T, dt=dt)
mmff.setSwitches(CheckInvariants=-1, PBC_nonBond=-1, PBC_evalAtom=-1, NonBonded=1, MMFF=1, doBonds=-1, Angles=-1, PiSigma=-1, PiPiI=-1, bNonBondNeighs=-1, bSurfAtoms=-1, bGridFF=-1, bTricubic=-1, bConstrZ=-1, bConstrains=-1)
print("E=", E)


'''
natoms = 50
N = 100

#mmff.init( xyz_name="data/nHexadecan_dicarboxylic", bMMFF=True  )
#     
# E = 0.0
# mmff.run(omp=True, nstepMax=20)
# print( "E=", E )
#for i in range(200):
#mmff.addSnapshot()
# with open("gopt_trajectory.xyz", "w") as file:
#     pass
a = np.full( natoms*3, -10.0 )
b = np.full( natoms*3, 10.0 )
boundaryRules = np.zeros(natoms*3, dtype=int)

RMSD_array = np.zeros(N)
Fs_array = np.zeros(N) 
nevals_array = np.zeros(N, dtype=np.int32)

for k in range (N):
    RMSD,Fs,nevals=mmff.runGlobalOptimization(maxeval=1000000, bSave=1, bShow=0, par_mut=np.array([0.1]), a=a, b=b, boundaryRules=boundaryRules)
    RMSD_array[k] = RMSD
    Fs_array[k] = Fs
    nevals_array[k] = nevals

data = np.column_stack((RMSD_array, Fs_array, nevals_array))
np.savetxt('results.txt', data, delimiter=',', header='RMSD,Fs,nevals', comments='')
'''

# # Histogram RMSD
# plt.figure(figsize=(10, 5))
# plt.hist(RMSD_array, bins=20, edgecolor='black')
# plt.xlabel('RMSD')
# plt.ylabel('Frequency')
# plt.title('Histogram of RMSD')
# plt.grid(True)
# plt.show()

# # Histogram Fs
# plt.figure(figsize=(10, 5))
# plt.hist(Fs_array, bins=20, edgecolor='black')
# plt.xlabel('Fs')
# plt.ylabel('Frequency')
# plt.title('Histogram of Fs')
# plt.grid(True)
# plt.show()

# # Histogram nevals
# plt.figure(figsize=(10, 5))
# plt.hist(nevals_array, bins=20, edgecolor='black')
# plt.xlabel('nevals')
# plt.ylabel('Frequency')
# plt.title('Histogram of nevals')
# plt.grid(True)
# plt.show()


#mmff.init( xyz_name="data/pyridine", surf_name="data/NaCl_1x1_L2", bMMFF=False, gridStep=-1 )  # without gridFF
#mmff.getBuffs()
#mmff.eval()
#mmff.relax(1000)
#print( "FORCES:\n mmff.fapos:\n ", mmff.fapos )
#mmff.plot(bForce=True, Fscale=10.0 )
#plt.show()
exit(0)

'''
n = 10
D = []
N = 10000
T = np.zeros(N)

for _ in range(n):
    mmff.init(xyz_name="data/polymer-2_new", surf_name="data/NaCl_1x1_L2", bMMFF=True)

    start_time = time.time() 
    mmff.run(nstepMax=N, outE=T, omp=True)
    end_time = time.time() 

    elapsed_time = end_time - start_time 
    
    average_T = T[8000:10000].mean()

    D.append((average_T, elapsed_time))
    print(f"{average_T:<20} K {elapsed_time:<20.4f} second")
    mmff.clear()

for average_T, elapsed_time in D:
    print(f"{average_T:<10.4f} K {elapsed_time:<10.4f} second")

plt.plot(T)
plt.xlabel('Iterace')
plt.ylabel('T')
plt.title('Závislost T na iteraci')
plt.grid(True)
plt.show()
'''


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
'''