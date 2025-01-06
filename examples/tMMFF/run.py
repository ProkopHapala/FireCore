import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

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

###### Mexican hat potential

nbStep = 100
nMDsteps = 100000
nEQsteps = 10000
t_damp = 20
T = 10
dt = 0.05

mmff.setVerbosity(verbosity=3)
mmff.init( xyz_name="data/enthropic_spring_10", bMMFF=True )
collectiveVariable = np.array([0], dtype=np.int32)
E = mmff.compute_Free_energy(0.5, 4.0, collectiveVariable, nbStep=nbStep, nMDsteps=nMDsteps, nEQsteps=nEQsteps, t_damp=t_damp, T=T, dt=dt)
print("E=", E)
print("Konec Milane")



'''# entropic spring
nbStep = 100
nMDsteps = 100000
t_damp = 100
T = 500
dt = 0.05
N = [5, 10, 20, 30]

import argparse

# Create parser for command-line arguments
parser = argparse.ArgumentParser(description='Run simulation with specified parameters')
parser.add_argument('--nMDsteps', type=int, help='Number of MD steps')
parser.add_argument('--nEQsteps', type=int, help='Number of equilibration steps')
parser.add_argument('--t_damp', type=float, help='Damping time')
parser.add_argument('--T', type=float, help='Temperature')
parser.add_argument('--dt', type=float, help='Time step')
parser.add_argument('--nSteps', type=int, help='Number of steps')
parser.add_argument('--N', type=int, help='System size')
parser.add_argument('--k', type=float, help='Spring constant')
parser.add_argument('--K', type=float, help='Spring constant')

# Parse arguments
args = parser.parse_args()

# Assign values from command line arguments
MY_N = args.N
MY_K = args.K  # Kept constant as it wasn't in command line args
MY_k = args.k  # Kept constant as it wasn't in command line args
MY_nbStep = args.nSteps
MY_nMDsteps = args.nMDsteps
MY_nEQsteps = args.nEQsteps
MY_t_damp = args.t_damp
MY_T = args.T
MY_dt = args.dt

xyz_name = "data/enthropic_spring_"+str(MY_N)
constr_name = "data/enthropic_spring_"+str(MY_N)+".cons"


with open(constr_name, "r") as file:
    line = file.readline().strip()
parts = line.split()
for i in range(5, 7):
    parts[i] = str(MY_K)

# Join the parts back together
new_line = ' '.join(parts)

# Write the modified content back to the file
with open(constr_name, "w") as file:
    file.write(new_line)

with open("stiffness.txt", "w") as file:
    file.write(str(MY_K))


mmff.setVerbosity( verbosity=1, idebug=1 )

mmff.init( xyz_name=xyz_name, constr_name=constr_name ,bMMFF=True)
mmff.setSwitches(CheckInvariants=-1, PBC_nonBond=-1, PBC_evalAtom=-1, NonBonded=-1, MMFF=1, doBonds=1, Angles=-1, PiSigma=-1, PiPiI=-1, bNonBondNeighs=-1, bSurfAtoms=-1, bGridFF=-1, bTricubic=-1, bConstrZ=-1, bConstrains=-1)
colectiveVariable = np.array([0], dtype=np.int32)
E = mmff.compute_Free_energy(1.0, 3.0, colectiveVariable, nbStep=MY_nbStep, nMDsteps=MY_nMDsteps, nEQsteps=MY_nbStep, t_damp=MY_t_damp, T=MY_T, dt=MY_dt)
print("E=", E)
'''



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

'''
natoms = 50
N = 100

#mmff.init( xyz_name="data/nHexadecan_dicarboxylic", bMMFF=True  )
#     


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
>>>>>>> Milan1:examples/tMMFF/run.py
#mmff.getBuffs()
#mmff.eval()
#mmff.relax(1000)
#print( "FORCES:\n mmff.fapos:\n ", mmff.fapos )
#mmff.plot(bForce=True, Fscale=10.0 )
#plt.show()
exit(0)


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

