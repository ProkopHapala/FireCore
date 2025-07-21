import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF as mmff


########## 1D Scan #####################################################
def scanPlot(nscan=1000, span=(0.0,8.0), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="E", saveFig=None, saveData=None):
    
    ts = np.linspace(span[0], span[1], nscan, endpoint=False)
  
    poss = np.zeros((nscan, 3))
    poss[:, 0] = p0[0] + ts * dir[0]
    poss[:, 1] = p0[1] + ts * dir[1]
    poss[:, 2] = p0[2] + ts * dir[2]

   
    Es, Fs, Ps = mmff.scan(poss, bF=True, bP=True)
    

    if saveData is not None:
        np.savetxt(saveData, np.column_stack((ts, Es)), header="ts\tEnergy", comments="# ")

    plt.title(label)
    plt.plot(ts, Es, '-', lw=0.5, label=label)
    plt.xlabel(f"Scaned along ({dir[0]}_{dir[1]}_{dir[2]}) direction ")
    plt.ylabel(f"Scaned Energy (eV)")
    
    # Optionally, save the figure to a file.
    if saveFig is not None:
        plt.savefig(saveFig)
    plt.show()



#======== Body
'''
###### Mexican hat potential
nbStep = 100
nMDsteps = 100000
nEQsteps = 10000
t_damp = 100
T = 100
dt = 0.5

mmff.setVerbosity(verbosity=3)
mmff.init( xyz_name="data/xyz/pyridine", bMMFF=True )
collectiveVariable = np.array([0], dtype=np.int32)
E = mmff.compute_Free_energy(0.5, 4.0, collectiveVariable, nbStep=nbStep, nMDsteps=nMDsteps, nEQsteps=nEQsteps, t_damp=t_damp, T=T, dt=dt)
print("E=", E)
'''

'''
###### Three body problem
mmff.init( xyz_name="data/xyz/H2O", bMMFF=True, nPBC=(2, 2, 2))
colectiveVariable = np.array([0], dtype=np.int32)
nbStep = 100000         # JE -> nProdSteps
nMDsteps = 100         # JE -> nrealization
nEQsteps = 10000
t_damp = 100
T = 300
dt = 0.01
lamda1 = 0.0
lamda2 = 5.0
mmff.setVerbosity(verbosity=3)
mmff.setSwitches(CheckInvariants=-1, PBC_nonBond=-1, PBC_evalAtom=-1, NonBonded=1, MMFF=1, doBonds=-1, Angles=-1, PiSigma=-1, PiPiI=-1, bNonBondNeighs=-1, bSurfAtoms=-1, bGridFF=-1, bTricubic=-1, bConstrZ=-1, bConstrains=-1)
E = mmff.compute_Free_energy(lamda1, lamda2, colectiveVariable, nbStep=nbStep, nMDsteps=nMDsteps, nEQsteps=nEQsteps, t_damp=t_damp, T=T, dt=dt)
print("E=", E)
'''

###### Entropic spring (call with run_TI.sh)
import argparse

# Create parser for command-line arguments
parser = argparse.ArgumentParser(description='Run simulation with specified parameters')
parser.add_argument('--nMDsteps', type=int, help='Number of MD steps')
parser.add_argument('--nEQsteps', type=int, help='Number of equilibration steps')
parser.add_argument('--t_damp', type=float, help='Damping time')
parser.add_argument('--T', type=float, help='Temperature')
parser.add_argument('--dt', type=float, help='Time step')
parser.add_argument('--nbSteps', type=int, help='Number of steps')
parser.add_argument('--N', type=int, help='System size')
parser.add_argument('--K', type=float, help='Spring constant')
parser.add_argument('--lamda1', type=float, help='Lambda 1')
parser.add_argument('--lamda2', type=float, help='Lambda 2')

# Parse arguments
args = parser.parse_args()

# Assign values from command line arguments
MY_N = args.N
MY_K = args.K 
MY_nbStep = args.nbSteps
MY_nMDsteps = args.nMDsteps
MY_nEQsteps = args.nEQsteps
MY_t_damp = args.t_damp
MY_T = args.T
MY_dt = args.dt
MY_lamda1 = args.lamda1
MY_lamda2 = args.lamda2

# xyz_name = "data/entropic_spring_"+str(MY_N)
# constr_name = "data/entropic_spring_"+str(MY_N)+".cons"
#colectiveVariable = np.array([0], dtype=np.int32)

xyz_name = "data/nucleobasis_AT"
constr_name = "data/nucleobasis_AT.cons"
colectiveVariable = np.array([1,2], dtype=np.int32)

mmff.setVerbosity( verbosity=3, idebug=1 )

mmff.init( xyz_name=xyz_name, constr_name=constr_name ,bMMFF=True)
mmff.setSwitches(CheckInvariants=-1, PBC_nonBond=-1, PBC_evalAtom=-1, NonBonded=-1, MMFF=1, doBonds=1, Angles=-1, PiSigma=-1, PiPiI=-1, bNonBondNeighs=-1, bSurfAtoms=-1, bGridFF=-1, bTricubic=-1, bConstrZ=-1, bConstrains=-1)

E = mmff.compute_Free_energy(MY_lamda1, MY_lamda2, colectiveVariable, nbStep=MY_nbStep, nMDsteps=MY_nMDsteps, nEQsteps=MY_nEQsteps, t_damp=MY_t_damp, T=MY_T, dt=MY_dt)
print("E=", E)
exit()

natoms = 50
N = 100

#mmff.init( xyz_name="data/nHexadecan_dicarboxylic", bMMFF=True  )
#     


#mmff.init( xyz_name="data/xyz/pyridine", surf_name="data/NaCl_1x1_L2" )    
#mmff.init( xyz_name="data/xyz/nHexadecan_dicarboxylic", bMMFF=True  )     
# mmff.init( xyz_name="data/xyz/O", surf_name="data/xyz/NaCl_1x1_L3" )  
#mmff.init( xyz_name="data/xyz/H2O", surf_name="data/xyz/NaCl_1x1_L3" )    
# mmff.init( xyz_name="data/xyz/PTCDA", surf_name="data/xyz/NaCl_8x8_L3" )    
# mmff.init( xyz_name="data/xyz/PTCDA", surf_name="data/xyz/NaCl_8x8_L3_copy" )

# mmff.init( xyz_name="data/xyz/C.iz0", surf_name="data/xyz/C.iz0" )
# mmff.init( xyz_name="data/xyz/C.iz0", surf_name="data/xyz/Na.iz0" )
# mmff.init( xyz_name="data/xyz/C.iz0", surf_name="data/xyz/Cl.iz0" )
# mmff.init( xyz_name="data/xyz/C.iz0", surf_name="data/xyz/NaCl.iz0" )
# mmff.init( xyz_name="data/xyz/C_Cl_iz0", surf_name="data/xyz/NaCl.iz0" )

# mmff.init( xyz_name="data/xyz/Na+Cl-", surf_name="data/xyz/NaCl_8x8_L3_Coulumb" )

# mmff.init( xyz_name="data/xyz/PTCDA", surf_name="data/xyz/NaCl_coulomb.iz0" )

mmff.init( xyz_name="data/xyz/PTCDA_paolo", surf_name="data/xyz/NaCl.ptcda" )


# mmff.init( xyz_name="data/xyz/molNaCl_Na.iz0", surf_name="data/xyz/NaCl_coulomb.iz0" )
# mmff.init( xyz_name="data/xyz/molNaCl_Cl.iz0", surf_name="data/xyz/NaCl_coulomb.iz0" )
# mmff.init( xyz_name="data/xyz/molNaCl.ix0.iy0", surf_name="data/xyz/NaCl_coulomb.iz0" )

# mmff.init( xyz_name="data/xyz/2_atom_NaCl_mol", surf_name="data/xyz/2_atom_NaCl" )




'''
Debug :: REQs = [[1.4915     0.03606831 0.         0.        ]] For Na single atom substrate  Na.iz0.xyz
Debug :: REQs = [[1.9735     0.09921518 0.         0.        ]] For Cl single atom substrate  Cl.iz0.xyz


MolWorld_sp3::scan_rigid()[ia=0] pos(  0.0000,  0.0000,  7.1569) REQ(  1.9255,      0.06747763,  0.0000,  0.0000) PLQd(     21.77017471,      1.21202304,      0.00000000,      0.00000000)
MolWorld_sp3::scan_rigid()[ia=0] pos(  0.0000,  0.0000,  7.1569) REQ(  1.9255,      0.06747763,  0.0000,  0.0000) PLQd(     21.77017471,      1.21202304,      0.00000000,      0.00000000)

'''

mmff.getBuffs()

#print( "ffflags ", mmff.ffflags )

mmff.setSwitches( NonBonded=-1, MMFF=-1, SurfAtoms=0, GridFF=1 )


mmff.PLQs[:,0] = 0.0  # delete Pauli
mmff.PLQs[:,1] = 0.0  # delete London
# mmff.PLQs[:,2 ] = 0.0 # delete Coulomb (charges)

# scanPlot( nscan=10, span=(0.0,4.0), dir=(1.0,0.0,0.0), p0=(0.0,0.0,0.0),  label="E_x", saveFig="E_x_scan.png", saveData="E_x_scan.dat")
# scanPlot( nscan=1000, span=(0.0,4.0), dir=(0.0,1.0,0.0), p0=(0.0,0.0,0.0),  label="E_y", saveFig="E_y_scan.png", saveData="E_y_scan.dat" )

# scanPlot( nscan=100, span=(3.3,8.3), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="Atop_C", saveFig="E_z_scan_Atop_C.png", saveData="E_z_scan_Atop_C.dat" )
# scanPlot( nscan=100, span=(2.9,8.3), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="Atop_C on Na", saveFig="E_z_scan_Atop_Na.png", saveData="E_z_scan_Atop_Na.dat" )
# scanPlot( nscan=100, span=(3.3,8.3), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="Atop_C on Cl ", saveFig="E_z_scan_Atop_Cl.png", saveData="E_z_scan_Atop_Cl.dat" )
# scanPlot( nscan=100, span=(0.8,6.3), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="Atop_Na", saveFig="E_z_scan_Atop_Na_SYM.png", saveData="E_z_scan_Atop_Na_SYM.dat" )
# scanPlot( nscan=100, span=(2.3,8.3), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="Atop_Na", saveFig="E_z_scan_Atop_Na_Ewald.png", saveData="E_z_scan_Atop_Na_Ewald.dat" )
# scanPlot( nscan=100, span=(7.9,18.3), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="Atop_Na", saveFig="E_z_scan_Atop_Na_Ewald.png", saveData="E_z_scan_Atop_Na_Ewald.dat" )
# scanPlot( nscan=100, span=(3.1,8), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="Atop_Cl", saveFig="E_z_scan_Atop_Cl_No_SYM.png", saveData="E_z_scan_Atop_Cl_No_SYM.dat" )

# scanPlot( nscan=100, span=(3,10), dir=(0.0,0.0,1.0), p0=(0.0,0.0,5.656854), label="PTCDA z0=0", saveFig="E_z_scan_PTCDA_Ewald.png", saveData="E_z_scan_PTCDA_Ewald.dat" )

# scanPlot( nscan=100, span=(3,10), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0), label="PTCDA", saveFig="E_z_scan_PTCDA_Ewald_new_Direct.png", saveData="E_z_scan_PTCDA_Ewald_new_Direct.dat" )
scanPlot( nscan=100, span=(3,10), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0), label="PTCDA", saveFig="E_z_scan_PTCDA_Ewald_new_trial.png", saveData="E_z_scan_PTCDA_Ewald_new_trial.dat" )
# scanPlot( nscan=100, span=(3,10), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0), label="PTCDA", saveFig="E_z_scan_PTCDA_Morse_new.png", saveData="E_z_scan_PTCDA_Morse_new.dat" )

# scanPlot( nscan=101, span=(3,8.3), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="Atop_Na", saveFig="E_z_scan_Atop_Na_Ewald.png", saveData="E_z_scan_Atop_Na_Ewald.dat" )
# scanPlot( nscan=101, span=(3,8.3), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="Atop_Cl", saveFig="E_z_scan_Atop_Cl_Ewald.png", saveData="E_z_scan_Atop_Cl_Ewald.dat" )

# scanPlot( nscan=100, span=(3.1,8.1), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="Atop_Na", saveFig="E_z_scan_Atop_Na_Morse.png", saveData="E_z_scan_Atop_Na_Morse.dat" )
# scanPlot( nscan=100, span=(2.5,8.5), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="Atop_Cl", saveFig="E_z_scan_Atop_Cl_Morse.png", saveData="E_z_scan_Atop_Cl_Morse.dat" )

# scanPlot( nscan=100, span=(0,20), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="Atop_Na", saveFig="E_z_scan_2atom_NaCl_Ewald.png", saveData="E_z_scan_2atom_NaCl_Ewald.dat" )

# # For add pos0+ adjust(1.5) -4.9,4.5, 
# 2.35,8.5  
# # For sub pos0 only 8,12.157
#add pos0 only -3.5,4.5
# 1.9255
##################
# scanPlot2D(nscan1=41, nscan2=41,
#            span1=(0.0, 4.1), span2=(0.0, 4.1),
#            dir1=(1.0,0.0,0.0), dir2=(0.0,1.0,0.0),
#            p0=(0.0,0.0,3.4), label="E_xy for z=3.4",
#            saveFig="E_xy_scan_mol_NaCl.png", saveData="E_xy_scan_mol_NaCl.dat")

# scanPlot2D(nscan1=41, nscan2=61,
#            span1=(0.0, 4.1), span2=(2.3, 8.3),
#            dir1=(1.0,0.0,0.0), dir2=(0.0,0.0,1.0),
#            p0=(0.0,0.0,0), label="E_xz for Y=3.4",
#            saveFig="E_xz_scan_C.png", saveData="E_xz_scan_C.dat")

# scanPlot2D(nscan1=41, nscan2=61,
#            span1=(0.0, 4.1), span2=(2.3, 8.3),
#            dir1=(0.0,1.0,0.0), dir2=(0.0,0.0,1.0),
#            p0=(0.0,0.0,0), label="E_xz for Y=3.4",
#            saveFig="E_yz_scan_C.png", saveData="E_yz_scan_C.dat")           

# plt.legend()
# plt.grid()
# plt.show()
# print(mmff.PLQs.shape)



# Initialize the system. Make sure you call mmff.init() and then mmff.getBuffs()
# mmff.setVerbosity(verbosity=1, idebug=1)
# mmff.init(xyz_name="data/xyz/PTCDA", surf_name="data/xyz/NaCl_8x8_L3_copy")
# mmff.getBuffs()

# (Optionally) set up any switches, etc.
# mmff.setSwitches(NonBonded=-1, MMFF=-1, SurfAtoms=0, GridFF=1)

# Now, call our new saveXSF_geometry function.
# We assume that: 
#   - The grid data is stored in an array registered in the C++ global dictionary with key "BsplinePaul_pbc"
#   - mmff.natoms holds the total number of atoms of the substrate plus PTCDA.
#   - mmff.AtomType is a pointer or array of atom types (or use mmff.params.atypes if that’s your interface)
#   - mmff.fapos holds the geometry (atomic positions)

# It might be that the name "BsplinePaul_pbc" is not registered; adjust the key if needed.
# grid_array = mmff.getArrayPointer("BsplinePaul_pbc")
# if grid_array is None:
#     raise RuntimeError('Grid array "BsplinePaul_pbc" not found. Ensure the grid is built.')

# # Save to XSF file.
# ret = mmff.saveXSF_geometry("rigid_scan.xsf", grid_array, pitch=1, offset=0, 
#                               natoms=mmff.natoms,
#                               atypes=mmff.params.atypes if hasattr(mmff, "params") else None,
#                               fapos=mmff.fapos, bPrimCoord=True)
# if ret != 0:
#     print("XSF file saved successfully!")
# else:
#     print("Error saving XSF file.")

# # Optionally, run the dynamics or further analysis
# mmff.run()



# Save the grid and geometry to an XSF file:
# mmff.saveXSF("rigid_scan.xsf", mmff.getArrayPointer("BsplinePaul_pbc"), 1, 0, mmff.natoms, mmff.AtomType, mmff.fapos)
# mmff.saveXSF("rigid_scan.xsf", mmff.getArrayPointer("BsplinePaul_pbc"))


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
# exit(0)


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


# #mmff.initWithSMILES( "C=C1NC#CC1CO" )
# mmff.initWithSMILES( "C=C" )
# mmff.getBuffs()
# mmff.relax(1000)
# mmff.plot()
# #mmff.plot_selection( mmff.selection[:nsel] )
# #mmff.plot_selection( [1,2,3] )


# plt.show()

