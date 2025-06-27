#!/usr/bin/env python3
import sys
import os
import numpy as np
np.set_printoptions(linewidth=np.inf, threshold=np.inf)

sys.path.append("../../")
from pyBall import eFF as eff

eff.setVerbosity(2)


# Initialize eFF (if needed)

# Set atom parameters (mode=2 for double8 format)
atomParams2 = np.array([
    # Z_nuc, R_eff, Zcore_eff,   PA,        PB,        PC,        PD,        PE
    [0.0,   1.0,      0.0,      0.0,       0.0,       0.0,       0.0,       0.0], # 0: Dummy
    [1.0,   0.0,      0.0,      0.0,       0.0,       0.0,       0.0,       0.0], # 1: H
    [2.0,   0.1,      2.0,      0.0,       0.0,       0.0,       0.0,       0.0], # 2: He
    [3.0,   0.1,      2.0,      0.0,       0.0,       0.0,       0.0,       0.0], # 3: Li
    [4.0,   0.1,      2.0,      0.0,       0.0,       0.0,       0.0,       0.0], # 4: Be
    [5.0,   0.1,      2.0,      0.0,       0.0,       0.0,       0.0,       0.0], # 5: B
    [6.0,   0.621427, 2.0,     22.721015, 0.728733,  1.103199, 17.695345, 6.693621], # 6: C
    [7.0,   0.0,      2.0,      0.0,       0.0,       0.0,       0.0,       0.0], # 7: N
    [8.0,   0.167813, 2.0,     25.080199, 0.331574,  1.276183, 12.910142, 3.189333], # 8: O
    [9.0,   0.3,      2.0,      0.0,       0.0,       0.0,       0.0,       0.0]  # 9: F
], dtype=np.float64)

eff.setAtomParams(atomParams2, mode=1)

# Process XYZ file
#xyz_file = os.path.join(os.path.dirname(__file__), 'export/scan_data/distscan_Oe.xyz')
#outEs, apos, epos = eff.processXYZ_e( xyz_file, nstepMax=0,  dt=0.001, Fconv=1e-3, optAlg=2, bOutputs=(1,1,1) )

trjname = 'relaxation.xyz'
print("trjname set")
eff.setTrjName( trjname, savePerNsteps=10)
#outEs, apos, epos = eff.processXYZ_e( "H2O_pairs.xyz", nstepMax=10000,  dt=0.001, Fconv=1e-3, optAlg=2, bOutputs=(1,1,1) )
#outEs, apos, epos = eff.processXYZ_e( "H2O_pairs_fc.xyz", nstepMax=1000,  dt=0.01, Fconv=1e-3, optAlg=2, bOutputs=(1,1,1) )
print("trj name set2 ")
outEs, apos, epos = eff.processXYZ_e( "H2O_spins.xyz", nstepMax=1000,  dt=0.001, Fconv=1e-3, optAlg=2, bOutputs=(1,1,1) )
#outEs, apos, epos = eff.processXYZ_e( "H2O_spins_fc.xyz", nstepMax=1000,  dt=0.01, Fconv=1e-3, optAlg=2, bOutputs=(1,1,1) )

print("==== Results:")
print("Energies:", outEs)
print("Atom positions:\n", apos)
print("Electron positions and sizes:\n", epos)

# ----- Visualization
from xyz_view_new import MolViewer
from pyBall import atomicUtils as au
trj = au.load_xyz_movie(trjname)
trj = au.trj_to_ename(trj)
trj = au.trj_fill_radius(trj, bVdw=True, rFactor=0.001, rmin=0.05) # Use a larger rFactor
#trj = au.trj_fill_radius(trj, bVdw=False, rFactor=1.0)
#print( "trj.enames", trj[0])
MolViewer.launch(trj=trj)
