import numpy as np
import sys
sys.path.append("../../")
from pyBall import MMFF as mmff

# Verbose to see trajectory debug prints
#mmff.setVerbosity(verbosity=1, idebug=1)
mmff.setVerbosity(verbosity=0, idebug=0)

mmff.setSwitches(NonBonded=1, MMFF=1)

# 1) Init molecule only (no GridFF substrate)
# /home/prokop/git/FireCore/cpp/common_resources/mol/xylitol.mol
mmff.init(xyz_name="data/mol/xylitol.mol", surf_name=None, bMMFF=True)
#mmff.init(xyz_name="data/mol/xylitol.mol2", surf_name=None, bMMFF=True)
#mmff.init(xyz_name="data/xyz/H2O", surf_name=None, bMMFF=True)

# 2) Define plane (optional if you’re happy with z-plane through origin)
mmff.setSurfFlatPlane(pos0=(0.0, 0.0, -10.0), normal=(0.0, 0.0, 1.0))

# 3) Choose potential and parameters
# Hamaker LJ 9-3 (mode=1): REQ.x = z0 (effective R), REQ.y = epsilon
mmff.setSurfFlatParams(mode=1, REQ=(1.5, 1.000, 0.0, 0.0), K=1.6)
# For Morse instead, use mode=2 and set K accordingly:
# mmff.setSurfFlatParams(mode=2, REQ=(1.5, 0.010, 0.0, 0.0), K=1.6)

mmff.setSwitches( GridFF=0, NonBonded=1, MMFF=1)

mmff.setTrjName("/home/prokop/git/FireCore/tests/tMMFF/trj_flat_surf.xyz", savePerNsteps=1, bDel=True, nPBC=(1,1,1))

# 5) Optionally relax
mmff.run(nstepMax=5000, dt=0.05, Fconv=1e-3, ialg=2, omp=True)

# 6) Save final geometry
#mmff.saveXYZ("h2o_flat_relaxed.xyz", bN=True, bComment=True)
mmff.saveXYZ("relaxed_flat_surf.xyz", "relaxed flat surf", 1)