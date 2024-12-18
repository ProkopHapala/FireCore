import sys
import os

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import FireCore as fc

# ======== Body

# Get molecule name from command-line argument
if len(sys.argv) < 3:
    print("Usage: python script.py <input_path.xyz> <output_dir>")
    sys.exit(1)

name     = sys.argv[1]
outdir   = sys.argv[2]
basename = os.path.basename(name)
mol = au.AtomicSystem(name)
#print("Before relaxation: apos=\n", mol.apos)
fc.setVerbosity(0)
fc.initialize(atomType=mol.atypes, atomPos=mol.apos, verbosity=3)
fc.relax(mol.apos, nstepf=1000, nmax_scf=100)
mol.saveXYZ( os.path.join(outdir, basename ) )