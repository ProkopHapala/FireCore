import sys
import os
sys.path.append("/home/niko/work/GRIDFF/FireCore/")
from pyBall import MMFF as mmff

mmff.init( xyz_name="mol", surf_name="sub", bUFF=True, bSimple=True, gridStep=-1 )
print("mmff.init DONE")
#exit()

mmff.test_UFF( int(sys.argv[1]) )
print("mmff.test_UFF DONE")

exit()
