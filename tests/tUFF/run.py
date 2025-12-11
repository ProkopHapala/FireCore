import sys
import os
sys.path.append("../../../../")
from pyBall import MMFF as mmff

mmff.init( xyz_name="PTCDA", bUFF=True, bSimple=True )

mmff.test_UFF( int(sys.argv[1]) )

exit()
