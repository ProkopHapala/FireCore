import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

sys.path.append("../../")
from pyBall import eFF as eff

with open("processXYZ.xyz", "w") as f: f.write("")
eff.processXYZ( "export/scan_data/distscan_H2O.xyz", bOutXYZ=True );


