

# import sys
# import os

# sys.path.append('../../')
# from pyBall             import plotUtils   as plu
# from pyBall             import atomicUtils as au
# from pyBall.atomicUtils import AtomicSystem

# import numpy as np
# import matplotlib.pyplot as plt
# #from functools import partial



# B = AtomicSystem(fname='B2.mol2' )
# G = AtomicSystem(fname='G1.mol2' )

# print("B.print()"); B.print()
# print("G.print()"); G.print()

# B.attach_group_by_marker( G, markerX="Se", markerY="F", forward_default=np.array([1.0, 0.0, 0.0]), _0=1, pre="X")
# B.save_mol2('backbone_with_G.mol2')


import sys
import numpy as np
sys.path.append('../../')
from pyBall.atomicUtils import AtomicSystem
# (Make sure the free functions are imported as well.)

# Load systems from mol2 files.
backbone = AtomicSystem(fname='B2.mol2')
endgroup = AtomicSystem(fname='G1.mol2')

backbone.preinitialize_atomic_properties()
endgroup.preinitialize_atomic_properties()

print("Backbone:")
backbone.print()
print("Endgroup:")
endgroup.print()

# Attach endgroup copies to the backbone using marker elements.
# Here, the backbone is expected to contain marker atoms with element "Se" (attachment)
# and "F" (up direction), and the endgroup must contain exactly one such pair.
#backbone.attach_group_by_marker(endgroup, markerX="Se", markerY="F", forward_default=np.array([1.0, 0.0, 0.0]), _0=1, pre="X")

backbone.attach_group_by_marker(endgroup, markerX="Se", markerY="F", _0=1, pre="X")

# Save the modified backbone.
#backbone.save_mol2('backbone_with_endgroup.mol2')
backbone.saveXYZ( fname="backbone_with_endgroup.xyz" )