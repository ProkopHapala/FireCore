
import sys
import numpy as np
sys.path.append('../../')

from pyBall.AtomicSystem import AtomicSystem

#backbone  = AtomicSystem(fname='./backbones/PNA_sat-cis-.mol2')
backbone  = AtomicSystem(fname='./backbones/diacetylene-.mol2')
endgroup1 = AtomicSystem(fname='./endgroups/guanine-.mol2')
endgroup2 = AtomicSystem(fname='./endgroups/cytosine-.mol2')

#backbone.save_mol2('backbone.mol2')
#endgroup1.save_mol2('endgroup1.mol2')

backbone.attach_group_by_marker(endgroup1, markerX="S_3", markerY="F",  _0=1)

print( "saving subs 1 to files: backbone_subs_1.mol2, backbone_subs_1.xyz" )
backbone.save_mol2('backbone_subs_1.mol2')
backbone.saveXYZ('backbone_subs_1.xyz')

backbone.attach_group_by_marker(endgroup2, markerX="Se" , markerY="Cl", _0=1)

print( "saving subs 2 to files: backbone_subs_2.mol2, backbone_subs_2.xyz" )
backbone.save_mol2('backbone_subs_2.mol2')
backbone.saveXYZ('backbone_subs_2.xyz')