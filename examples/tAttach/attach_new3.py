
import sys
import numpy as np
sys.path.append('../../')
from pyBall.AtomicSystem import AtomicSystem
# (Make sure the free functions are imported as well.)


backbones=[
    #("porph2s_",'./porphironoids/porphirin_subs2.mol2'),
    #("porph3s_",'./porphironoids/porphirin_3sites_subs.mol2'),
    #("diphenatroline_CO_3s",'./porphironoids/diphenatroline_CO_subs.mol2'),
    #("hex3s_",'./porphironoids/tisite_subs.mol2')
    #("PorhQuad_","./porphironoids/PorhQuad_4SeCl.mol2")
    ("PorhQuad_","./porphironoids/PorhQuad_4SeCl_30deg.mol2")
]

endgroups=[
#("C",'./endgroups/cytosine-.mol2'),
#("T",'./endgroups/thymine-.mol2'),
#("A",'./endgroups/adenine-.mol2'),
("G",'./endgroups/guanine-SeCl.mol2'),
]



for backbone_fname in backbones:
    bname, bfname = backbone_fname
    for endgroup in endgroups:
        ename, efname = endgroup
        outname = bname+ename
        backbone  = AtomicSystem(fname=bfname)
        endgroup1 = AtomicSystem(fname=efname)

        backbone.preinitialize_atomic_properties()
        endgroup1.preinitialize_atomic_properties()

        backbone.attach_group_by_marker(endgroup1, markerX="Se", markerY="Cl", _0=1)

        
        backbone.save_mol2(outname+'.mol2')
        backbone.saveXYZ(outname+'.xyz')



#backbone.save_mol('backbone_with_endgroup.mol')