
import sys
import numpy as np
sys.path.append('../../')
from pyBall.AtomicSystem import AtomicSystem


monomers ={
#"T":  ("./backbones/TNA.mol2"          ,(5,6  )  ),
#"P":  ("./backbones/PNA_sat-cis-.mol2" ,(2,23 ) ),
#"D":  ("./backbones/diacetylene-.mol2" ,(1,4  ) ),
"D":  ("./DANA_CG.mol2" ,(1,4  ) ),
"T":  ("./TNA_CG-.mol2" ,(6,5  ) ),
"P":  ("./PNA_CG.mol2" ,(2,19  ) ),
}

#sequence = "TPTTDD"
sequence = "PPPP"
#sequence = "DDDD"
#sequence = "TTTT"

n0=0
for i,letter in enumerate(sequence):
    print( "Adding %s" % letter )
    if i==0:
        orec = monomers[letter]
        A    = AtomicSystem(fname=orec[0])
        pos  = np.zeros(3)
    else:
        rec  = monomers[letter]
        B    = AtomicSystem(fname=rec[0])
        pos += B.lvec[1]
        n0_ = len( A.atypes )
        A.addSystems(B, pos=pos, added_bonds=[ (n0+rec[1][0],orec[1][1])], _0=1 )
        orec = rec
        n0 = n0_

A.saveXYZ  (f'sequence_{sequence}.xyz')
A.save_mol2(f'sequence_{sequence}.mol2')

