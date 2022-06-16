#!/usr/bin/python

import numpy as np
#import os
import sys
#from . import oclfft as ocl

def countOrbs(atypes):
    norbs = [ 1 if (x == 1) else 4 for x in atypes ]
    return np.cumsum(norbs)

def convCoefs( atypes, oCs, oatoms, typeDict ):
    na = len( atypes)
    norbs = [ 1 if (x == 1) else 4 for x in atypes ]
    norb = sum(norbs)
    atoms = np.zeros( (na,4), dtype=np.float32)
    coefs = np.zeros( (na,4), dtype=np.float32)
    #print( "atoms.shape ", atoms.shape, "oatoms.shape ", oatoms.shape, " oCs.shape ", oCs.shape  )
    atoms[:,:3] = oatoms[:,:3]
    atoms[:,3]=0.1
    io=0
    for ia,no in enumerate(norbs):
        coefs[ia,3]=oCs[io]; io+=1
        if(no>1):
            # yzx : Fireball p-orbitals are in order  y,z,x;   see https://nanosurf.fzu.cz/wiki/doku.php?id=fireball
            coefs[ia,0]=oCs[io+2]
            coefs[ia,1]=oCs[io+0]
            coefs[ia,2]=oCs[io+1]
            io+=3
        atoms[ia,3] += typeDict[atypes[ia]] + 0.1
    #print( "atoms ", atoms )
    #print( "coefs ", coefs )
    print( " atomTypes CPU ", atypes[ia]  )
    print( " atomTypes GPU ", atoms[:,3]  )
    return atoms, coefs