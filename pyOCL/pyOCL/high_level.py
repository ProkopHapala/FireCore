#!/usr/bin/python

import numpy as np
#import os
import sys
from . import oclfft as ocl

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

def convolve( A, B, C=None, iA=0, iB=1, iC=2 ):
    ocl.tryInitFFT( A.shape )
    if C is None:
        C = np.empty(A.shape)
        #C = np.zeros(A.shape)
    ocl.upload ( iA, A )   # ;plt.figure(); plt.imshow( arrA.real )
    ocl.upload ( iB, B )   # ;plt.figure(); plt.imshow( arrB.real )
    ocl.convolve( iA,iB,iC )
    ocl.download( iC, C)   # ;plt.figure(); plt.imshow( arrC.real )
    return C

def poisson( A, dcell, C=None, iA=0, iC=1 ):
    ocl.tryInitFFT( A.shape )     ;print( "DEBUG poisson 1 " )
    ocl.upload_d(  iA, A )        ;print( "DEBUG poisson 2 " )
    ocl.poisson(  iA,iC, dcell )  ;print( "DEBUG poisson 3 " )
    C = ocl.download( iC, C, Ns=A.shape ) ;print( "DEBUG poisson 4 " )
    return C

def projectAtoms__( atoms, acoefs):
    import matplotlib.pyplot as plt
    import time
    print( "# ========= TEST   projectAtoms__()  " )

    #initFireBall( atypes, atoms )

    ocl.init()
    Ns=(100,100,100)
    ocl.initFFT( Ns  )
    ocl.loadWfBasis( [1,6], Rcuts=[4.5,4.5] )
    ocl.initAtoms( len(atoms) )
    #initBasisTable( basis.shape[0], basis.shape[1], basis )

    ocl.setGridShape( )
    t0 = time.clock()
    ocl.projectAtoms( atoms, acoefs, 0 )
    ocl.saveToXsf( "test.xsf", 0 )
    arrA = np.zeros(Ns,dtype=np.csingle  )
    ocl.download ( 0, arrA )    
    t = time.clock()-t0; print( "projection time ", t )
    plt.figure(); 
    #print( arrA[10,10].real )
    plt.imshow( arrA[10].real ) 
    #plt.imshow( np.log( np.abs(arrA[10])) ) 
    plt.grid()
    plt.show(); 
    ocl.cleanup()