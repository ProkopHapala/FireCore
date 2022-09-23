
import numpy as np
import os
import sys
from . import FireCore as fc
#from . import MMFF as mmff

def makeLinearScan_firecore( nstep, selection, d, apos, nmax_scf=200 ):
    forces = np.zeros(apos.shape)
    Etemp  = np.zeros(8)
    Es     = np.zeros(nstep)
    for i in range(nstep):
        #fc.SCF( positions, iforce=iforce, nmax_scf=nmax_scf )
        fc.evalForce( apos, forces=forces, nmax_scf=nmax_scf, Es=Etemp, ixyz=i )
        print( "makeLinearScan_firecore() step# ", i," E[eV]= ", Etemp[0] )
        Es[i] = Etemp[0]
        apos[selection,:] += d
    return Es

def makeRotMat( ang ):
    ca=np.cos(ang)
    sa=np.sin(ang)
    return np.array([
        [1, 0,  0],
        [0,ca,-sa],
        [0,sa, ca]])

def makeRotationScan_firecore( nstep, selection, rot, p0, apos, nmax_scf=200 ):
    p0=np.array(p0)
    #ax=np.array(ax)
    #up=np.
    forces = np.zeros(apos.shape)
    Etemp  = np.zeros(8)
    Es     = np.zeros(nstep)
    for i in range(nstep):
        #fc.SCF( positions, iforce=iforce, nmax_scf=nmax_scf )
        fc.evalForce( apos, forces=forces, nmax_scf=nmax_scf, Es=Etemp, ixyz=i )
        print( "makeRotationScan_firecore() step# ", i," E[eV]= ", Etemp[0] )
        Es[i] = Etemp[0]
        apos[selection,:]-=p0
        for ia in selection:
            apos[ia,:] = np.dot( rot, apos[ia,:] )
        apos[selection,:]+=p0
    return Es

