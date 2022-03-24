#!/usr/bin/python

import sys, os
import numpy as np
import matplotlib.pyplot as plt

up1=os.path.dirname(sys.path[0])   # ../
up2=os.path.dirname(up1 ) ; print( "up1 :", up1 )  # ../../
up3=os.path.dirname(up2)  ; print( "up2 :", up2 )  # ../../../

sys.path.append(up2)
sys.path.append(up3)

#import pyProbeParticle as ppm
import pyBall.FireCore as fc
from pyProbeParticle import atomicUtils as au
from pyProbeParticle import fieldFFT    as fFFT

xyzs,Zs,enames,qs = au.loadAtomsNP('Coronene.xyz')
natoms=len(Zs)

# =======   FFT For  PyCUDA and PyOpenCL    https://pythonhosted.org/pyfft/

# Fx,Fy,Fz,E = fFFT.potential2forces_mem( rho1, lvec1, nDim1, rho=rho2, doForce=True, doPot=True, deleteV=True )

print ("atomType ", Zs)
print ("atomPos  ", xyzs)
fc.preinit()
fc.init( natoms,  Zs, xyzs )
# =========== Electron Density
fc.assembleH( xyzs )
fc.solveH()
sigma= fc.updateCharges() ; print( sigma )
ngrid, dCell, lvs = fc.setupGrid()
ewfaux = fc.getGridDens( ngrid=ngrid )

print( "lvs", lvs )

Fx,Fy,Fz,E, rhoTip = fFFT.potential2forces_mem( ewfaux, lvs, ngrid, rho=None, multipole={'s':1.0}, doForce=True, doPot=True, deleteV=True )
ewfaux = E
#ewfaux = Fz
#ewfaux = Fx
#ewfaux = rhoTip

print( ewfaux.min(),ewfaux.max() )
import matplotlib.pyplot as plt
sh = ewfaux.shape
plt.figure(); plt.imshow( ewfaux[ sh[0]//2+5,:,: ] )
plt.figure(); plt.imshow( ewfaux[ sh[0]//2  ,:,: ] )
plt.figure(); plt.imshow( ewfaux[ sh[0]//2-5,:,: ] )
plt.show()


 