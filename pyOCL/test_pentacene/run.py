#!/usr/bin/python


# INTEROPERABLE GLOBAL VARIABLES
# https://gcc.gnu.org/onlinedocs/gfortran/Interoperable-Global-Variables.html#Interoperable-Global-Variables


import sys

sys.path.append("../")
from pyOCL import oclfft as ocl
from pyOCL import utils as oclu
from pyOCL import jobs
from pyOCL import atomicUtils as au


sys.path.append("../../")
import pyBall as pb
from   pyBall import FireCore as fc

atomPos,atomType,enames,qs = au.loadAtomsNP( "input.bas")

iMO = 1

fc.run_nonSCF( atomType, atomPos )

ngrid, dCell, lvs = fc.setupGrid()
#fc.orb2xsf(iMO); #exit()
#fc.dens2xsf(); #exit()

#wf = fc.getGridMO( iMO,ngrid=ngrid)
rho = fc.getGridDens( ngrid=ngrid)

dcell = [0.1,0.1,0.1,0.1]
V     = oclu.poisson( rho, dcell )

#ewfaux = fc.getGridMO( iMO,ngrid=ngrid)   ;print( "ewfaux.min(),ewfaux.max() ", ewfaux.min(),ewfaux.max() )
#sh = ewfaux.shape                         ;print( "ewfaux.shape ", sh )

'''
#i0orb  = jobs.countOrbs( atomType )           ;print("i0orb ", i0orb)  
#wfcoef = fc.get_wfcoef(norb=i0orb[-1])


#print( xyzs )

#jobs.Test_projectDens( )
'''
