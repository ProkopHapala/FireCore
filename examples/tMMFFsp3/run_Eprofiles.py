import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import plotUtils   as plu
from pyBall import MMFFsp3     as mmff

rad2ang = 180.0/np.pi

# ================ Functions

def sample_PiPi():
    angs = np.linspace(0.0, np.pi, 100 )
    Es,Fs = mmff.sample_evalPiAling(angs)
    plt.plot(angs*rad2ang, Es )
    plt.plot(angs*rad2ang, Fs )
    plt.grid()

def H2O_angle():
    mmff.buildMolecule_xyz( xyz_name="data/H2O", bEpairs=False, fAutoCharges=-1 )
    mmff.makeFFs()
    mmff.getBuffs()   

    theta0 = mmff.measureAngle(0,1,2)    ;print( theta0 )
    n=100
    theta = np.arange( 0.0, np.pi, np.pi/n )-theta0
    Es    = mmff.scanRotation_ax( np.pi, n, [1], [0,0,0],  [0.,0.,1.], trjName="H2O_angle.xyz" )
    #plt.plot(theta, Es); plt.xlabel( "theta[rad]" ); plt.ylabel( "Energy[eV]" ); plt.grid()
    plt.plot(theta*rad2ang, Es); plt.xlabel( "theta[˚]" ); plt.ylabel( "Energy[eV]" ); plt.grid()
    plt.savefig( "H2O_angle.png" )  


def C2H4_torsion():
    mmff.buildMolecule_xyz( xyz_name="data/C2H4", bEpairs=False, fAutoCharges=-1 )
    mmff.makeFFs()
    mmff.getBuffs() 

    theta0 = mmff.measureAnglePiPi(0,1)  ;print( theta0 )
    n=100
    theta = np.arange( 0.0, 2*np.pi, 2*np.pi/n )-theta0
    #Es    = mmff.scanRotation_ax( 2*np.pi, n, [1,4,5], [0,0,0],  [0.,1.,0.], trjName="C2H4_torsion.xyz" )
    Es    = mmff.scanRotation_ax( 2*np.pi, n, [1], [0,0,0],  [0.,1.,0.], trjName="C2H4_torsion.xyz" )
    
    plt.plot(theta*rad2ang, Es); plt.xlabel( "theta[˚]" ); plt.ylabel( "Energy[eV]" ); plt.grid()
    plt.savefig( "C2H4_torsion.png" )  

def NH3_inversion():
    #mmff.buildMolecule_xyz( xyz_name="data/NH3", bEpairs=False, fAutoCharges=-1 )
    mmff.buildMolecule_xyz( xyz_name="data/NH3_R", bEpairs=False, fAutoCharges=-1 )
    mmff.makeFFs()
    mmff.getBuffs()

    #mmff.setTrjName("NH3_relax.xyz",1); nsteps = mmff.run( 1000 )    # relax     

    #Es  = mmff.scanTranslation_ax( 2*np.pi, n, [1], [0,0,0],  [0.,1.,0.],  )
    angs = np.arange( -np.pi/4, np.pi/4, 0.5*np.pi/100 )
    Es = mmff.scanAngleToAxis_ax( angs, [1,2,3], r=1.019, p0=[0.,0.,0.], ax=[0.,0.,1.], R=0, Es=None, trjName="NH3_inversion.xyz" )
    #mmff.scanTranslation_ax( 10, Es=None, trjName=None, sel=[0], vec=[0.0,0.0,-0.1], trjName="NH3_inversion.xyz", bAddjustCaps=False )

    plt.plot(angs*rad2ang, Es); plt.xlabel( "theta[˚]" ); plt.ylabel( "Energy[eV]" ); plt.grid()
    plt.savefig( "NH3_inversion.png" ) 
    

#def do_inversion( name ):
#    #mmff.buildMolecule_xyz( xyz_name="data/NH3", bEpairs=False, fAutoCharges=-1 )
#    mmff.buildMolecule_xyz( xyz_name="data/"+name, bEpairs=False, fAutoCharges=-1 )
#    mmff.makeFFs()
#    mmff.getBuffs()
#
#    #mmff.setTrjName("NH3_relax.xyz",1); nsteps = mmff.run( 1000 )    # relax     
#
#    #Es  = mmff.scanTranslation_ax( 2*np.pi, n, [1], [0,0,0],  [0.,1.,0.],  )
#    angs = np.arange( -np.pi/4, np.pi/4, 0.5*np.pi/100 )
#    Es = mmff.scanAngleToAxis_ax( angs, [1,2,3], r=1.019, p0=[0.,0.,0.], ax=[0.,0.,1.], R=0, Es=None, trjName=name+"_inversion.xyz" )
#    #mmff.scanTranslation_ax( 10, Es=None, trjName=None, sel=[0], vec=[0.0,0.0,-0.1], trjName="NH3_inversion.xyz", bAddjustCaps=False )
#
#    plt.plot(angs*rad2ang, Es); plt.xlabel( "theta[˚]" ); plt.ylabel( "Energy[eV]" ); plt.grid()
#    plt.savefig( name+"_inversion.png" ) 


#======== Body

#name = sys.argv[1]

#mmff.setVerbosity( verbosity=1, idebug=0 )
mmff.setVerbosity( verbosity=0, idebug=0 )


#mmff.setSwitches(doAngles=0, doPiPiT=0, doPiSigma=0, doPiPiI=0, doBonded=0, PBC=0, CheckInvariants=0)
#mmff.setSwitches( NonBonded=-1, PiSigma=-1, Angles=-1 )
#mmff.setSwitches( NonBonded=-1, PiSigma=+1, Angles=-1 )
#mmff.setSwitches( NonBonded=-1, PiSigma=-1, Angles=+1 )
mmff.printSwitches()

#------ Long Initialization
mmff.initParams()

#sample_PiPi()
#H2O_angle()
C2H4_torsion()
#NH3_inversion()

#do_inversion(name)

#names = ['']

#for name in names:
#    do_inversion(name)

mmff.printAtomTypes()
mmff.printAtomParams()
mmff.printSwitches()

plt.show()






