import sys
import numpy as np
import os

exec_path = os.getcwd()
sys.path.insert(0, exec_path)

data_dir = exec_path+"/cpp/common_resources/"

from pyBall import atomicUtils as au
from pyBall import MMFF as mmff

def numDeriv( xs, Es): 
    dx = xs[1]-xs[0]
    Fs = (Es[2:]-Es[:-2])/(2*dx)
    return -Fs





# Bonding interaction
def test_sample_evalBond():
    xs    = np.linspace(0.1,3.0,100)
    mmff.setVerbosity(-1)
    Es,Fs = mmff.sample_evalBond( xs)
    Fnum = numDeriv(xs,Es)
    assert np.all(np.abs(Fs[1:-1]-Fnum)<0.1), "!Check Forces::evalBond! sample_evalBond force calculation failed: analytical forces do not match numerical derivatives"

def test_sample_evalAtom():
    N = 10
    xs    = np.linspace(2.1,3.0,N)
    mmff.setVerbosity(-1)
    mmff.init( xyz_name=data_dir+"xyz/polymer-2_new", sElementTypes=data_dir+"ElementTypes.dat", sAtomTypes=data_dir+"AtomTypes.dat", sBondTypes=data_dir+"BondTypes.dat", sAngleTypes=data_dir+"AngleTypes.dat", sDihedralTypes=data_dir+"DihedralTypes.dat")
    Es,Fs = mmff.sample_evalAtom( xs, ia=1)
    mmff.clear()
    Fnum = numDeriv(xs,Es)
    assert np.all(np.abs(Fs[1:-1]-Fnum)<1.0), "!Check Forces::evalAtom! sample_evalAtom force calculation failed: analytical forces do not match numerical derivatives"



# Nonbonding interaction
def test_sample_getLJQH():
    xs    = np.linspace(0.8,3.0,1000)
    mmff.setVerbosity(-1)
    Es,Fs = mmff.sample_getLJQH( xs)
    Fnum = numDeriv(xs,Es)
    assert np.all(np.abs(Fs[1:-1]-Fnum)<0.1), "!Check Forces::getLJQH! sample_getLJQH force calculation failed: analytical forces do not match numerical derivatives"

def test_sample_evalLJQs_ng4_PBC_atom_omp():
    xs    = np.linspace(-3.0,-2.7,1000)
    mmff.setVerbosity(-1)
    mmff.init( xyz_name=data_dir+"xyz/pyridine", sElementTypes=data_dir+"ElementTypes.dat", sAtomTypes=data_dir+"AtomTypes.dat", sBondTypes=data_dir+"BondTypes.dat", sAngleTypes=data_dir+"AngleTypes.dat", sDihedralTypes=data_dir+"DihedralTypes.dat")
    Es,Fs = mmff.sample_evalLJQs_ng4_PBC_atom_omp( xs)
    mmff.clear()
    Fnum = numDeriv(xs,Es)
    assert np.all(np.abs(Fs[1:-1]-Fnum)<0.1), "!Check Forces::evalLJQs_ng4_PBC_atom_omp! sample_evalLJQs_ng4_PBC_atom_omp force calculation failed: analytical forces do not match numerical derivatives"

def test_sample_evalLJQs_ng4_atom_omp():
    xs    = np.linspace(-3.0,-2.7,1000)
    mmff.setVerbosity(-1)
    mmff.init( xyz_name=data_dir+"xyz/pyridine", sElementTypes=data_dir+"ElementTypes.dat", sAtomTypes=data_dir+"AtomTypes.dat", sBondTypes=data_dir+"BondTypes.dat", sAngleTypes=data_dir+"AngleTypes.dat", sDihedralTypes=data_dir+"DihedralTypes.dat")
    Es,Fs = mmff.sample_evalLJQs_ng4_atom_omp      ( xs)
    mmff.clear()
    Fnum = numDeriv(xs,Es)
    assert np.all(np.abs(Fs[1:-1]-Fnum)<0.1), "!Check Forces::evalLJQs_ng4_atom_omp! sample_evalLJQs_ng4_atom_omp force calculation failed: analytical forces do not match numerical derivatives"

def test_sample_evalLJQs_PBC_atom_omp():
    xs    = np.linspace(-3.0,-2.7,1000)
    mmff.setVerbosity(-1)
    mmff.init( xyz_name=data_dir+"xyz/pyridine", sElementTypes=data_dir+"ElementTypes.dat", sAtomTypes=data_dir+"AtomTypes.dat", sBondTypes=data_dir+"BondTypes.dat", sAngleTypes=data_dir+"AngleTypes.dat", sDihedralTypes=data_dir+"DihedralTypes.dat")
    Es,Fs = mmff.sample_evalLJQs_PBC_atom_omp      ( xs)
    mmff.clear()
    Fnum = numDeriv(xs,Es)
    assert np.all(np.abs(Fs[1:-1]-Fnum)<0.1), "!Check Forces::evalLJQs_PBC_atom_omp!"

def test_sample_evalLJQs_atom_omp():
    xs    = np.linspace(-3.0,-2.7,1000)
    mmff.setVerbosity(-1)
    mmff.init( xyz_name=data_dir+"xyz/pyridine", sElementTypes=data_dir+"ElementTypes.dat", sAtomTypes=data_dir+"AtomTypes.dat", sBondTypes=data_dir+"BondTypes.dat", sAngleTypes=data_dir+"AngleTypes.dat", sDihedralTypes=data_dir+"DihedralTypes.dat")
    Es,Fs = mmff.sample_evalLJQs_atom_omp      ( xs)
    mmff.clear()
    Fnum = numDeriv(xs,Es)
    assert np.all(np.abs(Fs[1:-1]-Fnum)<0.1), "!Check Forces::evalLJQs_atom_omp! sample_evalLJQs_atom_omp force calculation failed: analytical forces do not match numerical derivatives"





# Interaction with surface
def test_sample_evalMorsePBC_sym():
    xs    = np.linspace(0.0,3.0,100)
    #mmff.setVerbosity(-1)
    mmff.init( xyz_name=data_dir+"xyz/pyridine", surf_name=data_dir+"xyz/NaCl_1x1_L2", sElementTypes=data_dir+"ElementTypes.dat", sAtomTypes=data_dir+"AtomTypes.dat", sBondTypes=data_dir+"BondTypes.dat", sAngleTypes=data_dir+"AngleTypes.dat", sDihedralTypes=data_dir+"DihedralTypes.dat")
    mmff.setSwitches(bGridFF=-1)
    Es,Fs = mmff.sample_evalMorsePBC_sym( xs)
    Fnum = numDeriv(xs,Es)
    mmff.clear()
    assert np.all(np.abs(Fs[1:-1]-Fnum)<0.1), "!Check Forces::evalMorsePBC_sym! sample_evalMorsePBC_sym force calculation failed: analytical forces do not match numerical derivatives"



def test_sample_springbound():
    xs    = np.linspace(0.0,3.0,100)
    mmff.setVerbosity(-1)
    mmff.init( xyz_name=data_dir+"xyz/pyridine", surf_name=data_dir+"xyz/NaCl_1x1_L2", sElementTypes=data_dir+"ElementTypes.dat", sAtomTypes=data_dir+"AtomTypes.dat", sBondTypes=data_dir+"BondTypes.dat", sAngleTypes=data_dir+"AngleTypes.dat", sDihedralTypes=data_dir+"DihedralTypes.dat")
    Es,Fs = mmff.sample_springbound( xs)
    Fnum = numDeriv(xs,Es)
    mmff.clear()
    assert np.all(np.abs(Fs[1:-1]-Fnum)<0.1), "!Check Forces::springbound! sample_springbound force calculation failed: analytical forces do not match numerical derivatives"


# Constrains
def test_sample_applyConstr():
    xs    = np.linspace(0.0,3.0,100)
    mmff.setVerbosity(-1)
    mmff.init( xyz_name=data_dir+"xyz/dicarboxylic_acid", constr_name=data_dir+"dicarboxylic_acid.cons", sElementTypes=data_dir+"ElementTypes.dat", sAtomTypes=data_dir+"AtomTypes.dat", sBondTypes=data_dir+"BondTypes.dat", sAngleTypes=data_dir+"AngleTypes.dat", sDihedralTypes=data_dir+"DihedralTypes.dat")
    Es,Fs = mmff.sample_applyConstr( xs )
    Fnum = numDeriv(xs,Es)
    mmff.clear()
    assert np.all(np.abs(Fs[1:-1]-Fnum)<0.1), "!Check Forces::applyConstr! sample_applyConstr force calculation failed: analytical forces do not match numerical derivatives"



# SetSwitches+run_omp --> bonding
def test_run_omp_bonds():
    N = 1000
    xs    = np.linspace(-3.0,3.0,N)
    mmff.setVerbosity(-1)
    mmff.init( xyz_name=data_dir+"xyz/pyridine", sElementTypes=data_dir+"ElementTypes.dat", sAtomTypes=data_dir+"AtomTypes.dat", sBondTypes=data_dir+"BondTypes.dat", sAngleTypes=data_dir+"AngleTypes.dat", sDihedralTypes=data_dir+"DihedralTypes.dat")
    mmff.setSwitches(CheckInvariants=-1, PBC_nonBond=-1, PBC_evalAtom=-1, NonBonded=-1, MMFF=1, doBonds=1, Angles=-1, PiSigma=1, PiPiI=-1, bNonBondNeighs=1, bSurfAtoms=-1, bGridFF=-1, bTricubic=-1, bConstrZ=-1, bConstrains=-1, bMoving=-1)
    Es,Fs = mmff.sample_movementOfAtom(xs, ia=0)
    Fnum = numDeriv(xs,Es)
    mmff.clear()
    assert np.all(np.abs(Fs[1:-1]-Fnum)<0.1), "!Check Forces::run_omp_bonds or sample_movementOfAtom! run_omp_bonds force calculation failed: analytical forces do not match numerical derivatives"

def test_run_omp_bonds_angles():
    N = 1000
    xs    = np.linspace(-3.0,3.0,N)
    mmff.setVerbosity(-1)
    mmff.init( xyz_name=data_dir+"xyz/pyridine", sElementTypes=data_dir+"ElementTypes.dat", sAtomTypes=data_dir+"AtomTypes.dat", sBondTypes=data_dir+"BondTypes.dat", sAngleTypes=data_dir+"AngleTypes.dat", sDihedralTypes=data_dir+"DihedralTypes.dat")
    mmff.setSwitches(CheckInvariants=-1, PBC_nonBond=-1, PBC_evalAtom=-1, NonBonded=-1, MMFF=1, doBonds=-1, Angles=1, PiSigma=1, PiPiI=-1, bNonBondNeighs=1, bSurfAtoms=-1, bGridFF=-1, bTricubic=-1, bConstrZ=-1, bConstrains=-1, bMoving=-1)
    Es,Fs = mmff.sample_movementOfAtom(xs, ia=0)
    Fnum = numDeriv(xs,Es)
    mmff.clear()
    assert np.all(np.abs(Fs[1:-1]-Fnum)<0.1), "!Check Forces::run_omp_bonds_angles! run_omp_bonds_angles force calculation failed: analytical forces do not match numerical derivatives"

def test_run_omp_pisigma():
    N = 1000
    xs    = np.linspace(-3.0,3.0,N)
    mmff.setVerbosity(-1)
    mmff.init( xyz_name=data_dir+"xyz/pyridine", sElementTypes=data_dir+"ElementTypes.dat", sAtomTypes=data_dir+"AtomTypes.dat", sBondTypes=data_dir+"BondTypes.dat", sAngleTypes=data_dir+"AngleTypes.dat", sDihedralTypes=data_dir+"DihedralTypes.dat")
    mmff.setSwitches(CheckInvariants=-1, PBC_nonBond=-1, PBC_evalAtom=-1, NonBonded=-1, MMFF=1, doBonds=-1, Angles=-1, PiSigma=1, PiPiI=-1, bNonBondNeighs=1, bSurfAtoms=-1, bGridFF=-1, bTricubic=-1, bConstrZ=-1, bConstrains=-1, bMoving=-1)
    Es,Fs = mmff.sample_movementOfAtom(xs, ia=0)
    Fnum = numDeriv(xs,Es)
    mmff.clear()
    assert np.all(np.abs(Fs[1:-1]-Fnum)<0.01), "!Check Forces::run_omp_pisigma! run_omp_pisigma force calculation failed: analytical forces do not match numerical derivatives"

def test_run_omp_pipi():
    N = 1000
    xs    = np.linspace(-3.0,3.0,N)
    mmff.setVerbosity(-1)
    mmff.init( xyz_name=data_dir+"xyz/pyridine", sElementTypes=data_dir+"ElementTypes.dat", sAtomTypes=data_dir+"AtomTypes.dat", sBondTypes=data_dir+"BondTypes.dat", sAngleTypes=data_dir+"AngleTypes.dat", sDihedralTypes=data_dir+"DihedralTypes.dat")
    mmff.setSwitches(CheckInvariants=-1, PBC_nonBond=-1, PBC_evalAtom=-1, NonBonded=-1, MMFF=1, doBonds=-1, Angles=-1, PiSigma=-1, PiPiI=1, bNonBondNeighs=1, bSurfAtoms=-1, bGridFF=-1, bTricubic=-1, bConstrZ=-1, bConstrains=-1, bMoving=-1)
    Es,Fs = mmff.sample_movementOfAtom(xs, ia=0)
    Fnum = numDeriv(xs,Es)
    mmff.clear()
    assert np.all(np.abs(Fs[1:-1]-Fnum)<0.01), "!Check Forces::run_omp_pipi! run_omp_pipi force calculation failed: analytical forces do not match numerical derivatives"

def test_run_omp_bonds_angles_pisigma_pipi_surf():
    N = 1000
    xs    = np.linspace(-3.0,3.0,N)
    mmff.setVerbosity(-1)
    mmff.init( xyz_name=data_dir+"xyz/pyridine", sElementTypes=data_dir+"ElementTypes.dat", sAtomTypes=data_dir+"AtomTypes.dat", sBondTypes=data_dir+"BondTypes.dat", sAngleTypes=data_dir+"AngleTypes.dat", sDihedralTypes=data_dir+"DihedralTypes.dat")
    mmff.setSwitches(CheckInvariants=-1, PBC_nonBond=-1, PBC_evalAtom=-1, NonBonded=-1, MMFF=1, doBonds=1, Angles=1, PiSigma=1, PiPiI=1, bNonBondNeighs=1, bSurfAtoms=1, bGridFF=-1, bTricubic=-1, bConstrZ=-1, bConstrains=-1, bMoving=-1)
    Es,Fs = mmff.sample_movementOfAtom(xs, ia=0)
    Fnum = numDeriv(xs,Es)
    mmff.clear()
    assert np.all(np.abs(Fs[1:-1]-Fnum)<0.1), "!Check Forces::run_omp::eval_atom! run_omp force calculation failed: analytical forces do not match numerical derivatives"


# SetSwitches+run_omp --> non bonding
def test_run_omp_nonBond():
    N=1000
    xs    = np.linspace(-3.0,-2.1,N)
    mmff.setVerbosity(-1)
    mmff.init( xyz_name=data_dir+"xyz/pyridine", sElementTypes=data_dir+"ElementTypes.dat", sAtomTypes=data_dir+"AtomTypes.dat", sBondTypes=data_dir+"BondTypes.dat", sAngleTypes=data_dir+"AngleTypes.dat", sDihedralTypes=data_dir+"DihedralTypes.dat")
    mmff.setSwitches(CheckInvariants=-1, PBC_nonBond=-1, PBC_evalAtom=-1, NonBonded=1, MMFF=1, doBonds=-1, Angles=-1, PiSigma=-1, PiPiI=-1, bNonBondNeighs=-1, bSurfAtoms=-1, bGridFF=-1, bTricubic=-1, bConstrZ=-1, bConstrains=-1, bMoving=-1)
    Es,Fs = mmff.sample_movementOfAtom(xs, ia=0)
    Fnum = numDeriv(xs,Es)
    mmff.clear()
    assert np.all(np.abs(Fs[1:-1]-Fnum)<0.1), "!Check Forces::run_nonBond! run_nonBond force calculation failed: analytical forces do not match numerical derivatives"

def test_run_omp_nonBond_NonBondNeighs():
    N = 1000
    xs    = np.linspace(-3.0,3.0,N)
    mmff.setVerbosity(-1)
    mmff.init( xyz_name=data_dir+"xyz/pyridine", sElementTypes=data_dir+"ElementTypes.dat", sAtomTypes=data_dir+"AtomTypes.dat", sBondTypes=data_dir+"BondTypes.dat", sAngleTypes=data_dir+"AngleTypes.dat", sDihedralTypes=data_dir+"DihedralTypes.dat")
    mmff.setSwitches(CheckInvariants=-1, PBC_nonBond=-1, PBC_evalAtom=-1, NonBonded=1, MMFF=1, doBonds=-1, Angles=-1, PiSigma=-1, PiPiI=-1, bNonBondNeighs=1, bSurfAtoms=-1, bGridFF=-1, bTricubic=-1, bConstrZ=-1, bConstrains=-1, bMoving=-1)
    Es,Fs = mmff.sample_movementOfAtom(xs, ia=0)
    Fnum = numDeriv(xs,Es)
    mmff.clear()
    assert np.all(np.abs(Fs[1:-1]-Fnum)<0.1), "!Check Forces::run_omp_nonBond_NonBondNeighs! run_omp force calculation failed: analytical forces do not match numerical derivatives"

def test_run_omp_nonBond_PBC():
    N = 1000
    xs    = np.linspace(-3.0,-2.1,N)
    mmff.setVerbosity(-1)
    mmff.init( xyz_name=data_dir+"xyz/pyridine", sElementTypes=data_dir+"ElementTypes.dat", sAtomTypes=data_dir+"AtomTypes.dat", sBondTypes=data_dir+"BondTypes.dat", sAngleTypes=data_dir+"AngleTypes.dat", sDihedralTypes=data_dir+"DihedralTypes.dat")
    mmff.setSwitches(CheckInvariants=-1, PBC_nonBond=1, PBC_evalAtom=-1, NonBonded=1, MMFF=1, doBonds=-1, Angles=-1, PiSigma=-1, PiPiI=-1, bNonBondNeighs=-1, bSurfAtoms=-1, bGridFF=-1, bTricubic=-1, bConstrZ=-1, bConstrains=-1, bMoving=-1)
    Es,Fs = mmff.sample_movementOfAtom(xs, ia=0)
    Fnum = numDeriv(xs,Es)
    mmff.clear()
    assert np.all(np.abs(Fs[1:-1]-Fnum)<0.1), "!Check Forces::run_omp_nonBond_PBC! run_omp force calculation failed: analytical forces do not match numerical derivatives"

def test_run_omp_nonBond_PBC_NonBondNeighs():
    N = 1000
    xs    = np.linspace(-3.0,3.0,N)
    mmff.setVerbosity(-1)
    mmff.init( xyz_name=data_dir+"xyz/pyridine", sElementTypes=data_dir+"ElementTypes.dat", sAtomTypes=data_dir+"AtomTypes.dat", sBondTypes=data_dir+"BondTypes.dat", sAngleTypes=data_dir+"AngleTypes.dat", sDihedralTypes=data_dir+"DihedralTypes.dat")
    mmff.setSwitches(CheckInvariants=-1, PBC_nonBond=1, PBC_evalAtom=-1, NonBonded=1, MMFF=1, doBonds=-1, Angles=-1, PiSigma=-1, PiPiI=-1, bNonBondNeighs=1, bSurfAtoms=-1, bGridFF=-1, bTricubic=-1, bConstrZ=-1, bConstrains=-1, bMoving=-1)
    Es,Fs = mmff.sample_movementOfAtom(xs, ia=0)
    Fnum = numDeriv(xs,Es)
    mmff.clear()
    assert np.all(np.abs(Fs[1:-1]-Fnum)<0.1), "!Check Forces::run_omp_nonBond_PBC_NonBondNeighs! run_omp force calculation failed: analytical forces do not match numerical derivatives"



# SetSwitches+run_omp --> complete
def test_run_omp_complete():
    N = 1000
    xs    = np.linspace(2.3,3.0,N)
    mmff.setVerbosity(-1)
    mmff.init( xyz_name=data_dir+"xyz/dicarboxylic_acid", constr_name=data_dir+"dicarboxylic_acid.cons",  surf_name=data_dir+"xyz/NaCl_1x1_L2", sElementTypes=data_dir+"ElementTypes.dat", sAtomTypes=data_dir+"AtomTypes.dat", sBondTypes=data_dir+"BondTypes.dat", sAngleTypes=data_dir+"AngleTypes.dat", sDihedralTypes=data_dir+"DihedralTypes.dat")
    mmff.setSwitches(CheckInvariants=-1, PBC_nonBond=1, PBC_evalAtom=1, NonBonded=1, MMFF=1, doBonds=1, Angles=1, PiSigma=1, PiPiI=1, bNonBondNeighs=1, bSurfAtoms=1, bGridFF=-1, bTricubic=1, bConstrZ=-1, bConstrains=1, bMoving=-1)
    Es,Fs = mmff.sample_movementOfAtom(xs, ia=0)
    Fnum = numDeriv(xs,Es)
    mmff.clear()
    assert np.all(np.abs(Fs[1:-1]-Fnum)<0.1), "!Check Forces::run_omp_complete! run_omp force calculation failed: analytical forces do not match numerical derivatives"





# Constrains
def test_sample_DistConstr():
    N = 10
    xs    = np.linspace(0.0,3.0,N)
    Es,Fs = mmff.sample_DistConstr( xs, lmin=0.9, lmax=1.2, flim=0.5 )
    Fnum = numDeriv(xs,Es)
    assert np.all(np.abs(Fs[1:-1]-Fnum)<0.1), "!Check Constrains::DistConstr! sample_DistConstr force calculation failed: analytical forces do not match numerical derivatives"

def test_sample_SplineConstr():
    dx    = 1.5
    x0    = 0.5 
    Eps   = np.array( [1.0, 0.0,-1.0,-0.5,-0.2,-0.1] )
    xs    = np.linspace(0.0,6.0,100)
    Es,Fs = mmff.sample_SplineConstr( xs, Eps, x0=x0, dx=dx )
    Fnum = numDeriv(xs,Es)
    assert np.all(np.abs(Fs[1:-1]-Fnum)<0.01), "!Check Constrains::SplineConstr! sample_SplineConstr force calculation failed: analytical forces do not match numerical derivatives"




'''
# ----- lvecs
N = 100
xs    = np.linspace(1.0,3.0,N)
mmff.init( xyz_name=data_dir+"xyz/polymer-2_new",  surf_name=data_dir+"xyz/NaCl_1x1_L2", sElementTypes=data_dir+"ElementTypes.dat", sAtomTypes=data_dir+"AtomTypes.dat", sBondTypes=data_dir+"BondTypes.dat", sAngleTypes=data_dir+"AngleTypes.dat", sDihedralTypes=data_dir+"DihedralTypes.dat")
mmff.setSwitches(CheckInvariants=-1, PBC_nonBond=1, PBC_evalAtom=1, NonBonded=1, MMFF=1, doBonds=1, Angles=1, PiSigma=1, PiPiI=1, bNonBondNeighs=1, bSurfAtoms=-1, bGridFF=-1, bTricubic=-1, bConstrZ=-1, bConstrains=-1, bMoving=-1)
Es, Fs = mmff.sample_lvecs( xs )
Fnum = numDeriv(xs,Es)     #ratios = Fs[1:-1]/Fnum    ;print("ratios", ratios)
plt.figure(); plt.plot(xs, Es, label="E"); plt.plot(xs, Fs, label="F_ana");  plt.plot(xs[1:-1], Fnum, label="F_num"); plt.grid(); plt.legend()
#print(Fs)
plt.show()
'''

'''

# gridFF
xs    = np.linspace(-3.0,-2.7,1000)
mmff.init( xyz_name=data_dir+"xyz/pyridine", sElementTypes=data_dir+"ElementTypes.dat", sAtomTypes=data_dir+"AtomTypes.dat", sBondTypes=data_dir+"BondTypes.dat", sAngleTypes=data_dir+"AngleTypes.dat", sDihedralTypes=data_dir+"DihedralTypes.dat")
mmff.setVerbosity(3)
Es,Fs  = mmff.sample_addForce_Tricubic         ( xs)
Fnum = numDeriv(xs,Es)     #ratios = Fs[1:-1]/Fnum    ;print("ratios", ratios)
plt.figure(); plt.plot(xs, Es, label="E"); plt.plot(xs, Fs, label="F_ana");  plt.plot(xs[1:-1], Fnum, label="F_num"); plt.grid(); plt.legend()
plt.show()
'''

'''
# ----- Surf
xs    = np.linspace(0.0,3.0,100)
mmff.init( xyz_name="data/xyz/pyridine", surf_name="data/xyz/NaCl_1x1_L2" )
#Es,Fs = mmff.sample_addForce_Tricubic( xs)
#Es,Fs = mmff.sample_addForce(xs)
Fnum = numDeriv(xs,Es)     #ratios = Fs[1:-1]/Fnum    ;print("ratios", ratios)
plt.figure(); plt.plot(xs, Es, label="E"); plt.plot(xs, Fs, label="F_ana");  plt.plot(xs[1:-1], Fnum, label="F_num"); plt.grid(); plt.legend()
plt.show()
'''

'''
# ----- Angle
#mmff.sample_evalAngleCos( xs, lmin=1, lmax=1, kmin=1, kmax=1, flim=1e+300, Es=None, Fs=None)
xs    = np.linspace(-np.pi,np.pi,100)
#Es,Fs = mmff.sample_evalAngleCos( xs, ang0=np.pi*0.00 )
#print("Fs",Fs)
#Es,Fs = mmff.sample_evalAngleCos( xs, ang0=np.pi*0.25 )
#Es,Fs = mmff.sample_evalAngleCos( xs, ang0=np.pi*0.50 )
#Es,Fs = mmff.sample_evalAngleCos( xs, ang0=np.pi*0.75 )
#Es,Fs = mmff.sample_evalAngleCos( xs, ang0=np.pi*1.00 )

#Es,Fs = mmff.sample_evalAngleCosHalf( xs, ang0=np.pi*0.00 )
#Es,Fs = mmff.sample_evalAngleCosHalf( xs, ang0=np.pi*0.1 )
#Es,Fs = mmff.sample_evalAngleCosHalf( xs, ang0=np.pi*0.25 )
#Es,Fs = mmff.sample_evalAngleCosHalf( xs, ang0=np.pi*0.50 )
#Es,Fs = mmff.sample_evalAngleCosHalf ( xs, ang0=np.pi*0.75 )
#Es,Fs = mmff.sample_evalAngleCosHalf( xs, ang0=np.pi*0.9 )
#Es,Fs = mmff.sample_evalAngleCosHalf( xs, ang0=np.pi*1.00 )

Es,Fs = cos_half( xs)
Fnum = numDeriv(xs,Es)     #ratios = Fs[1:-1]/Fnum    ;print("ratios", ratios)
xs/=np.pi
plt.figure(); plt.plot(xs, Es, label="E"); plt.plot(xs, Fs, label="F_ana");  plt.plot(xs[1:-1], Fnum, label="F_num"); plt.grid(); plt.legend()
plt.show()




# ----- PiPiAlignment
xs    = np.linspace(-np.pi,np.pi,1000)
Es,Fs = mmff.sample_evalPiAling( xs, ang0=np.pi*0.00 )
Fnum = numDeriv(xs,Es)     #ratios = Fs[1:-1]/Fnum    ;print("ratios", ratios)
xs/=np.pi
plt.figure(); plt.plot(xs, Es, label="E"); plt.plot(xs, Fs, label="F_ana");  plt.plot(xs[1:-1], Fnum, label="F_num"); plt.grid(); plt.legend()
plt.show()
'''