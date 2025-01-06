import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    force=False  # This ensures our configuration takes precedence
)

LOGGER = logging.getLogger(__name__)

current_dir = os.path.dirname(__file__)
src_path = os.path.join(current_dir, "../")
sys.path.insert(0, os.path.abspath(src_path))

from pyBall import atomicUtils as au
from pyBall import MMFF as mmff

def numDeriv( xs, Es): 
    dx = xs[1]-xs[0]
    Fs = (Es[2:]-Es[:-2])/(2*dx)
    return -Fs

def cos_half( a ):
    E = (1-np.cos(a*0.5))*0.5
    F = np.sin(a*0.5)*0.25
    return E,F

#======== Body

'''
# ----- Angle
#mmff.sample_evalAngleCos( xs, lmin=1, lmax=1, kmin=1, kmax=1, flim=1e+300, Es=None, Fs=None)
xs    = np.linspace(-np.pi,np.pi,100)
Es,Fs = mmff.sample_evalAngleCos( xs, ang0=np.pi*0.00 )
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
Fnum = numDeriv(xs,Es)     #ratios = Fs[1:-1]/Fnum    ;print("ratios", ratios)
#Es,Fs = cos_half( xs + np.pi*0.5 )
xs/=np.pi
plt.figure(); plt.plot(xs, Es, label="E"); plt.plot(xs, Fs, label="F_ana");  plt.plot(xs[1:-1], Fnum, label="F_num"); plt.grid(); plt.legend()
plt.show()
'''


'''
# ----- PiPiAlignment
xs    = np.linspace(-np.pi,np.pi,1000)
Es,Fs = mmff.sample_evalPiAling( xs, ang0=np.pi*0.00 )
Fnum = numDeriv(xs,Es)     #ratios = Fs[1:-1]/Fnum    ;print("ratios", ratios)
xs/=np.pi
plt.figure(); plt.plot(xs, Es, label="E"); plt.plot(xs, Fs, label="F_ana");  plt.plot(xs[1:-1], Fnum, label="F_num"); plt.grid(); plt.legend()
plt.show()
'''  

# Bonding interaction
def test_sample_evalBond():
    xs    = np.linspace(0.1,3.0,100)
    mmff.setVerbosity(-1)
    Es,Fs = mmff.sample_evalBond( xs)
    Fnum = numDeriv(xs,Es)
    assert np.all(np.abs(Fs[1:-1]-Fnum)<0.1), "!Check Forces::evalBond! sample_evalBond force calculation failed: analytical forces do not match numerical derivatives"
 
def test_sample_evalAtom():
    xs    = np.linspace(-3.0,3.0,1000)
    mmff.setVerbosity(-1)
    mmff.init( xyz_name=current_dir+"/common_resources/xyz/polymer-2_new", sElementTypes=current_dir+"/common_resources/ElementTypes.dat", sAtomTypes=current_dir+"/common_resources/AtomTypes.dat", sBondTypes=current_dir+"/common_resources/BondTypes.dat", sAngleTypes=current_dir+"/common_resources/AngleTypes.dat", sDihedralTypes=current_dir+"/common_resources/DihedralTypes.dat")
    Es,Fs = mmff.sample_evalAtom( xs, ia=1)
    mmff.clear()
    Fnum = numDeriv(xs,Es)
    assert np.all(np.abs(Fs[1:-1]-Fnum)<1.0), "!Check Forces::evalAtom! sample_evalAtom force calculation failed: analytical forces do not match numerical derivatives"
xs    = np.linspace(-3.0,3.0,1000)
mmff.setVerbosity(-1)
mmff.init( xyz_name=current_dir+"/common_resources/xyz/polymer-2_new", sElementTypes=current_dir+"/common_resources/ElementTypes.dat", sAtomTypes=current_dir+"/common_resources/AtomTypes.dat", sBondTypes=current_dir+"/common_resources/BondTypes.dat", sAngleTypes=current_dir+"/common_resources/AngleTypes.dat", sDihedralTypes=current_dir+"/common_resources/DihedralTypes.dat")
mmff.setSwitches(CheckInvariants=-1, PBC_nonBond=-1, PBC_evalAtom=-1, NonBonded=-1, MMFF=1, doBonds=1, Angles=1, PiSigma=1, PiPiI=1, bNonBondNeighs=1, bSurfAtoms=1, bGridFF=-1, bTricubic=-1, bConstrZ=-1, bConstrains=-1)
Es,Fs = mmff.sample_evalAtom( xs, ia=0)
mmff.clear()
Fnum = numDeriv(xs,Es)
plt.figure(); plt.plot(xs, Es, label="E"); plt.plot(xs, Fs, label="F_ana");  plt.plot(xs[1:-1], Fnum, label="F_num"); plt.grid(); plt.legend()
plt.show()

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
    mmff.init( xyz_name=current_dir+"/data/xyz/pyridine", sElementTypes=current_dir+"/common_resources/ElementTypes.dat", sAtomTypes=current_dir+"/common_resources/AtomTypes.dat", sBondTypes=current_dir+"/common_resources/BondTypes.dat", sAngleTypes=current_dir+"/common_resources/AngleTypes.dat", sDihedralTypes=current_dir+"/common_resources/DihedralTypes.dat")
    Es,Fs = mmff.sample_evalLJQs_ng4_PBC_atom_omp( xs)
    mmff.clear()
    Fnum = numDeriv(xs,Es)
    assert np.all(np.abs(Fs[1:-1]-Fnum)<0.1), "!Check Forces::evalLJQs_ng4_PBC_atom_omp! sample_evalLJQs_ng4_PBC_atom_omp force calculation failed: analytical forces do not match numerical derivatives"

def test_sample_evalLJQs_ng4_atom_omp():
    xs    = np.linspace(-3.0,-2.7,1000)
    mmff.setVerbosity(-1)
    mmff.init( xyz_name=current_dir+"/data/xyz/pyridine", sElementTypes=current_dir+"/common_resources/ElementTypes.dat", sAtomTypes=current_dir+"/common_resources/AtomTypes.dat", sBondTypes=current_dir+"/common_resources/BondTypes.dat", sAngleTypes=current_dir+"/common_resources/AngleTypes.dat", sDihedralTypes=current_dir+"/common_resources/DihedralTypes.dat")
    Es,Fs = mmff.sample_evalLJQs_ng4_atom_omp      ( xs)
    mmff.clear()
    Fnum = numDeriv(xs,Es)
    assert np.all(np.abs(Fs[1:-1]-Fnum)<0.1), "!Check Forces::evalLJQs_ng4_atom_omp! sample_evalLJQs_ng4_atom_omp force calculation failed: analytical forces do not match numerical derivatives"

def test_sample_evalLJQs_PBC_atom_omp():
    xs    = np.linspace(-3.0,-2.7,1000)
    mmff.setVerbosity(-1)
    mmff.init( xyz_name=current_dir+"/data/xyz/pyridine", sElementTypes=current_dir+"/common_resources/ElementTypes.dat", sAtomTypes=current_dir+"/common_resources/AtomTypes.dat", sBondTypes=current_dir+"/common_resources/BondTypes.dat", sAngleTypes=current_dir+"/common_resources/AngleTypes.dat", sDihedralTypes=current_dir+"/common_resources/DihedralTypes.dat")
    Es,Fs = mmff.sample_evalLJQs_PBC_atom_omp      ( xs)
    mmff.clear()
    Fnum = numDeriv(xs,Es)
    assert np.all(np.abs(Fs[1:-1]-Fnum)<0.1), "!Check Forces::evalLJQs_PBC_atom_omp!"

def test_sample_evalLJQs_atom_omp():
    xs    = np.linspace(-3.0,-2.7,1000)
    mmff.setVerbosity(-1)
    mmff.init( xyz_name=current_dir+"/data/xyz/pyridine", sElementTypes=current_dir+"/common_resources/ElementTypes.dat", sAtomTypes=current_dir+"/common_resources/AtomTypes.dat", sBondTypes=current_dir+"/common_resources/BondTypes.dat", sAngleTypes=current_dir+"/common_resources/AngleTypes.dat", sDihedralTypes=current_dir+"/common_resources/DihedralTypes.dat")
    Es,Fs = mmff.sample_evalLJQs_atom_omp      ( xs)
    mmff.clear()
    Fnum = numDeriv(xs,Es)
    assert np.all(np.abs(Fs[1:-1]-Fnum)<0.1), "!Check Forces::evalLJQs_atom_omp! sample_evalLJQs_atom_omp force calculation failed: analytical forces do not match numerical derivatives"


'''
# ----- evalLJQs_ng4_PBC_atom_omp, evalLJQs_ng4_atom_omp,...
xs    = np.linspace(-3.0,-2.7,1000)
mmff.init( xyz_name="data/xyz/pyridine" )
mmff.setVerbosity(3)
Es,Fs  = mmff.sample_addForce_Tricubic         ( xs)
# Fnum = numDeriv(xs,Es)     #ratios = Fs[1:-1]/Fnum    ;print("ratios", ratios)
# plt.figure(); plt.plot(xs, Es, label="E"); plt.plot(xs, Fs, label="F_ana");  plt.plot(xs[1:-1], Fnum, label="F_num"); plt.grid(); plt.legend()
plt.show()
'''

# Interaction with surface
def test_sample_evalMorsePBC_sym():
    xs    = np.linspace(0.0,3.0,100)
    mmff.setVerbosity(-1)
    mmff.init( xyz_name=current_dir+"/data/xyz/pyridine", surf_name=current_dir+"/data/xyz/NaCl_1x1_L2", sElementTypes=current_dir+"/common_resources/ElementTypes.dat", sAtomTypes=current_dir+"/common_resources/AtomTypes.dat", sBondTypes=current_dir+"/common_resources/BondTypes.dat", sAngleTypes=current_dir+"/common_resources/AngleTypes.dat", sDihedralTypes=current_dir+"/common_resources/DihedralTypes.dat")
    mmff.setSwitches(bGridFF=-1)
    Es,Fs = mmff.sample_evalMorsePBC_sym( xs)
    Fnum = numDeriv(xs,Es)
    mmff.clear()
    assert np.all(np.abs(Fs[1:-1]-Fnum)<0.1), "!Check Forces::evalMorsePBC_sym! sample_evalMorsePBC_sym force calculation failed: analytical forces do not match numerical derivatives"

def test_sample_springbound():
    xs    = np.linspace(0.0,3.0,100)
    mmff.setVerbosity(-1)
    mmff.init( xyz_name=current_dir+"/data/xyz/pyridine", surf_name=current_dir+"/data/xyz/NaCl_1x1_L2", sElementTypes=current_dir+"/common_resources/ElementTypes.dat", sAtomTypes=current_dir+"/common_resources/AtomTypes.dat", sBondTypes=current_dir+"/common_resources/BondTypes.dat", sAngleTypes=current_dir+"/common_resources/AngleTypes.dat", sDihedralTypes=current_dir+"/common_resources/DihedralTypes.dat")
    Es,Fs = mmff.sample_springbound( xs)
    Fnum = numDeriv(xs,Es)
    mmff.clear()
    assert np.all(np.abs(Fs[1:-1]-Fnum)<0.1), "!Check Forces::springbound! sample_springbound force calculation failed: analytical forces do not match numerical derivatives"
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

# Constrains
def test_sample_applyConstr():
    xs    = np.linspace(0.0,3.0,100)
    mmff.setVerbosity(-1)
    mmff.init( xyz_name=current_dir+"/data/xyz/dicarboxylic_acid", constr_name=current_dir+"/data/dicarboxylic_acid.cons", sElementTypes=current_dir+"/common_resources/ElementTypes.dat", sAtomTypes=current_dir+"/common_resources/AtomTypes.dat", sBondTypes=current_dir+"/common_resources/BondTypes.dat", sAngleTypes=current_dir+"/common_resources/AngleTypes.dat", sDihedralTypes=current_dir+"/common_resources/DihedralTypes.dat")
    Es,Fs = mmff.sample_applyConstr( xs )
    Fnum = numDeriv(xs,Es)
    mmff.clear()
    assert np.all(np.abs(Fs[1:-1]-Fnum)<0.1), "!Check Forces::applyConstr! sample_applyConstr force calculation failed: analytical forces do not match numerical derivatives"



# ----- SetSwitches+run_omp
def test_run_omp_bonds_angles_pisigma_pipi_surf():
    N = 1000
    xs    = np.linspace(-3.0,3.0,N)
    mmff.setVerbosity(-1)
    mmff.init( xyz_name=current_dir+"/data/xyz/pyridine", surf_name=current_dir+"/data/xyz/NaCl_1x1_L2", sElementTypes=current_dir+"/common_resources/ElementTypes.dat", sAtomTypes=current_dir+"/common_resources/AtomTypes.dat", sBondTypes=current_dir+"/common_resources/BondTypes.dat", sAngleTypes=current_dir+"/common_resources/AngleTypes.dat", sDihedralTypes=current_dir+"/common_resources/DihedralTypes.dat")
    mmff.setSwitches(CheckInvariants=-1, PBC_nonBond=-1, PBC_evalAtom=-1, NonBonded=-1, MMFF=1, doBonds=1, Angles=1, PiSigma=1, PiPiI=1, bNonBondNeighs=1, bSurfAtoms=1, bGridFF=-1, bTricubic=-1, bConstrZ=-1, bConstrains=-1)
    Es,Fs = mmff.sample_movementOfAtom(xs, ia=0)
    Fnum = numDeriv(xs,Es)
    mmff.clear()
    assert np.all(np.abs(Fs[1:-1]-Fnum)<0.1)

'''
N = 1000
xs    = np.linspace(-3.0,3.0,N)
mmff.init( xyz_name="data/xyz/polymer-2_new",  surf_name="data/xyz/NaCl_1x1_L2" )
#mmff.init( xyz_name="data/O2")# )
mmff.setSwitches(CheckInvariants=-1, PBC_nonBond=-1, PBC_evalAtom=-1, NonBonded=1, MMFF=1, doBonds=1, Angles=1, PiSigma=1, PiPiI=1, bNonBondNeighs=1, bSurfAtoms=1, bGridFF=-1, bTricubic=-1, bConstrZ=-1, bConstrains=-1)
Es,Fs = mmff.sample_movementOfAtom(xs, ia=0)
Fnum = numDeriv(xs,Es)     #ratios = Fs[1:-1]/Fnum    ;print("ratios", ratios)
plt.figure(); plt.plot(xs, Es, label="E"); plt.plot(xs, Fs, label="F_ana");  plt.plot(xs[1:-1], Fnum, label="F_num"); plt.grid(); plt.legend()
plt.show()


# ----- lvecs
N = 10
xs    = np.linspace(1.0,3.0,N)
mmff.init( xyz_name="common_resources/xyz/polymer-2_new",  surf_name="common_resources/xyz/NaCl_1x1_L2" )
mmff.setSwitches(CheckInvariants=-1, PBC_nonBond=1, PBC_evalAtom=1, NonBonded=1, MMFF=1, doBonds=1, Angles=1, PiSigma=1, PiPiI=1, bNonBondNeighs=1, bSurfAtoms=-1, bGridFF=-1, bTricubic=-1, bConstrZ=-1, bConstrains=-1)
Es, Fs = mmff.sample_lvecs( xs )
Fnum = numDeriv(xs,Es)     #ratios = Fs[1:-1]/Fnum    ;print("ratios", ratios)
plt.figure(); plt.plot(xs, Es, label="E"); plt.plot(xs, Fs, label="F_ana");  plt.plot(xs[1:-1], Fnum, label="F_num"); plt.grid(); plt.legend()
#print(Fs)
plt.show()
'''

'''constr_name="hexan-dicarboxylic.cons",
# ----- DistConstr
#mmff.sample_DistConstr( xs, lmin=1, lmax=1, kmin=1, kmax=1, flim=1e+300, Es=None, Fs=None)
xs    = np.linspace(0.0,3.0,100)
Es,Fs = mmff.sample_DistConstr( xs, lmin=0.9, lmax=1.2, flim=0.5 )  # ;print("Fs",Fs)
plt.figure(); plt.plot(xs, Es,'.-', label="E"); plt.plot(xs, Fs, label="F_ana");  plt.plot(xs[1:-1], numDeriv(xs,Es), label="F_num"); plt.grid(); plt.legend()
'''

'''
# ----- SplineConstr
dx    = 1.5
x0    = 0.5 
Eps   = np.array( [1.0, 0.0,-1.0,-0.5,-0.2,-0.1] )
xp    = (np.array(range(len(Eps)))-1)*dx + x0
xs    = np.linspace(0.0,6.0,100)
Es,Fs = mmff.sample_SplineConstr( xs, Eps, x0=x0, dx=dx )  # ;print("Fs",Fs)
plt.figure(); 
plt.plot(xp, Eps, 'o-k', lw=0.2, label="Eps"); 
plt.plot(xs, Es,'.-',    label="E"); 
plt.plot(xs, Fs,'-',     label="F_ana");  
plt.plot(xs[1:-1],numDeriv(xs,Es),':', label="F_num"); 
plt.grid(); plt.legend()

plt.show()'''