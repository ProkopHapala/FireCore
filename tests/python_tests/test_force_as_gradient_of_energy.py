import sys
#import numpy as np
import os
#import matplotlib.pyplot as plt

sys.path.append("/../")
#from pyBall import MMFF as mmff


def numDeriv( xs, Es): 
    dx = xs[1]-xs[0]
    Fs = (Es[2:]-Es[:-2])/(2*dx)
    return -Fs

'''
# ----- Distance
#mmff.sample_evalDist( xs, lmin=1, lmax=1, kmin=1, kmax=1, flim=1e+300, Es=None, Fs=None)
xs    = np.linspace(0,10,100)
#Es,Fs = mmff.sample_evalDist( xs, lmin=1, lmax=1, kmin=1, kmax=1, flim=1e+300, Es=None, Fs=None)
#Es,Fs = mmff.sample_evalDist( xs, lmin=1, lmax=1, kmin=1, kmax=2, flim=1e+300, Es=None, Fs=None)

# ----- Angle
#mmff.sample_evalAngleCos( xs, lmin=1, lmax=1, kmin=1, kmax=1, flim=1e+300, Es=None, Fs=None)
xs    = np.linspace(-np.pi,np.pi,100)
#Es,Fs = mmff.sample_evalAngleCos( xs, ang0=np.pi*0.00 )
#print("Fs",Fs)
#Es,Fs = mmff.sample_evalAngleCos( xs, ang0=np.pi*0.25 )
#Es,Fs = mmff.sample_evalAngleCos( xs, ang0=np.pi*0.50 )
#Es,Fs = mmff.sample_evalAngleCos( xs, ang0=np.pi*0.75 )
#Es,Fs = mmff.sample_evalAngleCos( xs, ang0=np.pi*1.00 )

Es,Fs = mmff.sample_evalAngleCosHalf( xs, ang0=np.pi*0.00 )
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

# ----- PiPiAlignment
xs    = np.linspace(-np.pi,np.pi,1000)
Es,Fs = mmff.sample_evalPiAling( xs, ang0=np.pi*0.00 )
Fnum = numDeriv(xs,Es)     #ratios = Fs[1:-1]/Fnum    ;print("ratios", ratios)
xs/=np.pi
plt.figure(); plt.plot(xs, Es, label="E"); plt.plot(xs, Fs, label="F_ana");  plt.plot(xs[1:-1], Fnum, label="F_num"); plt.grid(); plt.legend()
plt.show()



# ----- Bond
xs    = np.linspace(0.1,3.0,100)
Es,Fs = mmff.sample_evalBond( xs)
Fnum = numDeriv(xs,Es)     #ratios = Fs[1:-1]/Fnum    ;print("ratios", ratios)
plt.figure(); plt.plot(xs, Es, label="E"); plt.plot(xs, Fs, label="F_ana");  plt.plot(xs[1:-1], Fnum, label="F_num"); plt.grid(); plt.legend()
plt.show()

'''
'''
# ----- evalAtom
xs    = np.linspace(-3.0,3.0,1000)
mmff.init( xyz_name="data/polymer-2_new" )
Es,Fs = mmff.sample_evalAtom( xs, ia=1)
Fnum = numDeriv(xs,Es)     #ratios = Fs[1:-1]/Fnum    ;print("ratios", ratios)
plt.figure(); plt.plot(xs, Es, label="E"); plt.plot(xs, Fs, label="F_ana");  plt.plot(xs[1:-1], Fnum, label="F_num"); plt.grid(); plt.legend()
plt.show()

# ----- getLJQH
xs    = np.linspace(0.1,3.0,1000)
Es,Fs = mmff.sample_getLJQH( xs)
Fnum = numDeriv(xs,Es)     #ratios = Fs[1:-1]/Fnum    ;print("ratios", ratios)
plt.figure(); plt.plot(xs, Es, label="E"); plt.plot(xs, Fs, label="F_ana");  plt.plot(xs[1:-1], Fnum, label="F_num"); plt.grid(); plt.legend()
plt.show()

# ----- evalLJQs_ng4_PBC_atom_omp, evalLJQs_ng4_atom_omp,...
xs    = np.linspace(-3.0,3.0,100)
mmff.init( xyz_name="data/pyridine" )
#Es,Fs = mmff.sample_evalLJQs_ng4_PBC_atom_omp  ( xs)
#Es,Fs = mmff.sample_evalLJQs_ng4_atom_omp      ( xs)
Es,Fs  = mmff.sample_evalLJQs_PBC_atom_omp     ( xs)
#Es,Fs  = mmff.sample_evalLJQs_atom_omp         ( xs)
#Es,Fs  = mmff.sample_addForce_Tricubic         ( xs)
Fnum = numDeriv(xs,Es)     #ratios = Fs[1:-1]/Fnum    ;print("ratios", ratios)
plt.figure(); plt.plot(xs, Es, label="E"); plt.plot(xs, Fs, label="F_ana");  plt.plot(xs[1:-1], Fnum, label="F_num"); plt.grid(); plt.legend()
plt.show()

# ----- Surf
xs    = np.linspace(0.0,3.0,100)
mmff.init( xyz_name="data/pyridine", surf_name="data/NaCl_1x1_L2" )
#Es,Fs = mmff.sample_addForce_Tricubic( xs)
#Es,Fs = mmff.sample_addForce(xs)
#Es,Fs = mmff.sample_evalMorsePBC_sym( xs)
#Es,Fs = mmff.sample_springbound( xs)
Fnum = numDeriv(xs,Es)     #ratios = Fs[1:-1]/Fnum    ;print("ratios", ratios)
plt.figure(); plt.plot(xs, Es, label="E"); plt.plot(xs, Fs, label="F_ana");  plt.plot(xs[1:-1], Fnum, label="F_num"); plt.grid(); plt.legend()
plt.show()

# ----- Constrains
xs    = np.linspace(0.0,3.0,100)
mmff.init( xyz_name="data/nHexadecan_dicarboxylic", constr_name="hexan-dicarboxylic.cons" )
Es,Fs = mmff.sample_applyConstr( xs )
Fnum = numDeriv(xs,Es)     #ratios = Fs[1:-1]/Fnum    ;print("ratios", ratios)
plt.figure(); plt.plot(xs, Es, label="E"); plt.plot(xs, Fs, label="F_ana");  plt.plot(xs[1:-1], Fnum, label="F_num"); plt.grid(); plt.legend()
plt.show()

# ----- SetSwitches+run_omp
N = 1000
xs    = np.linspace(-3.0,3.0,N)
mmff.init( xyz_name="data/polymer-2_new",  surf_name="data/NaCl_1x1_L2" )
#mmff.init( xyz_name="data/O2")# )
mmff.setSwitches(CheckInvariants=-1, PBC_nonBond=1, PBC_evalAtom=1, NonBonded=1, MMFF=1, doBonds=1, Angles=1, PiSigma=1, PiPiI=1, bNonBondNeighs=1, bSurfAtoms=1, bGridFF=-1, bTricubic=1, bConstrZ=-1, bConstrains=-1)
Es,Fs = mmff.sample_movementOfAtom(xs, ia=1)
Fnum = numDeriv(xs,Es)     #ratios = Fs[1:-1]/Fnum    ;print("ratios", ratios)
plt.figure(); plt.plot(xs, Es, label="E"); plt.plot(xs, Fs, label="F_ana");  plt.plot(xs[1:-1], Fnum, label="F_num"); plt.grid(); plt.legend()
plt.show()

# ----- lvecs
N = 10
xs    = np.linspace(1.0,3.0,N)
mmff.init( xyz_name="data/polymer-2_new",  surf_name="data/NaCl_1x1_L2" )
mmff.setSwitches(CheckInvariants=-1, PBC_nonBond=1, PBC_evalAtom=1, NonBonded=1, MMFF=1, doBonds=1, Angles=1, PiSigma=1, PiPiI=1, bNonBondNeighs=1, bSurfAtoms=-1, bGridFF=-1, bTricubic=-1, bConstrZ=-1, bConstrains=-1)
Es, Fs = mmff.sample_lvecs( xs )
Fnum = numDeriv(xs,Es)     #ratios = Fs[1:-1]/Fnum    ;print("ratios", ratios)
plt.figure(); plt.plot(xs, Es, label="E"); plt.plot(xs, Fs, label="F_ana");  plt.plot(xs[1:-1], Fnum, label="F_num"); plt.grid(); plt.legend()
#print(Fs)
plt.show()
'''