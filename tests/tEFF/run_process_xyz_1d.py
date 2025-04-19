import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

sys.path.append("../../")
from pyBall import eFF as eff
from run_process_xyz import extract_blocks, extract_nae


print("#=========== RUN /home/prokophapala/git/FireCore/tests/tEFF/run_process_xyz_1d.py ")
fname = "export/scan_data/distscan_H2O.xyz"
#fname = "export/scan_data/distscan_H2.xyz"

params, nrec = extract_blocks(fname)

eff.setVerbosity(1,0)
eff.initOpt( dt=0.005, damping=0.005, f_limit=1000.0)

bCoreElectrons = False
eff.setSwitches( coreCoul=1 )

atomParams = np.array([
#  Q   sQ   sP   cP
[ 0.,  1.0, 1.00, 0.0 ], # 0
[ 1.,  0.0, 0.00, 0.0 ], # 1 H
[ 0.,  1.0, 1.00, 1.0 ], # 2 He
[ 1.,  0.0, 0.10, 1.0 ], # 3 Li
[ 2.,  0.0, 0.10, 1.0 ], # 4 Be
[ 3.,  0.0, 0.10, 1.0 ], # 5 B
[ 4.,  0.0, 0.10, 1.0 ], # 6 C
[ 5.,  0.0, 0.10, 1.0 ], # 7 N
[ 6.,  0.0, 0.15, 1.0 ], # 8 O
[ 7.,  0.0, 0.10, 1.0 ], # 9 F
], dtype=np.float64)
eff.setAtomParams( atomParams )


with open("processXYZ.xyz", "w") as f: f.write("")
eff.preAllocateXYZ(fname, Rfac=-1.35, bCoreElectrons=bCoreElectrons)
eff.getBuffs()

#na,ne = extract_nae(fname)
outEs = np.zeros((nrec,5))
out_apos = np.zeros((nrec,eff.na,3))  ;print( "out_apos.shape ", out_apos.shape )
out_epos = np.zeros((nrec,eff.ne,4))  ;print( "out_epos.shape ", out_epos.shape )

#eff.processXYZ("export/scan_data/distscan_H2O.xyz", bOutXYZ=False, outEs=outEs, bCoreElectrons=True, bChangeCore=True, nstepMax=0)
eff.processXYZ( fname, bOutXYZ=True, outEs=outEs, apos=out_apos, epos=out_epos, bCoreElectrons=bCoreElectrons, bChangeCore=False, bChangeEsize=True, nstepMax=100000, dt=0.005, Fconv=1e-4, ialg=2 );

bSubstractLast=True
#bSubstractLast=False

plt.figure(figsize=(5,10))
plt.subplot(2,1,1)
if bSubstractLast:
    params['Etot'] -= params['Etot'][-1]
    outEs[:,0]     -= outEs[-1,0]
plt.plot(params['dist'], params['Etot'], '.-r', label='DFT (ref)')
plt.plot(params['dist'], outEs[:,0],     '.-k', label='eFF')
plt.xlabel('Distance (Å)')
plt.ylabel('Total Energy (eV)')
plt.title('1D Energy Scan H2O vs Distance')
plt.legend()
plt.grid()
plt.subplot(2,1,2)
plt.plot(params['dist'], out_epos[:,0,0], '.-r', label='e1 x')
plt.plot(params['dist'], out_epos[:,1,0], '.-b', label='e2 x')
plt.plot(params['dist'], out_epos[:,0,3], '.-m', label='e1 s')
plt.plot(params['dist'], out_epos[:,1,3], '.-c', label='e2 s')
plt.xlabel('Distance (Å)')
plt.ylabel('Position (Å)')
plt.title('1D Energy Scan H2O vs Distance')
plt.legend()
plt.grid()
plt.savefig('scan1d_eFF.png')

print("#=========== DONE /home/prokophapala/git/FireCore/tests/tEFF/run_process_xyz_1d.py ")
plt.show()