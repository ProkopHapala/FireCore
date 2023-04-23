import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import plotUtils   as plu
from pyBall import MMFFsp3     as mmff



# ======= Setup


'''
Molecules of interest:

H2O
HCOOH
CH2O
HCN
HCCH
CO
NH3


CH3NO2


HF,HCl,HBr

'''



molecules = [
"HCCH",    
"HCN",
#"NH3",
#"H2O",
"CH2O",
#"HCOOH",
"CHONH2",

#"pyridine",
#"pyrrole",
#"furan",
#"oxazine",
#"dihydropyrazine",
#"quinone", 
#"maleic_anhydride",
]






#======== Body

#mmff.setVerbosity( verbosity=1, idebug=0 )
mmff.setVerbosity( verbosity=0, idebug=0 )
mmff.initParams()

#------ Job1: Orient
# mmff.buildMolecule_xyz( xyz_name="data/HCOOH", bEpairs=True )
# #mmff.orient( [1,3], [0,4] )
# mmff.orient( [3,4], [0,4] )

#f1 = mmff.buildMolecule_xyz( xyz_name="data/pyridine", bEpairs=True )
#f1 = mmff.buildMolecule_xyz( xyz_name="data/OCH2", bEpairs=True )
f1 = mmff.buildMolecule_xyz( xyz_name="data/PTCDA", bEpairs=True )
#f1 = mmff.buildMolecule_xyz( xyz_name="data/HCOOH", bEpairs=True )
#f1 = mmff.buildMolecule_xyz( xyz_name="data/H2O", bEpairs=True )
rot = mmff.findMainAxes( ifrag=f1 );  # print(" rot ", rot)
sym = mmff.findSymmetry( ifrag=f1 ); print("sym=",sym)

mmff.saveXYZ( "auto_rot.xyz","", 0 )


#mmff.scanAllHBonds( "data/HCOOH", "data/H2O" )

'''
#------ Job2: HBond-scan
f1 = mmff.buildMolecule_xyz( xyz_name="data/HCOOH", bEpairs=True )
f2 = mmff.buildMolecule_xyz( xyz_name="data/H2O",   bEpairs=True )
b1,_,_ = mmff.getFrament(f1);    print( "bonds_1", b1[2] )
b2,_,_ = mmff.getFrament(f2);    print( "bonds_2", b2[2] )
#mmff.buildMolecule_xyz( xyz_name="data/HCOOH" )
#mmff.saveXYZ( "builder_1.xyz","", 0 )
#mmff.printTypes();
#mmff.printAtomConfs()
# mmff.printAtomTypes()
#mmff.printBonds()
donors_1    = mmff.selectBondsBetweenTypes( b1[2,0], b1[2,1], 8, 1,   True, True ); print( "donors_1    ", donors_1    )
acceptors_1 = mmff.selectBondsBetweenTypes( b1[2,0], b1[2,1], 8, 200, True, True ); print( "acceptors_1 ", acceptors_1 )
donors_2    = mmff.selectBondsBetweenTypes( b2[2,0], b2[2,1], 8, 1,   True, True ); print( "donors_2    ", donors_2    )
acceptors_2 = mmff.selectBondsBetweenTypes( b2[2,0], b2[2,1], 8, 200, True, True ); print( "acceptors_2 ", acceptors_2 )
#mmff.scanHBond( donors_1[0], acceptors_2[0], l0=2.0, ups=[ (0.,0.,1.), (1.,0.,0.) ] )
mmff.scanHBond( acceptors_1[0], donors_2[0], l0=2.0, isDonor=(False,True), ups=[ (0.,0.,1.), (1.,0.,0.) ] )
#mmff.scanHBond( [3,4], [13,10],                     ups=[ (0.,0.,1.), (1.,0.,0.) ] )
'''


# mmff.makeFFs()
# mmff.getBuffs()   
# mmff.eval()
# nmaxiter= 10000
# cos,f,v,dt,damp = mmff.setOptLog( nmaxiter )
# print( "py DEBUG 6" )
# mmff.setTrjName("relax.xyz",1)
# print( "py DEBUG 7" )
# nsteps = mmff.run(nmaxiter, ialg=3 )  # run with FIRE_smooth     
# #nsteps = mmff.run(nmaxiter, ialg=2 )  # run with FIRE 
# print( "py DEBUG 8" )
# print( "Relaxed in ", nsteps )   


# print("f\n",f)
# print("v\n",v)
# print("damp\n",damp)
# print("cos\n",cos)


# plt.figure()
# plt.subplot(2,1,1)
# plt.plot( cos [:nsteps-1], label="cos(v,f)" );
# plt.legend(); plt.grid(); 
# plt.subplot(2,1,2)
# plt.plot( f   [:nsteps-1], label="|f|"  );
# plt.plot( v   [:nsteps-1], label="|v|"  );
# plt.plot( dt  [:nsteps-1], label="dt"   );
# plt.plot( damp[:nsteps-1], label="damp" );
# plt.legend(); plt.grid(); plt.yscale('log');  plt.ylim( 1e-8, 1e+4 );
# plt.show()
# print(f,v,damp,cos,dt)
#print( "FORCES:\n mmff.fapos:\n ", mmff.fapos )
#mmff.plot(bForce=True, Fscale=10.0 )
#plt.show()
#exit(0)

