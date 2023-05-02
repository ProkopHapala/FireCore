import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import plotUtils   as plu
from pyBall import MMFFsp3     as mmff


'''

pointer axes violation - try to use Mudflap
https://gcc.gnu.org/wiki/Mudflap_Pointer_Debugging ?

see: https://stackoverflow.com/questions/19534051/glibc-detect-smallbin-linked-list-corrupted

'''


#======== Body

#mmff.setVerbosity( verbosity=1, idebug=0 )
mmff.setVerbosity( verbosity=0, idebug=0 )
#mmff.init( xyz_name="data/pyridine"  ) 

#------ Short Initialization
#mmff.init( xyz_name="data/HCOOH"  ) 


'''
#------ Long Initialization
mmff.initParams()
print("=========== HCOOH ")
mmff.buildMolecule_xyz( xyz_name="data/HCOOH"  )
mmff.makeFFs()
'''

#------ Long Initialization
mmff.initParams()
print( "py DEBUG 1" )
#mmff.buildMolecule_xyz( xyz_name="data/HCOOH"  )
#mmff.buildMolecule_xyz( xyz_name="data/C2H4"  );    mmff.setSwitches( NonBonded=-1 )
#mmff.buildMolecule_xyz( xyz_name="data/CHOCHO"  );    mmff.setSwitches( NonBonded=-1 )
mmff.buildMolecule_xyz( xyz_name="data/CHONH2"  );    mmff.setSwitches( NonBonded=-1 )


#mmff.buildMolecule_xyz( xyz_name="data/propandiol"  )
print( "py DEBUG 2" )
mmff.makeFFs()
print( "py DEBUG 3" )
mmff.getBuffs()   
print( "py DEBUG 4" )
mmff.eval()
print( "py DEBUG 5" )

print( "mmff.ndoms \n", mmff.ndims )
print( "mmff.Es    \n", mmff.Es )
print( "mmff.DOFs  \n", mmff.DOFs )
print( "mmff.fDOFs \n", mmff.fDOFs )
print( "mmff.apos  \n", mmff.apos )
print( "mmff.pipos \n", mmff.pipos )
print( "mmff.fapos \n", mmff.fapos )
print( "mmff.fpipos\n", mmff.fpipos )



'''
mmff.getBuffs() 
labels = [ mmff.getType(i) for i in range(mmff.natoms) ]   ;print( labels ); print( mmff.apos )
plu.plotAtoms( mmff.apos, labels=labels, axes=(0,2) )
plt.show()

print("=========== pyridine ")
mmff.clear()
#mmff.buildMolecule_xyz( xyz_name="data/HCOOH"   )
mmff.buildMolecule_xyz( xyz_name="data/pyridine"  )
mmff.makeFFs()
print("=========== HCOOH ")
mmff.clear()
#mmff.buildMolecule_xyz( xyz_name="data/HCOOH"   )
mmff.buildMolecule_xyz( xyz_name="data/pyridine"  )
mmff.makeFFs()
'''



nmaxiter= 10000

cos,f,v,dt,damp = mmff.setOptLog( nmaxiter )
print( "py DEBUG 6" )
mmff.setTrjName("relax.xyz",1)
print( "py DEBUG 7" )
nsteps = mmff.run(nmaxiter, ialg=3 )  # run with FIRE_smooth     
#nsteps = mmff.run(nmaxiter, ialg=2 )  # run with FIRE 
print( "py DEBUG 8" )
print( "Relaxed in ", nsteps )   


print("f\n",f)
print("v\n",v)
print("damp\n",damp)
print("cos\n",cos)


plt.figure()
plt.subplot(2,1,1)
plt.plot( cos [:nsteps-1], label="cos(v,f)" );
plt.legend(); plt.grid(); 
plt.subplot(2,1,2)
plt.plot( f   [:nsteps-1], label="|f|"  );
plt.plot( v   [:nsteps-1], label="|v|"  );
plt.plot( dt  [:nsteps-1], label="dt"   );
plt.plot( damp[:nsteps-1], label="damp" );
plt.legend(); plt.grid(); plt.yscale('log');  plt.ylim( 1e-8, 1e+4 );
plt.show()
print(f,v,damp,cos,dt)



#print( "FORCES:\n mmff.fapos:\n ", mmff.fapos )
#mmff.plot(bForce=True, Fscale=10.0 )
#plt.show()
#exit(0)

