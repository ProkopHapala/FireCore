import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

sys.path.append("../../")
from pyBall import FitREQ as fit
from pyBall import atomicUtils as au

# ============== Setup

Hartree2eV   = 27.2114079527
kcal2eV      =  0.0433634

# --------- Load Reference energy
# #fname = "scan_H2O_b3lyp_cc-pvdz.xyz"
# r0 = 1.91
# Eref,xs=fit.EnergyFromXYZ(fname)    # load reference energy and position from xyz file
# xs+=r0
# Eref*=Hartree2eV


dat   = np.genfromtxt( "angScan_HBond_OCH2_vs_H2O_Emap.dat", comments='@' )   ;print( "dat[1,:] ", dat[1,:] )
Eref  =  dat[:,6]
rs    =  dat[:,2]
angs  =  dat[:,4]
Eref/=0.0433634

# ------ Set strength of hydrogen bond correction
#Hfactor = 0.99
#Hfactor = 0.9
Hfactor = 0.95
#Hfactor = 0.99
#Hfactor = 1.0
1.908   

# ------ Set atom types
types = [
    ([     1.487 ,    np.sqrt(0.0006808),     +0.2,    +np.sqrt(0.0006808)*Hfactor   ], [0,0,0,1] ),    # 0 H H2O
    ([     1.487 ,    np.sqrt(0.0006808),     +0.18,   +np.sqrt(0.0006808)*Hfactor*0 ], [0,0,0,0] ),    # 1 H CH2=O
    ([     1.661 ,    np.sqrt(0.0091063),     -0.4 ,   -np.sqrt(0.0091063)*Hfactor*0 ], [0,0,0,0] ),    # 2 O H2O
    ([     1.661 ,    np.sqrt(0.0091063),     -0.36,   -np.sqrt(0.0091063)*Hfactor   ], [0,0,0,1] ),    # 3 O CH2=O 
    ([     1.908 ,    np.sqrt(0.0037292),     -0.0 ,   -np.sqrt(0.0037292)*Hfactor*0 ], [0,0,0,0] ),    # 4 C CH2=O
]
typeMask = np.array( [ t[1] for t in types ], dtype=np.int32  )
typREQs  = np.array( [ t[0] for t in types ], dtype=np.double )
fit.init_types( typeMask, typREQs, bCopy=True ) 
#fit.loadXYZ( fname, [0,1,2], [3,4,5] )
#fit.loadXYZ( fname+"_scan.xyz", [3,4,5], [0,1,2], types0=[1,0,0], testtypes=[1,0,0]  )     # load reference geometry
fit.loadXYZ( "angScan_HBond_OCH2_vs_H2O.xyz", [0,1,2,3], [4,5,6], types0=[3,4,1,1], testtypes=[2,0,0]  )     # load reference geometry


fit.getBuffs()
print( "typToREQ" , fit.typToREQ )
print( "typeREQs" , fit.typeREQs )
print( "types1"   , fit.types1   )
#print( "types2"   , fit.types2   )
#print( "types3"   , fit.types3   )

# ------ obtain energy profile from classical model (fit)
Es     = fit.getEs( imodel=2, bRigid=False)     #;print( "Es_noH", Es     )

# ------ obtain energy progile with HBond correction set to zero
typREQs[0,3] = 0.0              # set HBond correction to zero
fit.setType(0, typREQs[0,:] )   # set atom type 0
Es_noH = fit.getEs( imodel=1, bRigid=False)     #;print( "Es_noH", Es_noH )


#exit()


n2=6

rs    =  np.reshape( rs    ,(-1,n2)).transpose()
angs  =  np.reshape( angs  ,(-1,n2)).transpose()
Eref   = np.reshape( Eref  ,(-1,n2)).transpose()
Es     = np.reshape( Es    ,(-1,n2)).transpose()
Es_noH = np.reshape( Es_noH,(-1,n2)).transpose()

iang = 0

plt.plot( rs[iang], Eref  [iang], '-k', label="E_ref", lw=3 )
plt.plot( rs[iang], Es    [iang], '-g', label="E",     lw=2 )
plt.plot( rs[iang], Es_noH[iang], '-r', label="E_noH", lw=2 )
plt.axhline(0.0, c='k', ls='--')
plt.grid()
plt.legend()
plt.ylim( -0.3, 0.3 )


'''
plt.figure(figsize=(15,5))
plt.subplot( 3,2,1 ); plt.imshow( Eref   )
plt.subplot( 3,2,3 ); plt.imshow( Es     )
plt.subplot( 3,2,4 ); plt.imshow( Es     - Eref, cmap='bwr' )
plt.subplot( 3,2,5 ); plt.imshow( Es_noH )
plt.subplot( 3,2,6 ); plt.imshow( Es_noH - Eref, cmap='bwr' )
'''

plt.show()


'''
# -------- Plot graphs
Emin = Eref.min()
plt.plot(xs,Eref  , '-',  label="E_ref", lw=3, c='grey' )
plt.plot(xs,Es_noH, '-r', label="E_noH" )
plt.plot(xs,Es    , '-g', label="E_fit" )
plt.ylim( Emin*1.2, -Emin )

plt.xlabel("r(O-H) [A]"); plt.ylabel("E[eV]");
plt.legend(); plt.grid()
plt.savefig("Hbond_correction.png", bbox_inches='tight')
plt.show()
'''



'''''
examples of AtomType parameters

*      *     H       E    1  0   0   0   0.0  0.0   0.0006808  0.0     0.0
#1     2     3       4    5  6   7   8    9       10      11       12      13       15     16      17  18  19  20
#name parent element epair nv ne npi sym  Ruff   RvdW    EvdW     Qbase    Hb      Ass    Asp     Kss Ksp Kep Kpp
E      *    E       *    0  0   0   0   0.5    0.5     0.0010000  0.0     0.0    
H      *    H       E    1  0   0   0   0.354  1.187   0.0006808  0.0     0.0       
He     *    He      E    0  0   0   0   0.849  1.4815  0.0009450  0.0     0.0        
Li     *    Li      E    1  0   0   0   1.336  1.2255  0.0010841  0.0     0.0      
Be     *    Be      E    2  0   0   0   1.074  1.3725  0.0036859  0.0     0.0       
B      *    B       E    3 -1   0   0   0.838  2.0415  0.0078054  0.0     0.0       
C      *    C       E    4  0   0   0   0.757  1.908   0.0037292  0.0     0.0      109.5  -19.5   6.32  0.0 0.0 0.0           
N      *    N       E    3  1   0   0   0.700  1.780   0.0073719  0.0     0.0      109.5  -19.5   5.60  1.0 0.0 0.5 
O      *    O       E    2  2   0   0   0.658  1.661   0.0091063  0.0     0.0      109.5  -19.5   5.99  1.0 0.0 0.5 
F      *    F       E    1  1   0   0   0.668  1.661   0.0091063  0.0     0.0      
Ne     *    Ne      E    0  0   0   0   0.920  2.0000  0.0010000  0.0     0.0       
Na     *    Na      E    1  0   0   0   1.539  1.300   0.0020000  0.0     0.0       
Mg     *    Mg      E    2  0   0   0   1.421  1.5105  0.0048133  0.0     0.0       
Al     *    Al      E    3 -1   0   0   1.244  2.2495  0.0218985  0.0     0.0      
Si     *    Si      E    4  0   0   0   1.117  2.1475  0.0174321  0.0     0.0      109.5  -19.5   6.32  0.0 0.0 0.0 
P      *    P       E    3  1   0   0   1.101  2.0735  0.0132258  0.0     0.0      109.5  -19.5   5.60  1.0 0.0 0.5
S      *    S       E    2  2   0   0   1.064  2.0175  0.0118816  0.0     0.0      109.5  -19.5   5.99  1.0 0.0 0.5    
Cl     *    Cl      E    1  1   0   0   1.044  1.948   0.0200000  0.0     0.0  
#
H_     H      H       E    1  0   0   0   0.354  1.487   0.0006808  0.0     0.0       
#
C_3    C      C       E    4  0   0   0   0.757  1.908   0.0037292  0.0     0.0      109.47 -19.47  6.32  1.6 1.0 0.0
C_2    C      C       E    4  0   0   0   0.732  1.908   0.0037292  0.0     0.0      120.0    0.0   4.21  1.6 1.0 2.0
C_1    C      C       E    4  0   0   0   0.706  1.908   0.0037292  0.0     0.0      180.0    0.0   3.88  1.6 1.0 0.0
C_R    C      C       E    4  0   0   0   0.729  1.908   0.0037292  0.0     0.0      120.0    0.0   6.34  1.6 1.0 1.3
C_CH3  C_3    C       E    4  0   0   0   0.757  1.908   0.0037292  0.0     0.0      109.5  -19.5   6.32  1.0 0.0 0.0
C_ene  C_2    C       E    4  0   0   0   0.757  1.908   0.0037292  0.0     0.0      120.0    0.0   5.94  2.2 1.0 1.0
C_yne  C_1    C       E    4  0   0   0   0.757  1.908   0.0037292  0.0     0.0      180.0    0.0   3.88  0.0 0.0 0.0
C_ald  C_2    C       E    4  0   0   0   0.757  1.908   0.0037292  0.0     0.0      120.0    0.0   5.94  2.2 1.0 1.0
C_COO  C_2    C       E    4  0   0   0   0.757  1.908   0.0037292  0.0     0.0      120.0    0.0   5.94  2.2 1.0 1.0
#
N_3    N      N       E    3  1   0   0   0.700  1.780   0.0073719  0.0     0.0      106.7  -22.1   5.60  1.6 1.0 0.4
N_2    N      N       E    3  1   0   0   0.685  1.780   0.0073719  0.0     0.0      111.2    0.0   7.43  1.6 1.0 0.9
N_1    N      N       E    3  1   0   0   0.656  1.780   0.0073719  0.0     0.0      180.0    0.0   4.82  1.6 1.0 0.0
N_R    N      N       E    3  1   0   0   0.699  1.780   0.0073719  0.0     0.0      120.0    0.0   7.43  1.6 1.0 1.3
#
O_3    O      O       E    2  2   0   0   0.658  1.661   0.0091063  0.0     0.0      104.51 -37.745 5.99  1.6 1.0 0.4
O_2    O      O       E    2  2   0   0   0.634  1.661   0.0091063  0.0     0.0      120.0    0.0   8.02  1.6 1.0 0.0
O_1    O      O       E    2  2   0   0   0.639  1.661   0.0091063  0.0     0.0      180.0    0.0   8.02  1.6 1.0 0.0
O_R    O      O       E    2  2   0   0   0.680  1.661   0.0091063  0.0     0.0      110.0    0.0   8.02  1.6 1.0 0.9
O_OH   O_3    O       E    2  2   0   0   0.658  1.661   0.0091063  0.0     0.0      109.5  -19.5   5.99  1.0 0.0 0.5
O_ald  O_2    O       E    2  2   0   0   0.658  1.661   0.0091063  0.0     0.0      120.0    0.5   8.02  2.2 0.0 1.0
O_sCOO O_3    O       E    2  2   0   0   0.658  1.661   0.0091063  0.0     0.0      109.5  -19.5   5.99  1.0 0.0 0.5
O_pCOO O_2    O       E    2  2   0   0   0.658  1.661   0.0091063  0.0     0.0      120.0    0.0   8.02  2.2 0.0 1.0
#
H_OH   H      H       E    1  0   0   0   0.354  1.487   0.0006808  0.0     0.0      
H_NH2  H      H       E    1  0   0   0   0.354  1.487   0.0006808  0.0     0.0   
H_CH3  H      H       E    1  0   0   0   0.354  1.487   0.0006808  0.0     0.0  
H_ene  H      H       E    1  0   0   0   0.354  1.487   0.0006808  0.0     0.0     
H_yne  H      H       E    1  0   0   0   0.354  1.487   0.0006808  0.0     0.0     
H_ald  H      H       E    1  0   0   0   0.354  1.487   0.0006808  0.0     0.0     
H_COO  H      H       E    1  0   0   0   0.354  1.487   0.0006808  0.0     0.0 

'''