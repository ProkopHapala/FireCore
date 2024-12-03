import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

sys.path.append("../../")
from pyBall import FitREQ as fit
from pyBall import atomicUtils as au

# ============== Setup

fname = "scan_OH"

imodel = 3
Hartree2eV   = 27.2114079527
kcal2eV      =  0.0433634

# ------ Set strength of hydrogen bond correction
#Hfactor = 0.99
Hfactor = 0.9
#Hfactor = 0.95  

# ------ Set atom types
types = [
    ([     1.4,    np.sqrt(0.0006808),     0,   +np.sqrt(0.0006808)*Hfactor ], [0,0,0,1] ),    # 0 H H2O
    ([     1.6,    np.sqrt(0.0091063),     0,   -np.sqrt(0.0091063)*Hfactor ], [0,0,0,0] ),    # 2 O H2O
]
typeMask = np.array( [ t[1] for t in types ], dtype=np.int32  )
typREQs  = np.array( [ t[0] for t in types ], dtype=np.double )

# ============== Functions

def getE_BuckinghamQ_fake( r, Rij, Eij, qq=0.0, H=0.0, kMorse=1.6 ):
    C6 = 2*Eij * (Rij**6)
    B  =   Eij*np.exp(kMorse*Rij)
    E = B*np.exp(-kMorse*r) - C6/r**6 + 14.399644*qq/r
    return E

def getE_BuckinghamQ_true( r, Rij, Eij, qq=0.0, H=0.0, kMorse=1.6 ):
    # https://en.wikipedia.org/wiki/Buckingham_potential#Modified_Buckingham_(Exp-Six)_potential
    a = kMorse*Rij
    C6 = (a/(a-6))*Eij * (Rij**6)
    B  = Eij*(6/(a-6))*np.exp(kMorse*Rij)
    E = B*np.exp(-kMorse*r) - C6/r**6 + 14.399644*qq/r
    return E

def getE_LJQ( r, Rij, Eij, qq=0.0, H=0.0 ):
    C6  = 2* Eij * (Rij**6)
    C12 =    Eij * (Rij**12)
    E   = C12/r**12 - C6/r**6 + 14.399644*qq/r
    return E

def genTestFile( fname, n=100,  Rmin=1.5, Rmax=4.0, Rij=2.0, Eij=np.sqrt(0.0006808*0.0091063), qij=0, H=0, kMorse=6.0 ):
    # np.sqrt(0.0006808*0.0091063) = 0.001648
    x = np.linspace( Rmin, Rmax, n )
    f  = open(fname+".xyz",'w')
    f2 = open(fname+".dat",'w')
    Es = getE_BuckinghamQ_true( x, Rij, Eij, qij, H, kMorse )
    #Es = getE_LJQ( x, Rij, Eij, qij, H )
    for i in range(n):
        comment = "# r   %f ang    0.00000 E_tot         %f" %(x[i], Es[i] )
        #comment = "R %f E= %f" %(x[i], Es[i] )
        f.write( "%i\n" %(2) )
        f.write( comment+"\n" )
        f.write( "H %f %f %f\n" %(0.0,0.0,0.0)  )
        f.write( "O %f %f %f\n" %(x[i],0.0,0.0) )
        f2.write( comment+"\n" )
    f.close()
    f2.close()

def plot_LJ_vs_BK(  n=100,  Rmin=1.5, Rmax=4.0, Rij=2.0, Eij=1, qij=0, H=0, kMorse=6.0 ):
    x = np.linspace( Rmin, Rmax, n )
    #Es_Bk = getE_BuckinghamQ_fake( x, Rij, Eij, qij, H, kMorse )
    Es_Bk = getE_BuckinghamQ_true( x, Rij, Eij, qij, H, kMorse )
    Es_LJ = getE_LJQ( x, Rij, Eij, qij, H )
    plt.plot( x, Es_Bk, '-k', label="BuckinghamQ" )
    plt.plot( x, Es_LJ, '-r', label="LJQ" )

def load_dat( fname ):
    dat   = np.genfromtxt( fname, comments='@' )   ;print( "dat[1,:] ", dat[1,:] )
    Es  =  dat[:,6]
    rs    =  dat[:,2]
    angs  =  dat[:,4]
    #Es/=0.0433634
    return Es,rs,angs

# ============== Main

'''
plot_LJ_vs_BK(  n=100,  Rmin=1.5, Rmax=4.0, Rij=2.0, Eij=1, qij=0, H=0, kMorse=4.5 )
plt.grid()
plt.legend()
plt.ylim( -2.0, 2.0 )
plt.show()
'''




genTestFile( fname )
Eref,rs,angs = load_dat( fname+".dat" )

fit.setVerbosity(1)
fit.init_types( typeMask, typREQs, bCopy=True ) 
fit.loadXYZ( "scan_OH.xyz", [0], [1], types0=[0], testtypes=[1]  )     # load reference geometry

fit.getBuffs()
print( "typToREQ\n" , fit.typToREQ )
print( "typeREQs\n" , fit.typeREQs )
print( "types1\n"   , fit.types1   )
#print( "types2"   , fit.types2   )
#print( "types3"   , fit.types3   )

# ------ obtain energy profile from classical model (fit)
#Es     = fit.getEs( imodel=imodel, bRigid=False)     #;print( "Es_noH", Es     )
Es     = fit.getEs( imodel=imodel)     #;print( "Es_noH", Es     )

# ------ obtain energy progile with HBond correction set to zero
#typREQs[0,3] = 0.0              # set HBond correction to zero
#fit.setType(0, typREQs[0,:] )   # set atom type 0
#Es_noH = fit.getEs( imodel=imodel, bRigid=False)     #;print( "Es_noH", Es_noH )

print( "rs.shape", rs.shape, "Es.shape", Es.shape, "Eref.shape", Eref.shape ) 
plt.plot( rs, Eref  , '-k', label="E_ref", lw=3 )
plt.plot( rs, Es    , '-g', label="E",     lw=2 )
#plt.plot( rs, Es_noH, '-r', label="E_noH", lw=2 )
plt.axhline(0.0, c='k', ls='--')
plt.grid()
plt.legend()
plt.ylim( -0.03, 0.03 )
#plt.ylim( -2.0, 2.0 )


plt.show()
