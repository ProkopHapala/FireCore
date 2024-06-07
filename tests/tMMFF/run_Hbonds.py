import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

import ctypes

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF as mmff


# =========== Setting

nmaxiter = 10000
#xyz_name = "butandiol-2"
#xyz_name = "polymer-2_new-OH"
xyz_name = "formic_dimer"

ref_path = "/home/prokop/Desktop/CARBSIS/PEOPLE/Mithun/HBond_Fit_Radial_Scans/energy_angle_0/"

#  Charges :
#  NH3    +0.28363  -0.84922      
#  H2O    +0.34068  -0.68135
#  HF     +0.41277  -0.41277    

# =========== Functions

def makeWeights( xs, Es, wfac=10.0, dimin=3, bRenorm=True ):
    imin = np.argmin(Es)
    i0 = imin - dimin
    if i0<0: i0=0 
    ws   = np.exp( -wfac*Es )
    ws[:i0] = 0.0
    if bRenorm:
        ws *= 1./ws.sum()
    return ws

def compare_profile( typstr, fname, qH, qX, xshift=-0.3, clr=None, bPlt=True, bPltW=True, errSc=10.0 ):
    ref   = np.genfromtxt(ref_path+fname ); 
    ref[:,0] += xshift
    rs    = ref[:,0].copy(); 
    Es,Fs = mmff.sampleNonBondTypes( typstr, rs, qH=qH, qX=qX );    
    Ediff =  Es -  ref[:,1]
    ws = makeWeights( ref[:,0], ref[:,1], wfac=10.0, dimin=3 );  
    error = Ediff * ws
    RMS = np.sqrt((error**2).sum())
    if bPlt:
        plt.plot( rs,             Es, '.-', c=clr,  label=typstr   )
        plt.plot( ref[:,0], ref[:,1], ".:", c=clr, label=fname   )
        if bPltW:
            plt.plot( ref[:,0], ws,":",    c='gray', label='weights' )
            plt.plot( ref[:,0], error*errSc,"-", c='gray', label='weights' )
    print( typstr, " RMS = ", RMS );
    return RMS  


#ialg  = 2
#ialg  = 3
#alg_names=["GD","MDdamp","FIRE","FIREsmooth"]


# print( "offset(RvdW): ", mmff.AtomType.__dict__['RvdW'].offset )
# print( "offset(RvdW): ", mmff.AtomType.__dict__['EvdW'].offset )
# print( "offset(RvdW): ", mmff.AtomType.__dict__['Qbase'].offset )
# print( "offset(RvdW): ", mmff.AtomType.__dict__['Hb'].offset )
#exit()


#ctypes.c_double.from_address(ctypes.addressof(atom) + field_offset).value
#ctypes.c_double.from_address(ctypes.addressof(atom) + field_offset).value = value

def mutateInRange( xs, ranges, lims ):
    dx = np.random.rand( len(xs) )
    dx = dx*( ranges[:,1] - ranges[:,0] ) + ranges[:,0]
    #print( "dx ", dx )
    xs += dx
    #print( "xs: ", xs, lims[:,0], lims[:,1] )
    mask = xs<lims[:,0];   xs[mask] = lims[mask,0]
    mask = xs>lims[:,1];   xs[mask] = lims[mask,1]
    return xs

def getVals( indexes, atyps ):
    n = len(indexes)
    xs = np.zeros(n)
    for i in range( n ):
        ityp,ioff = indexes[ i ]
        xs[i] = ctypes.c_double.from_address(ctypes.addressof(atyps[ityp]) + ioff).value
    return xs

def setVals( xs, indexes, atyps ):
    n = len(xs)
    for i in range(n):
        ityp,ioff = indexes[ i ]
        ctypes.c_double.from_address(ctypes.addressof(atyps[ityp]) + ioff).value = xs[i]

def getOffsets(what):
    #np.zeros( (len(what),2),dtype=np.int32)
    return [ (i,mmff.AtomType.__dict__[s].offset) for i,s in what  ]

lims    = np.array( [ [+0.6,+0.999], [-0.999,-0.6],  [0.0,0.1]  ] )
rndlims = np.array( [ [-0.01,0.01],  [-0.01,0.01],   [0.0,0.0]  ] )
what = [
    [42,'Hb'],
    [11,'Hb'],
    [11,'EvdW'],
]

def run_rnd_opt( typstr, fname, nitr = 100, qH=0.2, qX=-0.2, xshift=-0.3, clr=None, bPlt=True, bPltW=True, errSc=10.0 ):
    atyps, ntyp = mmff.getAtomTypes() 
    indexes     = getOffsets(what)
    ref         = np.genfromtxt(ref_path+fname ); 
    ref[:,0]   += xshift
    rs          = ref[:,0].copy(); 
    ws          = makeWeights( ref[:,0], ref[:,1], wfac=10.0, dimin=3 ); 
    plt.plot( ref[:,0], ref[:,1], ".:", c=clr, label=fname   )
    if bPltW:
        plt.plot( ref[:,0], ws,":",    c='gray', label='weights' )
    RMS_old = 1e+300
    xs = getVals(  indexes, atyps, )
    print( "xs.0:", xs )
    for itr in range(nitr):
        if itr>0:
            x_new = xs
            x_new = mutateInRange( xs, rndlims, lims )
            #print( "x_new ", x_new )
            setVals( x_new, indexes, atyps )
        Es,Fs = mmff.sampleNonBondTypes( typstr, rs, qH=qH, qX=qX );    
        Ediff =  Es -  ref[:,1]
        error = Ediff * ws
        RMS  = np.sqrt((error**2).sum())
        if (itr>0) and (RMS < RMS_old ):
            xs = x_new
            print("ACCEPT  RMS=%g <(%g)" %(RMS, RMS_old), xs )
            RMS_old = RMS
            if bPlt:
                plt.plot( rs,             Es, '.-', c=clr,  label=("itr=%i RMS(%g)" %(itr,RMS) )   )
        else:
            #print("REJECT  RMS=%g >(%g)" %(RMS, RMS_old) )
            pass
        # if bPlt:
        #     plt.plot( rs,             Es, '.-', c=clr,  label=("itr=%i RMS(%g)" %(itr,RMS) )   )
        #     #if bPltW:
        #     #    plt.plot( ref[:,0], error*errSc,"-", c='gray', label='err' )



#======== Body

#alg_name = alg_names[ialg]
mmff.setVerbosity( verbosity=2, idebug=0 )
mmff.init( xyz_name="data/"+xyz_name, bMMFF=True )     
#mmff.getBuffs()
#mmff.eval()
#rs = np.linspace(1.5,20.0, 100 )

atyps, ntyp = mmff.getAtomTypes() #  ;print( ntyp, atyps[5].name.decode('utf-8'), atyps[5].RvdW  )

# ---- print known types ( check if mapping works
for i in range( ntyp ):
    a = atyps[i]
    print( i,a.name.decode('utf-8'), a.RvdW, a.EvdW, a.Qbase,  a.Hb  )

#exit()

run_rnd_opt(  "H_F   F",   'ref-HF_HF-000.dat',  nitr = 100,  qH=+0.41277, qX=-0.41277 )

#compare_profile( "H_F   F",   'ref-HF_HF-000.dat',   qH=+0.41277, qX=-0.41277, clr='g' ) 
#atyps[11].Hb = 0.8
#compare_profile( "H_F   F",   'ref-HF_HF-000.dat',   qH=+0.41277, qX=-0.41277, clr='cyan' ) 

#compare_profile( "H_OH  O_3", 'ref-H2O_H2O-000.dat', qH=+0.34068, qX=-0.34068, clr='r' ) 
#compare_profile( "H_NH2 N_3", 'ref-NH3_NH3-000.dat', qH=+0.28363, qX=-0.28363, clr='b' ) 


plt.legend()
plt.axhline(0,ls='--',c='k')
plt.ylim(1.5,20.0);   # plt.xticks(np.arange(1,5, 20.5, 0.2)) 
plt.ylim(-0.3,0.3)
plt.minorticks_on()
plt.grid(which='both')
#plt.ylim(-0.1,0.1)
#plt.ylim(-0.03,0.03)

plt.show()


print("ALL DONE")
#plt.show()



'''
clrs = [ 'r','g','b','m','c' ]

nb = mmff.findHbonds( Rcut=4.0, angMax=30.0 );
for ib in range( nb )[:1]:
    c = clrs[ib]
    Es,Fs,s = mmff.sampleHbond( ib, rs                   , kind=1            );    plt.plot( rs, Es - Es[-1], '-k',  label="Hb #"+str(ib) + s         )
    Es,Fs,s = mmff.sampleHbond( ib, rs                   , kind=2, dcomp=3.0 );    plt.plot( rs, Es - Es[-1], '-k',  label="Hb #"+str(ib) + s         )
    Es,Fs,s = mmff.sampleHbond( ib, rs                   , kind=2, dcomp=1.5 );    plt.plot( rs, Es - Es[-1], '--k', label="Hb #"+str(ib) + s         )
    #Es,Fs,s = mmff.sampleHbond( ib, rs,           maskH=0.0  );    plt.plot( rs, Es - Es[-1], c=c, ls='--k', label="Hb #"+str(ib)+" noHb" )
    Es,Fs,s = mmff.sampleHbond( ib, rs, maskQ=0.0            );    plt.plot( rs, Es - Es[-1], ':k',  label="Hb #"+str(ib)+" noQ |"+s  )
    #Es,Fs,s = mmff.sampleHbond( ib, rs, maskQ=0.0, maskH=0.0 );    plt.plot( rs, Es - Es[-1], c=c, ls=':k',  label="Hb #"+str(ib)+" noQ |"+s  )
    #print ,s
'''