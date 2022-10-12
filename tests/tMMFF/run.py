import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF as mmff

def deriv( xs,Es ):
    xs_ =(xs[2:]+xs[:-2])*0.5
    dxs =(xs[2:]-xs[:-2])
    dEs =(Es[2:]-Es[:-2])
    return dEs/dxs,xs_

#======== Body

mmff.init()
# -------- Non-covalent Surace (GridFF)
Q =1.00; K=-1.0; kind=1; Rdamp=1.0
sc=1.2
rs = np.linspace(-5.0,10.0,150)
#EsM,fsM  = mmff.sampleNonBond( rs, kind=kind, REQi=(R,e  ,Q*0), K=K, Rdamp=Rdamp)  
#EsC,fsC  = mmff.sampleNonBond( rs, kind=kind, REQi=(R,e*0,Q  ), K=K, Rdamp=Rdamp)   
#Es ,fs   = mmff.sampleSurf( "data/NaCl_sym-center", rs, kind=kind, atyp=0, Q=Q, K=K, Rdamp=Rdamp)  
Es,fs   = mmff.sampleSurf( "data/NaCl_sym-center", rs, kind=12, atyp=1, Q=Q, K=K, Rdamp=Rdamp, pos0=(2.0,0.0,0.0), bSave=True)  

#print("Es \n", Es);

fnum,xfs = deriv( rs,Es )

#Emin=Es.min()
Emax=Es.max()
Emin=Es.min()
plt.figure(figsize=(5,10))
plt.subplot(2,1,1); plt.plot(rs,Es,'k',label='PLQ');  plt.grid(); plt.axhline(0,c='k',ls='--');                                           #plt.ylim(Emin*sc,Emax*sc)  ;plt.legend()
plt.subplot(2,1,2); plt.plot(rs,fs,    label='Fana'); plt.plot(xfs,-fnum, ':r',label='Fnum')    ;plt.grid(); plt.axhline(0,c='k',ls='--'); #plt.ylim(Emin*sc,Emax*sc)  ;plt.legend()
plt.show()
exit()


'''
# -------- Non-covalent
Q =-0.05; R=1.908; e=0.0037292; K=-1.0; kind=1; Rdamp=1.0
sc=1.2
rs = np.linspace(0.0,10.0,100)
EsM,fsM  = mmff.sampleNonBond( rs, kind=kind, REQi=(R,e  ,Q*0), REQj=(R,e,abs(Q)), K=K, Rdamp=Rdamp)  
EsC,fsC  = mmff.sampleNonBond( rs, kind=kind, REQi=(R,e*0,Q  ), REQj=(R,e,abs(Q)), K=K, Rdamp=Rdamp)   
Es ,fs   = mmff.sampleNonBond( rs, kind=kind, REQi=(R,e  ,Q  ), REQj=(R,e,abs(Q)), K=K, Rdamp=Rdamp)  

fnum,xfs = deriv( rs,Es )

#Emin=Es.min()
Emin=e*-3.0
plt.figure(figsize=(5,10))
plt.subplot(2,1,1); plt.plot(rs,Es,'k',label='PLQ'); plt.plot(rs,EsM,label="PL"); plt.plot(rs,EsC,label="Q")   ;plt.grid(); plt.axhline(0,c='k',ls='--'); plt.ylim(Emin*sc,-Emin*sc)  ;plt.legend()
plt.subplot(2,1,2); plt.plot(rs,fs,    label='Fana'); plt.plot(xfs,fnum, ':r',label='Fnum')                    ;plt.grid(); plt.axhline(0,c='k',ls='--'); plt.ylim(Emin*sc,-Emin*sc)  ;plt.legend()
plt.show()
exit()
'''


#mmff.init_params( "data/AtomTypes.dat", "data/BondTypes.dat", "data/AngleTypes.dat" )
#mmff.insertSMILES("CC");
#mmff.insertSMILES("C=C");
#mmff.insertSMILES("C#C");
#mmff.insertSMILES("C#CCN=C", True );
#mmff.insertSMILES("C1#CCN=C1", True );
#mmff.insertSMILES("C=C1NC#CC1CO", True, True );


#mmff.initWithSMILES( "C=C1NC#CC1CO" )
mmff.initWithSMILES( "C=C" )
mmff.getBuffs()
mmff.relax(1000)
mmff.plot()
plt.show()

exit()


'''
# ======== Oritent Molecule
xyzs,Zs,enames,qs = au.loadAtomsNP( "data/Benzene_deriv.xyz" )
au.orient( 2, (5,2), (1,3), xyzs, bFlipXZ=True )
au.saveXYZ( enames, xyzs, "data/Benzene_deriv_.xyz", qs=qs, Rs=None )
plt.plot( xyzs[:,0],xyzs[:,1], "o" )
plt.axis('equal')
plt.show()
exit()
'''

'''
# ============== C2H4,xyz
#mmff.initWithMolFile( "C2H4.xyz", bNonBonded=False, bOptimizer=True)
#mmff.printBuffNames()
#mmff.getBuffs() #;print( mmff.ndims )
#mmff.eval()
#mmff.relax(1000, bWriteTrj=True )
#Es=mmff.scanRotation( [1,4,5], 0, 0,1, np.pi*2, 100, bWriteTrj=True)   ;print("Es=", Es)
#plt.plot(Es)
#print( "Es(Etot,Eb,Ea,Eps,EppT,EppI):", mmff.Es )
#nsel = mmff.splitAtBond(6-1)  ;print( "split to:\n", mmff.selection[:nsel],"\n", mmff.selection[nsel:] )
'''

# ============== Benzene_deriv.xyz
mmff.initWithMolFile( "data/Benzene_deriv.xyz", bNonBonded=False, bOptimizer=True)
mmff.getBuffs() #;print( mmff.ndims )

#nsel = mmff.splitAtBond(5)  ;print( "split to:\n", mmff.selection[:nsel],"\n", mmff.selection[nsel:] )
#nsel = mmff.splitAtBond(6)  ;print( "split to:\n", mmff.selection[:nsel],"\n", mmff.selection[nsel:] )
#nsel = mmff.splitAtBond(10)  ;print( "split to:\n", mmff.selection[:nsel],"\n", mmff.selection[nsel:] )
#nsel = mmff.splitAtBond(2)  ;print( "split to:\n", mmff.selection[:nsel],"\n", mmff.selection[nsel:] )
#nsel = mmff.splitAtBond(4)  ;print( "split to:\n", mmff.selection[:nsel],"\n", mmff.selection[nsel:] )
#print( "nsel ", nsel, len(mmff.selection)-nsel )
#Es=mmff.scanRotation( 1, 1,11, np.pi*2, 100, bWriteTrj=True, _0=1) ;plt.plot(Es) ;print("Es=", Es)

#Es = mmff.scanBondRotation( 6, np.pi*2, 100, bWriteTrj=True );  plt.plot(Es)
Es = mmff.scanBondRotation( 2, np.pi*2, 100, bWriteTrj=True );  plt.plot(Es); plt.grid()


#mmff.eval()
#mmff.relax(1000, Ftol=1e-4, bWriteTrj=True )
#Es=mmff.scanRotation( 1, 1,11, np.pi*2, 100, sel=[11,13,14,20]+[29,30,31,32], bWriteTrj=True, _0=1) ;plt.plot(Es) ;print("Es=", Es)

plt.figure()
mmff.plot()
#mmff.plot_selection( mmff.selection[:nsel] )
#mmff.plot_selection( [1,2,3] )

plt.show()