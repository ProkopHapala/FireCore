
import numpy as np
import os
import sys
#from . import FireCore as fc
from . import atomicUtils as au
#from . import MMFF as mmff

verbosity=0

def makeLinearScan_firecore( nstep, selection, d, apos, nmax_scf=200 ):
    forces = np.zeros(apos.shape)
    Etemp  = np.zeros(8)
    Es     = np.zeros(nstep)
    for i in range(nstep):
        #fc.SCF( positions, iforce=iforce, nmax_scf=nmax_scf )
        fc.evalForce( apos, forces=forces, nmax_scf=nmax_scf, Es=Etemp, ixyz=i )
        print( "makeLinearScan_firecore() step# ", i," E[eV]= ", Etemp[0] )
        Es[i] = Etemp[0]
        apos[selection,:] += d
    return Es

def makeRotMat( ang ):
    ca=np.cos(ang)
    sa=np.sin(ang)
    return np.array([
        [1, 0,  0],
        [0,ca,-sa],
        [0,sa, ca]])

def makeRotationScan_firecore( nstep, selection, rot, p0, apos, nmax_scf=200 ):
    p0=np.array(p0)
    #ax=np.array(ax)
    #up=np.
    forces = np.zeros(apos.shape)
    Etemp  = np.zeros(8)
    Es     = np.zeros(nstep)
    for i in range(nstep):
        #fc.SCF( positions, iforce=iforce, nmax_scf=nmax_scf )
        fc.evalForce( apos, forces=forces, nmax_scf=nmax_scf, Es=Etemp, ixyz=i )
        print( "makeRotationScan_firecore() step# ", i," E[eV]= ", Etemp[0] )
        Es[i] = Etemp[0]
        apos[selection,:]-=p0
        for ia in selection:
            apos[ia,:] = np.dot( rot, apos[ia,:] )
        apos[selection,:]+=p0
    return Es

def linearScan( molFix, molMov, xs, dir=(1.,0,0), xyz_file="scan.xyz", Eout_file=None, callback=None ):
    dir=np.array(dir); dir/=np.linalg.norm(dir)
    if xyz_file is not None: fout = open(xyz_file,'w')
    Es = np.zeros(xs.shape)
    apos1,es1=molFix
    apos2,es2=molMov
    na1=len(es1); na2=len(es2); # na=na1+na2
    apos = np.concatenate( (apos1,apos2) )  #;print( apos )
    es   = np.concatenate( (es1,es2) )      #;print( es )
    for i,x in enumerate(xs):
        apos[na1:,:] = apos2 + dir[None,:]*x
        if callback is not None: Es[i] = callback( (apos,es) )
        #print( " out.mol._atom  \n", out.mol._atom )
        if xyz_file is not None:
            au.writeToXYZ( fout, es, apos, commet=(" x %10.5f E_tot %20.10f " %(x,Es[i]) ) )
    if xyz_file  is not None: fout.close()
    if Eout_file is not None: np.savetxt( Eout_file, np.array( [xs, Es] ).transpose() )
    return Es

def angularScan( n, dang, molFix, molMov, ang0=None, drot=None, ax1=0,ax2=1, xyz_file="scan.xyz", Eout_file=None, callback=None ):
    if xyz_file is not None: fout = open(xyz_file,'w')
    if drot is None: drot = au.makeRotMatAng( dang, ax1=ax1, ax2=ax2 ).transpose()
    Es   = np.zeros(n)
    angs = np.zeros(n)
    apos1,es1=molFix
    apos2,es2=molMov
    na1=len(es1); na2=len(es2); # na=na1+na2
    apos = np.concatenate( (apos1,apos2) )  #;print( apos )
    es   = np.concatenate( (es1,es2) )      #;print( es )
    apos_ = apos2.copy()
    if ang0 is not None:
        rot0 = au.makeRotMatAng( ang0, ax1=1, ax2=2 ).transpose()
        au.mulpos(apos_,rot0)
    else:
        ang0=0
    for i in range(n):
        angs[i]=dang*i+ang0
        apos[na1:,:]=apos_[:,:]
        if callback is not None: Es[i] = callback( (apos,es) )
        #print( " out.mol._atom  \n", out.mol._atom )
        if xyz_file is not None:
            au.writeToXYZ( fout, es, apos, commet=(" ang %10.5f E_tot %20.10f " %(angs[i],Es[i]) ) )
        au.mulpos(apos_,drot)
    if xyz_file is not None: fout.close()
    if Eout_file is not None: np.savetxt( Eout_file, np.array( [angs, Es] ).transpose() )
    return Es,angs

def scan_xyz( fxyzin, fxyzout="out.xyz", Eout=None, callback=None ):
    fin =open(fxyzin,'r')
    if fxyzout is not None: fout=open(fxyzout,'w')
    Es=[]
    i=0
    while True:
        apos,Zs,es,qs = au.loadAtomsNP( fin=fin, bReadN=True )   #;print(apos) ;print(es)
        if(len(es)==0): break
        if callback is not None: Es.append(  callback((apos,es)) )
        comment =  "i %i E_tot %20.10f " %( i, Es[i] )
        if verbosity>0: print(comment)
        if fxyzout is not None:
            au.writeToXYZ( fout, es, apos, comment=comment )
        i+=1
    if fxyzout  is not None: fout.close()
    if Eout     is not None: np.savetxt( Eout, np.array( [i, Es] ).transpose() )
    return Es, range(len(Es))