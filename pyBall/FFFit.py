
import numpy as np
import os
import sys
#from . import FireCore as fc
from . import atomicUtils as au
#from . import MMFF as mmff

verbosity=0

def exprange( xmax, n=10 ):
    return ((1.+xmax)**(1./(n-1)))**np.arange(0,n)-1.0

def scan_ranges( xs=[-0.4,0.4,2.0,6.0,10.0], dxs=[0.1,0.2,0.4,1.0], sigma=1e-9 ):
    ars = []
    n = len(dxs)
    for i in range(n):
        d = -sigma
        if i==n-1: d=sigma
        ars.append( np.arange(xs[i],xs[i+1]+d,dxs[i]) )
    return np.concatenate( ars )

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
            au.writeToXYZ( fout, es, apos, comment=(" x %10.5f E_tot %20.10f " %(x,Es[i]) ) )
    if xyz_file  is not None: fout.close()
    if Eout_file is not None: np.savetxt( Eout_file, np.array( [xs, Es] ).transpose() )
    return Es

def linearScan_1mol( mol, selection, xs, dir=(1.,0,0), xyz_file="scan.xyz", Eout_file=None, callback=None ):
    dir=np.array(dir); dir/=np.linalg.norm(dir)
    if xyz_file is not None: fout = open(xyz_file,'w')
    Es = np.zeros(xs.shape)
    apos,es=mol
    apos0 = apos.copy()
    na=len(es)
    for i,x in enumerate(xs):
        apos[:,:] = apos0 
        apos[selection,:] += dir[None,:]*x
        if callback is not None: Es[i] = callback( (apos,es) )
        if xyz_file is not None:
            au.writeToXYZ( fout, es, apos, comment=(" x %10.5f E_tot %20.10f " %(x,Es[i]) ) )
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
            au.writeToXYZ( fout, es, apos, comment=(" ang %10.5f E_tot %20.10f " %(angs[i],Es[i]) ) )
        au.mulpos(apos_,drot)
    if xyz_file is not None: fout.close()
    if Eout_file is not None: np.savetxt( Eout_file, np.array( [angs, Es] ).transpose() )
    return Es,angs

def angularScan_1mol(  mol, selection, rs, angs, ax1=0,ax2=1, dir=None, xyz_file="scan.xyz", Eout_file=None, callback=None ):
    if xyz_file is not None: fout = open(xyz_file,'w')
    Es   = np.zeros( (len(rs),len(angs)) )
    if(dir is None): 
        dir=[0.,0.,0.]; dir[ax1]=1.0; dir=np.array(dir)
    else:
        dir=np.array(dir); dir*=(1/np.sqrt(np.dot(dir,dir)))
    apos,es = mol 
    apos0 = apos.copy()
    apos  = apos.copy()
    for i,r in enumerate(rs):
        #apos_i = apos0.copy()
        #apos_i[selection,:] += dir[None,:]*r
        aposi = apos0[selection,:] + dir[None,:]*r
        for j,ang in enumerate(angs):
            aposj = aposi.copy()
            #rot = au.makeRotMatAng( ang, ax1=1, ax2=2 ).transpose()
            rot = au.makeRotMatAng( ang, ax1=ax1, ax2=ax2 ).transpose()
            au.mulpos(aposj,rot)
            apos[selection,:] = aposj
            if callback is not None: Es[i] = callback( (apos,es) )
            #print( " out.mol._atom  \n", out.mol._atom )
            if xyz_file is not None:
                comment = ("# r %10.5f ang %10.5f " %(rs[i],angs[j]) )
                if callback is not None:
                    comment = " E_tot %20.10f" %Es[i,j]
                au.writeToXYZ( fout, es, apos, comment=comment )
    if xyz_file is not None: fout.close()
    if Eout_file is not None: np.savetxt( Eout_file, np.array( [angs, Es] ).transpose() )
    return Es,angs

def angularScan_1mol_vecs(  mol, selection, rs, angs, dir=(1.,0.,0.), up=(0.0,1.0,0.0), xyz_file="scan.xyz", Eout_file=None, callback=None ):
    if xyz_file is not None: fout = open(xyz_file,'w')
    Es   = np.zeros( (len(rs),len(angs)) )
    #dir=[0.,0.,0.]; dir[ax2]=1.0; dir=np.array(dir)
    dir=np.array(dir); up=np.array(up)  #dir[ax2]=1.0; dir=np.array(dir)
    apos,es = mol 
    apos0 = apos.copy()
    apos  = apos.copy()
    print( selection )
    for i,r in enumerate(rs):
        #apos_i = apos0.copy()
        #apos_i[selection,:] += dir[None,:]*r
        aposi = apos0[selection,:] + dir[None,:]*r
        for j,ang in enumerate(angs):
            aposj = aposi.copy()
            #rot = au.makeRotMatAng( ang, ax1=1, ax2=2 ).transpose()
            rot = au.makeRotMatAng2( dir, up, ang ).transpose()
            #rot = au.makeRotMatAng( ang, ax1=ax1, ax2=ax2 ).transpose()
            au.mulpos(aposj,rot)
            apos[selection,:] = aposj
            if callback is not None: Es[i] = callback( (apos,es) )
            #print( " out.mol._atom  \n", out.mol._atom )
            if xyz_file is not None:
                au.writeToXYZ( fout, es, apos, comment=("r %10.0f ang %10.5f E_tot %20.10f " %(rs[i],angs[j],Es[i,j]) ) )
    if xyz_file is not None: fout.close()
    if Eout_file is not None: np.savetxt( Eout_file, np.array( [angs, Es] ).transpose() )
    return Es,angs

def scan_xyz( fxyzin, fxyzout="out.xyz", Eout=None, callback=None, params=None, xs=None ):
    fin =open(fxyzin,'r')
    if fxyzout is not None: fout=open(fxyzout,'w')
    Es=[]
    i=0
    if xs is None: xs=range(i)
    while True:
        comments=[]
        apos,Zs,es,qs = au.loadAtomsNP( fin=fin, bReadN=True, comments=comments )  #;print(apos) ;print(es)
        if(len(es)==0): break
        E = 0.0
        if callback is not None:  E = callback( (apos,es), params=params, id=i )
        #comment = ( "i %i E_tot %20.10f x " %( i, Es[i] ) ) + str(xs[i])
        Es.append( E )
        #print( "comments "+str(i)+" "+ str(comments) )
        comment = comments[0].strip() + ( " E_tot %20.10f" %E )
        if verbosity>0: print(comment)
        if fxyzout is not None:
            au.writeToXYZ( fout, es, apos, comment=comment )
        i+=1
    if fxyzout  is not None: fout.close()
    if Eout     is not None: np.savetxt( Eout, np.array( [i, Es] ).transpose() )
    return Es, xs