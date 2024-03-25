import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF as mmff

def getLJ( x,y,z, R0, E0 ):
    r = np.sqrt(x**2 + y**2 + z**2)
    E  = E0*    ( (R0/r)**12 - 2*(R0/r)**6 )
    fr = E0*-12*( (R0/r)**12 -   (R0/r)**6 )/(r*r)
    return E,fr*x,fr*y,fr*z

def getLJ_2( x,y,z, R0, E0 ):
    R6   = R0**6
    C6   = E0*R6
    C#12  = E0*R6*R6
    ir2  = 1/(x**2 + y**2 + z**2)
    ir   = np.sqrt(ir2) 
    E    = C6*   (      R6*ir**6 -   2*ir2**3 )
    dE   = C6*12*(      R6*ir**6 -     ir2**3 )*ir2
    ddE  = C6*12*(  168*R6*ir**6 -  96*ir2**3 )*ir2*ir
    dddE = C6*12*( 2688*R6*ir**6 - 960*ir2**3 )*ir2*ir2
    return E,dE*x,dE*y,dE*z, ddE*x*y,ddE*x*z,ddE*y*z, dddE*x*y*z

def numDeriv( xs, Es): 
    dx = xs[1]-xs[0]
    Fs = (Es[2:]-Es[:-2])/(2*dx)
    return Fs

def makeGrid( atoms, ng, g0, dg ):
    xs = np.arange(ng[0])*dg[0] + g0[0]   #;print("xs: ",xs)
    ys = np.arange(ng[1])*dg[1] + g0[1] 
    zs = np.arange(ng[2])*dg[2] + g0[2] 
    Xs,Ys,Zs = np.meshgrid( xs,ys,zs )
    Eg = np.zeros( ng )
    for i,a in enumerate(atoms):
        X = Xs-a[0]
        Y = Ys-a[1]
        Z = Zs-a[2]
        Eg += a[4]*( 1/( X**2 + Y**2 + Z**2 + a[3]**2 ) )
    return Eg

def makeGrid_deriv_dir( atoms, ng, g0=(0.0,0.0,0.0), dg=(0.0,0.0,0.1) ):
    g0=np.array(g0)
    dg=np.array(dg)
    lg =np.sqrt(np.dot(dg,dg))
    Xs = g0[0] + np.arange(0, ng) * dg[0]
    Ys = g0[1] + np.arange(0, ng) * dg[1]
    Zs = g0[2] + np.arange(0, ng) * dg[2]
    ls =  np.sqrt( (Xs-g0[0])**2 + (Ys-g0[0])**2 + (Zs-g0[0])**2 )
    FE = np.zeros( (ng,2) )
    #print( FE.shape, Xs.shape, Ys.shape, Zs.shape )
    for i,a in enumerate(atoms):
        X = Xs-a[0]
        Y = Ys-a[1]
        Z = Zs-a[2]
        R2  = X**2 + Y**2 + Z**2 + a[3]**2;
        #iR2 = 1/R2
        #iR4 = iR2*iR2

        E,fx,fy,fz = getLJ( X,Y,Z, 3.0, 1.0 )
        FE[:,0]  = E
        FE[:,1]  = ( fx*dg[0] + fy*dg[1] + fz*dg[2] )/lg 

        #FE[:,0] +=    a[4]*iR2                                        # energy
        #FE[:,1] += -2*a[4]*iR4*( X*dg[0] + Y*dg[1] + Z*dg[2] )/lg   # directional derivative


        #FE[:,0] += Z**2
        #FE[:,1] += 2*Z
    return FE, ls

def makeGrid_deriv( atoms, ng, g0, dg ):
    xs = np.arange(ng[0])*dg[0] + g0[0]   #;print("xs: ",xs)
    ys = np.arange(ng[1])*dg[1] + g0[1] 
    zs = np.arange(ng[2])*dg[2] + g0[2] 
    Xs,Ys,Zs = np.meshgrid( xs,ys,zs )
    Eg = np.zeros( ng+[4,] )
    for i,a in enumerate(atoms):
        X = Xs-a[0]
        Y = Ys-a[1]
        Z = Zs-a[2]
        R2  = X**2 + Y**2 + Z**2 + a[3]**2;
        iR2 = 1/R2
        iR4 = iR2*iR2
        Eg[:,:,:,0] += a[4]*X*iR4
        Eg[:,:,:,1] += a[4]*Y*iR4
        Eg[:,:,:,1] += a[4]*Z*iR4
        Eg[:,:,:,3] += a[4]  *iR2
    return Eg

def make2dDeriv( FE, dg ):
    i3 = 1./3.
    dFE = np.zeros( FE.shape )
    dFE[:,:,:,0] =  ( np.roll(FE[:,:,:,2], -1, axis=1) - FE[:,:,:,2] )*(0.5/dg[1])   + ( np.roll(FE[:,:,:,1], -1, axis=2) - FE[:,:,:,1] )*(0.5/dg[2])        # yz
    dFE[:,:,:,1] =  ( np.roll(FE[:,:,:,0], -1, axis=1) - FE[:,:,:,0] )*(0.5/dg[1])   + ( np.roll(FE[:,:,:,1], -1, axis=0) - FE[:,:,:,1] )*(0.5/dg[0])        # xy
    dFE[:,:,:,2] =  ( np.roll(FE[:,:,:,0], -1, axis=2) - FE[:,:,:,0] )*(0.5/dg[2])   + ( np.roll(FE[:,:,:,2], -1, axis=0) - FE[:,:,:,2] )*(0.5/dg[0])        # xz
    dFE[:,:,:,2] = (       # xyz
        ( np.roll(dFE[:,:,:,0], -1, axis=0) - dFE[:,:,:,0] )*(i3/dg[0])  +
        ( np.roll(dFE[:,:,:,1], -1, axis=2) - dFE[:,:,:,1] )*(i3/dg[2])  +
        ( np.roll(dFE[:,:,:,2], -1, axis=1) - dFE[:,:,:,2] )*(i3/dg[1])  )
    return dFE

def getPoss( nsamp, extent ):
    nsamp=100
    extent = [0.0,4.2,0.0,4.2]
    xs,ys = np.meshgrid( np.linspace(extent[0],extent[1],nsamp), np.linspace(extent[2],extent[3],nsamp) )
    zs    =  xs*0.0
    ps = np.vstack([xs.flatten(), ys.flatten(), zs.flatten() ]).T.copy()    #;print("ps.shape: ",ps.shape)
    return ps

# # ----- Test 1: sample_SplineHermite 1D 
# dx    =  1.5
# x0    = -0.1 
# Eps   = np.array( [1.0, 0.0,-1.0,-0.5,-0.4,+0.1] )
# xp    = (np.array(range(len(Eps)))-1)*dx + x0
# xs    = np.linspace(0.0,4.2,100)
# Es,Fs = mmff.sample_SplineHermite( xs, Eps, x0=x0, dx=dx )  # ;print("Fs",Fs)
# plt.figure(); 
# plt.plot(xp, Eps, 'o-k', lw=0.2, label="Eps"); 
# plt.plot(xs, Es,'.-',    label="E"); 
# plt.plot(xs, Fs,'-',     label="F_ana");  
# plt.plot(xs[1:-1],numDeriv(xs,Es),':', label="F_num"); 
# plt.grid(); plt.legend()


# ----- Test 2: sample_SplineHermite 2D

# Eps   = np.array([ 
#     [1.0, 0.0,-1.0,-0.5,-0.4,+0.1], 
#     [1.0, 0.0,-1.0,-0.5,-0.4,+0.1],
#     [1.0, 0.0,-2.0,-1.5,-0.4,+0.1], 
#     [1.0, 0.0,-3.0,-2.5,-1.4,+0.1], 
#     [1.0, 0.0,-1.0,-0.5,-0.4,+0.1],
#     [1.0, 0.0,-1.0,-0.5,-0.4,+0.1],
# ])


# nsamp=100
# xs = np.linspace(0.0,4.2,nsamp)
# ps = np.zeros((nsamp,2))

# # ----- Test 2.1: sample_SplineHermite2D --- line along x
# ps[:,1] = 1.0
# ps[:,0] = xs
# FEs = mmff.sample_SplineHermite2D( ps, Eps, g0=[-0.1,-0.1], dg=[1.5,1.5] )  # ;print("Fs",Fs)
# plt.figure();
# plt.plot( xs, FEs[:,2], '.-', label="E" )
# plt.plot( xs, FEs[:,0], '.-', label="Fx" )
# plt.plot( xs, FEs[:,1], '.-', label="Fy" )
# plt.plot( xs[1:-1],numDeriv(xs,FEs[:,2]), ':',  label="Fx_num" )
# plt.grid(); plt.legend()

# # ----- Test 2.1: sample_SplineHermite2D --- line along y
# ps[:,1] = xs
# ps[:,0] = 1.0
# FEs = mmff.sample_SplineHermite2D( ps, Eps, g0=[-0.1,-0.1], dg=[1.5,1.5] )  # ;print("Fs",Fs)
# plt.figure();
# plt.plot( xs, FEs[:,2], '.-', label="E" )
# plt.plot( xs, FEs[:,0], '.-', label="Fx" )
# plt.plot( xs, FEs[:,1], '.-', label="Fy" )
# plt.plot( xs[1:-1],numDeriv(xs,FEs[:,2]), ':',  label="Fy_num" )
# plt.grid(); plt.legend()

# # ----- Test 2.3: sample_SplineHermite2D ---   imshow ( 2D color map )
# #nsamp=200
# nsamp=100
# extent = [0.0,4.2,0.0,4.2]
# xs,ys = np.meshgrid( np.linspace(extent[0],extent[1],nsamp), np.linspace(extent[2],extent[3],nsamp) )
# ps = np.vstack([xs.flatten(), ys.flatten()]).T.copy()
# #print("ps",ps)
# print("ps.shape: ",ps.shape)

# FEs = mmff.sample_SplineHermite2D( ps, Eps, g0=[-0.1,-0.1], dg=[1.5,1.5] )  # ;print("Fs",Fs)
# FEs = FEs.reshape(nsamp,nsamp,3)

# cmap='plasma'
# plt.figure(figsize=(15,5));
# plt.subplot(1,3,1); plt.imshow( FEs[:,:,2], extent=extent, origin='lower', cmap=cmap ); plt.colorbar(); plt.title("E")
# plt.subplot(1,3,2); plt.imshow( FEs[:,:,0], extent=extent, origin='lower', cmap=cmap ); plt.colorbar(); plt.title("Fx")
# plt.subplot(1,3,3); plt.imshow( FEs[:,:,1], extent=extent, origin='lower', cmap=cmap ); plt.colorbar(); plt.title("Fy")

# ----- Test 3: sample_SplineHermite 3D

ng = [6,6,6]
g0=np.array( [0.1, 0.1, 0.1 ] )
#dg=np.array( [ 1.5,  1.5, 1.5 ] )
dg=np.array( [ 0.2,  0.2, 0.2 ] )
atoms = np.array([
 [0.0, 0.0, 0.0,  1.0,   1.0],
 #[5.0, 1.5, 2.0,  0.5,  -1.0],
 #[2.5, 3.5, 4.0,  0.8,   0.5],
])

g0 = 3.0-0.3
dg = 0.1
FEg, xg = makeGrid_deriv_dir( atoms, 60, g0=(0.0,0.0,g0), dg=(0.0,0.0,dg) )
Eg = FEg[:,0].copy()
plt.plot( xg, FEg[:,0], 'ok', label="Eg   " )
plt.plot( xg, FEg[:,1], 'or',label="Fg_an" )
plt.plot( xg[1:-1], numDeriv( xg, FEg[:,0] ), '+g', label="Fg_num" )
#plt.grid()
#plt.legend()

# xs  = np.linspace( g0, 6.0, 1000 )
# FEs =  mmff.sample_SplineHermite1D_deriv( xs, FEg, g0, dg )
# plt.plot( xs, FEs[:,0], '-k', label="Es_deriv    ", lw=3 )
# plt.plot( xs, FEs[:,1], '-r', label="Fs_deriv_an ", lw=3 )
# plt.plot( xs[1:-1], numDeriv( xs, FEs[:,0] ), ':g', lw=3, label="Fs_deriv_num" )
# # plt.grid()
# # plt.legend()


# xs  = np.linspace( g0, 6.0, 1000 )
# FEs =  mmff.sample_SplineHermite( xs, Eg, g0, dg )
# plt.plot( xs, FEs[:,0], '--', label="Es_findif    ", c='gray' )
# plt.plot( xs, FEs[:,1], '--', label="Fs_findif_an ", c='orange')
# #plt.plot( xs[1:-1], numDeriv( xs, FEs[:,0] ), ':', label="Fs_findif_num" )

#FEg[:,1]*=-1
#Gs, Ws = mmff.fit_Bspline( FEg[:,0].copy(), dt=0.4, nmaxiter=1000, Ftol=1e-5 )

plt.figure()
Eg = FEg[:,0].copy()
Gs = FEg[:,0].copy()
#fmax=FEg[:,1].max()
Emin=FEg[:,0].min()
plt.plot( FEg[:,0], "ok", label="E_ref" )
plt.plot( FEg[:,1], "or", label="F_ref" )

Gs, Ws =  mmff.fitEF_Bspline( dg, FEg, Gs=Gs, nmaxiter=1, dt=1.0 )
plt.plot( Ws[:,0], "-", label=("F_fit  " ),  c='r')
plt.plot( (Ws[:,0]-FEg[:,1])*-1, "-", label=("Err W.x" ),  c='g' )
plt.plot( Ws[:,1], "-", label=("dErr/dF W.y" ), c='m' )
plt.legend()

colors = [  'k', 'r', 'g', 'b', 'm' ]
for i in range(1):
    Gs, Ws =  mmff.fitEF_Bspline( dg, FEg, Gs=Gs, nmaxiter=1, dt=1.0 )
    c=colors[i]
    #plt.plot( Ws[:,0], "-", label=("E_fit[%i]" %i ) )
    #plt.plot( Ws[:,1], "-", label=("E_fit[%i]" %i ) )
    plt.plot( Ws[:,0], "-", label=("F_fit  [%i]" %i ),  c=c)
    plt.plot( Ws[:,0]-FEg[:,1], "--", label=("Err   [%i]" %i ),  c=c)
    plt.plot( Ws[:,1], ":", label=("dErr/dF[%i]" %i ), c=c)
plt.legend()

plt.grid()
#plt.ylim(-fmax, fmax*1.5 )
plt.ylim(Emin*1.2, -Emin )

#plt.figure()
# xs  = np.linspace( g0, 6.0, 1000 )
# FEs =  mmff.sample_Bspline( xs, Gs, g0, dg )
# #FEs =  mmff.sample_Bspline( xs, Eg, g0, dg )
# plt.plot( xs, FEs[:,0], '-', lw=3, label="Es_findif    ", c='gray' )
# plt.plot( xs, FEs[:,1], '-', lw=3, label="Fs_findif_an ", c='orange')
# #plt.plot( xs[1:-1], numDeriv( xs, FEs[:,0] ), ':', label="Fs_findif_num" )
# plt.xlim(2,5)
# plt.ylim(-2,1)
# plt.grid()
# plt.legend()





# Eg  = makeGrid( atoms, ng, g0, dg )
# FE  = makeGrid_deriv( atoms, ng, g0, dg )
# dFE = make2dDeriv( FE, dg )



# ----- Test 3: sample_SplineHermite3D_deriv   1D-plot

# nsamp = 100
# xs    = np.linspace(0.1,0.5,nsamp)
# ps    = np.zeros((nsamp,3))
# ps[:,0],ps[:,1],ps[:,2] = xs,0.2,0.2
# FEs = mmff.sample_SplineHermite3D_deriv( ps, FE, dFE, g0=g0, dg=dg )
# #FEs = mmff.sample_SplineHermite3D_deriv( ps, FE, dFE*0, g0=g0, dg=dg )
# #FEs = mmff.sample_SplineHermite3D( ps, Eg, g0=g0, dg=dg )
# plt.figure();
# plt.plot( xs, FEs[:,3], '.-', label="E"  )
# plt.plot( xs, FEs[:,0], '.-', label="Fx" )
# # plt.plot( xs, FEs[:,1], '.-', label="Fy" )
# # plt.plot( xs, FEs[:,2], '.-', label="Fz" )
# plt.plot( xs[1:-1],numDeriv(xs,FEs[:,3]), ':',lw=2, label="Fx_num" )
# plt.grid(); plt.legend()
# FEs = mmff.sample_SplineHermite3D_deriv( ps, FE, dFE*0, g0=g0, dg=dg )
# plt.figure();
# plt.plot( xs, FEs[:,3], '.-', label="E"  )
# plt.plot( xs, FEs[:,0], '.-', label="Fx" )
# # plt.plot( xs, FEs[:,1], '.-', label="Fy" )
# # plt.plot( xs, FEs[:,2], '.-', label="Fz" )
# plt.plot( xs[1:-1],numDeriv(xs,FEs[:,3]), ':',lw=2, label="Fx_num" )
# plt.grid(); plt.legend()



# print(  "FE.shape, dFE.shape ",  FE.shape, dFE.shape )
# nsamp=100
# extent = [0.0,4.2,0.0,4.2]
# ps = getPoss( nsamp, extent )
# #FEs = mmff.sample_SplineHermite3D( ps, Eg, g0=g0, dg=dg )
# FEs = mmff.sample_SplineHermite3D_deriv( ps, FE, dFE, g0=g0, dg=dg )
# FEs = FEs.reshape(nsamp,nsamp,4)
# cmap='plasma'
# plt.figure(figsize=(15,5));
# plt.subplot(1,4,1); plt.imshow( FEs[:,:,3], extent=extent, origin='lower', cmap=cmap ); plt.colorbar(); plt.title("E" )
# plt.subplot(1,4,2); plt.imshow( FEs[:,:,0], extent=extent, origin='lower', cmap=cmap ); plt.colorbar(); plt.title("Fx")
# plt.subplot(1,4,3); plt.imshow( FEs[:,:,1], extent=extent, origin='lower', cmap=cmap ); plt.colorbar(); plt.title("Fy")
# plt.subplot(1,4,4); plt.imshow( FEs[:,:,2], extent=extent, origin='lower', cmap=cmap ); plt.colorbar(); plt.title("Fz")












# # ----- Test 3.1: sample_SplineHermite2D --- line along x
# nsamp = 100
# xs    = np.linspace(0.0,4.2,nsamp)
# ps    = np.zeros((nsamp,3))
# ps[:,0],ps[:,1],ps[:,2] = xs,1.0,1.0
# FEs = mmff.sample_SplineHermite3D( ps, Eg, g0=g0, dg=dg )
# plt.figure();
# plt.plot( xs, FEs[:,3], '.-', label="E"  )
# plt.plot( xs, FEs[:,0], '.-', label="Fx" )
# plt.plot( xs, FEs[:,1], '.-', label="Fy" )
# plt.plot( xs, FEs[:,2], '.-', label="Fz" )
# plt.plot( xs[1:-1],numDeriv(xs,FEs[:,3]), ':',lw=2, label="Fx_num" )
# plt.grid(); plt.legend()

# # ----- Test 3.2: sample_SplineHermite2D --- line along y
# ps[:,0],ps[:,1],ps[:,2] = 1.0,xs,1.0
# FEs = mmff.sample_SplineHermite3D( ps, Eg, g0=g0, dg=dg )
# plt.figure();
# plt.plot( xs, FEs[:,3], '.-', label="E"  )
# plt.plot( xs, FEs[:,0], '.-', label="Fx" )
# plt.plot( xs, FEs[:,1], '.-', label="Fy" )
# plt.plot( xs, FEs[:,2], '.-', label="Fz" )
# plt.plot( xs[1:-1],numDeriv(xs,FEs[:,3]), ':', lw=3, zorder=3,  label="Fy_num" )
# plt.grid(); plt.legend()

# ----- Test 3.3: sample_SplineHermite3D ---  imshow ( 2D color map )

'''
nsamp=100
extent = [0.0,4.2,0.0,4.2]
xs,ys = np.meshgrid( np.linspace(extent[0],extent[1],nsamp), np.linspace(extent[2],extent[3],nsamp) )
zs    =  xs*0.5
ps = np.vstack([xs.flatten(), ys.flatten(), zs.flatten() ]).T.copy()    #;print("ps.shape: ",ps.shape)

FEs = mmff.sample_SplineHermite3D( ps, Eg, g0=g0, dg=dg )
#FEs = mmff.sample_SplineHermite3D_f( ps.astype(np.float32), Eg.astype(np.float32), g0=g0, dg=dg )
FEs = FEs.reshape(nsamp,nsamp,4)

cmap='plasma'
plt.figure(figsize=(15,5));
plt.subplot(1,4,1); plt.imshow( FEs[:,:,3], extent=extent, origin='lower', cmap=cmap ); plt.colorbar(); plt.title("E" )
plt.subplot(1,4,2); plt.imshow( FEs[:,:,0], extent=extent, origin='lower', cmap=cmap ); plt.colorbar(); plt.title("Fx")
plt.subplot(1,4,3); plt.imshow( FEs[:,:,1], extent=extent, origin='lower', cmap=cmap ); plt.colorbar(); plt.title("Fy")
plt.subplot(1,4,4); plt.imshow( FEs[:,:,2], extent=extent, origin='lower', cmap=cmap ); plt.colorbar(); plt.title("Fz")
'''

plt.show()