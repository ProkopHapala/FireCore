import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF as mmff

def numDeriv( xs, Es): 
    dx = xs[1]-xs[0]
    Fs = (Es[2:]-Es[:-2])/(2*dx)
    return Fs

def makeGrid( ng, g0, dg ):
    xs = np.arange(ng[0])*dg[0] + g0[0]   #;print("xs: ",xs)
    ys = np.arange(ng[1])*dg[1] + g0[1] 
    zs = np.arange(ng[2])*dg[2] + g0[2] 
    Xs,Ys,Zs = np.meshgrid( xs,ys,zs )
    atoms = np.array([
    [2.0, 5.0, 3.0,  1.0,   1.0],
    [5.0, 1.5, 2.0,  0.5,  -1.0],
    [2.5, 3.5, 4.0,  0.8,   0.5],
    ])
    Eg = np.zeros( ng )
    for i,a in enumerate(atoms):
        X = Xs-a[0]
        Y = Ys-a[1]
        Z = Zs-a[2]
        Eg += a[4]*( 1/( X**2 + Y**2 + Z**2 + a[3]**2 ) )
    return Eg




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
g0=np.array( [-0.1, -0.1, 0.1 ] )
dg=np.array( [ 1.5,  1.5, 1.5 ] )
atoms = np.array([
 [2.0, 5.0, 3.0,  1.0,   1.0],
 [5.0, 1.5, 2.0,  0.5,  -1.0],
 [2.5, 3.5, 4.0,  0.8,   0.5],
])

Eg = makeGrid( ng, g0, dg )

nsamp = 100
xs    = np.linspace(0.0,4.2,nsamp)
ps    = np.zeros((nsamp,3))

# # ----- Test 3.1: sample_SplineHermite2D --- line along x
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

plt.show()