#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import argparse

# ======== Functions

def getLJ( r, R0=3.0, E0=1.0 ):
    u = R0/r
    u6 = u**6
    return E0*u6*(u6-2.0)

def getMorse( r, R0=3.0, E0=1.0, k=1.6 ):
    e = np.exp( -k*(r-R0) )
    return E0*e*(e-2.0)

def getSR( r,  E0=1.0, k=1.0 ):
    u = R0/r
    return E0*np.exp(-k*r)

def getSR2( r,  E0=1.0, w=1.0 ):
    u = r/w
    return E0*np.exp(-r*r)

def makeCircle( R=3.0 ):
    angs = np.linspace(0.0,2.0*np.pi,100)
    ps = np.zeros( (len(angs),2) )
    ps[:,0] = R*np.cos(angs)
    ps[:,1] = R*np.sin(angs)
    return ps

def plot( Es, Emax, Rcirc=None, apos=None, ps_circ=None ):
    plt.imshow( Es, origin='lower', vmin=-Emax, vmax=Emax, cmap='bwr', extent=[xs[0],xs[-1],ys[0],ys[-1]] )
    plt.colorbar()
    if apos is not None:
        plt.plot( apos[:,0], apos[:,1], '+k'  )
    if Rcirc is not None:
        ps_circ = makeCircle(Rcirc)    
    if ps_circ is not None:
        plt.plot( ps_circ[:,0], ps_circ[:,1], '--', lw=0.5, c='k' )
    plt.xlim( xs[0], xs[-1] )
    plt.ylim( ys[0], ys[-1] )
    plt.xlabel("x [A]")
    plt.ylabel("y [A]")
    
def plot3( Euncorr, Ecorr, Emax1, Emax2, Rcirc=3.0 ):
    ps_circ = makeCircle(Rcirc)
    apos = np.array( [[0.0,0.0],[Le,0.0]] )
    plt.figure(figsize=(15,5))
    plt.subplot(1,3,1); plot( Euncorr , Emax1, apos=apos, ps_circ=ps_circ )
    plt.subplot(1,3,2); plot( Ecorr   , Emax1, apos=apos, ps_circ=ps_circ )
    plt.subplot(1,3,3); plot( Ecorr   , Emax2, apos=apos, ps_circ=ps_circ )

# ======== Main

# make agrparse for  Morse (short and long agrument name)

parser = argparse.ArgumentParser()
parser.add_argument( "-L", type=float, default=1.0  )
parser.add_argument( "-R", type=float, default=3.0  )
parser.add_argument( "-E", type=float, default=0.01 )
parser.add_argument( "-k", type=float, default=1.6  )
parser.add_argument( "-M", "--bMorse", action='store_true' ) # if -M is given, bMorse=True
args = parser.parse_args()

# make arrays
xs = np.linspace( 0.0,8.0,300)
ys = np.linspace(-3.0,3.0,300)
Xs,Ys = np.meshgrid(xs,ys)
Rs = np.sqrt(Xs**2 + Ys**2)

E0  = args.E
R0  = args.R
Le  = args.L

bMorse = args.bMorse
if bMorse:
    EHb = 1.0
    E_uncorr = getMorse(Rs, R0, E0)
else:
    EHb = 10.0
    E_uncorr = getLJ(Rs, R0, E0)
Rse    = np.sqrt((Xs-Le)**2 + Ys**2)
E_corr = E_uncorr + getSR2(Rse, -EHb, w=0.1 )


plot3( E_uncorr, E_corr, E0, EHb*0.2, Rcirc=R0 )


plt.grid()

plt.show()