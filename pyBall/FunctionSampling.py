
import numpy as np



# ========= Morse Potential =========

def getMorse( r, R0, E0, a=1.7 ):
    #r = np.sqrt(x**2 + y**2 + z**2)
    e = np.exp( -a*(r-R0) )
    E  =     E0*(e*e - 2*e )
    fr = 2*a*E0*(e*e -   e ) 
    return E,fr*r

def getMorseSplit( r, R0, E0, a=1.7 ):
    #r = np.sqrt(x**2 + y**2 + z**2)
    e = np.exp( -a*(r-R0) )
    Ep  =  E0*e*e
    El  =  E0*-2*e
    Fp  =  E0*a*2*e*e
    Fl  =  E0*a*-2*e 
    return Ep,El,Fp,Fl

# ========= Lennard-Jones Potential =========

def getLJ( r, R0, E0, a=0 ):
    #r = np.sqrt(x**2 + y**2 + z**2)
    E  = E0*    ( (R0/r)**12 - 2*(R0/r)**6 )
    fr = E0*-12*( (R0/r)**12 -   (R0/r)**6 )/(r*r)
    return E,fr*r

def getLJSplit( r, R0, E0, a=0 ):
    e = (R0/r)**6
    Ep = E0*    e*e
    El = E0* -2*e
    Fp = E0*-12*e*e/r
    Fl = E0* 12*e  /r
    return Ep,El,Fp,Fl

def getLJ_atoms( apos, REs, Xs,Ys,Zs ):
    ng  = len(Xs)
    #Es = np.zeros( ng )
    #Fx = np.zeros( ng )
    #Fy = np.zeros( ng )
    #Fz = np.zeros( ng )
    Es = Xs*0.0
    Fx = Xs*0.0
    Fy = Xs*0.0
    Fz = Xs*0.0
    for i,p in enumerate(apos):
        print( p )
        dx = Xs-p[0]
        dy = Ys-p[1]
        dz = Zs-p[2]
        r = np.sqrt( dx**2 + dy**2 + dz**2 )
        R0,E0 = REs[i]
        E, fr = getLJ( r, R0, E0 )
        Es += E
        Fx += fr*dx/r
        Fy += fr*dy/r
        Fz += fr*dz/r
    return Es, Fx,Fy,Fz

# ========= Cosine Potential =========

def getCos( xs, freq ):
    E =       np.cos( freq*xs )
    F = -freq*np.sin( freq*xs )
    return E,F

def getCos2D( Xs, Ys, freq=(np.pi,np.pi) ):
    E =           np.cos( freq[0]*Xs )*np.cos( freq[1]*Ys )
    Fx = -freq[0]*np.sin( freq[0]*Xs )*np.cos( freq[1]*Ys )
    Fy = -freq[1]*np.cos( freq[0]*Xs )*np.sin( freq[1]*Ys )
    return E,Fx,Fy

def getCos3D( Xs, Ys, Zs, freq=(np.pi,np.pi,np.pi) ):
    fx,fy,fz = freq
    cx = np.cos( fx*Xs );  sx = np.sin( fx*Xs )
    cy = np.cos( fy*Ys );  sy = np.sin( fy*Ys )
    cz = np.cos( fz*Zs );  sz = np.sin( fz*Zs )
    E =           cx*cy*cz
    Fx = -fx*sx*cy*cz
    Fy = -fy*cx*sy*cz
    Fz = -fz*cx*cy*sz
    return E,Fx,Fy,Fz


# ========= Grid sampling generator in 2D and 3D =========

def make2Dsampling(  g0=(-5.0,2.0), gmax=(5.0,10.0), dg=(0.1,0.1) ):
    xs  = np.arange(g0[0], gmax[0], dg[0])
    ys  = np.arange(g0[1], gmax[1], dg[0])
    Xs,Ys = np.meshgrid(xs,ys)
    return Xs,Ys

def pack_ps2D( Xs, Ys):
    ps = np.zeros( ( len(Xs.flat), 2) )
    ps[:,0] = Xs.flat
    ps[:,1] = Ys.flat
    return ps

def make2Dsampling_ps(  g0=(-5.0,2.0), gmax=(5.0,10.0), dg=(0.1,0.1) ):
    Xs,Ys = make2Dsampling(  g0=g0, gmax=gmax, dg=dg )
    sh = Xs.shape
    ps = pack_ps2D( Xs, Ys)
    return ps, sh

def make3Dsampling(  g0=(-5.0,2.0), gmax=(5.0,10.0), dg=(0.1,0.1) ):
    xs  = np.arange(g0[0], gmax[0], dg[0])
    ys  = np.arange(g0[1], gmax[1], dg[0])
    zs  = np.arange(g0[2], gmax[2], dg[2])
    Xs,Ys,Zs = np.meshgrid(xs,ys,zs)
    return Xs,Ys,Zs

def pack_ps3D( Xs, Ys, Zs):
    ps = np.zeros( ( len(Xs.flat), 3) )
    ps[:,0] = Xs.flat
    ps[:,1] = Ys.flat
    ps[:,2] = Zs.flat
    return ps