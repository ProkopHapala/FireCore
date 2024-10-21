import numpy as np
import matplotlib.pyplot as plt

from . import utils as ut
from .. import atomicUtils as au
from .. import MMFF as mmff
from .. import FunctionSampling as fu

# =============  Functions

def test_gridFF_npy_Ewald( name="data/NaCl_1x1_L2",  ps_xy=[(0.0,0.0),(0.0,0.5),(0.5,0.0),(0.5,0.5)], mode=6, title="", bSaveFig=False, Vname= "VCoul" ):
    print( "py======= test_gridFF_npy() START" );
    
    # ----- 1D Ewald Debug
    fname="debug_"+Vname+".npy"
    fname="debug_"+Vname+".npy"
    mmff.makeGridFF( name=name, mode=mode, bSymmetrize=False, bFit=False )
    # ----- 2D "debug_Ewald_dens.npy"
    dens = np.load(name+"/"+"debug_Ewald_dens.npy")
    V    = np.load(name+"/"+"debug_Ewald_V.npy")
    Qtot=dens.sum();
    Qabs=np.abs(dens).sum();
    print( "Qtot= ",Qtot," Qabs=", Qabs );
    plt.figure(figsize=(15,10))
    plt.subplot(2,2,1); plt.imshow( dens[:,:,0].transpose() ); plt.title("dens(z,y)");
    plt.subplot(2,2,3); plt.imshow( V   [:,:,0].transpose() ); plt.title("V(z,y)");
    plt.subplot(2,2,2); plt.imshow( dens[10,:,:]            ); plt.title("dens(x,y)");
    plt.subplot(2,2,4); plt.imshow( V   [10,:,:]            ); plt.title("V(x,y)");
    plt.grid()
    plt.legend()
    plt.title(title+" "+Vname)
    
    if bSaveFig:
        plt.savefig(name+"/"+fname+".png")
    print( "py======= test_gridFF() DONE" );

#def test_gridFF_npy( name="data/NaCl_1x1_L2",  ps_xy=[(0.0,0.0),(0.0,0.5),(0.5,0.0),(0.5,0.5)], mode=6, title="", bSaveFig=False, Vname= "Coul" ):
def test_gridFF_npy( name="data/NaCl_1x1_L2",  ps_xy=[(0.0,0.0)], mode=6, title="", bSaveFig=False, Vname= "Coul", Bsc=2./3. ):
    print( "py======= test_gridFF_npy() START" );
    # ----- 1D "debug_"+Vname+".npy"
    plt.figure()
    fname="debug_V"+Vname+".npy"

    dat  = np.load(name+"/"+fname)
    datB = np.load(name+"/"+"debug_Bspline"+Vname+"_pbc.npy")
    nz,ny,nx = dat.shape
    for i,p in enumerate(ps_xy):
        plt.plot( dat [:,int(p[1]*ny),int(p[0]*nx)],       label='z-cut '+str(p) );
        plt.plot( datB[:,int(p[1]*ny),int(p[0]*nx)]*Bsc, '--', label='z-cut '+str(p) );
    plt.grid()
    plt.legend()
    plt.title(title+" "+Vname)
    if bSaveFig:
        plt.savefig(name+"/"+fname+".png")
    #plt.show()
    print( "py======= test_gridFF() DONE" );


#def test_gridFF_npy( name="data/NaCl_1x1_L2",  ps_xy=[(0.0,0.0),(0.0,0.5),(0.5,0.0),(0.5,0.5)], mode=6, title="", bSaveFig=False, Vname= "Coul" ):
def test_gridFF_npy_lat( name="data/NaCl_1x1_L2", ps_zy=[(0.0,0.0)], mode=6, title="", bSaveFig=False, Vname= "Coul", Bsc=2./3. ):
    print( "py======= test_gridFF_npy() START" );
    # ----- 1D "debug_"+Vname+".npy"
    plt.figure()
    fname="debug_V"+Vname+".npy"

    dat  = np.load(name+"/"+fname)
    datB = np.load(name+"/"+"debug_Bspline"+Vname+"_pbc.npy")
    nz,ny,nx = dat.shape
    for i,p in enumerate(ps_zy):
        plt.plot( dat [int(p[1]*ny),int(p[0]*nx),:],           label='Vref z-cut '+str(p) );
        plt.plot( datB[int(p[1]*ny),int(p[0]*nx),:]*Bsc, '--', label='Bspl z-cut '+str(p) );
    plt.grid()
    plt.legend()
    plt.title(title+" "+Vname)
    if bSaveFig:
        plt.savefig(name+"/"+fname+".png")
    #plt.show()
    print( "py======= test_gridFF() DONE" );

def test_gridFF( name="data/xyz/NaCl_1x1_L2", mode=6, dsamp=0.02, tmin=0.0,tmax=10.0, p0=[0.0,0.0,2.0], R0=3.5, E0=0.1, a=1.6, Q=0.4, H=0.0, scErr=100.0, Emax=None, Fmax=None, maxSc=5.0, title=None, bSaveFig=True, bRefine=True, nPBC=None ):
    print( "py======= test_gridFF() START" );
    #print( "test_gridFF() START" )
    #mode = 4
    #mode = 1
    mmff.makeGridFF( name=name, mode=mode, bRefine=bRefine )
    #plotGridFF_1D( EFg, ix=20,iy=20 )
    PLQH = ut.getPLQH( R0, E0, a, Q, H )

    ps,ts = ut.make_sample_points( p0, t0=tmin, tmax=tmax, dsamp=dsamp, iax=2 )
    
    #FF_ref = mmff.evalGridFFAtPoints( ps, PLQH=PLQH, bSplit=False, nPBC=nPBC )
    FF_ref = mmff.evalGridFFAtPoints( ps, PLQH=PLQH, bSplit=True, nPBC=nPBC  )
    
    #FFout  = mmff.sample_SplineHermite3D_comb3( ps, EFg, g0=[0.0,0.0,0.0], dg=[0.1,0.1,0.1], fes=None, Cs=PLQH )
    # Es,Fs = sampleSurf( name, zs, Es=None, fs=None, kind=1, atyp=0, Q=0.0, K=-1.0, Rdamp=1.0, pos0=(0.,0.,0.), bSave=False )

    ps_ = ps.copy();
    
    FFout = mmff.sampleSurf_new( ps_, PLQH, mode=mode, Rdamp=1.0 )  # * 15.1
    
    #Emin = FFout[:,3].min();
    #Fmin = FFout[:,2].min();

    if Emax is None: Emax = -FF_ref[:,3].min()*maxSc; 
    if Fmax is None: Fmax = -FF_ref[:,2].min()*maxSc;

    plt.figure(figsize=(5,10))
    plt.subplot(2,1,1);
    plt.plot( ts, FFout[:,3],  "-g", lw=0.5, label="Etot_fit" )
    plt.plot( ts, FF_ref[:,3], ":k", lw=2.0, label="Etot_ref" )
    plt.plot( ts, (FFout[:,3]-FF_ref[:,3])*scErr, "-r", lw=0.5, label=("Etot_err*%.2f" %scErr) )
    #plt.plot( ts, FFout[:,0], "-r", lw=0.5, label="Ftot_x" )
    #plt.plot( ts, FFout[:,1], "-g", lw=0.5, label="Ftot_y" )
    #plt.plot( ts, FFout[:,2], "-b", lw=0.5, label="Etot_z" )
    plt.axhline(0.0, c="k", ls='--', lw=0.5)
    plt.ylim( -Emax, Emax )
    plt.legend()
    plt.subplot(2,1,2);
    plt.plot( ts, FFout[:,2],  "-g", lw=0.5, label="Ftot_fit" )
    plt.plot( ts, FF_ref[:,2], ":k", lw=2.0, label="Ftot_ref" )
    plt.plot( ts, (FFout[:,2]-FF_ref[:,2])*scErr, "-r", lw=0.5, label=("Ftot_err*%.2f" %scErr) )
    plt.axhline(0.0, c="k", ls='--', lw=0.5)
    plt.ylim( -Fmax, Fmax )
    plt.legend()

    title_ = "p0="+str(p0)+"\nnPBC="+str(nPBC)
    if ( title is not None ): title_=title+"\n"+title_
    plt.suptitle( title_ )
    if ( bSaveFig ): plt.savefig( "test_gridFF_zcut.png", bbox_inches='tight' )
    
    #print( "ff.shape ", EFg.shape )
    #print( "test_gridFF() DONE" )
    print( "py======= test_gridFF() DONE" );
    #return EFg

def test_gridFF_lat( name="data/xyz/NaCl_1x1_L2", iax=0, tmin=0.0,tmax=10.0, p0=[1.05,1.05,2.0], mode=6, dsamp=0.02,  R0=3.5, E0=0.1, a=1.6, Q=0.4, H=0.0, scErr=100.0, title=None, Emax=None, Fmax=None, maxSc=5.0, bSaveFig=True, bRefine=True, nPBC=None ):
    print( "py======= test_gridFF_lat() START" );
    #print( "test_gridFF() START" )
    #mode = 4
    #mode = 1
    mmff.makeGridFF( name=name, mode=mode, bRefine=bRefine )

    #plotGridFF_1D( EFg, ix=20,iy=20 )
    PLQH = ut.getPLQH( R0, E0, a, Q, H )

    ps,ts = ut.make_sample_points( p0, t0=tmin, tmax=tmax, dsamp=dsamp, iax=iax )    

    FF_ref = mmff.evalGridFFAtPoints( ps, PLQH=PLQH, nPBC=nPBC )

    print( "ps ", ps)
    print( "FF_ref ", FF_ref)
    
    #FFout  = mmff.sample_SplineHermite3D_comb3( ps, EFg, g0=[0.0,0.0,0.0], dg=[0.1,0.1,0.1], fes=None, Cs=PLQH )
    # Es,Fs = sampleSurf( name, zs, Es=None, fs=None, kind=1, atyp=0, Q=0.0, K=-1.0, Rdamp=1.0, pos0=(0.,0.,0.), bSave=False )

    ps_ = ps.copy(); 

    FFout = mmff.sampleSurf_new( ps_, PLQH, mode=mode, Rdamp=1.0 )
    
    # Emin = FFout[:,3].min();
    # Fmin = FFout[:,2].min();

    #Emin = FF_ref[:,3].min(); 
    #Fmin = FF_ref[:,2].min()

    if Emax is None: Emax = -FF_ref[:,3].min()*maxSc; 
    if Fmax is None: Fmax = -FF_ref[:,2].min()*maxSc;

    plt.figure(figsize=(5,10))
    plt.subplot(2,1,1);
    plt.plot( ts, FFout[:,3], "-g", lw=0.5, label="Etot_fit" )
    plt.plot( ts, FF_ref[:,3], ":k", lw=2.0, label="Etot_ref" )
    plt.plot( ts, (FFout[:,3]-FF_ref[:,3])*scErr, "-r", lw=0.5, label=("Etot_err*%.2f" %scErr) )
    #plt.plot( zs, FFout[:,0], "-r", lw=0.5, label="Ftot_x" )
    #plt.plot( zs, FFout[:,1], "-g", lw=0.5, label="Ftot_y" )
    #plt.plot( zs, FFout[:,2], "-b", lw=0.5, label="Etot_z" )
    plt.axhline(0.0, c="k", ls='--', lw=0.5)
    plt.ylim( -Emax, Emax )
    plt.legend()
    plt.grid()
    plt.subplot(2,1,2);
    plt.plot( ts, FFout[:,iax],  "-g", lw=0.5, label="Ftot_fit" )
    plt.plot( ts, FF_ref[:,iax], ":k", lw=2.0, label="Ftot_ref" )
    plt.plot( ts, (FFout[:,iax]-FF_ref[:,iax])*scErr, "-r", lw=0.5, label=("Ftot_err*%.2f" %scErr) )
    plt.axhline(0.0, c="k", ls='--', lw=0.5)
    plt.ylim( -Fmax, Fmax )
    plt.legend()
    plt.grid()


    title_ = "p0="+str(p0)+"\nnPBC="+str(nPBC)
    if ( title is not None ): title_=title+"\n"+title_
    plt.suptitle( title_ )
    if ( bSaveFig ): plt.savefig( "test_gridFF_lat.png", bbox_inches='tight' )
    
    #print( "ff.shape ", EFg.shape )
    #print( "test_gridFF() DONE" )
    print( "py======= test_gridFF_lat() DONE" );
    #return EFg


def test_gridFF_2D( name="data/NaCl_1x1_L2", axs=(0,1), tmin=[0.0,0.0],tmax=[10.0,10.0], p0=[1.05,1.05,2.0], mode=6, dsamp=0.1,  R0=3.5, E0=0.1, a=1.6, Q=0.4, H=0.0, scErr=100.0, title=None, bSaveFig=True ):
    print( "py======= test_gridFF_lat() START" );
    #print( "test_gridFF() START" )
    #mode = 4
    #mode = 1
    mmff.makeGridFF( name=name, mode=mode )

    #plotGridFF_1D( EFg, ix=20,iy=20 )
    PLQH = ut.getPLQH( R0, E0, a, Q, H )

    xs = np.arange( tmin[0], tmax[0], dsamp)
    ys = np.arange( tmin[0], tmax[0], dsamp)
    Xs,Ys = np.meshgrid( xs, ys )
    ps = np.zeros( (len(ts), 3) )
    ps[:,0] = p0[0]
    ps[:,1] = p0[1]
    ps[:,2] = p0[2]
    ps[:,axs[0]] = Xs.flat
    ps[:,axs[1]] = Ys.flat

    FF_ref = mmff.evalGridFFAtPoints( ps, PLQH=PLQH )
    
    # ps_ = ps.copy(); 
    # ps_[:,2]+=-2.0+3.25;
    # ps_[:,0]+=2.0;
    # ps_[:,1]+=2.0;

    FFout = mmff.sampleSurf_new( ps_, PLQH, mode=mode, Rdamp=1.0 )
    
    Emin = FFout[:,3].min();
    Fmin = FFout[:,2].min();

    if ( title is not None ): plt.suptitle( title )
    if ( bSaveFig ): plt.savefig( "test_gridFF_2d.png", bbox_inches='tight' )
    
    #print( "ff.shape ", EFg.shape )
    #print( "test_gridFF() DONE" )
    print( "py======= test_gridFF_lat() DONE" );
    #return EFg

# ============= BODY

if __name__ == "__main__":
    # apos = [
    #     [-2.0,0.0,0.0],
    #     [ 2.0,0.0,0.0],
    # ]
    # REs=[
    #     [3.5,1.0],
    #     [3.5,1.0],
    # ]

    R0 = 3.5
    E0 = 1.0
    a  = 1.8

    #Q = 0.0
    Q = 0.4
    #p0 = [1.0,-5.05,2.0]
    #p0 = [0.0,0.0,2.0]
    p0 = [-2.0,-2.0,0.0]

    mmff.initParams()

    # test_gridFF    ( mode=6, title="Bspline_o3 \n(z-cut)" ,    Q=0.0, )
    # test_gridFF_lat( mode=6, title="Bspline_o3 \n(lat iax=0)", Q=0.0, )
    #test_gridFF_lat( mode=6, title="Bspline_o3 \n(lat iax=0)", Q=0.0, p0=p0, iax=0 )

    #test_gridFF    ( mode=1, title="tri-linar force \n(z-cut)"          )
    #test_gridFF_lat( mode=1, title="tri-Linear Force", Q=0.0, p0=p0, iax=0 )


    #Emax=0.01 
    #Fmax=0.01

    Emax=0.00001 
    Fmax=0.00001

    #test_gridFF    ( mode=6, title="Bspline_o3 \n(z-cut)" , p0=[0.0,0.0,2.0],  Q=0.4, E0=0, bRefine=False, nPBC=[50 ,50 ,0], Emax=Emax, Fmax=Fmax )
    #test_gridFF    ( mode=6, title="Bspline_o3 \n(z-cut)" , p0=[0.0,0.0,2.0],  Q=0.4, E0=0, bRefine=False, nPBC=[100,100,0], Emax=Emax, Fmax=Fmax )
    #test_gridFF    ( mode=6, title="Bspline_o3 \n(z-cut)" , p0=[0.0,0.0,2.0],  Q=0.4, E0=0, bRefine=False, nPBC=[150,150,0], Emax=Emax, Fmax=Fmax )

    #test_gridFF    ( mode=6, title="Bspline_o3 \n(z-cut)" , p0=[0.0,0.0,2.0],  Q=0.4, E0=0, bRefine=False, nPBC=[100,100,0], Emax=Emax, Fmax=Fmax )
    #test_gridFF    ( mode=6, title="Bspline_o3 \n(z-cut)" , p0=[0.0,0.0,2.0],  Q=0.4, E0=0, bRefine=False, nPBC=[300,300,0], Emax=Emax, Fmax=Fmax )
    test_gridFF    ( mode=6, title="Bspline_o3 \n(z-cut)" , p0=[0.0,0.0,2.0],  Q=0.4, E0=0, bRefine=False, nPBC=[400,400,0], Emax=Emax, Fmax=Fmax )
    test_gridFF    ( mode=6, title="Bspline_o3 \n(z-cut)" , p0=[2.0,2.0,2.0],  Q=0.4, E0=0, bRefine=False, nPBC=[400,400,0], Emax=Emax, Fmax=Fmax )

    Emax=0.01 
    Fmax=0.01

    test_gridFF_lat( mode=6, title="Bspline_o3 \n(lat iax=0)", p0=[0.0,0.0,2.0], nPBC=[400,400,0],  Q=0.4, E0=0.0, Emax=Emax, Fmax=Fmax  )
    test_gridFF_lat( mode=6, title="Bspline_o3 \n(lat iax=0)", p0=[2.0,2.0,2.0], nPBC=[400,400,0],  Q=0.4, E0=0.0, Emax=Emax, Fmax=Fmax  )

    Emax=0.1 
    Fmax=0.1

    # test_gridFF    ( mode=6, title="Bspline_o3 \n(z-cut)" , p0=[0.0,0.0,2.0],  Q=0.0, E0=0.1, bRefine=False, nPBC=[5,5,0], Emax=Emax, Fmax=Fmax )
    # test_gridFF    ( mode=6, title="Bspline_o3 \n(z-cut)" , p0=[2.0,2.0,2.0],  Q=0.0, E0=0.1, bRefine=False, nPBC=[5,5,0], Emax=Emax, Fmax=Fmax )
    # test_gridFF_lat( mode=6, title="Bspline_o3 \n(lat iax=0)", p0=[0.0,0.0,2.0], Q=0.0, E0=0.1 )
    # test_gridFF_lat( mode=6, title="Bspline_o3 \n(lat iax=0)", p0=[2.0,2.0,2.0], Q=0.0, E0=0.1 )

    #test_gridFF_npy( ps_xy=[(0.0,0.0),(0.0,0.5),(0.5,0.0),(0.5,0.5)], mode=6, title="" )

    #test_gridFF_npy( ps_xy=[(0.0,0.0)], mode=6, title="" )
    #test_gridFF_npy_lat( ps_zy=[(0.0,0.0)], mode=6, title="" )
    #test_gridFF_npy_lat( ps_zy=[(0.1,0.1)], mode=6, title="" )

    plt.show()