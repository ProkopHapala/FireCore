import numpy as np
import os
#import matplotlib.pyplot as plt
#import matplotlib.patches as patches

def getPLQH( R0, E0, a, Q, H ):
    e  = np.exp(a*R0);
    cL = e*E0;
    cP = e*cL;
    cH = e*e*H;
    return np.array([ cP, cL, Q, cH ])

def make_sample_points( p0, t0=0.0, tmax=10.0, dsamp=0.02, iax=2 ):
    ts = np.arange( t0, tmax, dsamp)
    ps = np.zeros( (len(ts), 3) )
    ps[:,0] = p0[0]
    ps[:,1] = p0[1]
    ps[:,2] = p0[2]
    ps[:,iax] = ts
    return ps, ts

def make_sample_points_f4( p0, t0=0.0, tmax=10.0, dsamp=0.02, iax=2 ):
    ts = np.arange( t0, tmax, dsamp)
    ps = np.zeros( (len(ts), 4), dtype=np.float32 )
    ps[:,0] = p0[0]
    ps[:,1] = p0[1]
    ps[:,2] = p0[2]
    ps[:,3] = 0.0
    ps[:,iax] = ts
    return ps, ts

def points_to_ocl(ps_):
    ps = np.zeros( ( len(ps_), 4), dtype=np.float32 )
    ps[:,0] = ps_[:,0]
    ps[:,1] = ps_[:,1]
    ps[:,2] = ps_[:,2]
    ps[:,3] = 0.0
    return ps


def load_potential_comb( path ):
    VPaul = np.load( os.path.join( path, "debug_BsplinePaul_pbc.npy" ) )
    VLond = np.load( os.path.join( path, "debug_BsplineLond_pbc.npy" ) )
    VCoul = np.load( os.path.join( path, "debug_BsplineCoul_pbc.npy" ) )
    #print("VCoul.shape = ", VCoul.shape)
    #print("VLond.shape = ", VLond.shape)
    #print("VPaul.shape = ", VPaul.shape)
    VPLQ = np.zeros(VCoul.shape+(4,), dtype=np.float32)
    VPLQ[:,:,:,0] = VPaul
    VPLQ[:,:,:,1] = VLond
    VPLQ[:,:,:,2] = VCoul
    VPLQ[:,:,:,3] = 0.0
    print("VPLQ.shape(BEFORE) = ", VPLQ.shape)
    VPLQ = VPLQ.transpose( (2,1,0,3) ).copy()
    sh = VPLQ.shape; print("sh = ", sh)
    print("VPLQ.shape(AFTER) = ", VPLQ.shape)
    return VPLQ