#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def plot_fname(fname='gridFF_EFprofile.log', iax=2):
    dat = np.genfromtxt(fname, skip_header=1 )
    Emin=dat[:,7].min()
    Fmin=dat[:,6].min()
    plt.subplot(2,1,1); plt.plot( dat[:,iax+1], dat[:,    7], lw=0.5, label="E "+fname  ); plt.ylim(Emin, -Emin); plt.legend(); plt.grid()
    plt.subplot(2,1,2); plt.plot( dat[:,iax+1], dat[:,iax+4], lw=0.5, label="F_z"+fname ); plt.ylim(Fmin, -Fmin); plt.legend(); plt.grid()
    return Emin, Fmin

plt.figure(figsize=(5,10))
plot_fname( 'gridFF_EFprofile_mod1.log' );
plot_fname( 'gridFF_EFprofile_mod2.log' );
plot_fname( 'gridFF_EFprofile_mod6.log' );

plt.figure(figsize=(5,10))
plot_fname( 'gridFF_EFprofile_mod1_x.log', iax=0 );
plot_fname( 'gridFF_EFprofile_mod2_x.log', iax=0 );
plot_fname( 'gridFF_EFprofile_mod6_x.log', iax=0 );

plt.figure(figsize=(5,10))
plot_fname( 'gridFF_EFprofile_mod1_y.log', iax=1 );
plot_fname( 'gridFF_EFprofile_mod2_y.log', iax=1 );
plot_fname( 'gridFF_EFprofile_mod6_y.log', iax=1 );


plt.show()