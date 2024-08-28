#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def getEF( fname, iax=2 ):
    dat = np.genfromtxt(fname, skip_header=1 )
    return dat[:,iax+1], dat[:,iax+4],  dat[:,7]

def plot_fname(fname='gridFF_EFprofile.log', iax=2, ls='-', c=None, lw=0.5, Eref=None, Fref=None, scErr=100.0 ):
    #dat = np.genfromtxt(fname, skip_header=1 )
    xs,F,E = getEF( fname, iax=iax )

    Emin=E.min()
    Fmin=F.min()

    plt.subplot(2,1,1); 
    plt.plot(                      xs,  E,             ls=ls,   lw=lw, c=c, label="E "+fname  ); 
    if Eref is not None: plt.plot( xs, (E-Eref)*scErr, ls='--', lw=lw, c=c, label=("E_err*%3.0f" %scErr )+fname  ); 
    plt.ylim(Emin, -Emin); 
    plt.legend(); 
    plt.grid()
    
    plt.subplot(2,1,2); 
    plt.plot(                          xs, F,              ls=ls, lw=lw, c=c, label="Fz "+fname ); 
    if Fref is not None: plt.plot( xs, (F-Fref)*scErr, ls='--', lw=lw, c=c, label=("Fz_err*%3.0f" %scErr )+fname  ); 
    plt.ylim(Fmin, -Fmin); 
    plt.legend(); 
    plt.grid()

    return Emin, Fmin




plt.figure(figsize=(5,10))
_,Fr,Er = getEF( 'gridFF_EFprofile_mod0.log' )
plot_fname( 'gridFF_EFprofile_mod0.log', ls=':', lw='1',  c='k' );
#plot_fname( 'gridFF_EFprofile_mod1.log', Eref=Er,Fref=Fr, c='r' );
plot_fname( 'gridFF_EFprofile_mod2.log', Eref=Er,Fref=Fr, c='g' );
plot_fname( 'gridFF_EFprofile_mod6.log', Eref=Er,Fref=Fr, c='b' );


plt.figure(figsize=(5,10))
_,Fr,Er = getEF( 'gridFF_EFprofile_mod0_x.log' )
plot_fname( 'gridFF_EFprofile_mod0_x.log', ls=':', lw='1',  c='k', iax=0 );
#plot_fname( 'gridFF_EFprofile_mod1_x.log', Eref=Er,Fref=Fr, c='r', iax=0 );
plot_fname( 'gridFF_EFprofile_mod2_x.log', Eref=Er,Fref=Fr, c='g', iax=0 );
plot_fname( 'gridFF_EFprofile_mod6_x.log', Eref=Er,Fref=Fr, c='b', iax=0 );


plt.figure(figsize=(5,10))
_,Fr,Er = getEF( 'gridFF_EFprofile_mod0_y.log' )
plot_fname( 'gridFF_EFprofile_mod0_y.log', ls=':', lw='1',  c='k', iax=1 );
#plot_fname( 'gridFF_EFprofile_mod1_y.log', Eref=Er,Fref=Fr, c='r', iax=1 );
plot_fname( 'gridFF_EFprofile_mod2_y.log', Eref=Er,Fref=Fr, c='g', iax=1 );
plot_fname( 'gridFF_EFprofile_mod6_y.log', Eref=Er,Fref=Fr, c='b', iax=1 );


plt.show()