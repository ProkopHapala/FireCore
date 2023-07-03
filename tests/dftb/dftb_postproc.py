#!/usr/bin/python

import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt

from pyBall import atomicUtils as au
#from pyBall import dftb_utils as dftbu
#from pyBall import Process as pc
from pyBall import plotUtils as plu

names=[
"2_ammonia_ammonia",
"2_formamide_formamide",
"2_formicAcid_formicAcid", 
"2_HCN_HCN",
"2_water_water",

"2_carbonicAcid_carbonicAcid",
"2_cyanidenitrate_urea",
"2_formicAcid_formamide",
"2_nitro_diamine",
"2_urea_urea",
"2_water_ammonia",

"2_water_carbonicAcid",
"2_water_urea",

"ammonia",
"formamide",
"formicAcid",
"HCN",
"water",

"carbonicAcid",
"cyanidenitrate",
"diamine",
"nitro",
"urea",
]

# ========= Setup

#in_dir   = "/home/prokophapala/git/FireCore/tests/dftb/input_xyz"
#workdir = "/home/prokophapala/git/FireCore/tests/dftb/HBsmall_relax"
#workdir = cwd = os.getcwd()

in_dir   = "/home/prokophapala/git/FireCore/tests/dftb/HBsmall_relax"
out_dir  = "/home/prokophapala/git/FireCore/tests/dftb/HBsmall_relax/results"

# ========= Functions

def plotBondLenghts( mol, fname, axes=(0,1) ):
    plt.figure(figsize=(5,5))
    hbs,rhbs = mol.findHBonds( bPrint=True )     #;print( hbs, rhbs )
    bs,  rbs = mol.findBonds ( )                 #;print( bs, rbs )
    rh_labs = [ ("%3.2f" %r) for r in rhbs ] 
    rb_labs = [ ("%3.2f" %r) for r in rbs ] 
    plu.plotSystem( mol,                    axes=axes, bBonds=False, bLabels=False )
    plu.plotBonds ( ps=mol.apos, links=hbs, axes=axes, colors="g",       labels=rh_labs )
    plu.plotBonds ( ps=mol.apos, links=bs,  axes=axes, colors="#808080", labels=rb_labs )
    plt.savefig(fname+".png", bbox_inches='tight' )
    plt.close( plt.gca() )

def getE(fname):
    hosts = subprocess.check_output( 'grep "Total Energy" %s | tail -1' %fname, shell=True )
    ws = hosts.split()
    return float(ws[4])
    
def getBindingE(name):
    ws=name.split('_')
    E12 = getE( in_dir+'/'+name+"/stdout.log" )
    E1  = getE( in_dir+'/'+ws[1]+"/stdout.log" )
    E2  = getE( in_dir+'/'+ws[2]+"/stdout.log" )
    return E12-E1-E2

# ========= Main

for name in names:
    try:
        dirname=in_dir+'/'+name        
        # if( name[0]=='2' ):
        #     E = getBindingE( name )
        #     print( E*23.0609 , " [kcal/mol] ", name  )
        # else:
        #     E = getE( dirname+"/stdout.log" )
        #     print( E, " [eV]",  name )

        #os.system("cp %s/geom.out.xyz %s/%s.xyz" %(dirname,out_dir,name)   )

        atoms = au.AtomicSystem( dirname+"/geom.out.xyz" )
        plotBondLenghts( atoms, fname=out_dir+"/"+name, axes=(0,1) )
    except Exception as e:
        print( name, e )
    
print( "DONE " )

