import sys
import os
import subprocess
import time
import numpy as np
import matplotlib.pyplot as plt

# https://stackoverflow.com/questions/16450788/python-running-subprocess-in-parallel

sys.path.append("../../")
from pyBall import dftb_utils as dftbu
from pyBall import atomicUtils as au
from pyBall import plotUtils   as plu


# ========= Setup

#workdir="/home/prokop/git/FireCore/tests/dftb/inputs/"
#workdir="/home/prokop/git/FireCore/tests/dftb/inputs2/"
#workdir="/home/prokop/git/FireCore/tests/dftb/inputs3/"
#workdir="/home/prokop/git/FireCore/tests/dftb/input_rigid/"
#workdir="/home/prokop/git/FireCore/tests/dftb/inputs_relaxed"
workdir="/home/prokop/git/FireCore/tests/dftb/input_relaxed_tight"

# ====== Functions

def plotBondLenghts( mol, fname="geom.png", axes=(0,1) ):
    hbs,rhbs = mol.findHBonds( bPrint=True )     #;print( hbs, rhbs )
    bs,  rbs = mol.findBonds ( )                 #;print( bs, rbs )
    rh_labs = [ ("%3.2f" %r) for r in rhbs ] 
    rb_labs = [ ("%3.2f" %r) for r in rbs ] 
    plu.plotSystem( mol,                    axes=axes, bBonds=False, bLabels=False )
    plu.plotBonds ( ps=mol.apos, links=bs,  axes=axes, colors="#808080", labels=rb_labs )
    if len( hbs )>0:
        plu.plotBonds ( ps=mol.apos, links=hbs, axes=axes, colors="g",       labels=rh_labs )
    plt.savefig( fname, bbox_inches='tight' )
    plt.close( plt.gcf() )

def getBindingEnergy( name ):
    outf =  "E_%s.log" %name
    os.system( "> "+outf  )
    s = 'grep "Total Energy:" %s/stdout.log | tail -1 >> ' + outf
    os.system( s %name )
    if bFrags:
        name1 = name+"_1"
        name2 = name+"_2"   
        os.system( s %name1 )
        os.system( s %name2 )
    data = np.genfromtxt( outf )
    dE_eV =  data[0,-2]- data[1,-2]- data[2,-2]
    print( "Binding Eenergy(%s): %g[eV] %g[kcal/mol] " %(fname,dE_eV,dE_eV*23.0609 ) )

    try:
        atoms = au.AtomicSystem("%s/geom.out.xyz" %name )
        plotBondLenghts( atoms, fname="opt/"+name+".png", axes=(1,0) )
        os.system( "cp %s/geom.out.xyz opt/%s.xyz " %(name,name) )
    except Exception as e:
        print( e )

    return dE_eV

# ======= Main

fnames = os.listdir(workdir)
fnames.sort()
print( fnames )

bFrags = True
all_dirs  = [] 
processes = []
os.chdir( workdir )
for fname in fnames:
    name, ename = os.path.splitext(fname)
    if ename == '.xyz':
        try:
            getBindingEnergy( name )
        except:
            print( "ERROR in ", name  )
        