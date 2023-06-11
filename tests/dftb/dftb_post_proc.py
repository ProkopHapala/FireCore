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

workdir="/home/prokop/git/FireCore/tests/dftb/inputs/"

# ====== Functions

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
    return dE_eV

# ======= Main

fnames = os.listdir(workdir)
print( fnames )

bFrags = True
all_dirs  = [] 
processes = []
os.chdir( workdir )
for fname in fnames:
    name, ename = os.path.splitext(fname)
    if ename == '.xyz':
        getBindingEnergy( name )
        