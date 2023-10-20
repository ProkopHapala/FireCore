import sys
import os
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import plotUtils   as plu

import BasePairUtils as bpu

# ============ Setup 

#Econtrast_min  = 10.0
#Econtrast_min  = 5.0

dEcontrast_min = 10.0


dirname=sys.argv[1]

files = os.listdir( dirname )
files.sort()
lss=['-', '--', ':', '-.']

# ======== Functions

def read_dat( fname, ni=2, nf=5 ):
    '''
    Read a .dat file with columns of data
    '''
    fin = open( fname, 'r' )
    names = []
    fdata = []
    for line in fin:
        ws = line.split()
        if len(ws)<nf: continue
        fdata.append([ float(w) for w in ws[ni:nf] ])      
        names.append( ws[0] )  
    fin.close()
    return names, np.array( fdata )


# ======== Main

xs=None
ifile=0
for fname in files:
    # if extensiton is not .dat, skip
    if not fname.endswith(".dat"): continue
    print( "ifile fname ", ifile, fname )
    # read file
    names, fdata = read_dat( dirname+"/"+fname, ni=2, nf=5 )
    #print( "fdata.shape ", fdata.shape )
    nnames = len(names)

    if ifile==0:
        plt.figure( figsize=( nnames*0.20+0.5,5) )
        xs = list( range( fdata.shape[0] ) )
    
    #print( "fdta ", fdata )
    #print( "len(xs), fdata.shape[0] ",  len(xs), fdata.shape[0] )
    #print( "xs ", xs )

    i_ = fname.find("_")
    label = fname.split(".")[0][i_+1:]

    # plot 2d data as imshow
    plt.imshow( fdata.T, origin='lower', cmap='bwr', vmin=-Econtrast_min, vmax=Econtrast_min, extent=[-0.5,nnames-0.5,-0.5,2.5] )

    plt.plot( xs,fdata[:,0], '.b', ls=lss[ifile], label=label )
    plt.plot( xs,fdata[:,1], '.r', ls=lss[ifile] )
    plt.plot( xs,fdata[:,2], '.g', ls=lss[ifile] )
    ifile+=1


font_prop = FontProperties(family='monospace', size=12 , weight='bold')
plt.xticks( xs, [ f"{n:<8}" for n in names], rotation='vertical', fontproperties=font_prop )  # set xticks to names of pairs in this class

plt.grid()
plt.ylabel( "Energy [kcal/mol]", fontproperties=font_prop ) 
plt.ylim(0.,+60)
plt.xlim(-0.5,nnames-0.5)
plt.yticks(np.arange( 0.0, 45.0+0.001, 5.0 ))


plt.legend()
plt.tight_layout()
plt.savefig( dirname+"/Econtrast_2endgroup.png", bbox_inches='tight' )
plt.savefig( dirname+"/Econtrast_2endgroup.svg", bbox_inches='tight' )

plt.show()