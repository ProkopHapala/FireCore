import sys
import os
import matplotlib.pyplot as plt
from   matplotlib import collections  as mc
sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import plotUtils   as plu


# ============ Setup

sz      =5
dirname ='.'
rec_lelve = 0
if len(sys.argv)>1: dirname=sys.argv[1]
if len(sys.argv)>2: rec_lelve=int(sys.argv[2])

# ============ Functions

def plotBondLenghts( mol, axes=(0,1), HBs=None, Bs=None, bBls=False ):
    plu.plotSystem( mol,                    axes=axes, bBonds=False, bLabels=False )
  
    # Bonds
    if Bs is None:
        bs,rbs = mol.findBonds ( )                 #;print( bs, rbs )
    else:
        bs,rbs = Bs
    if bBls:
        rb_labs = [ ("%3.2f" %r) for r in rbs ] 
    else:
        rb_labs = None
    plu.plotBonds ( ps=mol.apos, links=bs,  axes=axes, colors="#808080", labels=rb_labs )
    
    # H-Bonds
    if HBs is None:
        hbs,rhbs = mol.findHBonds( bPrint=True )     #;print( hbs, rhbs )
    else:
        hbs,rhbs = HBs
    if len(hbs)>0:
        rh_labs = [ ("%3.2f" %r) for r in rhbs ] 
        plu.plotBonds ( ps=mol.apos, links=hbs, axes=axes, colors="g",       labels=rh_labs )
    return (bs,rbs),(hbs,rhbs)

def plot_fname_list(names):
    for name in names:
        name_,ext = os.path.splitext(name)
        if ext != '.xyz': continue
        print(name)
        mol = au.AtomicSystem( name ) 

        plt.figure(figsize=(2*sz,sz))

        plt.subplot(1,2,1); Bs,HBs = plotBondLenghts( mol, axes=axes );  plt.xlabel(ax_names[axes[0]]);plt.ylabel(ax_names[axes[1]])
        mol.orientPCA()
        #mol.saveXYZ("post_PCA."+name  ,   other_lines=other_lines)
        plt.subplot(1,2,2);          plotBondLenghts( mol, axes=axes, Bs=Bs, HBs=HBs );  plt.xlabel(ax_names2[axes[0]]);plt.ylabel(ax_names2[axes[1]])
        
        plt.savefig(name+".png", bbox_inches='tight' )
        plt.close(plt.gcf())

# ============ Main

ax_names =['x','y','z']
ax_names2=['PCA.a','PCA.b','PBC.x']
axes=(0,1)

bak_dir = os.getcwd()
if rec_lelve>0:
    # only sub-directories
    dirs = os.listdir( dirname )   ;print( "dirs ",  dirs )
    #dirs = [ f for f in dirs if os.path.isdir(f) ] 
    print( "names ",  dirs )
    for d in dirs:
        print( "SUBDIR: ",  d )
        os.chdir( os.path.join(dirname,d) )
        plot_fname_list( os.listdir( '.' ) )
        os.chdir( bak_dir )
else:
    plot_fname_list( os.listdir( dirname ) )




#plt.show()
