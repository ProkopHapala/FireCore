import sys
import os
import subprocess
import time
import matplotlib.pyplot as plt

# https://stackoverflow.com/questions/16450788/python-running-subprocess-in-parallel

sys.path.append("../../")
from pyBall import dftb_utils as dftbu
from pyBall import atomicUtils as au
from pyBall import plotUtils   as plu


# ========= Setup

workdir="/home/prokop/git/FireCore/tests/dftb/inputs/"

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
        print(  fname )
        atoms= au.AtomicSystem(fname)
        atoms.findBonds()
        dirs=[name]

        if bFrags:
            name1 = name+"_1"
            name2 = name+"_2"
            dirs += [name1,name2]      
            ins,outs = atoms.selectBondedCluster( {0} )
            A = atoms.selectSubset( ins  )
            B = atoms.selectSubset( outs )
            os.mkdir(name1)
            os.mkdir(name2)
            A.saveXYZ(name1+"/input.xyz", bQs=False)
            B.saveXYZ(name2+"/input.xyz", bQs=False)
            dftbu.makeDFTBjob( atoms=A, fname=name1+"/dftb_in.hsd" )
            dftbu.makeDFTBjob( atoms=B, fname=name2+"/dftb_in.hsd" )

        os.mkdir(name)
        os.system('cp %s %s/input.xyz' %(fname,name) )  
        dftbu.makeDFTBjob( atoms=atoms, fname=name+"/dftb_in.hsd" )
        all_dirs += dirs
        for dir in dirs:
            print( dir )
            os.chdir(dir)
            
            f = open("stdout.log",'w')
            p = subprocess.Popen(['dftb+',"> OUT"],stdout=f)
            processes.append((p, f))
            os.chdir( workdir )

for p, f in processes:
    p.wait()
    f.close()

print( all_dirs ) 
for dir in all_dirs:
    os.chdir(dir)
    os.system('cat input.xyz geom.out.xyz > relax.xyz' ) 
    atoms = au.AtomicSystem("geom.out.xyz")
    plotBondLenghts( atoms, fname=dir+".png", axes=(1,0) )
    os.chdir( workdir )