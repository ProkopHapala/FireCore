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

#workdir="/home/prokop/git/FireCore/tests/dftb/inputs/"
#workdir="/home/prokop/git/FireCore/tests/dftb/inputs2/"
#workdir="/home/prokop/git/FireCore/tests/dftb/inputs3/"
#workdir="/home/prokop/git/FireCore/tests/dftb/input_rigid/"
workdir="/home/prokop/git/FireCore/tests/dftb/input_relaxed_tight"


params = dftbu.default_params.copy()
params.update({
    "Optimizer"   : "Rational{}",
    "GradAMax"    : 1E-6,
    'Temperature' : 300
})

nproc_max = 6


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


def wait_proc_finis( procs ):
    for p, f in procs:
        p.wait()
        f.close()

# ======= Main

fnames = os.listdir(workdir)
fnames.sort()
print( fnames )


bFrags = True
all_dirs  = [] 
all_procs = []
procs     = []


os.chdir( workdir )
for fname in fnames:
    name, ename = os.path.splitext(fname)
    if ename == '.xyz':
        #print(  fname )
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
            dftbu.makeDFTBjob( enames=A.enames, fname=name1+"/dftb_in.hsd", params=params, opt=False )
            dftbu.makeDFTBjob( enames=B.enames, fname=name2+"/dftb_in.hsd", params=params, opt=False )

        os.mkdir(name)
        os.system('cp %s %s/input.xyz' %(fname,name) )  
        dftbu.makeDFTBjob( enames=atoms.enames, fname=name+"/dftb_in.hsd", params=params )
        all_dirs += dirs
        for dir in dirs:
            print( dir )
            os.chdir(dir)
            f = open("stdout.log",'w')
            p = subprocess.Popen(['dftb+',"> OUT"],stdout=f)
            procs.append((p, f))
            os.chdir( workdir )

            if len(procs) >= nproc_max:
                wait_proc_finis( procs )
                all_procs += procs
                print( "DONE %i jobs " %len(all_procs) )
                procs = []

wait_proc_finis( procs )
all_procs += procs
print( "DONE %i jobs " %len(all_procs) )
procs = []

'''
print( all_dirs ) 
for dir in all_dirs:
    os.chdir(dir)
    os.system('cat input.xyz geom.out.xyz > relax.xyz' ) 
    atoms = au.AtomicSystem("geom.out.xyz")
    #plotBondLenghts( atoms, fname=dir+".png", axes=(1,0) )
    os.chdir( workdir )
'''