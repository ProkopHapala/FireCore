import sys
import os
import subprocess
import time

# https://stackoverflow.com/questions/16450788/python-running-subprocess-in-parallel

sys.path.append("../../")
from pyBall import dftb_utils as dftbu
from pyBall import atomicUtils as au

# ========= Setup

workdir="/home/prokop/git/FireCore/tests/dftb/inputs/"
fnames = os.listdir(workdir)
print( fnames )

processes = []
os.chdir( workdir )
for fname in fnames:
    name, ename = os.path.splitext(fname)
    print(  name, ename )
    if ename == '.xyz':
        atoms= au.AtomicSystem(fname)
        os.mkdir(name )
        os.chdir(name )
        os.system('cp ../%s input.xyz' %fname )  
        dftbu.makeDFTBjob( atoms=atoms )
        f = open("stdout.log",'w')
        p = subprocess.Popen(['dftb+',"> OUT"],stdout=f)
        processes.append((p, f))
        os.chdir( workdir )

for p, f in processes:
    p.wait()
    f.close()