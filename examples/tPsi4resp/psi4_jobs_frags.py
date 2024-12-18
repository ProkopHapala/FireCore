#!/usr/bin/python


import os
import sys
import subprocess
import psutil
import time
from pyBall import  psi4_utils  as psi4u
from pyBall import  atomicUtils as au

names=[
"2_ammonia_ammonia",
"2_formamide_formamide",
"2_formicAcid_formicAcid", 
"2_HCN_HCN",
"2_water_ammonia", 
"2_water_water",
"ammonia",
"formamide",
"formicAcid",  
"HCN",
"water",
]

methods=[
#"scf",
#"pbe",
#"pbe-d3",
#"b3lyp",
"b3lyp-d3",
#'m06',
#'m06-d3',
#"mp2",
#"mp2.5",
#"omp2.5",
#"ccsd",
]

basiss=[
#"cc-pVDZ",
"aug-cc-pVDZ",
#"aug-cc-pVTZ",
"def2-SVP",
"def2-SVPD",
"def2-TZVPP",
#"def2-TZVPPD",
#"aug-cc-pVTZ",
#"def2-TZVPPD",
#"def2-QZVPPD",
]

 # ========= Setup

nproc_max = 10

params_block='''
    scf_type df 
    
    opt_type MIN
    geom_maxiter 1000
    g_convergence INTERFRAG_TIGHT
    print_trajectory_xyz_file true
    DYNAMIC_LEVEL 1
''' 

# ========= Functions

class Proces:
    def __init__(self, name, args ):
        self.name    = name
        self.rss_max = 0
        self.args    = args
        self.fout    = open('stdout.log','w')
        self.t0      = time.time()
        self.proc    = subprocess.Popen( args, stdout=self.fout )
        print( "START(%s)" %self.name )

    def live(self):
        #print( self.name+".live()" )
        if( self.proc.poll() == None ):   # running
            procs        = psutil.Process( self.proc.pid )
            mem          = procs.memory_info()
            self.rss_max = max(self.rss_max,mem.rss)
            return True 
        else:
            print( "FINISHED(%s) mem_max %g time %g " %(self.name, self.rss_max/1073741824.0, time.time()-self.t0)  )
            return False

    def __del__(self):
        self.fout.close()

def wait_proc_finish( procs ):
    while len(procs)>0:
        procs = [ p for p in procs if p.live() ]
        
def addJobWhenFree( name, args, procs, nproc_max=10, wait_sec=0.1 ):
    while True:
        if len(procs)<nproc_max:
            procs.append( Proces( name, args ) )
            break
        else:
            time.sleep( wait_sec )
            proc_new = []
            for proc in procs:
                if proc.live():
                    proc_new.append( proc )
            procs = proc_new
    return procs
            
# ========= Main

all_procs = []
procs     = []

workdir = cwd = os.getcwd()

bFrags = False
bOpt   = True
mem = '15GB'




for method in methods:
    for bas in basiss:
        out_prefix = method+"-"+bas     
        try:  
            os.mkdir( out_prefix )
        except:
            pass  
        os.chdir( workdir )
        for name in names:
            atoms = au.AtomicSystem('input_xyz/'+name+".xyz")
            
            if( name[0]=='2' ):
                ins,outs = atoms.selectBondedCluster( {0} )
                nhyphen=len(ins)-1
            else:
                nhyphen=None
            lines = atoms.toLines()    
            dirname = out_prefix+"/"+name
            #print( dirname )
            os.mkdir( dirname )   

            # ----- run jobs locally
           
            if bFrags:
                dirname1 = name+"_1"
                dirname2 = name+"_2"
                dirs += [name1,name2]      
                os.mkdir(dirname)
                os.mkdir(dirname)
                #dftbu.makeDFTBjob( enames=A.enames, fname=name1+"/dftb_in.hsd", params=params, opt=opt_frag )
                #dftbu.makeDFTBjob( enames=B.enames, fname=name2+"/dftb_in.hsd", params=params, opt=opt_frag )
                psi4u.write_psi4_in( lines[ins] ,  fname=dirname+"/psi.in", mem=mem, method=method, basis=bas, bsse=None, params_block=params_block, opt=bOpt )
                psi4u.write_psi4_in( lines[outs],  fname=dirname+"/psi.in", mem=mem, method=method, basis=bas, bsse=None, params_block=params_block, opt=bOpt )
                
                os.chdir( dirname1 )
                procs = addJobWhenFree( dirname1, ['psi4',"psi.in"], procs, nproc_max=nproc_max, wait_sec=0.5 )
                os.chdir( workdir )
                
                os.chdir( dirname2 )
                procs = addJobWhenFree( dirname2, ['psi4',"psi.in"], procs, nproc_max=nproc_max, wait_sec=0.5 )
                os.chdir( workdir )
            
            psi4u.write_psi4_in( lines, fname=dirname+"/psi.in", nhyphen=nhyphen,  mem=mem, method=method, basis=bas, bsse="'cp'", params_block=params_block, opt=bOpt )
            os.chdir( dirname )
            procs = addJobWhenFree( dirname, ['psi4',"psi.in"], procs, nproc_max=nproc_max, wait_sec=0.5 )
            os.chdir( workdir )
            
wait_proc_finish( procs )
print( "DONE " )

