import os
import sys
import subprocess
import psutil
import time



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