
import os
import sys

def parseFile( fname ):
    fin = open(fname)
    bKer     = False
    kernels  = []
    ker_args = []
    for line in fin:
        wds = line.split()
        if len(wds)<1: continue
        if(bKer):
            if( wds[0] == "){" ):
                kernels.append( (ker_name,ker_args) )
                ker_args = []
                bKer     = False
            else:
                l = line.strip()
                if( l[:2]!="//" ):
                    ker_args.append(l)
                    
        elif( wds[0]== "__kernel" ):
            ker_name = wds[2]
            bKer     = True
    return kernels


def kernels2buffers( kernels ):
    buffers  = {}
    for ker in kernels:
        for l in ker[1]:
            ws = l.split()
            if(ws[0] == "__global"):
                buffers[ ws[2] ] = ws[1] 
    return buffers
                

if __name__ == "__main__":
    kernels = parseFile( sys.argv[1] ) 
    print("\n######### KERNELS ########\n")
    for ker in kernels:
        print("\n",ker[0])
        for l in ker[1]:
            print( "\t", l )
    buffers = kernels2buffers( kernels )
    print("\n\n######### BUFFERS ########\n")
    for key,val in buffers.items():
        print(val, key )