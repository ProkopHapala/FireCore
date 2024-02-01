'''
This will generate headers from OpenCL
'''

import os
import sys


class Parser_OpenCL():
    """
    A class for parsing OpenCL files and extracting kernels and buffers.

    Methods:
    - parseFile(fname): Parses the specified OpenCL file and returns a list of kernels.
    - kernels2buffers(kernels): Converts a list of kernels into a dictionary of buffers.
    """

    def __init__(self, fname):
        self.kernels = self.parseFile( fname )
        self.buffers = self.kernels2buffers( self.kernels )

    def parseFile(self, fname ):
        """
        Parses the specified OpenCL file and returns a list of kernels.

        Args:
        - fname: The path to the OpenCL file.

        Returns:
        - A list of tuples, where each tuple contains the name of a kernel and its arguments.
        """
        fin = open(fname)
        bKer     = False   # Am I in kernell header right now ?
        kernels   = []
        ker_args  = []
        for line in fin:
            wds = line.split()
            if len(wds)<1: continue
            if(bKer):
                if( wds[0] == "){" ):   # seek for end of kernel header
                    kernels.append( (ker_name,ker_args) )
                    ker_args = []
                    bKer     = False
                else:
                    l = line.strip()
                    if( l[:2]!="//" ):  # if not commented append 
                        ker_args.append(l)
            elif( wds[0]== "__kernel" ):  # if line starts with kernell
                ker_name = wds[2].split('(')[0]
                bKer     = True
        return kernels

    def kernels2buffers(self, kernels ):
        """
        Converts a list of kernels into a dictionary of buffers.

        Args:
        - kernels: A list of tuples, where each tuple contains the name of a kernel and its arguments.

        Returns:
        - A dictionary where the keys are buffer names and the values are buffer types.
        """
        buffers  = {}
        for ker in kernels:
            for l in ker[1]:
                ws = l.split()
                if(ws[0] == "__global"):
                    buffers[ ws[2] ] = ws[1]
        return buffers
    
    def printKernels(self, outf=sys.stdout ):
        for kernel in self.kernels:
            outf.write("\n"+kernel[0]+"\n"   )
            outf.write("(\n")
            for arg in kernel[1]:
                outf.write(arg+"\n")  
            outf.write(")\n")

    # writeBuffers writes out buffers nicely to screen
    def printBuffers(self, outf=sys.stdout ):
        for key,val in self.buffers.items():
            outf.write(val +" "+ key+"\n" )
              
if __name__ == "__main__":

    fname = "/home/prokop/git/FireCore/cpp/common_resources/cl/myprog.cl"
    if( len(sys.argv)>1 ):
        fname = sys.argv[1]

    parser = Parser_OpenCL(fname)
    #kernels = parser.parseFile( fname ) 

    print("\n######### KERNELS ########\n")
    parser.printKernels()
    print("\n\n######### BUFFERS ########\n")
    parser.printBuffers()