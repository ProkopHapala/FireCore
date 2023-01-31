import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

def process( fname, dct=None, prefix="#Epiece", bPrintSame=True, bPrintMissing=True, maxErr=1e-6 ):
    maxErr2 = maxErr**2
    bFillDct=False
    if dct is None:
        dct={}
        bFillDct=True
    for line in open(fname):
        wds=line.split()
        if( wds[0]==prefix ):
            key=wds[1]+" "+wds[2]
            val=float(wds[3])
            if(bFillDct):
                dct[key]=val
            else:
                if key in dct.keys():
                    vref = dct[key]
                    dval = val - vref
                    if( (dval*dval)>maxErr2 ):
                        print( key, val, dct[key], " IS DIFFERENT" )
                    else:
                        if(bPrintSame):print( key, val, dct[key] )
                else:
                    if(bPrintMissing):print( key, val, " NOT FOUND" )
    return dct



dct = process( "/home/prokop/git/pyeff/tests/CH4/output.txt"  )
#process( "./output.txt", dct  )
#process( "./output.txt", dct,  bPrintSame=False, bPrintMissing=False, maxErr=1e-6 )
process( "./output.txt", dct,  bPrintSame=False, bPrintMissing=False, maxErr=1e-4 )