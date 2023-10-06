import sys
import numpy as np
import matplotlib.pyplot as plt
import re
sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import plotUtils   as plu


# ============ Functions

# list all .xyz files in directory
#files  = [ f for f in os.listdir('.') if os.path.isfile(f) and f.endswith(".xyz") ]

#re_name_to_class = re.compile(r'([ON])[^-]*')
re_name_to_class = re.compile(r'([NO])' )

re_remove_Ss     = re.compile(r'_S-?\d'     )

def split_pair_with_S( name, bJoin=False ):
    '''
    Convert names format from HHH-hhS1_NNO-hpS1  to (HHH-hh,NNO-hp)
    '''
    ws = name.split('S')
    n1=ws[0]
    n2=  "_".join( ws[1].split('_')[1:] )

    #if 'NNO-hh' in name:
    #    print( name, n1, n2 )
    #    exit()
    if bJoin:
        return n1+"_"+n2
    else:
        return n1,n2

def name_to_class( name ):
    return re.sub( re_name_to_class, 'e', name.split('-')[0] )

def read_dat( fname, ni=0, nf=1, iname=0, toRemove=None ):
    #format:            1  -96.294471702523595       -251.76919147019100       -48.443292828581697       # HHH-hhS1_NNO-hpS1 HHH-hhS1_NNO-hpS1 
    f=open( fname, 'r' )
    ints  =[]
    floats=[]
    names =[] 
    pair_names = []

    nc = ni+nf+2 # number of columns in the file
    for l in f:
        ws = l.split()
        nw = len(ws)
        if(nw<nc):
            ints_i    = [-1    ]*ni
            floats_i  = [np.nan]*nf
            name      = ws[nw-1]
            pname     = split_pair_with_S( name, bJoin=False ) 
        else:
            ints_i   = [ int(ws[i])   for i in range(0    ,ni           ) ]
            floats_i = [ float(ws[i]) for i in range(ni   ,ni+nf        ) ]
            name     = ws[ni+nf+1+iname]
            pname    = split_pair_with_S( name, bJoin=False ) 

        if toRemove is not None:
            if (pname[0] in toRemove)or(pname[1] in toRemove):
                continue

        ints      .append( ints_i   )
        floats    .append( floats_i )
        names     .append( name )
        pair_names.append( pname )
    return ints,floats,names,pair_names

def try_nickname( name, nicknames, sep='\n' ):
    if name in nicknames:
        return sep+nicknames[name]
    else:
        return ""


def find_minim_energy_confs( Es, pair_names, Emax=-1000, ipivot=0 ):
    #Emins = { n:(Emax,-1) for n in pair_names }   # dictionary of minimum energies for each pair
    Emins         = {}
    unique_pnames = []
    for i,n in enumerate(pair_names):             # find minimum energy for each pair
        e = Es[i][ipivot]  
        if n not in Emins: 
            unique_pnames.append(n)
            Emins[n]=(e,i)
        else:
            if Emins[n][0]>e:                         # if current energy is lower than the minimum energy for this pair
                Emins[n]=(e,i)                        # update minimum energy for this pair, and index of the minimum energy
    return Emins, unique_pnames