import sys
import matplotlib.pyplot as plt
import re
sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import plotUtils   as plu


# ============ Functions

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

def read_dat( fname, ni=0, nf=1, nn=0 ):
    #format:            1  -96.294471702523595       -251.76919147019100       -48.443292828581697       # HHH-hhS1_NNO-hpS1 HHH-hhS1_NNO-hpS1 
    f=open( fname, 'r' )
    ints  =[]
    floats=[]
    names =[] 
    for l in f:
        ws = l.split()
        ints  .append( [ int(ws[i])   for i in range(0    ,ni           ) ] )
        floats.append( [ float(ws[i]) for i in range(ni   ,ni+nf        ) ] )
        names .append( [ ws[i]        for i in range(ni+nf+1,ni+nf+1+nn ) ] )
    return ints,floats,names


def find_minim_energy_confs( Es, pair_names, Emax=1000, ipivot=0 ):
    Emins = { n:(Emax,-1) for n in pair_names }   # dictionary of minimum energies for each pair
    for i,n in enumerate(pair_names):             # find minimum energy for each pair
        e = Es[i][ipivot]                         # Es[i] is a list of energies for each pair (DFTB, B3LYP, DFTB+disp)
        if Emins[n][0]>e:                         # if current energy is lower than the minimum energy for this pair
            Emins[n]=(e,i)                        # update minimum energy for this pair, and index of the minimum energy
    return Emins