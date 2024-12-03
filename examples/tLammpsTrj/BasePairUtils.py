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
        if len(ws)<1: continue
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
    unique_inds   = {}
    for i,n in enumerate(pair_names):             # find minimum energy for each pair
        e = Es[i][ipivot]  
        if n not in Emins: 
            unique_pnames.append(n)
            Emins[n]=(e,i)
            unique_inds[n] = i
        else:
            if Emins[n][0]>e:                         # if current energy is lower than the minimum energy for this pair
                Emins[n]=(e,i)                        # update minimum energy for this pair, and index of the minimum energy
                unique_inds[n] = i
    # print( " find_minim_energy_confs().check() " )
    # col_set = set()
    # for i,(n,e) in enumerate(Emins.items() ): 
    #     if n in col_set:
    #         print( i, "Already in the set ",  n )
    #     else:
    #         col_set.add( n )
    #     #print( i, n, e )
    unique_inds = [ unique_inds[n] for n in unique_pnames ]
    return Emins, unique_pnames, unique_inds

def check_list_uniqueness( lst ):
    col_set = set()
    for i,n in lst :   
        if n in col_set:
            print( "check_list_uniqueness[%i]" %i, " already in the set ", n )
        else:
            col_set.add( n )

# =============  From /home/prokop/Desktop/CARBSIS/Paolo/correlations/utils.py


#!/usr/bin/python

import numpy as np

def find_E_min( E_dct ):
    #Emin_dct={ name: ('',1e+8) for name in name_set }
    Emin_dct = {}
    for key,val in E_dct.items():
        name = pair_strip(key)
        if name not in Emin_dct:
            Emin_dct[name] = (key,val)
        else:
            val0 = Emin_dct[name]
            if val < val0[1]:
                #print( "replace ", val0, (key,val)  )
                Emin_dct[name] = (key,val)
    return Emin_dct

def makeAllPairs( Xs,Ys, dct ):
    pairs=[]
    for x in Xs:
        for y in Ys:
            name=x+"_"+y
            if name in dct:
                pairs.append( name )
    return pairs            
        
def name_to_class(name, trans={ 'N':'e','O':'e' } ):
    return name.split('-')[0].replace('N','e').replace('O','e')

def pair_strip(name):
    ws = name.split("S")
    x=ws[0]
    y="_".join( ws[1].split('_')[1:] )
    name_ = x+"_"+y
    #if(x == 'H-h_1'):
    #    print( name, " -> ", name_ )
    return name_


#for name in names:
#    print( name,"    ", name_to_class(name) )

def load_dict( fname, dct=None ):
    f = open(fname, 'r')
    if dct is None: dct={}
    for line in f:
        ws=line.split()
        dct[ws[2]] = float(ws[0])
    f.close()
    return dct
    

def makeEmaps( names, Emin_dct ):
    n = len(names)
    Exx  = np.zeros( n     )
    Emap = np.zeros( (n,n) )
    Emap[:,:] = np.nan
    for i in range(n):
        name1=names[i]
        for j in range(n):
            name2=names[j]
            pair_name = name1+"_"+name2
            E = np.nan
            try:
                E = Emin_dct[pair_name][1]
            except:
                try:
                    E = Emin_dct[name2+"_"+name1][1]
                except:
                    pass
            Emap[i,j] = E
            if i==j:
                Exx[i] = E
            if E == np.nan:
                print( pair_name, E )
    dT = Emap - Emap.transpose()
    print( " dT(Emap).max, min ", dT.max(), dT.min() )
    return Emap, Exx

def makeEcontrast( Emap ):
    E_contrast = Emap.copy()
    n=len(Emap)
    for i in range(n):
        for j in range(n):
            E_contrast[i,j] = min( Emap[i,i], Emap[j,j] ) - E_contrast[i,j]        
    return E_contrast
    
'''
def try_Combine( Emap, pairs, new_pair ):
    E    = Emap[ new_pair[0],new_pair[1] ]
    E_Xa = Emap[ pairs.flat, new_pair[0] ]
    E_Yb = Emap[ pairs.flat, new_pair[1] ]
    if np.any( E_Xa<E ) or np.any( E_Yb<E ):
        return False
    else:
        return True
'''

'''
def try_Combine( Emap, pairs, new_pair, dEmin=4.0 ):
    E    = Emap[ new_pair[0],new_pair[1] ] + dEmin
    E_Xa = Emap[ pairs.flat, new_pair[0] ]
    E_Yb = Emap[ pairs.flat, new_pair[1] ]
    if np.any( E_Xa<E ) or np.any( E_Yb<E ):
        return False
    else:
        return True
'''

def min_cross_E( Emap, pairs, p, Epairs ):
    E    = Emap[ p[0],p[1] ]
    
    Eps = np.maximum( Epairs, E )
    E_Xa = Emap[ pairs[:,0], p[0] ] - Eps 
    E_Xb = Emap[ pairs[:,1], p[0] ] - Eps
    E_Ya = Emap[ pairs[:,0], p[1] ] - Eps
    E_Yb = Emap[ pairs[:,1], p[1] ] - Eps
    
    #print( "E_Xa, E_Xb, E_Ya, E_Yb      ", E_Xa, E_Xb, E_Ya, E_Yb )
    #E_Xa -= Eps 
    #E_Xb -= Eps
    #E_Ya -= Eps
    #E_Yb -= Eps
    #print( "E_Xa, E_Xb, E_Ya, E_Yb -Eps ", E_Xa, E_Xb, E_Ya, E_Yb )
    
    #p1=(43,0); p2=(30, 27)
    #p1=(29, 25); p2=(0, 43) 
    #if ((p1==new_pair) and (pairs[0,0]==p2[0]) and (pairs[0,1]==p2[1]) ) or ((p2==new_pair) and (pairs[0,0]==p1[0]) and (pairs[0,1]==p1[1]) ):
    #    print( "new_pair ", new_pair, "pairs ",  list(pairs) )
    #    print( E, Epairs, Eps )
    #    print( "E_X ", new_pair[0], pairs[:], E_Xa, E_Xb )
    #    print( "E_Y ", new_pair[1], pairs[:], E_Ya, E_Yb )
    #    print( "dE ",  min( E_Xa.min(), E_Xb.min() ,   min( E_Ya.min(), E_Yb.min() ) )  )
    #    print( "E[%i,%i] %g " %(p1[0],p1[0],Emap[p1[0],p1[0]]) )
    #    print( "E[%i,%i] %g " %(p1[1],p1[1],Emap[p1[1],p1[1]]) )
    #    print( "E[%i,%i] %g " %(p2[0],p2[0],Emap[p2[0],p2[0]]) )
    #    print( "E[%i,%i] %g " %(p2[1],p2[1],Emap[p2[1],p2[1]]) )

    
    return min( E_Xa.min(), E_Xb.min() ,   min( E_Ya.min(), E_Yb.min() ) )
    
    
    
    #E_Xa = Emap[ pairs.flat, new_pair[0] ] 
    #E_Yb = Emap[ pairs.flat, new_pair[1] ]
    
    '''
    p1=(43,0)
    p2=(30, 27)
    if ((p1==new_pair) and (pairs[0,0]==p2[0]) and (pairs[0,1]==p2[1]) ) or ((p2==new_pair) and (pairs[0,0]==p1[0]) and (pairs[0,1]==p1[1]) ):
        print( new_pair, E, list(pairs.flat) )
        print( "E_Xa ", E_Xa-E )
        print( "E_Yb ", E_Yb-E )
        print( "dE ", min( E_Xa.min(), E_Yb.min() ) - E  )
    #print( "min_cross_E ", min( E_Xa.min(), E_Yb.min() ) - E, "   | ",  E_Xa.min(), E_Yb.min(), E  )
    '''
    
    #return min( E_Xa.min(), E_Yb.min() ) - E  

def check_min_gap_old( Emap, pairs ):
    E_mins = np.zeros(len(pairs))
    for i,p in enumerate(pairs):
        ps = np.delete( pairs, [i], axis=0 )
        print( ps )
        E      = Emap[ p[0],p[1] ]
        E_Xa   = Emap[ ps.flat, p[0] ]
        E_Yb   = Emap[ ps.flat, p[1] ]
        E_min  = min( E_Xa.min(), E_Yb.min() )
        #print( E, E_Xa, E_Yb )
        E_mins[i] = E_min 
    return E_mins.min() - E


def check_min_gap(Emap, pairs, bPrint=False):
    #print( "pairs: \n",  pairs )
    E_mins = np.zeros(len(pairs))
    all_indices = np.arange(len(pairs))
    #print( "all_indices: \n",  all_indices )
    for i, p in enumerate(pairs):
              
        E = Emap[p[0], p[1]]   
  
        if bPrint:
            emin = 1e+8
            for j, pj in enumerate(pairs):
                
                if (j!=i):
                    e=Emap[p[0],pj[0]]-E; emin=min(emin,e); print( p[0],pj[0], e )
                    e=Emap[p[0],pj[1]]-E; emin=min(emin,e); print( p[0],pj[1], e )
                    e=Emap[p[1],pj[0]]-E; emin=min(emin,e); print( p[1],pj[0], e )
                    e=Emap[p[1],pj[1]]-E; emin=min(emin,e); print( p[1],pj[1], e )
                else:
                    e=Emap[p[0],p[0]]-E; emin=min(emin,e); print( p[0],p[0], e )
                    e=Emap[p[1],p[1]]-E; emin=min(emin,e); print( p[1],p[1], e )
            print( "-- E[%i,%i] " %(p[0],p[1]), E, " emin ", emin )        
        
        other_indices = all_indices[all_indices != i]
        #print( "other_indices \n", other_indices )
        E_min = 1e+8
        if len(other_indices)>0:
            E_Xa = Emap[pairs[other_indices, 0], p[1]]
            E_Yb = Emap[pairs[other_indices, 1], p[0]]
            E_min = min(E_Xa.min(), E_Yb.min()) 
        
        E_Xb = Emap[pairs[:,1],p[1]]
        E_Ya = Emap[pairs[:,0],p[0]]
        
        E_min = min( E_min, min(E_Xb.min(), E_Ya.min()) ) - E
        E_mins[i] = E_min
    #print( E_mins )
    return E_mins.min()
   

def try_all_combs( Emap, alphabets, pairs, dEmin=4.0, bCheck=True ):
    alphabets_new = set()
    for alphabet in alphabets:
        #alph = set( alphabet.flat )
        np_alph  = np.array(list(alphabet))
        
        Epairs = Emap[ np_alph[:,0], np_alph[:,1] ]   #;print(Epairs.shape)
        letters = set( np_alph.flat )
        #letters = { p[0],p[1] for p in alphabet  }
        #print( "np_alph ", np_alph )
        for pair in pairs:
            if (pair[0] not in letters) and (pair[1] not in letters):
            
            
                #a = alphabet | frozenset( [pair] )
                #if( a not in aplhabets_new ):
                #    dE = check_min_gap( Emap, np_alph )
                #    if ( dE >  dEmin ):
                #        aplhabets_new.add( a ) 
            
                #if try_Combine( Emap, np_alph, pair ):
                #if ( dE >  dEmin ):
                dE =  min_cross_E( Emap, np_alph, pair, Epairs )
                #print( dE, np_alph, pair )
                if dE>dEmin:
                    #a = alphabet + pair
                    #a = set(alphabet)
                    #a.add( pair )
                    a = alphabet | frozenset( [pair] )
                    
                    if( bCheck ):
                        dE_ = check_min_gap(Emap, np.array(list(a)) )
                        #if( np.abs(dE_-dE)>0.0001 ):
                        if( dE_<dEmin ):
                            print( "ERROR: min_cross_E(%g) != check_min_gap(%g) for " %(dE,dE_), pair, np_alph )
                            check_min_gap(Emap, np.array(list(a)), bPrint=True )
                            exit()
                    
                    if( a not in alphabets_new ):
                        alphabets_new.add( a )

           
    return  alphabets_new

def inds2names( inds, names ):
    return (names[inds[0]],names[inds[1]])

def alphabet2names( alph, names ):
    return [ (names[a[0]],names[a[1]]) for a in alph ]

def print_all_alphabets( alphabets, dEmin=4.0, Emap=None,  bCheckE=True,  bSort=True ):
    for a in alphabets:
        print( "alphabet_%03i " %len(a), end =" ")
        for p in a: print(p,             end =" ")
        if bCheckE:
            #print( a )
            dE = check_min_gap( Emap, np.array(list(a)) )
            print("dE: ", dE )
            if( dE<dEmin ): 
                print( "dE(%g) > dEmin(%g) => Exit()" %(dE,dEmin) )
                exit()
        else:
            print("")

def print_sorted_alphabets( Emap, alphabets, dEmin=4.0, nMax=1000000, names=None ):
    alphabets = [  (check_min_gap(Emap, np.array(list(a))), a )  for a in alphabets ]
    alphabets.sort( key = lambda a: -a[0] )
    for i,(dE,a) in enumerate(alphabets):
        print( "alphabet_%03i %10.2f " %(len(a), dE), end =" ")
        if( names is not None ):
            aa = alphabet2names( a, names )
            e = [ Emap[n[0],n[1]] for n in a ]
        for j,p in enumerate(aa): print(p, "%10.2f" %-e[j],             end =" ")
        print("")            
        if(i+1>=nMax): break

def findAlphabetsForRange( Emap, EbindRange, dEmin=4.0, nPairMax=4,E_contrast=None, verbosity=1, nMax=1000000, nMax2=5, names=None, excluded_names=None ):
    
    if E_contrast is None: E_contrast = uu.makeEcontrast( Emap )
    
    # --- select basepairs within EbindRange with minumum contrast dEmin
    mask   = np.logical_and( E_contrast > dEmin, np.logical_and( Emap>EbindRange[0], Emap<EbindRange[1] ) )   #;print("mask ", mask)
    inds   = np.where( mask );                                       #; print( "inds:\n", inds )

    #print( "inds", inds )
    #if excluded_names is not None:
    #    inds = np.array( [ p for p in inds if (names[p[0]] not in excluded_names) and (names[p[1]] not in excluded_names) ])
    
    #ind1d  = np.union1d( inds[0], inds[1] )                          #; print( "\nind1d:\n", ind1d )
    #pairs  = [ (i[0],i[1]) for i in np.array(inds).transpose() if i[0]!=i[1] ]    #; print( "pairs ", pairs ) 
    pairs  = [ (i[0],i[1]) for i in np.array(inds).transpose() if (i[0]<i[1]) ]    #; print( "pairs ", pairs ) 

    if excluded_names is not None:
        pairs = [ p for p in pairs if (names[p[0]] not in excluded_names) and (names[p[1]] not in excluded_names) ]

    levels = []

    alphs = { frozenset( [ p ] ) for p in pairs }                    #;print( "alphs1 ", alphs1 ) 
    if len(alphs)>0: levels.append( alphs )
    if (verbosity>0) and (nMax2>0):
        print("#### Alphabets(m_pairs=%i) found(%i) EB(%g,%g) dE>%g " %(1,len(alphs),EbindRange[0],EbindRange[1],dEmin) ); 
        #print_all_alphabets( alphs, Emap=Emap )
        print_sorted_alphabets( Emap, alphs, dEmin=dEmin, nMax=nMax2, names=names )

    for npair in range(2,nPairMax+1):
        alphs = try_all_combs( Emap, alphs, pairs, dEmin=dEmin ) 
        if len(alphs)>0: 
            levels.append( alphs )
        else: 
            break
        if verbosity>0:                         
            print("#### Alphabets(m_pairs=%i) found(%i) EB(%g,%g) dE>%g " %(npair,len(alphs),EbindRange[0],EbindRange[1],dEmin) ); 
            #print_all_alphabets( alphs, Emap=Emap )
            print_sorted_alphabets( Emap, alphs, dEmin=dEmin, nMax=nMax, names=names )
    
    return levels, inds
        


