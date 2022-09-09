import numpy as np

# load atoms species parameters form a file ( currently used to load Lenard-Jones parameters )
def loadSpeciesLines( lines ):
    #params = [ ( float(l[0]), float(l[1]), float(l[2]), int(l[3]), l[4] ) for l in ( l.split() for l in lines )]  
    params = []
    for l in lines:
        l = l.split()
        if len(l) >= 5:
            # print l
            params.append( ( float(l[0]), float(l[1]), float(l[2]), int(l[3]), l[4] ) )
    return np.array( params, dtype=[('rmin',np.float64),('epsilon',np.float64),('alpha',np.float64),('atom',np.int),('symbol', '|S10')])

def loadSpecies(fname):
    try:
        with open(fname, 'r') as f:  
            str_Species = f.read(); 
    except:
        #if(verbose>0): print("defaul atomtypes.ini")
        fpath = os.path.dirname( os.path.realpath( __file__ ) ) + '/defaults/atomtypes.ini'
        print("loadSpecies from : ", fpath)
        with open(fpath, 'r') as f:
            str_Species = f.read();
    str_Species = "\n".join( "\t".join( l.split()[:5] )  for l in str_Species.split('\n')  )
    return loadSpeciesLines( str_Species.split('\n') )

def get_C612( i, j, FFparams ):
    '''
    compute Lenard-Jones coefitioens C6 and C12 pair of atoms i,j
    '''
#	print i, j, FFparams[i], FFparams[j]
    Rij = FFparams[i][0] + FFparams[j][0]
    Eij = np.sqrt( FFparams[i][1] * FFparams[j][1] )
    return 2*Eij*(Rij**6), Eij*(Rij**12)

def getAtomsLJ( iZprobe, iZs,  FFparams ):
    '''
    compute Lenard-Jones coefitioens C6 and C12 for interaction between atoms in list "iZs" and probe-particle "iZprobe"
    '''
    n   = len(iZs)
    cLJs  = np.zeros((n,2))
    for i in range(n):
        cLJs[i,0],cLJs[i,1] = get_C612( iZprobe-1, iZs[i]-1, FFparams )
    return cLJs

def getAtomsREA(  iZprobe, iZs,  FFparams, alphaFac=-1.0 ):
    '''
    compute Lenard-Jones coefitioens C6 and C12 for interaction between atoms in list "iZs" and probe-particle "iZprobe"
    '''
    n   = len(iZs)
    REAs  = np.zeros( (n,4) )
    i = iZprobe-1
    for ii in range(n):
        j = iZs[ii]-1
        #print ii, i, j
        REAs[ii,0] = FFparams[i][0] + FFparams[j][0]
        REAs[ii,1] = -np.sqrt( FFparams[i][1] * FFparams[j][1] )
        REAs[ii,2] = FFparams[j][2] * alphaFac
    return REAs     #np.array( REAs, dtype=np.float32 )