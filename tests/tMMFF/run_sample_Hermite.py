import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF as mmff
from pyBall import FunctionSampling as fu


# =============

apos = [
    [-2.0,0.0,0.0],
    [ 2.0,0.0,0.0],
]
REs=[
    [3.5,1.0],
    [3.5,1.0],
]

# =============

### ==== ToDo: Extract these functions to GridUtils.py (?)


def test_1D( func, funcS, g0=2.0, gmax=10.0, dg=0.2, dsamp=0.02, scErr=100.0, title=None, mode=1 ):
    xs  = np.arange(g0, gmax+1e-8, dg)     ; ng=len(xs)
    xs_ = np.arange(g0, gmax+1e-8 -3*dg, dsamp)  ; nsamp=len(xs_)
    
    #E0 = 1.0
    #R0 = 3.5
    #a  = 1.7
    E_,F_ = func( xs_ )
    if   ( mode==1 ):   # simple Hermite using just values
        E,F   = func( xs)
        Gs = E
        print( "xs ",xs )
        FEout = mmff.sample_SplineHermite( xs_, Gs, g0=g0, dg=dg )
    elif ( mode==2 ):   # simple Hermite using derivatives
        E,F   = func( xs )
        Gs    = np.zeros( (ng,2) )
        Gs[:,0] = E
        Gs[:,1] = F
        FEout = mmff.sample_SplineHermite1D_deriv( xs_, Gs, g0=g0, dg=dg )
    elif ( mode==3 ):   # hermite using linear combination of splited function
        Ep,El,Fp,Fl = funcS( xs )
        Gs = np.zeros( (ng,4) )
        Gs[:,0] = Ep
        Gs[:,1] = Fp*dg
        Gs[:,2] = El
        Gs[:,3] = Fl*dg
        FEout = mmff.sample_SplineHermite_comb( xs_, Gs, [1.0,1.0], ncomb=2, g0=g0, dg=dg )

    # ------ plotting

    Emin =  E_.min()
    Fmin = -F_.max()

    plt.figure(figsize=(5,10))
    plt.subplot(2,1,1)    
    plt.plot( xs_, E_,                    "-k", lw=0.5,  label="E_ref" )
    plt.plot( xs_,  FEout[:,0]          , "-b", lw=0.5, label="E_fit" )
    plt.plot( xs_, (FEout[:,0]-E_)*scErr, "-r", lw=0.5, label=("error*%g" % scErr) )
    #print( "Gs: ", Gs )
    plt.ylim(Emin*1.2,-Emin*1.2)
    plt.grid()
    plt.legend()
    plt.title("Energy")
    plt.subplot(2,1,2)   
    plt.plot( xs_, -F_ ,                  "-k", lw=0.5, label="F_ref" )    
    plt.plot( xs_, -FEout[:,1],           "-b", lw=0.5, label="F_fit" )
    plt.plot( xs_, (FEout[:,1]-F_)*scErr, "-r", lw=0.5, label=("error*%g" % scErr) )
    plt.ylim(Fmin*1.2,-Fmin*1.2)
    plt.legend()
    plt.title("Force")
    plt.grid()
    if title is not None: plt.suptitle(title)
    #plt.show()


R0 = 3.5
E0 = 1.0
a  = 1.8


#test_1D( lambda x:fu.getLJ(x,R0,E0), lambda x:fu.getLJSplit(x,R0,E0), title="LJ simple",      mode=1 )
#test_1D( lambda x:fu.getLJ(x,R0,E0), lambda x:fu.getLJSplit(x,R0,E0), title="LJ deriv",       mode=2 )
#test_1D( lambda x:fu.getLJ(x,R0,E0), lambda x:fu.getLJSplit(x,R0,E0), title="LJ deriv split", mode=3 )

test_1D( lambda x:fu.getMorse(x,R0,E0), lambda x:fu.getMorseSplit(x,R0,E0), title="Morse simple",      mode=1 )
test_1D( lambda x:fu.getMorse(x,R0,E0), lambda x:fu.getMorseSplit(x,R0,E0), title="Morse deriv",       mode=2 )
test_1D( lambda x:fu.getMorse(x,R0,E0), lambda x:fu.getMorseSplit(x,R0,E0), title="Morse deriv split", mode=3 )

#fu.checkNumDeriv( lambda x: fu.getLJ(x,3.5,1.0),        2.0, 6.0, tol=1e-6, n=400, errSc=100.0, plt=plt, label="LJ"   ,c='r' )
#fu.checkNumDeriv( lambda x: fu.getMorse(x,3.5,1.0,1.8), 2.0, 6.0, tol=1e-6, n=400, errSc=100.0, plt=plt, label="Morse",c='b' )


#test_1D( title="No-Half", mode=1)
#test_1D( title="No-Half", mode=1)

plt.show()