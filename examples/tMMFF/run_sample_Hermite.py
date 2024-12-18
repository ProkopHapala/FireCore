import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
#from pyBall import atomicUtils as au
# from pyBall import FunctionSampling as fu
from pyBall import MMFF as mmff
from pyBall.tests import HermiteSpline as hsp

R0 = 3.5
E0 = 1.0
a  = 1.8

#hsp.test_1D( lambda x:fu.getLJ(x,R0,E0), lambda x:fu.getLJSplit(x,R0,E0), title="LJ simple",      mode=1 )
#hsp.test_1D( lambda x:fu.getLJ(x,R0,E0), lambda x:fu.getLJSplit(x,R0,E0), title="LJ deriv",       mode=2 )
#hsp.test_1D( lambda x:fu.getLJ(x,R0,E0), lambda x:fu.getLJSplit(x,R0,E0), title="LJ deriv split", mode=3 )

#hsp.test_1D( lambda x:fu.getMorse(x,R0,E0), lambda x:fu.getMorseSplit(x,R0,E0), title="Morse simple",      mode=1 )
#hsp.test_1D( lambda x:fu.getMorse(x,R0,E0), lambda x:fu.getMorseSplit(x,R0,E0), title="Morse deriv",       mode=2 )
#hsp.test_1D( lambda x:fu.getMorse(x,R0,E0), lambda x:fu.getMorseSplit(x,R0,E0), title="Morse deriv split", mode=3 )

#fu.checkNumDeriv( lambda x: fu.getLJ(x,3.5,1.0),        2.0, 6.0, tol=1e-6, n=400, errSc=100.0, plt=plt, label="LJ"   ,c='r' )
#fu.checkNumDeriv( lambda x: fu.getMorse(x,3.5,1.0,1.8), 2.0, 6.0, tol=1e-6, n=400, errSc=100.0, plt=plt, label="Morse",c='b' )

#hsp.test_1D( title="No-Half", mode=1)
#hsp.test_1D( title="No-Half", mode=1)


#hsp.test_fit_2D( title="test mode=1", mode=1 )
#hsp.test_fit_2D( title="test mode=2", mode=3 )

mmff.initParams()

#hsp.test_gridFF( mode=4, title="Hybrid Hermite tri-cubic\n(z-cut)" )
hsp.test_gridFF( mode=6, title="Bspline (from HH)\n(z-cut)" )

#hsp.test_gridFF_lat( mode=1, title="tri-linar force"          )
#hsp.test_gridFF_lat( mode=4, title="Hybrid Hermite tri-cubic" )

#Q = 0.0
Q = 0.4
#p0 = [1.0,-5.05,2.0]
#p0 = [0.0,0.0,2.0]
p0 = [-2.0,-2.0,0.0]

#hsp.test_gridFF_lat( mode=1, title="tri-linar force \n(y-cut)"         , Q=Q, p0=p0, iax=1, tmin=-10, tmax=10 )
#hsp.test_gridFF_lat( mode=4, title="Hybrid Hermite tri-cubic \n(y-cut)", Q=Q, p0=p0, iax=1, tmin=-10, tmax=10. )

#hsp.test_gridFF_lat( mode=1, title="tri-linar force \n(x-cut)"         , Q=Q, p0=p0, iax=0 )
#hsp.test_gridFF_lat( mode=4, title="Hybrid Hermite tri-cubic \n(x-cut)", Q=Q, p0=p0, iax=0 )

#hsp.test_gridFF_lat( mode=4, title="Hybrid Hermite tri-cubic", Q=0.0, p0=[1.0,1.05,2.0], iax=1 )
#hsp.test_gridFF_lat( mode=4, title="Hybrid Hermite tri-cubic", Q=0.0, p0=[1.0,1.05,2.0], iax=0 )

hsp.test_gridFF_lat( mode=6, title="Bspline tri-cubic", Q=0.0, p0=p0, iax=0 )

hsp.test_gridFF( mode=1, title="tri-linar force \n(z-cut)"          )
hsp.test_gridFF_lat( mode=1, title="tri-Linear Force", Q=0.0, p0=p0, iax=0 )

plt.show()