import numpy             as np
import matplotlib.pyplot as plt

'''

We should use the relation e^x = lim_{n->inf} (1+x/n)^n   to apporximate the exponential function with a polynomial.

slope of 1-x^2 at point x=1 is -2 

'''

def exp_approx( x,  Ri=1.4,  Rcut=6.0, n=8 ):
    x_ = 1-x/Rcut
    k = 1/(1-Ri/Rcut)
    y = (x_*k)**n
    mask = x>Rcut
    y[mask] = 0
    return y

def exp_approx_2( x, Ri=1.4, Rcut=5.0, n=8 ):
    x_ = 1-(x/Rcut)**2
    k  = 1/(1-(Ri/Rcut)**2)
    y = (x_*k)**n
    mask = x>Rcut
    y[mask] = 0
    return y

def numDeriv( x, y ):
    dx = x[2:]-x[:-2]
    dy = y[2:]-y[:-2]
    x_ = x[1:-1]
    return -dy/dx, x_

def Morse( r, Ei=0.1, Ri=1.4, b=1.6 ):
    expr = np.exp( -b*( r - Ri ) )
    E    =      Ei*( expr*expr - 2*expr )
    F    =  2*b*Ei*( expr*expr -   expr ) 
    return E,F, expr

def LennardJones( r, Ei=0.1, Ri=1.4 ):
    ir6  = (Ri/r)**6
    ir12 = ir6*ir6
    E    = Ei*   ( ir12 - 2*ir6 )
    F    = Ei*12*( ir12 -  ir6  )/r
    return E,F

def MorseCut_1( r, Ei=0.1, Ri=1.4, Rc=6.0, n=8 ):
    x_  = 1-r/Rc
    k   = 1/(1-Ri/Rc)
    e   = (x_*k)**n
    mask    = r>Rc
    e[mask] = 0
    E =    Ei*(e*e - 2*e)
    F =  2*Ei*(e*e -   e)
    return E,F, e

def MorseCut_2( r, Ei=0.1, Ri=1.4, Rc=6.0, n=8 ):
    r2 = r*r
    x_ = 1-(r2/(Rc*Rc))
    k  = 1/(1-(Ri/Rc)**2)
    e = (x_*k)**n
    mask = r2>(Rc*Rc)
    e[mask] = 0
    E =    Ei*(e*e - 2*e)
    F =  2*Ei*(e*e -   e)
    return E,F, e

def MorseCut_3( r, Ei=0.1, Ri=1.4, Rc=6.0 ):
    r2 = r*r
    x1 = 1-(r2/(Rc*Rc))
    x2 = 1-(r2/(0.7*Ri*Ri))
    x3 = 1-(r2/(2.0*Rc*Rc))
    mask = r2>(Rc*Rc)
    x1[mask] = 0
    #E =  5*Ei* x1*x1*x1*x1*( -1 + x2 )
    E =  35*Ei*x2*(x1**2)*(x3**16)
    #E =  10*Ei*x2*(x1**2)*(x3**4)
    F =  0
    return E,F

def MorseCut_4( r, Ei=0.1, Ri=1.4, Rc=6.0 ):
    r2 = r*r
    x1 = 1-(r2/(Rc*Rc))
    ir2 = (Ri*Ri*1.3)/r2 
    ir4 = ( (Ri*Ri*1.05)/r2 )**2
    #x2 = ir2*ir2 - 2*ir2
    x2 = ir4*ir4 - 2*ir4
    mask = r2>(Rc*Rc)
    x1[mask] = 0
    E =  1.8*Ei*(x1**2)*x2
    #E =  2.0*Ei*(x1**2)*x2*x2
    F =  0
    return E,F

xs = np.linspace(0,6,1000)

# e_n  = exp_approx( xs, Ri=1.4+1.6, Rcut=6.0, n=4 )
# plt.plot(xs,e_n,'-',label='e_n')

#e_n2  = exp_approx_2( xs, Ri=1.4+1.6, Rcut=6.0, n=4 )
# e_n2  = exp_approx_2( xs, Ri=1.4+1.6, Rcut=5.5, n=4 )
# plt.plot(xs,e_n2,'-',label='e_n2')



E_LJ,F_LJ  = LennardJones( xs, Ei=0.1, Ri=1.4+1.6 )
F_LJnm, x_ = numDeriv( xs, E_LJ )
plt.plot(xs,E_LJ,'-',label='E_LJ')
#plt.plot(xs,F_LJ,'-',label='F_LJ')
#plt.plot(x_,F_LJnm,':',label='Fnum_LJ')
vmin = np.nanmin(E_LJ)*1.2

E_M,F_M, e_M   = Morse( xs, Ei=0.1, Ri=1.4+1.6, b=1.6 )
# F_Mnm, x_ = numDeriv( xs, E_M )
#plt.plot(xs,e_M,'-',label='e_Morse')
plt.plot(xs,E_M,'-',label='E_Morse')
# plt.plot(xs,F_M,'-',label='F_Morse')
# plt.plot(x_,F_Mnm,':',label='Fnum_Morse')
vmin = np.min(E_M)*1.2

E_C1,F_C1, e_C1 =  MorseCut_1( xs, Ei=0.1, Ri=1.4+1.6, Rc=6.0, n=4 )
#plt.plot(xs,e_C1,'-',label='e_C1')
plt.plot(xs,E_C1,'-',label='E_C1')

E_C2,F_C2, e_C2 =  MorseCut_2( xs, Ei=0.1, Ri=1.4+1.6, Rc=5.2,  n=4 )
#plt.plot(xs,e_C2,'-',label='e_C1')
plt.plot(xs,E_C2,'-',label='E_C2')

#E_3,F_3  = MorseCut_3( xs, Ei=0.1, Ri=1.4+1.6, Rc=6.0 )
#plt.plot(xs,E_3,'-k', lw=3, label='E_3')

E_4,F_4  = MorseCut_4( xs, Ei=0.1, Ri=1.4+1.6, Rc=6.0 )
plt.plot(xs,E_4,'-k', lw=3, label='E_3')

#plt.ylim(0,2)

plt.ylim(vmin,-vmin)
#plt.xlim(0,6.0)
plt.grid()
plt.axhline(0,ls='--',c='k')

plt.legend()

plt.show()