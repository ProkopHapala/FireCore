import numpy             as np
import matplotlib.pyplot as plt

def R8func( r2, Rpeak, Rnod ):
    '''
    This functions should go smoothly from 1 to 0 in the interval [R1,R2]
    '''
    x0 = Rnod/Rpeak
    a = (x0**2 - 1)**2
    y1 = 1 - r2/(Rpeak*Rpeak)
    y2 = a - y1**2
    return (y2/a)**2

def R8fast( r2,  R, Rnod ):
    R2   = R*R;
    R2n  = Rnod*Rnod;
    y1   = R2  - r2;
    y2   = R2n +  R2 - y1*y1
    y2  *= R2/( R2 + R2n )
    return y2*y2;

def R8func_down( r2, Rpeak, Rnod ):
    y = R8func( r2, Rpeak, Rnod )
    y[ r2<(Rpeak*Rpeak) ]=1
    y[ r2>(Rnod*Rnod)   ]=0
    return y

def R8func_up  ( r2, Rpeak, Rnod ):
    y = R8func( r2, Rpeak, Rnod )
    y[ r2>(Rpeak*Rpeak) ]=1
    y[ r2<(Rnod*Rnod)   ]=0
    return y

def R8func_peak( r2, Rpeak, R1 ):
    R2_2 = 2*Rpeak**2 - R1**2
    y  = R8func( r2, Rpeak, R1 )
    y[ r2<R2_2 ]=0
    y[ r2>R1**2 ]=0
    return y

x = np.linspace(-3,3,1000)

R =1.4
R1=1.8

#  R1**2 + R2**2 = 2*R**2   


R2 = np.sqrt( 2*R**2 - R1**2 )  #;print(R3)


y_fast = R8func     ( x**2, R, R1 )  

y      = R8func     ( x**2, R, R1 )   
y_down = R8func_down( x**2, R, R1 )   
y_up   = R8func_up  ( x**2, R, R2 )   
y_peak = R8func_peak( x**2, R, R1 )

plt.plot(x,y_fast, 'g', lw=3 )

plt.plot(x,y, 'k')
plt.plot(x,y_down, 'b', lw=2 )
plt.plot(x,y_up  , 'r', lw=2 )
plt.plot(x,y_peak, '--m', lw=1 )
plt.ylim(-1,2)
plt.xlim(0,3)
plt.grid()
#plt.axhline(0,ls='-',c='k')
plt.axhline(1,ls=':',c='k')


plt.axvline(R,ls='--',c='k')
plt.axvline(R1,ls=':',c='b')
plt.axvline(R2,ls=':',c='r')
plt.show()