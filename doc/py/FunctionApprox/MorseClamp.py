import numpy             as np
import matplotlib.pyplot as plt

'''
Fmax = 12*Ei*( (Ri/Rc)^12 - (Ri/Rc)^6 )/Rc
Fmax - 12*Ei*(Ri/Rc)^12/Rc = 12*Ei*(Ri/Rc)^6 )/Rc
# solve for Rc
'''

COULOMB_CONST = 14.3996448915


def numDeriv( x, y ):
    dx = x[2:]-x[:-2]
    dy = y[2:]-y[:-2]
    x_ = x[1:-1]
    return -dy/dx, x_

#def Rc_LJ(  Fmax, R, E ):
#    Rc = np.pow( ( -2*np.sqrt(3.*E * (Fmax+3*E)*(R**12) ) + 3*E*(R**6) )/Fmax , 1./6 )

def Rc_Morse(  Fmax, E, b ):
    FE = Fmax/E
    dr = np.log( ( b + np.sqrt( b*(2*FE+b) ) )/(2*b) )/b
    return dr

def LennardJones( r, Ei=0.1, Ri=1.4 ):
    e  = (Ri/r)**6
    E  =    Ei*( e*e - 2*e )
    F  = 12*Ei*( e*e -   e )/r   # the /r makes problem here, because we cannot solve equation for Rc   
    return E,F,e

def Coulomb( r, Qij, Rc ):
    E = COULOMB_CONST*Qij/r
    F = COULOMB_CONST*Qij/(r*r)
    return E,F

def Coulomb_Clamp( r, Qij, Rc ):
    r2 = r*r + Rc*Rc
    ir = 1/np.sqrt(r2)
    E = COULOMB_CONST*Qij*ir
    F = COULOMB_CONST*Qij*ir*ir
    return E,F

def LennardJones_Rcut( r, Ei=0.1, Ri=1.4, Rc=1.0 ):
    e  = (Ri/r)**6
    E  =    Ei*( e*e - 2*e )
    F  = 12*Ei*( e*e -   e )/r
    # ---- Clamp the potential and force
    ec   =  (Ri/Rc)**6
    Fmax = 12*Ei*( ec*ec -   ec )/Rc
    Emax =    Ei*( ec*ec - 2*ec ) + Fmax*( Rc - r )
    mask = r<Rc
    E[mask] = Emax[mask] 
    F[mask] = Fmax
    return E,F, e

def Morse( r, Ei=0.1, Ri=1.4, b=1.6 ):
    e = np.exp( -b*( r - Ri ) )
    E =      Ei*( e*e - 2*e )
    F =  2*b*Ei*( e*e -   e ) 
    return E,F, e

def Morse_Rcut( r, Ei=0.1, Ri=1.4, b=1.6, Rc=1.0 ):
    e = np.exp( -b*( r - Ri ) )
    E    =      Ei*( e*e - 2*e )
    F    =  2*b*Ei*( e*e -   e ) 
    # ---- Clamp the potential and force
    ec   = np.exp( -b*( Rc - Ri ) )
    Fmax = 2*b*Ei*( ec*ec -   ec )
    Emax =     Ei*( ec*ec - 2*ec ) + Fmax*( Rc - r )
    mask = r<Rc
    E[mask] = Emax[mask] 
    F[mask] = Fmax
    return E,F, e

# def Morse_Fclamp( r, Ei=0.1, Ri=1.4, b=1.6 ):
#     expr = np.exp( -b*( r - Ri ) )
#     E    =      Ei*( expr*expr - 2*expr )
#     F    =  2*b*Ei*( expr*expr -   expr ) 
#     return E,F, expr

# ======= Main
par = {
    'H': [ 1.187, 0.0006808 ],
    'C': [ 1.908, 0.0037292 ],
    'N': [ 1.780, 0.0073719 ],
    'O': [ 1.661, 0.0091063 ],
}

el1 = 'O'
el2 = 'H'

#Eij = np.sqrt( par['H'][1]*par['H'][1] );  Rij = par['H'][0]+par['H'][0]; 
Eij = np.sqrt( par[el1][1]*par[el2][1] );  Rij = par[el1][0]+par[el2][0]; 


b   =  1.6
Rmax = 6.0
Fmax = 10.0
#dRc =  0.5
dRc = Rc_Morse(  Fmax, Eij, b )
print( "Rij=",Rij," Eij=", Eij," dRc ", dRc ," Elements ", el1," ", el2  )
xs = np.linspace(0.1,Rmax,1000)

E_LJ,F_LJ,e_LJ  = LennardJones( xs, Ei=Eij, Ri=Rij )
#F_LJnm, x_ = numDeriv( xs, E_LJ )
# plt.plot(xs,E_LJ,'-',label='E_LJ')
# plt.plot(xs,F_LJ,'-',label='F_LJ')
#plt.plot(x_,F_LJnm,':',label='Fnum_LJ')
vmin = np.nanmin(E_LJ)*1.2

E_LJC,F_LJC,e_LJC  = LennardJones_Rcut( xs, Ei=Eij, Ri=Rij, Rc=Rij-dRc )
#F_LJCnm, x_ = numDeriv( xs,  E_LJC )
# plt.plot(xs,E_LJC,'-',label='E_LJC')
# plt.plot(xs,F_LJC,'-',label='F_LJC')
#plt.plot(x_,F_LJnm,':',label='Fnum_LJ')
#vmin = np.nanmin(E_LJ)*1.2

E_M,F_M, e_M   = Morse( xs, Ei=Eij, Ri=Rij, b=b )
# # F_Mnm, x_ = numDeriv( xs, E_M )
# #plt.plot(xs,e_M,'-',label='e_Morse')
# plt.plot(xs,E_M,'-b',label='E_Morse')
# plt.plot(xs,F_M,':b',label='F_Morse')
# # plt.plot(x_,F_Mnm,':',label='Fnum_Morse')
# vmin = np.min(E_M)*1.2

E_Mc,F_Mc, e_Mc   = Morse_Rcut( xs, Ei=Eij, Ri=Rij, b=b, Rc=Rij-dRc )
# # F_Mnm, x_ = numDeriv( xs, E_M )
# #plt.plot(xs,e_M,'-',label='e_MorseCut')
plt.plot(xs,E_Mc,'-',label='E_MorseCut')
plt.plot(xs,F_Mc,'-',label='F_MorseCut')
# # plt.plot(x_,F_Mnm,':',label='Fnum_MorseCut')
# #vmin = np.min(E_M)*1.2

E_C,F_C = Coulomb( xs, Qij=-0.5*0.5, Rc=Rij )
plt.plot(xs,E_C,'-r',label='E_Coulomb')
plt.plot(xs,F_C,':r',label='F_Coulomb')
# plt.plot(xs,E_C+E_M,'-k',lw=2,label='E_tot')
# plt.plot(xs,F_C+F_M,':k',lw=2,label='F_tot')

E_Cc,F_Cc = Coulomb_Clamp( xs, Qij=-0.5*0.5, Rc=Rij*0.25 )
plt.plot(xs,E_Cc,'-m',label='E_CoulombClamp')
plt.plot(xs,F_Cc,':m',label='F_CoulombClamp')
# plt.plot(xs,E_Cc+E_M,'-k',lw=2,label='E_tot')
# plt.plot(xs,F_Cc+F_M,':k',lw=2,label='F_tot')
plt.plot(xs,E_Cc+E_Mc,'-k',lw=2,label='E_tot')
plt.plot(xs,F_Cc+F_Mc,':k',lw=2,label='F_tot')

#Fmax = 2.0
#fCut = np.power(  ( -3 + np.sqrt(3*Fmax/Eij + 9 ))/6.0 ,   1/6.0  )
#print( fCut )

# Fmax  = 2.0
# dRcut = Rc_Morse( Fmax, 0.0, Eij, b ); print( dRcut )
# #fCut = 0.85
# p = 1/2.5
# for rij in np.linspace(Rij-1.0,Rij+1.0,5):
#     RcLJ = rij-0.5  + np.power(Rij, p ) - np.power(rij, p )
#     E_LJC,F_LJC,e_LJC  = LennardJones_Rcut( xs, Ei=Eij, Ri=rij, Rc=RcLJ )
#     #plt.plot(xs,E_LJC,'-',label='E_LJC')
#     plt.plot(xs,F_LJC,'-',label='F_LJC')

#     E_MC,F_MC, e_MC = Morse_Rcut( xs, Ei=Eij, Ri=rij, b=b, Rc=rij-dRcut )
#     #plt.plot(xs,E_MC,'-',label='E_MorseCut')
#     #plt.plot(xs,F_MC,'-',label='F_MorseCut')
    
plt.axvline(Rij-dRc,ls='--',c='k')

#plt.ylim(vmin*5.0,-vmin*30)
plt.ylim(-5,5)
#plt.xlim(0,6.0)
plt.grid()
plt.axhline(0,ls='--',c='k')

plt.legend()

plt.show()