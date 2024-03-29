#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

import scipy.special as spc

'''

Note: It seems that H2 Molecule cannot be sable without varying Kinetic Energy

see:  

[1] https://link.aps.org/doi/10.1103/PhysRevLett.99.185003
Excited Electron Dynamics Modeling of Warm Dense Matter
Julius T. Su, William A. Goddard,

[2] http://aip.scitation.org/doi/10.1063/1.3272671
The dynamics of highly excited electronic systems: Applications of the electron force field
Julius T. Su, William A. Goddard

[3] http://dx.doi.org/10.1016/j.mechmat.2015.02.008
Non-adiabatic dynamics modeling framework for materials in extreme conditions
Hai Xiao, Andres Jaramillo-Botero, Patrick L. Theofanis, William A. Goddard


To check and obtain constants:

https://en.wikipedia.org/wiki/Hydrogen_atom#Bohr%E2%80%93Sommerfeld_Model
https://en.wikipedia.org/wiki/Fine-structure_constant

'''


# ==== constants in  SI Units


# see https://en.wikipedia.org/wiki/Fine-structure_constant

const_bohr = 0.5291772105638411
Hartree2eV = 27.211386245988
A2bohr     = 1/0.5291772105638411

const_hbar  = 1.054571817e-34 # [J.s]  #6.582119569e-16 # [eV/s]
const_Me    = 9.10938356e-31  # [kg]
const_e     = 1.602176620898e-19  # [Coulomb]
const_eps0  = 8.854187812813e-12 # [F.m = Coulomb/(Volt*m)]
const_eV    = 1.602176620898e-19 # [J]
const_Angstroem = 1.0e-10 

const_K     =  const_hbar**2/const_Me
const_El    =  const_e**2/(4.*np.pi*const_eps0)

const_Ry     = 0.5 * const_El**2/const_K
const_Ry_eV  = 13.6056925944
const_El_eVA = const_El/( const_e*const_Angstroem )

const_K_eVA  = (const_El_eVA**2)/(2*const_Ry_eV)

#const_K  =   const_hbar**2/const_Me   #   [ eV * A^2 ]
#const_K  = 0.1* 30.0824137226  # [eV*A^2] hbar[J.s]^2/(Me [kg])   /  (  eV[J]*A^2[m])    # (6.62607015e-34^2/9.10938356e-31)/1.602176620898e-19/10e-20
#const_Ke =  1.5*const_K

const_Ke_eVA = const_K_eVA*1.5

#print("const_El, const_El_eVA ", const_El, const_El_eVA)
#print("const_Ry  const_Ry_eV  ", const_Ry, const_Ry/const_eV)
#print("const_K,  const_K_eVA  ", const_K, const_K_eVA)
#print("const_Ke_eVA           ", const_Ke_eVA)

sqrt2 = np.sqrt(2.)

def Kinetic( s ):
    '''
    Ek = (hbar^2/me) (3./2.) 1/s^2
    '''
    return const_Ke_eVA/(s**2)

def El( r, qq, si=0, sj=0 ):
    if si>0:
        if sj>0:
            si = np.sqrt( si**2 + sj**2 )
        return const_El_eVA * (qq/r) * spc.erf( sqrt2 * r/s )
    else:
        return const_El_eVA * (qq/r)

def El_aa( r, qq ):
    return const_El_eVA * (qq/r)

def El_ae( r, qq, s ):
    return const_El_eVA * (qq/r) * spc.erf( sqrt2 * r/s )

def El_ee( r, qq, si, sj ):
    s = np.sqrt( si**2 + sj**2 )
    return const_El_eVA * (qq/r) * spc.erf( sqrt2 * r/s )

def getT( r, si, sj ):
    #print "getT si, sj ", si, sj
    #   r = r * 1.125
    #   s = s*0.9
    si2 = si**2
    sj2 = sj**2
    r2  = r**2
    #return const_K * ( 1.5*( (si2+sj2)/(si2*sj2) )   - 2.*( 3.*(si2+sj2) - 2.*r2 )/( si2 + sj2 )**2 )
    return const_K_eVA * ( 1.5*( 1./si2 + 1./sj2 )   - 2.*( 3.*(si2+sj2) - 2.*r2 )/( si2 + sj2 )**2 )

def getAmp( si, sj ):
    si2 = si**2
    sj2 = sj**2
    #return const_K_eVA * ( 1.5*( 1./si2 + 1./sj2 )   - 2.*( 3.*(si2+sj2) - 2.*0 )/( si2 + sj2 )**2 )
    #return const_K_eVA * ( 1.5*( 1./si2 + 1./sj2 )   - 2.*( 1.*(si2+sj2) )/( si2 + sj2 )**2 )
    #return const_K_eVA * 2.2*( 1.5*( 1/si2 + 1/sj2 ) - 4.9/( si2 + sj2 ) )
    #return const_K_eVA * 2.2*( 1.5*( (si2 + sj2)/(si2*sj2) ) - 4.9/( si2 + sj2 ) )
    #return const_K_eVA * 2.2*( 1.5*(si2*si2 + sj2*sj2) - 1.9*(si2*sj2)  )/((si2*sj2)*(si2+sj2))
    #return const_K_eVA * 2.2*( 1.5*(si2*si2 + sj2*sj2) - 1.9*(si2*sj2)  )/((si2*sj2)*(si2+sj2))
    return const_K_eVA * 3.3*( si2*si2 + sj2*sj2 - 1.25*(si2*sj2) )/((si2*sj2)*(si2+sj2))
    #return const_K_eVA * 3.14*( si2*si2 + sj2*sj2 - 1.25*(si2*sj2) )/((si2*sj2)*(si2+sj2))
    #return const_K_eVA * ( 1.5*( 1./si2 + 1./sj2 ) ) 
    #return const_K_eVA * ( 1.5*( 1./si2 + 1./sj2 )  - 2.*3./( si2 + sj2 ) )



def getS( r, si, sj ):
    #print "getS si, sj ", si, sj
    #   r = r * 1.125
    #   s = s*0.9
    si2 = si**2
    sj2 = sj**2
    r2  = r**2
    return ( 2.*(si*sj)/(si2+sj2) )**1.5 * np.exp(  -r2/( si2 + sj2 ) )

'''
def EPauli( r, si, sj, rho=0.2 ):
    T = getT( r, si, sj )
    S = getS( r, si, sj )
    S2 = S**2
    # ( S2*(1+S2) + (1-rho)* S2*(1-S2) )  / (1-S2*S2 )
    # ( S2+S2*S2 + (1-rho)*(S2-S2*S2) )  / (1-S2*S2 )
    # ( ( (2-rho)*S2 +rho*S2*S2 )  / (1-S2*S2 )
    return T * ( (S2/(1.-S2))   + ( 1.-rho )*(S2/(1.+S2))     )

def EPauli_pair( r, si, sj, rho=0.2 ):
    T  = getT( r, si, sj )
    S  = getS( r, si, sj )
    S2 = S**2
    return T * ( rho*S2/(1.+S2) )
'''

def EPauli( r, si, sj, anti=False, rho=0.2, kr=1.125, ks=0.9 ):
    r  = r*kr
    si = si*ks
    sj = sj*ks
    T = getT( r, si, sj )
    S = getS( r, si, sj )
    S2 = S**2
    if anti:
        return T * ( rho*S2/(1.+S2) )
    else:
        return T * ( (S2/(1.-S2))   + ( 1.-rho )*(S2/(1.+S2)) )

def DensOverlap( r, si, sj, amp=10 ):
    s2 = si**2+sj**2
    #amp *= 1.4/s2
    #amp *= 0.7/(si*sj)
    #amp *= (1/si**2 + 1/sj**2)
    #amp  *= (si**2+sj**2)/(si*sj)**2
    #amp  *= (si+sj)**2/(si*sj)**2
    #amp  *= (1+(si-sj)**2)/min(si,sj)**2
    #amp  *= 0.5*(1+4*(si-sj)**2) *( 1/si**2 + 1/sj**2 )
    a  = 2*si*sj/s2
    e1 = amp*a**3
    e2 = np.exp( -2*(r**2/s2) )
    return e1*e2

def Hatom( s ):
    Ek  = Kinetic( s )
    Eae = El_ae( 0.01, -1., s )
    #Etot = Ek+Eae
    return Ek,Eae

def H2cation( rHH, s, cr=0.5 ):
    Ek   = Kinetic( s )                  # kinetic energy of electron
    Eaa  = El_aa( rHH,  1. )             # Coulomb repulsion  nuclei_1 + nuclei_2
    Eae  = El_ae( rHH*(cr   ), -1., s )  # Coulomb attraction electron + nuclei_1
    Eae += El_ae( rHH*(1.-cr), -1., s )  # Coulomb attraction electron + nuclei_2
    return Ek, Eae, Eaa

def H2molecule( r, s, cr=0.5 ):
    Ek    = 2*Kinetic( s )                      # kinetic energy of electron_1 and electron_2
    Eaa   =   El_aa( r,  +1. )                  # Coulomb repulsion   nuclei_1 * nuclei_2
    Eae   = 2*El_ae( r*(cr   ), -1., s )        # Coulomb attraction (electron_1 * nuclei_1)   +   (electron_2 * nuclei_2)
    Eae  += 2*El_ae( r*(1.-cr), -1., s )        # Coulomb attraction (electron_1 * nuclei_2)   +   (electron_2 * nuclei_1)
    Eee   =   El_ee( r*(1.-2.*cr), +1., s, s )  # Coulomb repulsion   electron_1 * electron_2
    EPaul =   EPauli( r*(1.-2.*cr), s, s, anti=True )   # Pauli repulsion electron_1 * electron_2
    return Ek, Eae, Eaa, Eee, EPaul




def getExmin1D(Es,xs):
    imin = np.argmin( Es )
    return imin,Es[imin],xs[imin]
    #print("H-atom Rmin Emin(Ek,Eel) ", ys[imin], Etot[imin], Ek[imin], Eae[imin]) 

    def S_ij(r, si, sj):
        # energy
        S_ij = ((2./((si/sj)+(sj/si)))**(3./2.)) * \
            np.exp((-1.*(r**2))/(si**2+sj**2))
        return S_ij


def E_up_up__( r, si, sj, rho=-0.2):
    # from sympy calculations
    print( "DEBUG rho",rho, "r",r, "si",si, "sj",sj )
    E_up_up = (8.0*(-rho + 1.0)*(1/(si/sj + sj/si))**3.0*np.exp(-2.0*r**2/(si**2 + sj**2))/(8.0*(1/(si/sj + sj/si))**3.0*np.exp(-2.0*r**2/(si**2 + sj**2)) + 1.0) + 8.0*(1/(si/sj + sj/si))**3.0 *
                np.exp(-2.0*r**2/(si**2 + sj**2))/(-8.0*(1/(si/sj + sj/si))**3.0*np.exp(-2.0*r**2/(si**2 + sj**2)) + 1.0))*(-(-4.0*r**2 + 6.0*si**2 + 6.0*sj**2)/(si**2 + sj**2)**2 + 1.5/sj**2 + 1.5/si**2)
    print( "DEBUG E_up_up [Ht]", E_up_up            , "rho",rho, "r[au]",r,           "si[au]",si,           "sj[au]",sj )
    print( "DEBUG E_up_up [eV]", E_up_up*Hartree2eV , "rho",rho, "r[A]",r*const_bohr, "si[A]",si*const_bohr, "sj[A]",sj*const_bohr )
    # from sympy calculations
    dE_up_updr = -1*1.125*r*(-8.0*si**2*sj**2*(si*sj/(si**2 + sj**2))**3.0*(si**2 + sj**2)*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(64.0*(si*sj/(si**2 + sj**2))**3.0 + 8.0*(rho - 1.0)*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2))) + 8.0*np.exp(2.0*r**2/(si**2 + sj**2))) + (si**2*sj**2*(4.0*r**2 - 6.0*si**2 - 6.0*sj**2) + 1.5*si**2*(si**2 + sj**2)**2 + 1.5*sj**2*(si**2 + sj**2)**2)*(32.0*(si*sj/(si**2 + sj**2))**3.0*(rho - 1.0)*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2*(8.0*(si*sj/(si**2 + sj**2))
                                ** 3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2))) + (si*sj/(si**2 + sj**2))**3.0*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2*(256.0*(si*sj/(si**2 + sj**2))**3.0 - 32.0*np.exp(2.0*r**2/(si**2 + sj**2))) + (si*sj/(si**2 + sj**2))**6.0*(-256.0*rho + 256.0)*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2 - 256.0*(si*sj/(si**2 + sj**2))**6.0*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2))/(si**2*sj**2*(si**2 + sj**2)**3*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2)
    if np.isnan(dE_up_updr) == True:
        dE_up_updr = 0

    # from sympy calculations
    dE_up_upds1 = -1*0.9*(sj**2*(si*sj/(si**2 + sj**2))**3.0*(si**2 + sj**2)*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(12.0*si**4*(si**2 + sj**2) + si**4*(16.0*r**2 - 24.0*si**2 - 24.0*sj**2) + 3.0*(si**2 + sj**2)**3)*(64.0*(si*sj/(si**2 + sj**2))**3.0 + 8.0*(rho - 1.0)*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2))) + 8.0*np.exp(2.0*r**2/(si**2 + sj**2))) + (si**2*sj**2*(4.0*r**2 - 6.0*si**2 - 6.0*sj**2) + 1.5*si**2*(si**2 + sj**2)**2 + 1.5*sj**2*(si**2 + sj**2)**2)*(32.0*r**2*si**2*(si*sj/(si**2 + sj**2))**3.0*(-rho + 1)*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2))) - 32.0*r**2*si**2*(si*sj/(si**2 + sj**2))**3.0*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))
                            ** 2 + 24.0*(si*sj/(si**2 + sj**2))**3.0*(rho - 1.0)*(si**2 - sj**2)*(si**2 + sj**2)*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2))) + 24.0*(si*sj/(si**2 + sj**2))**3.0*(si**2 - sj**2)*(si**2 + sj**2)*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2 + 8.0*(si*sj/(si**2 + sj**2))**6.0*(rho - 1.0)*(32.0*r**2*si**2 + 24.0*(-si**2 + sj**2)*(si**2 + sj**2))*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2 + (si*sj/(si**2 + sj**2))**6.0*(256.0*r**2*si**2 + 192.0*(-si**2 + sj**2)*(si**2 + sj**2))*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2))/(si**3*sj**2*(si**2 + sj**2)**4*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2)
    if np.isnan(dE_up_upds1) == True:
        dE_up_upds1 = 0

    # from sympy calculations
    dE_up_upds2 = -1*0.9*(si**2*(si*sj/(si**2 + sj**2))**3.0*(si**2 + sj**2)*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(12.0*sj**4*(si**2 + sj**2) + sj**4*(16.0*r**2 - 24.0*si**2 - 24.0*sj**2) + 3.0*(si**2 + sj**2)**3)*(64.0*(si*sj/(si**2 + sj**2))**3.0 + 8.0*(rho - 1.0)*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2))) + 8.0*np.exp(2.0*r**2/(si**2 + sj**2))) - (si**2*sj**2*(4.0*r**2 - 6.0*si**2 - 6.0*sj**2) + 1.5*si**2*(si**2 + sj**2)**2 + 1.5*sj**2*(si**2 + sj**2)**2)*(32.0*r**2*sj**2*(si*sj/(si**2 + sj**2))**3.0*(rho - 1)*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2))) + 32.0*r**2*sj**2*(si*sj/(si**2 + sj**2))**3.0*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))
                            ** 2 + 24.0*(si*sj/(si**2 + sj**2))**3.0*(rho - 1.0)*(si**2 - sj**2)*(si**2 + sj**2)*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2))) + 24.0*(si*sj/(si**2 + sj**2))**3.0*(si**2 - sj**2)*(si**2 + sj**2)*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2 - 8.0*(si*sj/(si**2 + sj**2))**6.0*(rho - 1.0)*(32.0*r**2*sj**2 + 24.0*(si**2 - sj**2)*(si**2 + sj**2))*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2 + (si*sj/(si**2 + sj**2))**6.0*(-256.0*r**2*sj**2 + 192.0*(-si**2 + sj**2)*(si**2 + sj**2))*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2))/(si**2*sj**3*(si**2 + sj**2)**4*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2)
    if np.isnan(dE_up_upds2) == True:
        dE_up_upds2 = 0

    return E_up_up, dE_up_updr, dE_up_upds1, dE_up_upds2

def pyeff_E_up_up( r, si, sj, rho=-0.2):
    if r<1.125e-7: r=1.125e-7
    # from sympy calculations
    #E_up_up = (8.0*(-rho + 1.0)*(1/(si/sj + sj/si))**3.0*np.exp(-2.0*r**2/(si**2 + sj**2))/(8.0*(1/(si/sj + sj/si))**3.0*np.exp(-2.0*r**2/(si**2 + sj**2)) + 1.0) + 8.0*(1/(si/sj + sj/si))**3.0 *
    #            np.exp(-2.0*r**2/(si**2 + sj**2))/(-8.0*(1/(si/sj + sj/si))**3.0*np.exp(-2.0*r**2/(si**2 + sj**2)) + 1.0))*(-(-4.0*r**2 + 6.0*si**2 + 6.0*sj**2)/(si**2 + sj**2)**2 + 1.5/sj**2 + 1.5/si**2)
    print( "DEBUG rho",rho, "r",r, "si",si, "sj",sj )
    E_up_up = (8.0*(-rho + 1.0)*(1/(si/sj + sj/si))**3.0*np.exp(-2.0*r**2/(si**2 + sj**2))/(8.0*(1/(si/sj + sj/si))**3.0*np.exp(-2.0*r**2/(si**2 + sj**2)) + 1.0) + 8.0*(1/(si/sj + sj/si))**3.0 *
                   np.exp(-2.0*r**2/(si**2 + sj**2))/(-8.0*(1/(si/sj + sj/si))**3.0*np.exp(-2.0*r**2/(si**2 + sj**2)) + 1.0))*(-(-4.0*r**2 + 6.0*si**2 + 6.0*sj**2)/(si**2 + sj**2)**2 + 1.5/sj**2 + 1.5/si**2)
        
    print( "DEBUG E_up_up [Ht]", E_up_up           , "rho",rho, "r[au]",r,           "si[au]",si,           "sj[au]",sj )
    print( "DEBUG E_up_up [eV]", E_up_up*Hartree2eV, "rho",rho, "r[A]",r*const_bohr, "si[A]",si*const_bohr, "sj[A]",sj*const_bohr )
    
    DT = (3./2.)*((1./(si**2))+(1./(sj**2)))-2. * (3.*((si**2)+(sj**2))-2.*(r**2))/(((si**2)+(sj**2))**2)
    S  = ((2./((si/sj)+(sj/si)))**(3./2.)) * np.exp((-1.*(r**2))/(si**2+sj**2))
    # dE_up_updr = -1*1.125*r*(-8.0*si**2*sj**2*(si*sj/(si**2 + sj**2))**3.0*(si**2 + sj**2)*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(64.0*(si*sj/(si**2 + sj**2))**3.0 + 8.0*(rho - 1.0)*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2))) + 8.0*np.exp(2.0*r**2/(si**2 + sj**2))) + (si**2*sj**2*(4.0*r**2 - 6.0*si**2 - 6.0*sj**2) + 1.5*si**2*(si**2 + sj**2)**2 + 1.5*sj**2*(si**2 + sj**2)**2)*(32.0*(si*sj/(si**2 + sj**2))**3.0*(rho - 1.0)*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2*(8.0*(si*sj/(si**2 + sj**2))
    #                             ** 3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2))) + (si*sj/(si**2 + sj**2))**3.0*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2*(256.0*(si*sj/(si**2 + sj**2))**3.0 - 32.0*np.exp(2.0*r**2/(si**2 + sj**2))) + (si*sj/(si**2 + sj**2))**6.0*(-256.0*rho + 256.0)*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2 - 256.0*(si*sj/(si**2 + sj**2))**6.0*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2))/(si**2*sj**2*(si**2 + sj**2)**3*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2)
    # if np.isnan(dE_up_updr) == True:dE_up_updr = 0

    # dE_up_upds1 = -1*0.9*(sj**2*(si*sj/(si**2 + sj**2))**3.0*(si**2 + sj**2)*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(12.0*si**4*(si**2 + sj**2) + si**4*(16.0*r**2 - 24.0*si**2 - 24.0*sj**2) + 3.0*(si**2 + sj**2)**3)*(64.0*(si*sj/(si**2 + sj**2))**3.0 + 8.0*(rho - 1.0)*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2))) + 8.0*np.exp(2.0*r**2/(si**2 + sj**2))) + (si**2*sj**2*(4.0*r**2 - 6.0*si**2 - 6.0*sj**2) + 1.5*si**2*(si**2 + sj**2)**2 + 1.5*sj**2*(si**2 + sj**2)**2)*(32.0*r**2*si**2*(si*sj/(si**2 + sj**2))**3.0*(-rho + 1)*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2))) - 32.0*r**2*si**2*(si*sj/(si**2 + sj**2))**3.0*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))
    #                         ** 2 + 24.0*(si*sj/(si**2 + sj**2))**3.0*(rho - 1.0)*(si**2 - sj**2)*(si**2 + sj**2)*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2))) + 24.0*(si*sj/(si**2 + sj**2))**3.0*(si**2 - sj**2)*(si**2 + sj**2)*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2 + 8.0*(si*sj/(si**2 + sj**2))**6.0*(rho - 1.0)*(32.0*r**2*si**2 + 24.0*(-si**2 + sj**2)*(si**2 + sj**2))*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2 + (si*sj/(si**2 + sj**2))**6.0*(256.0*r**2*si**2 + 192.0*(-si**2 + sj**2)*(si**2 + sj**2))*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2))/(si**3*sj**2*(si**2 + sj**2)**4*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2)
    # if np.isnan(dE_up_upds1) == True:dE_up_upds1 = 0

    # dE_up_upds2 = -1*0.9*(si**2*(si*sj/(si**2 + sj**2))**3.0*(si**2 + sj**2)*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(12.0*sj**4*(si**2 + sj**2) + sj**4*(16.0*r**2 - 24.0*si**2 - 24.0*sj**2) + 3.0*(si**2 + sj**2)**3)*(64.0*(si*sj/(si**2 + sj**2))**3.0 + 8.0*(rho - 1.0)*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2))) + 8.0*np.exp(2.0*r**2/(si**2 + sj**2))) - (si**2*sj**2*(4.0*r**2 - 6.0*si**2 - 6.0*sj**2) + 1.5*si**2*(si**2 + sj**2)**2 + 1.5*sj**2*(si**2 + sj**2)**2)*(32.0*r**2*sj**2*(si*sj/(si**2 + sj**2))**3.0*(rho - 1)*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2))) + 32.0*r**2*sj**2*(si*sj/(si**2 + sj**2))**3.0*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))
    #                         ** 2 + 24.0*(si*sj/(si**2 + sj**2))**3.0*(rho - 1.0)*(si**2 - sj**2)*(si**2 + sj**2)*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2))) + 24.0*(si*sj/(si**2 + sj**2))**3.0*(si**2 - sj**2)*(si**2 + sj**2)*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2 - 8.0*(si*sj/(si**2 + sj**2))**6.0*(rho - 1.0)*(32.0*r**2*sj**2 + 24.0*(si**2 - sj**2)*(si**2 + sj**2))*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2 + (si*sj/(si**2 + sj**2))**6.0*(-256.0*r**2*sj**2 + 192.0*(-si**2 + sj**2)*(si**2 + sj**2))*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2))/(si**2*sj**3*(si**2 + sj**2)**4*(8.0*(si*sj/(si**2 + sj**2))**3.0 - 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2)
    # if np.isnan(dE_up_upds2) == True: dE_up_upds2 = 0

    return E_up_up, DT, S # dE_up_updr, dE_up_upds1, dE_up_upds2

def pyeff_E_up_down( r, si, sj, rho=-0.2):
    # from sympy calculation
    E_up_down = ( -8.0*rho*(-(-4.0*r**2 + 6.0*si**2 + 6.0*sj**2)/(si**2 + sj**2)**2 + 1.5/sj**2 + 1.5/si**2)*(1/(si/sj + sj/si))**3.0
                 *np.exp(-2.0*r**2/(si**2 + sj**2))/(8.0*(1/(si/sj + sj/si))**3.0*np.exp(-2.0*r**2/(si**2 + sj**2)) + 1.0) )
    
    # dE_up_downdr = -1*1.125*r*rho*(-64.0*si**2*sj**2*(si*sj/(si**2 + sj**2))**3.0*(si**2 + sj**2)*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2))) + 32.0*(si*sj/(si**2 + sj**2))**3.0*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(si**2*sj**2*(4.0*r**2 - 6.0*si**2 - 6.0*sj**2) + 1.5*si**2*(
    #     si**2 + sj**2)**2 + 1.5*sj**2*(si**2 + sj**2)**2) + (si*sj/(si**2 + sj**2))**6.0*(si**2*sj**2*(-1024.0*r**2 + 1536.0*si**2 + 1536.0*sj**2) - 384.0*si**2*(si**2 + sj**2)**2 - 384.0*sj**2*(si**2 + sj**2)**2))/(si**2*sj**2*(si**2 + sj**2)**3*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2)
    # if np.isnan(dE_up_downdr) == True: dE_up_downdr = 0

    # dE_up_downds1 = -1*0.9*rho*(-32.0*r**2*si**2*(si*sj/(si**2 + sj**2))**3.0*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(si**2*sj**2*(4.0*r**2 - 6.0*si**2 - 6.0*sj**2) + 1.5*si**2*(si**2 + sj**2)**2 + 1.5*sj**2*(si**2 + sj**2)**2) + 8.0*sj**2*(si*sj/(si**2 + sj**2))**3.0*(si**2 + sj**2)*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(12.0*si**4*(si**2 + sj**2) + si**4*(16.0*r**2 - 24.0*si**2 - 24.0*sj**2) + 3.0*(si**2 + sj**2)**3) + 24.0*(si*sj/(si**2 + sj**2))**3.0*(
    #     si**2 - sj**2)*(si**2 + sj**2)*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(si**2*sj**2*(4.0*r**2 - 6.0*si**2 - 6.0*sj**2) + 1.5*si**2*(si**2 + sj**2)**2 + 1.5*sj**2*(si**2 + sj**2)**2) + 8.0*(si*sj/(si**2 + sj**2))**6.0*(32.0*r**2*si**2 + 24.0*(-si**2 + sj**2)*(si**2 + sj**2))*(si**2*sj**2*(4.0*r**2 - 6.0*si**2 - 6.0*sj**2) + 1.5*si**2*(si**2 + sj**2)**2 + 1.5*sj**2*(si**2 + sj**2)**2))/(si**3*sj**2*(si**2 + sj**2)**4*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2)
    # if np.isnan(dE_up_downds1) == True:  dE_up_downds1 = 0

    # dE_up_downds2 = -1.0*0.9*rho*(-32.0*r**2*sj**2.0*(si*sj/(si**2.0 + sj**2.0))**3.0*(8.0*(si*sj/(si**2.0 + sj**2.0))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(si**2*sj**2*(4.0*r**2 - 6.0*si**2 - 6.0*sj**2) + 1.5*si**2*(si**2 + sj**2)**2 + 1.5*sj**2*(si**2 + sj**2)**2) + 8.0*si**2*(si*sj/(si**2 + sj**2))**3.0*(si**2 + sj**2)*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(12.0*sj**4*(si**2 + sj**2) + sj**4*(16.0*r**2 - 24.0*si**2 - 24.0*sj**2) + 3.0*(si**2 + sj**2)**3) - 24.0*(si*sj/(si**2 + sj**2))**3.0*(
    #     si**2 - sj**2)*(si**2 + sj**2)*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))*(si**2*sj**2*(4.0*r**2 - 6.0*si**2 - 6.0*sj**2) + 1.5*si**2*(si**2 + sj**2)**2 + 1.5*sj**2*(si**2 + sj**2)**2) + 8.0*(si*sj/(si**2 + sj**2))**6.0*(32.0*r**2*sj**2 + 24.0*(si**2 - sj**2)*(si**2 + sj**2))*(si**2*sj**2*(4.0*r**2 - 6.0*si**2 - 6.0*sj**2) + 1.5*si**2*(si**2 + sj**2)**2 + 1.5*sj**2*(si**2 + sj**2)**2))/(si**2*sj**3*(si**2 + sj**2)**4*(8.0*(si*sj/(si**2 + sj**2))**3.0 + 1.0*np.exp(2.0*r**2/(si**2 + sj**2)))**2)
    # if np.isnan(dE_up_downds2) == True: dE_up_downds2 = 0

    return E_up_down #, dE_up_downdr, dE_up_downds1, dE_up_downds2
    
def pyeff_E_up_up_( r, si, sj, rho=-0.2):
    #print("pyeff_E_up_up_ !!!!!!!!!!!!!!!!!!!!")
    si2  = si**2
    sj2  = sj**2
    sij  = si*sj
    si2sj2 = si2 + sj2
    #denom_sij  = 1./(si/sj + sj/si)
    denom_sij  = sij/si2sj2
    denom_sij3 = denom_sij**3

    expr  = np.exp(-1.0*r**2/si2sj2)
    expr2 = expr**2

    expr_denom3 = denom_sij3*expr2

    #DT = (3./2.)*( (1./(si**2))+(1./(sj**2) ))-2. * (3.*((si**2)+(sj**2))-2.*(r**2))/(((si**2)+(sj**2))**2)
    #DT = (3./2.)*( sij2/(si*sj)**2  )-2. * ( 3.*( sij2 )-2.*(r**2) )/(sij2)**2
    DT = 1.5*( si2sj2/sij**2  )   -   ( 6.*si2sj2-4.*r**2 )/(si2sj2)**2
    #DT = ( (3./2.)*sij2  -2. * ( 3.*( sij2 )-2.*(r**2) )      )/(sij2)**2
    S  = ((2.*denom_sij)**(3./2.)) * expr
    S2 = 4*expr_denom3
    
    #E_up_up = 8*expr_denom3 * (   (-rho + 1.0)/(  8.0*expr_denom3 + 1.0 )   -  1/(  8.0*expr_denom3 - 1.0 ) )*(-(-4.0*r**2 + 6.0*sij2 )/sij2**2 + 1.5/sj**2 + 1.5/si**2)
    #E_up_up = 8*expr_denom3 *  ( -rho*(8.0*expr_denom3) + (rho - 2)  )/( (8.0*expr_denom3)**2 - 1  )*(-(-4.0*r**2 + 6.0*sij2 )/sij2**2 + 1.5/sj**2 + 1.5/si**2)
    #E_up_up = 8*expr_denom3 *  ( -rho*(8.0*expr_denom3) + (rho - 2)  )/( (8.0*expr_denom3)**2 - 1  )*(-(-4*r**2 + 6*sij2 )/sij2**2 + 1.5*sij2/(sj*si)**2 )

    #E_up_up = 8*expr_denom3 *  ( -rho*(8.0*expr_denom3) + (rho - 2)  )/( (8.0*expr_denom3)**2 - 1  )*DT

    E_up_up = 2*S2 *  ( -rho*( 2*S2 ) + (rho - 2)  )/( (2*S2)**2 - 1  )*DT

    return E_up_up, DT, S #, dE_up_updr, dE_up_upds1, dE_up_upds2

    E_up_up = ( (  8.0*(-rho + 1.0)* (1/(si/sj + sj/si))**3.0* np.exp(-2.0*r**2/(si**2 + sj**2))/( 8.0*(1/(si/sj + sj/si))**3.0*np.exp(-2.0*r**2/(si**2 + sj**2)) + 1.0  ) 
                 + 8.0             * (1/(si/sj + sj/si))**3.0* np.exp(-2.0*r**2/(si**2 + sj**2))/(-8.0*(1/(si/sj + sj/si))**3.0*np.exp(-2.0*r**2/(si**2 + sj**2)) + 1.0  )   )
    *(-(-4.0*r**2 + 6.0*si**2 + 6.0*sj**2)/(si**2 + sj**2)**2 + 1.5/sj**2 + 1.5/si**2) )

def pyeff_E_up_down_( r, si, sj, rho=-0.2 ):
    si2  = si**2
    sj2  = sj**2
    sij2 = si2 + sj2
    denom_sij  = 1./(si/sj + sj/si)
    denom_sij3 = denom_sij**3
    expr  = np.exp(-1.0*r**2/sij2)
    expr2 = expr**2
    expr_denom3 = denom_sij3*expr2

    S  = ((2.*denom_sij)**(3./2.)) * expr
    S2 = 4*expr_denom3
    DT = 1.5*( sij2/(si*sj)**2  )   -   ( 6.*sij2-4.*r**2 )/(sij2)**2

    #E_up_down = ( -8*rho*(-(-4*r**2 + 6*si**2 + 6*sj**2)/(si**2 + sj**2)**2 + 1.5/sj**2 + 1.5/si**2)*(1/(si/sj + sj/si))**3
    #             *np.exp(-2.0*r**2/(si**2 + sj**2))/(8.0*(1/(si/sj + sj/si))**3.0*np.exp(-2.0*r**2/(si**2 + sj**2)) + 1.0) )
    
    #E_up_down = ( 8*rho*(   (-4*r**2 + 6*sij2)/sij2**2 - 1.5/sj**2 + 1.5/si**2 )*  expr_denom3 / (8.0*expr_denom3  + 1.0) )

    E_up_down = ( -2*rho*DT*S2/ (2*S2  + 1) )

    return E_up_down #, dE_up_downdr, dE_up_downds1, dE_up_downds2


def pyeff_EPaul_(r, si, sj, rho=-0.2):
    #print("pyeff_E_up_up_ !!!!!!!!!!!!!!!!!!!!")
    si2    = si**2
    sj2    = sj**2
    sij    = si*sj
    si2sj2 = si2 + sj2
    denom_sij  = sij/si2sj2
    denom_sij3 = denom_sij**3

    expr  = np.exp(-r**2/si2sj2)
    expr2 = expr**2

    expr_denom3 = denom_sij3*expr2

    DT = 1.5*( si2sj2/sij**2  )   -   ( 6.*si2sj2-4.*r**2 )/(si2sj2)**2
    S  = ((2.*denom_sij)**(3./2.)) * expr
    S2 = 4*expr_denom3
    
    E_up_up   = 2*S2 *  ( -rho*( 2*S2 ) + (rho - 2)  )/( (2*S2)**2 - 1  )*DT
    E_up_down = ( -2*rho*DT*S2/ (2*S2  + 1) )

    return E_up_up, E_up_down, DT, S

def pyeff_EPaul( r, si, sj, rho=-0.2, kr=1.125, ks=0.9, bPrint=True ):

    r *=A2bohr*kr
    si*=A2bohr*ks
    sj*=A2bohr*ks
    r2     = r*r + 1e-8
    
    #print("pyeff_E_up_up_ !!!!!!!!!!!!!!!!!!!!")
    si2    = si*si
    sj2    = sj*sj
    sij    = si*sj
    si2sj2 = si2 + sj2
    invsi2sj2  =  1/si2sj2
    denom_sij  = sij*invsi2sj2
    denom_sij3 = denom_sij*denom_sij*denom_sij

    expr       = np.exp( -r2 * invsi2sj2 )
    expr2      = expr**2

    DT = 1.5*( si2sj2/(sij*sij)  )   -   ( 6.*si2sj2-4.*r2 )*invsi2sj2*invsi2sj2
    S  = ((2.*denom_sij)**(3./2.)) * expr
    #S2 = 4*denom_sij3*expr2
    S22 = 8*denom_sij3*expr2
    
    Euu =       DT*S22*( -rho*S22 + (rho - 2)  )/( S22*S22 - 1 )
    Eud =  -rho*DT*S22                          /( S22     + 1 )

    Euu *= Hartree2eV
    Eud *= Hartree2eV
    if bPrint:
        print( "DEBUG pyeff_EPaul(): Euu",Euu, "Eud",Eud, "DT",DT, "S22",S22, "rho",rho, "r",np.sqrt(r2), "si",si, "sj",sj )
    return Euu, Eud, DT, S


# ================== Composed functions

def run_ee_onsite():
    # ============= e-e onsite 
    r0 = 0.01
    ss = np.arange( 0.25, 5.0, 0.1 )
    rho=0.2; kr=1.125; ks=0.9
    r_ = r0*kr
    s_ = ss*ks
    T = getT( r_, s_, s_ )
    S = getS( r_, s_, s_ )
    S2 = S**2
    EPminus  =  T * ( rho*S2/(1.+S2) )
    EPplus   =  T * ( (S2/(1.-S2))  + ( 1.-rho )*(S2/(1.+S2)) )
    plt.figure()
    plt.title( 'Onsite (R= %g [A])' %r0 )
    plt.xlabel('sigma[A]')
    plt.plot( ss, S,   ':', label="S" )
    plt.plot( ss, T,   ':', label="T" )
    plt.plot( ss, EPplus, 'b', label="EP+" )
    plt.plot( ss, EPminus,'c', label="EP-" )
    plt.legend()
    plt.grid()
    #plt.show(); exit()


def run_ee_offsite():
    # ============= e-e
    rs = np.arange( 0.1, 6.0, 0.05 )
    #ss = [0.5, 1.0, 1.5 ]
    ss = [0.25, 1.0, 2.5 ]
    rho=0.2; kr=1.125; ks=0.9
    plt.figure(figsize=(13,10))
    for i,si in enumerate(ss):
        for j,sj in enumerate(ss):
            Eee = El_ee( rs, +1., si, sj )
            r_ = rs*kr
            #s_ = s*ks
            T = getT( r_, si*ks, sj*ks )
            S = getS( r_, si*ks, sj*ks )
            S2 = S**2
            EPminus  =  T * ( rho*S2/(1.+S2) )
            EPplus   =  T * ( (S2/(1.-S2))  + ( 1.-rho )*(S2/(1.+S2)) )

            #amp  = 10*(1+(si-sj)**2)/min(si,sj)**2
            #amp  = 10/min(si,sj)**2
            #amp  = 10*(1+0.6*abs(si-sj))/min(si,sj)**2
            #amp  = 10*(si/sj+sj/si)
            #amp  = 10
            #amp  = T*1.8
            amp   = getAmp( si, sj )

            EPdens   =   DensOverlap( rs, si, sj, amp=amp )
            plt.subplot(3,3,3*j+i+1)

            #plt.plot( rs, S,   ':', label="S" )
            #plt.plot( xs, T,   ':', label="T" )
            #plt.plot( rs, Eee ,   'r', label="Eee" )
            plt.plot( rs, EPplus, 'b', label="EP+" )
            #plt.plot( rs, EPminus,'c', label="EP-" )
            plt.plot( rs, EPdens,  'm', label="EPdens"  )
            plt.title( 'sigma (%g,%g)' %(si,sj) )
            plt.legend()
            plt.grid()
            #plt.plot( ys, Etot, 'k', label="Etot" )

def run_Hatom():
    # ============= H-atom
    '''
    from http://aip.scitation.org/doi/10.1063/1.3272671
    s = (3/2)*sqrt(pi/2)/Z =  1.87997120597 /Z   [bohr   ] =   0.9948379191595168 [A]
    E = -(4/(3pi))Z^2      = -0.42441318157 *Z^2 [Hartree] = -11.5488764245607    [eV]
    '''
    #Ek  = Kinetic( ys )
    #Eae = El_ae( 0.01, -1., ys )
    Ek,Eae = Hatom( ys )
    Etot = Ek+Eae
    plt.figure()
    plt.plot( ys, Ek ,  'r', label="Ek" )
    plt.plot( ys, Eae,  'b', label="Eae" )
    plt.plot( ys, Etot, 'k', label="Etot" )
    imin = np.argmin( Etot )
    print("H-atom Rmin Emin(Ek,Eel) ", ys[imin], Etot[imin], Ek[imin], Eae[imin]) 

    EHatom = Etot[imin]
    plt.legend()
    plt.grid()


def run_H2_cation():
    # ============= H2-cation
    Xs,Ys = np.meshgrid( xs,ys )
    Ek, Eae, Eaa = H2cation( Xs, Ys, cr=0.5 )
    Etot = Ek + Eaa + Eae
    #Emin = Etot.min()
    imin = np.unravel_index( np.argmin(Etot), Etot.shape )
    Emin = Etot[imin]
    Rmin = xs[imin[0]]
    Smin = ys[imin[1]]
    print("H2cation Rmin, Smin Emin Ebond ", Rmin, Smin, Emin, Emin-EHatom)
    vmin=-20.0 # [eV]
    vmax=-vmin
    plt.figure(figsize=(20,5))
    plt.subplot(1,4,1); plt.imshow( Etot, origin='image', extent=extent, vmin=vmin,vmax=vmax ) ;plt.title('Etot')
    plt.subplot(1,4,2); plt.imshow( Ek  , origin='image', extent=extent, vmin=vmin,vmax=vmax ) ;plt.title('Ek'  )
    plt.subplot(1,4,3); plt.imshow( Eaa , origin='image', extent=extent, vmin=vmin,vmax=vmax ) ;plt.title('Eaa' )
    plt.subplot(1,4,4); plt.imshow( Eae , origin='image', extent=extent, vmin=vmin,vmax=vmax ) ;plt.title('Eel' )
    #plt.subplot(1,4,2); plt.imshow( Ek  , origin='image', extent=extent ) ;plt.colorbar()  ;plt.title('Ek'  )
    #plt.subplot(1,4,3); plt.imshow( Eaa , origin='image', extent=extent ) ;plt.colorbar()  ;plt.title('Eaa' )
    #plt.subplot(1,4,4); plt.imshow( Eae , origin='image', extent=extent ) ;plt.colorbar()  ;plt.title('Eel' )


def run_H2_molecule():
    # article seem to be wrong 
    #   original: bond_lenght(H-H): 1.47 [A] ( vs 1.4 [A] exact SE solution, 1.38 [A] HF solution )
    #   they probably means bohr
    '''
    from http://aip.scitation.org/doi/10.1063/1.3272671
    s = 1.77 [bohr] =  0.9366437 [A]
    bond_lenght(H-H): 1.47 [bohr] = 0.7778905 [A] ( vs 1.4 [bohr] exact SE solution, 1.38 [bohr] HF solution )
    binding energy: 67 [kcal/mol] = 2.91 [eV]  (vs. 109 [kcal/mol] exact 86 [kcal/mol] HF )
        => Etot = 2*-11.5488764245607 - 2.91 = -26.0077528491 eV
    '''
    # ============= H2-molecule
    Ek, Eae, Eaa, Eee, EPaul = H2molecule( Xs, Ys, cr=0.49 )
    Etot = Ek + Eae + Eaa + Eee + EPaul
    #Emin = Etot.min()
    imin = np.unravel_index( np.argmin(Etot), Etot.shape )
    Emin = Etot[imin]
    Rmin = xs[imin[0]]
    Smin = ys[imin[1]]
    print("H2molecule Rmin, Smin Emin Ebond ", Rmin, Smin, Emin, Emin - 2*EHatom)
    vmin=-50.0 # [eV]
    vmax= 0.0 # [eV]
    plt.figure( figsize=(18,3) )
    plt.subplot(1,6,1); plt.imshow( Etot, origin='image', extent=extent, vmin=vmin,vmax=vmax ) ;plt.colorbar()  ;plt.title('Etot')
    #plt.subplot(1,6,2); plt.imshow( Ek  , origin='image', extent=extent, vmin=vmin,vmax=vmax ) ;plt.colorbar()  ;plt.title('Ek'  )
    #plt.subplot(1,6,3); plt.imshow( Eaa , origin='image', extent=extent, vmin=vmin,vmax=vmax ) ;plt.colorbar()  ;plt.title('Eaa' )
    #plt.subplot(1,6,4); plt.imshow( Eae , origin='image', extent=extent, vmin=vmin,vmax=vmax ) ;plt.colorbar()  ;plt.title('Eea' )
    #plt.subplot(1,6,5); plt.imshow( Eee , origin='image', extent=extent, vmin=vmin,vmax=vmax ) ;plt.colorbar()  ;plt.title('Eee' )
    #plt.subplot(1,6,6); plt.imshow( EPaul, origin='image', extent=extent, vmin=vmin,vmax=vmax ) ;plt.colorbar()  ;plt.title('EPaul')
    plt.subplot(1,6,2); plt.imshow( Ek  , origin='image', extent=extent  ) ;plt.colorbar()  ;plt.title('Ek'  )
    plt.subplot(1,6,3); plt.imshow( Eaa , origin='image', extent=extent ) ;plt.colorbar()  ;plt.title('Eaa' )
    plt.subplot(1,6,4); plt.imshow( Eae , origin='image', extent=extent ) ;plt.colorbar()  ;plt.title('Eea' )
    plt.subplot(1,6,5); plt.imshow( Eee , origin='image', extent=extent  ) ;plt.colorbar()  ;plt.title('Eee' )
    plt.subplot(1,6,6); plt.imshow( EPaul, origin='image', extent=extent ) ;plt.colorbar()  ;plt.title('EPaul')
    plt.show()


if __name__ == "__main__":

    extent=( 0.5,8.0,  0.5,4.5 )
    xs = np.arange( extent[0], extent[1], 0.05 )
    ys = np.arange( extent[2], extent[3], 0.1  )    
    
    #run_ee_onsite()
    #run_ee_offsite()
    run_Hatom()
    #run_H2_cation()
    #run_H2_molecule()
    plt.show()





