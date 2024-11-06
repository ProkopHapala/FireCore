
#ifndef  constants_h
#define  constants_h

#define GOLDEN_RATIO    1.618033988749894
#define deg2rad         0.01745329251994329547
#define rad2deg        57.29577951308232286465

#define COULOMB_CONST      14.3996448915      // e^2/(4*pi*eps0) in eV*Angstrom
#define const_Bohr_Radius  0.529177210903     // Angstrom

#define const_kB         8.617333262145e-5    // eV/K

#define const_Angstrom   1.0e-10
#define const_eV         1.602176634e-19
#define const_amu        1.66053906660e-27
#define const_hbar       1.0545718176461565e-34

#define const_eV_Angstrom  1.602176634e-9

#define const_timeu      1.0180505710774743e-14   // F = m*a = m* l/t^2 => [eV/A] = [amu * A / dt^2 ] => dt^2= (A^2*amu) / (eV) = sqrt( 1.6605390666050e-27 * (1e-10)^2 /1.602176634e-19 ) = 1.0180506e-14 s
#define const_timefs     10.180505710774743       
#define const_fs_timeu   0.09822694750238489      // timeu/fs = 1/10.18050571077474 

#endif
