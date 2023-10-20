#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
#import utils as uu
import BasePairUtils as bpu

keys=[
'H'  ,'e'  ,
'HH' ,'ee' ,
'HHH','eee',
'HeH','eHe',
'HHe','Hee',
]

#names = []
#for key in keys:
#    names+=name_dct[key]

fname = "/home/prokop/Desktop/CARBSIS/Paolo/correlations/binding_energy.dat"

#Es         = uu.load_dict( "binding_energy.dat" )

_,Es,names,pair_names = bpu.read_dat( fname,   ni=0, nf=1, iname=0 )

#print( "Es", Es )
#print( "names" ,names ) 
#print( "pair_names" ,pair_names ) 

Emin_dct   = bpu.find_E_min( Es )                       ; print( "Emin_dct", Emin_dct )
Emap, Exx  = bpu.makeEmaps( names, Emin_dct )           ; print( "Emap", Emap )
E_contrast = bpu.makeEcontrast( Emap )                  ; print( "E_contrast", E_contrast )

np.savetxt(  "Emap_XY.txt", Emap )
np.savetxt(  "Emap_XX.txt", Exx )


print("\n======== NEW RUN ========== \n")

#print( E_contrast )

p=(46,46)
print("pair", p, " is ", names[46],names[46], E_contrast[p[0],p[1]] )


E_max  = -20.
E_min  = -30.
dEmin  = 5.

Espan  = 5.0

#uu.findAlphabetsForRange( Emap, EbindRange=(-40.,-30.), dEmin=dEmin, nPairMax=4, E_contrast=E_contrast, verbosity=1, nMax=1000 )

Eranges = np.arange( -60, -10, 5 )
inds=[]
for E0 in Eranges:
    print( "\n################## Erange (%g,%g) dE %g" %(E0-Espan,E0+Espan,dEmin) )
    alphs, ind = bpu.findAlphabetsForRange( Emap, EbindRange=(E0-Espan,E0+Espan), dEmin=dEmin, nPairMax=4, E_contrast=E_contrast, verbosity=1, nMax=5, names=names )
    inds.append(ind)
#exit()



#Emap_s = Emap[ind1d,ind1d].copy()   ;print(Emap_s,"\n", Emap_s.shape)
#m1 = len(Emap_s)
#plt.imshow( -Emap_s, origin='lower', cmap='plasma' )

plt.figure()
Emap_       = np.reshape( Emap,       -1)
E_contrast_ = np.reshape( E_contrast, -1)
#plt.scatter( Emap_*-1., E_contrast_, c='k', s=0.5 )
plt.scatter( Emap_*-1., E_contrast_, c='k', s=1.0 )
#for ind in inds:
#    plt.plot( Emap[ind[0],ind[1]]*-1., E_contrast[ind[0],ind[1]], '.' )

plt.axhline(dEmin,ls='--',color='k')
plt.xlabel('E_Bind [kcal/mol]')
plt.ylabel('E_contrastB [kcal/mol]')

plt.xlim(0,50)   ;plt.xticks(np.arange(0, 50+1, step=5.0))
plt.ylim(0,40)   ;plt.yticks(np.arange(0, 40+1, step=5.0))

#plt.axis('equal')

plt.tight_layout()
plt.savefig("Alphabet.png", bbox_inches='tight')
plt.savefig("Alphabet.svg", bbox_inches='tight')

'''
plt.figure()
#plt.imshow( Emap )
#plt.imshow( -Emap, origin='lower', cmap='plasma', extent=(0,n-1,0,n-1) )
plt.imshow( -E_contrast, origin='lower', cmap='plasma', extent=(0,n-1,0,n-1), vmin=0 )
ax=plt.gca()
ax.set_xticks(range(n))
ax.set_yticks(range(n))
ax.set_xticklabels(names, rotation = 90)
ax.set_yticklabels(names )
plt.colorbar()
'''

plt.show()


