#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import BasePairUtils as bpu

#=== 1Hb
name_dct = {
'H' :[
'H-h_1',
'H-h_2',
'H-hh',
'H-p',
],
#--- e
'e': [
'N-hh',
'N-h',
'O-h',
'O-p',
],
#=== 2Hb
#--- HH
'HH' :[
'HH-h_1',
'HH-h_2',
'HH-hh-p',
'HH-hh',
'HH-h-p',
'HH-hp',
'HH-p_1',
'HH-p_2',
'HH-pp',
],
#--- ee
'ee' :[
'NN-hh',
'NN-hp',
'NN-pp',
'NO-h-p',
'NO-h',
'NO-p',
'OO-h',
],
'He' :[
'HN-hh',
'HN-hp_1',
'HN-hp_2',
'HN-h-p',
'HN-h',
'HN-pp',
'HO-h',
],
#=== 3Hb


'HHH': [
'HHH-hh',
'HHH-h-p',
'HHH-hp',
'HHH-h',
'HHH-p',
],
'eee':[
'NNN-hhh',
'NNO-hh_1',
'NNO-hh_2',
'NNO-hp',
'ONO-p',
],
'HeH' :[
'HNH-hh',
'HNH-h-p',
'HNH-hp',
'HNH-h',
'HNH-p',
],
'eHe':[
'NHO-hh',
'NHO-hp',
'OHO-h_1',
'OHO-h_2',
'OHO-h-p',
'OHO-p',
],
'HHe':[
'HHN-hh',
'HHN-hp',
'HHO-hh',
'HHO-h-p_1',
'HHO-h-p_2',
'HHO-hp',
'HHO-h',
'HHO-p',
],
'Hee': [
'HNN-hh',
'HNN-hp',
'HNO-hh',
'HNO-h-p',
'HNO-hp',
'HNO-h',
'HNO-p',
],
}

keys=[
'H'  ,'e'  ,
'HH' ,'ee' ,
'HHH','eee',
'HeH','eHe',
'HHe','Hee',
]

names = []
for key in keys:
    names+=name_dct[key]


    
Es         = bpu.load_dict( "/home/prokop/Desktop/CARBSIS/Paolo/correlations/binding_energy.dat" )
Emin_dct   = bpu.find_E_min( Es )
Emap, Exx  = bpu.makeEmaps( names, Emin_dct )
E_contrast = bpu.makeEcontrast( Emap )

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

#bpu.findAlphabetsForRange( Emap, EbindRange=(-40.,-30.), dEmin=dEmin, nPairMax=4, E_contrast=E_contrast, verbosity=1, nMax=1000 )

# find 'HH-hh-p' in names
i_HH_hh_p =  names.index('HH-hh-p'); print( "ind1d", i_HH_hh_p, "names[ind1d]", names[i_HH_hh_p] )   #  ;exit()

L4s = []

contained_groups = {}
conainded_pairs  = {}


Eranges = np.arange( -60, -10, 5 )
inds=[]
for E0 in Eranges:
    print( "\n################## Erange (%g,%g) dE %g" %(E0-Espan,E0+Espan,dEmin) )
    #alphs, ind = bpu.findAlphabetsForRange( Emap, EbindRange=(E0-Espan,E0+Espan), dEmin=dEmin, nPairMax=4, E_contrast=E_contrast, verbosity=1, names=names )
    levels, ind = bpu.findAlphabetsForRange( Emap, EbindRange=(E0-Espan,E0+Espan), dEmin=dEmin, nPairMax=4, E_contrast=E_contrast, verbosity=1, names=names, nMax2=0 )

    l4  = None
    if len(levels)>1:
        l4 = levels[1]
        #print( "!!!!!!!!! levels[1]", l4 )
        for alph in l4:
            # count how many times each group is contained in the alphabet
            for p in alph:
                conainded_pairs[p] = conainded_pairs.get(p,0)+1
                n1,n2 = p
                contained_groups[n1] = contained_groups.get(n1,0)+1
                contained_groups[n2] = contained_groups.get(n2,0)+1

            # check if 'HH-hh-p' is in the alphabet
            if any(i_HH_hh_p in pair for pair in alph):
                e = [ ( (names[n1],names[n2]),Emap[n1,n2] ) for n1,n2 in alph ]
                print( "HH-hh-p found in ", e )
    L4s .append(l4)
    inds.append(ind)
#exit()


# print sorted conainded_pairs
print( "##### conainded_pairs ", len(conainded_pairs) )
conainded_pairs_sorted = sorted(conainded_pairs.items(), key=lambda kv: kv[1], reverse=True)
for k,v in conainded_pairs_sorted:
    print( "%3i" %v, bpu.inds2names( k, names ) )
# print sorted contained_groups
print( "##### contained_groups ", len(contained_groups) )
contained_groups_sorted = sorted(contained_groups.items(), key=lambda kv: kv[1], reverse=True)
for k,v in contained_groups_sorted:
    print( "%3i" %v, names[k] )
#for k,v in contained_groups.items():
#    print( "contained_groups ", names[k],v )

#Emap_s = Emap[ind1d,ind1d].copy()   ;print(Emap_s,"\n", Emap_s.shape)
#m1 = len(Emap_s)
#plt.imshow( -Emap_s, origin='lower', cmap='plasma' )

plt.figure()
Emap_       = np.reshape( Emap,       -1)
E_contrast_ = np.reshape( E_contrast, -1)
#plt.scatter( Emap_*-1., E_contrast_, c='k', s=0.5 )
plt.scatter( Emap_*-1., E_contrast_, c='k', s=1.0 )

#colors = [ 'b', 'r', 'g', 'b', 'c', 'k', 'k', 'k', 'k', 'k', 'k' ]
colors = [ 'k', 'k', 'k','k','k','k', 'c', 'b', 'g', 'r' ]

for i,E0 in enumerate( Eranges):
    # plt lines between indexes stored in l4
    l4 = L4s[i]
    if l4 is not None:
        l4=list(l4)
        print( "E0 ", E0, "len(l4) ", len(l4) );
        for p1,p2 in l4:
            #print( "p1,p2", p1,p2 )
            c = colors[i]
            plt.plot( [Emap[p1[0],p1[1]]*-1., Emap[p2[0],p2[1]]*-1.], [E_contrast[p1[0],p1[1]], E_contrast[p2[0],p2[1]]], 'o-', ms=3, color=c, lw=0.5 )


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


