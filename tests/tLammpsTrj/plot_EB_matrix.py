import sys
import os
import re
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import plotUtils   as plu
import BasePairUtils as bpu


# ============ Setup

mols_gs=[
["H-p","H-h_1","H-h_2","H-hh",],
["N-h","N-hh","O-p","O-h",],
["HH-p_1","HH-p_2","HH-h_1","HH-h_2","HH-pp","HH-hp","HH-h-p","HH-hh","HH-hh-p",],
["HN-h","HN-pp","HN-hp_1","HN-hp_2","HN-h-p","HN-hh","HO-h",],
["NN-pp","NN-hp","NN-hh","NO-p","NO-h","OO-h",],   #"NO-h-p",
["HHH-p","HHH-h","HHH-hp","HHH-h-p","HHH-hh",],
["HHN-hp","HHN-hh","HHO-p","HHO-h","HHO-hp","HHO-h-p_1","HHO-h-p_2","HHO-hh",],
["HNH-p","HNH-h","HNH-hp","HNH-h-p","HNH-hh",],
["HNN-hp","HNN-hh","HNO-p","HNO-h","HNO-hp","HNO-h-p","HNO-hh",],
["NHO-hp","NHO-hh","OHO-p","OHO-h_1","OHO-h_2","OHO-h-p"],
["NNN-hhh","NNO-hh_1","NNO-hh_2",], #"NNO-hp","ONO-p"
]
mols=[]
sparators=[]
isep=0
for x in mols_gs:
    mols += x
    isep+=len(x)
    sparators.append( isep )

#print( mols )
#print( sparators )

EBmax = 40.0

#fname = "/home/prokop/Desktop/CARBSIS/Clanke_1/Data_for_Figs/Energy_B3LYP/binding_energy_big_B3LYP.dat"
fname = "/home/prokop/Desktop/CARBSIS/Clanke_1/Data_for_Figs/Energy_B3LYP/binding_energy_big_B3LYP-mithun_correct.dat"

# ============ Functions

def check_file(fname):
    with open(fname,'r') as f:
        i=0
        for line in f:
            ws = line.split()
            E = float(ws[0])
            if(E>EBmax):
                print( "WARNING: line#", i, line )
            i+=1

# ============ Main

check_file(fname)

_,Es,names,pair_names = bpu.read_dat( fname,   ni=0, nf=1, iname=0 )
Emins, uique_pnames   = bpu.find_minim_energy_confs( Es, pair_names, Emax=1000, ipivot=0 )   

n= len(mols)

Es  = np.empty((n,n)); Es [:,:] = np.NaN
Es_ = np.empty((n,n)); Es_[:,:] = np.NaN



for i,iname in enumerate(mols):
    for j,jname in enumerate(mols):
        if (iname,jname) in Emins:
            Es[i,j] = Emins[ (iname,jname) ][0]
            #print( "found: ", iname, jname, Es[i,j]  )
        elif (jname,iname) in Emins:
            Es[i,j] = Emins[ (jname,iname) ][0]
            #print( "found: ", jname, iname, Es[i,j]  )
        if Es[i,j]>EBmax:
            print( "WARNING: ", iname, jname, Es[i,j]  )
        #else: 
        #    print( "found: ", iname, jname, Es[i,j]  ) 

#EBmax = np.max( np.abs( Es ) )

print( "EBmax ", EBmax, " n " ,n  )
print( "Es.isnan() ", np.isnan( Es ).any() )
print( "Es.min,max ", Es.min(), Es.max() )

Es[Es>EBmax] = np.nan

for i in range(n):
    for j in range(i+1):
        Es_[i,j] = Es[i,j]
        if (j<i):
            Es_[j,i] = -( Es[j,i]-0.5*(Es[i,i]+Es[j,j] ) )
            #if Es_[j,i]<0: Es_[j,i]=np.nan
            if Es_[j,i]<1: Es_[j,i]=np.nan


#Es_[Es>0] = np.nan


#print( pair_names )

cmap = matplotlib.cm.seismic
#cmap.set_bad('#00FF00',1.)
cmap.set_bad('#FFFFFF',1.)

sz=12
fig, ax = plt.subplots(figsize=(sz,sz))


plt.imshow( Es_, interpolation='nearest', cmap=cmap, vmin=-EBmax, vmax=EBmax, origin='lower' )

ax.set_xticks(np.arange(n))
ax.set_yticks(np.arange(n))
ax.set_xticklabels(mols, rotation=90)
ax.set_yticklabels(mols)


for sep in sparators:
    plt.axhline( sep-0.5, ls='--', c='k' )
    plt.axvline( sep-0.5, ls='--', c='k' )


for i in range(n):
    for j in range(n):
        value = Es_[i,j]
        #value  = Es_[i,j]*50.0
        if not np.isnan(value):
            text = ax.text(j,i, "%2i" %np.abs(value), ha="center", va="center", color="gray", fontsize=8, fontweight="bold", fontname="Arial")

plt.tight_layout()

plt.savefig( "Energy_B3LYP.png", bbox_inches='tight' )
plt.savefig( "Energy_B3LYP.svg", bbox_inches='tight' )
#plt.savefig( "Energy_B3LYP.pdf", bbox_inches='tight' )

plt.show()




'''
for t in dat_ECs:
    #print( t )
    i,j=int(t[0]-0.5),int(t[1]-0.5)
    Es [i,j] = -t[2]
    Es_[i,j] = max(  t[2]/ECmax, 0.0 )
    ECs[i,j]=Es[i,j]; ECs[j,i]=-Es[i,j]*0;
    if( Es[i,j]<0 ): Es[i,j]=np.NaN
    
for t in dat_EBs:
    #print( t )
    i,j=int(t[0]-0.5),int(t[1]-0.5)
    Es [j,i] = min(  t[2], 0.0 )
    Es_[j,i] = min( -t[2]/EBmax, 0.0 )
    EBs[i,j]=Es[j,i]; EBs[j,i]=-Es[j,i]*0;
    if Es[j,i]>0: Es[j,i]=np.NaN
'''