#! /bin/python3

import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

sys.path.append("../../")
#from pyBall import Kekule as kek
#from pyBall import atomicUtils as au
from pyBall import plotUtils   as pu




'''
      Donor         Acceptor
-----------------------------
1    -OH,-NH2             
2        -NH-         -O-,=O
3                     =N-
'''


groups={
'-CH3':1,'=CH2':2,'-CH2-':2,'=CH-':3,'C':4,'C/N':(3,4),
'-NH2':1,'=NH' :2,'-NH-' :2,'=N-' :3,
'-OH' :1,'=O'  :2,'-O-'  :2,
}

BO_groups={
'-CH3':1,'=CH2':2,'-CH2-':1,
'-NH2':1,'=NH' :2,'-NH-' :1,
'-OH' :1,'=O'  :2,'-O-'  :1, 
}

non=[None]
cdon =[ '-NH-',      ]    # central donors
cacc =[ '=N-', '-O-' ]
sdon =[ '-NH2'       ]
sacc =[ '=O'         ]
Cs=['=CH-']
Cc=['C']
CN=['C/N']

central = cdon+cacc
side    = sdon+sacc 
don     = cdon+sdon
acc     = cacc+sacc 




# ------- Sites
sites={
"A1": [(-2.4, 0.0),  side+central+non   ],
"A2":[( 0.0, 0.0),       central       ], 
"A3":[( 2.4, 0.0),  side+central+non   ],
"B1":[(-3.6,-0.7), Cs ],
"B2":[(-1.2,-0.7), Cc ], 
"B3":[( 1.2,-0.7), Cc ], 
"B4":[( 3.6,-0.7), Cs ],
"C1":[(-3.6,-2.1), Cs ], 
"C2":[(-1.2,-2.1), Cc ], 
"C3":[( 1.2,-2.1), Cc ], 
"C4":[( 3.6,-2.1), Cs ],
"D1":[(-2.4,-2.8), CN+non ], 
"D2":[( 0.0,-2.8), CN+non ], 
"D3":[( 2.4,-2.8), CN+non ],
}


sites_can={
"A1":[(-2.4, 0.0), 0      ],
"A2":[( 0.0, 0.0), 1      ], 
"A3":[( 2.4, 0.0), 2      ],
"B1":[(-3.6,-0.7), "=CH-" ],
"B2":[(-1.2,-0.7), 3      ], 
"B3":[( 1.2,-0.7), 4      ], 
"B4":[( 3.6,-0.7), "=CH-" ],
"C1":[(-3.6,-2.1), "C/N"  ], 
"C2":[(-1.2,-2.1), 5      ], 
"C3":[( 1.2,-2.1), 6      ], 
"C4":[( 3.6,-2.1), "C/N"  ],
"D1":[(-2.4,-2.8), "C/N"  ], 
"D2":[( 0.0,-2.8), "C/N"  ], 
"D3":[( 2.4,-2.8), "C/N"  ],
}


# ====== FUNCTIONS

def makeNeighs( bonds, na ):
    aneighs=[ [] for a in range(na) ]
    bneighs=[ [] for a in range(na) ]
    for ib,b in enumerate(bonds):
        i,j=b
        aneighs[i].append(j)
        aneighs[j].append(i)
        bneighs[i].append(ib)
        bneighs[j].append(ib)
    return aneighs,bneighs

def determineBondsOfAtoms( atypes, bneighs, BO, ibonds ):
    #print( "BO    ", BO )
    #print( "ibonds", ibonds )
    n_new = 0
    for ia,a2b in enumerate(bneighs):
        typ  = atypes[ia]
        bo   = 0
        n_un = 0
        i_un = -1
        i_uno = -1 
        for ib in a2b:
            boi = BO[ib]
            if( boi > 0 ):  # determined bond order
                bo+=boi
            else:           # undetermined bond order
                i_uno = i_un
                i_un  = ib
                n_un += 1
        ao = groups[typ]
        if isinstance(ao,int):
            #print( "atom[%i]%s nng=%i n_un=%i ao=%i bo=%i " %(ia, typ, len(a2b),n_un, ao, bo), [BO[j] for j in a2b], [ibonds[j] for j in a2b ]  )
            if bo>ao: 
                print( "ERROR sum of bond-orders for atom[%i]%s is %i but atom can have max %i bonds" %(ia,typ,bo,ao)  )
                #exit(0)
                return -1
            dbo = ao - bo
            if(n_un == 1):
                if bo>ao: 
                    print( "ERROR bond[%i] of atom[%i] determined bo(%i)<1 | bo(%i) ao(%i)" %(i_un,ia,dbo, bo, ao )  )
                    return -1
                BO[i_un] = dbo
                n_new+=1
            elif(n_un==2):
                if( dbo<2 ):
                    print( "ERROR  atom[%i]%s ao(%i) insufficient for bo(%i)+n_un(%i)" %(ia,typ,ao,bo,n_un)  )
                elif( dbo==2 ):
                    BO[i_un ] = 1
                    BO[i_uno] = 1
                    n_new+=1
        else:
            ao_min,ao_max = ao
            ao_req = bo+n_un
            #print( "atom[%i]%s: bo(%i)+n_un(%i)=%i "  %(ia,typ,bo,n_un,ao_req ) )
            if ao_req > ao_max:
                print( "ERROR atom[%i]%s: bo(%i)+n_un(%i)=%i > ao_max(%i) " %(ia,typ, bo,n_un,ao_req, ao_max)  )
            elif ao_req==3:
                #print( "replaced atom[%i] %s -> %s " %(ia,typ,"=CH-") )
                atypes[ia]="=CH-"
                n_new+=1
            elif (bo==2) and ( n_un==0 ) and (typ=="C/N"):
                #print( "replaced atom[%i] %s -> %s " %(ia,typ,"-NH-") )
                atypes[ia]="-NH-"
                n_new+=1
    
    return n_new

def determineBondOrder( atoms, bonds, BO ): 
    for i,b in enumerate(bonds):
        A = atoms[b[0]]
        B = atoms[b[1]]
        boA = BO_groups.get( A, -1 )
        boB = BO_groups.get( B, -1 )
        bo  = None
        bomin = min(boA,boB)
        bomax = max(boA,boB)
        #if(bomax > 1 ): print( bomax, A,B, boA,boB  )
        if ( bomin > 0 ):  # bond determined from both sides
            if boA!=boB: 
                print( "ERROR bond[%i](%s,%s) cannot be bonded ", i, A,B ); 
                return
            bo = boA
        elif bomax>-1:  # bond determined from one side
            bo=bomax
        if bo is not None:
            BO[i]=bo
        #print( i, bo, A,B, boA, boB )
    return BO


def settleBonds( bonds, atypes, BO, bneighs ):
    for ib,b in enumerate(bonds):
        bo = BO[ib]
        if bo == -1:
            i,j=b
            A = atypes[i]
            B = atypes[j]
            if( (A=="C/N") and (B=="C/N") ):

                ngA = bneighs[i]
                boA=0
                for ja in ngA:
                    if(ja!=ib):
                        boiA = BO[i]
                        if(boiA>0):
                            boA+=boiA

                ngB = bneighs[j]
                boB=0
                for jb in ngB:
                    if(jb!=ib):
                        boiB = BO[i]
                        if(boiB>0):
                            boB+=boiB

                if( (boA==1) and (boB==1) ): 
                    atypes[i]="=CH-"
                    atypes[j]="=CH-"
                    BO   [ib]=2



def makeBondOrder(mol):
    #atypes, _, ibonds = molecules[0]
    atypes, apos, ibonds, BO = mol
    #BO = [-1 for i in ibonds]
    BO = determineBondOrder( atypes, ibonds, BO ) 
    #print( BO )
    
    aneighs,bneighs = makeNeighs( ibonds, len(atypes) )
    #print(ibonds)
    #print(bneighs)
    #print(BO)
    ret = 1
    for itr in range(10):
        #print("#========= iter ", itr )
        ret = determineBondsOfAtoms( atypes, bneighs, BO, ibonds )
        if ret <=0: break
    
    mol = ( atypes, apos, ibonds, BO )
    if ret<0:
        plt.figure( figsize=(5,5))
        plotMolecule( mol )
        plt.show()
        plt.close()

    settleBonds( ibonds, atypes, BO, bneighs )

    return mol

def plotMolecule( mol ):
    atypes, apos, ibonds, BO = mol

    choice=( '#808080', '#8080f0' )
    colors = [ choice[bo==-1] for bo in BO ]

    BO = np.array(BO,dtype=np.float); BO[BO<1.0]=1.5      #print( BO )
    lws = (BO-0.8)*8
    
    pu.plotBonds( links=ibonds, ps=apos, lws=lws, axes=(0,1), colors=colors )
    pu.plotAtoms( apos, labels=atypes, axes=(0,1) )
    #plt.axis('equal')
    plt.gca().set_aspect('equal',adjustable='box')
    plt.ylim(-3.5,0.5)
    plt.xlim(-4.0,4.0)

def plotMolecules( molecules, sz=20 ):
    plt.figure( figsize=(sz,sz) )
    iimg = 0
    nmol = len(molecules)
    for i,mol in enumerate(molecules):
        print( "#### Molecule[%i] " %i, mol[0] )
        plt.subplot(6,4,iimg+1)
        plotMolecule( mol )
        plt.title("Mol[%i]" %i )
        iimg+=1
        if( (iimg >= 24) or (i>=nmol-1)):
            plt.savefig( "endgroups_%03i_%03i.png" %(i-24+1,i+1) , bbox_inches='tight')  
            iimg=0  
            plt.figure( figsize=(sz,sz) )

# ====== Generate skeletons (no-types assignment)

atoms0 = ["A2","B2","B3","C2","C3"]
bonds0 = [ ("A2","B2"),("A2","B3"), ("B2","C2"), ("B3","C3") ]

nmol=0
skeletons = [] 
for h1 in [0,1,5,6]:
    atoms1 = []
    bonds1 = []
    if   h1==0: # missing
        pass
    elif h1==1: # terminating
        atoms1.append( "A1" )
        bonds1.append( ("A1","B2") )
    elif h1==5: # pentagon
        atoms1 += [ "A1","B1","C1" ]
        bonds1 += [ ("A1","B2"),("A1","B1"),("B1","C1"),("C1","C2"), ]
    elif h1==6: # hexagon
        atoms1 += [ "A1","B1","C1","D1" ]
        bonds1 += [ ("A1","B2"),("A1","B1"),("B1","C1"),("C1","D1"),("D1","C2"), ]
            
    for h3 in [0,1,5,6]:

        if(h3>h1): continue    # prevent symmetric variants

        atoms3 = []
        bonds3 = []
        if   h3==0: # missing
            pass
        elif h3==1: # end
            atoms3.append( "A3" )
            bonds3.append( ("A3","B3") )
        elif h3==5: # pentagon
            atoms3 += [ "A3","B4","C4" ]
            bonds3 += [ ("A3","B3"),("A3","B4"),("B4","C4"),("C4","C3"), ]
        elif h3==6: # hexagon
            atoms3 += [ "A3","B4","C4","D3" ]
            bonds3 += [ ("A3","B3"),("A3","B4"),("B4","C4"),("C4","D3"),("D3","C3"), ]

        for h2 in [5,6]:
            atoms2 = []
            bonds2 = []
            if   h2==5:  # pentagon
                bonds2 += [ ("C2","C3") ]
            elif h2==6: # hexagon
                atoms2 += [ "D2" ]
                bonds2 += [ ("D2","C2"),("D2","C3") ]

            nmol+=1
            print( nmol, h1,h3,h2 )
            
            molecule = ( atoms0+atoms1+atoms2+atoms3,  bonds0+bonds1+bonds2+bonds3, (h1,h2,h3),  )
            skeletons.append( molecule )

print( len(skeletons) )

def assing_type( val, key, aaa ):
    if isinstance(val, int):
        return aaa[val]
    else:
        return val

molecules = []
for i,mol in enumerate(skeletons):
    atoms,bonds,variant = mol
    v1,v2,v3=variant
    dct={ k:i for i,k in enumerate(atoms) }
    ibonds = [ (dct[b[0]],dct[b[1]]) for b in bonds ]
    apos   = np.array([ sites[k][0] for k in atoms ])

    b2="=CH-"
    b3="=CH-"
    c2="C/N"
    c3="C/N"
    if   v1==0:  # missing 
        groups1=[None] 
    elif v1==1:  # end
        groups1=side
        b2="C"
    else:  # cycle
        groups1=central
        b2="C"
        c2="C"

    if   v3==0:  # missing 
        groups3=[None] 
    elif v3==1:  # end
        groups3=side
        b3="C"
    else:  # cycle
        groups3=central
        b3="C"
        c3="C"

    group2 = central

    for a1 in groups1:
        for a3 in groups3:
            for a2 in group2:
                aaa = (a1,a2,a3,b2,b3,c2,c3)
                atypes= [ assing_type( sites_can[k][1], k, aaa ) for k in atoms ] 

                BO = [-1 for i in ibonds]
                mol = ( atypes, apos, ibonds, BO )
                molecules.append(mol)


print(len(molecules))

'''
sz=20
plt.figure( figsize=(sz,sz) )
iimg = 0
for i,mol in enumerate(molecules):
    atypes, apos, ibonds = mol
    print( "#### Molecule[%i] " %i, atypes )

    plt.subplot(5,5,iimg+1)
    pu.plotAtoms( apos, labels=atypes, axes=(0,1) )
    pu.plotBonds( links=ibonds, ps=apos, axes=(0,1) )
    plt.axis('equal')
    plt.xlim(-4.,4.)
    plt.ylim(-7.,1.)
    plt.title("Mol[%i]" %i )
    iimg+=1
    if(iimg >= 25):
        plt.savefig( "endgroups_%03i_%03i.png" %(i-25+1,i+1) , bbox_inches='tight')  
        iimg=0  
        plt.figure( figsize=(sz,sz) )
'''

'''
plt.figure( figsize=(10,10) )
for i,mol in enumerate(skeletons):
    atoms,bonds,variant = mol
    print( "#### Molecule[%i] " %i, variant, atoms )
    #print( "atoms: ", atoms )
    #print( "bonds"  , bonds )

    plt.subplot(5,5,i+1)
    
    na   = len(atoms)
    apos = np.array([ sites[k][0] for k in atoms ])
    es   = ['C']*na
    dct={ k:i for i,k in enumerate(atoms) }
    #print( "dct: ", dct )

    ibonds = [ (dct[b[0]],dct[b[1]]) for b in bonds ]
    
    #print("ibonds", ibonds)
    pu.plotAtoms( apos, es=es, axes=(0,1) )
    pu.plotBonds( links=ibonds, ps=apos, axes=(0,1) )
    plt.axis('equal')
    plt.xlim(-4.,4.)
    plt.ylim(-7.,1.)
plt.show()
'''

molecules_new = []
for i,mol in enumerate(molecules):
    print( "########## mol", i )
    #if(i>10):break
    mol_new = makeBondOrder(mol)
    molecules_new.append(mol_new) 

plotMolecules( molecules_new, sz=20 )
#print( BO )







