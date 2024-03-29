import sys
import os
import matplotlib.pyplot as plt
sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import plotUtils   as plu


path=sys.argv[1]

# ============ Functions

def getFileName( path ):
    lst=[]
    # that directory
    for filename in os.listdir(path):
        f = os.path.join(path, filename)
        # checking if it is a file
        if os.path.isfile(f):
            #print(f)
            lst.append(f)
    return lst

def plotCharges( path, fname ):
    mol = au.AtomicSystem(path+"/"+fname)        #;print( mol.qs )
    hbs,rbs = mol.findHBonds( Rh=2.5 )           # ;print( hbs, rbs )
    r_labesl = [ ("%3.3fA" %r) for r in rbs ] 
    q_labels = [ ("%3.3f" %q) for q in mol.qs ]  #;print( q_labels )
    #plu.plotSystem( mol, labels=q_labels, axes=(1,2) )
    plt.subplot(1,3,1); plu.plotSystem( mol, labels=q_labels, axes=(0,1) )    ;plu.plotBonds( ps=mol.apos, links=hbs, axes=(0,1), colors="g", labels=r_labesl )
    plt.subplot(1,3,2); plu.plotSystem( mol, labels=q_labels, axes=(1,2) )    ;plu.plotBonds( ps=mol.apos, links=hbs, axes=(1,2), colors="g", labels=r_labesl )
    plt.subplot(1,3,3); plu.plotSystem( mol, labels=q_labels, axes=(2,0) )    ;plu.plotBonds( ps=mol.apos, links=hbs, axes=(2,0), colors="g", labels=r_labesl )
    plt.title(fname)

# ============ Main


#plotCharges( "/home/prokop/Desktop/CARBSIS/Paolo/res/relax/charges/CHONH2-e1_vs_CHONH2-H0", "dft.xyz"  )
#plt.show()
#exit()


#fnames = getFileName( "/home/prokop/Desktop/CARBSIS/Mithun/resp-Hbond_close_small" )
#path = "/home/prokop/Desktop/CARBSIS/Mithun/resp-Hbond_close_small"
fnames = os.listdir(path)

for fname in fnames:
    try:
        #fname = "/home/prokop/Desktop/CARBSIS/Mithun/resp-Hbond_close_small/cc-opt-CHONH2-e1_vs_CHONH2-H0.xyz"
        plt.figure( figsize=(15,5) )
        plotCharges( path, fname )
        plt.savefig( path+fname+".png", bbox_inches='tight' )
    except Exception as e:
        print( "ERROR in plotting :", fname )
        print( e )
        plt.close()



#plt.show()



