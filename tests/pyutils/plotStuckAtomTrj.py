#!/python3

import numpy as np
import matplotlib.pyplot as plt

dat = np.loadtxt('StuckAtomTrj.log' )

# atoms = set( dat[:,0] )
# for a in atoms:
#     print( "atom# ", a )
#     mask = dat[:,0] == a
#     dati = dat[mask]
#     print( dati )
#     x0 = dati[0,1]
#     y0 = dati[0,2]
#     plt.plot( dati[1]-x0,dati[2]-y0, ".-", label='atom#%i x,y' %a  )
#     #plt.plot( dat[:,1],dat[:,2], ".", label='x,y')

plt.figure( figsize=(15,5) )
x0=dat[0,2]
y0=dat[0,3]
z0=dat[0,4]
plt.subplot(1,3,1); plt.plot( dat[:,2]-x0, dat[:,3]-y0, ".-", label='x,y' ); plt.legend()
plt.subplot(1,3,2); plt.plot( dat[:,2]-x0, dat[:,4]-z0, ".-", label='x,z' ); plt.legend()
plt.subplot(1,3,3); plt.plot( dat[:,3]-y0, dat[:,4]-z0, ".-", label='y,z' ); plt.legend()


# plt.subplot(3,1,1); plt.plot( dat[:,-1], "-", label='|F|' ); plt.legend()
# plt.subplot(3,1,2); plt.plot( dat[:,-2], "-", label='|v|' ); plt.legend()
# plt.subplot(3,1,3); plt.plot( dat[:,-3], "-", label='cos(v,f)' ); plt.legend()

plt.grid()
plt.show()