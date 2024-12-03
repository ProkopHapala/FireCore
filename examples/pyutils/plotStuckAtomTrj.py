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

lw=0.5
ms=1.5


mask = (dat[:,1] == 0 )

plt.figure( figsize=(15,15) )
x0=dat[0,2]
y0=dat[0,3]
z0=dat[0,4]
# plt.subplot(3,3,1); plt.plot( dat[:,2]-x0, dat[:,3]-y0, ".-",lw=lw,ms=ms, label='x,y' ); plt.legend()
# plt.subplot(3,3,2); plt.plot( dat[:,2]-x0, dat[:,4]-z0, ".-",lw=lw,ms=ms, label='x,z' ); plt.legend()
# plt.subplot(3,3,3); plt.plot( dat[:,3]-y0, dat[:,4]-z0, ".-",lw=lw,ms=ms, label='y,z' ); plt.legend()

# plt.subplot(3,3,4); plt.plot( dat[:,2+3], dat[:,3+3], ".-",lw=lw,ms=ms, label='v x,y' ); plt.legend()
# plt.subplot(3,3,5); plt.plot( dat[:,2+3], dat[:,4+3], ".-",lw=lw,ms=ms, label='v x,z' ); plt.legend()
# plt.subplot(3,3,6); plt.plot( dat[:,3+3], dat[:,4+3], ".-",lw=lw,ms=ms, label='v y,z' ); plt.legend()

# plt.subplot(3,3,7); plt.plot( dat[:,2+6], dat[:,3+6], ".-",lw=lw,ms=ms, label='F x,y' ); plt.legend()
# plt.subplot(3,3,8); plt.plot( dat[:,2+6], dat[:,4+6], ".-",lw=lw,ms=ms, label='F x,z' ); plt.legend()
# plt.subplot(3,3,9); plt.plot( dat[:,3+6], dat[:,4+6], ".-",lw=lw,ms=ms, label='F y,z' ); plt.legend()

plt.subplot(3,3,1); plt.plot( dat[mask,2]-x0, dat[mask,3]-y0, ".-",lw=lw,ms=ms, label='x,y' ); plt.legend()
plt.subplot(3,3,2); plt.plot( dat[mask,2]-x0, dat[mask,4]-z0, ".-",lw=lw,ms=ms, label='x,z' ); plt.legend()
plt.subplot(3,3,3); plt.plot( dat[mask,3]-y0, dat[mask,4]-z0, ".-",lw=lw,ms=ms, label='y,z' ); plt.legend()

plt.subplot(3,3,4); plt.plot( dat[mask,2+3], dat[mask,3+3], ".-",lw=lw,ms=ms, label='v x,y' ); plt.legend()
plt.subplot(3,3,5); plt.plot( dat[mask,2+3], dat[mask,4+3], ".-",lw=lw,ms=ms, label='v x,z' ); plt.legend()
plt.subplot(3,3,6); plt.plot( dat[mask,3+3], dat[mask,4+3], ".-",lw=lw,ms=ms, label='v y,z' ); plt.legend()

plt.subplot(3,3,7); plt.plot( dat[mask,2+6], dat[mask,3+6], ".-",lw=lw,ms=ms, label='F x,y' ); plt.legend()
plt.subplot(3,3,8); plt.plot( dat[mask,2+6], dat[mask,4+6], ".-",lw=lw,ms=ms, label='F x,z' ); plt.legend()
plt.subplot(3,3,9); plt.plot( dat[mask,3+6], dat[mask,4+6], ".-",lw=lw,ms=ms, label='F y,z' ); plt.legend()


# plt.subplot(3,3,1); plt.plot( dat[:,2]-x0, dat[:,3]-y0, ".",lw=lw,ms=ms, label='x,y' ); plt.legend()
# plt.subplot(3,3,2); plt.plot( dat[:,2]-x0, dat[:,4]-z0, ".",lw=lw,ms=ms, label='x,z' ); plt.legend()
# plt.subplot(3,3,3); plt.plot( dat[:,3]-y0, dat[:,4]-z0, ".",lw=lw,ms=ms, label='y,z' ); plt.legend()

# plt.subplot(3,3,4); plt.plot( dat[:,2+3], dat[:,3+3], ".",lw=lw,ms=ms, label='v x,y' ); plt.legend()
# plt.subplot(3,3,5); plt.plot( dat[:,2+3], dat[:,4+3], ".",lw=lw,ms=ms, label='v x,z' ); plt.legend()
# plt.subplot(3,3,6); plt.plot( dat[:,3+3], dat[:,4+3], ".",lw=lw,ms=ms, label='v y,z' ); plt.legend()

# plt.subplot(3,3,7); plt.plot( dat[:,2+6], dat[:,3+6], ".",lw=lw,ms=ms, label='F x,y' ); plt.legend()
# plt.subplot(3,3,8); plt.plot( dat[:,2+6], dat[:,4+6], ".",lw=lw,ms=ms, label='F x,z' ); plt.legend()
# plt.subplot(3,3,9); plt.plot( dat[:,3+6], dat[:,4+6], ".",lw=lw,ms=ms, label='F y,z' ); plt.legend()


plt.savefig('StuckAtomTrj_AtomTrj.png')

plt.figure( figsize=(15,5) )
plt.subplot(3,1,1); plt.plot( dat[:,-1], "-",lw=lw,ms=ms, label='|F|' ); plt.legend()
plt.subplot(3,1,2); plt.plot( dat[:,-2], "-",lw=lw,ms=ms, label='|v|' ); plt.legend()
plt.subplot(3,1,3); plt.plot( dat[:,-3], "-",lw=lw,ms=ms, label='cos(v,f)' ); plt.legend()

plt.savefig('StuckAtomTrj_Conv.png')

#plt.grid()
plt.show()