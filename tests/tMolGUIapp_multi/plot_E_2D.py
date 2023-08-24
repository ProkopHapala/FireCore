#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

dat = np.genfromtxt( 'lattice_scan_2d_multi.dat' ).transpose()


axs=dat[1]     ;print("axs", axs )
ays=dat[2]     ;print("ays", ays )
Es =dat[11]    ;print("Es" , Es  )
ixs =dat[15]   ;print("ixs", ixs )
iys =dat[17]   ;print("iys", iys )

print(Es.shape)

extent = ( axs.min(), axs.max(),  ays.min(), ays.max() )
#extent = (  ays.min(), ays.max(), axs.min(), axs.max(), )

Emap = np.reshape( Es, (40,30)  )

Emin = Emap.min()
dEmax = 5.0

plt.imshow(Emap.transpose(), vmin=Emin, vmax=Emin+dEmax, origin='lower', extent=extent, cmap='plasma')
plt.colorbar()

plt.show()





