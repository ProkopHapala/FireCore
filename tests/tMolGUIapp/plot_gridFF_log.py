#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt( 'gridFF_vs_NBFF.log', skiprows=1 ).transpose()

x     = data[3]
E     = data[4]
Eref  = data[5]
Fz    = data[10]
FzRef = data[11]

# ---- Data
plt.figure(figsize=(15,5))

plt.subplot(1,3,1); plt.title('GridFF vs NBFF'); plt.plot(x,E     ,label='E'); plt.plot(x,Fz,label='Fz'); plt.plot(x,Eref,':', label='Eref');  plt.plot(x,FzRef,':',label='FzRef'); plt.grid()
# ---- Errr
plt.subplot(1,3,2); plt.title('Absolute Error'); plt.plot(x,E-Eref,'b',label='dE_abs');  plt.plot(x,Fz-FzRef,label='dFz_abs'); plt.grid()
# ---- Relative error
plt.subplot(1,3,3); plt.title('Relative Error'); plt.plot(x,np.abs((E-Eref)/(np.abs(Eref)+1e-6)),label='dE_rel');       plt.plot(x,np.abs((Fz-FzRef)/(np.abs(FzRef)+1e-6)),label='dFz_rel'); plt.grid()

plt.show()