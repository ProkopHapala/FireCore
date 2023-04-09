#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt

dat = np.loadtxt( 'GPU_makeGridFF.log',        skiprows=0 ).transpose()
ref = np.loadtxt( '../tMolGUIapp/initGridFF_iz_ix0_iy0.log', skiprows=1 ).transpose()



# ---- Data
#plt.figure(figsize=(15,5))

#plt.subplot(1,3,1); 

plt.title('GridFF vs NBFF'); 
plt.plot(dat[1],dat[2] ,'-r',label='E_Paul');
plt.plot(ref[1],ref[2] ,':r',label='E_Paul_ref'); 

plt.plot(dat[1],dat[4] ,'-b',label='E_Lond');
plt.plot(ref[1],ref[4] ,':b',label='E_Long_ref'); 

plt.plot(dat[1],dat[4] ,'-g',label='E_Coul');
plt.plot(ref[1],ref[4] ,':c',label='E_Coul_ref'); 

#plt.plot(x,Fz        ,label='Fz'); 
#plt.plot(x,Eref,':'  ,label='Eref');  
#plt.plot(x,FzRef,':' ,label='FzRef'); 
#plt.ylim(-0.2,0.2); 
plt.legend(); plt.grid()
plt.show()