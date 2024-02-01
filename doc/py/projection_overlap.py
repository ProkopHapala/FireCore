import numpy             as np
import matplotlib.pyplot as plt

# ========== Functions

def slater( x, b  ):
    return np.exp( -b* np.abs(x) )

x = np.linspace(-4,4,1000)

y1 = slater(x-1.0,1)
y2 = slater(x+1.0,1) 

ytot = y1 + y2

rho1 = y1**2
rho2 = y2**2

rho_tot = ytot**2

rho12 = rho_tot - rho1 - rho2

plt.figure(figsize=(6,12))

plt.subplot(2,1,1)
plt.title("Wavefunction")
plt.plot(x,y1,label='y1')
plt.plot(x,y2,label='y2')
plt.plot(x,ytot,label='ytot')
plt.grid()

plt.subplot(2,1,2)
plt.title("Density")
plt.plot(x,rho1,     label='rho1')
plt.plot(x,rho2,     label='rho2')
plt.plot(x,rho_tot,  label='rho_tot')
plt.plot(x,rho12,    label='rho12')
plt.grid()

plt.show()

