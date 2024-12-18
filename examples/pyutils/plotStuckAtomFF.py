#!/python3

import numpy as np
import matplotlib.pyplot as plt

datx = np.loadtxt('Stuck_Fx.txt', skiprows=1 )
daty = np.loadtxt('Stuck_Fy.txt', skiprows=1 )
datz = np.loadtxt('Stuck_Fz.txt', skiprows=1 )

def diff(x,y):
    x_ = x[2:] - x[:-2]
    y_ = y[2:] - y[:-2]
    return x_, y_

def plot_diff(x,y,label=None):
    x_, y_ = diff(x,y)
    plt.plot(x_,y_,label=label)

plt.plot(datx[:,1]-datx[0,1],datx[:,4],label='Fx')
plt.plot(daty[:,2]-daty[0,2],daty[:,5],label='Fy')
plt.plot(datz[:,3]-datz[0,3],datz[:,6],label='Fz')

plt.plot(datx[:,1]-datx[0,1],datx[:,5],label='Fx')
plt.plot(daty[:,2]-daty[0,2],daty[:,6],label='Fy')
plt.plot(datz[:,3]-datz[0,3],datz[:,4],label='Fz')

plt.plot(datx[:,1]-datx[0,1],datx[:,6],label='Fx')
plt.plot(daty[:,2]-daty[0,2],daty[:,4],label='Fy')
plt.plot(datz[:,3]-datz[0,3],datz[:,5],label='Fz')


# plot_diff(datx[:,1]-datx[0,1],datx[:,4],label='Fx')
# plot_diff(daty[:,2]-daty[0,2],daty[:,5],label='Fy')
# plot_diff(datz[:,3]-datz[0,3],datz[:,6],label='Fz')

plt.grid()
plt.show()