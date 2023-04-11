import numpy             as np
import matplotlib.pyplot as plt

# ========== Functions

def KvaziFIREdamp( c, clim, damps ):
    # ----- velocity & force ~ perpendicular
    t = (c-clim[0] )/( clim[1] - clim[0]  )
    cv = damps[0] + (damps[1]-damps[0])*t
    #cf =     damps[1] *t*(1-t)*4
    cf =     damps[1]*t*(1-t)*2
    # ----- velocity & force ~ against each other
    mask_lo     =  c < clim[0]
    cv[mask_lo] = damps[0]  # v    // like 0.5 (strong damping)
    cf[mask_lo] = 0             # f
    # ----- velocity & force ~ alligned
    mask_hi     =  c > clim[1]
    cv[mask_hi] = damps[1]  # v    // like 0.99 (weak dampong damping)
    cf[mask_hi] = 0           # f
    return cv,cf

# ============ Main

phi = np.linspace(-np.pi,np.pi,1000)
s   = np.sin(phi)
c   = np.cos(phi)
s2 = 1-c*c

cv,cf = KvaziFIREdamp( c, [-0.1,+0.1], [0.1,0.9] )
cv,cf = KvaziFIREdamp( c, [-0.3,+0.3], [0.1,0.9] )

plt.title("FIRE velocity damping")
plt.plot(phi,c ,label='cos(x)')
plt.plot(phi,cv,label='cv (damp)')
plt.plot(phi,cf,label='cf (turn)')

plt.legend()
plt.grid()

plt.show()