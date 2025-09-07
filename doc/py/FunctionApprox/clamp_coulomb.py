import numpy as np
import matplotlib.pyplot as plt
from func_utils import plot1d, plot_with_deriv, plot1d_zip

# Coulomb potential function
def getCoulomb(r, Q):
    Coulomb_const_eVA = 1.0  # Replace with the actual constant if needed
    Q *= Coulomb_const_eVA
    E = Q / r
    F = Q / r**2
    return E, F

# Soft clamp function using NumPy masks
def soft_clamp( x, y1, y2):
    y,dy = getCoulomb(x, 1.0)
    mask   = y > y1
    y12    = y2 - y1
    invdy  = 1.0 / y12
    z = (y[mask] - y1) * invdy
    y[mask] = y1 + y12 * (1 - 1 / (1 + z))
    dy[mask] *= 1.0 / (1.0 + z)**2
    return y, dy

# Smooth clamp function using NumPy masks
def smooth_clamp(x, y1, y2):
    y,dy = getCoulomb(x, 1.0)
    mask = y > y1
    y21  = y2 - y1
    z = (y[mask] - y1) / y21
    denom = 1 / (1 + z + z**2)
    y  [mask] = y2 - y21 * denom
    dy [mask] = dy[mask] * (1 + 2 * z) * denom**2
    return y, dy

# Soft clamp function for negative values using NumPy masks
def soft_clamp_neg(x, y1, y2):
    y,dy = getCoulomb(x, -1.0)
    mask  = y < y1
    y12   = y2 - y1
    invdy = 1.0 / y12
    z = (y[mask] - y1) * invdy
    y[mask] = y1 + y12 * (1 - 1 / (1 + z))
    dy[mask] *= 1.0 / (1.0 + z)**2
    return y, dy


# Smooth clamp function for negative values using NumPy masks
def smooth_clamp_neg(x, y1, y2):
    y,dy = getCoulomb(x, -1.0)
    mask = y < y1
    y21  = y2 - y1
    z = (y[mask] - y1) / y21
    denom = 1 / (1 + z + z**2)
    y [mask] = y2 - y21 * denom
    dy[mask] = dy[mask] * (1 + 2 * z) * denom**2
    return y, dy

# Example usage
if __name__ == "__main__":
    x      = np.linspace(0.0001, 4.0, 1000)
    y, dy  = getCoulomb(x, 1.0)
    y1, y2 = 1.0, 2.0
    ylim, dylim = 3.0, 4.0
            
    fig, (ax1, ax2) = plot1d_zip(x, [
        ( "Coulomb"      , getCoulomb  (x, 1.0   )  ),
        ( "soft_clamp"   , soft_clamp  (x, y1, y2)  ),
        ( "smooth_clamp" , smooth_clamp(x, y1, y2)  ),
    ])
    ax1.set_ylim(0.0, ylim  )
    ax2.set_ylim(0.0, dylim )

    fig, (ax1, ax2) = plot1d_zip(x, [
        ( "Coulomb"          , getCoulomb       ( x,  -1.0      ) ),
        ( "soft_clamp_neg"   , soft_clamp_neg    (x, -y1, -y2) ),
        ( "smooth_clamp_neg" , smooth_clamp_neg  (x, -y1, -y2) ),
    ])
    ax1.set_ylim(-ylim, 0.0)
    ax2.set_ylim(-dylim, 0.0)

    plt.tight_layout()
    plt.show()