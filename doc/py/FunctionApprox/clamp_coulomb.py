import numpy as np
import matplotlib.pyplot as plt
from func_utils import plot_with_deriv

# Coulomb potential function
def getCoulomb(r, Q):
    Coulomb_const_eVA = 1.0  # Replace with the actual constant if needed
    Q *= Coulomb_const_eVA
    E = Q / r
    F = Q / r**2
    return E, F

# Soft clamp function using NumPy masks
def soft_clamp(y, dy, y1, y2):
    y_     = y.copy()
    dy_    = dy.copy()
    mask   = y_ > y1
    y12    = y2 - y1
    invdy  = 1.0 / y12
    z = (y_[mask] - y1) * invdy
    y_ [mask] = y1 + y12 * (1 - 1 / (1 + z))
    dy_[mask] *= 1.0 / (1.0 + z)**2
    return y_, dy_

# Soft clamp function for negative values using NumPy masks
def soft_clamp_neg(y, dy, y1, y2):
    y_    = y.copy()
    dy_   = dy.copy()
    mask  = y_ < y1
    y12   = y2 - y1
    invdy = 1.0 / y12
    z = (y_[mask] - y1) * invdy
    y_  [mask] = y1 + y12 * (1 - 1 / (1 + z))
    dy_ [mask] *= 1.0 / (1.0 + z)**2
    return y_, dy_

# Smooth clamp function using NumPy masks
def smooth_clamp(y, dy, y1, y2):
    y_   = y.copy()
    dy_  = dy.copy()
    mask = y_ > y1
    y21  = y2 - y1
    z = (y_[mask] - y1) / y21
    denom = 1 / (1 + z + z**2)
    y_  [mask] = y2 - y21 * denom
    dy_ [mask] = dy[mask] * (1 + 2 * z) * denom**2
    return y_, dy_

# Smooth clamp function for negative values using NumPy masks
def smooth_clamp_neg(y, dy, y1, y2):
    y_   = y.copy()
    dy_  = dy.copy()
    mask = y_ < y1
    y21  = y2 - y1
    z = (y_[mask] - y1) / y21
    denom = 1 / (1 + z + z**2)
    y_ [mask] = y2 - y21 * denom
    dy_[mask] = dy[mask] * (1 + 2 * z) * denom**2
    return y_, dy_

# Example usage:
# y(x) = 1 / (x - 2)^6, derivative dy/dx = -6 / (x - 2)^7

# Define the domain, avoiding x = 2 to prevent division by zero
x = np.linspace(0.0001, 4.0, 1000)  # Start slightly above 2
#y =  1.0  / (x - 2.0)**6
#dy = -6.0 / (x - 2.0)**7

y,dy = getCoulomb(x, 1.0 )

y1 = 1.0
y2 = 2.0

ylim=3.0
dylim=8.0

y_clamped, dy_clamped = soft_clamp(y, dy, y1, y2)
y_clamped_neg, dy_clamped_neg = soft_clamp_neg(-y, -dy, -y1, -y2)

y_smooth_clamped, dy_smooth_clamped         = smooth_clamp(y, dy, y1, y2)
y_smooth_clamped_neg, dy_smooth_clamped_neg = smooth_clamp_neg(-y, -dy, -y1, -y2)

# Plotting the results
plt.figure(figsize=(18, 12))

# Plot y and its clamped version
plt.subplot(2, 2, 1)
plt.plot(x, y, label=r"$y(x)$", color='blue', alpha=0.5)
plt.plot(x, y_clamped, label="Soft Clamped y", color='red')
plt.plot(x, y_smooth_clamped, label="Smooth Clamped y", color='green')
plt.axhline(y1, ls='--', c='k', label="y1 and y2")
plt.axhline(y2, ls='--', c='k')
plt.title("Soft and Smooth Clamping Functions")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.grid(True)
plt.ylim(0.0, ylim)

# Plot dy and its clamped derivative
plt.subplot(2, 2, 2)
plt.plot(x, dy, label=r"$dy/dx$", color='blue', alpha=0.5)
plt.plot(x, dy_clamped, label="Soft Clamped dy/dx", color='red')
plt.plot(x, dy_smooth_clamped, label="Smooth Clamped dy/dx", color='green')
plt.axhline(0, ls='--', c='k')  # Reference line for zero derivative
plt.title("Derivative of Soft and Smooth Clamped Functions")
plt.xlabel("x")
plt.ylabel("dy/dx")
plt.legend()
plt.grid(True)
plt.ylim(-dylim, dylim)

# Plot smooth clamp for negative values
plt.subplot(2, 2, 3)
plt.plot(x, -y, label=r"$y(x)$", color='blue', alpha=0.5)
plt.plot(x, y_clamped_neg, label="Soft Clamon 'red')
plt.plot(x, y_smooth_clamped_neg, label="Smooth Clamped (Negative) y", color='purple')
plt.axhline(y1, ls='--', c='k', label="y1 and y2")
plt.axhline(y2, ls='--', c='k')
plt.title("Smooth Clamping Function (Negative Values)")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.grid(True)
plt.ylim(-ylim, 0.0)

# Plot derivative of smooth clamp for negative values
plt.subplot(2, 2, 4)
plt.plot(x, -dy, label=r"$dy/dx$", color='blue', alpha=0.5)
plt.plot(x, dy_clamped_neg, label="Soft Clamped y", color='red')
plt.plot(x, dy_smooth_clamped_neg, label="Smooth Clamped (Negative) dy/dx", color='purple')
plt.axhline(0, ls='--', c='k')  # Reference line for zero derivative
plt.title("Derivative of Smooth Clamped Function (Negative Values)")
plt.xlabel("x")
plt.ylabel("dy/dx")
plt.legend()
plt.grid(True)
plt.ylim(-dylim, dylim)

plt.tight_layout()
plt.show()