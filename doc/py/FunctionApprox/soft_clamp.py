import numpy as np
import matplotlib.pyplot as plt

def soft_clamp(y, dy, y1, y2):
    """
    Applies a soft clamp to y, smoothly transitioning values above y1 towards y2.
    Also computes the derivative dy accordingly using the chain rule.

    Parameters:
    - y: np.ndarray, input values to clamp
    - dy: np.ndarray, derivatives of y with respect to some variable x
    - y1: float, lower threshold for clamping
    - y2: float, upper threshold for clamping

    Returns:
    - y_new: np.ndarray, clamped y values
    - dy_new: np.ndarray, updated derivatives
    """
    y_new  = y.copy()
    dy_new = dy.copy()
    mask   = y_new > y1
    y12    = y2 - y1
    invdy  = 1.0 / y12
    z = (y_new[mask] - y1) * invdy
    y_new[mask]   = y1 + y12 * (1 - 1 / (1 + z))
    dy_new[mask] *= 1.0 / (1.0 + z)**2
    return y_new, dy_new

def soft_clamp_exp(array_in, min_value, max_value, width=10.0):
    """
    Soft clamp function using NumPy masks, smoothly clamping values.

    Parameters:
    - array_in: np.ndarray, input array
    - min_value: float, minimum clamp value
    - max_value: float, maximum clamp value
    - width: float, width of the transition region

    Returns:
    - array_out: np.ndarray, output array after applying the soft clamp
    """
    array_out = array_in.copy()

    # Mask for values greater than max_value
    if max_value < np.inf:
        mask_max = array_out > max_value
        array_out[mask_max] = (
            (array_out[mask_max] - max_value) /
            (1 + np.exp((array_out[mask_max] - max_value) / width)) + max_value
        )

    # Mask for values less than min_value
    if min_value > -np.inf:
        mask_min = array_out < min_value
        array_out[mask_min] = (
            (array_out[mask_min] - min_value) /
            (1 + np.exp((array_out[mask_min] - min_value) / -width)) + min_value
        )

    return array_out

def plot_with_deriv(ax1, ax2, x, y, dy, label, color, alpha=1.0):
    ax1.plot(x, y, label=label, color=color, alpha=alpha)
    if dy is not None:
        ax2.plot(x, dy, label=label, color=color, alpha=alpha)

# Example usage:
# y(x) = 1 / (x - 2)^6, derivative dy/dx = -6 / (x - 2)^7

# Define the domain, avoiding x = 2 to prevent division by zero
x  = np.linspace(0.0, 4.0, 1000)  # Start slightly above 2
y  = 1. / (x - 2.)**6
dy = -6. / (x - 2.)**7

y1 = 1.0
y2 = 3.0

# Apply the soft clamp
y_clamped, dy_clamped = soft_clamp(y, dy, y1, y2 )

y_clamped_exp = soft_clamp_exp(y, -y1, y1, width=(y2-y1) )

# Plotting the results
plt.figure(figsize=(12, 6))

# Create subplots
ax1 = plt.subplot(1, 2, 1)
ax2 = plt.subplot(1, 2, 2)

# Plot functions and derivatives
plot_with_deriv(ax1, ax2, x, y, dy,
               label=r"$y(x) = \frac{1}{(x-2)^6}$", color='blue', alpha=0.5)
plot_with_deriv(ax1, ax2, x, y_clamped, dy_clamped,
               label="Clamped y", color='red')
plot_with_deriv(ax1, ax2, x, y_clamped_exp, None,
               label="Clamped_exp y", color='g')

# Add reference lines and titles
ax1.axhline(y1, ls='--', c='k', label="y1 and y2")
ax1.axhline(y2, ls='--', c='k')
ax1.set_title("Soft Clamping Function")
ax1.set_xlabel("x")
ax1.set_ylabel("y")
ax1.legend()
ax1.grid(True)
ax1.set_ylim(0.0, 5.0)

ax2.axhline(0, ls='--', c='k')
ax2.set_title("Derivative of Soft Clamped Function")
ax2.set_xlabel("x")
ax2.set_ylabel("dy/dx")
ax2.legend()
ax2.grid(True)
ax2.set_ylim(-15.0, 15.0)

plt.tight_layout()
plt.show()