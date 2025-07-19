import numpy as np
import matplotlib.pyplot as plt

def numDeriv(x, y):
    """Numerical derivative using central difference"""
    dx = x[2:]-x[:-2]
    dy = y[2:]-y[:-2]
    x_ = x[1:-1]
    return -dy/dx, x_

def plot_with_deriv(x, y, deriv=None, label='Function', color='blue', bNumDeriv=True, linestyle='-', linewidth=2, ax1=None, ax2=None ):
    """
    Plots a function and its derivative on given axes
    
    Args:
        ax_func: axis for function plot
        ax_deriv: axis for derivative plot
        x: x values
        y: y values (function)
        deriv: analytical derivative values (optional)
        label: legend label
        color: plot color
        show_num_deriv: whether to show numerical derivative
        linestyle: line style
        linewidth: line width
    """
    if ax1 is None:
        fig,(ax1,ax2) = plt.subplots(2,1)
    ax1.plot(x, y, label=label, color=color,  linestyle=linestyle, linewidth=linewidth)
    if deriv is not None:
        ax2.plot(x, deriv, label=f'{label} (analytical)',  color=color, linestyle=linestyle, linewidth=linewidth)
        if bNumDeriv:
            num_deriv, num_x = numDeriv(x, y)
            ax2.plot(num_x, num_deriv, label=f'{label} (numerical)', color=color, linestyle=':', linewidth=1.5, alpha=0.7)
    return fig, (ax1, ax2)

# Example usage:
"""
Basic usage example for plot_with_deriv:

# Create figure
fig, (ax1, ax2) = plt.subplots(2, 1)

# Plot function and derivative
plot_with_deriv(ax1, ax2, x, y, dy, label='Function 1', color='blue')
plot_with_deriv(ax1, ax2, x, y2, dy2, label='Function 2', color='red', linestyle='--')

# Customize axes
ax1.set_title('Function Comparison')
ax2.set_title('Derivative Comparison')
ax1.legend()
ax2.legend()
plt.show()

Advanced usage with automatic figure creation:

# This will create a new figure automatically
fig, (ax1, ax2) = plot_with_deriv(x, y, dy, label='Function', color='green')

# Add reference lines
ax1.axhline(0, ls='--', c='k')
ax2.axhline(0, ls='--', c='k')
plt.show()
"""